# Targets, drugs and association evidence

This page documents **62 newly added public functions**. [Complete index](index.md) · [Conventions and examples](usage.md)

Each function lists its current Python signature, parameters, returns and registered execution requirements. Provider descriptions retain source wording. In generated database wrappers, None normally omits a parameter; actual required fields and conditional constraints come from the input schema. Follow each operation's coordinate and assembly conventions.

## Function index

- [`dnakit.databases.ch_embl_get_drug_mechanisms`](#dnakit.databases.ch_embl_get_drug_mechanisms)
- [`dnakit.databases.ch_embl_get_molecule`](#dnakit.databases.ch_embl_get_molecule)
- [`dnakit.databases.ch_embl_get_target_activities`](#dnakit.databases.ch_embl_get_target_activities)
- [`dnakit.databases.ch_embl_search_drugs`](#dnakit.databases.ch_embl_search_drugs)
- [`dnakit.databases.ch_embl_search_targets`](#dnakit.databases.ch_embl_search_targets)
- [`dnakit.databases.ctd_get_chemical_gene_interactions`](#dnakit.databases.ctd_get_chemical_gene_interactions)
- [`dnakit.databases.dgidb_get_drug_gene_interactions`](#dnakit.databases.dgidb_get_drug_gene_interactions)
- [`dnakit.databases.dgidb_get_gene_druggability`](#dnakit.databases.dgidb_get_gene_druggability)
- [`dnakit.databases.drugbank_get_drug_basic_info_by_drug_name_or_id`](#dnakit.databases.drugbank_get_drug_basic_info_by_drug_name_or_id)
- [`dnakit.databases.drugbank_get_pharmacology_by_drug_name_or_drugbank_id`](#dnakit.databases.drugbank_get_pharmacology_by_drug_name_or_drugbank_id)
- [`dnakit.databases.drugbank_get_targets_by_drug_name_or_drugbank_id`](#dnakit.databases.drugbank_get_targets_by_drug_name_or_drugbank_id)
- [`dnakit.databases.drugbank_vocab_search`](#dnakit.databases.drugbank_vocab_search)
- [`dnakit.databases.europepmc_disease_target_score`](#dnakit.databases.europepmc_disease_target_score)
- [`dnakit.databases.expression_atlas_disease_target_score`](#dnakit.databases.expression_atlas_disease_target_score)
- [`dnakit.databases.fda_get_active_ingredient_info_by_drug_name`](#dnakit.databases.fda_get_active_ingredient_info_by_drug_name)
- [`dnakit.databases.fda_get_adverse_reactions_by_drug_name`](#dnakit.databases.fda_get_adverse_reactions_by_drug_name)
- [`dnakit.databases.fda_get_boxed_warning_info_by_drug_name`](#dnakit.databases.fda_get_boxed_warning_info_by_drug_name)
- [`dnakit.databases.fda_get_clinical_studies_info_by_drug_name`](#dnakit.databases.fda_get_clinical_studies_info_by_drug_name)
- [`dnakit.databases.fda_get_indications_by_drug_name`](#dnakit.databases.fda_get_indications_by_drug_name)
- [`dnakit.databases.fda_get_mechanism_of_action_by_drug_name`](#dnakit.databases.fda_get_mechanism_of_action_by_drug_name)
- [`dnakit.databases.hmdb_get_diseases`](#dnakit.databases.hmdb_get_diseases)
- [`dnakit.databases.hmdb_get_metabolite`](#dnakit.databases.hmdb_get_metabolite)
- [`dnakit.databases.hmdb_search`](#dnakit.databases.hmdb_search)
- [`dnakit.databases.open_targets_get_associated_drugs_by_disease_efo_id`](#dnakit.databases.open_targets_get_associated_drugs_by_disease_efo_id)
- [`dnakit.databases.open_targets_get_associated_drugs_by_target_ensembl_id`](#dnakit.databases.open_targets_get_associated_drugs_by_target_ensembl_id)
- [`dnakit.databases.open_targets_get_associated_phenotypes_by_disease_efo_id`](#dnakit.databases.open_targets_get_associated_phenotypes_by_disease_efo_id)
- [`dnakit.databases.open_targets_get_associated_targets_by_disease_efo_id`](#dnakit.databases.open_targets_get_associated_targets_by_disease_efo_id)
- [`dnakit.databases.open_targets_get_credible_set_detail`](#dnakit.databases.open_targets_get_credible_set_detail)
- [`dnakit.databases.open_targets_get_disease_ancestors_parents_by_efo_id`](#dnakit.databases.open_targets_get_disease_ancestors_parents_by_efo_id)
- [`dnakit.databases.open_targets_get_disease_descendants_children_by_efo_id`](#dnakit.databases.open_targets_get_disease_descendants_children_by_efo_id)
- [`dnakit.databases.open_targets_get_disease_description_by_efo_id`](#dnakit.databases.open_targets_get_disease_description_by_efo_id)
- [`dnakit.databases.open_targets_get_disease_id_description_by_name`](#dnakit.databases.open_targets_get_disease_id_description_by_name)
- [`dnakit.databases.open_targets_get_disease_synonyms_by_efo_id`](#dnakit.databases.open_targets_get_disease_synonyms_by_efo_id)
- [`dnakit.databases.open_targets_get_disease_therapeutic_areas_by_efo_id`](#dnakit.databases.open_targets_get_disease_therapeutic_areas_by_efo_id)
- [`dnakit.databases.open_targets_get_diseases_phenotypes_by_target_ensembl`](#dnakit.databases.open_targets_get_diseases_phenotypes_by_target_ensembl)
- [`dnakit.databases.open_targets_get_drug_adverse_events_by_chembl_id`](#dnakit.databases.open_targets_get_drug_adverse_events_by_chembl_id)
- [`dnakit.databases.open_targets_get_drug_chemb_id_by_generic_name`](#dnakit.databases.open_targets_get_drug_chemb_id_by_generic_name)
- [`dnakit.databases.open_targets_get_drug_indications_by_chembl_id`](#dnakit.databases.open_targets_get_drug_indications_by_chembl_id)
- [`dnakit.databases.open_targets_get_drug_mechanisms_of_action_by_chembl_id`](#dnakit.databases.open_targets_get_drug_mechanisms_of_action_by_chembl_id)
- [`dnakit.databases.open_targets_get_drug_warnings_by_chembl_id`](#dnakit.databases.open_targets_get_drug_warnings_by_chembl_id)
- [`dnakit.databases.open_targets_get_evidence_by_datasource`](#dnakit.databases.open_targets_get_evidence_by_datasource)
- [`dnakit.databases.open_targets_get_gwas_study`](#dnakit.databases.open_targets_get_gwas_study)
- [`dnakit.databases.open_targets_get_publications_by_target_ensembl_id`](#dnakit.databases.open_targets_get_publications_by_target_ensembl_id)
- [`dnakit.databases.open_targets_get_similar_entities_by_disease_efo_id`](#dnakit.databases.open_targets_get_similar_entities_by_disease_efo_id)
- [`dnakit.databases.open_targets_get_study_credible_sets`](#dnakit.databases.open_targets_get_study_credible_sets)
- [`dnakit.databases.open_targets_get_target_classes_by_ensembl_id`](#dnakit.databases.open_targets_get_target_classes_by_ensembl_id)
- [`dnakit.databases.open_targets_get_target_gene_ontology_by_ensembl_id`](#dnakit.databases.open_targets_get_target_gene_ontology_by_ensembl_id)
- [`dnakit.databases.open_targets_get_target_genomic_location_by_ensembl_id`](#dnakit.databases.open_targets_get_target_genomic_location_by_ensembl_id)
- [`dnakit.databases.open_targets_get_target_homologues_by_ensembl_id`](#dnakit.databases.open_targets_get_target_homologues_by_ensembl_id)
- [`dnakit.databases.open_targets_get_target_id_description_by_name`](#dnakit.databases.open_targets_get_target_id_description_by_name)
- [`dnakit.databases.open_targets_get_target_interactions_by_ensembl_id`](#dnakit.databases.open_targets_get_target_interactions_by_ensembl_id)
- [`dnakit.databases.open_targets_get_target_safety_profile_by_ensembl_id`](#dnakit.databases.open_targets_get_target_safety_profile_by_ensembl_id)
- [`dnakit.databases.open_targets_get_target_tractability_by_ensembl_id`](#dnakit.databases.open_targets_get_target_tractability_by_ensembl_id)
- [`dnakit.databases.open_targets_get_variant_credible_sets`](#dnakit.databases.open_targets_get_variant_credible_sets)
- [`dnakit.databases.open_targets_get_variant_info`](#dnakit.databases.open_targets_get_variant_info)
- [`dnakit.databases.open_targets_map_any_disease_id_to_all_other_ids`](#dnakit.databases.open_targets_map_any_disease_id_to_all_other_ids)
- [`dnakit.databases.open_targets_multi_entity_search_by_query_string`](#dnakit.databases.open_targets_multi_entity_search_by_query_string)
- [`dnakit.databases.open_targets_search_gwas_studies_by_disease`](#dnakit.databases.open_targets_search_gwas_studies_by_disease)
- [`dnakit.databases.open_targets_target_disease_evidence`](#dnakit.databases.open_targets_target_disease_evidence)
- [`dnakit.databases.pub_chem_get_cid_by_compound_name`](#dnakit.databases.pub_chem_get_cid_by_compound_name)
- [`dnakit.databases.pub_chem_get_compound_properties_by_cid`](#dnakit.databases.pub_chem_get_compound_properties_by_cid)
- [`dnakit.databases.reactome_disease_target_score`](#dnakit.databases.reactome_disease_target_score)

## `dnakit.databases.ch_embl_get_drug_mechanisms` {#dnakit.databases.ch_embl_get_drug_mechanisms}

Get mechanisms of action for a drug by ChEMBL drug ID or drug name. Accepts drug_chembl_id (e.g., 'CHEMBL1201581') or drug_name (e.g., 'trastuzumab') — drug_name triggers automatic ChEMBL ID lookup.

```text
dnakit.databases.ch_embl_get_drug_mechanisms(*, limit=None, offset=None, drug_chembl_id=None, drug_name=None, molecule_chembl_id=None, chembl_id=None, drug_chembl_id__exact=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `limit` | integer | None | Maximum results for the current request, not a database-wide total. {&quot;default&quot;: 20, &quot;maximum&quot;: 1000} |
| `offset` | integer | None | Pagination offset within the current query. {&quot;default&quot;: 0} |
| `drug_chembl_id` | string | None | ChEMBL drug/molecule ID (e.g., &quot;CHEMBL1201580&quot; for adalimumab, &quot;CHEMBL4535757&quot; for sotorasib). |
| `drug_name` | string | None | Drug name for automatic ChEMBL ID lookup (e.g., &quot;trastuzumab&quot;, &quot;lapatinib&quot;, &quot;aspirin&quot;). Case-insensitive. Triggers internal ChEMBL molecule search — use drug_chembl_id if you already know the ID. |
| `molecule_chembl_id` | string | None | Alias for drug_chembl_id. ChEMBL molecule ID (e.g., &quot;CHEMBL25&quot; for aspirin). |
| `chembl_id` | string | None | Alias for drug_chembl_id. ChEMBL ID (e.g., &quot;CHEMBL3301622&quot;). |
| `drug_chembl_id__exact` | string | None | DNAKit optional alias for drug_chembl_id |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "limit": {
          "type": "integer",
          "default": 20,
          "maximum": 1000
        },
        "offset": {
          "type": "integer",
          "default": 0
        },
        "drug_chembl_id": {
          "type": "string",
          "description": "ChEMBL drug/molecule ID (e.g., \"CHEMBL1201580\" for adalimumab, \"CHEMBL4535757\" for sotorasib)."
        },
        "drug_name": {
          "type": "string",
          "description": "Drug name for automatic ChEMBL ID lookup (e.g., \"trastuzumab\", \"lapatinib\", \"aspirin\"). Case-insensitive. Triggers internal ChEMBL molecule search — use drug_chembl_id if you already know the ID."
        },
        "molecule_chembl_id": {
          "type": "string",
          "description": "Alias for drug_chembl_id. ChEMBL molecule ID (e.g., \"CHEMBL25\" for aspirin)."
        },
        "chembl_id": {
          "type": "string",
          "description": "Alias for drug_chembl_id. ChEMBL ID (e.g., \"CHEMBL3301622\")."
        },
        "drug_chembl_id__exact": {
          "type": "string",
          "description": "DNAKit optional alias for drug_chembl_id"
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
            "mechanisms": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "mechanism_of_action": {
                    "type": "string"
                  },
                  "target_chembl_id": {
                    "type": "string"
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

??? info "DNAKit argument extension rules"

    ```json
    {
      "aliases": {
        "drug_chembl_id__exact": "drug_chembl_id"
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `ChEMBL_get_drug_mechanisms`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("ch_embl_get_drug_mechanisms", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ch_embl_get_molecule` {#dnakit.databases.ch_embl_get_molecule}

Get detailed information about a molecule by its ChEMBL ID. Returns molecule properties, structures, synonyms, and related data.

```text
dnakit.databases.ch_embl_get_molecule(*, chembl_id=None, format=None, molecule_chembl_id=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `chembl_id` | string | None | ChEMBL molecule ID, e.g., &#x27;CHEMBL25&#x27; |
| `format` | string | None | Response format {&quot;enum&quot;: &#91;&quot;json&quot;, &quot;xml&quot;, &quot;yaml&quot;&#93;, &quot;default&quot;: &quot;json&quot;} |
| `molecule_chembl_id` | string | None | Alias for chembl_id. ChEMBL molecule ID, e.g., CHEMBL1229517. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "chembl_id": {
          "type": "string",
          "description": "ChEMBL molecule ID, e.g., 'CHEMBL25'"
        },
        "format": {
          "type": "string",
          "enum": [
            "json",
            "xml",
            "yaml"
          ],
          "default": "json",
          "description": "Response format"
        },
        "molecule_chembl_id": {
          "type": "string",
          "description": "Alias for chembl_id. ChEMBL molecule ID, e.g., CHEMBL1229517."
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
          "type": "string",
          "description": "Request status"
        },
        "data": {
          "type": "object",
          "description": "Molecule data including structures, properties, synonyms"
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `ChEMBL_get_molecule`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("ch_embl_get_molecule", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ch_embl_get_target_activities` {#dnakit.databases.ch_embl_get_target_activities}

Get all activity data for a target by ChEMBL target ID. Returns bioactivity measurements (IC50, Ki, etc.) for molecules tested against this target. Note: activity values may be in different units (nM or uM) depending on the assay — check the 'standard_units' field in each result. To find a target ID, use ChEMBL_search_targets with a target name or gene symbol.

```text
dnakit.databases.ch_embl_get_target_activities(*, target_chembl_id__exact=None, limit=None, offset=None, target_chembl_id=None, chembl_id=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `target_chembl_id__exact` | string | None | ChEMBL target ID (e.g., &#x27;CHEMBL2074&#x27;). To find a target ID, use ChEMBL_search_targets with a target name or gene symbol. |
| `limit` | integer | None | Maximum results for the current request, not a database-wide total. {&quot;default&quot;: 20, &quot;maximum&quot;: 1000} |
| `offset` | integer | None | Pagination offset within the current query. {&quot;default&quot;: 0} |
| `target_chembl_id` | string | None | Alias for target_chembl_id__exact. ChEMBL target ID (e.g., CHEMBL213). |
| `chembl_id` | string | None | DNAKit optional alias for target_chembl_id |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "target_chembl_id__exact": {
          "type": "string",
          "description": "ChEMBL target ID (e.g., 'CHEMBL2074'). To find a target ID, use ChEMBL_search_targets with a target name or gene symbol."
        },
        "limit": {
          "type": "integer",
          "default": 20,
          "maximum": 1000
        },
        "offset": {
          "type": "integer",
          "default": 0
        },
        "target_chembl_id": {
          "type": "string",
          "description": "Alias for target_chembl_id__exact. ChEMBL target ID (e.g., CHEMBL213)."
        },
        "chembl_id": {
          "type": "string",
          "description": "DNAKit optional alias for target_chembl_id"
        }
      },
      "oneOf": [
        {
          "required": [
            "target_chembl_id__exact"
          ]
        },
        {
          "required": [
            "target_chembl_id"
          ]
        }
      ],
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
            "activities": {
              "type": "array",
              "items": {
                "type": "object"
              }
            }
          }
        },
        "count": {
          "type": "integer"
        }
      }
    }
    ```

??? info "DNAKit argument extension rules"

    ```json
    {
      "aliases": {
        "chembl_id": "target_chembl_id"
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `ChEMBL_get_target_activities`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("ch_embl_get_target_activities", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ch_embl_search_drugs` {#dnakit.databases.ch_embl_search_drugs}

Search drugs by name, approval status, or other criteria. Use the `query` parameter to search by drug name (e.g., {"query": "sotorasib"}). Use `max_phase` to filter by development phase.

```text
dnakit.databases.ch_embl_search_drugs(*, query=None, max_phase=None, limit=None, offset=None, molecule_chembl_id=None, q=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | None | Drug name to search for (partial match, case-insensitive). E.g., &quot;sotorasib&quot;, &quot;olaparib&quot;, &quot;imatinib&quot;. |
| `max_phase` | integer | None | Filter by maximum development phase (0-4) |
| `limit` | integer | None | Maximum results for the current request, not a database-wide total. {&quot;default&quot;: 20, &quot;maximum&quot;: 1000} |
| `offset` | integer | None | Pagination offset within the current query. {&quot;default&quot;: 0} |
| `molecule_chembl_id` | string | None | Filter by ChEMBL molecule ID (e.g., &quot;CHEMBL1201580&quot; for adalimumab). |
| `q` | string | None | DNAKit optional alias for query |
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
          "description": "Drug name to search for (partial match, case-insensitive). E.g., \"sotorasib\", \"olaparib\", \"imatinib\"."
        },
        "max_phase": {
          "type": "integer",
          "description": "Filter by maximum development phase (0-4)"
        },
        "limit": {
          "type": "integer",
          "default": 20,
          "maximum": 1000
        },
        "offset": {
          "type": "integer",
          "default": 0
        },
        "molecule_chembl_id": {
          "type": "string",
          "description": "Filter by ChEMBL molecule ID (e.g., \"CHEMBL1201580\" for adalimumab)."
        },
        "q": {
          "type": "string",
          "description": "DNAKit optional alias for query"
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
            "drugs": {
              "type": "array",
              "items": {
                "type": "object"
              }
            }
          }
        },
        "count": {
          "type": "integer"
        }
      }
    }
    ```

??? info "DNAKit argument extension rules"

    ```json
    {
      "aliases": {
        "q": "query"
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `ChEMBL_search_drugs`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("ch_embl_search_drugs", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ch_embl_search_targets` {#dnakit.databases.ch_embl_search_targets}

Search targets by name, organism, target type, or other criteria. Use this tool to find target ChEMBL IDs needed for other tools.

```text
dnakit.databases.ch_embl_search_targets(*, target_chembl_id=None, pref_name__contains=None, organism=None, target_type=None, fields=None, limit=None, offset=None, q=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `target_chembl_id` | string | None | Filter by target ChEMBL ID |
| `pref_name__contains` | string | None | Filter by target name (contains) |
| `organism` | string | None | Filter by organism (e.g., &#x27;Homo sapiens&#x27;) |
| `target_type` | string | None | Filter by target type (e.g., &#x27;SINGLE PROTEIN&#x27;, &#x27;PROTEIN COMPLEX&#x27;) |
| `fields` | array | None | Optional list of ChEMBL target fields to include in each returned target object (projection). ToolUniverse maps this to ChEMBL&#x27;s `only` query parameter (comma-separated). Supported fields: target_chembl_id, pref_name, organism, target_type, target_components. {&quot;uniqueItems&quot;: true, &quot;items&quot;: {&quot;type&quot;: &quot;string&quot;, &quot;enum&quot;: &#91;&quot;target_chembl_id&quot;, &quot;pref_name&quot;, &quot;organism&quot;, &quot;target_type&quot;, &quot;target_components&quot;&#93;}} |
| `limit` | integer | None | Maximum results for the current request, not a database-wide total. {&quot;default&quot;: 20, &quot;maximum&quot;: 1000} |
| `offset` | integer | None | Pagination offset within the current query. {&quot;default&quot;: 0} |
| `q` | string | None | DNAKit optional alias for pref_name__contains |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "target_chembl_id": {
          "type": "string",
          "description": "Filter by target ChEMBL ID"
        },
        "pref_name__contains": {
          "type": "string",
          "description": "Filter by target name (contains)"
        },
        "organism": {
          "type": "string",
          "description": "Filter by organism (e.g., 'Homo sapiens')"
        },
        "target_type": {
          "type": "string",
          "description": "Filter by target type (e.g., 'SINGLE PROTEIN', 'PROTEIN COMPLEX')"
        },
        "fields": {
          "type": "array",
          "uniqueItems": true,
          "items": {
            "type": "string",
            "enum": [
              "target_chembl_id",
              "pref_name",
              "organism",
              "target_type",
              "target_components"
            ]
          },
          "description": "Optional list of ChEMBL target fields to include in each returned target object (projection). ToolUniverse maps this to ChEMBL's `only` query parameter (comma-separated). Supported fields: target_chembl_id, pref_name, organism, target_type, target_components."
        },
        "limit": {
          "type": "integer",
          "default": 20,
          "maximum": 1000
        },
        "offset": {
          "type": "integer",
          "default": 0
        },
        "q": {
          "type": "string",
          "description": "DNAKit optional alias for pref_name__contains"
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
            "targets": {
              "type": "array",
              "items": {
                "type": "object",
                "description": "Target objects returned by ChEMBL. When 'fields' is provided, each object will contain only the selected fields."
              }
            }
          }
        },
        "count": {
          "type": "integer"
        }
      }
    }
    ```

??? info "DNAKit argument extension rules"

    ```json
    {
      "aliases": {
        "q": "pref_name__contains"
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `ChEMBL_search_targets`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("ch_embl_search_targets", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ctd_get_chemical_gene_interactions` {#dnakit.databases.ctd_get_chemical_gene_interactions}

PERMANENTLY UNAVAILABLE: this tool always returns an error. CTD's native batchQuery.go is CAPTCHA-blocked for programmatic clients, and the RENCI Automat mirror this tool used to fall back on has been fully decommissioned (its registry no longer lists a 'ctd' backend at all). No other free, live source of CTD's curated chemical-gene interactions was found. Use DGIdb_get_drug_gene_interactions or ChEMBL/OpenTargets tools for target-compound bioactivity data instead.

```text
dnakit.databases.ctd_get_chemical_gene_interactions(*, input_terms, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `input_terms` | string | Required by input schema | Chemical name, CAS RN, or MeSH ID. Unused -- this tool always returns an error; see the tool description. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "input_terms": {
          "type": "string",
          "description": "Chemical name, CAS RN, or MeSH ID. Unused -- this tool always returns an error; see the tool description."
        }
      },
      "required": [
        "input_terms"
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
              "source_id": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "CURIE of the subject node (e.g., CHEBI:15365)"
              },
              "source_name": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Human-readable name of the subject"
              },
              "target_id": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "CURIE of the object node (e.g., NCBIGene:7076)"
              },
              "target_name": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Human-readable name of the object"
              },
              "predicate": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "biolink predicate (e.g., biolink:affects)"
              },
              "qualified_predicate": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Causal flavour of the predicate (e.g., biolink:causes)"
              },
              "object_direction_qualifier": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Direction of effect on the object (e.g., increased, decreased)"
              },
              "knowledge_level": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "How the edge was obtained (e.g., knowledge_assertion)"
              },
              "agent_type": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Who created the edge (e.g., manual_agent)"
              },
              "primary_knowledge_source": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "InfoRes id of the source database (e.g., infores:ctd)"
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

**Execution and provenance:** `network`；upstream name `CTD_get_chemical_gene_interactions`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("ctd_get_chemical_gene_interactions", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.dgidb_get_drug_gene_interactions` {#dnakit.databases.dgidb_get_drug_gene_interactions}

Get drug-gene interactions for a list of genes from DGIdb. Returns drugs targeting each gene, interaction types, and evidence sources. Gene symbols DGIdb does not recognize are omitted from the results and named in metadata.unmatched_genes, so a missing gene is distinguishable from a gene with no interactions.

```text
dnakit.databases.dgidb_get_drug_gene_interactions(*, genes=None, gene_name=None, gene=None, interaction_sources=None, interaction_types=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `genes` | array | None | List of gene symbols (e.g., &#91;&#x27;EGFR&#x27;, &#x27;BRAF&#x27;&#93;). Also accepts a single gene as string. Aliases: gene_name, gene. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `gene_name` | string | None | Alias for genes. Single gene symbol (e.g., &#x27;EGFR&#x27;). |
| `gene` | string | None | Alias for genes. Single gene symbol (e.g., &#x27;EGFR&#x27;). |
| `interaction_sources` | array | None | Optional filter by data sources (e.g., &#91;&#x27;DrugBank&#x27;, &#x27;ChEMBL&#x27;&#93;). {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `interaction_types` | array | None | Optional filter by interaction types (e.g., &#91;&#x27;inhibitor&#x27;, &#x27;antagonist&#x27;&#93;). {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "genes": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "List of gene symbols (e.g., ['EGFR', 'BRAF']). Also accepts a single gene as string. Aliases: gene_name, gene."
        },
        "gene_name": {
          "type": "string",
          "description": "Alias for genes. Single gene symbol (e.g., 'EGFR')."
        },
        "gene": {
          "type": "string",
          "description": "Alias for genes. Single gene symbol (e.g., 'EGFR')."
        },
        "interaction_sources": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "Optional filter by data sources (e.g., ['DrugBank', 'ChEMBL'])."
        },
        "interaction_types": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "Optional filter by interaction types (e.g., ['inhibitor', 'antagonist'])."
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
        "data": {
          "type": "object",
          "properties": {
            "genes": {
              "type": "object",
              "properties": {
                "nodes": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "name": {
                        "type": "string"
                      },
                      "longName": {
                        "type": "string"
                      },
                      "interactions": {
                        "type": "array",
                        "items": {
                          "type": "object",
                          "properties": {
                            "drug": {
                              "type": "object"
                            },
                            "interactionTypes": {
                              "type": "array"
                            },
                            "sources": {
                              "type": "array"
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `DGIdb_get_drug_gene_interactions`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("dgidb_get_drug_gene_interactions", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.dgidb_get_gene_druggability` {#dnakit.databases.dgidb_get_gene_druggability}

Get druggability information for genes. Returns gene categories indicating if a gene is druggable (kinase, GPCR, ion channel, etc.). Gene symbols DGIdb does not recognize are omitted from the results and named in metadata.unmatched_genes.

```text
dnakit.databases.dgidb_get_gene_druggability(*, genes=None, gene_name=None, gene=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `genes` | array | None | List of gene symbols to check druggability. Aliases: gene_name, gene. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `gene_name` | string | None | Alias for genes. Single gene symbol (e.g., &#x27;EGFR&#x27;). |
| `gene` | string | None | Alias for genes. Single gene symbol (e.g., &#x27;EGFR&#x27;). |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "genes": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "List of gene symbols to check druggability. Aliases: gene_name, gene."
        },
        "gene_name": {
          "type": "string",
          "description": "Alias for genes. Single gene symbol (e.g., 'EGFR')."
        },
        "gene": {
          "type": "string",
          "description": "Alias for genes. Single gene symbol (e.g., 'EGFR')."
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
        "data": {
          "type": "object",
          "properties": {
            "genes": {
              "type": "object",
              "properties": {
                "nodes": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "name": {
                        "type": "string"
                      },
                      "geneCategories": {
                        "type": "array"
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `DGIdb_get_gene_druggability`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("dgidb_get_gene_druggability", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.drugbank_get_drug_basic_info_by_drug_name_or_id` {#dnakit.databases.drugbank_get_drug_basic_info_by_drug_name_or_id}

Get basic drug information including name, description, CAS number, and approval status by drug name or DrugBank ID.

```text
dnakit.databases.drugbank_get_drug_basic_info_by_drug_name_or_id(*, query, case_sensitive=None, exact_match=None, limit=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required by input schema | Drug name or DrugBank ID to search for {&quot;minLength&quot;: 1} |
| `case_sensitive` | boolean | None | Select True to perform a case-sensitive search {&quot;default&quot;: false} |
| `exact_match` | boolean | None | Select True to require an exact match with the queried name or ID {&quot;default&quot;: false} |
| `limit` | integer | None | Maximum number of results to return {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
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
          "description": "Drug name or DrugBank ID to search for",
          "minLength": 1
        },
        "case_sensitive": {
          "type": "boolean",
          "description": "Select True to perform a case-sensitive search",
          "default": false
        },
        "exact_match": {
          "type": "boolean",
          "description": "Select True to require an exact match with the queried name or ID",
          "default": false
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of results to return",
          "default": 10,
          "minimum": 1,
          "maximum": 100
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
        "error": {
          "type": "string",
          "description": "Error message (if any)"
        },
        "query": {
          "type": "string",
          "description": "The original search query"
        },
        "total_matches": {
          "type": "integer",
          "description": "Total number of matching records found"
        },
        "total_returned_results": {
          "type": "integer",
          "description": "Number of results returned in this response"
        },
        "results": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "_tag": {
                "type": "string",
                "description": "XML tag name of the record element"
              },
              "_text": {
                "type": "string",
                "description": "Text content of the record element"
              },
              "_attributes": {
                "type": "object",
                "description": "Attributes of the record element"
              },
              "drug_name": {
                "type": "string",
                "description": "Name of the drug"
              },
              "drugbank_id": {
                "type": "string",
                "description": "Primary DrugBank ID"
              },
              "secondary_ids": {
                "type": "string",
                "description": "Secondary DrugBank IDs, separated by ' | '"
              },
              "description": {
                "type": "string",
                "description": "Drug description"
              },
              "cas_number": {
                "type": "string",
                "description": "Chemical Abstracts Service registry number"
              },
              "unii": {
                "type": "string",
                "description": "Unique Ingredient Identifier"
              },
              "state": {
                "type": "string",
                "description": "Physical state of the drug"
              },
              "drug_type": {
                "type": "string",
                "description": "Type of drug (e.g., small molecule, biotech)"
              },
              "approval_groups": {
                "type": "string",
                "description": "Drug approval status groups, separated by ' | '"
              },
              "synonyms": {
                "type": "string",
                "description": "Drug synonyms, separated by ' | '"
              },
              "entry_creation_date": {
                "type": "string",
                "description": "Date when the entry was created"
              },
              "entry_last_updated": {
                "type": "string",
                "description": "Date when the entry was last updated"
              },
              "matched_fields": {
                "type": "array",
                "items": {
                  "type": "string"
                },
                "description": "List of fields that matched the search query"
              }
            }
          }
        },
        "search_parameters": {
          "type": "object",
          "properties": {
            "case_sensitive": {
              "type": "boolean",
              "description": "Whether the search was case-sensitive"
            },
            "exact_match": {
              "type": "boolean",
              "description": "Whether exact matching was used"
            },
            "limit": {
              "type": "integer",
              "description": "Maximum number of results requested"
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `drugbank_get_drug_basic_info_by_drug_name_or_id`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("drugbank_get_drug_basic_info_by_drug_name_or_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.drugbank_get_pharmacology_by_drug_name_or_drugbank_id` {#dnakit.databases.drugbank_get_pharmacology_by_drug_name_or_drugbank_id}

Get drug pharmacodynamics, mechanism of action, and pharmacokinetics by drug name or Drugbank ID.

```text
dnakit.databases.drugbank_get_pharmacology_by_drug_name_or_drugbank_id(*, query, case_sensitive=None, exact_match=None, limit=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required by input schema | Drug name or Drugbank ID to search for {&quot;minLength&quot;: 1} |
| `case_sensitive` | boolean | None | Select True to perform a case-sensitive search {&quot;default&quot;: false} |
| `exact_match` | boolean | None | Select True to require an exact match {&quot;default&quot;: false} |
| `limit` | integer | None | Maximum number of results to return {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
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
          "description": "Drug name or Drugbank ID to search for",
          "minLength": 1
        },
        "case_sensitive": {
          "type": "boolean",
          "description": "Select True to perform a case-sensitive search",
          "default": false
        },
        "exact_match": {
          "type": "boolean",
          "description": "Select True to require an exact match",
          "default": false
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of results to return",
          "default": 10,
          "minimum": 1,
          "maximum": 100
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
        "error": {
          "type": "string",
          "description": "Error message (if any)"
        },
        "query": {
          "type": "string",
          "description": "The original search query"
        },
        "total_matches": {
          "type": "integer",
          "description": "Total number of matching records found"
        },
        "total_returned_results": {
          "type": "integer",
          "description": "Number of results returned in this response"
        },
        "results": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "_tag": {
                "type": "string",
                "description": "XML tag name of the record element"
              },
              "_text": {
                "type": "string",
                "description": "Text content of the record element"
              },
              "_attributes": {
                "type": "object",
                "description": "Attributes of the record element"
              },
              "drug_name": {
                "type": "string",
                "description": "Name of the drug"
              },
              "drugbank_id": {
                "type": "string",
                "description": "Primary DrugBank ID"
              },
              "synonyms": {
                "type": "string",
                "description": "Drug synonyms, separated by ' | '"
              },
              "pharmacodynamics": {
                "type": "string",
                "description": "Pharmacodynamics information describing the drug's effects on the body"
              },
              "mechanism_of_action": {
                "type": "string",
                "description": "Mechanism of action describing how the drug works at the molecular level"
              },
              "absorption": {
                "type": "string",
                "description": "Information about drug absorption and bioavailability"
              },
              "metabolism": {
                "type": "string",
                "description": "Information about drug metabolism and metabolic pathways"
              },
              "half_life": {
                "type": "string",
                "description": "Plasma elimination half-life of the drug"
              },
              "protein_binding": {
                "type": "string",
                "description": "Information about protein binding properties"
              },
              "route_of_elimination": {
                "type": "string",
                "description": "Information about how the drug is eliminated from the body"
              },
              "volume_of_distribution": {
                "type": "string",
                "description": "Volume of distribution of the drug"
              },
              "clearance": {
                "type": "string",
                "description": "Information about drug clearance"
              },
              "matched_fields": {
                "type": "array",
                "items": {
                  "type": "string"
                },
                "description": "List of fields that matched the search query"
              }
            }
          }
        },
        "search_parameters": {
          "type": "object",
          "properties": {
            "case_sensitive": {
              "type": "boolean",
              "description": "Whether the search was case-sensitive"
            },
            "exact_match": {
              "type": "boolean",
              "description": "Whether exact matching was used"
            },
            "limit": {
              "type": "integer",
              "description": "Maximum number of results requested"
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `drugbank_get_pharmacology_by_drug_name_or_drugbank_id`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("drugbank_get_pharmacology_by_drug_name_or_drugbank_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.drugbank_get_targets_by_drug_name_or_drugbank_id` {#dnakit.databases.drugbank_get_targets_by_drug_name_or_drugbank_id}

Get drug targets, enzymes, carriers, and transporters by drug name or DrugBank ID.

```text
dnakit.databases.drugbank_get_targets_by_drug_name_or_drugbank_id(*, query, case_sensitive=None, exact_match=None, limit=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required by input schema | Drug name or Drugbank ID to search for {&quot;minLength&quot;: 1} |
| `case_sensitive` | boolean | None | Select True to perform a case-sensitive search {&quot;default&quot;: false} |
| `exact_match` | boolean | None | Select True to require an exact match {&quot;default&quot;: false} |
| `limit` | integer | None | Maximum number of results to return {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
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
          "description": "Drug name or Drugbank ID to search for",
          "minLength": 1
        },
        "case_sensitive": {
          "type": "boolean",
          "description": "Select True to perform a case-sensitive search",
          "default": false
        },
        "exact_match": {
          "type": "boolean",
          "description": "Select True to require an exact match",
          "default": false
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of results to return",
          "default": 10,
          "minimum": 1,
          "maximum": 100
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
        "error": {
          "type": "string",
          "description": "Error message (if any)"
        },
        "query": {
          "type": "string",
          "description": "The original search query"
        },
        "total_matches": {
          "type": "integer",
          "description": "Total number of matching records found"
        },
        "total_returned_results": {
          "type": "integer",
          "description": "Number of results returned in this response"
        },
        "results": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "_tag": {
                "type": "string",
                "description": "XML tag name of the record element"
              },
              "_text": {
                "type": "string",
                "description": "Text content of the record element"
              },
              "_attributes": {
                "type": "object",
                "description": "Attributes of the record element"
              },
              "drug_name": {
                "type": "string",
                "description": "Name of the drug"
              },
              "drugbank_id": {
                "type": "string",
                "description": "Primary DrugBank ID"
              },
              "synonyms": {
                "type": "string",
                "description": "Drug synonyms, separated by ' | '"
              },
              "targets": {
                "type": "array",
                "items": {
                  "type": "object",
                  "properties": {
                    "id": {
                      "type": "string",
                      "description": "Target ID"
                    },
                    "name": {
                      "type": "string",
                      "description": "Target name"
                    },
                    "organism": {
                      "type": "string",
                      "description": "Target organism"
                    },
                    "actions": {
                      "type": "string",
                      "description": "Target actions, separated by ' | '"
                    }
                  }
                },
                "description": "Drug targets"
              },
              "enzymes": {
                "type": "array",
                "items": {
                  "type": "object",
                  "properties": {
                    "id": {
                      "type": "string",
                      "description": "Enzyme ID"
                    },
                    "name": {
                      "type": "string",
                      "description": "Enzyme name"
                    },
                    "organism": {
                      "type": "string",
                      "description": "Enzyme organism"
                    },
                    "actions": {
                      "type": "string",
                      "description": "Enzyme actions, separated by ' | '"
                    }
                  }
                },
                "description": "Drug enzymes"
              },
              "carriers": {
                "type": "array",
                "items": {
                  "type": "object",
                  "properties": {
                    "id": {
                      "type": "string",
                      "description": "Carrier ID"
                    },
                    "name": {
                      "type": "string",
                      "description": "Carrier name"
                    },
                    "organism": {
                      "type": "string",
                      "description": "Carrier organism"
                    },
                    "actions": {
                      "type": "string",
                      "description": "Carrier actions, separated by ' | '"
                    }
                  }
                },
                "description": "Drug carriers"
              },
              "transporters": {
                "type": "array",
                "items": {
                  "type": "object",
                  "properties": {
                    "id": {
                      "type": "string",
                      "description": "Transporter ID"
                    },
                    "name": {
                      "type": "string",
                      "description": "Transporter name"
                    },
                    "organism": {
                      "type": "string",
                      "description": "Transporter organism"
                    },
                    "actions": {
                      "type": "string",
                      "description": "Transporter actions, separated by ' | '"
                    }
                  }
                },
                "description": "Drug transporters"
              },
              "matched_fields": {
                "type": "array",
                "items": {
                  "type": "string"
                },
                "description": "List of fields that matched the search query"
              }
            }
          }
        },
        "search_parameters": {
          "type": "object",
          "properties": {
            "case_sensitive": {
              "type": "boolean",
              "description": "Whether the search was case-sensitive"
            },
            "exact_match": {
              "type": "boolean",
              "description": "Whether exact matching was used"
            },
            "limit": {
              "type": "integer",
              "description": "Maximum number of results requested"
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `drugbank_get_targets_by_drug_name_or_drugbank_id`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("drugbank_get_targets_by_drug_name_or_drugbank_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.drugbank_vocab_search` {#dnakit.databases.drugbank_vocab_search}

Search the DrugBank vocabulary dataset for drugs by name, ID, synonyms, or other fields using text-based queries. Returns detailed drug information including DrugBank ID, common name, CAS number, UNII, and synonyms.

```text
dnakit.databases.drugbank_vocab_search(*, query, search_fields=None, case_sensitive=None, exact_match=None, limit=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required by input schema | Search query string. Can be drug name, synonym, DrugBank ID, or any text to search for. {&quot;minLength&quot;: 1} |
| `search_fields` | array | None | Fields to search in. Available fields: &#x27;DrugBank ID&#x27;, &#x27;Accession Numbers&#x27;, &#x27;Common name&#x27;, &#x27;CAS&#x27;, &#x27;UNII&#x27;, &#x27;Synonyms&#x27;, &#x27;Standard InChI Key&#x27;. Defaults to &#91;&#x27;Common name&#x27;, &#x27;Synonyms&#x27;, &#x27;DrugBank ID&#x27;&#93;. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `case_sensitive` | boolean | None | Whether the search should be case sensitive. {&quot;default&quot;: false} |
| `exact_match` | boolean | None | Whether to perform exact matching instead of substring matching. {&quot;default&quot;: false} |
| `limit` | integer | None | Maximum number of results to return. {&quot;minimum&quot;: 1, &quot;maximum&quot;: 1000, &quot;default&quot;: 10} |
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
          "description": "Search query string. Can be drug name, synonym, DrugBank ID, or any text to search for.",
          "minLength": 1
        },
        "search_fields": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "Fields to search in. Available fields: 'DrugBank ID', 'Accession Numbers', 'Common name', 'CAS', 'UNII', 'Synonyms', 'Standard InChI Key'. Defaults to ['Common name', 'Synonyms', 'DrugBank ID']."
        },
        "case_sensitive": {
          "type": "boolean",
          "description": "Whether the search should be case sensitive.",
          "default": false
        },
        "exact_match": {
          "type": "boolean",
          "description": "Whether to perform exact matching instead of substring matching.",
          "default": false
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of results to return.",
          "minimum": 1,
          "maximum": 1000,
          "default": 10
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
      "description": "Search results with matched drug records and search metadata",
      "properties": {
        "query": {
          "type": "string",
          "description": "The original search query string"
        },
        "total_results": {
          "type": "integer",
          "description": "Total number of matching results found"
        },
        "results": {
          "type": "array",
          "description": "Array of matching drug records",
          "items": {
            "type": "object",
            "properties": {
              "DrugBank ID": {
                "type": "string",
                "description": "Drug identifier (e.g., 'DB00945')"
              },
              "Accession Numbers": {
                "type": "string",
                "description": "Alternative drug identifiers"
              },
              "Common name": {
                "type": "string",
                "description": "Drug common name (e.g., 'Acetylsalicylic acid')"
              },
              "CAS": {
                "type": "string",
                "description": "Chemical Abstracts Service number"
              },
              "UNII": {
                "type": "string",
                "description": "Unique Ingredient Identifier"
              },
              "Synonyms": {
                "type": "string",
                "description": "Alternative drug names and synonyms"
              },
              "Standard InChI Key": {
                "type": "string",
                "description": "International Chemical Identifier key"
              },
              "matched_fields": {
                "type": "array",
                "items": {
                  "type": "string"
                },
                "description": "Fields that matched the search query (e.g., ['Synonyms', 'Common name'])"
              }
            }
          }
        },
        "search_parameters": {
          "type": "object",
          "description": "Parameters used for the search",
          "properties": {
            "search_fields": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "case_sensitive": {
              "type": "boolean"
            },
            "exact_match": {
              "type": "boolean"
            },
            "limit": {
              "type": "integer"
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `drugbank_vocab_search`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("drugbank_vocab_search", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.europepmc_disease_target_score` {#dnakit.databases.europepmc_disease_target_score}

Extract disease-target association scores from Europe PMC literature. This includes literature-based evidence. Returns EVERY target with a Europe PMC score for the disease, sorted strongest-first -- it is not a paged or top-N view. Scanning is capped at 25 seconds; when the cap is hit the response sets 'truncated' and a 'note' saying how many of the disease's associated targets were scanned.

```text
dnakit.databases.europepmc_disease_target_score(*, efoId, pageSize=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | Required by input schema | The EFO (Experimental Factor Ontology) ID of the disease, e.g., &#x27;MONDO_0011996&#x27; for chronic myelogenous leukemia |
| `pageSize` | integer | None | Upstream fetch batch size (max 100), NOT a limit on results returned. The tool pages through all of the disease&#x27;s associated targets regardless of this value -- &#x27;pageSize&#x27;: 10 on MONDO_0005011 (Crohn disease) still returns 390 scored targets. Lower it only to make each upstream request smaller. {&quot;default&quot;: 100} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "efoId": {
          "type": "string",
          "description": "The EFO (Experimental Factor Ontology) ID of the disease, e.g., 'MONDO_0011996' for chronic myelogenous leukemia"
        },
        "pageSize": {
          "type": "integer",
          "description": "Upstream fetch batch size (max 100), NOT a limit on results returned. The tool pages through all of the disease's associated targets regardless of this value -- 'pageSize': 10 on MONDO_0005011 (Crohn disease) still returns 390 scored targets. Lower it only to make each upstream request smaller.",
          "default": 100
        }
      },
      "required": [
        "efoId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "description": "Europe PMC literature disease-target association scores",
      "properties": {
        "disease_info": {
          "type": "object",
          "properties": {
            "disease_id": {
              "type": "string"
            },
            "disease_name": {
              "type": "string"
            }
          }
        },
        "datasource": {
          "type": "string",
          "description": "Always 'europepmc'"
        },
        "total_targets_with_scores": {
          "type": "integer"
        },
        "target_scores": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "target_symbol": {
                "type": "string"
              },
              "target_id": {
                "type": "string"
              },
              "datasource": {
                "type": "string"
              },
              "score": {
                "type": "number",
                "description": "Europe PMC literature score (0-1)"
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `europepmc_disease_target_score`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("europepmc_disease_target_score", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.expression_atlas_disease_target_score` {#dnakit.databases.expression_atlas_disease_target_score}

Extract disease-target association scores from Expression Atlas. This provides gene expression data.

```text
dnakit.databases.expression_atlas_disease_target_score(*, efoId, pageSize=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | Required by input schema | The EFO (Experimental Factor Ontology) ID of the disease, e.g., &#x27;MONDO_0011996&#x27; for chronic myelogenous leukemia |
| `pageSize` | integer | None | Number of results per page (default: 100, max: 100) {&quot;default&quot;: 100} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "efoId": {
          "type": "string",
          "description": "The EFO (Experimental Factor Ontology) ID of the disease, e.g., 'MONDO_0011996' for chronic myelogenous leukemia"
        },
        "pageSize": {
          "type": "integer",
          "description": "Number of results per page (default: 100, max: 100)",
          "default": 100
        }
      },
      "required": [
        "efoId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "description": "Expression Atlas gene expression disease-target association scores",
      "properties": {
        "disease_info": {
          "type": "object",
          "properties": {
            "disease_id": {
              "type": "string"
            },
            "disease_name": {
              "type": "string"
            }
          }
        },
        "datasource": {
          "type": "string",
          "description": "Always 'expression_atlas'"
        },
        "total_targets_with_scores": {
          "type": "integer"
        },
        "target_scores": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "target_symbol": {
                "type": "string"
              },
              "target_id": {
                "type": "string"
              },
              "datasource": {
                "type": "string"
              },
              "score": {
                "type": "number",
                "description": "Expression Atlas score (0-1)"
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `expression_atlas_disease_target_score`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("expression_atlas_disease_target_score", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.fda_get_active_ingredient_info_by_drug_name` {#dnakit.databases.fda_get_active_ingredient_info_by_drug_name}

Fetch a list of active ingredients in a specific drug product.

```text
dnakit.databases.fda_get_active_ingredient_info_by_drug_name(*, drug_name, limit=None, skip=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `drug_name` | string | Required by input schema | The name of the drug. |
| `limit` | integer | None | The number of records to return. |
| `skip` | integer | None | The number of records to skip. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "drug_name": {
          "type": "string",
          "description": "The name of the drug."
        },
        "limit": {
          "type": "integer",
          "description": "The number of records to return."
        },
        "skip": {
          "type": "integer",
          "description": "The number of records to skip."
        }
      },
      "required": [
        "drug_name"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "meta": {
          "type": "object",
          "properties": {
            "skip": {
              "type": "integer",
              "description": "Number of upstream records skipped before this page."
            },
            "limit": {
              "type": "integer",
              "description": "Maximum number of records requested from openFDA for this call (page size)."
            },
            "total": {
              "type": "integer",
              "description": "Total number of label records openFDA matched UPSTREAM, before any local processing. This is NOT the length of 'results': local deduplication and the applied limit both shrink the list. Use 'result_count' for the number of rows actually returned and 'duplicates_removed' to reconcile the two."
            }
          }
        },
        "results": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "openfda.brand_name": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                }
              },
              "openfda.generic_name": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                }
              },
              "active_ingredient": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                }
              },
              "spl_product_data_elements": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                },
                "description": "Product/ingredient identity string from the SPL record (e.g. 'Ethyol amifostine AMIFOSTINE AMIFOSTINE'). Populated on records whose 'openfda' block is empty and therefore whose openfda.brand_name / openfda.generic_name come back null; use it to confirm which product a row describes when those two fields are unavailable."
              }
            }
          }
        },
        "result_count": {
          "type": "integer",
          "description": "Number of rows in 'results' -- i.e. after byte-identical duplicate labels were removed and the caller's limit was applied. Compare against meta.total, which is the upstream pre-deduplication hit count."
        },
        "duplicates_removed": {
          "type": "integer",
          "description": "Number of fetched records dropped because another returned record was byte-identical to them. 0 when nothing was deduplicated."
        },
        "dedup_note": {
          "type": "string",
          "description": "Present only when local deduplication ran. Explains how meta.total, 'result_count' and 'duplicates_removed' relate, and warns that paging with 'skip' deduplicates each page independently, so consecutive pages may overlap or omit records and 'result_count' must not be summed across pages."
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `FDA_get_active_ingredient_info_by_drug_name`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("fda_get_active_ingredient_info_by_drug_name", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.fda_get_adverse_reactions_by_drug_name` {#dnakit.databases.fda_get_adverse_reactions_by_drug_name}

Retrieve adverse reactions information based on the drug name. IMPORTANT -- FDA labels file safety content under different section names depending on the label's FORMAT, not the drug: modern PLR-format labels (Physician Labeling Rule, 2006 onward) use 'warnings_and_cautions' and usually carry NO 'warnings' section, while legacy and OTC-monograph labels use 'warnings'/'precautions' and usually carry no 'warnings_and_cautions'. Measured against openFDA drug/label on 2026-08-11 (search=_exists_:<field>): 261,639 labels total; 208,073 have 'warnings'; 46,925 have 'warnings_and_cautions'; 43,486 have 'warnings_and_cautions' but no 'warnings' at all; 33,054 have 'boxed_warning'. A null section therefore does NOT mean the drug lacks that content. When every requested section comes back null for a matched label but that same label carries a related section, the related section's content is added to the row under ITS OWN key (never rewritten into the requested key), its name is listed in 'related_sections_present', and a top-level 'section_note' explains the split and names the sibling tool that returns it. Example: GIAPREZA (angiotensin II) has warnings=null and boxed_warning=null but documents its thrombosis risk in 'warnings_and_cautions'. Related sections consulted for this tool: warnings, boxed_warning, precautions, general_precautions (tools: FDA_get_warnings_by_drug_name, FDA_get_boxed_warning_info_by_drug_name, FDA_get_precautions_by_drug_name, FDA_get_general_precautions_by_drug_name).

```text
dnakit.databases.fda_get_adverse_reactions_by_drug_name(*, drug_name, limit=None, skip=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `drug_name` | string | Required by input schema | The name of the drug. |
| `limit` | integer | None | The number of records to return. |
| `skip` | integer | None | The number of records to skip. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "drug_name": {
          "type": "string",
          "description": "The name of the drug."
        },
        "limit": {
          "type": "integer",
          "description": "The number of records to return."
        },
        "skip": {
          "type": "integer",
          "description": "The number of records to skip."
        }
      },
      "required": [
        "drug_name"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "meta": {
          "type": "object",
          "properties": {
            "skip": {
              "type": "integer",
              "description": "Number of upstream records skipped before this page."
            },
            "limit": {
              "type": "integer",
              "description": "Maximum number of records requested from openFDA for this call (page size)."
            },
            "total": {
              "type": "integer",
              "description": "Total number of label records openFDA matched UPSTREAM, before any local processing. This is NOT the length of 'results': local deduplication and the applied limit both shrink the list. Use 'result_count' for the number of rows actually returned and 'duplicates_removed' to reconcile the two."
            }
          }
        },
        "results": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "openfda.brand_name": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                }
              },
              "openfda.generic_name": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                }
              },
              "adverse_reactions": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                }
              },
              "warnings_and_cautions": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                }
              },
              "spl_product_data_elements": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                },
                "description": "Product/ingredient identity string from the SPL record (e.g. 'Ethyol amifostine AMIFOSTINE AMIFOSTINE'). Populated on records whose 'openfda' block is empty and therefore whose openfda.brand_name / openfda.generic_name come back null; use it to confirm which product a row describes when those two fields are unavailable."
              },
              "warnings": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                },
                "description": "Related label section, included only when every requested section came back null for this label and this section is populated instead. Content is verbatim from that section under its own name, so provenance is never misrepresented. See 'section_note' and 'related_sections_present'. Retrievable directly via FDA_get_warnings_by_drug_name."
              },
              "boxed_warning": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                },
                "description": "Related label section, included only when every requested section came back null for this label and this section is populated instead. Content is verbatim from that section under its own name, so provenance is never misrepresented. See 'section_note' and 'related_sections_present'. Retrievable directly via FDA_get_boxed_warning_info_by_drug_name."
              },
              "precautions": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                },
                "description": "Related label section, included only when every requested section came back null for this label and this section is populated instead. Content is verbatim from that section under its own name, so provenance is never misrepresented. See 'section_note' and 'related_sections_present'. Retrievable directly via FDA_get_precautions_by_drug_name."
              },
              "general_precautions": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                },
                "description": "Related label section, included only when every requested section came back null for this label and this section is populated instead. Content is verbatim from that section under its own name, so provenance is never misrepresented. See 'section_note' and 'related_sections_present'. Retrievable directly via FDA_get_general_precautions_by_drug_name."
              },
              "related_sections_present": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                },
                "description": "Names of related label sections that ARE populated on this label. Present only when every requested section came back null for this row. Each named section is also attached to this row under its own key. Its presence means 'the requested section is missing from this label FORMAT', not 'the drug has no such content'."
              }
            }
          }
        },
        "result_count": {
          "type": "integer",
          "description": "Number of rows in 'results' -- i.e. after byte-identical duplicate labels were removed and the caller's limit was applied. Compare against meta.total, which is the upstream pre-deduplication hit count."
        },
        "duplicates_removed": {
          "type": "integer",
          "description": "Number of fetched records dropped because another returned record was byte-identical to them. 0 when nothing was deduplicated."
        },
        "dedup_note": {
          "type": "string",
          "description": "Present only when local deduplication ran. Explains how meta.total, 'result_count' and 'duplicates_removed' relate, and warns that paging with 'skip' deduplicates each page independently, so consecutive pages may overlap or omit records and 'result_count' must not be summed across pages."
        },
        "section_note": {
          "type": "string",
          "description": "Present only when at least one returned label lacks every requested section while carrying a related one. Explains that FDA splits this content by label format (PLR labels use 'warnings_and_cautions'; legacy/OTC labels use 'warnings'/'precautions'), names the sections actually present, and names the sibling tool that returns them. Never read a null section as 'this drug has no warnings' without checking this field."
        },
        "note": {
          "type": "string",
          "description": "Present only when the query had to be broadened to return anything -- closest-spelling candidate, term-based match, or full-label-text match. Describes exactly how the match was weakened; verify the returned openfda.brand_name / openfda.generic_name is the product you intended."
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `FDA_get_adverse_reactions_by_drug_name`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("fda_get_adverse_reactions_by_drug_name", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.fda_get_boxed_warning_info_by_drug_name` {#dnakit.databases.fda_get_boxed_warning_info_by_drug_name}

Retrieve boxed warning and adverse effects information for a specific drug. IMPORTANT -- FDA labels file safety content under different section names depending on the label's FORMAT, not the drug: modern PLR-format labels (Physician Labeling Rule, 2006 onward) use 'warnings_and_cautions' and usually carry NO 'warnings' section, while legacy and OTC-monograph labels use 'warnings'/'precautions' and usually carry no 'warnings_and_cautions'. Measured against openFDA drug/label on 2026-08-11 (search=_exists_:<field>): 261,639 labels total; 208,073 have 'warnings'; 46,925 have 'warnings_and_cautions'; 43,486 have 'warnings_and_cautions' but no 'warnings' at all; 33,054 have 'boxed_warning'. A null section therefore does NOT mean the drug lacks that content. When every requested section comes back null for a matched label but that same label carries a related section, the related section's content is added to the row under ITS OWN key (never rewritten into the requested key), its name is listed in 'related_sections_present', and a top-level 'section_note' explains the split and names the sibling tool that returns it. Example: GIAPREZA (angiotensin II) has warnings=null and boxed_warning=null but documents its thrombosis risk in 'warnings_and_cautions'. Related sections consulted for this tool: warnings_and_cautions, warnings (tools: FDA_get_warnings_and_cautions_by_drug_name, FDA_get_warnings_by_drug_name).

```text
dnakit.databases.fda_get_boxed_warning_info_by_drug_name(*, drug_name, limit=None, skip=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `drug_name` | string | Required by input schema | The name of the drug. |
| `limit` | integer | None | The number of records to return. |
| `skip` | integer | None | The number of records to skip. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "drug_name": {
          "type": "string",
          "description": "The name of the drug."
        },
        "limit": {
          "type": "integer",
          "description": "The number of records to return."
        },
        "skip": {
          "type": "integer",
          "description": "The number of records to skip."
        }
      },
      "required": [
        "drug_name"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "meta": {
          "type": "object",
          "properties": {
            "skip": {
              "type": "integer",
              "description": "Number of upstream records skipped before this page."
            },
            "limit": {
              "type": "integer",
              "description": "Maximum number of records requested from openFDA for this call (page size)."
            },
            "total": {
              "type": "integer",
              "description": "Total number of label records openFDA matched UPSTREAM, before any local processing. This is NOT the length of 'results': local deduplication and the applied limit both shrink the list. Use 'result_count' for the number of rows actually returned and 'duplicates_removed' to reconcile the two."
            }
          }
        },
        "results": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "openfda.brand_name": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                }
              },
              "openfda.generic_name": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                }
              },
              "boxed_warning": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                }
              },
              "spl_product_data_elements": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                },
                "description": "Product/ingredient identity string from the SPL record (e.g. 'Ethyol amifostine AMIFOSTINE AMIFOSTINE'). Populated on records whose 'openfda' block is empty and therefore whose openfda.brand_name / openfda.generic_name come back null; use it to confirm which product a row describes when those two fields are unavailable."
              },
              "warnings_and_cautions": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                },
                "description": "Related label section, included only when every requested section came back null for this label and this section is populated instead. Content is verbatim from that section under its own name, so provenance is never misrepresented. See 'section_note' and 'related_sections_present'. Retrievable directly via FDA_get_warnings_and_cautions_by_drug_name."
              },
              "warnings": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                },
                "description": "Related label section, included only when every requested section came back null for this label and this section is populated instead. Content is verbatim from that section under its own name, so provenance is never misrepresented. See 'section_note' and 'related_sections_present'. Retrievable directly via FDA_get_warnings_by_drug_name."
              },
              "related_sections_present": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                },
                "description": "Names of related label sections that ARE populated on this label. Present only when every requested section came back null for this row. Each named section is also attached to this row under its own key. Its presence means 'the requested section is missing from this label FORMAT', not 'the drug has no such content'."
              }
            }
          }
        },
        "result_count": {
          "type": "integer",
          "description": "Number of rows in 'results' -- i.e. after byte-identical duplicate labels were removed and the caller's limit was applied. Compare against meta.total, which is the upstream pre-deduplication hit count."
        },
        "duplicates_removed": {
          "type": "integer",
          "description": "Number of fetched records dropped because another returned record was byte-identical to them. 0 when nothing was deduplicated."
        },
        "dedup_note": {
          "type": "string",
          "description": "Present only when local deduplication ran. Explains how meta.total, 'result_count' and 'duplicates_removed' relate, and warns that paging with 'skip' deduplicates each page independently, so consecutive pages may overlap or omit records and 'result_count' must not be summed across pages."
        },
        "section_note": {
          "type": "string",
          "description": "Present only when at least one returned label lacks every requested section while carrying a related one. Explains that FDA splits this content by label format (PLR labels use 'warnings_and_cautions'; legacy/OTC labels use 'warnings'/'precautions'), names the sections actually present, and names the sibling tool that returns them. Never read a null section as 'this drug has no warnings' without checking this field."
        },
        "note": {
          "type": "string",
          "description": "Present only when the query had to be broadened to return anything -- closest-spelling candidate, term-based match, or full-label-text match. Describes exactly how the match was weakened; verify the returned openfda.brand_name / openfda.generic_name is the product you intended."
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `FDA_get_boxed_warning_info_by_drug_name`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("fda_get_boxed_warning_info_by_drug_name", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.fda_get_clinical_studies_info_by_drug_name` {#dnakit.databases.fda_get_clinical_studies_info_by_drug_name}

Retrieve clinical studies information based on the drug name.

```text
dnakit.databases.fda_get_clinical_studies_info_by_drug_name(*, drug_name, limit=None, skip=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `drug_name` | string | Required by input schema | The name of the drug. |
| `limit` | integer | None | The number of records to return. |
| `skip` | integer | None | The number of records to skip. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "drug_name": {
          "type": "string",
          "description": "The name of the drug."
        },
        "limit": {
          "type": "integer",
          "description": "The number of records to return."
        },
        "skip": {
          "type": "integer",
          "description": "The number of records to skip."
        }
      },
      "required": [
        "drug_name"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "meta": {
          "type": "object",
          "properties": {
            "skip": {
              "type": "integer",
              "description": "Number of upstream records skipped before this page."
            },
            "limit": {
              "type": "integer",
              "description": "Maximum number of records requested from openFDA for this call (page size)."
            },
            "total": {
              "type": "integer",
              "description": "Total number of label records openFDA matched UPSTREAM, before any local processing. This is NOT the length of 'results': local deduplication and the applied limit both shrink the list. Use 'result_count' for the number of rows actually returned and 'duplicates_removed' to reconcile the two."
            }
          }
        },
        "results": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "openfda.brand_name": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                }
              },
              "openfda.generic_name": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                }
              },
              "clinical_studies": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                }
              },
              "spl_product_data_elements": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                },
                "description": "Product/ingredient identity string from the SPL record (e.g. 'Ethyol amifostine AMIFOSTINE AMIFOSTINE'). Populated on records whose 'openfda' block is empty and therefore whose openfda.brand_name / openfda.generic_name come back null; use it to confirm which product a row describes when those two fields are unavailable."
              }
            }
          }
        },
        "result_count": {
          "type": "integer",
          "description": "Number of rows in 'results' -- i.e. after byte-identical duplicate labels were removed and the caller's limit was applied. Compare against meta.total, which is the upstream pre-deduplication hit count."
        },
        "duplicates_removed": {
          "type": "integer",
          "description": "Number of fetched records dropped because another returned record was byte-identical to them. 0 when nothing was deduplicated."
        },
        "dedup_note": {
          "type": "string",
          "description": "Present only when local deduplication ran. Explains how meta.total, 'result_count' and 'duplicates_removed' relate, and warns that paging with 'skip' deduplicates each page independently, so consecutive pages may overlap or omit records and 'result_count' must not be summed across pages."
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `FDA_get_clinical_studies_info_by_drug_name`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("fda_get_clinical_studies_info_by_drug_name", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.fda_get_indications_by_drug_name` {#dnakit.databases.fda_get_indications_by_drug_name}

Retrieve indications and usage information based on a specific drug name.

```text
dnakit.databases.fda_get_indications_by_drug_name(*, drug_name, limit=None, skip=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `drug_name` | string | Required by input schema | The name of the drug. |
| `limit` | integer | None | The number of records to return. |
| `skip` | integer | None | The number of records to skip. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "drug_name": {
          "type": "string",
          "description": "The name of the drug."
        },
        "limit": {
          "type": "integer",
          "description": "The number of records to return."
        },
        "skip": {
          "type": "integer",
          "description": "The number of records to skip."
        }
      },
      "required": [
        "drug_name"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "meta": {
          "type": "object",
          "properties": {
            "skip": {
              "type": "integer",
              "description": "Number of upstream records skipped before this page."
            },
            "limit": {
              "type": "integer",
              "description": "Maximum number of records requested from openFDA for this call (page size)."
            },
            "total": {
              "type": "integer",
              "description": "Total number of label records openFDA matched UPSTREAM, before any local processing. This is NOT the length of 'results': local deduplication and the applied limit both shrink the list. Use 'result_count' for the number of rows actually returned and 'duplicates_removed' to reconcile the two."
            }
          }
        },
        "results": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "openfda.brand_name": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                }
              },
              "openfda.generic_name": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                }
              },
              "indications_and_usage": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                }
              },
              "spl_product_data_elements": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                },
                "description": "Product/ingredient identity string from the SPL record (e.g. 'Ethyol amifostine AMIFOSTINE AMIFOSTINE'). Populated on records whose 'openfda' block is empty and therefore whose openfda.brand_name / openfda.generic_name come back null; use it to confirm which product a row describes when those two fields are unavailable."
              }
            }
          }
        },
        "result_count": {
          "type": "integer",
          "description": "Number of rows in 'results' -- i.e. after byte-identical duplicate labels were removed and the caller's limit was applied. Compare against meta.total, which is the upstream pre-deduplication hit count."
        },
        "duplicates_removed": {
          "type": "integer",
          "description": "Number of fetched records dropped because another returned record was byte-identical to them. 0 when nothing was deduplicated."
        },
        "dedup_note": {
          "type": "string",
          "description": "Present only when local deduplication ran. Explains how meta.total, 'result_count' and 'duplicates_removed' relate, and warns that paging with 'skip' deduplicates each page independently, so consecutive pages may overlap or omit records and 'result_count' must not be summed across pages."
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `FDA_get_indications_by_drug_name`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("fda_get_indications_by_drug_name", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.fda_get_mechanism_of_action_by_drug_name` {#dnakit.databases.fda_get_mechanism_of_action_by_drug_name}

Retrieve the mechanism of action information for a specific drug.

```text
dnakit.databases.fda_get_mechanism_of_action_by_drug_name(*, drug_name, limit=None, skip=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `drug_name` | string | Required by input schema | The name of the drug. |
| `limit` | integer | None | The number of records to return. |
| `skip` | integer | None | The number of records to skip. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "drug_name": {
          "type": "string",
          "description": "The name of the drug."
        },
        "limit": {
          "type": "integer",
          "description": "The number of records to return."
        },
        "skip": {
          "type": "integer",
          "description": "The number of records to skip."
        }
      },
      "required": [
        "drug_name"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "meta": {
          "type": "object",
          "properties": {
            "skip": {
              "type": "integer",
              "description": "Number of upstream records skipped before this page."
            },
            "limit": {
              "type": "integer",
              "description": "Maximum number of records requested from openFDA for this call (page size)."
            },
            "total": {
              "type": "integer",
              "description": "Total number of label records openFDA matched UPSTREAM, before any local processing. This is NOT the length of 'results': local deduplication and the applied limit both shrink the list. Use 'result_count' for the number of rows actually returned and 'duplicates_removed' to reconcile the two."
            }
          }
        },
        "results": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "openfda.brand_name": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                }
              },
              "openfda.generic_name": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                }
              },
              "mechanism_of_action": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                }
              },
              "spl_product_data_elements": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                },
                "description": "Product/ingredient identity string from the SPL record (e.g. 'Ethyol amifostine AMIFOSTINE AMIFOSTINE'). Populated on records whose 'openfda' block is empty and therefore whose openfda.brand_name / openfda.generic_name come back null; use it to confirm which product a row describes when those two fields are unavailable."
              }
            }
          }
        },
        "result_count": {
          "type": "integer",
          "description": "Number of rows in 'results' -- i.e. after byte-identical duplicate labels were removed and the caller's limit was applied. Compare against meta.total, which is the upstream pre-deduplication hit count."
        },
        "duplicates_removed": {
          "type": "integer",
          "description": "Number of fetched records dropped because another returned record was byte-identical to them. 0 when nothing was deduplicated."
        },
        "dedup_note": {
          "type": "string",
          "description": "Present only when local deduplication ran. Explains how meta.total, 'result_count' and 'duplicates_removed' relate, and warns that paging with 'skip' deduplicates each page independently, so consecutive pages may overlap or omit records and 'result_count' must not be summed across pages."
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `FDA_get_mechanism_of_action_by_drug_name`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("fda_get_mechanism_of_action_by_drug_name", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.hmdb_get_diseases` {#dnakit.databases.hmdb_get_diseases}

Alias for Metabolite_get_diseases. Get disease associations for a metabolite via CTD database.

```text
dnakit.databases.hmdb_get_diseases(*, operation=None, hmdb_id=None, compound_name=None, pubchem_cid=None, limit=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Provider operation mode; see this row&#x27;s default and allowed values. {&quot;const&quot;: &quot;get_diseases&quot;, &quot;default&quot;: &quot;get_diseases&quot;} |
| `hmdb_id` | string | None | HMDB ID (e.g., HMDB0000122) |
| `compound_name` | string | None | Compound name (e.g., glucose) |
| `pubchem_cid` | string | None | PubChem CID as an alternative identifier to hmdb_id or compound_name. |
| `limit` | integer | None | Maximum number of disease associations to return (default 50). {&quot;default&quot;: 50} |
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
          "const": "get_diseases",
          "default": "get_diseases"
        },
        "hmdb_id": {
          "type": "string",
          "description": "HMDB ID (e.g., HMDB0000122)"
        },
        "compound_name": {
          "type": "string",
          "description": "Compound name (e.g., glucose)"
        },
        "pubchem_cid": {
          "type": "string",
          "description": "PubChem CID as an alternative identifier to hmdb_id or compound_name."
        },
        "limit": {
          "type": "integer",
          "default": 50,
          "description": "Maximum number of disease associations to return (default 50)."
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
            "identifier": {
              "type": "string",
              "description": "The identifier that was submitted"
            },
            "compound_name": {
              "type": "string",
              "description": "Canonical compound name resolved via PubChem"
            },
            "pubchem_cid": {
              "type": [
                "integer",
                "null"
              ]
            },
            "ctd_query_term": {
              "type": "string",
              "description": "The name variant used to query CTD"
            },
            "disease_count": {
              "type": "integer"
            },
            "diseases": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "disease_name": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "disease_id": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "disease_categories": {
                    "type": [
                      "string",
                      "array",
                      "null"
                    ]
                  },
                  "direct_evidence": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "pubmed_ids": {
                    "type": "array",
                    "items": {
                      "type": "string"
                    }
                  }
                }
              }
            }
          },
          "required": [
            "identifier",
            "compound_name",
            "disease_count",
            "diseases"
          ]
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

**Execution and provenance:** `network`；upstream name `HMDB_get_diseases`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("hmdb_get_diseases", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.hmdb_get_metabolite` {#dnakit.databases.hmdb_get_metabolite}

Alias for Metabolite_get_info. Get compound information by HMDB ID, compound name, or PubChem CID.

```text
dnakit.databases.hmdb_get_metabolite(*, operation=None, hmdb_id=None, compound_name=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Provider operation mode; see this row&#x27;s default and allowed values. {&quot;const&quot;: &quot;get_info&quot;, &quot;default&quot;: &quot;get_info&quot;} |
| `hmdb_id` | string | None | HMDB ID (e.g., HMDB0000122) |
| `compound_name` | string | None | Compound name (e.g., glucose) |
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
          "const": "get_info",
          "default": "get_info"
        },
        "hmdb_id": {
          "type": "string",
          "description": "HMDB ID (e.g., HMDB0000122)"
        },
        "compound_name": {
          "type": "string",
          "description": "Compound name (e.g., glucose)"
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
            "pubchem_cid": {
              "type": [
                "integer",
                "null"
              ]
            },
            "name": {
              "type": [
                "string",
                "null"
              ]
            },
            "iupac_name": {
              "type": [
                "string",
                "null"
              ]
            },
            "formula": {
              "type": [
                "string",
                "null"
              ]
            },
            "molecular_weight": {
              "type": [
                "string",
                "null"
              ]
            },
            "smiles": {
              "type": [
                "string",
                "null"
              ]
            },
            "inchikey": {
              "type": [
                "string",
                "null"
              ]
            }
          },
          "required": [
            "pubchem_cid"
          ]
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

**Execution and provenance:** `network`；upstream name `HMDB_get_metabolite`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("hmdb_get_metabolite", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.hmdb_search` {#dnakit.databases.hmdb_search}

Alias for Metabolite_search. Search for metabolites by compound name or molecular formula via PubChem. Accepts common names, synonyms (e.g., 'glucosylceramide', 'glucose', 'ATP').

```text
dnakit.databases.hmdb_search(*, operation=None, query, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Provider operation mode; see this row&#x27;s default and allowed values. {&quot;const&quot;: &quot;search&quot;, &quot;default&quot;: &quot;search&quot;} |
| `query` | string | Required by input schema | Compound name or molecular formula to search |
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
          "default": "search"
        },
        "query": {
          "type": "string",
          "description": "Compound name or molecular formula to search"
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
          "type": "object",
          "properties": {
            "query": {
              "type": "string",
              "description": "The search query that was submitted"
            },
            "results": {
              "type": "array",
              "description": "Matching compounds from PubChem",
              "items": {
                "type": "object",
                "properties": {
                  "pubchem_cid": {
                    "type": [
                      "integer",
                      "null"
                    ]
                  },
                  "name": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "formula": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "molecular_weight": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "smiles": {
                    "type": [
                      "string",
                      "null"
                    ]
                  }
                }
              }
            },
            "count": {
              "type": "integer",
              "description": "Number of compounds returned"
            }
          },
          "required": [
            "query",
            "results",
            "count"
          ]
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

**Execution and provenance:** `network`；upstream name `HMDB_search`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("hmdb_search", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_associated_drugs_by_disease_efo_id` {#dnakit.databases.open_targets_get_associated_drugs_by_disease_efo_id}

Retrieve known drugs associated with a specific disease by disease efoId.

```text
dnakit.databases.open_targets_get_associated_drugs_by_disease_efo_id(*, efoId, size=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | Required by input schema | The EFO ID of the disease. |
| `size` | integer | None | Limit for the returned page; does not automatically retrieve later pages. {&quot;minimum&quot;: 1, &quot;maximum&quot;: 1000} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "efoId": {
          "type": "string",
          "description": "The EFO ID of the disease."
        },
        "size": {
          "type": "integer",
          "minimum": 1,
          "maximum": 1000
        }
      },
      "required": [
        "efoId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "disease": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "name": {
                  "type": "string"
                },
                "drugAndClinicalCandidates": {
                  "type": "object",
                  "properties": {
                    "count": {
                      "type": "integer"
                    },
                    "returned": {
                      "type": "integer",
                      "description": "Number of drugs and clinical candidates actually present in `rows` for this page."
                    },
                    "truncated": {
                      "type": "boolean",
                      "description": "True when `returned` is less than `count`, i.e. rows remain unfetched; `metadata.truncation_note` then names the exact argument that fetches the rest."
                    },
                    "rows": {
                      "type": "array",
                      "items": {
                        "type": "object",
                        "properties": {
                          "drug": {
                            "type": "object",
                            "properties": {
                              "id": {
                                "type": "string"
                              },
                              "name": {
                                "type": "string"
                              },
                              "tradeNames": {
                                "type": "array",
                                "items": {
                                  "type": "string"
                                }
                              },
                              "maximumClinicalTrialPhase": {
                                "type": "integer"
                              },
                              "isApproved": {
                                "type": "boolean"
                              },
                              "hasBeenWithdrawn": {
                                "type": "boolean"
                              }
                            }
                          },
                          "maxClinicalStage": {
                            "type": [
                              "integer",
                              "string",
                              "null"
                            ]
                          },
                          "mechanismOfAction": {
                            "type": "string"
                          },
                          "target": {
                            "type": "object",
                            "properties": {
                              "id": {
                                "type": "string"
                              },
                              "approvedSymbol": {
                                "type": "string"
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        },
        "metadata": {
          "type": "object",
          "properties": {
            "truncation_note": {
              "type": "string",
              "description": "Present only when a collection came back truncated: states how many of how many rows were returned and names the exact argument that fetches the rest."
            }
          }
        }
      }
    }
    ```

??? info "DNAKit argument extension rules"

    ```json
    {
      "parameters": {
        "size": {
          "type": "integer",
          "minimum": 1,
          "maximum": 1000
        }
      },
      "handler": "page_size"
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_associated_drugs_by_disease_efoId`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_associated_drugs_by_disease_efo_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_associated_drugs_by_target_ensembl_id` {#dnakit.databases.open_targets_get_associated_drugs_by_target_ensembl_id}

Get known drugs and information (e.g. id, name, MoA) associated with a specific target ensemblID, including clinical trial phase and mechanism of action of the drugs.

```text
dnakit.databases.open_targets_get_associated_drugs_by_target_ensembl_id(*, ensemblId, size=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `ensemblId` | string | Required by input schema | The Ensembl ID of the target. |
| `size` | integer | None | Limit for the returned page; does not automatically retrieve later pages. {&quot;minimum&quot;: 1, &quot;maximum&quot;: 1000} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "ensemblId": {
          "type": "string",
          "description": "The Ensembl ID of the target."
        },
        "size": {
          "type": "integer",
          "minimum": 1,
          "maximum": 1000
        }
      },
      "required": [
        "ensemblId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "target": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "approvedSymbol": {
                  "type": "string"
                },
                "drugAndClinicalCandidates": {
                  "type": "object",
                  "properties": {
                    "count": {
                      "type": "integer"
                    },
                    "returned": {
                      "type": "integer",
                      "description": "Number of drugs and clinical candidates actually present in `rows` for this page."
                    },
                    "truncated": {
                      "type": "boolean",
                      "description": "True when `returned` is less than `count`, i.e. rows remain unfetched; `metadata.truncation_note` then names the exact argument that fetches the rest."
                    },
                    "rows": {
                      "type": "array",
                      "items": {
                        "type": "object",
                        "properties": {
                          "drug": {
                            "type": "object",
                            "properties": {
                              "id": {
                                "type": "string"
                              },
                              "name": {
                                "type": "string"
                              },
                              "maximumClinicalTrialPhase": {
                                "type": "integer"
                              },
                              "isApproved": {
                                "type": "boolean"
                              },
                              "hasBeenWithdrawn": {
                                "type": "boolean"
                              }
                            }
                          },
                          "maxClinicalStage": {
                            "type": "string"
                          },
                          "mechanismOfAction": {
                            "type": "string"
                          },
                          "disease": {
                            "type": "object",
                            "properties": {
                              "id": {
                                "type": "string"
                              },
                              "name": {
                                "type": "string"
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        },
        "metadata": {
          "type": "object",
          "properties": {
            "truncation_note": {
              "type": "string",
              "description": "Present only when a collection came back truncated: states how many of how many rows were returned and names the exact argument that fetches the rest."
            }
          }
        }
      }
    }
    ```

??? info "DNAKit argument extension rules"

    ```json
    {
      "parameters": {
        "size": {
          "type": "integer",
          "minimum": 1,
          "maximum": 1000
        }
      },
      "handler": "page_size"
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_associated_drugs_by_target_ensemblID`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_associated_drugs_by_target_ensembl_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_associated_phenotypes_by_disease_efo_id` {#dnakit.databases.open_targets_get_associated_phenotypes_by_disease_efo_id}

Find HPO phenotypes asosciated with the specified disease efoId.

```text
dnakit.databases.open_targets_get_associated_phenotypes_by_disease_efo_id(*, efoId, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | Required by input schema | The efoId of a disease or phenotype. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "efoId": {
          "type": "string",
          "description": "The efoId of a disease or phenotype."
        }
      },
      "required": [
        "efoId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "disease": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "name": {
                  "type": "string"
                },
                "phenotypes": {
                  "type": "object",
                  "properties": {
                    "rows": {
                      "type": "array",
                      "items": {
                        "type": "object",
                        "properties": {
                          "phenotypeHPO": {
                            "type": "object",
                            "properties": {
                              "id": {
                                "type": "string"
                              },
                              "name": {
                                "type": "string"
                              },
                              "description": {
                                "type": "string"
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_associated_phenotypes_by_disease_efoId`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_associated_phenotypes_by_disease_efo_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_associated_targets_by_disease_efo_id` {#dnakit.databases.open_targets_get_associated_targets_by_disease_efo_id}

Find targets associated with a specific disease or phenotype based on efoId. Results are PAGINATED (50 per call by default): `associatedTargets.count` is the total number of associated targets, `associatedTargets.returned` is how many are in `rows`, and `associatedTargets.truncated` is true when targets remain unfetched, in which case `metadata.truncation_note` names the `size`/`index` values to pass. Raise `size` (Open Targets rejects a page above 3000) and step `index` to page through; a common disease has thousands of associated targets, so the first page is NOT the full list.

```text
dnakit.databases.open_targets_get_associated_targets_by_disease_efo_id(*, efoId, size=None, index=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | Required by input schema | The efoId of a disease or phenotype. |
| `size` | integer | None | Number of top-scored associated targets to return (default 50). A disease can have thousands of associated targets; raise this to retrieve more than the first page. {&quot;default&quot;: 50} |
| `index` | integer | None | Zero-based page index for paginating beyond the first `size` targets (default 0). |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "efoId": {
          "type": "string",
          "description": "The efoId of a disease or phenotype."
        },
        "size": {
          "type": "integer",
          "description": "Number of top-scored associated targets to return (default 50). A disease can have thousands of associated targets; raise this to retrieve more than the first page.",
          "default": 50
        },
        "index": {
          "type": "integer",
          "description": "Zero-based page index for paginating beyond the first `size` targets (default 0)."
        }
      },
      "required": [
        "efoId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "disease": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "name": {
                  "type": "string"
                },
                "associatedTargets": {
                  "type": "object",
                  "properties": {
                    "count": {
                      "type": "integer"
                    },
                    "returned": {
                      "type": "integer",
                      "description": "Number of associated targets actually present in `rows` for this page."
                    },
                    "truncated": {
                      "type": "boolean",
                      "description": "True when `returned` is less than `count`, i.e. targets remain unfetched; `metadata.truncation_note` then names the `size`/`index` values to pass."
                    },
                    "rows": {
                      "type": "array",
                      "items": {
                        "type": "object",
                        "properties": {
                          "target": {
                            "type": "object",
                            "properties": {
                              "id": {
                                "type": "string"
                              },
                              "approvedSymbol": {
                                "type": "string"
                              }
                            }
                          },
                          "score": {
                            "type": "number"
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        },
        "metadata": {
          "type": "object",
          "properties": {
            "truncation_note": {
              "type": "string",
              "description": "Present only when a collection came back truncated: states how many of how many rows were returned and names the exact argument that fetches the rest."
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_associated_targets_by_disease_efoId`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_associated_targets_by_disease_efo_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_credible_set_detail` {#dnakit.databases.open_targets_get_credible_set_detail}

Get detailed information about a specific GWAS credible set (fine-mapped locus) from Open Targets. Returns the lead variant, associated study, fine-mapping method, statistical significance, locus-to-gene (L2G) predictions mapping the signal to likely causal genes, and colocalization evidence with QTLs.

```text
dnakit.databases.open_targets_get_credible_set_detail(*, studyLocusId, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `studyLocusId` | string | Required by input schema | Credible set identifier (studyLocusId), a 32-character hash (e.g., &#x27;b758d8fb10924f5338cbad8d27c7dee8&#x27;). Obtain from variant credible sets or study credible sets queries. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "studyLocusId": {
          "type": "string",
          "description": "Credible set identifier (studyLocusId), a 32-character hash (e.g., 'b758d8fb10924f5338cbad8d27c7dee8'). Obtain from variant credible sets or study credible sets queries."
        }
      },
      "required": [
        "studyLocusId"
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
            "credibleSet": {
              "type": "object",
              "properties": {
                "studyLocusId": {
                  "type": "string"
                },
                "studyId": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "region": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "chromosome": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "position": {
                  "type": [
                    "integer",
                    "null"
                  ]
                },
                "finemappingMethod": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "pValueMantissa": {
                  "type": [
                    "number",
                    "null"
                  ]
                },
                "pValueExponent": {
                  "type": [
                    "integer",
                    "null"
                  ]
                },
                "beta": {
                  "type": [
                    "number",
                    "null"
                  ]
                },
                "standardError": {
                  "type": [
                    "number",
                    "null"
                  ]
                },
                "confidence": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "credibleSetlog10BF": {
                  "type": [
                    "number",
                    "null"
                  ]
                },
                "sampleSize": {
                  "type": [
                    "integer",
                    "null"
                  ]
                },
                "variant": {
                  "type": "object",
                  "properties": {
                    "id": {
                      "type": "string"
                    },
                    "rsIds": {
                      "type": "array",
                      "items": {
                        "type": "string"
                      }
                    },
                    "chromosome": {
                      "type": "string"
                    },
                    "position": {
                      "type": "integer"
                    },
                    "referenceAllele": {
                      "type": "string"
                    },
                    "alternateAllele": {
                      "type": "string"
                    },
                    "mostSevereConsequence": {
                      "type": "object",
                      "properties": {
                        "id": {
                          "type": "string"
                        },
                        "label": {
                          "type": "string"
                        }
                      }
                    }
                  }
                },
                "study": {
                  "type": "object",
                  "properties": {
                    "id": {
                      "type": "string"
                    },
                    "traitFromSource": {
                      "type": [
                        "string",
                        "null"
                      ]
                    },
                    "publicationFirstAuthor": {
                      "type": [
                        "string",
                        "null"
                      ]
                    },
                    "publicationTitle": {
                      "type": [
                        "string",
                        "null"
                      ]
                    },
                    "pubmedId": {
                      "type": [
                        "string",
                        "null"
                      ]
                    },
                    "nSamples": {
                      "type": [
                        "integer",
                        "null"
                      ]
                    },
                    "diseases": {
                      "type": "array",
                      "items": {
                        "type": "object",
                        "properties": {
                          "id": {
                            "type": "string"
                          },
                          "name": {
                            "type": "string"
                          }
                        }
                      }
                    }
                  }
                },
                "l2GPredictions": {
                  "type": "object",
                  "properties": {
                    "rows": {
                      "type": "array",
                      "items": {
                        "type": "object",
                        "properties": {
                          "target": {
                            "type": "object",
                            "properties": {
                              "id": {
                                "type": "string"
                              },
                              "approvedSymbol": {
                                "type": "string"
                              }
                            }
                          },
                          "score": {
                            "type": "number"
                          }
                        }
                      }
                    }
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
              "type": "object"
            }
          },
          "required": [
            "error"
          ]
        }
      ]
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_credible_set_detail`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_credible_set_detail", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_disease_ancestors_parents_by_efo_id` {#dnakit.databases.open_targets_get_disease_ancestors_parents_by_efo_id}

Retrieve the disease ancestors and parents in the ontology using the disease EFO ID.

```text
dnakit.databases.open_targets_get_disease_ancestors_parents_by_efo_id(*, efoId, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | Required by input schema | The EFO ID of the disease. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "efoId": {
          "type": "string",
          "description": "The EFO ID of the disease."
        }
      },
      "required": [
        "efoId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "disease": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "name": {
                  "type": "string"
                },
                "ancestors": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  }
                },
                "parents": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "id": {
                        "type": "string"
                      },
                      "name": {
                        "type": "string"
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_disease_ancestors_parents_by_efoId`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_disease_ancestors_parents_by_efo_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_disease_descendants_children_by_efo_id` {#dnakit.databases.open_targets_get_disease_descendants_children_by_efo_id}

Retrieve the disease descendants and children in the ontology using the disease EFO ID.

```text
dnakit.databases.open_targets_get_disease_descendants_children_by_efo_id(*, efoId, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | Required by input schema | The EFO ID of the disease. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "efoId": {
          "type": "string",
          "description": "The EFO ID of the disease."
        }
      },
      "required": [
        "efoId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "disease": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "name": {
                  "type": "string"
                },
                "descendants": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  }
                },
                "children": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "id": {
                        "type": "string"
                      },
                      "name": {
                        "type": "string"
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_disease_descendants_children_by_efoId`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_disease_descendants_children_by_efo_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_disease_description_by_efo_id` {#dnakit.databases.open_targets_get_disease_description_by_efo_id}

Retrieve disease description, name, database cros references, obsolete terms, and whether it's a therapeutic area, all using the specified efoId.

```text
dnakit.databases.open_targets_get_disease_description_by_efo_id(*, efoId, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | Required by input schema | The EFO ID of the disease. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "efoId": {
          "type": "string",
          "description": "The EFO ID of the disease."
        }
      },
      "required": [
        "efoId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "disease": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "name": {
                  "type": "string"
                },
                "description": {
                  "type": "string"
                },
                "dbXRefs": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  }
                },
                "isTherapeuticArea": {
                  "type": "boolean"
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_disease_description_by_efoId`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_disease_description_by_efo_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_disease_id_description_by_name` {#dnakit.databases.open_targets_get_disease_id_description_by_name}

Retrieve the efoId and additional details of a disease based on its name.

```text
dnakit.databases.open_targets_get_disease_id_description_by_name(*, diseaseName=None, query=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `diseaseName` | string | None | The name of the disease to search for. |
| `query` | string | None | DNAKit optional alias for diseaseName |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "diseaseName": {
          "type": "string",
          "description": "The name of the disease to search for."
        },
        "query": {
          "type": "string",
          "description": "DNAKit optional alias for diseaseName"
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
        "data": {
          "type": "object",
          "properties": {
            "search": {
              "type": "object",
              "properties": {
                "hits": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "id": {
                        "type": "string"
                      },
                      "name": {
                        "type": "string"
                      },
                      "description": {
                        "type": "string"
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

??? info "DNAKit argument extension rules"

    ```json
    {
      "aliases": {
        "query": "diseaseName"
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_disease_id_description_by_name`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_disease_id_description_by_name", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_disease_synonyms_by_efo_id` {#dnakit.databases.open_targets_get_disease_synonyms_by_efo_id}

Retrieve disease synonyms by its EFO ID.

```text
dnakit.databases.open_targets_get_disease_synonyms_by_efo_id(*, efoId, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | Required by input schema | The EFO ID of the disease. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "efoId": {
          "type": "string",
          "description": "The EFO ID of the disease."
        }
      },
      "required": [
        "efoId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "disease": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "name": {
                  "type": "string"
                },
                "synonyms": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "relation": {
                        "type": "string"
                      },
                      "terms": {
                        "type": "array",
                        "items": {
                          "type": "string"
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_disease_synonyms_by_efoId`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_disease_synonyms_by_efo_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_disease_therapeutic_areas_by_efo_id` {#dnakit.databases.open_targets_get_disease_therapeutic_areas_by_efo_id}

Retrieve the therapeutic areas associated with a specific disease efoId.

```text
dnakit.databases.open_targets_get_disease_therapeutic_areas_by_efo_id(*, efoId, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | Required by input schema | The EFO ID of the disease. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "efoId": {
          "type": "string",
          "description": "The EFO ID of the disease."
        }
      },
      "required": [
        "efoId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "disease": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "name": {
                  "type": "string"
                },
                "therapeuticAreas": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "id": {
                        "type": "string"
                      },
                      "name": {
                        "type": "string"
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_disease_therapeutic_areas_by_efoId`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_disease_therapeutic_areas_by_efo_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_diseases_phenotypes_by_target_ensembl` {#dnakit.databases.open_targets_get_diseases_phenotypes_by_target_ensembl}

Find diseases or phenotypes associated with a specific target using ensemblId. `associatedDiseases.count` is the total number of associations; `rows` returns only one page of them (25 by default) -- pass `page` to retrieve more (e.g. page 2 of a target with thousands of associations like SIRT1's 3432). `associatedDiseases.returned` reports how many rows this response actually contains and `associatedDiseases.truncated` is true when associations remain unfetched, in which case `metadata.truncation_note` names the exact `page` argument to pass.

```text
dnakit.databases.open_targets_get_diseases_phenotypes_by_target_ensembl(*, ensemblId, page=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `ensemblId` | string | Required by input schema | The ensemblId of a target. |
| `page` | object | None | Pagination parameters. `count` in the response is the total across all pages; omit this to get the first page only. {&quot;properties&quot;: {&quot;index&quot;: {&quot;type&quot;: &quot;integer&quot;, &quot;description&quot;: &quot;The index of the page to retrieve (0-based).&quot;}, &quot;size&quot;: {&quot;type&quot;: &quot;integer&quot;, &quot;description&quot;: &quot;The number of associations per page.&quot;}}} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "ensemblId": {
          "type": "string",
          "description": "The ensemblId of a target."
        },
        "page": {
          "type": "object",
          "properties": {
            "index": {
              "type": "integer",
              "description": "The index of the page to retrieve (0-based)."
            },
            "size": {
              "type": "integer",
              "description": "The number of associations per page."
            }
          },
          "description": "Pagination parameters. `count` in the response is the total across all pages; omit this to get the first page only."
        }
      },
      "required": [
        "ensemblId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "target": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "approvedSymbol": {
                  "type": "string"
                },
                "associatedDiseases": {
                  "type": "object",
                  "properties": {
                    "count": {
                      "type": "integer"
                    },
                    "returned": {
                      "type": "integer",
                      "description": "Number of associations actually present in `rows` for this page."
                    },
                    "truncated": {
                      "type": "boolean",
                      "description": "True when `returned` is less than `count`, i.e. associations remain unfetched; `metadata.truncation_note` then names the exact `page` argument to pass."
                    },
                    "rows": {
                      "type": "array",
                      "items": {
                        "type": "object",
                        "properties": {
                          "disease": {
                            "type": "object",
                            "properties": {
                              "id": {
                                "type": "string"
                              },
                              "name": {
                                "type": "string"
                              }
                            }
                          },
                          "datasourceScores": {
                            "type": "array",
                            "items": {
                              "type": "object",
                              "properties": {
                                "id": {
                                  "type": "string"
                                },
                                "score": {
                                  "type": "number"
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        },
        "metadata": {
          "type": "object",
          "properties": {
            "truncation_note": {
              "type": "string",
              "description": "Present only when a collection came back truncated: states how many of how many rows were returned and names the exact argument that fetches the rest."
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_diseases_phenotypes_by_target_ensembl`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_diseases_phenotypes_by_target_ensembl", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_drug_adverse_events_by_chembl_id` {#dnakit.databases.open_targets_get_drug_adverse_events_by_chembl_id}

Retrieve significant adverse events reported for a specific drug chemblId, ranked by log-likelihood ratio. Results are PAGINATED and Open Targets returns only 25 rows when no `page` is given: `adverseEvents.count` is the total number of signals, `adverseEvents.returned` is how many are present in `rows`, and `adverseEvents.truncated` is true when signals remain unfetched, in which case a `truncation_note` naming the exact `page` argument to pass is added to `metadata`. Pass page: {"index": 0, "size": <count>} to retrieve every signal before drawing any safety conclusion — a default page is NOT the full list (ibuprofen, CHEMBL521, has 55 signals; the default 25 omit toxic epidermal necrolysis, urticaria, systemic lupus erythematosus and hepatic enzyme increased).

```text
dnakit.databases.open_targets_get_drug_adverse_events_by_chembl_id(*, chemblId, page=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `chemblId` | string | Required by input schema | The ChEMBL ID of the drug. |
| `page` | object | None | Pagination settings. {&quot;properties&quot;: {&quot;index&quot;: {&quot;type&quot;: &quot;integer&quot;, &quot;description&quot;: &quot;Index of the page to fetch, starting from 0.&quot;}, &quot;size&quot;: {&quot;type&quot;: &quot;integer&quot;, &quot;description&quot;: &quot;Number of entries per page.&quot;}}} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "chemblId": {
          "type": "string",
          "description": "The ChEMBL ID of the drug."
        },
        "page": {
          "type": "object",
          "properties": {
            "index": {
              "type": "integer",
              "description": "Index of the page to fetch, starting from 0."
            },
            "size": {
              "type": "integer",
              "description": "Number of entries per page."
            }
          },
          "description": "Pagination settings."
        }
      },
      "required": [
        "chemblId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "drug": {
          "type": "object",
          "properties": {
            "id": {
              "type": "string"
            },
            "name": {
              "type": "string"
            },
            "adverseEvents": {
              "type": "object",
              "properties": {
                "count": {
                  "type": "integer",
                  "description": "Total number of significant adverse-event signals for this drug across all pages — NOT the number of rows in this response."
                },
                "returned": {
                  "type": "integer",
                  "description": "Number of signals actually present in `rows` for this page."
                },
                "truncated": {
                  "type": "boolean",
                  "description": "True when `returned` is less than `count`, i.e. signals remain unfetched; `metadata.truncation_note` then names the exact `page` argument to pass."
                },
                "criticalValue": {
                  "type": "number"
                },
                "rows": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "name": {
                        "type": "string"
                      },
                      "meddraCode": {
                        "type": "string"
                      },
                      "count": {
                        "type": "integer"
                      },
                      "logLR": {
                        "type": "number"
                      }
                    }
                  }
                }
              }
            }
          }
        },
        "metadata": {
          "type": "object",
          "properties": {
            "truncation_note": {
              "type": "string",
              "description": "Present only when a collection came back truncated: states how many of how many rows were returned and names the exact `page` argument that fetches the rest."
            }
          }
        },
        "error": {
          "type": "string"
        },
        "error_details": {
          "type": "object",
          "properties": {
            "type": {
              "type": "string"
            },
            "message": {
              "type": "string"
            },
            "retriable": {
              "type": "boolean"
            },
            "next_steps": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "details": {
              "type": "object",
              "properties": {
                "validation_error": {
                  "type": "string"
                },
                "path": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  }
                },
                "schema": {
                  "type": "object",
                  "properties": {
                    "type": {
                      "type": "string"
                    },
                    "properties": {
                      "type": "object",
                      "properties": {
                        "chemblId": {
                          "type": "object",
                          "properties": {
                            "type": {
                              "type": "string"
                            },
                            "description": {
                              "type": "string"
                            }
                          }
                        },
                        "page": {
                          "type": "object",
                          "properties": {
                            "type": {
                              "type": "string"
                            },
                            "properties": {
                              "type": "object",
                              "properties": {
                                "index": {
                                  "type": "object",
                                  "properties": {
                                    "type": {
                                      "type": "string"
                                    },
                                    "description": {
                                      "type": "string"
                                    }
                                  }
                                },
                                "size": {
                                  "type": "object",
                                  "properties": {
                                    "type": {
                                      "type": "string"
                                    },
                                    "description": {
                                      "type": "string"
                                    }
                                  }
                                }
                              }
                            },
                            "description": {
                              "type": "string"
                            }
                          }
                        }
                      }
                    },
                    "required": {
                      "type": "array",
                      "items": {
                        "type": "string"
                      }
                    }
                  }
                },
                "parameter": {
                  "type": "string"
                },
                "expected": {
                  "type": "string"
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_drug_adverse_events_by_chemblId`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_drug_adverse_events_by_chembl_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_drug_chemb_id_by_generic_name` {#dnakit.databases.open_targets_get_drug_chemb_id_by_generic_name}

Fetch the drug chemblId and description based on the drug generic name.

```text
dnakit.databases.open_targets_get_drug_chemb_id_by_generic_name(*, drugName=None, query=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `drugName` | string | None | The generic name of the drug for which the ID is required. |
| `query` | string | None | DNAKit optional alias for drugName |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "drugName": {
          "type": "string",
          "description": "The generic name of the drug for which the ID is required."
        },
        "query": {
          "type": "string",
          "description": "DNAKit optional alias for drugName"
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
        "data": {
          "type": "object",
          "properties": {
            "search": {
              "type": "object",
              "properties": {
                "hits": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "id": {
                        "type": "string"
                      },
                      "name": {
                        "type": "string"
                      },
                      "description": {
                        "type": "string"
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

??? info "DNAKit argument extension rules"

    ```json
    {
      "aliases": {
        "query": "drugName"
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_drug_chembId_by_generic_name`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_drug_chemb_id_by_generic_name", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_drug_indications_by_chembl_id` {#dnakit.databases.open_targets_get_drug_indications_by_chembl_id}

Fetch indications (treatable phenotypes/diseases) for a given drug chemblId.

```text
dnakit.databases.open_targets_get_drug_indications_by_chembl_id(*, chemblId, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `chemblId` | string | Required by input schema | The chemblId of the drug for which to retrieve treatable phenotypes information. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "chemblId": {
          "type": "string",
          "description": "The chemblId of the drug for which to retrieve treatable phenotypes information."
        }
      },
      "required": [
        "chemblId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "drug": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "name": {
                  "type": "string"
                },
                "indications": {
                  "type": "object",
                  "properties": {
                    "rows": {
                      "type": "array",
                      "items": {
                        "type": "object",
                        "properties": {
                          "disease": {
                            "type": "object",
                            "properties": {
                              "id": {
                                "type": "string"
                              },
                              "name": {
                                "type": "string"
                              }
                            }
                          },
                          "maxPhaseForIndication": {
                            "type": "number"
                          },
                          "references": {
                            "type": "array",
                            "items": {
                              "type": "object",
                              "properties": {
                                "source": {
                                  "type": "string"
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_drug_indications_by_chemblId`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_drug_indications_by_chembl_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_drug_mechanisms_of_action_by_chembl_id` {#dnakit.databases.open_targets_get_drug_mechanisms_of_action_by_chembl_id}

Retrieve the mechanisms of action associated with a specific drug using chemblId.

```text
dnakit.databases.open_targets_get_drug_mechanisms_of_action_by_chembl_id(*, chemblId, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `chemblId` | string | Required by input schema | The ChEMBL ID of the drug. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "chemblId": {
          "type": "string",
          "description": "The ChEMBL ID of the drug."
        }
      },
      "required": [
        "chemblId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "drug": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "name": {
                  "type": "string"
                },
                "mechanismsOfAction": {
                  "type": "object",
                  "properties": {
                    "rows": {
                      "type": "array",
                      "items": {
                        "type": "object",
                        "properties": {
                          "mechanismOfAction": {
                            "type": "string"
                          },
                          "actionType": {
                            "type": "string"
                          },
                          "targetName": {
                            "type": "string"
                          },
                          "targets": {
                            "type": "array",
                            "items": {
                              "type": "object",
                              "properties": {
                                "id": {
                                  "type": "string"
                                },
                                "approvedSymbol": {
                                  "type": "string"
                                }
                              }
                            }
                          },
                          "references": {
                            "type": "array",
                            "items": {
                              "type": "object",
                              "properties": {
                                "source": {
                                  "type": "string"
                                },
                                "urls": {
                                  "type": "array",
                                  "items": {
                                    "type": "string"
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_drug_mechanisms_of_action_by_chemblId`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_drug_mechanisms_of_action_by_chembl_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_drug_warnings_by_chembl_id` {#dnakit.databases.open_targets_get_drug_warnings_by_chembl_id}

Retrieve warnings for a specific drug using ChEMBL ID.

```text
dnakit.databases.open_targets_get_drug_warnings_by_chembl_id(*, chemblId, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `chemblId` | string | Required by input schema | The ChEMBL ID of the drug. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "chemblId": {
          "type": "string",
          "description": "The ChEMBL ID of the drug."
        }
      },
      "required": [
        "chemblId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "drug": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "name": {
                  "type": "string"
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_drug_warnings_by_chemblId`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_drug_warnings_by_chembl_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_evidence_by_datasource` {#dnakit.databases.open_targets_get_evidence_by_datasource}

Get target-disease evidence from Open Targets filtered by configurable data sources. Unlike OpenTargets_target_disease_evidence (intogen-only), this tool accepts any datasourceIds. Valid sources include: gwas_credible_sets, eva, eva_somatic, gene_burden, genomics_england, intogen, cancer_gene_census, cancer_biomarkers, clinical_precedence, crispr, europepmc, expression_atlas, gene2phenotype, impc, orphanet, uniprot_literature, uniprot_variants, clingen, phewas_catalog, reactome, progeny, slapenrich, sysbio. Note: two source IDs were renamed upstream and the retired names now match nothing rather than erroring — 'ot_genetics_portal' became 'gwas_credible_sets' (confirmed live: IL23R/Crohn disease returns 0 rows for the old name and 43 for the new one) and 'chembl' became 'clinical_precedence'. This tool remaps both retired names automatically and reports the substitution in metadata.datasource_rename_note. Returns evidence rows with scores, literature references, and source attribution.

```text
dnakit.databases.open_targets_get_evidence_by_datasource(*, efoId=None, ensemblId=None, gene_symbol=None, disease_name=None, datasourceIds=None, size=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | None | Disease EFO/MONDO ID (e.g., &#x27;MONDO_0005011&#x27; for Crohn disease). Alternative to disease_name. |
| `ensemblId` | string | None | Target Ensembl gene ID (e.g., &#x27;ENSG00000141510&#x27; for TP53). Alternative to gene_symbol. |
| `gene_symbol` | string | None | HGNC gene symbol (e.g., &#x27;TP53&#x27;, &#x27;BRCA1&#x27;). Auto-resolved to ensemblId. |
| `disease_name` | string | None | Disease or phenotype name (e.g., &#x27;Crohn disease&#x27;). Auto-resolved to efoId. |
| `datasourceIds` | array | None | List of datasource IDs to filter evidence. Examples: &#91;&#x27;gwas_credible_sets&#x27;, &#x27;europepmc&#x27;&#93;, &#91;&#x27;eva&#x27;, &#x27;clinical_precedence&#x27;&#93;, &#91;&#x27;intogen&#x27;, &#x27;cancer_gene_census&#x27;&#93;. For genetic (GWAS) association evidence use &#x27;gwas_credible_sets&#x27; — the retired &#x27;ot_genetics_portal&#x27; name is remapped automatically. Omit or pass empty array for all sources. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `size` | integer | None | Maximum evidence rows to return (default: 50) {&quot;default&quot;: 50} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "efoId": {
          "type": "string",
          "description": "Disease EFO/MONDO ID (e.g., 'MONDO_0005011' for Crohn disease). Alternative to disease_name."
        },
        "ensemblId": {
          "type": "string",
          "description": "Target Ensembl gene ID (e.g., 'ENSG00000141510' for TP53). Alternative to gene_symbol."
        },
        "gene_symbol": {
          "type": "string",
          "description": "HGNC gene symbol (e.g., 'TP53', 'BRCA1'). Auto-resolved to ensemblId."
        },
        "disease_name": {
          "type": "string",
          "description": "Disease or phenotype name (e.g., 'Crohn disease'). Auto-resolved to efoId."
        },
        "datasourceIds": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "List of datasource IDs to filter evidence. Examples: ['gwas_credible_sets', 'europepmc'], ['eva', 'clinical_precedence'], ['intogen', 'cancer_gene_census']. For genetic (GWAS) association evidence use 'gwas_credible_sets' — the retired 'ot_genetics_portal' name is remapped automatically. Omit or pass empty array for all sources."
        },
        "size": {
          "type": "integer",
          "description": "Maximum evidence rows to return (default: 50)",
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
        "data": {
          "type": "object",
          "properties": {
            "disease": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "name": {
                  "type": "string"
                },
                "evidences": {
                  "type": "object",
                  "properties": {
                    "count": {
                      "type": "integer"
                    },
                    "returned": {
                      "type": "integer",
                      "description": "Number of evidence rows actually present in `rows` for this page."
                    },
                    "truncated": {
                      "type": "boolean",
                      "description": "True when `returned` is less than `count`, i.e. evidence rows remain unfetched; `metadata.truncation_note` then names the `size` value to pass."
                    },
                    "rows": {
                      "type": "array",
                      "items": {
                        "type": "object",
                        "properties": {
                          "datasourceId": {
                            "type": "string"
                          },
                          "datatypeId": {
                            "type": "string"
                          },
                          "score": {
                            "type": "number"
                          },
                          "literature": {
                            "type": "array"
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        },
        "metadata": {
          "type": "object",
          "properties": {
            "truncation_note": {
              "type": "string",
              "description": "Present only when a collection came back truncated: states how many of how many rows were returned and names the exact argument that fetches the rest."
            },
            "datasource_rename_note": {
              "type": "string",
              "description": "Present only when a retired datasource ID was remapped to its current name."
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_evidence_by_datasource`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_evidence_by_datasource", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_gwas_study` {#dnakit.databases.open_targets_get_gwas_study}

Get detailed metadata for a GWAS study from Open Targets Platform, including trait information, sample sizes, publication details, disease mappings, cohorts, LD population structure, and availability of summary statistics. Use GCST accession IDs.

```text
dnakit.databases.open_targets_get_gwas_study(*, studyId, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `studyId` | string | Required by input schema | GWAS study accession ID (e.g., &#x27;GCST000392&#x27; for T1D Barrett et al., &#x27;GCST90029024&#x27; for T2D Loh et al.) |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "studyId": {
          "type": "string",
          "description": "GWAS study accession ID (e.g., 'GCST000392' for T1D Barrett et al., 'GCST90029024' for T2D Loh et al.)"
        }
      },
      "required": [
        "studyId"
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
            "study": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "studyType": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "traitFromSource": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "publicationFirstAuthor": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "publicationTitle": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "publicationDate": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "publicationJournal": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "pubmedId": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "nSamples": {
                  "type": [
                    "integer",
                    "null"
                  ]
                },
                "nCases": {
                  "type": [
                    "integer",
                    "null"
                  ]
                },
                "nControls": {
                  "type": [
                    "integer",
                    "null"
                  ]
                },
                "initialSampleSize": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "hasSumstats": {
                  "type": [
                    "boolean",
                    "null"
                  ]
                },
                "summarystatsLocation": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "cohorts": {
                  "type": [
                    "array",
                    "null"
                  ],
                  "items": {
                    "type": "string"
                  }
                },
                "diseases": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "id": {
                        "type": "string"
                      },
                      "name": {
                        "type": "string"
                      }
                    }
                  }
                },
                "ldPopulationStructure": {
                  "type": [
                    "array",
                    "null"
                  ],
                  "items": {
                    "type": "object",
                    "properties": {
                      "ldPopulation": {
                        "type": "string"
                      },
                      "relativeSampleSize": {
                        "type": "number"
                      }
                    }
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
              "type": "object"
            }
          },
          "required": [
            "error"
          ]
        }
      ]
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_gwas_study`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_gwas_study", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_publications_by_target_ensembl_id` {#dnakit.databases.open_targets_get_publications_by_target_ensembl_id}

Retrieve publications related to a target ensemblID, including PubMed IDs and publication dates. Only the first 25 publications are returned and this tool exposes no pagination parameter: `literatureOcurrences.count` is the total (often hundreds of thousands), `literatureOcurrences.returned` is how many are in `rows`, and `literatureOcurrences.truncated` is true when publications remain unfetched, in which case `metadata.truncation_note` says so. Narrow the query with startYear/endYear/additionalIds rather than treating `rows` as the literature for this target.

```text
dnakit.databases.open_targets_get_publications_by_target_ensembl_id(*, entityId, additionalIds=None, startYear=None, startMonth=None, endYear=None, endMonth=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `entityId` | string | Required by input schema | The ID of the entity (ensemblID). |
| `additionalIds` | array | None | List of additional IDs to include in the search. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `startYear` | integer | None | Year at the lower end of the filter. |
| `startMonth` | integer | None | Month at the lower end of the filter. |
| `endYear` | integer | None | Year at the higher end of the filter. |
| `endMonth` | integer | None | Month at the higher end of the filter. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "entityId": {
          "type": "string",
          "description": "The ID of the entity (ensemblID)."
        },
        "additionalIds": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "List of additional IDs to include in the search."
        },
        "startYear": {
          "type": "integer",
          "description": "Year at the lower end of the filter."
        },
        "startMonth": {
          "type": "integer",
          "description": "Month at the lower end of the filter."
        },
        "endYear": {
          "type": "integer",
          "description": "Year at the higher end of the filter."
        },
        "endMonth": {
          "type": "integer",
          "description": "Month at the higher end of the filter."
        }
      },
      "required": [
        "entityId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "target": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "approvedSymbol": {
                  "type": "string"
                },
                "literatureOcurrences": {
                  "type": "object",
                  "properties": {
                    "count": {
                      "type": "integer",
                      "description": "Total number of publications for this target — NOT the number of rows in this response."
                    },
                    "filteredCount": {
                      "type": "integer",
                      "description": "Total number of publications left after the startYear/startMonth/endYear/endMonth filter."
                    },
                    "earliestPubYear": {
                      "type": "integer"
                    },
                    "returned": {
                      "type": "integer",
                      "description": "Number of publications actually present in `rows` for this page."
                    },
                    "truncated": {
                      "type": "boolean",
                      "description": "True when `returned` is less than `count`, i.e. publications remain unfetched; `metadata.truncation_note` then says so. This tool exposes no pagination parameter, so narrow the query with startYear/endYear/additionalIds instead."
                    },
                    "rows": {
                      "type": "array",
                      "items": {
                        "type": "object",
                        "properties": {
                          "pmid": {
                            "type": "string"
                          },
                          "pmcid": {
                            "type": [
                              "string",
                              "null"
                            ]
                          },
                          "publicationDate": {
                            "type": [
                              "string",
                              "null"
                            ]
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        },
        "metadata": {
          "type": "object",
          "properties": {
            "truncation_note": {
              "type": "string",
              "description": "Present only when a collection came back truncated: states how many of how many rows were returned and why the withheld rows cannot be fetched through this tool."
            }
          }
        },
        "error": {
          "type": "string"
        },
        "error_details": {
          "type": "object",
          "properties": {
            "type": {
              "type": "string"
            },
            "message": {
              "type": "string"
            },
            "retriable": {
              "type": "boolean"
            },
            "next_steps": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "details": {
              "type": "object",
              "properties": {
                "validation_error": {
                  "type": "string"
                },
                "path": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  }
                },
                "schema": {
                  "type": "object",
                  "properties": {
                    "type": {
                      "type": "string"
                    },
                    "properties": {
                      "type": "object",
                      "properties": {
                        "entityId": {
                          "type": "object",
                          "properties": {
                            "type": {
                              "type": "string"
                            },
                            "description": {
                              "type": "string"
                            }
                          }
                        },
                        "additionalIds": {
                          "type": "object",
                          "properties": {
                            "type": {
                              "type": "string"
                            },
                            "items": {
                              "type": "object",
                              "properties": {
                                "type": {
                                  "type": "string"
                                }
                              }
                            },
                            "description": {
                              "type": "string"
                            }
                          }
                        },
                        "startYear": {
                          "type": "object",
                          "properties": {
                            "type": {
                              "type": "string"
                            },
                            "description": {
                              "type": "string"
                            }
                          }
                        },
                        "startMonth": {
                          "type": "object",
                          "properties": {
                            "type": {
                              "type": "string"
                            },
                            "description": {
                              "type": "string"
                            }
                          }
                        },
                        "endYear": {
                          "type": "object",
                          "properties": {
                            "type": {
                              "type": "string"
                            },
                            "description": {
                              "type": "string"
                            }
                          }
                        },
                        "endMonth": {
                          "type": "object",
                          "properties": {
                            "type": {
                              "type": "string"
                            },
                            "description": {
                              "type": "string"
                            }
                          }
                        }
                      }
                    },
                    "required": {
                      "type": "array",
                      "items": {
                        "type": "string"
                      }
                    }
                  }
                },
                "parameter": {
                  "type": "string"
                },
                "expected": {
                  "type": "string"
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_publications_by_target_ensemblID`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_publications_by_target_ensembl_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_similar_entities_by_disease_efo_id` {#dnakit.databases.open_targets_get_similar_entities_by_disease_efo_id}

Retrieve similar entities for a given disease efoId using a model trained with PubMed.

```text
dnakit.databases.open_targets_get_similar_entities_by_disease_efo_id(*, efoId, threshold, size, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | Required by input schema | The EFO ID of the disease. |
| `threshold` | number | Required by input schema | Threshold similarity between 0 and 1. Only results above threshold are returned. |
| `size` | integer | Required by input schema | Number of similar entities to fetch. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "efoId": {
          "type": "string",
          "description": "The EFO ID of the disease."
        },
        "threshold": {
          "type": "number",
          "description": "Threshold similarity between 0 and 1. Only results above threshold are returned."
        },
        "size": {
          "type": "integer",
          "description": "Number of similar entities to fetch."
        }
      },
      "required": [
        "efoId",
        "size",
        "threshold"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "disease": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "name": {
                  "type": "string"
                },
                "similarEntities": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "id": {
                        "type": "string"
                      },
                      "category": {
                        "type": "string"
                      },
                      "score": {
                        "type": "number"
                      },
                      "object": {
                        "type": "object",
                        "properties": {
                          "id": {
                            "type": "string"
                          },
                          "name": {
                            "type": "string"
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_similar_entities_by_disease_efoId`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_similar_entities_by_disease_efo_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_study_credible_sets` {#dnakit.databases.open_targets_get_study_credible_sets}

Get all credible sets (fine-mapped GWAS loci) for a specific GWAS study from Open Targets. Returns the lead variants, fine-mapping details, p-values, effect sizes, and locus-to-gene (L2G) predictions for each locus in the study. Use this to explore all significant loci discovered in a GWAS.

```text
dnakit.databases.open_targets_get_study_credible_sets(*, studyIds, size=None, index=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `studyIds` | array | Required by input schema | GWAS study accession IDs (e.g., &#91;&#x27;GCST000392&#x27;&#93;) {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `size` | integer | None | Number of credible sets to return (default 20) {&quot;default&quot;: 20} |
| `index` | integer | None | Page index for pagination (0-based, default 0) {&quot;default&quot;: 0} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "studyIds": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "GWAS study accession IDs (e.g., ['GCST000392'])"
        },
        "size": {
          "type": "integer",
          "description": "Number of credible sets to return (default 20)",
          "default": 20
        },
        "index": {
          "type": "integer",
          "description": "Page index for pagination (0-based, default 0)",
          "default": 0
        }
      },
      "required": [
        "studyIds"
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
            "credibleSets": {
              "type": "object",
              "properties": {
                "count": {
                  "type": "integer"
                },
                "rows": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "studyLocusId": {
                        "type": "string"
                      },
                      "studyId": {
                        "type": "string"
                      },
                      "region": {
                        "type": [
                          "string",
                          "null"
                        ]
                      },
                      "chromosome": {
                        "type": [
                          "string",
                          "null"
                        ]
                      },
                      "position": {
                        "type": [
                          "integer",
                          "null"
                        ]
                      },
                      "finemappingMethod": {
                        "type": [
                          "string",
                          "null"
                        ]
                      },
                      "pValueMantissa": {
                        "type": [
                          "number",
                          "null"
                        ]
                      },
                      "pValueExponent": {
                        "type": [
                          "integer",
                          "null"
                        ]
                      },
                      "beta": {
                        "type": [
                          "number",
                          "null"
                        ]
                      },
                      "confidence": {
                        "type": [
                          "string",
                          "null"
                        ]
                      },
                      "variant": {
                        "type": "object",
                        "properties": {
                          "id": {
                            "type": "string"
                          },
                          "rsIds": {
                            "type": "array",
                            "items": {
                              "type": "string"
                            }
                          }
                        }
                      },
                      "study": {
                        "type": "object",
                        "properties": {
                          "traitFromSource": {
                            "type": [
                              "string",
                              "null"
                            ]
                          },
                          "diseases": {
                            "type": "array",
                            "items": {
                              "type": "object",
                              "properties": {
                                "id": {
                                  "type": "string"
                                },
                                "name": {
                                  "type": "string"
                                }
                              }
                            }
                          }
                        }
                      },
                      "l2GPredictions": {
                        "type": "object",
                        "properties": {
                          "rows": {
                            "type": "array",
                            "items": {
                              "type": "object",
                              "properties": {
                                "target": {
                                  "type": "object",
                                  "properties": {
                                    "id": {
                                      "type": "string"
                                    },
                                    "approvedSymbol": {
                                      "type": "string"
                                    }
                                  }
                                },
                                "score": {
                                  "type": "number"
                                }
                              }
                            }
                          }
                        }
                      }
                    }
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
              "type": "object"
            }
          },
          "required": [
            "error"
          ]
        }
      ]
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_study_credible_sets`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_study_credible_sets", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_target_classes_by_ensembl_id` {#dnakit.databases.open_targets_get_target_classes_by_ensembl_id}

Retrieve the target classes associated with a specific target ensemblID.

```text
dnakit.databases.open_targets_get_target_classes_by_ensembl_id(*, ensemblId, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `ensemblId` | string | Required by input schema | The Ensembl ID of the target. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "ensemblId": {
          "type": "string",
          "description": "The Ensembl ID of the target."
        }
      },
      "required": [
        "ensemblId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "target": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "approvedSymbol": {
                  "type": "string"
                },
                "targetClass": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "id": {
                        "type": "integer"
                      },
                      "label": {
                        "type": "string"
                      },
                      "level": {
                        "type": "string"
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_target_classes_by_ensemblID`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_target_classes_by_ensembl_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_target_gene_ontology_by_ensembl_id` {#dnakit.databases.open_targets_get_target_gene_ontology_by_ensembl_id}

Retrieve Gene Ontology annotations for a specific target by Ensembl ID.

```text
dnakit.databases.open_targets_get_target_gene_ontology_by_ensembl_id(*, ensemblId, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `ensemblId` | string | Required by input schema | The Ensembl ID of the target for which to retrieve Gene Ontology annotations. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "ensemblId": {
          "type": "string",
          "description": "The Ensembl ID of the target for which to retrieve Gene Ontology annotations."
        }
      },
      "required": [
        "ensemblId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "target": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "approvedSymbol": {
                  "type": "string"
                },
                "geneOntology": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "aspect": {
                        "type": "string"
                      },
                      "evidence": {
                        "type": "string"
                      },
                      "geneProduct": {
                        "type": "string"
                      },
                      "source": {
                        "type": "string"
                      },
                      "term": {
                        "type": "object",
                        "properties": {
                          "id": {
                            "type": "string"
                          },
                          "name": {
                            "type": "string"
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_target_gene_ontology_by_ensemblID`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_target_gene_ontology_by_ensembl_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_target_genomic_location_by_ensembl_id` {#dnakit.databases.open_targets_get_target_genomic_location_by_ensembl_id}

Retrieve genomic location data for a specific target, including chromosome, start, end, and strand.

```text
dnakit.databases.open_targets_get_target_genomic_location_by_ensembl_id(*, ensemblId, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `ensemblId` | string | Required by input schema | The Ensembl ID of the target for which to retrieve genomic location information. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "ensemblId": {
          "type": "string",
          "description": "The Ensembl ID of the target for which to retrieve genomic location information."
        }
      },
      "required": [
        "ensemblId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "target": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "approvedSymbol": {
                  "type": "string"
                },
                "genomicLocation": {
                  "type": "object",
                  "properties": {
                    "chromosome": {
                      "type": "string"
                    },
                    "start": {
                      "type": "integer"
                    },
                    "end": {
                      "type": "integer"
                    },
                    "strand": {
                      "type": "integer"
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_target_genomic_location_by_ensemblID`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_target_genomic_location_by_ensembl_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_target_homologues_by_ensembl_id` {#dnakit.databases.open_targets_get_target_homologues_by_ensembl_id}

Fetch homologues for a specific target by Ensembl ID.

```text
dnakit.databases.open_targets_get_target_homologues_by_ensembl_id(*, ensemblId, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `ensemblId` | string | Required by input schema | The Ensembl ID of the target for which to retrieve homologues. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "ensemblId": {
          "type": "string",
          "description": "The Ensembl ID of the target for which to retrieve homologues."
        }
      },
      "required": [
        "ensemblId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "target": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "approvedSymbol": {
                  "type": "string"
                },
                "homologues": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "homologyType": {
                        "type": "string"
                      },
                      "queryPercentageIdentity": {
                        "type": "number"
                      },
                      "speciesId": {
                        "type": "string"
                      },
                      "speciesName": {
                        "type": "string"
                      },
                      "targetGeneId": {
                        "type": "string"
                      },
                      "targetGeneSymbol": {
                        "type": "string"
                      },
                      "targetPercentageIdentity": {
                        "type": "number"
                      },
                      "isHighConfidence": {
                        "type": "string"
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_target_homologues_by_ensemblID`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_target_homologues_by_ensembl_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_target_id_description_by_name` {#dnakit.databases.open_targets_get_target_id_description_by_name}

Get the ensemblId and description based on the target name.

```text
dnakit.databases.open_targets_get_target_id_description_by_name(*, targetName, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `targetName` | string | Required by input schema | The name of the target for which the ID is required. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "targetName": {
          "type": "string",
          "description": "The name of the target for which the ID is required."
        }
      },
      "required": [
        "targetName"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "search": {
              "type": "object",
              "properties": {
                "hits": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "id": {
                        "type": "string"
                      },
                      "name": {
                        "type": "string"
                      },
                      "description": {
                        "type": "string"
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_target_id_description_by_name`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_target_id_description_by_name", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_target_interactions_by_ensembl_id` {#dnakit.databases.open_targets_get_target_interactions_by_ensembl_id}

Retrieve interaction data for a specific target ensemblID, including interaction partners and evidence. Results are PAGINATED and Open Targets returns only 25 rows when no `page` is given: `interactions.count` is the total number of interactions (a well-studied target has thousands), `interactions.returned` is how many are in `rows`, and `interactions.truncated` is true when interactions remain unfetched, in which case `metadata.truncation_note` names the exact `page` argument to pass. Pass page: {"index": 0, "size": N} (Open Targets rejects a size above 3000) and step `index` to page through.

```text
dnakit.databases.open_targets_get_target_interactions_by_ensembl_id(*, ensemblId, page=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `ensemblId` | string | Required by input schema | The Ensembl ID of the target. |
| `page` | object | None | Pagination parameters. {&quot;properties&quot;: {&quot;index&quot;: {&quot;type&quot;: &quot;integer&quot;, &quot;description&quot;: &quot;The index of the page to retrieve.&quot;}, &quot;size&quot;: {&quot;type&quot;: &quot;integer&quot;, &quot;description&quot;: &quot;The number of items per page.&quot;}}} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "ensemblId": {
          "type": "string",
          "description": "The Ensembl ID of the target."
        },
        "page": {
          "type": "object",
          "properties": {
            "index": {
              "type": "integer",
              "description": "The index of the page to retrieve."
            },
            "size": {
              "type": "integer",
              "description": "The number of items per page."
            }
          },
          "description": "Pagination parameters."
        }
      },
      "required": [
        "ensemblId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "target": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "approvedSymbol": {
                  "type": "string"
                },
                "interactions": {
                  "type": "object",
                  "properties": {
                    "count": {
                      "type": "integer",
                      "description": "Total number of interactions for this target across all pages — NOT the number of rows in this response."
                    },
                    "returned": {
                      "type": "integer",
                      "description": "Number of interactions actually present in `rows` for this page."
                    },
                    "truncated": {
                      "type": "boolean",
                      "description": "True when `returned` is less than `count`, i.e. interactions remain unfetched; `metadata.truncation_note` then names the exact `page` argument to pass."
                    },
                    "rows": {
                      "type": "array",
                      "items": {
                        "type": "object",
                        "properties": {
                          "intA": {
                            "type": "string"
                          },
                          "targetA": {
                            "type": [
                              "object",
                              "null"
                            ],
                            "properties": {
                              "id": {
                                "type": "string"
                              },
                              "approvedSymbol": {
                                "type": "string"
                              }
                            }
                          },
                          "intB": {
                            "type": "string"
                          },
                          "targetB": {
                            "type": [
                              "object",
                              "null"
                            ],
                            "properties": {
                              "id": {
                                "type": "string"
                              },
                              "approvedSymbol": {
                                "type": "string"
                              }
                            }
                          },
                          "score": {
                            "type": [
                              "number",
                              "null"
                            ]
                          },
                          "sourceDatabase": {
                            "type": "string"
                          },
                          "evidences": {
                            "type": "array",
                            "items": {
                              "type": "object",
                              "properties": {
                                "evidenceScore": {
                                  "type": [
                                    "number",
                                    "null"
                                  ]
                                },
                                "pubmedId": {
                                  "type": [
                                    "string",
                                    "null"
                                  ]
                                },
                                "interactionDetectionMethodShortName": {
                                  "type": [
                                    "string",
                                    "null"
                                  ]
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        },
        "metadata": {
          "type": "object",
          "properties": {
            "truncation_note": {
              "type": "string",
              "description": "Present only when a collection came back truncated: states how many of how many rows were returned and names the exact `page` argument that fetches the rest."
            }
          }
        },
        "error": {
          "type": "string"
        },
        "error_details": {
          "type": "object",
          "properties": {
            "type": {
              "type": "string"
            },
            "message": {
              "type": "string"
            },
            "retriable": {
              "type": "boolean"
            },
            "next_steps": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "details": {
              "type": "object",
              "properties": {
                "validation_error": {
                  "type": "string"
                },
                "path": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  }
                },
                "schema": {
                  "type": "object",
                  "properties": {
                    "type": {
                      "type": "string"
                    },
                    "properties": {
                      "type": "object",
                      "properties": {
                        "ensemblId": {
                          "type": "object",
                          "properties": {
                            "type": {
                              "type": "string"
                            },
                            "description": {
                              "type": "string"
                            }
                          }
                        },
                        "page": {
                          "type": "object",
                          "properties": {
                            "type": {
                              "type": "string"
                            },
                            "properties": {
                              "type": "object",
                              "properties": {
                                "index": {
                                  "type": "object",
                                  "properties": {
                                    "type": {
                                      "type": "string"
                                    },
                                    "description": {
                                      "type": "string"
                                    }
                                  }
                                },
                                "size": {
                                  "type": "object",
                                  "properties": {
                                    "type": {
                                      "type": "string"
                                    },
                                    "description": {
                                      "type": "string"
                                    }
                                  }
                                }
                              }
                            },
                            "description": {
                              "type": "string"
                            }
                          }
                        }
                      }
                    },
                    "required": {
                      "type": "array",
                      "items": {
                        "type": "string"
                      }
                    }
                  }
                },
                "parameter": {
                  "type": "string"
                },
                "expected": {
                  "type": "string"
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_target_interactions_by_ensemblID`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_target_interactions_by_ensembl_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_target_safety_profile_by_ensembl_id` {#dnakit.databases.open_targets_get_target_safety_profile_by_ensembl_id}

Retrieve known target safety liabilities for a specific target Ensembl ID.

```text
dnakit.databases.open_targets_get_target_safety_profile_by_ensembl_id(*, ensemblId, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `ensemblId` | string | Required by input schema | The Ensembl ID of the target for which to retrieve safety liabilities. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "ensemblId": {
          "type": "string",
          "description": "The Ensembl ID of the target for which to retrieve safety liabilities."
        }
      },
      "required": [
        "ensemblId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "target": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "approvedSymbol": {
                  "type": "string"
                },
                "safetyLiabilities": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "event": {
                        "type": "string"
                      },
                      "eventId": {
                        "type": "string"
                      },
                      "studies": {
                        "type": "array",
                        "items": {
                          "type": "object",
                          "properties": {
                            "name": {
                              "type": "string"
                            },
                            "type": {
                              "type": "string"
                            },
                            "description": {
                              "type": "string"
                            }
                          }
                        }
                      },
                      "datasource": {
                        "type": "string"
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_target_safety_profile_by_ensemblID`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_target_safety_profile_by_ensembl_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_target_tractability_by_ensembl_id` {#dnakit.databases.open_targets_get_target_tractability_by_ensembl_id}

Retrieve tractability assessments, including modality and values, for a specific target ensembl ID.

```text
dnakit.databases.open_targets_get_target_tractability_by_ensembl_id(*, ensemblId, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `ensemblId` | string | Required by input schema | The Ensembl ID of the target. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "ensemblId": {
          "type": "string",
          "description": "The Ensembl ID of the target."
        }
      },
      "required": [
        "ensemblId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "target": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "approvedSymbol": {
                  "type": "string"
                },
                "tractability": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "label": {
                        "type": "string"
                      },
                      "modality": {
                        "type": "string"
                      },
                      "value": {
                        "type": "boolean"
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_target_tractability_by_ensemblID`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_target_tractability_by_ensembl_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_variant_credible_sets` {#dnakit.databases.open_targets_get_variant_credible_sets}

Get GWAS credible sets (fine-mapped loci) associated with a genetic variant from Open Targets. Returns studies where this variant appears in a fine-mapped credible set, with locus-to-gene (L2G) predictions, p-values, effect sizes, and the fine-mapping method used. Essential for understanding which GWAS signals a variant is part of and which genes it likely affects.

```text
dnakit.databases.open_targets_get_variant_credible_sets(*, variantId=None, size=None, variant_id=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `variantId` | string | None | Open Targets variant ID in format chr_position_ref_alt (e.g., &#x27;10_112998590_C_T&#x27; for rs7903146) |
| `size` | integer | None | Number of credible sets to return (default 10) {&quot;default&quot;: 10} |
| `variant_id` | string | None | DNAKit optional alias for variantId |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "variantId": {
          "type": "string",
          "description": "Open Targets variant ID in format chr_position_ref_alt (e.g., '10_112998590_C_T' for rs7903146)"
        },
        "size": {
          "type": "integer",
          "description": "Number of credible sets to return (default 10)",
          "default": 10
        },
        "variant_id": {
          "type": "string",
          "description": "DNAKit optional alias for variantId"
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
            "variant": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "rsIds": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  }
                },
                "credibleSets": {
                  "type": "object",
                  "properties": {
                    "count": {
                      "type": "integer"
                    },
                    "rows": {
                      "type": "array",
                      "items": {
                        "type": "object",
                        "properties": {
                          "studyLocusId": {
                            "type": "string"
                          },
                          "studyId": {
                            "type": "string"
                          },
                          "region": {
                            "type": [
                              "string",
                              "null"
                            ]
                          },
                          "finemappingMethod": {
                            "type": [
                              "string",
                              "null"
                            ]
                          },
                          "pValueMantissa": {
                            "type": [
                              "number",
                              "null"
                            ]
                          },
                          "pValueExponent": {
                            "type": [
                              "integer",
                              "null"
                            ]
                          },
                          "beta": {
                            "type": [
                              "number",
                              "null"
                            ]
                          },
                          "confidence": {
                            "type": [
                              "string",
                              "null"
                            ]
                          },
                          "study": {
                            "type": "object",
                            "properties": {
                              "id": {
                                "type": "string"
                              },
                              "traitFromSource": {
                                "type": [
                                  "string",
                                  "null"
                                ]
                              },
                              "publicationFirstAuthor": {
                                "type": [
                                  "string",
                                  "null"
                                ]
                              },
                              "pubmedId": {
                                "type": [
                                  "string",
                                  "null"
                                ]
                              },
                              "diseases": {
                                "type": "array",
                                "items": {
                                  "type": "object",
                                  "properties": {
                                    "id": {
                                      "type": "string"
                                    },
                                    "name": {
                                      "type": "string"
                                    }
                                  }
                                }
                              }
                            }
                          },
                          "l2GPredictions": {
                            "type": "object",
                            "properties": {
                              "rows": {
                                "type": "array",
                                "items": {
                                  "type": "object",
                                  "properties": {
                                    "target": {
                                      "type": "object",
                                      "properties": {
                                        "id": {
                                          "type": "string"
                                        },
                                        "approvedSymbol": {
                                          "type": "string"
                                        }
                                      }
                                    },
                                    "score": {
                                      "type": "number"
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
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
              "type": "object"
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
        "variant_id": "variantId"
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_variant_credible_sets`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_variant_credible_sets", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_get_variant_info` {#dnakit.databases.open_targets_get_variant_info}

Get detailed information about a genetic variant from Open Targets Platform, including genomic coordinates, rsIDs, allele frequencies across populations (gnomAD), functional consequence, transcript consequences, and HGVS notation. Useful for GWAS variant annotation and population genetics.

```text
dnakit.databases.open_targets_get_variant_info(*, variantId, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `variantId` | string | Required by input schema | Open Targets variant ID in format chr_position_ref_alt (e.g., &#x27;10_112998590_C_T&#x27; for rs7903146, &#x27;19_44908684_T_C&#x27; for rs429358). Use the search tool to find variant IDs from rsIDs. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "variantId": {
          "type": "string",
          "description": "Open Targets variant ID in format chr_position_ref_alt (e.g., '10_112998590_C_T' for rs7903146, '19_44908684_T_C' for rs429358). Use the search tool to find variant IDs from rsIDs."
        }
      },
      "required": [
        "variantId"
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
            "variant": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string",
                  "description": "Variant identifier in chr_pos_ref_alt format"
                },
                "rsIds": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  },
                  "description": "dbSNP rs identifiers"
                },
                "chromosome": {
                  "type": "string"
                },
                "position": {
                  "type": "integer"
                },
                "referenceAllele": {
                  "type": "string"
                },
                "alternateAllele": {
                  "type": "string"
                },
                "hgvsId": {
                  "type": [
                    "string",
                    "null"
                  ],
                  "description": "HGVS notation"
                },
                "variantDescription": {
                  "type": "string",
                  "description": "Human-readable variant description"
                },
                "mostSevereConsequence": {
                  "type": "object",
                  "properties": {
                    "id": {
                      "type": "string"
                    },
                    "label": {
                      "type": "string"
                    }
                  }
                },
                "alleleFrequencies": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "populationName": {
                        "type": "string"
                      },
                      "alleleFrequency": {
                        "type": "number"
                      }
                    }
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
              "type": "object"
            }
          },
          "required": [
            "error"
          ]
        }
      ]
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_get_variant_info`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_get_variant_info", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_map_any_disease_id_to_all_other_ids` {#dnakit.databases.open_targets_map_any_disease_id_to_all_other_ids}

Given any known disease or phenotype ID (EFO, OMIM, MONDO, UMLS, ICD10, MedDRA, etc.), return all known cross-referenced IDs including the EFO ID.

```text
dnakit.databases.open_targets_map_any_disease_id_to_all_other_ids(*, inputId, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `inputId` | string | Required by input schema | Any known disease ID (e.g. OMIM:604302, UMLS:C0003873, ICD10:M05, MONDO_0008383, etc.) |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "inputId": {
          "type": "string",
          "description": "Any known disease ID (e.g. OMIM:604302, UMLS:C0003873, ICD10:M05, MONDO_0008383, etc.)"
        }
      },
      "required": [
        "inputId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "mapIds": {
              "type": "object",
              "properties": {
                "mappings": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "term": {
                        "type": "string"
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_map_any_disease_id_to_all_other_ids`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_map_any_disease_id_to_all_other_ids", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_multi_entity_search_by_query_string` {#dnakit.databases.open_targets_multi_entity_search_by_query_string}

Perform a multi-entity search based on a query string, filtering by entity names and pagination settings.

```text
dnakit.databases.open_targets_multi_entity_search_by_query_string(*, queryString, entityNames=None, page=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `queryString` | string | Required by input schema | The search string for querying information. |
| `entityNames` | array | None | List of entity names to search for (e.g., target, disease, drug). {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `page` | object | None | Pagination settings with index and size. {&quot;properties&quot;: {&quot;index&quot;: {&quot;type&quot;: &quot;integer&quot;, &quot;description&quot;: &quot;Pagination index.&quot;}, &quot;size&quot;: {&quot;type&quot;: &quot;integer&quot;, &quot;description&quot;: &quot;Pagination size.&quot;}}} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "queryString": {
          "type": "string",
          "description": "The search string for querying information."
        },
        "entityNames": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "List of entity names to search for (e.g., target, disease, drug)."
        },
        "page": {
          "type": "object",
          "properties": {
            "index": {
              "type": "integer",
              "description": "Pagination index."
            },
            "size": {
              "type": "integer",
              "description": "Pagination size."
            }
          },
          "description": "Pagination settings with index and size."
        }
      },
      "required": [
        "queryString"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "error": {
          "type": "string"
        },
        "error_details": {
          "type": "object",
          "properties": {
            "type": {
              "type": "string"
            },
            "message": {
              "type": "string"
            },
            "retriable": {
              "type": "boolean"
            },
            "next_steps": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "details": {
              "type": "object",
              "properties": {
                "validation_error": {
                  "type": "string"
                },
                "path": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  }
                },
                "schema": {
                  "type": "object",
                  "properties": {
                    "type": {
                      "type": "string"
                    },
                    "properties": {
                      "type": "object",
                      "properties": {
                        "queryString": {
                          "type": "object",
                          "properties": {
                            "type": {
                              "type": "string"
                            },
                            "description": {
                              "type": "string"
                            }
                          }
                        },
                        "entityNames": {
                          "type": "object",
                          "properties": {
                            "type": {
                              "type": "string"
                            },
                            "items": {
                              "type": "object",
                              "properties": {
                                "type": {
                                  "type": "string"
                                }
                              }
                            },
                            "description": {
                              "type": "string"
                            }
                          }
                        },
                        "page": {
                          "type": "object",
                          "properties": {
                            "type": {
                              "type": "string"
                            },
                            "properties": {
                              "type": "object",
                              "properties": {
                                "index": {
                                  "type": "object",
                                  "properties": {
                                    "type": {
                                      "type": "string"
                                    },
                                    "description": {
                                      "type": "string"
                                    }
                                  }
                                },
                                "size": {
                                  "type": "object",
                                  "properties": {
                                    "type": {
                                      "type": "string"
                                    },
                                    "description": {
                                      "type": "string"
                                    }
                                  }
                                }
                              }
                            },
                            "description": {
                              "type": "string"
                            }
                          }
                        }
                      }
                    },
                    "required": {
                      "type": "array",
                      "items": {
                        "type": "string"
                      }
                    }
                  }
                },
                "parameter": {
                  "type": "string"
                },
                "expected": {
                  "type": "string"
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_multi_entity_search_by_query_string`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_multi_entity_search_by_query_string", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_search_gwas_studies_by_disease` {#dnakit.databases.open_targets_search_gwas_studies_by_disease}

Search for GWAS studies associated with a disease or trait in Open Targets. Takes disease ontology IDs (EFO, MONDO) and returns matching GWAS studies with their metadata including sample sizes, traits, publications, and associated diseases. Use enableIndirect=true to include child disease terms.

```text
dnakit.databases.open_targets_search_gwas_studies_by_disease(*, diseaseIds=None, enableIndirect=None, disease_name=None, disease=None, trait=None, size=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `diseaseIds` | array | None | Disease ontology IDs (e.g., &#91;&#x27;MONDO_0005148&#x27;&#93; for type 2 diabetes, &#91;&#x27;MONDO_0004975&#x27;&#93; for Alzheimer disease) {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `enableIndirect` | boolean | None | Include studies for child disease terms (default true) {&quot;default&quot;: true} |
| `disease_name` | string | None | Disease or trait name for auto-resolution (e.g., &#x27;type 2 diabetes&#x27;). Resolved to diseaseIds automatically. |
| `disease` | string | None | Alias for disease_name. |
| `trait` | string | None | Alias for disease_name. |
| `size` | integer | None | Limit for the returned page; does not automatically retrieve later pages. {&quot;minimum&quot;: 1, &quot;maximum&quot;: 1000} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "diseaseIds": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "Disease ontology IDs (e.g., ['MONDO_0005148'] for type 2 diabetes, ['MONDO_0004975'] for Alzheimer disease)"
        },
        "enableIndirect": {
          "type": "boolean",
          "description": "Include studies for child disease terms (default true)",
          "default": true
        },
        "disease_name": {
          "type": "string",
          "description": "Disease or trait name for auto-resolution (e.g., 'type 2 diabetes'). Resolved to diseaseIds automatically."
        },
        "disease": {
          "type": "string",
          "description": "Alias for disease_name."
        },
        "trait": {
          "type": "string",
          "description": "Alias for disease_name."
        },
        "size": {
          "type": "integer",
          "minimum": 1,
          "maximum": 1000
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
            "studies": {
              "type": "object",
              "properties": {
                "count": {
                  "type": "integer"
                },
                "rows": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "id": {
                        "type": "string"
                      },
                      "studyType": {
                        "type": [
                          "string",
                          "null"
                        ]
                      },
                      "traitFromSource": {
                        "type": [
                          "string",
                          "null"
                        ]
                      },
                      "publicationFirstAuthor": {
                        "type": [
                          "string",
                          "null"
                        ]
                      },
                      "publicationTitle": {
                        "type": [
                          "string",
                          "null"
                        ]
                      },
                      "publicationDate": {
                        "type": [
                          "string",
                          "null"
                        ]
                      },
                      "pubmedId": {
                        "type": [
                          "string",
                          "null"
                        ]
                      },
                      "nSamples": {
                        "type": [
                          "integer",
                          "null"
                        ]
                      },
                      "nCases": {
                        "type": [
                          "integer",
                          "null"
                        ]
                      },
                      "nControls": {
                        "type": [
                          "integer",
                          "null"
                        ]
                      },
                      "hasSumstats": {
                        "type": [
                          "boolean",
                          "null"
                        ]
                      },
                      "diseases": {
                        "type": "array",
                        "items": {
                          "type": "object",
                          "properties": {
                            "id": {
                              "type": "string"
                            },
                            "name": {
                              "type": "string"
                            }
                          }
                        }
                      }
                    }
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
              "type": "object"
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
      "parameters": {
        "size": {
          "type": "integer",
          "minimum": 1,
          "maximum": 1000
        }
      },
      "handler": "page_size"
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_search_gwas_studies_by_disease`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_search_gwas_studies_by_disease", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.open_targets_target_disease_evidence` {#dnakit.databases.open_targets_target_disease_evidence}

Explore IntOGen somatic driver evidence for a target-disease association. IMPORTANT: IntOGen only covers somatic tumor driver mutations — this tool returns empty data ({}) for non-cancer diseases (e.g., type 2 diabetes, hypertension). Empty results are expected and normal when querying non-oncology phenotypes. For non-cancer diseases or to query other data sources, use OpenTargets_get_evidence_by_datasource instead. Accepts gene_symbol (e.g. 'TP53') and disease_name (e.g. 'Crohn disease') as alternatives to ensemblId/efoId — IDs are auto-resolved via OpenTargets search. Results are PAGINATED and Open Targets returns only 25 evidence rows: `evidences.count` is the total, `evidences.returned` is how many are in `rows`, and `evidences.truncated` is true when evidence remains unfetched, in which case `metadata.truncation_note` says so. This tool exposes no pagination parameter, so use OpenTargets_get_evidence_by_datasource (which takes `size`) when the full evidence set is needed.

```text
dnakit.databases.open_targets_target_disease_evidence(*, efoId=None, ensemblId=None, gene_symbol=None, disease_name=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | None | EFO/MONDO disease ID (e.g., MONDO_0005011). Alternative to disease_name. |
| `ensemblId` | string | None | Ensembl gene ID (e.g., ENSG00000141510). Alternative to gene_symbol. |
| `gene_symbol` | string | None | HGNC gene symbol (e.g., &#x27;TP53&#x27;, &#x27;BRCA1&#x27;). Auto-resolved to ensemblId. |
| `disease_name` | string | None | Disease or phenotype name (e.g., &#x27;Crohn disease&#x27;, &#x27;breast carcinoma&#x27;). Auto-resolved to efoId. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "efoId": {
          "type": "string",
          "description": "EFO/MONDO disease ID (e.g., MONDO_0005011). Alternative to disease_name."
        },
        "ensemblId": {
          "type": "string",
          "description": "Ensembl gene ID (e.g., ENSG00000141510). Alternative to gene_symbol."
        },
        "gene_symbol": {
          "type": "string",
          "description": "HGNC gene symbol (e.g., 'TP53', 'BRCA1'). Auto-resolved to ensemblId."
        },
        "disease_name": {
          "type": "string",
          "description": "Disease or phenotype name (e.g., 'Crohn disease', 'breast carcinoma'). Auto-resolved to efoId."
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
        "data": {
          "type": "object",
          "properties": {
            "disease": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "string"
                },
                "name": {
                  "type": "string"
                },
                "evidences": {
                  "type": "object",
                  "properties": {
                    "count": {
                      "type": "integer"
                    },
                    "returned": {
                      "type": "integer",
                      "description": "Number of evidence rows actually present in `rows` for this page."
                    },
                    "truncated": {
                      "type": "boolean",
                      "description": "True when `returned` is less than `count`, i.e. evidence rows remain unfetched; `metadata.truncation_note` then says so. This tool exposes no pagination parameter, so use OpenTargets_get_evidence_by_datasource (which takes `size`) to reach the rest."
                    }
                  }
                }
              }
            }
          }
        },
        "metadata": {
          "type": "object",
          "properties": {
            "truncation_note": {
              "type": "string",
              "description": "Present only when a collection came back truncated: states how many of how many rows were returned and why the withheld rows cannot be fetched through this tool."
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OpenTargets_target_disease_evidence`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("open_targets_target_disease_evidence", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.pub_chem_get_cid_by_compound_name` {#dnakit.databases.pub_chem_get_cid_by_compound_name}

Search PubChem by compound name to get CID(s). This is the PRIMARY tool for looking up a drug or chemical by name (e.g., 'cytarabine', 'aspirin', 'metformin', 'venetoclax') to obtain PubChem Compound IDs. Use the returned CID with other PubChem tools for properties, patents, etc. Only works with chemical/drug names, not disease names.

```text
dnakit.databases.pub_chem_get_cid_by_compound_name(*, name=None, compound_name=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `name` | string | None | Chemical compound name (e.g., &quot;Aspirin&quot;, &quot;Acetaminophen&quot;) or IUPAC name. Do not use disease names or medical conditions. |
| `compound_name` | string | None | Alias for name. The compound name to look up. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "name": {
          "type": "string",
          "description": "Chemical compound name (e.g., \"Aspirin\", \"Acetaminophen\") or IUPAC name. Do not use disease names or medical conditions."
        },
        "compound_name": {
          "type": "string",
          "description": "Alias for name. The compound name to look up."
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
        "IdentifierList": {
          "type": "object",
          "properties": {
            "CID": {
              "type": "array",
              "items": {
                "type": "integer"
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `PubChem_get_CID_by_compound_name`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("pub_chem_get_cid_by_compound_name", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.pub_chem_get_compound_properties_by_cid` {#dnakit.databases.pub_chem_get_compound_properties_by_cid}

Get a set of specified molecular properties through CID (Compound ID), such as molecular weight, IUPAC name, and SMILES. PubChem's PUG-REST renamed its SMILES properties: ask for 'ConnectivitySMILES' (formerly 'CanonicalSMILES') and 'SMILES' (formerly 'IsomericSMILES'). Requesting a legacy name returns the value under the CURRENT key, so a caller that looks up the key it asked for will not find it.

```text
dnakit.databases.pub_chem_get_compound_properties_by_cid(*, cid, properties=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `cid` | integer | Required by input schema | PubChem compound ID to query, e.g., 2244 (Aspirin). |
| `properties` | array &#124; string &#124; null | None | Optional list of PubChem property names to retrieve (e.g. &#91;&#x27;MolecularFormula&#x27;,&#x27;InChIKey&#x27;,&#x27;XLogP&#x27;,&#x27;IUPACName&#x27;,&#x27;MolecularWeight&#x27;&#93;). If omitted, returns a default set (MolecularWeight, IUPACName, ConnectivitySMILES). For SMILES use the current names &#x27;ConnectivitySMILES&#x27; (was &#x27;CanonicalSMILES&#x27;) or &#x27;SMILES&#x27; (was &#x27;IsomericSMILES&#x27;) — a legacy name still returns data, but under the current key, not the one requested. See PubChem PUG-REST property tables for valid names. Accepts a list or a comma-separated string. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "cid": {
          "type": "integer",
          "description": "PubChem compound ID to query, e.g., 2244 (Aspirin)."
        },
        "properties": {
          "type": [
            "array",
            "string",
            "null"
          ],
          "items": {
            "type": "string"
          },
          "description": "Optional list of PubChem property names to retrieve (e.g. ['MolecularFormula','InChIKey','XLogP','IUPACName','MolecularWeight']). If omitted, returns a default set (MolecularWeight, IUPACName, ConnectivitySMILES). For SMILES use the current names 'ConnectivitySMILES' (was 'CanonicalSMILES') or 'SMILES' (was 'IsomericSMILES') — a legacy name still returns data, but under the current key, not the one requested. See PubChem PUG-REST property tables for valid names. Accepts a list or a comma-separated string."
        }
      },
      "required": [
        "cid"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "PropertyTable": {
          "type": "object",
          "properties": {
            "Properties": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "CID": {
                    "type": "integer"
                  },
                  "MolecularWeight": {
                    "type": "string"
                  },
                  "ConnectivitySMILES": {
                    "type": "string"
                  },
                  "IUPACName": {
                    "type": "string"
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `PubChem_get_compound_properties_by_CID`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("pub_chem_get_compound_properties_by_cid", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.reactome_disease_target_score` {#dnakit.databases.reactome_disease_target_score}

Extract disease-target association scores from Reactome pathway data. This includes pathway-based evidence.

```text
dnakit.databases.reactome_disease_target_score(*, efoId, pageSize=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | Required by input schema | The EFO (Experimental Factor Ontology) ID of the disease, e.g., &#x27;MONDO_0011996&#x27; for chronic myelogenous leukemia |
| `pageSize` | integer | None | Number of results per page (default: 100, max: 100) {&quot;default&quot;: 100} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "efoId": {
          "type": "string",
          "description": "The EFO (Experimental Factor Ontology) ID of the disease, e.g., 'MONDO_0011996' for chronic myelogenous leukemia"
        },
        "pageSize": {
          "type": "integer",
          "description": "Number of results per page (default: 100, max: 100)",
          "default": 100
        }
      },
      "required": [
        "efoId"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "description": "Reactome pathway disease-target association scores",
      "properties": {
        "disease_info": {
          "type": "object",
          "properties": {
            "disease_id": {
              "type": "string"
            },
            "disease_name": {
              "type": "string"
            }
          }
        },
        "datasource": {
          "type": "string",
          "description": "Always 'reactome'"
        },
        "total_targets_with_scores": {
          "type": "integer"
        },
        "target_scores": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "target_symbol": {
                "type": "string"
              },
              "target_id": {
                "type": "string"
              },
              "datasource": {
                "type": "string"
              },
              "score": {
                "type": "number",
                "description": "Reactome pathway score (0-1)"
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `reactome_disease_target_score`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("reactome_disease_target_score", verbose=False)
print(info["input_schema"])
```

