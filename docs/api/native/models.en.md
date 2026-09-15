# Prediction with trained models

This page documents **14 newly added public functions**. [Complete index](index.md) · [Conventions and examples](usage.md)

Each function lists its current Python signature, parameters, returns and registered execution requirements. Provider descriptions retain source wording. In generated database wrappers, None normally omits a parameter; actual required fields and conditional constraints come from the input schema. Follow each operation's coordinate and assembly conventions.

These functions call existing models or prediction services. Weights, SDKs, credentials and service terms remain provider requirements; DNAKit did not retrain these models.

## Function index

- [`dnakit.databases.alpha_genome_predict_interval`](#dnakit.databases.alpha_genome_predict_interval)
- [`dnakit.databases.alpha_genome_score_variant`](#dnakit.databases.alpha_genome_score_variant)
- [`dnakit.databases.dyna_mut2_predict_stability`](#dnakit.databases.dyna_mut2_predict_stability)
- [`dnakit.databases.esm2_score_missense_variant`](#dnakit.databases.esm2_score_missense_variant)
- [`dnakit.databases.esm_describe_sae_feature`](#dnakit.databases.esm_describe_sae_feature)
- [`dnakit.databases.esm_explain_variant_mechanism`](#dnakit.databases.esm_explain_variant_mechanism)
- [`dnakit.databases.esm_get_region_sae_features`](#dnakit.databases.esm_get_region_sae_features)
- [`dnakit.databases.esm_get_sae_features`](#dnakit.databases.esm_get_sae_features)
- [`dnakit.databases.esm_score_sequence`](#dnakit.databases.esm_score_sequence)
- [`dnakit.databases.esm_score_variant_sae_batch`](#dnakit.databases.esm_score_variant_sae_batch)
- [`dnakit.databases.esm_score_variant_sae_disruption`](#dnakit.databases.esm_score_variant_sae_disruption)
- [`dnakit.databases.evo2_score_variant`](#dnakit.databases.evo2_score_variant)
- [`dnakit.databases.nvidia_nim_alphafold2`](#dnakit.databases.nvidia_nim_alphafold2)
- [`dnakit.databases.nvidia_nim_esmfold`](#dnakit.databases.nvidia_nim_esmfold)

## `dnakit.databases.alpha_genome_predict_interval` {#dnakit.databases.alpha_genome_predict_interval}

Predict multimodal genomic tracks for a genomic interval with DeepMind AlphaGenome (single DNA-sequence model; up to 1 Mb at single-base resolution). Returns a compact per-modality summary (track counts and shapes) for the requested outputs. Hosted API via the alphagenome SDK; free for non-commercial use (set ALPHA_GENOME_API_KEY).

```text
dnakit.databases.alpha_genome_predict_interval(*, chromosome, start, end, output_types=None, ontology_terms=None, organism=None, sequence_length=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `chromosome` | string | Required by input schema | Chromosome, e.g. &#x27;chr19&#x27;. |
| `start` | integer | Required by input schema | Interval start (0-based). |
| `end` | integer | Required by input schema | Interval end. |
| `output_types` | array &#124; null | None | Modalities to predict, e.g. &#91;&#x27;RNA_SEQ&#x27;,&#x27;ATAC&#x27;&#93; (default &#91;&#x27;RNA_SEQ&#x27;&#93;). {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `ontology_terms` | array &#124; null | None | Optional tissue/cell ontology terms (e.g. &#91;&#x27;UBERON:0001114&#x27;&#93; = liver) to restrict tracks. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `organism` | string &#124; null | None | &#x27;human&#x27; (default) or &#x27;mouse&#x27;. |
| `sequence_length` | string &#124; null | None | Context window: 16KB, 100KB, 500KB, or 1MB (default). |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "required": [
        "chromosome",
        "end",
        "start"
      ],
      "properties": {
        "chromosome": {
          "type": "string",
          "description": "Chromosome, e.g. 'chr19'."
        },
        "start": {
          "type": "integer",
          "description": "Interval start (0-based)."
        },
        "end": {
          "type": "integer",
          "description": "Interval end."
        },
        "output_types": {
          "type": [
            "array",
            "null"
          ],
          "items": {
            "type": "string"
          },
          "description": "Modalities to predict, e.g. ['RNA_SEQ','ATAC'] (default ['RNA_SEQ'])."
        },
        "ontology_terms": {
          "type": [
            "array",
            "null"
          ],
          "items": {
            "type": "string"
          },
          "description": "Optional tissue/cell ontology terms (e.g. ['UBERON:0001114'] = liver) to restrict tracks."
        },
        "organism": {
          "type": [
            "string",
            "null"
          ],
          "description": "'human' (default) or 'mouse'."
        },
        "sequence_length": {
          "type": [
            "string",
            "null"
          ],
          "description": "Context window: 16KB, 100KB, 500KB, or 1MB (default)."
        }
      }
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "interval": {
              "type": "string"
            },
            "tracks": {
              "type": "array",
              "items": {
                "type": "object"
              }
            }
          }
        },
        {
          "type": "object",
          "properties": {
            "status": {
              "const": "error"
            },
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

**Execution and provenance:** `network`；upstream name `AlphaGenome_predict_interval`。

**Provider credentials:** `ALPHA_GENOME_API_KEY`

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("alpha_genome_predict_interval", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.alpha_genome_score_variant` {#dnakit.databases.alpha_genome_score_variant}

Score a regulatory variant's effect with DeepMind AlphaGenome (Avsec, Nature 2026), the hosted successor to Enformer/Borzoi. Predicts ref vs alt over up to 1 Mb at single-base resolution and returns the recommended per-track effect scores (sorted by |effect|). Hosted API via the alphagenome SDK; free for non-commercial use (set ALPHA_GENOME_API_KEY).

```text
dnakit.databases.alpha_genome_score_variant(*, chromosome, position, reference_bases, alternate_bases, output_type=None, organism=None, sequence_length=None, top_n=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `chromosome` | string | Required by input schema | Chromosome, e.g. &#x27;chr22&#x27;. |
| `position` | integer | Required by input schema | 1-based variant position. |
| `reference_bases` | string | Required by input schema | Reference allele, e.g. &#x27;A&#x27;. |
| `alternate_bases` | string | Required by input schema | Alternate allele, e.g. &#x27;C&#x27;. |
| `output_type` | string &#124; null | None | Modality to score: RNA_SEQ (default), ATAC, DNASE, CAGE, CHIP_HISTONE, CHIP_TF, SPLICE_SITES, CONTACT_MAPS. |
| `organism` | string &#124; null | None | &#x27;human&#x27; (default) or &#x27;mouse&#x27;. |
| `sequence_length` | string &#124; null | None | Context window: 16KB, 100KB, 500KB, or 1MB (default). |
| `top_n` | integer &#124; null | None | Number of top &#124;effect&#124; tracks to return (default 20). |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "required": [
        "alternate_bases",
        "chromosome",
        "position",
        "reference_bases"
      ],
      "properties": {
        "chromosome": {
          "type": "string",
          "description": "Chromosome, e.g. 'chr22'."
        },
        "position": {
          "type": "integer",
          "description": "1-based variant position."
        },
        "reference_bases": {
          "type": "string",
          "description": "Reference allele, e.g. 'A'."
        },
        "alternate_bases": {
          "type": "string",
          "description": "Alternate allele, e.g. 'C'."
        },
        "output_type": {
          "type": [
            "string",
            "null"
          ],
          "description": "Modality to score: RNA_SEQ (default), ATAC, DNASE, CAGE, CHIP_HISTONE, CHIP_TF, SPLICE_SITES, CONTACT_MAPS."
        },
        "organism": {
          "type": [
            "string",
            "null"
          ],
          "description": "'human' (default) or 'mouse'."
        },
        "sequence_length": {
          "type": [
            "string",
            "null"
          ],
          "description": "Context window: 16KB, 100KB, 500KB, or 1MB (default)."
        },
        "top_n": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Number of top |effect| tracks to return (default 20)."
        }
      }
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
              "type": "string"
            },
            "output_type": {
              "type": "string"
            },
            "scores": {
              "type": "array",
              "items": {
                "type": "object"
              }
            }
          }
        },
        {
          "type": "object",
          "properties": {
            "status": {
              "const": "error"
            },
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

**Execution and provenance:** `network`；upstream name `AlphaGenome_score_variant`。

**Provider credentials:** `ALPHA_GENOME_API_KEY`

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("alpha_genome_score_variant", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.dyna_mut2_predict_stability` {#dnakit.databases.dyna_mut2_predict_stability}

Predict the effect of a single-point amino acid mutation on protein stability and dynamics using DynaMut2. Provide a PDB ID (4-character code from RCSB PDB), chain identifier, and mutation in the format WtPosNew (e.g., V1A means Valine at position 1 mutated to Alanine). Returns the predicted change in Gibbs free energy (ddG in kcal/mol): positive = stabilizing, negative = destabilizing. The tool downloads the PDB structure from RCSB and submits it to DynaMut2 for analysis. Computation typically takes 15-30 seconds. Based on Normal Mode Analysis and graph-based signatures.

```text
dnakit.databases.dyna_mut2_predict_stability(*, operation=None, pdb_id, chain, mutation, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;predict_stability&quot;&#93;, &quot;default&quot;: &quot;predict_stability&quot;} |
| `pdb_id` | string | Required by input schema | 4-character PDB accession code from RCSB PDB (e.g., &#x27;4HHB&#x27; for hemoglobin, &#x27;1U46&#x27; for ACK1 kinase). The structure will be downloaded automatically. |
| `chain` | string | Required by input schema | Chain identifier in the PDB structure (e.g., &#x27;A&#x27;). Must be a valid chain in the specified PDB entry. |
| `mutation` | string | Required by input schema | Single-point mutation in the format WtPosNew where Wt is the wild-type amino acid (1-letter code), Pos is the residue number, and New is the mutant amino acid (1-letter code). Examples: &#x27;V1A&#x27; (Val1Ala), &#x27;E346K&#x27; (Glu346Lys), &#x27;G12D&#x27; (Gly12Asp). |
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
          "enum": [
            "predict_stability"
          ],
          "description": "Operation type",
          "default": "predict_stability"
        },
        "pdb_id": {
          "type": "string",
          "description": "4-character PDB accession code from RCSB PDB (e.g., '4HHB' for hemoglobin, '1U46' for ACK1 kinase). The structure will be downloaded automatically."
        },
        "chain": {
          "type": "string",
          "description": "Chain identifier in the PDB structure (e.g., 'A'). Must be a valid chain in the specified PDB entry."
        },
        "mutation": {
          "type": "string",
          "description": "Single-point mutation in the format WtPosNew where Wt is the wild-type amino acid (1-letter code), Pos is the residue number, and New is the mutant amino acid (1-letter code). Examples: 'V1A' (Val1Ala), 'E346K' (Glu346Lys), 'G12D' (Gly12Asp)."
        }
      },
      "required": [
        "chain",
        "mutation",
        "pdb_id"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "pdb_id": {
          "type": "string",
          "description": "PDB accession code used"
        },
        "chain": {
          "type": "string",
          "description": "Chain identifier"
        },
        "mutation": {
          "type": "string",
          "description": "Mutation in WtPosNew format"
        },
        "wild_type": {
          "type": "string",
          "description": "Wild-type residue (3-letter code)"
        },
        "mutant": {
          "type": "string",
          "description": "Mutant residue (3-letter code)"
        },
        "position": {
          "type": "string",
          "description": "Residue position number"
        },
        "ddg_prediction": {
          "type": "number",
          "description": "Predicted change in Gibbs free energy (ddG, kcal/mol). Positive = stabilizing, negative = destabilizing."
        },
        "effect": {
          "type": "string",
          "description": "Predicted effect: 'stabilizing' or 'destabilizing'"
        },
        "job_id": {
          "type": "string",
          "description": "DynaMut2 job identifier for retrieving results later"
        },
        "results_page": {
          "type": "string",
          "description": "URL to the full DynaMut2 results page with visualizations"
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `DynaMut2_predict_stability`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("dyna_mut2_predict_stability", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.esm2_score_missense_variant` {#dnakit.databases.esm2_score_missense_variant}

Score a missense protein variant with ESM-2 masked-marginal log-likelihood ratio (Meier 2021) over HuggingFace's keyless hf-inference. Input wild-type 1-letter sequence + 1-based position + mutant residue; returns LLR = logP(mut) - logP(wt). Negative => mutant disfavored (candidate deleterious). No API key. For ESMC/SAE scoring with a key use ESM_score_sequence.

```text
dnakit.databases.esm2_score_missense_variant(*, sequence, position, mutant, wild_type=None, model_id=None, wait_for_model=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `sequence` | string | Required by input schema | Wild-type protein sequence in 1-letter amino-acid code (whitespace is ignored). |
| `position` | integer | Required by input schema | 1-based residue position of the variant within the sequence. |
| `mutant` | string | Required by input schema | Mutant amino acid (single 1-letter code, one of the 20 standard residues). |
| `wild_type` | string &#124; null | None | Optional wild-type residue (1-letter). If given, it is validated against the residue at `position` to catch coordinate/sequence mismatches. |
| `model_id` | string &#124; null | None | ESM-2 model to use. Default &#x27;facebook/esm2_t33_650M_UR50D&#x27;. Other sizes: esm2_t30_150M_UR50D, esm2_t12_35M_UR50D, esm2_t6_8M_UR50D. |
| `wait_for_model` | boolean &#124; null | None | If true, block until the model finishes loading on the HF servers instead of returning a retryable &#x27;loading&#x27; status. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "required": [
        "mutant",
        "position",
        "sequence"
      ],
      "properties": {
        "sequence": {
          "type": "string",
          "description": "Wild-type protein sequence in 1-letter amino-acid code (whitespace is ignored)."
        },
        "position": {
          "type": "integer",
          "description": "1-based residue position of the variant within the sequence."
        },
        "mutant": {
          "type": "string",
          "description": "Mutant amino acid (single 1-letter code, one of the 20 standard residues)."
        },
        "wild_type": {
          "type": [
            "string",
            "null"
          ],
          "description": "Optional wild-type residue (1-letter). If given, it is validated against the residue at `position` to catch coordinate/sequence mismatches."
        },
        "model_id": {
          "type": [
            "string",
            "null"
          ],
          "description": "ESM-2 model to use. Default 'facebook/esm2_t33_650M_UR50D'. Other sizes: esm2_t30_150M_UR50D, esm2_t12_35M_UR50D, esm2_t6_8M_UR50D."
        },
        "wait_for_model": {
          "type": [
            "boolean",
            "null"
          ],
          "description": "If true, block until the model finishes loading on the HF servers instead of returning a retryable 'loading' status."
        }
      }
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
              "type": "string"
            },
            "wild_type": {
              "type": "string"
            },
            "mutant": {
              "type": "string"
            },
            "p_wild_type": {
              "type": "number"
            },
            "p_mutant": {
              "type": "number"
            },
            "log_likelihood_ratio": {
              "type": "number"
            },
            "direction": {
              "type": "string"
            }
          }
        },
        {
          "type": "object",
          "properties": {
            "status": {
              "const": "error"
            },
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

**Execution and provenance:** `network`；upstream name `ESM2_score_missense_variant`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("esm2_score_missense_variant", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.esm_describe_sae_feature` {#dnakit.databases.esm_describe_sae_feature}

Label a single SAE feature_id with its dominant biological category by aggregating UniProt feature-type overlap at SAE-activating residues across a curated 10-protein panel. Returns one of: catalytic, ligand-binding, ptm, domain, motif, structural-stability, secondary-structure, transmembrane, signal-peptide, propeptide, or 'uncategorized' (no informative UniProt overlap found). Use when ESM_get_sae_features / ESM_score_variant_sae_disruption returned a feature_id that needs biological interpretation. Cost: 10 Forge credits FIRST time (1 per panel protein). Subsequent calls for the same feature_id are FREE (filesystem cache at ~/.cache/tooluniverse/sae_labels/{sae_model}/feature_{id}.json). Latency: ~30-60s first call (panel of 10 SAE inferences), <1s on cache hit. Prerequisites: same as ESM_get_sae_features. License: Cambrian Inference License — non-commercial.

```text
dnakit.databases.esm_describe_sae_feature(*, feature_id, sae_model=None, model=None, n_proteins=None, top_residues_per_protein=None, use_cache=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `feature_id` | integer | Required by input schema | SAE feature index in &#91;0, 16383&#93;. Use the values returned by ESM_get_sae_features.active_features&#91;&#93;.feature_id. |
| `sae_model` | string | None | SAE checkpoint to label. Cache is keyed on this — different SAEs produce different labels for the same id. {&quot;default&quot;: &quot;esmc-6b-2024-12_k64_codebook16384_layer60&quot;} |
| `model` | string | None | ESMC backbone. {&quot;enum&quot;: &#91;&quot;esmc-6b-2024-12&quot;&#93;, &quot;default&quot;: &quot;esmc-6b-2024-12&quot;} |
| `n_proteins` | integer | None | Number of panel proteins to run SAE on. Default 10 (the full curated panel). Lower for cheaper exploration at the cost of confidence. {&quot;default&quot;: 10} |
| `top_residues_per_protein` | integer | None | For each protein, take the top-K residues where the target feature activates most strongly, then check their UniProt annotations. Default 3. {&quot;default&quot;: 3} |
| `use_cache` | boolean | None | If true (default) and a cached label exists at ~/.cache/tooluniverse/sae_labels/.../feature_{id}.json, return it instead of rerunning the panel. Set false to force a recompute. {&quot;default&quot;: true} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "feature_id": {
          "type": "integer",
          "description": "SAE feature index in [0, 16383]. Use the values returned by ESM_get_sae_features.active_features[].feature_id."
        },
        "sae_model": {
          "type": "string",
          "description": "SAE checkpoint to label. Cache is keyed on this — different SAEs produce different labels for the same id.",
          "default": "esmc-6b-2024-12_k64_codebook16384_layer60"
        },
        "model": {
          "type": "string",
          "description": "ESMC backbone.",
          "enum": [
            "esmc-6b-2024-12"
          ],
          "default": "esmc-6b-2024-12"
        },
        "n_proteins": {
          "type": "integer",
          "description": "Number of panel proteins to run SAE on. Default 10 (the full curated panel). Lower for cheaper exploration at the cost of confidence.",
          "default": 10
        },
        "top_residues_per_protein": {
          "type": "integer",
          "description": "For each protein, take the top-K residues where the target feature activates most strongly, then check their UniProt annotations. Default 3.",
          "default": 3
        },
        "use_cache": {
          "type": "boolean",
          "description": "If true (default) and a cached label exists at ~/.cache/tooluniverse/sae_labels/.../feature_{id}.json, return it instead of rerunning the panel. Set false to force a recompute.",
          "default": true
        }
      },
      "required": [
        "feature_id"
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
            "feature_id": {
              "type": "integer"
            },
            "sae_model": {
              "type": "string"
            },
            "category": {
              "type": "string",
              "description": "Dominant high-level category derived from UniProt feature-type voting across the panel"
            },
            "confidence": {
              "type": "number",
              "description": "Fraction of votes for the dominant category (0-1)"
            },
            "n_proteins_with_activation": {
              "type": "integer"
            },
            "n_proteins_analyzed": {
              "type": "integer"
            },
            "category_vote_counts": {
              "type": "object",
              "description": "Raw vote counts per category, useful for spotting ties / weak signals"
            },
            "supporting_evidence": {
              "type": "array",
              "description": "Per-residue evidence rows: which panel protein, which position, activation magnitude, which UniProt feature types overlapped",
              "items": {
                "type": "object",
                "properties": {
                  "protein": {
                    "type": "string"
                  },
                  "position_1based": {
                    "type": "integer"
                  },
                  "activation": {
                    "type": "number"
                  },
                  "uniprot_types": {
                    "type": "array",
                    "items": {
                      "type": "string"
                    }
                  },
                  "uniprot_categories": {
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
            "feature_id"
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

**Execution and provenance:** `model`；upstream name `ESM_describe_sae_feature`。

**Provider credentials:** `ESM_API_KEY`

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("esm_describe_sae_feature", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.esm_explain_variant_mechanism` {#dnakit.databases.esm_explain_variant_mechanism}

One-call composite for variant mechanism: runs ESMC-6B SAE variant disruption + describe_sae_feature on each top affected feature + composes a category-level mechanism summary (e.g. 'Disrupted feature categories (lost): catalytic=2, ligand-binding=1'). Use when a variant-interpretation skill wants the mechanism in a single tool call rather than orchestrating disruption + N description calls. Set include_descriptions=false to skip the category labeling step (2 Forge calls only). Otherwise: 2 disruption calls + up to (2 * top_k_features) describe calls (most cached after first use). Prerequisites: pip install 'esm @ git+https://github.com/evolutionaryscale/esm@ee891c52'; set ESM_API_KEY. Outputs governed by EvolutionaryScale Cambrian Inference License — non-commercial use only unless covered by a separate commercial agreement.

```text
dnakit.databases.esm_explain_variant_mechanism(*, sequence, position, ref_aa, alt_aa, window=None, top_k_features=None, include_descriptions=None, model=None, sae_model=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `sequence` | string | Required by input schema | Reference (wild-type) protein amino acid sequence in single-letter code. Up to ~2700 AA. |
| `position` | integer | Required by input schema | 1-indexed residue position of the variant |
| `ref_aa` | string | Required by input schema | Single-letter wild-type amino acid (must match sequence&#91;position-1&#93;) |
| `alt_aa` | string | Required by input schema | Single-letter substituted amino acid |
| `window` | integer | None | Residue window centered on the mutation for activation summation. {&quot;default&quot;: 8} |
| `top_k_features` | integer | None | Number of top lost / top gained features to describe and include in the category summary. {&quot;default&quot;: 5} |
| `include_descriptions` | boolean | None | If true (default), call ESM_describe_sae_feature on each top feature to get category labels. If false, return raw feature IDs without categories (cheaper). {&quot;default&quot;: true} |
| `model` | string | None | ESMC base model {&quot;enum&quot;: &#91;&quot;esmc-6b-2024-12&quot;&#93;, &quot;default&quot;: &quot;esmc-6b-2024-12&quot;} |
| `sae_model` | string | None | SAE codebook identifier {&quot;default&quot;: &quot;esmc-6b-2024-12_k64_codebook16384_layer60&quot;} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "sequence": {
          "type": "string",
          "description": "Reference (wild-type) protein amino acid sequence in single-letter code. Up to ~2700 AA."
        },
        "position": {
          "type": "integer",
          "description": "1-indexed residue position of the variant"
        },
        "ref_aa": {
          "type": "string",
          "description": "Single-letter wild-type amino acid (must match sequence[position-1])"
        },
        "alt_aa": {
          "type": "string",
          "description": "Single-letter substituted amino acid"
        },
        "window": {
          "type": "integer",
          "description": "Residue window centered on the mutation for activation summation.",
          "default": 8
        },
        "top_k_features": {
          "type": "integer",
          "description": "Number of top lost / top gained features to describe and include in the category summary.",
          "default": 5
        },
        "include_descriptions": {
          "type": "boolean",
          "description": "If true (default), call ESM_describe_sae_feature on each top feature to get category labels. If false, return raw feature IDs without categories (cheaper).",
          "default": true
        },
        "model": {
          "type": "string",
          "enum": [
            "esmc-6b-2024-12"
          ],
          "default": "esmc-6b-2024-12",
          "description": "ESMC base model"
        },
        "sae_model": {
          "type": "string",
          "description": "SAE codebook identifier",
          "default": "esmc-6b-2024-12_k64_codebook16384_layer60"
        }
      },
      "required": [
        "alt_aa",
        "position",
        "ref_aa",
        "sequence"
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
              "type": "string"
            },
            "position": {
              "type": "integer"
            },
            "window": {
              "type": "integer"
            },
            "mechanism_summary": {
              "type": "string"
            },
            "lost_feature_categories": {
              "type": "object"
            },
            "gained_feature_categories": {
              "type": "object"
            },
            "top_features_lost": {
              "type": "array"
            },
            "top_features_gained": {
              "type": "array"
            }
          },
          "required": [
            "variant"
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

**Execution and provenance:** `model`；upstream name `ESM_explain_variant_mechanism`。

**Provider credentials:** `ESM_API_KEY`

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("esm_explain_variant_mechanism", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.esm_get_region_sae_features` {#dnakit.databases.esm_get_region_sae_features}

Aggregate ESMC-6B SAE features over a contiguous residue range to characterize the region's biological identity. Returns features ranked by total |activation| across the region, with each feature's per-residue hit pattern. Use for: interpreting a known domain (residues X-Y), profiling an antibody epitope, characterizing a binding pocket, comparing two regions of the same protein, or seeding ESM_describe_sae_feature calls on the top-K to get category labels. Max region length 500 residues. Prerequisites: pip install 'esm @ git+https://github.com/evolutionaryscale/esm@ee891c52'; set ESM_API_KEY. Outputs governed by EvolutionaryScale Cambrian Inference License — non-commercial use only unless covered by a separate commercial agreement.

```text
dnakit.databases.esm_get_region_sae_features(*, sequence, start_position, end_position, top_k_features=None, model=None, sae_model=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `sequence` | string | Required by input schema | Protein amino acid sequence in single-letter code. Up to ~2700 AA. |
| `start_position` | integer | Required by input schema | 1-indexed inclusive start of the region of interest. |
| `end_position` | integer | Required by input schema | 1-indexed inclusive end of the region of interest. Must be &gt;= start_position and within sequence length. Region length cap: 500. |
| `top_k_features` | integer | None | Number of top features (by total &#124;activation&#124; over region) to return. {&quot;default&quot;: 20} |
| `model` | string | None | ESMC base model {&quot;enum&quot;: &#91;&quot;esmc-6b-2024-12&quot;&#93;, &quot;default&quot;: &quot;esmc-6b-2024-12&quot;} |
| `sae_model` | string | None | SAE codebook identifier {&quot;default&quot;: &quot;esmc-6b-2024-12_k64_codebook16384_layer60&quot;} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "sequence": {
          "type": "string",
          "description": "Protein amino acid sequence in single-letter code. Up to ~2700 AA."
        },
        "start_position": {
          "type": "integer",
          "description": "1-indexed inclusive start of the region of interest."
        },
        "end_position": {
          "type": "integer",
          "description": "1-indexed inclusive end of the region of interest. Must be >= start_position and within sequence length. Region length cap: 500."
        },
        "top_k_features": {
          "type": "integer",
          "description": "Number of top features (by total |activation| over region) to return.",
          "default": 20
        },
        "model": {
          "type": "string",
          "enum": [
            "esmc-6b-2024-12"
          ],
          "default": "esmc-6b-2024-12",
          "description": "ESMC base model"
        },
        "sae_model": {
          "type": "string",
          "description": "SAE codebook identifier",
          "default": "esmc-6b-2024-12_k64_codebook16384_layer60"
        }
      },
      "required": [
        "end_position",
        "sequence",
        "start_position"
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
            "sequence_length": {
              "type": "integer"
            },
            "region": {
              "type": "array"
            },
            "region_length": {
              "type": "integer"
            },
            "n_features_active_in_region": {
              "type": "integer"
            },
            "top_features": {
              "type": "array"
            }
          },
          "required": [
            "sequence_length"
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

**Execution and provenance:** `model`；upstream name `ESM_get_region_sae_features`。

**Provider credentials:** `ESM_API_KEY`

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("esm_get_region_sae_features", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.esm_get_sae_features` {#dnakit.databases.esm_get_sae_features}

Run a protein sequence through an ESMC Sparse Autoencoder (SAE) and return sparse feature activations per residue. SAEs decompose the protein language model's hidden states into a 16,384-feature sparse codebook (top-k=64 active per residue). Each feature is an interpretable latent — many activate on specific biological categories (catalytic site, ligand-binding region, PTM sequon, structural motif, etc.) that can be looked up via ESM_describe_sae_feature. Use for: variant interpretation (compare ref vs mutant activations), mechanism inference (which functional features the protein engages), protein-LM interpretability research. Prerequisites: (1) pip install 'esm @ git+https://github.com/evolutionaryscale/esm@ee891c52' — the PyPI esm release does NOT yet include SAE support, which lives on an upstream feature branch. (2) ESM_API_KEY env var with a Forge token (https://forge.evolutionaryscale.ai). License note: SAE outputs are governed by the EvolutionaryScale Cambrian Inference Clickthrough License — non-commercial / academic use only.

```text
dnakit.databases.esm_get_sae_features(*, sequence, model=None, sae_model=None, position=None, window=None, top_k_per_residue=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `sequence` | string | Required by input schema | Protein amino acid sequence in single-letter code (e.g. &#x27;MEEPQSDPSVEPPLSQETFSDLWKLLPENN&#x27;). Standard 20 amino acids; avoid gaps or non-standard characters. Practical length cap ~2,700 AA. |
| `model` | string | None | ESMC backbone model. Currently SAE is only released for esmc-6b-2024-12 (default). {&quot;enum&quot;: &#91;&quot;esmc-6b-2024-12&quot;&#93;, &quot;default&quot;: &quot;esmc-6b-2024-12&quot;} |
| `sae_model` | string | None | SAE checkpoint name. Default is the layer-60 SAE matching esmc-6b backbone. {&quot;default&quot;: &quot;esmc-6b-2024-12_k64_codebook16384_layer60&quot;} |
| `position` | integer &#124; null | None | Optional 1-indexed residue position. If set, only activations within +/- window residues of this position are returned (use for variant analysis to keep response compact). If null, all residues are returned (large response for long sequences). {&quot;default&quot;: null} |
| `window` | integer | None | Residue window radius around position (only used when position is set). Default 8 = +/- 8 residues = 17-residue window centered on the variant. {&quot;default&quot;: 8} |
| `top_k_per_residue` | integer | None | Cap features returned per residue, sorted by absolute activation. Default 64 = full sparsity (no cap). Set lower (e.g. 16) to reduce response size when no position filter is used. {&quot;default&quot;: 64} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "sequence": {
          "type": "string",
          "description": "Protein amino acid sequence in single-letter code (e.g. 'MEEPQSDPSVEPPLSQETFSDLWKLLPENN'). Standard 20 amino acids; avoid gaps or non-standard characters. Practical length cap ~2,700 AA."
        },
        "model": {
          "type": "string",
          "description": "ESMC backbone model. Currently SAE is only released for esmc-6b-2024-12 (default).",
          "enum": [
            "esmc-6b-2024-12"
          ],
          "default": "esmc-6b-2024-12"
        },
        "sae_model": {
          "type": "string",
          "description": "SAE checkpoint name. Default is the layer-60 SAE matching esmc-6b backbone.",
          "default": "esmc-6b-2024-12_k64_codebook16384_layer60"
        },
        "position": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Optional 1-indexed residue position. If set, only activations within +/- window residues of this position are returned (use for variant analysis to keep response compact). If null, all residues are returned (large response for long sequences).",
          "default": null
        },
        "window": {
          "type": "integer",
          "description": "Residue window radius around position (only used when position is set). Default 8 = +/- 8 residues = 17-residue window centered on the variant.",
          "default": 8
        },
        "top_k_per_residue": {
          "type": "integer",
          "description": "Cap features returned per residue, sorted by absolute activation. Default 64 = full sparsity (no cap). Set lower (e.g. 16) to reduce response size when no position filter is used.",
          "default": 64
        }
      },
      "required": [
        "sequence"
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
            "sequence_length": {
              "type": "integer"
            },
            "model": {
              "type": "string"
            },
            "sae_model": {
              "type": "string"
            },
            "residues_returned": {
              "type": "integer"
            },
            "position": {
              "type": [
                "integer",
                "null"
              ]
            },
            "window": {
              "type": [
                "integer",
                "null"
              ]
            },
            "top_k_per_residue": {
              "type": "integer"
            },
            "activations": {
              "type": "array",
              "description": "Per-residue active SAE features. Sparse format: each entry has the residue position (1-indexed) and its top-K non-zero features.",
              "items": {
                "type": "object",
                "properties": {
                  "residue_idx_1based": {
                    "type": "integer"
                  },
                  "active_features": {
                    "type": "array",
                    "items": {
                      "type": "object",
                      "properties": {
                        "feature_id": {
                          "type": "integer",
                          "description": "SAE feature index 0-16383"
                        },
                        "activation": {
                          "type": "number",
                          "description": "Activation magnitude (normalized)"
                        }
                      }
                    }
                  }
                }
              }
            }
          },
          "required": [
            "sequence_length"
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

**Execution and provenance:** `model`；upstream name `ESM_get_sae_features`。

**Provider credentials:** `ESM_API_KEY`

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("esm_get_sae_features", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.esm_score_sequence` {#dnakit.databases.esm_score_sequence}

Score a protein sequence using ESMC logits to compute per-residue log-probabilities and mean pseudo-log-likelihood. Higher (less negative) mean log-probability indicates the sequence is more natural/likely according to the model's learned protein distribution. Use for: evaluating mutant fitness, comparing sequence variants, identifying low-confidence residues. Mean log-probability near 0 is highest confidence; values below -4 suggest unusual or low-quality residues. Prerequisites: pip install esm; set ESM_API_KEY (from https://forge.evolutionaryscale.ai).

```text
dnakit.databases.esm_score_sequence(*, sequence, model=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `sequence` | string | Required by input schema | Protein amino acid sequence in single-letter code to score (e.g. &#x27;MKTAYIAKQRQISFVKSHFSRQ&#x27;). |
| `model` | string | None | ESMC model to use for scoring. {&quot;enum&quot;: &#91;&quot;esmc-300m-2024-12&quot;, &quot;esmc-600m-2024-12&quot;&#93;, &quot;default&quot;: &quot;esmc-300m-2024-12&quot;} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "sequence": {
          "type": "string",
          "description": "Protein amino acid sequence in single-letter code to score (e.g. 'MKTAYIAKQRQISFVKSHFSRQ')."
        },
        "model": {
          "type": "string",
          "description": "ESMC model to use for scoring.",
          "enum": [
            "esmc-300m-2024-12",
            "esmc-600m-2024-12"
          ],
          "default": "esmc-300m-2024-12"
        }
      },
      "required": [
        "sequence"
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
            "status": {
              "type": "string",
              "enum": [
                "success"
              ]
            },
            "model": {
              "type": "string"
            },
            "sequence": {
              "type": "string"
            },
            "sequence_length": {
              "type": "integer"
            },
            "mean_log_probability": {
              "type": [
                "number",
                "null"
              ],
              "description": "Mean per-residue log-probability (pseudo-log-likelihood). Higher values indicate more natural sequences."
            },
            "per_residue_log_probabilities": {
              "type": "array",
              "items": {
                "type": "number"
              },
              "description": "Per-residue log-probabilities for each amino acid position"
            }
          },
          "required": [
            "status"
          ]
        },
        {
          "type": "object",
          "properties": {
            "status": {
              "type": "string",
              "enum": [
                "error"
              ]
            },
            "error": {
              "type": "string"
            }
          },
          "required": [
            "status",
            "error"
          ]
        }
      ]
    }
    ```

**Execution and provenance:** `model`；upstream name `ESM_score_sequence`。

**Provider credentials:** `ESM_API_KEY`

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("esm_score_sequence", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.esm_score_variant_sae_batch` {#dnakit.databases.esm_score_variant_sae_batch}

Score many missense variants against one reference protein using ESMC-6B SAE. More Forge-efficient than calling ESM_score_variant_sae_disruption N times: runs the reference SAE once and reuses it across every variant, giving N+1 Forge calls instead of 2N. Use for saturation mutagenesis at a single position (all 19 alternative residues), DMS-style sweeps, clinical variant panels, or any workflow scoring multiple variants on the same protein. Max 100 variants per call. Prerequisites: pip install 'esm @ git+https://github.com/evolutionaryscale/esm@ee891c52' (SAE support is on an unmerged feature branch — PyPI esm 3.2.x does NOT include SAEConfig). Set ESM_API_KEY environment variable. Outputs governed by EvolutionaryScale Cambrian Inference License — non-commercial use only unless covered by a separate commercial agreement.

```text
dnakit.databases.esm_score_variant_sae_batch(*, sequence, variants, window=None, top_k_features=None, model=None, sae_model=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `sequence` | string | Required by input schema | Reference (wild-type) protein amino acid sequence in single-letter code. Up to ~2700 AA. |
| `variants` | array | Required by input schema | List of variants to score. Each variant: {position: 1-indexed int, ref_aa: single letter matching sequence&#91;position-1&#93;, alt_aa: single letter substitution}. {&quot;maxItems&quot;: 100, &quot;items&quot;: {&quot;type&quot;: &quot;object&quot;, &quot;properties&quot;: {&quot;position&quot;: {&quot;type&quot;: &quot;integer&quot;, &quot;description&quot;: &quot;1-indexed residue position&quot;}, &quot;ref_aa&quot;: {&quot;type&quot;: &quot;string&quot;, &quot;description&quot;: &quot;Single-letter wild-type amino acid (must match sequence&#91;position-1&#93;)&quot;}, &quot;alt_aa&quot;: {&quot;type&quot;: &quot;string&quot;, &quot;description&quot;: &quot;Single-letter substituted amino acid&quot;}}, &quot;required&quot;: &#91;&quot;position&quot;, &quot;ref_aa&quot;, &quot;alt_aa&quot;&#93;}} |
| `window` | integer | None | Residue window centered on each variant for activation summation. Default 8 (i.e. positions &#91;pos-8, pos+8&#93;). {&quot;default&quot;: 8} |
| `top_k_features` | integer | None | Number of top features to return per variant (top-K lost + top-K gained). {&quot;default&quot;: 10} |
| `model` | string | None | ESMC base model (must be the 6B model — SAEs are trained against it) {&quot;enum&quot;: &#91;&quot;esmc-6b-2024-12&quot;&#93;, &quot;default&quot;: &quot;esmc-6b-2024-12&quot;} |
| `sae_model` | string | None | SAE codebook identifier {&quot;default&quot;: &quot;esmc-6b-2024-12_k64_codebook16384_layer60&quot;} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "sequence": {
          "type": "string",
          "description": "Reference (wild-type) protein amino acid sequence in single-letter code. Up to ~2700 AA."
        },
        "variants": {
          "type": "array",
          "description": "List of variants to score. Each variant: {position: 1-indexed int, ref_aa: single letter matching sequence[position-1], alt_aa: single letter substitution}.",
          "maxItems": 100,
          "items": {
            "type": "object",
            "properties": {
              "position": {
                "type": "integer",
                "description": "1-indexed residue position"
              },
              "ref_aa": {
                "type": "string",
                "description": "Single-letter wild-type amino acid (must match sequence[position-1])"
              },
              "alt_aa": {
                "type": "string",
                "description": "Single-letter substituted amino acid"
              }
            },
            "required": [
              "position",
              "ref_aa",
              "alt_aa"
            ]
          }
        },
        "window": {
          "type": "integer",
          "description": "Residue window centered on each variant for activation summation. Default 8 (i.e. positions [pos-8, pos+8]).",
          "default": 8
        },
        "top_k_features": {
          "type": "integer",
          "description": "Number of top features to return per variant (top-K lost + top-K gained).",
          "default": 10
        },
        "model": {
          "type": "string",
          "enum": [
            "esmc-6b-2024-12"
          ],
          "default": "esmc-6b-2024-12",
          "description": "ESMC base model (must be the 6B model — SAEs are trained against it)"
        },
        "sae_model": {
          "type": "string",
          "description": "SAE codebook identifier",
          "default": "esmc-6b-2024-12_k64_codebook16384_layer60"
        }
      },
      "required": [
        "sequence",
        "variants"
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
            "sequence_length": {
              "type": "integer"
            },
            "n_variants": {
              "type": "integer"
            },
            "n_succeeded": {
              "type": "integer"
            },
            "window": {
              "type": "integer"
            },
            "results": {
              "type": "array"
            }
          },
          "required": [
            "sequence_length"
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

**Execution and provenance:** `model`；upstream name `ESM_score_variant_sae_batch`。

**Provider credentials:** `ESM_API_KEY`

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("esm_score_variant_sae_batch", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.esm_score_variant_sae_disruption` {#dnakit.databases.esm_score_variant_sae_disruption}

Composite SAE-based variant scoring. Given a protein sequence and a missense variant (position + ref_aa + alt_aa), runs ESMC-6B SAE on both reference and mutant sequences, computes per-feature activation deltas summed over a window centered on the mutation site, and returns ranked top-K features LOST and GAINED in the mutant. This is the convenience layer over ESM_get_sae_features — use this for one-shot variant interpretation; use ESM_get_sae_features directly only when you need raw per-residue features for non-variant analyses. Validates ref_aa matches the position in the supplied sequence (returns clear error if not, so you can detect transcript / isoform mismatches). Cost: 2 Forge API credits (1 ref + 1 mut). Latency: ~3-6s for typical-length human proteins. Prerequisites: same as ESM_get_sae_features (pip install 'esm @ git+https://github.com/evolutionaryscale/esm@ee891c52', ESM_API_KEY env var). License: SAE outputs governed by Cambrian Inference License — non-commercial / academic only.

```text
dnakit.databases.esm_score_variant_sae_disruption(*, sequence, position, ref_aa, alt_aa, window=None, top_k_features=None, model=None, sae_model=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `sequence` | string | Required by input schema | Reference protein sequence (canonical isoform). Single-letter codes, no gaps. The mutant sequence is built internally by substituting alt_aa at position. |
| `position` | integer | Required by input schema | 1-indexed mutation position. The amino acid at sequence&#91;position-1&#93; must equal ref_aa; otherwise the tool returns an error (the wrong sequence was supplied or the variant notation references a different isoform). |
| `ref_aa` | string | Required by input schema | Reference amino acid at the mutation position, single-letter code (e.g. &#x27;R&#x27;). Tool validates this matches the supplied sequence. |
| `alt_aa` | string | Required by input schema | Mutant amino acid, single-letter code (e.g. &#x27;H&#x27; for R175H). |
| `window` | integer | None | Residue window radius around the mutation. Per-feature activations are summed across this window before computing the delta. Default 8 = +/- 8 residues = 17-residue window. {&quot;default&quot;: 8} |
| `top_k_features` | integer | None | Number of top LOST and top GAINED features to return. Default 10. {&quot;default&quot;: 10} |
| `model` | string | None | ESMC backbone (default esmc-6b-2024-12). {&quot;enum&quot;: &#91;&quot;esmc-6b-2024-12&quot;&#93;, &quot;default&quot;: &quot;esmc-6b-2024-12&quot;} |
| `sae_model` | string | None | SAE checkpoint (default layer-60 6B SAE). {&quot;default&quot;: &quot;esmc-6b-2024-12_k64_codebook16384_layer60&quot;} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "sequence": {
          "type": "string",
          "description": "Reference protein sequence (canonical isoform). Single-letter codes, no gaps. The mutant sequence is built internally by substituting alt_aa at position."
        },
        "position": {
          "type": "integer",
          "description": "1-indexed mutation position. The amino acid at sequence[position-1] must equal ref_aa; otherwise the tool returns an error (the wrong sequence was supplied or the variant notation references a different isoform)."
        },
        "ref_aa": {
          "type": "string",
          "description": "Reference amino acid at the mutation position, single-letter code (e.g. 'R'). Tool validates this matches the supplied sequence."
        },
        "alt_aa": {
          "type": "string",
          "description": "Mutant amino acid, single-letter code (e.g. 'H' for R175H)."
        },
        "window": {
          "type": "integer",
          "description": "Residue window radius around the mutation. Per-feature activations are summed across this window before computing the delta. Default 8 = +/- 8 residues = 17-residue window.",
          "default": 8
        },
        "top_k_features": {
          "type": "integer",
          "description": "Number of top LOST and top GAINED features to return. Default 10.",
          "default": 10
        },
        "model": {
          "type": "string",
          "description": "ESMC backbone (default esmc-6b-2024-12).",
          "enum": [
            "esmc-6b-2024-12"
          ],
          "default": "esmc-6b-2024-12"
        },
        "sae_model": {
          "type": "string",
          "description": "SAE checkpoint (default layer-60 6B SAE).",
          "default": "esmc-6b-2024-12_k64_codebook16384_layer60"
        }
      },
      "required": [
        "alt_aa",
        "position",
        "ref_aa",
        "sequence"
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
              "type": "string",
              "description": "e.g. 'R175H'"
            },
            "position": {
              "type": "integer"
            },
            "window": {
              "type": "integer"
            },
            "n_unique_features_touched": {
              "type": "integer"
            },
            "top_features_lost": {
              "type": "array",
              "description": "Features with most negative delta (strong in ref, weak/absent in mutant)",
              "items": {
                "type": "object",
                "properties": {
                  "feature_id": {
                    "type": "integer"
                  },
                  "delta": {
                    "type": "number",
                    "description": "mut_activation_sum - ref_activation_sum (negative)"
                  },
                  "ref_activation_sum": {
                    "type": "number"
                  },
                  "mut_activation_sum": {
                    "type": "number"
                  }
                }
              }
            },
            "top_features_gained": {
              "type": "array",
              "description": "Features with most positive delta (absent/weak in ref, strong in mutant)",
              "items": {
                "type": "object",
                "properties": {
                  "feature_id": {
                    "type": "integer"
                  },
                  "delta": {
                    "type": "number"
                  },
                  "ref_activation_sum": {
                    "type": "number"
                  },
                  "mut_activation_sum": {
                    "type": "number"
                  }
                }
              }
            }
          },
          "required": [
            "variant"
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

**Execution and provenance:** `model`；upstream name `ESM_score_variant_sae_disruption`。

**Provider credentials:** `ESM_API_KEY`

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("esm_score_variant_sae_disruption", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.evo2_score_variant` {#dnakit.databases.evo2_score_variant}

Zero-shot variant-effect scoring with NVIDIA-hosted Evo 2 (Arc Institute genome foundation model). Runs the model's forward pass on the reference and alternate DNA windows and returns the delta log-likelihood (alt - ref). Negative = variant disfavored by the genome model (candidate deleterious). Two hosted sizes are selectable via `model`: evo2-40b (default, most accurate) and evo2-7b (smaller, faster). Hosted NIM API; requires NVIDIA_API_KEY. Complements NvidiaNIM_evo2 (which only generates sequences).

```text
dnakit.databases.evo2_score_variant(*, ref_sequence=None, alt_sequence=None, sequence=None, position=None, reference=None, alternate=None, model=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `ref_sequence` | string &#124; null | None | Reference DNA window (A/C/G/T/N). Use with alt_sequence (same length). |
| `alt_sequence` | string &#124; null | None | Alternate DNA window, same length/centering as ref_sequence. |
| `sequence` | string &#124; null | None | Reference DNA window for point-substitution mode (use with position + alternate). |
| `position` | integer &#124; null | None | 1-based position of the substituted base within `sequence`. |
| `reference` | string &#124; null | None | Optional reference base (single letter) at `position`, validated against `sequence`. |
| `alternate` | string &#124; null | None | Alternate base (single letter) substituted at `position`. |
| `model` | string | None | Hosted Evo 2 model size: evo2-40b (default, most accurate) or evo2-7b (smaller, faster). {&quot;enum&quot;: &#91;&quot;evo2-40b&quot;, &quot;evo2-7b&quot;&#93;, &quot;default&quot;: &quot;evo2-40b&quot;} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "ref_sequence": {
          "type": [
            "string",
            "null"
          ],
          "description": "Reference DNA window (A/C/G/T/N). Use with alt_sequence (same length)."
        },
        "alt_sequence": {
          "type": [
            "string",
            "null"
          ],
          "description": "Alternate DNA window, same length/centering as ref_sequence."
        },
        "sequence": {
          "type": [
            "string",
            "null"
          ],
          "description": "Reference DNA window for point-substitution mode (use with position + alternate)."
        },
        "position": {
          "type": [
            "integer",
            "null"
          ],
          "description": "1-based position of the substituted base within `sequence`."
        },
        "reference": {
          "type": [
            "string",
            "null"
          ],
          "description": "Optional reference base (single letter) at `position`, validated against `sequence`."
        },
        "alternate": {
          "type": [
            "string",
            "null"
          ],
          "description": "Alternate base (single letter) substituted at `position`."
        },
        "model": {
          "type": "string",
          "enum": [
            "evo2-40b",
            "evo2-7b"
          ],
          "default": "evo2-40b",
          "description": "Hosted Evo 2 model size: evo2-40b (default, most accurate) or evo2-7b (smaller, faster)."
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
            "delta_loglik": {
              "type": "number"
            },
            "ref_loglik": {
              "type": "number"
            },
            "alt_loglik": {
              "type": "number"
            },
            "direction": {
              "type": "string"
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

**Execution and provenance:** `network`；upstream name `Evo2_score_variant`。

**Provider credentials:** `NVIDIA_API_KEY`

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("evo2_score_variant", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.nvidia_nim_alphafold2` {#dnakit.databases.nvidia_nim_alphafold2}

Predict protein 3D structure from amino acid sequence using DeepMind's AlphaFold2 via NVIDIA NIM. Returns PDB format structure. Requires NVIDIA_API_KEY. Input: sequence (amino acids), algorithm (jackhmmer/mmseqs2), databases (uniref90, mgnify, small_bfd). Long-running async operation.

```text
dnakit.databases.nvidia_nim_alphafold2(*, sequence, algorithm=None, databases=None, e_value=None, iterations=None, relax_prediction=None, skip_template_search=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `sequence` | string | Required by input schema | Amino acid sequence to predict structure for (single letter codes) |
| `algorithm` | string | None | MSA search algorithm. mmseqs2 is faster, jackhmmer is more sensitive {&quot;enum&quot;: &#91;&quot;jackhmmer&quot;, &quot;mmseqs2&quot;&#93;, &quot;default&quot;: &quot;mmseqs2&quot;} |
| `databases` | array | None | Sequence databases for MSA search: uniref90, mgnify, small_bfd {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}, &quot;default&quot;: &#91;&quot;small_bfd&quot;&#93;} |
| `e_value` | number | None | E-value threshold for MSA search {&quot;default&quot;: 0.0001} |
| `iterations` | integer | None | Number of search iterations {&quot;default&quot;: 1} |
| `relax_prediction` | boolean | None | Whether to perform structure relaxation {&quot;default&quot;: false} |
| `skip_template_search` | boolean | None | Skip template-based prediction {&quot;default&quot;: true} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "sequence": {
          "type": "string",
          "description": "Amino acid sequence to predict structure for (single letter codes)"
        },
        "algorithm": {
          "type": "string",
          "enum": [
            "jackhmmer",
            "mmseqs2"
          ],
          "description": "MSA search algorithm. mmseqs2 is faster, jackhmmer is more sensitive",
          "default": "mmseqs2"
        },
        "databases": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "Sequence databases for MSA search: uniref90, mgnify, small_bfd",
          "default": [
            "small_bfd"
          ]
        },
        "e_value": {
          "type": "number",
          "description": "E-value threshold for MSA search",
          "default": 0.0001
        },
        "iterations": {
          "type": "integer",
          "description": "Number of search iterations",
          "default": 1
        },
        "relax_prediction": {
          "type": "boolean",
          "description": "Whether to perform structure relaxation",
          "default": false
        },
        "skip_template_search": {
          "type": "boolean",
          "description": "Skip template-based prediction",
          "default": true
        }
      },
      "required": [
        "sequence"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "description": "AlphaFold2 structure prediction result with PDB coordinates, confidence scores, and model metrics",
      "additionalProperties": true
    }
    ```

**Execution and provenance:** `network`；upstream name `NvidiaNIM_alphafold2`。

**Provider credentials:** `NVIDIA_API_KEY`

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("nvidia_nim_alphafold2", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.nvidia_nim_esmfold` {#dnakit.databases.nvidia_nim_esmfold}

Fast alignment-free protein structure prediction using ESMFold via NVIDIA NIM. No MSA required - predicts directly from sequence. Input: sequence (max 1024 amino acids). Returns PDB structure. Fast synchronous response.

```text
dnakit.databases.nvidia_nim_esmfold(*, sequence, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `sequence` | string | Required by input schema | Amino acid sequence (max 1024 chars, valid: ARNDCQEGHILKMFPSTWYVXBOU) {&quot;maxLength&quot;: 1024} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "sequence": {
          "type": "string",
          "description": "Amino acid sequence (max 1024 chars, valid: ARNDCQEGHILKMFPSTWYVXBOU)",
          "maxLength": 1024
        }
      },
      "required": [
        "sequence"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": [
        "string",
        "null"
      ],
      "description": "PDB format structure text (returned via 'structure' field, not 'data'). Contains ATOM records with coordinates and pLDDT confidence in B-factor column."
    }
    ```

**Execution and provenance:** `network`；upstream name `NvidiaNIM_esmfold`。

**Provider credentials:** `NVIDIA_API_KEY`

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("nvidia_nim_esmfold", verbose=False)
print(info["input_schema"])
```

