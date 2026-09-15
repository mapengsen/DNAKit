# Sequence analysis, VCF, cloning and inheritance

This page documents **15 newly added public functions**. [Complete index](index.md) · [Conventions and examples](usage.md)

Each function lists its current Python signature, parameters, returns and registered execution requirements. Provider descriptions retain source wording. In generated database wrappers, None normally omits a parameter; actual required fields and conditional constraints come from the input schema. Follow each operation's coordinate and assembly conventions.

## Function index

- [`dnakit.databases.coding_variant_fraction`](#dnakit.databases.coding_variant_fraction)
- [`dnakit.databases.dna_gibson_design`](#dnakit.databases.dna_gibson_design)
- [`dnakit.databases.dna_golden_gate_assemble`](#dnakit.databases.dna_golden_gate_assemble)
- [`dnakit.databases.dna_golden_gate_design`](#dnakit.databases.dna_golden_gate_design)
- [`dnakit.databases.dna_translate_reading_frames`](#dnakit.databases.dna_translate_reading_frames)
- [`dnakit.databases.dna_virtual_digest`](#dnakit.databases.dna_virtual_digest)
- [`dnakit.databases.sequence_count_residues`](#dnakit.databases.sequence_count_residues)
- [`dnakit.databases.sequence_dn_ds`](#dnakit.databases.sequence_dn_ds)
- [`dnakit.databases.sequence_gc_content`](#dnakit.databases.sequence_gc_content)
- [`dnakit.databases.sequence_reverse_complement`](#dnakit.databases.sequence_reverse_complement)
- [`dnakit.databases.sequence_stats`](#dnakit.databases.sequence_stats)
- [`dnakit.databases.vcf_count_variants`](#dnakit.databases.vcf_count_variants)
- [`dnakit.databases.vcf_normalize`](#dnakit.databases.vcf_normalize)
- [`dnakit.databases.vcf_summary_stats`](#dnakit.databases.vcf_summary_stats)
- [`dnakit.workflows.simulate_inheritance`](#dnakit.workflows.simulate_inheritance)

## `dnakit.databases.coding_variant_fraction` {#dnakit.databases.coding_variant_fraction}

Compute the fraction of CODING variants in an Excel/CSV variant table that match a given Sequence Ontology annotation and fall below a VAF threshold. Uses the CODING-only denominator convention (synonymous_variant, missense_variant, splice_region_variant, stop_gained, stop_lost, start_lost, frameshift_variant, inframe_insertion, inframe_deletion) - intronic, UTR, intergenic, and up/downstream records are excluded. Auto-detects 'variant allele frequency'/'VAF' and 'Sequence Ontology' columns if not specified. Supports flat headers (header_rows=1) and MultiIndex headers (header_rows=2) common in cancer variant reports. Returns the fraction with coding denominator plus the 'naive' fraction (over all VAF-filtered records) for comparison.

```text
dnakit.databases.coding_variant_fraction(*, file, vaf_threshold=None, annotation=None, header_rows=None, vaf_col=None, so_col=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `file` | string | Required by input schema | Path to variant Excel (.xlsx) or CSV file. |
| `vaf_threshold` | number | None | Keep variants with VAF strictly less than this value. Default 0.3. |
| `annotation` | string | None | Sequence Ontology term to count in numerator. Default &#x27;synonymous_variant&#x27;. |
| `header_rows` | integer | None | 1 for flat headers, 2 for MultiIndex (common in Excel variant reports). Default 1. {&quot;enum&quot;: &#91;1, 2&#93;} |
| `vaf_col` | string &#124; null | None | Explicit VAF column name. If empty, auto-detect from &#x27;variant allele freq&#x27; or &#x27;vaf&#x27;. |
| `so_col` | string &#124; null | None | Explicit Sequence Ontology column name. If empty, auto-detect from &#x27;sequence ontology&#x27;. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "file": {
          "type": "string",
          "description": "Path to variant Excel (.xlsx) or CSV file."
        },
        "vaf_threshold": {
          "type": "number",
          "description": "Keep variants with VAF strictly less than this value. Default 0.3."
        },
        "annotation": {
          "type": "string",
          "description": "Sequence Ontology term to count in numerator. Default 'synonymous_variant'."
        },
        "header_rows": {
          "type": "integer",
          "enum": [
            1,
            2
          ],
          "description": "1 for flat headers, 2 for MultiIndex (common in Excel variant reports). Default 1."
        },
        "vaf_col": {
          "type": [
            "string",
            "null"
          ],
          "description": "Explicit VAF column name. If empty, auto-detect from 'variant allele freq' or 'vaf'."
        },
        "so_col": {
          "type": [
            "string",
            "null"
          ],
          "description": "Explicit Sequence Ontology column name. If empty, auto-detect from 'sequence ontology'."
        }
      },
      "required": [
        "file"
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
            "total_variants": {
              "type": "integer"
            },
            "vaf_filtered": {
              "type": "integer"
            },
            "coding_variants_below_vaf": {
              "type": "integer"
            },
            "matching_annotation": {
              "type": "integer"
            },
            "fraction": {
              "type": "number"
            },
            "naive_fraction": {
              "type": "number"
            },
            "vaf_column": {
              "type": "string"
            },
            "so_column": {
              "type": "string"
            },
            "vaf_threshold": {
              "type": "number"
            },
            "annotation": {
              "type": "string"
            }
          },
          "required": [
            "total_variants",
            "coding_variants_below_vaf",
            "matching_annotation",
            "fraction"
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

**Execution and provenance:** `network`；upstream name `coding_variant_fraction`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("coding_variant_fraction", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.dna_gibson_design` {#dnakit.databases.dna_gibson_design}

Design Gibson Assembly overlaps for seamless DNA fragment assembly. For each fragment, computes the 5' overhang (from the previous fragment's 3' end) and 3' overhang (from the next fragment's 5' end) based on the specified overlap length (default 20 bp). Returns the modified fragment sequences ready for Gibson Assembly cloning. Use for: seamless cloning, multi-fragment assembly, vector backbone insertion.

```text
dnakit.databases.dna_gibson_design(*, operation=None, fragments, overlap_length=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation (optional; defaults to &#x27;gibson_design&#x27; for this tool). {&quot;enum&quot;: &#91;&quot;gibson_design&quot;&#93;, &quot;default&quot;: &quot;gibson_design&quot;} |
| `fragments` | array | Required by input schema | List of DNA fragment sequences (at least 2). Each must be longer than overlap_length. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `overlap_length` | integer &#124; null | None | Overlap length in bp for Gibson Assembly (default: 20 bp). Must be at least 1. {&quot;default&quot;: 20} |
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
            "gibson_design"
          ],
          "description": "Operation (optional; defaults to 'gibson_design' for this tool).",
          "default": "gibson_design"
        },
        "fragments": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "List of DNA fragment sequences (at least 2). Each must be longer than overlap_length."
        },
        "overlap_length": {
          "type": [
            "integer",
            "null"
          ],
          "default": 20,
          "description": "Overlap length in bp for Gibson Assembly (default: 20 bp). Must be at least 1."
        }
      },
      "required": [
        "fragments"
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
            "assembly_fragments": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "name": {
                    "type": "string"
                  },
                  "original_sequence": {
                    "type": "string"
                  },
                  "with_overlaps": {
                    "type": "string"
                  },
                  "left_overlap": {
                    "type": "string"
                  },
                  "right_overlap": {
                    "type": "string"
                  }
                }
              }
            },
            "n_fragments": {
              "type": "integer"
            },
            "overlap_length": {
              "type": "integer"
            }
          },
          "required": [
            "assembly_fragments",
            "n_fragments",
            "overlap_length"
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

**Execution and provenance:** `network`；upstream name `DNA_gibson_design`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("dna_gibson_design", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.dna_golden_gate_assemble` {#dnakit.databases.dna_golden_gate_assemble}

Work out the product of a Golden Gate reaction that has already been performed: given the plasmids/parts that were combined and the Type IIS enzyme used (BsaI, BbsI, Esp3I/BsmBI, SapI), digest each input, discard fragments that keep a recognition site (they are re-cut in the reaction), and ligate the rest by matching 4-bp overhangs into the circular product. Returns the assembled sequence, its length, and the fragment order with the overhang at each junction. Use to answer 'I combined these plasmids with <enzyme> — what does the resulting plasmid contain/express?'. This is the reverse of DNA_golden_gate_design, which designs parts rather than analysing a finished reaction.

```text
dnakit.databases.dna_golden_gate_assemble(*, operation=None, fragments, enzyme=None, circular=None, labels=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation (optional; defaults to &#x27;golden_gate_assemble&#x27; for this tool). {&quot;enum&quot;: &#91;&quot;golden_gate_assemble&quot;&#93;, &quot;default&quot;: &quot;golden_gate_assemble&quot;} |
| `fragments` | array | Required by input schema | The DNA sequences combined in the reaction (plasmids or linear parts), A/T/G/C only. At least 2. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `enzyme` | string &#124; null | None | Type IIS enzyme used in the reaction, e.g. &#x27;BsaI&#x27; (default), &#x27;BbsI&#x27;, &#x27;Esp3I&#x27;/&#x27;BsmBI&#x27;, &#x27;SapI&#x27;. {&quot;default&quot;: &quot;BsaI&quot;} |
| `circular` | boolean &#124; null | None | Whether the inputs are circular plasmids (default true). Set false for linear parts. {&quot;default&quot;: true} |
| `labels` | array &#124; null | None | Optional names for the inputs (e.g. plasmid names), echoed in the assembly order. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
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
            "golden_gate_assemble"
          ],
          "description": "Operation (optional; defaults to 'golden_gate_assemble' for this tool).",
          "default": "golden_gate_assemble"
        },
        "fragments": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "The DNA sequences combined in the reaction (plasmids or linear parts), A/T/G/C only. At least 2."
        },
        "enzyme": {
          "type": [
            "string",
            "null"
          ],
          "default": "BsaI",
          "description": "Type IIS enzyme used in the reaction, e.g. 'BsaI' (default), 'BbsI', 'Esp3I'/'BsmBI', 'SapI'."
        },
        "circular": {
          "type": [
            "boolean",
            "null"
          ],
          "default": true,
          "description": "Whether the inputs are circular plasmids (default true). Set false for linear parts."
        },
        "labels": {
          "type": [
            "array",
            "null"
          ],
          "items": {
            "type": "string"
          },
          "description": "Optional names for the inputs (e.g. plasmid names), echoed in the assembly order."
        }
      },
      "required": [
        "fragments"
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
            "enzyme": {
              "type": "string"
            },
            "product_length": {
              "type": "integer",
              "description": "Length of the assembled circular product in bp."
            },
            "num_fragments_assembled": {
              "type": "integer"
            },
            "assembly_order": {
              "type": "array",
              "description": "Fragments in assembly order, each with source, length and its two overhangs."
            },
            "unused_fragments": {
              "type": "array",
              "description": "Site-free fragments that did not enter the circular product."
            },
            "product_sequence": {
              "type": "string"
            }
          },
          "required": [
            "enzyme",
            "product_length",
            "num_fragments_assembled",
            "assembly_order",
            "product_sequence"
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

**Execution and provenance:** `network`；upstream name `DNA_golden_gate_assemble`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("dna_golden_gate_assemble", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.dna_golden_gate_design` {#dnakit.databases.dna_golden_gate_design}

Design Golden Gate Assembly parts with Type IIS restriction enzyme sites (BsaI or BbsI). Assigns unique non-palindromic 4-bp overhangs to each junction and builds the full synthetic part sequences including enzyme recognition sites. Use for: scarless multi-part assembly, modular cloning (MoClo), Golden Gate cloning workflows.

```text
dnakit.databases.dna_golden_gate_design(*, operation=None, parts, enzyme=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation (optional; defaults to &#x27;golden_gate_design&#x27; for this tool). {&quot;enum&quot;: &#91;&quot;golden_gate_design&quot;&#93;, &quot;default&quot;: &quot;golden_gate_design&quot;} |
| `parts` | array | Required by input schema | List of DNA part sequences to assemble (at least 2). A, T, G, C only. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `enzyme` | string &#124; null | None | Type IIS restriction enzyme to use. Options: &#x27;BsaI&#x27; (default) or &#x27;BbsI&#x27;. {&quot;enum&quot;: &#91;&quot;BsaI&quot;, &quot;BbsI&quot;, null&#93;, &quot;default&quot;: &quot;BsaI&quot;} |
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
            "golden_gate_design"
          ],
          "description": "Operation (optional; defaults to 'golden_gate_design' for this tool).",
          "default": "golden_gate_design"
        },
        "parts": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "List of DNA part sequences to assemble (at least 2). A, T, G, C only."
        },
        "enzyme": {
          "type": [
            "string",
            "null"
          ],
          "enum": [
            "BsaI",
            "BbsI",
            null
          ],
          "default": "BsaI",
          "description": "Type IIS restriction enzyme to use. Options: 'BsaI' (default) or 'BbsI'."
        }
      },
      "required": [
        "parts"
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
            "parts_with_overhangs": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "name": {
                    "type": "string"
                  },
                  "sequence": {
                    "type": "string"
                  },
                  "left_overhang": {
                    "type": "string"
                  },
                  "right_overhang": {
                    "type": "string"
                  },
                  "full_sequence": {
                    "type": "string"
                  }
                }
              }
            },
            "overhangs": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "enzyme": {
              "type": "string"
            },
            "protocol_note": {
              "type": "string"
            }
          },
          "required": [
            "parts_with_overhangs",
            "overhangs",
            "enzyme",
            "protocol_note"
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

**Execution and provenance:** `network`；upstream name `DNA_golden_gate_design`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("dna_golden_gate_design", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.dna_translate_reading_frames` {#dnakit.databases.dna_translate_reading_frames}

Translate a DNA sequence to protein in all 3 forward reading frames and identify the longest ORF. Uses the standard genetic code (NCBI Code 1). Returns protein sequences for each frame, longest ORF per frame, and the best frame with longest uninterrupted protein. Use for: reading frame analysis, finding the correct ORF in an unknown sequence, comparing all three frames.

```text
dnakit.databases.dna_translate_reading_frames(*, operation=None, sequence, frame=None, genetic_code=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;translate_reading_frames&quot;&#93;, &quot;default&quot;: &quot;translate_reading_frames&quot;} |
| `sequence` | string | Required by input schema | DNA sequence (A, T, G, C only). Case insensitive. Spaces and newlines are stripped. |
| `frame` | string &#124; null | None | Reading frame to translate: &#x27;1&#x27;, &#x27;2&#x27;, &#x27;3&#x27;, or &#x27;all&#x27; (default). Frame 1 starts at position 0, frame 2 at position 1, frame 3 at position 2. {&quot;enum&quot;: &#91;&quot;1&quot;, &quot;2&quot;, &quot;3&quot;, &quot;all&quot;, null&#93;, &quot;default&quot;: &quot;all&quot;} |
| `genetic_code` | integer &#124; null | None | NCBI genetic code number. Only 1 (the standard code) is implemented; any other value is rejected with an error rather than translated with the standard table. Alternative code tables (2 = vertebrate mitochondrial, 3 = yeast mitochondrial, 11 = bacterial, etc.) are not supported by this tool. {&quot;enum&quot;: &#91;1, null&#93;, &quot;default&quot;: 1} |
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
            "translate_reading_frames"
          ],
          "description": "Operation type",
          "default": "translate_reading_frames"
        },
        "sequence": {
          "type": "string",
          "description": "DNA sequence (A, T, G, C only). Case insensitive. Spaces and newlines are stripped."
        },
        "frame": {
          "type": [
            "string",
            "null"
          ],
          "enum": [
            "1",
            "2",
            "3",
            "all",
            null
          ],
          "default": "all",
          "description": "Reading frame to translate: '1', '2', '3', or 'all' (default). Frame 1 starts at position 0, frame 2 at position 1, frame 3 at position 2."
        },
        "genetic_code": {
          "type": [
            "integer",
            "null"
          ],
          "enum": [
            1,
            null
          ],
          "default": 1,
          "description": "NCBI genetic code number. Only 1 (the standard code) is implemented; any other value is rejected with an error rather than translated with the standard table. Alternative code tables (2 = vertebrate mitochondrial, 3 = yeast mitochondrial, 11 = bacterial, etc.) are not supported by this tool."
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
            "frames": {
              "type": "object",
              "description": "Translation results per frame"
            },
            "best_frame": {
              "type": "integer"
            },
            "best_orf": {
              "type": "string"
            },
            "best_orf_length_aa": {
              "type": "integer"
            }
          },
          "required": [
            "sequence_length",
            "frames",
            "best_frame"
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

**Execution and provenance:** `network`；upstream name `DNA_translate_reading_frames`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("dna_translate_reading_frames", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.dna_virtual_digest` {#dnakit.databases.dna_virtual_digest}

Perform a virtual restriction enzyme digest on a DNA sequence using the NEB enzyme library (25 curated enzymes plus Biopython's ~600-enzyme library, including Type IIS enzymes such as BsaI, BbsI, Esp3I/BsmBI and SapI, which cut outside their recognition site). Returns all resulting fragments with their sequences, lengths, and positions, plus a list of cut sites. Supports both linear and circular DNA. Use for: cloning planning, gel electrophoresis prediction, verifying construct integrity.

```text
dnakit.databases.dna_virtual_digest(*, operation=None, sequence, enzymes=None, circular=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation (optional; defaults to &#x27;virtual_digest&#x27; for this tool). {&quot;enum&quot;: &#91;&quot;virtual_digest&quot;&#93;, &quot;default&quot;: &quot;virtual_digest&quot;} |
| `sequence` | string | Required by input schema | DNA sequence (A, T, G, C, N only). Case insensitive. |
| `enzymes` | array &#124; null | None | List of enzyme names to use for digestion (e.g., &#91;&#x27;EcoRI&#x27;, &#x27;BamHI&#x27;&#93;). If null, uses all 24 NEB enzymes. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `circular` | boolean &#124; null | None | Treat sequence as circular DNA (e.g., plasmid). Default: false (linear). {&quot;default&quot;: false} |
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
            "virtual_digest"
          ],
          "description": "Operation (optional; defaults to 'virtual_digest' for this tool).",
          "default": "virtual_digest"
        },
        "sequence": {
          "type": "string",
          "description": "DNA sequence (A, T, G, C, N only). Case insensitive."
        },
        "enzymes": {
          "type": [
            "array",
            "null"
          ],
          "items": {
            "type": "string"
          },
          "description": "List of enzyme names to use for digestion (e.g., ['EcoRI', 'BamHI']). If null, uses all 24 NEB enzymes."
        },
        "circular": {
          "type": [
            "boolean",
            "null"
          ],
          "default": false,
          "description": "Treat sequence as circular DNA (e.g., plasmid). Default: false (linear)."
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
            "fragments": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "sequence": {
                    "type": "string"
                  },
                  "length": {
                    "type": "integer"
                  },
                  "start": {
                    "type": "integer"
                  },
                  "end": {
                    "type": "integer"
                  }
                }
              }
            },
            "n_fragments": {
              "type": "integer"
            },
            "enzymes_used": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "cut_sites": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "enzyme": {
                    "type": "string"
                  },
                  "position": {
                    "type": "integer"
                  }
                }
              }
            }
          },
          "required": [
            "fragments",
            "n_fragments",
            "enzymes_used",
            "cut_sites"
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

**Execution and provenance:** `network`；upstream name `DNA_virtual_digest`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("dna_virtual_digest", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.sequence_count_residues` {#dnakit.databases.sequence_count_residues}

Count occurrences of a specific amino acid or nucleotide residue in a sequence. Accepts either a raw sequence string or a UniProt accession (fetches protein sequence from UniProt REST API). Reports count, fraction, percentage, and 1-based positions. Useful for cysteine counting (disulfide bond estimation), charged residue analysis, or nucleotide composition.

```text
dnakit.databases.sequence_count_residues(*, operation=None, sequence=None, uniprot_id=None, residue, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;count_residues&quot;&#93;, &quot;default&quot;: &quot;count_residues&quot;} |
| `sequence` | string | None | Sequence string (DNA/RNA/protein). Optional if uniprot_id is provided. |
| `uniprot_id` | string | None | UniProt accession to fetch sequence from (e.g., &#x27;P00533&#x27; for EGFR). Optional if sequence is provided. |
| `residue` | string | Required by input schema | Single character residue to count (e.g., &#x27;C&#x27; for cysteine) |
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
            "count_residues"
          ],
          "description": "Operation type",
          "default": "count_residues"
        },
        "sequence": {
          "type": "string",
          "description": "Sequence string (DNA/RNA/protein). Optional if uniprot_id is provided."
        },
        "uniprot_id": {
          "type": "string",
          "description": "UniProt accession to fetch sequence from (e.g., 'P00533' for EGFR). Optional if sequence is provided."
        },
        "residue": {
          "type": "string",
          "description": "Single character residue to count (e.g., 'C' for cysteine)"
        }
      },
      "required": [
        "residue"
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
            "data": {
              "type": "object",
              "properties": {
                "count": {
                  "type": "integer"
                },
                "sequence_length": {
                  "type": "integer"
                },
                "fraction": {
                  "type": "number"
                },
                "percent": {
                  "type": "number"
                },
                "positions_1based": {
                  "type": "array"
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

**Execution and provenance:** `network`；upstream name `Sequence_count_residues`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("sequence_count_residues", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.sequence_dn_ds` {#dnakit.databases.sequence_dn_ds}

Compute dN/dS (Ka/Ks) between two coding sequences via the Nei-Gojobori (1986) estimator with Jukes-Cantor correction — the standard test for selection: dN/dS > 1 = positive/diversifying, << 1 = purifying, ~ 1 = near-neutral. Give two in-frame, codon-aligned coding sequences inline ('seq1'/'seq2') or as single-record FASTA files. Returns dN, dS, the ratio, per-site counts, and an interpretation. Pure Python, no API key. Works for any organisms/genes.

```text
dnakit.databases.sequence_dn_ds(*, seq1=None, seq2=None, fasta1_path=None, fasta2_path=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `seq1` | string &#124; null | None | First coding sequence (in-frame, codon-aligned to seq2). DNA/RNA letters. |
| `seq2` | string &#124; null | None | Second coding sequence (same length/frame as seq1). |
| `fasta1_path` | string &#124; null | None | Alternative to seq1: path to a single-record FASTA |
| `fasta2_path` | string &#124; null | None | Alternative to seq2: path to a single-record FASTA |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "required": [],
      "properties": {
        "seq1": {
          "type": [
            "string",
            "null"
          ],
          "description": "First coding sequence (in-frame, codon-aligned to seq2). DNA/RNA letters."
        },
        "seq2": {
          "type": [
            "string",
            "null"
          ],
          "description": "Second coding sequence (same length/frame as seq1)."
        },
        "fasta1_path": {
          "type": [
            "string",
            "null"
          ],
          "description": "Alternative to seq1: path to a single-record FASTA"
        },
        "fasta2_path": {
          "type": [
            "string",
            "null"
          ],
          "description": "Alternative to seq2: path to a single-record FASTA"
        }
      }
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "description": "Nei-Gojobori (1986) dN/dS estimate between two coding sequences.",
      "properties": {
        "dN_dS": {
          "type": "number",
          "description": "dN/dS (Ka/Ks) ratio."
        },
        "dN": {
          "type": "number",
          "description": "Nonsynonymous substitution rate."
        },
        "dS": {
          "type": "number",
          "description": "Synonymous substitution rate."
        },
        "pN": {
          "type": "number",
          "description": "Proportion of nonsynonymous differences."
        },
        "pS": {
          "type": "number",
          "description": "Proportion of synonymous differences."
        },
        "N_sites": {
          "type": "number",
          "description": "Number of nonsynonymous sites."
        },
        "S_sites": {
          "type": "number",
          "description": "Number of synonymous sites."
        },
        "Nd": {
          "type": "number",
          "description": "Observed number of nonsynonymous differences."
        },
        "Sd": {
          "type": "number",
          "description": "Observed number of synonymous differences."
        },
        "codons_compared": {
          "type": "integer",
          "description": "Number of codons compared."
        },
        "interpretation": {
          "type": "string",
          "description": "Human-readable interpretation of the selection signal."
        }
      },
      "required": [
        "dN_dS",
        "dN",
        "dS",
        "codons_compared"
      ]
    }
    ```

**Execution and provenance:** `network`；upstream name `Sequence_dn_ds`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("sequence_dn_ds", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.sequence_gc_content` {#dnakit.databases.sequence_gc_content}

Calculate GC content of a DNA or RNA sequence. Reports GC count, GC fraction, GC percentage, and full nucleotide composition. Ambiguous bases (N) are excluded from the denominator. Useful for primer design, genome characterization, and codon usage analysis.

```text
dnakit.databases.sequence_gc_content(*, operation=None, sequence, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;gc_content&quot;&#93;, &quot;default&quot;: &quot;gc_content&quot;} |
| `sequence` | string | Required by input schema | DNA or RNA sequence string |
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
            "gc_content"
          ],
          "description": "Operation type",
          "default": "gc_content"
        },
        "sequence": {
          "type": "string",
          "description": "DNA or RNA sequence string"
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
            "data": {
              "type": "object",
              "properties": {
                "gc_fraction": {
                  "type": "number"
                },
                "gc_percent": {
                  "type": "number"
                },
                "length": {
                  "type": "integer"
                },
                "composition": {
                  "type": "object"
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

**Execution and provenance:** `network`；upstream name `Sequence_gc_content`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("sequence_gc_content", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.sequence_reverse_complement` {#dnakit.databases.sequence_reverse_complement}

Return the reverse complement of a DNA sequence. Supports A, T, C, G, and N bases. Useful for primer design, probe design, and strand orientation analysis.

```text
dnakit.databases.sequence_reverse_complement(*, operation=None, sequence, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;reverse_complement&quot;&#93;, &quot;default&quot;: &quot;reverse_complement&quot;} |
| `sequence` | string | Required by input schema | DNA sequence (A/T/C/G/N only) |
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
            "reverse_complement"
          ],
          "description": "Operation type",
          "default": "reverse_complement"
        },
        "sequence": {
          "type": "string",
          "description": "DNA sequence (A/T/C/G/N only)"
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
            "data": {
              "type": "object",
              "properties": {
                "original": {
                  "type": "string"
                },
                "reverse_complement": {
                  "type": "string"
                },
                "length": {
                  "type": "integer"
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

**Execution and provenance:** `network`；upstream name `Sequence_reverse_complement`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("sequence_reverse_complement", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.sequence_stats` {#dnakit.databases.sequence_stats}

Compute basic statistics for a DNA, RNA, or protein sequence. Auto-detects sequence type, reports length, full composition, and type-specific metrics (GC% for nucleotides, estimated monoisotopic molecular weight for proteins -- for the conventional average mass reported by ExPASy ProtParam/UniProt, use ProtParam_calculate instead). Accepts sequence string or UniProt accession.

```text
dnakit.databases.sequence_stats(*, operation=None, sequence=None, uniprot_id=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;stats&quot;&#93;, &quot;default&quot;: &quot;stats&quot;} |
| `sequence` | string | None | Sequence string. Optional if uniprot_id is provided. |
| `uniprot_id` | string | None | UniProt accession to fetch protein sequence from. Optional if sequence is provided. |
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
            "stats"
          ],
          "description": "Operation type",
          "default": "stats"
        },
        "sequence": {
          "type": "string",
          "description": "Sequence string. Optional if uniprot_id is provided."
        },
        "uniprot_id": {
          "type": "string",
          "description": "UniProt accession to fetch protein sequence from. Optional if sequence is provided."
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
            "data": {
              "type": "object",
              "properties": {
                "sequence_type": {
                  "type": "string"
                },
                "length": {
                  "type": "integer"
                },
                "composition": {
                  "type": "object"
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

**Execution and provenance:** `network`；upstream name `Sequence_stats`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("sequence_stats", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.vcf_count_variants` {#dnakit.databases.vcf_count_variants}

Count variants in a local VCF/BCF after applying filters, via bcftools. Returns total records plus SNP/indel/MNP breakdown. Supports PASS-only, region, min QUAL, and a custom bcftools include expression so SNP/indel counts under a filter are deterministic. Requires bcftools on PATH.

```text
dnakit.databases.vcf_count_variants(*, operation=None, vcf_path, pass_only=None, regions=None, min_qual=None, include_expr=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | object &#124; None | None | Operation (fixed) {&quot;const&quot;: &quot;count_variants&quot;} |
| `vcf_path` | string | Required by input schema | Path to a local .vcf, .vcf.gz, or .bcf file |
| `pass_only` | boolean &#124; null | None | If true, count only FILTER=PASS or &#x27;.&#x27; records |
| `regions` | string &#124; null | None | Restrict to a region, e.g. &#x27;chr1&#x27; or &#x27;chr1:1-1000000&#x27; |
| `min_qual` | number &#124; null | None | Keep only records with QUAL &gt;= this value |
| `include_expr` | string &#124; null | None | Extra bcftools -i include expression, e.g. &#x27;INFO/AF&gt;0.01&#x27; |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "required": [
        "vcf_path"
      ],
      "properties": {
        "operation": {
          "const": "count_variants",
          "description": "Operation (fixed)"
        },
        "vcf_path": {
          "type": "string",
          "description": "Path to a local .vcf, .vcf.gz, or .bcf file"
        },
        "pass_only": {
          "type": [
            "boolean",
            "null"
          ],
          "description": "If true, count only FILTER=PASS or '.' records"
        },
        "regions": {
          "type": [
            "string",
            "null"
          ],
          "description": "Restrict to a region, e.g. 'chr1' or 'chr1:1-1000000'"
        },
        "min_qual": {
          "type": [
            "number",
            "null"
          ],
          "description": "Keep only records with QUAL >= this value"
        },
        "include_expr": {
          "type": [
            "string",
            "null"
          ],
          "description": "Extra bcftools -i include expression, e.g. 'INFO/AF>0.01'"
        }
      }
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

**Execution and provenance:** `external_process`；upstream name `VCF_count_variants`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("vcf_count_variants", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.vcf_normalize` {#dnakit.databases.vcf_normalize}

Normalize a local VCF/BCF with bcftools norm (split or join multiallelics; left-align indels against a reference FASTA when given) and report SNP/indel/record counts before vs after. Resolves the common cause of disagreeing indel counts: un-split multiallelics and non-left-aligned indels. Requires bcftools on PATH.

```text
dnakit.databases.vcf_normalize(*, operation=None, vcf_path, multiallelics=None, reference_fasta=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | object &#124; None | None | Operation (fixed) {&quot;const&quot;: &quot;normalize&quot;} |
| `vcf_path` | string | Required by input schema | Path to a local .vcf, .vcf.gz, or .bcf file |
| `multiallelics` | string &#124; null | None | &#x27;split&#x27; (default) breaks multiallelic records into biallelic; &#x27;join&#x27; merges; &#x27;none&#x27; leaves as-is {&quot;enum&quot;: &#91;&quot;split&quot;, &quot;join&quot;, &quot;none&quot;, null&#93;} |
| `reference_fasta` | string &#124; null | None | Optional reference FASTA path; enables indel left-alignment (-f) |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "required": [
        "vcf_path"
      ],
      "properties": {
        "operation": {
          "const": "normalize",
          "description": "Operation (fixed)"
        },
        "vcf_path": {
          "type": "string",
          "description": "Path to a local .vcf, .vcf.gz, or .bcf file"
        },
        "multiallelics": {
          "type": [
            "string",
            "null"
          ],
          "enum": [
            "split",
            "join",
            "none",
            null
          ],
          "description": "'split' (default) breaks multiallelic records into biallelic; 'join' merges; 'none' leaves as-is"
        },
        "reference_fasta": {
          "type": [
            "string",
            "null"
          ],
          "description": "Optional reference FASTA path; enables indel left-alignment (-f)"
        }
      }
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

**Execution and provenance:** `external_process`；upstream name `VCF_normalize`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("vcf_normalize", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.vcf_summary_stats` {#dnakit.databases.vcf_summary_stats}

Canonical variant statistics for a local VCF/BCF file via bcftools stats. Returns records, SNPs, indels, MNPs, multiallelic sites, ts/tv ratio, and per-sample het/hom/missing counts. Optional PASS-only, region, QUAL, and bcftools-expression filters. Use this instead of writing an ad-hoc parser so counts are reproducible. Requires bcftools on PATH.

```text
dnakit.databases.vcf_summary_stats(*, operation=None, vcf_path, pass_only=None, regions=None, min_qual=None, include_expr=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | object &#124; None | None | Operation (fixed) {&quot;const&quot;: &quot;summary_stats&quot;} |
| `vcf_path` | string | Required by input schema | Path to a local .vcf, .vcf.gz, or .bcf file |
| `pass_only` | boolean &#124; null | None | If true, count only FILTER=PASS or &#x27;.&#x27; records |
| `regions` | string &#124; null | None | Restrict to a region, e.g. &#x27;chr1&#x27; or &#x27;chr1:1-1000000&#x27; (streamed, no index needed) |
| `min_qual` | number &#124; null | None | Keep only records with QUAL &gt;= this value |
| `include_expr` | string &#124; null | None | Extra bcftools -i include expression, e.g. &#x27;INFO/DP&gt;10&#x27; |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "required": [
        "vcf_path"
      ],
      "properties": {
        "operation": {
          "const": "summary_stats",
          "description": "Operation (fixed)"
        },
        "vcf_path": {
          "type": "string",
          "description": "Path to a local .vcf, .vcf.gz, or .bcf file"
        },
        "pass_only": {
          "type": [
            "boolean",
            "null"
          ],
          "description": "If true, count only FILTER=PASS or '.' records"
        },
        "regions": {
          "type": [
            "string",
            "null"
          ],
          "description": "Restrict to a region, e.g. 'chr1' or 'chr1:1-1000000' (streamed, no index needed)"
        },
        "min_qual": {
          "type": [
            "number",
            "null"
          ],
          "description": "Keep only records with QUAL >= this value"
        },
        "include_expr": {
          "type": [
            "string",
            "null"
          ],
          "description": "Extra bcftools -i include expression, e.g. 'INFO/DP>10'"
        }
      }
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

**Execution and provenance:** `external_process`；upstream name `VCF_summary_stats`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("vcf_summary_stats", verbose=False)
print(info["input_schema"])
```

## `dnakit.workflows.simulate_inheritance` {#dnakit.workflows.simulate_inheritance}

Simulate independent allele inheritance using explicit synthetic or supplied data.

The attributed Genomebook engine implements per-locus segregation and mutation;
it does not simulate chromosome-scale crossover or predict real clinical outcomes.
A local random generator avoids changing randomness elsewhere in the process.

```text
dnakit.workflows.simulate_inheritance(father, mother, trait_registry, disease_registry, *, offspring_count=3, generation=1, mutation_rate=0.001, seed=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `father` | dict&#91;str, Any&#93; | Required | Father record with sex=Male and Genomebook genotype fields. |
| `mother` | dict&#91;str, Any&#93; | Required | Mother record with sex=Female and Genomebook genotype fields. |
| `trait_registry` | dict&#91;str, Any&#93; | Required | Explicit trait and locus registry. |
| `disease_registry` | dict&#91;str, Any&#93; | Required | Explicit disease locus registry. |
| `offspring_count` | int | 3 | Number of simulated offspring, 1–10000. |
| `generation` | int | 1 | Positive integer generation label. |
| `mutation_rate` | float | 0.001 | Per-locus simulated mutation rate, between 0 and 1. |
| `seed` | int &#124; None | None | Local random seed; does not modify global RNG state. |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Mapping with offspring, seed, mutation_rate, scope and upstream; independent-locus simulation, not chromosome crossover or clinical prediction.

