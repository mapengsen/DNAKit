# 遗传变异、GWAS 与变异评分

本页收录 **80 个新增公开函数**。 [完整索引](index.md) · [调用约定与完整示例](usage.md)

每个函数均列出当前 Python 签名、参数、返回值和已登记的执行要求。提供方说明及字段描述保留源码原文。`None` 默认值在生成的数据库包装中通常表示不向服务传该参数；实际必填项与组合约束以输入 Schema 为准。坐标系和基因组版本按各接口说明使用。

## 函数目录

- [`dnakit.databases.alpha_missense_get_protein_scores`](#dnakit.databases.alpha_missense_get_protein_scores)
- [`dnakit.databases.alpha_missense_get_residue_scores`](#dnakit.databases.alpha_missense_get_residue_scores)
- [`dnakit.databases.alpha_missense_get_variant_score`](#dnakit.databases.alpha_missense_get_variant_score)
- [`dnakit.databases.bio_bank_japan_phewas_by_variant`](#dnakit.databases.bio_bank_japan_phewas_by_variant)
- [`dnakit.databases.cadd_get_position_scores`](#dnakit.databases.cadd_get_position_scores)
- [`dnakit.databases.cadd_get_range_scores`](#dnakit.databases.cadd_get_range_scores)
- [`dnakit.databases.cadd_get_scores`](#dnakit.databases.cadd_get_scores)
- [`dnakit.databases.cadd_get_variant_score`](#dnakit.databases.cadd_get_variant_score)
- [`dnakit.databases.civic_get_molecular_profile`](#dnakit.databases.civic_get_molecular_profile)
- [`dnakit.databases.civic_get_variant`](#dnakit.databases.civic_get_variant)
- [`dnakit.databases.civic_get_variants_by_gene`](#dnakit.databases.civic_get_variants_by_gene)
- [`dnakit.databases.civic_search_assertions`](#dnakit.databases.civic_search_assertions)
- [`dnakit.databases.civic_search_evidence_items`](#dnakit.databases.civic_search_evidence_items)
- [`dnakit.databases.civic_search_genes`](#dnakit.databases.civic_search_genes)
- [`dnakit.databases.civic_search_therapies`](#dnakit.databases.civic_search_therapies)
- [`dnakit.databases.clin_gen_dosage_by_gene`](#dnakit.databases.clin_gen_dosage_by_gene)
- [`dnakit.databases.clin_gen_dosage_region_search`](#dnakit.databases.clin_gen_dosage_region_search)
- [`dnakit.databases.clin_gen_get_actionability_adult`](#dnakit.databases.clin_gen_get_actionability_adult)
- [`dnakit.databases.clin_gen_get_actionability_pediatric`](#dnakit.databases.clin_gen_get_actionability_pediatric)
- [`dnakit.databases.clin_gen_get_dosage_sensitivity`](#dnakit.databases.clin_gen_get_dosage_sensitivity)
- [`dnakit.databases.clin_gen_get_gene_validity`](#dnakit.databases.clin_gen_get_gene_validity)
- [`dnakit.databases.clin_gen_get_variant_classifications`](#dnakit.databases.clin_gen_get_variant_classifications)
- [`dnakit.databases.clin_gen_search_actionability`](#dnakit.databases.clin_gen_search_actionability)
- [`dnakit.databases.clin_gen_search_dosage_sensitivity`](#dnakit.databases.clin_gen_search_dosage_sensitivity)
- [`dnakit.databases.clin_gen_search_gene_validity`](#dnakit.databases.clin_gen_search_gene_validity)
- [`dnakit.databases.clin_var_get_clinical_significance`](#dnakit.databases.clin_var_get_clinical_significance)
- [`dnakit.databases.clin_var_get_variant_details`](#dnakit.databases.clin_var_get_variant_details)
- [`dnakit.databases.clin_var_search_variants`](#dnakit.databases.clin_var_search_variants)
- [`dnakit.databases.cosmic_get_mutations_by_gene`](#dnakit.databases.cosmic_get_mutations_by_gene)
- [`dnakit.databases.cosmic_search_mutations`](#dnakit.databases.cosmic_search_mutations)
- [`dnakit.databases.dbsnp_get_frequencies`](#dnakit.databases.dbsnp_get_frequencies)
- [`dnakit.databases.dbsnp_get_variant_by_rsid`](#dnakit.databases.dbsnp_get_variant_by_rsid)
- [`dnakit.databases.ensembl_var_get_variant_consequences`](#dnakit.databases.ensembl_var_get_variant_consequences)
- [`dnakit.databases.ensembl_vep_annotate_hgvs`](#dnakit.databases.ensembl_vep_annotate_hgvs)
- [`dnakit.databases.ensembl_vep_annotate_rsid`](#dnakit.databases.ensembl_vep_annotate_rsid)
- [`dnakit.databases.ensembl_vep_variant_recoder`](#dnakit.databases.ensembl_vep_variant_recoder)
- [`dnakit.databases.eve_get_gene_info`](#dnakit.databases.eve_get_gene_info)
- [`dnakit.databases.eve_get_variant_score`](#dnakit.databases.eve_get_variant_score)
- [`dnakit.databases.favor_annotate_variant`](#dnakit.databases.favor_annotate_variant)
- [`dnakit.databases.finn_gen_get_variant_finemapping`](#dnakit.databases.finn_gen_get_variant_finemapping)
- [`dnakit.databases.genebass_gene_burden_phewas`](#dnakit.databases.genebass_gene_burden_phewas)
- [`dnakit.databases.gnomad_get_gene_constraints`](#dnakit.databases.gnomad_get_gene_constraints)
- [`dnakit.databases.gnomad_get_sv_by_gene`](#dnakit.databases.gnomad_get_sv_by_gene)
- [`dnakit.databases.gnomad_get_sv_by_region`](#dnakit.databases.gnomad_get_sv_by_region)
- [`dnakit.databases.gnomad_get_sv_detail`](#dnakit.databases.gnomad_get_sv_detail)
- [`dnakit.databases.gnomad_get_variant`](#dnakit.databases.gnomad_get_variant)
- [`dnakit.databases.gnomad_get_variant_populations`](#dnakit.databases.gnomad_get_variant_populations)
- [`dnakit.databases.gnomad_search_variants`](#dnakit.databases.gnomad_search_variants)
- [`dnakit.databases.gwas_get_associations_for_snp`](#dnakit.databases.gwas_get_associations_for_snp)
- [`dnakit.databases.gwas_get_associations_for_study`](#dnakit.databases.gwas_get_associations_for_study)
- [`dnakit.databases.gwas_get_associations_for_trait`](#dnakit.databases.gwas_get_associations_for_trait)
- [`dnakit.databases.gwas_get_snp_by_id`](#dnakit.databases.gwas_get_snp_by_id)
- [`dnakit.databases.gwas_get_snps_for_gene`](#dnakit.databases.gwas_get_snps_for_gene)
- [`dnakit.databases.gwas_get_studies_for_trait`](#dnakit.databases.gwas_get_studies_for_trait)
- [`dnakit.databases.gwas_get_study_by_id`](#dnakit.databases.gwas_get_study_by_id)
- [`dnakit.databases.gwas_get_variants_for_trait`](#dnakit.databases.gwas_get_variants_for_trait)
- [`dnakit.databases.gwas_search_associations`](#dnakit.databases.gwas_search_associations)
- [`dnakit.databases.gwas_search_associations_by_gene`](#dnakit.databases.gwas_search_associations_by_gene)
- [`dnakit.databases.gwas_search_snps`](#dnakit.databases.gwas_search_snps)
- [`dnakit.databases.gwas_search_studies`](#dnakit.databases.gwas_search_studies)
- [`dnakit.databases.igsr_list_data_collections`](#dnakit.databases.igsr_list_data_collections)
- [`dnakit.databases.igsr_search_populations`](#dnakit.databases.igsr_search_populations)
- [`dnakit.databases.igsr_search_samples`](#dnakit.databases.igsr_search_samples)
- [`dnakit.databases.ldlink_get_proxies`](#dnakit.databases.ldlink_get_proxies)
- [`dnakit.databases.mave_db_get_effect_matrix`](#dnakit.databases.mave_db_get_effect_matrix)
- [`dnakit.databases.my_variant_get_pathogenicity_scores`](#dnakit.databases.my_variant_get_pathogenicity_scores)
- [`dnakit.databases.my_variant_get_variant_annotation`](#dnakit.databases.my_variant_get_variant_annotation)
- [`dnakit.databases.my_variant_query_variants`](#dnakit.databases.my_variant_query_variants)
- [`dnakit.databases.onco_kb_annotate_variant`](#dnakit.databases.onco_kb_annotate_variant)
- [`dnakit.databases.progenetix_cnv_search`](#dnakit.databases.progenetix_cnv_search)
- [`dnakit.databases.progenetix_list_cohorts`](#dnakit.databases.progenetix_list_cohorts)
- [`dnakit.databases.progenetix_list_filtering_terms`](#dnakit.databases.progenetix_list_filtering_terms)
- [`dnakit.databases.progenetix_search_biosamples`](#dnakit.databases.progenetix_search_biosamples)
- [`dnakit.databases.splice_ai_get_max_delta`](#dnakit.databases.splice_ai_get_max_delta)
- [`dnakit.databases.splice_ai_predict_pangolin`](#dnakit.databases.splice_ai_predict_pangolin)
- [`dnakit.databases.splice_ai_predict_splice`](#dnakit.databases.splice_ai_predict_splice)
- [`dnakit.databases.tpmi_phewas_by_variant`](#dnakit.databases.tpmi_phewas_by_variant)
- [`dnakit.databases.ukbtopmed_phewas_by_variant`](#dnakit.databases.ukbtopmed_phewas_by_variant)
- [`dnakit.databases.variant_validator_gene2transcripts`](#dnakit.databases.variant_validator_gene2transcripts)
- [`dnakit.databases.variant_validator_validate_variant`](#dnakit.databases.variant_validator_validate_variant)

## `dnakit.databases.alpha_missense_get_protein_scores` {#dnakit.databases.alpha_missense_get_protein_scores}

Get AlphaMissense pathogenicity scores for EVERY residue of a protein in one call. The upstream hegelab API only exposes per-residue queries; this tool hides that by looking up the protein length from UniProt and concurrently fetching all residues via a thread pool. Returns a list of per-position records (one per residue) plus protein length, count of successful fetches, and the pathogenic/ambiguous/benign thresholds. Use this when you want whole-protein vulnerability assessment, hotspot detection, or to build a (20 × n_positions) AlphaMissense score matrix for benchmarking. Set max_residues > 0 to cap the loop for fast diagnostic calls.

```text
dnakit.databases.alpha_missense_get_protein_scores(*, uniprot_id, max_residues=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `uniprot_id` | string | 必填（输入 Schema） | UniProt accession ID (e.g., &#x27;P00533&#x27; for EGFR, &#x27;P01116&#x27; for KRAS, &#x27;P04637&#x27; for TP53). |
| `max_residues` | integer &#124; null | None | Optional cap on the number of residues fetched. 0 / null / omitted = no cap (fetch every residue of the protein, recommended for whole-protein analysis). Positive integer = fetch only the first N residues (useful for fast diagnostic calls on very long proteins). {&quot;default&quot;: 0} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "uniprot_id": {
          "type": "string",
          "description": "UniProt accession ID (e.g., 'P00533' for EGFR, 'P01116' for KRAS, 'P04637' for TP53)."
        },
        "max_residues": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Optional cap on the number of residues fetched. 0 / null / omitted = no cap (fetch every residue of the protein, recommended for whole-protein analysis). Positive integer = fetch only the first N residues (useful for fast diagnostic calls on very long proteins).",
          "default": 0
        }
      },
      "required": [
        "uniprot_id"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "uniprot_id": {
              "type": "string"
            },
            "protein_length": {
              "type": "integer"
            },
            "scores": {
              "type": "array",
              "description": "Per-residue records, one per position with successful fetch. Each record carries position, uid, aa, resi, benign/ambiguous/pathogenic bin strings (e.g. '6:A,C,D,R,S,V'), benign_all/ambiguous_all/pathogenic_all (covering all 19 substitutions, not just SNV-reachable), and mean/mean_all numeric scores.",
              "items": {
                "type": "object"
              }
            },
            "n_positions_returned": {
              "type": "integer"
            },
            "n_positions_attempted": {
              "type": "integer"
            },
            "max_residues_cap": {
              "type": [
                "integer",
                "null"
              ]
            },
            "thresholds": {
              "type": "object"
            },
            "pdb_download": {
              "type": "string"
            }
          },
          "required": [
            "uniprot_id",
            "protein_length",
            "scores",
            "n_positions_returned"
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

**执行与来源：** `network`；上游接口名 `AlphaMissense_get_protein_scores`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("alpha_missense_get_protein_scores", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.alpha_missense_get_residue_scores` {#dnakit.databases.alpha_missense_get_residue_scores}

Get the AlphaMissense class breakdown for substitutions at a specific protein position: which residues fall in each class (benign/ambiguous/pathogenic), reported separately for the SNV-accessible substitutions and for all 19, plus the residue mean scores. NOTE: the upstream hotspot API does not expose a numeric score per individual substitution, so this returns class membership rather than 20 scores.

```text
dnakit.databases.alpha_missense_get_residue_scores(*, uniprot_id, position, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `uniprot_id` | string | 必填（输入 Schema） | UniProt accession ID (e.g., &#x27;P00533&#x27; for EGFR) |
| `position` | integer | 必填（输入 Schema） | Amino acid position in the protein (1-indexed) |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "uniprot_id": {
          "type": "string",
          "description": "UniProt accession ID (e.g., 'P00533' for EGFR)"
        },
        "position": {
          "type": "integer",
          "description": "Amino acid position in the protein (1-indexed)"
        }
      },
      "required": [
        "position",
        "uniprot_id"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

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
            "uniprot_id": {
              "type": "string"
            },
            "position": {
              "type": "integer"
            },
            "scores": {
              "type": "object"
            },
            "thresholds": {
              "type": "object"
            }
          }
        },
        "error": {
          "type": "string"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `AlphaMissense_get_residue_scores`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("alpha_missense_get_residue_scores", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.alpha_missense_get_variant_score` {#dnakit.databases.alpha_missense_get_variant_score}

Classify a specific missense variant with AlphaMissense. Input: UniProt ID and variant (e.g., 'p.R123H'). Returns the variant's AlphaMissense class (pathogenic/ambiguous/benign) and whether it came from the SNV-accessible or the all-substitutions set. NOTE: the upstream hotspot API publishes per-residue class membership, not a numeric score per substitution, so `pathogenicity_score` is null and `score_available` is false; the `residue_mean_score_*` fields are residue-level averages, NOT this variant's score. For a numeric per-variant score use GenomeNexus_annotate_variant.

```text
dnakit.databases.alpha_missense_get_variant_score(*, uniprot_id, variant, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `uniprot_id` | string | 必填（输入 Schema） | UniProt accession ID (e.g., &#x27;P00533&#x27; for EGFR) |
| `variant` | string | 必填（输入 Schema） | Variant in protein notation: p.X123Y or X123Y where X is reference amino acid, 123 is position, Y is variant (e.g., &#x27;p.R123H&#x27;, &#x27;V600E&#x27;) |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "uniprot_id": {
          "type": "string",
          "description": "UniProt accession ID (e.g., 'P00533' for EGFR)"
        },
        "variant": {
          "type": "string",
          "description": "Variant in protein notation: p.X123Y or X123Y where X is reference amino acid, 123 is position, Y is variant (e.g., 'p.R123H', 'V600E')"
        }
      },
      "required": [
        "uniprot_id",
        "variant"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

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
            "uniprot_id": {
              "type": "string"
            },
            "variant": {
              "type": "string"
            },
            "position": {
              "type": "integer"
            },
            "reference_aa": {
              "type": "string"
            },
            "variant_aa": {
              "type": "string"
            },
            "pathogenicity_score": {
              "type": [
                "number",
                "null"
              ],
              "description": "Always null: the upstream hotspot API exposes no per-substitution score."
            },
            "score_available": {
              "type": "boolean"
            },
            "classification": {
              "type": [
                "string",
                "null"
              ],
              "enum": [
                "pathogenic",
                "ambiguous",
                "benign",
                null
              ]
            },
            "substitution_set": {
              "type": [
                "string",
                "null"
              ],
              "enum": [
                "snv_accessible",
                "all_substitutions",
                null
              ]
            },
            "residue_mean_score_snv_accessible": {
              "type": [
                "number",
                "null"
              ]
            },
            "residue_mean_score_all_substitutions": {
              "type": [
                "number",
                "null"
              ]
            },
            "residue_class_counts": {
              "type": "object"
            },
            "note": {
              "type": "string"
            }
          }
        },
        "error": {
          "type": "string"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `AlphaMissense_get_variant_score`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("alpha_missense_get_variant_score", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.bio_bank_japan_phewas_by_variant` {#dnakit.databases.bio_bank_japan_phewas_by_variant}

Phenome-wide association study (PheWAS) lookup for a single variant in BioBank Japan (BBJ), the largest Japanese-ancestry biobank. Returns every phenotype associated with the variant (the inverse of a GWAS), sorted by p-value. Accepts an rsID (auto-resolved to GRCh37 via Ensembl) or an explicit chr:pos:ref:alt variant in GRCh37. Use for cross-ancestry replication of a genetic association in an East Asian (Japanese) population.

```text
dnakit.databases.bio_bank_japan_phewas_by_variant(*, rsid=None, variant=None, limit=None, max_pval=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `rsid` | string | None | dbSNP rsID, e.g. &#x27;rs7903146&#x27;. Auto-resolved to GRCh37 coordinates via Ensembl. Provide this OR &#x27;variant&#x27;. |
| `variant` | string | None | Variant in chr:pos:ref:alt or chr-pos-ref-alt format, in GRCh37, e.g. &#x27;10:114758349:C:T&#x27;. Provide this OR &#x27;rsid&#x27;. A bare rsID is also accepted here. |
| `limit` | integer | None | Maximum number of phenotype associations to return (sorted by ascending p-value). {&quot;default&quot;: 20} |
| `max_pval` | number &#124; null | None | Optional: only return associations with p-value &lt;= this threshold (e.g. 5e-8 for genome-wide significance). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "rsid": {
          "type": "string",
          "description": "dbSNP rsID, e.g. 'rs7903146'. Auto-resolved to GRCh37 coordinates via Ensembl. Provide this OR 'variant'."
        },
        "variant": {
          "type": "string",
          "description": "Variant in chr:pos:ref:alt or chr-pos-ref-alt format, in GRCh37, e.g. '10:114758349:C:T'. Provide this OR 'rsid'. A bare rsID is also accepted here."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of phenotype associations to return (sorted by ascending p-value).",
          "default": 20
        },
        "max_pval": {
          "type": [
            "number",
            "null"
          ],
          "description": "Optional: only return associations with p-value <= this threshold (e.g. 5e-8 for genome-wide significance)."
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "variant": {
              "type": "string"
            },
            "rsids": {
              "type": [
                "string",
                "null"
              ]
            },
            "nearest_genes": {
              "type": [
                "string",
                "null"
              ]
            },
            "associations": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "phenocode": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "phenostring": {
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
                  "pval": {
                    "type": [
                      "number",
                      "null"
                    ]
                  },
                  "beta": {
                    "type": [
                      "number",
                      "null"
                    ]
                  },
                  "sebeta": {
                    "type": [
                      "number",
                      "null"
                    ]
                  },
                  "af": {
                    "type": [
                      "number",
                      "null"
                    ]
                  },
                  "num_cases": {
                    "type": [
                      "number",
                      "null"
                    ]
                  },
                  "num_controls": {
                    "type": [
                      "number",
                      "null"
                    ]
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

**执行与来源：** `network`；上游接口名 `BioBankJapan_phewas_by_variant`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("bio_bank_japan_phewas_by_variant", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.cadd_get_position_scores` {#dnakit.databases.cadd_get_position_scores}

Get CADD scores for all possible substitutions at a genomic position. Returns PHRED scores for A, C, G, T variants at the position. Useful for saturation mutagenesis analysis.

```text
dnakit.databases.cadd_get_position_scores(*, chrom, pos, version=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `chrom` | string | 必填（输入 Schema） | Chromosome (1-22, X, Y, MT) |
| `pos` | integer | 必填（输入 Schema） | Genomic position (1-based) |
| `version` | string | None | CADD version and genome build {&quot;default&quot;: &quot;GRCh38-v1.7&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "chrom": {
          "type": "string",
          "description": "Chromosome (1-22, X, Y, MT)"
        },
        "pos": {
          "type": "integer",
          "description": "Genomic position (1-based)"
        },
        "version": {
          "type": "string",
          "default": "GRCh38-v1.7",
          "description": "CADD version and genome build"
        }
      },
      "required": [
        "chrom",
        "pos"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "chrom": {
          "type": "string"
        },
        "pos": {
          "type": "integer"
        },
        "variants": {
          "type": "array"
        },
        "version": {
          "type": "string"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `CADD_get_position_scores`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("cadd_get_position_scores", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.cadd_get_range_scores` {#dnakit.databases.cadd_get_range_scores}

Get CADD scores for all variants in a genomic range (max 100bp). Returns all pre-computed deleteriousness scores in the region. Use for hotspot analysis or exon scanning.

```text
dnakit.databases.cadd_get_range_scores(*, chrom, start, end, version=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `chrom` | string | 必填（输入 Schema） | Chromosome (1-22, X, Y, MT) |
| `start` | integer | 必填（输入 Schema） | Start position (1-based, inclusive) |
| `end` | integer | 必填（输入 Schema） | End position (1-based, inclusive). Max range: 100bp. |
| `version` | string | None | CADD version and genome build {&quot;default&quot;: &quot;GRCh38-v1.7&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "chrom": {
          "type": "string",
          "description": "Chromosome (1-22, X, Y, MT)"
        },
        "start": {
          "type": "integer",
          "description": "Start position (1-based, inclusive)"
        },
        "end": {
          "type": "integer",
          "description": "End position (1-based, inclusive). Max range: 100bp."
        },
        "version": {
          "type": "string",
          "default": "GRCh38-v1.7",
          "description": "CADD version and genome build"
        }
      },
      "required": [
        "chrom",
        "end",
        "start"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "chrom": {
          "type": "string"
        },
        "start": {
          "type": "integer"
        },
        "end": {
          "type": "integer"
        },
        "variants": {
          "type": "array"
        },
        "count": {
          "type": "integer"
        },
        "version": {
          "type": "string"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `CADD_get_range_scores`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("cadd_get_range_scores", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.cadd_get_scores` {#dnakit.databases.cadd_get_scores}

Query indexed CADD annotations through MyVariant; preserve missing scores as missing.

```text
dnakit.databases.cadd_get_scores(*, variant, assembly='hg38', verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant` | string | 必填（输入 Schema） | 用于 MyVariant CADD 注释检索的变异查询条件。 |
| `assembly` | string | &#x27;hg38&#x27; | 参考基因组版本，例如 hg38。 {&quot;default&quot;: &quot;hg38&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回 MyVariant 检索结果；字典结果另含 dnakit_score_source，表明分数来自已索引的 CADD 注释。缺失分数保持缺失，不执行新的模型预测。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "variant": {
          "type": "string"
        },
        "assembly": {
          "type": "string",
          "default": "hg38"
        }
      },
      "required": [
        "variant"
      ]
    }
    ```

**执行与来源：** `network`；上游接口名 `CADD_get_scores`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("cadd_get_scores", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.cadd_get_variant_score` {#dnakit.databases.cadd_get_variant_score}

Get CADD deleteriousness score for a specific variant. PHRED scores: >=20 top 1% deleterious (pathogenic threshold ~15-20), >=30 top 0.1%. Input: chromosome, position, ref, alt alleles. Supports GRCh37 and GRCh38.

```text
dnakit.databases.cadd_get_variant_score(*, chrom, pos, ref, alt, version=None, include_annotations=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `chrom` | string | 必填（输入 Schema） | Chromosome (1-22, X, Y, MT). Can include &#x27;chr&#x27; prefix. |
| `pos` | integer | 必填（输入 Schema） | Genomic position (1-based) |
| `ref` | string | 必填（输入 Schema） | Reference allele (e.g., &#x27;A&#x27;, &#x27;G&#x27;) |
| `alt` | string | 必填（输入 Schema） | Alternate allele (e.g., &#x27;T&#x27;, &#x27;C&#x27;) |
| `version` | string | None | CADD version and genome build {&quot;default&quot;: &quot;GRCh38-v1.7&quot;, &quot;enum&quot;: &#91;&quot;GRCh38-v1.7&quot;, &quot;GRCh37-v1.7&quot;, &quot;GRCh38-v1.6&quot;, &quot;GRCh37-v1.6&quot;&#93;} |
| `include_annotations` | boolean | None | Include full annotation details in response {&quot;default&quot;: false} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "chrom": {
          "type": "string",
          "description": "Chromosome (1-22, X, Y, MT). Can include 'chr' prefix."
        },
        "pos": {
          "type": "integer",
          "description": "Genomic position (1-based)"
        },
        "ref": {
          "type": "string",
          "description": "Reference allele (e.g., 'A', 'G')"
        },
        "alt": {
          "type": "string",
          "description": "Alternate allele (e.g., 'T', 'C')"
        },
        "version": {
          "type": "string",
          "default": "GRCh38-v1.7",
          "enum": [
            "GRCh38-v1.7",
            "GRCh37-v1.7",
            "GRCh38-v1.6",
            "GRCh37-v1.6"
          ],
          "description": "CADD version and genome build"
        },
        "include_annotations": {
          "type": "boolean",
          "default": false,
          "description": "Include full annotation details in response"
        }
      },
      "required": [
        "alt",
        "chrom",
        "pos",
        "ref"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": [
        "object",
        "null"
      ],
      "properties": {
        "chrom": {
          "type": "string"
        },
        "pos": {
          "type": "integer"
        },
        "ref": {
          "type": "string"
        },
        "alt": {
          "type": "string"
        },
        "phred_score": {
          "type": [
            "number",
            "null"
          ]
        },
        "raw_score": {
          "type": [
            "number",
            "null"
          ]
        },
        "interpretation": {
          "type": [
            "string",
            "null"
          ]
        },
        "version": {
          "type": "string"
        },
        "thresholds": {
          "type": "object"
        },
        "annotations": {
          "type": [
            "object",
            "null"
          ]
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `CADD_get_variant_score`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("cadd_get_variant_score", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.civic_get_molecular_profile` {#dnakit.databases.civic_get_molecular_profile}

Get detailed information about a specific molecular profile in CIViC database by molecular profile ID. Molecular profiles represent combinations of variants or features (e.g., 'BRAF V600E') that serve as biomarkers for clinical interpretation. Returns profile name, description, clinical significance score, and associated evidence items with diseases, therapies, and evidence levels.

```text
dnakit.databases.civic_get_molecular_profile(*, molecular_profile_id, evidence_limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `molecular_profile_id` | integer | 必填（输入 Schema） | CIViC molecular profile ID (e.g., 12 for BRAF V600E) |
| `evidence_limit` | integer | None | Maximum number of evidence items to return (default: 10) {&quot;default&quot;: 10} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "molecular_profile_id": {
          "type": "integer",
          "description": "CIViC molecular profile ID (e.g., 12 for BRAF V600E)"
        },
        "evidence_limit": {
          "type": "integer",
          "description": "Maximum number of evidence items to return (default: 10)",
          "default": 10
        }
      },
      "required": [
        "molecular_profile_id"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "molecularProfile": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "integer"
                },
                "name": {
                  "type": "string"
                },
                "molecularProfileScore": {
                  "type": "number",
                  "description": "CIViC clinical actionability score"
                },
                "description": {
                  "type": "string",
                  "description": "Clinical description of the molecular profile"
                },
                "evidenceItems": {
                  "type": "object",
                  "properties": {
                    "totalCount": {
                      "type": "integer",
                      "description": "Total number of records matching the query in CIViC, independent of how many were returned in this page. Use this as the denominator; len(nodes) is only the page size."
                    },
                    "pageInfo": {
                      "type": "object",
                      "description": "Pagination state for this connection.",
                      "properties": {
                        "hasNextPage": {
                          "type": "boolean",
                          "description": "True when more records match than were returned; the returned nodes are an incomplete sample of totalCount."
                        }
                      }
                    },
                    "nodes": {
                      "type": "array",
                      "items": {
                        "type": "object",
                        "properties": {
                          "id": {
                            "type": "integer"
                          },
                          "name": {
                            "type": "string"
                          },
                          "evidenceLevel": {
                            "type": "string",
                            "description": "A-E evidence level"
                          },
                          "evidenceDirection": {
                            "type": "string"
                          },
                          "significance": {
                            "type": "string"
                          },
                          "status": {
                            "type": "string"
                          },
                          "disease": {
                            "type": "object",
                            "properties": {
                              "name": {
                                "type": "string"
                              }
                            }
                          },
                          "therapies": {
                            "type": "array",
                            "items": {
                              "type": "object",
                              "properties": {
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
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `civic_get_molecular_profile`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("civic_get_molecular_profile", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.civic_get_variant` {#dnakit.databases.civic_get_variant}

Get detailed information about a specific variant in CIViC database by variant ID. Variants represent specific genetic alterations (SNVs, indels, structural variants, etc.) with clinical significance in cancer. Returns variant name, coordinates, gene association, and related information.

```text
dnakit.databases.civic_get_variant(*, variant_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant_id` | integer | 必填（输入 Schema） | CIViC variant ID (e.g., 4170) |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "variant_id": {
          "type": "integer",
          "description": "CIViC variant ID (e.g., 4170)"
        }
      },
      "required": [
        "variant_id"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "variant": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "integer",
                  "description": "CIViC variant ID"
                },
                "name": {
                  "type": "string",
                  "description": "Variant name"
                },
                "molecularProfiles": {
                  "type": "object",
                  "description": "Molecular profiles this variant participates in.",
                  "properties": {
                    "totalCount": {
                      "type": "integer",
                      "description": "Total number of molecular profiles for this variant in CIViC, independent of how many were returned in this page. Use this as the denominator; len(nodes) is only the page size."
                    },
                    "pageInfo": {
                      "type": "object",
                      "description": "Pagination state for this connection.",
                      "properties": {
                        "hasNextPage": {
                          "type": "boolean",
                          "description": "True when more molecular profiles match than were returned; the returned nodes are an incomplete sample of totalCount."
                        }
                      }
                    },
                    "nodes": {
                      "type": "array",
                      "items": {
                        "type": "object",
                        "properties": {
                          "id": {
                            "type": "integer"
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

**执行与来源：** `network`；上游接口名 `civic_get_variant`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("civic_get_variant", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.civic_get_variants_by_gene` {#dnakit.databases.civic_get_variants_by_gene}

Get all variants associated with a specific gene in CIViC database. Returns variant information including names, coordinates, and associated evidence. Use gene_name or gene_symbol (e.g., 'EGFR', 'BRAF') for convenience, or gene_id if already known. Uses cursor-based pagination to return up to 500 variants per call.

```text
dnakit.databases.civic_get_variants_by_gene(*, gene_id=None, gene_name=None, gene_symbol=None, gene=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | integer | None | CIViC gene ID (e.g., 19 for EGFR, 12 for BRAF). Find gene IDs using civic_search_genes. |
| `gene_name` | string | None | Gene symbol (e.g., &#x27;EGFR&#x27;, &#x27;BRAF&#x27;, &#x27;TP53&#x27;). Will be looked up automatically. Aliases: gene, gene_symbol, query. |
| `gene_symbol` | string | None | Alias for gene_name. Standard gene symbol (e.g., &#x27;KRAS&#x27;, &#x27;BRCA1&#x27;, &#x27;EGFR&#x27;). |
| `gene` | string | None | Alias for gene_name. Gene symbol (e.g., &#x27;EGFR&#x27;, &#x27;KRAS&#x27;). |
| `limit` | integer | None | Maximum number of variants to return (default: 500, uses cursor pagination to bypass CIViC&#x27;s 100/page server cap) {&quot;default&quot;: 500} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene_id": {
          "type": "integer",
          "description": "CIViC gene ID (e.g., 19 for EGFR, 12 for BRAF). Find gene IDs using civic_search_genes."
        },
        "gene_name": {
          "type": "string",
          "description": "Gene symbol (e.g., 'EGFR', 'BRAF', 'TP53'). Will be looked up automatically. Aliases: gene, gene_symbol, query."
        },
        "gene_symbol": {
          "type": "string",
          "description": "Alias for gene_name. Standard gene symbol (e.g., 'KRAS', 'BRCA1', 'EGFR')."
        },
        "gene": {
          "type": "string",
          "description": "Alias for gene_name. Gene symbol (e.g., 'EGFR', 'KRAS')."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of variants to return (default: 500, uses cursor pagination to bypass CIViC's 100/page server cap)",
          "default": 500
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "gene": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "integer"
                },
                "name": {
                  "type": "string"
                },
                "variants": {
                  "type": "object",
                  "properties": {
                    "totalCount": {
                      "type": "integer",
                      "description": "Total number of records matching the query in CIViC, independent of how many were returned in this page. Use this as the denominator; len(nodes) is only the page size."
                    },
                    "pageInfo": {
                      "type": "object",
                      "description": "Pagination state for this connection.",
                      "properties": {
                        "hasNextPage": {
                          "type": "boolean",
                          "description": "True when more records match than were returned; the returned nodes are an incomplete sample of totalCount."
                        }
                      }
                    },
                    "nodes": {
                      "type": "array",
                      "items": {
                        "type": "object",
                        "properties": {
                          "id": {
                            "type": "integer",
                            "description": "CIViC variant ID"
                          },
                          "name": {
                            "type": "string",
                            "description": "Variant name (e.g., 'G12C', 'L858R', 'V600E')"
                          },
                          "feature": {
                            "type": "object",
                            "description": "Parent gene feature",
                            "properties": {
                              "id": {
                                "type": "integer"
                              },
                              "name": {
                                "type": "string",
                                "description": "Gene name"
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
            "source": {
              "type": "string"
            },
            "format": {
              "type": "string"
            },
            "note": {
              "type": "string",
              "description": "Optional note about duplicate variant names or other data quality observations"
            }
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `civic_get_variants_by_gene`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("civic_get_variants_by_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.civic_search_assertions` {#dnakit.databases.civic_search_assertions}

Search for assertions in CIViC database. Assertions are higher-level clinical interpretations that integrate multiple evidence items into formal statements about clinical actionability. Returns a list of assertions with descriptions and associated molecular profiles or variants. Filter by therapy (e.g., 'ponatinib'), disease (e.g., 'leukemia'), or variant name.

```text
dnakit.databases.civic_search_assertions(*, limit=None, therapy=None, therapy_name=None, disease=None, disease_name=None, variant_name=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `limit` | integer | None | Maximum number of assertions to return (default: 20, recommended max: 100) {&quot;default&quot;: 20} |
| `therapy` | string | None | Filter by therapy/drug name (e.g., &#x27;imatinib&#x27;, &#x27;ponatinib&#x27;). Alias: therapy_name. |
| `therapy_name` | string | None | Alias for therapy. Filter by therapy/drug name. |
| `disease` | string | None | Filter by disease name (e.g., &#x27;leukemia&#x27;, &#x27;melanoma&#x27;). Alias: disease_name. |
| `disease_name` | string | None | Alias for disease. Filter by disease name. |
| `variant_name` | string | None | Filter by variant name (e.g., &#x27;V600E&#x27;, &#x27;T315I&#x27;). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "limit": {
          "type": "integer",
          "description": "Maximum number of assertions to return (default: 20, recommended max: 100)",
          "default": 20
        },
        "therapy": {
          "type": "string",
          "description": "Filter by therapy/drug name (e.g., 'imatinib', 'ponatinib'). Alias: therapy_name."
        },
        "therapy_name": {
          "type": "string",
          "description": "Alias for therapy. Filter by therapy/drug name."
        },
        "disease": {
          "type": "string",
          "description": "Filter by disease name (e.g., 'leukemia', 'melanoma'). Alias: disease_name."
        },
        "disease_name": {
          "type": "string",
          "description": "Alias for disease. Filter by disease name."
        },
        "variant_name": {
          "type": "string",
          "description": "Filter by variant name (e.g., 'V600E', 'T315I')."
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "assertions": {
              "type": "object",
              "properties": {
                "totalCount": {
                  "type": "integer",
                  "description": "Total number of records matching the query in CIViC, independent of how many were returned in this page. Use this as the denominator; len(nodes) is only the page size."
                },
                "pageInfo": {
                  "type": "object",
                  "description": "Pagination state for this connection.",
                  "properties": {
                    "hasNextPage": {
                      "type": "boolean",
                      "description": "True when more records match than were returned; the returned nodes are an incomplete sample of totalCount."
                    }
                  }
                },
                "nodes": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "id": {
                        "type": "integer"
                      },
                      "description": {
                        "type": "string"
                      },
                      "status": {
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

**执行与来源：** `network`；上游接口名 `civic_search_assertions`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("civic_search_assertions", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.civic_search_evidence_items` {#dnakit.databases.civic_search_evidence_items}

Search for evidence items in CIViC database. Evidence items are curated statements linking variants or molecular profiles to clinical outcomes, therapies, or disease states. Returns evidence with disease, therapies, significance, direction, and citation. Filter by therapy name (e.g., 'imatinib'), disease name (e.g., 'leukemia'), or molecular profile name (e.g., 'BRAF V600E').

```text
dnakit.databases.civic_search_evidence_items(*, limit=None, status=None, therapy=None, therapy_name=None, disease=None, disease_name=None, molecular_profile=None, evidence_type=None, significance=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `limit` | integer | None | Maximum number of evidence items to return (default: 20, recommended max: 100) {&quot;default&quot;: 20} |
| `status` | string &#124; null | None | Filter by curation status. Default: ACCEPTED (peer-reviewed). Options: ACCEPTED, SUBMITTED, REJECTED, ALL (returns all statuses combined). |
| `therapy` | string | None | Filter by therapy/drug name (e.g., &#x27;imatinib&#x27;, &#x27;pembrolizumab&#x27;). Alias: therapy_name. Note: matches any evidence item where the therapy appears, including combination regimens — results may include multi-drug combinations. |
| `therapy_name` | string | None | Alias for therapy. Filter by therapy/drug name. |
| `disease` | string | None | Filter by disease name (e.g., &#x27;leukemia&#x27;, &#x27;melanoma&#x27;, &#x27;lung cancer&#x27;). Alias: disease_name. Note: CIViC uses specific disease names (e.g., &#x27;Lung Non-small Cell Carcinoma&#x27;, not &#x27;NSCLC&#x27;); try partial names or multiple searches if results are empty. |
| `disease_name` | string | None | Alias for disease. Filter by disease name. |
| `molecular_profile` | string | None | Filter by molecular profile name (e.g., &#x27;BRAF V600E&#x27;, &#x27;EGFR T790M&#x27;, &#x27;KRAS G12C&#x27;). Uses substring matching — &#x27;FLT3 ITD&#x27; will also match &#x27;FLT3 ITD AND FLT3 D835Y&#x27;. For gene fusions, CIViC uses double-colon notation: &#x27;GENE::PARTNER Fusion&#x27; (e.g., &#x27;FGFR2::BICC1 Fusion&#x27;, &#x27;ALK::EML4 Fusion&#x27;). Use civic_search_molecular_profiles to discover exact profile names. Gene-prefixed clinical spellings (&#x27;EGFRvIII&#x27;, &#x27;BRAFV600E&#x27;) are also handled: CIViC names the profile &#x27;&lt;gene&gt; &lt;designation&gt;&#x27; (&#x27;EGFR VIII&#x27;), so that separated form is queried too and both result sets are merged, reported in normalization_note. |
| `evidence_type` | string &#124; null | None | Filter by evidence type. Values: PREDICTIVE (drug response), DIAGNOSTIC (disease diagnosis), PROGNOSTIC (patient outcomes), PREDISPOSING (disease risk), ONCOGENIC (variant pathogenicity), FUNCTIONAL (molecular function). |
| `significance` | string &#124; null | None | Filter by clinical significance. Values depend on evidence_type, e.g. SENSITIVITYRESPONSE, RESISTANCE (predictive); BETTER_OUTCOME, POOR_OUTCOME (prognostic); POSITIVE, NEGATIVE (diagnostic); ONCOGENIC, PROTECTIVE. Use to separate, e.g., resistance from sensitivity evidence for a variant. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "limit": {
          "type": "integer",
          "description": "Maximum number of evidence items to return (default: 20, recommended max: 100)",
          "default": 20
        },
        "status": {
          "type": [
            "string",
            "null"
          ],
          "description": "Filter by curation status. Default: ACCEPTED (peer-reviewed). Options: ACCEPTED, SUBMITTED, REJECTED, ALL (returns all statuses combined)."
        },
        "therapy": {
          "type": "string",
          "description": "Filter by therapy/drug name (e.g., 'imatinib', 'pembrolizumab'). Alias: therapy_name. Note: matches any evidence item where the therapy appears, including combination regimens — results may include multi-drug combinations."
        },
        "therapy_name": {
          "type": "string",
          "description": "Alias for therapy. Filter by therapy/drug name."
        },
        "disease": {
          "type": "string",
          "description": "Filter by disease name (e.g., 'leukemia', 'melanoma', 'lung cancer'). Alias: disease_name. Note: CIViC uses specific disease names (e.g., 'Lung Non-small Cell Carcinoma', not 'NSCLC'); try partial names or multiple searches if results are empty."
        },
        "disease_name": {
          "type": "string",
          "description": "Alias for disease. Filter by disease name."
        },
        "molecular_profile": {
          "type": "string",
          "description": "Filter by molecular profile name (e.g., 'BRAF V600E', 'EGFR T790M', 'KRAS G12C'). Uses substring matching — 'FLT3 ITD' will also match 'FLT3 ITD AND FLT3 D835Y'. For gene fusions, CIViC uses double-colon notation: 'GENE::PARTNER Fusion' (e.g., 'FGFR2::BICC1 Fusion', 'ALK::EML4 Fusion'). Use civic_search_molecular_profiles to discover exact profile names. Gene-prefixed clinical spellings ('EGFRvIII', 'BRAFV600E') are also handled: CIViC names the profile '<gene> <designation>' ('EGFR VIII'), so that separated form is queried too and both result sets are merged, reported in normalization_note."
        },
        "evidence_type": {
          "type": [
            "string",
            "null"
          ],
          "description": "Filter by evidence type. Values: PREDICTIVE (drug response), DIAGNOSTIC (disease diagnosis), PROGNOSTIC (patient outcomes), PREDISPOSING (disease risk), ONCOGENIC (variant pathogenicity), FUNCTIONAL (molecular function)."
        },
        "significance": {
          "type": [
            "string",
            "null"
          ],
          "description": "Filter by clinical significance. Values depend on evidence_type, e.g. SENSITIVITYRESPONSE, RESISTANCE (predictive); BETTER_OUTCOME, POOR_OUTCOME (prognostic); POSITIVE, NEGATIVE (diagnostic); ONCOGENIC, PROTECTIVE. Use to separate, e.g., resistance from sensitivity evidence for a variant."
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "evidenceItems": {
              "type": "object",
              "properties": {
                "totalCount": {
                  "type": "integer",
                  "description": "Total number of records matching the query in CIViC, independent of how many were returned in this page. Use this as the denominator; len(nodes) is only the page size."
                },
                "pageInfo": {
                  "type": "object",
                  "description": "Pagination state for this connection.",
                  "properties": {
                    "hasNextPage": {
                      "type": "boolean",
                      "description": "True when more records match than were returned; the returned nodes are an incomplete sample of totalCount."
                    }
                  }
                },
                "nodes": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "id": {
                        "type": "integer"
                      },
                      "description": {
                        "type": "string"
                      },
                      "evidenceLevel": {
                        "type": "string"
                      },
                      "evidenceType": {
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

**执行与来源：** `network`；上游接口名 `civic_search_evidence_items`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("civic_search_evidence_items", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.civic_search_genes` {#dnakit.databases.civic_search_genes}

Search for genes in CIViC (Clinical Interpretation of Variants in Cancer) database. CIViC is a community knowledgebase for expert-curated interpretations of variants in cancer. Returns genes with their IDs, names, and descriptions. Use this to find genes of interest before querying variants or evidence.

```text
dnakit.databases.civic_search_genes(*, limit=None, name=None, query=None, gene_name=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `limit` | integer | None | Maximum number of genes to return (default: 10, recommended max: 100) {&quot;default&quot;: 10} |
| `name` | string | None | Gene symbol to search for (e.g., &quot;EGFR&quot;, &quot;BRAF&quot;, &quot;BRCA1&quot;). Alias: use &#x27;query&#x27; or &#x27;gene_name&#x27; instead. |
| `query` | string | None | Gene symbol to search for (e.g., &quot;FLT3&quot;, &quot;KRAS&quot;, &quot;TP53&quot;). Alias for &#x27;name&#x27;. |
| `gene_name` | string | None | Gene symbol to search for. Alias for &#x27;name&#x27;. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "limit": {
          "type": "integer",
          "description": "Maximum number of genes to return (default: 10, recommended max: 100)",
          "default": 10
        },
        "name": {
          "type": "string",
          "description": "Gene symbol to search for (e.g., \"EGFR\", \"BRAF\", \"BRCA1\"). Alias: use 'query' or 'gene_name' instead."
        },
        "query": {
          "type": "string",
          "description": "Gene symbol to search for (e.g., \"FLT3\", \"KRAS\", \"TP53\"). Alias for 'name'."
        },
        "gene_name": {
          "type": "string",
          "description": "Gene symbol to search for. Alias for 'name'."
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

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
                "totalCount": {
                  "type": "integer",
                  "description": "Total number of records matching the query in CIViC, independent of how many were returned in this page. Use this as the denominator; len(nodes) is only the page size."
                },
                "pageInfo": {
                  "type": "object",
                  "description": "Pagination state for this connection.",
                  "properties": {
                    "hasNextPage": {
                      "type": "boolean",
                      "description": "True when more records match than were returned; the returned nodes are an incomplete sample of totalCount."
                    }
                  }
                },
                "nodes": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "id": {
                        "type": "integer",
                        "description": "CIViC gene ID"
                      },
                      "name": {
                        "type": "string",
                        "description": "Gene symbol/name"
                      },
                      "description": {
                        "type": [
                          "string",
                          "null"
                        ],
                        "description": "Gene description"
                      },
                      "entrezId": {
                        "type": [
                          "string",
                          "integer"
                        ],
                        "description": "Entrez Gene ID"
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
            "source": {
              "type": "string"
            },
            "format": {
              "type": "string"
            }
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `civic_search_genes`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("civic_search_genes", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.civic_search_therapies` {#dnakit.databases.civic_search_therapies}

Search for therapies (drugs/treatments) in CIViC database. Returns a list of cancer therapies and drugs (with IDs and names) that are associated with variants and evidence in CIViC. Use 'query' or 'name' to filter by therapy name (e.g., 'imatinib', 'pembrolizumab').

```text
dnakit.databases.civic_search_therapies(*, limit=None, name=None, query=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `limit` | integer | None | Maximum number of therapies to return (default: 20, recommended max: 100) {&quot;default&quot;: 20} |
| `name` | string | None | Filter by therapy/drug name (e.g., &#x27;imatinib&#x27;, &#x27;ponatinib&#x27;, &#x27;pembrolizumab&#x27;). Alias: query. |
| `query` | string | None | Alias for name. Filter by therapy/drug name. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "limit": {
          "type": "integer",
          "description": "Maximum number of therapies to return (default: 20, recommended max: 100)",
          "default": 20
        },
        "name": {
          "type": "string",
          "description": "Filter by therapy/drug name (e.g., 'imatinib', 'ponatinib', 'pembrolizumab'). Alias: query."
        },
        "query": {
          "type": "string",
          "description": "Alias for name. Filter by therapy/drug name."
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "browseTherapies": {
              "type": "object",
              "properties": {
                "pageInfo": {
                  "type": "object",
                  "description": "Pagination state for this connection. This connection reports no filter-aware total, so hasNextPage is the only completeness signal available.",
                  "properties": {
                    "hasNextPage": {
                      "type": "boolean",
                      "description": "True when more records match than were returned. This connection publishes no filter-aware total, so this is the only completeness signal available."
                    }
                  }
                },
                "nodes": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "id": {
                        "type": "integer"
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

**执行与来源：** `network`；上游接口名 `civic_search_therapies`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("civic_search_therapies", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.clin_gen_dosage_by_gene` {#dnakit.databases.clin_gen_dosage_by_gene}

Get ClinGen dosage sensitivity curation for a specific gene using the JSON API. Returns haploinsufficiency (HI) and triplosensitivity (TS) assertion scores with genomic coordinates in both GRCh37 and GRCh38. Scores: 0=No evidence, 1=Little evidence, 2=Emerging evidence, 3=Sufficient evidence, 30=Gene associated with autosomal recessive phenotype, 40=Dosage sensitivity unlikely. Also includes pLI score, HI prediction score, and associated diseases with MONDO IDs.

```text
dnakit.databases.clin_gen_dosage_by_gene(*, gene=None, gene_symbol=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | None | Gene symbol to search for (e.g., &#x27;BRCA1&#x27;, &#x27;MECP2&#x27;, &#x27;TP53&#x27;). Case-insensitive. |
| `gene_symbol` | string | None | DNAKit optional alias for gene |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene": {
          "type": "string",
          "description": "Gene symbol to search for (e.g., 'BRCA1', 'MECP2', 'TP53'). Case-insensitive."
        },
        "gene_symbol": {
          "type": "string",
          "description": "DNAKit optional alias for gene"
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "array",
      "items": {
        "type": "object",
        "properties": {
          "symbol": {
            "type": [
              "string",
              "null"
            ],
            "description": "Gene symbol"
          },
          "hgnc_id": {
            "type": [
              "string",
              "null"
            ],
            "description": "HGNC ID (e.g., HGNC:1100)"
          },
          "location": {
            "type": [
              "string",
              "null"
            ],
            "description": "Cytogenetic band (e.g., 17q21.31)"
          },
          "grch37": {
            "type": [
              "string",
              "null"
            ],
            "description": "GRCh37 coordinates"
          },
          "grch38": {
            "type": [
              "string",
              "null"
            ],
            "description": "GRCh38 coordinates"
          },
          "haplo_assertion": {
            "type": [
              "integer",
              "string",
              "null"
            ],
            "description": "Haploinsufficiency score: 0=None, 1=Little, 2=Emerging, 3=Sufficient, 30=AR, 40=Unlikely"
          },
          "triplo_assertion": {
            "type": [
              "integer",
              "string",
              "null"
            ],
            "description": "Triplosensitivity score: 0=None, 1=Little, 2=Emerging, 3=Sufficient, 40=Unlikely"
          },
          "pli": {
            "type": [
              "number",
              "null"
            ],
            "description": "gnomAD pLI score"
          },
          "hi": {
            "type": [
              "number",
              "null"
            ],
            "description": "Haploinsufficiency prediction score (%)"
          },
          "plof": {
            "type": [
              "number",
              "null"
            ],
            "description": "pLoF score"
          },
          "haplo_disease": {
            "type": [
              "string",
              "null"
            ],
            "description": "Disease associated with haploinsufficiency"
          },
          "haplo_disease_id": {
            "type": [
              "string",
              "null"
            ],
            "description": "MONDO ID for haplo disease"
          },
          "triplo_disease": {
            "type": [
              "string",
              "null"
            ],
            "description": "Disease associated with triplosensitivity"
          },
          "triplo_disease_id": {
            "type": [
              "string",
              "null"
            ],
            "description": "MONDO ID for triplo disease"
          },
          "omim": {
            "type": [
              "string",
              "null"
            ],
            "description": "Whether gene has OMIM entry (Yes/No)"
          },
          "omimlink": {
            "type": [
              "string",
              "null"
            ],
            "description": "OMIM gene ID"
          },
          "morbid": {
            "type": [
              "string",
              "null"
            ],
            "description": "Whether gene is in OMIM morbid map"
          },
          "acmgsf": {
            "type": [
              "integer",
              "null"
            ],
            "description": "ACMG secondary findings gene (1=yes, 0=no)"
          },
          "locus": {
            "type": [
              "string",
              "null"
            ],
            "description": "Locus type (e.g., protein-coding gene)"
          },
          "date": {
            "type": [
              "string",
              "null"
            ],
            "description": "Last curation date"
          }
        },
        "additionalProperties": true
      }
    }
    ```

??? info "DNAKit 参数扩展规则"

    ```json
    {
      "aliases": {
        "gene_symbol": "gene"
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `ClinGen_dosage_by_gene`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("clin_gen_dosage_by_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.clin_gen_dosage_region_search` {#dnakit.databases.clin_gen_dosage_region_search}

Search ClinGen dosage sensitivity curations by genomic region. Returns all genes and recurrent CNV regions with dosage sensitivity curations overlapping the specified chromosomal coordinates. Supports both GRCh37 and GRCh38. Results include haploinsufficiency and triplosensitivity scores. Useful for interpreting the clinical significance of deletions and duplications detected by chromosomal microarray.

```text
dnakit.databases.clin_gen_dosage_region_search(*, chromosome, start, end, assembly=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `chromosome` | string | 必填（输入 Schema） | Chromosome (e.g., &#x27;17&#x27;, &#x27;X&#x27;). |
| `start` | integer | 必填（输入 Schema） | Start position. |
| `end` | integer | 必填（输入 Schema） | End position. |
| `assembly` | string | None | Genome assembly version. {&quot;enum&quot;: &#91;&quot;GRCh37&quot;, &quot;GRCh38&quot;&#93;, &quot;default&quot;: &quot;GRCh38&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "chromosome": {
          "type": "string",
          "description": "Chromosome (e.g., '17', 'X')."
        },
        "start": {
          "type": "integer",
          "description": "Start position."
        },
        "end": {
          "type": "integer",
          "description": "End position."
        },
        "assembly": {
          "type": "string",
          "enum": [
            "GRCh37",
            "GRCh38"
          ],
          "default": "GRCh38",
          "description": "Genome assembly version."
        }
      },
      "required": [
        "chromosome",
        "end",
        "start"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "array",
      "items": {
        "type": "object",
        "properties": {
          "symbol": {
            "type": [
              "string",
              "null"
            ]
          },
          "hgnc_id": {
            "type": [
              "string",
              "null"
            ]
          },
          "type": {
            "type": [
              "integer",
              "null"
            ],
            "description": "0=gene, 1=region"
          },
          "location": {
            "type": [
              "string",
              "null"
            ]
          },
          "grch37": {
            "type": [
              "string",
              "null"
            ]
          },
          "grch38": {
            "type": [
              "string",
              "null"
            ]
          },
          "haplo_assertion": {
            "type": [
              "integer",
              "string",
              "null"
            ]
          },
          "triplo_assertion": {
            "type": [
              "integer",
              "string",
              "null"
            ]
          },
          "pli": {
            "type": [
              "number",
              "null"
            ]
          },
          "hi": {
            "type": [
              "number",
              "null"
            ]
          },
          "haplo_disease": {
            "type": [
              "string",
              "null"
            ]
          },
          "triplo_disease": {
            "type": [
              "string",
              "null"
            ]
          }
        },
        "additionalProperties": true
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `ClinGen_dosage_region_search`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("clin_gen_dosage_region_search", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.clin_gen_get_actionability_adult` {#dnakit.databases.clin_gen_get_actionability_adult}

Get ClinGen clinical actionability curations for adult context. Returns genes with actionable findings, intervention scores, and outcome scores. Actionability scores help prioritize which incidental findings to return. Example: BRCA1/2 have high actionability for cancer screening.

```text
dnakit.databases.clin_gen_get_actionability_adult(*, gene=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | None | Optional: Filter by gene symbol |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene": {
          "type": "string",
          "description": "Optional: Filter by gene symbol"
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "array",
      "items": {
        "type": "object",
        "properties": {
          "docId": {
            "type": "string"
          },
          "topicIri": {
            "type": "string"
          },
          "curationType": {
            "type": "string"
          },
          "latestSearchDate": {
            "type": "string"
          },
          "lastUpdated": {
            "type": "string"
          },
          "lastAuthor": {
            "type": "string"
          },
          "context": {
            "type": "string"
          },
          "contextIri": {
            "type": "string"
          },
          "release": {
            "type": "string"
          },
          "releaseDate": {
            "type": "string"
          },
          "geneOrVariant": {
            "type": "string"
          },
          "geneOmim": {
            "type": "string"
          },
          "disease": {
            "type": "string"
          },
          "omim": {
            "type": "string"
          },
          "status-overall": {
            "type": "string"
          },
          "status-stg1": {
            "type": "string"
          },
          "status-stg2": {
            "type": "string"
          },
          "status-scoring": {
            "type": "string"
          },
          "outcome": {
            "type": [
              "string",
              "null"
            ]
          },
          "outcomeScoringGroup": {
            "type": "string"
          },
          "intervention": {
            "type": [
              "string",
              "null"
            ]
          },
          "interventionScoringGroup": {
            "type": "string"
          },
          "severity": {
            "type": [
              "string",
              "null"
            ]
          },
          "likelihood": {
            "type": [
              "string",
              "null"
            ]
          },
          "natureOfIntervention": {
            "type": [
              "string",
              "null"
            ]
          },
          "effectiveness": {
            "type": [
              "string",
              "null"
            ]
          },
          "overall": {
            "type": [
              "string",
              "null"
            ]
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `ClinGen_get_actionability_adult`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("clin_gen_get_actionability_adult", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.clin_gen_get_actionability_pediatric` {#dnakit.databases.clin_gen_get_actionability_pediatric}

Get ClinGen clinical actionability curations for pediatric context. Returns genes with actionable findings in children, intervention recommendations, and evidence scores. Pediatric curations may differ from adult context due to age-specific interventions.

```text
dnakit.databases.clin_gen_get_actionability_pediatric(*, gene=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | None | Optional: Filter by gene symbol |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene": {
          "type": "string",
          "description": "Optional: Filter by gene symbol"
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "array",
      "items": {
        "type": "object",
        "properties": {
          "docId": {
            "type": "string"
          },
          "topicIri": {
            "type": "string"
          },
          "curationType": {
            "type": "string"
          },
          "latestSearchDate": {
            "type": "string"
          },
          "lastUpdated": {
            "type": "string"
          },
          "lastAuthor": {
            "type": "string"
          },
          "context": {
            "type": "string"
          },
          "contextIri": {
            "type": "string"
          },
          "release": {
            "type": "string"
          },
          "releaseDate": {
            "type": "string"
          },
          "geneOrVariant": {
            "type": "string"
          },
          "geneOmim": {
            "type": "string"
          },
          "disease": {
            "type": "string"
          },
          "omim": {
            "type": "string"
          },
          "status-overall": {
            "type": "string"
          },
          "status-stg1": {
            "type": "string"
          },
          "status-stg2": {
            "type": "string"
          },
          "status-scoring": {
            "type": "string"
          },
          "outcome": {
            "type": [
              "string",
              "null"
            ]
          },
          "outcomeScoringGroup": {
            "type": "string"
          },
          "intervention": {
            "type": [
              "string",
              "null"
            ]
          },
          "interventionScoringGroup": {
            "type": "string"
          },
          "severity": {
            "type": [
              "string",
              "null"
            ]
          },
          "likelihood": {
            "type": [
              "string",
              "null"
            ]
          },
          "natureOfIntervention": {
            "type": [
              "string",
              "null"
            ]
          },
          "effectiveness": {
            "type": [
              "string",
              "null"
            ]
          },
          "overall": {
            "type": [
              "string",
              "null"
            ]
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `ClinGen_get_actionability_pediatric`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("clin_gen_get_actionability_pediatric", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.clin_gen_get_dosage_sensitivity` {#dnakit.databases.clin_gen_get_dosage_sensitivity}

Get all ClinGen dosage sensitivity curations. Returns haploinsufficiency and triplosensitivity scores for curated genes. Set include_regions=true to also get recurrent CNV regions. Scores: 0=No evidence, 1=Little, 2=Emerging, 3=Sufficient evidence.

```text
dnakit.databases.clin_gen_get_dosage_sensitivity(*, gene=None, include_regions=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | None | Optional: Filter by gene symbol |
| `include_regions` | boolean | None | Include recurrent CNV regions (default: false) {&quot;default&quot;: false} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene": {
          "type": "string",
          "description": "Optional: Filter by gene symbol"
        },
        "include_regions": {
          "type": "boolean",
          "description": "Include recurrent CNV regions (default: false)",
          "default": false
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "array",
      "items": {
        "type": "object"
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `ClinGen_get_dosage_sensitivity`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("clin_gen_get_dosage_sensitivity", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.clin_gen_get_gene_validity` {#dnakit.databases.clin_gen_get_gene_validity}

Get all ClinGen gene-disease validity curations. Returns comprehensive list of gene-disease relationships with classifications (Definitive, Strong, Moderate, Limited). Useful for batch analysis or exploring curated genes. Optional gene filter available.

```text
dnakit.databases.clin_gen_get_gene_validity(*, gene=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | None | Optional: Filter by gene symbol |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene": {
          "type": "string",
          "description": "Optional: Filter by gene symbol"
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "array",
      "items": {
        "type": "object"
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `ClinGen_get_gene_validity`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("clin_gen_get_gene_validity", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.clin_gen_get_variant_classifications` {#dnakit.databases.clin_gen_get_variant_classifications}

Get variant pathogenicity classifications from ClinGen Evidence Repository. Returns expert-curated variant classifications with supporting evidence. Optional filters by gene or variant. Classifications follow ACMG/AMP guidelines with expert review.

```text
dnakit.databases.clin_gen_get_variant_classifications(*, gene=None, variant=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | None | Optional: Filter by gene symbol |
| `variant` | string | None | Optional: filter to one variant. Accepts a ClinGen canonical allele ID (&#x27;CA16020993&#x27;, with or without a &#x27;CAR:&#x27; prefix), a ClinVar VariationID (&#x27;370074&#x27;), full HGVS (&#x27;NM_000277.3:c.1315+1G&gt;T&#x27;), or a protein change (&#x27;p.Arg408Trp&#x27;). The form supplied selects the matching server-side filter. HGVS and protein-change strings match on substring, so a short fragment can return many rows -- combine with &#x27;gene&#x27; to scope them. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene": {
          "type": "string",
          "description": "Optional: Filter by gene symbol"
        },
        "variant": {
          "type": "string",
          "description": "Optional: filter to one variant. Accepts a ClinGen canonical allele ID ('CA16020993', with or without a 'CAR:' prefix), a ClinVar VariationID ('370074'), full HGVS ('NM_000277.3:c.1315+1G>T'), or a protein change ('p.Arg408Trp'). The form supplied selects the matching server-side filter. HGVS and protein-change strings match on substring, so a short fragment can return many rows -- combine with 'gene' to scope them."
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "array",
      "items": {
        "type": "object",
        "properties": {
          "Variation": {
            "type": [
              "string",
              "null"
            ]
          },
          "ClinVar Variation Id": {
            "type": [
              "string",
              "null"
            ]
          },
          "HGNC Gene Symbol": {
            "type": [
              "string",
              "null"
            ]
          },
          "Disease": {
            "type": [
              "string",
              "null"
            ]
          },
          "Mondo Id": {
            "type": [
              "string",
              "null"
            ]
          },
          "Assertion": {
            "type": [
              "string",
              "null"
            ]
          },
          "Expert Panel": {
            "type": [
              "string",
              "null"
            ]
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `ClinGen_get_variant_classifications`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("clin_gen_get_variant_classifications", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.clin_gen_search_actionability` {#dnakit.databases.clin_gen_search_actionability}

Search ClinGen clinical actionability across both adult and pediatric contexts by gene. Returns actionability data from both contexts for comparison. Useful for comprehensive actionability assessment of a gene across age groups.

```text
dnakit.databases.clin_gen_search_actionability(*, gene, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | 必填（输入 Schema） | Gene symbol to search (e.g., BRCA1, MLH1, LDLR) |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene": {
          "type": "string",
          "description": "Gene symbol to search (e.g., BRCA1, MLH1, LDLR)"
        }
      },
      "required": [
        "gene"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "Adult": {
              "type": "array"
            },
            "Pediatric": {
              "type": "array"
            }
          }
        },
        "gene_searched": {
          "type": "string"
        },
        "adult_count": {
          "type": "integer"
        },
        "pediatric_count": {
          "type": "integer"
        },
        "failed_contexts": {
          "type": "object",
          "description": "Maps context name (Adult/Pediatric) to the fetch error. Present only on partial or total failure."
        },
        "note": {
          "type": "string",
          "description": "Present when some but not all contexts failed: warns that the failed context's count of 0 means 'not retrieved', not 'no curation exists'."
        },
        "source": {
          "type": "string"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `ClinGen_search_actionability`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("clin_gen_search_actionability", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.clin_gen_search_dosage_sensitivity` {#dnakit.databases.clin_gen_search_dosage_sensitivity}

Search ClinGen dosage sensitivity curations by gene symbol. Returns haploinsufficiency (HI) and triplosensitivity (TS) scores (0-3 scale). Score 3 = Sufficient evidence for dosage pathogenicity. Critical for CNV interpretation. Example: MECP2 has HI score of 3.

```text
dnakit.databases.clin_gen_search_dosage_sensitivity(*, gene, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | 必填（输入 Schema） | Gene symbol to search (e.g., MECP2, PMP22, RAI1) |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene": {
          "type": "string",
          "description": "Gene symbol to search (e.g., MECP2, PMP22, RAI1)"
        }
      },
      "required": [
        "gene"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "array",
      "items": {
        "type": "object",
        "properties": {
          "GENE SYMBOL": {
            "type": [
              "string",
              "null"
            ]
          },
          "HGNC ID": {
            "type": [
              "string",
              "null"
            ]
          },
          "HAPLOINSUFFICIENCY": {
            "type": [
              "string",
              "null"
            ]
          },
          "TRIPLOSENSITIVITY": {
            "type": [
              "string",
              "null"
            ]
          },
          "ONLINE REPORT": {
            "type": [
              "string",
              "null"
            ]
          },
          "DATE": {
            "type": [
              "string",
              "null"
            ]
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `ClinGen_search_dosage_sensitivity`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("clin_gen_search_dosage_sensitivity", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.clin_gen_search_gene_validity` {#dnakit.databases.clin_gen_search_gene_validity}

Search ClinGen gene-disease validity curations by gene symbol. Returns classification (Definitive, Strong, Moderate, Limited, Disputed, Refuted), disease associations, and inheritance patterns. Essential for ACMG variant classification. Example: BRCA1 returns validity scores for hereditary breast/ovarian cancer.

```text
dnakit.databases.clin_gen_search_gene_validity(*, gene, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | 必填（输入 Schema） | Gene symbol to search (e.g., BRCA1, TP53, CFTR) |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene": {
          "type": "string",
          "description": "Gene symbol to search (e.g., BRCA1, TP53, CFTR)"
        }
      },
      "required": [
        "gene"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "array",
      "items": {
        "type": "object",
        "properties": {
          "GENE SYMBOL": {
            "type": [
              "string",
              "null"
            ]
          },
          "GENE ID (HGNC)": {
            "type": [
              "string",
              "null"
            ]
          },
          "DISEASE LABEL": {
            "type": [
              "string",
              "null"
            ]
          },
          "DISEASE ID (MONDO)": {
            "type": [
              "string",
              "null"
            ]
          },
          "CLASSIFICATION": {
            "type": [
              "string",
              "null"
            ]
          },
          "MOI": {
            "type": [
              "string",
              "null"
            ]
          },
          "SOP": {
            "type": [
              "string",
              "null"
            ]
          },
          "GCEP": {
            "type": [
              "string",
              "null"
            ]
          },
          "ONLINE REPORT": {
            "type": [
              "string",
              "null"
            ]
          },
          "CLASSIFICATION DATE": {
            "type": [
              "string",
              "null"
            ]
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `ClinGen_search_gene_validity`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("clin_gen_search_gene_validity", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.clin_var_get_clinical_significance` {#dnakit.databases.clin_var_get_clinical_significance}

Get clinical significance information for a variant from ClinVar. Returns pathogenicity classification and clinical interpretations.

```text
dnakit.databases.clin_var_get_clinical_significance(*, variant_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant_id` | string | 必填（输入 Schema） | ClinVar variant ID (e.g., &#x27;12345&#x27;, &#x27;123456&#x27;) |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "variant_id": {
          "type": "string",
          "description": "ClinVar variant ID (e.g., '12345', '123456')"
        }
      },
      "required": [
        "variant_id"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "description": "ClinVar clinical significance data",
      "properties": {
        "status": {
          "type": "string",
          "description": "Status of the request"
        },
        "data": {
          "type": "object",
          "description": "Formatted clinical significance. Replaces the raw esummary envelope; the unprocessed record for this variant is preserved at raw_data.",
          "properties": {
            "variant_id": {
              "type": "string",
              "description": "ClinVar Variation ID that was requested"
            },
            "germline_classification": {
              "type": "object",
              "description": "Germline classification for the variant",
              "properties": {
                "description": {
                  "type": "string",
                  "description": "Classification text"
                },
                "review_status": {
                  "type": "string",
                  "description": "ClinVar review status for this classification"
                },
                "last_evaluated": {
                  "type": "string",
                  "description": "Date the classification was last evaluated"
                },
                "fda_recognized": {
                  "type": "string",
                  "description": "FDA-recognized database designation, when present"
                },
                "traits": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  },
                  "description": "Trait/condition names associated with the classification"
                }
              }
            },
            "clinical_impact": {
              "type": "object",
              "description": "Clinical impact classification (somatic variants)",
              "properties": {
                "description": {
                  "type": "string",
                  "description": "Classification text"
                },
                "review_status": {
                  "type": "string",
                  "description": "ClinVar review status for this classification"
                },
                "last_evaluated": {
                  "type": "string",
                  "description": "Date the classification was last evaluated"
                }
              }
            },
            "oncogenicity": {
              "type": "object",
              "description": "Oncogenicity classification",
              "properties": {
                "description": {
                  "type": "string",
                  "description": "Classification text"
                },
                "review_status": {
                  "type": "string",
                  "description": "ClinVar review status for this classification"
                },
                "last_evaluated": {
                  "type": "string",
                  "description": "Date the classification was last evaluated"
                }
              }
            },
            "raw_data": {
              "type": "object",
              "description": "Unprocessed NCBI esummary record for this variant"
            }
          }
        },
        "variant_id": {
          "type": "string",
          "description": "Input variant ID"
        },
        "url": {
          "type": "string",
          "description": "API URL used"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `ClinVar_get_clinical_significance`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("clin_var_get_clinical_significance", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.clin_var_get_variant_details` {#dnakit.databases.clin_var_get_variant_details}

Get variant summary information from ClinVar by variant ID. Returns accession, title, genes, clinical significance, review status, chromosome location, and variation name. For the same variant, this provides a different field view than ClinVar_get_clinical_significance.

```text
dnakit.databases.clin_var_get_variant_details(*, variant_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant_id` | string | 必填（输入 Schema） | ClinVar variant ID (e.g., &#x27;12345&#x27;, &#x27;123456&#x27;) |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "variant_id": {
          "type": "string",
          "description": "ClinVar variant ID (e.g., '12345', '123456')"
        }
      },
      "required": [
        "variant_id"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "description": "ClinVar variant details",
      "properties": {
        "status": {
          "type": "string",
          "description": "Status of the request"
        },
        "data": {
          "type": "object",
          "description": "Formatted variant details. Replaces the raw esummary envelope; the unprocessed record for this variant is preserved at raw_data.",
          "properties": {
            "variant_id": {
              "type": "string",
              "description": "ClinVar Variation ID that was requested"
            },
            "accession": {
              "type": "string",
              "description": "ClinVar accession (e.g. VCV000012345)"
            },
            "obj_type": {
              "type": "string",
              "description": "Variant object type (e.g. single nucleotide variant, Deletion)"
            },
            "chromosome": {
              "type": "string",
              "description": "Chromosome the variant maps to"
            },
            "location": {
              "type": "string",
              "description": "Cytogenetic band of the variant's first variation location"
            },
            "variation_name": {
              "type": "string",
              "description": "HGVS-style variation name (e.g. NM_000059.4(BRCA2):c.1234A>G)"
            },
            "title": {
              "type": "string",
              "description": "ClinVar record title"
            },
            "genes": {
              "type": "array",
              "items": {
                "type": "string"
              },
              "description": "HGNC symbols of genes overlapping the variant"
            },
            "clinical_significance": {
              "type": "string",
              "description": "Germline classification description (e.g. Pathogenic, Likely benign, Uncertain significance)"
            },
            "review_status": {
              "type": "string",
              "description": "ClinVar review status for the germline classification"
            },
            "raw_data": {
              "type": "object",
              "description": "Unprocessed NCBI esummary record for this variant"
            }
          }
        },
        "variant_id": {
          "type": "string",
          "description": "Input variant ID"
        },
        "url": {
          "type": "string",
          "description": "API URL used"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `ClinVar_get_variant_details`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("clin_var_get_variant_details", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.clin_var_search_variants` {#dnakit.databases.clin_var_search_variants}

Search for variants in ClinVar database by gene name, condition, variant ID, or variant name (protein change / HGVS notation). Returns variant identifiers and basic information. At least one of gene, condition, variant_id, or variant_name must be provided. Combine gene + variant_name for an exact-match lookup of one specific known variant instead of browsing all of a gene's ClinVar entries.

```text
dnakit.databases.clin_var_search_variants(*, gene=None, condition=None, variant_id=None, variant_name=None, max_results=None, limit=None, clinical_significance=None, gene_symbol=None, significance=None, query=None, variant=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | None | Gene name or symbol (e.g., &#x27;BRCA1&#x27;, &#x27;BRCA2&#x27;) At least one of gene, condition, or variant_id must be provided. |
| `condition` | string | None | Disease or condition name (e.g., &#x27;breast cancer&#x27;, &#x27;diabetes&#x27;) At least one of gene, condition, or variant_id must be provided. Matched against ClinVar&#x27;s disease index, which carries MedGen&#x27;s own condition names exactly -- a synonym, a locus name or a misspelling (&#x27;Lebers congenital amaurosis&#x27; for the indexed &#x27;Leber congenital amaurosis&#x27;) matches nothing there, so the search falls back to an unrestricted free-text match. When that happens the response says so: &#x27;condition_filter_applied&#x27;: false plus a &#x27;condition_filter_warning&#x27;, meaning total_count is NOT a disease-filtered count. Use MedGen_search_conditions to get the indexed name. |
| `variant_id` | string | None | ClinVar variant ID (e.g., &#x27;12345&#x27;) At least one of gene, condition, or variant_id must be provided. |
| `variant_name` | string &#124; array &#124; null | None | Protein change or HGVS notation to search for (e.g., &#x27;Glu6Val&#x27;, &#x27;V600E&#x27;, &#x27;c.20A&gt;T&#x27;), matched against ClinVar&#x27;s own &#91;Variant name&#93; index -- combine with &#x27;gene&#x27; for an exact-match lookup of a specific known variant instead of browsing all of a gene&#x27;s ClinVar entries. A gene-level browse can silently miss well-known but old/low-ID records once a gene has more entries than max_results (e.g. HBB&#x27;s canonical sickle-cell record, NM_000518.5(HBB):c.20A&gt;T (p.Glu7Val), isn&#x27;t within the first 100-425 gene-level rows). Accepts a list to try multiple candidate names in one call (e.g. every allele at a multi-allelic site) -- combined with OR. |
| `max_results` | integer | None | Maximum number of results to return (default: 20). Alias: limit. {&quot;default&quot;: 20, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `limit` | integer | None | Alias for max_results: maximum number of results to return. {&quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `clinical_significance` | string | None | Filter by clinical significance. Must be one of the classes ClinVar&#x27;s index actually recognizes: &#x27;Pathogenic&#x27;, &#x27;Likely pathogenic&#x27;, &#x27;Uncertain significance&#x27;, &#x27;VUS&#x27;, &#x27;Likely benign&#x27;, &#x27;Benign&#x27;, &#x27;Established risk allele&#x27;, &#x27;Likely risk allele&#x27;, &#x27;Uncertain risk allele&#x27;, &#x27;drug response&#x27;, &#x27;risk factor&#x27;, &#x27;not provided&#x27;, &#x27;other&#x27;. Case-insensitive; underscores/hyphens are treated as spaces. Combine two classes with &#x27;/&#x27; to search their union (e.g. &#x27;Pathogenic/Likely pathogenic&#x27;). Any other value is rejected with an error rather than silently filtering the results to nothing. |
| `gene_symbol` | string &#124; null | None | Alias for gene. HGNC gene symbol (e.g., &quot;DPYD&quot;, &quot;CYP2C19&quot;). |
| `significance` | string &#124; null | None | Alias for clinical_significance (e.g., &quot;pathogenic&quot;, &quot;benign&quot;, &quot;uncertain_significance&quot;). Same accepted value set and same validation as clinical_significance. |
| `query` | string &#124; null | None | Alias for condition. Free-text search mapped to condition/disease field. |
| `variant` | string &#124; array &#124; null | None | DNAKit optional alias for variant_name |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene": {
          "type": "string",
          "description": "Gene name or symbol (e.g., 'BRCA1', 'BRCA2') At least one of gene, condition, or variant_id must be provided."
        },
        "condition": {
          "type": "string",
          "description": "Disease or condition name (e.g., 'breast cancer', 'diabetes') At least one of gene, condition, or variant_id must be provided. Matched against ClinVar's disease index, which carries MedGen's own condition names exactly -- a synonym, a locus name or a misspelling ('Lebers congenital amaurosis' for the indexed 'Leber congenital amaurosis') matches nothing there, so the search falls back to an unrestricted free-text match. When that happens the response says so: 'condition_filter_applied': false plus a 'condition_filter_warning', meaning total_count is NOT a disease-filtered count. Use MedGen_search_conditions to get the indexed name."
        },
        "variant_id": {
          "type": "string",
          "description": "ClinVar variant ID (e.g., '12345') At least one of gene, condition, or variant_id must be provided."
        },
        "variant_name": {
          "type": [
            "string",
            "array",
            "null"
          ],
          "description": "Protein change or HGVS notation to search for (e.g., 'Glu6Val', 'V600E', 'c.20A>T'), matched against ClinVar's own [Variant name] index -- combine with 'gene' for an exact-match lookup of a specific known variant instead of browsing all of a gene's ClinVar entries. A gene-level browse can silently miss well-known but old/low-ID records once a gene has more entries than max_results (e.g. HBB's canonical sickle-cell record, NM_000518.5(HBB):c.20A>T (p.Glu7Val), isn't within the first 100-425 gene-level rows). Accepts a list to try multiple candidate names in one call (e.g. every allele at a multi-allelic site) -- combined with OR."
        },
        "max_results": {
          "type": "integer",
          "description": "Maximum number of results to return (default: 20). Alias: limit.",
          "default": 20,
          "minimum": 1,
          "maximum": 100
        },
        "limit": {
          "type": "integer",
          "description": "Alias for max_results: maximum number of results to return.",
          "minimum": 1,
          "maximum": 100
        },
        "clinical_significance": {
          "type": "string",
          "description": "Filter by clinical significance. Must be one of the classes ClinVar's index actually recognizes: 'Pathogenic', 'Likely pathogenic', 'Uncertain significance', 'VUS', 'Likely benign', 'Benign', 'Established risk allele', 'Likely risk allele', 'Uncertain risk allele', 'drug response', 'risk factor', 'not provided', 'other'. Case-insensitive; underscores/hyphens are treated as spaces. Combine two classes with '/' to search their union (e.g. 'Pathogenic/Likely pathogenic'). Any other value is rejected with an error rather than silently filtering the results to nothing."
        },
        "gene_symbol": {
          "type": [
            "string",
            "null"
          ],
          "description": "Alias for gene. HGNC gene symbol (e.g., \"DPYD\", \"CYP2C19\")."
        },
        "significance": {
          "type": [
            "string",
            "null"
          ],
          "description": "Alias for clinical_significance (e.g., \"pathogenic\", \"benign\", \"uncertain_significance\"). Same accepted value set and same validation as clinical_significance."
        },
        "query": {
          "type": [
            "string",
            "null"
          ],
          "description": "Alias for condition. Free-text search mapped to condition/disease field."
        },
        "variant": {
          "type": [
            "string",
            "array",
            "null"
          ],
          "description": "DNAKit optional alias for variant_name"
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "description": "ClinVar search results",
      "properties": {
        "status": {
          "type": "string",
          "description": "Status of the request"
        },
        "data": {
          "type": "object",
          "description": "Search results from ClinVar",
          "properties": {
            "esearchresult": {
              "type": "object",
              "properties": {
                "count": {
                  "type": "string",
                  "description": "Number of results found"
                },
                "idlist": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  },
                  "description": "List of ClinVar variant IDs"
                },
                "querytranslation": {
                  "type": "string",
                  "description": "Translated query"
                }
              }
            }
          }
        },
        "search_params": {
          "type": "object",
          "description": "Search parameters used",
          "properties": {
            "gene": {
              "type": "string"
            },
            "condition": {
              "type": "string"
            },
            "variant_id": {
              "type": "string"
            }
          }
        },
        "url": {
          "type": "string",
          "description": "API URL used"
        }
      }
    }
    ```

??? info "DNAKit 参数扩展规则"

    ```json
    {
      "aliases": {
        "variant": "variant_name"
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `ClinVar_search_variants`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("clin_var_search_variants", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.cosmic_get_mutations_by_gene` {#dnakit.databases.cosmic_get_mutations_by_gene}

Get all somatic mutations for a specific gene from COSMIC database. Returns comprehensive mutation data for a gene symbol including CDS changes, amino acid changes, cancer types, and genomic positions. Useful for exploring the mutational landscape of cancer genes like BRAF, TP53, EGFR, KRAS, etc.

```text
dnakit.databases.cosmic_get_mutations_by_gene(*, operation=None, gene=None, gene_name=None, max_results=None, genome_build=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_by_gene) {&quot;const&quot;: &quot;get_by_gene&quot;, &quot;default&quot;: &quot;get_by_gene&quot;} |
| `gene` | string | None | Gene symbol (e.g., BRAF, TP53, EGFR, KRAS, PIK3CA). Alias: gene_name also accepted. |
| `gene_name` | string | None | Alias for gene parameter. Gene symbol (e.g., FLT3, BRAF, TP53). |
| `max_results` | integer | None | Maximum number of mutations to return (default: 100, max: 500) {&quot;default&quot;: 100, &quot;minimum&quot;: 1, &quot;maximum&quot;: 500} |
| `genome_build` | integer | None | Genome build version: 37 (GRCh37/hg19) or 38 (GRCh38/hg38). Default: 37 {&quot;default&quot;: 37, &quot;enum&quot;: &#91;37, 38&#93;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "get_by_gene",
          "description": "Operation type (fixed: get_by_gene)",
          "default": "get_by_gene"
        },
        "gene": {
          "type": "string",
          "description": "Gene symbol (e.g., BRAF, TP53, EGFR, KRAS, PIK3CA). Alias: gene_name also accepted."
        },
        "gene_name": {
          "type": "string",
          "description": "Alias for gene parameter. Gene symbol (e.g., FLT3, BRAF, TP53)."
        },
        "max_results": {
          "type": "integer",
          "description": "Maximum number of mutations to return (default: 100, max: 500)",
          "default": 100,
          "minimum": 1,
          "maximum": 500
        },
        "genome_build": {
          "type": "integer",
          "description": "Genome build version: 37 (GRCh37/hg19) or 38 (GRCh38/hg38). Default: 37",
          "default": 37,
          "enum": [
            37,
            38
          ]
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "gene": {
          "type": "string",
          "description": "Gene symbol queried"
        },
        "total_count": {
          "type": "integer",
          "description": "Total number of mutations found"
        },
        "results": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "mutation_id": {
                "type": "string",
                "description": "COSMIC mutation ID"
              },
              "GeneName": {
                "type": "string"
              },
              "MutationCDS": {
                "type": "string"
              },
              "MutationAA": {
                "type": "string"
              },
              "MutationGenomePosition": {
                "type": "string",
                "description": "Genomic coordinates"
              },
              "PrimarySite": {
                "type": "string"
              },
              "PrimaryHistology": {
                "type": "string"
              },
              "Fathmm": {
                "type": "string",
                "description": "FATHMM pathogenicity prediction"
              }
            }
          }
        },
        "genome_build": {
          "type": "string"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `COSMIC_get_mutations_by_gene`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("cosmic_get_mutations_by_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.cosmic_search_mutations` {#dnakit.databases.cosmic_search_mutations}

Search COSMIC database for somatic mutations in cancer. COSMIC is the world's largest expert-curated database of somatic mutations in human cancer. Search by gene name, mutation ID, or mutation description. Returns mutation IDs, gene names, CDS changes, and amino acid changes. No authentication required.

```text
dnakit.databases.cosmic_search_mutations(*, operation=None, terms=None, query=None, max_results=None, genome_build=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search) {&quot;const&quot;: &quot;search&quot;, &quot;default&quot;: &quot;search&quot;} |
| `terms` | string | None | Search query - gene name (e.g., BRAF), mutation (e.g., V600E), or mutation ID (e.g., COSM476). Aliases: query or gene also accepted. |
| `query` | string | None | Alias for terms. Search query - gene name, mutation, or COSMIC ID. |
| `max_results` | integer | None | Maximum number of results to return (default: 20, max: 500) {&quot;default&quot;: 20, &quot;minimum&quot;: 1, &quot;maximum&quot;: 500} |
| `genome_build` | integer | None | Genome build version: 37 (GRCh37/hg19) or 38 (GRCh38/hg38). Default: 37 {&quot;default&quot;: 37, &quot;enum&quot;: &#91;37, 38&#93;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

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
        "terms": {
          "type": "string",
          "description": "Search query - gene name (e.g., BRAF), mutation (e.g., V600E), or mutation ID (e.g., COSM476). Aliases: query or gene also accepted."
        },
        "query": {
          "type": "string",
          "description": "Alias for terms. Search query - gene name, mutation, or COSMIC ID."
        },
        "max_results": {
          "type": "integer",
          "description": "Maximum number of results to return (default: 20, max: 500)",
          "default": 20,
          "minimum": 1,
          "maximum": 500
        },
        "genome_build": {
          "type": "integer",
          "description": "Genome build version: 37 (GRCh37/hg19) or 38 (GRCh38/hg38). Default: 37",
          "default": 37,
          "enum": [
            37,
            38
          ]
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "total_count": {
          "type": "integer",
          "description": "Total number of matching mutations"
        },
        "results": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "mutation_id": {
                "type": "string",
                "description": "COSMIC mutation ID (e.g., COSM476)"
              },
              "display": {
                "type": "array",
                "description": "Display array with mutation details [mutation_id, gene, CDS, AA]",
                "items": {
                  "type": "string"
                }
              },
              "GeneName": {
                "type": "string",
                "description": "Gene symbol"
              },
              "MutationCDS": {
                "type": "string",
                "description": "CDS-level mutation notation"
              },
              "MutationAA": {
                "type": "string",
                "description": "Amino acid change notation"
              },
              "PrimarySite": {
                "type": "string",
                "description": "Primary tissue/organ site"
              },
              "PrimaryHistology": {
                "type": "string",
                "description": "Primary histology/cancer type"
              }
            }
          }
        },
        "genome_build": {
          "type": "string",
          "description": "Genome build used (GRCh37 or GRCh38)"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `COSMIC_search_mutations`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("cosmic_search_mutations", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.dbsnp_get_frequencies` {#dnakit.databases.dbsnp_get_frequencies}

Get allele frequencies for a variant from dbSNP. Returns population-specific allele frequency data. Each study reports an 'allele_count' -- the number of ALLELES (chromosomes) observed carrying the allele, NOT a number of individuals. Do not read it as study size: for a diploid cohort the number of people is at most half of it, and the denominator is the cohort's total allele count, which dbSNP does not report here.

```text
dnakit.databases.dbsnp_get_frequencies(*, rsid, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `rsid` | string | 必填（输入 Schema） | dbSNP rsID (e.g., &#x27;rs12345&#x27;, &#x27;12345&#x27;) |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "rsid": {
          "type": "string",
          "description": "dbSNP rsID (e.g., 'rs12345', '12345')"
        }
      },
      "required": [
        "rsid"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "description": "dbSNP allele frequency data",
      "properties": {
        "status": {
          "type": "string",
          "description": "Status of the request"
        },
        "data": {
          "type": "object",
          "description": "Frequency data from dbSNP",
          "properties": {
            "frequencies": {
              "type": "array",
              "description": "Array of allele frequency data",
              "items": {
                "type": "object",
                "properties": {
                  "study": {
                    "type": "string",
                    "description": "Study/population name (e.g. '1000Genomes', 'TOPMED', 'GnomAD_genomes')"
                  },
                  "allele": {
                    "type": "string",
                    "description": "The allele this frequency refers to. Note this is NOT always the minor or alternate allele -- dbSNP reports whichever allele the study submitted, so check it before comparing studies."
                  },
                  "frequency": {
                    "type": "number",
                    "description": "Frequency of 'allele' in this study, between 0 and 1"
                  },
                  "allele_count": {
                    "type": "integer",
                    "description": "Number of ALLELES (chromosomes) observed carrying 'allele' in this study -- the numerator of 'frequency'. This is the value after the slash in dbSNP's 'global_mafs' encoding 'allele=frequency/count'. It is NOT a count of individuals: for 1000 Genomes phase 3 (2504 samples = 5008 alleles), rs429358 gives 0.1505591 x 5008 = 754 and rs80338939 gives 0.002396 x 5008 = 12, so a study reporting 12 here observed 12 carrier chromosomes out of 5008, not a 12-person cohort. dbSNP does not report the denominator, so this number alone does not give study size."
                  },
                  "sample_count": {
                    "type": "integer",
                    "description": "LEGACY MISNOMER retained for backwards compatibility -- identical in value to 'allele_count'. Despite the name this is the allele/chromosome count from dbSNP's 'global_mafs', NOT a number of individuals or a sample size. Prefer 'allele_count'."
                  }
                }
              }
            },
            "total_studies": {
              "type": "integer",
              "description": "Number of studies with frequency data"
            }
          }
        },
        "rsid": {
          "type": "string",
          "description": "Input rsID"
        },
        "url": {
          "type": "string",
          "description": "API URL used"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `dbsnp_get_frequencies`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("dbsnp_get_frequencies", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.dbsnp_get_variant_by_rsid` {#dnakit.databases.dbsnp_get_variant_by_rsid}

Get variant information from dbSNP by rsID. Returns genomic coordinates, alleles, and basic variant information. Coordinates are reported on BOTH assemblies: 'position'/'position_grch38' are on GRCh38 (dbSNP's current assembly, named by the 'assembly' field) and 'position_grch37' is the same variant on GRCh37. Always pick the field matching your pipeline's build -- the two can differ by hundreds of kilobases.

```text
dnakit.databases.dbsnp_get_variant_by_rsid(*, rsid, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `rsid` | string | 必填（输入 Schema） | dbSNP rsID (e.g., &#x27;rs12345&#x27;, &#x27;12345&#x27;) |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "rsid": {
          "type": "string",
          "description": "dbSNP rsID (e.g., 'rs12345', '12345')"
        }
      },
      "required": [
        "rsid"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "description": "dbSNP variant information",
      "properties": {
        "status": {
          "type": "string",
          "description": "Status of the request"
        },
        "data": {
          "type": "object",
          "description": "Variant information from dbSNP, flattened from the NCBI E-utilities esummary document summary",
          "properties": {
            "refsnp_id": {
              "type": "string",
              "description": "dbSNP reference ID, with the 'rs' prefix (e.g. 'rs80338939')"
            },
            "snp_id": {
              "type": [
                "integer",
                "null"
              ],
              "description": "Numeric dbSNP ID without the 'rs' prefix (e.g. 80338939)"
            },
            "chromosome": {
              "type": [
                "string",
                "null"
              ],
              "description": "Chromosome name (e.g. '13'), without a 'chr' prefix"
            },
            "position": {
              "type": [
                "string",
                "null"
              ],
              "description": "Genomic coordinate on GRCh38 as a 'chromosome:position' STRING (e.g. '13:20189547'), not an integer. Identical to 'position_grch38'; see the 'assembly' field for the build label. Retained under this name for backwards compatibility."
            },
            "assembly": {
              "type": "string",
              "description": "Genome assembly that 'position' and 'position_grch38' are expressed on. Always 'GRCh38', dbSNP's current assembly."
            },
            "position_grch38": {
              "type": [
                "string",
                "null"
              ],
              "description": "Genomic coordinate on GRCh38 as a 'chromosome:position' string, from esummary 'chrpos'. Same value as 'position', explicitly labelled."
            },
            "position_grch37": {
              "type": [
                "string",
                "null"
              ],
              "description": "Genomic coordinate for the same variant on GRCh37 (hg19) as a 'chromosome:position' string, from esummary 'chrpos_prev_assm'. Null when dbSNP has no previous-assembly mapping. Use this, not 'position', for GRCh37/hg19 pipelines."
            },
            "allele": {
              "type": [
                "string",
                "null"
              ],
              "description": "Variant allele as reported by dbSNP"
            },
            "snp_class": {
              "type": [
                "string",
                "null"
              ],
              "description": "Variant class (e.g. 'snv', 'delins', 'del')"
            },
            "clinical_significance": {
              "type": "array",
              "description": "ClinVar clinical significance terms (e.g. 'pathogenic'). A single-element list containing an empty string when dbSNP reports none.",
              "items": {
                "type": "string"
              }
            },
            "genes": {
              "type": "array",
              "description": "Symbols of genes overlapping the variant",
              "items": {
                "type": "string"
              }
            },
            "allele_frequencies": {
              "type": "array",
              "description": "Raw per-study frequency entries from esummary 'global_mafs', UNPARSED. Use dbsnp_get_frequencies for parsed numeric values.",
              "items": {
                "type": "object",
                "properties": {
                  "study": {
                    "type": "string",
                    "description": "Study/population name (e.g. '1000Genomes')"
                  },
                  "freq": {
                    "type": "string",
                    "description": "Raw dbSNP encoding 'allele=frequency/allele_count', e.g. '-=0.002396/12'. The number after the slash is an ALLELE (chromosome) count, not a number of individuals."
                  }
                }
              }
            },
            "hgvs_notation": {
              "type": [
                "string",
                "null"
              ],
              "description": "dbSNP 'docsum' string. Not bare HGVS: a pipe-delimited record whose first segment is 'HGVS=' followed by a comma-separated list of HGVS expressions, then SEQ=/LEN=/GENE= segments. Parse before use."
            },
            "spdi_notation": {
              "type": [
                "string",
                "null"
              ],
              "description": "Comma-separated SPDI expressions (sequence:position:deleted:inserted)"
            },
            "function_class": {
              "type": "array",
              "description": "Molecular consequence terms (e.g. 'frameshift_variant')",
              "items": {
                "type": "string"
              }
            },
            "validated": {
              "type": "array",
              "description": "dbSNP validation status terms (e.g. 'by-frequency')",
              "items": {
                "type": "string"
              }
            },
            "created_date": {
              "type": [
                "string",
                "null"
              ],
              "description": "Record creation date, 'YYYY/MM/DD HH:MM'"
            },
            "updated_date": {
              "type": [
                "string",
                "null"
              ],
              "description": "Record last-update date, 'YYYY/MM/DD HH:MM'"
            }
          }
        },
        "rsid": {
          "type": "string",
          "description": "Input rsID, normalized to include the 'rs' prefix"
        },
        "url": {
          "type": "string",
          "description": "API URL used (NCBI E-utilities esummary.fcgi)"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `dbsnp_get_variant_by_rsid`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("dbsnp_get_variant_by_rsid", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ensembl_var_get_variant_consequences` {#dnakit.databases.ensembl_var_get_variant_consequences}

Get detailed variant information including consequence type, clinical significance, synonyms, evidence attributes, and allele mappings from Ensembl. More comprehensive than VEP - returns the full variant record with all known synonyms, mapped positions, ancestral allele, and minor allele data. Example: rs429358 returns most_severe_consequence='missense_variant', 27 synonyms, mapped to chr19:44908684 (GRCh38), and is classified as a missense_variant affecting APOE.

```text
dnakit.databases.ensembl_var_get_variant_consequences(*, variant_id, species=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant_id` | string | 必填（输入 Schema） | dbSNP rsID of the variant. Examples: &#x27;rs429358&#x27; (APOE), &#x27;rs7903146&#x27; (TCF7L2), &#x27;rs1042779&#x27;. |
| `species` | string | None | Species name. Default &#x27;human&#x27;. Examples: &#x27;human&#x27;, &#x27;homo_sapiens&#x27;. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "variant_id": {
          "type": "string",
          "description": "dbSNP rsID of the variant. Examples: 'rs429358' (APOE), 'rs7903146' (TCF7L2), 'rs1042779'."
        },
        "species": {
          "type": "string",
          "description": "Species name. Default 'human'. Examples: 'human', 'homo_sapiens'."
        }
      },
      "required": [
        "variant_id"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "name": {
              "type": "string"
            },
            "source": {
              "type": [
                "string",
                "null"
              ]
            },
            "most_severe_consequence": {
              "type": [
                "string",
                "null"
              ]
            },
            "ancestral_allele": {
              "type": [
                "string",
                "null"
              ]
            },
            "minor_allele": {
              "type": [
                "string",
                "null"
              ]
            },
            "MAF": {
              "type": [
                "number",
                "null"
              ]
            },
            "synonyms": {
              "type": "array",
              "items": {
                "type": "string"
              },
              "description": "Alternative variant identifiers"
            },
            "mappings": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "seq_region_name": {
                    "type": "string"
                  },
                  "start": {
                    "type": "integer"
                  },
                  "end": {
                    "type": "integer"
                  },
                  "allele_string": {
                    "type": "string"
                  },
                  "strand": {
                    "type": "integer"
                  },
                  "assembly_name": {
                    "type": "string"
                  }
                }
              }
            },
            "evidence": {
              "type": "array",
              "items": {
                "type": "string"
              },
              "description": "Evidence tags (e.g., 1000Genomes, gnomAD, Cited)"
            },
            "clinical_significance": {
              "type": "array",
              "items": {
                "type": "string"
              },
              "description": "Clinical significance annotations"
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

**执行与来源：** `network`；上游接口名 `EnsemblVar_get_variant_consequences`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ensembl_var_get_variant_consequences", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ensembl_vep_annotate_hgvs` {#dnakit.databases.ensembl_vep_annotate_hgvs}

Predict functional consequences of a genetic variant using HGVS notation via the Ensembl Variant Effect Predictor (VEP). Returns impact on genes, transcripts, and proteins including SIFT/PolyPhen pathogenicity scores, consequence types, amino acid changes, and colocated known variants. HGVS notation examples: 'BRAF:p.Val600Glu' (protein), 'ENST00000366667:c.803C>T' (coding DNA). This complements the existing region-based VEP tool by accepting standard clinical variant nomenclature.

```text
dnakit.databases.ensembl_vep_annotate_hgvs(*, hgvs_notation, species=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `hgvs_notation` | string | 必填（输入 Schema） | Variant in HGVS notation. Supports protein (p.), coding DNA (c.), and genomic (g.) notation. Examples: &#x27;BRAF:p.Val600Glu&#x27;, &#x27;ENST00000366667:c.803C&gt;T&#x27;, &#x27;NC_000007.14:g.140753336A&gt;T&#x27;. |
| `species` | string | None | Species name. Default: &#x27;human&#x27;. Other options: &#x27;mouse&#x27;, &#x27;rat&#x27;, &#x27;zebrafish&#x27;, etc. {&quot;default&quot;: &quot;human&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "hgvs_notation": {
          "type": "string",
          "description": "Variant in HGVS notation. Supports protein (p.), coding DNA (c.), and genomic (g.) notation. Examples: 'BRAF:p.Val600Glu', 'ENST00000366667:c.803C>T', 'NC_000007.14:g.140753336A>T'."
        },
        "species": {
          "type": "string",
          "description": "Species name. Default: 'human'. Other options: 'mouse', 'rat', 'zebrafish', etc.",
          "default": "human"
        }
      },
      "required": [
        "hgvs_notation"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "input": {
              "type": [
                "string",
                "null"
              ],
              "description": "Input variant identifier"
            },
            "assembly_name": {
              "type": [
                "string",
                "null"
              ],
              "description": "Reference genome assembly"
            },
            "seq_region_name": {
              "type": [
                "string",
                "null"
              ],
              "description": "Chromosome or contig"
            },
            "start": {
              "type": [
                "integer",
                "null"
              ],
              "description": "Genomic start position"
            },
            "end": {
              "type": [
                "integer",
                "null"
              ],
              "description": "Genomic end position"
            },
            "strand": {
              "type": [
                "integer",
                "null"
              ],
              "description": "Strand (1 or -1)"
            },
            "allele_string": {
              "type": [
                "string",
                "null"
              ],
              "description": "Allele string (ref/alt)"
            },
            "most_severe_consequence": {
              "type": [
                "string",
                "null"
              ],
              "description": "Most severe consequence (e.g., missense_variant)"
            },
            "transcript_consequences": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "gene_symbol": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "gene_id": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "transcript_id": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "biotype": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "variant_allele": {
                    "type": [
                      "string",
                      "null"
                    ],
                    "description": "The alternate allele this consequence is predicted for. For multi-allelic variants the same transcript appears once per allele, so this field is required to attribute a consequence to the correct allele."
                  },
                  "consequence_terms": {
                    "type": "array",
                    "items": {
                      "type": "string"
                    }
                  },
                  "impact": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "amino_acids": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "codons": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "sift_prediction": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "sift_score": {
                    "type": [
                      "number",
                      "null"
                    ]
                  },
                  "polyphen_prediction": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "polyphen_score": {
                    "type": [
                      "number",
                      "null"
                    ]
                  }
                }
              }
            },
            "colocated_variants": {
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
                  "allele_string": {
                    "type": [
                      "string",
                      "null"
                    ]
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

**执行与来源：** `network`；上游接口名 `EnsemblVEP_annotate_hgvs`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ensembl_vep_annotate_hgvs", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ensembl_vep_annotate_rsid` {#dnakit.databases.ensembl_vep_annotate_rsid}

Predict functional consequences of a genetic variant using its dbSNP rs identifier via the Ensembl Variant Effect Predictor (VEP). Returns impact on genes, transcripts, and proteins including SIFT/PolyPhen scores, consequence types, amino acid changes, and colocated variants. Useful when you have a variant's rs number from GWAS, clinical reports, or literature. Example: 'rs7903146' (TCF7L2 T2D variant), 'rs429358' (APOE Alzheimer variant).

```text
dnakit.databases.ensembl_vep_annotate_rsid(*, variant_id=None, rsid=None, rs_id=None, species=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant_id` | string | None | dbSNP rs identifier. Examples: &#x27;rs7903146&#x27; (TCF7L2, type 2 diabetes), &#x27;rs429358&#x27; (APOE, Alzheimer&#x27;s), &#x27;rs1801133&#x27; (MTHFR C677T). Aliases: rsid, rs_id. |
| `rsid` | string | None | Alias for variant_id. dbSNP rs identifier (e.g., &#x27;rs7903146&#x27;). |
| `rs_id` | string | None | Alias for variant_id. dbSNP rs identifier (e.g., &#x27;rs7903146&#x27;). |
| `species` | string | None | Species name. Default: &#x27;human&#x27;. {&quot;default&quot;: &quot;human&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "variant_id": {
          "type": "string",
          "description": "dbSNP rs identifier. Examples: 'rs7903146' (TCF7L2, type 2 diabetes), 'rs429358' (APOE, Alzheimer's), 'rs1801133' (MTHFR C677T). Aliases: rsid, rs_id."
        },
        "rsid": {
          "type": "string",
          "description": "Alias for variant_id. dbSNP rs identifier (e.g., 'rs7903146')."
        },
        "rs_id": {
          "type": "string",
          "description": "Alias for variant_id. dbSNP rs identifier (e.g., 'rs7903146')."
        },
        "species": {
          "type": "string",
          "description": "Species name. Default: 'human'.",
          "default": "human"
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "input": {
              "type": [
                "string",
                "null"
              ]
            },
            "assembly_name": {
              "type": [
                "string",
                "null"
              ]
            },
            "seq_region_name": {
              "type": [
                "string",
                "null"
              ]
            },
            "start": {
              "type": [
                "integer",
                "null"
              ]
            },
            "end": {
              "type": [
                "integer",
                "null"
              ]
            },
            "allele_string": {
              "type": [
                "string",
                "null"
              ]
            },
            "most_severe_consequence": {
              "type": [
                "string",
                "null"
              ]
            },
            "transcript_consequences": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "gene_symbol": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "variant_allele": {
                    "type": [
                      "string",
                      "null"
                    ],
                    "description": "The alternate allele this consequence is predicted for. For multi-allelic variants the same transcript appears once per allele, so this field is required to attribute a consequence to the correct allele."
                  },
                  "consequence_terms": {
                    "type": "array",
                    "items": {
                      "type": "string"
                    }
                  },
                  "impact": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "amino_acids": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "sift_prediction": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "polyphen_prediction": {
                    "type": [
                      "string",
                      "null"
                    ]
                  }
                }
              }
            },
            "colocated_variants": {
              "type": "array"
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

**执行与来源：** `network`；上游接口名 `EnsemblVEP_annotate_rsid`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ensembl_vep_annotate_rsid", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ensembl_vep_variant_recoder` {#dnakit.databases.ensembl_vep_variant_recoder}

Convert genetic variant identifiers between different nomenclature formats using the Ensembl Variant Recoder. Given a variant ID (e.g., rs number), returns equivalent representations in HGVS genomic (hgvsg), HGVS coding (hgvsc), HGVS protein (hgvsp), and SPDI format. Essential for cross-referencing variants across databases that use different notation systems. Example: 'rs429358' returns HGVS coding notation for all APOE transcripts.

```text
dnakit.databases.ensembl_vep_variant_recoder(*, variant_id, species=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant_id` | string | 必填（输入 Schema） | Variant identifier to recode. Accepts: rsID (e.g., &#x27;rs429358&#x27;), HGVS notation (e.g., &#x27;AGT:c.803T&gt;C&#x27;), SPDI (e.g., &#x27;NC_000019.10:44908683:T:C&#x27;). |
| `species` | string | None | Species name. Default: &#x27;human&#x27;. {&quot;default&quot;: &quot;human&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "variant_id": {
          "type": "string",
          "description": "Variant identifier to recode. Accepts: rsID (e.g., 'rs429358'), HGVS notation (e.g., 'AGT:c.803T>C'), SPDI (e.g., 'NC_000019.10:44908683:T:C')."
        },
        "species": {
          "type": "string",
          "description": "Species name. Default: 'human'.",
          "default": "human"
        }
      },
      "required": [
        "variant_id"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "allele": {
                "type": "string",
                "description": "Alternative allele"
              },
              "input": {
                "type": "string",
                "description": "Original input variant ID"
              },
              "id": {
                "type": "array",
                "items": {
                  "type": "string"
                },
                "description": "Known variant IDs (rsIDs, COSMIC, etc.)"
              },
              "hgvsg": {
                "type": "array",
                "items": {
                  "type": "string"
                },
                "description": "HGVS genomic notation(s)"
              },
              "hgvsc": {
                "type": "array",
                "items": {
                  "type": "string"
                },
                "description": "HGVS coding notation(s) for each transcript"
              },
              "hgvsp": {
                "type": "array",
                "items": {
                  "type": "string"
                },
                "description": "HGVS protein notation(s)"
              },
              "spdi": {
                "type": "array",
                "items": {
                  "type": "string"
                },
                "description": "SPDI notation(s)"
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

**执行与来源：** `network`；上游接口名 `EnsemblVEP_variant_recoder`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ensembl_vep_variant_recoder", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.eve_get_gene_info` {#dnakit.databases.eve_get_gene_info}

Check if EVE scores are available for a gene and get gene information. EVE covers ~3,000 disease-related human genes. Returns Ensembl ID and link to EVE website.

```text
dnakit.databases.eve_get_gene_info(*, gene_symbol, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 必填（输入 Schema） | Gene symbol (e.g., &#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;EGFR&#x27;) |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene_symbol": {
          "type": "string",
          "description": "Gene symbol (e.g., 'TP53', 'BRCA1', 'EGFR')"
        }
      },
      "required": [
        "gene_symbol"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

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
            "gene_symbol": {
              "type": "string"
            },
            "ensembl_id": {
              "type": "string"
            },
            "eve_website": {
              "type": "string"
            }
          }
        },
        "error": {
          "type": "string"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `EVE_get_gene_info`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("eve_get_gene_info", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.eve_get_variant_score` {#dnakit.databases.eve_get_variant_score}

Get EVE (Evolutionary Variant Effect) pathogenicity score for a variant. EVE uses unsupervised deep learning on evolutionary data. Scores: 0-1 where >0.5 is likely pathogenic. Accessed via Ensembl VEP API.

```text
dnakit.databases.eve_get_variant_score(*, variant=None, chrom=None, pos=None, ref=None, alt=None, species=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant` | string | None | Variant in HGVS format (e.g., &#x27;ENST00000269305.4:c.100G&gt;A&#x27;, &#x27;NM_000546.5:c.215C&gt;G&#x27;) |
| `chrom` | string | None | Chromosome (1-22, X, Y). Use with pos, ref, alt instead of variant. |
| `pos` | integer | None | Genomic position (GRCh38) |
| `ref` | string | None | Reference allele |
| `alt` | string | None | Alternate allele |
| `species` | string | None | Species (default: human) {&quot;default&quot;: &quot;human&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "variant": {
          "type": "string",
          "description": "Variant in HGVS format (e.g., 'ENST00000269305.4:c.100G>A', 'NM_000546.5:c.215C>G')"
        },
        "chrom": {
          "type": "string",
          "description": "Chromosome (1-22, X, Y). Use with pos, ref, alt instead of variant."
        },
        "pos": {
          "type": "integer",
          "description": "Genomic position (GRCh38)"
        },
        "ref": {
          "type": "string",
          "description": "Reference allele"
        },
        "alt": {
          "type": "string",
          "description": "Alternate allele"
        },
        "species": {
          "type": "string",
          "default": "human",
          "description": "Species (default: human)"
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

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
            "variant": {
              "type": "string"
            },
            "eve_scores": {
              "type": "array"
            },
            "threshold": {
              "type": "object"
            }
          }
        },
        "error": {
          "type": "string"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `EVE_get_variant_score`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("eve_get_variant_score", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.favor_annotate_variant` {#dnakit.databases.favor_annotate_variant}

Comprehensive functional annotation of a single GRCh38 variant via FAVOR (Functional Annotation of Variants Online Resource, Harvard). One call returns: allele frequencies (BRAVO/TOPMed, gnomAD by ancestry, 1000 Genomes), gene and consequence (GENCODE category, exonic effect, protein variant, HGVS), deleteriousness scores (CADD, SIFT, PolyPhen-2, AlphaMissense, MetaSVM, MutationTaster/Assessor, FATHMM-XF, Grantham, LINSIGHT, FunSeq), conservation (GERP, phyloP mammalian/primate/vertebrate, phastCons), ClinVar clinical significance/disease, and regulatory/epigenomic context (CAGE promoter/enhancer, GeneHancer, ENCODE DNase, ReMap TF overlap). Use for variant interpretation across both coding and non-coding regions in one step. Input is a GRCh38 (hg38) coordinate string; resolve rsIDs to GRCh38 coordinates first if needed. No API key required.

```text
dnakit.databases.favor_annotate_variant(*, variant, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant` | string | 必填（输入 Schema） | GRCh38 SNV/indel as chr-pos-ref-alt. Accepts &#x27;19-44908822-C-T&#x27;, &#x27;chr19:44908822:C:T&#x27;, etc. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "variant": {
          "type": "string",
          "description": "GRCh38 SNV/indel as chr-pos-ref-alt. Accepts '19-44908822-C-T', 'chr19:44908822:C:T', etc."
        }
      },
      "required": [
        "variant"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "found": {
              "type": "boolean"
            },
            "variant": {
              "type": "object",
              "description": "Identity: variant_vcf, rsid, chromosome, position, HGVS genomic"
            },
            "gene_consequence": {
              "type": "object",
              "description": "Gene, GENCODE category, exonic effect, protein variant, HGVS c/p"
            },
            "allele_frequency": {
              "type": "object",
              "description": "BRAVO/TOPMed, gnomAD by ancestry, 1000 Genomes frequencies"
            },
            "deleteriousness": {
              "type": "object",
              "description": "CADD, SIFT, PolyPhen, AlphaMissense, MetaSVM and other pathogenicity scores"
            },
            "conservation": {
              "type": "object",
              "description": "GERP, phyloP, phastCons conservation scores"
            },
            "clinical": {
              "type": "object",
              "description": "ClinVar significance, disease, review status"
            },
            "regulatory": {
              "type": "object",
              "description": "CAGE, GeneHancer, ENCODE, ReMap regulatory context"
            },
            "all_annotations": {
              "type": "object",
              "description": "Full ~180-field FAVOR record"
            }
          },
          "required": [
            "found"
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

**执行与来源：** `network`；上游接口名 `FAVOR_annotate_variant`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("favor_annotate_variant", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.finn_gen_get_variant_finemapping` {#dnakit.databases.finn_gen_get_variant_finemapping}

Get fine-mapping regions associated with a genomic variant in FinnGen. Returns phenotypes where this variant appears in a credible set from statistical fine-mapping. Input variant in chr:pos:ref:alt format (GRCh38). Useful for understanding which disease associations a variant contributes to in the Finnish population.

```text
dnakit.databases.finn_gen_get_variant_finemapping(*, variant, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant` | string | 必填（输入 Schema） | Genomic variant in chr:pos:ref:alt format (GRCh38). Examples: &#x27;19:44908684:T:C&#x27; (APOE rs429358), &#x27;10:112998590:C:T&#x27; (TCF7L2 rs7903146). Also accepts chr-pos-ref-alt format. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "variant": {
          "type": "string",
          "description": "Genomic variant in chr:pos:ref:alt format (GRCh38). Examples: '19:44908684:T:C' (APOE rs429358), '10:112998590:C:T' (TCF7L2 rs7903146). Also accepts chr-pos-ref-alt format."
        }
      },
      "required": [
        "variant"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "array",
      "items": {
        "type": "object",
        "properties": {
          "phenocode": {
            "type": "string"
          },
          "chromosome": {
            "type": "integer"
          },
          "start": {
            "type": "integer"
          },
          "end": {
            "type": "integer"
          },
          "type": {
            "type": "string"
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `FinnGen_get_variant_finemapping`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("finn_gen_get_variant_finemapping", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.genebass_gene_burden_phewas` {#dnakit.databases.genebass_gene_burden_phewas}

Exome-wide gene-level burden PheWAS from Genebass (UK Biobank ~394K exomes, SAIGE-GENE+). Given one gene, returns its rare-variant burden association across ~4,500 phenotypes (GRCh38), sorted by p-value. Aggregates rare coding variants into a single gene-level test, complementing single-variant PheWAS. Accepts an Ensembl gene ID (ENSG...) or a gene symbol (auto-resolved via Ensembl). Choose a burden_set: pLoF (predicted loss-of-function, default), missense|LC, or synonymous.

```text
dnakit.databases.genebass_gene_burden_phewas(*, gene, burden_set=None, limit=None, max_pval=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | 必填（输入 Schema） | Ensembl gene ID (e.g. &#x27;ENSG00000148737&#x27;) or gene symbol (e.g. &#x27;TCF7L2&#x27;). Symbols are auto-resolved via Ensembl. |
| `burden_set` | string | None | Variant category aggregated in the burden test. One of: &#x27;pLoF&#x27; (default), &#x27;missense&#124;LC&#x27;, &#x27;synonymous&#x27;. Aliases accepted: lof/plof -&gt; pLoF, missense -&gt; missense&#124;LC. {&quot;default&quot;: &quot;pLoF&quot;} |
| `limit` | integer | None | Maximum number of phenotype associations to return (sorted by ascending burden p-value). {&quot;default&quot;: 25} |
| `max_pval` | number &#124; null | None | Optional: only return associations with burden p-value &lt;= this threshold (e.g. 2.5e-6 for exome-wide significance). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene": {
          "type": "string",
          "description": "Ensembl gene ID (e.g. 'ENSG00000148737') or gene symbol (e.g. 'TCF7L2'). Symbols are auto-resolved via Ensembl."
        },
        "burden_set": {
          "type": "string",
          "description": "Variant category aggregated in the burden test. One of: 'pLoF' (default), 'missense|LC', 'synonymous'. Aliases accepted: lof/plof -> pLoF, missense -> missense|LC.",
          "default": "pLoF"
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of phenotype associations to return (sorted by ascending burden p-value).",
          "default": 25
        },
        "max_pval": {
          "type": [
            "number",
            "null"
          ],
          "description": "Optional: only return associations with burden p-value <= this threshold (e.g. 2.5e-6 for exome-wide significance)."
        }
      },
      "required": [
        "gene"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "gene_id": {
              "type": "string"
            },
            "gene_symbol": {
              "type": [
                "string",
                "null"
              ]
            },
            "burden_set": {
              "type": "string"
            },
            "associations": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "phenocode": {
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
                  "trait_type": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "pval": {
                    "type": [
                      "number",
                      "null"
                    ]
                  },
                  "pval_skat": {
                    "type": [
                      "number",
                      "null"
                    ]
                  },
                  "beta": {
                    "type": [
                      "number",
                      "null"
                    ]
                  },
                  "total_variants": {
                    "type": [
                      "number",
                      "null"
                    ]
                  }
                }
              }
            }
          },
          "required": [
            "gene_id"
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

**执行与来源：** `network`；上游接口名 `Genebass_gene_burden_phewas`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("genebass_gene_burden_phewas", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gnomad_get_gene_constraints` {#dnakit.databases.gnomad_get_gene_constraints}

Get gene constraint metrics from gnomAD for a gene symbol (GRCh38). Returns constraint fields under `exac_constraint` and `gnomad_constraint` (e.g., pLI, oe_lof/oe_mis/oe_syn, exp/obs counts). Mitochondrially-encoded genes (MT-*) are scored under `mitochondrial_constraint` instead -- their `gnomad_constraint` and `exac_constraint` are always null, which does not mean unconstrained. Note: does NOT include LOEUF, missense Z-score, or synonymous Z-score -- for the full constraint panel (pLI + LOEUF + mis_z + syn_z), use gnomad_get_constraint instead.

```text
dnakit.databases.gnomad_get_gene_constraints(*, gene_symbol, reference_genome=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 必填（输入 Schema） | Gene symbol (e.g., &#x27;BRCA1&#x27;, &#x27;TP53&#x27;) |
| `reference_genome` | string | None | Reference genome. {&quot;enum&quot;: &#91;&quot;GRCh37&quot;, &quot;GRCh38&quot;&#93;, &quot;default&quot;: &quot;GRCh38&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene_symbol": {
          "type": "string",
          "description": "Gene symbol (e.g., 'BRCA1', 'TP53')"
        },
        "reference_genome": {
          "type": "string",
          "enum": [
            "GRCh37",
            "GRCh38"
          ],
          "default": "GRCh38",
          "description": "Reference genome."
        }
      },
      "required": [
        "gene_symbol"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "gene": {
              "type": "object",
              "properties": {
                "symbol": {
                  "type": "string"
                },
                "gene_id": {
                  "type": "string"
                },
                "exac_constraint": {
                  "type": [
                    "object",
                    "null"
                  ]
                },
                "gnomad_constraint": {
                  "type": [
                    "object",
                    "null"
                  ],
                  "properties": {
                    "exp_lof": {
                      "type": "number"
                    },
                    "obs_lof": {
                      "type": [
                        "integer",
                        "number"
                      ]
                    },
                    "oe_lof": {
                      "type": "number"
                    },
                    "pLI": {
                      "type": "number"
                    },
                    "exp_mis": {
                      "type": "number"
                    },
                    "obs_mis": {
                      "type": [
                        "integer",
                        "number"
                      ]
                    },
                    "oe_mis": {
                      "type": "number"
                    },
                    "exp_syn": {
                      "type": "number"
                    },
                    "obs_syn": {
                      "type": [
                        "integer",
                        "number"
                      ]
                    },
                    "oe_syn": {
                      "type": "number"
                    }
                  }
                },
                "mitochondrial_constraint": {
                  "type": [
                    "object",
                    "null"
                  ],
                  "description": "Constraint for mitochondrially-encoded genes, which gnomAD scores here and never under gnomad_constraint/exac_constraint. __typename is ProteinMitochondrialGeneConstraint (exp_lof, obs_lof, oe_lof, oe_lof_upper) for protein-coding mtDNA genes such as MT-ND4, or RNAMitochondrialGeneConstraint (observed, expected, oe, oe_upper) for mt-tRNA/mt-rRNA genes such as MT-TL1. Null for nuclear genes."
                }
              }
            },
            "note": {
              "type": "string",
              "description": "Present only when gnomAD scored the gene under none of gnomad_constraint, exac_constraint or mitochondrial_constraint -- an absence of data rather than evidence that the gene is unconstrained."
            }
          },
          "required": [
            "gene"
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

**执行与来源：** `network`；上游接口名 `gnomad_get_gene_constraints`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gnomad_get_gene_constraints", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gnomad_get_sv_by_gene` {#dnakit.databases.gnomad_get_sv_by_gene}

Get structural variants (SVs) from gnomAD for a gene. The SV callset follows the requested assembly: GRCh38 (default) uses the gnomAD v4 SV callset (gnomad_sv_r4), GRCh37 uses the gnomAD v2.1 SV callset (gnomad_sv_r2_1), so SV coordinates always match the gene coordinates returned. Returns deletions, duplications, inversions, translocations, and complex SVs with population allele frequencies, SV type, consequence, and homozygote counts (null for multi-allelic CNVs). Essential for interpreting CNVs in clinical genomics.

```text
dnakit.databases.gnomad_get_sv_by_gene(*, gene_symbol, reference_genome=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 必填（输入 Schema） | Gene symbol (e.g., &#x27;BRCA1&#x27;, &#x27;TP53&#x27;, &#x27;MECP2&#x27;) |
| `reference_genome` | string | None | Reference genome assembly. {&quot;enum&quot;: &#91;&quot;GRCh37&quot;, &quot;GRCh38&quot;&#93;, &quot;default&quot;: &quot;GRCh38&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene_symbol": {
          "type": "string",
          "description": "Gene symbol (e.g., 'BRCA1', 'TP53', 'MECP2')"
        },
        "reference_genome": {
          "type": "string",
          "enum": [
            "GRCh37",
            "GRCh38"
          ],
          "default": "GRCh38",
          "description": "Reference genome assembly."
        }
      },
      "required": [
        "gene_symbol"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "gene": {
              "type": "object",
              "properties": {
                "symbol": {
                  "type": "string"
                },
                "gene_id": {
                  "type": "string"
                },
                "chrom": {
                  "type": "string"
                },
                "start": {
                  "type": [
                    "integer",
                    "number"
                  ]
                },
                "stop": {
                  "type": [
                    "integer",
                    "number"
                  ]
                },
                "structural_variants": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "variant_id": {
                        "type": "string"
                      },
                      "chrom": {
                        "type": "string"
                      },
                      "pos": {
                        "type": [
                          "integer",
                          "number"
                        ]
                      },
                      "end": {
                        "type": [
                          "integer",
                          "number"
                        ]
                      },
                      "type": {
                        "type": "string"
                      },
                      "length": {
                        "type": [
                          "integer",
                          "number"
                        ]
                      },
                      "ac": {
                        "type": [
                          "integer",
                          "number"
                        ]
                      },
                      "an": {
                        "type": [
                          "integer",
                          "number"
                        ]
                      },
                      "af": {
                        "type": "number"
                      },
                      "homozygote_count": {
                        "type": [
                          "integer",
                          "number",
                          "null"
                        ]
                      },
                      "major_consequence": {
                        "type": [
                          "string",
                          "null"
                        ]
                      },
                      "filters": {
                        "type": "array",
                        "items": {
                          "type": "string"
                        }
                      }
                    }
                  }
                }
              }
            },
            "dataset": {
              "type": "string",
              "description": "The gnomAD dataset that actually answered, including when it came from this tool's default."
            },
            "reference_genome": {
              "type": "string",
              "description": "Reference assembly the dataset is built on (gnomad_r2_1*/exac = GRCh37; gnomad_r3*/gnomad_r4* = GRCh38)."
            },
            "dataset_note": {
              "type": "string",
              "description": "Present only when the dataset structurally lacks a callset the payload reports as null (e.g. gnomad_r3 has no exomes), so `null` is not misread as 'not observed'."
            }
          },
          "required": [
            "gene"
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

**执行与来源：** `network`；上游接口名 `gnomad_get_sv_by_gene`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gnomad_get_sv_by_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gnomad_get_sv_by_region` {#dnakit.databases.gnomad_get_sv_by_region}

Get structural variants (SVs) from gnomAD overlapping a genomic region. The SV callset follows the requested assembly: GRCh38 (default) uses the gnomAD v4 SV callset (gnomad_sv_r4), GRCh37 uses the gnomAD v2.1 SV callset (gnomad_sv_r2_1), so SV coordinates are always in the same assembly as the requested region. Returns DEL/DUP/INV/BND/CPX variants with allele frequencies and consequences. Use for CNV interpretation in specific chromosomal regions. Region should be < 5Mb.

```text
dnakit.databases.gnomad_get_sv_by_region(*, chrom, start, stop=None, reference_genome=None, end=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `chrom` | string | 必填（输入 Schema） | Chromosome (e.g., &#x27;17&#x27;, &#x27;X&#x27;). |
| `start` | integer | 必填（输入 Schema） | Start position (1-based, GRCh38 by default). |
| `stop` | integer | None | Stop position. |
| `reference_genome` | string | None | Reference genome assembly. {&quot;enum&quot;: &#91;&quot;GRCh37&quot;, &quot;GRCh38&quot;&#93;, &quot;default&quot;: &quot;GRCh38&quot;} |
| `end` | integer | None | DNAKit optional alias for stop |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "chrom": {
          "type": "string",
          "description": "Chromosome (e.g., '17', 'X')."
        },
        "start": {
          "type": "integer",
          "description": "Start position (1-based, GRCh38 by default)."
        },
        "stop": {
          "type": "integer",
          "description": "Stop position."
        },
        "reference_genome": {
          "type": "string",
          "enum": [
            "GRCh37",
            "GRCh38"
          ],
          "default": "GRCh38",
          "description": "Reference genome assembly."
        },
        "end": {
          "type": "integer",
          "description": "DNAKit optional alias for stop"
        }
      },
      "required": [
        "chrom",
        "start"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "region": {
              "type": "object",
              "properties": {
                "chrom": {
                  "type": "string"
                },
                "start": {
                  "type": [
                    "integer",
                    "number"
                  ]
                },
                "stop": {
                  "type": [
                    "integer",
                    "number"
                  ]
                },
                "genes": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "gene_id": {
                        "type": "string"
                      },
                      "symbol": {
                        "type": "string"
                      }
                    }
                  }
                },
                "structural_variants": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "variant_id": {
                        "type": "string"
                      },
                      "chrom": {
                        "type": "string"
                      },
                      "pos": {
                        "type": [
                          "integer",
                          "number"
                        ]
                      },
                      "end": {
                        "type": [
                          "integer",
                          "number"
                        ]
                      },
                      "type": {
                        "type": "string"
                      },
                      "length": {
                        "type": [
                          "integer",
                          "number"
                        ]
                      },
                      "ac": {
                        "type": [
                          "integer",
                          "number"
                        ]
                      },
                      "an": {
                        "type": [
                          "integer",
                          "number"
                        ]
                      },
                      "af": {
                        "type": "number"
                      },
                      "homozygote_count": {
                        "type": [
                          "integer",
                          "number",
                          "null"
                        ]
                      },
                      "major_consequence": {
                        "type": [
                          "string",
                          "null"
                        ]
                      },
                      "filters": {
                        "type": "array",
                        "items": {
                          "type": "string"
                        }
                      }
                    }
                  }
                }
              }
            },
            "dataset": {
              "type": "string",
              "description": "The gnomAD dataset that actually answered, including when it came from this tool's default."
            },
            "reference_genome": {
              "type": "string",
              "description": "Reference assembly the dataset is built on (gnomad_r2_1*/exac = GRCh37; gnomad_r3*/gnomad_r4* = GRCh38)."
            },
            "dataset_note": {
              "type": "string",
              "description": "Present only when the dataset structurally lacks a callset the payload reports as null (e.g. gnomad_r3 has no exomes), so `null` is not misread as 'not observed'."
            }
          },
          "required": [
            "region"
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

??? info "DNAKit 参数扩展规则"

    ```json
    {
      "aliases": {
        "end": "stop"
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `gnomad_get_sv_by_region`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gnomad_get_sv_by_region", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gnomad_get_sv_detail` {#dnakit.databases.gnomad_get_sv_detail}

Get detailed information for a specific gnomAD structural variant by its ID. Resolves IDs from both SV callsets, selected automatically from the ID format: GRCh38 IDs such as 'DEL_chr17_24e4872b' resolve against the gnomAD v4 SV callset (gnomad_sv_r4), while GRCh37 IDs such as 'DEL_1_1937' resolve against the gnomAD v2.1 SV callset (gnomad_sv_r2_1). Returns SV coordinates, type, length, population allele frequency, consequence, and filter status (consequence and homozygote/hemizygote counts are null for multi-allelic CNVs and intergenic SVs). Use gnomad_get_sv_by_gene or gnomad_get_sv_by_region to discover variant IDs first.

```text
dnakit.databases.gnomad_get_sv_detail(*, variant_id=None, dataset=None, sv_id=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant_id` | string | None | gnomAD structural variant ID. GRCh38/v4 IDs look like &#x27;DEL_chr17_24e4872b&#x27;; GRCh37/v2.1 IDs look like &#x27;DEL_1_1937&#x27;. |
| `dataset` | string | None | Optional SV callset override. Defaults to the callset auto-detected from the variant ID format (gnomad_sv_r2_1 for GRCh37-style IDs, gnomad_sv_r4 otherwise). Set explicitly only to query a v2.1 subset such as gnomad_sv_r2_1_controls. {&quot;enum&quot;: &#91;&quot;gnomad_sv_r2_1&quot;, &quot;gnomad_sv_r2_1_controls&quot;, &quot;gnomad_sv_r2_1_non_neuro&quot;, &quot;gnomad_sv_r4&quot;&#93;} |
| `sv_id` | string | None | DNAKit optional alias for variant_id |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "variant_id": {
          "type": "string",
          "description": "gnomAD structural variant ID. GRCh38/v4 IDs look like 'DEL_chr17_24e4872b'; GRCh37/v2.1 IDs look like 'DEL_1_1937'."
        },
        "dataset": {
          "type": "string",
          "enum": [
            "gnomad_sv_r2_1",
            "gnomad_sv_r2_1_controls",
            "gnomad_sv_r2_1_non_neuro",
            "gnomad_sv_r4"
          ],
          "description": "Optional SV callset override. Defaults to the callset auto-detected from the variant ID format (gnomad_sv_r2_1 for GRCh37-style IDs, gnomad_sv_r4 otherwise). Set explicitly only to query a v2.1 subset such as gnomad_sv_r2_1_controls."
        },
        "sv_id": {
          "type": "string",
          "description": "DNAKit optional alias for variant_id"
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "structural_variant": {
              "type": "object",
              "properties": {
                "variant_id": {
                  "type": "string"
                },
                "chrom": {
                  "type": "string"
                },
                "pos": {
                  "type": [
                    "integer",
                    "number"
                  ]
                },
                "end": {
                  "type": [
                    "integer",
                    "number"
                  ]
                },
                "type": {
                  "type": "string"
                },
                "length": {
                  "type": [
                    "integer",
                    "number"
                  ]
                },
                "ac": {
                  "type": [
                    "integer",
                    "number"
                  ]
                },
                "an": {
                  "type": [
                    "integer",
                    "number"
                  ]
                },
                "af": {
                  "type": "number"
                },
                "homozygote_count": {
                  "type": [
                    "integer",
                    "number",
                    "null"
                  ]
                },
                "hemizygote_count": {
                  "type": [
                    "integer",
                    "number",
                    "null"
                  ]
                },
                "major_consequence": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "consequence": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "reference_genome": {
                  "type": "string"
                },
                "filters": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  }
                },
                "chrom2": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "pos2": {
                  "type": [
                    "integer",
                    "number",
                    "null"
                  ]
                },
                "end2": {
                  "type": [
                    "integer",
                    "number",
                    "null"
                  ]
                },
                "ac_hom": {
                  "type": [
                    "integer",
                    "number",
                    "null"
                  ]
                },
                "ac_hemi": {
                  "type": [
                    "integer",
                    "number",
                    "null"
                  ]
                }
              }
            },
            "dataset": {
              "type": "string",
              "description": "The gnomAD dataset that actually answered, including when it came from this tool's default."
            },
            "reference_genome": {
              "type": "string",
              "description": "Reference assembly the dataset is built on (gnomad_r2_1*/exac = GRCh37; gnomad_r3*/gnomad_r4* = GRCh38)."
            },
            "dataset_note": {
              "type": "string",
              "description": "Present only when the dataset structurally lacks a callset the payload reports as null (e.g. gnomad_r3 has no exomes), so `null` is not misread as 'not observed'."
            }
          },
          "required": [
            "structural_variant"
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

??? info "DNAKit 参数扩展规则"

    ```json
    {
      "aliases": {
        "sv_id": "variant_id"
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `gnomad_get_sv_detail`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gnomad_get_sv_detail", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gnomad_get_variant` {#dnakit.databases.gnomad_get_variant}

Get basic variant metadata from gnomAD by `variant_id` (format like '19-44908822-C-T'). Use `gnomad_search_variants` to discover valid `variant_id` values from an rsID query. `dataset` is required by gnomAD and defaults to `gnomad_r3`. `gnomad_r3` is a GENOMES-ONLY callset with no exome component, so under that default `exome` is null for every variant -- that means 'not present in this callset', NOT 'not observed in exomes'; responses carry a `dataset_note` saying so. Pass dataset='gnomad_r4' for exome allele counts (on GJB2 c.35delG that is a further 10,423 allele observations out of 1,461,600 alleles). Defaults deliberately differ across this gnomAD family: gnomad_get_variant, gnomad_search_variants and gnomad_get_region default to `gnomad_r3` while gnomad_get_variant_populations and gnomad_get_constraint default to `gnomad_r4`, so gnomad_get_variant and gnomad_get_variant_populations report different numbers for the same variant unless you pass the same `dataset` to both. Every response echoes `dataset` and `reference_genome` (gnomad_r2_1* and exac = GRCh37; gnomad_r3* and gnomad_r4* = GRCh38) so the callset that answered is never a guess.

```text
dnakit.databases.gnomad_get_variant(*, variant_id, dataset=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant_id` | string | 必填（输入 Schema） | Variant ID (e.g., &#x27;19-44908822-C-T&#x27;). |
| `dataset` | string | None | gnomAD dataset ID. Allowed values: gnomad_r4, gnomad_r4_non_ukb, gnomad_r3, gnomad_r3_controls_and_biobanks, gnomad_r3_non_cancer, gnomad_r3_non_neuro, gnomad_r3_non_topmed, gnomad_r3_non_v2, gnomad_r2_1, gnomad_r2_1_controls, gnomad_r2_1_non_neuro, gnomad_r2_1_non_cancer, gnomad_r2_1_non_topmed, exac. Coverage: gnomad_r4* = GRCh38 with genomes + exomes; gnomad_r3* = GRCh38 with genomes ONLY (no exome component, so `exome` is always null); gnomad_r2_1* = GRCh37 with genomes + exomes; exac = GRCh37 with exomes only. {&quot;default&quot;: &quot;gnomad_r3&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "variant_id": {
          "type": "string",
          "description": "Variant ID (e.g., '19-44908822-C-T')."
        },
        "dataset": {
          "type": "string",
          "default": "gnomad_r3",
          "description": "gnomAD dataset ID. Allowed values: gnomad_r4, gnomad_r4_non_ukb, gnomad_r3, gnomad_r3_controls_and_biobanks, gnomad_r3_non_cancer, gnomad_r3_non_neuro, gnomad_r3_non_topmed, gnomad_r3_non_v2, gnomad_r2_1, gnomad_r2_1_controls, gnomad_r2_1_non_neuro, gnomad_r2_1_non_cancer, gnomad_r2_1_non_topmed, exac. Coverage: gnomad_r4* = GRCh38 with genomes + exomes; gnomad_r3* = GRCh38 with genomes ONLY (no exome component, so `exome` is always null); gnomad_r2_1* = GRCh37 with genomes + exomes; exac = GRCh37 with exomes only."
        }
      },
      "required": [
        "variant_id"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "variant": {
              "type": "object",
              "properties": {
                "variant_id": {
                  "type": "string"
                },
                "chrom": {
                  "type": "string"
                },
                "pos": {
                  "type": [
                    "integer",
                    "number"
                  ]
                },
                "ref": {
                  "type": "string"
                },
                "alt": {
                  "type": "string"
                },
                "rsid": {
                  "type": "string"
                },
                "genome": {
                  "type": "object",
                  "properties": {
                    "af": {
                      "type": "number"
                    },
                    "an": {
                      "type": [
                        "integer",
                        "number"
                      ]
                    },
                    "ac": {
                      "type": [
                        "integer",
                        "number"
                      ]
                    }
                  }
                },
                "exome": {
                  "type": [
                    "object",
                    "null"
                  ],
                  "description": "Exome callset frequencies. Null when the variant is absent from the exome callset OR when the selected dataset has no exome component at all (gnomad_r3*) -- see `dataset_note`.",
                  "properties": {
                    "af": {
                      "type": "number"
                    },
                    "an": {
                      "type": [
                        "integer",
                        "number"
                      ]
                    },
                    "ac": {
                      "type": [
                        "integer",
                        "number"
                      ]
                    }
                  }
                }
              }
            },
            "dataset": {
              "type": "string",
              "description": "The gnomAD dataset that actually answered, including when it came from this tool's default."
            },
            "reference_genome": {
              "type": "string",
              "description": "Reference assembly the dataset is built on (gnomad_r2_1*/exac = GRCh37; gnomad_r3*/gnomad_r4* = GRCh38)."
            },
            "dataset_note": {
              "type": "string",
              "description": "Present only when the dataset structurally lacks a callset the payload reports as null (e.g. gnomad_r3 has no exomes), so `null` is not misread as 'not observed'."
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

**执行与来源：** `network`；上游接口名 `gnomad_get_variant`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gnomad_get_variant", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gnomad_get_variant_populations` {#dnakit.databases.gnomad_get_variant_populations}

Get per-ancestry (population-stratified) allele frequencies for a variant from gnomAD by `variant_id` (format like '1-55051215-G-GA'). Unlike `gnomad_get_variant` (aggregate AF only), this returns per-ancestry rows (`amr`, `nfe`, `fin`, `afr`, `eas`, `sas`, `asj`, `ami`, `mid`, `remaining`, plus `_XX`/`_XY` sex splits) with computed `af = ac / an` (af is null when an = 0), separated into `genome` and `exome` callsets. `dataset` defaults to `gnomad_r4`. Use `gnomad_search_variants` to find a `variant_id` from an rsID. `gnomad_r4` carries both genomes and exomes; `gnomad_r3` is genomes-only, so selecting it makes `exome` null for every variant (a property of the callset, not of the variant) and the response says so in `dataset_note`. Defaults deliberately differ across this gnomAD family: gnomad_get_variant, gnomad_search_variants and gnomad_get_region default to `gnomad_r3` while gnomad_get_variant_populations and gnomad_get_constraint default to `gnomad_r4`, so gnomad_get_variant and gnomad_get_variant_populations report different numbers for the same variant unless you pass the same `dataset` to both. Every response echoes `dataset` and `reference_genome` (gnomad_r2_1* and exac = GRCh37; gnomad_r3* and gnomad_r4* = GRCh38) so the callset that answered is never a guess.

```text
dnakit.databases.gnomad_get_variant_populations(*, variant_id, dataset=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant_id` | string | 必填（输入 Schema） | Variant ID (format &#x27;chrom-pos-ref-alt&#x27;, e.g., &#x27;1-55051215-G-GA&#x27;). |
| `dataset` | string | None | gnomAD dataset ID. Allowed values: gnomad_r4, gnomad_r4_non_ukb, gnomad_r3, gnomad_r3_controls_and_biobanks, gnomad_r3_non_cancer, gnomad_r3_non_neuro, gnomad_r3_non_topmed, gnomad_r3_non_v2, gnomad_r2_1, gnomad_r2_1_controls, gnomad_r2_1_non_neuro, gnomad_r2_1_non_cancer, gnomad_r2_1_non_topmed, exac. Coverage: gnomad_r4* = GRCh38 with genomes + exomes; gnomad_r3* = GRCh38 with genomes ONLY (no exome component, so `exome` is always null); gnomad_r2_1* = GRCh37 with genomes + exomes; exac = GRCh37 with exomes only. {&quot;default&quot;: &quot;gnomad_r4&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "variant_id": {
          "type": "string",
          "description": "Variant ID (format 'chrom-pos-ref-alt', e.g., '1-55051215-G-GA')."
        },
        "dataset": {
          "type": "string",
          "default": "gnomad_r4",
          "description": "gnomAD dataset ID. Allowed values: gnomad_r4, gnomad_r4_non_ukb, gnomad_r3, gnomad_r3_controls_and_biobanks, gnomad_r3_non_cancer, gnomad_r3_non_neuro, gnomad_r3_non_topmed, gnomad_r3_non_v2, gnomad_r2_1, gnomad_r2_1_controls, gnomad_r2_1_non_neuro, gnomad_r2_1_non_cancer, gnomad_r2_1_non_topmed, exac. Coverage: gnomad_r4* = GRCh38 with genomes + exomes; gnomad_r3* = GRCh38 with genomes ONLY (no exome component, so `exome` is always null); gnomad_r2_1* = GRCh37 with genomes + exomes; exac = GRCh37 with exomes only."
        }
      },
      "required": [
        "variant_id"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "variant_id": {
              "type": "string"
            },
            "chrom": {
              "type": "string"
            },
            "pos": {
              "type": [
                "integer",
                "number"
              ]
            },
            "ref": {
              "type": "string"
            },
            "alt": {
              "type": "string"
            },
            "rsid": {
              "type": [
                "string",
                "null"
              ]
            },
            "dataset": {
              "type": "string",
              "description": "The gnomAD dataset that actually answered, including when it came from this tool's default."
            },
            "genome": {
              "type": [
                "object",
                "null"
              ],
              "properties": {
                "ac": {
                  "type": [
                    "integer",
                    "number",
                    "null"
                  ]
                },
                "an": {
                  "type": [
                    "integer",
                    "number",
                    "null"
                  ]
                },
                "af": {
                  "type": [
                    "number",
                    "null"
                  ]
                },
                "populations": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "id": {
                        "type": "string"
                      },
                      "ac": {
                        "type": [
                          "integer",
                          "number",
                          "null"
                        ]
                      },
                      "an": {
                        "type": [
                          "integer",
                          "number",
                          "null"
                        ]
                      },
                      "af": {
                        "type": [
                          "number",
                          "null"
                        ]
                      }
                    }
                  }
                }
              }
            },
            "exome": {
              "type": [
                "object",
                "null"
              ],
              "properties": {
                "ac": {
                  "type": [
                    "integer",
                    "number",
                    "null"
                  ]
                },
                "an": {
                  "type": [
                    "integer",
                    "number",
                    "null"
                  ]
                },
                "af": {
                  "type": [
                    "number",
                    "null"
                  ]
                },
                "populations": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "id": {
                        "type": "string"
                      },
                      "ac": {
                        "type": [
                          "integer",
                          "number",
                          "null"
                        ]
                      },
                      "an": {
                        "type": [
                          "integer",
                          "number",
                          "null"
                        ]
                      },
                      "af": {
                        "type": [
                          "number",
                          "null"
                        ]
                      }
                    }
                  }
                }
              }
            },
            "reference_genome": {
              "type": "string",
              "description": "Reference assembly the dataset is built on (gnomad_r2_1*/exac = GRCh37; gnomad_r3*/gnomad_r4* = GRCh38)."
            },
            "dataset_note": {
              "type": "string",
              "description": "Present only when the dataset structurally lacks a callset the payload reports as null (e.g. gnomad_r3 has no exomes), so `null` is not misread as 'not observed'."
            }
          },
          "required": [
            "variant_id"
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

**执行与来源：** `network`；上游接口名 `gnomad_get_variant_populations`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gnomad_get_variant_populations", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gnomad_search_variants` {#dnakit.databases.gnomad_search_variants}

Search for variants in gnomAD by free-text query (commonly an rsID like 'rs7412'). Returns matching `variant_id` strings. Note: gnomAD's `variant_search` response only exposes `variant_id`. Use the returned `variant_id` with `gnomad_get_variant`. `dataset` defaults to `gnomad_r3` (GRCh38, genomes only); which variants exist depends on the callset searched. Defaults deliberately differ across this gnomAD family: gnomad_get_variant, gnomad_search_variants and gnomad_get_region default to `gnomad_r3` while gnomad_get_variant_populations and gnomad_get_constraint default to `gnomad_r4`, so gnomad_get_variant and gnomad_get_variant_populations report different numbers for the same variant unless you pass the same `dataset` to both. Every response echoes `dataset` and `reference_genome` (gnomad_r2_1* and exac = GRCh37; gnomad_r3* and gnomad_r4* = GRCh38) so the callset that answered is never a guess.

```text
dnakit.databases.gnomad_search_variants(*, query=None, dataset=None, variant=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | None | Variant search query (e.g., &#x27;rs7412&#x27;). |
| `dataset` | string | None | gnomAD dataset ID. Allowed values: gnomad_r4, gnomad_r4_non_ukb, gnomad_r3, gnomad_r3_controls_and_biobanks, gnomad_r3_non_cancer, gnomad_r3_non_neuro, gnomad_r3_non_topmed, gnomad_r3_non_v2, gnomad_r2_1, gnomad_r2_1_controls, gnomad_r2_1_non_neuro, gnomad_r2_1_non_cancer, gnomad_r2_1_non_topmed, exac. Coverage: gnomad_r4* = GRCh38 with genomes + exomes; gnomad_r3* = GRCh38 with genomes ONLY (no exome component, so `exome` is always null); gnomad_r2_1* = GRCh37 with genomes + exomes; exac = GRCh37 with exomes only. {&quot;default&quot;: &quot;gnomad_r3&quot;} |
| `variant` | string | None | DNAKit optional alias for query |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "query": {
          "type": "string",
          "description": "Variant search query (e.g., 'rs7412')."
        },
        "dataset": {
          "type": "string",
          "default": "gnomad_r3",
          "description": "gnomAD dataset ID. Allowed values: gnomad_r4, gnomad_r4_non_ukb, gnomad_r3, gnomad_r3_controls_and_biobanks, gnomad_r3_non_cancer, gnomad_r3_non_neuro, gnomad_r3_non_topmed, gnomad_r3_non_v2, gnomad_r2_1, gnomad_r2_1_controls, gnomad_r2_1_non_neuro, gnomad_r2_1_non_cancer, gnomad_r2_1_non_topmed, exac. Coverage: gnomad_r4* = GRCh38 with genomes + exomes; gnomad_r3* = GRCh38 with genomes ONLY (no exome component, so `exome` is always null); gnomad_r2_1* = GRCh37 with genomes + exomes; exac = GRCh37 with exomes only."
        },
        "variant": {
          "type": "string",
          "description": "DNAKit optional alias for query"
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "variant_search": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "variant_id": {
                    "type": "string"
                  }
                }
              }
            },
            "dataset": {
              "type": "string",
              "description": "The gnomAD dataset that actually answered, including when it came from this tool's default."
            },
            "reference_genome": {
              "type": "string",
              "description": "Reference assembly the dataset is built on (gnomad_r2_1*/exac = GRCh37; gnomad_r3*/gnomad_r4* = GRCh38)."
            },
            "dataset_note": {
              "type": "string",
              "description": "Present only when the dataset structurally lacks a callset the payload reports as null (e.g. gnomad_r3 has no exomes), so `null` is not misread as 'not observed'."
            }
          },
          "required": [
            "variant_search"
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

??? info "DNAKit 参数扩展规则"

    ```json
    {
      "aliases": {
        "variant": "query"
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `gnomad_search_variants`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gnomad_search_variants", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gwas_get_associations_for_snp` {#dnakit.databases.gwas_get_associations_for_snp}

Get all associations for a specific SNP with optional sorting.

```text
dnakit.databases.gwas_get_associations_for_snp(*, rs_id, sort=None, direction=None, size=None, page=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `rs_id` | string | 必填（输入 Schema） | dbSNP rs identifier |
| `sort` | string | None | Sort field (e.g., &#x27;p_value&#x27;, &#x27;or_value&#x27;) |
| `direction` | string | None | Sort direction (&#x27;asc&#x27; or &#x27;desc&#x27;) |
| `size` | integer | None | Number of results to return per page |
| `page` | integer | None | Page number for pagination |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "rs_id": {
          "type": "string",
          "description": "dbSNP rs identifier"
        },
        "sort": {
          "type": "string",
          "description": "Sort field (e.g., 'p_value', 'or_value')"
        },
        "direction": {
          "type": "string",
          "description": "Sort direction ('asc' or 'desc')"
        },
        "size": {
          "type": "integer",
          "description": "Number of results to return per page"
        },
        "page": {
          "type": "integer",
          "description": "Page number for pagination"
        }
      },
      "required": [
        "rs_id"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "array",
          "description": "Array of GWAS association objects",
          "items": {
            "type": "object",
            "properties": {
              "association_id": {
                "type": "integer",
                "description": "Unique association identifier"
              },
              "p_value": {
                "type": "number",
                "description": "Statistical p-value"
              },
              "beta": {
                "type": "string",
                "description": "Effect size (beta coefficient)"
              },
              "efo_traits": {
                "type": "array",
                "items": {
                  "type": "object",
                  "properties": {
                    "efo_id": {
                      "type": "string"
                    },
                    "efo_trait": {
                      "type": "string"
                    }
                  }
                }
              },
              "reported_trait": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "accession_id": {
                "type": "string",
                "description": "Study accession ID"
              },
              "pubmed_id": {
                "type": [
                  "string",
                  "integer"
                ],
                "description": "PubMed reference (can be string or integer)"
              },
              "first_author": {
                "type": "string",
                "description": "First author name"
              },
              "snp_effect_allele": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "snp_allele": {
                "type": "array",
                "items": {
                  "type": "object",
                  "properties": {
                    "rs_id": {
                      "type": "string"
                    },
                    "effect_allele": {
                      "type": "string"
                    }
                  }
                }
              },
              "locations": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "mapped_genes": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "ci_lower": {
                "type": "number",
                "description": "Lower confidence interval"
              },
              "ci_upper": {
                "type": "number",
                "description": "Upper confidence interval"
              },
              "risk_frequency": {
                "type": "string",
                "description": "Risk allele frequency"
              },
              "range": {
                "type": "string",
                "description": "Effect size range"
              },
              "_links": {
                "type": "object",
                "description": "API links"
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

**执行与来源：** `network`；上游接口名 `gwas_get_associations_for_snp`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gwas_get_associations_for_snp", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gwas_get_associations_for_study` {#dnakit.databases.gwas_get_associations_for_study}

Get all associations for a specific study, sorted by p-value.

```text
dnakit.databases.gwas_get_associations_for_study(*, accession_id, size=None, page=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession_id` | string | 必填（输入 Schema） | Study accession identifier |
| `size` | integer | None | Number of results to return per page |
| `page` | integer | None | Page number for pagination |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "accession_id": {
          "type": "string",
          "description": "Study accession identifier"
        },
        "size": {
          "type": "integer",
          "description": "Number of results to return per page"
        },
        "page": {
          "type": "integer",
          "description": "Page number for pagination"
        }
      },
      "required": [
        "accession_id"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "array",
          "description": "Array of GWAS association objects",
          "items": {
            "type": "object",
            "properties": {
              "association_id": {
                "type": "integer",
                "description": "Unique association identifier"
              },
              "p_value": {
                "type": "number",
                "description": "Statistical p-value"
              },
              "beta": {
                "type": "string",
                "description": "Effect size (beta coefficient)"
              },
              "efo_traits": {
                "type": "array",
                "items": {
                  "type": "object",
                  "properties": {
                    "efo_id": {
                      "type": "string"
                    },
                    "efo_trait": {
                      "type": "string"
                    }
                  }
                }
              },
              "reported_trait": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "accession_id": {
                "type": "string",
                "description": "Study accession ID"
              },
              "pubmed_id": {
                "type": [
                  "string",
                  "integer"
                ],
                "description": "PubMed reference (can be string or integer)"
              },
              "first_author": {
                "type": "string",
                "description": "First author name"
              },
              "snp_effect_allele": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "snp_allele": {
                "type": "array",
                "items": {
                  "type": "object",
                  "properties": {
                    "rs_id": {
                      "type": "string"
                    },
                    "effect_allele": {
                      "type": "string"
                    }
                  }
                }
              },
              "locations": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "mapped_genes": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "ci_lower": {
                "type": "number",
                "description": "Lower confidence interval"
              },
              "ci_upper": {
                "type": "number",
                "description": "Upper confidence interval"
              },
              "risk_frequency": {
                "type": "string",
                "description": "Risk allele frequency"
              },
              "range": {
                "type": "string",
                "description": "Effect size range"
              },
              "_links": {
                "type": "object",
                "description": "API links"
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

**执行与来源：** `network`；上游接口名 `gwas_get_associations_for_study`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gwas_get_associations_for_study", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gwas_get_associations_for_trait` {#dnakit.databases.gwas_get_associations_for_trait}

Get all associations for a specific trait, sorted by p-value (most significant first). A free-text disease_trait is resolved to an EFO/MONDO term first, and that resolution is disclosed in the response: 'query_trait' (what you asked for), 'resolved_efo_id', 'resolved_efo_label' (the term's label in GWAS Catalog) and 'trait_resolution_source', plus a 'trait_resolution_note' whenever the resolved label is not an exact match for your query. Check those before trusting the results: e.g. 'cardiorespiratory fitness' resolves to EFO_0009184 'heart rate response to exercise', which is a different phenotype from VO2max (EFO_0004887 'maximal oxygen uptake measurement'). Pass efo_id explicitly to bypass the resolver.

```text
dnakit.databases.gwas_get_associations_for_trait(*, disease_trait=None, efo_uri=None, efo_id=None, efo_trait=None, size=None, page=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `disease_trait` | string | None | Disease or trait name for text-based search (e.g., &#x27;diabetes&#x27;, &#x27;hypertension&#x27;) |
| `efo_uri` | string | None | Full EFO ontology URI (e.g., &#x27;http://www.ebi.ac.uk/efo/EFO_0001645&#x27;) |
| `efo_id` | string | None | EFO/OBA term ID (e.g., &#x27;EFO_0001645&#x27;, &#x27;OBA_2050062&#x27;). Recommended for reliable trait filtering. |
| `efo_trait` | string | None | Exact EFO trait label. Use when you know the canonical trait string. |
| `size` | integer | None | Number of results to return per page |
| `page` | integer | None | Page number for pagination |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "disease_trait": {
          "type": "string",
          "description": "Disease or trait name for text-based search (e.g., 'diabetes', 'hypertension')"
        },
        "efo_uri": {
          "type": "string",
          "description": "Full EFO ontology URI (e.g., 'http://www.ebi.ac.uk/efo/EFO_0001645')"
        },
        "efo_id": {
          "type": "string",
          "description": "EFO/OBA term ID (e.g., 'EFO_0001645', 'OBA_2050062'). Recommended for reliable trait filtering."
        },
        "efo_trait": {
          "type": "string",
          "description": "Exact EFO trait label. Use when you know the canonical trait string."
        },
        "size": {
          "type": "integer",
          "description": "Number of results to return per page"
        },
        "page": {
          "type": "integer",
          "description": "Page number for pagination"
        }
      },
      "oneOf": [
        {
          "required": [
            "disease_trait"
          ]
        },
        {
          "required": [
            "efo_uri"
          ]
        },
        {
          "required": [
            "efo_id"
          ]
        },
        {
          "required": [
            "efo_trait"
          ]
        }
      ],
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "array",
          "description": "Array of GWAS association objects",
          "items": {
            "type": "object",
            "properties": {
              "association_id": {
                "type": "integer",
                "description": "Unique association identifier"
              },
              "p_value": {
                "type": "number",
                "description": "Statistical p-value"
              },
              "beta": {
                "type": "string",
                "description": "Effect size (beta coefficient)"
              },
              "efo_traits": {
                "type": "array",
                "items": {
                  "type": "object",
                  "properties": {
                    "efo_id": {
                      "type": "string"
                    },
                    "efo_trait": {
                      "type": "string"
                    }
                  }
                }
              },
              "reported_trait": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "accession_id": {
                "type": "string",
                "description": "Study accession ID"
              },
              "pubmed_id": {
                "type": [
                  "string",
                  "integer"
                ],
                "description": "PubMed reference (can be string or integer)"
              },
              "first_author": {
                "type": "string",
                "description": "First author name"
              },
              "snp_effect_allele": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "snp_allele": {
                "type": "array",
                "items": {
                  "type": "object",
                  "properties": {
                    "rs_id": {
                      "type": "string"
                    },
                    "effect_allele": {
                      "type": "string"
                    }
                  }
                }
              },
              "locations": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "mapped_genes": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "ci_lower": {
                "type": "number",
                "description": "Lower confidence interval"
              },
              "ci_upper": {
                "type": "number",
                "description": "Upper confidence interval"
              },
              "risk_frequency": {
                "type": "string",
                "description": "Risk allele frequency"
              },
              "range": {
                "type": "string",
                "description": "Effect size range"
              },
              "_links": {
                "type": "object",
                "description": "API links"
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

**执行与来源：** `network`；上游接口名 `gwas_get_associations_for_trait`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gwas_get_associations_for_trait", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gwas_get_snp_by_id` {#dnakit.databases.gwas_get_snp_by_id}

Get a specific GWAS SNP by its rs ID.

```text
dnakit.databases.gwas_get_snp_by_id(*, rs_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `rs_id` | string | 必填（输入 Schema） | dbSNP rs identifier |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "rs_id": {
          "type": "string",
          "description": "dbSNP rs identifier"
        }
      },
      "required": [
        "rs_id"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "description": "Single SNP object",
      "properties": {
        "rs_id": {
          "type": "string",
          "description": "dbSNP rs identifier"
        },
        "merged": {
          "type": "integer",
          "description": "Merged status"
        },
        "functional_class": {
          "type": "string",
          "description": "Functional classification"
        },
        "last_update_date": {
          "type": "string",
          "description": "Last update timestamp"
        },
        "locations": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "chromosome_name": {
                "type": "string"
              },
              "chromosome_position": {
                "type": "integer"
              },
              "region": {
                "type": "object",
                "properties": {
                  "name": {
                    "type": "string"
                  }
                }
              }
            }
          }
        },
        "alleles": {
          "type": "string",
          "description": "Allele information"
        },
        "most_severe_consequence": {
          "type": "string",
          "description": "Most severe functional consequence"
        },
        "mapped_genes": {
          "type": "array",
          "items": {
            "type": "string"
          }
        },
        "maf": {
          "type": "number",
          "description": "Minor allele frequency"
        },
        "minor_allele": {
          "type": "string",
          "description": "Minor allele"
        },
        "_links": {
          "type": "object",
          "description": "API links"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `gwas_get_snp_by_id`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gwas_get_snp_by_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gwas_get_snps_for_gene` {#dnakit.databases.gwas_get_snps_for_gene}

Get all GWAS Catalog SNPs mapped to a specific gene by gene name (e.g., 'BRCA1', 'TP53'). Returns variant details including rs IDs, chromosomal positions, and functional annotations. Use this when you want to find all catalogued variants within or near a specific gene.

```text
dnakit.databases.gwas_get_snps_for_gene(*, size=None, page=None, gene_symbol=None, gene=None, mapped_gene=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `size` | integer | None | Number of results to return per page |
| `page` | integer | None | Page number for pagination |
| `gene_symbol` | string | None | Gene name or symbol (e.g., &#x27;BRCA1&#x27;, &#x27;TP53&#x27;, &#x27;GBA&#x27;) |
| `gene` | string | None | Alias for gene_symbol |
| `mapped_gene` | string | None | Alias for gene_symbol |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "size": {
          "type": "integer",
          "description": "Number of results to return per page"
        },
        "page": {
          "type": "integer",
          "description": "Page number for pagination"
        },
        "gene_symbol": {
          "type": "string",
          "description": "Gene name or symbol (e.g., 'BRCA1', 'TP53', 'GBA')"
        },
        "gene": {
          "type": "string",
          "description": "Alias for gene_symbol"
        },
        "mapped_gene": {
          "type": "string",
          "description": "Alias for gene_symbol"
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "array",
          "description": "Array of SNP objects",
          "items": {
            "type": "object",
            "properties": {
              "rs_id": {
                "type": "string",
                "description": "dbSNP rs identifier"
              },
              "merged": {
                "type": "integer",
                "description": "Merged status"
              },
              "functional_class": {
                "type": "string",
                "description": "Functional classification"
              },
              "last_update_date": {
                "type": "string",
                "description": "Last update timestamp"
              },
              "locations": {
                "type": "array",
                "items": {
                  "type": "object",
                  "properties": {
                    "chromosome_name": {
                      "type": "string"
                    },
                    "chromosome_position": {
                      "type": "integer"
                    },
                    "region": {
                      "type": "object",
                      "properties": {
                        "name": {
                          "type": "string"
                        }
                      }
                    }
                  }
                }
              },
              "alleles": {
                "type": "string",
                "description": "Allele information"
              },
              "most_severe_consequence": {
                "type": "string",
                "description": "Most severe functional consequence"
              },
              "mapped_genes": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "maf": {
                "type": "number",
                "description": "Minor allele frequency"
              },
              "minor_allele": {
                "type": "string",
                "description": "Minor allele"
              },
              "_links": {
                "type": "object",
                "description": "API links"
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

**执行与来源：** `network`；上游接口名 `gwas_get_snps_for_gene`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gwas_get_snps_for_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gwas_get_studies_for_trait` {#dnakit.databases.gwas_get_studies_for_trait}

Get studies for a specific trait with optional filters for cohort, GxE interactions, and summary statistics.

```text
dnakit.databases.gwas_get_studies_for_trait(*, disease_trait=None, efo_uri=None, efo_id=None, efo_trait=None, cohort=None, gxe=None, full_pvalue_set=None, size=None, page=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `disease_trait` | string | None | Disease or trait name for text-based search (e.g., &#x27;diabetes&#x27;, &#x27;alzheimer disease&#x27;) |
| `efo_uri` | string | None | Full EFO ontology URI (e.g., &#x27;http://www.ebi.ac.uk/efo/EFO_0001645&#x27;) |
| `efo_id` | string | None | EFO/OBA term ID (e.g., &#x27;EFO_0001645&#x27;, &#x27;OBA_2050062&#x27;). Recommended for reliable trait filtering. |
| `efo_trait` | string | None | Exact EFO trait label. Use when you know the canonical trait string. |
| `cohort` | string | None | Cohort name (e.g., &#x27;UKB&#x27; for UK Biobank) |
| `gxe` | boolean | None | Filter for Gene-by-Environment interaction studies |
| `full_pvalue_set` | boolean | None | Filter for studies with full summary statistics |
| `size` | integer | None | Number of results to return per page |
| `page` | integer | None | Page number for pagination |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "disease_trait": {
          "type": "string",
          "description": "Disease or trait name for text-based search (e.g., 'diabetes', 'alzheimer disease')"
        },
        "efo_uri": {
          "type": "string",
          "description": "Full EFO ontology URI (e.g., 'http://www.ebi.ac.uk/efo/EFO_0001645')"
        },
        "efo_id": {
          "type": "string",
          "description": "EFO/OBA term ID (e.g., 'EFO_0001645', 'OBA_2050062'). Recommended for reliable trait filtering."
        },
        "efo_trait": {
          "type": "string",
          "description": "Exact EFO trait label. Use when you know the canonical trait string."
        },
        "cohort": {
          "type": "string",
          "description": "Cohort name (e.g., 'UKB' for UK Biobank)"
        },
        "gxe": {
          "type": "boolean",
          "description": "Filter for Gene-by-Environment interaction studies"
        },
        "full_pvalue_set": {
          "type": "boolean",
          "description": "Filter for studies with full summary statistics"
        },
        "size": {
          "type": "integer",
          "description": "Number of results to return per page"
        },
        "page": {
          "type": "integer",
          "description": "Page number for pagination"
        }
      },
      "oneOf": [
        {
          "required": [
            "disease_trait"
          ]
        },
        {
          "required": [
            "efo_uri"
          ]
        },
        {
          "required": [
            "efo_id"
          ]
        },
        {
          "required": [
            "efo_trait"
          ]
        }
      ],
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "array",
          "description": "Array of GWAS study objects",
          "items": {
            "type": "object",
            "properties": {
              "accession_id": {
                "type": "string",
                "description": "Study accession identifier"
              },
              "disease_trait": {
                "type": "string",
                "description": "Disease or trait name"
              },
              "efo_traits": {
                "type": "array",
                "items": {
                  "type": "object",
                  "properties": {
                    "efo_id": {
                      "type": "string"
                    },
                    "efo_trait": {
                      "type": "string"
                    }
                  }
                }
              },
              "initial_sample_size": {
                "type": "string",
                "description": "Initial sample size description"
              },
              "replication_sample_size": {
                "type": "string",
                "description": "Replication sample size"
              },
              "discovery_ancestry": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "replication_ancestry": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "cohort": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "genotyping_technologies": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "platforms": {
                "type": "string",
                "description": "Genotyping platforms used"
              },
              "snp_count": {
                "type": "integer",
                "description": "Number of SNPs analyzed"
              },
              "imputed": {
                "type": "boolean",
                "description": "Whether data was imputed"
              },
              "pooled": {
                "type": "boolean",
                "description": "Whether data was pooled"
              },
              "gxe": {
                "type": "boolean",
                "description": "Gene-environment interaction study"
              },
              "gxg": {
                "type": "boolean",
                "description": "Gene-gene interaction study"
              },
              "full_summary_stats_available": {
                "type": "boolean",
                "description": "Full summary statistics available"
              },
              "full_summary_stats": {
                "type": "string",
                "description": "URL to full summary statistics"
              },
              "pubmed_id": {
                "type": [
                  "string",
                  "integer"
                ],
                "description": "PubMed reference (can be string or integer)"
              },
              "terms_of_license": {
                "type": "string",
                "description": "Data usage license"
              },
              "array_manufacturer": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "_links": {
                "type": "object",
                "description": "API links"
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

**执行与来源：** `network`；上游接口名 `gwas_get_studies_for_trait`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gwas_get_studies_for_trait", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gwas_get_study_by_id` {#dnakit.databases.gwas_get_study_by_id}

Get a specific GWAS study by its unique identifier. Accepts either 'study_id' or 'accession_id' (the same GCST identifier field name used by GWASAssociationsForStudy and returned on every association record) so IDs can be chained directly between gwas_* tools without renaming.

```text
dnakit.databases.gwas_get_study_by_id(*, study_id=None, accession_id=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `study_id` | string | None | GWAS study identifier (GCST accession, e.g. &#x27;GCST000392&#x27;) |
| `accession_id` | string | None | Alias for study_id -- the same GCST accession identifier, accepted so results from other gwas_* tools (which use the field name &#x27;accession_id&#x27;) can be passed straight in. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "study_id": {
          "type": "string",
          "description": "GWAS study identifier (GCST accession, e.g. 'GCST000392')"
        },
        "accession_id": {
          "type": "string",
          "description": "Alias for study_id -- the same GCST accession identifier, accepted so results from other gwas_* tools (which use the field name 'accession_id') can be passed straight in."
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "description": "Single GWAS study object",
      "properties": {
        "accession_id": {
          "type": "string",
          "description": "Study accession identifier"
        },
        "disease_trait": {
          "type": "string",
          "description": "Disease or trait name"
        },
        "efo_traits": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "efo_id": {
                "type": "string"
              },
              "efo_trait": {
                "type": "string"
              }
            }
          }
        },
        "initial_sample_size": {
          "type": "string",
          "description": "Initial sample size description"
        },
        "replication_sample_size": {
          "type": "string",
          "description": "Replication sample size"
        },
        "discovery_ancestry": {
          "type": "array",
          "items": {
            "type": "string"
          }
        },
        "replication_ancestry": {
          "type": "array",
          "items": {
            "type": "string"
          }
        },
        "cohort": {
          "type": "array",
          "items": {
            "type": "string"
          }
        },
        "genotyping_technologies": {
          "type": "array",
          "items": {
            "type": "string"
          }
        },
        "platforms": {
          "type": "string",
          "description": "Genotyping platforms used"
        },
        "snp_count": {
          "type": "integer",
          "description": "Number of SNPs analyzed"
        },
        "imputed": {
          "type": "boolean",
          "description": "Whether data was imputed"
        },
        "pooled": {
          "type": "boolean",
          "description": "Whether data was pooled"
        },
        "gxe": {
          "type": "boolean",
          "description": "Gene-environment interaction study"
        },
        "gxg": {
          "type": "boolean",
          "description": "Gene-gene interaction study"
        },
        "full_summary_stats_available": {
          "type": "boolean",
          "description": "Full summary statistics available"
        },
        "full_summary_stats": {
          "type": "string",
          "description": "URL to full summary statistics"
        },
        "pubmed_id": {
          "type": [
            "string",
            "integer"
          ],
          "description": "PubMed reference (can be string or integer)"
        },
        "terms_of_license": {
          "type": "string",
          "description": "Data usage license"
        },
        "array_manufacturer": {
          "type": "array",
          "items": {
            "type": "string"
          }
        },
        "_links": {
          "type": "object",
          "description": "API links"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `gwas_get_study_by_id`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gwas_get_study_by_id", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gwas_get_variants_for_trait` {#dnakit.databases.gwas_get_variants_for_trait}

Search the GWAS Catalog for all genetic variants (SNPs) linked to a specific disease or trait. Accepts free-text disease/trait names (e.g., 'diabetes', 'breast cancer'), EFO ontology IDs, or EFO trait labels. Returns variant details including rsIDs, p-values, mapped genes, and genomic locations. Results are sorted by p-value ascending (most significant first) by default, so the first page holds the strongest loci; override with the 'sort' and 'direction' parameters. This is the best tool for finding all GWAS variants for a disease when you do NOT already have a specific variant ID. A free-text disease_trait is resolved to an EFO/MONDO term first, and that resolution is disclosed: 'query_trait', 'resolved_efo_id', 'resolved_efo_label' and 'trait_resolution_source' are returned, plus a 'trait_resolution_note' whenever the resolved label is not an exact match for your query - check them before trusting the results, and pass efo_id explicitly to bypass the resolver.

```text
dnakit.databases.gwas_get_variants_for_trait(*, disease_trait=None, trait=None, efo_uri=None, efo_id=None, efo_trait=None, size=None, limit=None, page=None, sort=None, direction=None, p_value_threshold=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `disease_trait` | string | None | Disease or trait name for text-based search (e.g., &#x27;diabetes&#x27;, &#x27;breast cancer&#x27;) |
| `trait` | string | None | Alias for disease_trait. Disease or trait name (e.g., &#x27;type 2 diabetes&#x27;, &#x27;breast cancer&#x27;) |
| `efo_uri` | string | None | Full EFO ontology URI (e.g., &#x27;http://www.ebi.ac.uk/efo/EFO_0001645&#x27;) |
| `efo_id` | string | None | EFO/OBA term ID (e.g., &#x27;EFO_0001645&#x27;, &#x27;OBA_2050062&#x27;). Recommended for reliable trait filtering. |
| `efo_trait` | string | None | Exact EFO trait label. Use when you know the canonical trait string. |
| `size` | integer | None | Number of results to return per page |
| `limit` | integer | None | Alias for size. Number of results to return per page |
| `page` | integer | None | Page number for pagination |
| `sort` | string | None | Sort field (e.g., &#x27;p_value&#x27;, &#x27;or_value&#x27;). Defaults to &#x27;p_value&#x27; so the most significant associations appear first. |
| `direction` | string | None | Sort direction (&#x27;asc&#x27; or &#x27;desc&#x27;). Defaults to &#x27;asc&#x27;. |
| `p_value_threshold` | number | None | 对当前返回页应用的 P 值过滤阈值；不代表已遍历全库。 {&quot;minimum&quot;: 0, &quot;maximum&quot;: 1} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "disease_trait": {
          "type": "string",
          "description": "Disease or trait name for text-based search (e.g., 'diabetes', 'breast cancer')"
        },
        "trait": {
          "type": "string",
          "description": "Alias for disease_trait. Disease or trait name (e.g., 'type 2 diabetes', 'breast cancer')"
        },
        "efo_uri": {
          "type": "string",
          "description": "Full EFO ontology URI (e.g., 'http://www.ebi.ac.uk/efo/EFO_0001645')"
        },
        "efo_id": {
          "type": "string",
          "description": "EFO/OBA term ID (e.g., 'EFO_0001645', 'OBA_2050062'). Recommended for reliable trait filtering."
        },
        "efo_trait": {
          "type": "string",
          "description": "Exact EFO trait label. Use when you know the canonical trait string."
        },
        "size": {
          "type": "integer",
          "description": "Number of results to return per page"
        },
        "limit": {
          "type": "integer",
          "description": "Alias for size. Number of results to return per page"
        },
        "page": {
          "type": "integer",
          "description": "Page number for pagination"
        },
        "sort": {
          "type": "string",
          "description": "Sort field (e.g., 'p_value', 'or_value'). Defaults to 'p_value' so the most significant associations appear first."
        },
        "direction": {
          "type": "string",
          "description": "Sort direction ('asc' or 'desc'). Defaults to 'asc'."
        },
        "p_value_threshold": {
          "type": "number",
          "minimum": 0,
          "maximum": 1
        }
      },
      "oneOf": [
        {
          "required": [
            "disease_trait"
          ]
        },
        {
          "required": [
            "trait"
          ]
        },
        {
          "required": [
            "efo_uri"
          ]
        },
        {
          "required": [
            "efo_id"
          ]
        },
        {
          "required": [
            "efo_trait"
          ]
        }
      ],
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "array",
          "description": "Array of variant objects",
          "items": {
            "type": "object",
            "properties": {
              "association_id": {
                "type": "integer"
              },
              "risk_frequency": {
                "type": "string"
              },
              "pvalue_description": {
                "type": "string"
              },
              "range": {
                "type": "string"
              },
              "locations": {
                "type": "array",
                "items": {
                  "type": [
                    "object",
                    "string"
                  ]
                },
                "description": "Array of location objects or location strings"
              },
              "rs_id": {
                "type": "string",
                "description": "dbSNP rs identifier"
              },
              "merged": {
                "type": "integer",
                "description": "Merged status"
              },
              "functional_class": {
                "type": "string",
                "description": "Functional classification"
              },
              "last_update_date": {
                "type": "string",
                "description": "Last update timestamp"
              },
              "alleles": {
                "type": "string",
                "description": "Allele information"
              },
              "most_severe_consequence": {
                "type": "string",
                "description": "Most severe functional consequence"
              },
              "mapped_genes": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "maf": {
                "type": "number",
                "description": "Minor allele frequency"
              },
              "minor_allele": {
                "type": "string",
                "description": "Minor allele"
              },
              "_links": {
                "type": "object",
                "description": "API links"
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

??? info "DNAKit 参数扩展规则"

    ```json
    {
      "parameters": {
        "p_value_threshold": {
          "type": "number",
          "minimum": 0,
          "maximum": 1
        }
      },
      "handler": "p_value_filter"
    }
    ```

**执行与来源：** `network`；上游接口名 `gwas_get_variants_for_trait`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gwas_get_variants_for_trait", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gwas_search_associations` {#dnakit.databases.gwas_search_associations}

Search the GWAS Catalog for genome-wide association study results by keyword. Accepts free-text disease/trait names (e.g., 'melanoma', 'diabetes', 'breast cancer'), SNP rs IDs, EFO ontology IDs, or study accession IDs. This is the primary tool for discovering GWAS associations when you do NOT already have a specific association ID. Returns SNPs, p-values, odds ratios, mapped genes, and study metadata with pagination support. Response envelope keys beyond 'data': 'metadata.pagination' (upstream page block - its totalElements/totalPages always describe the UNFILTERED result set); when a p_value threshold is given, 'metadata.p_value_filter', 'metadata.p_value_filter_scope' and 'metadata.p_value_filter_note' state that the threshold was applied CLIENT-SIDE (GWAS Catalog REST v2 has no server-side p-value filter), 'metadata.pagination_totals_are_prefilter' is true, 'metadata.filtered_total' / 'filtered_total_is_exact' / 'filtered_total_at_least' / 'filtered_total_note' give the post-filter denominator, and 'metadata.sort_override' records that the page was fetched p_value-ascending so the threshold could match anything. When a disease_trait is resolved to an ontology term, 'query_trait', 'resolved_efo_id', 'resolved_efo_label' and 'trait_resolution_source' are returned, plus a 'trait_resolution_note' whenever the resolved label is not an exact match for the trait you asked for. An empty 'data' comes with a 'note' that names the real cause: a client-side filter that removed every row is reported as such, and only a genuinely empty unfiltered result set suggests trying a trait synonym or a gene-based tool.

```text
dnakit.databases.gwas_search_associations(*, disease_trait=None, query=None, efo_uri=None, efo_id=None, efo_trait=None, rs_id=None, accession_id=None, p_value=None, p_value_threshold=None, sort=None, direction=None, size=None, page=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `disease_trait` | string | None | Disease or trait name for text-based search (e.g., &#x27;diabetes&#x27;, &#x27;coronary artery disease&#x27;, &#x27;breast cancer&#x27;). NOTE: GWAS associations are indexed by trait/disease, not by gene or drug name — passing a gene symbol or drug name here will fail. Use efo_id for precise filtering. |
| `query` | string | None | Alias for disease_trait. Disease or trait name (e.g., &#x27;type 2 diabetes&#x27;). NOTE: gene symbols and drug names are not valid GWAS trait queries. |
| `efo_uri` | string | None | Full EFO ontology URI (e.g., &#x27;http://www.ebi.ac.uk/efo/EFO_0001645&#x27;) |
| `efo_id` | string | None | EFO/OBA term ID (e.g., &#x27;EFO_0001645&#x27;, &#x27;OBA_2050062&#x27;). Recommended for reliable trait filtering. |
| `efo_trait` | string | None | Exact EFO trait label (e.g., &#x27;serum alanine aminotransferase amount&#x27;). Use when you know the canonical trait string. |
| `rs_id` | string | None | dbSNP rs identifier |
| `accession_id` | string | None | Study accession identifier |
| `p_value` | number | None | Maximum p-value threshold for filtering (e.g., 5e-8 for genome-wide significance). IMPORTANT: the GWAS Catalog REST API v2 has NO server-side p-value filter (its OpenAPI document lists only pubmed_id, rs_id, full_pvalue_set, accession_id, efo_trait, efo_id, show_child_trait, mapped_gene, extended_geneset, sort, direction, page and size on GET /v2/associations), so this threshold is applied CLIENT-SIDE to the page of associations fetched. To make it meaningful the page is fetched p_value-ascending (most-significant-first) whenever sort is &#x27;p_value&#x27; (the default), and the surviving rows are then re-ordered into the direction you asked for. metadata.p_value_filter_scope, metadata.filtered_total and metadata.p_value_filter_note report exactly what was filtered. |
| `p_value_threshold` | number | None | Alias for p_value: maximum p-value threshold (also applied client-side; see p_value). |
| `sort` | string | None | Sort field (e.g., &#x27;p_value&#x27;, &#x27;or_value&#x27;). Defaults to &#x27;p_value&#x27;. When a p_value threshold is supplied and sort is &#x27;p_value&#x27;, the upstream page is always fetched ascending so the threshold sees the significant end; metadata.sort_override records this and the returned rows are ordered per your &#x27;direction&#x27;. |
| `direction` | string | None | Sort direction (&#x27;asc&#x27; or &#x27;desc&#x27;). Defaults to &#x27;asc&#x27;. Combined with a p_value threshold this orders the RETURNED rows; see metadata.sort_override. |
| `size` | integer | None | Number of results to return |
| `page` | integer | None | Page number for pagination |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "disease_trait": {
          "type": "string",
          "description": "Disease or trait name for text-based search (e.g., 'diabetes', 'coronary artery disease', 'breast cancer'). NOTE: GWAS associations are indexed by trait/disease, not by gene or drug name — passing a gene symbol or drug name here will fail. Use efo_id for precise filtering."
        },
        "query": {
          "type": "string",
          "description": "Alias for disease_trait. Disease or trait name (e.g., 'type 2 diabetes'). NOTE: gene symbols and drug names are not valid GWAS trait queries."
        },
        "efo_uri": {
          "type": "string",
          "description": "Full EFO ontology URI (e.g., 'http://www.ebi.ac.uk/efo/EFO_0001645')"
        },
        "efo_id": {
          "type": "string",
          "description": "EFO/OBA term ID (e.g., 'EFO_0001645', 'OBA_2050062'). Recommended for reliable trait filtering."
        },
        "efo_trait": {
          "type": "string",
          "description": "Exact EFO trait label (e.g., 'serum alanine aminotransferase amount'). Use when you know the canonical trait string."
        },
        "rs_id": {
          "type": "string",
          "description": "dbSNP rs identifier"
        },
        "accession_id": {
          "type": "string",
          "description": "Study accession identifier"
        },
        "p_value": {
          "type": "number",
          "description": "Maximum p-value threshold for filtering (e.g., 5e-8 for genome-wide significance). IMPORTANT: the GWAS Catalog REST API v2 has NO server-side p-value filter (its OpenAPI document lists only pubmed_id, rs_id, full_pvalue_set, accession_id, efo_trait, efo_id, show_child_trait, mapped_gene, extended_geneset, sort, direction, page and size on GET /v2/associations), so this threshold is applied CLIENT-SIDE to the page of associations fetched. To make it meaningful the page is fetched p_value-ascending (most-significant-first) whenever sort is 'p_value' (the default), and the surviving rows are then re-ordered into the direction you asked for. metadata.p_value_filter_scope, metadata.filtered_total and metadata.p_value_filter_note report exactly what was filtered."
        },
        "p_value_threshold": {
          "type": "number",
          "description": "Alias for p_value: maximum p-value threshold (also applied client-side; see p_value)."
        },
        "sort": {
          "type": "string",
          "description": "Sort field (e.g., 'p_value', 'or_value'). Defaults to 'p_value'. When a p_value threshold is supplied and sort is 'p_value', the upstream page is always fetched ascending so the threshold sees the significant end; metadata.sort_override records this and the returned rows are ordered per your 'direction'."
        },
        "direction": {
          "type": "string",
          "description": "Sort direction ('asc' or 'desc'). Defaults to 'asc'. Combined with a p_value threshold this orders the RETURNED rows; see metadata.sort_override."
        },
        "size": {
          "type": "integer",
          "description": "Number of results to return"
        },
        "page": {
          "type": "integer",
          "description": "Page number for pagination"
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "array",
          "description": "Array of GWAS association objects",
          "items": {
            "type": "object",
            "properties": {
              "association_id": {
                "type": "integer",
                "description": "Unique association identifier"
              },
              "p_value": {
                "type": "number",
                "description": "Statistical p-value"
              },
              "beta": {
                "type": "string",
                "description": "Effect size (beta coefficient)"
              },
              "efo_traits": {
                "type": "array",
                "items": {
                  "type": "object",
                  "properties": {
                    "efo_id": {
                      "type": "string"
                    },
                    "efo_trait": {
                      "type": "string"
                    }
                  }
                }
              },
              "reported_trait": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "accession_id": {
                "type": "string",
                "description": "Study accession ID"
              },
              "pubmed_id": {
                "type": [
                  "string",
                  "integer"
                ],
                "description": "PubMed reference (can be string or integer)"
              },
              "first_author": {
                "type": "string",
                "description": "First author name"
              },
              "snp_effect_allele": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "snp_allele": {
                "type": "array",
                "items": {
                  "type": "object",
                  "properties": {
                    "rs_id": {
                      "type": "string"
                    },
                    "effect_allele": {
                      "type": "string"
                    }
                  }
                }
              },
              "locations": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "mapped_genes": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "ci_lower": {
                "type": "number",
                "description": "Lower confidence interval"
              },
              "ci_upper": {
                "type": "number",
                "description": "Upper confidence interval"
              },
              "risk_frequency": {
                "type": "string",
                "description": "Risk allele frequency"
              },
              "range": {
                "type": "string",
                "description": "Effect size range"
              },
              "_links": {
                "type": "object",
                "description": "API links"
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

**执行与来源：** `network`；上游接口名 `gwas_search_associations`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gwas_search_associations", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gwas_search_associations_by_gene` {#dnakit.databases.gwas_search_associations_by_gene}

Search GWAS Catalog associations by gene name (returns strongest risk allele and p-value fields).

```text
dnakit.databases.gwas_search_associations_by_gene(*, gene_name, size=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_name` | string | 必填（输入 Schema） | Gene symbol (e.g., BRCA1). |
| `size` | integer | None | Max associations to return (default 100; results are sorted by p-value, most significant first). {&quot;default&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene_name": {
          "type": "string",
          "description": "Gene symbol (e.g., BRCA1)."
        },
        "size": {
          "type": "integer",
          "description": "Max associations to return (default 100; results are sorted by p-value, most significant first).",
          "default": 100
        }
      },
      "required": [
        "gene_name"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "gene_name": {
          "type": "string"
        },
        "association_count": {
          "type": "integer"
        },
        "associations": {
          "type": "array",
          "items": {
            "type": "object"
          },
          "description": "Raw GWAS Catalog association objects (v2 API). Each may include fields such as association_id, beta, standard_error, range, p_value, p_value_exponent, p_value_mantissa, efo_traits, mapped_genes, strongest_risk_alleles, study_accession, risk_allele_frequency."
        },
        "total_found": {
          "type": "integer"
        },
        "queried_gene_name": {
          "type": "string",
          "description": "Present only when the requested symbol matched nothing and an unhyphenated spelling was used instead; holds the symbol originally requested."
        },
        "gene_name_note": {
          "type": "string",
          "description": "Present only when the requested symbol matched nothing and an unhyphenated spelling was used instead; explains which spelling produced these results."
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `GWAS_search_associations_by_gene`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gwas_search_associations_by_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gwas_search_snps` {#dnakit.databases.gwas_search_snps}

Search the GWAS Catalog for single nucleotide polymorphisms (SNPs) by rs ID (e.g., 'rs7903146') or mapped gene name (e.g., 'TCF7L2'). Returns SNP details including chromosomal location, mapped genes, and functional class. Use this to look up specific variants or find all GWAS-catalogued SNPs for a gene.

```text
dnakit.databases.gwas_search_snps(*, rs_id=None, mapped_gene=None, size=None, page=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `rs_id` | string | None | dbSNP rs identifier |
| `mapped_gene` | string | None | Gene name or symbol |
| `size` | integer | None | Number of results to return |
| `page` | integer | None | Page number for pagination |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "rs_id": {
          "type": "string",
          "description": "dbSNP rs identifier"
        },
        "mapped_gene": {
          "type": "string",
          "description": "Gene name or symbol"
        },
        "size": {
          "type": "integer",
          "description": "Number of results to return"
        },
        "page": {
          "type": "integer",
          "description": "Page number for pagination"
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "array",
          "description": "Array of SNP objects",
          "items": {
            "type": "object",
            "properties": {
              "rs_id": {
                "type": "string",
                "description": "dbSNP rs identifier"
              },
              "merged": {
                "type": "integer",
                "description": "Merged status"
              },
              "functional_class": {
                "type": "string",
                "description": "Functional classification"
              },
              "last_update_date": {
                "type": "string",
                "description": "Last update timestamp"
              },
              "locations": {
                "type": "array",
                "items": {
                  "type": "object",
                  "properties": {
                    "chromosome_name": {
                      "type": "string"
                    },
                    "chromosome_position": {
                      "type": "integer"
                    },
                    "region": {
                      "type": "object",
                      "properties": {
                        "name": {
                          "type": "string"
                        }
                      }
                    }
                  }
                }
              },
              "alleles": {
                "type": "string",
                "description": "Allele information"
              },
              "most_severe_consequence": {
                "type": "string",
                "description": "Most severe functional consequence"
              },
              "mapped_genes": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "maf": {
                "type": "number",
                "description": "Minor allele frequency"
              },
              "minor_allele": {
                "type": "string",
                "description": "Minor allele"
              },
              "_links": {
                "type": "object",
                "description": "API links"
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

**执行与来源：** `network`；上游接口名 `gwas_search_snps`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gwas_search_snps", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gwas_search_studies` {#dnakit.databases.gwas_search_studies}

Search for GWAS studies by various criteria including disease trait, cohort, GxE interactions, and summary statistics availability.

```text
dnakit.databases.gwas_search_studies(*, disease_trait=None, efo_uri=None, efo_id=None, efo_trait=None, cohort=None, gxe=None, full_pvalue_set=None, size=None, page=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `disease_trait` | string | None | Disease or trait name for text-based search (e.g., &#x27;diabetes&#x27;, &#x27;breast cancer&#x27;) |
| `efo_uri` | string | None | Full EFO ontology URI (e.g., &#x27;http://www.ebi.ac.uk/efo/EFO_0001645&#x27;) |
| `efo_id` | string | None | EFO/OBA term ID (e.g., &#x27;EFO_0001645&#x27;, &#x27;OBA_2050062&#x27;). Recommended for reliable trait filtering. |
| `efo_trait` | string | None | Exact EFO trait label. Use when you know the canonical trait string. |
| `cohort` | string | None | Cohort name (e.g., &#x27;UKB&#x27; for UK Biobank) |
| `gxe` | boolean | None | Filter for Gene-by-Environment interaction studies |
| `full_pvalue_set` | boolean | None | Filter for studies with full summary statistics |
| `size` | integer | None | Number of results to return |
| `page` | integer | None | Page number for pagination |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "disease_trait": {
          "type": "string",
          "description": "Disease or trait name for text-based search (e.g., 'diabetes', 'breast cancer')"
        },
        "efo_uri": {
          "type": "string",
          "description": "Full EFO ontology URI (e.g., 'http://www.ebi.ac.uk/efo/EFO_0001645')"
        },
        "efo_id": {
          "type": "string",
          "description": "EFO/OBA term ID (e.g., 'EFO_0001645', 'OBA_2050062'). Recommended for reliable trait filtering."
        },
        "efo_trait": {
          "type": "string",
          "description": "Exact EFO trait label. Use when you know the canonical trait string."
        },
        "cohort": {
          "type": "string",
          "description": "Cohort name (e.g., 'UKB' for UK Biobank)"
        },
        "gxe": {
          "type": "boolean",
          "description": "Filter for Gene-by-Environment interaction studies"
        },
        "full_pvalue_set": {
          "type": "boolean",
          "description": "Filter for studies with full summary statistics"
        },
        "size": {
          "type": "integer",
          "description": "Number of results to return"
        },
        "page": {
          "type": "integer",
          "description": "Page number for pagination"
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "array",
          "description": "Array of GWAS study objects",
          "items": {
            "type": "object",
            "properties": {
              "accession_id": {
                "type": "string",
                "description": "Study accession identifier"
              },
              "disease_trait": {
                "type": "string",
                "description": "Disease or trait name"
              },
              "efo_traits": {
                "type": "array",
                "items": {
                  "type": "object",
                  "properties": {
                    "efo_id": {
                      "type": "string"
                    },
                    "efo_trait": {
                      "type": "string"
                    }
                  }
                }
              },
              "initial_sample_size": {
                "type": "string",
                "description": "Initial sample size description"
              },
              "replication_sample_size": {
                "type": "string",
                "description": "Replication sample size"
              },
              "discovery_ancestry": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "replication_ancestry": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "cohort": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "genotyping_technologies": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "platforms": {
                "type": "string",
                "description": "Genotyping platforms used"
              },
              "snp_count": {
                "type": "integer",
                "description": "Number of SNPs analyzed"
              },
              "imputed": {
                "type": "boolean",
                "description": "Whether data was imputed"
              },
              "pooled": {
                "type": "boolean",
                "description": "Whether data was pooled"
              },
              "gxe": {
                "type": "boolean",
                "description": "Gene-environment interaction study"
              },
              "gxg": {
                "type": "boolean",
                "description": "Gene-gene interaction study"
              },
              "full_summary_stats_available": {
                "type": "boolean",
                "description": "Full summary statistics available"
              },
              "full_summary_stats": {
                "type": "string",
                "description": "URL to full summary statistics"
              },
              "pubmed_id": {
                "type": [
                  "string",
                  "integer"
                ],
                "description": "PubMed reference (can be string or integer)"
              },
              "terms_of_license": {
                "type": "string",
                "description": "Data usage license"
              },
              "array_manufacturer": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "_links": {
                "type": "object",
                "description": "API links"
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

**执行与来源：** `network`；上游接口名 `gwas_search_studies`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gwas_search_studies", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.igsr_list_data_collections` {#dnakit.databases.igsr_list_data_collections}

List available data collections from the International Genome Sample Resource (IGSR). Data collections include the original 1000 Genomes Project phases, the 30x high-coverage resequencing, Human Genome Structural Variation Consortium (HGSVC), Human Pangenome Reference Consortium (HPRC), and other related studies. Returns collection codes, titles, sample/population counts, and available data types (sequence, alignment, variants). Use this to discover what genomics datasets are available through IGSR.

```text
dnakit.databases.igsr_list_data_collections(*, operation=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: list_data_collections). {&quot;enum&quot;: &#91;&quot;list_data_collections&quot;&#93;, &quot;default&quot;: &quot;list_data_collections&quot;} |
| `limit` | integer | None | Maximum number of results to return. {&quot;default&quot;: 50, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "enum": [
            "list_data_collections"
          ],
          "default": "list_data_collections",
          "description": "Operation type (fixed: list_data_collections)."
        },
        "limit": {
          "type": "integer",
          "default": 50,
          "minimum": 1,
          "maximum": 100,
          "description": "Maximum number of results to return."
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "total": {
              "type": "integer",
              "description": "Total number of matching records in IGSR, independent of the `limit` page size."
            },
            "returned": {
              "type": "integer",
              "description": "Number of records included in this response (at most `limit`)."
            },
            "collections": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "code": {
                    "type": "string"
                  },
                  "title": {
                    "type": "string"
                  },
                  "short_title": {
                    "type": "string"
                  },
                  "sample_count": {
                    "type": "integer"
                  },
                  "population_count": {
                    "type": "integer"
                  },
                  "data_types": {
                    "type": "array",
                    "items": {
                      "type": "string"
                    }
                  },
                  "website": {
                    "type": [
                      "string",
                      "null"
                    ]
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

**执行与来源：** `network`；上游接口名 `IGSR_list_data_collections`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("igsr_list_data_collections", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.igsr_search_populations` {#dnakit.databases.igsr_search_populations}

Search 1000 Genomes Project populations from the International Genome Sample Resource (IGSR). The 1000 Genomes Project characterized genetic variation across 26 populations grouped into 5 superpopulations: AFR (African), AMR (Admixed American), EAS (East Asian), EUR (European), SAS (South Asian). Returns population codes, names, geographic coordinates, sample counts, and superpopulation assignments. Use this to find populations for allele frequency comparisons, population stratification analysis, or ancestry-specific variant studies.

```text
dnakit.databases.igsr_search_populations(*, operation=None, superpopulation=None, query=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search_populations). {&quot;enum&quot;: &#91;&quot;search_populations&quot;&#93;, &quot;default&quot;: &quot;search_populations&quot;} |
| `superpopulation` | string &#124; null | None | Filter by superpopulation code: AFR (African), AMR (Admixed American), EAS (East Asian), EUR (European), SAS (South Asian). Leave empty for all. |
| `query` | string &#124; null | None | Free-text search query to find populations by name or description (e.g., &#x27;Chinese&#x27;, &#x27;Yoruba&#x27;, &#x27;British&#x27;). |
| `limit` | integer | None | Maximum number of results to return. {&quot;default&quot;: 25, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "enum": [
            "search_populations"
          ],
          "default": "search_populations",
          "description": "Operation type (fixed: search_populations)."
        },
        "superpopulation": {
          "type": [
            "string",
            "null"
          ],
          "description": "Filter by superpopulation code: AFR (African), AMR (Admixed American), EAS (East Asian), EUR (European), SAS (South Asian). Leave empty for all."
        },
        "query": {
          "type": [
            "string",
            "null"
          ],
          "description": "Free-text search query to find populations by name or description (e.g., 'Chinese', 'Yoruba', 'British')."
        },
        "limit": {
          "type": "integer",
          "default": 25,
          "minimum": 1,
          "maximum": 100,
          "description": "Maximum number of results to return."
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "total": {
              "type": "integer",
              "description": "Total number of matching records in IGSR, independent of the `limit` page size."
            },
            "returned": {
              "type": "integer",
              "description": "Number of records included in this response (at most `limit`)."
            },
            "populations": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "code": {
                    "type": "string"
                  },
                  "name": {
                    "type": "string"
                  },
                  "description": {
                    "type": "string"
                  },
                  "sample_count": {
                    "type": "integer"
                  },
                  "superpopulation_code": {
                    "type": "string"
                  },
                  "superpopulation_name": {
                    "type": "string"
                  },
                  "latitude": {
                    "type": [
                      "number",
                      "null"
                    ]
                  },
                  "longitude": {
                    "type": [
                      "number",
                      "null"
                    ]
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

**执行与来源：** `network`；上游接口名 `IGSR_search_populations`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("igsr_search_populations", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.igsr_search_samples` {#dnakit.databases.igsr_search_samples}

Search samples from the 1000 Genomes Project by population, data collection, or sample name. Returns sample names (e.g., HG00096, NA12878), sex, BioSample IDs, population assignments, and available data collections. Use this to find samples from specific populations for variant analysis, to identify samples with particular sequencing data types (WGS, exome, RNA-seq), or to look up specific well-characterized samples like NA12878 (reference genome sample).

```text
dnakit.databases.igsr_search_samples(*, operation=None, population=None, data_collection=None, sample_name=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search_samples). {&quot;enum&quot;: &#91;&quot;search_samples&quot;&#93;, &quot;default&quot;: &quot;search_samples&quot;} |
| `population` | string &#124; null | None | Filter by population code (e.g., &#x27;CHB&#x27;, &#x27;YRI&#x27;, &#x27;GBR&#x27;, &#x27;CEU&#x27;). Use IGSR_search_populations to find valid codes. |
| `data_collection` | string &#124; null | None | Filter by data collection title (e.g., &#x27;1000 Genomes 30x on GRCh38&#x27;, &#x27;1000 Genomes on GRCh38&#x27;). |
| `sample_name` | string &#124; null | None | Look up specific sample by name (e.g., &#x27;HG00096&#x27;, &#x27;NA12878&#x27;). |
| `limit` | integer | None | Maximum number of results to return. {&quot;default&quot;: 25, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "enum": [
            "search_samples"
          ],
          "default": "search_samples",
          "description": "Operation type (fixed: search_samples)."
        },
        "population": {
          "type": [
            "string",
            "null"
          ],
          "description": "Filter by population code (e.g., 'CHB', 'YRI', 'GBR', 'CEU'). Use IGSR_search_populations to find valid codes."
        },
        "data_collection": {
          "type": [
            "string",
            "null"
          ],
          "description": "Filter by data collection title (e.g., '1000 Genomes 30x on GRCh38', '1000 Genomes on GRCh38')."
        },
        "sample_name": {
          "type": [
            "string",
            "null"
          ],
          "description": "Look up specific sample by name (e.g., 'HG00096', 'NA12878')."
        },
        "limit": {
          "type": "integer",
          "default": 25,
          "minimum": 1,
          "maximum": 100,
          "description": "Maximum number of results to return."
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "total": {
              "type": "integer",
              "description": "Total number of matching records in IGSR, independent of the `limit` page size."
            },
            "returned": {
              "type": "integer",
              "description": "Number of records included in this response (at most `limit`)."
            },
            "samples": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "name": {
                    "type": "string"
                  },
                  "sex": {
                    "type": "string"
                  },
                  "biosample_id": {
                    "type": "string"
                  },
                  "populations": {
                    "type": "array",
                    "items": {
                      "type": "object",
                      "properties": {
                        "code": {
                          "type": "string"
                        },
                        "name": {
                          "type": "string"
                        },
                        "superpopulation": {
                          "type": "string"
                        }
                      }
                    }
                  },
                  "data_collections": {
                    "type": "array",
                    "items": {
                      "type": "string"
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

**执行与来源：** `network`；上游接口名 `IGSR_search_samples`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("igsr_search_samples", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ldlink_get_proxies` {#dnakit.databases.ldlink_get_proxies}

Find variants in linkage disequilibrium (LD proxies) with a query SNP, population-specific, via NIH LDlink's LDproxy module. Essential for GWAS interpretation: a lead SNP is often just a tag for the true causal variant, so high-R2 proxies (especially in coding/regulatory regions) are better mechanistic candidates. Returns proxy rsIDs with R2, coordinates, alleles, MAF, and distance. Requires a free LDlink token (LDLINK_TOKEN env var; register at ldlink.nih.gov).

```text
dnakit.databases.ldlink_get_proxies(*, variant, population=None, r2_threshold=None, genome_build=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant` | string | 必填（输入 Schema） | Query SNP rsID, e.g. &#x27;rs7903146&#x27;. |
| `population` | string &#124; null | None | 1000 Genomes population code (default &#x27;CEU&#x27;). Examples: &#x27;CEU&#x27; (European), &#x27;YRI&#x27; (African), &#x27;CHB&#x27; (East Asian), or super-populations &#x27;EUR&#x27;/&#x27;AFR&#x27;/&#x27;EAS&#x27;/&#x27;SAS&#x27;/&#x27;AMR&#x27;. |
| `r2_threshold` | number &#124; null | None | Minimum R2 to report (default 0.8; high LD). |
| `genome_build` | string &#124; null | None | &#x27;grch38&#x27; (default) or &#x27;grch37&#x27;. |
| `limit` | integer &#124; null | None | Max proxies to return (1-500, default 50). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "variant": {
          "type": "string",
          "description": "Query SNP rsID, e.g. 'rs7903146'."
        },
        "population": {
          "type": [
            "string",
            "null"
          ],
          "description": "1000 Genomes population code (default 'CEU'). Examples: 'CEU' (European), 'YRI' (African), 'CHB' (East Asian), or super-populations 'EUR'/'AFR'/'EAS'/'SAS'/'AMR'."
        },
        "r2_threshold": {
          "type": [
            "number",
            "null"
          ],
          "description": "Minimum R2 to report (default 0.8; high LD)."
        },
        "genome_build": {
          "type": [
            "string",
            "null"
          ],
          "description": "'grch38' (default) or 'grch37'."
        },
        "limit": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Max proxies to return (1-500, default 50)."
        }
      },
      "required": [
        "variant"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "rsid": {
                "type": "string"
              },
              "r2": {
                "type": [
                  "number",
                  "null"
                ]
              },
              "coord": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "alleles": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "maf": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "distance": {
                "type": [
                  "string",
                  "null"
                ]
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

**执行与来源：** `network`；上游接口名 `LDlink_get_proxies`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ldlink_get_proxies", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.mave_db_get_effect_matrix` {#dnakit.databases.mave_db_get_effect_matrix}

One-shot DMS data loader: fetch a MaveDB score set, parse HGVS to (ref_aa, position, alt_aa), filter to single missense, optionally verify position numbering against a UniProt canonical sequence, and reshape into a (20 amino acids × n_positions) effect matrix ready for downstream analysis. Hides the boilerplate that variant-effect-prediction benchmarking and residue-mechanism interpretation skills used to inline. Returns the matrix, the position list, the standard amino acid order, and explicit numbering / drop / score-field metadata so the caller can audit what was loaded.

```text
dnakit.databases.mave_db_get_effect_matrix(*, urn, uniprot_id=None, score_field=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `urn` | string | 必填（输入 Schema） | MaveDB score set URN (e.g., &#x27;urn:mavedb:00000115-a-7&#x27;). Obtain from MaveDB_search_score_sets. |
| `uniprot_id` | string &#124; null | None | Optional UniProt accession (e.g., &#x27;P01116&#x27; for KRAS). If supplied, the tool fetches the canonical sequence and verifies the MaveDB position numbering against it via a landmark check; sets numbering_offset and surfaces a warning in numbering_check if the positions don&#x27;t match (silent off-by-N is the #1 silent bug in DMS analysis). |
| `score_field` | string &#124; null | None | Optional explicit column name to use as the variant score (e.g., &#x27;score&#x27;, &#x27;ddG&#x27;, &#x27;fitness&#x27;). If omitted, tries &#x27;score&#x27; → &#x27;ddG&#x27; → &#x27;fitness&#x27; → first numeric non-HGVS field in that order. The actual field used is reported in data.score_field_used. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "urn": {
          "type": "string",
          "description": "MaveDB score set URN (e.g., 'urn:mavedb:00000115-a-7'). Obtain from MaveDB_search_score_sets."
        },
        "uniprot_id": {
          "type": [
            "string",
            "null"
          ],
          "description": "Optional UniProt accession (e.g., 'P01116' for KRAS). If supplied, the tool fetches the canonical sequence and verifies the MaveDB position numbering against it via a landmark check; sets numbering_offset and surfaces a warning in numbering_check if the positions don't match (silent off-by-N is the #1 silent bug in DMS analysis)."
        },
        "score_field": {
          "type": [
            "string",
            "null"
          ],
          "description": "Optional explicit column name to use as the variant score (e.g., 'score', 'ddG', 'fitness'). If omitted, tries 'score' → 'ddG' → 'fitness' → first numeric non-HGVS field in that order. The actual field used is reported in data.score_field_used."
        }
      },
      "required": [
        "urn"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "urn": {
              "type": "string"
            },
            "matrix": {
              "type": "array",
              "description": "(20 × n_positions) effect matrix. Rows = 20 standard amino acids in 'ACDEFGHIKLMNPQRSTVWY' order. Cells = score for that (alt_aa, position) variant, or null for unmeasured / WT-diagonal.",
              "items": {
                "type": "array",
                "items": {
                  "type": [
                    "number",
                    "null"
                  ]
                }
              }
            },
            "positions": {
              "type": "array",
              "description": "Sorted list of 1-based residue positions (length = matrix's column count).",
              "items": {
                "type": "integer"
              }
            },
            "amino_acid_order": {
              "type": "string",
              "description": "Always 'ACDEFGHIKLMNPQRSTVWY' for consistency across downstream skills."
            },
            "shape": {
              "type": "array",
              "items": {
                "type": "integer"
              }
            },
            "n_parsed_single_missense": {
              "type": "integer"
            },
            "n_dropped": {
              "type": "integer"
            },
            "n_total_rows": {
              "type": "integer"
            },
            "score_field_used": {
              "type": [
                "string",
                "null"
              ]
            },
            "numbering_check": {
              "type": [
                "object",
                "null"
              ]
            },
            "numbering_offset": {
              "type": [
                "integer",
                "null"
              ]
            }
          },
          "required": [
            "matrix",
            "positions",
            "amino_acid_order"
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

**执行与来源：** `network`；上游接口名 `MaveDB_get_effect_matrix`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("mave_db_get_effect_matrix", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.my_variant_get_pathogenicity_scores` {#dnakit.databases.my_variant_get_pathogenicity_scores}

Get comprehensive pathogenicity prediction scores for a variant from dbNSFP. Returns REVEL, CADD, AlphaMissense, SIFT, PolyPhen-2, MetaRNN, GERP, PhyloP, and PhastCons scores in a single call. Use this when MyVariant_query_variants does not return dbnsfp data, or when you specifically need REVEL/AlphaMissense scores for ACMG PP3/BP4 classification. Coordinates default to hg19/GRCh37 -- set assembly='hg38' for hg38/GRCh38 (e.g. rs4244285 is chr10:g.96541616G>C on hg19 and chr10:g.94781859G>C on hg38). An HGVS variant_id must be written in the SAME assembly you request: this endpoint looks ids up verbatim and never lifts coordinates over, so an hg19 id with assembly='hg38' is a 404, not a translation. An rsID works in either assembly. The response reports the frame in a top-level 'coordinate_assembly' key. MyVariant sometimes files several records under one rsID and resolves the rsID to the highest-scoring one, which can be an identifier-only stub; when that happens the sibling records that do carry the requested fields are reported in the top-level 'sibling_records_with_requested_fields', 'sibling_variant_ids_with_requested_fields' and 'sibling_record_note' keys, and 'data' still holds the record MyVariant resolved.

```text
dnakit.databases.my_variant_get_pathogenicity_scores(*, variant_id, fields=None, assembly=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant_id` | string | 必填（输入 Schema） | Variant ID: rsID (e.g., rs45478192) or HGVS genomic (e.g., chr16:g.23635348A&gt;C). Both formats work; an HGVS id must be expressed in the assembly given by the &#x27;assembly&#x27; parameter. |
| `fields` | string | None | Fields to return (pre-configured for pathogenicity scores) {&quot;default&quot;: &quot;dbnsfp.revel.score,cadd.phred,dbnsfp.alphamissense.score,dbnsfp.alphamissense.pred,dbnsfp.sift.score,dbnsfp.sift.pred,dbnsfp.polyphen2.hdiv.score,dbnsfp.polyphen2.hdiv.pred,dbnsfp.metarnn.score,dbnsfp.metarnn.pred,cadd.gerp.rs,dbnsfp.phylop.100way_vertebrate.rankscore,dbnsfp.phastcons.100way_vertebrate.rankscore,dbnsfp.vest4.score,dbnsfp.mutationtaster.pred,clinvar.rcv.clinical_significance,dbsnp.rsid&quot;} |
| `assembly` | string | None | Reference genome assembly for both the HGVS variant_id and the returned coordinates: &#x27;hg19&#x27; (GRCh37, the default) or &#x27;hg38&#x27; (GRCh38). {&quot;enum&quot;: &#91;&quot;hg19&quot;, &quot;hg38&quot;&#93;, &quot;default&quot;: &quot;hg19&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "variant_id": {
          "type": "string",
          "description": "Variant ID: rsID (e.g., rs45478192) or HGVS genomic (e.g., chr16:g.23635348A>C). Both formats work; an HGVS id must be expressed in the assembly given by the 'assembly' parameter."
        },
        "fields": {
          "type": "string",
          "description": "Fields to return (pre-configured for pathogenicity scores)",
          "default": "dbnsfp.revel.score,cadd.phred,dbnsfp.alphamissense.score,dbnsfp.alphamissense.pred,dbnsfp.sift.score,dbnsfp.sift.pred,dbnsfp.polyphen2.hdiv.score,dbnsfp.polyphen2.hdiv.pred,dbnsfp.metarnn.score,dbnsfp.metarnn.pred,cadd.gerp.rs,dbnsfp.phylop.100way_vertebrate.rankscore,dbnsfp.phastcons.100way_vertebrate.rankscore,dbnsfp.vest4.score,dbnsfp.mutationtaster.pred,clinvar.rcv.clinical_significance,dbsnp.rsid"
        },
        "assembly": {
          "type": "string",
          "description": "Reference genome assembly for both the HGVS variant_id and the returned coordinates: 'hg19' (GRCh37, the default) or 'hg38' (GRCh38).",
          "enum": [
            "hg19",
            "hg38"
          ],
          "default": "hg19"
        }
      },
      "required": [
        "variant_id"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

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

**执行与来源：** `network`；上游接口名 `MyVariant_get_pathogenicity_scores`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("my_variant_get_pathogenicity_scores", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.my_variant_get_variant_annotation` {#dnakit.databases.my_variant_get_variant_annotation}

Get detailed annotations for a specific variant by its HGVS ID (e.g., 'chr1:g.2186318G>A') or rsID. Returns information from dbSNP, ClinVar, CADD, and other sources. Coordinates default to hg19/GRCh37 -- set assembly='hg38' for hg38/GRCh38 (e.g. rs4244285 is chr10:g.96541616G>C on hg19 and chr10:g.94781859G>C on hg38). An HGVS variant_id must be written in the SAME assembly you request: this endpoint looks ids up verbatim and never lifts coordinates over, so an hg19 id with assembly='hg38' is a 404, not a translation. An rsID works in either assembly. The response reports the frame in a top-level 'coordinate_assembly' key. MyVariant sometimes files several records under one rsID and resolves the rsID to the highest-scoring one, which can be an identifier-only stub; when that happens the sibling records that do carry the requested fields are reported in the top-level 'sibling_records_with_requested_fields', 'sibling_variant_ids_with_requested_fields' and 'sibling_record_note' keys, and 'data' still holds the record MyVariant resolved.

```text
dnakit.databases.my_variant_get_variant_annotation(*, variant_id, fields=None, assembly=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant_id` | string | 必填（输入 Schema） | HGVS variant ID (e.g., &#x27;chr1:g.2186318G&gt;A&#x27;) or rsID (e.g., &#x27;rs123&#x27;). An HGVS id must be expressed in the assembly given by the &#x27;assembly&#x27; parameter. |
| `fields` | string | None | Fields to return. Useful fields: dbsnp, clinvar, cadd, gnomad_genome, gnomad_exome, dbnsfp (functional predictions). {&quot;default&quot;: &quot;dbsnp,clinvar,cadd,gnomad_genome,dbnsfp&quot;} |
| `assembly` | string | None | Reference genome assembly for both the HGVS variant_id and the returned coordinates: &#x27;hg19&#x27; (GRCh37, the default) or &#x27;hg38&#x27; (GRCh38). {&quot;enum&quot;: &#91;&quot;hg19&quot;, &quot;hg38&quot;&#93;, &quot;default&quot;: &quot;hg19&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "variant_id": {
          "type": "string",
          "description": "HGVS variant ID (e.g., 'chr1:g.2186318G>A') or rsID (e.g., 'rs123'). An HGVS id must be expressed in the assembly given by the 'assembly' parameter."
        },
        "fields": {
          "type": "string",
          "description": "Fields to return. Useful fields: dbsnp, clinvar, cadd, gnomad_genome, gnomad_exome, dbnsfp (functional predictions).",
          "default": "dbsnp,clinvar,cadd,gnomad_genome,dbnsfp"
        },
        "assembly": {
          "type": "string",
          "description": "Reference genome assembly for both the HGVS variant_id and the returned coordinates: 'hg19' (GRCh37, the default) or 'hg38' (GRCh38).",
          "enum": [
            "hg19",
            "hg38"
          ],
          "default": "hg19"
        }
      },
      "required": [
        "variant_id"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "_id": {
          "type": "string"
        },
        "dbsnp": {
          "type": "object",
          "properties": {
            "rsid": {
              "type": "string"
            },
            "gene": {
              "type": "object"
            },
            "alleles": {
              "type": "array"
            }
          }
        },
        "clinvar": {
          "type": "object",
          "description": "ClinVar clinical significance and condition data"
        },
        "cadd": {
          "type": "object",
          "properties": {
            "phred": {
              "type": "number",
              "description": "CADD Phred-scaled score (>20 = top 1% damaging)"
            },
            "consequence": {
              "type": "string"
            }
          }
        },
        "gnomad_genome": {
          "type": "object",
          "description": "gnomAD genome allele frequencies"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `MyVariant_get_variant_annotation`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("my_variant_get_variant_annotation", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.my_variant_query_variants` {#dnakit.databases.my_variant_query_variants}

Search for genetic variants by various criteria including rsID, gene, chromosome position, or clinical significance. Returns variant annotations from sources like dbSNP, ClinVar, gnomAD, and CADD. Coordinates default to hg19/GRCh37 -- set assembly='hg38' for hg38/GRCh38 (e.g. rs4244285 is chr10:g.96541616G>C on hg19 and chr10:g.94781859G>C on hg38). Any HGVS coordinate in the query must be written in the assembly you request; an rsID works in either. The response reports the frame in a top-level 'coordinate_assembly' key.

```text
dnakit.databases.my_variant_query_variants(*, query=None, fields=None, size=None, assembly=None, variant_id=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | None | Search query. Examples: &#x27;rs58991260&#x27; (rsID), &#x27;chr7:g.55249071G&gt;A&#x27; (HGVS), &#x27;clinvar.gene.symbol:BRCA1&#x27; (gene in ClinVar), &#x27;cadd.phred:&gt;30&#x27; (high CADD score). |
| `fields` | string | None | Comma-separated fields to return. Common: dbsnp, clinvar, cadd, gnomad_genome, dbnsfp. {&quot;default&quot;: &quot;dbsnp.rsid,clinvar.rcv.clinical_significance,cadd.phred,gnomad_genome.af.af&quot;} |
| `size` | integer | None | Maximum number of results (1-100). {&quot;default&quot;: 10} |
| `assembly` | string | None | Reference genome assembly for the coordinates in the query and in the results: &#x27;hg19&#x27; (GRCh37, the default) or &#x27;hg38&#x27; (GRCh38). {&quot;enum&quot;: &#91;&quot;hg19&quot;, &quot;hg38&quot;&#93;, &quot;default&quot;: &quot;hg19&quot;} |
| `variant_id` | string | None | DNAKit optional alias for query |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "query": {
          "type": "string",
          "description": "Search query. Examples: 'rs58991260' (rsID), 'chr7:g.55249071G>A' (HGVS), 'clinvar.gene.symbol:BRCA1' (gene in ClinVar), 'cadd.phred:>30' (high CADD score)."
        },
        "fields": {
          "type": "string",
          "description": "Comma-separated fields to return. Common: dbsnp, clinvar, cadd, gnomad_genome, dbnsfp.",
          "default": "dbsnp.rsid,clinvar.rcv.clinical_significance,cadd.phred,gnomad_genome.af.af"
        },
        "size": {
          "type": "integer",
          "description": "Maximum number of results (1-100).",
          "default": 10
        },
        "assembly": {
          "type": "string",
          "description": "Reference genome assembly for the coordinates in the query and in the results: 'hg19' (GRCh37, the default) or 'hg38' (GRCh38).",
          "enum": [
            "hg19",
            "hg38"
          ],
          "default": "hg19"
        },
        "variant_id": {
          "type": "string",
          "description": "DNAKit optional alias for query"
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "took": {
          "type": "integer"
        },
        "total": {
          "type": "integer"
        },
        "hits": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "_id": {
                "type": "string",
                "description": "Variant ID in HGVS format"
              },
              "dbsnp": {
                "type": "object",
                "properties": {
                  "rsid": {
                    "type": "string"
                  }
                }
              },
              "clinvar": {
                "type": "object",
                "properties": {
                  "rcv": {
                    "type": "array",
                    "description": "ClinVar RCV records"
                  }
                }
              },
              "cadd": {
                "type": "object",
                "properties": {
                  "phred": {
                    "type": "number",
                    "description": "CADD Phred score"
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

??? info "DNAKit 参数扩展规则"

    ```json
    {
      "aliases": {
        "variant_id": "query"
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `MyVariant_query_variants`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("my_variant_query_variants", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.onco_kb_annotate_variant` {#dnakit.databases.onco_kb_annotate_variant}

Annotate a cancer variant for oncogenic potential and treatment implications using OncoKB. Returns evidence level, oncogenicity, mutation effect, and FDA-approved therapies. Requires ONCOKB_API_TOKEN env var (demo mode available for BRAF, TP53, ROS1). Example: BRAF V600E in melanoma.

```text
dnakit.databases.onco_kb_annotate_variant(*, operation=None, gene=None, gene_symbol=None, variant=None, tumor_type=None, alteration=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: annotate_variant) {&quot;const&quot;: &quot;annotate_variant&quot;, &quot;default&quot;: &quot;annotate_variant&quot;} |
| `gene` | string | None | Gene symbol (e.g., BRAF, EGFR, TP53, KRAS) |
| `gene_symbol` | string | None | Alias for gene. Gene symbol (e.g., BRAF, EGFR, TP53, KRAS). |
| `variant` | string | None | Variant notation - protein change (e.g., V600E, T790M, G12D) |
| `tumor_type` | string | None | Optional OncoTree tumor type code (e.g., MEL for melanoma, LUAD for lung adenocarcinoma) |
| `alteration` | string | None | Alias for variant. Protein change notation (e.g., &quot;V600E&quot;, &quot;T790M&quot;, &quot;G12D&quot;). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "annotate_variant",
          "description": "Operation type (fixed: annotate_variant)",
          "default": "annotate_variant"
        },
        "gene": {
          "type": "string",
          "description": "Gene symbol (e.g., BRAF, EGFR, TP53, KRAS)"
        },
        "gene_symbol": {
          "type": "string",
          "description": "Alias for gene. Gene symbol (e.g., BRAF, EGFR, TP53, KRAS)."
        },
        "variant": {
          "type": "string",
          "description": "Variant notation - protein change (e.g., V600E, T790M, G12D)"
        },
        "tumor_type": {
          "type": "string",
          "description": "Optional OncoTree tumor type code (e.g., MEL for melanoma, LUAD for lung adenocarcinoma)"
        },
        "alteration": {
          "type": "string",
          "description": "Alias for variant. Protein change notation (e.g., \"V600E\", \"T790M\", \"G12D\")."
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "query": {
              "type": "object",
              "properties": {
                "id": {
                  "type": "null"
                },
                "referenceGenome": {
                  "type": "string"
                },
                "hugoSymbol": {
                  "type": "string"
                },
                "entrezGeneId": {
                  "type": [
                    "integer",
                    "number"
                  ]
                },
                "alteration": {
                  "type": "string"
                },
                "alterationType": {
                  "type": "null"
                },
                "svType": {
                  "type": "null"
                },
                "tumorType": {
                  "type": [
                    "null",
                    "string"
                  ]
                },
                "consequence": {
                  "type": "null"
                },
                "proteinStart": {
                  "type": "null"
                },
                "proteinEnd": {
                  "type": "null"
                },
                "hgvs": {
                  "type": "null"
                },
                "hgvsInfo": {
                  "type": "null"
                },
                "canonicalTranscript": {
                  "type": "null"
                },
                "germline": {
                  "type": "boolean"
                }
              }
            },
            "geneExist": {
              "type": "boolean"
            },
            "variantExist": {
              "type": "boolean"
            },
            "mutationEffect": {
              "type": "object",
              "properties": {
                "knownEffect": {
                  "type": "string"
                },
                "description": {
                  "type": "string"
                },
                "citations": {
                  "type": "object",
                  "properties": {
                    "pmids": {
                      "type": "array",
                      "items": {
                        "type": "string"
                      }
                    },
                    "abstracts": {
                      "type": "array"
                    }
                  }
                }
              }
            },
            "highestSensitiveLevel": {
              "type": "string"
            },
            "highestResistanceLevel": {
              "type": "null"
            },
            "highestDiagnosticImplicationLevel": {
              "type": [
                "null",
                "string"
              ]
            },
            "highestPrognosticImplicationLevel": {
              "type": "null"
            },
            "geneSummary": {
              "type": "string"
            },
            "variantSummary": {
              "type": "string"
            },
            "tumorTypeSummary": {
              "type": "string"
            },
            "prognosticSummary": {
              "type": "string"
            },
            "diagnosticSummary": {
              "type": "string"
            },
            "diagnosticImplications": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "levelOfEvidence": {
                    "type": "string"
                  },
                  "alterations": {
                    "type": "array",
                    "items": {
                      "type": "string"
                    }
                  },
                  "tumorType": {
                    "type": "object",
                    "properties": {
                      "id": {
                        "type": [
                          "integer",
                          "number"
                        ]
                      },
                      "code": {
                        "type": "string"
                      },
                      "color": {
                        "type": "string"
                      },
                      "name": {
                        "type": "string"
                      },
                      "mainType": {
                        "type": "object",
                        "properties": {
                          "id": {
                            "type": "null"
                          },
                          "name": {
                            "type": "string"
                          },
                          "tumorForm": {
                            "type": "string"
                          }
                        }
                      },
                      "tissue": {
                        "type": "string"
                      },
                      "children": {
                        "type": "object"
                      },
                      "parent": {
                        "type": "string"
                      },
                      "level": {
                        "type": [
                          "integer",
                          "number"
                        ]
                      },
                      "tumorForm": {
                        "type": "string"
                      }
                    }
                  },
                  "pmids": {
                    "type": "array",
                    "items": {
                      "type": "string"
                    }
                  },
                  "abstracts": {
                    "type": "array"
                  },
                  "description": {
                    "type": "string"
                  }
                }
              }
            },
            "prognosticImplications": {
              "type": "array"
            },
            "treatments": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "alterations": {
                    "type": "array",
                    "items": {
                      "type": "string"
                    }
                  },
                  "drugs": {
                    "type": "array",
                    "items": {
                      "type": "object",
                      "properties": {
                        "ncitCode": {
                          "type": "string"
                        },
                        "drugName": {
                          "type": "string"
                        }
                      }
                    }
                  },
                  "approvedIndications": {
                    "type": "array"
                  },
                  "level": {
                    "type": "string"
                  },
                  "fdaLevel": {
                    "type": "string"
                  },
                  "levelAssociatedCancerType": {
                    "type": "object",
                    "properties": {
                      "id": {
                        "type": [
                          "integer",
                          "number"
                        ]
                      },
                      "code": {
                        "type": "string"
                      },
                      "color": {
                        "type": "string"
                      },
                      "name": {
                        "type": "string"
                      },
                      "mainType": {
                        "type": "object",
                        "properties": {
                          "id": {
                            "type": "null"
                          },
                          "name": {
                            "type": "string"
                          },
                          "tumorForm": {
                            "type": "string"
                          }
                        }
                      },
                      "tissue": {
                        "type": "string"
                      },
                      "children": {
                        "type": "object"
                      },
                      "parent": {
                        "type": [
                          "null",
                          "string"
                        ]
                      },
                      "level": {
                        "type": [
                          "integer",
                          "number"
                        ]
                      },
                      "tumorForm": {
                        "type": "string"
                      }
                    }
                  },
                  "levelExcludedCancerTypes": {
                    "type": "array",
                    "items": {
                      "type": "object",
                      "properties": {
                        "id": {
                          "type": [
                            "integer",
                            "number"
                          ]
                        },
                        "code": {
                          "type": "string"
                        },
                        "color": {
                          "type": "string"
                        },
                        "name": {
                          "type": "string"
                        },
                        "mainType": {
                          "type": "object",
                          "properties": {
                            "id": {
                              "type": "null"
                            },
                            "name": {
                              "type": "string"
                            },
                            "tumorForm": {
                              "type": "string"
                            }
                          }
                        },
                        "tissue": {
                          "type": "string"
                        },
                        "children": {
                          "type": "object"
                        },
                        "parent": {
                          "type": [
                            "null",
                            "string"
                          ]
                        },
                        "level": {
                          "type": [
                            "integer",
                            "number"
                          ]
                        },
                        "tumorForm": {
                          "type": "string"
                        }
                      }
                    }
                  },
                  "pmids": {
                    "type": "array",
                    "items": {
                      "type": "string"
                    }
                  },
                  "abstracts": {
                    "type": "array",
                    "items": {
                      "type": "object",
                      "properties": {
                        "link": {
                          "type": "string"
                        },
                        "abstract": {
                          "type": "string"
                        }
                      }
                    }
                  },
                  "description": {
                    "type": "string"
                  }
                }
              }
            },
            "dataVersion": {
              "type": "string"
            },
            "lastUpdate": {
              "type": "string"
            },
            "alleleExist": {
              "type": "boolean"
            },
            "oncogenic": {
              "type": "string"
            },
            "hotspot": {
              "type": "boolean"
            },
            "highestFdaLevel": {
              "type": "string"
            },
            "otherSignificantSensitiveLevels": {
              "type": "array"
            },
            "otherSignificantResistanceLevels": {
              "type": "array"
            },
            "exon": {
              "type": "null"
            },
            "vus": {
              "type": "boolean"
            }
          },
          "required": [
            "query"
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

**执行与来源：** `network`；上游接口名 `OncoKB_annotate_variant`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("onco_kb_annotate_variant", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.progenetix_cnv_search` {#dnakit.databases.progenetix_cnv_search}

Search for cancer biosamples with copy number variations (CNVs) in a specific genomic region. Combines disease filter with genomic coordinates to find samples with amplifications or deletions in genes of interest. Uses RefSeq chromosome accessions for coordinates. Example: searching for EGFR amplifications (refseq:NC_000007.14:55019017-55211628, DUP) in breast cancer (NCIT:C4017) returns 9 matching biosamples.

```text
dnakit.databases.progenetix_cnv_search(*, reference_name, start, end, variant_type=None, filters=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `reference_name` | string | 必填（输入 Schema） | RefSeq chromosome accession. Examples: &#x27;refseq:NC_000007.14&#x27; (chr7/GRCh38), &#x27;refseq:NC_000017.11&#x27; (chr17/GRCh38), &#x27;refseq:NC_000001.11&#x27; (chr1/GRCh38). |
| `start` | integer | 必填（输入 Schema） | Start position (1-based, GRCh38). Example: 55019017 for EGFR start. |
| `end` | integer | 必填（输入 Schema） | End position (1-based, GRCh38). Example: 55211628 for EGFR end. |
| `variant_type` | string | None | CNV type: &#x27;DUP&#x27; for amplification/duplication, &#x27;DEL&#x27; for deletion. Leave empty for any CNV type. {&quot;default&quot;: &quot;&quot;} |
| `filters` | string | None | Optional NCIt ontology code to filter by cancer type. Example: &#x27;NCIT:C4017&#x27; for breast cancer. {&quot;default&quot;: &quot;&quot;} |
| `limit` | integer | None | Maximum number of biosamples to return (default: 10). {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "reference_name": {
          "type": "string",
          "description": "RefSeq chromosome accession. Examples: 'refseq:NC_000007.14' (chr7/GRCh38), 'refseq:NC_000017.11' (chr17/GRCh38), 'refseq:NC_000001.11' (chr1/GRCh38)."
        },
        "start": {
          "type": "integer",
          "description": "Start position (1-based, GRCh38). Example: 55019017 for EGFR start."
        },
        "end": {
          "type": "integer",
          "description": "End position (1-based, GRCh38). Example: 55211628 for EGFR end."
        },
        "variant_type": {
          "type": "string",
          "description": "CNV type: 'DUP' for amplification/duplication, 'DEL' for deletion. Leave empty for any CNV type.",
          "default": ""
        },
        "filters": {
          "type": "string",
          "description": "Optional NCIt ontology code to filter by cancer type. Example: 'NCIT:C4017' for breast cancer.",
          "default": ""
        },
        "limit": {
          "type": "integer",
          "default": 10,
          "minimum": 1,
          "maximum": 100,
          "description": "Maximum number of biosamples to return (default: 10)."
        }
      },
      "required": [
        "end",
        "reference_name",
        "start"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "region": {
          "type": "string"
        },
        "variant_type": {
          "type": "string"
        },
        "filters": {
          "type": "string"
        },
        "total_count": {
          "type": "integer"
        },
        "returned_count": {
          "type": "integer"
        },
        "biosamples": {
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
              "biosample_status": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "histological_diagnosis": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "histological_diagnosis_id": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "external_references": {
                "type": "array",
                "items": {
                  "type": "object"
                }
              }
            }
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `Progenetix_cnv_search`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("progenetix_cnv_search", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.progenetix_list_cohorts` {#dnakit.databases.progenetix_list_cohorts}

List available cohorts in the Progenetix cancer CNV database. Returns cohort IDs, names, sizes, and data types. Progenetix organizes samples into cohorts like TCGA, DIPG studies, and curated cancer sets. Example: returns 6 cohorts including TCGA (517 samples) and DIPG (16 samples).

```text
dnakit.databases.progenetix_list_cohorts(*, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `limit` | integer | None | Maximum number of cohorts to return (default: 10). {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 50} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "limit": {
          "type": "integer",
          "default": 10,
          "minimum": 1,
          "maximum": 50,
          "description": "Maximum number of cohorts to return (default: 10)."
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "total_cohorts": {
          "type": "integer"
        },
        "cohorts": {
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
              "cohort_size": {
                "type": [
                  "integer",
                  "null"
                ]
              },
              "cohort_type": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "data_types": {
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
    ```

**执行与来源：** `network`；上游接口名 `Progenetix_list_cohorts`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("progenetix_list_cohorts", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.progenetix_list_filtering_terms` {#dnakit.databases.progenetix_list_filtering_terms}

List available filtering terms (ontology codes) in the Progenetix cancer CNV database. Returns NCIt cancer type codes, EDAM operation codes, or other ontology terms with the number of matching samples. Useful for discovering which cancer types are available before querying biosamples or individuals. Example: prefixes='NCIT' returns 890 NCIt terms across all cancer types.

```text
dnakit.databases.progenetix_list_filtering_terms(*, prefixes=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `prefixes` | string | None | Ontology prefix to filter terms by. Use &#x27;NCIT&#x27; for NCI Thesaurus cancer types (default), &#x27;EDAM&#x27; for operation types, &#x27;HP&#x27; for HPO phenotypes. {&quot;default&quot;: &quot;NCIT&quot;} |
| `limit` | integer | None | Maximum number of terms to return (default: 25). {&quot;default&quot;: 25, &quot;minimum&quot;: 1, &quot;maximum&quot;: 200} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "prefixes": {
          "type": "string",
          "default": "NCIT",
          "description": "Ontology prefix to filter terms by. Use 'NCIT' for NCI Thesaurus cancer types (default), 'EDAM' for operation types, 'HP' for HPO phenotypes."
        },
        "limit": {
          "type": "integer",
          "default": 25,
          "minimum": 1,
          "maximum": 200,
          "description": "Maximum number of terms to return (default: 25)."
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "prefixes": {
          "type": "string"
        },
        "total_terms": {
          "type": "integer"
        },
        "terms": {
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
              "label": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "count": {
                "type": [
                  "integer",
                  "null"
                ]
              },
              "type": {
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
    ```

**执行与来源：** `network`；上游接口名 `Progenetix_list_filtering_terms`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("progenetix_list_filtering_terms", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.progenetix_search_biosamples` {#dnakit.databases.progenetix_search_biosamples}

Search cancer tumor biosamples in the Progenetix database by NCIt disease ontology code. Progenetix contains genome-wide copy number variation (CNV) profiles from over 100,000 tumor samples. Returns biosample IDs, histological diagnosis, pathological stage, and external references. Uses the GA4GH Beacon v2 protocol. Example: NCIT:C3058 returns 28 glioblastoma biosamples with histological classification and study references.

```text
dnakit.databases.progenetix_search_biosamples(*, filters, limit=None, skip=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `filters` | string | 必填（输入 Schema） | NCIt ontology code to filter biosamples by cancer type. Examples: &#x27;NCIT:C3058&#x27; (Glioblastoma), &#x27;NCIT:C4017&#x27; (Breast Carcinoma), &#x27;NCIT:C9145&#x27; (Acute Myeloid Leukemia), &#x27;NCIT:C3224&#x27; (Melanoma), &#x27;NCIT:C4872&#x27; (Lung Adenocarcinoma), &#x27;NCIT:C3167&#x27; (Neuroblastoma). |
| `limit` | integer | None | Maximum number of biosamples to return (default: 10). {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `skip` | integer | None | Number of results to skip for pagination (default: 0). {&quot;default&quot;: 0} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "filters": {
          "type": "string",
          "description": "NCIt ontology code to filter biosamples by cancer type. Examples: 'NCIT:C3058' (Glioblastoma), 'NCIT:C4017' (Breast Carcinoma), 'NCIT:C9145' (Acute Myeloid Leukemia), 'NCIT:C3224' (Melanoma), 'NCIT:C4872' (Lung Adenocarcinoma), 'NCIT:C3167' (Neuroblastoma)."
        },
        "limit": {
          "type": "integer",
          "default": 10,
          "minimum": 1,
          "maximum": 100,
          "description": "Maximum number of biosamples to return (default: 10)."
        },
        "skip": {
          "type": "integer",
          "default": 0,
          "description": "Number of results to skip for pagination (default: 0)."
        }
      },
      "required": [
        "filters"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "filters": {
          "type": "string"
        },
        "total_count": {
          "type": "integer"
        },
        "returned_count": {
          "type": "integer"
        },
        "biosamples": {
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
              "biosample_status": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "histological_diagnosis": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "histological_diagnosis_id": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "pathological_stage": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "pathological_tnm": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "sample_origin_type": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "external_references": {
                "type": "array",
                "items": {
                  "type": "object"
                }
              }
            }
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `Progenetix_search_biosamples`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("progenetix_search_biosamples", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.splice_ai_get_max_delta` {#dnakit.databases.splice_ai_get_max_delta}

Get simplified maximum SpliceAI delta score and pathogenicity interpretation for a variant. Returns single max score across all splice effects (AG/AL/DG/DL) with interpretation: benign (<0.2), low (0.2-0.5), moderate (0.5-0.8), high (≥0.8). Best for quick variant triage.

```text
dnakit.databases.splice_ai_get_max_delta(*, variant, genome=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant` | string | 必填（输入 Schema） | Variant in chr-pos-ref-alt format (e.g., chr8-140300616-T-G) |
| `genome` | string | None | Genome build: 37 or 38 (default: 38) {&quot;enum&quot;: &#91;&quot;37&quot;, &quot;38&quot;&#93;, &quot;default&quot;: &quot;38&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "variant": {
          "type": "string",
          "description": "Variant in chr-pos-ref-alt format (e.g., chr8-140300616-T-G)"
        },
        "genome": {
          "type": "string",
          "enum": [
            "37",
            "38"
          ],
          "description": "Genome build: 37 or 38 (default: 38)",
          "default": "38"
        }
      },
      "required": [
        "variant"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "variant": {
              "type": "string"
            },
            "genome": {
              "type": "string"
            },
            "max_delta_score": {
              "type": "number"
            },
            "interpretation": {
              "type": "string"
            },
            "pathogenicity_threshold": {
              "type": "string"
            }
          }
        },
        "source": {
          "type": "string"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `SpliceAI_get_max_delta`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("splice_ai_get_max_delta", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.splice_ai_predict_pangolin` {#dnakit.databases.splice_ai_predict_pangolin}

Predict splice-altering effects using Pangolin model (alternative to SpliceAI). Pangolin uses a different architecture and may provide complementary predictions. Returns splice effect scores per transcript. Useful for validation or when SpliceAI gives borderline results.

```text
dnakit.databases.splice_ai_predict_pangolin(*, variant, genome=None, distance=None, mask=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant` | string | 必填（输入 Schema） | Variant in chr-pos-ref-alt format (e.g., chr8-140300616-T-G) |
| `genome` | string | None | Genome build: 37 or 38 (default: 38) {&quot;enum&quot;: &#91;&quot;37&quot;, &quot;38&quot;&#93;, &quot;default&quot;: &quot;38&quot;} |
| `distance` | integer | None | Distance parameter for model (default: 1000) {&quot;default&quot;: 1000} |
| `mask` | boolean | None | Use masked scores {&quot;default&quot;: false} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "variant": {
          "type": "string",
          "description": "Variant in chr-pos-ref-alt format (e.g., chr8-140300616-T-G)"
        },
        "genome": {
          "type": "string",
          "enum": [
            "37",
            "38"
          ],
          "description": "Genome build: 37 or 38 (default: 38)",
          "default": "38"
        },
        "distance": {
          "type": "integer",
          "description": "Distance parameter for model (default: 1000)",
          "default": 1000
        },
        "mask": {
          "type": "boolean",
          "description": "Use masked scores",
          "default": false
        }
      },
      "required": [
        "variant"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "variant": {
              "type": "string"
            },
            "genome": {
              "type": "string"
            },
            "scores": {
              "type": "array"
            }
          }
        },
        "source": {
          "type": "string"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `SpliceAI_predict_pangolin`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("splice_ai_predict_pangolin", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.splice_ai_predict_splice` {#dnakit.databases.splice_ai_predict_splice}

Predict splice-altering effects using SpliceAI deep learning model. Returns delta scores for acceptor gain (DS_AG), acceptor loss (DS_AL), donor gain (DS_DG), donor loss (DS_DL). Scores ≥0.2 suggest splice impact, ≥0.5 moderate, ≥0.8 high pathogenicity. Critical for intronic/exonic variant interpretation.

```text
dnakit.databases.splice_ai_predict_splice(*, variant, genome=None, distance=None, mask=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant` | string | 必填（输入 Schema） | Variant in chr-pos-ref-alt format (e.g., chr8-140300616-T-G) or colon-separated |
| `genome` | string | None | Genome build: 37 (GRCh37/hg19) or 38 (GRCh38/hg38). Default: 38 {&quot;enum&quot;: &#91;&quot;37&quot;, &quot;38&quot;&#93;, &quot;default&quot;: &quot;38&quot;} |
| `distance` | integer | None | Distance parameter for model (default: 50). Larger values capture more distal effects. {&quot;default&quot;: 50} |
| `mask` | boolean | None | Use masked scores (recommended for variant interpretation). Raw scores better for alternative splicing analysis. {&quot;default&quot;: false} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "variant": {
          "type": "string",
          "description": "Variant in chr-pos-ref-alt format (e.g., chr8-140300616-T-G) or colon-separated"
        },
        "genome": {
          "type": "string",
          "enum": [
            "37",
            "38"
          ],
          "description": "Genome build: 37 (GRCh37/hg19) or 38 (GRCh38/hg38). Default: 38",
          "default": "38"
        },
        "distance": {
          "type": "integer",
          "description": "Distance parameter for model (default: 50). Larger values capture more distal effects.",
          "default": 50
        },
        "mask": {
          "type": "boolean",
          "description": "Use masked scores (recommended for variant interpretation). Raw scores better for alternative splicing analysis.",
          "default": false
        }
      },
      "required": [
        "variant"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "data": {
          "type": "object",
          "properties": {
            "variant": {
              "type": "string"
            },
            "genome": {
              "type": "string"
            },
            "scores": {
              "type": "array"
            },
            "max_delta_score": {
              "type": "number"
            },
            "interpretation": {
              "type": "string"
            }
          }
        },
        "source": {
          "type": "string"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `SpliceAI_predict_splice`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("splice_ai_predict_splice", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.tpmi_phewas_by_variant` {#dnakit.databases.tpmi_phewas_by_variant}

Phenome-wide association study (PheWAS) lookup for a single variant in the Taiwan Precision Medicine Initiative (TPMI) PheWeb, a large Taiwanese (Han Chinese) biobank. Returns every phenotype associated with the variant, sorted by p-value. Accepts an rsID (auto-resolved to GRCh38 via Ensembl) or an explicit chr:pos:ref:alt variant in GRCh38. Use for cross-ancestry replication in an East Asian (Han Chinese) population.

```text
dnakit.databases.tpmi_phewas_by_variant(*, rsid=None, variant=None, limit=None, max_pval=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `rsid` | string | None | dbSNP rsID, e.g. &#x27;rs7903146&#x27;. Auto-resolved to GRCh38 coordinates via Ensembl. Provide this OR &#x27;variant&#x27;. |
| `variant` | string | None | Variant in chr:pos:ref:alt or chr-pos-ref-alt format, in GRCh38, e.g. &#x27;10:112998590:C:T&#x27;. Provide this OR &#x27;rsid&#x27;. A bare rsID is also accepted here. |
| `limit` | integer | None | Maximum number of phenotype associations to return (sorted by ascending p-value). {&quot;default&quot;: 20} |
| `max_pval` | number &#124; null | None | Optional: only return associations with p-value &lt;= this threshold (e.g. 5e-8 for genome-wide significance). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "rsid": {
          "type": "string",
          "description": "dbSNP rsID, e.g. 'rs7903146'. Auto-resolved to GRCh38 coordinates via Ensembl. Provide this OR 'variant'."
        },
        "variant": {
          "type": "string",
          "description": "Variant in chr:pos:ref:alt or chr-pos-ref-alt format, in GRCh38, e.g. '10:112998590:C:T'. Provide this OR 'rsid'. A bare rsID is also accepted here."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of phenotype associations to return (sorted by ascending p-value).",
          "default": 20
        },
        "max_pval": {
          "type": [
            "number",
            "null"
          ],
          "description": "Optional: only return associations with p-value <= this threshold (e.g. 5e-8 for genome-wide significance)."
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "variant": {
              "type": "string"
            },
            "rsids": {
              "type": [
                "string",
                "null"
              ]
            },
            "nearest_genes": {
              "type": [
                "string",
                "null"
              ]
            },
            "associations": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "phenocode": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "phenostring": {
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
                  "pval": {
                    "type": [
                      "number",
                      "null"
                    ]
                  },
                  "beta": {
                    "type": [
                      "number",
                      "null"
                    ]
                  },
                  "sebeta": {
                    "type": [
                      "number",
                      "null"
                    ]
                  },
                  "af": {
                    "type": [
                      "number",
                      "null"
                    ]
                  },
                  "num_cases": {
                    "type": [
                      "number",
                      "null"
                    ]
                  },
                  "num_controls": {
                    "type": [
                      "number",
                      "null"
                    ]
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

**执行与来源：** `network`；上游接口名 `TPMI_phewas_by_variant`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("tpmi_phewas_by_variant", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ukbtopmed_phewas_by_variant` {#dnakit.databases.ukbtopmed_phewas_by_variant}

Phenome-wide association study (PheWAS) lookup for a single variant in the UKB-TOPMed PheWeb (UK Biobank, ~400K European-ancestry samples, 1,400+ phecode phenotypes). Returns every phenotype associated with the variant, sorted by p-value. Accepts an rsID (auto-resolved to GRCh38 via Ensembl) or an explicit chr:pos:ref:alt variant in GRCh38. Use for variant pleiotropy assessment and replication in a European population.

```text
dnakit.databases.ukbtopmed_phewas_by_variant(*, rsid=None, variant=None, limit=None, max_pval=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `rsid` | string | None | dbSNP rsID, e.g. &#x27;rs7903146&#x27;. Auto-resolved to GRCh38 coordinates via Ensembl. Provide this OR &#x27;variant&#x27;. |
| `variant` | string | None | Variant in chr:pos:ref:alt or chr-pos-ref-alt format, in GRCh38, e.g. &#x27;10:112998590:C:T&#x27;. Provide this OR &#x27;rsid&#x27;. A bare rsID is also accepted here. |
| `limit` | integer | None | Maximum number of phenotype associations to return (sorted by ascending p-value). {&quot;default&quot;: 20} |
| `max_pval` | number &#124; null | None | Optional: only return associations with p-value &lt;= this threshold (e.g. 5e-8 for genome-wide significance). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "rsid": {
          "type": "string",
          "description": "dbSNP rsID, e.g. 'rs7903146'. Auto-resolved to GRCh38 coordinates via Ensembl. Provide this OR 'variant'."
        },
        "variant": {
          "type": "string",
          "description": "Variant in chr:pos:ref:alt or chr-pos-ref-alt format, in GRCh38, e.g. '10:112998590:C:T'. Provide this OR 'rsid'. A bare rsID is also accepted here."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of phenotype associations to return (sorted by ascending p-value).",
          "default": 20
        },
        "max_pval": {
          "type": [
            "number",
            "null"
          ],
          "description": "Optional: only return associations with p-value <= this threshold (e.g. 5e-8 for genome-wide significance)."
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "variant": {
              "type": "string"
            },
            "rsids": {
              "type": [
                "string",
                "null"
              ]
            },
            "nearest_genes": {
              "type": [
                "string",
                "null"
              ]
            },
            "associations": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "phenocode": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "phenostring": {
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
                  "pval": {
                    "type": [
                      "number",
                      "null"
                    ]
                  },
                  "beta": {
                    "type": [
                      "number",
                      "null"
                    ]
                  },
                  "sebeta": {
                    "type": [
                      "number",
                      "null"
                    ]
                  },
                  "af": {
                    "type": [
                      "number",
                      "null"
                    ]
                  },
                  "num_cases": {
                    "type": [
                      "number",
                      "null"
                    ]
                  },
                  "num_controls": {
                    "type": [
                      "number",
                      "null"
                    ]
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

**执行与来源：** `network`；上游接口名 `UKBTOPMed_phewas_by_variant`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ukbtopmed_phewas_by_variant", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.variant_validator_gene2transcripts` {#dnakit.databases.variant_validator_gene2transcripts}

Get all transcripts for a gene using VariantValidator. Returns MANE Select and MANE Plus Clinical transcript annotations, Ensembl and RefSeq transcript identifiers, genomic coordinates, and cross-references. Useful for identifying the correct transcript for variant annotation. Supports filter for MANE transcripts only or all transcripts.

```text
dnakit.databases.variant_validator_gene2transcripts(*, gene_symbol, gene=None, gene_name=None, transcript_set=None, genome_build=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 必填（输入 Schema） | HGNC gene symbol (e.g., &#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;EGFR&#x27;). Aliases: gene, gene_name. |
| `gene` | string | None | Alias for gene_symbol. |
| `gene_name` | string | None | Alias for gene_symbol. |
| `transcript_set` | string | None | Transcript filter: &#x27;mane&#x27; for MANE Select/Plus Clinical only, &#x27;refseq&#x27; for RefSeq transcripts, &#x27;ensembl&#x27; for Ensembl, &#x27;all&#x27; for everything {&quot;default&quot;: &quot;mane&quot;} |
| `genome_build` | string | None | Reference genome assembly: &#x27;GRCh37&#x27; or &#x27;GRCh38&#x27; (default: GRCh38) {&quot;default&quot;: &quot;GRCh38&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene_symbol": {
          "type": "string",
          "description": "HGNC gene symbol (e.g., 'TP53', 'BRCA1', 'EGFR'). Aliases: gene, gene_name."
        },
        "gene": {
          "type": "string",
          "description": "Alias for gene_symbol."
        },
        "gene_name": {
          "type": "string",
          "description": "Alias for gene_symbol."
        },
        "transcript_set": {
          "type": "string",
          "description": "Transcript filter: 'mane' for MANE Select/Plus Clinical only, 'refseq' for RefSeq transcripts, 'ensembl' for Ensembl, 'all' for everything",
          "default": "mane"
        },
        "genome_build": {
          "type": "string",
          "description": "Reference genome assembly: 'GRCh37' or 'GRCh38' (default: GRCh38)",
          "default": "GRCh38"
        }
      },
      "required": [
        "gene_symbol"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "current_name": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "current_symbol": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "hgnc": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "previous_symbol": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "requested_symbol": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "transcripts": {
                "type": "array",
                "items": {
                  "type": "object",
                  "properties": {
                    "reference": {
                      "type": [
                        "string",
                        "null"
                      ],
                      "description": "Transcript accession (e.g., NM_000546.6 or ENST00000269305.9)"
                    },
                    "description": {
                      "type": [
                        "string",
                        "null"
                      ]
                    },
                    "coding_start": {
                      "type": [
                        "integer",
                        "null"
                      ]
                    },
                    "coding_end": {
                      "type": [
                        "integer",
                        "null"
                      ]
                    },
                    "length": {
                      "type": [
                        "integer",
                        "null"
                      ]
                    },
                    "annotations": {
                      "type": "object",
                      "properties": {
                        "chromosome": {
                          "type": [
                            "string",
                            "null"
                          ]
                        },
                        "mane_select": {
                          "type": [
                            "boolean",
                            "null"
                          ]
                        },
                        "mane_plus_clinical": {
                          "type": [
                            "boolean",
                            "null"
                          ]
                        },
                        "ensembl_select": {
                          "type": [
                            "boolean",
                            "null"
                          ]
                        },
                        "db_xref": {
                          "type": [
                            "object",
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

**执行与来源：** `network`；上游接口名 `VariantValidator_gene2transcripts`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("variant_validator_gene2transcripts", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.variant_validator_validate_variant` {#dnakit.databases.variant_validator_validate_variant}

Validate and convert genetic variant descriptions using VariantValidator. Takes HGVS-formatted variant notation and returns validated variant descriptions with genomic context, transcript context, protein consequence, and cross-references. Validates variant nomenclature against current RefSeq sequences and LOVD syntax checker. Supports substitutions, deletions, insertions, duplications, inversions. Reference assemblies: GRCh37 (hg19) or GRCh38 (hg38).

```text
dnakit.databases.variant_validator_validate_variant(*, genome_build, variant_description, select_transcripts, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `genome_build` | string | 必填（输入 Schema） | Reference genome assembly: &#x27;GRCh37&#x27; (hg19) or &#x27;GRCh38&#x27; (hg38) |
| `variant_description` | string | 必填（输入 Schema） | HGVS variant description (e.g., &#x27;NM_007294.4:c.5266dup&#x27; for BRCA1 c.5266dupC, &#x27;NM_000179.3:c.943C&gt;T&#x27;, &#x27;chr17:g.43092919del&#x27;) |
| `select_transcripts` | string | 必填（输入 Schema） | Transcript to validate against (e.g., &#x27;NM_007294.4&#x27;). Use &#x27;all&#x27; to get all transcripts for a genomic variant. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "genome_build": {
          "type": "string",
          "description": "Reference genome assembly: 'GRCh37' (hg19) or 'GRCh38' (hg38)"
        },
        "variant_description": {
          "type": "string",
          "description": "HGVS variant description (e.g., 'NM_007294.4:c.5266dup' for BRCA1 c.5266dupC, 'NM_000179.3:c.943C>T', 'chr17:g.43092919del')"
        },
        "select_transcripts": {
          "type": "string",
          "description": "Transcript to validate against (e.g., 'NM_007294.4'). Use 'all' to get all transcripts for a genomic variant."
        }
      },
      "required": [
        "genome_build",
        "select_transcripts",
        "variant_description"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "description": "Validated variant information. Contains 'flag' field, 'metadata' field, and dynamic keys for each validated variant.",
          "properties": {
            "flag": {
              "type": [
                "string",
                "null"
              ],
              "description": "Validation status: 'gene_variant', 'warning', 'error'"
            },
            "metadata": {
              "type": [
                "object",
                "null"
              ],
              "description": "VariantValidator version information",
              "properties": {
                "variantvalidator_version": {
                  "type": "string"
                },
                "variantvalidator_hgvs_version": {
                  "type": "string"
                },
                "vvdb_version": {
                  "type": "string"
                }
              }
            }
          },
          "additionalProperties": {
            "type": "object",
            "description": "HGVS variant key with validated transcript/genomic coordinates, protein consequence, and warnings",
            "properties": {
              "hgvs_transcript_variant": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Validated HGVS transcript notation"
              },
              "hgvs_refseqgene_variant": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "HGVS RefSeqGene notation"
              },
              "hgvs_predicted_protein_consequence": {
                "type": [
                  "object",
                  "null"
                ],
                "properties": {
                  "tlr": {
                    "type": [
                      "string",
                      "null"
                    ],
                    "description": "Three-letter amino acid code"
                  },
                  "slr": {
                    "type": [
                      "string",
                      "null"
                    ],
                    "description": "Single-letter amino acid code"
                  }
                }
              },
              "primary_assembly_loci": {
                "type": [
                  "object",
                  "null"
                ],
                "description": "Genomic coordinates in GRCh37/GRCh38",
                "properties": {
                  "grch38": {
                    "type": [
                      "object",
                      "null"
                    ]
                  },
                  "grch37": {
                    "type": [
                      "object",
                      "null"
                    ]
                  }
                }
              },
              "gene_ids": {
                "type": [
                  "object",
                  "null"
                ],
                "properties": {
                  "hgnc_id": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "entrez_gene_id": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "ncbi_id": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "ensembl_gene_id": {
                    "type": [
                      "string",
                      "null"
                    ]
                  }
                }
              },
              "gene_symbol": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "validation_warnings": {
                "type": [
                  "array",
                  "null"
                ]
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

**执行与来源：** `network`；上游接口名 `VariantValidator_validate_variant`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("variant_validator_validate_variant", verbose=False)
print(info["input_schema"])
```

