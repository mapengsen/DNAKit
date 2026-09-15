# 通路、互作、蛋白与结构

本页收录 **64 个新增公开函数**。 [完整索引](index.md) · [调用约定与完整示例](usage.md)

每个函数均列出当前 Python 签名、参数、返回值和已登记的执行要求。提供方说明及字段描述保留源码原文。`None` 默认值在生成的数据库包装中通常表示不向服务传该参数；实际必填项与组合约束以输入 Schema 为准。坐标系和基因组版本按各接口说明使用。

## 函数目录

- [`dnakit.databases.alphafold_get_prediction`](#dnakit.databases.alphafold_get_prediction)
- [`dnakit.databases.bio_grid_get_interactions`](#dnakit.databases.bio_grid_get_interactions)
- [`dnakit.databases.blast_protein_search`](#dnakit.databases.blast_protein_search)
- [`dnakit.databases.bvbrc_search_epitopes`](#dnakit.databases.bvbrc_search_epitopes)
- [`dnakit.databases.enrichr_enrich`](#dnakit.databases.enrichr_enrich)
- [`dnakit.databases.enrichr_gene_enrichment_analysis`](#dnakit.databases.enrichr_gene_enrichment_analysis)
- [`dnakit.databases.go_get_annotations_for_gene`](#dnakit.databases.go_get_annotations_for_gene)
- [`dnakit.databases.go_search_terms`](#dnakit.databases.go_search_terms)
- [`dnakit.databases.humanbase_ppi_analysis`](#dnakit.databases.humanbase_ppi_analysis)
- [`dnakit.databases.iedb_get_epitope_antigens`](#dnakit.databases.iedb_get_epitope_antigens)
- [`dnakit.databases.iedb_get_epitope_mhc`](#dnakit.databases.iedb_get_epitope_mhc)
- [`dnakit.databases.iedb_search_bcell`](#dnakit.databases.iedb_search_bcell)
- [`dnakit.databases.iedb_search_epitopes`](#dnakit.databases.iedb_search_epitopes)
- [`dnakit.databases.iedb_search_mhc`](#dnakit.databases.iedb_search_mhc)
- [`dnakit.databases.iedb_search_tcell_assays`](#dnakit.databases.iedb_search_tcell_assays)
- [`dnakit.databases.imgt_get_gene_info`](#dnakit.databases.imgt_get_gene_info)
- [`dnakit.databases.imgt_search_genes`](#dnakit.databases.imgt_search_genes)
- [`dnakit.databases.imm_port_search_studies`](#dnakit.databases.imm_port_search_studies)
- [`dnakit.databases.intact_get_complex_details`](#dnakit.databases.intact_get_complex_details)
- [`dnakit.databases.intact_get_interaction_network`](#dnakit.databases.intact_get_interaction_network)
- [`dnakit.databases.intact_get_interactions`](#dnakit.databases.intact_get_interactions)
- [`dnakit.databases.intact_search_interactions`](#dnakit.databases.intact_search_interactions)
- [`dnakit.databases.inter_pro_get_entries_for_protein`](#dnakit.databases.inter_pro_get_entries_for_protein)
- [`dnakit.databases.inter_pro_get_protein_domains`](#dnakit.databases.inter_pro_get_protein_domains)
- [`dnakit.databases.kegg_find_genes`](#dnakit.databases.kegg_find_genes)
- [`dnakit.databases.kegg_get_gene_info`](#dnakit.databases.kegg_get_gene_info)
- [`dnakit.databases.kegg_get_gene_pathways`](#dnakit.databases.kegg_get_gene_pathways)
- [`dnakit.databases.kegg_get_pathway_genes`](#dnakit.databases.kegg_get_pathway_genes)
- [`dnakit.databases.kegg_get_pathway_info`](#dnakit.databases.kegg_get_pathway_info)
- [`dnakit.databases.kegg_search_pathway`](#dnakit.databases.kegg_search_pathway)
- [`dnakit.databases.ncbiprotein_get_summary`](#dnakit.databases.ncbiprotein_get_summary)
- [`dnakit.databases.panther_enrichment`](#dnakit.databases.panther_enrichment)
- [`dnakit.databases.panther_ortholog`](#dnakit.databases.panther_ortholog)
- [`dnakit.databases.pdb_get_structure`](#dnakit.databases.pdb_get_structure)
- [`dnakit.databases.pdbe_get_uniprot_mappings`](#dnakit.databases.pdbe_get_uniprot_mappings)
- [`dnakit.databases.pdbe_sifts_get_best_structures`](#dnakit.databases.pdbe_sifts_get_best_structures)
- [`dnakit.databases.pfam_get_protein_annotations`](#dnakit.databases.pfam_get_protein_annotations)
- [`dnakit.databases.plant_reactome_get_pathway`](#dnakit.databases.plant_reactome_get_pathway)
- [`dnakit.databases.plant_reactome_list_species`](#dnakit.databases.plant_reactome_list_species)
- [`dnakit.databases.plant_reactome_search_pathways`](#dnakit.databases.plant_reactome_search_pathways)
- [`dnakit.databases.proteins_api_get_features`](#dnakit.databases.proteins_api_get_features)
- [`dnakit.databases.proteins_api_search`](#dnakit.databases.proteins_api_search)
- [`dnakit.databases.quick_go_annotations_by_gene`](#dnakit.databases.quick_go_annotations_by_gene)
- [`dnakit.databases.rcsbdata_get_entry`](#dnakit.databases.rcsbdata_get_entry)
- [`dnakit.databases.reactome_analysis_pathway_enrichment`](#dnakit.databases.reactome_analysis_pathway_enrichment)
- [`dnakit.databases.reactome_content_search`](#dnakit.databases.reactome_content_search)
- [`dnakit.databases.reactome_get_participants`](#dnakit.databases.reactome_get_participants)
- [`dnakit.databases.reactome_get_pathway`](#dnakit.databases.reactome_get_pathway)
- [`dnakit.databases.reactome_get_pathway_reactions`](#dnakit.databases.reactome_get_pathway_reactions)
- [`dnakit.databases.reactome_map_uniprot_to_pathways`](#dnakit.databases.reactome_map_uniprot_to_pathways)
- [`dnakit.databases.string_functional_enrichment`](#dnakit.databases.string_functional_enrichment)
- [`dnakit.databases.string_get_functional_annotations`](#dnakit.databases.string_get_functional_annotations)
- [`dnakit.databases.string_get_functional_enrichment`](#dnakit.databases.string_get_functional_enrichment)
- [`dnakit.databases.string_get_interaction_partners`](#dnakit.databases.string_get_interaction_partners)
- [`dnakit.databases.string_get_network`](#dnakit.databases.string_get_network)
- [`dnakit.databases.string_ppi_enrichment`](#dnakit.databases.string_ppi_enrichment)
- [`dnakit.databases.structure_annotate_per_residue`](#dnakit.databases.structure_annotate_per_residue)
- [`dnakit.databases.uni_prot_get_disease_variants_by_accession`](#dnakit.databases.uni_prot_get_disease_variants_by_accession)
- [`dnakit.databases.uni_prot_get_entry_by_accession`](#dnakit.databases.uni_prot_get_entry_by_accession)
- [`dnakit.databases.uni_prot_get_features_by_accession`](#dnakit.databases.uni_prot_get_features_by_accession)
- [`dnakit.databases.uni_prot_get_function_by_accession`](#dnakit.databases.uni_prot_get_function_by_accession)
- [`dnakit.databases.uni_prot_get_sequence_by_accession`](#dnakit.databases.uni_prot_get_sequence_by_accession)
- [`dnakit.databases.uni_prot_search`](#dnakit.databases.uni_prot_search)
- [`dnakit.databases.wiki_pathways_search`](#dnakit.databases.wiki_pathways_search)

## `dnakit.databases.alphafold_get_prediction` {#dnakit.databases.alphafold_get_prediction}

Retrieve full AlphaFold 3D structure predictions for a given protein. IMPORTANT: The qualifier must be a UniProt ACCESSION (e.g., 'P69905' for HBA_HUMAN). Do NOT use UniProt entry names like 'HBA_HUMAN' - they will cause API errors. To find UniProt accession from a gene/protein name, use `UniProt_search` (e.g., query='gene:HBA' organism='human') or `UniProt_id_mapping` for ID conversion. Returns residue-level metadata including sequence, per-residue confidence scores (pLDDT), and structure download links (PDB, CIF, PAE). If you already have the accession and want UniProt details, call `UniProt_get_entry_by_accession`. For a quick overview, use `alphafold_get_summary`. For mutation/variant impact, see `alphafold_get_annotations`.

```text
dnakit.databases.alphafold_get_prediction(*, qualifier=None, uniprot_id=None, uniprot_accession=None, accession=None, sequence_checksum=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `qualifier` | string | None | UniProt ACCESSION (e.g., &#x27;P69905&#x27;). Do NOT use entry names like &#x27;HBA_HUMAN&#x27;. Aliases: uniprot_id, uniprot_accession. |
| `uniprot_id` | string | None | Alias for qualifier: UniProt accession (e.g., &#x27;P69905&#x27;). |
| `uniprot_accession` | string | None | Alias for qualifier: UniProt accession (e.g., &#x27;P69905&#x27;). |
| `accession` | string | None | Alias for qualifier. UniProt accession (e.g., &#x27;P04637&#x27;) -- the same parameter name the UniProt_*_by_accession tools use. |
| `sequence_checksum` | string | None | Optional CRC64 checksum of the UniProt sequence. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "qualifier": {
          "type": "string",
          "description": "UniProt ACCESSION (e.g., 'P69905'). Do NOT use entry names like 'HBA_HUMAN'. Aliases: uniprot_id, uniprot_accession."
        },
        "uniprot_id": {
          "type": "string",
          "description": "Alias for qualifier: UniProt accession (e.g., 'P69905')."
        },
        "uniprot_accession": {
          "type": "string",
          "description": "Alias for qualifier: UniProt accession (e.g., 'P69905')."
        },
        "accession": {
          "type": "string",
          "description": "Alias for qualifier. UniProt accession (e.g., 'P04637') -- the same parameter name the UniProt_*_by_accession tools use."
        },
        "sequence_checksum": {
          "type": "string",
          "description": "Optional CRC64 checksum of the UniProt sequence."
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
          "items": {
            "type": "object",
            "properties": {
              "toolUsed": {
                "type": "string"
              },
              "providerId": {
                "type": "string"
              },
              "entityType": {
                "type": "string"
              },
              "isUniProt": {
                "type": "boolean"
              },
              "modelEntityId": {
                "type": "string"
              },
              "modelCreatedDate": {
                "type": "string"
              },
              "sequenceVersionDate": {
                "type": "string"
              },
              "globalMetricValue": {
                "type": "number"
              },
              "fractionPlddtVeryLow": {
                "type": "number"
              },
              "fractionPlddtLow": {
                "type": "number"
              },
              "fractionPlddtConfident": {
                "type": "number"
              },
              "fractionPlddtVeryHigh": {
                "type": "number"
              },
              "latestVersion": {
                "type": "integer"
              },
              "allVersions": {
                "type": "array",
                "items": {
                  "type": "integer"
                }
              },
              "sequence": {
                "type": "string"
              },
              "sequenceStart": {
                "type": "integer"
              },
              "sequenceEnd": {
                "type": "integer"
              },
              "sequenceChecksum": {
                "type": "string"
              },
              "isUniProtReviewed": {
                "type": "boolean"
              },
              "gene": {
                "type": "string"
              },
              "uniprotAccession": {
                "type": "string"
              },
              "uniprotId": {
                "type": "string"
              },
              "uniprotDescription": {
                "type": "string"
              },
              "taxId": {
                "type": "integer"
              },
              "organismScientificName": {
                "type": "string"
              },
              "isUniProtReferenceProteome": {
                "type": "boolean"
              },
              "organismCommonNames": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "organismSynonyms": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "geneSynonyms": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "proteinFullNames": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "proteinShortNames": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "complexName": {
                "type": "string"
              },
              "stoichiometry": {
                "type": "number"
              },
              "ipTM": {
                "type": "number"
              },
              "ipSAE": {
                "type": "number"
              },
              "keywords": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "taxonomyLineage": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "functions": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "alternativeNames": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "catalyticActivities": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "bcifUrl": {
                "type": "string"
              },
              "cifUrl": {
                "type": "string"
              },
              "pdbUrl": {
                "type": "string"
              },
              "paeImageUrl": {
                "type": "string"
              },
              "paeDocUrl": {
                "type": "string"
              },
              "amAnnotationsUrl": {
                "type": "string"
              },
              "amAnnotationsHg19Url": {
                "type": "string"
              },
              "amAnnotationsHg38Url": {
                "type": "string"
              },
              "entryId": {
                "type": "string"
              },
              "isReviewed": {
                "type": "boolean"
              },
              "isReferenceProteome": {
                "type": "boolean"
              },
              "uniprotStart": {
                "type": "integer"
              },
              "uniprotEnd": {
                "type": "integer"
              },
              "uniprotSequence": {
                "type": "string"
              }
            }
          }
        },
        {
          "type": "object",
          "properties": {
            "error": {
              "type": "string"
            },
            "detail": {
              "type": "string"
            },
            "endpoint": {
              "type": "string"
            },
            "query": {
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

**执行与来源：** `network`；上游接口名 `alphafold_get_prediction`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("alphafold_get_prediction", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.bio_grid_get_interactions` {#dnakit.databases.bio_grid_get_interactions}

Query experimentally validated protein and genetic interactions from BioGRID (Biological General Repository for Interaction Datasets). Returns curated interactions from published studies with evidence methods, PubMed citations, and throughput information. **Prerequisites**: Requires BIOGRID_API_KEY environment variable (free academic API key at https://webservice.thebiogrid.org/). BioGRID contains 2.3M+ interactions from 80+ organisms, all experimentally validated (no predictions). Use for: finding experimentally proven interactions, getting literature evidence for interactions, validating predicted interactions, finding interaction methods used, accessing high-confidence curated data.

```text
dnakit.databases.bio_grid_get_interactions(*, gene_names=None, gene_name=None, organism=None, interaction_type=None, evidence_types=None, limit=None, throughput=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_names` | array | None | List of gene names or protein identifiers (e.g., &#91;&#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;MYC&#x27;&#93;). Accepts official gene symbols. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}, &quot;minItems&quot;: 1} |
| `gene_name` | string | None | Alias for gene_names: single gene symbol (e.g., &#x27;TP53&#x27;). Converted to gene_names list automatically. |
| `organism` | string | None | Organism name (e.g., &#x27;Homo sapiens&#x27;, &#x27;Mus musculus&#x27;) or NCBI taxonomy ID (e.g., &#x27;9606&#x27; for human, &#x27;10090&#x27; for mouse). Default: 9606 (human) {&quot;default&quot;: &quot;9606&quot;} |
| `interaction_type` | string | None | Type of interaction: &#x27;physical&#x27; (protein-protein), &#x27;genetic&#x27; (epistasis, synthetic lethality), &#x27;both&#x27; (all interactions) {&quot;enum&quot;: &#91;&quot;physical&quot;, &quot;genetic&quot;, &quot;both&quot;&#93;, &quot;default&quot;: &quot;both&quot;} |
| `evidence_types` | array | None | Filter by evidence types (e.g., &#91;&#x27;Affinity Capture-MS&#x27;, &#x27;Two-hybrid&#x27;&#93; for physical, &#91;&#x27;Synthetic Lethality&#x27;, &#x27;Dosage Rescue&#x27;&#93; for genetic). Leave empty for all evidence types. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `limit` | integer | None | Maximum number of interactions to return (default: 100, max: 10000) {&quot;minimum&quot;: 1, &quot;maximum&quot;: 10000, &quot;default&quot;: 100} |
| `throughput` | string | None | Filter by throughput: &#x27;low&#x27; (low-throughput studies), &#x27;high&#x27; (high-throughput screens), or leave empty for all {&quot;enum&quot;: &#91;&quot;low&quot;, &quot;high&quot;, &quot;&quot;&#93;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene_names": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "List of gene names or protein identifiers (e.g., ['TP53', 'BRCA1', 'MYC']). Accepts official gene symbols.",
          "minItems": 1
        },
        "gene_name": {
          "type": "string",
          "description": "Alias for gene_names: single gene symbol (e.g., 'TP53'). Converted to gene_names list automatically."
        },
        "organism": {
          "type": "string",
          "description": "Organism name (e.g., 'Homo sapiens', 'Mus musculus') or NCBI taxonomy ID (e.g., '9606' for human, '10090' for mouse). Default: 9606 (human)",
          "default": "9606"
        },
        "interaction_type": {
          "type": "string",
          "description": "Type of interaction: 'physical' (protein-protein), 'genetic' (epistasis, synthetic lethality), 'both' (all interactions)",
          "enum": [
            "physical",
            "genetic",
            "both"
          ],
          "default": "both"
        },
        "evidence_types": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "Filter by evidence types (e.g., ['Affinity Capture-MS', 'Two-hybrid'] for physical, ['Synthetic Lethality', 'Dosage Rescue'] for genetic). Leave empty for all evidence types."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of interactions to return (default: 100, max: 10000)",
          "minimum": 1,
          "maximum": 10000,
          "default": 100
        },
        "throughput": {
          "type": "string",
          "description": "Filter by throughput: 'low' (low-throughput studies), 'high' (high-throughput screens), or leave empty for all",
          "enum": [
            "low",
            "high",
            ""
          ]
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
          "description": "BioGRID API response keyed by interaction ID. Each value is an interaction object.",
          "additionalProperties": {
            "type": "object",
            "properties": {
              "BIOGRID_INTERACTION_ID": {
                "type": "integer"
              },
              "ENTREZ_GENE_A": {
                "type": "string"
              },
              "ENTREZ_GENE_B": {
                "type": "string"
              },
              "OFFICIAL_SYMBOL_A": {
                "type": "string"
              },
              "OFFICIAL_SYMBOL_B": {
                "type": "string"
              },
              "SYNONYMS_A": {
                "type": "string"
              },
              "SYNONYMS_B": {
                "type": "string"
              },
              "EXPERIMENTAL_SYSTEM": {
                "type": "string"
              },
              "EXPERIMENTAL_SYSTEM_TYPE": {
                "type": "string"
              },
              "PUBMED_ID": {
                "type": "integer"
              },
              "THROUGHPUT": {
                "type": "string"
              },
              "ORGANISM_A_ID": {
                "type": "integer"
              },
              "ORGANISM_B_ID": {
                "type": "integer"
              }
            }
          }
        },
        {
          "type": "object",
          "description": "Error response",
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

**执行与来源：** `network`；上游接口名 `BioGRID_get_interactions`。

**提供方凭据：** `BIOGRID_ACCESS_KEY`

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("bio_grid_get_interactions", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.blast_protein_search` {#dnakit.databases.blast_protein_search}

Search protein sequences using NCBI BLAST blastp against protein databases. Requires Biopython (pip install biopython). IMPORTANT: This tool uses NCBI's remote BLAST service which is inherently SLOW - expect 5-30 minutes per query depending on server load. This is normal NCBI behavior, not a bug. For faster results, use 'swissprot' or 'pdb' databases instead of 'nr'.

```text
dnakit.databases.blast_protein_search(*, sequence, database=None, expect=None, hitlist_size=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `sequence` | string | 必填（输入 Schema） | Protein sequence to search (minimum 10 amino acids) |
| `database` | string | None | Database to search. Options: &#x27;nr&#x27; (comprehensive but slow), &#x27;swissprot&#x27; (curated, faster), &#x27;pdb&#x27; (protein structures only) {&quot;default&quot;: &quot;nr&quot;} |
| `expect` | number | None | E-value threshold for reporting hits {&quot;default&quot;: 10.0} |
| `hitlist_size` | integer | None | Maximum number of hits to return {&quot;default&quot;: 50} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "sequence": {
          "type": "string",
          "description": "Protein sequence to search (minimum 10 amino acids)"
        },
        "database": {
          "type": "string",
          "default": "nr",
          "description": "Database to search. Options: 'nr' (comprehensive but slow), 'swissprot' (curated, faster), 'pdb' (protein structures only)"
        },
        "expect": {
          "type": "number",
          "default": 10.0,
          "description": "E-value threshold for reporting hits"
        },
        "hitlist_size": {
          "type": "integer",
          "default": 50,
          "description": "Maximum number of hits to return"
        }
      },
      "required": [
        "sequence"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "query_id": {
          "type": [
            "string",
            "null"
          ]
        },
        "query_length": {
          "type": [
            "integer",
            "null"
          ]
        },
        "database": {
          "type": [
            "string",
            "null"
          ]
        },
        "algorithm": {
          "type": [
            "string",
            "null"
          ]
        },
        "alignments": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "hit_id": {
                "type": "string"
              },
              "hit_def": {
                "type": "string"
              },
              "hit_length": {
                "type": "integer"
              },
              "hsps": {
                "type": "array",
                "items": {
                  "type": "object",
                  "properties": {
                    "score": {
                      "type": "number"
                    },
                    "bits": {
                      "type": "number"
                    },
                    "expect": {
                      "type": "number"
                    },
                    "identities": {
                      "type": "integer"
                    },
                    "positives": {
                      "type": [
                        "integer",
                        "null"
                      ]
                    },
                    "gaps": {
                      "type": "integer"
                    },
                    "align_length": {
                      "type": "integer"
                    },
                    "query_start": {
                      "type": "integer"
                    },
                    "query_end": {
                      "type": "integer"
                    },
                    "hit_start": {
                      "type": "integer"
                    },
                    "hit_end": {
                      "type": "integer"
                    },
                    "query": {
                      "type": "string"
                    },
                    "match": {
                      "type": "string"
                    },
                    "sbjct": {
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
    ```

**执行与来源：** `network`；上游接口名 `BLAST_protein_search`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("blast_protein_search", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.bvbrc_search_epitopes` {#dnakit.databases.bvbrc_search_epitopes}

Search for pathogen epitopes (B-cell and T-cell) in BV-BRC. Returns epitope sequences, types, protein associations, and assay results for vaccine and immunology research. Covers epitopes from IEDB for bacteria and viruses including SARS-CoV-2, influenza, tuberculosis, and other pathogens. Example: search taxon_id 2697049 for SARS-CoV-2 spike glycoprotein epitopes.

```text
dnakit.databases.bvbrc_search_epitopes(*, taxon_id=None, protein_name=None, epitope_type=None, organism=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `taxon_id` | string &#124; null | None | NCBI Taxonomy ID for the pathogen organism. Examples: &#x27;2697049&#x27; (SARS-CoV-2), &#x27;11320&#x27; (Influenza A), &#x27;1773&#x27; (M. tuberculosis). |
| `protein_name` | string &#124; null | None | Target protein name to filter epitopes. Examples: &#x27;Spike glycoprotein&#x27;, &#x27;Nucleocapsid protein&#x27;, &#x27;Replicase polyprotein 1ab&#x27;. |
| `epitope_type` | string &#124; null | None | Type of epitope. Options: &#x27;Linear peptide&#x27;, &#x27;Discontinuous peptide&#x27;. Most epitopes are linear peptides. |
| `organism` | string &#124; null | None | Organism name keyword to search across epitope records. Example: &#x27;coronavirus&#x27;, &#x27;influenza&#x27;. |
| `limit` | integer &#124; null | None | Maximum number of results. Default: 25. Max: 100. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "taxon_id": {
          "type": [
            "string",
            "null"
          ],
          "description": "NCBI Taxonomy ID for the pathogen organism. Examples: '2697049' (SARS-CoV-2), '11320' (Influenza A), '1773' (M. tuberculosis)."
        },
        "protein_name": {
          "type": [
            "string",
            "null"
          ],
          "description": "Target protein name to filter epitopes. Examples: 'Spike glycoprotein', 'Nucleocapsid protein', 'Replicase polyprotein 1ab'."
        },
        "epitope_type": {
          "type": [
            "string",
            "null"
          ],
          "description": "Type of epitope. Options: 'Linear peptide', 'Discontinuous peptide'. Most epitopes are linear peptides."
        },
        "organism": {
          "type": [
            "string",
            "null"
          ],
          "description": "Organism name keyword to search across epitope records. Example: 'coronavirus', 'influenza'."
        },
        "limit": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Maximum number of results. Default: 25. Max: 100."
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
          "items": {
            "type": "object",
            "properties": {
              "epitope_id": {
                "type": "string"
              },
              "epitope_type": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "epitope_sequence": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "organism": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "protein_name": {
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
              "bcell_assays": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "tcell_assays": {
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

**执行与来源：** `network`；上游接口名 `BVBRC_search_epitopes`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("bvbrc_search_epitopes", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.enrichr_enrich` {#dnakit.databases.enrichr_enrich}

Perform gene set enrichment analysis using Enrichr. Submit a list of genes and get enriched terms from a specific library. Returns ranked terms with p-values, z-scores, combined scores, and overlapping genes. No authentication required. Set species to query the fly, worm, yeast, or fish Enrichr instances; libraries and gene symbols differ per organism.

```text
dnakit.databases.enrichr_enrich(*, operation=None, gene_list, library=None, top_n=None, species=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: enrich) {&quot;const&quot;: &quot;enrich&quot;, &quot;default&quot;: &quot;enrich&quot;} |
| `gene_list` | array | 必填（输入 Schema） | List of gene symbols (e.g., &#91;&#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;EGFR&#x27;&#93;). Minimum 1 gene. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `library` | string | None | Enrichr library name (default: GO_Biological_Process_2023). Common: GO_Biological_Process_2023, KEGG_2021_Human, Reactome_2022, WikiPathways_2024_Human, MSigDB_Hallmark_2020 {&quot;default&quot;: &quot;GO_Biological_Process_2023&quot;} |
| `top_n` | integer | None | Number of top enriched terms to return (default: 10) {&quot;default&quot;: 10} |
| `species` | string &#124; null | None | Organism instance to query: &#x27;human&#x27; (default), &#x27;mouse&#x27;, &#x27;fly&#x27; (FlyEnrichr), &#x27;worm&#x27; (WormEnrichr), &#x27;yeast&#x27; (YeastEnrichr), or &#x27;fish&#x27; (FishEnrichr). Gene set libraries differ per organism, and gene symbols must match that organism&#x27;s nomenclature. {&quot;enum&quot;: &#91;null, &quot;human&quot;, &quot;mouse&quot;, &quot;fly&quot;, &quot;worm&quot;, &quot;yeast&quot;, &quot;fish&quot;&#93;} |
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
          "const": "enrich",
          "description": "Operation type (fixed: enrich)",
          "default": "enrich"
        },
        "gene_list": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "List of gene symbols (e.g., ['TP53', 'BRCA1', 'EGFR']). Minimum 1 gene."
        },
        "library": {
          "type": "string",
          "description": "Enrichr library name (default: GO_Biological_Process_2023). Common: GO_Biological_Process_2023, KEGG_2021_Human, Reactome_2022, WikiPathways_2024_Human, MSigDB_Hallmark_2020",
          "default": "GO_Biological_Process_2023"
        },
        "top_n": {
          "type": "integer",
          "description": "Number of top enriched terms to return (default: 10)",
          "default": 10
        },
        "species": {
          "type": [
            "string",
            "null"
          ],
          "enum": [
            null,
            "human",
            "mouse",
            "fly",
            "worm",
            "yeast",
            "fish"
          ],
          "description": "Organism instance to query: 'human' (default), 'mouse', 'fly' (FlyEnrichr), 'worm' (WormEnrichr), 'yeast' (YeastEnrichr), or 'fish' (FishEnrichr). Gene set libraries differ per organism, and gene symbols must match that organism's nomenclature."
        }
      },
      "required": [
        "gene_list"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object"
    }
    ```

**执行与来源：** `network`；上游接口名 `Enrichr_enrich`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("enrichr_enrich", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.enrichr_gene_enrichment_analysis` {#dnakit.databases.enrichr_gene_enrichment_analysis}

Perform gene enrichment analysis using Enrichr to find biological pathways, processes, and molecular functions associated with a list of genes. Returns connectivity paths between genes and enrichment terms.

```text
dnakit.databases.enrichr_gene_enrichment_analysis(*, gene_list, libs=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_list` | array | 必填（输入 Schema） | List of gene names or symbols to analyze. At least 2 genes are required for path ranking analysis. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}, &quot;minItems&quot;: 2} |
| `libs` | array | None | List of Enrichr enrichment libraries to query (e.g. &#x27;KEGG_2021_Human&#x27;, &#x27;GO_Biological_Process_2023&#x27;). Omitting this or passing an empty list uses a broad default set of pathway/ontology libraries. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}, &quot;default&quot;: &#91;&quot;WikiPathways_2024_Human&quot;, &quot;Reactome_Pathways_2024&quot;, &quot;MSigDB_Hallmark_2020&quot;, &quot;GO_Molecular_Function_2023&quot;, &quot;GO_Biological_Process_2023&quot;&#93;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene_list": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "List of gene names or symbols to analyze. At least 2 genes are required for path ranking analysis.",
          "minItems": 2
        },
        "libs": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "List of Enrichr enrichment libraries to query (e.g. 'KEGG_2021_Human', 'GO_Biological_Process_2023'). Omitting this or passing an empty list uses a broad default set of pathway/ontology libraries.",
          "default": [
            "WikiPathways_2024_Human",
            "Reactome_Pathways_2024",
            "MSigDB_Hallmark_2020",
            "GO_Molecular_Function_2023",
            "GO_Biological_Process_2023"
          ]
        }
      },
      "required": [
        "gene_list"
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
        "connected_paths": {
          "type": "object",
          "description": "Ranked weighted paths between the first two genes routed through shared enrichment terms. Each key is 'Path: [gene, term, ..., gene]' and its value is 'Total Weight: <sum>', where the weight of each gene-term edge is that term's Enrichr combined score (higher = stronger enrichment support for that gene through that pathway/term)."
        },
        "connections": {
          "type": "object",
          "description": "Per gene-to-enrichment-term connectivity paths found in the enrichment graph (only real pathway/ontology terms are considered as targets, not other genes). Each key is 'Connectivity: <gene> - <term>' and its value is a list of (path, weight) pairs, where weight is the sum of Enrichr combined scores along that path; higher weight means stronger combined enrichment evidence linking the gene to that term."
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `enrichr_gene_enrichment_analysis`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("enrichr_gene_enrichment_analysis", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.go_get_annotations_for_gene` {#dnakit.databases.go_get_annotations_for_gene}

Finds GO annotations for a specific gene/protein using GOlr search. Returns one row per annotation, so a gene with several evidence codes for the same GO term appears several times. This returns a PAGE, not the complete set: heavily annotated genes have far more annotations than the default page holds (measured: 'TP53' matches 1812 annotation rows), and the response carries no total, so a short list is NOT evidence the gene has no further annotations -- raise 'rows' if completeness matters.

```text
dnakit.databases.go_get_annotations_for_gene(*, gene_id, rows=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填（输入 Schema） | A gene identifier such as gene symbol (e.g., &#x27;TP53&#x27;) or database ID. |
| `rows` | integer | None | Maximum number of annotation rows to return. Default: 100. Raise it (e.g. 2000) for a heavily annotated gene such as TP53, which has 1812 annotation rows; lower it to keep responses small. {&quot;minimum&quot;: 1, &quot;default&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene_id": {
          "type": "string",
          "description": "A gene identifier such as gene symbol (e.g., 'TP53') or database ID."
        },
        "rows": {
          "type": "integer",
          "description": "Maximum number of annotation rows to return. Default: 100. Raise it (e.g. 2000) for a heavily annotated gene such as TP53, which has 1812 annotation rows; lower it to keep responses small.",
          "minimum": 1,
          "default": 100
        }
      },
      "required": [
        "gene_id"
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
          "document_category": {
            "type": "string"
          },
          "source": {
            "type": "string"
          },
          "bioentity": {
            "type": "string"
          },
          "bioentity_internal_id": {
            "type": "string"
          },
          "bioentity_label": {
            "type": "string"
          },
          "bioentity_label_searchable": {
            "type": "string"
          },
          "annotation_class": {
            "type": "string"
          },
          "annotation_class_label": {
            "type": "string"
          },
          "annotation_class_label_searchable": {
            "type": "string"
          },
          "aspect": {
            "type": "string"
          },
          "bioentity_name": {
            "type": "string"
          },
          "bioentity_name_searchable": {
            "type": "string"
          },
          "type": {
            "type": "string"
          },
          "date": {
            "type": "string"
          },
          "assigned_by": {
            "type": "string"
          },
          "bioentity_isoform": {
            "type": "string"
          },
          "taxon": {
            "type": "string"
          },
          "taxon_label": {
            "type": "string"
          },
          "taxon_label_searchable": {
            "type": "string"
          },
          "taxon_closure_map": {
            "type": "string"
          },
          "taxon_subset_closure_map": {
            "type": "string"
          },
          "panther_family": {
            "type": "string"
          },
          "panther_family_searchable": {
            "type": "string"
          },
          "panther_family_label": {
            "type": "string"
          },
          "panther_family_label_searchable": {
            "type": "string"
          },
          "evidence_type": {
            "type": "string"
          },
          "evidence": {
            "type": "string"
          },
          "evidence_label": {
            "type": "string"
          },
          "evidence_label_searchable": {
            "type": "string"
          },
          "evidence_closure_map": {
            "type": "string"
          },
          "evidence_subset_closure_map": {
            "type": "string"
          },
          "isa_partof_closure_map": {
            "type": "string"
          },
          "regulates_closure_map": {
            "type": "string"
          },
          "id": {
            "type": "string"
          },
          "isa_partof_closure_label": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "taxon_subset_closure_label_searchable": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "taxon_subset_closure": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "evidence_with": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "synonym_searchable": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "isa_partof_closure_label_searchable": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "evidence_with_searchable": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "reference": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "synonym": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "evidence_closure_label_searchable": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "regulates_closure_label_searchable": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "taxon_closure_label_searchable": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "isa_partof_closure": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "taxon_closure": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "evidence_closure": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "evidence_subset_closure_label_searchable": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "evidence_closure_label": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "evidence_subset_closure": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "evidence_subset_closure_label": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "taxon_closure_label": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "regulates_closure": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "evidence_type_closure": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "regulates_closure_label": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "taxon_subset_closure_label": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "reference_searchable": {
            "type": "array",
            "items": {
              "type": "string"
            }
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `GO_get_annotations_for_gene`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("go_get_annotations_for_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.go_search_terms` {#dnakit.databases.go_search_terms}

Searches for Gene Ontology (GO) terms by a keyword using the GOlr search engine. Returns GO terms and related biological entities.

```text
dnakit.databases.go_search_terms(*, query, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填（输入 Schema） | The keyword to search for, e.g., &#x27;apoptosis&#x27; or &#x27;kinase activity&#x27;. |
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
          "description": "The keyword to search for, e.g., 'apoptosis' or 'kinase activity'."
        }
      },
      "required": [
        "query"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "responseHeader": {
          "type": "object",
          "properties": {
            "status": {
              "type": "integer"
            },
            "QTime": {
              "type": "integer"
            },
            "params": {
              "type": "object",
              "properties": {
                "q": {
                  "type": "string"
                },
                "fq": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  }
                },
                "rows": {
                  "type": "string"
                },
                "wt": {
                  "type": "string"
                }
              }
            }
          }
        },
        "response": {
          "type": "object",
          "properties": {
            "numFound": {
              "type": "integer"
            },
            "start": {
              "type": "integer"
            },
            "docs": {
              "type": "array"
            }
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `GO_search_terms`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("go_search_terms", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.humanbase_ppi_analysis` {#dnakit.databases.humanbase_ppi_analysis}

Retrieve tissue-specific protein-protein interactions and biological processes from HumanBase. Returns a NetworkX graph of tissue specific protein-protein interactions and a list of associated biological processes involeed by the given genes from Gene Ontology.

```text
dnakit.databases.humanbase_ppi_analysis(*, gene_list=None, tissue=None, max_node=None, interaction=None, string_mode=None, genes=None, top_n=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_list` | array | None | List of gene names or symbols to analyze for protein-protein interactions. The gene name should be the official gene symbol, not the synonym. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}, &quot;minItems&quot;: 1} |
| `tissue` | string | None | Tissue type for tissue-specific interactions. Use simple names like &#x27;brain&#x27;, &#x27;heart&#x27;, &#x27;liver&#x27;, &#x27;kidney&#x27;, &#x27;lung&#x27;, &#x27;blood&#x27;. For v3 networks with more cell-type specificity, append &#x27;-v3&#x27; (e.g., &#x27;brain-v3&#x27;, &#x27;liver-v3&#x27;). {&quot;default&quot;: &quot;brain&quot;} |
| `max_node` | integer | None | Maximum number of nodes to retrieve in the interaction network. Warning: the more nodes, the more time it takes to retrieve the data. Default is 10 (~30 seconds). {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `interaction` | string | None | Specific interaction type to filter by. Available types: &#x27;co-expression&#x27;, &#x27;interaction&#x27;, &#x27;tf-binding&#x27;, &#x27;gsea-microrna-targets&#x27;, &#x27;gsea-perturbations&#x27;. If not specified, all types will be included. {&quot;enum&quot;: &#91;&quot;co-expression&quot;, &quot;interaction&quot;, &quot;tf-binding&quot;, &quot;gsea-microrna-targets&quot;, &quot;gsea-perturbations&quot;&#93;} |
| `string_mode` | boolean | None | Whether to return the result in string mode. If True, the result will be a string of the network graph and the biological processes. If False, the result will be a NetworkX graph and a list of biological processes. {&quot;default&quot;: true} |
| `genes` | array | None | Alias for gene_list. List of gene symbols to analyze. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `top_n` | integer | None | Alias for max_node. Maximum number of network nodes to retrieve. {&quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene_list": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "List of gene names or symbols to analyze for protein-protein interactions. The gene name should be the official gene symbol, not the synonym.",
          "minItems": 1
        },
        "tissue": {
          "type": "string",
          "description": "Tissue type for tissue-specific interactions. Use simple names like 'brain', 'heart', 'liver', 'kidney', 'lung', 'blood'. For v3 networks with more cell-type specificity, append '-v3' (e.g., 'brain-v3', 'liver-v3').",
          "default": "brain"
        },
        "max_node": {
          "type": "integer",
          "description": "Maximum number of nodes to retrieve in the interaction network. Warning: the more nodes, the more time it takes to retrieve the data. Default is 10 (~30 seconds).",
          "default": 10,
          "minimum": 1,
          "maximum": 100
        },
        "interaction": {
          "type": "string",
          "description": "Specific interaction type to filter by. Available types: 'co-expression', 'interaction', 'tf-binding', 'gsea-microrna-targets', 'gsea-perturbations'. If not specified, all types will be included.",
          "enum": [
            "co-expression",
            "interaction",
            "tf-binding",
            "gsea-microrna-targets",
            "gsea-perturbations"
          ]
        },
        "string_mode": {
          "type": "boolean",
          "description": "Whether to return the result in string mode. If True, the result will be a string of the network graph and the biological processes. If False, the result will be a NetworkX graph and a list of biological processes.",
          "default": true
        },
        "genes": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "Alias for gene_list. List of gene symbols to analyze."
        },
        "top_n": {
          "type": "integer",
          "description": "Alias for max_node. Maximum number of network nodes to retrieve.",
          "minimum": 1,
          "maximum": 100
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "string",
      "description": "Formatted string representation of the protein-protein interaction network, including network summary, protein nodes with descriptions, interaction details with evidence types, associated biological processes from Gene Ontology, and network analysis metrics."
    }
    ```

**执行与来源：** `network`；上游接口名 `humanbase_ppi_analysis`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("humanbase_ppi_analysis", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.iedb_get_epitope_antigens` {#dnakit.databases.iedb_get_epitope_antigens}

Given an epitope `structure_id`, get linked antigens (IEDB Query API). Returned rows contain `parent_source_antigen_iri` (an antigen identifier) that you can use in other filters/searches.

```text
dnakit.databases.iedb_get_epitope_antigens(*, structure_id, limit=None, offset=None, order=None, select=None, filters=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `structure_id` | integer | 必填（输入 Schema） | IEDB epitope structure ID (from `iedb_search_epitopes`). |
| `limit` | integer | None | Maximum number of rows to return. {&quot;default&quot;: 10} |
| `offset` | integer | None | Pagination offset. {&quot;default&quot;: 0} |
| `order` | string | None | Optional PostgREST order clause. |
| `select` | object &#124; None | None | Optional projection list (comma-separated string or array of field names). {&quot;anyOf&quot;: &#91;{&quot;type&quot;: &quot;string&quot;}, {&quot;type&quot;: &quot;array&quot;, &quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}}&#93;} |
| `filters` | object | None | Advanced PostgREST filters mapping column -&gt; filter expression. {&quot;additionalProperties&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "structure_id": {
          "type": "integer",
          "description": "IEDB epitope structure ID (from `iedb_search_epitopes`)."
        },
        "limit": {
          "type": "integer",
          "default": 10,
          "description": "Maximum number of rows to return."
        },
        "offset": {
          "type": "integer",
          "default": 0,
          "description": "Pagination offset."
        },
        "order": {
          "type": "string",
          "description": "Optional PostgREST order clause."
        },
        "select": {
          "description": "Optional projection list (comma-separated string or array of field names).",
          "anyOf": [
            {
              "type": "string"
            },
            {
              "type": "array",
              "items": {
                "type": "string"
              }
            }
          ]
        },
        "filters": {
          "type": "object",
          "description": "Advanced PostgREST filters mapping column -> filter expression.",
          "additionalProperties": {
            "type": "string"
          }
        }
      },
      "required": [
        "structure_id"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "array",
      "description": "List of epitope-to-antigen linkage rows",
      "items": {
        "type": "object",
        "properties": {
          "structure_id": {
            "type": [
              "integer",
              "null"
            ]
          },
          "parent_source_antigen_iri": {
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

**执行与来源：** `network`；上游接口名 `iedb_get_epitope_antigens`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("iedb_get_epitope_antigens", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.iedb_get_epitope_mhc` {#dnakit.databases.iedb_get_epitope_mhc}

Given an epitope `structure_id`, get linked MHC restriction association rows (IEDB Query API). Returned rows include `elution_id` for the association.

```text
dnakit.databases.iedb_get_epitope_mhc(*, structure_id, limit=None, offset=None, order=None, select=None, filters=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `structure_id` | integer | 必填（输入 Schema） | IEDB epitope structure ID (from `iedb_search_epitopes`). |
| `limit` | integer | None | Maximum number of rows to return. {&quot;default&quot;: 10} |
| `offset` | integer | None | Pagination offset. {&quot;default&quot;: 0} |
| `order` | string | None | Optional PostgREST order clause. |
| `select` | object &#124; None | None | Optional projection list (comma-separated string or array of field names). {&quot;anyOf&quot;: &#91;{&quot;type&quot;: &quot;string&quot;}, {&quot;type&quot;: &quot;array&quot;, &quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}}&#93;} |
| `filters` | object | None | Advanced PostgREST filters mapping column -&gt; filter expression. {&quot;additionalProperties&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "structure_id": {
          "type": "integer",
          "description": "IEDB epitope structure ID (from `iedb_search_epitopes`)."
        },
        "limit": {
          "type": "integer",
          "default": 10,
          "description": "Maximum number of rows to return."
        },
        "offset": {
          "type": "integer",
          "default": 0,
          "description": "Pagination offset."
        },
        "order": {
          "type": "string",
          "description": "Optional PostgREST order clause."
        },
        "select": {
          "description": "Optional projection list (comma-separated string or array of field names).",
          "anyOf": [
            {
              "type": "string"
            },
            {
              "type": "array",
              "items": {
                "type": "string"
              }
            }
          ]
        },
        "filters": {
          "type": "object",
          "description": "Advanced PostgREST filters mapping column -> filter expression.",
          "additionalProperties": {
            "type": "string"
          }
        }
      },
      "required": [
        "structure_id"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "array",
      "description": "List of epitope-to-MHC linkage rows",
      "items": {
        "type": "object",
        "properties": {
          "structure_id": {
            "type": [
              "integer",
              "null"
            ]
          },
          "elution_id": {
            "type": [
              "integer",
              "null"
            ]
          }
        },
        "additionalProperties": true
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `iedb_get_epitope_mhc`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("iedb_get_epitope_mhc", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.iedb_search_bcell` {#dnakit.databases.iedb_search_bcell}

Search B-cell assay rows (IEDB Query API). Use this to discover `bcell_id` values and related epitope `structure_id` values.

Tip: use `select` to request only a few fields (smaller + faster responses).

```text
dnakit.databases.iedb_search_bcell(*, limit=None, offset=None, order=None, select=None, filters=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `limit` | integer | None | Maximum number of rows to return. {&quot;default&quot;: 10} |
| `offset` | integer | None | Pagination offset. {&quot;default&quot;: 0} |
| `order` | string | None | Optional PostgREST order clause. |
| `select` | object &#124; None | None | Optional projection list (comma-separated string or array of field names). {&quot;anyOf&quot;: &#91;{&quot;type&quot;: &quot;string&quot;}, {&quot;type&quot;: &quot;array&quot;, &quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}}&#93;} |
| `filters` | object | None | Advanced PostgREST filters mapping column -&gt; filter expression. {&quot;additionalProperties&quot;: {&quot;type&quot;: &quot;string&quot;}} |
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
          "description": "Maximum number of rows to return."
        },
        "offset": {
          "type": "integer",
          "default": 0,
          "description": "Pagination offset."
        },
        "order": {
          "type": "string",
          "description": "Optional PostgREST order clause."
        },
        "select": {
          "description": "Optional projection list (comma-separated string or array of field names).",
          "anyOf": [
            {
              "type": "string"
            },
            {
              "type": "array",
              "items": {
                "type": "string"
              }
            }
          ]
        },
        "filters": {
          "type": "object",
          "description": "Advanced PostgREST filters mapping column -> filter expression.",
          "additionalProperties": {
            "type": "string"
          }
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "array",
      "description": "List of B-cell assay rows",
      "items": {
        "type": "object",
        "properties": {
          "bcell_id": {
            "type": [
              "integer",
              "null"
            ]
          },
          "structure_id": {
            "type": [
              "integer",
              "null"
            ]
          },
          "assay_description": {
            "type": [
              "string",
              "null"
            ]
          },
          "reference_id": {
            "type": [
              "integer",
              "null"
            ]
          }
        },
        "additionalProperties": true
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `iedb_search_bcell`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("iedb_search_bcell", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.iedb_search_epitopes` {#dnakit.databases.iedb_search_epitopes}

Search immune epitopes (IEDB Query API). Use this to discover epitope `structure_id` values, then use the epitope linking tools (e.g., `iedb_get_epitope_references`) to retrieve related references, MHC, antigens, etc.

Tips:
- Use `select` to request only a few fields (smaller + faster responses).
- For advanced filtering, pass `filters` with values like `eq.123` or `ilike.*KVF*`.

```text
dnakit.databases.iedb_search_epitopes(*, sequence_contains=None, structure_type=None, limit=None, offset=None, order=None, select=None, filters=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `sequence_contains` | string | None | Sequence fragment to search within `linear_sequence` using a contains match. Example: &#x27;KVF&#x27; finds epitopes containing that motif. |
| `structure_type` | string | None | Filter by epitope structure type. Example: &#x27;Linear peptide&#x27;. {&quot;enum&quot;: &#91;&quot;Linear peptide&quot;, &quot;Discontinuous peptide&quot;, &quot;Non-peptidic&quot;&#93;} |
| `limit` | integer | None | Maximum number of rows to return. {&quot;default&quot;: 10} |
| `offset` | integer | None | Pagination offset. {&quot;default&quot;: 0} |
| `order` | string | None | Optional PostgREST order clause (e.g., &#x27;structure_id.asc&#x27;). |
| `select` | object &#124; None | None | Optional projection list to reduce payload size. Provide a comma-separated string or an array of field names (e.g., &#91;&#x27;structure_id&#x27;,&#x27;linear_sequence&#x27;&#93;). {&quot;anyOf&quot;: &#91;{&quot;type&quot;: &quot;string&quot;}, {&quot;type&quot;: &quot;array&quot;, &quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}}&#93;} |
| `filters` | object | None | Advanced PostgREST filters mapping column -&gt; filter expression (e.g., {&quot;linear_sequence&quot;:&quot;ilike.*KVF*&quot;}). {&quot;additionalProperties&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "sequence_contains": {
          "type": "string",
          "description": "Sequence fragment to search within `linear_sequence` using a contains match. Example: 'KVF' finds epitopes containing that motif."
        },
        "structure_type": {
          "type": "string",
          "description": "Filter by epitope structure type. Example: 'Linear peptide'.",
          "enum": [
            "Linear peptide",
            "Discontinuous peptide",
            "Non-peptidic"
          ]
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of rows to return.",
          "default": 10
        },
        "offset": {
          "type": "integer",
          "description": "Pagination offset.",
          "default": 0
        },
        "order": {
          "type": "string",
          "description": "Optional PostgREST order clause (e.g., 'structure_id.asc')."
        },
        "select": {
          "description": "Optional projection list to reduce payload size. Provide a comma-separated string or an array of field names (e.g., ['structure_id','linear_sequence']).",
          "anyOf": [
            {
              "type": "string"
            },
            {
              "type": "array",
              "items": {
                "type": "string"
              }
            }
          ]
        },
        "filters": {
          "type": "object",
          "description": "Advanced PostgREST filters mapping column -> filter expression (e.g., {\"linear_sequence\":\"ilike.*KVF*\"}).",
          "additionalProperties": {
            "type": "string"
          }
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "array",
      "description": "List of matching epitope rows (projected via `select` if provided).",
      "items": {
        "type": "object",
        "properties": {
          "structure_id": {
            "type": [
              "integer",
              "string",
              "null"
            ]
          },
          "structure_type": {
            "type": [
              "string",
              "null"
            ]
          },
          "linear_sequence": {
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

**执行与来源：** `network`；上游接口名 `iedb_search_epitopes`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("iedb_search_epitopes", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.iedb_search_mhc` {#dnakit.databases.iedb_search_mhc}

Search MHC-related assay rows (IEDB Query API). Use this to discover epitope `structure_id` values involved in MHC restriction.

Tip: use `select` to request only a few fields (smaller + faster responses).

```text
dnakit.databases.iedb_search_mhc(*, limit=None, offset=None, order=None, select=None, filters=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `limit` | integer | None | Maximum number of rows to return. {&quot;default&quot;: 10} |
| `offset` | integer | None | Pagination offset. {&quot;default&quot;: 0} |
| `order` | string | None | Optional PostgREST order clause. |
| `select` | object &#124; None | None | Optional projection list (comma-separated string or array of field names). {&quot;anyOf&quot;: &#91;{&quot;type&quot;: &quot;string&quot;}, {&quot;type&quot;: &quot;array&quot;, &quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}}&#93;} |
| `filters` | object | None | Advanced PostgREST filters mapping column -&gt; filter expression. {&quot;additionalProperties&quot;: {&quot;type&quot;: &quot;string&quot;}} |
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
          "description": "Maximum number of rows to return."
        },
        "offset": {
          "type": "integer",
          "default": 0,
          "description": "Pagination offset."
        },
        "order": {
          "type": "string",
          "description": "Optional PostgREST order clause."
        },
        "select": {
          "description": "Optional projection list (comma-separated string or array of field names).",
          "anyOf": [
            {
              "type": "string"
            },
            {
              "type": "array",
              "items": {
                "type": "string"
              }
            }
          ]
        },
        "filters": {
          "type": "object",
          "description": "Advanced PostgREST filters mapping column -> filter expression.",
          "additionalProperties": {
            "type": "string"
          }
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "array",
      "description": "List of MHC-related assay rows",
      "items": {
        "type": "object",
        "properties": {
          "structure_id": {
            "type": [
              "integer",
              "null"
            ]
          },
          "mhc_restriction": {
            "type": [
              "string",
              "null"
            ]
          },
          "mhc_class": {
            "type": [
              "string",
              "null"
            ]
          },
          "mhc_allele_resolution": {
            "type": [
              "string",
              "null"
            ]
          },
          "reference_id": {
            "type": [
              "integer",
              "null"
            ]
          }
        },
        "additionalProperties": true
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `iedb_search_mhc`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("iedb_search_mhc", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.iedb_search_tcell_assays` {#dnakit.databases.iedb_search_tcell_assays}

Search T-cell assay data from the IEDB. Returns T-cell response assays with epitope sequences, MHC restriction, source organism, and assay results.

Useful for finding T-cell epitopes for vaccine design, immunotherapy, and understanding cellular immunity. Filter by peptide sequence, MHC class, organism, or qualitative result.

Tip: use `select` to request only key fields for faster responses.

```text
dnakit.databases.iedb_search_tcell_assays(*, sequence=None, sequence_contains=None, mhc_class=None, qualitative_measure=None, limit=None, offset=None, select=None, filters=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `sequence` | string | None | Exact peptide sequence to search (e.g., GILGFVFTL for influenza M1 epitope). |
| `sequence_contains` | string | None | Partial peptide sequence (substring match). Example: &#x27;SIINFEKL&#x27;. |
| `mhc_class` | string | None | MHC class restriction (I or II). {&quot;enum&quot;: &#91;&quot;I&quot;, &quot;II&quot;&#93;} |
| `qualitative_measure` | string | None | Filter by assay result. {&quot;enum&quot;: &#91;&quot;Positive&quot;, &quot;Negative&quot;, &quot;Positive-High&quot;, &quot;Positive-Intermediate&quot;, &quot;Positive-Low&quot;&#93;} |
| `limit` | integer | None | Maximum rows to return. {&quot;default&quot;: 10} |
| `offset` | integer | None | Pagination offset. {&quot;default&quot;: 0} |
| `select` | object &#124; None | None | Columns to return (array of strings). Example: &#91;&#x27;tcell_id&#x27;,&#x27;linear_sequence&#x27;,&#x27;mhc_restriction&#x27;,&#x27;qualitative_measure&#x27;&#93;. {&quot;oneOf&quot;: &#91;{&quot;type&quot;: &quot;string&quot;}, {&quot;type&quot;: &quot;array&quot;, &quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}}&#93;} |
| `filters` | object | None | Advanced PostgREST filters (e.g., {&quot;source_organism_iri&quot;:&quot;eq.NCBITaxon:11320&quot;}). {&quot;additionalProperties&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "sequence": {
          "type": "string",
          "description": "Exact peptide sequence to search (e.g., GILGFVFTL for influenza M1 epitope)."
        },
        "sequence_contains": {
          "type": "string",
          "description": "Partial peptide sequence (substring match). Example: 'SIINFEKL'."
        },
        "mhc_class": {
          "type": "string",
          "description": "MHC class restriction (I or II).",
          "enum": [
            "I",
            "II"
          ]
        },
        "qualitative_measure": {
          "type": "string",
          "description": "Filter by assay result.",
          "enum": [
            "Positive",
            "Negative",
            "Positive-High",
            "Positive-Intermediate",
            "Positive-Low"
          ]
        },
        "limit": {
          "type": "integer",
          "description": "Maximum rows to return.",
          "default": 10
        },
        "offset": {
          "type": "integer",
          "description": "Pagination offset.",
          "default": 0
        },
        "select": {
          "description": "Columns to return (array of strings). Example: ['tcell_id','linear_sequence','mhc_restriction','qualitative_measure'].",
          "oneOf": [
            {
              "type": "string"
            },
            {
              "type": "array",
              "items": {
                "type": "string"
              }
            }
          ]
        },
        "filters": {
          "type": "object",
          "description": "Advanced PostgREST filters (e.g., {\"source_organism_iri\":\"eq.NCBITaxon:11320\"}).",
          "additionalProperties": {
            "type": "string"
          }
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "array",
      "description": "List of T-cell assay records.",
      "items": {
        "type": "object",
        "properties": {
          "tcell_id": {
            "type": [
              "integer",
              "null"
            ]
          },
          "linear_sequence": {
            "type": [
              "string",
              "null"
            ]
          },
          "structure_type": {
            "type": [
              "string",
              "null"
            ]
          },
          "mhc_restriction": {
            "type": [
              "string",
              "null"
            ]
          },
          "mhc_class": {
            "type": [
              "string",
              "null"
            ]
          },
          "qualitative_measure": {
            "type": [
              "string",
              "null"
            ]
          },
          "pubmed_id": {
            "type": [
              "string",
              "null"
            ]
          },
          "journal_name": {
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

**执行与来源：** `network`；上游接口名 `iedb_search_tcell_assays`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("iedb_search_tcell_assays", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.imgt_get_gene_info` {#dnakit.databases.imgt_get_gene_info}

Get information about IMGT databases, gene nomenclature, and analysis tools. Returns descriptions of IMGT resources (LIGM-DB, GENE-DB, V-QUEST) and gene naming conventions. Good starting point for immunogenetics research.

```text
dnakit.databases.imgt_get_gene_info(*, operation=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | 提供方操作模式；默认值与允许值见本行约束。 {&quot;const&quot;: &quot;get_gene_info&quot;, &quot;default&quot;: &quot;get_gene_info&quot;} |
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
          "const": "get_gene_info",
          "default": "get_gene_info"
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
          "type": "string"
        },
        "data": {
          "type": "object"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `IMGT_get_gene_info`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("imgt_get_gene_info", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.imgt_search_genes` {#dnakit.databases.imgt_search_genes}

Search IMGT for immunoglobulin/TCR germline genes. Filter by gene type (IGHV, IGKV, TRAV, etc.) and species. Returns search URLs and gene type information for navigating IMGT GENE-DB.

```text
dnakit.databases.imgt_search_genes(*, operation=None, query=None, gene_type=None, species=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | 提供方操作模式；默认值与允许值见本行约束。 {&quot;const&quot;: &quot;search_genes&quot;, &quot;default&quot;: &quot;search_genes&quot;} |
| `query` | string | None | Search query (gene name) |
| `gene_type` | string | None | Gene type: IGHV, IGHD, IGHJ, IGKV, IGLV, TRAV, TRBV, etc. |
| `species` | string | None | Species name (default: Homo sapiens) {&quot;default&quot;: &quot;Homo sapiens&quot;} |
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
          "const": "search_genes",
          "default": "search_genes"
        },
        "query": {
          "type": "string",
          "description": "Search query (gene name)"
        },
        "gene_type": {
          "type": "string",
          "description": "Gene type: IGHV, IGHD, IGHJ, IGKV, IGLV, TRAV, TRBV, etc."
        },
        "species": {
          "type": "string",
          "description": "Species name (default: Homo sapiens)",
          "default": "Homo sapiens"
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
          "type": "string"
        },
        "data": {
          "type": "object"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `IMGT_search_genes`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("imgt_search_genes", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.imm_port_search_studies` {#dnakit.databases.imm_port_search_studies}

Search the ImmPort immunology database for studies by keyword, disease condition, assay method, or research focus. ImmPort is a NIAID-funded repository with 900+ open studies covering vaccine trials, flow cytometry, ELISPOT, RNA-seq, and clinical immunology data. Assay methods include: Flow Cytometry, ELISPOT, ELISA, RNA sequencing, Transcription profiling by array, Hemagglutination Inhibition, Luminex xMAP, Q-PCR, CyTOF, and more. Research focuses include: Vaccine Response, Infection Response, Immune Response, Transplantation, Autoimmune, Allergy, Preterm Birth. The species filter is applied server-side; condition_or_disease, assay_method, and research_focus filters are applied client-side.

```text
dnakit.databases.imm_port_search_studies(*, query, condition_or_disease=None, assay_method=None, research_focus=None, species=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填（输入 Schema） | Search keywords (e.g., &#x27;influenza vaccine&#x27;, &#x27;COVID-19&#x27;, &#x27;T cell response&#x27;, &#x27;malaria&#x27;) |
| `condition_or_disease` | string | None | Filter by disease/condition (e.g., &#x27;influenza&#x27;, &#x27;malaria&#x27;, &#x27;HIV&#x27;, &#x27;dengue&#x27;, &#x27;tuberculosis&#x27;). Client-side substring match. |
| `assay_method` | string | None | Filter by assay method (e.g., &#x27;Flow Cytometry&#x27;, &#x27;ELISPOT&#x27;, &#x27;RNA sequencing&#x27;, &#x27;ELISA&#x27;, &#x27;CyTOF&#x27;, &#x27;Luminex xMAP&#x27;). Client-side substring match. |
| `research_focus` | string | None | Filter by research focus (e.g., &#x27;Vaccine Response&#x27;, &#x27;Infection Response&#x27;, &#x27;Immune Response&#x27;, &#x27;Transplantation&#x27;, &#x27;Autoimmune&#x27;). Client-side substring match. |
| `species` | string | None | Filter by species (e.g., &#x27;Homo sapiens&#x27;, &#x27;Mus musculus&#x27;). Applied server-side. |
| `limit` | integer | None | Maximum number of studies to return (1-100, default 10) {&quot;default&quot;: 10} |
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
          "description": "Search keywords (e.g., 'influenza vaccine', 'COVID-19', 'T cell response', 'malaria')"
        },
        "condition_or_disease": {
          "type": "string",
          "description": "Filter by disease/condition (e.g., 'influenza', 'malaria', 'HIV', 'dengue', 'tuberculosis'). Client-side substring match."
        },
        "assay_method": {
          "type": "string",
          "description": "Filter by assay method (e.g., 'Flow Cytometry', 'ELISPOT', 'RNA sequencing', 'ELISA', 'CyTOF', 'Luminex xMAP'). Client-side substring match."
        },
        "research_focus": {
          "type": "string",
          "description": "Filter by research focus (e.g., 'Vaccine Response', 'Infection Response', 'Immune Response', 'Transplantation', 'Autoimmune'). Client-side substring match."
        },
        "species": {
          "type": "string",
          "description": "Filter by species (e.g., 'Homo sapiens', 'Mus musculus'). Applied server-side."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of studies to return (1-100, default 10)",
          "default": 10
        }
      },
      "required": [
        "query"
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
          "study_accession": {
            "type": [
              "string",
              "null"
            ]
          },
          "title": {
            "type": [
              "string",
              "null"
            ]
          },
          "brief_description": {
            "type": [
              "string",
              "null"
            ]
          },
          "condition_or_disease": {
            "type": [
              "array",
              "null"
            ]
          },
          "research_focus": {
            "type": [
              "array",
              "null"
            ]
          },
          "species": {
            "type": [
              "array",
              "null"
            ]
          },
          "assay_methods": {
            "type": [
              "array",
              "null"
            ]
          },
          "assay_method_counts": {
            "type": [
              "array",
              "null"
            ]
          },
          "biosample_types": {
            "type": [
              "array",
              "null"
            ]
          },
          "actual_enrollment": {
            "type": [
              "integer",
              "null"
            ]
          },
          "age_range": {
            "type": [
              "string",
              "null"
            ]
          },
          "gender_included": {
            "type": [
              "string",
              "null"
            ]
          },
          "pubmed_ids": {
            "type": [
              "array",
              "null"
            ]
          },
          "doi": {
            "type": [
              "string",
              "null"
            ]
          },
          "program_name": {
            "type": [
              "array",
              "null"
            ]
          },
          "study_pi": {
            "type": [
              "array",
              "null"
            ]
          },
          "clinical_trial": {
            "type": [
              "string",
              "null"
            ]
          },
          "latest_data_release_date": {
            "type": [
              "string",
              "null"
            ]
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `ImmPort_search_studies`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("imm_port_search_studies", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.intact_get_complex_details` {#dnakit.databases.intact_get_complex_details}

Get detailed information about a specific protein complex by its complex AC (Complex Accession, e.g., 'CPX-915'). Returns comprehensive complex information including: complex AC, name, synonyms, systematic name, functions, properties, species, ligands, complex assemblies, participants (with interactor details), and associated diseases. To find complex ACs, use intact_get_interactions_by_complex to search for complexes first.

```text
dnakit.databases.intact_get_complex_details(*, complex_ac, format=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `complex_ac` | string | 必填（输入 Schema） | Complex AC (Complex Accession) in format &#x27;CPX-XXXXX&#x27; (e.g., &#x27;CPX-915&#x27;, &#x27;CPX-10081&#x27;). Get complex ACs from intact_get_interactions_by_complex search results. |
| `format` | string | None | 返回格式；允许值与服务默认值见当前接口说明。 {&quot;enum&quot;: &#91;&quot;json&quot;, &quot;xml&quot;&#93;, &quot;default&quot;: &quot;json&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "complex_ac": {
          "type": "string",
          "description": "Complex AC (Complex Accession) in format 'CPX-XXXXX' (e.g., 'CPX-915', 'CPX-10081'). Get complex ACs from intact_get_interactions_by_complex search results."
        },
        "format": {
          "type": "string",
          "enum": [
            "json",
            "xml"
          ],
          "default": "json"
        }
      },
      "required": [
        "complex_ac"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "description": "Complete complex details including complexAC, name, synonyms, systematicName, functions, properties, species, ligands, complexAssemblies, participants, diseases, and other metadata",
      "properties": {
        "complexAc": {
          "type": "string",
          "description": "Complex AC (e.g., 'CPX-915')"
        },
        "name": {
          "type": "string",
          "description": "Complex name"
        },
        "systematicName": {
          "type": "string",
          "description": "Systematic name of the complex"
        },
        "synonyms": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "Alternative names for the complex"
        },
        "functions": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "Functional descriptions of the complex"
        },
        "properties": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "Properties and characteristics of the complex"
        },
        "ac": {
          "type": "string",
          "description": "Additional accession identifier"
        },
        "species": {
          "type": "string",
          "description": "Organism information"
        },
        "ligands": {
          "type": "array",
          "description": "Associated ligands"
        },
        "complexAssemblies": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "Complex assembly information"
        },
        "participants": {
          "type": "array",
          "description": "Complex participants with interactor details, stoichiometry, and features"
        },
        "diseases": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "Associated diseases"
        },
        "crossReferences": {
          "type": "array",
          "description": "Cross-references to other databases (GO, PubMed, PDB, etc.)"
        },
        "institution": {
          "type": "string",
          "description": "Institution that curated the complex data"
        },
        "institutionURL": {
          "type": "string",
          "description": "URL of the institution"
        },
        "agonists": {
          "type": "array",
          "description": "Agonists associated with the complex"
        },
        "antagonists": {
          "type": "array",
          "description": "Antagonists associated with the complex"
        },
        "comments": {
          "type": "array",
          "description": "Additional comments about the complex"
        },
        "predictedComplex": {
          "type": "boolean",
          "description": "Whether the complex is predicted or experimentally verified"
        },
        "evidenceType": {
          "type": "object",
          "description": "Evidence type information including identifier, description, and confidence score",
          "properties": {
            "identifier": {
              "type": "string",
              "description": "Evidence type identifier (e.g., ECO code)"
            },
            "description": {
              "type": "string",
              "description": "Description of the evidence type"
            },
            "searchURL": {
              "type": "string",
              "description": "URL to search for evidence type details"
            },
            "confidenceScore": {
              "type": "integer",
              "description": "Confidence score for the evidence"
            }
          }
        },
        "releaseDates": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "Release dates for the complex data"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `intact_get_complex_details`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("intact_get_complex_details", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.intact_get_interaction_network` {#dnakit.databases.intact_get_interaction_network}

Get interaction network centered on a specific interactor. Uses EBI Search API (IntAct domain) for reliable access. Accepts IntAct identifiers (e.g., 'EBI-1004115'), UniProt IDs (e.g., 'P04637'), or gene names (e.g., 'BRCA1'). Returns interactions with interaction IDs in the 'interaction_ids' field. These IDs can be used with intact_get_interaction_details for detailed information.

```text
dnakit.databases.intact_get_interaction_network(*, identifier=None, uniprot_id=None, protein_id=None, gene_symbol=None, gene_name=None, depth=None, limit=None, size=None, format=None, protein_name=None, protein=None, gene=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `identifier` | string | None | IntAct identifier, UniProt ID, or gene name. Aliases: uniprot_id, protein_id, gene_symbol, gene_name. |
| `uniprot_id` | string | None | Alias for identifier: UniProt accession (e.g., &#x27;P04637&#x27;). |
| `protein_id` | string | None | Alias for identifier: protein identifier. |
| `gene_symbol` | string | None | Alias for identifier: gene symbol (e.g., &#x27;BRCA1&#x27;). |
| `gene_name` | string | None | Alias for identifier: gene name. |
| `depth` | integer | None | Network depth: 1 for direct interactions only, 2 for 2-hop network, etc. (default: 1) {&quot;default&quot;: 1, &quot;minimum&quot;: 1, &quot;maximum&quot;: 3} |
| `limit` | integer | None | Maximum number of interactions to return (default: 50, max: 200). Alias: size. {&quot;default&quot;: 50} |
| `size` | integer | None | Alias for limit. Maximum number of interactions to return (default: 50). {&quot;default&quot;: 50} |
| `format` | string | None | 返回格式；允许值与服务默认值见当前接口说明。 {&quot;enum&quot;: &#91;&quot;json&quot;, &quot;xml&quot;&#93;, &quot;default&quot;: &quot;json&quot;} |
| `protein_name` | string | None | Alias for gene_symbol/identifier. Common protein name (e.g., MDM2, TP53). |
| `protein` | string | None | Alias for identifier. Gene symbol or protein name (e.g., &#x27;TP53&#x27;, &#x27;BRCA1&#x27;). |
| `gene` | string | None | DNAKit optional alias for gene_symbol |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "identifier": {
          "type": "string",
          "description": "IntAct identifier, UniProt ID, or gene name. Aliases: uniprot_id, protein_id, gene_symbol, gene_name."
        },
        "uniprot_id": {
          "type": "string",
          "description": "Alias for identifier: UniProt accession (e.g., 'P04637')."
        },
        "protein_id": {
          "type": "string",
          "description": "Alias for identifier: protein identifier."
        },
        "gene_symbol": {
          "type": "string",
          "description": "Alias for identifier: gene symbol (e.g., 'BRCA1')."
        },
        "gene_name": {
          "type": "string",
          "description": "Alias for identifier: gene name."
        },
        "depth": {
          "type": "integer",
          "description": "Network depth: 1 for direct interactions only, 2 for 2-hop network, etc. (default: 1)",
          "default": 1,
          "minimum": 1,
          "maximum": 3
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of interactions to return (default: 50, max: 200). Alias: size.",
          "default": 50
        },
        "size": {
          "type": "integer",
          "description": "Alias for limit. Maximum number of interactions to return (default: 50).",
          "default": 50
        },
        "format": {
          "type": "string",
          "enum": [
            "json",
            "xml"
          ],
          "default": "json"
        },
        "protein_name": {
          "type": "string",
          "description": "Alias for gene_symbol/identifier. Common protein name (e.g., MDM2, TP53)."
        },
        "protein": {
          "type": "string",
          "description": "Alias for identifier. Gene symbol or protein name (e.g., 'TP53', 'BRCA1')."
        },
        "gene": {
          "type": "string",
          "description": "DNAKit optional alias for gene_symbol"
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
          "type": "array"
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
        "gene": "gene_symbol"
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `intact_get_interaction_network`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("intact_get_interaction_network", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.intact_get_interactions` {#dnakit.databases.intact_get_interactions}

Get all interactions involving a specific interactor (protein or gene) by identifier. Uses EBI Search API (IntAct domain) for reliable access. Returns list of interactions with interaction IDs (e.g., 'EBI-366083-EBI-366083') in the 'interaction_ids' field. These IDs can be used with intact_get_interaction_details or intact_get_interaction_network to get detailed information. Accepts UniProt IDs, gene names, or IntAct identifiers. Returns 25 interactions per call by default; compare 'count' with 'hitCount' in the response metadata and raise 'size' to retrieve the rest.

```text
dnakit.databases.intact_get_interactions(*, identifier, format=None, size=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `identifier` | string | 必填（输入 Schema） | IntAct identifier, UniProt ID, or gene name |
| `format` | string | None | Response format. Only &#x27;json&#x27; is produced by this endpoint. {&quot;enum&quot;: &#91;&quot;json&quot;&#93;, &quot;default&quot;: &quot;json&quot;} |
| `size` | integer &#124; null | None | Maximum number of interactions to return (default 25). IntAct entries commonly have far more -- the response reports the true total as &#x27;hitCount&#x27;, so raise &#x27;size&#x27; to that value to retrieve them all. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "identifier": {
          "type": "string",
          "description": "IntAct identifier, UniProt ID, or gene name"
        },
        "format": {
          "type": "string",
          "enum": [
            "json"
          ],
          "default": "json",
          "description": "Response format. Only 'json' is produced by this endpoint."
        },
        "size": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Maximum number of interactions to return (default 25). IntAct entries commonly have far more -- the response reports the true total as 'hitCount', so raise 'size' to that value to retrieve them all."
        }
      },
      "required": [
        "identifier"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "array",
      "description": "List of interactions involving the interactor"
    }
    ```

**执行与来源：** `network`；上游接口名 `intact_get_interactions`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("intact_get_interactions", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.intact_search_interactions` {#dnakit.databases.intact_search_interactions}

Search for interactions in IntAct database using EBI Search API (IntAct domain). Supports searching by protein names, gene names, interaction types, or other criteria. Use '*' as query to get all interactions (may be slow). Returns interaction entries with interaction IDs (e.g., 'EBI-366083-EBI-366083') in the 'interaction_ids' field. These IDs can be used with intact_get_interaction_details or intact_get_interaction_network to get detailed information.

```text
dnakit.databases.intact_search_interactions(*, query, first=None, max=None, format=None, species=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填（输入 Schema） | Search query string (e.g., protein name, gene name, interaction type) |
| `first` | integer | None | First result index for pagination (default: 0) {&quot;default&quot;: 0, &quot;minimum&quot;: 0} |
| `max` | integer | None | Maximum number of results to return (default: 25, max: 100) {&quot;default&quot;: 25, &quot;maximum&quot;: 100} |
| `format` | string | None | 返回格式；允许值与服务默认值见当前接口说明。 {&quot;enum&quot;: &#91;&quot;json&quot;, &quot;xml&quot;&#93;, &quot;default&quot;: &quot;json&quot;} |
| `species` | string | None | 物种名或 Taxonomy ID；类型与支持值取决于当前接口。 |
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
          "description": "Search query string (e.g., protein name, gene name, interaction type)"
        },
        "first": {
          "type": "integer",
          "description": "First result index for pagination (default: 0)",
          "default": 0,
          "minimum": 0
        },
        "max": {
          "type": "integer",
          "description": "Maximum number of results to return (default: 25, max: 100)",
          "default": 25,
          "maximum": 100
        },
        "format": {
          "type": "string",
          "enum": [
            "json",
            "xml"
          ],
          "default": "json"
        },
        "species": {
          "type": "string"
        }
      },
      "required": [
        "query"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "array",
      "description": "Search results with interactions"
    }
    ```

??? info "DNAKit 参数扩展规则"

    ```json
    {
      "parameters": {
        "species": {
          "type": "string"
        }
      },
      "handler": "species_filter"
    }
    ```

**执行与来源：** `network`；上游接口名 `intact_search_interactions`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("intact_search_interactions", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.inter_pro_get_entries_for_protein` {#dnakit.databases.inter_pro_get_entries_for_protein}

Get all InterPro domain and family entries annotated on a specific protein. Performs a reverse lookup from a UniProt protein accession to all InterPro entries (domains, families, homologous superfamilies, conserved sites) that have been detected on that protein. Example: P04637 (TP53) has 9 InterPro entries including IPR002117 (p53 family), IPR011615 (p53 DNA-binding domain), IPR010991 (p53 tetramerisation domain), IPR013872 (p53 transactivation domain).

```text
dnakit.databases.inter_pro_get_entries_for_protein(*, accession, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string | 必填（输入 Schema） | UniProt protein accession. Examples: &#x27;P04637&#x27; (TP53, 9 entries), &#x27;P00533&#x27; (EGFR), &#x27;P01308&#x27; (Insulin), &#x27;P38398&#x27; (BRCA1). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "accession": {
          "type": "string",
          "description": "UniProt protein accession. Examples: 'P04637' (TP53, 9 entries), 'P00533' (EGFR), 'P01308' (Insulin), 'P38398' (BRCA1)."
        }
      },
      "required": [
        "accession"
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
            "protein_accession": {
              "type": "string"
            },
            "total_entries": {
              "type": [
                "integer",
                "number"
              ]
            },
            "entries": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "accession": {
                    "type": "string"
                  },
                  "name": {
                    "type": "string"
                  },
                  "type": {
                    "type": "string"
                  },
                  "source_database": {
                    "type": "string"
                  }
                }
              }
            }
          },
          "required": [
            "protein_accession"
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

**执行与来源：** `network`；上游接口名 `InterPro_get_entries_for_protein`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("inter_pro_get_entries_for_protein", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.inter_pro_get_protein_domains` {#dnakit.databases.inter_pro_get_protein_domains}

Get protein domain annotations from InterPro database using UniProt protein ID. Returns domain families, signatures, and functional annotations.

```text
dnakit.databases.inter_pro_get_protein_domains(*, protein_id=None, accession=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `protein_id` | string | None | UniProt protein ID (e.g., P05067, Q9Y6K9) |
| `accession` | string | None | DNAKit optional alias for protein_id |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "protein_id": {
          "type": "string",
          "description": "UniProt protein ID (e.g., P05067, Q9Y6K9)"
        },
        "accession": {
          "type": "string",
          "description": "DNAKit optional alias for protein_id"
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
          "items": {
            "type": "object",
            "properties": {
              "metadata": {
                "type": "object",
                "properties": {
                  "accession": {
                    "type": "string"
                  },
                  "name": {
                    "type": "string"
                  },
                  "source_database": {
                    "type": "string"
                  },
                  "type": {
                    "type": "string"
                  },
                  "integrated": {
                    "type": "null"
                  },
                  "member_databases": {
                    "type": "object",
                    "properties": {
                      "panther": {
                        "type": "object",
                        "properties": {
                          "PTHR11447": {
                            "type": "string"
                          }
                        }
                      },
                      "prints": {
                        "type": "object",
                        "properties": {
                          "PR00386": {
                            "type": "string"
                          }
                        }
                      },
                      "ssf": {
                        "type": "object",
                        "properties": {
                          "SSF49417": {
                            "type": "string"
                          },
                          "SSF47719": {
                            "type": "string"
                          }
                        }
                      },
                      "pfam": {
                        "type": "object",
                        "properties": {
                          "PF07710": {
                            "type": "string"
                          },
                          "PF00870": {
                            "type": "string"
                          },
                          "PF08563": {
                            "type": "string"
                          },
                          "PF18521": {
                            "type": "string"
                          }
                        }
                      },
                      "cdd": {
                        "type": "object",
                        "properties": {
                          "cd08367": {
                            "type": "string"
                          }
                        }
                      },
                      "cathgene3d": {
                        "type": "object",
                        "properties": {
                          "G3DSA:2.60.40.720": {
                            "type": "string"
                          },
                          "G3DSA:4.10.170.10": {
                            "type": "string"
                          }
                        }
                      },
                      "prosite": {
                        "type": "object",
                        "properties": {
                          "PS00348": {
                            "type": "string"
                          }
                        }
                      }
                    }
                  },
                  "go_terms": {
                    "type": [
                      "array",
                      "null"
                    ],
                    "items": {
                      "type": "object",
                      "properties": {
                        "identifier": {
                          "type": "string"
                        },
                        "name": {
                          "type": "string"
                        },
                        "category": {
                          "type": "object",
                          "properties": {
                            "code": {
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
              },
              "proteins": {
                "type": "array",
                "items": {
                  "type": "object",
                  "properties": {
                    "accession": {
                      "type": "string"
                    },
                    "protein_length": {
                      "type": [
                        "integer",
                        "number"
                      ]
                    },
                    "source_database": {
                      "type": "string"
                    },
                    "organism": {
                      "type": "string"
                    },
                    "in_alphafold": {
                      "type": "boolean"
                    },
                    "in_bfvd": {
                      "type": "boolean"
                    },
                    "entry_protein_locations": {
                      "type": "array",
                      "items": {
                        "type": "object",
                        "properties": {
                          "fragments": {
                            "type": "array",
                            "items": {
                              "type": "object",
                              "properties": {
                                "start": {
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
                                "dc-status": {
                                  "type": "string"
                                }
                              }
                            }
                          },
                          "representative": {
                            "type": "boolean"
                          },
                          "model": {
                            "type": "null"
                          },
                          "score": {
                            "type": "null"
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
        "accession": "protein_id"
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `InterPro_get_protein_domains`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("inter_pro_get_protein_domains", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.kegg_find_genes` {#dnakit.databases.kegg_find_genes}

Find genes in KEGG database by keyword. Can search across all organisms or within a specific organism.

```text
dnakit.databases.kegg_find_genes(*, keyword=None, organism=None, query=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `keyword` | string | None | Search keyword for gene names or descriptions |
| `organism` | string | None | Organism code (e.g., &#x27;hsa&#x27; for human, &#x27;mmu&#x27; for mouse). Optional - searches all organisms if not specified {&quot;default&quot;: &quot;&quot;} |
| `query` | string | None | DNAKit optional alias for keyword |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "keyword": {
          "type": "string",
          "description": "Search keyword for gene names or descriptions"
        },
        "organism": {
          "type": "string",
          "description": "Organism code (e.g., 'hsa' for human, 'mmu' for mouse). Optional - searches all organisms if not specified",
          "default": ""
        },
        "query": {
          "type": "string",
          "description": "DNAKit optional alias for keyword"
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
          "items": {
            "type": "object",
            "properties": {
              "gene_id": {
                "type": "string",
                "description": "KEGG gene identifier"
              },
              "description": {
                "type": "string",
                "description": "Gene description"
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
      "aliases": {
        "query": "keyword"
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `kegg_find_genes`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("kegg_find_genes", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.kegg_get_gene_info` {#dnakit.databases.kegg_get_gene_info}

Get detailed gene information from KEGG by gene ID. Returns gene data including sequence, function, and pathway associations.

```text
dnakit.databases.kegg_get_gene_info(*, gene_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填（输入 Schema） | KEGG gene identifier (e.g., &#x27;hsa:348&#x27;, &#x27;hsa:3480&#x27;) |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene_id": {
          "type": "string",
          "description": "KEGG gene identifier (e.g., 'hsa:348', 'hsa:3480')"
        }
      },
      "required": [
        "gene_id"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "description": "KEGG gene information",
      "properties": {
        "status": {
          "type": "string",
          "description": "Status of the request"
        },
        "data": {
          "type": "object",
          "description": "Gene information",
          "properties": {
            "gene_id": {
              "type": "string",
              "description": "KEGG gene identifier"
            },
            "raw_data": {
              "type": "string",
              "description": "Raw gene data from KEGG"
            },
            "line_count": {
              "type": "integer",
              "description": "Number of data lines"
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

**执行与来源：** `network`；上游接口名 `kegg_get_gene_info`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("kegg_get_gene_info", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.kegg_get_gene_pathways` {#dnakit.databases.kegg_get_gene_pathways}

Get all KEGG pathways that a gene participates in. Returns pathway IDs and names for a given KEGG gene identifier. Example: hsa:7157 (human TP53) participates in 51 pathways including p53 signaling pathway (hsa04115), apoptosis (hsa04210), cell cycle (hsa04110), and 48 disease pathways (various cancers, viral infections). Requires KEGG organism:id format.

```text
dnakit.databases.kegg_get_gene_pathways(*, gene_id=None, gene_symbol=None, gene=None, organism=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | None | KEGG gene identifier in organism:id format. Examples: &#x27;hsa:7157&#x27; (human TP53), &#x27;hsa:672&#x27; (human BRCA1), &#x27;mmu:22059&#x27; (mouse Trp53). Use organism prefix: hsa=human, mmu=mouse, rno=rat, dme=fly, sce=yeast. |
| `gene_symbol` | string | None | Alias for gene_id: human gene symbol (e.g. &#x27;TP53&#x27;, &#x27;BRCA1&#x27;). Auto-resolved to KEGG gene ID using organism=&#x27;hsa&#x27;. |
| `gene` | string | None | Alias for gene_id: human gene symbol (e.g. &#x27;TP53&#x27;, &#x27;EGFR&#x27;). |
| `organism` | string | None | KEGG organism code for gene symbol resolution (default: &#x27;hsa&#x27; for human). Used only when gene_symbol is provided instead of gene_id. {&quot;default&quot;: &quot;hsa&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene_id": {
          "type": "string",
          "description": "KEGG gene identifier in organism:id format. Examples: 'hsa:7157' (human TP53), 'hsa:672' (human BRCA1), 'mmu:22059' (mouse Trp53). Use organism prefix: hsa=human, mmu=mouse, rno=rat, dme=fly, sce=yeast."
        },
        "gene_symbol": {
          "type": "string",
          "description": "Alias for gene_id: human gene symbol (e.g. 'TP53', 'BRCA1'). Auto-resolved to KEGG gene ID using organism='hsa'."
        },
        "gene": {
          "type": "string",
          "description": "Alias for gene_id: human gene symbol (e.g. 'TP53', 'EGFR')."
        },
        "organism": {
          "type": "string",
          "description": "KEGG organism code for gene symbol resolution (default: 'hsa' for human). Used only when gene_symbol is provided instead of gene_id.",
          "default": "hsa"
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
            "gene_id": {
              "type": "string",
              "description": "Queried KEGG gene ID"
            },
            "pathway_count": {
              "type": "integer",
              "description": "Number of pathways the gene participates in"
            },
            "pathways": {
              "type": "array",
              "description": "List of KEGG pathways with IDs and names",
              "items": {
                "type": "object",
                "properties": {
                  "pathway_id": {
                    "type": "string"
                  },
                  "pathway_name": {
                    "type": "string"
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

**执行与来源：** `network`；上游接口名 `KEGG_get_gene_pathways`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("kegg_get_gene_pathways", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.kegg_get_pathway_genes` {#dnakit.databases.kegg_get_pathway_genes}

Get all genes in a KEGG pathway. Returns gene IDs for all member genes of a specific pathway. Example: hsa04115 (p53 signaling pathway) contains genes including CDK2 (hsa:1017), CDK4 (hsa:1019), CDK6 (hsa:1021), CDKN1A (hsa:1026), CDKN2A (hsa:1029), and many more involved in p53-mediated cell cycle arrest, apoptosis, and senescence.

```text
dnakit.databases.kegg_get_pathway_genes(*, pathway_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `pathway_id` | string | 必填（输入 Schema） | KEGG pathway identifier (e.g., &#x27;hsa04115&#x27; for p53 signaling pathway, &#x27;hsa04110&#x27; for cell cycle, &#x27;hsa04210&#x27; for apoptosis). Use organism prefix: hsa=human, mmu=mouse. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "pathway_id": {
          "type": "string",
          "description": "KEGG pathway identifier (e.g., 'hsa04115' for p53 signaling pathway, 'hsa04110' for cell cycle, 'hsa04210' for apoptosis). Use organism prefix: hsa=human, mmu=mouse."
        }
      },
      "required": [
        "pathway_id"
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
            "pathway_id": {
              "type": "string",
              "description": "KEGG pathway ID"
            },
            "pathway_name": {
              "type": "string",
              "description": "Pathway name and species"
            },
            "gene_count": {
              "type": "integer",
              "description": "Number of genes in the pathway"
            },
            "genes": {
              "type": "array",
              "description": "List of KEGG gene IDs in the pathway",
              "items": {
                "type": "string"
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

**执行与来源：** `network`；上游接口名 `KEGG_get_pathway_genes`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("kegg_get_pathway_genes", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.kegg_get_pathway_info` {#dnakit.databases.kegg_get_pathway_info}

Get detailed pathway information from KEGG by pathway ID. Returns pathway data including genes, compounds, and reactions.

```text
dnakit.databases.kegg_get_pathway_info(*, pathway_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `pathway_id` | string | 必填（输入 Schema） | KEGG pathway identifier (e.g., &#x27;hsa00010&#x27;, &#x27;path:hsa00010&#x27;) |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "pathway_id": {
          "type": "string",
          "description": "KEGG pathway identifier (e.g., 'hsa00010', 'path:hsa00010')"
        }
      },
      "required": [
        "pathway_id"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "description": "KEGG pathway information",
      "properties": {
        "status": {
          "type": "string",
          "description": "Status of the request"
        },
        "data": {
          "type": "object",
          "description": "Pathway information",
          "properties": {
            "pathway_id": {
              "type": "string",
              "description": "KEGG pathway identifier"
            },
            "raw_data": {
              "type": "string",
              "description": "Raw pathway data from KEGG"
            },
            "line_count": {
              "type": "integer",
              "description": "Number of data lines"
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

**执行与来源：** `network`；上游接口名 `kegg_get_pathway_info`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("kegg_get_pathway_info", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.kegg_search_pathway` {#dnakit.databases.kegg_search_pathway}

Search KEGG pathways by keyword. Returns pathway IDs and descriptions matching the search term.

```text
dnakit.databases.kegg_search_pathway(*, keyword, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `keyword` | string | 必填（输入 Schema） | Search keyword for pathway names or descriptions (e.g., &#x27;diabetes&#x27;, &#x27;metabolism&#x27;) |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "keyword": {
          "type": "string",
          "description": "Search keyword for pathway names or descriptions (e.g., 'diabetes', 'metabolism')"
        }
      },
      "required": [
        "keyword"
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
              "pathway_id": {
                "type": "string",
                "description": "KEGG pathway identifier"
              },
              "description": {
                "type": "string",
                "description": "Pathway description"
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

**执行与来源：** `network`；上游接口名 `kegg_search_pathway`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("kegg_search_pathway", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ncbiprotein_get_summary` {#dnakit.databases.ncbiprotein_get_summary}

Get protein record summary from the NCBI Protein database by GI number or accession using the E-utilities ESummary API. Returns protein title, organism, accession version, sequence length, molecular type, source database, and creation/update dates. NCBI Protein aggregates sequences from GenBank, RefSeq, Swiss-Prot, PDB, and other sources. Use for resolving NCBI protein GI numbers to names, accessions, and species. Example: GI 119395750 returns 'keratin, type II cytoskeletal 1 [Homo sapiens]' with accession NP_006112.3 and length 644 aa.

```text
dnakit.databases.ncbiprotein_get_summary(*, id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `id` | string | 必填（输入 Schema） | NCBI Protein ID (GI number or accession). Can be a single ID or comma-separated list. Examples: &#x27;119395750&#x27; (keratin 1, NP_006112), &#x27;4507879&#x27; (TP53, NP_000537), &#x27;4503747&#x27; (EGFR, NP_005219). GI numbers can be obtained from NCBI Gene genomic info or from BLAST results. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "id": {
          "type": "string",
          "description": "NCBI Protein ID (GI number or accession). Can be a single ID or comma-separated list. Examples: '119395750' (keratin 1, NP_006112), '4507879' (TP53, NP_000537), '4503747' (EGFR, NP_005219). GI numbers can be obtained from NCBI Gene genomic info or from BLAST results."
        }
      },
      "required": [
        "id"
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
            "result": {
              "type": "object",
              "properties": {
                "uids": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  },
                  "description": "List of returned protein IDs (GI numbers)"
                }
              },
              "additionalProperties": {
                "type": "object",
                "properties": {
                  "uid": {
                    "type": "string",
                    "description": "GI number"
                  },
                  "caption": {
                    "type": "string",
                    "description": "Accession number (e.g., NP_006112)"
                  },
                  "title": {
                    "type": "string",
                    "description": "Protein name with organism"
                  },
                  "extra": {
                    "type": "string",
                    "description": "Additional identifier info"
                  },
                  "gi": {
                    "type": "integer",
                    "description": "GI number (integer form)"
                  },
                  "taxid": {
                    "type": "integer",
                    "description": "NCBI Taxonomy ID"
                  },
                  "slen": {
                    "type": "integer",
                    "description": "Sequence length in amino acids"
                  },
                  "moltype": {
                    "type": "string",
                    "description": "Molecule type (aa for amino acid)"
                  },
                  "sourcedb": {
                    "type": "string",
                    "description": "Source database (refseq, insd, etc.)"
                  },
                  "organism": {
                    "type": "string",
                    "description": "Scientific name of the organism"
                  },
                  "accessionversion": {
                    "type": "string",
                    "description": "Accession with version (e.g., NP_006112.3)"
                  },
                  "createdate": {
                    "type": "string",
                    "description": "Record creation date"
                  },
                  "updatedate": {
                    "type": "string",
                    "description": "Record last update date"
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

**执行与来源：** `network`；上游接口名 `NCBIProtein_get_summary`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ncbiprotein_get_summary", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.panther_enrichment` {#dnakit.databases.panther_enrichment}

Perform gene set enrichment (overrepresentation) analysis using PANTHER. Given a gene list, identifies statistically enriched GO terms, PANTHER pathways, Reactome pathways, or protein classes using Fisher's exact test with FDR correction. PANTHER is one of the most widely used tools for functional enrichment in bioinformatics. Example: gene list 'TP53,BRCA1,EGFR,KRAS,PTEN' enriched for GO biological_process returns 'regulation of mitotic cell cycle phase transition' (FDR=0.003).

```text
dnakit.databases.panther_enrichment(*, gene_list, organism=None, annotation_dataset=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_list` | string | 必填（输入 Schema） | Comma-separated list of gene symbols, UniProt IDs, or Ensembl IDs. Example: &#x27;TP53,BRCA1,EGFR,KRAS,PTEN,APC,RB1&#x27;. |
| `organism` | integer &#124; null | None | NCBI taxonomy ID. Default: 9606 (human). Others: 10090 (mouse), 7227 (fly), 559292 (yeast). |
| `annotation_dataset` | string &#124; null | None | Annotation dataset for enrichment. Options: &#x27;GO:0008150&#x27; (biological_process, default), &#x27;GO:0003674&#x27; (molecular_function), &#x27;GO:0005575&#x27; (cellular_component), &#x27;ANNOT_TYPE_ID_PANTHER_PATHWAY&#x27; (PANTHER pathway), &#x27;ANNOT_TYPE_ID_REACTOME_PATHWAY&#x27; (Reactome pathway), &#x27;ANNOT_TYPE_ID_PANTHER_PC&#x27; (protein class). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene_list": {
          "type": "string",
          "description": "Comma-separated list of gene symbols, UniProt IDs, or Ensembl IDs. Example: 'TP53,BRCA1,EGFR,KRAS,PTEN,APC,RB1'."
        },
        "organism": {
          "type": [
            "integer",
            "null"
          ],
          "description": "NCBI taxonomy ID. Default: 9606 (human). Others: 10090 (mouse), 7227 (fly), 559292 (yeast)."
        },
        "annotation_dataset": {
          "type": [
            "string",
            "null"
          ],
          "description": "Annotation dataset for enrichment. Options: 'GO:0008150' (biological_process, default), 'GO:0003674' (molecular_function), 'GO:0005575' (cellular_component), 'ANNOT_TYPE_ID_PANTHER_PATHWAY' (PANTHER pathway), 'ANNOT_TYPE_ID_REACTOME_PATHWAY' (Reactome pathway), 'ANNOT_TYPE_ID_PANTHER_PC' (protein class)."
        }
      },
      "required": [
        "gene_list"
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
            "gene_list": {
              "type": "string"
            },
            "organism": {
              "type": "integer"
            },
            "annotation_dataset": {
              "type": "string"
            },
            "result_count": {
              "type": "integer"
            },
            "enriched_terms": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "term_id": {
                    "type": "string"
                  },
                  "term_label": {
                    "type": "string"
                  },
                  "number_in_list": {
                    "type": "integer"
                  },
                  "number_in_reference": {
                    "type": "integer"
                  },
                  "expected": {
                    "type": "number"
                  },
                  "fold_enrichment": {
                    "type": "number"
                  },
                  "pvalue": {
                    "type": "number"
                  },
                  "fdr": {
                    "type": "number"
                  },
                  "direction": {
                    "type": "string"
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

**执行与来源：** `network`；上游接口名 `PANTHER_enrichment`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("panther_enrichment", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.panther_ortholog` {#dnakit.databases.panther_ortholog}

Find orthologs of a gene across species using PANTHER evolutionary classification. Returns EVERY matching ortholog in 'mappings' (count in 'total_mappings'), each with its persistent ID, gene symbol, and ortholog type (LDO=Least Diverged Ortholog, O=ortholog, P=paralog, X=xenolog, LDX=least diverged xenolog); 'mapping' repeats the first entry for backward compatibility and is not the whole answer. A gene often has several orthologs in one target species -- human ABCB1 maps to both mouse Abcb1a (LDO) and Abcb1b (O) -- so read 'mappings', especially with ortholog_type='O'. PANTHER classifies orthologs by evolutionary divergence across 144 fully sequenced genomes. Example: P04637 (human TP53) maps to Tp53 (mouse P02340) as LDO.

```text
dnakit.databases.panther_ortholog(*, gene_id, organism=None, target_organism=None, ortholog_type=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填（输入 Schema） | Gene identifier - UniProt accession, Ensembl ID, or gene symbol. Examples: &#x27;P04637&#x27; (TP53), &#x27;P38398&#x27; (BRCA1). |
| `organism` | integer &#124; null | None | Source organism NCBI taxonomy ID. Default: 9606 (human). |
| `target_organism` | integer &#124; null | None | Target organism NCBI taxonomy ID. Default: 10090 (mouse). Others: 10116 (rat), 7227 (fly), 6239 (worm). |
| `ortholog_type` | string &#124; null | None | Ortholog type filter. Options: &#x27;LDO&#x27; (least diverged ortholog, default), &#x27;O&#x27; (all orthologs), &#x27;P&#x27; (paralogs), &#x27;X&#x27; (xenologs), &#x27;LDX&#x27; (least diverged xenolog), &#x27;all&#x27; (all types). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene_id": {
          "type": "string",
          "description": "Gene identifier - UniProt accession, Ensembl ID, or gene symbol. Examples: 'P04637' (TP53), 'P38398' (BRCA1)."
        },
        "organism": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Source organism NCBI taxonomy ID. Default: 9606 (human)."
        },
        "target_organism": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Target organism NCBI taxonomy ID. Default: 10090 (mouse). Others: 10116 (rat), 7227 (fly), 6239 (worm)."
        },
        "ortholog_type": {
          "type": [
            "string",
            "null"
          ],
          "description": "Ortholog type filter. Options: 'LDO' (least diverged ortholog, default), 'O' (all orthologs), 'P' (paralogs), 'X' (xenologs), 'LDX' (least diverged xenolog), 'all' (all types)."
        }
      },
      "required": [
        "gene_id"
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
            "source_organism": {
              "type": "integer"
            },
            "target_organism": {
              "type": "integer"
            },
            "ortholog_type": {
              "type": "string"
            },
            "mapping": {
              "type": [
                "object",
                "null"
              ],
              "properties": {
                "source_gene": {
                  "type": "string"
                },
                "target_gene": {
                  "type": "string"
                },
                "target_gene_symbol": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "ortholog_type": {
                  "type": "string"
                },
                "persistent_id": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "target_persistent_id": {
                  "type": [
                    "string",
                    "null"
                  ]
                }
              }
            },
            "mappings": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "source_gene": {
                    "type": "string"
                  },
                  "target_gene": {
                    "type": "string"
                  },
                  "target_gene_symbol": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "ortholog_type": {
                    "type": "string"
                  },
                  "persistent_id": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "target_persistent_id": {
                    "type": [
                      "string",
                      "null"
                    ]
                  }
                }
              }
            },
            "total_mappings": {
              "type": "integer"
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

**执行与来源：** `network`；上游接口名 `PANTHER_ortholog`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("panther_ortholog", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.pdb_get_structure` {#dnakit.databases.pdb_get_structure}

Retrieve a bounded RCSB PDB coordinate file as text, without writing it to disk.

    Service: https://www.rcsb.org/docs/programmatic-access/file-download-services
    

```text
dnakit.databases.pdb_get_structure(*, pdb_id, format='cif', max_bytes=10000000, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `pdb_id` | string | 必填（输入 Schema） | PDB 结构编号，支持 4 字符编号或 pdb_ 加 8 字符编号。 |
| `format` | string | &#x27;cif&#x27; | 返回格式；允许值与服务默认值见当前接口说明。 {&quot;default&quot;: &quot;cif&quot;} |
| `max_bytes` | integer | 10000000 | 坐标文件大小上限，1–100000000 字节。 {&quot;default&quot;: 10000000} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

字典含 pdb_id、format、data（坐标文件文本）、url、sha256；不写入磁盘。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "pdb_id": {
          "type": "string"
        },
        "format": {
          "type": "string",
          "default": "cif"
        },
        "max_bytes": {
          "type": "integer",
          "default": 10000000
        }
      },
      "required": [
        "pdb_id"
      ]
    }
    ```

**执行与来源：** `network`；上游接口名 `PDB_get_structure`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("pdb_get_structure", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.pdbe_get_uniprot_mappings` {#dnakit.databases.pdbe_get_uniprot_mappings}

Get UniProt-to-PDB chain mappings for a PDB structure from the PDBe Graph API. For each UniProt accession represented in the structure, returns all chain mappings with residue-level alignment (UniProt start/end, PDB start/end, identity, coverage). Essential for understanding which proteins are in a crystal structure and how they map to UniProt reference sequences. Example: 4HHB returns P69905 (HBA_HUMAN) mapped to chains A,C and P68871 (HBB_HUMAN) mapped to chains B,D.

```text
dnakit.databases.pdbe_get_uniprot_mappings(*, pdb_id=None, uniprot_id=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `pdb_id` | string | None | PDB identifier (4-character code). Examples: &#x27;4hhb&#x27; (hemoglobin), &#x27;1cbs&#x27; (CRABP2), &#x27;6lu7&#x27; (SARS-CoV-2 main protease). |
| `uniprot_id` | string | None | UniProt accession；反向结构映射时不要同时传 pdb_id。 |
| `limit` | integer | None | 当前请求的返回数量上限；不是全库命中总数。 {&quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "pdb_id": {
          "type": "string",
          "description": "PDB identifier (4-character code). Examples: '4hhb' (hemoglobin), '1cbs' (CRABP2), '6lu7' (SARS-CoV-2 main protease)."
        },
        "uniprot_id": {
          "type": "string"
        },
        "limit": {
          "type": "integer",
          "minimum": 1,
          "maximum": 100
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
          "description": "Keyed by PDB ID, then 'UniProt', then UniProt accession",
          "additionalProperties": {
            "type": "object",
            "properties": {
              "UniProt": {
                "type": "object",
                "additionalProperties": {
                  "type": "object",
                  "properties": {
                    "name": {
                      "type": "string",
                      "description": "UniProt entry name"
                    },
                    "identifier": {
                      "type": "string"
                    },
                    "mappings": {
                      "type": "array",
                      "items": {
                        "type": "object",
                        "properties": {
                          "entity_id": {
                            "type": "integer"
                          },
                          "chain_id": {
                            "type": "string"
                          },
                          "struct_asym_id": {
                            "type": "string"
                          },
                          "unp_start": {
                            "type": "integer",
                            "description": "UniProt residue start"
                          },
                          "unp_end": {
                            "type": "integer",
                            "description": "UniProt residue end"
                          },
                          "identity": {
                            "type": "number",
                            "description": "Sequence identity (0-1)"
                          },
                          "coverage": {
                            "type": "number",
                            "description": "Fraction of UniProt sequence covered"
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
        "uniprot_id": {
          "type": "string"
        },
        "limit": {
          "type": "integer",
          "minimum": 1,
          "maximum": 100
        }
      },
      "optional": [
        "pdb_id"
      ]
    }
    ```

**执行与来源：** `network`；上游接口名 `PDBe_get_uniprot_mappings`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("pdbe_get_uniprot_mappings", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.pdbe_sifts_get_best_structures` {#dnakit.databases.pdbe_sifts_get_best_structures}

Get the best available PDB structures for a UniProt protein, ranked by coverage and resolution from PDBe SIFTS. Returns PDB IDs with chain mapping, resolution, experimental method (X-ray, cryo-EM, NMR), and UniProt sequence coverage. Essential for structural biologists to find the most suitable structure for a protein of interest. Example: P04637 (TP53) returns 641 structure entries, P00533 (EGFR) returns 633 entries, with X-ray and cryo-EM structures ranked by quality. Returns 50 entries per call by default; use 'limit' (max 500) and 'offset' to page through the full ranked list, whose size is reported in the response.

```text
dnakit.databases.pdbe_sifts_get_best_structures(*, uniprot_accession, limit=None, offset=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `uniprot_accession` | string | 必填（输入 Schema） | UniProt accession for the protein. Examples: &#x27;P04637&#x27; (TP53), &#x27;P00533&#x27; (EGFR), &#x27;P38398&#x27; (BRCA1), &#x27;P01308&#x27; (Insulin). |
| `limit` | integer &#124; null | None | Maximum number of entries to return (default 50, max 500). Proteins such as P04637 (TP53) have hundreds of structure-chain entries; raise this or use &#x27;offset&#x27; to reach past the first page. {&quot;default&quot;: 50} |
| `offset` | integer &#124; null | None | Zero-based index of the first entry to return, for paging through the ranked list (e.g. offset=50 returns entries 51-100). Default: 0. {&quot;default&quot;: 0} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "uniprot_accession": {
          "type": "string",
          "description": "UniProt accession for the protein. Examples: 'P04637' (TP53), 'P00533' (EGFR), 'P38398' (BRCA1), 'P01308' (Insulin)."
        },
        "limit": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Maximum number of entries to return (default 50, max 500). Proteins such as P04637 (TP53) have hundreds of structure-chain entries; raise this or use 'offset' to reach past the first page.",
          "default": 50
        },
        "offset": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Zero-based index of the first entry to return, for paging through the ranked list (e.g. offset=50 returns entries 51-100). Default: 0.",
          "default": 0
        }
      },
      "required": [
        "uniprot_accession"
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
            "uniprot_accession": {
              "type": "string",
              "description": "UniProt accession queried"
            },
            "structures": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "pdb_id": {
                    "type": [
                      "string",
                      "null"
                    ],
                    "description": "PDB entry ID (4 characters)"
                  },
                  "chain_id": {
                    "type": [
                      "string",
                      "null"
                    ],
                    "description": "Chain ID in PDB entry"
                  },
                  "uniprot_start": {
                    "type": [
                      "integer",
                      "null"
                    ],
                    "description": "UniProt sequence start covered"
                  },
                  "uniprot_end": {
                    "type": [
                      "integer",
                      "null"
                    ],
                    "description": "UniProt sequence end covered"
                  },
                  "resolution": {
                    "type": [
                      "number",
                      "null"
                    ],
                    "description": "Structure resolution in Angstroms"
                  },
                  "experimental_method": {
                    "type": [
                      "string",
                      "null"
                    ],
                    "description": "Method (X-ray diffraction, Electron Microscopy, NMR)"
                  },
                  "coverage": {
                    "type": [
                      "number",
                      "null"
                    ],
                    "description": "Fraction of UniProt sequence covered"
                  },
                  "tax_id": {
                    "type": [
                      "integer",
                      "null"
                    ],
                    "description": "Taxonomy ID of source organism"
                  }
                }
              },
              "description": "Ranked PDB structures (this page)"
            },
            "total_structures": {
              "type": "integer",
              "description": "Total available structure-chain entries"
            },
            "returned_structures": {
              "type": "integer",
              "description": "Number of entries in this page of results"
            },
            "offset": {
              "type": "integer",
              "description": "Zero-based index of the first entry returned"
            },
            "note": {
              "type": "string",
              "description": "Present when more entries remain; explains how to page to them"
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

**执行与来源：** `network`；上游接口名 `PDBeSIFTS_get_best_structures`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("pdbe_sifts_get_best_structures", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.pfam_get_protein_annotations` {#dnakit.databases.pfam_get_protein_annotations}

Get all Pfam domain annotations for a specific protein by its UniProt accession. Returns every Pfam domain found in the protein with exact start/end residue positions, E-values (scores), short names, and integrated InterPro entries. This provides the complete Pfam domain architecture of a protein. For example, P04637 (TP53) has 4 Pfam domains: PF08563 (P53_TAD, 6-30), PF18521 (P53_TAD2, 35-59), PF00870 (P53, 100-288), PF07710 (P53_tetramer, 319-357). Use this to understand the domain composition of any protein.

```text
dnakit.databases.pfam_get_protein_annotations(*, accession, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string | 必填（输入 Schema） | UniProt protein accession. Examples: &#x27;P04637&#x27; (TP53, 4 domains), &#x27;P00533&#x27; (EGFR), &#x27;P68431&#x27; (Histone H3.1), &#x27;P38398&#x27; (BRCA1), &#x27;P01308&#x27; (Insulin). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "accession": {
          "type": "string",
          "description": "UniProt protein accession. Examples: 'P04637' (TP53, 4 domains), 'P00533' (EGFR), 'P68431' (Histone H3.1), 'P38398' (BRCA1), 'P01308' (Insulin)."
        }
      },
      "required": [
        "accession"
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
            "accession": {
              "type": "string",
              "description": "UniProt accession queried"
            },
            "protein_length": {
              "type": [
                "integer",
                "null"
              ],
              "description": "Total protein length"
            },
            "domain_count": {
              "type": "integer",
              "description": "Number of Pfam domains found"
            },
            "domains": {
              "type": "array",
              "description": "Pfam domains with positions, sorted by start position",
              "items": {
                "type": "object",
                "properties": {
                  "pfam_accession": {
                    "type": "string"
                  },
                  "name": {
                    "type": "string"
                  },
                  "short_name": {
                    "type": "string"
                  },
                  "type": {
                    "type": "string"
                  },
                  "integrated_interpro": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "start": {
                    "type": "integer"
                  },
                  "end": {
                    "type": "integer"
                  },
                  "score": {
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

**执行与来源：** `network`；上游接口名 `Pfam_get_protein_annotations`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("pfam_get_protein_annotations", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.plant_reactome_get_pathway` {#dnakit.databases.plant_reactome_get_pathway}

Get detailed information about a specific Plant Reactome pathway by stable ID. Returns pathway name, species, sub-events (reactions), compartments, GO annotations, and cross-references. Example: 'R-OSA-1119312' returns the Photorespiration pathway in Oryza sativa.

```text
dnakit.databases.plant_reactome_get_pathway(*, pathway_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `pathway_id` | string | 必填（输入 Schema） | Plant Reactome stable identifier. Format: R-{species_code}-{number}. Examples: &#x27;R-OSA-1119312&#x27; (Photorespiration in rice), &#x27;R-OSA-1119519&#x27; (Calvin cycle in rice). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "pathway_id": {
          "type": "string",
          "description": "Plant Reactome stable identifier. Format: R-{species_code}-{number}. Examples: 'R-OSA-1119312' (Photorespiration in rice), 'R-OSA-1119519' (Calvin cycle in rice)."
        }
      },
      "required": [
        "pathway_id"
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
            "stId": {
              "type": "string"
            },
            "displayName": {
              "type": "string"
            },
            "speciesName": {
              "type": [
                "string",
                "null"
              ]
            },
            "isInDisease": {
              "type": [
                "boolean",
                "null"
              ]
            },
            "isInferred": {
              "type": [
                "boolean",
                "null"
              ]
            },
            "hasEvent": {
              "type": [
                "array",
                "null"
              ],
              "description": "Sub-events (reactions and sub-pathways)"
            },
            "compartment": {
              "type": [
                "array",
                "null"
              ]
            },
            "goBiologicalProcess": {
              "type": [
                "object",
                "null"
              ]
            },
            "crossReference": {
              "type": [
                "array",
                "null"
              ]
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

**执行与来源：** `network`；上游接口名 `PlantReactome_get_pathway`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("plant_reactome_get_pathway", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.plant_reactome_list_species` {#dnakit.databases.plant_reactome_list_species}

List all plant species available in Plant Reactome with their taxonomy IDs. Returns the complete catalog of 140+ plant species including model organisms and major crop species. Useful for discovering which species have curated pathway data available.

```text
dnakit.databases.plant_reactome_list_species(*, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {},
      "required": []
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
              "dbId": {
                "type": "integer"
              },
              "displayName": {
                "type": "string",
                "description": "Species name"
              },
              "taxId": {
                "type": [
                  "string",
                  "integer"
                ]
              },
              "abbreviation": {
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

**执行与来源：** `network`；上游接口名 `PlantReactome_list_species`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("plant_reactome_list_species", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.plant_reactome_search_pathways` {#dnakit.databases.plant_reactome_search_pathways}

Search for plant biological pathways in Plant Reactome (Gramene). Plant Reactome is a curated resource for plant metabolic and regulatory pathways, covering model and crop species including rice (Oryza sativa), Arabidopsis thaliana, maize, wheat, and 140+ plant species. Returns pathways matching the query with their stable IDs and species. Example: 'photosynthesis' finds photosynthesis-related pathways.

```text
dnakit.databases.plant_reactome_search_pathways(*, query, species=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填（输入 Schema） | Search query for plant pathways. Examples: &#x27;photosynthesis&#x27;, &#x27;Calvin cycle&#x27;, &#x27;flavonoid biosynthesis&#x27;, &#x27;nitrogen fixation&#x27;, &#x27;auxin signaling&#x27;. |
| `species` | string &#124; null | None | Species name to filter results. Examples: &#x27;Oryza sativa&#x27;, &#x27;Arabidopsis thaliana&#x27;, &#x27;Zea mays&#x27;. Default: all species. |
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
          "description": "Search query for plant pathways. Examples: 'photosynthesis', 'Calvin cycle', 'flavonoid biosynthesis', 'nitrogen fixation', 'auxin signaling'."
        },
        "species": {
          "type": [
            "string",
            "null"
          ],
          "description": "Species name to filter results. Examples: 'Oryza sativa', 'Arabidopsis thaliana', 'Zea mays'. Default: all species."
        }
      },
      "required": [
        "query"
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
              "stId": {
                "type": "string",
                "description": "Stable identifier (e.g., R-OSA-1119312)"
              },
              "name": {
                "type": "string"
              },
              "species": {
                "type": [
                  "string",
                  "array",
                  "null"
                ]
              },
              "exact_type": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Entity type (Pathway, Reaction, etc.)"
              },
              "compartment_names": {
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

**执行与来源：** `network`；上游接口名 `PlantReactome_search_pathways`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("plant_reactome_search_pathways", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.proteins_api_get_features` {#dnakit.databases.proteins_api_get_features}

Get detailed sequence features for a protein, including domains, regions, motifs, modified residues, and other annotated features. Supports batch operations: pass a list of accessions or comma-separated string to retrieve multiple proteins in parallel. Note: This tool extracts features from the main protein endpoint (the separate features endpoint is not available). Returns all sequence feature annotations for the protein.

```text
dnakit.databases.proteins_api_get_features(*, accession, format=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string &#124; array | 必填（输入 Schema） | UniProt protein accession(s). Can be a single accession (e.g., &#x27;P05067&#x27;), comma-separated string (e.g., &#x27;P05067,P04637&#x27;), or list (e.g., &#91;&#x27;P05067&#x27;, &#x27;P04637&#x27;&#93;). Use UniProt search or ebi_search_domain with domain=&#x27;uniprot&#x27; to find valid accessions. |
| `format` | string | None | Response format. JSON is recommended for most use cases. {&quot;enum&quot;: &#91;&quot;json&quot;, &quot;xml&quot;&#93;, &quot;default&quot;: &quot;json&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "accession": {
          "type": [
            "string",
            "array"
          ],
          "description": "UniProt protein accession(s). Can be a single accession (e.g., 'P05067'), comma-separated string (e.g., 'P05067,P04637'), or list (e.g., ['P05067', 'P04637']). Use UniProt search or ebi_search_domain with domain='uniprot' to find valid accessions."
        },
        "format": {
          "type": "string",
          "enum": [
            "json",
            "xml"
          ],
          "default": "json",
          "description": "Response format. JSON is recommended for most use cases."
        }
      },
      "required": [
        "accession"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "array",
      "description": "List of sequence feature annotations for the requested protein."
    }
    ```

**执行与来源：** `network`；上游接口名 `proteins_api_get_features`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("proteins_api_get_features", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.proteins_api_search` {#dnakit.databases.proteins_api_search}

Search proteins in Proteins API by gene name, protein name, or accession. Defaults to human proteins (taxid=9606) for gene/name searches; use organism param to override. Gene names (e.g., 'BRCA1') search by gene parameter, accessions (e.g., 'P05067') search by accession parameter, and longer queries (e.g., 'insulin') search by protein parameter. Returns list of matching proteins with full annotations.

```text
dnakit.databases.proteins_api_search(*, query=None, size=None, offset=None, format=None, gene_symbol=None, gene=None, reviewed=None, organism=None, taxid=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | None | Search query. Can be: gene name (e.g., &#x27;BRCA1&#x27;, &#x27;TP53&#x27;), protein name (e.g., &#x27;insulin&#x27;, &#x27;hemoglobin&#x27;), or UniProt accession (e.g., &#x27;P05067&#x27;). The tool automatically determines the appropriate search parameter. |
| `size` | integer | None | Maximum number of results to return (default: 25, max: 100) {&quot;default&quot;: 25, &quot;maximum&quot;: 100} |
| `offset` | integer | None | Offset for pagination (default: 0) {&quot;default&quot;: 0} |
| `format` | string | None | Response format. JSON is recommended for most use cases. {&quot;enum&quot;: &#91;&quot;json&quot;, &quot;xml&quot;&#93;, &quot;default&quot;: &quot;json&quot;} |
| `gene_symbol` | string | None | Alias for query. Gene symbol to search (e.g., TP53, BRCA1). |
| `gene` | string | None | Alias for query. Gene name to search. |
| `reviewed` | boolean | None | If true, return only reviewed Swiss-Prot entries. If false or omitted, includes unreviewed TrEMBL entries. |
| `organism` | string | None | Restrict the search to an organism by NAME (e.g. &#x27;Klebsiella pneumoniae&#x27;, &#x27;Escherichia coli&#x27;); use &#x27;taxid&#x27; for numeric NCBI taxonomy ids. A purely numeric value here (e.g. &#x27;9606&#x27;) is auto-corrected to a taxid lookup and the response reports it in &#x27;normalization_note&#x27;. By default a gene/protein-name search is biased toward human (taxid 9606); set this for non-human proteins. If omitted and no human match is found, the search auto-widens to all organisms. |
| `taxid` | string | None | Restrict the search to an NCBI taxonomy id (e.g. &#x27;573&#x27; for Klebsiella pneumoniae, &#x27;562&#x27; for E. coli). Alternative to &#x27;organism&#x27; for precise scoping. |
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
          "description": "Search query. Can be: gene name (e.g., 'BRCA1', 'TP53'), protein name (e.g., 'insulin', 'hemoglobin'), or UniProt accession (e.g., 'P05067'). The tool automatically determines the appropriate search parameter."
        },
        "size": {
          "type": "integer",
          "description": "Maximum number of results to return (default: 25, max: 100)",
          "default": 25,
          "maximum": 100
        },
        "offset": {
          "type": "integer",
          "description": "Offset for pagination (default: 0)",
          "default": 0
        },
        "format": {
          "type": "string",
          "enum": [
            "json",
            "xml"
          ],
          "default": "json",
          "description": "Response format. JSON is recommended for most use cases."
        },
        "gene_symbol": {
          "type": "string",
          "description": "Alias for query. Gene symbol to search (e.g., TP53, BRCA1)."
        },
        "gene": {
          "type": "string",
          "description": "Alias for query. Gene name to search."
        },
        "reviewed": {
          "type": "boolean",
          "description": "If true, return only reviewed Swiss-Prot entries. If false or omitted, includes unreviewed TrEMBL entries."
        },
        "organism": {
          "type": "string",
          "description": "Restrict the search to an organism by NAME (e.g. 'Klebsiella pneumoniae', 'Escherichia coli'); use 'taxid' for numeric NCBI taxonomy ids. A purely numeric value here (e.g. '9606') is auto-corrected to a taxid lookup and the response reports it in 'normalization_note'. By default a gene/protein-name search is biased toward human (taxid 9606); set this for non-human proteins. If omitted and no human match is found, the search auto-widens to all organisms."
        },
        "taxid": {
          "type": "string",
          "description": "Restrict the search to an NCBI taxonomy id (e.g. '573' for Klebsiella pneumoniae, '562' for E. coli). Alternative to 'organism' for precise scoping."
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "array",
      "description": "Array of protein objects matching the search query. Each protein contains accession, id, organism, protein names, gene names, features, and other annotations."
    }
    ```

**执行与来源：** `network`；上游接口名 `proteins_api_search`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("proteins_api_search", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.quick_go_annotations_by_gene` {#dnakit.databases.quick_go_annotations_by_gene}

Search Gene Ontology (GO) annotations for a specific gene product using the EBI QuickGO browser. Returns GO terms annotated to the protein/gene with evidence codes, qualifiers, and annotation sources. Accepts UniProtKB accessions (e.g., 'UniProtKB:P04637' for TP53). Provides more detailed annotation data than the standard GO API, including evidence provenance and assigned-by information. Useful for understanding gene function, localization, and biological processes.

```text
dnakit.databases.quick_go_annotations_by_gene(*, gene_product_id, aspect=None, taxon_id=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_product_id` | string | 必填（输入 Schema） | Gene product identifier, typically a UniProtKB accession. Format: &#x27;UniProtKB:P04637&#x27; (TP53), &#x27;UniProtKB:P38398&#x27; (BRCA1), &#x27;UniProtKB:P00533&#x27; (EGFR). |
| `aspect` | string &#124; null | None | Filter by GO aspect: &#x27;biological_process&#x27;, &#x27;molecular_function&#x27;, or &#x27;cellular_component&#x27;. Default: all aspects. |
| `taxon_id` | integer &#124; null | None | Filter by NCBI taxonomy ID. Example: 9606 (human), 10090 (mouse). |
| `limit` | integer | None | Maximum number of annotations to return. Default: 25, max: 100. {&quot;default&quot;: 25} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene_product_id": {
          "type": "string",
          "description": "Gene product identifier, typically a UniProtKB accession. Format: 'UniProtKB:P04637' (TP53), 'UniProtKB:P38398' (BRCA1), 'UniProtKB:P00533' (EGFR)."
        },
        "aspect": {
          "type": [
            "string",
            "null"
          ],
          "description": "Filter by GO aspect: 'biological_process', 'molecular_function', or 'cellular_component'. Default: all aspects."
        },
        "taxon_id": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Filter by NCBI taxonomy ID. Example: 9606 (human), 10090 (mouse)."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of annotations to return. Default: 25, max: 100.",
          "default": 25
        }
      },
      "required": [
        "gene_product_id"
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
              "go_id": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "GO term accession"
              },
              "go_name": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "GO term name"
              },
              "go_aspect": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "GO aspect (biological_process, etc.)"
              },
              "gene_product_id": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Gene product accession"
              },
              "symbol": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Gene symbol"
              },
              "qualifier": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Annotation qualifier"
              },
              "go_evidence": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Evidence code (IDA, IMP, ISS, etc.)"
              },
              "assigned_by": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Database that made the annotation"
              },
              "reference": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Literature reference"
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

**执行与来源：** `network`；上游接口名 `QuickGO_annotations_by_gene`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("quick_go_annotations_by_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.rcsbdata_get_entry` {#dnakit.databases.rcsbdata_get_entry}

Get comprehensive structure metadata for a PDB entry from the RCSB Data API. Returns experimental method, resolution, deposition/release dates, polymer/nonpolymer entity counts, molecular weight, assembly count, space group, and unit cell dimensions. This is the direct REST endpoint for PDB entry details, complementing GraphQL-based tools. Example: '4HHB' returns hemoglobin at 1.74A by X-ray diffraction, deposited 1984-03-07, 2 polymer entities.

```text
dnakit.databases.rcsbdata_get_entry(*, pdb_id=None, entry_id=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `pdb_id` | string | None | PDB entry ID (4 characters). Examples: &#x27;4HHB&#x27; (hemoglobin), &#x27;1TUP&#x27; (p53-DNA complex), &#x27;1M17&#x27; (EGFR kinase), &#x27;6LU7&#x27; (SARS-CoV-2 main protease). |
| `entry_id` | string | None | Alias for pdb_id (the RCSB-native term; the REST endpoint is /core/entry/{id}). Provide either pdb_id or entry_id. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "pdb_id": {
          "type": "string",
          "description": "PDB entry ID (4 characters). Examples: '4HHB' (hemoglobin), '1TUP' (p53-DNA complex), '1M17' (EGFR kinase), '6LU7' (SARS-CoV-2 main protease)."
        },
        "entry_id": {
          "type": "string",
          "description": "Alias for pdb_id (the RCSB-native term; the REST endpoint is /core/entry/{id}). Provide either pdb_id or entry_id."
        }
      },
      "anyOf": [
        {
          "required": [
            "pdb_id"
          ]
        },
        {
          "required": [
            "entry_id"
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
          "type": "object",
          "properties": {
            "pdb_id": {
              "type": [
                "string",
                "null"
              ],
              "description": "PDB entry ID"
            },
            "title": {
              "type": [
                "string",
                "null"
              ],
              "description": "Structure title"
            },
            "method": {
              "type": [
                "string",
                "null"
              ],
              "description": "Experimental method (X-RAY DIFFRACTION, ELECTRON MICROSCOPY, NMR, etc.)"
            },
            "resolution": {
              "type": [
                "number",
                "null"
              ],
              "description": "Resolution in Angstroms"
            },
            "deposit_date": {
              "type": [
                "string",
                "null"
              ],
              "description": "Deposition date"
            },
            "release_date": {
              "type": [
                "string",
                "null"
              ],
              "description": "Initial release date"
            },
            "revision_date": {
              "type": [
                "string",
                "null"
              ],
              "description": "Last revision date"
            },
            "polymer_entity_count": {
              "type": [
                "integer",
                "null"
              ],
              "description": "Number of polymer entities (proteins, nucleic acids)"
            },
            "nonpolymer_entity_count": {
              "type": [
                "integer",
                "null"
              ],
              "description": "Number of non-polymer entities (ligands, ions)"
            },
            "deposited_atom_count": {
              "type": [
                "integer",
                "null"
              ],
              "description": "Total deposited atom count"
            },
            "deposited_model_count": {
              "type": [
                "integer",
                "null"
              ],
              "description": "Deposited modeled polymer monomer count"
            },
            "molecular_weight": {
              "type": [
                "number",
                "null"
              ],
              "description": "Total molecular weight in kDa"
            },
            "assembly_count": {
              "type": [
                "integer",
                "null"
              ],
              "description": "Number of biological assemblies"
            },
            "space_group": {
              "type": [
                "string",
                "null"
              ],
              "description": "Crystal space group"
            },
            "unit_cell": {
              "type": [
                "object",
                "null"
              ],
              "description": "Unit cell dimensions (a, b, c, alpha, beta, gamma)"
            }
          },
          "required": [
            "pdb_id"
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

**执行与来源：** `network`；上游接口名 `RCSBData_get_entry`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("rcsbdata_get_entry", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.reactome_analysis_pathway_enrichment` {#dnakit.databases.reactome_analysis_pathway_enrichment}

Perform pathway overrepresentation (enrichment) analysis using the Reactome Analysis Service. Submit a list of gene symbols, UniProt IDs, or Ensembl IDs and identify which Reactome pathways are statistically enriched. Returns pathways ranked by p-value with FDR (False Discovery Rate) correction, entity counts, and reaction coverage. Supports projection to human pathways for cross-species analysis. INTERPRETING RESULTS: rank pathways by 'fdr' or 'p_value' - those are the enrichment statistics. 'entities_coverage' (= entities_found / entities_total) is the fraction of a pathway your identifiers cover. 'pathway_size_fraction_of_reactome' (deprecated alias 'entities_ratio') is Reactome's raw 'ratio': pathway size over Reactome's entity universe for the species, independent of the identifiers you submit - never rank by it. Example: submitting TP53, BRCA1, CDH1, EGFR, KRAS returns 394 enriched pathways including 'Developmental Lineage of Mammary Gland Myoepithelial Cells' (p=3.94e-08) and 'Signaling by ERBB2 KD Mutants'.

```text
dnakit.databases.reactome_analysis_pathway_enrichment(*, identifiers, page_size=None, include_disease=None, projection=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `identifiers` | string | 必填（输入 Schema） | Newline-separated list of gene/protein identifiers. Supports gene symbols (TP53, BRCA1), UniProt IDs (P04637, P38398), or Ensembl IDs. Example: &#x27;TP53\nBRCA1\nCDH1\nEGFR\nKRAS&#x27;. |
| `page_size` | integer &#124; null | None | Number of pathways to return (default 20, max 50). |
| `include_disease` | boolean &#124; null | None | Include disease pathways in results (default true). |
| `projection` | boolean &#124; null | None | Project identifiers to human Reactome pathways for cross-species analysis (default true). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "identifiers": {
          "type": "string",
          "description": "Newline-separated list of gene/protein identifiers. Supports gene symbols (TP53, BRCA1), UniProt IDs (P04637, P38398), or Ensembl IDs. Example: 'TP53\\nBRCA1\\nCDH1\\nEGFR\\nKRAS'."
        },
        "page_size": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Number of pathways to return (default 20, max 50)."
        },
        "include_disease": {
          "type": [
            "boolean",
            "null"
          ],
          "description": "Include disease pathways in results (default true)."
        },
        "projection": {
          "type": [
            "boolean",
            "null"
          ],
          "description": "Project identifiers to human Reactome pathways for cross-species analysis (default true)."
        }
      },
      "required": [
        "identifiers"
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
            "token": {
              "type": [
                "string",
                "null"
              ],
              "description": "Analysis token for retrieving results later"
            },
            "analysis_type": {
              "type": [
                "string",
                "null"
              ],
              "description": "OVERREPRESENTATION or EXPRESSION"
            },
            "projection": {
              "type": [
                "boolean",
                "null"
              ],
              "description": "Whether projection was applied"
            },
            "identifiers_not_found": {
              "type": "integer",
              "description": "Count of unmatched identifiers"
            },
            "pathways_found": {
              "type": "integer",
              "description": "Total enriched pathways"
            },
            "pathways": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "pathway_id": {
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
                  "species": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "is_disease": {
                    "type": [
                      "boolean",
                      "null"
                    ]
                  },
                  "is_lowest_level": {
                    "type": [
                      "boolean",
                      "null"
                    ]
                  },
                  "entities_found": {
                    "type": [
                      "integer",
                      "null"
                    ],
                    "description": "Number of submitted identifiers that mapped into this pathway."
                  },
                  "entities_total": {
                    "type": [
                      "integer",
                      "null"
                    ],
                    "description": "Total number of Reactome entities in this pathway, i.e. the pathway's size."
                  },
                  "entities_coverage": {
                    "type": [
                      "number",
                      "null"
                    ],
                    "description": "entities_found / entities_total: the fraction of this pathway covered by your submitted identifiers. Null when entities_total is 0 or unavailable."
                  },
                  "p_value": {
                    "type": [
                      "number",
                      "null"
                    ],
                    "description": "Overrepresentation p-value for this pathway. One of the two enrichment statistics - rank pathways by this or by fdr."
                  },
                  "fdr": {
                    "type": [
                      "number",
                      "null"
                    ],
                    "description": "Benjamini-Hochberg FDR-corrected p-value. Preferred statistic for ranking enriched pathways."
                  },
                  "reactions_found": {
                    "type": [
                      "integer",
                      "null"
                    ]
                  },
                  "reactions_total": {
                    "type": [
                      "integer",
                      "null"
                    ]
                  },
                  "pathway_size_fraction_of_reactome": {
                    "type": [
                      "number",
                      "null"
                    ],
                    "description": "Reactome's raw 'ratio' field: entities_total divided by the total entity count of Reactome for this pathway's species (measured at 16251 entities for Homo sapiens in Reactome release 97; the divisor differs per species - e.g. 8896 for Gallus gallus, 10752 for Bos taurus - and changes between releases, so never hard-code it). It measures pathway SIZE only and is completely independent of the submitted identifiers: the same pathway returns the same value for any gene list. It is NOT entities_found/entities_total - use entities_coverage for that - and it must not be used to rank pathways by enrichment."
                  },
                  "entities_ratio": {
                    "type": [
                      "number",
                      "null"
                    ],
                    "description": "DEPRECATED alias of pathway_size_fraction_of_reactome, kept for backward compatibility and carrying an identical value. Despite its name it is NOT entities_found/entities_total; use entities_coverage for pathway coverage and p_value/fdr for enrichment."
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

**执行与来源：** `network`；上游接口名 `ReactomeAnalysis_pathway_enrichment`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("reactome_analysis_pathway_enrichment", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.reactome_content_search` {#dnakit.databases.reactome_content_search}

Search Reactome for pathways, reactions, and biological entities by keyword. Performs free-text search across Reactome's knowledge base of biological pathways. Supports filtering by species and entity type. Example: searching 'apoptosis' in Homo sapiens returns pathways like R-HSA-109581 (Apoptosis), R-HSA-169911 (Regulation of Apoptosis), R-HSA-9635465 (Suppression of apoptosis). Use this to discover relevant pathways before querying their details. Results are paged: 'count' is how many entries came back, 'total_results' is how many matched in Reactome overall, and 'has_more' flags that the list is a slice -- use 'rows' and 'start' to page through the rest.

```text
dnakit.databases.reactome_content_search(*, query, species=None, types=None, rows=None, start=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填（输入 Schema） | Search query. Examples: &#x27;apoptosis&#x27;, &#x27;TP53&#x27;, &#x27;cell cycle&#x27;, &#x27;DNA repair&#x27;, &#x27;insulin signaling&#x27;. |
| `species` | string | None | Species name for filtering results. Default: &#x27;Homo sapiens&#x27;. Other options: &#x27;Mus musculus&#x27;, &#x27;Rattus norvegicus&#x27;, &#x27;Danio rerio&#x27;, &#x27;Drosophila melanogaster&#x27;. {&quot;default&quot;: &quot;Homo sapiens&quot;} |
| `types` | string | None | Entity type to search. Options: &#x27;Pathway&#x27;, &#x27;Reaction&#x27;, &#x27;Complex&#x27;, &#x27;Protein&#x27;, &#x27;SmallMolecule&#x27;. Default: &#x27;Pathway&#x27;. {&quot;default&quot;: &quot;Pathway&quot;} |
| `rows` | integer &#124; null | None | Number of results to return per page (1-1000). Omit to use Reactome&#x27;s default page of 10. Queries often match far more than that (&#x27;DNA repair&#x27; in Homo sapiens matches 455 pathways), so raise this to see beyond the first page. Compare the returned &#x27;count&#x27; against &#x27;total_results&#x27; to know how much you are holding. |
| `start` | integer &#124; null | None | Zero-based offset of the first result to return, for paging through large result sets. Use together with &#x27;rows&#x27; (e.g. rows=50 with start=0, then start=50). Omit to start at the first result. |
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
          "description": "Search query. Examples: 'apoptosis', 'TP53', 'cell cycle', 'DNA repair', 'insulin signaling'."
        },
        "species": {
          "type": "string",
          "description": "Species name for filtering results. Default: 'Homo sapiens'. Other options: 'Mus musculus', 'Rattus norvegicus', 'Danio rerio', 'Drosophila melanogaster'.",
          "default": "Homo sapiens"
        },
        "types": {
          "type": "string",
          "description": "Entity type to search. Options: 'Pathway', 'Reaction', 'Complex', 'Protein', 'SmallMolecule'. Default: 'Pathway'.",
          "default": "Pathway"
        },
        "rows": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Number of results to return per page (1-1000). Omit to use Reactome's default page of 10. Queries often match far more than that ('DNA repair' in Homo sapiens matches 455 pathways), so raise this to see beyond the first page. Compare the returned 'count' against 'total_results' to know how much you are holding."
        },
        "start": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Zero-based offset of the first result to return, for paging through large result sets. Use together with 'rows' (e.g. rows=50 with start=0, then start=50). Omit to start at the first result."
        }
      },
      "required": [
        "query"
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
            "query": {
              "type": "string",
              "description": "Search query used"
            },
            "species": {
              "type": "string",
              "description": "Species filter applied"
            },
            "count": {
              "type": "integer",
              "description": "Number of entries actually returned in 'results' on this page"
            },
            "total_results": {
              "type": "integer",
              "description": "Total number of results matching the query in Reactome (upstream 'numberOfMatches'), which is usually far larger than 'count'"
            },
            "start": {
              "type": "integer",
              "description": "Zero-based offset of the first returned entry within the full result set"
            },
            "has_more": {
              "type": "boolean",
              "description": "True when more matching entries exist beyond this page; raise 'rows' or advance 'start' to reach them"
            },
            "results": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "type": {
                    "type": "string",
                    "description": "Entity type (Pathway, Reaction, etc.)"
                  },
                  "stId": {
                    "type": [
                      "string",
                      "null"
                    ],
                    "description": "Reactome stable identifier"
                  },
                  "name": {
                    "type": [
                      "string",
                      "null"
                    ],
                    "description": "Entity name"
                  },
                  "species": {
                    "type": "array",
                    "description": "Species list"
                  },
                  "is_disease": {
                    "type": "boolean",
                    "description": "Whether this is a disease-related pathway"
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

**执行与来源：** `network`；上游接口名 `ReactomeContent_search`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("reactome_content_search", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.reactome_get_participants` {#dnakit.databases.reactome_get_participants}

Get all participants (physical entities) of a reaction or event. Returns list of physical entities (proteins, complexes, small molecules) participating in the event, including their reference entities.

```text
dnakit.databases.reactome_get_participants(*, stId, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `stId` | string | 必填（输入 Schema） | Event Stable ID (reaction or pathway, e.g., &#x27;R-HSA-111289&#x27;). To find reaction IDs, use Reactome_get_pathway_reactions with a pathway ID. To find pathway IDs, use Reactome_list_top_pathways. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "stId": {
          "type": "string",
          "description": "Event Stable ID (reaction or pathway, e.g., 'R-HSA-111289'). To find reaction IDs, use Reactome_get_pathway_reactions with a pathway ID. To find pathway IDs, use Reactome_list_top_pathways."
        }
      },
      "required": [
        "stId"
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
          "peDbId": {
            "type": "integer",
            "description": "Physical entity database ID"
          },
          "displayName": {
            "type": "string",
            "description": "Participant display name"
          },
          "schemaClass": {
            "type": "string",
            "description": "Schema class (e.g., 'Complex', 'SimpleEntity')"
          },
          "refEntities": {
            "type": "array",
            "items": {
              "type": "object"
            },
            "description": "Reference entities (UniProt, ChEBI, etc.) associated with this participant"
          }
        },
        "description": "Physical entity participating in the event"
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `Reactome_get_participants`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("reactome_get_participants", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.reactome_get_pathway` {#dnakit.databases.reactome_get_pathway}

Get detailed information about a specific pathway by its Stable ID. Returns comprehensive pathway data including metadata, compartments, literature references, orthologous events, and related pathways. Uses the unified query endpoint that works for all entity types.

```text
dnakit.databases.reactome_get_pathway(*, stId, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `stId` | string | 必填（输入 Schema） | Pathway Stable ID (e.g., &#x27;R-HSA-73817&#x27;). To find pathway IDs, use Reactome_list_top_pathways to get top-level pathways, or Reactome_map_uniprot_to_pathways if you have a UniProt protein ID. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "stId": {
          "type": "string",
          "description": "Pathway Stable ID (e.g., 'R-HSA-73817'). To find pathway IDs, use Reactome_list_top_pathways to get top-level pathways, or Reactome_map_uniprot_to_pathways if you have a UniProt protein ID."
        }
      },
      "required": [
        "stId"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "dbId": {
          "type": "integer",
          "description": "Database ID"
        },
        "displayName": {
          "type": "string",
          "description": "Pathway display name"
        },
        "stId": {
          "type": "string",
          "description": "Stable identifier (e.g., 'R-HSA-73817')"
        },
        "stIdVersion": {
          "type": "string",
          "description": "Versioned stable ID"
        },
        "isInDisease": {
          "type": "boolean",
          "description": "Whether pathway is associated with disease"
        },
        "isInferred": {
          "type": "boolean",
          "description": "Whether pathway is inferred"
        },
        "maxDepth": {
          "type": "integer",
          "description": "Maximum depth of pathway hierarchy"
        },
        "name": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "Pathway names"
        },
        "releaseDate": {
          "type": "string",
          "description": "Release date"
        },
        "speciesName": {
          "type": "string",
          "description": "Species name"
        },
        "compartment": {
          "type": "array",
          "items": {
            "type": "object"
          },
          "description": "Cellular compartments"
        },
        "figure": {
          "type": "array",
          "items": {
            "type": "object"
          },
          "description": "Pathway figures"
        },
        "goBiologicalProcess": {
          "type": "object",
          "description": "GO biological process annotation"
        },
        "literatureReference": {
          "type": "array",
          "items": {
            "type": "object"
          },
          "description": "Literature references"
        },
        "orthologousEvent": {
          "type": "array",
          "items": {
            "type": "object"
          },
          "description": "Orthologous events in other species"
        },
        "hasEvent": {
          "type": "array",
          "items": {
            "type": "object"
          },
          "description": "Events (reactions/subpathways) in this pathway"
        },
        "schemaClass": {
          "type": "string",
          "description": "Schema class (e.g., 'Pathway')"
        },
        "className": {
          "type": "string",
          "description": "Class name"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `Reactome_get_pathway`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("reactome_get_pathway", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.reactome_get_pathway_reactions` {#dnakit.databases.reactome_get_pathway_reactions}

Query all reactions and subpathways contained under a Pathway using Pathway Stable ID. Returns list of events (reactions and subpathways) contained in the pathway, including their metadata and relationships.

```text
dnakit.databases.reactome_get_pathway_reactions(*, stId, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `stId` | string | 必填（输入 Schema） | Pathway Stable ID (e.g., &#x27;R-HSA-73817&#x27;). To find pathway IDs, use Reactome_list_top_pathways to get top-level pathways, or Reactome_map_uniprot_to_pathways if you have a UniProt protein ID. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "stId": {
          "type": "string",
          "description": "Pathway Stable ID (e.g., 'R-HSA-73817'). To find pathway IDs, use Reactome_list_top_pathways to get top-level pathways, or Reactome_map_uniprot_to_pathways if you have a UniProt protein ID."
        }
      },
      "required": [
        "stId"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "array",
      "items": {
        "oneOf": [
          {
            "type": "integer",
            "description": "Event ID reference"
          },
          {
            "type": "object",
            "properties": {
              "dbId": {
                "type": "integer",
                "description": "Database ID"
              },
              "displayName": {
                "type": "string",
                "description": "Event display name"
              },
              "stId": {
                "type": "string",
                "description": "Stable identifier"
              },
              "stIdVersion": {
                "type": "string",
                "description": "Versioned stable ID"
              },
              "isInDisease": {
                "type": "boolean",
                "description": "Whether event is associated with disease"
              },
              "isInferred": {
                "type": "boolean",
                "description": "Whether event is inferred"
              },
              "maxDepth": {
                "type": "integer",
                "description": "Maximum depth"
              },
              "name": {
                "type": "array",
                "items": {
                  "type": "string"
                },
                "description": "Event names"
              },
              "releaseDate": {
                "type": "string",
                "description": "Release date"
              },
              "speciesName": {
                "type": "string",
                "description": "Species name"
              },
              "eventOf": {
                "type": "array",
                "items": {
                  "type": [
                    "integer",
                    "object"
                  ]
                },
                "description": "Parent events (can be IDs or objects)"
              },
              "category": {
                "type": "string",
                "description": "Event category"
              },
              "schemaClass": {
                "type": "string",
                "description": "Schema class (e.g., 'Reaction', 'Pathway')"
              },
              "className": {
                "type": "string",
                "description": "Class name"
              }
            },
            "description": "Event object with full details"
          }
        ]
      },
      "description": "List of events (reactions or subpathways), can be IDs or full objects"
    }
    ```

**执行与来源：** `network`；上游接口名 `Reactome_get_pathway_reactions`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("reactome_get_pathway_reactions", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.reactome_map_uniprot_to_pathways` {#dnakit.databases.reactome_map_uniprot_to_pathways}

Map a UniProt protein identifier to Reactome pathways. Returns all pathways that contain this protein.

```text
dnakit.databases.reactome_map_uniprot_to_pathways(*, uniprot_id=None, id=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `uniprot_id` | string | None | UniProt protein accession (e.g., &#x27;P04637&#x27; for TP53, &#x27;P00533&#x27; for EGFR) |
| `id` | string | None | DNAKit optional alias for uniprot_id |
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
          "description": "UniProt protein accession (e.g., 'P04637' for TP53, 'P00533' for EGFR)"
        },
        "id": {
          "type": "string",
          "description": "DNAKit optional alias for uniprot_id"
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
          "dbId": {
            "type": "integer",
            "description": "Database ID"
          },
          "displayName": {
            "type": "string",
            "description": "Pathway display name"
          },
          "stId": {
            "type": "string",
            "description": "Stable identifier"
          },
          "stIdVersion": {
            "type": "string",
            "description": "Versioned stable ID"
          },
          "isInDisease": {
            "type": "boolean",
            "description": "Whether pathway is associated with disease"
          },
          "isInferred": {
            "type": "boolean",
            "description": "Whether pathway is inferred"
          },
          "maxDepth": {
            "type": "integer",
            "description": "Maximum depth"
          },
          "name": {
            "type": "array",
            "items": {
              "type": "string"
            },
            "description": "Pathway names"
          },
          "releaseDate": {
            "type": "string",
            "description": "Release date"
          },
          "speciesName": {
            "type": "string",
            "description": "Species name"
          },
          "hasDiagram": {
            "type": "boolean",
            "description": "Whether pathway has diagram"
          },
          "hasEHLD": {
            "type": "boolean",
            "description": "Whether pathway has EHLD diagram"
          },
          "lastUpdatedDate": {
            "type": "string",
            "description": "Last update date"
          },
          "schemaClass": {
            "type": "string",
            "description": "Schema class (e.g., 'Pathway')"
          },
          "className": {
            "type": "string",
            "description": "Class name"
          }
        },
        "description": "Pathway objects containing the UniProt protein"
      }
    }
    ```

??? info "DNAKit 参数扩展规则"

    ```json
    {
      "aliases": {
        "id": "uniprot_id"
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `Reactome_map_uniprot_to_pathways`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("reactome_map_uniprot_to_pathways", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.string_functional_enrichment` {#dnakit.databases.string_functional_enrichment}

Identify enriched biological functions, pathways, and processes for a protein set using STRING (Search Tool for Retrieval of Interacting Genes/Proteins) database. Performs statistical enrichment analysis against GO (Gene Ontology) terms, KEGG pathways, Reactome pathways, and other annotation databases. Returns significantly enriched terms with FDR-corrected p-values. Minimum 3-4 proteins required for meaningful statistical analysis. No API key required. Use for: discovering pathways enriched in protein network, identifying biological processes, finding shared functions among proteins, interpreting omics data (proteomics, interactomics), hypothesis generation.

```text
dnakit.databases.string_functional_enrichment(*, protein_ids, species=None, category=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `protein_ids` | array | 必填（输入 Schema） | List of protein identifiers (UniProt IDs, gene names, Ensembl IDs). Minimum 3 proteins recommended for meaningful enrichment. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}, &quot;minItems&quot;: 2} |
| `species` | integer | None | NCBI taxonomy ID (default: 9606 for human) {&quot;default&quot;: 9606} |
| `category` | string | None | Enrichment category: &#x27;Process&#x27; (GO Biological Process), &#x27;Component&#x27; (GO Cellular Component), &#x27;Function&#x27; (GO Molecular Function), &#x27;KEGG&#x27;, &#x27;Reactome&#x27;, &#x27;WikiPathways&#x27;, &#x27;COMPARTMENTS&#x27;, &#x27;TISSUES&#x27;, &#x27;DISEASES&#x27; {&quot;enum&quot;: &#91;&quot;Process&quot;, &quot;Component&quot;, &quot;Function&quot;, &quot;KEGG&quot;, &quot;Reactome&quot;, &quot;WikiPathways&quot;, &quot;COMPARTMENTS&quot;, &quot;TISSUES&quot;, &quot;DISEASES&quot;&#93;, &quot;default&quot;: &quot;Process&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "protein_ids": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "List of protein identifiers (UniProt IDs, gene names, Ensembl IDs). Minimum 3 proteins recommended for meaningful enrichment.",
          "minItems": 2
        },
        "species": {
          "type": "integer",
          "description": "NCBI taxonomy ID (default: 9606 for human)",
          "default": 9606
        },
        "category": {
          "type": "string",
          "description": "Enrichment category: 'Process' (GO Biological Process), 'Component' (GO Cellular Component), 'Function' (GO Molecular Function), 'KEGG', 'Reactome', 'WikiPathways', 'COMPARTMENTS', 'TISSUES', 'DISEASES'",
          "enum": [
            "Process",
            "Component",
            "Function",
            "KEGG",
            "Reactome",
            "WikiPathways",
            "COMPARTMENTS",
            "TISSUES",
            "DISEASES"
          ],
          "default": "Process"
        }
      },
      "required": [
        "protein_ids"
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
          "category": {
            "type": "string"
          },
          "term": {
            "type": "string"
          },
          "number_of_genes": {
            "type": "integer"
          },
          "number_of_genes_in_background": {
            "type": "integer"
          },
          "ncbiTaxonId": {
            "type": "integer"
          },
          "inputGenes": {
            "type": [
              "string",
              "array",
              "null"
            ]
          },
          "preferredNames": {
            "type": [
              "string",
              "array",
              "null"
            ]
          },
          "p_value": {
            "type": "number"
          },
          "fdr": {
            "type": "number"
          },
          "description": {
            "type": "string"
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `STRING_functional_enrichment`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("string_functional_enrichment", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.string_get_functional_annotations` {#dnakit.databases.string_get_functional_annotations}

Get comprehensive functional annotations for a protein from STRING database. Returns per-protein annotations from Gene Ontology (Process, Function, Component), KEGG pathways, Reactome, WikiPathways, InterPro domains, Pfam, HPO phenotypes, DISEASES, TISSUES expression, and COMPARTMENTS subcellular localization. Different from STRING_functional_enrichment which tests for statistical enrichment of a protein SET - this tool returns all known annotations for individual proteins. Example: 'TP53' returns 1,877 annotations across 13 categories including 655 biological processes, 462 HPO phenotypes, 212 disease associations.

```text
dnakit.databases.string_get_functional_annotations(*, identifiers, species=None, category=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `identifiers` | string | 必填（输入 Schema） | Protein identifier (gene name or STRING ID). Examples: &#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;EGFR&#x27;, &#x27;9606.ENSP00000269305&#x27;. |
| `species` | integer | None | NCBI taxonomy ID (default 9606 for human). Examples: 9606 (human), 10090 (mouse), 7227 (fly), 6239 (worm). {&quot;default&quot;: 9606} |
| `category` | string | None | Filter by annotation category. Options: &#x27;Process&#x27; (GO BP), &#x27;Function&#x27; (GO MF), &#x27;Component&#x27; (GO CC), &#x27;KEGG&#x27;, &#x27;RCTM&#x27; (Reactome), &#x27;WikiPathways&#x27;, &#x27;InterPro&#x27;, &#x27;Pfam&#x27;, &#x27;HPO&#x27;, &#x27;DISEASES&#x27;, &#x27;TISSUES&#x27;, &#x27;COMPARTMENTS&#x27;, &#x27;Keyword&#x27;. Leave null for all categories. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "identifiers": {
          "type": "string",
          "description": "Protein identifier (gene name or STRING ID). Examples: 'TP53', 'BRCA1', 'EGFR', '9606.ENSP00000269305'."
        },
        "species": {
          "type": "integer",
          "description": "NCBI taxonomy ID (default 9606 for human). Examples: 9606 (human), 10090 (mouse), 7227 (fly), 6239 (worm).",
          "default": 9606
        },
        "category": {
          "type": "string",
          "description": "Filter by annotation category. Options: 'Process' (GO BP), 'Function' (GO MF), 'Component' (GO CC), 'KEGG', 'RCTM' (Reactome), 'WikiPathways', 'InterPro', 'Pfam', 'HPO', 'DISEASES', 'TISSUES', 'COMPARTMENTS', 'Keyword'. Leave null for all categories."
        }
      },
      "required": [
        "identifiers"
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
            "identifiers": {
              "type": "string",
              "description": "Queried protein identifier"
            },
            "species": {
              "type": "integer",
              "description": "NCBI taxonomy ID"
            },
            "total_annotations": {
              "type": "integer",
              "description": "Total number of annotations across all categories"
            },
            "category_summary": {
              "type": "object",
              "description": "Number of annotations per category"
            },
            "annotations": {
              "type": "object",
              "description": "Annotations grouped by category, each with term, description, number_of_genes, input_genes (max 15 per category)"
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

**执行与来源：** `network`；上游接口名 `STRING_get_functional_annotations`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("string_get_functional_annotations", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.string_get_functional_enrichment` {#dnakit.databases.string_get_functional_enrichment}

Run STRING enrichment for carriage-return or newline-separated protein identifiers.

```text
dnakit.databases.string_get_functional_enrichment(*, identifiers, species=9606, category=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `identifiers` | string | 必填（输入 Schema） | 以换行或回车分隔的蛋白标识；将转换成 STRING 标识列表。 |
| `species` | integer | 9606 | 物种名或 Taxonomy ID；类型与支持值取决于当前接口。 {&quot;default&quot;: 9606} |
| `category` | string | None | 分类或结果类别，按当前接口的枚举/来源定义填写。 |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回 STRING_functional_enrichment 的富集结果，见 [string_functional_enrichment](pathways.md#dnakit.databases.string_functional_enrichment)。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "identifiers": {
          "type": "string"
        },
        "species": {
          "type": "integer",
          "default": 9606
        },
        "category": {
          "type": "string"
        }
      },
      "required": [
        "identifiers"
      ]
    }
    ```

**执行与来源：** `network`；上游接口名 `STRING_get_functional_enrichment`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("string_get_functional_enrichment", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.string_get_interaction_partners` {#dnakit.databases.string_get_interaction_partners}

Get interaction partners for a single protein from STRING database. Returns a ranked list of proteins that interact with the query protein, along with detailed confidence scores by evidence type (experimental, database, text-mining, etc.). Different from STRING_get_network which returns all edges in a network - this returns a simple partner list focused on one protein. Useful for finding what proteins a given protein interacts with.

```text
dnakit.databases.string_get_interaction_partners(*, identifiers=None, species=None, limit=None, required_score=None, protein_ids=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `identifiers` | string | None | Protein identifier (gene name or STRING ID). Examples: &#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;EGFR&#x27;, &#x27;9606.ENSP00000269305&#x27; |
| `species` | integer | None | NCBI taxonomy ID (9606=human, 10090=mouse, 7227=fly, 6239=worm) {&quot;default&quot;: 9606} |
| `limit` | integer | None | Maximum number of interaction partners to return {&quot;default&quot;: 10} |
| `required_score` | integer &#124; null | None | Minimum combined STRING score (0-1000). 400=medium, 700=high, 900=highest confidence. |
| `protein_ids` | array &#124; string | None | 蛋白标识；输入列表或字符串类型见当前接口 Schema。 {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "identifiers": {
          "type": "string",
          "description": "Protein identifier (gene name or STRING ID). Examples: 'TP53', 'BRCA1', 'EGFR', '9606.ENSP00000269305'"
        },
        "species": {
          "type": "integer",
          "default": 9606,
          "description": "NCBI taxonomy ID (9606=human, 10090=mouse, 7227=fly, 6239=worm)"
        },
        "limit": {
          "type": "integer",
          "default": 10,
          "description": "Maximum number of interaction partners to return"
        },
        "required_score": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Minimum combined STRING score (0-1000). 400=medium, 700=high, 900=highest confidence."
        },
        "protein_ids": {
          "type": [
            "array",
            "string"
          ],
          "items": {
            "type": "string"
          }
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
          "description": "Array of interaction partner objects",
          "items": {
            "type": "object",
            "properties": {
              "stringId_A": {
                "type": "string"
              },
              "stringId_B": {
                "type": "string"
              },
              "preferredName_A": {
                "type": "string",
                "description": "Gene symbol of query protein"
              },
              "preferredName_B": {
                "type": "string",
                "description": "Gene symbol of interaction partner"
              },
              "ncbiTaxonId": {
                "type": [
                  "string",
                  "integer",
                  "null"
                ]
              },
              "score": {
                "type": "number",
                "description": "Combined confidence score (0-1)"
              },
              "nscore": {
                "type": "number"
              },
              "fscore": {
                "type": "number"
              },
              "pscore": {
                "type": "number"
              },
              "ascore": {
                "type": "number"
              },
              "escore": {
                "type": "number",
                "description": "Experimental evidence score"
              },
              "dscore": {
                "type": "number",
                "description": "Database evidence score"
              },
              "tscore": {
                "type": "number",
                "description": "Text-mining score"
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
      "aliases": {
        "protein_ids": "identifiers"
      },
      "parameters": {
        "protein_ids": {
          "type": [
            "array",
            "string"
          ],
          "items": {
            "type": "string"
          }
        }
      },
      "conversions": {
        "identifiers": "join_cr"
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `STRING_get_interaction_partners`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("string_get_interaction_partners", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.string_get_network` {#dnakit.databases.string_get_network}

Get the protein-protein interaction network for one or more proteins from STRING database. Returns all interactions (edges) between the query proteins and their neighbors, including confidence scores broken down by evidence type (experiments, co-expression, databases, text-mining, etc.). Different from STRING_get_functional_annotations (which returns GO/pathway annotations) and STRING_get_enrichment (which tests enrichment of a set). Use this to build or visualize PPI networks. Multiple proteins can be separated by carriage return characters (\r).

```text
dnakit.databases.string_get_network(*, identifiers=None, species=None, limit=None, required_score=None, protein_ids=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `identifiers` | string | None | Protein identifier(s). For multiple proteins, separate with &#x27;\r&#x27; (carriage return character). Examples: &#x27;TP53&#x27;, &#x27;BRCA1\rBRCA2\rTP53&#x27;, &#x27;9606.ENSP00000269305&#x27; |
| `species` | integer | None | NCBI taxonomy ID. Examples: 9606 (human), 10090 (mouse), 7227 (Drosophila), 6239 (C. elegans) {&quot;default&quot;: 9606} |
| `limit` | integer | None | Maximum number of interaction partners to include (per protein). 0 = only interactions among input proteins. {&quot;default&quot;: 10} |
| `required_score` | integer &#124; null | None | Minimum combined STRING score (0-1000). 400=medium, 700=high, 900=highest confidence. Default: 400. |
| `protein_ids` | array &#124; string | None | 蛋白标识；输入列表或字符串类型见当前接口 Schema。 {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "identifiers": {
          "type": "string",
          "description": "Protein identifier(s). For multiple proteins, separate with '\\r' (carriage return character). Examples: 'TP53', 'BRCA1\\rBRCA2\\rTP53', '9606.ENSP00000269305'"
        },
        "species": {
          "type": "integer",
          "default": 9606,
          "description": "NCBI taxonomy ID. Examples: 9606 (human), 10090 (mouse), 7227 (Drosophila), 6239 (C. elegans)"
        },
        "limit": {
          "type": "integer",
          "default": 10,
          "description": "Maximum number of interaction partners to include (per protein). 0 = only interactions among input proteins."
        },
        "required_score": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Minimum combined STRING score (0-1000). 400=medium, 700=high, 900=highest confidence. Default: 400."
        },
        "protein_ids": {
          "type": [
            "array",
            "string"
          ],
          "items": {
            "type": "string"
          }
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
          "description": "Array of interaction objects (edges in the network)",
          "items": {
            "type": "object",
            "properties": {
              "stringId_A": {
                "type": "string",
                "description": "STRING ID of protein A (e.g., '9606.ENSP00000269305')"
              },
              "stringId_B": {
                "type": "string",
                "description": "STRING ID of protein B"
              },
              "preferredName_A": {
                "type": "string",
                "description": "Gene symbol of protein A (e.g., 'TP53')"
              },
              "preferredName_B": {
                "type": "string",
                "description": "Gene symbol of protein B"
              },
              "ncbiTaxonId": {
                "type": [
                  "string",
                  "integer",
                  "null"
                ],
                "description": "NCBI taxonomy ID"
              },
              "score": {
                "type": "number",
                "description": "Combined confidence score (0-1)"
              },
              "nscore": {
                "type": "number",
                "description": "Neighborhood score"
              },
              "fscore": {
                "type": "number",
                "description": "Gene fusion score"
              },
              "pscore": {
                "type": "number",
                "description": "Phylogenetic profile score"
              },
              "ascore": {
                "type": "number",
                "description": "Co-expression score"
              },
              "escore": {
                "type": "number",
                "description": "Experimental evidence score"
              },
              "dscore": {
                "type": "number",
                "description": "Database evidence score"
              },
              "tscore": {
                "type": "number",
                "description": "Text-mining score"
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
      "aliases": {
        "protein_ids": "identifiers"
      },
      "parameters": {
        "protein_ids": {
          "type": [
            "array",
            "string"
          ],
          "items": {
            "type": "string"
          }
        }
      },
      "conversions": {
        "identifiers": "join_cr"
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `STRING_get_network`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("string_get_network", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.string_ppi_enrichment` {#dnakit.databases.string_ppi_enrichment}

Test if your protein set has more interactions than expected by chance (PPI = Protein-Protein Interaction enrichment). Compares observed interactions in your network to random expectation - significant enrichment suggests proteins work together functionally. Returns enrichment p-value and expected vs observed interaction counts. Useful for validating that your protein list forms a real functional module rather than random proteins. No API key required. Use for: validating protein complex predictions, testing if proteins form functional module, quality control for network analysis, distinguishing real complexes from random protein lists, assessing functional coherence.

```text
dnakit.databases.string_ppi_enrichment(*, protein_ids, species=None, confidence_score=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `protein_ids` | array | 必填（输入 Schema） | List of protein identifiers (UniProt IDs, gene names, Ensembl IDs). Minimum 3 proteins required. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}, &quot;minItems&quot;: 3} |
| `species` | integer | None | NCBI taxonomy ID (default: 9606 for human) {&quot;default&quot;: 9606} |
| `confidence_score` | number | None | Minimum confidence score for counting interactions (0-1, default: 0.4) {&quot;minimum&quot;: 0, &quot;maximum&quot;: 1, &quot;default&quot;: 0.4} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "protein_ids": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "List of protein identifiers (UniProt IDs, gene names, Ensembl IDs). Minimum 3 proteins required.",
          "minItems": 3
        },
        "species": {
          "type": "integer",
          "description": "NCBI taxonomy ID (default: 9606 for human)",
          "default": 9606
        },
        "confidence_score": {
          "type": "number",
          "description": "Minimum confidence score for counting interactions (0-1, default: 0.4)",
          "minimum": 0,
          "maximum": 1,
          "default": 0.4
        }
      },
      "required": [
        "protein_ids"
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
            "number_of_nodes": {
              "type": "integer"
            },
            "number_of_edges": {
              "type": "integer"
            },
            "expected_number_of_edges": {
              "type": "number"
            },
            "p_value": {
              "type": "number"
            },
            "avg_node_degree": {
              "type": "number"
            },
            "local_clustering_coefficient": {
              "type": "number"
            }
          }
        },
        {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "number_of_nodes": {
                "type": "integer"
              },
              "number_of_edges": {
                "type": "integer"
              },
              "expected_number_of_edges": {
                "type": "number"
              },
              "p_value": {
                "type": "number"
              },
              "avg_node_degree": {
                "type": "number"
              },
              "local_clustering_coefficient": {
                "type": "number"
              }
            }
          }
        }
      ]
    }
    ```

**执行与来源：** `network`；上游接口名 `STRING_ppi_enrichment`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("string_ppi_enrichment", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.structure_annotate_per_residue` {#dnakit.databases.structure_annotate_per_residue}

Per-residue structural annotation from a PDB structure: which residues sit at a binding interface (vs a partner chain), which line a ligand pocket, which are buried (core) vs solvent-exposed (surface). For each residue of the target chain returns: 1-letter amino acid, min side-chain-heavy-atom distance to the partner chain(s), min distance to ligand heavy atoms, relative SASA (computed with freesasa on the ISOLATED target chain — partner chains would bias the call by burying interface residues), a region label in {interface, ligand, both, other}, and an is_core boolean. Optionally fetches secondary structure (helix/strand/coil) from PDBe REST. This is the structural annotation track plotted under a DMS heatmap and the structural prior SAE feature drops are read against. Methodology adapted from an upstream research workflow.

```text
dnakit.databases.structure_annotate_per_residue(*, operation=None, pdb_id=None, pdb_content=None, target_chain=None, partner_chains=None, ligand_resnames=None, distance_cutoff=None, core_rsa_cutoff=None, include_secondary_structure=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;annotate_per_residue&quot;&#93;, &quot;default&quot;: &quot;annotate_per_residue&quot;} |
| `pdb_id` | string &#124; null | None | 4-character PDB ID (e.g. &#x27;6VJJ&#x27;). If provided, the structure is downloaded from RCSB. Exactly one of pdb_id or pdb_content must be supplied. |
| `pdb_content` | string &#124; null | None | Raw PDB file content (string). Use when the structure is local (e.g. AlphaFold-predicted, in-memory). |
| `target_chain` | string | None | Chain ID of the protein to annotate (e.g. &#x27;A&#x27;). {&quot;default&quot;: &quot;A&quot;} |
| `partner_chains` | array | None | Chain IDs that define the binding partner(s). Used to compute the interface. Empty list = no interface analysis (all dist_partner values null). {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}, &quot;default&quot;: &#91;&#93;} |
| `ligand_resnames` | array | None | Three-letter residue names of bound ligands to define the ligand pocket (e.g. &#91;&#x27;GNP&#x27;, &#x27;MG&#x27;&#93; for GTP-analogue + magnesium). Searched across all chains. Empty list = no ligand-pocket analysis. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}, &quot;default&quot;: &#91;&#93;} |
| `distance_cutoff` | number | None | Distance cutoff in angstroms for interface and ligand-pocket classification. {&quot;default&quot;: 5.0} |
| `core_rsa_cutoff` | number | None | Relative solvent accessibility threshold below which a residue is classified as &#x27;core&#x27;. {&quot;default&quot;: 0.25} |
| `include_secondary_structure` | boolean | None | If true and pdb_id is provided, fetch per-residue secondary structure (helix/strand/coil) from PDBe REST. Only meaningful for deposited PDB entries; ignored when pdb_content is supplied. {&quot;default&quot;: false} |
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
          "description": "Operation type",
          "enum": [
            "annotate_per_residue"
          ],
          "default": "annotate_per_residue"
        },
        "pdb_id": {
          "type": [
            "string",
            "null"
          ],
          "description": "4-character PDB ID (e.g. '6VJJ'). If provided, the structure is downloaded from RCSB. Exactly one of pdb_id or pdb_content must be supplied."
        },
        "pdb_content": {
          "type": [
            "string",
            "null"
          ],
          "description": "Raw PDB file content (string). Use when the structure is local (e.g. AlphaFold-predicted, in-memory)."
        },
        "target_chain": {
          "type": "string",
          "description": "Chain ID of the protein to annotate (e.g. 'A').",
          "default": "A"
        },
        "partner_chains": {
          "type": "array",
          "description": "Chain IDs that define the binding partner(s). Used to compute the interface. Empty list = no interface analysis (all dist_partner values null).",
          "items": {
            "type": "string"
          },
          "default": []
        },
        "ligand_resnames": {
          "type": "array",
          "description": "Three-letter residue names of bound ligands to define the ligand pocket (e.g. ['GNP', 'MG'] for GTP-analogue + magnesium). Searched across all chains. Empty list = no ligand-pocket analysis.",
          "items": {
            "type": "string"
          },
          "default": []
        },
        "distance_cutoff": {
          "type": "number",
          "description": "Distance cutoff in angstroms for interface and ligand-pocket classification.",
          "default": 5.0
        },
        "core_rsa_cutoff": {
          "type": "number",
          "description": "Relative solvent accessibility threshold below which a residue is classified as 'core'.",
          "default": 0.25
        },
        "include_secondary_structure": {
          "type": "boolean",
          "description": "If true and pdb_id is provided, fetch per-residue secondary structure (helix/strand/coil) from PDBe REST. Only meaningful for deposited PDB entries; ignored when pdb_content is supplied.",
          "default": false
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
            "status": {
              "type": "string",
              "enum": [
                "success"
              ]
            },
            "pdb_id": {
              "type": [
                "string",
                "null"
              ]
            },
            "target_chain": {
              "type": "string"
            },
            "partner_chains": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "ligand_resnames": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "distance_cutoff": {
              "type": "number"
            },
            "core_rsa_cutoff": {
              "type": "number"
            },
            "n_residues": {
              "type": "integer"
            },
            "annotations": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "position": {
                    "type": "integer"
                  },
                  "aa": {
                    "type": "string"
                  },
                  "dist_partner": {
                    "type": [
                      "number",
                      "null"
                    ]
                  },
                  "dist_ligand": {
                    "type": [
                      "number",
                      "null"
                    ]
                  },
                  "rsa": {
                    "type": [
                      "number",
                      "null"
                    ]
                  },
                  "region": {
                    "type": "string",
                    "enum": [
                      "interface",
                      "ligand",
                      "both",
                      "other"
                    ]
                  },
                  "is_core": {
                    "type": "boolean"
                  },
                  "ss_element": {
                    "type": "string"
                  }
                }
              }
            },
            "method": {
              "type": "object"
            },
            "provenance": {
              "type": "string"
            }
          },
          "required": [
            "status",
            "target_chain",
            "n_residues",
            "annotations"
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

**执行与来源：** `network`；上游接口名 `Structure_annotate_per_residue`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("structure_annotate_per_residue", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.uni_prot_get_disease_variants_by_accession` {#dnakit.databases.uni_prot_get_disease_variants_by_accession}

Extract all variants (feature type = VARIANT) and their related annotations from UniProtKB entry.

```text
dnakit.databases.uni_prot_get_disease_variants_by_accession(*, accession, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string | 必填（输入 Schema） | UniProtKB accession, e.g., P05067. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "accession": {
          "type": "string",
          "description": "UniProtKB accession, e.g., P05067."
        }
      },
      "required": [
        "accession"
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
          "type": {
            "type": "string"
          },
          "location": {
            "type": "object",
            "properties": {
              "start": {
                "type": "object",
                "properties": {
                  "value": {
                    "type": [
                      "integer",
                      "number"
                    ]
                  },
                  "modifier": {
                    "type": "string"
                  }
                }
              },
              "end": {
                "type": "object",
                "properties": {
                  "value": {
                    "type": [
                      "integer",
                      "number"
                    ]
                  },
                  "modifier": {
                    "type": "string"
                  }
                }
              }
            }
          },
          "description": {
            "type": "string"
          },
          "featureId": {
            "type": "string"
          },
          "alternativeSequence": {
            "type": "object",
            "properties": {
              "originalSequence": {
                "type": "string"
              },
              "alternativeSequences": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              }
            }
          },
          "featureCrossReferences": {
            "type": "array",
            "items": {
              "type": "object",
              "properties": {
                "database": {
                  "type": "string"
                },
                "id": {
                  "type": "string"
                }
              }
            }
          },
          "evidences": {
            "type": "array",
            "items": {
              "type": "object",
              "properties": {
                "evidenceCode": {
                  "type": "string"
                },
                "source": {
                  "type": "string"
                },
                "id": {
                  "type": "string"
                }
              }
            }
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `UniProt_get_disease_variants_by_accession`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("uni_prot_get_disease_variants_by_accession", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.uni_prot_get_entry_by_accession` {#dnakit.databases.uni_prot_get_entry_by_accession}

Get the complete JSON entry for a specified UniProtKB accession. WARNING: This tool returns the complete UniProtKB entry, which can be extremely large (40,000+ lines for some entries) and may cause LLM context limits to be exceeded or crashes. For most use cases, consider using more specific extraction tools instead: UniProt_get_function_by_accession, UniProt_get_sequence_by_accession, UniProt_get_recommended_name_by_accession, UniProt_get_organism_by_accession, UniProt_get_subcellular_location_by_accession, UniProt_get_disease_variants_by_accession, or UniProt_get_ptm_processing_by_accession.

```text
dnakit.databases.uni_prot_get_entry_by_accession(*, accession, compact=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string | 必填（输入 Schema） | UniProtKB entry accession, e.g., P05067. |
| `compact` | boolean | None | Return a bounded summary instead of the complete UniProtKB JSON entry. Defaults to true to avoid oversized LLM outputs. Set compact=false only when you explicitly need the raw UniProtKB JSON. {&quot;default&quot;: true} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "accession": {
          "type": "string",
          "description": "UniProtKB entry accession, e.g., P05067."
        },
        "compact": {
          "type": "boolean",
          "description": "Return a bounded summary instead of the complete UniProtKB JSON entry. Defaults to true to avoid oversized LLM outputs. Set compact=false only when you explicitly need the raw UniProtKB JSON.",
          "default": true
        }
      },
      "required": [
        "accession"
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
            "entryType": {
              "type": "string"
            },
            "primaryAccession": {
              "type": "string"
            },
            "uniProtkbId": {
              "type": "string"
            },
            "protein_name": {
              "type": "string"
            },
            "gene_names": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "organism": {
              "type": "object",
              "properties": {
                "scientificName": {
                  "type": "string"
                },
                "commonName": {
                  "type": "string"
                },
                "taxonId": {
                  "type": [
                    "integer",
                    "number"
                  ]
                },
                "lineage": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  }
                }
              }
            },
            "sequence": {
              "type": "object",
              "properties": {
                "value": {
                  "type": "string"
                },
                "length": {
                  "type": [
                    "integer",
                    "number"
                  ]
                },
                "molWeight": {
                  "type": [
                    "integer",
                    "number"
                  ]
                },
                "crc64": {
                  "type": "string"
                },
                "md5": {
                  "type": "string"
                }
              }
            },
            "annotationScore": {
              "type": "number"
            },
            "proteinExistence": {
              "type": "string"
            },
            "keywords": {
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
                  "name": {
                    "type": "string"
                  }
                }
              }
            },
            "comments": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "commentType": {
                    "type": "string"
                  },
                  "texts": {
                    "type": "array",
                    "items": {
                      "type": "string"
                    }
                  }
                }
              }
            },
            "features": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "type": {
                    "type": "string"
                  },
                  "description": {
                    "type": "string"
                  },
                  "location": {
                    "type": "object",
                    "properties": {
                      "start": {
                        "type": "object",
                        "properties": {
                          "value": {
                            "type": [
                              "integer",
                              "number"
                            ]
                          },
                          "modifier": {
                            "type": "string"
                          }
                        }
                      },
                      "end": {
                        "type": "object",
                        "properties": {
                          "value": {
                            "type": [
                              "integer",
                              "number"
                            ]
                          },
                          "modifier": {
                            "type": "string"
                          }
                        }
                      }
                    }
                  }
                }
              }
            },
            "xref_count": {
              "type": [
                "integer",
                "number"
              ]
            }
          },
          "required": [
            "primaryAccession"
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

**执行与来源：** `network`；上游接口名 `UniProt_get_entry_by_accession`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("uni_prot_get_entry_by_accession", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.uni_prot_get_features_by_accession` {#dnakit.databases.uni_prot_get_features_by_accession}

Extract all sequence features (domains, active/binding sites, motifs, modified residues, signal peptides, transmembrane regions, etc.) from a UniProtKB entry. Returns each feature's type, start/end positions, and description. Use to map functional/structural regions of a protein or to see what is annotated at a residue of interest. Example: P04637 (TP53) returns binding sites, the DNA-binding domain, modified residues, and more.

```text
dnakit.databases.uni_prot_get_features_by_accession(*, accession, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string | 必填（输入 Schema） | UniProtKB accession, e.g., &#x27;P04637&#x27; (TP53) or &#x27;P05067&#x27; (APP). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "accession": {
          "type": "string",
          "description": "UniProtKB accession, e.g., 'P04637' (TP53) or 'P05067' (APP)."
        }
      },
      "required": [
        "accession"
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
              "type": {
                "type": "string"
              },
              "location": {
                "type": "object"
              },
              "description": {
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

**执行与来源：** `network`；上游接口名 `UniProt_get_features_by_accession`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("uni_prot_get_features_by_accession", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.uni_prot_get_function_by_accession` {#dnakit.databases.uni_prot_get_function_by_accession}

Extract functional annotations from UniProtKB entry (Comment type = FUNCTION).

```text
dnakit.databases.uni_prot_get_function_by_accession(*, accession, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string | 必填（输入 Schema） | UniProtKB accession, e.g., P05067. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "accession": {
          "type": "string",
          "description": "UniProtKB accession, e.g., P05067."
        }
      },
      "required": [
        "accession"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "array",
      "items": {
        "type": "string"
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `UniProt_get_function_by_accession`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("uni_prot_get_function_by_accession", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.uni_prot_get_sequence_by_accession` {#dnakit.databases.uni_prot_get_sequence_by_accession}

Extract the canonical sequence from UniProtKB entry.

```text
dnakit.databases.uni_prot_get_sequence_by_accession(*, accession, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string | 必填（输入 Schema） | UniProtKB accession, e.g., P05067. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "accession": {
          "type": "string",
          "description": "UniProtKB accession, e.g., P05067."
        }
      },
      "required": [
        "accession"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "string"
    }
    ```

**执行与来源：** `network`；上游接口名 `UniProt_get_sequence_by_accession`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("uni_prot_get_sequence_by_accession", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.uni_prot_search` {#dnakit.databases.uni_prot_search}

Search UniProtKB database with flexible query syntax. Returns protein entries with accession numbers and metadata. Query syntax supports: field searches (gene:TP53, taxonomy_id:9606, reviewed:true), ranges (length:[100 TO 500], mass:[20000 TO 50000]), wildcards (gene:MEIOB*), boolean operators (AND/OR/NOT), and parentheses for grouping. Examples: 'gene:TP53 AND taxonomy_id:9606', 'length:[400 TO 500] AND reviewed:true', 'tissue:brain NOT taxonomy_id:10090'. Scope by taxon with taxonomy_id (lineage match, 'Taxonomy [OC]') rather than organism_id ('Organism [OS]', exact organism node only): most bacterial entries are annotated at strain level, so a species-level organism_id silently returns 0 or an undercount (gene:gyrA AND organism_id:1423 -> 0 hits; taxonomy_id:1423 -> 1). A caller-supplied organism_id:<taxid> is rewritten to taxonomy_id:<taxid> and the change is reported in data.normalization_note.

```text
dnakit.databases.uni_prot_search(*, query, organism=None, limit=None, min_length=None, max_length=None, fields=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填（输入 Schema） | Search query using UniProt syntax. Simple: &#x27;MEIOB&#x27;, &#x27;insulin&#x27;. Field searches: &#x27;gene:TP53&#x27;, &#x27;protein_name:insulin&#x27;, &#x27;taxonomy_id:9606&#x27;, &#x27;reviewed:true&#x27;. Ranges: &#x27;length:&#91;100 TO 500&#93;&#x27;, &#x27;mass:&#91;20000 TO 50000&#93;&#x27;. Wildcards: &#x27;gene:MEIOB*&#x27;. Boolean: &#x27;gene:TP53 AND taxonomy_id:9606&#x27;, &#x27;tissue:brain OR tissue:liver&#x27;, &#x27;reviewed:true NOT fragment:true&#x27;. Use parentheses for grouping: &#x27;(taxonomy_id:9606 OR taxonomy_id:10090) AND gene:TP53&#x27;. Taxon scoping: use &#x27;taxonomy_id:&lt;taxid&gt;&#x27; (matches the whole lineage, so strain-level entries under a species taxid are included) or &#x27;organism_name:&quot;...&quot;&#x27;. A bare &#x27;organism:&#x27; auto-converts, and &#x27;organism_id:&lt;taxid&gt;&#x27; is rewritten to &#x27;taxonomy_id:&lt;taxid&gt;&#x27; with a note in data.normalization_note, because organism_id matches only the exact organism node and misses strain-level entries. |
| `organism` | string | None | Optional organism filter. Use common names (&#x27;human&#x27;, &#x27;mouse&#x27;, &#x27;rat&#x27;, &#x27;yeast&#x27;) or taxonomy ID (&#x27;9606&#x27;). Automatically combined with query using AND. Will not duplicate if organism is already in query. |
| `limit` | integer | None | Maximum number of results to return (default: 25, max: 500). Accepts string or integer. |
| `min_length` | integer | None | Minimum sequence length. Auto-converts to an open-ended length range query (min to unbounded). |
| `max_length` | integer | None | Maximum sequence length. Auto-converts to an open-ended length range query (unbounded to max). |
| `fields` | array | None | List of field names to return (e.g., &#91;&#x27;accession&#x27;,&#x27;gene_primary&#x27;,&#x27;length&#x27;,&#x27;organism_name&#x27;&#93;). When specified, returns raw API response with requested fields. Common fields: accession, id, gene_names, gene_primary, protein_name, organism_name, organism_id, length, mass, sequence, reviewed, cc_function. See UniProt API docs for full list. Default (no fields): returns formatted response with accession, id, protein_name, gene_names, organism, length. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
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
          "description": "Search query using UniProt syntax. Simple: 'MEIOB', 'insulin'. Field searches: 'gene:TP53', 'protein_name:insulin', 'taxonomy_id:9606', 'reviewed:true'. Ranges: 'length:[100 TO 500]', 'mass:[20000 TO 50000]'. Wildcards: 'gene:MEIOB*'. Boolean: 'gene:TP53 AND taxonomy_id:9606', 'tissue:brain OR tissue:liver', 'reviewed:true NOT fragment:true'. Use parentheses for grouping: '(taxonomy_id:9606 OR taxonomy_id:10090) AND gene:TP53'. Taxon scoping: use 'taxonomy_id:<taxid>' (matches the whole lineage, so strain-level entries under a species taxid are included) or 'organism_name:\"...\"'. A bare 'organism:' auto-converts, and 'organism_id:<taxid>' is rewritten to 'taxonomy_id:<taxid>' with a note in data.normalization_note, because organism_id matches only the exact organism node and misses strain-level entries."
        },
        "organism": {
          "type": "string",
          "description": "Optional organism filter. Use common names ('human', 'mouse', 'rat', 'yeast') or taxonomy ID ('9606'). Automatically combined with query using AND. Will not duplicate if organism is already in query."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of results to return (default: 25, max: 500). Accepts string or integer."
        },
        "min_length": {
          "type": "integer",
          "description": "Minimum sequence length. Auto-converts to an open-ended length range query (min to unbounded)."
        },
        "max_length": {
          "type": "integer",
          "description": "Maximum sequence length. Auto-converts to an open-ended length range query (unbounded to max)."
        },
        "fields": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "List of field names to return (e.g., ['accession','gene_primary','length','organism_name']). When specified, returns raw API response with requested fields. Common fields: accession, id, gene_names, gene_primary, protein_name, organism_name, organism_id, length, mass, sequence, reviewed, cc_function. See UniProt API docs for full list. Default (no fields): returns formatted response with accession, id, protein_name, gene_names, organism, length."
        }
      },
      "required": [
        "query"
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
            "total_results": {
              "type": [
                "integer",
                "number"
              ]
            },
            "returned": {
              "type": [
                "integer",
                "number"
              ]
            },
            "results": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "accession": {
                    "type": "string"
                  },
                  "id": {
                    "type": "string"
                  },
                  "protein_name": {
                    "type": "string"
                  },
                  "gene_names": {
                    "type": "array",
                    "items": {
                      "type": "string"
                    }
                  },
                  "organism": {
                    "type": "string"
                  },
                  "length": {
                    "type": [
                      "integer",
                      "number"
                    ]
                  }
                }
              }
            },
            "normalization_note": {
              "type": "string",
              "description": "Present only when the submitted query was rewritten. Currently emitted when an 'organism_id:<taxid>' clause was rewritten to 'taxonomy_id:<taxid>', widening the scope to include descendant strains."
            }
          },
          "required": [
            "total_results"
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

**执行与来源：** `network`；上游接口名 `UniProt_search`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("uni_prot_search", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.wiki_pathways_search` {#dnakit.databases.wiki_pathways_search}

Search WikiPathways for pathways by free-text title match. Returns WPIDs, titles, and organisms. Backed by the sparql.wikipathways.org SPARQL endpoint.

```text
dnakit.databases.wiki_pathways_search(*, query, organism=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填（输入 Schema） | Text to search, e.g., p53 |
| `organism` | string | None | Optional organism filter, e.g., &#x27;Homo sapiens&#x27; |
| `limit` | integer | None | 当前请求的返回数量上限；不是全库命中总数。 {&quot;default&quot;: 20} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回类型为 object，具体结构由所调用操作决定；不保证统一的 data 层。提供方错误会抛出异常。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "query": {
          "type": "string",
          "description": "Text to search, e.g., p53"
        },
        "organism": {
          "type": "string",
          "description": "Optional organism filter, e.g., 'Homo sapiens'"
        },
        "limit": {
          "type": "integer",
          "default": 20
        }
      },
      "required": [
        "query"
      ]
    }
    ```

**执行与来源：** `network`；上游接口名 `WikiPathways_search`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("wiki_pathways_search", verbose=False)
print(info["input_schema"])
```

