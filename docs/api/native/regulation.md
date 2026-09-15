# 表达、表观调控与肿瘤数据

本页收录 **65 个新增公开函数**。 [完整索引](index.md) · [调用约定与完整示例](usage.md)

每个函数均列出当前 Python 签名、参数、返回值和已登记的执行要求。提供方说明及字段描述保留源码原文。`None` 默认值在生成的数据库包装中通常表示不向服务传该参数；实际必填项与组合约束以输入 Schema 为准。坐标系和基因组版本按各接口说明使用。

ChIP-Atlas 按 gene 查峰时，未提供 peak_file 会返回 requires_data 和检索计划；提供已核对组装的 BED 后才计算区间重叠。富集入口返回 requires_submission 提交说明，不能把它当成已计算的富集统计。参数扩展分支的结果可能不同于提供方原始 Schema。

## 函数目录

- [`dnakit.databases.c_bio_portal_get_cancer_studies`](#dnakit.databases.c_bio_portal_get_cancer_studies)
- [`dnakit.databases.c_bio_portal_get_gene_info`](#dnakit.databases.c_bio_portal_get_gene_info)
- [`dnakit.databases.c_bio_portal_get_molecular_profiles`](#dnakit.databases.c_bio_portal_get_molecular_profiles)
- [`dnakit.databases.c_bio_portal_get_mutations`](#dnakit.databases.c_bio_portal_get_mutations)
- [`dnakit.databases.c_bio_portal_get_samples`](#dnakit.databases.c_bio_portal_get_samples)
- [`dnakit.databases.cellx_gene_download_h5ad`](#dnakit.databases.cellx_gene_download_h5ad)
- [`dnakit.databases.cellx_gene_get_cell_metadata`](#dnakit.databases.cellx_gene_get_cell_metadata)
- [`dnakit.databases.cellx_gene_get_embeddings`](#dnakit.databases.cellx_gene_get_embeddings)
- [`dnakit.databases.cellx_gene_get_expression_data`](#dnakit.databases.cellx_gene_get_expression_data)
- [`dnakit.databases.ch_ipatlas_enrichment_analysis`](#dnakit.databases.ch_ipatlas_enrichment_analysis)
- [`dnakit.databases.ch_ipatlas_get_experiments`](#dnakit.databases.ch_ipatlas_get_experiments)
- [`dnakit.databases.ch_ipatlas_get_peak_data`](#dnakit.databases.ch_ipatlas_get_peak_data)
- [`dnakit.databases.ch_ipatlas_search_datasets`](#dnakit.databases.ch_ipatlas_search_datasets)
- [`dnakit.databases.dep_map_get_gene_dependencies`](#dnakit.databases.dep_map_get_gene_dependencies)
- [`dnakit.databases.dep_map_search_cell_lines`](#dnakit.databases.dep_map_search_cell_lines)
- [`dnakit.databases.e_qtl_get_associations`](#dnakit.databases.e_qtl_get_associations)
- [`dnakit.databases.e_qtl_list_datasets`](#dnakit.databases.e_qtl_list_datasets)
- [`dnakit.databases.encode_get_biosample`](#dnakit.databases.encode_get_biosample)
- [`dnakit.databases.encode_get_chromatin_state`](#dnakit.databases.encode_get_chromatin_state)
- [`dnakit.databases.encode_get_experiment`](#dnakit.databases.encode_get_experiment)
- [`dnakit.databases.encode_list_files`](#dnakit.databases.encode_list_files)
- [`dnakit.databases.encode_search_annotations`](#dnakit.databases.encode_search_annotations)
- [`dnakit.databases.encode_search_biosamples`](#dnakit.databases.encode_search_biosamples)
- [`dnakit.databases.encode_search_chromatin_accessibility`](#dnakit.databases.encode_search_chromatin_accessibility)
- [`dnakit.databases.encode_search_experiments`](#dnakit.databases.encode_search_experiments)
- [`dnakit.databases.encode_search_histone_experiments`](#dnakit.databases.encode_search_histone_experiments)
- [`dnakit.databases.encode_search_rnaseq_experiments`](#dnakit.databases.encode_search_rnaseq_experiments)
- [`dnakit.databases.expression_atlas_search_differential`](#dnakit.databases.expression_atlas_search_differential)
- [`dnakit.databases.expression_atlas_search_experiments`](#dnakit.databases.expression_atlas_search_experiments)
- [`dnakit.databases.four_dn_search_data`](#dnakit.databases.four_dn_search_data)
- [`dnakit.databases.gdc_get_clinical_data`](#dnakit.databases.gdc_get_clinical_data)
- [`dnakit.databases.gdc_get_mutation_frequency`](#dnakit.databases.gdc_get_mutation_frequency)
- [`dnakit.databases.gdc_get_ssm_by_gene`](#dnakit.databases.gdc_get_ssm_by_gene)
- [`dnakit.databases.gdc_get_survival`](#dnakit.databases.gdc_get_survival)
- [`dnakit.databases.gdc_list_files`](#dnakit.databases.gdc_list_files)
- [`dnakit.databases.gdc_list_projects`](#dnakit.databases.gdc_list_projects)
- [`dnakit.databases.gdc_search_cases`](#dnakit.databases.gdc_search_cases)
- [`dnakit.databases.geo_search_atacseq_datasets`](#dnakit.databases.geo_search_atacseq_datasets)
- [`dnakit.databases.geo_search_chipseq_datasets`](#dnakit.databases.geo_search_chipseq_datasets)
- [`dnakit.databases.geo_search_rnaseq_datasets`](#dnakit.databases.geo_search_rnaseq_datasets)
- [`dnakit.databases.gtex_calculate_eqtl`](#dnakit.databases.gtex_calculate_eqtl)
- [`dnakit.databases.gtex_get_expression_summary`](#dnakit.databases.gtex_get_expression_summary)
- [`dnakit.databases.gtex_get_median_gene_expression`](#dnakit.databases.gtex_get_median_gene_expression)
- [`dnakit.databases.gtex_get_multi_tissue_eqtls`](#dnakit.databases.gtex_get_multi_tissue_eqtls)
- [`dnakit.databases.gtex_get_single_tissue_eqtls`](#dnakit.databases.gtex_get_single_tissue_eqtls)
- [`dnakit.databases.gtex_get_tissue_sites`](#dnakit.databases.gtex_get_tissue_sites)
- [`dnakit.databases.gtex_query_eqtl`](#dnakit.databases.gtex_query_eqtl)
- [`dnakit.databases.harmonizome_get_gene`](#dnakit.databases.harmonizome_get_gene)
- [`dnakit.databases.hpa_get_cancer_prognostics_by_gene`](#dnakit.databases.hpa_get_cancer_prognostics_by_gene)
- [`dnakit.databases.hpa_get_protein_interactions_by_gene`](#dnakit.databases.hpa_get_protein_interactions_by_gene)
- [`dnakit.databases.hpa_get_rna_expression_by_source`](#dnakit.databases.hpa_get_rna_expression_by_source)
- [`dnakit.databases.hpa_get_rna_expression_in_specific_tissues`](#dnakit.databases.hpa_get_rna_expression_in_specific_tissues)
- [`dnakit.databases.hpa_get_subcellular_location`](#dnakit.databases.hpa_get_subcellular_location)
- [`dnakit.databases.hpa_search_genes_by_query`](#dnakit.databases.hpa_search_genes_by_query)
- [`dnakit.databases.jaspar_get_matrix`](#dnakit.databases.jaspar_get_matrix)
- [`dnakit.databases.jaspar_get_matrix_versions`](#dnakit.databases.jaspar_get_matrix_versions)
- [`dnakit.databases.jaspar_get_transcription_factors`](#dnakit.databases.jaspar_get_transcription_factors)
- [`dnakit.databases.jaspar_search_matrices`](#dnakit.databases.jaspar_search_matrices)
- [`dnakit.databases.re_map_get_transcription_factor_binding`](#dnakit.databases.re_map_get_transcription_factor_binding)
- [`dnakit.databases.regulome_db_query_variant`](#dnakit.databases.regulome_db_query_variant)
- [`dnakit.databases.screen_get_regulatory_elements`](#dnakit.databases.screen_get_regulatory_elements)
- [`dnakit.databases.ucsc_get_encode_c_cres`](#dnakit.databases.ucsc_get_encode_c_cres)
- [`dnakit.databases.uni_bind_get_dataset`](#dnakit.databases.uni_bind_get_dataset)
- [`dnakit.databases.uni_bind_list_tfs`](#dnakit.databases.uni_bind_list_tfs)
- [`dnakit.databases.uni_bind_search_datasets`](#dnakit.databases.uni_bind_search_datasets)

## `dnakit.databases.c_bio_portal_get_cancer_studies` {#dnakit.databases.c_bio_portal_get_cancer_studies}

Get one page of the cBioPortal study catalogue. This returns a SLICE, not the whole catalogue: `limit` defaults to 20 while cBioPortal currently holds several hundred studies. The response always reports `total_available` (the full catalogue size) alongside `count` (studies actually returned), and sets `truncated: true` with a `truncation_note` whenever studies were left behind. Never conclude a study is missing from cBioPortal from a truncated page — raise `limit` to `total_available`, or page with `offset`, before drawing that conclusion.

```text
dnakit.databases.c_bio_portal_get_cancer_studies(*, limit=None, offset=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `limit` | integer | None | Number of studies to return in this page. The default of 20 covers only a small fraction of the catalogue; set it to the `total_available` reported in the response to retrieve every study. {&quot;default&quot;: 20} |
| `offset` | integer | None | 0-based index of the first study to return, for paging through the catalogue. Applied client-side because the cBioPortal /studies endpoint ignores its own pageNumber parameter. {&quot;default&quot;: 0} |
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
          "default": 20,
          "description": "Number of studies to return in this page. The default of 20 covers only a small fraction of the catalogue; set it to the `total_available` reported in the response to retrieve every study."
        },
        "offset": {
          "type": "integer",
          "default": 0,
          "description": "0-based index of the first study to return, for paging through the catalogue. Applied client-side because the cBioPortal /studies endpoint ignores its own pageNumber parameter."
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
          "studyId": {
            "type": "string"
          },
          "name": {
            "type": "string"
          },
          "description": {
            "type": "string"
          },
          "cancerTypeId": {
            "type": "string"
          },
          "publicStudy": {
            "type": "boolean"
          },
          "pmid": {
            "type": "string"
          },
          "citation": {
            "type": "string"
          },
          "groups": {
            "type": [
              "array",
              "string"
            ]
          },
          "status": {
            "type": [
              "string",
              "integer"
            ]
          },
          "importDate": {
            "type": "string"
          },
          "allSampleCount": {
            "type": "integer"
          },
          "readPermission": {
            "type": "boolean"
          },
          "referenceGenome": {
            "type": "string"
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `cBioPortal_get_cancer_studies`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("c_bio_portal_get_cancer_studies", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.c_bio_portal_get_gene_info` {#dnakit.databases.c_bio_portal_get_gene_info}

Get basic information about a specific gene by Entrez Gene ID. Returns the HUGO gene symbol and gene type (e.g., 'protein-coding'); does not return aliases -- use HGNC_fetch_gene_by_symbol or NCBIGene_get_summary for gene aliases/chromosome location.

```text
dnakit.databases.c_bio_portal_get_gene_info(*, entrez_gene_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `entrez_gene_id` | integer | 必填（输入 Schema） | Entrez Gene ID (e.g., 672 for BRCA1) |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "entrez_gene_id": {
          "type": "integer",
          "description": "Entrez Gene ID (e.g., 672 for BRCA1)"
        }
      },
      "required": [
        "entrez_gene_id"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "entrezGeneId": {
          "type": "integer"
        },
        "hugoGeneSymbol": {
          "type": "string"
        },
        "type": {
          "type": "string"
        },
        "cytoband": {
          "type": "string"
        },
        "length": {
          "type": "integer"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `cBioPortal_get_gene_info`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("c_bio_portal_get_gene_info", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.c_bio_portal_get_molecular_profiles` {#dnakit.databases.c_bio_portal_get_molecular_profiles}

Get molecular profiles for a cancer study. Molecular profiles include mutation data, copy number alterations, mRNA expression, etc.

```text
dnakit.databases.c_bio_portal_get_molecular_profiles(*, study_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `study_id` | string | 必填（输入 Schema） | Cancer study ID |
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
          "description": "Cancer study ID"
        }
      },
      "required": [
        "study_id"
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
          "molecularProfileId": {
            "type": "string"
          },
          "molecularAlterationType": {
            "type": "string"
          },
          "datatype": {
            "type": "string"
          },
          "name": {
            "type": "string"
          },
          "description": {
            "type": "string"
          },
          "showProfileInAnalysisTab": {
            "type": "boolean"
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `cBioPortal_get_molecular_profiles`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("c_bio_portal_get_molecular_profiles", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.c_bio_portal_get_mutations` {#dnakit.databases.c_bio_portal_get_mutations}

Get mutation data for specific genes in a cancer study. This uses the updated cBioPortal API that requires molecular profiles.

```text
dnakit.databases.c_bio_portal_get_mutations(*, study_id, gene_list, sample_list_id=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `study_id` | string | 必填（输入 Schema） | Cancer study ID (e.g., &#x27;brca_tcga&#x27;) |
| `gene_list` | string | 必填（输入 Schema） | Comma-separated gene symbols (e.g., &#x27;BRCA1,BRCA2&#x27;) |
| `sample_list_id` | string | None | Optional sample list ID. If not provided, uses all samples in the study. |
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
          "description": "Cancer study ID (e.g., 'brca_tcga')"
        },
        "gene_list": {
          "type": "string",
          "description": "Comma-separated gene symbols (e.g., 'BRCA1,BRCA2')"
        },
        "sample_list_id": {
          "type": "string",
          "description": "Optional sample list ID. If not provided, uses all samples in the study."
        }
      },
      "required": [
        "gene_list",
        "study_id"
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
          "molecularProfileId": {
            "type": "string"
          },
          "sampleId": {
            "type": "string"
          },
          "patientId": {
            "type": "string"
          },
          "studyId": {
            "type": "string"
          },
          "entrezGeneId": {
            "type": "integer"
          },
          "gene": {
            "type": "object",
            "properties": {
              "entrezGeneId": {
                "type": "integer"
              },
              "hugoGeneSymbol": {
                "type": "string"
              },
              "type": {
                "type": "string"
              }
            }
          },
          "mutationType": {
            "type": "string"
          },
          "proteinChange": {
            "type": "string"
          },
          "mutationStatus": {
            "type": "string"
          },
          "validationStatus": {
            "type": "string"
          },
          "startPosition": {
            "type": "integer"
          },
          "endPosition": {
            "type": "integer"
          },
          "referenceAllele": {
            "type": "string"
          },
          "variantAllele": {
            "type": "string"
          },
          "keyword": {
            "type": "string"
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `cBioPortal_get_mutations`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("c_bio_portal_get_mutations", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.c_bio_portal_get_samples` {#dnakit.databases.c_bio_portal_get_samples}

Get all samples in a cancer study. Returns the full cohort by default; the response is a bare array with no truncation marker, so lower `page_size` only when you deliberately want a sample. Use this count as the denominator for mutation-frequency calculations.

```text
dnakit.databases.c_bio_portal_get_samples(*, study_id, page_size=None, page_number=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `study_id` | string | 必填（输入 Schema） | Cancer study ID |
| `page_size` | integer | None | Maximum number of samples to return. Defaults high enough to return every sample in a study; a smaller value silently truncates the cohort, which will corrupt any frequency denominator computed from it. {&quot;default&quot;: 100000} |
| `page_number` | integer | None | 0-based page index, used with `page_size` to page through very large studies. {&quot;default&quot;: 0} |
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
          "description": "Cancer study ID"
        },
        "page_size": {
          "type": "integer",
          "default": 100000,
          "description": "Maximum number of samples to return. Defaults high enough to return every sample in a study; a smaller value silently truncates the cohort, which will corrupt any frequency denominator computed from it."
        },
        "page_number": {
          "type": "integer",
          "default": 0,
          "description": "0-based page index, used with `page_size` to page through very large studies."
        }
      },
      "required": [
        "study_id"
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
          "sampleId": {
            "type": "string"
          },
          "sampleType": {
            "type": "string"
          },
          "patientId": {
            "type": "string"
          },
          "studyId": {
            "type": "string"
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `cBioPortal_get_samples`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("c_bio_portal_get_samples", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.cellx_gene_download_h5ad` {#dnakit.databases.cellx_gene_download_h5ad}

Download original H5AD (HDF5-based AnnData) files from CELLxGENE datasets or get their URIs. Access unprocessed data from individual studies with complete metadata and analysis results. Requires dataset_id from CELLxGENE Discover. Note: Files can be large (GBs). Use for: accessing raw data, reproducibility, offline analysis, custom processing pipelines. Note: requires the cellxgene-census Python package; install with `pip install cellxgene-census` for live access.

```text
dnakit.databases.cellx_gene_download_h5ad(*, operation=None, dataset_id, output_path=None, census_version=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;download_h5ad&quot;&#93;, &quot;default&quot;: &quot;download_h5ad&quot;} |
| `dataset_id` | string | 必填（输入 Schema） | CELLxGENE dataset identifier (required) |
| `output_path` | string | None | Local path to save H5AD file (omit to get URI only) |
| `census_version` | string | None | Census version to query. &#x27;stable&#x27; (recommended, Long-Term Support release), &#x27;latest&#x27; (newest data, may change), or specific date &#x27;YYYY-MM-DD&#x27; (for reproducibility). Default &#x27;stable&#x27; is best for production analyses. {&quot;default&quot;: &quot;stable&quot;} |
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
            "download_h5ad"
          ],
          "default": "download_h5ad"
        },
        "dataset_id": {
          "type": "string",
          "description": "CELLxGENE dataset identifier (required)"
        },
        "output_path": {
          "type": "string",
          "description": "Local path to save H5AD file (omit to get URI only)"
        },
        "census_version": {
          "type": "string",
          "description": "Census version to query. 'stable' (recommended, Long-Term Support release), 'latest' (newest data, may change), or specific date 'YYYY-MM-DD' (for reproducibility). Default 'stable' is best for production analyses.",
          "default": "stable"
        }
      },
      "required": [
        "dataset_id"
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
        "dataset_id": {
          "type": "string"
        },
        "output_path": {
          "type": [
            "string",
            "null"
          ]
        },
        "uri": {
          "type": [
            "string",
            "null"
          ]
        },
        "message": {
          "type": "string"
        }
      },
      "additionalProperties": true
    }
    ```

**执行与来源：** `filesystem`；上游接口名 `CELLxGENE_download_h5ad`。

**提供方凭据：** `CELLXGENE_CENSUS_PACKAGE_INSTALLED`

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("cellx_gene_download_h5ad", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.cellx_gene_get_cell_metadata` {#dnakit.databases.cellx_gene_get_cell_metadata}

Query cell metadata from CELLxGENE Census (50M+ human/mouse single cells). CRITICAL: obs_value_filter is REQUIRED - unfiltered queries will timeout. Common filters: tissue_general (brain, heart, lung, liver, blood, kidney, skin, eye, intestine), cell_type (T cell, B cell, neuron, macrophage, fibroblast), disease (normal, COVID-19, Alzheimer disease, cancer). Start with: obs_value_filter='tissue_general == "lung" and disease == "normal"'. Returns cell type, tissue, disease, donor info. Use for: finding cells by criteria, cohort selection, exploring cell populations. Note: requires the cellxgene-census Python package; install with `pip install cellxgene-census` for live access.

```text
dnakit.databases.cellx_gene_get_cell_metadata(*, operation=None, organism=None, obs_value_filter=None, column_names=None, census_version=None, gene=None, tissue=None, max_cells=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;get_obs_metadata&quot;&#93;, &quot;default&quot;: &quot;get_obs_metadata&quot;} |
| `organism` | string | None | Organism name {&quot;enum&quot;: &#91;&quot;Homo sapiens&quot;, &quot;Mus musculus&quot;&#93;, &quot;default&quot;: &quot;Homo sapiens&quot;} |
| `obs_value_filter` | string | None | REQUIRED - filter cells using SQL-like syntax. Unfiltered queries timeout (50M+ cells). Format: &#x27;field == &quot;value&quot;&#x27;. Operators: ==, !=, in, &lt;, &gt;, &lt;=, &gt;=. Combine with &#x27;and&#x27;/&#x27;or&#x27;. Common fields and values: tissue_general (brain, heart, lung, liver, blood, kidney, skin, eye, intestine, pancreas, bone marrow), cell_type (T cell, B cell, neuron, macrophage, fibroblast, epithelial cell, endothelial cell), disease (normal, COVID-19, Alzheimer disease, lung adenocarcinoma), assay (10x 3&#x27; v3, Smart-seq2), sex (male, female). Examples: &#x27;tissue_general == &quot;lung&quot;&#x27;, &#x27;disease == &quot;COVID-19&quot; and tissue_general == &quot;blood&quot;&#x27;, &#x27;cell_type in &#91;&quot;T cell&quot;, &quot;B cell&quot;&#93; and disease == &quot;normal&quot;&#x27;. |
| `column_names` | array | None | Specific columns to return (default: all columns) {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `census_version` | string | None | Census version to query. &#x27;stable&#x27; (recommended, Long-Term Support release), &#x27;latest&#x27; (newest data, may change), or specific date &#x27;YYYY-MM-DD&#x27; (for reproducibility). Default &#x27;stable&#x27; is best for production analyses. {&quot;default&quot;: &quot;stable&quot;} |
| `gene` | string | None | 基因标识或基因列表；命名要求和类型见当前接口。 |
| `tissue` | string | None | 目标组织名称或标识，按当前数据源填写。 |
| `max_cells` | integer | None | 本次返回或处理的细胞数上限。 {&quot;minimum&quot;: 1, &quot;maximum&quot;: 10000} |
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
            "get_obs_metadata"
          ],
          "default": "get_obs_metadata"
        },
        "organism": {
          "type": "string",
          "description": "Organism name",
          "enum": [
            "Homo sapiens",
            "Mus musculus"
          ],
          "default": "Homo sapiens"
        },
        "obs_value_filter": {
          "type": "string",
          "description": "REQUIRED - filter cells using SQL-like syntax. Unfiltered queries timeout (50M+ cells). Format: 'field == \"value\"'. Operators: ==, !=, in, <, >, <=, >=. Combine with 'and'/'or'. Common fields and values: tissue_general (brain, heart, lung, liver, blood, kidney, skin, eye, intestine, pancreas, bone marrow), cell_type (T cell, B cell, neuron, macrophage, fibroblast, epithelial cell, endothelial cell), disease (normal, COVID-19, Alzheimer disease, lung adenocarcinoma), assay (10x 3' v3, Smart-seq2), sex (male, female). Examples: 'tissue_general == \"lung\"', 'disease == \"COVID-19\" and tissue_general == \"blood\"', 'cell_type in [\"T cell\", \"B cell\"] and disease == \"normal\"'."
        },
        "column_names": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "Specific columns to return (default: all columns)"
        },
        "census_version": {
          "type": "string",
          "description": "Census version to query. 'stable' (recommended, Long-Term Support release), 'latest' (newest data, may change), or specific date 'YYYY-MM-DD' (for reproducibility). Default 'stable' is best for production analyses.",
          "default": "stable"
        },
        "gene": {
          "type": "string"
        },
        "tissue": {
          "type": "string"
        },
        "max_cells": {
          "type": "integer",
          "minimum": 1,
          "maximum": 10000
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": [
        "array",
        "null"
      ],
      "items": {
        "type": "object",
        "properties": {
          "cell_type": {
            "type": "string"
          },
          "tissue_general": {
            "type": "string"
          },
          "disease": {
            "type": "string"
          },
          "assay": {
            "type": "string"
          },
          "sex": {
            "type": "string"
          },
          "donor_id": {
            "type": "string"
          },
          "dataset_id": {
            "type": "string"
          },
          "soma_joinid": {
            "type": "integer"
          }
        },
        "additionalProperties": true
      }
    }
    ```

??? info "DNAKit 参数扩展规则"

    ```json
    {
      "parameters": {
        "gene": {
          "type": "string"
        },
        "tissue": {
          "type": "string"
        },
        "max_cells": {
          "type": "integer",
          "minimum": 1,
          "maximum": 10000
        }
      },
      "optional": [
        "obs_value_filter"
      ]
    }
    ```

**执行与来源：** `network`；上游接口名 `CELLxGENE_get_cell_metadata`。

**提供方凭据：** `CELLXGENE_CENSUS_PACKAGE_INSTALLED`

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("cellx_gene_get_cell_metadata", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.cellx_gene_get_embeddings` {#dnakit.databases.cellx_gene_get_embeddings}

Access pre-calculated cell embeddings (scVI, Geneformer) from CELLxGENE Census. Returns available embedding names or specific embedding data. Embeddings enable rapid visualization and analysis without recomputation. Use for: UMAP/t-SNE visualization, cell similarity analysis, transfer learning. Note: requires the cellxgene-census Python package; install with `pip install cellxgene-census` for live access.

```text
dnakit.databases.cellx_gene_get_embeddings(*, operation=None, organism=None, embedding_name=None, census_version=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;get_embeddings&quot;&#93;, &quot;default&quot;: &quot;get_embeddings&quot;} |
| `organism` | string | None | Organism name {&quot;enum&quot;: &#91;&quot;Homo sapiens&quot;, &quot;Mus musculus&quot;&#93;, &quot;default&quot;: &quot;Homo sapiens&quot;} |
| `embedding_name` | string | None | Name of specific embedding to retrieve (omit to list available embeddings) |
| `census_version` | string | None | Census version to query. &#x27;stable&#x27; (recommended, Long-Term Support release), &#x27;latest&#x27; (newest data, may change), or specific date &#x27;YYYY-MM-DD&#x27; (for reproducibility). Default &#x27;stable&#x27; is best for production analyses. {&quot;default&quot;: &quot;stable&quot;} |
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
            "get_embeddings"
          ],
          "default": "get_embeddings"
        },
        "organism": {
          "type": "string",
          "description": "Organism name",
          "enum": [
            "Homo sapiens",
            "Mus musculus"
          ],
          "default": "Homo sapiens"
        },
        "embedding_name": {
          "type": "string",
          "description": "Name of specific embedding to retrieve (omit to list available embeddings)"
        },
        "census_version": {
          "type": "string",
          "description": "Census version to query. 'stable' (recommended, Long-Term Support release), 'latest' (newest data, may change), or specific date 'YYYY-MM-DD' (for reproducibility). Default 'stable' is best for production analyses.",
          "default": "stable"
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
        "message": {
          "type": "string"
        },
        "organism": {
          "type": [
            "string",
            "null"
          ]
        },
        "embedding_name": {
          "type": [
            "string",
            "null"
          ]
        },
        "available": {
          "type": [
            "boolean",
            "null"
          ]
        },
        "shape": {
          "type": [
            "array",
            "null"
          ]
        },
        "available_embeddings": {
          "type": [
            "object",
            "null"
          ]
        }
      },
      "additionalProperties": true
    }
    ```

**执行与来源：** `network`；上游接口名 `CELLxGENE_get_embeddings`。

**提供方凭据：** `CELLXGENE_CENSUS_PACKAGE_INSTALLED`

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("cellx_gene_get_embeddings", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.cellx_gene_get_expression_data` {#dnakit.databases.cellx_gene_get_expression_data}

Query gene expression data from CELLxGENE Census (50M+ cells, 60K+ genes). CRITICAL: At least one filter REQUIRED - unfiltered queries timeout. Start with: obs_value_filter='tissue_general == "lung" and disease == "normal"' AND var_value_filter='feature_name in ["TP53", "EGFR"]'. Returns AnnData summary (dimensions, metadata). For full matrix data, use Python API with sufficient memory. Note: requires the cellxgene-census Python package; install with `pip install cellxgene-census` for live access.

```text
dnakit.databases.cellx_gene_get_expression_data(*, operation=None, organism=None, obs_value_filter=None, var_value_filter=None, obs_column_names=None, var_column_names=None, census_version=None, gene=None, tissue=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;get_anndata&quot;&#93;, &quot;default&quot;: &quot;get_anndata&quot;} |
| `organism` | string | None | Organism name {&quot;enum&quot;: &#91;&quot;Homo sapiens&quot;, &quot;Mus musculus&quot;&#93;, &quot;default&quot;: &quot;Homo sapiens&quot;} |
| `obs_value_filter` | string | None | REQUIRED (or use var_value_filter) - Cell filter. Common values: tissue_general (brain, lung, heart, blood, liver), cell_type (T cell, B cell, neuron, macrophage), disease (normal, COVID-19, Alzheimer disease). Example: &#x27;tissue_general == &quot;lung&quot; and disease == &quot;normal&quot;&#x27; |
| `var_value_filter` | string | None | REQUIRED (or use obs_value_filter) - Gene filter by symbol or Ensembl ID. Examples: &#x27;feature_name in &#91;&quot;TP53&quot;, &quot;BRCA1&quot;, &quot;EGFR&quot;&#93;&#x27;, &#x27;feature_name == &quot;CD4&quot;&#x27;. Common genes: TP53, EGFR, BRCA1, CD4, CD8A, IL6, TNF, GAPDH. |
| `obs_column_names` | array | None | Cell metadata columns to include {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `var_column_names` | array | None | Gene metadata columns to include {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `census_version` | string | None | Census version to query. &#x27;stable&#x27; (recommended, Long-Term Support release), &#x27;latest&#x27; (newest data, may change), or specific date &#x27;YYYY-MM-DD&#x27; (for reproducibility). Default &#x27;stable&#x27; is best for production analyses. {&quot;default&quot;: &quot;stable&quot;} |
| `gene` | string | None | 基因标识或基因列表；命名要求和类型见当前接口。 |
| `tissue` | string | None | 目标组织名称或标识，按当前数据源填写。 |
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
            "get_anndata"
          ],
          "default": "get_anndata"
        },
        "organism": {
          "type": "string",
          "description": "Organism name",
          "enum": [
            "Homo sapiens",
            "Mus musculus"
          ],
          "default": "Homo sapiens"
        },
        "obs_value_filter": {
          "type": "string",
          "description": "REQUIRED (or use var_value_filter) - Cell filter. Common values: tissue_general (brain, lung, heart, blood, liver), cell_type (T cell, B cell, neuron, macrophage), disease (normal, COVID-19, Alzheimer disease). Example: 'tissue_general == \"lung\" and disease == \"normal\"'"
        },
        "var_value_filter": {
          "type": "string",
          "description": "REQUIRED (or use obs_value_filter) - Gene filter by symbol or Ensembl ID. Examples: 'feature_name in [\"TP53\", \"BRCA1\", \"EGFR\"]', 'feature_name == \"CD4\"'. Common genes: TP53, EGFR, BRCA1, CD4, CD8A, IL6, TNF, GAPDH."
        },
        "obs_column_names": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "Cell metadata columns to include"
        },
        "var_column_names": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "Gene metadata columns to include"
        },
        "census_version": {
          "type": "string",
          "description": "Census version to query. 'stable' (recommended, Long-Term Support release), 'latest' (newest data, may change), or specific date 'YYYY-MM-DD' (for reproducibility). Default 'stable' is best for production analyses.",
          "default": "stable"
        },
        "gene": {
          "type": "string"
        },
        "tissue": {
          "type": "string"
        }
      },
      "required": []
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
        "n_obs": {
          "type": "integer"
        },
        "n_vars": {
          "type": "integer"
        },
        "obs_columns": {
          "type": "array"
        },
        "var_columns": {
          "type": "array"
        },
        "message": {
          "type": "string"
        }
      },
      "additionalProperties": true
    }
    ```

??? info "DNAKit 参数扩展规则"

    ```json
    {
      "parameters": {
        "gene": {
          "type": "string"
        },
        "tissue": {
          "type": "string"
        }
      },
      "handler": "census_filters"
    }
    ```

**执行与来源：** `network`；上游接口名 `CELLxGENE_get_expression_data`。

**提供方凭据：** `CELLXGENE_CENSUS_PACKAGE_INSTALLED`

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("cellx_gene_get_expression_data", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ch_ipatlas_enrichment_analysis` {#dnakit.databases.ch_ipatlas_enrichment_analysis}

根据 BED、motif 或基因列表生成 ChIP-Atlas 富集任务的提交说明；本函数不计算富集统计量。

```text
dnakit.databases.ch_ipatlas_enrichment_analysis(*, operation=None, bed_data=None, motif=None, gene_list=None, genome=None, antigen_class=None, cell_type_class=None, threshold=None, distance=None, gene=None, cell_type=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | 提供方操作模式；默认值与允许值见本行约束。 {&quot;enum&quot;: &#91;&quot;enrichment_analysis&quot;&#93;, &quot;default&quot;: &quot;enrichment_analysis&quot;} |
| `bed_data` | string | None | **Option 1**: BED format genomic regions (tab-separated: chr, start, end). For finding proteins bound to specific genomic regions. Example: &#x27;chr1\t1000\t2000\nchr2\t3000\t4000&#x27;. |
| `motif` | string | None | **Option 2**: DNA sequence motif in IUPAC notation. Use: A/T/G/C (bases), W=A&#124;T, S=G&#124;C, M=A&#124;C, K=G&#124;T, R=A&#124;G, Y=C&#124;T, B=C&#124;G&#124;T, D=A&#124;G&#124;T, H=A&#124;C&#124;T, V=A&#124;C&#124;G, N=any. For finding proteins bound to specific DNA sequences. Example: &#x27;CANNTG&#x27; (E-box motif). |
| `gene_list` | array &#124; string | None | **Option 3**: Gene symbols (HGNC for human, MGI for mouse, RGD for rat, FlyBase, WormBase, SGD for yeast). Provide as array or single gene. For finding transcription factors regulating genes. Example: &#91;&#x27;TP53&#x27;, &#x27;MDM2&#x27;, &#x27;CDKN1A&#x27;&#93;. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `genome` | string | None | Genome assembly {&quot;enum&quot;: &#91;&quot;hg38&quot;, &quot;hg19&quot;, &quot;mm10&quot;, &quot;mm9&quot;, &quot;rn6&quot;, &quot;dm6&quot;, &quot;dm3&quot;, &quot;ce11&quot;, &quot;ce10&quot;, &quot;sacCer3&quot;&#93;, &quot;default&quot;: &quot;hg38&quot;} |
| `antigen_class` | string | None | Filter by antigen class (e.g., &#x27;TFs and others&#x27;, &#x27;Histone&#x27;, &#x27;RNA polymerase&#x27;) |
| `cell_type_class` | string | None | Filter by cell type class (e.g., &#x27;Blood&#x27;, &#x27;Liver&#x27;, &#x27;Brain&#x27;) |
| `threshold` | string | None | Peak calling stringency (MACS2 Q-value). Options: &#x27;05&#x27;=1e-5 (permissive, more peaks, broader features), &#x27;10&#x27;=1e-10 (moderate, balanced), &#x27;20&#x27;=1e-20 (strict, high confidence only, narrow peaks). Default &#x27;05&#x27; suitable for most analyses. Higher values = fewer but more confident peaks. {&quot;enum&quot;: &#91;&quot;05&quot;, &quot;10&quot;, &quot;20&quot;&#93;, &quot;default&quot;: &quot;05&quot;} |
| `distance` | string | None | Distance from Transcription Start Site (TSS) in base pairs for gene-TF association. Defines promoter region. Default 5000 (±5kb, captures typical promoters). Use 1000-2000 for narrow promoters, 10000+ for enhancer regions. {&quot;default&quot;: &quot;5000&quot;} |
| `gene` | array &#124; string | None | DNAKit optional alias for gene_list {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `cell_type` | string | None | DNAKit optional alias for cell_type_class |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

字典包含 status=requires_submission、analysis_completed=False 和提交说明；不会返回已完成的富集检验。下面保留提供方原始结果 Schema 供追溯，不能据此认定分析已经完成。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "enum": [
            "enrichment_analysis"
          ],
          "default": "enrichment_analysis"
        },
        "bed_data": {
          "type": "string",
          "description": "**Option 1**: BED format genomic regions (tab-separated: chr, start, end). For finding proteins bound to specific genomic regions. Example: 'chr1\\t1000\\t2000\\nchr2\\t3000\\t4000'."
        },
        "motif": {
          "type": "string",
          "description": "**Option 2**: DNA sequence motif in IUPAC notation. Use: A/T/G/C (bases), W=A|T, S=G|C, M=A|C, K=G|T, R=A|G, Y=C|T, B=C|G|T, D=A|G|T, H=A|C|T, V=A|C|G, N=any. For finding proteins bound to specific DNA sequences. Example: 'CANNTG' (E-box motif)."
        },
        "gene_list": {
          "type": [
            "array",
            "string"
          ],
          "items": {
            "type": "string"
          },
          "description": "**Option 3**: Gene symbols (HGNC for human, MGI for mouse, RGD for rat, FlyBase, WormBase, SGD for yeast). Provide as array or single gene. For finding transcription factors regulating genes. Example: ['TP53', 'MDM2', 'CDKN1A']."
        },
        "genome": {
          "type": "string",
          "description": "Genome assembly",
          "enum": [
            "hg38",
            "hg19",
            "mm10",
            "mm9",
            "rn6",
            "dm6",
            "dm3",
            "ce11",
            "ce10",
            "sacCer3"
          ],
          "default": "hg38"
        },
        "antigen_class": {
          "type": "string",
          "description": "Filter by antigen class (e.g., 'TFs and others', 'Histone', 'RNA polymerase')"
        },
        "cell_type_class": {
          "type": "string",
          "description": "Filter by cell type class (e.g., 'Blood', 'Liver', 'Brain')"
        },
        "threshold": {
          "type": "string",
          "description": "Peak calling stringency (MACS2 Q-value). Options: '05'=1e-5 (permissive, more peaks, broader features), '10'=1e-10 (moderate, balanced), '20'=1e-20 (strict, high confidence only, narrow peaks). Default '05' suitable for most analyses. Higher values = fewer but more confident peaks.",
          "enum": [
            "05",
            "10",
            "20"
          ],
          "default": "05"
        },
        "distance": {
          "type": "string",
          "description": "Distance from Transcription Start Site (TSS) in base pairs for gene-TF association. Defines promoter region. Default 5000 (±5kb, captures typical promoters). Use 1000-2000 for narrow promoters, 10000+ for enhancer regions.",
          "default": "5000"
        },
        "gene": {
          "type": [
            "array",
            "string"
          ],
          "items": {
            "type": "string"
          },
          "description": "DNAKit optional alias for gene_list"
        },
        "cell_type": {
          "type": "string",
          "description": "DNAKit optional alias for cell_type_class"
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "description": "Response object with enrichment analysis URL and parameters",
      "additionalProperties": true,
      "properties": {
        "message": {
          "type": "string"
        },
        "genes": {
          "type": "array",
          "items": {
            "type": "string"
          }
        },
        "distance_from_tss": {
          "type": "string"
        },
        "url": {
          "type": "string"
        },
        "note": {
          "type": "string"
        }
      }
    }
    ```

??? info "DNAKit 参数扩展规则"

    ```json
    {
      "aliases": {
        "gene": "gene_list",
        "cell_type": "cell_type_class"
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `ChIPAtlas_enrichment_analysis`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ch_ipatlas_enrichment_analysis", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ch_ipatlas_get_experiments` {#dnakit.databases.ch_ipatlas_get_experiments}

Search ChIP-Atlas experiment metadata including experiment IDs, antigens, cell types, and processing statistics. Filter by genome, antigen (TF/histone), or cell type. Returns experiment details with SRX/ERX/DRX IDs. Use for: finding relevant experiments, getting experiment metadata, dataset discovery, batch analysis planning.

```text
dnakit.databases.ch_ipatlas_get_experiments(*, operation=None, genome=None, antigen=None, cell_type=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | 提供方操作模式；默认值与允许值见本行约束。 {&quot;enum&quot;: &#91;&quot;get_experiment_list&quot;&#93;, &quot;default&quot;: &quot;get_experiment_list&quot;} |
| `genome` | string | None | Filter by genome assembly {&quot;enum&quot;: &#91;&quot;hg38&quot;, &quot;hg19&quot;, &quot;mm10&quot;, &quot;mm9&quot;, &quot;rn6&quot;, &quot;dm6&quot;, &quot;dm3&quot;, &quot;ce11&quot;, &quot;ce10&quot;, &quot;sacCer3&quot;&#93;} |
| `antigen` | string | None | Filter by antigen/protein name (e.g., &#x27;CTCF&#x27;, &#x27;H3K4me3&#x27;) |
| `cell_type` | string | None | Filter by cell type (e.g., &#x27;K-562&#x27;, &#x27;HepG2&#x27;, &#x27;CD4&#x27;) |
| `limit` | integer | None | Maximum results to return {&quot;default&quot;: 100, &quot;minimum&quot;: 1, &quot;maximum&quot;: 1000} |
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
            "get_experiment_list"
          ],
          "default": "get_experiment_list"
        },
        "genome": {
          "type": "string",
          "description": "Filter by genome assembly",
          "enum": [
            "hg38",
            "hg19",
            "mm10",
            "mm9",
            "rn6",
            "dm6",
            "dm3",
            "ce11",
            "ce10",
            "sacCer3"
          ]
        },
        "antigen": {
          "type": "string",
          "description": "Filter by antigen/protein name (e.g., 'CTCF', 'H3K4me3')"
        },
        "cell_type": {
          "type": "string",
          "description": "Filter by cell type (e.g., 'K-562', 'HepG2', 'CD4')"
        },
        "limit": {
          "type": "integer",
          "description": "Maximum results to return",
          "default": 100,
          "minimum": 1,
          "maximum": 1000
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "description": "Response with experiment metadata including URLs and file info",
      "additionalProperties": true,
      "properties": {
        "message": {
          "type": "string"
        },
        "metadata_file_url": {
          "type": "string"
        },
        "web_search_url": {
          "type": "string"
        },
        "file_size": {
          "type": "string"
        },
        "total_experiments": {
          "type": "string"
        },
        "filters_requested": {
          "type": "object",
          "additionalProperties": true
        },
        "recommendation": {
          "type": "string"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `ChIPAtlas_get_experiments`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ch_ipatlas_get_experiments", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ch_ipatlas_get_peak_data` {#dnakit.databases.ch_ipatlas_get_peak_data}

查询 ChIP-Atlas 实验峰数据入口；gene 模式先定位 GRCh38 基因区域，提供 peak_file 后计算 BED 重叠，否则返回数据获取计划。

```text
dnakit.databases.ch_ipatlas_get_peak_data(*, operation=None, experiment_id=None, genome=None, format=None, threshold=None, gene=None, experiment_type=None, limit=None, peak_file=None, flank_bp=None, max_file_bytes=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | 提供方操作模式；默认值与允许值见本行约束。 {&quot;enum&quot;: &#91;&quot;get_peak_data&quot;&#93;, &quot;default&quot;: &quot;get_peak_data&quot;} |
| `experiment_id` | string | None | Experiment ID (SRX/ERX/DRX format, required) |
| `genome` | string | None | Genome assembly {&quot;enum&quot;: &#91;&quot;hg38&quot;, &quot;hg19&quot;, &quot;mm10&quot;, &quot;mm9&quot;, &quot;rn6&quot;, &quot;dm6&quot;, &quot;dm3&quot;, &quot;ce11&quot;, &quot;ce10&quot;, &quot;sacCer3&quot;&#93;, &quot;default&quot;: &quot;hg38&quot;} |
| `format` | string | None | Output format {&quot;enum&quot;: &#91;&quot;bigwig&quot;, &quot;bed&quot;, &quot;bigbed&quot;&#93;, &quot;default&quot;: &quot;bigwig&quot;} |
| `threshold` | string | None | Q-value threshold for BED/BigBed peak files. &#x27;05&#x27;=1e-5 (more peaks), &#x27;10&#x27;=1e-10 (moderate), &#x27;20&#x27;=1e-20 (high confidence only). Only applies to BED/BigBed formats. Default &#x27;05&#x27;. {&quot;enum&quot;: &#91;&quot;05&quot;, &quot;10&quot;, &quot;20&quot;&#93;, &quot;default&quot;: &quot;05&quot;} |
| `gene` | string | None | 基因标识或基因列表；命名要求和类型见当前接口。 |
| `experiment_type` | string | None | 实验类别，按数据源接受的类型筛选。 |
| `limit` | integer | None | 当前请求的返回数量上限；不是全库命中总数。 {&quot;minimum&quot;: 1, &quot;maximum&quot;: 10000} |
| `peak_file` | string | None | 本地 BED 峰文件；必须核对其基因组版本和轨道来源。 |
| `flank_bp` | integer | None | 基因区域两侧扩展的碱基数。 {&quot;minimum&quot;: 0, &quot;maximum&quot;: 1000000} |
| `max_file_bytes` | integer | None | 允许读取的峰文件字节上限。 {&quot;minimum&quot;: 1} |
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
            "get_peak_data"
          ],
          "default": "get_peak_data"
        },
        "experiment_id": {
          "type": "string",
          "description": "Experiment ID (SRX/ERX/DRX format, required)"
        },
        "genome": {
          "type": "string",
          "description": "Genome assembly",
          "enum": [
            "hg38",
            "hg19",
            "mm10",
            "mm9",
            "rn6",
            "dm6",
            "dm3",
            "ce11",
            "ce10",
            "sacCer3"
          ],
          "default": "hg38"
        },
        "format": {
          "type": "string",
          "description": "Output format",
          "enum": [
            "bigwig",
            "bed",
            "bigbed"
          ],
          "default": "bigwig"
        },
        "threshold": {
          "type": "string",
          "description": "Q-value threshold for BED/BigBed peak files. '05'=1e-5 (more peaks), '10'=1e-10 (moderate), '20'=1e-20 (high confidence only). Only applies to BED/BigBed formats. Default '05'.",
          "enum": [
            "05",
            "10",
            "20"
          ],
          "default": "05"
        },
        "gene": {
          "type": "string"
        },
        "experiment_type": {
          "type": "string"
        },
        "limit": {
          "type": "integer",
          "minimum": 1,
          "maximum": 10000
        },
        "peak_file": {
          "type": "string"
        },
        "flank_bp": {
          "type": "integer",
          "minimum": 0,
          "maximum": 1000000
        },
        "max_file_bytes": {
          "type": "integer",
          "minimum": 1
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "description": "Response with download URL for peak data",
      "additionalProperties": true,
      "properties": {
        "experiment_id": {
          "type": "string"
        },
        "genome": {
          "type": "string"
        },
        "format": {
          "type": "string"
        },
        "url": {
          "type": "string"
        },
        "message": {
          "type": "string"
        }
      }
    }
    ```

??? info "DNAKit 参数扩展规则"

    ```json
    {
      "parameters": {
        "gene": {
          "type": "string"
        },
        "experiment_type": {
          "type": "string"
        },
        "limit": {
          "type": "integer",
          "minimum": 1,
          "maximum": 10000
        },
        "peak_file": {
          "type": "string"
        },
        "flank_bp": {
          "type": "integer",
          "minimum": 0,
          "maximum": 1000000
        },
        "max_file_bytes": {
          "type": "integer",
          "minimum": 1
        }
      },
      "optional": [
        "experiment_id"
      ]
    }
    ```

**执行与来源：** `network`；上游接口名 `ChIPAtlas_get_peak_data`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ch_ipatlas_get_peak_data", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ch_ipatlas_search_datasets` {#dnakit.databases.ch_ipatlas_search_datasets}

Search ChIP-Atlas by antigen or cell type to find available datasets. Returns number of experiments and experiment IDs. Use antigenList for finding TF/histone datasets across cell types. Use celltypeList for finding all factors in specific cell type. Use for: dataset discovery, planning comparative analysis, checking data availability.

```text
dnakit.databases.ch_ipatlas_search_datasets(*, operation=None, antigen=None, cell_type=None, genome=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | 提供方操作模式；默认值与允许值见本行约束。 {&quot;enum&quot;: &#91;&quot;search_datasets&quot;&#93;, &quot;default&quot;: &quot;search_datasets&quot;} |
| `antigen` | string | None | Search by antigen/protein name (e.g., &#x27;CTCF&#x27;, &#x27;H3K27ac&#x27;) |
| `cell_type` | string | None | Search by cell type (e.g., &#x27;K-562&#x27;, &#x27;HeLa&#x27;, &#x27;MCF-7&#x27;) |
| `genome` | string | None | Genome assembly {&quot;enum&quot;: &#91;&quot;hg38&quot;, &quot;hg19&quot;, &quot;mm10&quot;, &quot;mm9&quot;, &quot;rn6&quot;, &quot;dm6&quot;, &quot;dm3&quot;, &quot;ce11&quot;, &quot;ce10&quot;, &quot;sacCer3&quot;&#93;, &quot;default&quot;: &quot;hg38&quot;} |
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
            "search_datasets"
          ],
          "default": "search_datasets"
        },
        "antigen": {
          "type": "string",
          "description": "Search by antigen/protein name (e.g., 'CTCF', 'H3K27ac')"
        },
        "cell_type": {
          "type": "string",
          "description": "Search by cell type (e.g., 'K-562', 'HeLa', 'MCF-7')"
        },
        "genome": {
          "type": "string",
          "description": "Genome assembly",
          "enum": [
            "hg38",
            "hg19",
            "mm10",
            "mm9",
            "rn6",
            "dm6",
            "dm3",
            "ce11",
            "ce10",
            "sacCer3"
          ],
          "default": "hg38"
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "description": "Response with search results for datasets",
      "additionalProperties": true,
      "properties": {
        "search_key": {
          "type": "string"
        },
        "search_value": {
          "type": "string"
        },
        "genome": {
          "type": "string"
        },
        "num_results": {
          "type": "integer"
        },
        "results": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "genome": {
                "type": "string"
              },
              "class": {
                "type": "string"
              },
              "name": {
                "type": "string"
              },
              "num_experiments": {
                "type": "string"
              },
              "experiment_ids": {
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

**执行与来源：** `network`；上游接口名 `ChIPAtlas_search_datasets`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ch_ipatlas_search_datasets", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.dep_map_get_gene_dependencies` {#dnakit.databases.dep_map_get_gene_dependencies}

Get gene dependency information from CRISPR knockout screens. Returns gene IDs and annotations. Negative gene effect scores indicate essentiality (cell death upon knockout).

```text
dnakit.databases.dep_map_get_gene_dependencies(*, gene_symbol, model_id=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 必填（输入 Schema） | Gene symbol (e.g., &#x27;EGFR&#x27;, &#x27;KRAS&#x27;, &#x27;TP53&#x27;) |
| `model_id` | string | None | Optional: Filter by specific cell line |
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
          "description": "Gene symbol (e.g., 'EGFR', 'KRAS', 'TP53')"
        },
        "model_id": {
          "type": "string",
          "description": "Optional: Filter by specific cell line"
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
            "genes": {
              "type": "array"
            }
          }
        },
        "error": {
          "type": "string"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `DepMap_get_gene_dependencies`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("dep_map_get_gene_dependencies", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.dep_map_search_cell_lines` {#dnakit.databases.dep_map_search_cell_lines}

Search cancer cell lines by name. Returns matching cell lines with IDs and cancer types. Use for finding cell models by common names (A549, HeLa, MCF7, etc.).

```text
dnakit.databases.dep_map_search_cell_lines(*, query, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填（输入 Schema） | Search query (cell line name like &#x27;A549&#x27;, &#x27;MCF7&#x27;) |
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
          "description": "Search query (cell line name like 'A549', 'MCF7')"
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
            "query": {
              "type": "string"
            },
            "cell_lines": {
              "type": "array"
            },
            "count": {
              "type": "integer"
            }
          }
        },
        "error": {
          "type": "string"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `DepMap_search_cell_lines`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("dep_map_search_cell_lines", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.e_qtl_get_associations` {#dnakit.databases.e_qtl_get_associations}

Get eQTL associations (variant-gene expression correlations) from a specific dataset in the EBI eQTL Catalogue. Returns statistical associations between genetic variants and molecular traits (gene/exon/transcript expression or splicing). Each result includes beta effect size, p-value, allele frequency, and variant coordinates. Use eQTL_list_datasets first to find available dataset IDs. Can filter by gene_id (Ensembl gene ID) or variant (chr_pos_ref_alt format).

```text
dnakit.databases.e_qtl_get_associations(*, dataset_id, size=None, gene_id=None, variant=None, molecular_trait_id=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `dataset_id` | string | 必填（输入 Schema） | Dataset identifier from eQTL_list_datasets (e.g., &#x27;QTD000001&#x27;, &#x27;QTD000584&#x27;) |
| `size` | integer | None | Number of associations to return (default 20) {&quot;default&quot;: 20} |
| `gene_id` | string &#124; null | None | Filter by Ensembl gene ID (e.g., &#x27;ENSG00000141510&#x27; for TP53). Returns all eQTLs affecting this gene. |
| `variant` | string &#124; null | None | Filter by variant in chr_pos_ref_alt format (e.g., &#x27;chr1_791100_G_GGGA&#x27;). Returns all genes affected by this variant. |
| `molecular_trait_id` | string &#124; null | None | Filter by molecular trait ID (usually same as gene_id for &#x27;ge&#x27; datasets) |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "dataset_id": {
          "type": "string",
          "description": "Dataset identifier from eQTL_list_datasets (e.g., 'QTD000001', 'QTD000584')"
        },
        "size": {
          "type": "integer",
          "default": 20,
          "description": "Number of associations to return (default 20)"
        },
        "gene_id": {
          "type": [
            "string",
            "null"
          ],
          "description": "Filter by Ensembl gene ID (e.g., 'ENSG00000141510' for TP53). Returns all eQTLs affecting this gene."
        },
        "variant": {
          "type": [
            "string",
            "null"
          ],
          "description": "Filter by variant in chr_pos_ref_alt format (e.g., 'chr1_791100_G_GGGA'). Returns all genes affected by this variant."
        },
        "molecular_trait_id": {
          "type": [
            "string",
            "null"
          ],
          "description": "Filter by molecular trait ID (usually same as gene_id for 'ge' datasets)"
        }
      },
      "required": [
        "dataset_id"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "array",
          "description": "Array of eQTL association objects",
          "items": {
            "type": "object",
            "properties": {
              "molecular_trait_id": {
                "type": "string",
                "description": "Molecular trait identifier"
              },
              "gene_id": {
                "type": "string",
                "description": "Ensembl gene ID"
              },
              "variant": {
                "type": "string",
                "description": "Variant ID in chr_pos_ref_alt format"
              },
              "rsid": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "dbSNP rsID if available"
              },
              "chromosome": {
                "type": "string",
                "description": "Chromosome"
              },
              "position": {
                "type": "integer",
                "description": "Genomic position"
              },
              "ref": {
                "type": "string"
              },
              "alt": {
                "type": "string"
              },
              "type": {
                "type": "string",
                "description": "Variant type (SNP, INDEL)"
              },
              "beta": {
                "type": "number",
                "description": "Effect size (beta coefficient)"
              },
              "se": {
                "type": "number",
                "description": "Standard error"
              },
              "pvalue": {
                "type": "number",
                "description": "P-value for the association"
              },
              "nlog10p": {
                "type": "number",
                "description": "Negative log10 of p-value"
              },
              "maf": {
                "type": "number",
                "description": "Minor allele frequency"
              },
              "ac": {
                "type": "integer",
                "description": "Allele count"
              },
              "an": {
                "type": "integer",
                "description": "Total allele number"
              },
              "median_tpm": {
                "type": [
                  "number",
                  "null"
                ],
                "description": "Median TPM expression of the gene"
              },
              "r2": {
                "type": "number",
                "description": "Imputation quality (R-squared)"
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

**执行与来源：** `network`；上游接口名 `eQTL_get_associations`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("e_qtl_get_associations", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.e_qtl_list_datasets` {#dnakit.databases.e_qtl_list_datasets}

List available eQTL (expression Quantitative Trait Loci) datasets from the EBI eQTL Catalogue. Each dataset represents a study-tissue-quantification combination with measured genetic associations to gene expression. Returns dataset IDs, study labels, tissue types, sample sizes, and quantification methods (ge=gene expression, exon, tx=transcript, txrev=transcript ratio, leafcutter=splicing). Use dataset IDs with eQTL_get_associations for specific variant-gene associations.

```text
dnakit.databases.e_qtl_list_datasets(*, size=None, quant_method=None, tissue_id=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `size` | integer | None | Number of datasets to return (default 20) {&quot;default&quot;: 20} |
| `quant_method` | string &#124; null | None | Filter by quantification method: &#x27;ge&#x27; (gene expression), &#x27;exon&#x27;, &#x27;tx&#x27; (transcript), &#x27;txrev&#x27; (transcript ratio), &#x27;leafcutter&#x27; (splicing) |
| `tissue_id` | string &#124; null | None | Filter by tissue ontology ID (e.g., &#x27;CL_0000235&#x27; for macrophage, &#x27;UBERON_0002107&#x27; for liver) |
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
          "default": 20,
          "description": "Number of datasets to return (default 20)"
        },
        "quant_method": {
          "type": [
            "string",
            "null"
          ],
          "description": "Filter by quantification method: 'ge' (gene expression), 'exon', 'tx' (transcript), 'txrev' (transcript ratio), 'leafcutter' (splicing)"
        },
        "tissue_id": {
          "type": [
            "string",
            "null"
          ],
          "description": "Filter by tissue ontology ID (e.g., 'CL_0000235' for macrophage, 'UBERON_0002107' for liver)"
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
          "description": "Array of eQTL dataset objects",
          "items": {
            "type": "object",
            "properties": {
              "dataset_id": {
                "type": "string",
                "description": "Unique dataset identifier (e.g., 'QTD000001')"
              },
              "study_id": {
                "type": "string",
                "description": "Study identifier (e.g., 'QTS000001')"
              },
              "study_label": {
                "type": "string",
                "description": "Human-readable study name (e.g., 'Alasoo_2018')"
              },
              "quant_method": {
                "type": "string",
                "description": "Quantification method (ge, exon, tx, txrev, leafcutter)"
              },
              "tissue_id": {
                "type": "string",
                "description": "Tissue ontology ID"
              },
              "tissue_label": {
                "type": "string",
                "description": "Human-readable tissue name"
              },
              "condition_label": {
                "type": "string",
                "description": "Experimental condition (e.g., 'naive', 'LPS')"
              },
              "sample_size": {
                "type": "integer",
                "description": "Number of samples in the dataset"
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

**执行与来源：** `network`；上游接口名 `eQTL_list_datasets`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("e_qtl_list_datasets", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.encode_get_biosample` {#dnakit.databases.encode_get_biosample}

Get detailed metadata for a specific ENCODE biosample by its accession ID. Returns comprehensive information including biosample type, organism, life stage, age, sex, tissue/cell line name, donor information, treatments, genetic modifications, culture conditions, and characterizations. Essential for understanding sample provenance, experimental conditions, and biological context of ENCODE datasets. Find biosample accessions using ENCODE_search_biosamples (format: ENCBS######).

```text
dnakit.databases.encode_get_biosample(*, accession, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string | 必填（输入 Schema） | ENCODE biosample accession identifier (format: ENCBS######, e.g., &#x27;ENCBS000AAA&#x27;, &#x27;ENCBS130ENC&#x27;). Find accessions using ENCODE_search_biosamples or from experiment metadata. |
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
          "description": "ENCODE biosample accession identifier (format: ENCBS######, e.g., 'ENCBS000AAA', 'ENCBS130ENC'). Find accessions using ENCODE_search_biosamples or from experiment metadata."
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
              "type": "string"
            },
            "url": {
              "type": "string"
            },
            "aliases": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "schema_version": {
              "type": "string"
            },
            "status": {
              "type": "string"
            },
            "lab": {
              "type": "string"
            },
            "award": {
              "type": "string"
            },
            "date_created": {
              "type": "string"
            },
            "submitted_by": {
              "type": "string"
            },
            "notes": {
              "type": "string"
            },
            "documents": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "references": {
              "type": "array"
            },
            "source": {
              "type": "string"
            },
            "product_id": {
              "type": "string"
            },
            "biosample_ontology": {
              "type": "string"
            },
            "genetic_modifications": {
              "type": "array"
            },
            "alternate_accessions": {
              "type": "array"
            },
            "description": {
              "type": "string"
            },
            "treatments": {
              "type": "array"
            },
            "dbxrefs": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "donor": {
              "type": "string"
            },
            "organism": {
              "type": "string"
            },
            "passage_number": {
              "type": [
                "integer",
                "number"
              ]
            },
            "internal_tags": {
              "type": "array"
            },
            "culture_harvest_date": {
              "type": "string"
            },
            "culture_start_date": {
              "type": "string"
            },
            "@id": {
              "type": "string"
            },
            "@type": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "uuid": {
              "type": "string"
            },
            "sex": {
              "type": "string"
            },
            "age": {
              "type": "string"
            },
            "age_units": {
              "type": "string"
            },
            "health_status": {
              "type": "string"
            },
            "life_stage": {
              "type": "string"
            },
            "applied_modifications": {
              "type": "array"
            },
            "characterizations": {
              "type": "array"
            },
            "parent_of": {
              "type": "array"
            },
            "age_display": {
              "type": "string"
            },
            "origin_batch": {
              "type": "string"
            },
            "perturbed": {
              "type": "boolean"
            },
            "simple_summary": {
              "type": "string"
            },
            "summary": {
              "type": "string"
            },
            "model_organism_sex": {
              "type": "string"
            },
            "mouse_life_stage": {
              "type": "string"
            },
            "model_organism_age": {
              "type": "string"
            },
            "model_organism_age_units": {
              "type": "string"
            },
            "model_organism_health_status": {
              "type": "string"
            },
            "model_organism_donor_modifications": {
              "type": "array"
            }
          },
          "required": [
            "accession"
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

**执行与来源：** `network`；上游接口名 `ENCODE_get_biosample`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("encode_get_biosample", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.encode_get_chromatin_state` {#dnakit.databases.encode_get_chromatin_state}

Search ENCODE chromatin state annotations (ChromHMM segmentations) for specific biosamples or tissues. Chromatin states are genome-wide annotations that classify each genomic position into a functional state (e.g., active TSS, flanking active TSS, strong transcription, weak transcription, genic enhancers, active enhancers, bivalent enhancer, bivalent/poised TSS, flanking bivalent TSS, repressed polycomb, weak repressed polycomb, quiescent, heterochromatin) based on combinatorial patterns of histone marks. Use this to understand the chromatin landscape of specific cell types.

```text
dnakit.databases.encode_get_chromatin_state(*, biosample_term_name=None, organism=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `biosample_term_name` | string &#124; null | None | Biosample name (e.g., &#x27;K562&#x27;, &#x27;HepG2&#x27;, &#x27;GM12878&#x27;, &#x27;liver&#x27;). Leave empty for all. |
| `organism` | string | None | Organism scientific name. {&quot;default&quot;: &quot;Homo sapiens&quot;} |
| `limit` | integer | None | Maximum number of results. {&quot;default&quot;: 25, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "biosample_term_name": {
          "type": [
            "string",
            "null"
          ],
          "description": "Biosample name (e.g., 'K562', 'HepG2', 'GM12878', 'liver'). Leave empty for all."
        },
        "organism": {
          "type": "string",
          "default": "Homo sapiens",
          "description": "Organism scientific name."
        },
        "limit": {
          "type": "integer",
          "default": 25,
          "minimum": 1,
          "maximum": 100,
          "description": "Maximum number of results."
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
              "type": "integer"
            },
            "annotations": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "accession": {
                    "type": "string"
                  },
                  "annotation_type": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "description": {
                    "type": "string"
                  },
                  "biosample_summary": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "status": {
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

**执行与来源：** `network`；上游接口名 `ENCODE_get_chromatin_state`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("encode_get_chromatin_state", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.encode_get_experiment` {#dnakit.databases.encode_get_experiment}

Get comprehensive metadata for a specific ENCODE experiment by its accession ID. Returns detailed information including assay type, target protein, biosample details, replicates, files, quality metrics, analysis pipeline, lab, award, and publication references. Use this to understand experimental design, data quality, and available files before downloading. Find experiment accessions using ENCODE_search_experiments (format: ENCSR######).

```text
dnakit.databases.encode_get_experiment(*, accession, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string | 必填（输入 Schema） | ENCODE experiment accession identifier (format: ENCSR######, e.g., &#x27;ENCSR000AKS&#x27;, &#x27;ENCSR000CAG&#x27;). Find accessions using ENCODE_search_experiments. |
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
          "description": "ENCODE experiment accession identifier (format: ENCSR######, e.g., 'ENCSR000AKS', 'ENCSR000CAG'). Find accessions using ENCODE_search_experiments."
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
            "assay_term_name": {
              "type": "string"
            },
            "biosample_ontology": {
              "type": "string"
            },
            "documents": {
              "type": "array"
            },
            "references": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "schema_version": {
              "type": "string"
            },
            "accession": {
              "type": "string"
            },
            "alternate_accessions": {
              "type": "array"
            },
            "analyses": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "description": {
              "type": "string"
            },
            "dbxrefs": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "date_released": {
              "type": "string"
            },
            "doi": {
              "type": "string"
            },
            "internal_tags": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "status": {
              "type": "string"
            },
            "date_created": {
              "type": "string"
            },
            "submitted_by": {
              "type": "string"
            },
            "lab": {
              "type": "string"
            },
            "award": {
              "type": "string"
            },
            "aliases": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "date_submitted": {
              "type": "string"
            },
            "target": {
              "type": "string"
            },
            "possible_controls": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "supersedes": {
              "type": "array"
            },
            "related_files": {
              "type": "array"
            },
            "internal_status": {
              "type": "string"
            },
            "bio_replicate_count": {
              "type": [
                "integer",
                "number"
              ]
            },
            "tech_replicate_count": {
              "type": [
                "integer",
                "number"
              ]
            },
            "replication_type": {
              "type": "string"
            },
            "objective_slims": {
              "type": "array"
            },
            "type_slims": {
              "type": "array"
            },
            "category_slims": {
              "type": "array"
            },
            "assay_title": {
              "type": "string"
            },
            "assay_slims": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "replicates": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "simple_biosample_summary": {
              "type": "string"
            },
            "biosample_summary": {
              "type": "string"
            },
            "assay_term_id": {
              "type": "string"
            },
            "@id": {
              "type": "string"
            },
            "@type": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "uuid": {
              "type": "string"
            },
            "original_files": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "contributing_files": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "files": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "revoked_files": {
              "type": "array"
            },
            "assembly": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "hub": {
              "type": "string"
            },
            "default_analysis": {
              "type": "string"
            },
            "related_series": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "superseded_by": {
              "type": "array"
            },
            "related_annotations": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "life_stage_age": {
              "type": "string"
            },
            "perturbed": {
              "type": "boolean"
            }
          },
          "required": [
            "assay_term_name"
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

**执行与来源：** `network`；上游接口名 `ENCODE_get_experiment`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("encode_get_experiment", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.encode_list_files` {#dnakit.databases.encode_list_files}

List ENCODE data files with filters by file format, output type, or assay. Returns file metadata including accessions, formats (FASTQ, BAM, bigWig), output types (alignments, peaks, signal), and file sizes. Use this to programmatically discover downloadable artifacts for specific experiments or data types. Find file accessions to use with ENCODE_get_file for detailed file information and download URLs.

```text
dnakit.databases.encode_list_files(*, file_type=None, assay_title=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `file_type` | string | None | File format filter (e.g., &#x27;fastq&#x27;, &#x27;bam&#x27;, &#x27;bigWig&#x27;, &#x27;bed&#x27;). Common formats: fastq (raw reads), bam (alignments), bigWig (coverage), bed (peaks/regions). |
| `assay_title` | string | None | Assay filter (e.g., &#x27;ChIP-seq&#x27;, &#x27;RNA-seq&#x27;). Filters files by the experimental assay that generated them. |
| `limit` | integer | None | Maximum number of results to return (1–100). {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "file_type": {
          "type": "string",
          "description": "File format filter (e.g., 'fastq', 'bam', 'bigWig', 'bed'). Common formats: fastq (raw reads), bam (alignments), bigWig (coverage), bed (peaks/regions)."
        },
        "assay_title": {
          "type": "string",
          "description": "Assay filter (e.g., 'ChIP-seq', 'RNA-seq'). Filters files by the experimental assay that generated them."
        },
        "limit": {
          "type": "integer",
          "default": 10,
          "minimum": 1,
          "maximum": 100,
          "description": "Maximum number of results to return (1–100)."
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
            "@graph": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "@id": {
                    "type": "string"
                  },
                  "@type": {
                    "type": "array",
                    "items": {
                      "type": "string"
                    }
                  },
                  "accession": {
                    "type": "string"
                  },
                  "assay_term_name": {
                    "type": "string"
                  },
                  "award": {
                    "type": "object",
                    "properties": {
                      "project": {
                        "type": "string"
                      }
                    }
                  },
                  "biological_replicates": {
                    "type": "array",
                    "items": {
                      "type": [
                        "integer",
                        "number"
                      ]
                    }
                  },
                  "biosample_ontology": {
                    "type": "object",
                    "properties": {
                      "organ_slims": {
                        "type": "array",
                        "items": {
                          "type": "string"
                        }
                      },
                      "term_name": {
                        "type": "string"
                      }
                    }
                  },
                  "dataset": {
                    "type": "string"
                  },
                  "date_created": {
                    "type": "string"
                  },
                  "file_format": {
                    "type": "string"
                  },
                  "file_size": {
                    "type": [
                      "integer",
                      "number"
                    ]
                  },
                  "file_type": {
                    "type": "string"
                  },
                  "href": {
                    "type": "string"
                  },
                  "lab": {
                    "type": "object",
                    "properties": {
                      "title": {
                        "type": "string"
                      }
                    }
                  },
                  "output_category": {
                    "type": "string"
                  },
                  "output_type": {
                    "type": "string"
                  },
                  "quality_metrics": {
                    "type": "array"
                  },
                  "replicate": {
                    "type": "object",
                    "properties": {
                      "library": {
                        "type": "object",
                        "properties": {
                          "accession": {
                            "type": "string"
                          }
                        }
                      }
                    }
                  },
                  "run_type": {
                    "type": "string"
                  },
                  "simple_biosample_summary": {
                    "type": "string"
                  },
                  "status": {
                    "type": "string"
                  },
                  "technical_replicates": {
                    "type": "array",
                    "items": {
                      "type": "string"
                    }
                  },
                  "title": {
                    "type": "string"
                  }
                }
              }
            },
            "@type": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "all": {
              "type": "string"
            },
            "notification": {
              "type": "string"
            },
            "title": {
              "type": "string"
            },
            "total": {
              "type": [
                "integer",
                "number"
              ]
            }
          },
          "required": [
            "@graph"
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

**执行与来源：** `network`；上游接口名 `ENCODE_list_files`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("encode_list_files", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.encode_search_annotations` {#dnakit.databases.encode_search_annotations}

Search ENCODE annotations including candidate cis-Regulatory Elements (cCREs), chromatin states, and imputed signals. cCREs are classified as promoter-like (PLS), proximal enhancer-like (pELS), distal enhancer-like (dELS), DNase-H3K4me3 (DNase-H3K4me3), and CTCF-only. Chromatin state annotations use ChromHMM models to segment the genome into functional states. Use this to find regulatory element annotations for specific biosamples or genome assemblies.

```text
dnakit.databases.encode_search_annotations(*, annotation_type=None, biosample_term_name=None, organism=None, assembly=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `annotation_type` | string | None | Annotation type filter. Options: &#x27;candidate Cis-Regulatory Elements&#x27; (cCREs), &#x27;chromatin state&#x27; (ChromHMM), &#x27;imputed signal&#x27; (predicted epigenetic marks). Default: cCREs. {&quot;default&quot;: &quot;candidate Cis-Regulatory Elements&quot;} |
| `biosample_term_name` | string &#124; null | None | Biosample filter (e.g., &#x27;K562&#x27;, &#x27;HepG2&#x27;). Leave empty for all. |
| `organism` | string | None | Organism scientific name. {&quot;default&quot;: &quot;Homo sapiens&quot;} |
| `assembly` | string | None | Genome assembly (e.g., &#x27;GRCh38&#x27;, &#x27;hg19&#x27;, &#x27;mm10&#x27;). {&quot;default&quot;: &quot;GRCh38&quot;} |
| `limit` | integer | None | Maximum number of results. {&quot;default&quot;: 25, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "annotation_type": {
          "type": "string",
          "default": "candidate Cis-Regulatory Elements",
          "description": "Annotation type filter. Options: 'candidate Cis-Regulatory Elements' (cCREs), 'chromatin state' (ChromHMM), 'imputed signal' (predicted epigenetic marks). Default: cCREs."
        },
        "biosample_term_name": {
          "type": [
            "string",
            "null"
          ],
          "description": "Biosample filter (e.g., 'K562', 'HepG2'). Leave empty for all."
        },
        "organism": {
          "type": "string",
          "default": "Homo sapiens",
          "description": "Organism scientific name."
        },
        "assembly": {
          "type": "string",
          "default": "GRCh38",
          "description": "Genome assembly (e.g., 'GRCh38', 'hg19', 'mm10')."
        },
        "limit": {
          "type": "integer",
          "default": 25,
          "minimum": 1,
          "maximum": 100,
          "description": "Maximum number of results."
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
              "type": "integer"
            },
            "annotations": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "accession": {
                    "type": "string"
                  },
                  "annotation_type": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "description": {
                    "type": "string"
                  },
                  "biosample_summary": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "status": {
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

**执行与来源：** `network`；上游接口名 `ENCODE_search_annotations`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("encode_search_annotations", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.encode_search_biosamples` {#dnakit.databases.encode_search_biosamples}

Search ENCODE biosamples (cell lines, tissues, primary cells) by organism, biosample type, or treatment. Returns biosample accessions, types, organism, life stage, and descriptions. Use this to discover biological samples used in ENCODE experiments, find cell lines or tissues of interest, or identify treated vs. control samples. Find biosample accessions to use with ENCODE_get_biosample for detailed characterization data.

```text
dnakit.databases.encode_search_biosamples(*, organism=None, biosample_type=None, treatment=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `organism` | string | None | Organism filter (e.g., &#x27;Homo sapiens&#x27;, &#x27;Mus musculus&#x27;). |
| `biosample_type` | string | None | Biosample classification filter — maps to ENCODE&#x27;s biosample_ontology.classification field (e.g., &#x27;cell line&#x27;, &#x27;tissue&#x27;, &#x27;primary cell&#x27;, &#x27;in vitro differentiated cells&#x27;). |
| `treatment` | string | None | Treatment filter (e.g., &#x27;interferon gamma&#x27;, &#x27;ethanol&#x27;). Use &#x27;*&#x27; to find any treated samples. |
| `limit` | integer | None | Maximum number of results to return (1–100). {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "organism": {
          "type": "string",
          "description": "Organism filter (e.g., 'Homo sapiens', 'Mus musculus')."
        },
        "biosample_type": {
          "type": "string",
          "description": "Biosample classification filter — maps to ENCODE's biosample_ontology.classification field (e.g., 'cell line', 'tissue', 'primary cell', 'in vitro differentiated cells')."
        },
        "treatment": {
          "type": "string",
          "description": "Treatment filter (e.g., 'interferon gamma', 'ethanol'). Use '*' to find any treated samples."
        },
        "limit": {
          "type": "integer",
          "default": 10,
          "minimum": 1,
          "maximum": 100,
          "description": "Maximum number of results to return (1–100)."
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
            "@graph": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "@id": {
                    "type": "string"
                  },
                  "@type": {
                    "type": "array",
                    "items": {
                      "type": "string"
                    }
                  },
                  "accession": {
                    "type": "string"
                  },
                  "age": {
                    "type": "string"
                  },
                  "applied_modifications": {
                    "type": "array",
                    "items": {
                      "type": "object",
                      "properties": {
                        "method": {
                          "type": "string"
                        },
                        "purpose": {
                          "type": "string"
                        },
                        "modified_site_by_target_id": {
                          "type": "object",
                          "properties": {
                            "organism": {
                              "type": "object",
                              "properties": {
                                "schema_version": {
                                  "type": "string"
                                },
                                "@type": {
                                  "type": "array",
                                  "items": {
                                    "type": "string"
                                  }
                                },
                                "name": {
                                  "type": "string"
                                },
                                "taxon_id": {
                                  "type": "string"
                                },
                                "scientific_name": {
                                  "type": "string"
                                },
                                "@id": {
                                  "type": "string"
                                },
                                "uuid": {
                                  "type": "string"
                                },
                                "status": {
                                  "type": "string"
                                }
                              }
                            },
                            "name": {
                              "type": "string"
                            }
                          }
                        },
                        "category": {
                          "type": "string"
                        }
                      }
                    }
                  },
                  "biosample_ontology": {
                    "type": "object",
                    "properties": {
                      "term_name": {
                        "type": "string"
                      },
                      "classification": {
                        "type": "string"
                      }
                    }
                  },
                  "lab": {
                    "type": "object",
                    "properties": {
                      "title": {
                        "type": "string"
                      }
                    }
                  },
                  "life_stage": {
                    "type": "string"
                  },
                  "organism": {
                    "type": "object",
                    "properties": {
                      "scientific_name": {
                        "type": "string"
                      }
                    }
                  },
                  "parent_of": {
                    "type": "array"
                  },
                  "source": {
                    "type": "object",
                    "properties": {
                      "title": {
                        "type": "string"
                      }
                    }
                  },
                  "status": {
                    "type": "string"
                  },
                  "summary": {
                    "type": "string"
                  },
                  "treatments": {
                    "type": "array",
                    "items": {
                      "type": "object",
                      "properties": {
                        "treatment_term_name": {
                          "type": "string"
                        },
                        "purpose": {
                          "type": "string"
                        }
                      }
                    }
                  }
                }
              }
            },
            "@type": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "notification": {
              "type": "string"
            },
            "title": {
              "type": "string"
            },
            "total": {
              "type": [
                "integer",
                "number"
              ]
            },
            "all": {
              "type": "string"
            }
          },
          "required": [
            "@graph"
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

**执行与来源：** `network`；上游接口名 `ENCODE_search_biosamples`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("encode_search_biosamples", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.encode_search_chromatin_accessibility` {#dnakit.databases.encode_search_chromatin_accessibility}

Search ENCODE chromatin accessibility experiments (ATAC-seq and DNase-seq). ATAC-seq and DNase-seq identify open chromatin regions including promoters, enhancers, and other regulatory elements. Returns experiment accessions, biosample details, and accessibility assay metadata. Open chromatin regions are sites of active gene regulation where transcription factors can bind. Use this to find chromatin accessibility data for specific tissues or cell types.

```text
dnakit.databases.encode_search_chromatin_accessibility(*, assay_type=None, biosample_term_name=None, organism=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `assay_type` | string | None | Chromatin accessibility assay: &#x27;ATAC-seq&#x27; (Assay for Transposase-Accessible Chromatin, newer/preferred) or &#x27;DNase-seq&#x27; (DNase I hypersensitive sites, larger dataset). {&quot;default&quot;: &quot;ATAC-seq&quot;} |
| `biosample_term_name` | string &#124; null | None | Biosample name filter (e.g., &#x27;K562&#x27;, &#x27;GM12878&#x27;, &#x27;liver&#x27;). Leave empty for all. |
| `organism` | string | None | Organism scientific name. {&quot;default&quot;: &quot;Homo sapiens&quot;} |
| `limit` | integer | None | Maximum number of results to return. {&quot;default&quot;: 25, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "assay_type": {
          "type": "string",
          "default": "ATAC-seq",
          "description": "Chromatin accessibility assay: 'ATAC-seq' (Assay for Transposase-Accessible Chromatin, newer/preferred) or 'DNase-seq' (DNase I hypersensitive sites, larger dataset)."
        },
        "biosample_term_name": {
          "type": [
            "string",
            "null"
          ],
          "description": "Biosample name filter (e.g., 'K562', 'GM12878', 'liver'). Leave empty for all."
        },
        "organism": {
          "type": "string",
          "default": "Homo sapiens",
          "description": "Organism scientific name."
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
              "type": "integer"
            },
            "experiments": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "accession": {
                    "type": "string"
                  },
                  "assay_title": {
                    "type": "string"
                  },
                  "biosample_summary": {
                    "type": "string"
                  },
                  "status": {
                    "type": "string"
                  },
                  "lab": {
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

**执行与来源：** `network`；上游接口名 `ENCODE_search_chromatin_accessibility`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("encode_search_chromatin_accessibility", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.encode_search_experiments` {#dnakit.databases.encode_search_experiments}

Search ENCODE experiments

```text
dnakit.databases.encode_search_experiments(*, assay_title=None, target=None, organism=None, status=None, limit=None, biosample=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `assay_title` | string | None | 实验方法名称，例如 ATAC-seq。 |
| `target` | string | None | 实验靶标或目标标识，按当前数据源填写。 |
| `organism` | string | None | 目标生物名称，按提供方命名规则填写。 |
| `status` | string | None | 记录的发布或处理状态筛选条件。 {&quot;default&quot;: &quot;released&quot;} |
| `limit` | integer | None | 当前请求的返回数量上限；不是全库命中总数。 {&quot;default&quot;: 10} |
| `biosample` | string | None | 生物样本名称或检索条件。 |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回类型为 object，具体结构由所调用操作决定；不保证统一的 data 层。提供方错误会抛出异常。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "assay_title": {
          "type": "string"
        },
        "target": {
          "type": "string"
        },
        "organism": {
          "type": "string"
        },
        "status": {
          "type": "string",
          "default": "released"
        },
        "limit": {
          "type": "integer",
          "default": 10
        },
        "biosample": {
          "type": "string"
        }
      },
      "required": []
    }
    ```

??? info "DNAKit 参数扩展规则"

    ```json
    {
      "parameters": {
        "biosample": {
          "type": "string"
        }
      },
      "optional": []
    }
    ```

**执行与来源：** `network`；上游接口名 `ENCODE_search_experiments`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("encode_search_experiments", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.encode_search_histone_experiments` {#dnakit.databases.encode_search_histone_experiments}

Search ENCODE histone ChIP-seq experiments by histone modification mark, biosample, or organism. Returns experiment accessions, histone marks (H3K4me3, H3K27ac, H3K27me3, H3K36me3, H3K4me1, H3K9me3), biosample summaries, and metadata. Use this to find histone modification profiling data for specific marks or tissues. Common marks: H3K4me3 (active promoters), H3K27ac (active enhancers), H3K27me3 (polycomb repression), H3K4me1 (poised enhancers), H3K36me3 (gene bodies), H3K9me3 (heterochromatin).

```text
dnakit.databases.encode_search_histone_experiments(*, histone_mark=None, target=None, biosample_term_name=None, biosample=None, organism=None, limit=None, biosample_term=None, cell_type=None, tissue=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `histone_mark` | string &#124; null | None | Histone modification mark to filter by (e.g., &#x27;H3K4me3&#x27;, &#x27;H3K27ac&#x27;, &#x27;H3K27me3&#x27;, &#x27;H3K4me1&#x27;, &#x27;H3K36me3&#x27;, &#x27;H3K9me3&#x27;). Leave empty to search all histone marks. |
| `target` | string &#124; null | None | Alias for histone_mark. Histone modification mark (e.g., &#x27;H3K27ac&#x27;, &#x27;H3K4me3&#x27;). |
| `biosample_term_name` | string &#124; null | None | Biosample name from ENCODE ontology (cell lines or tissues, NOT disease names). Examples: &#x27;K562&#x27;, &#x27;HepG2&#x27;, &#x27;GM12878&#x27;, &#x27;liver&#x27;, &#x27;brain&#x27;. For disease-specific data, search by relevant cell lines (e.g., K562/MOLM-13 for AML). Leave empty to search all. |
| `biosample` | string &#124; null | None | Alias for biosample_term_name. Biosample (tissue or cell line, e.g., &#x27;liver&#x27;, &#x27;K562&#x27;). |
| `organism` | string | None | Organism scientific name (e.g., &#x27;Homo sapiens&#x27;, &#x27;Mus musculus&#x27;). {&quot;default&quot;: &quot;Homo sapiens&quot;} |
| `limit` | integer | None | Maximum number of results to return (1-100). {&quot;default&quot;: 25, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `biosample_term` | string &#124; null | None | Alias for biosample_term_name. Biosample tissue or cell line (e.g., &quot;breast epithelium&quot;, &quot;liver&quot;, &quot;K562&quot;). |
| `cell_type` | string &#124; null | None | Alias for biosample_term_name. Cell type or tissue (e.g., &#x27;GM12878&#x27;, &#x27;K562&#x27;, &#x27;liver&#x27;, &#x27;brain&#x27;). |
| `tissue` | string &#124; null | None | Alias for biosample_term_name. Tissue or cell line (e.g., &#x27;liver&#x27;, &#x27;brain&#x27;, &#x27;breast epithelium&#x27;, &#x27;MCF-7&#x27;). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "histone_mark": {
          "type": [
            "string",
            "null"
          ],
          "description": "Histone modification mark to filter by (e.g., 'H3K4me3', 'H3K27ac', 'H3K27me3', 'H3K4me1', 'H3K36me3', 'H3K9me3'). Leave empty to search all histone marks."
        },
        "target": {
          "type": [
            "string",
            "null"
          ],
          "description": "Alias for histone_mark. Histone modification mark (e.g., 'H3K27ac', 'H3K4me3')."
        },
        "biosample_term_name": {
          "type": [
            "string",
            "null"
          ],
          "description": "Biosample name from ENCODE ontology (cell lines or tissues, NOT disease names). Examples: 'K562', 'HepG2', 'GM12878', 'liver', 'brain'. For disease-specific data, search by relevant cell lines (e.g., K562/MOLM-13 for AML). Leave empty to search all."
        },
        "biosample": {
          "type": [
            "string",
            "null"
          ],
          "description": "Alias for biosample_term_name. Biosample (tissue or cell line, e.g., 'liver', 'K562')."
        },
        "organism": {
          "type": "string",
          "default": "Homo sapiens",
          "description": "Organism scientific name (e.g., 'Homo sapiens', 'Mus musculus')."
        },
        "limit": {
          "type": "integer",
          "default": 25,
          "minimum": 1,
          "maximum": 100,
          "description": "Maximum number of results to return (1-100)."
        },
        "biosample_term": {
          "type": [
            "string",
            "null"
          ],
          "description": "Alias for biosample_term_name. Biosample tissue or cell line (e.g., \"breast epithelium\", \"liver\", \"K562\")."
        },
        "cell_type": {
          "type": [
            "string",
            "null"
          ],
          "description": "Alias for biosample_term_name. Cell type or tissue (e.g., 'GM12878', 'K562', 'liver', 'brain')."
        },
        "tissue": {
          "type": [
            "string",
            "null"
          ],
          "description": "Alias for biosample_term_name. Tissue or cell line (e.g., 'liver', 'brain', 'breast epithelium', 'MCF-7')."
        }
      },
      "additionalProperties": false,
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
              "type": "integer"
            },
            "experiments": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "accession": {
                    "type": "string"
                  },
                  "histone_mark": {
                    "type": "string"
                  },
                  "biosample_summary": {
                    "type": "string"
                  },
                  "status": {
                    "type": "string"
                  },
                  "lab": {
                    "type": "string"
                  },
                  "date_released": {
                    "type": [
                      "string",
                      "null"
                    ]
                  }
                }
              }
            }
          },
          "required": [
            "experiments"
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

**执行与来源：** `network`；上游接口名 `ENCODE_search_histone_experiments`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("encode_search_histone_experiments", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.encode_search_rnaseq_experiments` {#dnakit.databases.encode_search_rnaseq_experiments}

Search ENCODE RNA-seq experiments by biosample (cell type or tissue), organism, and assay type (total RNA-seq, polyA plus RNA-seq, small RNA-seq, microRNA-seq). Returns experiment accessions, biosample summaries, lab names, and metadata from the ENCODE project. Use this to find ENCODE-curated RNA-seq data for specific cell lines or tissues.

```text
dnakit.databases.encode_search_rnaseq_experiments(*, biosample_term_name=None, biosample=None, cell_type=None, tissue=None, organism=None, assay_type=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `biosample_term_name` | string &#124; null | None | Biosample name from ENCODE ontology (e.g., &#x27;K562&#x27;, &#x27;HepG2&#x27;, &#x27;GM12878&#x27;, &#x27;liver&#x27;, &#x27;brain&#x27;). Leave empty to search all biosamples. |
| `biosample` | string &#124; null | None | Alias for biosample_term_name. Cell type or tissue (e.g., &#x27;K562&#x27;, &#x27;liver&#x27;). |
| `cell_type` | string &#124; null | None | Alias for biosample_term_name. Cell type (e.g., &#x27;K562&#x27;, &#x27;HepG2&#x27;, &#x27;GM12878&#x27;). |
| `tissue` | string &#124; null | None | Alias for biosample_term_name. Tissue type (e.g., &#x27;liver&#x27;, &#x27;brain&#x27;, &#x27;heart&#x27;). |
| `organism` | string | None | Organism scientific name (e.g., &#x27;Homo sapiens&#x27;, &#x27;Mus musculus&#x27;). {&quot;default&quot;: &quot;Homo sapiens&quot;} |
| `assay_type` | string | None | RNA-seq assay type: &#x27;total RNA-seq&#x27; (default, all RNA), &#x27;polyA plus RNA-seq&#x27; (mRNA), &#x27;small RNA-seq&#x27;, &#x27;microRNA-seq&#x27;. Aliases: &#x27;total&#x27;, &#x27;polya&#x27;, &#x27;mirna&#x27;, &#x27;small&#x27;. {&quot;default&quot;: &quot;total RNA-seq&quot;} |
| `limit` | integer | None | Maximum number of results to return (1-100). {&quot;default&quot;: 25, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "biosample_term_name": {
          "type": [
            "string",
            "null"
          ],
          "description": "Biosample name from ENCODE ontology (e.g., 'K562', 'HepG2', 'GM12878', 'liver', 'brain'). Leave empty to search all biosamples."
        },
        "biosample": {
          "type": [
            "string",
            "null"
          ],
          "description": "Alias for biosample_term_name. Cell type or tissue (e.g., 'K562', 'liver')."
        },
        "cell_type": {
          "type": [
            "string",
            "null"
          ],
          "description": "Alias for biosample_term_name. Cell type (e.g., 'K562', 'HepG2', 'GM12878')."
        },
        "tissue": {
          "type": [
            "string",
            "null"
          ],
          "description": "Alias for biosample_term_name. Tissue type (e.g., 'liver', 'brain', 'heart')."
        },
        "organism": {
          "type": "string",
          "default": "Homo sapiens",
          "description": "Organism scientific name (e.g., 'Homo sapiens', 'Mus musculus')."
        },
        "assay_type": {
          "type": "string",
          "default": "total RNA-seq",
          "description": "RNA-seq assay type: 'total RNA-seq' (default, all RNA), 'polyA plus RNA-seq' (mRNA), 'small RNA-seq', 'microRNA-seq'. Aliases: 'total', 'polya', 'mirna', 'small'."
        },
        "limit": {
          "type": "integer",
          "default": 25,
          "minimum": 1,
          "maximum": 100,
          "description": "Maximum number of results to return (1-100)."
        }
      },
      "additionalProperties": false,
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
              "type": "integer"
            },
            "experiments": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "accession": {
                    "type": "string"
                  },
                  "assay_title": {
                    "type": "string"
                  },
                  "biosample_summary": {
                    "type": "string"
                  },
                  "status": {
                    "type": "string"
                  },
                  "lab": {
                    "type": "string"
                  },
                  "date_released": {
                    "type": [
                      "string",
                      "null"
                    ]
                  }
                }
              }
            }
          },
          "required": [
            "total"
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

**执行与来源：** `network`；上游接口名 `ENCODE_search_rnaseq_experiments`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("encode_search_rnaseq_experiments", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.expression_atlas_search_differential` {#dnakit.databases.expression_atlas_search_differential}

Search for differential expression EXPERIMENTS from EBI Expression Atlas by condition text and/or a gene mentioned in the experiment description. Supports disease-context expression analysis. Caveat: `condition` does real text filtering, but `gene` is a text-search coincidence (experiment descriptions rarely name individual genes) and does not return actual per-gene differential-expression values -- treat gene-tagged results as a starting point for manual review, not a validated gene-specific result set.

```text
dnakit.databases.expression_atlas_search_differential(*, gene=None, condition=None, species=None, limit=None, offset=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | None | Gene symbol or Ensembl ID (e.g., &#x27;TP53&#x27;, &#x27;ENSG00000141510&#x27;) |
| `condition` | string | None | Condition/disease to filter by (e.g., &#x27;cancer&#x27;, &#x27;inflammation&#x27;, &#x27;breast&#x27;) |
| `species` | string | None | Species name (default: &#x27;homo sapiens&#x27;) {&quot;default&quot;: &quot;homo sapiens&quot;} |
| `limit` | integer &#124; null | None | Maximum number of experiments to return (default 50, max 500). Compare with the total reported in the response to see how many exist. {&quot;default&quot;: 50} |
| `offset` | integer &#124; null | None | Zero-based index of the first experiment to return, for paging through results beyond the first page (e.g. offset=50 returns 51-100). Default: 0. {&quot;default&quot;: 0} |
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
          "description": "Gene symbol or Ensembl ID (e.g., 'TP53', 'ENSG00000141510')"
        },
        "condition": {
          "type": "string",
          "description": "Condition/disease to filter by (e.g., 'cancer', 'inflammation', 'breast')"
        },
        "species": {
          "type": "string",
          "description": "Species name (default: 'homo sapiens')",
          "default": "homo sapiens"
        },
        "limit": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Maximum number of experiments to return (default 50, max 500). Compare with the total reported in the response to see how many exist.",
          "default": 50
        },
        "offset": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Zero-based index of the first experiment to return, for paging through results beyond the first page (e.g. offset=50 returns 51-100). Default: 0.",
          "default": 0
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
            "gene": {
              "type": "string"
            },
            "condition": {
              "type": "string"
            },
            "species": {
              "type": "string"
            },
            "experiments": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "experiment_accession": {
                    "type": "string"
                  },
                  "experiment_type": {
                    "type": "string"
                  },
                  "experiment_description": {
                    "type": "string"
                  },
                  "species": {
                    "type": "string"
                  },
                  "num_assays": {
                    "type": [
                      "integer",
                      "number"
                    ]
                  },
                  "gene_mentioned": {
                    "type": "boolean"
                  },
                  "factors": {
                    "type": "array",
                    "items": {
                      "type": "string"
                    }
                  }
                }
              }
            },
            "experiment_count": {
              "type": [
                "integer",
                "number"
              ]
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

**执行与来源：** `network`；上游接口名 `ExpressionAtlas_search_differential`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("expression_atlas_search_differential", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.expression_atlas_search_experiments` {#dnakit.databases.expression_atlas_search_experiments}

Search all EBI Expression Atlas experiments (baseline + differential) by condition text and/or a gene mentioned in the experiment description. Returns experiment accessions, types, descriptions, and assay counts. Caveat: `condition` does real text filtering, but `gene` is a text-search coincidence (experiment descriptions rarely name individual genes), so a gene-only query effectively returns the unfiltered catalog -- use to find relevant expression datasets, not as a validated per-gene result set.

```text
dnakit.databases.expression_atlas_search_experiments(*, gene=None, condition=None, species=None, limit=None, offset=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | None | Gene symbol or Ensembl ID |
| `condition` | string | None | Biological condition or tissue (e.g., &#x27;brain&#x27;, &#x27;cancer&#x27;, &#x27;kidney&#x27;) |
| `species` | string | None | Species name (optional, e.g., &#x27;homo sapiens&#x27;) |
| `limit` | integer &#124; null | None | Maximum number of experiments to return (default 50, max 500). Compare with the total reported in the response to see how many exist. {&quot;default&quot;: 50} |
| `offset` | integer &#124; null | None | Zero-based index of the first experiment to return, for paging through results beyond the first page (e.g. offset=50 returns 51-100). Default: 0. {&quot;default&quot;: 0} |
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
          "description": "Gene symbol or Ensembl ID"
        },
        "condition": {
          "type": "string",
          "description": "Biological condition or tissue (e.g., 'brain', 'cancer', 'kidney')"
        },
        "species": {
          "type": "string",
          "description": "Species name (optional, e.g., 'homo sapiens')"
        },
        "limit": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Maximum number of experiments to return (default 50, max 500). Compare with the total reported in the response to see how many exist.",
          "default": 50
        },
        "offset": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Zero-based index of the first experiment to return, for paging through results beyond the first page (e.g. offset=50 returns 51-100). Default: 0.",
          "default": 0
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
            "gene": {
              "type": "string"
            },
            "condition": {
              "type": "string"
            },
            "species": {
              "type": "string"
            },
            "experiments": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "experiment_accession": {
                    "type": "string"
                  },
                  "experiment_type": {
                    "type": "string"
                  },
                  "experiment_description": {
                    "type": "string"
                  },
                  "species": {
                    "type": "string"
                  },
                  "num_assays": {
                    "type": [
                      "integer",
                      "number"
                    ]
                  },
                  "gene_mentioned": {
                    "type": [
                      "boolean",
                      "null"
                    ]
                  }
                }
              }
            },
            "total_count": {
              "type": [
                "integer",
                "number"
              ]
            },
            "gene_specific_count": {
              "type": [
                "integer",
                "number"
              ]
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

**执行与来源：** `network`；上游接口名 `ExpressionAtlas_search_experiments`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("expression_atlas_search_experiments", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.four_dn_search_data` {#dnakit.databases.four_dn_search_data}

Search 4DN Data Portal for Hi-C chromatin conformation files and experiments. Access 350+ uniformly processed Hi-C contact files from multiple technologies (in situ Hi-C, Capture Hi-C, Micro-C XL). Search by assay type, cell type, or keywords. Returns file metadata with download URLs. Note: Downloads require free 4DN account (create at data.4dnucleome.org). Use for: finding Hi-C datasets, discovering 3D genome data, accessing TAD (Topologically Associating Domain) boundaries, A/B compartments, chromatin loop calls.

```text
dnakit.databases.four_dn_search_data(*, operation=None, query=None, item_type=None, file_type=None, assay_title=None, biosource_name=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | 提供方操作模式；默认值与允许值见本行约束。 {&quot;enum&quot;: &#91;&quot;search&quot;&#93;, &quot;default&quot;: &quot;search&quot;} |
| `query` | string | None | Search query (default: &#x27;*&#x27; for all) {&quot;default&quot;: &quot;*&quot;} |
| `item_type` | string | None | Type of item to search {&quot;enum&quot;: &#91;&quot;File&quot;, &quot;Experiment&quot;, &quot;Biosource&quot;&#93;, &quot;default&quot;: &quot;File&quot;} |
| `file_type` | string | None | Filter by file type. Common types: &#x27;contact list&#x27; (processed contact matrices), &#x27;pairs&#x27; (aligned read pairs), &#x27;cooler&#x27; (multi-resolution contact matrices), &#x27;mcool&#x27; (multi-resolution cooler), &#x27;hic&#x27; (Juicer format). |
| `assay_title` | string | None | Filter by assay (e.g., &#x27;Hi-C&#x27;, &#x27;Micro-C&#x27;, &#x27;Capture Hi-C&#x27;) |
| `biosource_name` | string | None | Filter by cell type (e.g., &#x27;H1-hESC&#x27;, &#x27;HFFc6&#x27;, &#x27;GM12878&#x27;) |
| `limit` | integer | None | Maximum results to return {&quot;default&quot;: 25, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
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
            "search"
          ],
          "default": "search"
        },
        "query": {
          "type": "string",
          "description": "Search query (default: '*' for all)",
          "default": "*"
        },
        "item_type": {
          "type": "string",
          "description": "Type of item to search",
          "enum": [
            "File",
            "Experiment",
            "Biosource"
          ],
          "default": "File"
        },
        "file_type": {
          "type": "string",
          "description": "Filter by file type. Common types: 'contact list' (processed contact matrices), 'pairs' (aligned read pairs), 'cooler' (multi-resolution contact matrices), 'mcool' (multi-resolution cooler), 'hic' (Juicer format)."
        },
        "assay_title": {
          "type": "string",
          "description": "Filter by assay (e.g., 'Hi-C', 'Micro-C', 'Capture Hi-C')"
        },
        "biosource_name": {
          "type": "string",
          "description": "Filter by cell type (e.g., 'H1-hESC', 'HFFc6', 'GM12878')"
        },
        "limit": {
          "type": "integer",
          "description": "Maximum results to return",
          "default": 25,
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
      "type": "object",
      "properties": {
        "status": {
          "type": "string"
        },
        "num_results": {
          "type": "integer"
        },
        "total": {
          "type": "integer"
        },
        "results": {
          "type": "array"
        },
        "search_url": {
          "type": "string"
        },
        "error": {
          "type": "string"
        }
      },
      "additionalProperties": true
    }
    ```

**执行与来源：** `network`；上游接口名 `FourDN_search_data`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("four_dn_search_data", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gdc_get_clinical_data` {#dnakit.databases.gdc_get_clinical_data}

Get detailed clinical data for cancer cases from NCI GDC/TCGA. Returns demographics (gender, race, vital_status, age_at_index), diagnoses (primary_diagnosis, tumor_stage, age_at_diagnosis, days_to_last_follow_up), and treatments (therapeutic_agents, treatment_type). Filter by project, primary_site, disease_type, or vital_status.

```text
dnakit.databases.gdc_get_clinical_data(*, project_id=None, primary_site=None, disease_type=None, vital_status=None, gender=None, size=None, offset=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `project_id` | string | None | GDC project identifier (e.g., &#x27;TCGA-BRCA&#x27;, &#x27;TCGA-LUAD&#x27;, &#x27;TARGET-AML&#x27;) |
| `primary_site` | string | None | Primary anatomical site (e.g., &#x27;Breast&#x27;, &#x27;Lung&#x27;, &#x27;Brain&#x27;) |
| `disease_type` | string | None | Disease type filter (e.g., &#x27;Ductal and Lobular Neoplasms&#x27;) |
| `vital_status` | string | None | Vital status filter: &#x27;Alive&#x27; or &#x27;Dead&#x27; {&quot;enum&quot;: &#91;&quot;Alive&quot;, &quot;Dead&quot;&#93;} |
| `gender` | string | None | Gender filter: &#x27;female&#x27; or &#x27;male&#x27; {&quot;enum&quot;: &#91;&quot;female&quot;, &quot;male&quot;&#93;} |
| `size` | integer | None | Number of cases to return (1-100) {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `offset` | integer | None | Pagination offset (0-based) {&quot;default&quot;: 0, &quot;minimum&quot;: 0} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回类型为 object，具体结构由所调用操作决定；不保证统一的 data 层。提供方错误会抛出异常。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "project_id": {
          "type": "string",
          "description": "GDC project identifier (e.g., 'TCGA-BRCA', 'TCGA-LUAD', 'TARGET-AML')"
        },
        "primary_site": {
          "type": "string",
          "description": "Primary anatomical site (e.g., 'Breast', 'Lung', 'Brain')"
        },
        "disease_type": {
          "type": "string",
          "description": "Disease type filter (e.g., 'Ductal and Lobular Neoplasms')"
        },
        "vital_status": {
          "type": "string",
          "description": "Vital status filter: 'Alive' or 'Dead'",
          "enum": [
            "Alive",
            "Dead"
          ]
        },
        "gender": {
          "type": "string",
          "description": "Gender filter: 'female' or 'male'",
          "enum": [
            "female",
            "male"
          ]
        },
        "size": {
          "type": "integer",
          "default": 10,
          "minimum": 1,
          "maximum": 100,
          "description": "Number of cases to return (1-100)"
        },
        "offset": {
          "type": "integer",
          "default": 0,
          "minimum": 0,
          "description": "Pagination offset (0-based)"
        }
      },
      "required": []
    }
    ```

**执行与来源：** `network`；上游接口名 `GDC_get_clinical_data`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gdc_get_clinical_data", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gdc_get_mutation_frequency` {#dnakit.databases.gdc_get_mutation_frequency}

Get pan-cancer mutation frequency statistics for a gene across all TCGA projects. Returns overall and per-project mutation rates. Note: this tool is pan-cancer only and does not support filtering by cancer type.

```text
dnakit.databases.gdc_get_mutation_frequency(*, gene_symbol=None, gene=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_symbol` | string | None | Gene symbol (e.g., &#x27;TP53&#x27;, &#x27;KRAS&#x27;) |
| `gene` | string | None | Gene symbol alias — alternative to gene_symbol |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回类型为 object，具体结构由所调用操作决定；不保证统一的 data 层。提供方错误会抛出异常。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene_symbol": {
          "type": "string",
          "description": "Gene symbol (e.g., 'TP53', 'KRAS')"
        },
        "gene": {
          "type": "string",
          "description": "Gene symbol alias — alternative to gene_symbol"
        }
      },
      "required": []
    }
    ```

**执行与来源：** `network`；上游接口名 `GDC_get_mutation_frequency`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gdc_get_mutation_frequency", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gdc_get_ssm_by_gene` {#dnakit.databases.gdc_get_ssm_by_gene}

Get somatic mutations (SSMs) for a gene across TCGA/GDC projects

```text
dnakit.databases.gdc_get_ssm_by_gene(*, gene_symbol, project_id=None, size=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 必填（输入 Schema） | Gene symbol (e.g., &#x27;TP53&#x27;, &#x27;EGFR&#x27;, &#x27;BRAF&#x27;) |
| `project_id` | string | None | Optional: Filter by project (e.g., &#x27;TCGA-BRCA&#x27;) |
| `size` | integer | None | Number of results (1–100) {&quot;default&quot;: 20, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回类型为 object，具体结构由所调用操作决定；不保证统一的 data 层。提供方错误会抛出异常。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene_symbol": {
          "type": "string",
          "description": "Gene symbol (e.g., 'TP53', 'EGFR', 'BRAF')"
        },
        "project_id": {
          "type": "string",
          "description": "Optional: Filter by project (e.g., 'TCGA-BRCA')"
        },
        "size": {
          "type": "integer",
          "default": 20,
          "minimum": 1,
          "maximum": 100,
          "description": "Number of results (1–100)"
        }
      },
      "required": [
        "gene_symbol"
      ]
    }
    ```

**执行与来源：** `network`；上游接口名 `GDC_get_ssm_by_gene`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gdc_get_ssm_by_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gdc_get_survival` {#dnakit.databases.gdc_get_survival}

Get Kaplan-Meier survival data for a GDC/TCGA cancer cohort. Returns time-to-event data with censoring status and survival estimates for each patient. Filter by project and optionally by gene mutation status. Use for overall survival analysis of TCGA cancer types.

```text
dnakit.databases.gdc_get_survival(*, project_id, gene_symbol=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `project_id` | string | 必填（输入 Schema） | GDC project identifier (e.g., &#x27;TCGA-BRCA&#x27;, &#x27;TCGA-LUAD&#x27;, &#x27;TCGA-GBM&#x27;) |
| `gene_symbol` | string | None | Optional: gene symbol to filter cases with mutations in this gene (e.g., &#x27;TP53&#x27;, &#x27;KRAS&#x27;) |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回类型为 object，具体结构由所调用操作决定；不保证统一的 data 层。提供方错误会抛出异常。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "project_id": {
          "type": "string",
          "description": "GDC project identifier (e.g., 'TCGA-BRCA', 'TCGA-LUAD', 'TCGA-GBM')"
        },
        "gene_symbol": {
          "type": "string",
          "description": "Optional: gene symbol to filter cases with mutations in this gene (e.g., 'TP53', 'KRAS')"
        }
      },
      "required": [
        "project_id"
      ]
    }
    ```

**执行与来源：** `network`；上游接口名 `GDC_get_survival`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gdc_get_survival", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gdc_list_files` {#dnakit.databases.gdc_list_files}

List NCI GDC files via /files with optional data_type filter

```text
dnakit.databases.gdc_list_files(*, data_type=None, size=None, offset=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `data_type` | string | None | Data type filter (e.g., &#x27;Gene Expression Quantification&#x27;) |
| `size` | integer | None | Number of results (1–100) {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `offset` | integer | None | Offset for pagination (0-based) {&quot;default&quot;: 0, &quot;minimum&quot;: 0} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回类型为 object，具体结构由所调用操作决定；不保证统一的 data 层。提供方错误会抛出异常。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "data_type": {
          "type": "string",
          "description": "Data type filter (e.g., 'Gene Expression Quantification')"
        },
        "size": {
          "type": "integer",
          "default": 10,
          "minimum": 1,
          "maximum": 100,
          "description": "Number of results (1–100)"
        },
        "offset": {
          "type": "integer",
          "default": 0,
          "minimum": 0,
          "description": "Offset for pagination (0-based)"
        }
      },
      "required": []
    }
    ```

**执行与来源：** `network`；上游接口名 `GDC_list_files`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gdc_list_files", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gdc_list_projects` {#dnakit.databases.gdc_list_projects}

List GDC projects (TCGA, TARGET, etc.) with summary statistics

```text
dnakit.databases.gdc_list_projects(*, program=None, size=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `program` | string | None | Filter by program (e.g., &#x27;TCGA&#x27;, &#x27;TARGET&#x27;) |
| `size` | integer | None | Number of results (1–100) {&quot;default&quot;: 20, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回类型为 object，具体结构由所调用操作决定；不保证统一的 data 层。提供方错误会抛出异常。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "program": {
          "type": "string",
          "description": "Filter by program (e.g., 'TCGA', 'TARGET')"
        },
        "size": {
          "type": "integer",
          "default": 20,
          "minimum": 1,
          "maximum": 100,
          "description": "Number of results (1–100)"
        }
      },
      "required": []
    }
    ```

**执行与来源：** `network`；上游接口名 `GDC_list_projects`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gdc_list_projects", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gdc_search_cases` {#dnakit.databases.gdc_search_cases}

Search NCI GDC cases via /cases

```text
dnakit.databases.gdc_search_cases(*, project_id=None, size=None, offset=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `project_id` | string | None | GDC project identifier (e.g., &#x27;TCGA-BRCA&#x27;) |
| `size` | integer | None | Number of results (1–100) {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `offset` | integer | None | Offset for pagination (0-based) {&quot;default&quot;: 0, &quot;minimum&quot;: 0} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回类型为 object，具体结构由所调用操作决定；不保证统一的 data 层。提供方错误会抛出异常。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "project_id": {
          "type": "string",
          "description": "GDC project identifier (e.g., 'TCGA-BRCA')"
        },
        "size": {
          "type": "integer",
          "default": 10,
          "minimum": 1,
          "maximum": 100,
          "description": "Number of results (1–100)"
        },
        "offset": {
          "type": "integer",
          "default": 0,
          "minimum": 0,
          "description": "Offset for pagination (0-based)"
        }
      },
      "required": []
    }
    ```

**执行与来源：** `network`；上游接口名 `GDC_search_cases`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gdc_search_cases", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.geo_search_atacseq_datasets` {#dnakit.databases.geo_search_atacseq_datasets}

Search NCBI GEO for ATAC-seq chromatin accessibility datasets by query keyword, tissue, or cell type. Returns dataset accessions (GSE IDs), titles, summaries, organism, and sample counts. ATAC-seq measures open chromatin regions and is used to map regulatory elements, enhancers, and TF binding sites genome-wide.

```text
dnakit.databases.geo_search_atacseq_datasets(*, query=None, organism=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string &#124; null | None | Search query (e.g., &#x27;T cell differentiation&#x27;, &#x27;liver hepatocyte&#x27;, &#x27;hematopoiesis&#x27;). Leave empty to search broadly. |
| `organism` | string | None | Organism scientific name (e.g., &#x27;Homo sapiens&#x27;, &#x27;Mus musculus&#x27;). {&quot;default&quot;: &quot;Homo sapiens&quot;} |
| `limit` | integer | None | Maximum number of results to return (1-100). {&quot;default&quot;: 20, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "query": {
          "type": [
            "string",
            "null"
          ],
          "description": "Search query (e.g., 'T cell differentiation', 'liver hepatocyte', 'hematopoiesis'). Leave empty to search broadly."
        },
        "organism": {
          "type": "string",
          "default": "Homo sapiens",
          "description": "Organism scientific name (e.g., 'Homo sapiens', 'Mus musculus')."
        },
        "limit": {
          "type": "integer",
          "default": 20,
          "minimum": 1,
          "maximum": 100,
          "description": "Maximum number of results to return (1-100)."
        }
      },
      "additionalProperties": false,
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
              "type": "integer"
            },
            "datasets": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "accession": {
                    "type": "string"
                  },
                  "title": {
                    "type": "string"
                  },
                  "summary": {
                    "type": "string"
                  },
                  "organism": {
                    "type": "string"
                  },
                  "n_samples": {
                    "type": "integer"
                  },
                  "date_published": {
                    "type": [
                      "string",
                      "null"
                    ]
                  }
                }
              }
            }
          },
          "required": [
            "total"
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

**执行与来源：** `network`；上游接口名 `GEO_search_atacseq_datasets`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("geo_search_atacseq_datasets", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.geo_search_chipseq_datasets` {#dnakit.databases.geo_search_chipseq_datasets}

Search NCBI GEO for ChIP-seq (Chromatin Immunoprecipitation followed by sequencing) datasets. ChIP-seq profiles genome-wide binding of proteins (transcription factors, histones) to DNA. Returns dataset accessions, titles, summaries, and sample counts. Use this to find ChIP-seq studies for specific transcription factors, histone marks, or in specific cell types/diseases.

```text
dnakit.databases.geo_search_chipseq_datasets(*, query, organism=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填（输入 Schema） | Search terms for ChIP-seq datasets (e.g., &#x27;H3K27ac liver&#x27;, &#x27;CTCF cancer&#x27;, &#x27;p53 ChIP-seq&#x27;). |
| `organism` | string | None | Organism filter. {&quot;default&quot;: &quot;Homo sapiens&quot;} |
| `limit` | integer | None | Maximum number of dataset IDs to return. {&quot;default&quot;: 20, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
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
          "description": "Search terms for ChIP-seq datasets (e.g., 'H3K27ac liver', 'CTCF cancer', 'p53 ChIP-seq')."
        },
        "organism": {
          "type": "string",
          "default": "Homo sapiens",
          "description": "Organism filter."
        },
        "limit": {
          "type": "integer",
          "default": 20,
          "minimum": 1,
          "maximum": 100,
          "description": "Maximum number of dataset IDs to return."
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
            "total": {
              "type": "integer"
            },
            "datasets": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "accession": {
                    "type": "string"
                  },
                  "title": {
                    "type": "string"
                  },
                  "summary": {
                    "type": "string"
                  },
                  "organism": {
                    "type": "string"
                  },
                  "n_samples": {
                    "type": "integer"
                  },
                  "date_published": {
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

**执行与来源：** `network`；上游接口名 `GEO_search_chipseq_datasets`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("geo_search_chipseq_datasets", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.geo_search_rnaseq_datasets` {#dnakit.databases.geo_search_rnaseq_datasets}

Search NCBI GEO for RNA-seq gene expression datasets by query keyword, disease, tissue, or gene of interest. Returns dataset accessions (GSE IDs), titles, summaries, organism, and sample counts. Use this to find transcriptomics experiments relevant to a biological question. Covers bulk RNA-seq and total RNA-seq studies across all tissues and diseases.

```text
dnakit.databases.geo_search_rnaseq_datasets(*, query=None, organism=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string &#124; null | None | Search query (e.g., &#x27;BRCA1 breast cancer&#x27;, &#x27;T cell activation&#x27;, &#x27;Alzheimer disease&#x27;). Leave empty to search broadly. |
| `organism` | string | None | Organism scientific name (e.g., &#x27;Homo sapiens&#x27;, &#x27;Mus musculus&#x27;). {&quot;default&quot;: &quot;Homo sapiens&quot;} |
| `limit` | integer | None | Maximum number of results to return (1-100). {&quot;default&quot;: 20, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "query": {
          "type": [
            "string",
            "null"
          ],
          "description": "Search query (e.g., 'BRCA1 breast cancer', 'T cell activation', 'Alzheimer disease'). Leave empty to search broadly."
        },
        "organism": {
          "type": "string",
          "default": "Homo sapiens",
          "description": "Organism scientific name (e.g., 'Homo sapiens', 'Mus musculus')."
        },
        "limit": {
          "type": "integer",
          "default": 20,
          "minimum": 1,
          "maximum": 100,
          "description": "Maximum number of results to return (1-100)."
        }
      },
      "additionalProperties": false,
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
              "type": "integer"
            },
            "datasets": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "accession": {
                    "type": "string"
                  },
                  "title": {
                    "type": "string"
                  },
                  "summary": {
                    "type": "string"
                  },
                  "organism": {
                    "type": "string"
                  },
                  "n_samples": {
                    "type": "integer"
                  },
                  "date_published": {
                    "type": [
                      "string",
                      "null"
                    ]
                  }
                }
              }
            }
          },
          "required": [
            "total"
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

**执行与来源：** `network`；上游接口名 `GEO_search_rnaseq_datasets`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("geo_search_rnaseq_datasets", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gtex_calculate_eqtl` {#dnakit.databases.gtex_calculate_eqtl}

Calculate custom eQTL for any gene-variant pair in any tissue. Dynamically calculates gene-variant association statistics including p-value, normalized effect size (NES), t-statistic, and genotype data. Works for both significant and non-significant associations. Returns detailed statistics even if association not in precomputed results. Use for: testing specific hypotheses, custom variant-gene pairs, dynamic eQTL analysis, validating findings. Requires all three: gencode_id, variant_id, and tissue_site_detail_id.

```text
dnakit.databases.gtex_calculate_eqtl(*, operation=None, gencode_id, variant_id, tissue_site_detail_id, dataset_id=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;calculate_eqtl&quot;&#93;, &quot;default&quot;: &quot;calculate_eqtl&quot;} |
| `gencode_id` | string | 必填（输入 Schema） | Required: Versioned GENCODE ID (e.g. ENSG00000141510.18) |
| `variant_id` | string | 必填（输入 Schema） | Required: GTEx variant ID |
| `tissue_site_detail_id` | string | 必填（输入 Schema） | Required: Tissue ID (e.g. &#x27;Whole_Blood&#x27;) |
| `dataset_id` | string | None | 提供方的数据集标识。 {&quot;enum&quot;: &#91;&quot;gtex_v8&quot;, &quot;gtex_v10&quot;, &quot;gtex_snrnaseq_pilot&quot;&#93;, &quot;default&quot;: &quot;gtex_v8&quot;} |
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
            "calculate_eqtl"
          ],
          "description": "Operation type",
          "default": "calculate_eqtl"
        },
        "gencode_id": {
          "type": "string",
          "description": "Required: Versioned GENCODE ID (e.g. ENSG00000141510.18)"
        },
        "variant_id": {
          "type": "string",
          "description": "Required: GTEx variant ID"
        },
        "tissue_site_detail_id": {
          "type": "string",
          "description": "Required: Tissue ID (e.g. 'Whole_Blood')"
        },
        "dataset_id": {
          "type": "string",
          "enum": [
            "gtex_v8",
            "gtex_v10",
            "gtex_snrnaseq_pilot"
          ],
          "default": "gtex_v8"
        }
      },
      "required": [
        "gencode_id",
        "tissue_site_detail_id",
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
            "data": {
              "type": "array",
              "items": {
                "type": "number"
              }
            },
            "error": {
              "type": "number"
            },
            "gencodeId": {
              "type": "string"
            },
            "geneSymbol": {
              "type": "string"
            },
            "genotypes": {
              "type": "array",
              "items": {
                "type": [
                  "integer",
                  "number"
                ]
              }
            },
            "hetCount": {
              "type": [
                "integer",
                "number"
              ]
            },
            "homoAltCount": {
              "type": [
                "integer",
                "number"
              ]
            },
            "homoRefCount": {
              "type": [
                "integer",
                "number"
              ]
            },
            "maf": {
              "type": "number"
            },
            "nes": {
              "type": "number"
            },
            "pValue": {
              "type": "number"
            },
            "pValueThreshold": {
              "type": "number"
            },
            "tStatistic": {
              "type": "number"
            },
            "tissueSiteDetailId": {
              "type": "string"
            },
            "variantId": {
              "type": "string"
            }
          },
          "required": [
            "data"
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

**执行与来源：** `network`；上游接口名 `GTEx_calculate_eqtl`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gtex_calculate_eqtl", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gtex_get_expression_summary` {#dnakit.databases.gtex_get_expression_summary}

Get GTEx expression summary for a gene via /expression/geneExpression

```text
dnakit.databases.gtex_get_expression_summary(*, gene_symbol=None, ensembl_gene_id=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_symbol` | string | None | Gene symbol (e.g., TP53, BRCA1). Auto-resolved to GENCODE ID. |
| `ensembl_gene_id` | string | None | Ensembl gene ID, e.g., ENSG00000141510 |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回类型为 object，具体结构由所调用操作决定；不保证统一的 data 层。提供方错误会抛出异常。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene_symbol": {
          "type": "string",
          "description": "Gene symbol (e.g., TP53, BRCA1). Auto-resolved to GENCODE ID."
        },
        "ensembl_gene_id": {
          "type": "string",
          "description": "Ensembl gene ID, e.g., ENSG00000141510"
        }
      },
      "required": []
    }
    ```

**执行与来源：** `network`；上游接口名 `GTEx_get_expression_summary`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gtex_get_expression_summary", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gtex_get_median_gene_expression` {#dnakit.databases.gtex_get_median_gene_expression}

Get median gene expression levels across GTEx tissues. Returns median expression in TPM (Transcripts Per Million) for one or more genes across 54 tissue sites. Filter by specific tissues or get expression across all tissues. Defaults to GTEx V8 (2020) -- confirmed live no V11 dataset is queryable via this API (dataset_id accepts gtex_v7/gtex_v8/gtex_v10). Set dataset_id='gtex_v10' for the newer release; gene IDs are resolved to each dataset's own GENCODE version automatically, so V10 returns full results. Use for: tissue-specific expression profiling, identifying highly expressed genes, comparing expression patterns across tissues. Example: Get brain vs liver expression for TP53 gene.

```text
dnakit.databases.gtex_get_median_gene_expression(*, operation=None, gencode_id=None, tissue_site_detail_id=None, dataset_id=None, page=None, items_per_page=None, gene_symbol=None, gene=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;get_median_gene_expression&quot;&#93;, &quot;default&quot;: &quot;get_median_gene_expression&quot;} |
| `gencode_id` | string &#124; array | None | Gene identifier(s): gene symbol (e.g. &#x27;TP53&#x27;), unversioned Ensembl ID (e.g. &#x27;ENSG00000141510&#x27;), or versioned GENCODE ID (e.g. &#x27;ENSG00000141510.18&#x27;). Auto-resolved to versioned GENCODE ID. Can be single string or array. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `tissue_site_detail_id` | array | None | Optional: Tissue IDs to filter (e.g. &#91;&#x27;Liver&#x27;, &#x27;Brain_Cortex&#x27;&#93;). Omit for all tissues. See GTEx_get_tissue_sites for valid IDs {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `dataset_id` | string | None | GTEx dataset version (default: gtex_v8, the most widely cited release). gtex_v10 is the newer release (943 eQTL donors vs 838) and is fully queryable -- gene IDs are resolved to each dataset&#x27;s own GENCODE version automatically. {&quot;enum&quot;: &#91;&quot;gtex_v8&quot;, &quot;gtex_v10&quot;, &quot;gtex_snrnaseq_pilot&quot;&#93;, &quot;default&quot;: &quot;gtex_v8&quot;} |
| `page` | integer | None | Page number for pagination (0-based) {&quot;minimum&quot;: 0, &quot;default&quot;: 0} |
| `items_per_page` | integer | None | Results per page {&quot;minimum&quot;: 1, &quot;maximum&quot;: 100000, &quot;default&quot;: 250} |
| `gene_symbol` | string | None | Gene symbol alias for gencode_id (e.g., &quot;TP53&quot;, &quot;COL5A1&quot;) |
| `gene` | string | None | DNAKit optional alias for gene_symbol |
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
            "get_median_gene_expression"
          ],
          "description": "Operation type",
          "default": "get_median_gene_expression"
        },
        "gencode_id": {
          "type": [
            "string",
            "array"
          ],
          "items": {
            "type": "string"
          },
          "description": "Gene identifier(s): gene symbol (e.g. 'TP53'), unversioned Ensembl ID (e.g. 'ENSG00000141510'), or versioned GENCODE ID (e.g. 'ENSG00000141510.18'). Auto-resolved to versioned GENCODE ID. Can be single string or array."
        },
        "tissue_site_detail_id": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "Optional: Tissue IDs to filter (e.g. ['Liver', 'Brain_Cortex']). Omit for all tissues. See GTEx_get_tissue_sites for valid IDs"
        },
        "dataset_id": {
          "type": "string",
          "enum": [
            "gtex_v8",
            "gtex_v10",
            "gtex_snrnaseq_pilot"
          ],
          "default": "gtex_v8",
          "description": "GTEx dataset version (default: gtex_v8, the most widely cited release). gtex_v10 is the newer release (943 eQTL donors vs 838) and is fully queryable -- gene IDs are resolved to each dataset's own GENCODE version automatically."
        },
        "page": {
          "type": "integer",
          "minimum": 0,
          "default": 0,
          "description": "Page number for pagination (0-based)"
        },
        "items_per_page": {
          "type": "integer",
          "minimum": 1,
          "maximum": 100000,
          "default": 250,
          "description": "Results per page"
        },
        "gene_symbol": {
          "type": "string",
          "description": "Gene symbol alias for gencode_id (e.g., \"TP53\", \"COL5A1\")"
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
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "median": {
                "type": "number"
              },
              "tissueSiteDetailId": {
                "type": "string"
              },
              "ontologyId": {
                "type": "string"
              },
              "datasetId": {
                "type": "string"
              },
              "gencodeId": {
                "type": "string"
              },
              "geneSymbol": {
                "type": "string"
              },
              "unit": {
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

??? info "DNAKit 参数扩展规则"

    ```json
    {
      "aliases": {
        "gene": "gene_symbol"
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `GTEx_get_median_gene_expression`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gtex_get_median_gene_expression", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gtex_get_multi_tissue_eqtls` {#dnakit.databases.gtex_get_multi_tissue_eqtls}

Get multi-tissue eQTL meta-analysis results (Metasoft). Returns m-values (posterior probability of eQTL effect in each tissue), normalized effect sizes, p-values, and standard errors across all tissues for a gene-variant pair. Shows tissue-sharing patterns of eQTL effects. Use for: understanding tissue-specificity of eQTLs, finding shared vs tissue-specific effects, cross-tissue meta-analysis, regulatory mechanism studies. Requires gencode_id, optionally filter by variant_id.

```text
dnakit.databases.gtex_get_multi_tissue_eqtls(*, operation=None, gencode_id, variant_id=None, dataset_id=None, page=None, items_per_page=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;get_multi_tissue_eqtls&quot;&#93;, &quot;default&quot;: &quot;get_multi_tissue_eqtls&quot;} |
| `gencode_id` | string | 必填（输入 Schema） | Required: Versioned GENCODE ID |
| `variant_id` | string | None | Optional: GTEx variant ID to filter specific variant |
| `dataset_id` | string | None | 提供方的数据集标识。 {&quot;enum&quot;: &#91;&quot;gtex_v8&quot;, &quot;gtex_v10&quot;, &quot;gtex_snrnaseq_pilot&quot;&#93;, &quot;default&quot;: &quot;gtex_v8&quot;} |
| `page` | integer | None | 结果页码；起始页约定按提供方参数说明。 {&quot;minimum&quot;: 0, &quot;default&quot;: 0} |
| `items_per_page` | integer | None | 每页返回的记录数。 {&quot;minimum&quot;: 1, &quot;maximum&quot;: 100000, &quot;default&quot;: 250} |
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
            "get_multi_tissue_eqtls"
          ],
          "description": "Operation type",
          "default": "get_multi_tissue_eqtls"
        },
        "gencode_id": {
          "type": "string",
          "description": "Required: Versioned GENCODE ID"
        },
        "variant_id": {
          "type": "string",
          "description": "Optional: GTEx variant ID to filter specific variant"
        },
        "dataset_id": {
          "type": "string",
          "enum": [
            "gtex_v8",
            "gtex_v10",
            "gtex_snrnaseq_pilot"
          ],
          "default": "gtex_v8"
        },
        "page": {
          "type": "integer",
          "minimum": 0,
          "default": 0
        },
        "items_per_page": {
          "type": "integer",
          "minimum": 1,
          "maximum": 100000,
          "default": 250
        }
      },
      "required": [
        "gencode_id"
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
              "gencodeId": {
                "type": "string"
              },
              "datasetId": {
                "type": "string"
              },
              "metaP": {
                "type": "number"
              },
              "variantId": {
                "type": "string"
              },
              "tissues": {
                "type": "object",
                "properties": {
                  "Thyroid": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Testis": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Small_Intestine_Terminal_Ileum": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Brain_Frontal_Cortex_BA9": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Skin_Not_Sun_Exposed_Suprapubic": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Vagina": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Whole_Blood": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Breast_Mammary_Tissue": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Pituitary": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Minor_Salivary_Gland": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Adipose_Subcutaneous": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Adrenal_Gland": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Heart_Atrial_Appendage": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Brain_Amygdala": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Stomach": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Brain_Caudate_basal_ganglia": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Colon_Transverse": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Brain_Cerebellum": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Esophagus_Muscularis": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Liver": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Brain_Hypothalamus": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Prostate": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Pancreas": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Nerve_Tibial": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Spleen": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Colon_Sigmoid": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Brain_Anterior_cingulate_cortex_BA24": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Muscle_Skeletal": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Brain_Hippocampus": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Brain_Cortex": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Heart_Left_Ventricle": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Artery_Tibial": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Uterus": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Ovary": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Esophagus_Gastroesophageal_Junction": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Cells_EBV-transformed_lymphocytes": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Artery_Coronary": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Brain_Spinal_cord_cervical_c-1": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Adipose_Visceral_Omentum": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Brain_Nucleus_accumbens_basal_ganglia": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Brain_Cerebellar_Hemisphere": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Esophagus_Mucosa": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Cells_Cultured_fibroblasts": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Artery_Aorta": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Brain_Putamen_basal_ganglia": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Lung": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Kidney_Cortex": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": [
                          "null",
                          "number"
                        ]
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": [
                          "null",
                          "number"
                        ]
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Brain_Substantia_nigra": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
                        "type": "number"
                      }
                    }
                  },
                  "Skin_Sun_Exposed_Lower_leg": {
                    "type": "object",
                    "properties": {
                      "mValue": {
                        "type": "number"
                      },
                      "pValue": {
                        "type": "number"
                      },
                      "se": {
                        "type": "number"
                      },
                      "nes": {
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

**执行与来源：** `network`；上游接口名 `GTEx_get_multi_tissue_eqtls`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gtex_get_multi_tissue_eqtls", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gtex_get_single_tissue_eqtls` {#dnakit.databases.gtex_get_single_tissue_eqtls}

Get significant single-tissue eQTL associations. Returns precomputed gene-variant associations with p-values, effect sizes (NES), and genomic positions. Query by gene ID, variant ID, tissue, or combination. Results show which genetic variants affect gene expression in specific tissues. Use for: finding eQTLs for candidate genes, looking up variant effects, tissue-specific regulatory analysis, GWAS follow-up. Must provide at least one filter: gencode_id, variant_id, or tissue_site_detail_id.

```text
dnakit.databases.gtex_get_single_tissue_eqtls(*, operation=None, gencode_id=None, variant_id=None, tissue_site_detail_id=None, dataset_id=None, page=None, items_per_page=None, gene_symbol=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;get_single_tissue_eqtls&quot;&#93;, &quot;default&quot;: &quot;get_single_tissue_eqtls&quot;} |
| `gencode_id` | array | None | Optional: Versioned GENCODE ID(s) to query {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `variant_id` | array | None | Optional: GTEx variant ID(s) to query {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `tissue_site_detail_id` | array | None | Optional: Tissue ID(s) to filter. At least one of gencode_id, variant_id, or tissue_site_detail_id required {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `dataset_id` | string | None | 提供方的数据集标识。 {&quot;enum&quot;: &#91;&quot;gtex_v8&quot;, &quot;gtex_v10&quot;, &quot;gtex_snrnaseq_pilot&quot;&#93;, &quot;default&quot;: &quot;gtex_v8&quot;} |
| `page` | integer | None | 结果页码；起始页约定按提供方参数说明。 {&quot;minimum&quot;: 0, &quot;default&quot;: 0} |
| `items_per_page` | integer | None | 每页返回的记录数。 {&quot;minimum&quot;: 1, &quot;maximum&quot;: 100000, &quot;default&quot;: 250} |
| `gene_symbol` | string | None | 基因符号，例如 BRCA1。 |
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
            "get_single_tissue_eqtls"
          ],
          "description": "Operation type",
          "default": "get_single_tissue_eqtls"
        },
        "gencode_id": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "Optional: Versioned GENCODE ID(s) to query"
        },
        "variant_id": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "Optional: GTEx variant ID(s) to query"
        },
        "tissue_site_detail_id": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "Optional: Tissue ID(s) to filter. At least one of gencode_id, variant_id, or tissue_site_detail_id required"
        },
        "dataset_id": {
          "type": "string",
          "enum": [
            "gtex_v8",
            "gtex_v10",
            "gtex_snrnaseq_pilot"
          ],
          "default": "gtex_v8"
        },
        "page": {
          "type": "integer",
          "minimum": 0,
          "default": 0
        },
        "items_per_page": {
          "type": "integer",
          "minimum": 1,
          "maximum": 100000,
          "default": 250
        },
        "gene_symbol": {
          "type": "string"
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
              "snpId": {
                "type": "string"
              },
              "pos": {
                "type": [
                  "integer",
                  "number"
                ]
              },
              "snpIdUpper": {
                "type": "string"
              },
              "variantId": {
                "type": "string"
              },
              "geneSymbol": {
                "type": "string"
              },
              "pValue": {
                "type": "number"
              },
              "geneSymbolUpper": {
                "type": "string"
              },
              "datasetId": {
                "type": "string"
              },
              "tissueSiteDetailId": {
                "type": "string"
              },
              "ontologyId": {
                "type": "string"
              },
              "chromosome": {
                "type": "string"
              },
              "gencodeId": {
                "type": "string"
              },
              "nes": {
                "type": "number"
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
        "gene_symbol": "gencode_id"
      },
      "parameters": {
        "gene_symbol": {
          "type": "string"
        }
      },
      "conversions": {
        "gencode_id": "as_list"
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `GTEx_get_single_tissue_eqtls`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gtex_get_single_tissue_eqtls", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gtex_get_tissue_sites` {#dnakit.databases.gtex_get_tissue_sites}

Get all available GTEx tissue sites with metadata. Returns tissue IDs, names, sample counts, eGene counts, and color codes for 54 non-diseased tissue sites. Includes RNA-seq sample summaries with donor age and sex statistics. Use for: discovering available tissues, finding tissue IDs for queries, checking sample sizes, understanding tissue coverage. Essential first step before querying tissue-specific data.

```text
dnakit.databases.gtex_get_tissue_sites(*, operation=None, dataset_id=None, page=None, items_per_page=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;get_tissue_sites&quot;&#93;, &quot;default&quot;: &quot;get_tissue_sites&quot;} |
| `dataset_id` | string | None | GTEx dataset version {&quot;enum&quot;: &#91;&quot;gtex_v8&quot;, &quot;gtex_v10&quot;, &quot;gtex_snrnaseq_pilot&quot;&#93;, &quot;default&quot;: &quot;gtex_v8&quot;} |
| `page` | integer | None | Page number {&quot;minimum&quot;: 0, &quot;default&quot;: 0} |
| `items_per_page` | integer | None | Results per page {&quot;minimum&quot;: 1, &quot;maximum&quot;: 100000, &quot;default&quot;: 250} |
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
            "get_tissue_sites"
          ],
          "description": "Operation type",
          "default": "get_tissue_sites"
        },
        "dataset_id": {
          "type": "string",
          "enum": [
            "gtex_v8",
            "gtex_v10",
            "gtex_snrnaseq_pilot"
          ],
          "default": "gtex_v8",
          "description": "GTEx dataset version"
        },
        "page": {
          "type": "integer",
          "minimum": 0,
          "default": 0,
          "description": "Page number"
        },
        "items_per_page": {
          "type": "integer",
          "minimum": 1,
          "maximum": 100000,
          "default": 250,
          "description": "Results per page"
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
              "tissueSiteDetailId": {
                "type": "string"
              },
              "colorHex": {
                "type": "string"
              },
              "colorRgb": {
                "type": "string"
              },
              "datasetId": {
                "type": "string"
              },
              "eGeneCount": {
                "type": [
                  "integer",
                  "null",
                  "number"
                ]
              },
              "expressedGeneCount": {
                "type": [
                  "integer",
                  "number"
                ]
              },
              "hasEGenes": {
                "type": "boolean"
              },
              "hasSGenes": {
                "type": "boolean"
              },
              "mappedInHubmap": {
                "type": "boolean"
              },
              "eqtlSampleSummary": {
                "type": "object",
                "properties": {
                  "totalCount": {
                    "type": [
                      "integer",
                      "number"
                    ]
                  },
                  "female": {
                    "type": "object",
                    "properties": {
                      "ageMax": {
                        "type": [
                          "integer",
                          "null",
                          "number"
                        ]
                      },
                      "ageMin": {
                        "type": [
                          "integer",
                          "null",
                          "number"
                        ]
                      },
                      "ageMean": {
                        "type": [
                          "null",
                          "number"
                        ]
                      },
                      "count": {
                        "type": [
                          "integer",
                          "number"
                        ]
                      }
                    }
                  },
                  "male": {
                    "type": "object",
                    "properties": {
                      "ageMax": {
                        "type": [
                          "integer",
                          "null",
                          "number"
                        ]
                      },
                      "ageMin": {
                        "type": [
                          "integer",
                          "null",
                          "number"
                        ]
                      },
                      "ageMean": {
                        "type": [
                          "null",
                          "number"
                        ]
                      },
                      "count": {
                        "type": [
                          "integer",
                          "number"
                        ]
                      }
                    }
                  }
                }
              },
              "rnaSeqSampleSummary": {
                "type": "object",
                "properties": {
                  "totalCount": {
                    "type": [
                      "integer",
                      "number"
                    ]
                  },
                  "female": {
                    "type": "object",
                    "properties": {
                      "ageMax": {
                        "type": [
                          "integer",
                          "null",
                          "number"
                        ]
                      },
                      "ageMin": {
                        "type": [
                          "integer",
                          "null",
                          "number"
                        ]
                      },
                      "ageMean": {
                        "type": [
                          "null",
                          "number"
                        ]
                      },
                      "count": {
                        "type": [
                          "integer",
                          "number"
                        ]
                      }
                    }
                  },
                  "male": {
                    "type": "object",
                    "properties": {
                      "ageMax": {
                        "type": [
                          "integer",
                          "null",
                          "number"
                        ]
                      },
                      "ageMin": {
                        "type": [
                          "integer",
                          "null",
                          "number"
                        ]
                      },
                      "ageMean": {
                        "type": [
                          "null",
                          "number"
                        ]
                      },
                      "count": {
                        "type": [
                          "integer",
                          "number"
                        ]
                      }
                    }
                  }
                }
              },
              "sGeneCount": {
                "type": [
                  "integer",
                  "null",
                  "number"
                ]
              },
              "samplingSite": {
                "type": "string"
              },
              "tissueSite": {
                "type": "string"
              },
              "tissueSiteDetail": {
                "type": "string"
              },
              "tissueSiteDetailAbbr": {
                "type": "string"
              },
              "ontologyId": {
                "type": "string"
              },
              "ontologyIri": {
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

**执行与来源：** `network`；上游接口名 `GTEx_get_tissue_sites`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gtex_get_tissue_sites", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gtex_query_eqtl` {#dnakit.databases.gtex_query_eqtl}

Query GTEx single-tissue eQTL via /association/singleTissueEqtl

```text
dnakit.databases.gtex_query_eqtl(*, gene_symbol=None, ensembl_gene_id=None, page=None, size=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_symbol` | string | None | Gene symbol (e.g., TP53, BRCA1). Auto-resolved to GENCODE ID. |
| `ensembl_gene_id` | string | None | Ensembl gene ID, e.g., ENSG00000141510 |
| `page` | integer | None | Page number (1-based) {&quot;default&quot;: 1, &quot;minimum&quot;: 1} |
| `size` | integer | None | Page size (1–100) {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回类型为 object，具体结构由所调用操作决定；不保证统一的 data 层。提供方错误会抛出异常。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene_symbol": {
          "type": "string",
          "description": "Gene symbol (e.g., TP53, BRCA1). Auto-resolved to GENCODE ID."
        },
        "ensembl_gene_id": {
          "type": "string",
          "description": "Ensembl gene ID, e.g., ENSG00000141510"
        },
        "page": {
          "type": "integer",
          "default": 1,
          "minimum": 1,
          "description": "Page number (1-based)"
        },
        "size": {
          "type": "integer",
          "default": 10,
          "minimum": 1,
          "maximum": 100,
          "description": "Page size (1–100)"
        }
      },
      "required": []
    }
    ```

**执行与来源：** `network`；上游接口名 `GTEx_query_eqtl`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gtex_query_eqtl", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.harmonizome_get_gene` {#dnakit.databases.harmonizome_get_gene}

Get comprehensive gene information from Harmonizome (Ma'ayan Lab), which integrates data from 100+ genomics datasets covering gene expression, protein interactions, pathways, diseases, drug targets, and more. Returns gene symbol, name, NCBI Entrez ID, description, synonyms, and linked protein identifiers. Useful as a quick gene summary with curated descriptions from RefSeq. Example: TP53 returns description about tumor suppressor protein, NCBI ID 7157, synonyms [TRP53, LFS1, P53, BMFS5, BCC7], and linked UniProt protein P53_HUMAN.

```text
dnakit.databases.harmonizome_get_gene(*, gene_symbol, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 必填（输入 Schema） | Official gene symbol (case-insensitive). Examples: &#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;EGFR&#x27;, &#x27;KRAS&#x27;, &#x27;CDK2&#x27;. |
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
          "description": "Official gene symbol (case-insensitive). Examples: 'TP53', 'BRCA1', 'EGFR', 'KRAS', 'CDK2'."
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
        "symbol": {
          "type": "string"
        },
        "name": {
          "type": "string"
        },
        "ncbi_entrez_gene_id": {
          "type": "integer"
        },
        "ncbi_entrez_gene_url": {
          "type": "string"
        },
        "description": {
          "type": [
            "string",
            "null"
          ]
        },
        "synonyms": {
          "type": "array",
          "items": {
            "type": "string"
          }
        },
        "proteins": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "symbol": {
                "type": "string"
              },
              "href": {
                "type": "string"
              }
            }
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `Harmonizome_get_gene`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("harmonizome_get_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.hpa_get_cancer_prognostics_by_gene` {#dnakit.databases.hpa_get_cancer_prognostics_by_gene}

Retrieve prognostic value of a gene across various cancer types, indicating if its expression level correlates with patient survival outcomes. Reports both the significant results ('prognostic_summary') and the cancer types that were assessed but found NOT prognostic ('non_prognostic_summary'), plus 'cancers_assessed_count' as the denominator, so that a negative result can be distinguished from missing data.

```text
dnakit.databases.hpa_get_cancer_prognostics_by_gene(*, ensembl_id=None, gene_name=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `ensembl_id` | string | None | Ensembl Gene ID of the gene to check, e.g., &#x27;ENSG00000141510&#x27; for TP53, &#x27;ENSG00000012048&#x27; for BRCA1. |
| `gene_name` | string | None | DNAKit optional alias for ensembl_id |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "ensembl_id": {
          "type": "string",
          "description": "Ensembl Gene ID of the gene to check, e.g., 'ENSG00000141510' for TP53, 'ENSG00000012048' for BRCA1."
        },
        "gene_name": {
          "type": "string",
          "description": "DNAKit optional alias for ensembl_id"
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
        "ensembl_id": {
          "type": "string"
        },
        "gene": {
          "type": "string"
        },
        "gene_synonym": {
          "type": "array",
          "items": {
            "type": "string"
          }
        },
        "cancers_assessed_count": {
          "type": "integer",
          "description": "Total number of cancer types HPA assessed for this gene (the denominator). Equals prognostic_cancers_count + non_prognostic_cancers_count. A cancer type absent from both summary arrays was not assessed by HPA for this gene."
        },
        "prognostic_cancers_count": {
          "type": "integer",
          "description": "Number of assessed cancer types where expression significantly correlates with survival. This is a filtered subset of cancers_assessed_count, not a total."
        },
        "non_prognostic_cancers_count": {
          "type": "integer",
          "description": "Number of assessed cancer types where HPA found NO significant survival correlation."
        },
        "prognostic_summary": {
          "type": "array",
          "description": "Cancer types with SIGNIFICANT prognostic value only. Compare its length against cancers_assessed_count before concluding a cancer type has no data — non-significant assessed cancers are in non_prognostic_summary.",
          "items": {
            "type": "object",
            "properties": {
              "cancer_type": {
                "type": "string"
              },
              "prognostic_type": {
                "type": "string"
              },
              "p_value": {
                "type": "string"
              },
              "is_prognostic": {
                "type": "boolean"
              }
            }
          }
        },
        "non_prognostic_summary": {
          "type": "array",
          "description": "Cancer types HPA assessed and found NOT prognostic, with their p-values. These are citable negative results and are distinct from cancer types HPA never assessed (absent from both arrays).",
          "items": {
            "type": "object",
            "properties": {
              "cancer_type": {
                "type": "string"
              },
              "prognostic_type": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "p_value": {
                "type": "string"
              },
              "is_prognostic": {
                "type": "boolean"
              }
            }
          }
        },
        "note": {
          "type": "string"
        }
      }
    }
    ```

??? info "DNAKit 参数扩展规则"

    ```json
    {
      "aliases": {
        "gene_name": "ensembl_id"
      },
      "conversions": {
        "ensembl_id": "ensembl_gene"
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `HPA_get_cancer_prognostics_by_gene`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("hpa_get_cancer_prognostics_by_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.hpa_get_protein_interactions_by_gene` {#dnakit.databases.hpa_get_protein_interactions_by_gene}

DEPRECATED / non-functional: the Human Protein Atlas search API no longer exposes protein-protein interaction (ppi) data, so this tool returns an error. For PPI partners use STRING_get_protein_interactions (by gene symbol, parameter 'protein_ids') or EBIProteins_get_interactions (by UniProt accession) instead.

```text
dnakit.databases.hpa_get_protein_interactions_by_gene(*, gene_name, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_name` | string | 必填（输入 Schema） | Official gene symbol, e.g., &#x27;EGFR&#x27;, &#x27;TP53&#x27;, &#x27;BRCA1&#x27;, etc. |
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
          "description": "Official gene symbol, e.g., 'EGFR', 'TP53', 'BRCA1', etc."
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
        "gene": {
          "type": "string"
        },
        "gene_synonym": {
          "type": "array",
          "items": {
            "type": "string"
          }
        },
        "interactions": {
          "type": "string"
        },
        "interactor_count": {
          "type": "integer"
        },
        "interactors": {
          "type": "array"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `HPA_get_protein_interactions_by_gene`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("hpa_get_protein_interactions_by_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.hpa_get_rna_expression_by_source` {#dnakit.databases.hpa_get_rna_expression_by_source}

Get RNA expression for a gene in a specific biological source, read from HPA's dedicated per-source download columns (t_RNA_*, blood_RNA_*, brain_RNA_*, sc_RNA_*). Supports source_type 'tissue', 'blood', 'brain' and 'single_cell'. Units differ by source: single_cell values are nCPM, all others nTPM -- the returned expression_unit is read from HPA's own column header, never assumed.

```text
dnakit.databases.hpa_get_rna_expression_by_source(*, gene_name, source_type, source_name, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_name` | string | 必填（输入 Schema） | Gene name or gene symbol, e.g., &#x27;GFAP&#x27;, &#x27;TP53&#x27;, &#x27;BRCA1&#x27;, etc. |
| `source_type` | string | 必填（输入 Schema） | The type of biological source. Choose from: &#x27;tissue&#x27; (51 consensus tissues, nTPM), &#x27;blood&#x27; (19 immune cell subsets, nTPM), &#x27;brain&#x27; (13 brain regions, nTPM), &#x27;single_cell&#x27; (~160 single cell types, nCPM). |
| `source_name` | string | 必填（输入 Schema） | The specific source name, e.g. &#x27;liver&#x27;, &#x27;heart_muscle&#x27;, &#x27;skin&#x27;, &#x27;t_cell&#x27;, &#x27;hepatocyte&#x27;, &#x27;cerebellum&#x27;. Case-insensitive; spaces and hyphens are treated as underscores. Names HPA publishes only as subsets (blood &#x27;t_cell&#x27;, &#x27;b_cell&#x27;, &#x27;monocyte&#x27;, &#x27;dendritic_cell&#x27;) resolve to a representative subset column and the response reports which one in a &#x27;note&#x27;. An invalid name returns an error listing valid names -- it never returns a silent &#x27;N/A&#x27;. |
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
          "description": "Gene name or gene symbol, e.g., 'GFAP', 'TP53', 'BRCA1', etc."
        },
        "source_type": {
          "type": "string",
          "description": "The type of biological source. Choose from: 'tissue' (51 consensus tissues, nTPM), 'blood' (19 immune cell subsets, nTPM), 'brain' (13 brain regions, nTPM), 'single_cell' (~160 single cell types, nCPM)."
        },
        "source_name": {
          "type": "string",
          "description": "The specific source name, e.g. 'liver', 'heart_muscle', 'skin', 't_cell', 'hepatocyte', 'cerebellum'. Case-insensitive; spaces and hyphens are treated as underscores. Names HPA publishes only as subsets (blood 't_cell', 'b_cell', 'monocyte', 'dendritic_cell') resolve to a representative subset column and the response reports which one in a 'note'. An invalid name returns an error listing valid names -- it never returns a silent 'N/A'."
        }
      },
      "required": [
        "gene_name",
        "source_name",
        "source_type"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "error": {
          "type": "string"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `HPA_get_rna_expression_by_source`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("hpa_get_rna_expression_by_source", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.hpa_get_rna_expression_in_specific_tissues` {#dnakit.databases.hpa_get_rna_expression_in_specific_tissues}

Query RNA expression levels (nTPM) for a specific gene in one or more user-specified tissues. Reads HPA's full per-tissue nTPM panel (the 'Tissue RNA - <tissue> [nTPM]' columns), so any of HPA's 51 published consensus tissues and brain regions can be queried, e.g. 'cerebral cortex', 'basal ganglia', 'midbrain', 'liver'. Common names are mapped to HPA's actual column names, including the ones HPA suffixes with '1' ('skin' -> 'Skin 1', 'stomach' -> 'Stomach 1', 'endometrium' -> 'Endometrium 1'); when a name resolves to a differently-named HPA column the response says so in a 'note'. A tissue reported as 'Not found' means HPA publishes no RNA column for that name (it is a naming problem or a tissue HPA covers only with protein/IHC data), not that the gene lacks expression there.

```text
dnakit.databases.hpa_get_rna_expression_in_specific_tissues(*, ensembl_id=None, tissue_names=None, tissues=None, gene_name=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `ensembl_id` | string | None | Ensembl Gene ID for the gene, e.g., &#x27;ENSG00000141510&#x27; for TP53. |
| `tissue_names` | array | None | List of tissue names to query, e.g., &#91;&#x27;skin&#x27;, &#x27;liver&#x27;, &#x27;heart muscle&#x27;, &#x27;kidney&#x27;&#93;. Case-insensitive; underscores, spaces and hyphens are interchangeable. HPA has no RNA column for &#x27;bronchus&#x27;, &#x27;nasopharynx&#x27;, &#x27;oral mucosa&#x27; or &#x27;soft tissue&#x27; (protein/IHC only). {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `tissues` | array | None | DNAKit optional alias for tissue_names {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `gene_name` | string | None | DNAKit optional alias for ensembl_id |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "ensembl_id": {
          "type": "string",
          "description": "Ensembl Gene ID for the gene, e.g., 'ENSG00000141510' for TP53."
        },
        "tissue_names": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "List of tissue names to query, e.g., ['skin', 'liver', 'heart muscle', 'kidney']. Case-insensitive; underscores, spaces and hyphens are interchangeable. HPA has no RNA column for 'bronchus', 'nasopharynx', 'oral mucosa' or 'soft tissue' (protein/IHC only)."
        },
        "tissues": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "DNAKit optional alias for tissue_names"
        },
        "gene_name": {
          "type": "string",
          "description": "DNAKit optional alias for ensembl_id"
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
        "error": {
          "type": "string"
        }
      }
    }
    ```

??? info "DNAKit 参数扩展规则"

    ```json
    {
      "aliases": {
        "tissues": "tissue_names",
        "gene_name": "ensembl_id"
      },
      "conversions": {
        "ensembl_id": "ensembl_gene"
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `HPA_get_rna_expression_in_specific_tissues`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("hpa_get_rna_expression_in_specific_tissues", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.hpa_get_subcellular_location` {#dnakit.databases.hpa_get_subcellular_location}

Get annotated subcellular locations for a protein using optimized columns parameter. Retrieves both main and additional subcellular locations efficiently.

```text
dnakit.databases.hpa_get_subcellular_location(*, gene_name, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_name` | string | 必填（输入 Schema） | Gene name or gene symbol, e.g., &#x27;CCNB1&#x27;, &#x27;TP53&#x27;, &#x27;EGFR&#x27;, etc. |
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
          "description": "Gene name or gene symbol, e.g., 'CCNB1', 'TP53', 'EGFR', etc."
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
        "gene_synonym": {
          "type": "array",
          "items": {
            "type": "string"
          }
        },
        "main_locations": {
          "type": "array",
          "items": {
            "type": "string"
          }
        },
        "additional_locations": {
          "type": "array",
          "items": {
            "type": "string"
          }
        },
        "total_locations": {
          "type": "integer"
        },
        "location_summary": {
          "type": "string"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `HPA_get_subcellular_location`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("hpa_get_subcellular_location", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.hpa_search_genes_by_query` {#dnakit.databases.hpa_search_genes_by_query}

Search for matching genes by gene name, keywords, or cell line names and return Ensembl ID list. This is the entry tool for many HPA query workflows.

```text
dnakit.databases.hpa_search_genes_by_query(*, search_query, max_results=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `search_query` | string | 必填（输入 Schema） | Gene name, alias, keyword, or cell line name to search for, e.g., &#x27;EGFR&#x27;, &#x27;TP53&#x27;, or &#x27;MCF7&#x27;. |
| `max_results` | integer &#124; null | None | Maximum number of gene matches to return (default: 50). This search does a broad full-text match with no server-side limit, so a short/common query can otherwise return thousands of results; exact gene-symbol matches are ranked first. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "search_query": {
          "type": "string",
          "description": "Gene name, alias, keyword, or cell line name to search for, e.g., 'EGFR', 'TP53', or 'MCF7'."
        },
        "max_results": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Maximum number of gene matches to return (default: 50). This search does a broad full-text match with no server-side limit, so a short/common query can otherwise return thousands of results; exact gene-symbol matches are ranked first."
        }
      },
      "required": [
        "search_query"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "error": {
          "type": "string"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `HPA_search_genes_by_query`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("hpa_search_genes_by_query", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.jaspar_get_matrix` {#dnakit.databases.jaspar_get_matrix}

Get full details for a single JASPAR matrix by `matrix_id` (e.g., MA0002.2). Returns metadata plus PFM/PWM. For reproducible analyses, use `jaspar_get_matrix_versions` with the corresponding `base_id` (e.g., MA0002) to see all versions and choose/cite the exact version you used.

```text
dnakit.databases.jaspar_get_matrix(*, matrix_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `matrix_id` | string | 必填（输入 Schema） | JASPAR matrix ID (e.g., &#x27;MA0002.2&#x27;). You can discover IDs with `JASPAR_get_transcription_factors`. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "matrix_id": {
          "type": "string",
          "description": "JASPAR matrix ID (e.g., 'MA0002.2'). You can discover IDs with `JASPAR_get_transcription_factors`."
        }
      },
      "required": [
        "matrix_id"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "description": "Full matrix detail including motif counts (PFM) and metadata",
      "properties": {
        "matrix_id": {
          "type": "string"
        },
        "base_id": {
          "type": "string"
        },
        "name": {
          "type": "string"
        },
        "collection": {
          "type": "string"
        },
        "version": {
          "type": [
            "string",
            "integer"
          ]
        },
        "pfm": {
          "type": "object",
          "description": "Position frequency matrix",
          "properties": {
            "A": {
              "type": "array",
              "items": {
                "type": "number"
              }
            },
            "C": {
              "type": "array",
              "items": {
                "type": "number"
              }
            },
            "G": {
              "type": "array",
              "items": {
                "type": "number"
              }
            },
            "T": {
              "type": "array",
              "items": {
                "type": "number"
              }
            }
          }
        },
        "tax_group": {
          "type": [
            "string",
            "null"
          ]
        },
        "species": {
          "type": [
            "array",
            "null"
          ],
          "items": {
            "type": "object",
            "properties": {
              "tax_id": {
                "type": [
                  "integer",
                  "string"
                ]
              },
              "name": {
                "type": "string"
              }
            }
          }
        },
        "class": {
          "type": [
            "array",
            "null"
          ],
          "items": {
            "type": "string"
          }
        },
        "sequence_logo": {
          "type": [
            "string",
            "null"
          ]
        },
        "versions_url": {
          "type": [
            "string",
            "null"
          ]
        },
        "sites_url": {
          "type": [
            "string",
            "null"
          ]
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `jaspar_get_matrix`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("jaspar_get_matrix", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.jaspar_get_matrix_versions` {#dnakit.databases.jaspar_get_matrix_versions}

List all versions available for a JASPAR matrix base ID (e.g., base_id MA0002). This is useful for reproducibility (choose a specific matrix version like MA0002.2) and for tracking updates across releases.

```text
dnakit.databases.jaspar_get_matrix_versions(*, base_id, page=None, page_size=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `base_id` | string | 必填（输入 Schema） | JASPAR base matrix ID without version (e.g., &#x27;MA0002&#x27;). You can get `base_id` from `jaspar_get_matrix`. |
| `page` | integer | None | Page number (1-based). {&quot;default&quot;: 1} |
| `page_size` | integer | None | Results per page. {&quot;default&quot;: 20} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "base_id": {
          "type": "string",
          "description": "JASPAR base matrix ID without version (e.g., 'MA0002'). You can get `base_id` from `jaspar_get_matrix`."
        },
        "page": {
          "type": "integer",
          "default": 1,
          "description": "Page number (1-based)."
        },
        "page_size": {
          "type": "integer",
          "default": 20,
          "description": "Results per page."
        }
      },
      "required": [
        "base_id"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "description": "Available versions for a base matrix ID",
      "properties": {
        "count": {
          "type": "integer"
        },
        "next": {
          "type": [
            "string",
            "null"
          ]
        },
        "previous": {
          "type": [
            "string",
            "null"
          ]
        },
        "results": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "matrix_id": {
                "type": "string"
              },
              "base_id": {
                "type": "string"
              },
              "name": {
                "type": "string"
              },
              "collection": {
                "type": "string"
              },
              "version": {
                "type": [
                  "string",
                  "integer"
                ]
              },
              "url": {
                "type": "string"
              },
              "sequence_logo": {
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

**执行与来源：** `network`；上游接口名 `jaspar_get_matrix_versions`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("jaspar_get_matrix_versions", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.jaspar_get_transcription_factors` {#dnakit.databases.jaspar_get_transcription_factors}

List all JASPAR transcription factor binding profiles. Returns all matrices; does NOT support filtering by TF name. To search by name, use jaspar_search_matrices with {"name": "FOXA1"}. Supports pagination via page/page_size.

```text
dnakit.databases.jaspar_get_transcription_factors(*, collection=None, page=None, page_size=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `collection` | string | None | JASPAR collection (e.g., CORE). {&quot;default&quot;: &quot;CORE&quot;} |
| `page` | integer | None | Page number (1-based). {&quot;default&quot;: 1} |
| `page_size` | integer | None | Results per page (JASPAR `page_size`). {&quot;default&quot;: 20} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "collection": {
          "type": "string",
          "default": "CORE",
          "description": "JASPAR collection (e.g., CORE)."
        },
        "page": {
          "type": "integer",
          "default": 1,
          "description": "Page number (1-based)."
        },
        "page_size": {
          "type": "integer",
          "default": 20,
          "description": "Results per page (JASPAR `page_size`)."
        }
      },
      "additionalProperties": false,
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "description": "Paginated list of JASPAR matrices (summary metadata)",
      "properties": {
        "count": {
          "type": "integer"
        },
        "next": {
          "type": [
            "string",
            "null"
          ]
        },
        "previous": {
          "type": [
            "string",
            "null"
          ]
        },
        "results": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "matrix_id": {
                "type": "string"
              },
              "base_id": {
                "type": "string"
              },
              "name": {
                "type": "string"
              },
              "collection": {
                "type": "string"
              },
              "version": {
                "type": [
                  "string",
                  "integer"
                ]
              },
              "url": {
                "type": "string"
              },
              "sequence_logo": {
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

**执行与来源：** `network`；上游接口名 `JASPAR_get_transcription_factors`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("jaspar_get_transcription_factors", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.jaspar_search_matrices` {#dnakit.databases.jaspar_search_matrices}

Search/filter JASPAR matrices using common query parameters. Typical filters include `search` (free text), `name`, `collection`, `tax_group`, `species`, and paging via `page`/`page_size`.

```text
dnakit.databases.jaspar_search_matrices(*, search=None, name=None, collection=None, tax_group=None, species=None, page=None, page_size=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `search` | string | None | Free-text search across matrices (e.g., &#x27;CTCF&#x27;). |
| `name` | string | None | Filter by TF name (e.g., &#x27;CTCF&#x27;). |
| `collection` | string | None | Filter by collection (e.g., &#x27;CORE&#x27;). |
| `tax_group` | string | None | Filter by taxonomic group (e.g., &#x27;vertebrates&#x27;). |
| `species` | string | None | Filter by NCBI taxonomy ID (e.g., &#x27;9606&#x27; for human). Sent to JASPAR as its `tax_id` query parameter. |
| `page` | integer | None | Page number (1-based). {&quot;default&quot;: 1} |
| `page_size` | integer | None | Results per page. {&quot;default&quot;: 20} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "search": {
          "type": "string",
          "description": "Free-text search across matrices (e.g., 'CTCF')."
        },
        "name": {
          "type": "string",
          "description": "Filter by TF name (e.g., 'CTCF')."
        },
        "collection": {
          "type": "string",
          "description": "Filter by collection (e.g., 'CORE')."
        },
        "tax_group": {
          "type": "string",
          "description": "Filter by taxonomic group (e.g., 'vertebrates')."
        },
        "species": {
          "type": "string",
          "description": "Filter by NCBI taxonomy ID (e.g., '9606' for human). Sent to JASPAR as its `tax_id` query parameter."
        },
        "page": {
          "type": "integer",
          "default": 1,
          "description": "Page number (1-based)."
        },
        "page_size": {
          "type": "integer",
          "default": 20,
          "description": "Results per page."
        }
      },
      "additionalProperties": false,
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "description": "Paginated list of JASPAR matrices matching the provided filters",
      "properties": {
        "count": {
          "type": "integer"
        },
        "next": {
          "type": [
            "string",
            "null"
          ]
        },
        "previous": {
          "type": [
            "string",
            "null"
          ]
        },
        "results": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "matrix_id": {
                "type": "string"
              },
              "base_id": {
                "type": "string"
              },
              "name": {
                "type": "string"
              },
              "collection": {
                "type": "string"
              },
              "version": {
                "type": [
                  "string",
                  "integer"
                ]
              },
              "url": {
                "type": "string"
              },
              "sequence_logo": {
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

**执行与来源：** `network`；上游接口名 `jaspar_search_matrices`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("jaspar_search_matrices", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.re_map_get_transcription_factor_binding` {#dnakit.databases.re_map_get_transcription_factor_binding}

Get transcription factor binding sites from ReMap database for specific genes and cell types

```text
dnakit.databases.re_map_get_transcription_factor_binding(*, gene_name, cell_type=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_name` | string | 必填（输入 Schema） | Gene symbol (e.g., BRCA1, TP53, MYC) |
| `cell_type` | string | None | Cell type (e.g., HepG2, K562, MCF-7) {&quot;default&quot;: &quot;HepG2&quot;} |
| `limit` | integer | None | Number of results to return {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 50} |
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
          "description": "Gene symbol (e.g., BRCA1, TP53, MYC)"
        },
        "cell_type": {
          "type": "string",
          "description": "Cell type (e.g., HepG2, K562, MCF-7)",
          "default": "HepG2"
        },
        "limit": {
          "type": "integer",
          "description": "Number of results to return",
          "default": 10,
          "minimum": 1,
          "maximum": 50
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
      "oneOf": [
        {
          "type": "object",
          "required": [
            "experiments",
            "count",
            "gene_name",
            "cell_type"
          ],
          "properties": {
            "experiments": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "accession": {
                    "type": "string"
                  },
                  "assay_title": {
                    "type": "string"
                  },
                  "target": {
                    "type": [
                      "object",
                      "null"
                    ]
                  },
                  "biosample_ontology": {
                    "type": [
                      "object",
                      "null"
                    ]
                  },
                  "description": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "status": {
                    "type": [
                      "string",
                      "null"
                    ]
                  }
                }
              }
            },
            "gene_name": {
              "type": "string"
            },
            "cell_type": {
              "type": "string"
            },
            "count": {
              "type": "integer"
            },
            "url": {
              "type": "string"
            },
            "note": {
              "type": "string"
            }
          }
        },
        {
          "type": "object",
          "required": [
            "error"
          ],
          "properties": {
            "error": {
              "type": "string"
            }
          }
        }
      ]
    }
    ```

**执行与来源：** `network`；上游接口名 `ReMap_get_transcription_factor_binding`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("re_map_get_transcription_factor_binding", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.regulome_db_query_variant` {#dnakit.databases.regulome_db_query_variant}

Query RegulomeDB for regulatory variant annotations using rsID

```text
dnakit.databases.regulome_db_query_variant(*, rsid=None, variant=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `rsid` | string | None | dbSNP rsID (e.g., &#x27;rs4994&#x27;, &#x27;rs429358&#x27;) |
| `variant` | string | None | Alias for rsid. dbSNP rsID (e.g., &#x27;rs7903146&#x27;, &#x27;rs429358&#x27;). |
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
          "description": "dbSNP rsID (e.g., 'rs4994', 'rs429358')"
        },
        "variant": {
          "type": "string",
          "description": "Alias for rsid. dbSNP rsID (e.g., 'rs7903146', 'rs429358')."
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
          "description": "Successful regulatory annotation",
          "properties": {
            "rsid": {
              "type": "string"
            },
            "assembly": {
              "type": [
                "string",
                "null"
              ]
            },
            "query_coordinates": {
              "type": [
                "array",
                "null"
              ],
              "items": {
                "type": "string"
              }
            },
            "regulome_score": {
              "type": [
                "object",
                "null"
              ],
              "properties": {
                "probability": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "ranking": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "tissue_specific_scores": {
                  "type": [
                    "object",
                    "null"
                  ]
                }
              }
            },
            "variants": {
              "type": "array",
              "items": {
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
                  "rsids": {
                    "type": "array",
                    "items": {
                      "type": "string"
                    }
                  }
                }
              }
            },
            "features": {
              "type": [
                "object",
                "array",
                "null"
              ]
            },
            "nearby_snps": {
              "type": "array"
            },
            "notifications": {
              "type": [
                "object",
                "null"
              ]
            },
            "total_supporting_datasets": {
              "type": "integer"
            }
          },
          "required": [
            "rsid",
            "regulome_score",
            "variants"
          ]
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

**执行与来源：** `network`；上游接口名 `RegulomeDB_query_variant`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("regulome_db_query_variant", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.screen_get_regulatory_elements` {#dnakit.databases.screen_get_regulatory_elements}

Get candidate cis-regulatory elements (cCREs) from SCREEN database for specific genomic regions

```text
dnakit.databases.screen_get_regulatory_elements(*, gene_name=None, query=None, element_type=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_name` | string | None | Gene symbol to search for regulatory elements (e.g., BRCA1, TP53). Alias: &#x27;query&#x27;. |
| `query` | string | None | Gene symbol to search for regulatory elements (alias for gene_name). |
| `element_type` | string | None | Type of regulatory element (promoter, enhancer, insulator) {&quot;default&quot;: &quot;enhancer&quot;} |
| `limit` | integer | None | Number of results to return {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 50} |
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
          "description": "Gene symbol to search for regulatory elements (e.g., BRCA1, TP53). Alias: 'query'."
        },
        "query": {
          "type": "string",
          "description": "Gene symbol to search for regulatory elements (alias for gene_name)."
        },
        "element_type": {
          "type": "string",
          "description": "Type of regulatory element (promoter, enhancer, insulator)",
          "default": "enhancer"
        },
        "limit": {
          "type": "integer",
          "description": "Number of results to return",
          "default": 10,
          "minimum": 1,
          "maximum": 50
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
            "gene_name": {
              "type": "string"
            },
            "element_type": {
              "type": "string"
            },
            "count": {
              "type": [
                "integer",
                "number"
              ]
            },
            "regulatory_elements": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "accession": {
                    "type": "string"
                  },
                  "assay_title": {
                    "type": "string"
                  },
                  "target": {
                    "type": "string"
                  },
                  "biosample": {
                    "type": "string"
                  },
                  "status": {
                    "type": "string"
                  }
                }
              }
            }
          },
          "required": [
            "gene_name"
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

**执行与来源：** `network`；上游接口名 `SCREEN_get_regulatory_elements`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("screen_get_regulatory_elements", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ucsc_get_encode_c_cres` {#dnakit.databases.ucsc_get_encode_c_cres}

Get ENCODE4 candidate cis-Regulatory Elements (cCREs) for a genomic region from UCSC Genome Browser. cCREs are classified by epigenetic signatures into: Promoter-like (PLS, red), Proximal enhancer-like (pELS, orange), Distal enhancer-like (dELS, yellow), DNase-H3K4me3 (blue), CTCF-only/CTCF-bound (green). Each cCRE includes max Z-scores for DNase, H3K4me3, H3K27ac, and CTCF signals across biosamples. Use this to identify regulatory elements in a region of interest for functional genomics or GWAS variant interpretation.

```text
dnakit.databases.ucsc_get_encode_c_cres(*, genome=None, chrom, start, end, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `genome` | string | None | Genome assembly (e.g., &#x27;hg38&#x27;, &#x27;mm10&#x27;). {&quot;default&quot;: &quot;hg38&quot;} |
| `chrom` | string | 必填（输入 Schema） | Chromosome name (e.g., &#x27;chr17&#x27;). |
| `start` | integer | 必填（输入 Schema） | Start position (0-based). |
| `end` | integer | 必填（输入 Schema） | End position (0-based). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "genome": {
          "type": "string",
          "default": "hg38",
          "description": "Genome assembly (e.g., 'hg38', 'mm10')."
        },
        "chrom": {
          "type": "string",
          "description": "Chromosome name (e.g., 'chr17')."
        },
        "start": {
          "type": "integer",
          "description": "Start position (0-based)."
        },
        "end": {
          "type": "integer",
          "description": "End position (0-based)."
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
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "genome": {
              "type": "string"
            },
            "region": {
              "type": "string"
            },
            "ccre_count": {
              "type": "integer"
            },
            "ccres": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "name": {
                    "type": "string"
                  },
                  "chrom": {
                    "type": "string"
                  },
                  "chromStart": {
                    "type": "integer"
                  },
                  "chromEnd": {
                    "type": "integer"
                  },
                  "cCRE_class": {
                    "type": "string"
                  },
                  "DNase_maxZ": {
                    "type": "number"
                  },
                  "H3K4me3_maxZ": {
                    "type": "number"
                  },
                  "H3K27ac_maxZ": {
                    "type": "number"
                  },
                  "CTCF_maxZ": {
                    "type": "number"
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

**执行与来源：** `network`；上游接口名 `UCSC_get_encode_cCREs`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ucsc_get_encode_c_cres", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.uni_bind_get_dataset` {#dnakit.databases.uni_bind_get_dataset}

Retrieve full direct transcription factor (TF)-DNA binding-site detail for one UniBind dataset by its dataset_id (e.g. 'EXP030726.neural_stem_cells.SMAD3', obtained from UniBind_search_datasets). Returns the TF name, cell line/tissue, biological condition, source identifiers (GEO/ENCODE/PubMed), the JASPAR motif IDs used, and a list of TFBS model blocks. Each TFBS model gives the JASPAR motif (jaspar_id + version), total number of predicted binding sites (total_tfbs), the PWM score threshold, the centrality distance threshold, the CentriMo adjusted enrichment p-value, and direct download URLs for the binding sites in BED format, the underlying sequences in FASTA, and a summary plot. These are motif-anchored, base-pair-resolution, experimentally supported direct binding sites - distinct from raw ChIP-seq peaks. Public, no API key.

```text
dnakit.databases.uni_bind_get_dataset(*, operation=None, dataset_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation selector. Must be &#x27;get_dataset&#x27;. {&quot;enum&quot;: &#91;&quot;get_dataset&quot;&#93;, &quot;default&quot;: &quot;get_dataset&quot;} |
| `dataset_id` | string | 必填（输入 Schema） | UniBind dataset identifier, e.g. &#x27;EXP030726.neural_stem_cells.SMAD3&#x27;. Obtain from UniBind_search_datasets (the &#x27;dataset_id&#x27; field). |
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
            "get_dataset"
          ],
          "default": "get_dataset",
          "description": "Operation selector. Must be 'get_dataset'."
        },
        "dataset_id": {
          "type": "string",
          "description": "UniBind dataset identifier, e.g. 'EXP030726.neural_stem_cells.SMAD3'. Obtain from UniBind_search_datasets (the 'dataset_id' field)."
        }
      },
      "required": [
        "dataset_id"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "required": [
            "dataset_id",
            "tf_name",
            "tfbs_models"
          ],
          "properties": {
            "dataset_id": {
              "type": "string"
            },
            "tf_name": {
              "type": [
                "string",
                "null"
              ]
            },
            "cell_line": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "biological_condition": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "identifier": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "jaspar_id": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "prediction_models": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "tfbs_models": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "prediction_model": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "jaspar_id": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "jaspar_version": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "total_tfbs": {
                    "type": [
                      "string",
                      "integer",
                      "null"
                    ]
                  },
                  "score_threshold": {
                    "type": [
                      "string",
                      "number",
                      "null"
                    ]
                  },
                  "distance_threshold": {
                    "type": [
                      "string",
                      "number",
                      "null"
                    ]
                  },
                  "adj_centrimo_pvalue": {
                    "type": [
                      "number",
                      "string",
                      "null"
                    ]
                  },
                  "bed_url": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "fasta_url": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "summary_plot_url": {
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
          "required": [
            "error"
          ],
          "properties": {
            "error": {
              "type": "string"
            }
          }
        }
      ]
    }
    ```

**执行与来源：** `network`；上游接口名 `UniBind_get_dataset`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("uni_bind_get_dataset", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.uni_bind_list_tfs` {#dnakit.databases.uni_bind_list_tfs}

List the catalog of transcription factors (TFs) profiled in the UniBind database of direct TF-DNA binding sites. Optionally filter by a case-insensitive substring of the TF gene symbol via 'search' (e.g. 'SMAD' returns SMAD2, SMAD3, SMAD4...). Returns a sorted, de-duplicated list of TF gene symbols. Use this to discover valid tf_name values for UniBind_search_datasets. Public, no API key.

```text
dnakit.databases.uni_bind_list_tfs(*, operation=None, search=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation selector. Must be &#x27;list_tfs&#x27;. {&quot;enum&quot;: &#91;&quot;list_tfs&quot;&#93;, &quot;default&quot;: &quot;list_tfs&quot;} |
| `search` | string &#124; null | None | Case-insensitive substring filter applied client-side to TF gene symbols, e.g. &#x27;SMAD&#x27;, &#x27;FOX&#x27;, &#x27;GATA&#x27;. Optional; omit to list all TFs. |
| `limit` | integer | None | Maximum number of TF names to return after filtering. {&quot;default&quot;: 200, &quot;minimum&quot;: 1, &quot;maximum&quot;: 1000} |
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
            "list_tfs"
          ],
          "default": "list_tfs",
          "description": "Operation selector. Must be 'list_tfs'."
        },
        "search": {
          "type": [
            "string",
            "null"
          ],
          "description": "Case-insensitive substring filter applied client-side to TF gene symbols, e.g. 'SMAD', 'FOX', 'GATA'. Optional; omit to list all TFs."
        },
        "limit": {
          "type": "integer",
          "default": 200,
          "minimum": 1,
          "maximum": 1000,
          "description": "Maximum number of TF names to return after filtering."
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
            "type": "string"
          }
        },
        {
          "type": "object",
          "required": [
            "error"
          ],
          "properties": {
            "error": {
              "type": "string"
            }
          }
        }
      ]
    }
    ```

**执行与来源：** `network`；上游接口名 `UniBind_list_tfs`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("uni_bind_list_tfs", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.uni_bind_search_datasets` {#dnakit.databases.uni_bind_search_datasets}

Search the UniBind catalog of curated, experimentally derived direct transcription factor (TF)-DNA binding-site datasets (one dataset = one ChIP-seq experiment for one TF in one cell type/condition). Filter by TF gene symbol (tf_name), species scientific name (species, e.g. 'Homo sapiens'), cell type/tissue (cell_line), and collection ('Robust' or 'Permissive'). All filters are optional and compose. Returns a paginated list of matching datasets with their TF name, total peak count, and dataset_id (use UniBind_get_dataset for full binding-site detail). UniBind provides motif-anchored, base-pair-resolution direct binding sites derived from ChIP-seq via the DAMO/ChIP-eat pipeline with JASPAR motifs - distinct from JASPAR/HOCOMOCO motif matrices and from ReMap raw ChIP-seq peaks. Public, no API key.

```text
dnakit.databases.uni_bind_search_datasets(*, operation=None, tf_name=None, species=None, cell_line=None, collection=None, order=None, page=None, page_size=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation selector. Must be &#x27;search_datasets&#x27;. {&quot;enum&quot;: &#91;&quot;search_datasets&quot;&#93;, &quot;default&quot;: &quot;search_datasets&quot;} |
| `tf_name` | string &#124; null | None | Transcription factor gene symbol to filter by, e.g. &#x27;CTCF&#x27;, &#x27;SMAD3&#x27;, &#x27;ESR1&#x27;. Exact match on UniBind&#x27;s TF name. Optional. |
| `species` | string &#124; null | None | Species scientific name to filter by, e.g. &#x27;Homo sapiens&#x27;, &#x27;Mus musculus&#x27;, &#x27;Drosophila melanogaster&#x27;. Optional. |
| `cell_line` | string &#124; null | None | Cell type / tissue / condition name to filter by, e.g. &#x27;neural stem cells&#x27;, &#x27;MCF7&#x27;. Optional. |
| `collection` | string &#124; null | None | UniBind collection: &#x27;Robust&#x27; (high-confidence, metadata-curated) or &#x27;Permissive&#x27; (broader). Optional. {&quot;enum&quot;: &#91;&quot;Robust&quot;, &quot;Permissive&quot;, null&#93;} |
| `order` | string &#124; null | None | Field to order results by, e.g. &#x27;tf_name&#x27; or &#x27;-tf_name&#x27; (prefix &#x27;-&#x27; for descending). Optional. |
| `page` | integer | None | 1-based page number for pagination. {&quot;default&quot;: 1, &quot;minimum&quot;: 1} |
| `page_size` | integer | None | Number of datasets per page (max 1000). {&quot;default&quot;: 25, &quot;minimum&quot;: 1, &quot;maximum&quot;: 1000} |
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
            "search_datasets"
          ],
          "default": "search_datasets",
          "description": "Operation selector. Must be 'search_datasets'."
        },
        "tf_name": {
          "type": [
            "string",
            "null"
          ],
          "description": "Transcription factor gene symbol to filter by, e.g. 'CTCF', 'SMAD3', 'ESR1'. Exact match on UniBind's TF name. Optional."
        },
        "species": {
          "type": [
            "string",
            "null"
          ],
          "description": "Species scientific name to filter by, e.g. 'Homo sapiens', 'Mus musculus', 'Drosophila melanogaster'. Optional."
        },
        "cell_line": {
          "type": [
            "string",
            "null"
          ],
          "description": "Cell type / tissue / condition name to filter by, e.g. 'neural stem cells', 'MCF7'. Optional."
        },
        "collection": {
          "type": [
            "string",
            "null"
          ],
          "enum": [
            "Robust",
            "Permissive",
            null
          ],
          "description": "UniBind collection: 'Robust' (high-confidence, metadata-curated) or 'Permissive' (broader). Optional."
        },
        "order": {
          "type": [
            "string",
            "null"
          ],
          "description": "Field to order results by, e.g. 'tf_name' or '-tf_name' (prefix '-' for descending). Optional."
        },
        "page": {
          "type": "integer",
          "default": 1,
          "minimum": 1,
          "description": "1-based page number for pagination."
        },
        "page_size": {
          "type": "integer",
          "default": 25,
          "minimum": 1,
          "maximum": 1000,
          "description": "Number of datasets per page (max 1000)."
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
              "tf_name": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "total_peaks": {
                "type": [
                  "string",
                  "integer",
                  "null"
                ]
              },
              "dataset_url": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "dataset_id": {
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
          "required": [
            "error"
          ],
          "properties": {
            "error": {
              "type": "string"
            }
          }
        }
      ]
    }
    ```

**执行与来源：** `network`；上游接口名 `UniBind_search_datasets`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("uni_bind_search_datasets", verbose=False)
print(info["input_schema"])
```

