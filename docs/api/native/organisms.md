# 模式生物、物种与微生物组

本页收录 **49 个新增公开函数**。 [完整索引](index.md) · [调用约定与完整示例](usage.md)

每个函数均列出当前 Python 签名、参数、返回值和已登记的执行要求。提供方说明及字段描述保留源码原文。`None` 默认值在生成的数据库包装中通常表示不向服务传该参数；实际必填项与组合约束以输入 Schema 为准。坐标系和基因组版本按各接口说明使用。

## 函数目录

- [`dnakit.databases.fly_base_get_gene`](#dnakit.databases.fly_base_get_gene)
- [`dnakit.databases.fly_base_get_gene_alleles`](#dnakit.databases.fly_base_get_gene_alleles)
- [`dnakit.databases.fly_base_get_gene_disease_models`](#dnakit.databases.fly_base_get_gene_disease_models)
- [`dnakit.databases.fly_base_get_gene_expression`](#dnakit.databases.fly_base_get_gene_expression)
- [`dnakit.databases.fly_base_get_gene_interactions`](#dnakit.databases.fly_base_get_gene_interactions)
- [`dnakit.databases.fly_base_get_gene_orthologs`](#dnakit.databases.fly_base_get_gene_orthologs)
- [`dnakit.databases.fly_mine_search`](#dnakit.databases.fly_mine_search)
- [`dnakit.databases.gbif_get_species`](#dnakit.databases.gbif_get_species)
- [`dnakit.databases.gbif_search_species`](#dnakit.databases.gbif_search_species)
- [`dnakit.databases.gtdb_get_species`](#dnakit.databases.gtdb_get_species)
- [`dnakit.databases.gtdb_get_taxon_info`](#dnakit.databases.gtdb_get_taxon_info)
- [`dnakit.databases.gtdb_search_genomes`](#dnakit.databases.gtdb_search_genomes)
- [`dnakit.databases.gtdb_search_taxon`](#dnakit.databases.gtdb_search_taxon)
- [`dnakit.databases.mgi_get_gene`](#dnakit.databases.mgi_get_gene)
- [`dnakit.databases.mgi_get_phenotypes`](#dnakit.databases.mgi_get_phenotypes)
- [`dnakit.databases.mgi_search_genes`](#dnakit.databases.mgi_search_genes)
- [`dnakit.databases.mgnify_get_genome`](#dnakit.databases.mgnify_get_genome)
- [`dnakit.databases.mgnify_get_go_terms`](#dnakit.databases.mgnify_get_go_terms)
- [`dnakit.databases.mgnify_get_interpro`](#dnakit.databases.mgnify_get_interpro)
- [`dnakit.databases.mgnify_get_study_detail`](#dnakit.databases.mgnify_get_study_detail)
- [`dnakit.databases.mgnify_get_taxonomy`](#dnakit.databases.mgnify_get_taxonomy)
- [`dnakit.databases.mgnify_list_analyses`](#dnakit.databases.mgnify_list_analyses)
- [`dnakit.databases.mgnify_list_biomes`](#dnakit.databases.mgnify_list_biomes)
- [`dnakit.databases.mgnify_search_genomes`](#dnakit.databases.mgnify_search_genomes)
- [`dnakit.databases.mgnify_search_studies`](#dnakit.databases.mgnify_search_studies)
- [`dnakit.databases.powo_search_plants`](#dnakit.databases.powo_search_plants)
- [`dnakit.databases.rgd_get_annotations`](#dnakit.databases.rgd_get_annotations)
- [`dnakit.databases.rgd_get_gene`](#dnakit.databases.rgd_get_gene)
- [`dnakit.databases.rgd_get_orthologs`](#dnakit.databases.rgd_get_orthologs)
- [`dnakit.databases.rgd_get_qtls_in_region`](#dnakit.databases.rgd_get_qtls_in_region)
- [`dnakit.databases.rgd_get_strain`](#dnakit.databases.rgd_get_strain)
- [`dnakit.databases.rgd_get_strain_annotations`](#dnakit.databases.rgd_get_strain_annotations)
- [`dnakit.databases.rgd_search_genes`](#dnakit.databases.rgd_search_genes)
- [`dnakit.databases.rgd_search_strains`](#dnakit.databases.rgd_search_strains)
- [`dnakit.databases.sgd_get_gene`](#dnakit.databases.sgd_get_gene)
- [`dnakit.databases.sgd_get_go_annotations`](#dnakit.databases.sgd_get_go_annotations)
- [`dnakit.databases.sgd_get_interactions`](#dnakit.databases.sgd_get_interactions)
- [`dnakit.databases.sgd_get_phenotypes`](#dnakit.databases.sgd_get_phenotypes)
- [`dnakit.databases.sgd_search`](#dnakit.databases.sgd_search)
- [`dnakit.databases.usda_plants_get_characteristics`](#dnakit.databases.usda_plants_get_characteristics)
- [`dnakit.databases.usda_plants_get_profile`](#dnakit.databases.usda_plants_get_profile)
- [`dnakit.databases.worm_base_get_expression`](#dnakit.databases.worm_base_get_expression)
- [`dnakit.databases.worm_base_get_gene`](#dnakit.databases.worm_base_get_gene)
- [`dnakit.databases.worm_base_get_phenotypes`](#dnakit.databases.worm_base_get_phenotypes)
- [`dnakit.databases.xenbase_get_gene`](#dnakit.databases.xenbase_get_gene)
- [`dnakit.databases.xenbase_search_genes`](#dnakit.databases.xenbase_search_genes)
- [`dnakit.databases.zfin_get_gene`](#dnakit.databases.zfin_get_gene)
- [`dnakit.databases.zfin_get_gene_expression`](#dnakit.databases.zfin_get_gene_expression)
- [`dnakit.databases.zfin_get_gene_phenotypes`](#dnakit.databases.zfin_get_gene_phenotypes)

## `dnakit.databases.fly_base_get_gene` {#dnakit.databases.fly_base_get_gene}

Get detailed Drosophila melanogaster (fruit fly) gene information from FlyBase via the Alliance of Genome Resources. Returns gene symbol, name, synopsis, genomic location, synonyms, and cross-references. FlyBase is the primary database for Drosophila genetics. Example: 'FB:FBgn0003996' returns the white (w) gene with its function as an ABCG2-class transporter.

```text
dnakit.databases.fly_base_get_gene(*, gene_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填（输入 Schema） | FlyBase gene ID with &#x27;FB:&#x27; prefix. Examples: &#x27;FB:FBgn0003996&#x27; (white/w), &#x27;FB:FBgn0000490&#x27; (decapentaplegic/dpp), &#x27;FB:FBgn0004647&#x27; (Notch/N), &#x27;FB:FBgn0000546&#x27; (eyeless/ey), &#x27;FB:FBgn0003731&#x27; (hedgehog/hh). |
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
          "description": "FlyBase gene ID with 'FB:' prefix. Examples: 'FB:FBgn0003996' (white/w), 'FB:FBgn0000490' (decapentaplegic/dpp), 'FB:FBgn0004647' (Notch/N), 'FB:FBgn0000546' (eyeless/ey), 'FB:FBgn0003731' (hedgehog/hh)."
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
            "id": {
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
            "name": {
              "type": [
                "string",
                "null"
              ]
            },
            "species": {
              "type": "object",
              "properties": {
                "name": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "short_name": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "taxon_id": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "data_provider": {
                  "type": [
                    "string",
                    "null"
                  ]
                }
              }
            },
            "gene_synopsis": {
              "type": [
                "string",
                "null"
              ]
            },
            "automated_gene_synopsis": {
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
            "so_term": {
              "type": [
                "string",
                "null"
              ]
            },
            "genomic_location": {
              "type": "object",
              "properties": {
                "chromosome": {
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
                "assembly": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "strand": {
                  "type": [
                    "string",
                    "null"
                  ]
                }
              }
            },
            "cross_references": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "name": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "url": {
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

**执行与来源：** `network`；上游接口名 `FlyBase_get_gene`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("fly_base_get_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.fly_base_get_gene_alleles` {#dnakit.databases.fly_base_get_gene_alleles}

Get alleles and variants for a Drosophila gene from FlyBase via the Alliance of Genome Resources. Returns allele symbols, categories (insertion, point mutation, etc.), associated phenotypes, disease annotations, and variant details including genomic coordinates. Example: 'FB:FBgn0000490' (dpp) returns alleles like dpp[MI03752] with variant location and molecular consequence.

```text
dnakit.databases.fly_base_get_gene_alleles(*, gene_id, limit=None, page=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填（输入 Schema） | FlyBase gene ID with &#x27;FB:&#x27; prefix. Examples: &#x27;FB:FBgn0003996&#x27; (white/w), &#x27;FB:FBgn0000490&#x27; (dpp), &#x27;FB:FBgn0004647&#x27; (Notch/N). |
| `limit` | integer | None | Maximum number of alleles to return (1-100). Default: 20. |
| `page` | integer | None | Page number for pagination. Default: 1. |
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
          "description": "FlyBase gene ID with 'FB:' prefix. Examples: 'FB:FBgn0003996' (white/w), 'FB:FBgn0000490' (dpp), 'FB:FBgn0004647' (Notch/N)."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of alleles to return (1-100). Default: 20."
        },
        "page": {
          "type": "integer",
          "description": "Page number for pagination. Default: 1."
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
              "symbol": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "symbol_text": {
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
              "has_disease": {
                "type": [
                  "boolean",
                  "null"
                ]
              },
              "has_phenotype": {
                "type": [
                  "boolean",
                  "null"
                ]
              },
              "variants": {
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
                    "type": {
                      "type": [
                        "string",
                        "null"
                      ]
                    },
                    "location": {
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

**执行与来源：** `network`；上游接口名 `FlyBase_get_gene_alleles`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("fly_base_get_gene_alleles", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.fly_base_get_gene_disease_models` {#dnakit.databases.fly_base_get_gene_disease_models}

Get disease models involving a Drosophila gene from the Alliance of Genome Resources. Drosophila is a powerful model organism for human disease research. Returns genotype-based disease models with associated human disease names and Disease Ontology IDs. Example: 'FB:FBgn0000490' (dpp) has 946 disease model entries including models for ALS, Parkinson's disease, and other neurodegenerative conditions.

```text
dnakit.databases.fly_base_get_gene_disease_models(*, gene_id, limit=None, page=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填（输入 Schema） | FlyBase gene ID with &#x27;FB:&#x27; prefix. Examples: &#x27;FB:FBgn0000490&#x27; (dpp), &#x27;FB:FBgn0004647&#x27; (Notch/N), &#x27;FB:FBgn0003996&#x27; (white/w). |
| `limit` | integer | None | Maximum number of disease models to return (1-100). Default: 20. |
| `page` | integer | None | Page number for pagination. Default: 1. |
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
          "description": "FlyBase gene ID with 'FB:' prefix. Examples: 'FB:FBgn0000490' (dpp), 'FB:FBgn0004647' (Notch/N), 'FB:FBgn0003996' (white/w)."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of disease models to return (1-100). Default: 20."
        },
        "page": {
          "type": "integer",
          "description": "Page number for pagination. Default: 1."
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
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "model_id": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "model_name": {
                "type": [
                  "string",
                  "null"
                ]
              },
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
              "data_provider": {
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
                    "association_type": {
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

**执行与来源：** `network`；上游接口名 `FlyBase_get_gene_disease_models`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("fly_base_get_gene_disease_models", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.fly_base_get_gene_expression` {#dnakit.databases.fly_base_get_gene_expression}

Get expression annotations for a Drosophila gene from the Alliance of Genome Resources. Returns per-annotation records of developmental stage (whenExpressed) and anatomical location (whereExpressed) for the gene, plus the total annotation count. Example: 'FB:FBgn0000490' (dpp) has 317 expression annotations such as 'embryonic stage 14' in 'abdominal transverse muscle cell'. Useful for understanding when and where a fly gene is expressed.

```text
dnakit.databases.fly_base_get_gene_expression(*, gene_id, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填（输入 Schema） | FlyBase gene ID, with or without the &#x27;FB:&#x27; prefix. Examples: &#x27;FB:FBgn0000490&#x27; (dpp), &#x27;FB:FBgn0003996&#x27; (white/w), &#x27;FB:FBgn0004647&#x27; (Notch/N). |
| `limit` | integer | None | Max annotations to return (default 50, max 100). {&quot;default&quot;: 50} |
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
          "description": "FlyBase gene ID, with or without the 'FB:' prefix. Examples: 'FB:FBgn0000490' (dpp), 'FB:FBgn0003996' (white/w), 'FB:FBgn0004647' (Notch/N)."
        },
        "limit": {
          "type": "integer",
          "description": "Max annotations to return (default 50, max 100).",
          "default": 50
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
            "total_annotations": {
              "type": "integer"
            },
            "returned": {
              "type": "integer"
            },
            "annotations": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "stage": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "location": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "data_provider": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "reference": {
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
            "status": {
              "type": "string"
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

**执行与来源：** `network`；上游接口名 `FlyBase_get_gene_expression`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("fly_base_get_gene_expression", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.fly_base_get_gene_interactions` {#dnakit.databases.fly_base_get_gene_interactions}

Get molecular or genetic interactions for a Drosophila gene from the Alliance of Genome Resources. Molecular interactions are physical protein-protein interactions. Genetic interactions are phenotype-based (suppression, enhancement, etc.). Returns interacting gene pairs with species and interaction type. Example: 'FB:FBgn0000490' (dpp) has 169 genetic interactions including genetic suppression with anchor.

```text
dnakit.databases.fly_base_get_gene_interactions(*, gene_id, interaction_type=None, limit=None, page=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填（输入 Schema） | FlyBase gene ID with &#x27;FB:&#x27; prefix. Examples: &#x27;FB:FBgn0003996&#x27; (white/w), &#x27;FB:FBgn0000490&#x27; (dpp), &#x27;FB:FBgn0004647&#x27; (Notch/N). |
| `interaction_type` | string | None | Type of interaction: &#x27;molecular&#x27; (physical PPI) or &#x27;genetic&#x27; (phenotype-based). Default: &#x27;molecular&#x27;. |
| `limit` | integer | None | Maximum number of interactions to return (1-100). Default: 20. |
| `page` | integer | None | Page number for pagination. Default: 1. |
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
          "description": "FlyBase gene ID with 'FB:' prefix. Examples: 'FB:FBgn0003996' (white/w), 'FB:FBgn0000490' (dpp), 'FB:FBgn0004647' (Notch/N)."
        },
        "interaction_type": {
          "type": "string",
          "description": "Type of interaction: 'molecular' (physical PPI) or 'genetic' (phenotype-based). Default: 'molecular'."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of interactions to return (1-100). Default: 20."
        },
        "page": {
          "type": "integer",
          "description": "Page number for pagination. Default: 1."
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
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "subject_gene_id": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "subject_symbol": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "interactor_gene_id": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "interactor_symbol": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "interactor_species": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "interaction_type": {
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

**执行与来源：** `network`；上游接口名 `FlyBase_get_gene_interactions`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("fly_base_get_gene_interactions", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.fly_base_get_gene_orthologs` {#dnakit.databases.fly_base_get_gene_orthologs}

Get orthologs of a Drosophila gene across other species from the Alliance of Genome Resources. Returns orthologous genes in human, mouse, zebrafish, worm, yeast, and other model organisms with prediction method details and best-score indicators. Essential for translational research from fly models to human disease genes. Example: 'FB:FBgn0004647' (Notch) returns orthologs including human NOTCH1-4.

```text
dnakit.databases.fly_base_get_gene_orthologs(*, gene_id, limit=None, page=None, stringency=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填（输入 Schema） | FlyBase gene ID with &#x27;FB:&#x27; prefix. Examples: &#x27;FB:FBgn0004647&#x27; (Notch/N), &#x27;FB:FBgn0003996&#x27; (white/w), &#x27;FB:FBgn0000546&#x27; (eyeless/ey). |
| `limit` | integer | None | Maximum number of orthologs to return (1-100). Default: 20. |
| `page` | integer | None | Page number for pagination. Default: 1. |
| `stringency` | string | None | Stringency filter for orthology calls: &#x27;stringent&#x27; (high confidence) or &#x27;moderate&#x27;. Default: &#x27;stringent&#x27;. |
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
          "description": "FlyBase gene ID with 'FB:' prefix. Examples: 'FB:FBgn0004647' (Notch/N), 'FB:FBgn0003996' (white/w), 'FB:FBgn0000546' (eyeless/ey)."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of orthologs to return (1-100). Default: 20."
        },
        "page": {
          "type": "integer",
          "description": "Page number for pagination. Default: 1."
        },
        "stringency": {
          "type": "string",
          "description": "Stringency filter for orthology calls: 'stringent' (high confidence) or 'moderate'. Default: 'stringent'."
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
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "subject_gene_id": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "subject_symbol": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "subject_species": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "ortholog_gene_id": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "ortholog_symbol": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "ortholog_species": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "is_best_score": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "is_best_score_reverse": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "methods": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "method_count": {
                "type": "integer"
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

**执行与来源：** `network`；上游接口名 `FlyBase_get_gene_orthologs`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("fly_base_get_gene_orthologs", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.fly_mine_search` {#dnakit.databases.fly_mine_search}

Search FlyMine, the InterMine data warehouse for Drosophila melanogaster (fruit fly) genomics data from FlyBase. Searches across genes, proteins, publications, protein domains, gene ontology terms, and more. FlyMine integrates data from FlyBase, UniProt, Reactome, KEGG, TreeFam, and other model organism databases. Returns matching records with category facets and relevance scores. Use facet_Category parameter to filter by entity type (e.g., 'Gene', 'Protein', 'Publication', 'ProteinDomain', 'CDNAClone', 'GOTerm', 'InteractionExperiment').

```text
dnakit.databases.fly_mine_search(*, q, size=None, format=None, facet_Category=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `q` | string | 必填（输入 Schema） | Search query: gene symbol (e.g., &#x27;dpp&#x27;, &#x27;Notch&#x27;, &#x27;hedgehog&#x27;), protein name, domain name, or free text. |
| `size` | integer | None | Number of results to return (default: 10, max: 100). {&quot;default&quot;: 10} |
| `format` | string | None | Response format. Must be &#x27;json&#x27;. {&quot;default&quot;: &quot;json&quot;} |
| `facet_Category` | string | None | Filter by entity category. Options: &#x27;Gene&#x27;, &#x27;Protein&#x27;, &#x27;Publication&#x27;, &#x27;ProteinDomain&#x27;, &#x27;CDNAClone&#x27;, &#x27;GOTerm&#x27;, &#x27;InteractionExperiment&#x27;, &#x27;RNAiScreen&#x27;, &#x27;UniProtFeature&#x27;. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "q": {
          "type": "string",
          "description": "Search query: gene symbol (e.g., 'dpp', 'Notch', 'hedgehog'), protein name, domain name, or free text."
        },
        "size": {
          "type": "integer",
          "default": 10,
          "description": "Number of results to return (default: 10, max: 100)."
        },
        "format": {
          "type": "string",
          "default": "json",
          "description": "Response format. Must be 'json'."
        },
        "facet_Category": {
          "type": "string",
          "description": "Filter by entity category. Options: 'Gene', 'Protein', 'Publication', 'ProteinDomain', 'CDNAClone', 'GOTerm', 'InteractionExperiment', 'RNAiScreen', 'UniProtFeature'."
        }
      },
      "required": [
        "q"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "description": "FlyMine search results",
          "properties": {
            "totalHits": {
              "type": "integer",
              "description": "Total matching records"
            },
            "facets": {
              "type": "object",
              "description": "Facet counts by category, organism, and pathway"
            },
            "results": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "id": {
                    "type": "integer"
                  },
                  "type": {
                    "type": "string",
                    "description": "Entity type (Gene, Protein, etc.)"
                  },
                  "relevance": {
                    "type": "number"
                  },
                  "fields": {
                    "type": "object",
                    "description": "Entity-specific fields. For genes: symbol, organism.name, primaryIdentifier, secondaryIdentifier. For proteins: primaryAccession, organism.name. For domains: name, primaryIdentifier (InterPro ID), shortName, type."
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

**执行与来源：** `network`；上游接口名 `FlyMine_search`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("fly_mine_search", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gbif_get_species` {#dnakit.databases.gbif_get_species}

Get detailed species information from the Global Biodiversity Information Facility (GBIF) by its taxon key. Returns the full taxonomic record including scientific name, authorship, rank, full classification (kingdom through species), taxonomic status, basionym, vernacular name, and number of descendant taxa. Use GBIF_search_species first to find taxon keys. Example: key 2435099 returns Puma concolor (puma) with full classification from Animalia to Puma concolor.

```text
dnakit.databases.gbif_get_species(*, speciesKey, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `speciesKey` | integer | 必填（输入 Schema） | GBIF species/taxon key. Examples: 2435099 (Puma concolor / puma), 2878688 (Quercus robur / English oak), 5219404 (Canis lupus / wolf), 2480528 (Homo sapiens), 2495662 (Aquila chrysaetos / golden eagle). Use GBIF_search_species to find keys. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "speciesKey": {
          "type": "integer",
          "description": "GBIF species/taxon key. Examples: 2435099 (Puma concolor / puma), 2878688 (Quercus robur / English oak), 5219404 (Canis lupus / wolf), 2480528 (Homo sapiens), 2495662 (Aquila chrysaetos / golden eagle). Use GBIF_search_species to find keys."
        }
      },
      "required": [
        "speciesKey"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "description": "GBIF species record",
          "properties": {
            "key": {
              "type": "integer",
              "description": "GBIF taxon key"
            },
            "scientificName": {
              "type": "string",
              "description": "Full scientific name with authority"
            },
            "canonicalName": {
              "type": [
                "string",
                "null"
              ],
              "description": "Scientific name without authority"
            },
            "vernacularName": {
              "type": [
                "string",
                "null"
              ],
              "description": "Common name"
            },
            "authorship": {
              "type": [
                "string",
                "null"
              ],
              "description": "Taxonomic authority"
            },
            "rank": {
              "type": "string",
              "description": "Taxonomic rank (KINGDOM, PHYLUM, CLASS, ORDER, FAMILY, GENUS, SPECIES)"
            },
            "taxonomicStatus": {
              "type": "string",
              "description": "Status: ACCEPTED, SYNONYM, DOUBTFUL, etc."
            },
            "kingdom": {
              "type": [
                "string",
                "null"
              ],
              "description": "Kingdom name"
            },
            "phylum": {
              "type": [
                "string",
                "null"
              ],
              "description": "Phylum name"
            },
            "class": {
              "type": [
                "string",
                "null"
              ],
              "description": "Class name"
            },
            "order": {
              "type": [
                "string",
                "null"
              ],
              "description": "Order name"
            },
            "family": {
              "type": [
                "string",
                "null"
              ],
              "description": "Family name"
            },
            "genus": {
              "type": [
                "string",
                "null"
              ],
              "description": "Genus name"
            },
            "species": {
              "type": [
                "string",
                "null"
              ],
              "description": "Species binomial"
            },
            "kingdomKey": {
              "type": [
                "integer",
                "null"
              ],
              "description": "Kingdom taxon key"
            },
            "phylumKey": {
              "type": [
                "integer",
                "null"
              ],
              "description": "Phylum taxon key"
            },
            "classKey": {
              "type": [
                "integer",
                "null"
              ],
              "description": "Class taxon key"
            },
            "orderKey": {
              "type": [
                "integer",
                "null"
              ],
              "description": "Order taxon key"
            },
            "familyKey": {
              "type": [
                "integer",
                "null"
              ],
              "description": "Family taxon key"
            },
            "genusKey": {
              "type": [
                "integer",
                "null"
              ],
              "description": "Genus taxon key"
            },
            "speciesKey": {
              "type": [
                "integer",
                "null"
              ],
              "description": "Species taxon key"
            },
            "parent": {
              "type": [
                "string",
                "null"
              ],
              "description": "Parent taxon name"
            },
            "parentKey": {
              "type": [
                "integer",
                "null"
              ],
              "description": "Parent taxon key"
            },
            "basionym": {
              "type": [
                "string",
                "null"
              ],
              "description": "Original name (basionym)"
            },
            "numDescendants": {
              "type": "integer",
              "description": "Number of descendant taxa"
            },
            "origin": {
              "type": [
                "string",
                "null"
              ],
              "description": "Data origin (SOURCE, DENORMED_CLASSIFICATION)"
            },
            "publishedIn": {
              "type": [
                "string",
                "null"
              ],
              "description": "Publication reference"
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

**执行与来源：** `network`；上游接口名 `GBIF_get_species`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gbif_get_species", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gbif_search_species` {#dnakit.databases.gbif_search_species}

Search species via GBIF species/search

```text
dnakit.databases.gbif_search_species(*, query, limit=None, offset=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填（输入 Schema） | Query keyword, e.g., Homo |
| `limit` | integer | None | 当前请求的返回数量上限；不是全库命中总数。 {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 300} |
| `offset` | integer | None | 分页起始偏移；对应当前查询结果。 {&quot;default&quot;: 0, &quot;minimum&quot;: 0} |
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
          "description": "Query keyword, e.g., Homo"
        },
        "limit": {
          "type": "integer",
          "default": 10,
          "minimum": 1,
          "maximum": 300
        },
        "offset": {
          "type": "integer",
          "default": 0,
          "minimum": 0
        }
      },
      "required": [
        "query"
      ]
    }
    ```

**执行与来源：** `network`；上游接口名 `GBIF_search_species`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gbif_search_species", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gtdb_get_species` {#dnakit.databases.gtdb_get_species}

Get detailed information about a GTDB species cluster, including all genomes assigned to that species. A GTDB species cluster groups genomes by average nucleotide identity (ANI). Returns the full list of genomes with their accessions, NCBI organism names, NCBI taxonomy, and whether they are the GTDB species representative or NCBI type material. Useful for exploring species diversity, finding representative genomes, and comparing GTDB vs NCBI classifications.

```text
dnakit.databases.gtdb_get_species(*, operation=None, species, max_genomes=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_species) {&quot;enum&quot;: &#91;&quot;get_species&quot;&#93;, &quot;default&quot;: &quot;get_species&quot;} |
| `species` | string | 必填（输入 Schema） | Species name WITHOUT the s__ prefix (e.g., &#x27;Escherichia coli&#x27;, &#x27;Hydrogenothermus marinus&#x27;, &#x27;Bacillus subtilis&#x27;) |
| `max_genomes` | integer | None | Maximum number of genomes to return (default: 20). Species clusters can have thousands of genomes. |
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
            "get_species"
          ],
          "description": "Operation type (fixed: get_species)",
          "default": "get_species"
        },
        "species": {
          "type": "string",
          "description": "Species name WITHOUT the s__ prefix (e.g., 'Escherichia coli', 'Hydrogenothermus marinus', 'Bacillus subtilis')"
        },
        "max_genomes": {
          "type": "integer",
          "description": "Maximum number of genomes to return (default: 20). Species clusters can have thousands of genomes."
        }
      },
      "required": [
        "species"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "description": "Species cluster details",
          "properties": {
            "species_name": {
              "type": "string",
              "description": "The GTDB species name"
            },
            "total_genomes": {
              "type": "integer",
              "description": "Total number of genomes in this species cluster"
            },
            "genomes_shown": {
              "type": "integer",
              "description": "Number of genomes in the response (may be limited)"
            },
            "genomes": {
              "type": "array",
              "description": "Genome entries in the species cluster",
              "items": {
                "type": "object",
                "properties": {
                  "accession": {
                    "type": "string",
                    "description": "NCBI genome accession (GCA/GCF prefix)"
                  },
                  "ncbi_org_name": {
                    "type": "string",
                    "description": "NCBI organism name"
                  },
                  "ncbi_tax": {
                    "type": "string",
                    "description": "NCBI taxonomy string"
                  },
                  "gtdb_species_rep": {
                    "type": "boolean",
                    "description": "Whether this is the GTDB species representative genome"
                  },
                  "ncbi_type_material": {
                    "type": "string",
                    "description": "NCBI type material designation"
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
              "type": "string",
              "description": "Error message"
            }
          },
          "required": [
            "error"
          ]
        }
      ]
    }
    ```

**执行与来源：** `network`；上游接口名 `GTDB_get_species`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gtdb_get_species", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gtdb_get_taxon_info` {#dnakit.databases.gtdb_get_taxon_info}

Get information about a GTDB taxon including its rank, total genome count, higher-rank lineage, and release history. GTDB uses rank prefixes: d__ (domain), p__ (phylum), c__ (class), o__ (order), f__ (family), g__ (genus), s__ (species). For example, f__Lachnospiraceae returns its rank (Family), genome count, and lineage (Bacteria > Bacillota > Clostridia > Lachnospirales). Also returns the full lineage breakdown from domain to species. Useful for understanding where a taxon fits in the GTDB tree and how many genomes support it.

```text
dnakit.databases.gtdb_get_taxon_info(*, operation=None, taxon, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_taxon_info) {&quot;enum&quot;: &#91;&quot;get_taxon_info&quot;&#93;, &quot;default&quot;: &quot;get_taxon_info&quot;} |
| `taxon` | string | 必填（输入 Schema） | GTDB taxon name with rank prefix (e.g., &#x27;f__Lachnospiraceae&#x27;, &#x27;d__Bacteria&#x27;, &#x27;g__Escherichia&#x27;, &#x27;s__Escherichia coli&#x27;). Use GTDB_search_taxon to find valid names. |
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
            "get_taxon_info"
          ],
          "description": "Operation type (fixed: get_taxon_info)",
          "default": "get_taxon_info"
        },
        "taxon": {
          "type": "string",
          "description": "GTDB taxon name with rank prefix (e.g., 'f__Lachnospiraceae', 'd__Bacteria', 'g__Escherichia', 's__Escherichia coli'). Use GTDB_search_taxon to find valid names."
        }
      },
      "required": [
        "taxon"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "description": "Taxon card information",
          "properties": {
            "taxon": {
              "type": "string",
              "description": "The GTDB taxon name"
            },
            "rank": {
              "type": [
                "string",
                "null"
              ],
              "description": "Taxonomic rank (Domain, Phylum, Class, Order, Family, Genus, Species)"
            },
            "n_genomes": {
              "type": [
                "integer",
                "null"
              ],
              "description": "Total number of genomes classified under this taxon"
            },
            "higher_ranks": {
              "type": "array",
              "description": "Parent taxa from domain down to parent rank",
              "items": {
                "type": "string"
              }
            },
            "in_releases": {
              "type": "array",
              "description": "GTDB releases this taxon appears in",
              "items": {
                "type": "string"
              }
            },
            "lineage": {
              "type": [
                "object",
                "null"
              ],
              "description": "Full lineage from domain (d) to species (s), keys: d, p, c, o, f, g, s"
            }
          }
        },
        {
          "type": "object",
          "properties": {
            "error": {
              "type": "string",
              "description": "Error message"
            }
          },
          "required": [
            "error"
          ]
        }
      ]
    }
    ```

**执行与来源：** `network`；上游接口名 `GTDB_get_taxon_info`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gtdb_get_taxon_info", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gtdb_search_genomes` {#dnakit.databases.gtdb_search_genomes}

Search for prokaryotic genomes in GTDB by organism name. A multi-word query (e.g. a binomial like 'Mycobacterium tuberculosis') is first matched precisely, requiring EVERY term of the query to appear in the genome's name or taxonomy; 'match_type' is then 'exact'. If nothing matches the whole query, the tool falls back to GTDB's loose search, which matches ANY term of the query and therefore returns a broader taxon (often just the genus) or a different species; that case is reported as 'match_type': 'broad' with a top-level 'note', and 'total_results' then counts the broader match rather than the queried organism -- always read 'match_type' and 'total_results_scope' before quoting a count or trusting a result. Results also match against the whole taxonomic clade rather than organism-name text alone, so exact/prefix name matches are ranked first within a page. Returns genome accessions with both NCBI and GTDB taxonomy classifications, allowing comparison between the two taxonomic systems. Each result shows the genome accession (GCA/GCF), NCBI organism name, NCBI taxonomy, GTDB taxonomy, and whether the genome is a GTDB species representative. Supports pagination. Useful for finding genome accessions, exploring taxonomic reclassifications, and linking genomes to the GTDB taxonomy.

```text
dnakit.databases.gtdb_search_genomes(*, operation=None, query, page=None, items_per_page=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search_genomes) {&quot;enum&quot;: &#91;&quot;search_genomes&quot;&#93;, &quot;default&quot;: &quot;search_genomes&quot;} |
| `query` | string | 必填（输入 Schema） | Search query for organism name (e.g., &#x27;Escherichia&#x27;, &#x27;Staphylococcus aureus&#x27;, &#x27;Lactobacillus&#x27;) |
| `page` | integer | None | Page number for pagination (default: 1) |
| `items_per_page` | integer | None | Number of results per page (default: 10, max: 50) |
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
            "search_genomes"
          ],
          "description": "Operation type (fixed: search_genomes)",
          "default": "search_genomes"
        },
        "query": {
          "type": "string",
          "description": "Search query for organism name (e.g., 'Escherichia', 'Staphylococcus aureus', 'Lactobacillus')"
        },
        "page": {
          "type": "integer",
          "description": "Page number for pagination (default: 1)"
        },
        "items_per_page": {
          "type": "integer",
          "description": "Number of results per page (default: 10, max: 50)"
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
          "description": "Genome search results",
          "properties": {
            "query": {
              "type": "string",
              "description": "The search query"
            },
            "match_type": {
              "type": "string",
              "enum": [
                "exact",
                "broad",
                "none"
              ],
              "description": "How the results relate to the query. 'exact': every returned genome matches every term of the query. 'broad': nothing matched the whole query, so these results match only part of it (typically the genus) and are NOT the queried organism -- see 'note'. 'none': GTDB holds no matching genome."
            },
            "total_results": {
              "type": "integer",
              "description": "Total number of genomes in the match described by 'match_type'. When match_type is 'broad' this is the count for the broader match, NOT for the queried organism -- see 'total_results_scope'."
            },
            "total_results_scope": {
              "type": "string",
              "description": "Plain-language statement of exactly what 'total_results' counts"
            },
            "note": {
              "type": "string",
              "description": "Present when match_type is 'broad' or 'none': explains that the results are not the queried organism and that total_results must not be reported as its genome count"
            },
            "page": {
              "type": "integer",
              "description": "Current page number"
            },
            "results": {
              "type": "array",
              "description": "Genome search results",
              "items": {
                "type": "object",
                "properties": {
                  "gid": {
                    "type": "string",
                    "description": "Genome ID / accession"
                  },
                  "accession": {
                    "type": "string",
                    "description": "NCBI genome accession"
                  },
                  "ncbiOrgName": {
                    "type": "string",
                    "description": "NCBI organism name"
                  },
                  "ncbiTaxonomy": {
                    "type": "string",
                    "description": "Full NCBI taxonomy string"
                  },
                  "gtdbTaxonomy": {
                    "type": "string",
                    "description": "Full GTDB taxonomy string"
                  },
                  "isGtdbSpeciesRep": {
                    "type": "boolean",
                    "description": "Whether this is a GTDB species representative"
                  },
                  "isNcbiTypeMaterial": {
                    "type": "boolean",
                    "description": "Whether this is NCBI type material"
                  }
                }
              }
            },
            "count": {
              "type": "integer",
              "description": "Number of results on this page"
            }
          }
        },
        {
          "type": "object",
          "properties": {
            "error": {
              "type": "string",
              "description": "Error message"
            }
          },
          "required": [
            "error"
          ]
        }
      ]
    }
    ```

**执行与来源：** `network`；上游接口名 `GTDB_search_genomes`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gtdb_search_genomes", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gtdb_search_taxon` {#dnakit.databases.gtdb_search_taxon}

Search the Genome Taxonomy Database (GTDB) for prokaryotic taxa by partial name. GTDB provides a standardized genome-based taxonomy for Bacteria and Archaea, resolving polyphyletic groups in NCBI taxonomy. Returns matching taxon names with GTDB rank prefixes (d__=domain, p__=phylum, c__=class, o__=order, f__=family, g__=genus, s__=species). Search is case-insensitive and matches partial names. Use this for looking up standardized prokaryotic taxonomy names, finding correct GTDB identifiers, or exploring taxonomic groups like Lachnospiraceae, Enterobacteriaceae, or Clostridia.

```text
dnakit.databases.gtdb_search_taxon(*, operation=None, query, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search_taxon) {&quot;enum&quot;: &#91;&quot;search_taxon&quot;&#93;, &quot;default&quot;: &quot;search_taxon&quot;} |
| `query` | string | 必填（输入 Schema） | Partial taxon name to search for (e.g., &#x27;Lachnospiraceae&#x27;, &#x27;Escherichia&#x27;, &#x27;Bacillota&#x27;). Case-insensitive. |
| `limit` | integer | None | Maximum number of results per rank (default: 20, max: 100) |
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
            "search_taxon"
          ],
          "description": "Operation type (fixed: search_taxon)",
          "default": "search_taxon"
        },
        "query": {
          "type": "string",
          "description": "Partial taxon name to search for (e.g., 'Lachnospiraceae', 'Escherichia', 'Bacillota'). Case-insensitive."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of results per rank (default: 20, max: 100)"
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
          "description": "Taxon search results",
          "properties": {
            "query": {
              "type": "string",
              "description": "The search query"
            },
            "matches": {
              "type": "array",
              "description": "Matching GTDB taxon names with rank prefixes",
              "items": {
                "type": "string"
              }
            },
            "count": {
              "type": "integer",
              "description": "Number of matches"
            }
          }
        },
        {
          "type": "object",
          "properties": {
            "error": {
              "type": "string",
              "description": "Error message"
            }
          },
          "required": [
            "error"
          ]
        }
      ]
    }
    ```

**执行与来源：** `network`；上游接口名 `GTDB_search_taxon`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("gtdb_search_taxon", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.mgi_get_gene` {#dnakit.databases.mgi_get_gene}

Get detailed information about a mouse gene from MGI (Mouse Genome Informatics) by its MGI identifier. Returns gene symbol, full name, species info, gene synopsis (curated and automated), synonyms, genomic location (GRCm39 assembly), SO term, and cross-references to UniProt, NCBI, Ensembl, etc. Powered by the Alliance of Genome Resources API which integrates MGI curated data. Example: 'MGI:98834' returns Trp53 (transformation related protein 53) details.

```text
dnakit.databases.mgi_get_gene(*, gene_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填（输入 Schema） | MGI gene identifier with &#x27;MGI:&#x27; prefix. Examples: &#x27;MGI:98834&#x27; (Trp53/p53), &#x27;MGI:104537&#x27; (Brca1), &#x27;MGI:96677&#x27; (Kit), &#x27;MGI:88039&#x27; (Apoe), &#x27;MGI:88276&#x27; (Braf). |
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
          "description": "MGI gene identifier with 'MGI:' prefix. Examples: 'MGI:98834' (Trp53/p53), 'MGI:104537' (Brca1), 'MGI:96677' (Kit), 'MGI:88039' (Apoe), 'MGI:88276' (Braf)."
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
      "description": "Mouse gene details from MGI",
      "properties": {
        "id": {
          "type": [
            "string",
            "null"
          ],
          "description": "MGI identifier (e.g., 'MGI:98834')"
        },
        "symbol": {
          "type": [
            "string",
            "null"
          ],
          "description": "Official gene symbol (e.g., 'Trp53')"
        },
        "name": {
          "type": [
            "string",
            "null"
          ],
          "description": "Full gene name (e.g., 'transformation related protein 53')"
        },
        "species": {
          "type": "object",
          "properties": {
            "name": {
              "type": [
                "string",
                "null"
              ]
            },
            "short_name": {
              "type": [
                "string",
                "null"
              ]
            },
            "taxon_id": {
              "type": [
                "string",
                "null"
              ]
            },
            "data_provider": {
              "type": [
                "string",
                "null"
              ]
            }
          }
        },
        "gene_synopsis": {
          "type": [
            "string",
            "null"
          ],
          "description": "Curated gene synopsis from MGI"
        },
        "automated_gene_synopsis": {
          "type": [
            "string",
            "null"
          ],
          "description": "Automatically generated gene synopsis"
        },
        "synonyms": {
          "type": "array",
          "items": {
            "type": "string"
          }
        },
        "so_term": {
          "type": [
            "string",
            "null"
          ],
          "description": "Sequence Ontology term (e.g., 'protein_coding_gene')"
        },
        "genomic_location": {
          "type": "object",
          "properties": {
            "chromosome": {
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
            "assembly": {
              "type": [
                "string",
                "null"
              ]
            },
            "strand": {
              "type": [
                "string",
                "null"
              ]
            }
          }
        },
        "cross_references": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "name": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "url": {
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

**执行与来源：** `network`；上游接口名 `MGI_get_gene`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("mgi_get_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.mgi_get_phenotypes` {#dnakit.databases.mgi_get_phenotypes}

Get phenotype annotations for a mouse gene from MGI (Mouse Genome Informatics). Returns Mammalian Phenotype (MP) ontology annotations including phenotype statements derived from mutant allele studies. Uses Alliance of Genome Resources API which integrates MGI phenotype curation data. Useful for understanding mouse gene function through loss-of-function and gain-of-function phenotypes. Example: 'MGI:98834' returns Trp53 phenotypes including tumor susceptibility, apoptosis defects, etc.

```text
dnakit.databases.mgi_get_phenotypes(*, gene_id, limit=None, page=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填（输入 Schema） | MGI gene identifier with &#x27;MGI:&#x27; prefix. Examples: &#x27;MGI:98834&#x27; (Trp53), &#x27;MGI:104537&#x27; (Brca1), &#x27;MGI:96677&#x27; (Kit). |
| `limit` | integer | None | Number of phenotype annotations to return per page (default: 20, max: 100). {&quot;default&quot;: 20} |
| `page` | integer | None | Page number for pagination (default: 1). {&quot;default&quot;: 1} |
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
          "description": "MGI gene identifier with 'MGI:' prefix. Examples: 'MGI:98834' (Trp53), 'MGI:104537' (Brca1), 'MGI:96677' (Kit)."
        },
        "limit": {
          "type": "integer",
          "default": 20,
          "description": "Number of phenotype annotations to return per page (default: 20, max: 100)."
        },
        "page": {
          "type": "integer",
          "default": 1,
          "description": "Page number for pagination (default: 1)."
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
      "description": "List of phenotype annotations for the mouse gene",
      "items": {
        "type": "object",
        "properties": {
          "gene_symbol": {
            "type": [
              "string",
              "null"
            ],
            "description": "Gene symbol"
          },
          "gene_id": {
            "type": [
              "string",
              "null"
            ],
            "description": "MGI gene ID"
          },
          "phenotype_statement": {
            "type": [
              "string",
              "null"
            ],
            "description": "Phenotype description from MP ontology annotation"
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `MGI_get_phenotypes`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("mgi_get_phenotypes", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.mgi_search_genes` {#dnakit.databases.mgi_search_genes}

Search for mouse (Mus musculus) genes in MGI (Mouse Genome Informatics) by gene symbol, name, or keyword. Returns matching mouse gene symbols, names, and MGI identifiers. Powered by the Alliance of Genome Resources API which integrates MGI data. Use this to find mouse gene IDs before calling MGI_get_gene or MGI_get_phenotypes. Examples: 'Trp53' (tumor protein p53), 'Brca1', 'insulin', 'Notch'.

```text
dnakit.databases.mgi_search_genes(*, query, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填（输入 Schema） | Gene search query: mouse gene symbol (e.g., &#x27;Trp53&#x27;, &#x27;Brca1&#x27;, &#x27;Sox9&#x27;), gene name keyword (e.g., &#x27;insulin&#x27;, &#x27;tumor&#x27;), or MGI ID (e.g., &#x27;MGI:98834&#x27;). |
| `limit` | integer | None | Maximum number of results to return (default: 10, max: 50). {&quot;default&quot;: 10} |
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
          "description": "Gene search query: mouse gene symbol (e.g., 'Trp53', 'Brca1', 'Sox9'), gene name keyword (e.g., 'insulin', 'tumor'), or MGI ID (e.g., 'MGI:98834')."
        },
        "limit": {
          "type": "integer",
          "default": 10,
          "description": "Maximum number of results to return (default: 10, max: 50)."
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
      "description": "List of matching mouse genes",
      "items": {
        "type": "object",
        "properties": {
          "symbol": {
            "type": [
              "string",
              "null"
            ],
            "description": "Official gene symbol (e.g., 'Trp53')"
          },
          "name": {
            "type": [
              "string",
              "null"
            ],
            "description": "Gene name"
          },
          "gene_id": {
            "type": [
              "string",
              "null"
            ],
            "description": "MGI identifier (e.g., 'MGI:98834')"
          },
          "category": {
            "type": [
              "string",
              "null"
            ],
            "description": "Entity category (always 'gene')"
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `MGI_search_genes`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("mgi_search_genes", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.mgnify_get_genome` {#dnakit.databases.mgnify_get_genome}

Get detailed genome information from the MGnify genome catalog. Returns assembly quality metrics (completeness, contamination, N50), taxonomy, geographic origin, rRNA counts, protein count, pangenome statistics, and cross-references to ENA. Example: MGYG000000001 returns a high-quality Firmicutes isolate genome with 98.6% completeness.

```text
dnakit.databases.mgnify_get_genome(*, genome_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `genome_id` | string | 必填（输入 Schema） | MGnify genome accession. Examples: &#x27;MGYG000000001&#x27;, &#x27;MGYG000535630&#x27;. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "genome_id": {
          "type": "string",
          "description": "MGnify genome accession. Examples: 'MGYG000000001', 'MGYG000535630'."
        }
      },
      "required": [
        "genome_id"
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
            "genome_id": {
              "type": "string"
            },
            "type": {
              "type": "string",
              "description": "Isolate or MAG"
            },
            "taxonomy": {
              "type": "string",
              "description": "Full taxonomic lineage"
            },
            "completeness": {
              "type": "number"
            },
            "contamination": {
              "type": "number"
            },
            "length": {
              "type": "integer"
            },
            "num_proteins": {
              "type": "integer"
            },
            "gc_content": {
              "type": "number"
            },
            "n50": {
              "type": "integer"
            },
            "geographic_origin": {
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

**执行与来源：** `network`；上游接口名 `MGnify_get_genome`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("mgnify_get_genome", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.mgnify_get_go_terms` {#dnakit.databases.mgnify_get_go_terms}

Get Gene Ontology (GO) functional annotations from a MGnify metagenomics analysis. Returns GO terms with their accession, description, gene count, and category (biological_process, molecular_function, or cellular_component). Use to understand the functional capacity of a microbial community. Example: MGYA00585482 returns GO:0055114 oxidation-reduction process (count=1507).

```text
dnakit.databases.mgnify_get_go_terms(*, analysis_id, page_size=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `analysis_id` | string | 必填（输入 Schema） | MGnify analysis accession (e.g., &#x27;MGYA00585482&#x27;). Find IDs via MGnify_list_analyses. |
| `page_size` | integer | None | Number of GO terms to return per page (default 25, max 100). {&quot;default&quot;: 25} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "analysis_id": {
          "type": "string",
          "description": "MGnify analysis accession (e.g., 'MGYA00585482'). Find IDs via MGnify_list_analyses."
        },
        "page_size": {
          "type": "integer",
          "description": "Number of GO terms to return per page (default 25, max 100).",
          "default": 25
        }
      },
      "required": [
        "analysis_id"
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
          "go_id": {
            "type": "string"
          },
          "description": {
            "type": "string"
          },
          "count": {
            "type": "integer"
          },
          "category": {
            "type": "string"
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `MGnify_get_go_terms`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("mgnify_get_go_terms", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.mgnify_get_interpro` {#dnakit.databases.mgnify_get_interpro}

Get InterPro protein domain annotations from a MGnify metagenomics analysis. Returns InterPro identifiers with description and gene counts, ranked by abundance. Useful for understanding the protein family composition of a microbial community. Example: MGYA00585482 returns IPR036388 Winged helix-like DNA-binding domain (count=360).

```text
dnakit.databases.mgnify_get_interpro(*, analysis_id, page_size=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `analysis_id` | string | 必填（输入 Schema） | MGnify analysis accession (e.g., &#x27;MGYA00585482&#x27;). Find IDs via MGnify_list_analyses. |
| `page_size` | integer | None | Number of InterPro entries to return per page (default 25, max 100). {&quot;default&quot;: 25} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "analysis_id": {
          "type": "string",
          "description": "MGnify analysis accession (e.g., 'MGYA00585482'). Find IDs via MGnify_list_analyses."
        },
        "page_size": {
          "type": "integer",
          "description": "Number of InterPro entries to return per page (default 25, max 100).",
          "default": 25
        }
      },
      "required": [
        "analysis_id"
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
          "interpro_id": {
            "type": "string"
          },
          "description": {
            "type": "string"
          },
          "count": {
            "type": "integer"
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `MGnify_get_interpro`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("mgnify_get_interpro", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.mgnify_get_study_detail` {#dnakit.databases.mgnify_get_study_detail}

Get detailed information about a specific MGnify metagenomics study. Returns study name, abstract, associated BioProject, centre, public/private status, sample count, and biomes. Complementary to MGnify_search_studies which lists studies; use MGnify_list_analyses / MGnify_list_analysis_downloads for exact analyses/downloads counts (this endpoint does not provide them). Example: MGYS00002008 is a Tara Oceans marine metagenomics study.

```text
dnakit.databases.mgnify_get_study_detail(*, study_accession, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `study_accession` | string | 必填（输入 Schema） | MGnify study accession. Examples: &#x27;MGYS00002008&#x27;, &#x27;MGYS00005292&#x27;. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "study_accession": {
          "type": "string",
          "description": "MGnify study accession. Examples: 'MGYS00002008', 'MGYS00005292'."
        }
      },
      "required": [
        "study_accession"
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
            "study_id": {
              "type": "string"
            },
            "study_name": {
              "type": "string"
            },
            "study_abstract": {
              "type": "string"
            },
            "bioproject": {
              "type": "string"
            },
            "centre_name": {
              "type": "string"
            },
            "biomes": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "is_public": {
              "type": [
                "boolean",
                "null"
              ]
            },
            "samples_count": {
              "type": [
                "integer",
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

**执行与来源：** `network`；上游接口名 `MGnify_get_study_detail`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("mgnify_get_study_detail", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.mgnify_get_taxonomy` {#dnakit.databases.mgnify_get_taxonomy}

Get taxonomic composition from a MGnify metagenomics analysis. Returns organisms detected with their abundance counts, taxonomic lineage (domain/phylum/class/order/family/genus/species), and rank. Supports SSU (16S/18S) and LSU (23S/28S) rRNA marker genes. Use MGnify_list_analyses to find analysis IDs for a study. Example: MGYA00585482 returns Bacteria, Proteobacteria, Firmicutes lineages with counts.

```text
dnakit.databases.mgnify_get_taxonomy(*, analysis_id, rna_type=None, page_size=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `analysis_id` | string | 必填（输入 Schema） | MGnify analysis accession (e.g., &#x27;MGYA00585482&#x27;). Find IDs via MGnify_list_analyses. |
| `rna_type` | string | None | rRNA marker type: &#x27;ssu&#x27; for 16S/18S (default) or &#x27;lsu&#x27; for 23S/28S. {&quot;default&quot;: &quot;ssu&quot;} |
| `page_size` | integer | None | Number of taxa to return per page (default 25, max 100). {&quot;default&quot;: 25} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "analysis_id": {
          "type": "string",
          "description": "MGnify analysis accession (e.g., 'MGYA00585482'). Find IDs via MGnify_list_analyses."
        },
        "rna_type": {
          "type": "string",
          "description": "rRNA marker type: 'ssu' for 16S/18S (default) or 'lsu' for 23S/28S.",
          "default": "ssu"
        },
        "page_size": {
          "type": "integer",
          "description": "Number of taxa to return per page (default 25, max 100).",
          "default": 25
        }
      },
      "required": [
        "analysis_id"
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
          "organism_id": {
            "type": "string"
          },
          "name": {
            "type": "string"
          },
          "rank": {
            "type": "string"
          },
          "domain": {
            "type": [
              "string",
              "null"
            ]
          },
          "count": {
            "type": "integer"
          },
          "lineage": {
            "type": "string"
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `MGnify_get_taxonomy`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("mgnify_get_taxonomy", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.mgnify_list_analyses` {#dnakit.databases.mgnify_list_analyses}

List MGnify analyses via /analyses for a given study_accession

```text
dnakit.databases.mgnify_list_analyses(*, study_accession, size=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `study_accession` | string | 必填（输入 Schema） | MGnify study accession, e.g., &#x27;MGYS00000001&#x27; |
| `size` | integer | None | 当前返回页的数量上限；不自动遍历后续页。 {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回类型为 object，具体结构由所调用操作决定；不保证统一的 data 层。提供方错误会抛出异常。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "study_accession": {
          "type": "string",
          "description": "MGnify study accession, e.g., 'MGYS00000001'"
        },
        "size": {
          "type": "integer",
          "default": 10,
          "minimum": 1,
          "maximum": 100
        }
      },
      "required": [
        "study_accession"
      ]
    }
    ```

**执行与来源：** `network`；上游接口名 `MGnify_list_analyses`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("mgnify_list_analyses", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.mgnify_list_biomes` {#dnakit.databases.mgnify_list_biomes}

Browse the MGnify biome hierarchy. Returns biome identifiers, names, and sample counts. Biomes range from engineered environments (bioreactors) to natural habitats (ocean, soil, host-associated). Useful for discovering what metagenomic data is available for specific environments.

```text
dnakit.databases.mgnify_list_biomes(*, depth=None, page=None, page_size=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `depth` | integer | None | Hierarchy depth to filter by (1=root, 2=second level, etc.). |
| `page` | integer | None | Page number (default 1). |
| `page_size` | integer | None | Results per page (default 25, max 100). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "depth": {
          "type": "integer",
          "description": "Hierarchy depth to filter by (1=root, 2=second level, etc.)."
        },
        "page": {
          "type": "integer",
          "description": "Page number (default 1)."
        },
        "page_size": {
          "type": "integer",
          "description": "Results per page (default 25, max 100)."
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
              "biome_id": {
                "type": "string"
              },
              "biome_name": {
                "type": "string"
              },
              "samples_count": {
                "type": "integer"
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

**执行与来源：** `network`；上游接口名 `MGnify_list_biomes`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("mgnify_list_biomes", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.mgnify_search_genomes` {#dnakit.databases.mgnify_search_genomes}

Search the MGnify genome catalog for metagenome-assembled genomes (MAGs) and isolate genomes. Filter by taxonomy lineage. Returns quality metrics, taxonomy, and protein counts. Useful for finding reference genomes from specific microbial lineages. Taxonomy names follow GTDB, not NCBI: use 'Bacteroidota' (not 'Bacteroidetes'), 'Bacillota' (not 'Firmicutes'), 'Pseudomonadota' (not 'Proteobacteria') -- an NCBI-style name usually returns 0 hits rather than an error.

```text
dnakit.databases.mgnify_search_genomes(*, taxonomy=None, page=None, page_size=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `taxonomy` | string | None | GTDB taxonomy lineage filter, matched anywhere in the lineage. Examples: &#x27;Bacteroidota&#x27; (phylum), &#x27;Bacteroides&#x27; (genus), &#x27;Bacillota&#x27;, &#x27;Pseudomonadota&#x27;. NCBI-style phylum names such as &#x27;Bacteroidetes&#x27; or &#x27;Firmicutes&#x27; return few or no hits -- use the GTDB spelling. |
| `page` | integer | None | Page number (default 1). |
| `page_size` | integer | None | Results per page (default 25, max 100). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "taxonomy": {
          "type": "string",
          "description": "GTDB taxonomy lineage filter, matched anywhere in the lineage. Examples: 'Bacteroidota' (phylum), 'Bacteroides' (genus), 'Bacillota', 'Pseudomonadota'. NCBI-style phylum names such as 'Bacteroidetes' or 'Firmicutes' return few or no hits -- use the GTDB spelling."
        },
        "page": {
          "type": "integer",
          "description": "Page number (default 1)."
        },
        "page_size": {
          "type": "integer",
          "description": "Results per page (default 25, max 100)."
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
              "genome_id": {
                "type": "string"
              },
              "type": {
                "type": "string"
              },
              "taxonomy": {
                "type": "string"
              },
              "completeness": {
                "type": "number"
              },
              "contamination": {
                "type": "number"
              },
              "length": {
                "type": "integer"
              },
              "num_proteins": {
                "type": "integer"
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

**执行与来源：** `network`；上游接口名 `MGnify_search_genomes`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("mgnify_search_genomes", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.mgnify_search_studies` {#dnakit.databases.mgnify_search_studies}

Search MGnify studies via /studies with optional biome/search filters

```text
dnakit.databases.mgnify_search_studies(*, biome=None, search=None, size=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `biome` | string | None | Biome lineage identifier, e.g., &#x27;root:Host-associated&#x27; or &#x27;root:Host-associated:Human:Digestive system&#x27; |
| `search` | string | None | Keyword to search in study title/description |
| `size` | integer | None | 当前返回页的数量上限；不自动遍历后续页。 {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回类型为 object，具体结构由所调用操作决定；不保证统一的 data 层。提供方错误会抛出异常。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "biome": {
          "type": "string",
          "description": "Biome lineage identifier, e.g., 'root:Host-associated' or 'root:Host-associated:Human:Digestive system'"
        },
        "search": {
          "type": "string",
          "description": "Keyword to search in study title/description"
        },
        "size": {
          "type": "integer",
          "default": 10,
          "minimum": 1,
          "maximum": 100
        }
      },
      "required": []
    }
    ```

**执行与来源：** `network`；上游接口名 `MGnify_search_studies`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("mgnify_search_studies", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.powo_search_plants` {#dnakit.databases.powo_search_plants}

Search the Plants of the World Online (POWO) database maintained by the Royal Botanic Gardens, Kew. POWO is the most comprehensive global checklist of vascular plants with 1.3+ million names and 350,000+ accepted species. Search by plant name, genus, family, or common name. Returns taxonomic information including accepted names, synonyms, distribution, and images. Essential for plant taxonomy, botany, and biodiversity research.

```text
dnakit.databases.powo_search_plants(*, q, f=None, perPage=None, cursor=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `q` | string | 必填（输入 Schema） | Search query. Can be scientific name, common name, genus, or family. Examples: &#x27;rosa canina&#x27;, &#x27;orchid&#x27;, &#x27;Asteraceae&#x27;, &#x27;oak&#x27; |
| `f` | string &#124; null | None | Filter facets (comma-separated). Common filters: &#x27;species_f:true&#x27; (species only), &#x27;accepted_names:true&#x27; (only accepted names). Example: &#x27;species_f:true&#x27; |
| `perPage` | integer &#124; null | None | Results per page (default 24, max 500) |
| `cursor` | string &#124; null | None | Pagination cursor from previous response |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "q": {
          "type": "string",
          "description": "Search query. Can be scientific name, common name, genus, or family. Examples: 'rosa canina', 'orchid', 'Asteraceae', 'oak'"
        },
        "f": {
          "type": [
            "string",
            "null"
          ],
          "description": "Filter facets (comma-separated). Common filters: 'species_f:true' (species only), 'accepted_names:true' (only accepted names). Example: 'species_f:true'"
        },
        "perPage": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Results per page (default 24, max 500)"
        },
        "cursor": {
          "type": [
            "string",
            "null"
          ],
          "description": "Pagination cursor from previous response"
        }
      },
      "required": [
        "q"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "description": "Plant search results from POWO",
          "properties": {
            "totalResults": {
              "type": "integer",
              "description": "Total matching results"
            },
            "page": {
              "type": "integer"
            },
            "totalPages": {
              "type": "integer"
            },
            "perPage": {
              "type": "integer"
            },
            "cursor": {
              "type": [
                "string",
                "null"
              ],
              "description": "Cursor for next page"
            },
            "results": {
              "type": "array",
              "description": "List of plant taxa",
              "items": {
                "type": "object",
                "properties": {
                  "accepted": {
                    "type": "boolean",
                    "description": "Whether this is an accepted name"
                  },
                  "author": {
                    "type": [
                      "string",
                      "null"
                    ],
                    "description": "Taxonomic authority"
                  },
                  "kingdom": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "family": {
                    "type": [
                      "string",
                      "null"
                    ],
                    "description": "Plant family (e.g. 'Rosaceae')"
                  },
                  "name": {
                    "type": "string",
                    "description": "Scientific name (e.g. 'Rosa canina')"
                  },
                  "rank": {
                    "type": "string",
                    "description": "Taxonomic rank: Species, Genus, Family, etc."
                  },
                  "snippet": {
                    "type": [
                      "string",
                      "null"
                    ],
                    "description": "Brief description or common name"
                  },
                  "url": {
                    "type": "string",
                    "description": "POWO page URL"
                  },
                  "fqId": {
                    "type": "string",
                    "description": "Fully qualified ID for species lookup"
                  },
                  "images": {
                    "type": [
                      "array",
                      "null"
                    ],
                    "description": "Available images",
                    "items": {
                      "type": "object",
                      "properties": {
                        "fullsizeUrl": {
                          "type": "string"
                        },
                        "thumbnailUrl": {
                          "type": "string"
                        },
                        "caption": {
                          "type": "string"
                        },
                        "rights": {
                          "type": "string"
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

**执行与来源：** `network`；上游接口名 `POWO_search_plants`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("powo_search_plants", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.rgd_get_annotations` {#dnakit.databases.rgd_get_annotations}

Get disease, phenotype, pathway, and GO annotations for a rat gene from RGD. Returns curated annotations with ontology terms, evidence codes, and qualifiers. Filter by aspect: D=Disease, P=Pathway, W=Phenotype, F=Molecular Function, C=Cellular Component, E=Expression. Essential for rat disease model research and toxicogenomics.

```text
dnakit.databases.rgd_get_annotations(*, rgd_id, aspect=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `rgd_id` | string | 必填（输入 Schema） | RGD gene ID (e.g., &#x27;2218&#x27; for Brca1) |
| `aspect` | string | None | Filter by annotation type: &#x27;D&#x27; (Disease), &#x27;P&#x27; (Pathway), &#x27;W&#x27; (Phenotype), &#x27;F&#x27; (Function), &#x27;C&#x27; (Component), &#x27;E&#x27; (Expression). Leave empty for all. {&quot;default&quot;: &quot;&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "rgd_id": {
          "type": "string",
          "description": "RGD gene ID (e.g., '2218' for Brca1)"
        },
        "aspect": {
          "type": "string",
          "description": "Filter by annotation type: 'D' (Disease), 'P' (Pathway), 'W' (Phenotype), 'F' (Function), 'C' (Component), 'E' (Expression). Leave empty for all.",
          "default": ""
        }
      },
      "required": [
        "rgd_id"
      ]
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

**执行与来源：** `network`；上游接口名 `RGD_get_annotations`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("rgd_get_annotations", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.rgd_get_gene` {#dnakit.databases.rgd_get_gene}

Get detailed rat gene information from the Rat Genome Database by RGD ID. Returns gene symbol, name, description (including AGR-aggregated functional summary with disease associations, pathways, and interactions), gene type, and Ensembl cross-reference. Use when looking up specific rat genes for toxicology, pharmacology, or disease model research.

```text
dnakit.databases.rgd_get_gene(*, rgd_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `rgd_id` | string | 必填（输入 Schema） | RGD gene ID (numeric, e.g., &#x27;2218&#x27; for Brca1). Can include &#x27;RGD:&#x27; prefix. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "rgd_id": {
          "type": "string",
          "description": "RGD gene ID (numeric, e.g., '2218' for Brca1). Can include 'RGD:' prefix."
        }
      },
      "required": [
        "rgd_id"
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

**执行与来源：** `network`；上游接口名 `RGD_get_gene`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("rgd_get_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.rgd_get_orthologs` {#dnakit.databases.rgd_get_orthologs}

Get orthologs of a rat gene across species (human, mouse, chinchilla, dog, pig, etc.) from RGD. Returns ortholog gene symbols, RGD IDs, names, and species. Use for cross-species comparison in model organism genetics, finding human disease gene orthologs in rat.

```text
dnakit.databases.rgd_get_orthologs(*, rgd_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `rgd_id` | string | 必填（输入 Schema） | RGD gene ID (e.g., &#x27;2218&#x27; for Brca1) |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "rgd_id": {
          "type": "string",
          "description": "RGD gene ID (e.g., '2218' for Brca1)"
        }
      },
      "required": [
        "rgd_id"
      ]
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

**执行与来源：** `network`；上游接口名 `RGD_get_orthologs`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("rgd_get_orthologs", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.rgd_get_qtls_in_region` {#dnakit.databases.rgd_get_qtls_in_region}

Get QTLs (Quantitative Trait Loci) overlapping a genomic region from the Rat Genome Database (RGD's flagship data type). Returns QTL symbol, name, RGD ID, LOD score, p-value, and inheritance type. Requires a chromosome, start/stop coordinates, and a map_key identifying the assembly (e.g., 360 for rat rn7/GRCr8, 38 for human GRCh38, 35 for mouse GRCm39). Example: chromosome='10', start=1, stop=50000000, map_key=360 returns QTL Uae41 'Urinary albumin excretion QTL 41' lod=5.26 rgdId=7387235.

```text
dnakit.databases.rgd_get_qtls_in_region(*, chromosome, start, stop, map_key, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `chromosome` | string | 必填（输入 Schema） | Chromosome name, e.g. &#x27;10&#x27; or &#x27;chr10&#x27;. |
| `start` | integer | 必填（输入 Schema） | Region start coordinate (1-based). |
| `stop` | integer | 必填（输入 Schema） | Region stop coordinate. |
| `map_key` | integer | 必填（输入 Schema） | Assembly map key. 360=rat rn7/GRCr8, 372=rat mRatBN7.2, 38=human GRCh38, 17=human GRCh37, 35=mouse GRCm39, 18=mouse GRCm38. |
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
          "description": "Chromosome name, e.g. '10' or 'chr10'."
        },
        "start": {
          "type": "integer",
          "description": "Region start coordinate (1-based)."
        },
        "stop": {
          "type": "integer",
          "description": "Region stop coordinate."
        },
        "map_key": {
          "type": "integer",
          "description": "Assembly map key. 360=rat rn7/GRCr8, 372=rat mRatBN7.2, 38=human GRCh38, 17=human GRCh37, 35=mouse GRCm39, 18=mouse GRCm38."
        }
      },
      "required": [
        "chromosome",
        "map_key",
        "start",
        "stop"
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
              "rgd_id": {
                "type": [
                  "integer",
                  "null"
                ]
              },
              "symbol": {
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
              "chromosome": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "lod": {
                "type": [
                  "number",
                  "null"
                ]
              },
              "p_value": {
                "type": [
                  "number",
                  "null"
                ]
              },
              "variance": {
                "type": [
                  "number",
                  "null"
                ]
              },
              "inheritance_type": {
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

**执行与来源：** `network`；上游接口名 `RGD_get_qtls_in_region`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("rgd_get_qtls_in_region", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.rgd_get_strain` {#dnakit.databases.rgd_get_strain}

Get the full Rat Genome Database (RGD) record for a single rat strain by its exact strain symbol (e.g. 'SHR', 'BN', 'SS', 'WKY', 'GK'). Returns the strain RGD ID, name, full strain name, substrain, strain type, genetics formula, origin, background strain, and description. Use RGD_search_strains first if you only know a keyword.

```text
dnakit.databases.rgd_get_strain(*, symbol, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `symbol` | string | 必填（输入 Schema） | Exact rat strain symbol (case-insensitive), e.g. &#x27;SHR&#x27; (spontaneously hypertensive rat), &#x27;BN&#x27; (Brown Norway), &#x27;SS&#x27; (Dahl salt-sensitive). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "symbol": {
          "type": "string",
          "description": "Exact rat strain symbol (case-insensitive), e.g. 'SHR' (spontaneously hypertensive rat), 'BN' (Brown Norway), 'SS' (Dahl salt-sensitive)."
        }
      },
      "required": [
        "symbol"
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
            "rgd_id": {
              "type": [
                "integer",
                "null"
              ]
            },
            "symbol": {
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
            "full_strain_name": {
              "type": [
                "string",
                "null"
              ]
            },
            "substrain": {
              "type": [
                "string",
                "null"
              ]
            },
            "strain_type": {
              "type": [
                "string",
                "null"
              ]
            },
            "genetics": {
              "type": [
                "string",
                "null"
              ]
            },
            "origin": {
              "type": [
                "string",
                "null"
              ]
            },
            "source": {
              "type": [
                "string",
                "null"
              ]
            },
            "background_strain_rgd_id": {
              "type": [
                "integer",
                "null"
              ]
            },
            "research_use": {
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

**执行与来源：** `network`；上游接口名 `RGD_get_strain`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("rgd_get_strain", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.rgd_get_strain_annotations` {#dnakit.databases.rgd_get_strain_annotations}

Get curated disease, phenotype, and variant-trait ontology annotations for a rat strain from the Rat Genome Database (RGD), by strain symbol (e.g. 'SHR', 'BN', 'SS'). This exposes RGD's rat disease-MODEL annotations - for example SHR is annotated to 'Left Ventricular Hypertrophy' (disease, MODEL: spontaneous) and 'arterial blood pressure trait' (phenotype). Each annotation includes the ontology term, accession (DOID/VT/CMO/MP/RS), category, qualifier, and evidence code. Optional 'category' filter: disease, phenotype, strain, pathway, biological_process, molecular_function, cellular_component.

```text
dnakit.databases.rgd_get_strain_annotations(*, symbol, category=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `symbol` | string | 必填（输入 Schema） | Exact rat strain symbol (case-insensitive), e.g. &#x27;SHR&#x27;, &#x27;BN&#x27;, &#x27;SS&#x27;. |
| `category` | string &#124; null | None | Optional annotation category filter: &#x27;disease&#x27;, &#x27;phenotype&#x27;, &#x27;strain&#x27;, &#x27;pathway&#x27;, etc. Omit to return all categories. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "symbol": {
          "type": "string",
          "description": "Exact rat strain symbol (case-insensitive), e.g. 'SHR', 'BN', 'SS'."
        },
        "category": {
          "type": [
            "string",
            "null"
          ],
          "description": "Optional annotation category filter: 'disease', 'phenotype', 'strain', 'pathway', etc. Omit to return all categories."
        }
      },
      "required": [
        "symbol"
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
              "term": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "term_acc": {
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
              "aspect": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "qualifier": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "evidence": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "data_source": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "ref_rgd_id": {
                "type": [
                  "integer",
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

**执行与来源：** `network`；上游接口名 `RGD_get_strain_annotations`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("rgd_get_strain_annotations", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.rgd_search_genes` {#dnakit.databases.rgd_search_genes}

Search for rat genes by symbol or keyword in the Rat Genome Database. Returns matching genes with RGD IDs, symbols, names, types, and descriptions. Supports rat, human, and mouse species. Use for finding rat orthologs, toxicology studies, or rat disease model genes.

```text
dnakit.databases.rgd_search_genes(*, query, species=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填（输入 Schema） | Gene symbol or keyword (e.g., &#x27;Brca1&#x27;, &#x27;Tp53&#x27;, &#x27;insulin&#x27;) |
| `species` | string | None | Species to search: &#x27;rat&#x27; (default), &#x27;human&#x27;, &#x27;mouse&#x27; {&quot;default&quot;: &quot;rat&quot;} |
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
          "description": "Gene symbol or keyword (e.g., 'Brca1', 'Tp53', 'insulin')"
        },
        "species": {
          "type": "string",
          "description": "Species to search: 'rat' (default), 'human', 'mouse'",
          "default": "rat"
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
          "type": "array"
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

**执行与来源：** `network`；上游接口名 `RGD_search_genes`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("rgd_search_genes", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.rgd_search_strains` {#dnakit.databases.rgd_search_strains}

Search the Rat Genome Database (RGD) rat strain catalog by keyword and/or strain type. Rat inbred/congenic/transgenic strains are RGD's flagship rat-specific resource (not in Alliance of Genome Resources or MGI). Provide 'query' (a keyword matched against symbol, name, full strain name, genetics, and description) and/or 'strain_type' (one of: inbred, congenic, mutant, transgenic, recombinant_inbred, outbred, consomic, coisogenic, segregating_inbred, hybrid, wild, conplastic, advanced_intercross_line, transchromosomal). At least one of the two must be supplied. Returns matching strains with their RGD ID, symbol, name, type, genetics, and origin.

```text
dnakit.databases.rgd_search_strains(*, query=None, strain_type=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string &#124; null | None | Keyword to match against strain symbol, name, genetics, and description (case-insensitive substring). Example: &#x27;hypertensive&#x27; or &#x27;Dahl&#x27;. |
| `strain_type` | string &#124; null | None | Filter by strain type, e.g. &#x27;inbred&#x27;, &#x27;congenic&#x27;, &#x27;transgenic&#x27;, &#x27;recombinant_inbred&#x27;, &#x27;consomic&#x27;. |
| `limit` | integer &#124; null | None | Maximum number of strains to return (default 25, max 200). |
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
          "description": "Keyword to match against strain symbol, name, genetics, and description (case-insensitive substring). Example: 'hypertensive' or 'Dahl'."
        },
        "strain_type": {
          "type": [
            "string",
            "null"
          ],
          "description": "Filter by strain type, e.g. 'inbred', 'congenic', 'transgenic', 'recombinant_inbred', 'consomic'."
        },
        "limit": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Maximum number of strains to return (default 25, max 200)."
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
              "rgd_id": {
                "type": [
                  "integer",
                  "null"
                ]
              },
              "symbol": {
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
              "full_strain_name": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "substrain": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "strain_type": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "genetics": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "origin": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "source": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "background_strain_rgd_id": {
                "type": [
                  "integer",
                  "null"
                ]
              },
              "research_use": {
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

**执行与来源：** `network`；上游接口名 `RGD_search_strains`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("rgd_search_strains", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.sgd_get_gene` {#dnakit.databases.sgd_get_gene}

Get yeast gene/locus overview from SGD (Saccharomyces Genome Database) by SGD identifier. Returns gene name, systematic name, description, UniProt ID, aliases, and functional annotations. S. cerevisiae is the primary model organism for eukaryotic cell biology. Example: S000000259 returns GAL4, the galactose metabolism transcription factor.

```text
dnakit.databases.sgd_get_gene(*, sgd_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `sgd_id` | string | 必填（输入 Schema） | SGD locus identifier. Examples: &#x27;S000000259&#x27; (GAL4), &#x27;S000003219&#x27; (RMR1), &#x27;S000005564&#x27; (ACT1). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "sgd_id": {
          "type": "string",
          "description": "SGD locus identifier. Examples: 'S000000259' (GAL4), 'S000003219' (RMR1), 'S000005564' (ACT1)."
        }
      },
      "required": [
        "sgd_id"
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
            "sgd_id": {
              "type": "string",
              "description": "SGD identifier"
            },
            "display_name": {
              "type": "string",
              "description": "Gene name (e.g., GAL4)"
            },
            "gene_name": {
              "type": "string"
            },
            "systematic_name": {
              "type": "string",
              "description": "Systematic ORF name (e.g., YPL248C)"
            },
            "locus_type": {
              "type": "string"
            },
            "description": {
              "type": "string",
              "description": "Curated functional description"
            },
            "uniprot_id": {
              "type": "string"
            },
            "aliases": {
              "type": "array",
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

**执行与来源：** `network`；上游接口名 `SGD_get_gene`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("sgd_get_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.sgd_get_go_annotations` {#dnakit.databases.sgd_get_go_annotations}

Get Gene Ontology annotations for a yeast gene from SGD. Returns GO terms (biological process, molecular function, cellular component), evidence codes, qualifiers, and references. SGD GO annotations are among the most thoroughly curated in any model organism. Example: S000003219 (RMR1) has annotations for meiotic recombination and cytoplasm localization.

```text
dnakit.databases.sgd_get_go_annotations(*, sgd_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `sgd_id` | string | 必填（输入 Schema） | SGD locus identifier. Examples: &#x27;S000003219&#x27; (RMR1), &#x27;S000005564&#x27; (ACT1), &#x27;S000000259&#x27; (GAL4). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "sgd_id": {
          "type": "string",
          "description": "SGD locus identifier. Examples: 'S000003219' (RMR1), 'S000005564' (ACT1), 'S000000259' (GAL4)."
        }
      },
      "required": [
        "sgd_id"
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
                "type": "string",
                "description": "GO term identifier (e.g., GO:0007131)"
              },
              "go_term": {
                "type": "string",
                "description": "GO term name"
              },
              "go_aspect": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "GO namespace: biological_process, molecular_function, or cellular_component"
              },
              "evidence_code": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Evidence code (IDA, IMP, IGI, etc.)"
              },
              "qualifier": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "annotation_type": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "reference": {
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

**执行与来源：** `network`；上游接口名 `SGD_get_go_annotations`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("sgd_get_go_annotations", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.sgd_get_interactions` {#dnakit.databases.sgd_get_interactions}

Get genetic and physical interactions for a yeast gene from SGD. Returns interaction type (Physical/Genetic), experiment type (Two-hybrid, Affinity Capture-MS, Synthetic Lethality, etc.), interacting gene names, and references. Yeast interaction data is the most comprehensive of any organism. Example: S000003219 (RMR1) has 178 interactions.

```text
dnakit.databases.sgd_get_interactions(*, sgd_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `sgd_id` | string | 必填（输入 Schema） | SGD locus identifier. Examples: &#x27;S000003219&#x27; (RMR1), &#x27;S000005564&#x27; (ACT1), &#x27;S000000259&#x27; (GAL4). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "sgd_id": {
          "type": "string",
          "description": "SGD locus identifier. Examples: 'S000003219' (RMR1), 'S000005564' (ACT1), 'S000000259' (GAL4)."
        }
      },
      "required": [
        "sgd_id"
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
              "interaction_type": {
                "type": "string",
                "description": "Physical or Genetic"
              },
              "experiment_type": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Experimental method (Two-hybrid, Affinity Capture-MS, Synthetic Lethality, etc.)"
              },
              "bait_gene": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "hit_gene": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "annotation_type": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "reference": {
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

**执行与来源：** `network`；上游接口名 `SGD_get_interactions`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("sgd_get_interactions", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.sgd_get_phenotypes` {#dnakit.databases.sgd_get_phenotypes}

Get curated phenotype annotations for a yeast gene from SGD. Returns observable phenotypes (growth, morphology, fitness), mutant types, experimental methods, strain backgrounds, and chemical conditions. Yeast phenotypes are a key resource for functional genomics. Example: S000003219 (RMR1) has 29 phenotype annotations including competitive fitness, respiratory growth defects.

```text
dnakit.databases.sgd_get_phenotypes(*, sgd_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `sgd_id` | string | 必填（输入 Schema） | SGD locus identifier. Examples: &#x27;S000003219&#x27; (RMR1), &#x27;S000000259&#x27; (GAL4), &#x27;S000005564&#x27; (ACT1). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "sgd_id": {
          "type": "string",
          "description": "SGD locus identifier. Examples: 'S000003219' (RMR1), 'S000000259' (GAL4), 'S000005564' (ACT1)."
        }
      },
      "required": [
        "sgd_id"
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
              "phenotype": {
                "type": "string",
                "description": "Phenotype name with qualifier (e.g., 'competitive fitness: decreased')"
              },
              "mutant_type": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "experiment_type": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "strain_name": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "allele": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "chemical": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "reference": {
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

**执行与来源：** `network`；上游接口名 `SGD_get_phenotypes`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("sgd_get_phenotypes", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.sgd_search` {#dnakit.databases.sgd_search}

Search SGD across all data types (genes, GO terms, phenotypes, references, alleles, diseases, complexes, pathways). Returns matching entries with categories, descriptions, and links. Supports category filtering. Example: searching 'GAL4' finds the transcription factor gene, related references, alleles, and GO terms.

```text
dnakit.databases.sgd_search(*, query, category=None, limit=None, offset=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填（输入 Schema） | Search query string. Examples: &#x27;GAL4&#x27;, &#x27;cell cycle&#x27;, &#x27;sporulation&#x27;. |
| `category` | string | None | Filter by category. Options: &#x27;locus&#x27;, &#x27;reference&#x27;, &#x27;biological_process&#x27;, &#x27;molecular_function&#x27;, &#x27;cellular_component&#x27;, &#x27;phenotype&#x27;, &#x27;disease&#x27;, &#x27;allele&#x27;, &#x27;complex&#x27;, &#x27;pathway&#x27;, &#x27;chemical&#x27;, &#x27;strain&#x27;, &#x27;download&#x27;. |
| `limit` | integer | None | Maximum results to return (default 10, max 50). |
| `offset` | integer | None | Pagination offset (default 0). |
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
          "description": "Search query string. Examples: 'GAL4', 'cell cycle', 'sporulation'."
        },
        "category": {
          "type": "string",
          "description": "Filter by category. Options: 'locus', 'reference', 'biological_process', 'molecular_function', 'cellular_component', 'phenotype', 'disease', 'allele', 'complex', 'pathway', 'chemical', 'strain', 'download'."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum results to return (default 10, max 50)."
        },
        "offset": {
          "type": "integer",
          "description": "Pagination offset (default 0)."
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
              "name": {
                "type": "string"
              },
              "category": {
                "type": "string"
              },
              "description": {
                "type": "string"
              },
              "href": {
                "type": "string"
              },
              "aliases": {
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

**执行与来源：** `network`；上游接口名 `SGD_search`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("sgd_search", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.usda_plants_get_characteristics` {#dnakit.databases.usda_plants_get_characteristics}

Get the USDA PLANTS morphology/physiology, growth-requirement, and reproduction characteristics for a plant by its PLANTS symbol (e.g. 'ABBA'). Returns a list of characteristic name/value/category entries (e.g. active growth period, shade tolerance, drought tolerance, mature height, bloom period). Resolves the symbol to the plant id automatically. No API key required.

```text
dnakit.databases.usda_plants_get_characteristics(*, symbol, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `symbol` | string | 必填（输入 Schema） | USDA PLANTS symbol, e.g. &#x27;ABBA&#x27;, &#x27;ACRU&#x27;. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "symbol": {
          "type": "string",
          "description": "USDA PLANTS symbol, e.g. 'ABBA', 'ACRU'."
        }
      },
      "required": [
        "symbol"
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

**执行与来源：** `network`；上游接口名 `USDA_plants_get_characteristics`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("usda_plants_get_characteristics", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.usda_plants_get_profile` {#dnakit.databases.usda_plants_get_profile}

Get a USDA PLANTS Database profile for a plant by its PLANTS symbol (e.g. 'ABBA' = balsam fir): scientific name, common name, taxonomic group/rank, duration (annual/perennial), growth habits (tree/shrub/forb), and synonyms. 'native_statuses' gives region-level (e.g. 'L48', 'CAN') native/introduced status; fine-grained per-state/county FIPS distribution codes are not returned by this endpoint (fips_distribution is included for forward-compatibility but is null in the current USDA API response). Authoritative US plant taxonomy. No API key required.

```text
dnakit.databases.usda_plants_get_profile(*, symbol, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `symbol` | string | 必填（输入 Schema） | USDA PLANTS symbol, e.g. &#x27;ABBA&#x27; (balsam fir), &#x27;ACRU&#x27; (red maple), &#x27;QUAL&#x27; (white oak). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "symbol": {
          "type": "string",
          "description": "USDA PLANTS symbol, e.g. 'ABBA' (balsam fir), 'ACRU' (red maple), 'QUAL' (white oak)."
        }
      },
      "required": [
        "symbol"
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

**执行与来源：** `network`；上游接口名 `USDA_plants_get_profile`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("usda_plants_get_profile", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.worm_base_get_expression` {#dnakit.databases.worm_base_get_expression}

Get expression data for a C. elegans gene from WormBase. Returns tissues where the gene is expressed, developmental stages, subcellular localization, and expression clusters. Expression data is curated from GFP reporters, antibody staining, RNA in situ, and microarray experiments. Example: WBGene00006763 (unc-26) is expressed in 38 tissues including neurons, pharynx, and body wall muscle.

```text
dnakit.databases.worm_base_get_expression(*, gene_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填（输入 Schema） | WormBase gene identifier. Examples: &#x27;WBGene00006763&#x27; (unc-26), &#x27;WBGene00000898&#x27; (daf-2). |
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
          "description": "WormBase gene identifier. Examples: 'WBGene00006763' (unc-26), 'WBGene00000898' (daf-2)."
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
            "wormbase_id": {
              "type": "string"
            },
            "gene_name": {
              "type": [
                "string",
                "null"
              ]
            },
            "expressed_in_count": {
              "type": "integer"
            },
            "expressed_in": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "term_id": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "term_name": {
                    "type": [
                      "string",
                      "null"
                    ]
                  }
                }
              }
            },
            "expressed_during": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "term_id": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "term_name": {
                    "type": [
                      "string",
                      "null"
                    ]
                  }
                }
              }
            },
            "subcellular_localization": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "term_id": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "term_name": {
                    "type": [
                      "string",
                      "null"
                    ]
                  }
                }
              }
            },
            "expression_clusters_count": {
              "type": "integer"
            },
            "expression_clusters": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "cluster_id": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "cluster_label": {
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

**执行与来源：** `network`；上游接口名 `WormBase_get_expression`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("worm_base_get_expression", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.worm_base_get_gene` {#dnakit.databases.worm_base_get_gene}

Get detailed gene overview from WormBase for C. elegans (nematode worm) by WormBase gene ID. Returns gene name, sequence name, species, description, gene type, and status. C. elegans is one of the most important model organisms for genetics, neuroscience, and developmental biology. Example: WBGene00006763 returns unc-26 (synaptojanin, a polyphosphoinositide phosphatase required for synaptic vesicle recycling).

```text
dnakit.databases.worm_base_get_gene(*, gene_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填（输入 Schema） | WormBase gene identifier (WBGene format). Examples: &#x27;WBGene00006763&#x27; (unc-26), &#x27;WBGene00006789&#x27; (unc-52), &#x27;WBGene00000898&#x27; (daf-2). |
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
          "description": "WormBase gene identifier (WBGene format). Examples: 'WBGene00006763' (unc-26), 'WBGene00006789' (unc-52), 'WBGene00000898' (daf-2)."
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
            "wormbase_id": {
              "type": "string"
            },
            "gene_name": {
              "type": [
                "string",
                "null"
              ]
            },
            "sequence_name": {
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
            "description": {
              "type": [
                "string",
                "null"
              ]
            },
            "gene_type": {
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

**执行与来源：** `network`；上游接口名 `WormBase_get_gene`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("worm_base_get_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.worm_base_get_phenotypes` {#dnakit.databases.worm_base_get_phenotypes}

Get phenotype annotations for a C. elegans gene from WormBase. Returns observed phenotypes (from RNAi, mutation, transgene experiments), phenotypes not observed, and total counts. Phenotype data uses the WormBase Phenotype Ontology (WPO). Example: WBGene00006763 (unc-26) has 125+ observed phenotypes including locomotion defects, synaptic vesicle recycling abnormalities, and pharyngeal pumping defects.

```text
dnakit.databases.worm_base_get_phenotypes(*, gene_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填（输入 Schema） | WormBase gene identifier. Examples: &#x27;WBGene00006763&#x27; (unc-26), &#x27;WBGene00000898&#x27; (daf-2). |
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
          "description": "WormBase gene identifier. Examples: 'WBGene00006763' (unc-26), 'WBGene00000898' (daf-2)."
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
            "wormbase_id": {
              "type": "string"
            },
            "gene_name": {
              "type": [
                "string",
                "null"
              ]
            },
            "phenotype_count": {
              "type": "integer"
            },
            "phenotypes": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "phenotype_id": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "phenotype_name": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "evidence_type": {
                    "type": [
                      "string",
                      "null"
                    ]
                  }
                }
              }
            },
            "not_observed_count": {
              "type": "integer"
            },
            "phenotypes_not_observed": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "phenotype_id": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "phenotype_name": {
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

**执行与来源：** `network`；上游接口名 `WormBase_get_phenotypes`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("worm_base_get_phenotypes", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.xenbase_get_gene` {#dnakit.databases.xenbase_get_gene}

Get detailed Xenopus gene information from Xenbase via the Alliance of Genome Resources. Returns gene symbol, name, synopsis, genomic location, synonyms, cross-references, and species details. Supports both Xenopus tropicalis and Xenopus laevis genes. Gene IDs use Xenbase format (e.g., 'XB-GENE-484286' or 'Xenbase:XB-GENE-484286').

```text
dnakit.databases.xenbase_get_gene(*, gene_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填（输入 Schema） | Xenbase gene ID. Can be with or without &#x27;Xenbase:&#x27; prefix. Examples: &#x27;XB-GENE-484286&#x27; (tp53, X. tropicalis), &#x27;XB-GENE-6252168&#x27; (tp53.L, X. laevis), &#x27;XB-GENE-487192&#x27; (pax6, X. tropicalis). |
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
          "description": "Xenbase gene ID. Can be with or without 'Xenbase:' prefix. Examples: 'XB-GENE-484286' (tp53, X. tropicalis), 'XB-GENE-6252168' (tp53.L, X. laevis), 'XB-GENE-487192' (pax6, X. tropicalis)."
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
      "properties": {
        "gene_id": {
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
        "name": {
          "type": [
            "string",
            "null"
          ]
        },
        "species": {
          "type": "object",
          "properties": {
            "name": {
              "type": [
                "string",
                "null"
              ]
            },
            "short_name": {
              "type": [
                "string",
                "null"
              ]
            },
            "taxon_id": {
              "type": [
                "string",
                "null"
              ]
            },
            "data_provider": {
              "type": [
                "string",
                "null"
              ]
            }
          }
        },
        "synonyms": {
          "type": "array",
          "items": {
            "type": "string"
          }
        },
        "gene_synopsis": {
          "type": [
            "string",
            "null"
          ]
        },
        "automated_gene_synopsis": {
          "type": [
            "string",
            "null"
          ]
        },
        "so_term": {
          "type": [
            "string",
            "null"
          ]
        },
        "genomic_location": {
          "type": [
            "object",
            "null"
          ]
        },
        "xenbase_url": {
          "type": [
            "string",
            "null"
          ]
        },
        "cross_references": {
          "type": [
            "object",
            "null"
          ]
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `Xenbase_get_gene`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("xenbase_get_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.xenbase_search_genes` {#dnakit.databases.xenbase_search_genes}

Search for Xenopus (frog) genes in Xenbase via the Alliance of Genome Resources API. Returns gene IDs, symbols, names, and species for both Xenopus tropicalis and Xenopus laevis. Xenbase is the primary database for Xenopus genetics and developmental biology. Optionally filter by species (tropicalis or laevis). Useful for finding amphibian homologs of human disease genes.

```text
dnakit.databases.xenbase_search_genes(*, query, species=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填（输入 Schema） | Gene symbol or name to search. Examples: &#x27;tp53&#x27;, &#x27;pax6&#x27;, &#x27;sox2&#x27;, &#x27;bmp4&#x27;, &#x27;wnt3a&#x27;, &#x27;fgf8&#x27;. |
| `species` | string | None | Optional species filter: &#x27;tropicalis&#x27; for X. tropicalis, &#x27;laevis&#x27; for X. laevis. Leave empty for both species. {&quot;default&quot;: &quot;&quot;} |
| `limit` | integer | None | Maximum number of results to return. Default: 10. {&quot;default&quot;: 10} |
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
          "description": "Gene symbol or name to search. Examples: 'tp53', 'pax6', 'sox2', 'bmp4', 'wnt3a', 'fgf8'."
        },
        "species": {
          "type": "string",
          "description": "Optional species filter: 'tropicalis' for X. tropicalis, 'laevis' for X. laevis. Leave empty for both species.",
          "default": ""
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of results to return. Default: 10.",
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
          "gene_id": {
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
          "synonyms": {
            "type": "array",
            "items": {
              "type": "string"
            }
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `Xenbase_search_genes`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("xenbase_search_genes", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.zfin_get_gene` {#dnakit.databases.zfin_get_gene}

Get detailed Danio rerio (zebrafish) gene information from ZFIN via the Alliance of Genome Resources. Returns gene symbol, name, synopsis, genomic location, synonyms, and cross-references. ZFIN is the primary database for zebrafish genetics. Example: 'ZFIN:ZDB-GENE-990415-8' returns pax2a (paired box 2a) with its role in nervous system development.

```text
dnakit.databases.zfin_get_gene(*, gene_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填（输入 Schema） | ZFIN gene ID with &#x27;ZFIN:&#x27; prefix. Examples: &#x27;ZFIN:ZDB-GENE-990415-8&#x27; (pax2a), &#x27;ZFIN:ZDB-GENE-980526-166&#x27; (shha/sonic hedgehog a), &#x27;ZFIN:ZDB-GENE-990415-72&#x27; (tp53/tumor protein p53), &#x27;ZFIN:ZDB-GENE-030131-7740&#x27; (bmp4/bone morphogenetic protein 4). |
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
          "description": "ZFIN gene ID with 'ZFIN:' prefix. Examples: 'ZFIN:ZDB-GENE-990415-8' (pax2a), 'ZFIN:ZDB-GENE-980526-166' (shha/sonic hedgehog a), 'ZFIN:ZDB-GENE-990415-72' (tp53/tumor protein p53), 'ZFIN:ZDB-GENE-030131-7740' (bmp4/bone morphogenetic protein 4)."
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
            "id": {
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
            "name": {
              "type": [
                "string",
                "null"
              ]
            },
            "species": {
              "type": "object",
              "properties": {
                "name": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "short_name": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "taxon_id": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "data_provider": {
                  "type": [
                    "string",
                    "null"
                  ]
                }
              }
            },
            "gene_synopsis": {
              "type": [
                "string",
                "null"
              ]
            },
            "automated_gene_synopsis": {
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
            "so_term": {
              "type": [
                "string",
                "null"
              ]
            },
            "genomic_location": {
              "type": "object",
              "properties": {
                "chromosome": {
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
                "assembly": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "strand": {
                  "type": [
                    "string",
                    "null"
                  ]
                }
              }
            },
            "cross_references": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "name": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "url": {
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

**执行与来源：** `network`；上游接口名 `ZFIN_get_gene`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("zfin_get_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.zfin_get_gene_expression` {#dnakit.databases.zfin_get_gene_expression}

Get expression annotations for a zebrafish gene from the Alliance of Genome Resources. Returns per-annotation records of developmental stage (whenExpressed) and anatomical location (whereExpressed) for the gene, plus the total annotation count. Example: 'ZFIN:ZDB-GENE-980526-166' (shha) has hundreds of annotations across developmental stages. Zebrafish is a premier model for vertebrate development.

```text
dnakit.databases.zfin_get_gene_expression(*, gene_id, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填（输入 Schema） | ZFIN gene ID, with or without the &#x27;ZFIN:&#x27; prefix. Examples: &#x27;ZFIN:ZDB-GENE-990415-8&#x27; (pax2a), &#x27;ZFIN:ZDB-GENE-980526-166&#x27; (shha), &#x27;ZFIN:ZDB-GENE-990415-72&#x27; (tp53). |
| `limit` | integer | None | Max annotations to return (default 50, max 100). {&quot;default&quot;: 50} |
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
          "description": "ZFIN gene ID, with or without the 'ZFIN:' prefix. Examples: 'ZFIN:ZDB-GENE-990415-8' (pax2a), 'ZFIN:ZDB-GENE-980526-166' (shha), 'ZFIN:ZDB-GENE-990415-72' (tp53)."
        },
        "limit": {
          "type": "integer",
          "description": "Max annotations to return (default 50, max 100).",
          "default": 50
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
            "total_annotations": {
              "type": "integer"
            },
            "returned": {
              "type": "integer"
            },
            "annotations": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "stage": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "location": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "data_provider": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "reference": {
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
            "status": {
              "type": "string"
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

**执行与来源：** `network`；上游接口名 `ZFIN_get_gene_expression`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("zfin_get_gene_expression", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.zfin_get_gene_phenotypes` {#dnakit.databases.zfin_get_gene_phenotypes}

Get phenotype annotations for a zebrafish gene from ZFIN via the Alliance of Genome Resources. Returns phenotypes observed when the gene is mutated or knocked out in zebrafish, including affected anatomical structures and phenotype statements. Zebrafish is widely used for forward genetic screens and disease modeling. Example: 'ZFIN:ZDB-GENE-990415-8' (pax2a) returns phenotypes related to eye, kidney, and nervous system development.

```text
dnakit.databases.zfin_get_gene_phenotypes(*, gene_id, limit=None, page=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填（输入 Schema） | ZFIN gene ID with &#x27;ZFIN:&#x27; prefix. Examples: &#x27;ZFIN:ZDB-GENE-990415-8&#x27; (pax2a), &#x27;ZFIN:ZDB-GENE-990415-72&#x27; (tp53), &#x27;ZFIN:ZDB-GENE-980526-166&#x27; (shha). |
| `limit` | integer | None | Maximum number of phenotypes to return (1-100). Default: 20. |
| `page` | integer | None | Page number for pagination. Default: 1. |
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
          "description": "ZFIN gene ID with 'ZFIN:' prefix. Examples: 'ZFIN:ZDB-GENE-990415-8' (pax2a), 'ZFIN:ZDB-GENE-990415-72' (tp53), 'ZFIN:ZDB-GENE-980526-166' (shha)."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of phenotypes to return (1-100). Default: 20."
        },
        "page": {
          "type": "integer",
          "description": "Page number for pagination. Default: 1."
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
              "phenotype_statement": {
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

**执行与来源：** `network`；上游接口名 `ZFIN_get_gene_phenotypes`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("zfin_get_gene_phenotypes", verbose=False)
print(info["input_schema"])
```

