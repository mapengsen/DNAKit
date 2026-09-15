# 基因、序列、测序数据与下载

本页收录 **41 个新增公开函数**。 [完整索引](index.md) · [调用约定与完整示例](usage.md)

每个函数均列出当前 Python 签名、参数、返回值和已登记的执行要求。提供方说明及字段描述保留源码原文。`None` 默认值在生成的数据库包装中通常表示不向服务传该参数；实际必填项与组合约束以输入 Schema 为准。坐标系和基因组版本按各接口说明使用。

## 函数目录

- [`dnakit.databases.download_file`](#dnakit.databases.download_file)
- [`dnakit.databases.enaportal_search_samples`](#dnakit.databases.enaportal_search_samples)
- [`dnakit.databases.enaportal_search_studies`](#dnakit.databases.enaportal_search_studies)
- [`dnakit.databases.ensembl_compara_get_gene_tree`](#dnakit.databases.ensembl_compara_get_gene_tree)
- [`dnakit.databases.ensembl_compara_get_orthologues`](#dnakit.databases.ensembl_compara_get_orthologues)
- [`dnakit.databases.ensembl_compara_get_paralogues`](#dnakit.databases.ensembl_compara_get_paralogues)
- [`dnakit.databases.ensembl_get_homology`](#dnakit.databases.ensembl_get_homology)
- [`dnakit.databases.ensembl_get_overlap_features`](#dnakit.databases.ensembl_get_overlap_features)
- [`dnakit.databases.ensembl_get_regulatory_features`](#dnakit.databases.ensembl_get_regulatory_features)
- [`dnakit.databases.ensembl_get_sequence`](#dnakit.databases.ensembl_get_sequence)
- [`dnakit.databases.ensembl_get_structural_variants`](#dnakit.databases.ensembl_get_structural_variants)
- [`dnakit.databases.ensembl_get_sv_detail`](#dnakit.databases.ensembl_get_sv_detail)
- [`dnakit.databases.ensembl_lookup_gene`](#dnakit.databases.ensembl_lookup_gene)
- [`dnakit.databases.ensembl_seq_get_id_sequence`](#dnakit.databases.ensembl_seq_get_id_sequence)
- [`dnakit.databases.ensembl_seq_get_region_sequence`](#dnakit.databases.ensembl_seq_get_region_sequence)
- [`dnakit.databases.ensembl_vep_region`](#dnakit.databases.ensembl_vep_region)
- [`dnakit.databases.hgnc_fetch_gene_by_symbol`](#dnakit.databases.hgnc_fetch_gene_by_symbol)
- [`dnakit.databases.my_gene_batch_query`](#dnakit.databases.my_gene_batch_query)
- [`dnakit.databases.my_gene_get_gene_annotation`](#dnakit.databases.my_gene_get_gene_annotation)
- [`dnakit.databases.my_gene_query_genes`](#dnakit.databases.my_gene_query_genes)
- [`dnakit.databases.ncbi_fetch_accessions`](#dnakit.databases.ncbi_fetch_accessions)
- [`dnakit.databases.ncbi_get_sequence`](#dnakit.databases.ncbi_get_sequence)
- [`dnakit.databases.ncbi_search_nucleotide`](#dnakit.databases.ncbi_search_nucleotide)
- [`dnakit.databases.ncbi_sra_search_runs`](#dnakit.databases.ncbi_sra_search_runs)
- [`dnakit.databases.ncbidatasets_get_gene`](#dnakit.databases.ncbidatasets_get_gene)
- [`dnakit.databases.ncbidatasets_get_gene_by_symbol`](#dnakit.databases.ncbidatasets_get_gene_by_symbol)
- [`dnakit.databases.ncbidatasets_get_genome_assembly`](#dnakit.databases.ncbidatasets_get_genome_assembly)
- [`dnakit.databases.ncbidatasets_get_orthologs`](#dnakit.databases.ncbidatasets_get_orthologs)
- [`dnakit.databases.ncbidatasets_get_sequence_reports`](#dnakit.databases.ncbidatasets_get_sequence_reports)
- [`dnakit.databases.ncbidatasets_get_taxonomy`](#dnakit.databases.ncbidatasets_get_taxonomy)
- [`dnakit.databases.ncbidatasets_list_genomes_by_taxon`](#dnakit.databases.ncbidatasets_list_genomes_by_taxon)
- [`dnakit.databases.ncbidatasets_suggest_taxonomy`](#dnakit.databases.ncbidatasets_suggest_taxonomy)
- [`dnakit.databases.ncbigene_get_summary`](#dnakit.databases.ncbigene_get_summary)
- [`dnakit.databases.ncbigene_search`](#dnakit.databases.ncbigene_search)
- [`dnakit.databases.sra_get_experiment`](#dnakit.databases.sra_get_experiment)
- [`dnakit.databases.sra_search_experiments`](#dnakit.databases.sra_search_experiments)
- [`dnakit.databases.tark_get_mane_transcripts`](#dnakit.databases.tark_get_mane_transcripts)
- [`dnakit.databases.tark_get_transcript`](#dnakit.databases.tark_get_transcript)
- [`dnakit.download.ena.sequence`](#dnakit.download.ena.sequence)
- [`dnakit.download.ncbi.sequence`](#dnakit.download.ncbi.sequence)
- [`dnakit.search.ncbi.accessions`](#dnakit.search.ncbi.accessions)

## `dnakit.databases.download_file` {#dnakit.databases.download_file}

Download files from HTTP/HTTPS URLs with cross-platform support (Windows, Mac, Linux). Similar to curl but platform-independent. Can save to specified path or temporary directory.

```text
dnakit.databases.download_file(*, url, output_path=None, timeout=None, return_content=None, chunk_size=None, follow_redirects=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `url` | string | 必填（输入 Schema） | HTTP or HTTPS URL to download from (e.g., https://example.com/file.txt) {&quot;format&quot;: &quot;uri&quot;} |
| `output_path` | string | None | Optional path to save the file. If not specified, file will be saved to system temp directory. {&quot;example&quot;: &quot;/tmp/downloaded_file.txt or C:\\Users\\Downloads\\file.txt&quot;} |
| `timeout` | integer | None | Request timeout in seconds {&quot;minimum&quot;: 1, &quot;maximum&quot;: 300, &quot;default&quot;: 30} |
| `return_content` | boolean | None | If true, return file content as text instead of saving to disk (default: false) {&quot;default&quot;: false} |
| `chunk_size` | integer | None | Download chunk size in bytes (default: 8192) {&quot;minimum&quot;: 1024, &quot;maximum&quot;: 10485760, &quot;default&quot;: 8192} |
| `follow_redirects` | boolean | None | Follow HTTP redirects (default: true) {&quot;default&quot;: true} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "required": [
        "url"
      ],
      "additionalProperties": false,
      "properties": {
        "url": {
          "type": "string",
          "description": "HTTP or HTTPS URL to download from (e.g., https://example.com/file.txt)",
          "format": "uri"
        },
        "output_path": {
          "type": "string",
          "description": "Optional path to save the file. If not specified, file will be saved to system temp directory.",
          "example": "/tmp/downloaded_file.txt or C:\\Users\\Downloads\\file.txt"
        },
        "timeout": {
          "type": "integer",
          "description": "Request timeout in seconds",
          "minimum": 1,
          "maximum": 300,
          "default": 30
        },
        "return_content": {
          "type": "boolean",
          "description": "If true, return file content as text instead of saving to disk (default: false)",
          "default": false
        },
        "chunk_size": {
          "type": "integer",
          "description": "Download chunk size in bytes (default: 8192)",
          "minimum": 1024,
          "maximum": 10485760,
          "default": 8192
        },
        "follow_redirects": {
          "type": "boolean",
          "description": "Follow HTTP redirects (default: true)",
          "default": true
        }
      }
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "description": "Output of the tool execution.",
      "properties": {}
    }
    ```

**执行与来源：** `filesystem`；上游接口名 `download_file`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("download_file", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.enaportal_search_samples` {#dnakit.databases.enaportal_search_samples}

Search the European Nucleotide Archive (ENA) for biological samples by taxonomy or text query. Returns sample accessions, descriptions, taxonomy, and publication dates. Useful for finding source material metadata across the world's nucleotide archives (INSDC). Example: 'tax_tree(562)' finds all E. coli samples; 'description="liver"' finds samples with liver in the description.

```text
dnakit.databases.enaportal_search_samples(*, query, limit=None, fields=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填（输入 Schema） | ENA search query. Examples: &#x27;tax_tree(562)&#x27; (E. coli), &#x27;description=&quot;liver&quot;&#x27;, &#x27;tax_tree(9606) AND description=&quot;tumor&quot;&#x27;. |
| `limit` | integer &#124; null | None | Maximum results to return (1-100, default 10). |
| `fields` | string &#124; null | None | Comma-separated fields. Default: &#x27;sample_accession,sample_alias,description,tax_id,scientific_name,first_public&#x27;. |
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
          "description": "ENA search query. Examples: 'tax_tree(562)' (E. coli), 'description=\"liver\"', 'tax_tree(9606) AND description=\"tumor\"'."
        },
        "limit": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Maximum results to return (1-100, default 10)."
        },
        "fields": {
          "type": [
            "string",
            "null"
          ],
          "description": "Comma-separated fields. Default: 'sample_accession,sample_alias,description,tax_id,scientific_name,first_public'."
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
              "sample_accession": {
                "type": "string"
              },
              "sample_alias": {
                "type": "string"
              },
              "description": {
                "type": "string"
              },
              "tax_id": {
                "type": "string"
              },
              "scientific_name": {
                "type": "string"
              },
              "first_public": {
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

**执行与来源：** `network`；上游接口名 `ENAPortal_search_samples`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("enaportal_search_samples", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.enaportal_search_studies` {#dnakit.databases.enaportal_search_studies}

Search the European Nucleotide Archive (ENA) for sequencing studies using text queries or taxonomy filters. Returns study accessions, titles, submitting centers, and publication dates. ENA is the world's largest repository of nucleotide sequences (>3 billion records). Supports powerful query syntax: use 'description="keyword"' for text search, 'tax_tree(TAXID)' for taxonomy-based queries, or combine with AND/OR. Example: searching 'description="SARS-CoV-2"' returns COVID-19 sequencing studies.

```text
dnakit.databases.enaportal_search_studies(*, query, limit=None, fields=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填（输入 Schema） | ENA search query. IMPORTANT: Plain text keywords (e.g., &#x27;cancer&#x27;) will cause HTTP 400. Always wrap text keywords with description=&quot;...&quot;. Examples: &#x27;description=&quot;cancer&quot;&#x27;, &#x27;tax_tree(9606)&#x27; (human), &#x27;description=&quot;transcriptomics&quot; AND tax_tree(10090)&#x27; (mouse transcriptomics). |
| `limit` | integer &#124; null | None | Maximum results to return (1-100, default 10). |
| `fields` | string &#124; null | None | Comma-separated fields to return. Default: &#x27;study_accession,study_title,center_name,first_public,description&#x27;. Other fields: &#x27;last_updated,tax_id,scientific_name&#x27;. |
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
          "description": "ENA search query. IMPORTANT: Plain text keywords (e.g., 'cancer') will cause HTTP 400. Always wrap text keywords with description=\"...\". Examples: 'description=\"cancer\"', 'tax_tree(9606)' (human), 'description=\"transcriptomics\" AND tax_tree(10090)' (mouse transcriptomics)."
        },
        "limit": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Maximum results to return (1-100, default 10)."
        },
        "fields": {
          "type": [
            "string",
            "null"
          ],
          "description": "Comma-separated fields to return. Default: 'study_accession,study_title,center_name,first_public,description'. Other fields: 'last_updated,tax_id,scientific_name'."
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
              "study_accession": {
                "type": "string"
              },
              "study_title": {
                "type": "string"
              },
              "center_name": {
                "type": "string"
              },
              "first_public": {
                "type": "string"
              },
              "description": {
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

**执行与来源：** `network`；上游接口名 `ENAPortal_search_studies`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("enaportal_search_studies", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ensembl_compara_get_gene_tree` {#dnakit.databases.ensembl_compara_get_gene_tree}

Get the phylogenetic gene tree showing evolutionary relationships of a gene family across species from Ensembl Compara. Gene trees trace the history of speciation and duplication events for a gene family. Returns tree members (genes across species), tree structure, and optionally Newick format. Useful for evolutionary analysis and understanding gene family expansion. Example: TP53 gene tree shows conservation across vertebrates.

```text
dnakit.databases.ensembl_compara_get_gene_tree(*, gene, species=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | 必填（输入 Schema） | Gene symbol or Ensembl gene ID. Examples: &#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;ENSG00000141510&#x27;. Ensembl IDs are more reliable than symbols for gene tree retrieval. |
| `species` | string | None | Species (default: &#x27;human&#x27;). Used when providing a gene symbol. Examples: &#x27;human&#x27;, &#x27;mouse&#x27;. |
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
          "description": "Gene symbol or Ensembl gene ID. Examples: 'TP53', 'BRCA1', 'ENSG00000141510'. Ensembl IDs are more reliable than symbols for gene tree retrieval."
        },
        "species": {
          "type": "string",
          "description": "Species (default: 'human'). Used when providing a gene symbol. Examples: 'human', 'mouse'."
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
            "tree_id": {
              "type": [
                "string",
                "null"
              ],
              "description": "Gene tree stable ID"
            },
            "newick": {
              "type": [
                "string",
                "null"
              ],
              "description": "Tree in Newick format (if available)"
            },
            "rooted": {
              "type": [
                "boolean",
                "integer",
                "null"
              ],
              "description": "Whether tree is rooted (may be 0/1 integer or boolean)"
            },
            "members": {
              "type": "array",
              "items": {
                "type": "object"
              },
              "description": "Gene tree members (gene_id, species) - max 50 shown"
            },
            "total_members": {
              "type": "integer",
              "description": "Total members in the gene tree"
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

**执行与来源：** `network`；上游接口名 `EnsemblCompara_get_gene_tree`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ensembl_compara_get_gene_tree", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ensembl_compara_get_orthologues` {#dnakit.databases.ensembl_compara_get_orthologues}

Find orthologous genes (genes in different species that evolved from a common ancestor) using Ensembl Compara. Essential for comparative genomics - finding the mouse, zebrafish, or fly equivalent of a human gene. Returns orthologue gene IDs, species, homology type (one-to-one, one-to-many), and taxonomy level of divergence. Example: TP53 in human has a one-to-one orthologue (Trp53) in mouse.

```text
dnakit.databases.ensembl_compara_get_orthologues(*, gene, species=None, target_species=None, target_taxon=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | 必填（输入 Schema） | Gene symbol or Ensembl gene ID. Examples: &#x27;BRCA1&#x27;, &#x27;TP53&#x27;, &#x27;ENSG00000141510&#x27;, &#x27;ENSG00000012048&#x27;. |
| `species` | string | None | Source species (default: &#x27;human&#x27;). Examples: &#x27;human&#x27;, &#x27;mouse&#x27;, &#x27;zebrafish&#x27;, &#x27;drosophila_melanogaster&#x27;. |
| `target_species` | string | None | Target species to limit orthologue search. Examples: &#x27;mouse&#x27;, &#x27;zebrafish&#x27;, &#x27;rat&#x27;, &#x27;drosophila_melanogaster&#x27;. Default: all species. |
| `target_taxon` | integer | None | NCBI taxon ID to limit orthologue search. Examples: 10090 (mouse), 7955 (zebrafish), 7227 (fruit fly). Alternative to target_species. |
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
          "description": "Gene symbol or Ensembl gene ID. Examples: 'BRCA1', 'TP53', 'ENSG00000141510', 'ENSG00000012048'."
        },
        "species": {
          "type": "string",
          "description": "Source species (default: 'human'). Examples: 'human', 'mouse', 'zebrafish', 'drosophila_melanogaster'."
        },
        "target_species": {
          "type": "string",
          "description": "Target species to limit orthologue search. Examples: 'mouse', 'zebrafish', 'rat', 'drosophila_melanogaster'. Default: all species."
        },
        "target_taxon": {
          "type": "integer",
          "description": "NCBI taxon ID to limit orthologue search. Examples: 10090 (mouse), 7955 (zebrafish), 7227 (fruit fly). Alternative to target_species."
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
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "source_gene": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Source Ensembl gene ID"
              },
              "target_gene": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Target orthologue Ensembl gene ID"
              },
              "target_protein": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Target orthologue Ensembl protein ID"
              },
              "target_species": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Target species name"
              },
              "homology_type": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Type: ortholog_one2one, ortholog_one2many, ortholog_many2many"
              },
              "taxonomy_level": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Taxonomic level of divergence (e.g., Euarchontoglires, Mammalia)"
              },
              "method": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Method used (ENSEMBL_ORTHOLOGUES)"
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

**执行与来源：** `network`；上游接口名 `EnsemblCompara_get_orthologues`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ensembl_compara_get_orthologues", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ensembl_compara_get_paralogues` {#dnakit.databases.ensembl_compara_get_paralogues}

Find paralogous genes (genes within the same species that arose from gene duplication events) using Ensembl Compara. Paralogues often have related but diverged functions. Returns paralogue gene IDs, duplication type, and evolutionary context. Example: TP53 has paralogues TP63 (ENSG00000073282) and TP73 (ENSG00000078900) in human, from ancient vertebrate gene duplication.

```text
dnakit.databases.ensembl_compara_get_paralogues(*, gene, species=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | 必填（输入 Schema） | Gene symbol or Ensembl gene ID. Examples: &#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;ENSG00000141510&#x27;. |
| `species` | string | None | Species (default: &#x27;human&#x27;). Examples: &#x27;human&#x27;, &#x27;mouse&#x27;, &#x27;zebrafish&#x27;. |
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
          "description": "Gene symbol or Ensembl gene ID. Examples: 'TP53', 'BRCA1', 'ENSG00000141510'."
        },
        "species": {
          "type": "string",
          "description": "Species (default: 'human'). Examples: 'human', 'mouse', 'zebrafish'."
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
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "source_gene": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Source Ensembl gene ID"
              },
              "paralogue_gene": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Paralogue Ensembl gene ID"
              },
              "paralogue_protein": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Paralogue Ensembl protein ID"
              },
              "species": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Species name"
              },
              "paralogy_type": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Type of paralogy (within_species_paralog, etc.)"
              },
              "taxonomy_level": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Level of gene duplication event"
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

**执行与来源：** `network`；上游接口名 `EnsemblCompara_get_paralogues`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ensembl_compara_get_paralogues", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ensembl_get_homology` {#dnakit.databases.ensembl_get_homology}

Get homology (orthologues and paralogues) for a gene by symbol across species. Returns evolutionary relationships with gene IDs, protein IDs, chromosomal locations, and homology types (1-to-1, 1-to-many, etc). Essential for comparative genomics studies.

```text
dnakit.databases.ensembl_get_homology(*, species, symbol, target_species=None, target_taxon=None, type=None, sequence=None, aligned=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `species` | string | 必填（输入 Schema） | Source species (e.g., &#x27;human&#x27;, &#x27;homo_sapiens&#x27;) |
| `symbol` | string | 必填（输入 Schema） | Gene symbol (e.g., &#x27;BRCA2&#x27;, &#x27;TP53&#x27;) |
| `target_species` | string | None | Target species to find homologues in (optional, e.g., &#x27;mouse&#x27;, &#x27;mus_musculus&#x27;). Leave empty to search all species. |
| `target_taxon` | string | None | Target taxonomic group (optional, e.g., &#x27;Mammalia&#x27;, &#x27;Primates&#x27;) |
| `type` | string | None | Homology type: &#x27;orthologues&#x27; (different species), &#x27;paralogues&#x27; (same species), &#x27;all&#x27; (both) {&quot;enum&quot;: &#91;&quot;orthologues&quot;, &quot;paralogues&quot;, &quot;all&quot;&#93;, &quot;default&quot;: &quot;all&quot;} |
| `sequence` | string | None | Include sequences in response {&quot;enum&quot;: &#91;&quot;cdna&quot;, &quot;protein&quot;, &quot;none&quot;&#93;, &quot;default&quot;: &quot;none&quot;} |
| `aligned` | boolean | None | Return aligned sequences {&quot;default&quot;: false} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "species": {
          "type": "string",
          "description": "Source species (e.g., 'human', 'homo_sapiens')"
        },
        "symbol": {
          "type": "string",
          "description": "Gene symbol (e.g., 'BRCA2', 'TP53')"
        },
        "target_species": {
          "type": "string",
          "description": "Target species to find homologues in (optional, e.g., 'mouse', 'mus_musculus'). Leave empty to search all species."
        },
        "target_taxon": {
          "type": "string",
          "description": "Target taxonomic group (optional, e.g., 'Mammalia', 'Primates')"
        },
        "type": {
          "type": "string",
          "enum": [
            "orthologues",
            "paralogues",
            "all"
          ],
          "default": "all",
          "description": "Homology type: 'orthologues' (different species), 'paralogues' (same species), 'all' (both)"
        },
        "sequence": {
          "type": "string",
          "enum": [
            "cdna",
            "protein",
            "none"
          ],
          "default": "none",
          "description": "Include sequences in response"
        },
        "aligned": {
          "type": "boolean",
          "default": false,
          "description": "Return aligned sequences"
        }
      },
      "required": [
        "species",
        "symbol"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": [
        "array",
        "object"
      ],
      "items": {
        "type": "object",
        "properties": {
          "homologies": {
            "type": "array",
            "items": {
              "type": "object",
              "properties": {
                "type": {
                  "type": "string",
                  "description": "Homology type"
                },
                "target": {
                  "type": "object",
                  "properties": {
                    "id": {
                      "type": "string"
                    },
                    "species": {
                      "type": "string"
                    },
                    "protein_id": {
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

**执行与来源：** `network`；上游接口名 `ensembl_get_homology`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ensembl_get_homology", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ensembl_get_overlap_features` {#dnakit.databases.ensembl_get_overlap_features}

Get genomic features (genes, transcripts, repeats) overlapping a region. Returns comprehensive annotations including coding genes, ncRNAs, pseudogenes, and structural features. Use this for genome annotation and feature discovery.

```text
dnakit.databases.ensembl_get_overlap_features(*, species, region, feature=None, biotype=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `species` | string | 必填（输入 Schema） | Species (e.g., &#x27;human&#x27;, &#x27;homo_sapiens&#x27;) |
| `region` | string | 必填（输入 Schema） | Genomic region (e.g., &#x27;7:140424943-140624564&#x27;) |
| `feature` | string | None | Feature type to retrieve {&quot;enum&quot;: &#91;&quot;gene&quot;, &quot;transcript&quot;, &quot;cds&quot;, &quot;exon&quot;, &quot;repeat&quot;, &quot;simple_feature&quot;, &quot;misc_feature&quot;, &quot;constrained_element&quot;&#93;, &quot;default&quot;: &quot;gene&quot;} |
| `biotype` | string | None | Filter by biotype (optional, e.g., &#x27;protein_coding&#x27;, &#x27;lincRNA&#x27;, &#x27;miRNA&#x27;) |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "species": {
          "type": "string",
          "description": "Species (e.g., 'human', 'homo_sapiens')"
        },
        "region": {
          "type": "string",
          "description": "Genomic region (e.g., '7:140424943-140624564')"
        },
        "feature": {
          "type": "string",
          "enum": [
            "gene",
            "transcript",
            "cds",
            "exon",
            "repeat",
            "simple_feature",
            "misc_feature",
            "constrained_element"
          ],
          "default": "gene",
          "description": "Feature type to retrieve"
        },
        "biotype": {
          "type": "string",
          "description": "Filter by biotype (optional, e.g., 'protein_coding', 'lincRNA', 'miRNA')"
        }
      },
      "required": [
        "region",
        "species"
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
          "id": {
            "type": "string"
          },
          "external_name": {
            "type": "string"
          },
          "biotype": {
            "type": "string"
          },
          "start": {
            "type": "integer"
          },
          "end": {
            "type": "integer"
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `ensembl_get_overlap_features`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ensembl_get_overlap_features", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ensembl_get_regulatory_features` {#dnakit.databases.ensembl_get_regulatory_features}

Get regulatory features (promoters, enhancers, transcription factor binding sites) in a genomic region. Returns regulatory annotations including feature types, positions, and evidence. Use this to identify regulatory elements that may control gene expression.

```text
dnakit.databases.ensembl_get_regulatory_features(*, region, species=None, feature, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `region` | string | 必填（输入 Schema） | Genomic region in format &#x27;chromosome:start-end&#x27; (e.g., &#x27;1:1000000-2000000&#x27;) |
| `species` | string | None | Species name (default &#x27;human&#x27;) {&quot;default&quot;: &quot;human&quot;} |
| `feature` | string | 必填（输入 Schema） | Feature type (must be &#x27;regulatory&#x27; for this tool) |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "region": {
          "type": "string",
          "description": "Genomic region in format 'chromosome:start-end' (e.g., '1:1000000-2000000')"
        },
        "species": {
          "type": "string",
          "default": "human",
          "description": "Species name (default 'human')"
        },
        "feature": {
          "type": "string",
          "description": "Feature type (must be 'regulatory' for this tool)"
        }
      },
      "required": [
        "feature",
        "region"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "array",
      "description": "Array of regulatory feature objects",
      "items": {
        "type": "object",
        "properties": {
          "id": {
            "type": "string",
            "description": "Regulatory feature ID"
          },
          "feature_type": {
            "type": "string",
            "description": "Feature type (e.g., 'Promoter', 'Enhancer', 'TF_binding_site')"
          },
          "start": {
            "type": "integer",
            "description": "Feature start position"
          },
          "end": {
            "type": "integer",
            "description": "Feature end position"
          },
          "strand": {
            "type": "integer",
            "description": "Feature strand"
          },
          "seq_region_name": {
            "type": "string",
            "description": "Chromosome name"
          },
          "description": {
            "type": "string",
            "description": "Feature description"
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `ensembl_get_regulatory_features`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ensembl_get_regulatory_features", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ensembl_get_sequence` {#dnakit.databases.ensembl_get_sequence}

Get DNA or protein sequence for a gene, transcript, or genomic region. Returns sequence in JSON format with metadata. Supports multiple sequence types: genomic, CDS, cDNA, and peptide sequences.

```text
dnakit.databases.ensembl_get_sequence(*, id, type=None, species=None, multiple_sequences=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `id` | string | 必填（输入 Schema） | Ensembl gene/transcript ID (e.g., &#x27;ENSG00000139618&#x27; or &#x27;ENST00000380152&#x27;) |
| `type` | string | None | Sequence type: &#x27;genomic&#x27; (genomic DNA), &#x27;cds&#x27; (coding sequence), &#x27;cdna&#x27; (cDNA), &#x27;protein&#x27; (protein/peptide sequence -- &#x27;peptide&#x27; is also accepted and mapped to &#x27;protein&#x27;, the name Ensembl&#x27;s own REST API actually expects) {&quot;enum&quot;: &#91;&quot;genomic&quot;, &quot;cds&quot;, &quot;cdna&quot;, &quot;protein&quot;, &quot;peptide&quot;&#93;, &quot;default&quot;: &quot;genomic&quot;} |
| `species` | string | None | Species name (optional, defaults to human) |
| `multiple_sequences` | boolean | None | Return multiple sequences if available (e.g., all transcripts) {&quot;default&quot;: true} |
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
          "description": "Ensembl gene/transcript ID (e.g., 'ENSG00000139618' or 'ENST00000380152')"
        },
        "type": {
          "type": "string",
          "enum": [
            "genomic",
            "cds",
            "cdna",
            "protein",
            "peptide"
          ],
          "default": "genomic",
          "description": "Sequence type: 'genomic' (genomic DNA), 'cds' (coding sequence), 'cdna' (cDNA), 'protein' (protein/peptide sequence -- 'peptide' is also accepted and mapped to 'protein', the name Ensembl's own REST API actually expects)"
        },
        "species": {
          "type": "string",
          "description": "Species name (optional, defaults to human)"
        },
        "multiple_sequences": {
          "type": "boolean",
          "default": true,
          "description": "Return multiple sequences if available (e.g., all transcripts)"
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
      "type": [
        "array",
        "object"
      ],
      "description": "Array of sequence objects or single sequence object",
      "items": {
        "type": "object",
        "properties": {
          "id": {
            "type": "string",
            "description": "Sequence ID"
          },
          "seq": {
            "type": "string",
            "description": "DNA/protein sequence"
          },
          "desc": {
            "type": "string",
            "description": "Sequence description"
          },
          "molecule": {
            "type": "string",
            "description": "Molecule type (dna, rna, protein)"
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `ensembl_get_sequence`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ensembl_get_sequence", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ensembl_get_structural_variants` {#dnakit.databases.ensembl_get_structural_variants}

Get structural variants overlapping a genomic region from Ensembl. Returns known SVs from DGVa (Database of Genomic Variants archive), ClinGen, and other sources. Each SV includes its ID (nsv/esv accession), coordinates, and feature type. Use ensembl_get_sv_detail to get clinical significance, variant class (CNV, deletion, duplication, inversion), and supporting evidence for specific variants. Region must be < 5Mb.

```text
dnakit.databases.ensembl_get_structural_variants(*, species, region=None, feature=None, chrom=None, start=None, end=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `species` | string | 必填（输入 Schema） | Species name (e.g., &#x27;human&#x27;, &#x27;homo_sapiens&#x27;). |
| `region` | string | None | Genomic region in format &#x27;chr:start-end&#x27; (e.g., &#x27;17:43044295-43125370&#x27;). Maximum 5Mb. |
| `feature` | string | None | Feature type to retrieve. Fixed to &#x27;structural_variation&#x27; for this tool. {&quot;default&quot;: &quot;structural_variation&quot;} |
| `chrom` | string | None | 染色体标识；须与 start、end 一起传入。 |
| `start` | integer | None | 区域起点；本参数分支使用 Ensembl 的 1-based 闭区间。 {&quot;minimum&quot;: 1} |
| `end` | integer | None | 区域终点；本参数分支使用 Ensembl 的 1-based 闭区间。 {&quot;minimum&quot;: 1} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "species": {
          "type": "string",
          "description": "Species name (e.g., 'human', 'homo_sapiens')."
        },
        "region": {
          "type": "string",
          "description": "Genomic region in format 'chr:start-end' (e.g., '17:43044295-43125370'). Maximum 5Mb."
        },
        "feature": {
          "type": "string",
          "default": "structural_variation",
          "description": "Feature type to retrieve. Fixed to 'structural_variation' for this tool."
        },
        "chrom": {
          "type": "string"
        },
        "start": {
          "type": "integer",
          "minimum": 1
        },
        "end": {
          "type": "integer",
          "minimum": 1
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
      "type": "array",
      "description": "Array of structural variants overlapping the region.",
      "items": {
        "type": "object",
        "properties": {
          "id": {
            "type": "string",
            "description": "Structural variant accession (nsv or esv ID from DGVa/dbVar)"
          },
          "seq_region_name": {
            "type": "string",
            "description": "Chromosome number"
          },
          "start": {
            "type": "integer",
            "description": "Start coordinate (GRCh38)"
          },
          "end": {
            "type": "integer",
            "description": "End coordinate (GRCh38)"
          },
          "strand": {
            "type": "integer",
            "description": "Strand (1 or -1)"
          },
          "feature_type": {
            "type": "string",
            "description": "Always 'structural_variation'"
          },
          "assembly_name": {
            "type": "string",
            "description": "Genome assembly (e.g., GRCh38)"
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
        "chrom": {
          "type": "string"
        },
        "start": {
          "type": "integer",
          "minimum": 1
        },
        "end": {
          "type": "integer",
          "minimum": 1
        }
      },
      "optional": [
        "region"
      ],
      "handler": "region_coordinates"
    }
    ```

**执行与来源：** `network`；上游接口名 `ensembl_get_structural_variants`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ensembl_get_structural_variants", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ensembl_get_sv_detail` {#dnakit.databases.ensembl_get_sv_detail}

Get detailed information for a specific structural variant from Ensembl by its accession (nsv/esv ID). Returns variant class (CNV, deletion, duplication, inversion, gain, loss), clinical significance (pathogenic, benign, etc.), source database, supporting evidence with individual sample observations, study information, and genomic mappings. Essential for CNV clinical interpretation.

```text
dnakit.databases.ensembl_get_sv_detail(*, species, id=None, sv_id=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `species` | string | 必填（输入 Schema） | Species name (e.g., &#x27;human&#x27;, &#x27;homo_sapiens&#x27;). |
| `id` | string | None | Structural variant accession from Ensembl/DGVa (e.g., &#x27;nsv2769779&#x27;, &#x27;esv3647486&#x27;). Use ensembl_get_structural_variants to find IDs. |
| `sv_id` | string | None | DNAKit optional alias for id |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "species": {
          "type": "string",
          "description": "Species name (e.g., 'human', 'homo_sapiens')."
        },
        "id": {
          "type": "string",
          "description": "Structural variant accession from Ensembl/DGVa (e.g., 'nsv2769779', 'esv3647486'). Use ensembl_get_structural_variants to find IDs."
        },
        "sv_id": {
          "type": "string",
          "description": "DNAKit optional alias for id"
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
      "type": "object",
      "description": "Detailed structural variant information.",
      "properties": {
        "name": {
          "type": [
            "string",
            "null"
          ],
          "description": "Variant name/accession"
        },
        "var_class": {
          "type": "string",
          "description": "Variant class: CNV, deletion, duplication, inversion, insertion, gain, loss, tandem_duplication, complex_structural_alteration"
        },
        "source": {
          "type": "string",
          "description": "Source database (e.g., 'Database of Genomic Variants Archive')"
        },
        "clinical_significance": {
          "type": [
            "array",
            "null"
          ],
          "items": {
            "type": "string"
          },
          "description": "Clinical significance: pathogenic, likely pathogenic, benign, likely benign, uncertain significance"
        },
        "alias": {
          "type": [
            "string",
            "null"
          ],
          "description": "Alternative name or alias"
        },
        "supporting_evidence": {
          "type": [
            "array",
            "null"
          ],
          "description": "Individual sample-level supporting evidence",
          "items": {
            "type": "object",
            "properties": {
              "name": {
                "type": "string",
                "description": "Supporting variant accession (nssv/essv)"
              },
              "var_class": {
                "type": "string",
                "description": "Observed variant class (gain, loss, deletion, etc.)"
              },
              "source": {
                "type": "string"
              },
              "clinical_significance": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "string"
                }
              },
              "copy_number": {
                "type": [
                  "integer",
                  "null"
                ]
              },
              "study": {
                "type": [
                  "object",
                  "null"
                ],
                "properties": {
                  "name": {
                    "type": "string",
                    "description": "Study accession (e.g., estd214, nstd37)"
                  },
                  "description": {
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
                },
                "additionalProperties": true
              },
              "mappings": {
                "type": [
                  "array",
                  "null"
                ],
                "items": {
                  "type": "object",
                  "properties": {
                    "location": {
                      "type": "string"
                    },
                    "seq_region_name": {
                      "type": "string"
                    },
                    "start": {
                      "type": "integer"
                    },
                    "end": {
                      "type": "integer"
                    },
                    "assembly_name": {
                      "type": "string"
                    },
                    "genomic_size": {
                      "type": [
                        "integer",
                        "null"
                      ]
                    }
                  },
                  "additionalProperties": true
                }
              }
            },
            "additionalProperties": true
          }
        }
      },
      "additionalProperties": true
    }
    ```

??? info "DNAKit 参数扩展规则"

    ```json
    {
      "aliases": {
        "sv_id": "id"
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `ensembl_get_sv_detail`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ensembl_get_sv_detail", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ensembl_lookup_gene` {#dnakit.databases.ensembl_lookup_gene}

Lookup gene information by Ensembl gene ID or symbol. Returns gene details including location, description, and external references. Supports both stable IDs (e.g., ENSG00000139618) and gene symbols (e.g., BRCA1). For gene symbols, species parameter is required.

```text
dnakit.databases.ensembl_lookup_gene(*, gene_id=None, species=None, id=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | None | Ensembl gene ID or symbol (e.g., &#x27;ENSG00000139618&#x27; or &#x27;BRCA1&#x27;). If using a stable ID, the tool will automatically route to /lookup/id endpoint. |
| `species` | string | None | Species name required for gene symbols (default &#x27;homo_sapiens&#x27;). Examples: &#x27;homo_sapiens&#x27;, &#x27;mus_musculus&#x27;, &#x27;rattus_norvegicus&#x27; {&quot;default&quot;: &quot;homo_sapiens&quot;} |
| `id` | string | None | DNAKit optional alias for gene_id |
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
          "description": "Ensembl gene ID or symbol (e.g., 'ENSG00000139618' or 'BRCA1'). If using a stable ID, the tool will automatically route to /lookup/id endpoint."
        },
        "species": {
          "type": "string",
          "description": "Species name required for gene symbols (default 'homo_sapiens'). Examples: 'homo_sapiens', 'mus_musculus', 'rattus_norvegicus'",
          "default": "homo_sapiens"
        },
        "id": {
          "type": "string",
          "description": "DNAKit optional alias for gene_id"
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "description": "Gene information object",
      "properties": {
        "id": {
          "type": "string",
          "description": "Ensembl gene ID"
        },
        "display_name": {
          "type": "string",
          "description": "Gene symbol"
        },
        "description": {
          "type": "string",
          "description": "Gene description"
        },
        "biotype": {
          "type": "string",
          "description": "Gene biotype"
        },
        "start": {
          "type": "integer",
          "description": "Gene start position"
        },
        "end": {
          "type": "integer",
          "description": "Gene end position"
        },
        "strand": {
          "type": "integer",
          "description": "Gene strand (1 or -1)"
        },
        "seq_region_name": {
          "type": "string",
          "description": "Chromosome name"
        }
      }
    }
    ```

??? info "DNAKit 参数扩展规则"

    ```json
    {
      "aliases": {
        "id": "gene_id"
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `ensembl_lookup_gene`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ensembl_lookup_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ensembl_seq_get_id_sequence` {#dnakit.databases.ensembl_seq_get_id_sequence}

Get the amino acid or nucleotide sequence for an Ensembl ID (protein, transcript, or gene). Retrieves protein sequences from ENSP IDs, cDNA/CDS from ENST IDs, or genomic DNA from ENSG IDs. Example: ENSP00000269305 returns the 393 amino acid protein sequence of TP53 canonical isoform. ENST00000269305 with type='cdna' returns the 2512 bp cDNA sequence of TP53-201.

```text
dnakit.databases.ensembl_seq_get_id_sequence(*, ensembl_id, type=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `ensembl_id` | string | 必填（输入 Schema） | Ensembl stable ID. Examples: &#x27;ENSP00000269305&#x27; (TP53 protein, 393 aa), &#x27;ENST00000269305&#x27; (TP53-201 transcript), &#x27;ENSG00000141510&#x27; (TP53 gene). Protein IDs start with ENSP, transcript with ENST, gene with ENSG. |
| `type` | string | None | Sequence type to return. Options: &#x27;protein&#x27; (amino acid), &#x27;cdna&#x27; (spliced cDNA), &#x27;cds&#x27; (coding sequence), &#x27;genomic&#x27; (unspliced DNA). Default: &#x27;protein&#x27; for ENSP IDs, auto-detected for others. {&quot;default&quot;: &quot;protein&quot;} |
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
          "description": "Ensembl stable ID. Examples: 'ENSP00000269305' (TP53 protein, 393 aa), 'ENST00000269305' (TP53-201 transcript), 'ENSG00000141510' (TP53 gene). Protein IDs start with ENSP, transcript with ENST, gene with ENSG."
        },
        "type": {
          "type": "string",
          "description": "Sequence type to return. Options: 'protein' (amino acid), 'cdna' (spliced cDNA), 'cds' (coding sequence), 'genomic' (unspliced DNA). Default: 'protein' for ENSP IDs, auto-detected for others.",
          "default": "protein"
        }
      },
      "required": [
        "ensembl_id"
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
            "ensembl_id": {
              "type": "string",
              "description": "Ensembl stable ID"
            },
            "molecule": {
              "type": [
                "string",
                "null"
              ],
              "description": "Molecule type (protein or dna)"
            },
            "description": {
              "type": [
                "string",
                "null"
              ],
              "description": "Sequence description"
            },
            "sequence_length": {
              "type": "integer",
              "description": "Sequence length (aa for protein, bp for DNA)"
            },
            "sequence": {
              "type": "string",
              "description": "Amino acid or nucleotide sequence (truncated to 10kb)"
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

**执行与来源：** `network`；上游接口名 `EnsemblSeq_get_id_sequence`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ensembl_seq_get_id_sequence", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ensembl_seq_get_region_sequence` {#dnakit.databases.ensembl_seq_get_region_sequence}

Get the DNA nucleotide sequence for a specific genomic region from Ensembl. Retrieves the reference genome sequence for given chromosome coordinates. Useful for examining regulatory regions, intergenic sequences, or any arbitrary genomic interval. Example: region '17:7668421-7668520' in homo_sapiens returns 100 bp of the TP53 locus. Supports any Ensembl species and coordinate strand.

```text
dnakit.databases.ensembl_seq_get_region_sequence(*, region, species=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `region` | string | 必填（输入 Schema） | Genomic region in format &#x27;chr:start-end&#x27; or &#x27;chr:start..end:strand&#x27;. Examples: &#x27;17:7668421-7668520&#x27; (TP53 locus, 100bp), &#x27;7:140424943-140425043&#x27; (BRAF locus, 100bp), &#x27;X:15560138-15602945&#x27; (MECP2 gene). Strand: 1 for forward, -1 for reverse. |
| `species` | string | None | Ensembl species name. Default: &#x27;homo_sapiens&#x27;. Examples: &#x27;mus_musculus&#x27;, &#x27;danio_rerio&#x27;, &#x27;drosophila_melanogaster&#x27;. {&quot;default&quot;: &quot;homo_sapiens&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "region": {
          "type": "string",
          "description": "Genomic region in format 'chr:start-end' or 'chr:start..end:strand'. Examples: '17:7668421-7668520' (TP53 locus, 100bp), '7:140424943-140425043' (BRAF locus, 100bp), 'X:15560138-15602945' (MECP2 gene). Strand: 1 for forward, -1 for reverse."
        },
        "species": {
          "type": "string",
          "description": "Ensembl species name. Default: 'homo_sapiens'. Examples: 'mus_musculus', 'danio_rerio', 'drosophila_melanogaster'.",
          "default": "homo_sapiens"
        }
      },
      "required": [
        "region"
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
            "species": {
              "type": "string",
              "description": "Species name"
            },
            "region": {
              "type": "string",
              "description": "Genomic region queried"
            },
            "molecule": {
              "type": "string",
              "description": "Molecule type (dna)"
            },
            "description": {
              "type": [
                "string",
                "null"
              ],
              "description": "Region description from Ensembl"
            },
            "sequence_length": {
              "type": "integer",
              "description": "Length of returned sequence in bp"
            },
            "sequence": {
              "type": "string",
              "description": "DNA sequence (truncated to 10kb if longer)"
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

**执行与来源：** `network`；上游接口名 `EnsemblSeq_get_region_sequence`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ensembl_seq_get_region_sequence", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ensembl_vep_region` {#dnakit.databases.ensembl_vep_region}

Variant Effect Predictor (VEP) for genomic variants. Predicts functional consequences of variants including impact on genes, transcripts, and proteins. Returns SIFT/PolyPhen scores, consequence types, and affected transcripts. Essential for variant interpretation.

```text
dnakit.databases.ensembl_vep_region(*, species, region, allele=None, AncestralAllele=None, Blosum62=None, CADD=None, Conservation=None, DisGeNET=None, GeneSplicer=None, GO=None, LoF=None, MaxEntScan=None, Phenotypes=None, SIFT=None, PolyPhen=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `species` | string | 必填（输入 Schema） | Species (e.g., &#x27;human&#x27;, &#x27;homo_sapiens&#x27;) |
| `region` | string | 必填（输入 Schema） | Genomic region with variant (e.g., &#x27;21:25891796-25891796&#x27;) |
| `allele` | string | None | VEP 区域变异的等位基因参数，遵循 Ensembl VEP 格式。 |
| `AncestralAllele` | integer | None | Include ancestral allele information {&quot;default&quot;: 0} |
| `Blosum62` | integer | None | Include BLOSUM62 protein conservation score {&quot;default&quot;: 0} |
| `CADD` | integer | None | Include CADD scores {&quot;default&quot;: 0} |
| `Conservation` | integer | None | Include conservation scores {&quot;default&quot;: 0} |
| `DisGeNET` | integer | None | Include DisGeNET disease associations {&quot;default&quot;: 0} |
| `GeneSplicer` | integer | None | Include GeneSplicer predictions {&quot;default&quot;: 0} |
| `GO` | integer | None | Include Gene Ontology terms {&quot;default&quot;: 0} |
| `LoF` | integer | None | Include Loss-of-Function predictions {&quot;default&quot;: 0} |
| `MaxEntScan` | integer | None | Include MaxEntScan splice predictions {&quot;default&quot;: 0} |
| `Phenotypes` | integer | None | Include phenotype annotations {&quot;default&quot;: 0} |
| `SIFT` | string | None | SIFT prediction: p=prediction, s=score, b=both {&quot;enum&quot;: &#91;&quot;p&quot;, &quot;s&quot;, &quot;b&quot;&#93;} |
| `PolyPhen` | string | None | PolyPhen prediction: p=prediction, s=score, b=both {&quot;enum&quot;: &#91;&quot;p&quot;, &quot;s&quot;, &quot;b&quot;&#93;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "species": {
          "type": "string",
          "description": "Species (e.g., 'human', 'homo_sapiens')"
        },
        "region": {
          "type": "string",
          "description": "Genomic region with variant (e.g., '21:25891796-25891796')"
        },
        "allele": {
          "type": "string"
        },
        "AncestralAllele": {
          "type": "integer",
          "default": 0,
          "description": "Include ancestral allele information"
        },
        "Blosum62": {
          "type": "integer",
          "default": 0,
          "description": "Include BLOSUM62 protein conservation score"
        },
        "CADD": {
          "type": "integer",
          "default": 0,
          "description": "Include CADD scores"
        },
        "Conservation": {
          "type": "integer",
          "default": 0,
          "description": "Include conservation scores"
        },
        "DisGeNET": {
          "type": "integer",
          "default": 0,
          "description": "Include DisGeNET disease associations"
        },
        "GeneSplicer": {
          "type": "integer",
          "default": 0,
          "description": "Include GeneSplicer predictions"
        },
        "GO": {
          "type": "integer",
          "default": 0,
          "description": "Include Gene Ontology terms"
        },
        "LoF": {
          "type": "integer",
          "default": 0,
          "description": "Include Loss-of-Function predictions"
        },
        "MaxEntScan": {
          "type": "integer",
          "default": 0,
          "description": "Include MaxEntScan splice predictions"
        },
        "Phenotypes": {
          "type": "integer",
          "default": 0,
          "description": "Include phenotype annotations"
        },
        "SIFT": {
          "type": "string",
          "enum": [
            "p",
            "s",
            "b"
          ],
          "description": "SIFT prediction: p=prediction, s=score, b=both"
        },
        "PolyPhen": {
          "type": "string",
          "enum": [
            "p",
            "s",
            "b"
          ],
          "description": "PolyPhen prediction: p=prediction, s=score, b=both"
        }
      },
      "required": [
        "region",
        "species"
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
          "allele_string": {
            "type": "string"
          },
          "consequence_terms": {
            "type": "array"
          },
          "most_severe_consequence": {
            "type": "string"
          },
          "transcript_consequences": {
            "type": "array"
          }
        }
      }
    }
    ```

??? info "DNAKit 参数扩展规则"

    ```json
    {
      "parameters": {
        "allele": {
          "type": "string"
        }
      },
      "optional": [
        "allele"
      ],
      "handler": "vep_region"
    }
    ```

**执行与来源：** `network`；上游接口名 `ensembl_vep_region`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ensembl_vep_region", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.hgnc_fetch_gene_by_symbol` {#dnakit.databases.hgnc_fetch_gene_by_symbol}

Fetch detailed gene information from HGNC (HUGO Gene Nomenclature Committee) by official gene symbol. Returns comprehensive gene data including HGNC ID, full name, chromosomal location, cross-references to Ensembl, NCBI Entrez, UniProt, OMIM, and RefSeq. HGNC is the authoritative source for human gene nomenclature. Example: 'TP53' returns tumor protein p53 with all its identifiers. Retired symbols are resolved automatically: if the symbol is not a current approved symbol, it is looked up as a previous symbol and then as an alias symbol, and the current gene's full record is returned with the substitution disclosed in metadata ('resolved_from', 'resolved_symbol', 'resolution_relation', 'symbol_resolution_note'). Example: 'CYP2E' returns CYP2E1 (renamed in 2002), 'DIA4' returns NQO1. If the retired symbol is ambiguous (shared by several genes) an error naming every candidate is returned instead of an arbitrary pick; if the symbol is unknown to HGNC entirely, metadata carries 'no_results_note' pointing at HGNC_search_genes.

```text
dnakit.databases.hgnc_fetch_gene_by_symbol(*, symbol, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `symbol` | string | 必填（输入 Schema） | Official HGNC gene symbol. Examples: &#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;EGFR&#x27;, &#x27;INS&#x27;, &#x27;ACE2&#x27;. |
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
          "description": "Official HGNC gene symbol. Examples: 'TP53', 'BRCA1', 'EGFR', 'INS', 'ACE2'."
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
            "hgnc_id": {
              "type": "string",
              "description": "HGNC identifier (e.g., HGNC:11998)"
            },
            "symbol": {
              "type": "string",
              "description": "Approved gene symbol"
            },
            "name": {
              "type": "string",
              "description": "Approved gene name"
            },
            "locus_type": {
              "type": "string",
              "description": "Type of locus (e.g., gene with protein product)"
            },
            "locus_group": {
              "type": "string",
              "description": "Locus group classification"
            },
            "status": {
              "type": "string",
              "description": "Symbol status (Approved, etc.)"
            },
            "location": {
              "type": [
                "string",
                "null"
              ],
              "description": "Chromosomal location"
            },
            "entrez_id": {
              "type": [
                "string",
                "null"
              ],
              "description": "NCBI Entrez Gene ID"
            },
            "ensembl_gene_id": {
              "type": [
                "string",
                "null"
              ],
              "description": "Ensembl Gene ID"
            },
            "uniprot_ids": {
              "type": [
                "array",
                "null"
              ],
              "description": "UniProt accession(s)"
            },
            "refseq_accession": {
              "type": [
                "array",
                "null"
              ],
              "description": "RefSeq accession(s)"
            },
            "omim_id": {
              "type": [
                "array",
                "null"
              ],
              "description": "OMIM ID(s)"
            },
            "ccds_id": {
              "type": [
                "array",
                "null"
              ],
              "description": "CCDS ID(s)"
            },
            "prev_symbol": {
              "type": [
                "array",
                "null"
              ],
              "description": "Previous gene symbols"
            },
            "prev_name": {
              "type": [
                "array",
                "null"
              ],
              "description": "Previous gene names"
            },
            "alias_symbol": {
              "type": [
                "array",
                "null"
              ],
              "description": "Alias symbols"
            },
            "gene_group": {
              "type": [
                "array",
                "null"
              ],
              "description": "Gene group memberships"
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

**执行与来源：** `network`；上游接口名 `HGNC_fetch_gene_by_symbol`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("hgnc_fetch_gene_by_symbol", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.my_gene_batch_query` {#dnakit.databases.my_gene_batch_query}

Query multiple genes at once using a list of gene identifiers. Efficient for retrieving annotations for a gene set or performing ID mapping between different gene ID types.

```text
dnakit.databases.my_gene_batch_query(*, gene_ids, species=None, fields=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_ids` | array | 必填（输入 Schema） | List of gene IDs to query. Can mix different ID types (Entrez, Ensembl, symbols). {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `species` | string | None | Species filter for the query. {&quot;default&quot;: &quot;human&quot;} |
| `fields` | string | None | Comma-separated list of fields to return for each gene. {&quot;default&quot;: &quot;symbol,name,entrezgene,ensembl.gene&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene_ids": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "List of gene IDs to query. Can mix different ID types (Entrez, Ensembl, symbols)."
        },
        "species": {
          "type": "string",
          "description": "Species filter for the query.",
          "default": "human"
        },
        "fields": {
          "type": "string",
          "description": "Comma-separated list of fields to return for each gene.",
          "default": "symbol,name,entrezgene,ensembl.gene"
        }
      },
      "required": [
        "gene_ids"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "results": {
          "type": "array",
          "description": "Array of gene annotation results",
          "items": {
            "type": "object",
            "properties": {
              "_id": {
                "type": "string"
              },
              "query": {
                "type": "string",
                "description": "Original query term"
              },
              "symbol": {
                "type": "string"
              },
              "name": {
                "type": "string"
              },
              "entrezgene": {
                "type": [
                  "integer",
                  "string"
                ]
              },
              "notfound": {
                "type": "boolean",
                "description": "True if gene was not found"
              }
            }
          }
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `MyGene_batch_query`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("my_gene_batch_query", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.my_gene_get_gene_annotation` {#dnakit.databases.my_gene_get_gene_annotation}

Get detailed annotation for a specific gene by its ID (Entrez or Ensembl). Returns comprehensive information including gene ontology, pathways, protein domains, and literature references.

```text
dnakit.databases.my_gene_get_gene_annotation(*, gene_id, fields=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填（输入 Schema） | Gene ID to query. Can be Entrez Gene ID (e.g., &#x27;1017&#x27; for CDK2) or Ensembl ID (e.g., &#x27;ENSG00000123374&#x27;). |
| `fields` | string | None | Comma-separated list of fields to return. Available fields include: symbol, name, summary, go (gene ontology), pathway (KEGG/Reactome), interpro (protein domains), generif (gene references), pdb, uniprot. {&quot;default&quot;: &quot;symbol,name,entrezgene,ensembl,summary,go,pathway,interpro&quot;} |
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
          "description": "Gene ID to query. Can be Entrez Gene ID (e.g., '1017' for CDK2) or Ensembl ID (e.g., 'ENSG00000123374')."
        },
        "fields": {
          "type": "string",
          "description": "Comma-separated list of fields to return. Available fields include: symbol, name, summary, go (gene ontology), pathway (KEGG/Reactome), interpro (protein domains), generif (gene references), pdb, uniprot.",
          "default": "symbol,name,entrezgene,ensembl,summary,go,pathway,interpro"
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
        "_id": {
          "type": "string",
          "description": "Gene ID"
        },
        "symbol": {
          "type": "string",
          "description": "Official gene symbol"
        },
        "name": {
          "type": "string",
          "description": "Full gene name"
        },
        "entrezgene": {
          "type": [
            "integer",
            "string"
          ],
          "description": "Entrez Gene ID"
        },
        "ensembl": {
          "type": "object",
          "properties": {
            "gene": {
              "type": "string"
            },
            "transcript": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "protein": {
              "type": "array",
              "items": {
                "type": "string"
              }
            }
          }
        },
        "summary": {
          "type": "string",
          "description": "Gene function summary"
        },
        "go": {
          "type": "object",
          "description": "Gene Ontology annotations",
          "properties": {
            "BP": {
              "type": "array",
              "description": "Biological Process terms"
            },
            "CC": {
              "type": "array",
              "description": "Cellular Component terms"
            },
            "MF": {
              "type": "array",
              "description": "Molecular Function terms"
            }
          }
        },
        "pathway": {
          "type": "object",
          "description": "Pathway annotations from KEGG, Reactome, WikiPathways"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `MyGene_get_gene_annotation`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("my_gene_get_gene_annotation", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.my_gene_query_genes` {#dnakit.databases.my_gene_query_genes}

Search for genes by keyword, symbol, name, or other identifiers. Returns gene annotations from 30+ sources including Entrez Gene, Ensembl, UniProt, HGNC. Useful for finding gene IDs, symbols, descriptions, and basic annotations.

```text
dnakit.databases.my_gene_query_genes(*, query=None, species=None, fields=None, size=None, q=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | None | Search query. Can be gene symbol (e.g., &#x27;CDK2&#x27;), name (&#x27;cyclin dependent kinase&#x27;), Entrez ID (&#x27;1017&#x27;), or Ensembl ID (&#x27;ENSG00000123374&#x27;). Supports wildcards (*) and boolean operators (AND, OR, NOT). |
| `species` | string | None | Species filter. Use common name or NCBI taxonomy ID. Examples: &#x27;human&#x27;, &#x27;mouse&#x27;, &#x27;9606&#x27; (human), &#x27;all&#x27;. {&quot;default&quot;: &quot;human&quot;} |
| `fields` | string | None | Comma-separated list of fields to return. Common fields: symbol, name, entrezgene, ensembl.gene, summary, go, pathway, interpro. {&quot;default&quot;: &quot;symbol,name,entrezgene,ensembl.gene,summary&quot;} |
| `size` | integer | None | Maximum number of results to return (1-100). {&quot;default&quot;: 10} |
| `q` | string | None | DNAKit optional alias for query |
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
          "description": "Search query. Can be gene symbol (e.g., 'CDK2'), name ('cyclin dependent kinase'), Entrez ID ('1017'), or Ensembl ID ('ENSG00000123374'). Supports wildcards (*) and boolean operators (AND, OR, NOT)."
        },
        "species": {
          "type": "string",
          "description": "Species filter. Use common name or NCBI taxonomy ID. Examples: 'human', 'mouse', '9606' (human), 'all'.",
          "default": "human"
        },
        "fields": {
          "type": "string",
          "description": "Comma-separated list of fields to return. Common fields: symbol, name, entrezgene, ensembl.gene, summary, go, pathway, interpro.",
          "default": "symbol,name,entrezgene,ensembl.gene,summary"
        },
        "size": {
          "type": "integer",
          "description": "Maximum number of results to return (1-100).",
          "default": 10
        },
        "q": {
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
          "type": "integer",
          "description": "Query time in milliseconds"
        },
        "total": {
          "type": "integer",
          "description": "Total number of matching genes"
        },
        "max_score": {
          "type": "number",
          "description": "Maximum relevance score"
        },
        "hits": {
          "type": "array",
          "description": "Array of matching gene records",
          "items": {
            "type": "object",
            "properties": {
              "_id": {
                "type": "string",
                "description": "Gene ID (usually Entrez gene ID)"
              },
              "_score": {
                "type": "number",
                "description": "Relevance score"
              },
              "symbol": {
                "type": "string",
                "description": "Official gene symbol"
              },
              "name": {
                "type": "string",
                "description": "Full gene name"
              },
              "entrezgene": {
                "type": [
                  "integer",
                  "string"
                ],
                "description": "Entrez Gene ID"
              },
              "ensembl": {
                "type": "object",
                "properties": {
                  "gene": {
                    "type": "string",
                    "description": "Ensembl gene ID"
                  }
                }
              },
              "summary": {
                "type": "string",
                "description": "Gene summary/description"
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
        "q": "query"
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `MyGene_query_genes`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("my_gene_query_genes", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ncbi_fetch_accessions` {#dnakit.databases.ncbi_fetch_accessions}

Convert GenBank UIDs to accession numbers (U00096, NC_045512, etc.). Takes UIDs from NCBI_search_nucleotide and returns actual accession numbers that can be used with database tools. Returns both EMBL/GenBank accessions (U*, M*, AJ*) that work with ENA Browser tools, and RefSeq accessions (NC_*, NM_*) for NCBI tools. Use this as second step after NCBI_search_nucleotide to get usable accession numbers.

```text
dnakit.databases.ncbi_fetch_accessions(*, operation=None, uids, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | object &#124; None | None | Operation type (fixed: fetch_accession) {&quot;const&quot;: &quot;fetch_accession&quot;, &quot;default&quot;: &quot;fetch_accession&quot;} |
| `uids` | object | 必填（输入 Schema） | GenBank UID(s) from NCBI_search_nucleotide. Can be single UID or array of UIDs. {&quot;oneOf&quot;: &#91;{&quot;type&quot;: &quot;array&quot;, &quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}}, {&quot;type&quot;: &quot;string&quot;}&#93;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "required": [
        "uids"
      ],
      "properties": {
        "operation": {
          "const": "fetch_accession",
          "description": "Operation type (fixed: fetch_accession)",
          "default": "fetch_accession"
        },
        "uids": {
          "oneOf": [
            {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            {
              "type": "string"
            }
          ],
          "description": "GenBank UID(s) from NCBI_search_nucleotide. Can be single UID or array of UIDs."
        }
      }
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "array",
      "items": {
        "type": "string"
      },
      "description": "List of accession numbers"
    }
    ```

**执行与来源：** `network`；上游接口名 `NCBI_fetch_accessions`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ncbi_fetch_accessions", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ncbi_get_sequence` {#dnakit.databases.ncbi_get_sequence}

Retrieve DNA/RNA sequence data from NCBI by accession number. Returns sequences in specified format (fasta, gb, embl). Works with both EMBL/GenBank accessions (U00096, M15390) and RefSeq accessions (NC_045512, NM_000546). Use this to get actual sequence data after finding accessions with NCBI_search_nucleotide. Returns full sequence with annotations depending on format.

```text
dnakit.databases.ncbi_get_sequence(*, operation=None, accession, format=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | object &#124; None | None | Operation type (fixed: fetch_sequence) {&quot;const&quot;: &quot;fetch_sequence&quot;, &quot;default&quot;: &quot;fetch_sequence&quot;} |
| `accession` | string | 必填（输入 Schema） | NCBI accession number (e.g., &#x27;U00096&#x27;, &#x27;NC_045512&#x27;, &#x27;NM_000546&#x27;). Works with GenBank, EMBL, and RefSeq accessions. |
| `format` | string | None | Sequence format: &#x27;fasta&#x27; for FASTA, &#x27;gb&#x27; for GenBank, &#x27;embl&#x27; for EMBL, &#x27;gp&#x27; for GenPept protein. Default: fasta {&quot;enum&quot;: &#91;&quot;fasta&quot;, &quot;gb&quot;, &quot;embl&quot;, &quot;gp&quot;&#93;, &quot;default&quot;: &quot;fasta&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "required": [
        "accession"
      ],
      "properties": {
        "operation": {
          "const": "fetch_sequence",
          "description": "Operation type (fixed: fetch_sequence)",
          "default": "fetch_sequence"
        },
        "accession": {
          "type": "string",
          "description": "NCBI accession number (e.g., 'U00096', 'NC_045512', 'NM_000546'). Works with GenBank, EMBL, and RefSeq accessions."
        },
        "format": {
          "type": "string",
          "enum": [
            "fasta",
            "gb",
            "embl",
            "gp"
          ],
          "description": "Sequence format: 'fasta' for FASTA, 'gb' for GenBank, 'embl' for EMBL, 'gp' for GenPept protein. Default: fasta",
          "default": "fasta"
        }
      }
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "string",
      "description": "Sequence data in requested format (FASTA, GenBank, EMBL, or GenPept)"
    }
    ```

**执行与来源：** `network`；上游接口名 `NCBI_get_sequence`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ncbi_get_sequence", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ncbi_search_nucleotide` {#dnakit.databases.ncbi_search_nucleotide}

Search NCBI Nucleotide database (GenBank/EMBL/RefSeq) by organism, gene name, or keywords. Returns GenBank UIDs which can be used to fetch accession numbers. Supports filtering by organism, gene, strain, sequence type (complete genome, mRNA, RefSeq). Use this to discover DNA/RNA sequences for specific organisms or genes. Returns both EMBL/GenBank accessions (U00096, M15390) and RefSeq accessions (NC_*, NM_*). Follow up with NCBI_fetch_accessions to get actual accession numbers.

```text
dnakit.databases.ncbi_search_nucleotide(*, operation=None, organism=None, gene=None, strain=None, keywords=None, seq_type=None, query=None, limit=None, sort=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | object &#124; None | None | Operation type (fixed: search) {&quot;const&quot;: &quot;search&quot;, &quot;default&quot;: &quot;search&quot;} |
| `organism` | string | None | Organism name to search (e.g., &#x27;Escherichia coli&#x27;, &#x27;Homo sapiens&#x27;). Searches in &#91;Organism&#93; field. |
| `gene` | string | None | Gene name to search (e.g., &#x27;BRCA1&#x27;, &#x27;TP53&#x27;). Searches in &#91;Gene&#93; field. |
| `strain` | string | None | Strain name to filter (e.g., &#x27;K-12&#x27;, &#x27;MG1655&#x27;). Searches in &#91;Strain&#93; field. |
| `keywords` | string | None | Keywords to search in titles (e.g., &#x27;complete genome&#x27;, &#x27;mitochondrial&#x27;). Searches in &#91;Title&#93; field. |
| `seq_type` | string | None | Sequence type filter: &#x27;complete_genome&#x27; for complete genomes, &#x27;mrna&#x27; for mRNA sequences, &#x27;refseq&#x27; for RefSeq curated sequences {&quot;enum&quot;: &#91;&quot;complete_genome&quot;, &quot;mrna&quot;, &quot;refseq&quot;&#93;} |
| `query` | string | None | Free-form search query if not using structured filters. Use NCBI query syntax (e.g., &#x27;(Escherichia coli&#91;Organism&#93;) AND (K-12&#91;Strain&#93;)&#x27;) |
| `limit` | integer | None | Maximum number of results to return (default: 20, max: 100) {&quot;default&quot;: 20, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `sort` | string | None | Sort order for results (default: relevance) {&quot;enum&quot;: &#91;&quot;relevance&quot;, &quot;pub_date&quot;, &quot;title&quot;&#93;, &quot;default&quot;: &quot;relevance&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "required": [],
      "properties": {
        "operation": {
          "const": "search",
          "description": "Operation type (fixed: search)",
          "default": "search"
        },
        "organism": {
          "type": "string",
          "description": "Organism name to search (e.g., 'Escherichia coli', 'Homo sapiens'). Searches in [Organism] field."
        },
        "gene": {
          "type": "string",
          "description": "Gene name to search (e.g., 'BRCA1', 'TP53'). Searches in [Gene] field."
        },
        "strain": {
          "type": "string",
          "description": "Strain name to filter (e.g., 'K-12', 'MG1655'). Searches in [Strain] field."
        },
        "keywords": {
          "type": "string",
          "description": "Keywords to search in titles (e.g., 'complete genome', 'mitochondrial'). Searches in [Title] field."
        },
        "seq_type": {
          "type": "string",
          "enum": [
            "complete_genome",
            "mrna",
            "refseq"
          ],
          "description": "Sequence type filter: 'complete_genome' for complete genomes, 'mrna' for mRNA sequences, 'refseq' for RefSeq curated sequences"
        },
        "query": {
          "type": "string",
          "description": "Free-form search query if not using structured filters. Use NCBI query syntax (e.g., '(Escherichia coli[Organism]) AND (K-12[Strain])')"
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of results to return (default: 20, max: 100)",
          "default": 20,
          "minimum": 1,
          "maximum": 100
        },
        "sort": {
          "type": "string",
          "enum": [
            "relevance",
            "pub_date",
            "title"
          ],
          "description": "Sort order for results (default: relevance)",
          "default": "relevance"
        }
      }
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "uids": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "List of GenBank UIDs (internal IDs)"
        },
        "count": {
          "type": "integer",
          "description": "Total number of matching results"
        },
        "returned": {
          "type": "integer",
          "description": "Number of UIDs returned"
        },
        "search_term": {
          "type": "string",
          "description": "Generated NCBI search term used"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `NCBI_search_nucleotide`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ncbi_search_nucleotide", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ncbi_sra_search_runs` {#dnakit.databases.ncbi_sra_search_runs}

Search NCBI Sequence Read Archive (SRA) for NGS/RNA-seq sequencing runs by study accession, organism, sequencing strategy (RNA-Seq, WGS, ChIP-Seq), platform (ILLUMINA, OXFORD_NANOPORE), or library source. Returns SRA UIDs which can be used with other SRA tools to retrieve metadata and download links. Essential for discovering high-throughput sequencing data from public repositories. Use this to find relevant sequencing experiments before accessing the actual FASTQ data files.

```text
dnakit.databases.ncbi_sra_search_runs(*, operation=None, study=None, organism=None, strategy=None, platform=None, source=None, query=None, limit=None, sort=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | object &#124; None | None | Operation type (fixed: search) {&quot;const&quot;: &quot;search&quot;} |
| `study` | string | None | Study accession to search (e.g., &#x27;SRP000001&#x27;, &#x27;ERP000001&#x27;, &#x27;DRP000001&#x27;). Searches in &#91;Study&#93; field. |
| `organism` | string | None | Organism name to filter (e.g., &#x27;Homo sapiens&#x27;, &#x27;Mus musculus&#x27;, &#x27;SARS-CoV-2&#x27;). Searches in &#91;Organism&#93; field. |
| `strategy` | string | None | Sequencing strategy to filter. Common values: &#x27;RNA-Seq&#x27; (transcriptomics), &#x27;WGS&#x27; (whole genome), &#x27;WXS&#x27; (exome), &#x27;ChIP-Seq&#x27; (chromatin), &#x27;ATAC-Seq&#x27; (accessibility), &#x27;Bisulfite-Seq&#x27; (methylation), &#x27;AMPLICON&#x27; (targeted), &#x27;Hi-C&#x27; (3D genome) |
| `platform` | string | None | Sequencing platform to filter. Common values: &#x27;ILLUMINA&#x27;, &#x27;OXFORD_NANOPORE&#x27;, &#x27;PACBIO_SMRT&#x27;, &#x27;ION_TORRENT&#x27;, &#x27;BGISEQ&#x27; |
| `source` | string | None | Library source to filter. Common values: &#x27;GENOMIC&#x27; (DNA), &#x27;TRANSCRIPTOMIC&#x27; (RNA), &#x27;METAGENOMIC&#x27;, &#x27;METATRANSCRIPTOMIC&#x27;, &#x27;SYNTHETIC&#x27;, &#x27;VIRAL_RNA&#x27; |
| `query` | string | None | Free-form search query if not using structured filters. Use NCBI query syntax (e.g., &#x27;(Homo sapiens&#91;Organism&#93;) AND (RNA-Seq&#91;Strategy&#93;)&#x27;) |
| `limit` | integer | None | Maximum number of results to return (default: 20, max: 100) {&quot;default&quot;: 20, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `sort` | string | None | Sort order for results (default: relevance) {&quot;enum&quot;: &#91;&quot;relevance&quot;, &quot;pub_date&quot;&#93;, &quot;default&quot;: &quot;relevance&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "required": [],
      "properties": {
        "operation": {
          "const": "search",
          "description": "Operation type (fixed: search)"
        },
        "study": {
          "type": "string",
          "description": "Study accession to search (e.g., 'SRP000001', 'ERP000001', 'DRP000001'). Searches in [Study] field."
        },
        "organism": {
          "type": "string",
          "description": "Organism name to filter (e.g., 'Homo sapiens', 'Mus musculus', 'SARS-CoV-2'). Searches in [Organism] field."
        },
        "strategy": {
          "type": "string",
          "description": "Sequencing strategy to filter. Common values: 'RNA-Seq' (transcriptomics), 'WGS' (whole genome), 'WXS' (exome), 'ChIP-Seq' (chromatin), 'ATAC-Seq' (accessibility), 'Bisulfite-Seq' (methylation), 'AMPLICON' (targeted), 'Hi-C' (3D genome)"
        },
        "platform": {
          "type": "string",
          "description": "Sequencing platform to filter. Common values: 'ILLUMINA', 'OXFORD_NANOPORE', 'PACBIO_SMRT', 'ION_TORRENT', 'BGISEQ'"
        },
        "source": {
          "type": "string",
          "description": "Library source to filter. Common values: 'GENOMIC' (DNA), 'TRANSCRIPTOMIC' (RNA), 'METAGENOMIC', 'METATRANSCRIPTOMIC', 'SYNTHETIC', 'VIRAL_RNA'"
        },
        "query": {
          "type": "string",
          "description": "Free-form search query if not using structured filters. Use NCBI query syntax (e.g., '(Homo sapiens[Organism]) AND (RNA-Seq[Strategy])')"
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of results to return (default: 20, max: 100)",
          "default": 20,
          "minimum": 1,
          "maximum": 100
        },
        "sort": {
          "type": "string",
          "enum": [
            "relevance",
            "pub_date"
          ],
          "description": "Sort order for results (default: relevance)",
          "default": "relevance"
        }
      }
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "properties": {
        "uids": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "List of SRA UIDs"
        },
        "count": {
          "type": "integer",
          "description": "Total number of matching results"
        },
        "returned": {
          "type": "integer",
          "description": "Number of UIDs returned"
        },
        "search_term": {
          "type": "string",
          "description": "Generated NCBI search term used"
        }
      }
    }
    ```

**执行与来源：** `network`；上游接口名 `NCBI_SRA_search_runs`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ncbi_sra_search_runs", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ncbidatasets_get_gene` {#dnakit.databases.ncbidatasets_get_gene}

Get comprehensive gene information from the National Center for Biotechnology Information (NCBI) by Gene ID. Returns gene symbol, description, chromosomal location, genomic coordinates, cross-references to UniProt/Ensembl/OMIM, synonyms, and annotation details. Covers all organisms in the NCBI Gene database. No API key is required; optionally set NCBI_API_KEY for a higher request-rate limit. Example: gene_id 7157 returns TP53 (tumor protein p53) on chromosome 17.

```text
dnakit.databases.ncbidatasets_get_gene(*, gene_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填（输入 Schema） | NCBI Gene ID (numeric). Examples: &#x27;7157&#x27; (TP53), &#x27;672&#x27; (BRCA1), &#x27;3630&#x27; (INS), &#x27;22059&#x27; (mouse Trp53). {&quot;pattern&quot;: &quot;^&#91;1-9&#93;&#91;0-9&#93;*$&quot;} |
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
          "pattern": "^[1-9][0-9]*$",
          "description": "NCBI Gene ID (numeric). Examples: '7157' (TP53), '672' (BRCA1), '3630' (INS), '22059' (mouse Trp53)."
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
              "type": [
                "string",
                "null"
              ],
              "description": "NCBI Gene ID"
            },
            "symbol": {
              "type": [
                "string",
                "null"
              ],
              "description": "Gene symbol"
            },
            "description": {
              "type": [
                "string",
                "null"
              ],
              "description": "Gene description/full name"
            },
            "tax_id": {
              "type": [
                "string",
                "null"
              ],
              "description": "NCBI Taxonomy ID"
            },
            "taxname": {
              "type": [
                "string",
                "null"
              ],
              "description": "Species scientific name"
            },
            "common_name": {
              "type": [
                "string",
                "null"
              ],
              "description": "Species common name"
            },
            "type": {
              "type": [
                "string",
                "null"
              ],
              "description": "Gene biotype (e.g., PROTEIN_CODING)"
            },
            "chromosomes": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "orientation": {
              "type": [
                "string",
                "null"
              ]
            },
            "swiss_prot_accessions": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "ensembl_gene_ids": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "omim_ids": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "synonyms": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "nomenclature_authority": {
              "type": [
                "object",
                "null"
              ]
            },
            "genomic_locations": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "assembly": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "accession": {
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
                  "begin": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "end": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "orientation": {
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

**执行与来源：** `network`；上游接口名 `NCBIDatasets_get_gene`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ncbidatasets_get_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ncbidatasets_get_gene_by_symbol` {#dnakit.databases.ncbidatasets_get_gene_by_symbol}

Look up gene information by gene symbol and organism. Searches NCBI Gene database using official gene symbols. Specify the species with `taxon` (or its synonyms `organism`/`species`) as a common name (human, mouse, rat) or taxonomy ID. IMPORTANT -- the species defaults to HUMAN when none is given, so a gene symbol shared across species returns the human gene unless you say otherwise; `metadata.query_taxon_defaulted` is true whenever that default was applied, and `metadata.query_taxon` echoes the species actually used. Every returned row also carries its own `tax_id`/`taxname`, which is the authoritative check. Returns gene ID, description, chromosomal location, and database cross-references.

```text
dnakit.databases.ncbidatasets_get_gene_by_symbol(*, symbol, taxon=None, organism=None, species=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `symbol` | string | 必填（输入 Schema） | Gene symbol. Examples: &#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;INS&#x27;, &#x27;EGFR&#x27;. {&quot;minLength&quot;: 1, &quot;pattern&quot;: &quot;.*\\S.*&quot;} |
| `taxon` | string | None | Organism as common name or taxonomy ID. Examples: &#x27;human&#x27;, &#x27;mouse&#x27;, &#x27;rat&#x27;, &#x27;9606&#x27;, &#x27;10090&#x27;. Default: &#x27;human&#x27; (reported as metadata.query_taxon_defaulted). Synonyms `organism` and `species` are accepted. {&quot;minLength&quot;: 1, &quot;pattern&quot;: &quot;.*\\S.*&quot;} |
| `organism` | string | None | Synonym for `taxon`. Used when `taxon` is not supplied. {&quot;minLength&quot;: 1, &quot;pattern&quot;: &quot;.*\\S.*&quot;} |
| `species` | string | None | Synonym for `taxon`. Used when neither `taxon` nor `organism` is supplied. {&quot;minLength&quot;: 1, &quot;pattern&quot;: &quot;.*\\S.*&quot;} |
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
          "minLength": 1,
          "pattern": ".*\\S.*",
          "description": "Gene symbol. Examples: 'TP53', 'BRCA1', 'INS', 'EGFR'."
        },
        "taxon": {
          "type": "string",
          "minLength": 1,
          "pattern": ".*\\S.*",
          "description": "Organism as common name or taxonomy ID. Examples: 'human', 'mouse', 'rat', '9606', '10090'. Default: 'human' (reported as metadata.query_taxon_defaulted). Synonyms `organism` and `species` are accepted."
        },
        "organism": {
          "type": "string",
          "minLength": 1,
          "pattern": ".*\\S.*",
          "description": "Synonym for `taxon`. Used when `taxon` is not supplied."
        },
        "species": {
          "type": "string",
          "minLength": 1,
          "pattern": ".*\\S.*",
          "description": "Synonym for `taxon`. Used when neither `taxon` nor `organism` is supplied."
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
              "description": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "tax_id": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "taxname": {
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
              "chromosomes": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "swiss_prot_accessions": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "ensembl_gene_ids": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "synonyms": {
                "type": "array",
                "items": {
                  "type": "string"
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

**执行与来源：** `network`；上游接口名 `NCBIDatasets_get_gene_by_symbol`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ncbidatasets_get_gene_by_symbol", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ncbidatasets_get_genome_assembly` {#dnakit.databases.ncbidatasets_get_genome_assembly}

Get genome assembly metadata for an assembly accession (RefSeq GCF_ or GenBank GCA_) from the NCBI Datasets v2 API. Returns organism/strain, assembly name/level/status, RefSeq category (e.g. reference genome), release date, submitter, BioProject, and assembly statistics (total length, chromosome count, contig/scaffold N50, GC%, annotation provider). Use to characterize a specific genome assembly across any organism (bacteria, fungi, plant, animal). No API key required.

```text
dnakit.databases.ncbidatasets_get_genome_assembly(*, accession, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string | 必填（输入 Schema） | Assembly accession, e.g. &#x27;GCF_000005845.2&#x27; (E. coli K-12) or &#x27;GCF_000001405.40&#x27; (human GRCh38.p14). {&quot;pattern&quot;: &quot;^(GCF&#124;GCA)_&#91;0-9&#93;+(?:\\.&#91;0-9&#93;+)?$&quot;} |
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
          "pattern": "^(GCF|GCA)_[0-9]+(?:\\.[0-9]+)?$",
          "description": "Assembly accession, e.g. 'GCF_000005845.2' (E. coli K-12) or 'GCF_000001405.40' (human GRCh38.p14)."
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
          "description": "Curated assembly fields (accession, organism, level, stats)",
          "properties": {
            "accession": {
              "type": [
                "string",
                "null"
              ]
            },
            "paired_accession": {
              "type": [
                "string",
                "null"
              ]
            },
            "source_database": {
              "type": [
                "string",
                "null"
              ]
            },
            "organism_name": {
              "type": [
                "string",
                "null"
              ]
            },
            "tax_id": {
              "type": [
                "integer",
                "null"
              ]
            },
            "strain": {
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
            "assembly_level": {
              "type": [
                "string",
                "null"
              ]
            },
            "assembly_status": {
              "type": [
                "string",
                "null"
              ]
            },
            "refseq_category": {
              "type": [
                "string",
                "null"
              ]
            },
            "release_date": {
              "type": [
                "string",
                "null"
              ]
            },
            "submitter": {
              "type": [
                "string",
                "null"
              ]
            },
            "bioproject_accession": {
              "type": [
                "string",
                "null"
              ]
            },
            "total_sequence_length": {
              "type": [
                "string",
                "integer",
                "null"
              ]
            },
            "number_of_chromosomes": {
              "type": [
                "integer",
                "null"
              ]
            },
            "number_of_contigs": {
              "type": [
                "integer",
                "null"
              ]
            },
            "contig_n50": {
              "type": [
                "integer",
                "null"
              ]
            },
            "scaffold_n50": {
              "type": [
                "integer",
                "null"
              ]
            },
            "gc_percent": {
              "type": [
                "number",
                "null"
              ]
            },
            "annotation_provider": {
              "type": [
                "string",
                "null"
              ]
            },
            "annotation_name": {
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

**执行与来源：** `network`；上游接口名 `NCBIDatasets_get_genome_assembly`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ncbidatasets_get_genome_assembly", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ncbidatasets_get_orthologs` {#dnakit.databases.ncbidatasets_get_orthologs}

Get orthologous genes across species for a given NCBI Gene ID. Returns orthologs identified by NCBI's Ortholog pipeline, including gene symbols, species, and gene types for each ortholog. Useful for comparative genomics and cross-species gene function studies. Example: gene_id 7157 (TP53) returns mouse Trp53, rat Tp53, zebrafish tp53, etc.

```text
dnakit.databases.ncbidatasets_get_orthologs(*, gene_id, page_size=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填（输入 Schema） | NCBI Gene ID (numeric). Examples: &#x27;7157&#x27; (TP53), &#x27;672&#x27; (BRCA1), &#x27;3630&#x27; (INS). {&quot;pattern&quot;: &quot;^&#91;1-9&#93;&#91;0-9&#93;*$&quot;} |
| `page_size` | integer | None | Maximum number of orthologs to return (1-100). Default: 20. {&quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
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
          "pattern": "^[1-9][0-9]*$",
          "description": "NCBI Gene ID (numeric). Examples: '7157' (TP53), '672' (BRCA1), '3630' (INS)."
        },
        "page_size": {
          "type": "integer",
          "minimum": 1,
          "maximum": 100,
          "description": "Maximum number of orthologs to return (1-100). Default: 20."
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
              "description": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "tax_id": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "taxname": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "common_name": {
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
              "chromosomes": {
                "type": "array",
                "items": {
                  "type": "string"
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

**执行与来源：** `network`；上游接口名 `NCBIDatasets_get_orthologs`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ncbidatasets_get_orthologs", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ncbidatasets_get_sequence_reports` {#dnakit.databases.ncbidatasets_get_sequence_reports}

Get a page of per-sequence (chromosome/plasmid/scaffold) reports for a genome assembly via the NCBI Datasets v2 API. Returns each molecule's chromosome name, role, location type, RefSeq and GenBank accessions, length, and GC%, plus total_available and next_page_token metadata. Use next_page_token to retrieve additional pages for assemblies with many sequences. Use to map an assembly to its individual sequence accessions (e.g. resolve E. coli K-12 to NC_000913.3, or human chromosomes to their RefSeq NC_ accessions). No API key required.

```text
dnakit.databases.ncbidatasets_get_sequence_reports(*, accession, page_size=None, page_token=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string | 必填（输入 Schema） | Assembly accession, e.g. &#x27;GCF_000005845.2&#x27;. {&quot;pattern&quot;: &quot;^(GCF&#124;GCA)_&#91;0-9&#93;+(?:\\.&#91;0-9&#93;+)?$&quot;} |
| `page_size` | integer &#124; null | None | Maximum sequence records in this page (default 100, maximum 1000). {&quot;minimum&quot;: 1, &quot;maximum&quot;: 1000} |
| `page_token` | string &#124; null | None | Opaque next_page_token returned by a previous call. Omit for the first page. {&quot;minLength&quot;: 1, &quot;pattern&quot;: &quot;.*\\S.*&quot;} |
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
          "pattern": "^(GCF|GCA)_[0-9]+(?:\\.[0-9]+)?$",
          "description": "Assembly accession, e.g. 'GCF_000005845.2'."
        },
        "page_size": {
          "type": [
            "integer",
            "null"
          ],
          "minimum": 1,
          "maximum": 1000,
          "description": "Maximum sequence records in this page (default 100, maximum 1000)."
        },
        "page_token": {
          "type": [
            "string",
            "null"
          ],
          "minLength": 1,
          "pattern": ".*\\S.*",
          "description": "Opaque next_page_token returned by a previous call. Omit for the first page."
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
              "chr_name": {
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
              "role": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "location_type": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "refseq_accession": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "genbank_accession": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "length": {
                "type": [
                  "integer",
                  "null"
                ]
              },
              "gc_percent": {
                "type": [
                  "number",
                  "null"
                ]
              }
            }
          },
          "description": "Sequence records"
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

**执行与来源：** `network`；上游接口名 `NCBIDatasets_get_sequence_reports`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ncbidatasets_get_sequence_reports", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ncbidatasets_get_taxonomy` {#dnakit.databases.ncbidatasets_get_taxonomy}

Get detailed taxonomy information from NCBI by taxonomy ID. Returns organism name, rank, full lineage (as bare ancestor tax_ids in `lineage`, and as resolved organism_name/rank pairs in `lineage_names`), current counts of genes/assemblies/RNA types, and child taxa. Covers all organisms in NCBI Taxonomy including bacteria, archaea, eukaryotes, and viruses. Example: tax_id 9606 returns Homo sapiens and its root-to-parent lineage.

```text
dnakit.databases.ncbidatasets_get_taxonomy(*, tax_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `tax_id` | string | 必填（输入 Schema） | NCBI Taxonomy ID (numeric). Examples: &#x27;9606&#x27; (Homo sapiens), &#x27;10090&#x27; (Mus musculus), &#x27;7227&#x27; (Drosophila melanogaster), &#x27;562&#x27; (E. coli). {&quot;pattern&quot;: &quot;^&#91;1-9&#93;&#91;0-9&#93;*$&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "tax_id": {
          "type": "string",
          "pattern": "^[1-9][0-9]*$",
          "description": "NCBI Taxonomy ID (numeric). Examples: '9606' (Homo sapiens), '10090' (Mus musculus), '7227' (Drosophila melanogaster), '562' (E. coli)."
        }
      },
      "required": [
        "tax_id"
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
            "tax_id": {
              "type": [
                "integer",
                "null"
              ]
            },
            "organism_name": {
              "type": [
                "string",
                "null"
              ]
            },
            "genbank_common_name": {
              "type": [
                "string",
                "null"
              ]
            },
            "rank": {
              "type": [
                "string",
                "null"
              ]
            },
            "blast_name": {
              "type": [
                "string",
                "null"
              ]
            },
            "lineage": {
              "type": "array",
              "items": {
                "type": "integer"
              }
            },
            "lineage_names": {
              "type": "array",
              "description": "Ancestor tax_ids from `lineage`, resolved to organism_name/rank (root-first). Entries fall back to null names/ranks if the resolution lookup failed.",
              "items": {
                "type": "object",
                "properties": {
                  "tax_id": {
                    "type": "integer"
                  },
                  "organism_name": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "rank": {
                    "type": [
                      "string",
                      "null"
                    ]
                  }
                }
              }
            },
            "children": {
              "type": "array",
              "items": {
                "type": "integer"
              }
            },
            "counts": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "type": {
                    "type": "string"
                  },
                  "count": {
                    "type": "integer"
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

**执行与来源：** `network`；上游接口名 `NCBIDatasets_get_taxonomy`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ncbidatasets_get_taxonomy", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ncbidatasets_list_genomes_by_taxon` {#dnakit.databases.ncbidatasets_list_genomes_by_taxon}

List genome assemblies available for a taxon (NCBI tax id or scientific name) via the NCBI Datasets v2 API. Returns a summary list of assemblies (accession, assembly name/level, RefSeq category, release date, total length, N50, GC%) plus the total count available. Use to discover what genomes exist for an organism or clade, e.g. all Escherichia coli or Mycobacterium tuberculosis assemblies. Set reference_only to restrict to reference/representative genomes. No API key required.

```text
dnakit.databases.ncbidatasets_list_genomes_by_taxon(*, taxon, limit=None, reference_only=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `taxon` | string | 必填（输入 Schema） | NCBI tax id or scientific name, e.g. &#x27;562&#x27;, &#x27;Escherichia coli&#x27;, &#x27;Saccharomyces cerevisiae&#x27;. {&quot;minLength&quot;: 1, &quot;pattern&quot;: &quot;.*\\S.*&quot;} |
| `limit` | integer &#124; null | None | Max assemblies to return (default 20, max 100). {&quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `reference_only` | boolean &#124; null | None | If true, return only reference/representative genomes. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "taxon": {
          "type": "string",
          "minLength": 1,
          "pattern": ".*\\S.*",
          "description": "NCBI tax id or scientific name, e.g. '562', 'Escherichia coli', 'Saccharomyces cerevisiae'."
        },
        "limit": {
          "type": [
            "integer",
            "null"
          ],
          "minimum": 1,
          "maximum": 100,
          "description": "Max assemblies to return (default 20, max 100)."
        },
        "reference_only": {
          "type": [
            "boolean",
            "null"
          ],
          "description": "If true, return only reference/representative genomes."
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
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "accession": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "paired_accession": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "source_database": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "organism_name": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "tax_id": {
                "type": [
                  "integer",
                  "null"
                ]
              },
              "strain": {
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
              "assembly_level": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "assembly_status": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "refseq_category": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "release_date": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "submitter": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "bioproject_accession": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "total_sequence_length": {
                "type": [
                  "string",
                  "integer",
                  "null"
                ]
              },
              "number_of_chromosomes": {
                "type": [
                  "integer",
                  "null"
                ]
              },
              "number_of_contigs": {
                "type": [
                  "integer",
                  "null"
                ]
              },
              "contig_n50": {
                "type": [
                  "integer",
                  "null"
                ]
              },
              "scaffold_n50": {
                "type": [
                  "integer",
                  "null"
                ]
              },
              "gc_percent": {
                "type": [
                  "number",
                  "null"
                ]
              },
              "annotation_provider": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "annotation_name": {
                "type": [
                  "string",
                  "null"
                ]
              }
            }
          },
          "description": "Assembly summaries"
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

**执行与来源：** `network`；上游接口名 `NCBIDatasets_list_genomes_by_taxon`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ncbidatasets_list_genomes_by_taxon", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ncbidatasets_suggest_taxonomy` {#dnakit.databases.ncbidatasets_suggest_taxonomy}

Suggest taxonomic names matching a query string. Searches NCBI Taxonomy for organisms by partial name, returning matching scientific names, taxonomy IDs, common names, and ranks. Useful for finding taxonomy IDs from partial organism names. Example: 'drosophila' returns D. melanogaster (7227), D. simulans (7240), etc.

```text
dnakit.databases.ncbidatasets_suggest_taxonomy(*, query, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填（输入 Schema） | Partial organism name to search. Examples: &#x27;drosophila&#x27;, &#x27;escherichia&#x27;, &#x27;arabidopsis&#x27;, &#x27;saccharomyces&#x27;. {&quot;minLength&quot;: 1, &quot;pattern&quot;: &quot;.*\\S.*&quot;} |
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
          "minLength": 1,
          "pattern": ".*\\S.*",
          "description": "Partial organism name to search. Examples: 'drosophila', 'escherichia', 'arabidopsis', 'saccharomyces'."
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
              "scientific_name": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "tax_id": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "common_name": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "rank": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "group_name": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "matched_term": {
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

**执行与来源：** `network`；上游接口名 `NCBIDatasets_suggest_taxonomy`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ncbidatasets_suggest_taxonomy", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ncbigene_get_summary` {#dnakit.databases.ncbigene_get_summary}

Get a detailed gene summary from NCBI Gene by Entrez Gene ID using the E-utilities ESummary API. Returns gene name, symbol, description, function summary, chromosomal location, map position, aliases, genomic coordinates, exon count, and MIM (OMIM) IDs. The most authoritative source for human gene information. Use NCBIGene_search first to find Gene IDs. Example: Gene ID 7157 returns TP53 on chromosome 17p13.1 with a detailed functional summary.

```text
dnakit.databases.ncbigene_get_summary(*, id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `id` | string | 必填（输入 Schema） | NCBI Gene ID (Entrez Gene ID). Can be a single ID or comma-separated list. Examples: &#x27;7157&#x27; (TP53), &#x27;672&#x27; (BRCA1), &#x27;1956&#x27; (EGFR), &#x27;3845&#x27; (KRAS). Get IDs from NCBIGene_search. |
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
          "description": "NCBI Gene ID (Entrez Gene ID). Can be a single ID or comma-separated list. Examples: '7157' (TP53), '672' (BRCA1), '1956' (EGFR), '3845' (KRAS). Get IDs from NCBIGene_search."
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
                  "description": "List of returned Gene IDs"
                }
              },
              "additionalProperties": {
                "type": "object",
                "properties": {
                  "uid": {
                    "type": "string",
                    "description": "Gene ID"
                  },
                  "name": {
                    "type": "string",
                    "description": "Gene symbol (e.g., TP53)"
                  },
                  "description": {
                    "type": "string",
                    "description": "Gene full name"
                  },
                  "chromosome": {
                    "type": "string"
                  },
                  "maplocation": {
                    "type": "string",
                    "description": "Cytogenetic band (e.g., 17p13.1)"
                  },
                  "otheraliases": {
                    "type": "string",
                    "description": "Comma-separated gene aliases"
                  },
                  "otherdesignations": {
                    "type": "string",
                    "description": "Alternative protein names"
                  },
                  "nomenclaturesymbol": {
                    "type": "string",
                    "description": "Official gene symbol"
                  },
                  "nomenclaturename": {
                    "type": "string",
                    "description": "Official gene name"
                  },
                  "summary": {
                    "type": "string",
                    "description": "Detailed gene function summary from RefSeq"
                  },
                  "mim": {
                    "type": "array",
                    "items": {
                      "type": "string"
                    },
                    "description": "OMIM disease/phenotype IDs"
                  },
                  "genomicinfo": {
                    "type": "array",
                    "items": {
                      "type": "object",
                      "properties": {
                        "chrloc": {
                          "type": "string"
                        },
                        "chraccver": {
                          "type": "string",
                          "description": "Chromosome RefSeq accession"
                        },
                        "chrstart": {
                          "type": "integer"
                        },
                        "chrstop": {
                          "type": "integer"
                        },
                        "exoncount": {
                          "type": "integer"
                        }
                      }
                    }
                  },
                  "organism": {
                    "type": "object",
                    "properties": {
                      "scientificname": {
                        "type": "string"
                      },
                      "commonname": {
                        "type": "string"
                      },
                      "taxid": {
                        "type": "integer"
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

**执行与来源：** `network`；上游接口名 `NCBIGene_get_summary`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ncbigene_get_summary", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ncbigene_search` {#dnakit.databases.ncbigene_search}

Search the NCBI Gene database using the E-utilities API. Find gene IDs (Entrez Gene IDs) by gene symbol, name, or keyword combined with an organism filter. Returns a list of matching Gene IDs that can be used with NCBIGene_get_summary. Essential for resolving gene names to standardized NCBI Gene IDs. Example: searching 'TP53' in Homo sapiens returns Gene ID 7157; 'BRCA1' returns 672.

```text
dnakit.databases.ncbigene_search(*, term, retmax=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `term` | string | 必填（输入 Schema） | Search query. Use gene symbol with organism filter for best results. Format: &#x27;GENE&#91;Symbol&#93; AND Organism&#91;Organism&#93;&#x27;. Examples: &#x27;TP53&#91;Symbol&#93; AND Homo sapiens&#91;Organism&#93;&#x27;, &#x27;BRCA1&#91;Symbol&#93; AND Homo sapiens&#91;Organism&#93;&#x27;, &#x27;insulin AND Homo sapiens&#91;Organism&#93;&#x27;, &#x27;EGFR AND Mus musculus&#91;Organism&#93;&#x27;. |
| `retmax` | integer | None | Maximum number of results to return (default 10, max 100). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "term": {
          "type": "string",
          "description": "Search query. Use gene symbol with organism filter for best results. Format: 'GENE[Symbol] AND Organism[Organism]'. Examples: 'TP53[Symbol] AND Homo sapiens[Organism]', 'BRCA1[Symbol] AND Homo sapiens[Organism]', 'insulin AND Homo sapiens[Organism]', 'EGFR AND Mus musculus[Organism]'."
        },
        "retmax": {
          "type": "integer",
          "description": "Maximum number of results to return (default 10, max 100)."
        }
      },
      "required": [
        "term"
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
            "header": {
              "type": "object",
              "properties": {
                "type": {
                  "type": "string"
                },
                "version": {
                  "type": "string"
                }
              }
            },
            "esearchresult": {
              "type": "object",
              "properties": {
                "count": {
                  "type": "string",
                  "description": "Total number of matching genes"
                },
                "retmax": {
                  "type": "string"
                },
                "retstart": {
                  "type": "string"
                },
                "idlist": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  },
                  "description": "List of NCBI Gene IDs (Entrez Gene IDs)"
                },
                "querytranslation": {
                  "type": "string",
                  "description": "How NCBI interpreted the query"
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

**执行与来源：** `network`；上游接口名 `NCBIGene_search`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("ncbigene_search", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.sra_get_experiment` {#dnakit.databases.sra_get_experiment}

Get detailed metadata for a specific SRA experiment by accession number. Returns experiment title, organism, sequencing platform, library information, study details, and run statistics. Accepts SRA experiment accessions (SRX/ERX/DRX), study accessions (SRP/ERP/DRP), or sample accessions (SRS/ERS/DRS).

```text
dnakit.databases.sra_get_experiment(*, accession, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string | 必填（输入 Schema） | SRA accession number (e.g., &#x27;SRX26927637&#x27;, &#x27;SRP000001&#x27;, &#x27;SRS000001&#x27;) |
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
          "description": "SRA accession number (e.g., 'SRX26927637', 'SRP000001', 'SRS000001')"
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
            "data": {
              "type": "object",
              "properties": {
                "uid": {
                  "type": "string"
                },
                "title": {
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
                "taxid": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "platform": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "library_strategy": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "library_source": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "library_selection": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "library_layout": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "total_runs": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "total_spots": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "total_bases": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "experiment_accession": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "study_accession": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "study_name": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "bioproject": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "submitter": {
                  "type": [
                    "string",
                    "null"
                  ]
                },
                "runs": {
                  "type": [
                    "array",
                    "null"
                  ]
                }
              }
            },
            "metadata": {
              "type": "object"
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

**执行与来源：** `network`；上游接口名 `SRA_get_experiment`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("sra_get_experiment", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.sra_search_experiments` {#dnakit.databases.sra_search_experiments}

Search the NCBI Sequence Read Archive (SRA) for sequencing experiments by keyword, organism, library strategy, or sequencing platform. SRA is the largest repository of raw sequencing data, containing data from Illumina, PacBio, Oxford Nanopore, and other platforms. Library strategies include: RNA-Seq, WGS (whole genome), WXS (whole exome), ChIP-Seq, ATAC-seq, Bisulfite-Seq, AMPLICON, scRNA-Seq, and more.

```text
dnakit.databases.sra_search_experiments(*, query=None, organism=None, library_strategy=None, platform=None, limit=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | None | Search keywords (e.g., &#x27;human gut metagenome&#x27;, &#x27;BRCA1 RNA-Seq&#x27;, &#x27;single cell lung cancer&#x27;) |
| `organism` | string | None | Organism name to filter by (e.g., &#x27;Homo sapiens&#x27;, &#x27;Mus musculus&#x27;, &#x27;human gut metagenome&#x27;) |
| `library_strategy` | string | None | Sequencing strategy filter (e.g., &#x27;RNA-Seq&#x27;, &#x27;WGS&#x27;, &#x27;WXS&#x27;, &#x27;ChIP-Seq&#x27;, &#x27;ATAC-seq&#x27;, &#x27;Bisulfite-Seq&#x27;, &#x27;AMPLICON&#x27;) |
| `platform` | string | None | Sequencing platform filter (e.g., &#x27;ILLUMINA&#x27;, &#x27;PACBIO_SMRT&#x27;, &#x27;OXFORD_NANOPORE&#x27;, &#x27;ION_TORRENT&#x27;) |
| `limit` | integer | None | Maximum results to return (1-50, default 10) {&quot;default&quot;: 10} |
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
          "description": "Search keywords (e.g., 'human gut metagenome', 'BRCA1 RNA-Seq', 'single cell lung cancer')"
        },
        "organism": {
          "type": "string",
          "description": "Organism name to filter by (e.g., 'Homo sapiens', 'Mus musculus', 'human gut metagenome')"
        },
        "library_strategy": {
          "type": "string",
          "description": "Sequencing strategy filter (e.g., 'RNA-Seq', 'WGS', 'WXS', 'ChIP-Seq', 'ATAC-seq', 'Bisulfite-Seq', 'AMPLICON')"
        },
        "platform": {
          "type": "string",
          "description": "Sequencing platform filter (e.g., 'ILLUMINA', 'PACBIO_SMRT', 'OXFORD_NANOPORE', 'ION_TORRENT')"
        },
        "limit": {
          "type": "integer",
          "description": "Maximum results to return (1-50, default 10)",
          "default": 10
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
            "data": {
              "type": "object",
              "properties": {
                "total": {
                  "type": "integer"
                },
                "returned": {
                  "type": "integer"
                },
                "query_used": {
                  "type": "string"
                },
                "experiments": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "uid": {
                        "type": "string"
                      },
                      "title": {
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
                      "taxid": {
                        "type": [
                          "string",
                          "null"
                        ]
                      },
                      "platform": {
                        "type": [
                          "string",
                          "null"
                        ]
                      },
                      "library_strategy": {
                        "type": [
                          "string",
                          "null"
                        ]
                      },
                      "library_source": {
                        "type": [
                          "string",
                          "null"
                        ]
                      },
                      "library_layout": {
                        "type": [
                          "string",
                          "null"
                        ]
                      },
                      "total_runs": {
                        "type": [
                          "string",
                          "null"
                        ]
                      },
                      "total_spots": {
                        "type": [
                          "string",
                          "null"
                        ]
                      },
                      "total_bases": {
                        "type": [
                          "string",
                          "null"
                        ]
                      },
                      "experiment_accession": {
                        "type": [
                          "string",
                          "null"
                        ]
                      },
                      "study_accession": {
                        "type": [
                          "string",
                          "null"
                        ]
                      },
                      "bioproject": {
                        "type": [
                          "string",
                          "null"
                        ]
                      },
                      "submitter": {
                        "type": [
                          "string",
                          "null"
                        ]
                      },
                      "runs": {
                        "type": [
                          "array",
                          "null"
                        ]
                      }
                    }
                  }
                }
              }
            },
            "metadata": {
              "type": "object"
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

**执行与来源：** `network`；上游接口名 `SRA_search_experiments`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("sra_search_experiments", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.tark_get_mane_transcripts` {#dnakit.databases.tark_get_mane_transcripts}

Look up the MANE (Matched Annotation from NCBI and EMBL-EBI) transcript(s) for a gene, mapping the Ensembl transcript (ENST) to its equivalent RefSeq transcript (NM_). Returns MANE Select and, where defined, MANE Plus Clinical pairings from the Ensembl Tark transcript archive. Query by HGNC gene symbol, an Ensembl transcript id, or a RefSeq transcript id (versions optional — they are matched on the base id). Use to find the clinically-recommended representative transcript for a gene and to convert between RefSeq and Ensembl transcript namespaces. No API key required.

```text
dnakit.databases.tark_get_mane_transcripts(*, gene=None, ensembl_id=None, refseq_id=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string &#124; null | None | HGNC gene symbol, e.g. &#x27;BRCA2&#x27;, &#x27;CFTR&#x27;. |
| `ensembl_id` | string &#124; null | None | Ensembl transcript id, e.g. &#x27;ENST00000380152&#x27; (version optional). |
| `refseq_id` | string &#124; null | None | RefSeq transcript id, e.g. &#x27;NM_000059&#x27; (version optional). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "gene": {
          "type": [
            "string",
            "null"
          ],
          "description": "HGNC gene symbol, e.g. 'BRCA2', 'CFTR'."
        },
        "ensembl_id": {
          "type": [
            "string",
            "null"
          ],
          "description": "Ensembl transcript id, e.g. 'ENST00000380152' (version optional)."
        },
        "refseq_id": {
          "type": [
            "string",
            "null"
          ],
          "description": "RefSeq transcript id, e.g. 'NM_000059' (version optional)."
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
              "gene": {
                "type": "string"
              },
              "mane_type": {
                "type": "string",
                "description": "MANE SELECT or MANE PLUS CLINICAL"
              },
              "ensembl_transcript": {
                "type": "string",
                "description": "Versioned ENST id"
              },
              "refseq_transcript": {
                "type": "string",
                "description": "Versioned RefSeq NM id"
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

**执行与来源：** `network`；上游接口名 `Tark_get_mane_transcripts`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("tark_get_mane_transcripts", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.tark_get_transcript` {#dnakit.databases.tark_get_transcript}

Get the archived transcript record for an Ensembl transcript id (ENST) from the Ensembl Tark transcript archive: assembly, biotype, genomic location (region/start/end/strand), the transcript sequence checksum, and the Ensembl release versions in which the transcript appears. Use to resolve transcript versions and verify transcript identity/equivalence across releases. No API key required.

```text
dnakit.databases.tark_get_transcript(*, stable_id, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `stable_id` | string | 必填（输入 Schema） | Ensembl transcript stable id, e.g. &#x27;ENST00000380152&#x27; (version optional). |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "stable_id": {
          "type": "string",
          "description": "Ensembl transcript stable id, e.g. 'ENST00000380152' (version optional)."
        }
      },
      "required": [
        "stable_id"
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

**执行与来源：** `network`；上游接口名 `Tark_get_transcript`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("tark_get_transcript", verbose=False)
print(info["input_schema"])
```

## `dnakit.download.ena.sequence` {#dnakit.download.ena.sequence}

按 INSDC accession 下载 ENA 序列文件。

```text
dnakit.download.ena.sequence(accession, output_path, *, format='fasta', config=None, search_config=None, api_base_url='https://www.ebi.ac.uk/ena/browser/api', verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | str | 必填 | 版本化序列 accession；ENA 使用 INSDC 编号，NCBI 也支持 RefSeq。 |
| `output_path` | str &#124; Path | 必填 | 目标序列文件路径。 |
| `format` | Literal&#91;&#x27;fasta&#x27;, &#x27;embl&#x27;, &#x27;genbank&#x27;&#93; | &#x27;fasta&#x27; | fasta、genbank 或 embl；提供方及转换依赖见序列下载教程。 |
| `config` | DownloadConfig &#124; None | None | 查询使用 SearchConfig；下载使用 DownloadConfig。 |
| `search_config` | SearchConfig &#124; None | None | 控制查询超时、返回大小和 NCBI 凭据的 SearchConfig。 |
| `api_base_url` | str | &#x27;https://www.ebi.ac.uk/ena/browser/api&#x27; | 服务 API 根地址；默认值见签名。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |

**返回值**

DatasetDownloadResult：files 为文件记录，metadata 含来源、解析出的 accession 和序列长度，manifest_path 指向来源清单。详细格式支持与坐标约定见[序列下载教程](../features/15_download.md#accession-sequences)。

??? info "补充说明（源码原文）"

    Download an INSDC nucleotide accession from ENA with a checksummed manifest.

    FASTA and EMBL are native responses; GenBank conversion requires optional
    dnakit[bio-formats]. RefSeq accessions (including NZ_) require NCBI instead.
    No fallback silently changes the requested accession or source database.
    SearchConfig bounds the response (20 MB default, 100 MB maximum).

## `dnakit.download.ncbi.sequence` {#dnakit.download.ncbi.sequence}

按 accession 下载 NCBI 序列，可选择格式、核酸子区间和链方向。

```text
dnakit.download.ncbi.sequence(accession, output_path, *, database='nuccore', format='fasta', seq_start=None, seq_stop=None, strand=1, config=None, search_config=None, api_base_url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils', verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | str | 必填 | 版本化序列 accession；ENA 使用 INSDC 编号，NCBI 也支持 RefSeq。 |
| `output_path` | str &#124; Path | 必填 | 目标序列文件路径。 |
| `database` | Literal&#91;&#x27;nuccore&#x27;, &#x27;protein&#x27;&#93; | &#x27;nuccore&#x27; | NCBI 序列数据库 nuccore 或 protein。 |
| `format` | Literal&#91;&#x27;fasta&#x27;, &#x27;genbank&#x27;, &#x27;embl&#x27;&#93; | &#x27;fasta&#x27; | fasta、genbank 或 embl；提供方及转换依赖见序列下载教程。 |
| `seq_start` | int &#124; None | None | NCBI 核酸子区间起点，1-based 闭区间；须同时指定 seq_stop。 |
| `seq_stop` | int &#124; None | None | NCBI 核酸子区间终点，1-based 闭区间；须同时指定 seq_start。 |
| `strand` | Literal&#91;1, 2&#93; | 1 | NCBI 核酸序列链方向，1 正向、2 反向互补。 |
| `config` | DownloadConfig &#124; None | None | 查询使用 SearchConfig；下载使用 DownloadConfig。 |
| `search_config` | SearchConfig &#124; None | None | 控制查询超时、返回大小和 NCBI 凭据的 SearchConfig。 |
| `api_base_url` | str | &#x27;https://eutils.ncbi.nlm.nih.gov/entrez/eutils&#x27; | 服务 API 根地址；默认值见签名。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |

**返回值**

DatasetDownloadResult：files 为文件记录，metadata 含来源、解析出的 accession 和序列长度，manifest_path 指向来源清单。详细格式支持与坐标约定见[序列下载教程](../features/15_download.md#accession-sequences)。

??? info "补充说明（源码原文）"

    Download a nucleotide/RNA/protein accession via NCBI EFetch, with a manifest.

    FASTA and GenBank/GenPept are native responses. Nucleotide EMBL output is
    converted from GenBank using optional ``dnakit[bio-formats]`` (Biopython).
    Optional seq_start/seq_stop use NCBI's 1-based inclusive coordinates;
    strand=2 retrieves the reverse complement. Protein requests use full records.
    SearchConfig controls contact/API credentials and bounded response size
    (20 MB by default, at most 100 MB); DownloadConfig controls output limits
    and overwrite. DownloadConfig.api_key is used when search_config is absent.

## `dnakit.search.ncbi.accessions` {#dnakit.search.ncbi.accessions}

把 NCBI 序列 UID 转为版本化 accession，保留输入顺序和重复项。

```text
dnakit.search.ncbi.accessions(uids, *, database='nuccore', config=None, api_base_url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils', verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `uids` | Sequence&#91;str &#124; int&#93; | 必填 | NCBI 序列 UID 列表，按每批 200 个解析。 |
| `database` | Literal&#91;&#x27;nuccore&#x27;, &#x27;protein&#x27;&#93; | &#x27;nuccore&#x27; | NCBI 序列数据库 nuccore 或 protein。 |
| `config` | SearchConfig &#124; None | None | 查询使用 SearchConfig；下载使用 DownloadConfig。 |
| `api_base_url` | str | &#x27;https://eutils.ncbi.nlm.nih.gov/entrez/eutils&#x27; | 服务 API 根地址；默认值见签名。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |

**返回值**

tuple[str, ...]；空输入返回空元组，缺失 UID 或版本号抛出 QueryError。

??? info "补充说明（源码原文）"

    Resolve sequence UIDs to accession.version in input order, in batches of 200.

    Empty input returns an empty tuple without network access. Missing UIDs or
    accession versions raise QueryError instead of silently dropping results.

