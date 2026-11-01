# GrEBI (`grebi`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **15** 个 DNA/基因相关工具（服务器共列出 33 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [EBISPOT/GrEBI @ `74aef79`](https://github.com/EBISPOT/GrEBI/tree/74aef794aab9b3d0e7a52d01c15c48a86d043e09) |
| 许可证 | Apache-2.0 |
| 启动命令 | — |
| 传输方式 | http |
| 安装 | Select transport='http' for the author-documented development service. Local deployment is a separate Java/API plus graph-database and HPC data-pipeline task; no local startup recipe is claimed here. |
| 副作用 | Sends queries to the EMBL-EBI development knowledge-graph service. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`disease_to_genes`](#disease_to_genes) | `association` | Returns genes associated with a disease ID (e.g. MONDO) from Monarch and GWAS Catalog data in the GrEBI graph, with source datasets. |
| [`gene_to_diseases`](#gene_to_diseases) | `association` | Returns diseases associated with a gene ID from Monarch and GWAS Catalog data in the GrEBI graph, with source datasets. |
| [`gwas_by_cell_type`](#gwas_by_cell_type) | `association` | Returns GWAS Catalog associations (SNP, gene, allele, p-value, study) for traits linked to a Cell Ontology cell type. |
| [`gwas_by_disease_location`](#gwas_by_disease_location) | `association` | Returns GWAS Catalog associations (SNP, gene, allele, p-value, study) for diseases located in a given anatomical location. |
| [`gwas_by_drugs_indicated_for_disease`](#gwas_by_drugs_indicated_for_disease) | `association` | Returns GWAS Catalog associations (SNP, gene, allele, p-value) for traits involving drugs indicated for a disease. |
| [`gwas_by_gene_and_cell_type`](#gwas_by_gene_and_cell_type) | `association` | Returns GWAS Catalog associations mapped to a gene for traits linked to a given cell type. |
| [`gwas_by_gene_and_disease`](#gwas_by_gene_and_disease) | `association` | Returns GWAS Catalog associations mapped to a gene for traits within a given disease context. |
| [`gwas_by_gene_and_location`](#gwas_by_gene_and_location) | `association` | Returns GWAS Catalog associations mapped to a gene for traits linked to an anatomical location. |
| [`gwas_by_location`](#gwas_by_location) | `association` | Returns GWAS Catalog associations (SNP, gene, allele, p-value, study) for traits linked to an anatomical location or its parts. |
| [`gwas_by_pathway`](#gwas_by_pathway) | `association` | Returns GWAS Catalog SNP associations for genes in a Reactome pathway, with alleles, traits and p-values. |
| [`gwas_studies_by_cell_type`](#gwas_studies_by_cell_type) | `association` | Returns GWAS Catalog studies annotated with traits linked to a cell type, with sample sizes, association counts and summary statistics availability. |
| [`gwas_studies_by_disease_location`](#gwas_studies_by_disease_location) | `association` | Returns GWAS Catalog studies for diseases located in an anatomical location, with study metadata. |
| [`gwas_studies_by_drugs_indicated_for_disease`](#gwas_studies_by_drugs_indicated_for_disease) | `association` | Returns GWAS Catalog studies whose traits involve drugs indicated for a disease, with study metadata. |
| [`gwas_trait_to_mouse_models_via_embeddings`](#gwas_trait_to_mouse_models_via_embeddings) | `association` | Finds mouse gene models with phenotypes semantically similar to a GWAS trait whose human orthologs have GWAS associations with that trait. |
| [`mouse_gene_to_opentargets`](#mouse_gene_to_opentargets) | `association` | Returns Open Targets disease associations, with scores, for the human ortholog of a mouse gene. |

## `disease_to_genes` {#disease_to_genes}

Returns genes associated with a disease ID (e.g. MONDO) from Monarch and GWAS Catalog data in the GrEBI graph, with source datasets.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`graph`, `disease_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **关键词：**gene-disease association, mondo, monarch, gwas catalog, disease genes

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `graph` | any | 是 | 可选值：`ebi_monarch_xspecies` |
| `disease_id` | string | 是 | Disease |
| `sortBy` | any | 是 | 可选值：`disease`, `gene`, `from_datasources` |
| `sortDir` | any | 是 | 可选值：`asc`, `desc` |
| `pageNum` | integer | 是 | Page number (0-based) |
| `pageSize` | integer | 是 | Number of results per page |

??? note "上游说明（英文原文）"

    Search for gene disease associations

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'disease_to_genes',
    {'graph': 'ebi_monarch_xspecies',
     'disease_id': '...',
     'sortBy': 'disease',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gene_to_diseases` {#gene_to_diseases}

Returns diseases associated with a gene ID from Monarch and GWAS Catalog data in the GrEBI graph, with source datasets.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`graph`, `gene_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **关键词：**gene-disease association, hgnc, monarch, gwas catalog, disease

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `graph` | any | 是 | 可选值：`ebi_monarch_xspecies` |
| `gene_id` | string | 是 | Gene |
| `sortBy` | any | 是 | 可选值：`disease`, `gene`, `from_datasources` |
| `sortDir` | any | 是 | 可选值：`asc`, `desc` |
| `pageNum` | integer | 是 | Page number (0-based) |
| `pageSize` | integer | 是 | Number of results per page |

??? note "上游说明（英文原文）"

    Search for diseases associated with a gene

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gene_to_diseases',
    {'graph': 'ebi_monarch_xspecies',
     'gene_id': '...',
     'sortBy': 'disease',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gwas_by_cell_type` {#gwas_by_cell_type}

Returns GWAS Catalog associations (SNP, gene, allele, p-value, study) for traits linked to a Cell Ontology cell type.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`graph`, `cell_type_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **关键词：**gwas catalog, snp, cell ontology, trait, p-value, mapped gene

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `graph` | any | 是 | 可选值：`ebi_monarch_xspecies` |
| `cell_type_id` | string | 是 | Cell Type |
| `sortBy` | any | 是 | 可选值：`cell_type`, `trait`, `snp`, `effect_allele`, `gene`, `study_accession`, `p_value`, `pubmed_id`, `reported_trait`, `chromosome`, `position`, `risk_allele_frequency` |
| `sortDir` | any | 是 | 可选值：`asc`, `desc` |
| `pageNum` | integer | 是 | Page number (0-based) |
| `pageSize` | integer | 是 | Number of results per page |

??? note "上游说明（英文原文）"

    Search GWAS Catalog by cell type

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gwas_by_cell_type',
    {'graph': 'ebi_monarch_xspecies',
     'cell_type_id': '...',
     'sortBy': 'cell_type',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gwas_by_disease_location` {#gwas_by_disease_location}

Returns GWAS Catalog associations (SNP, gene, allele, p-value, study) for diseases located in a given anatomical location.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`graph`, `location_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **关键词：**gwas catalog, snp, anatomy, disease location, p-value

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `graph` | any | 是 | 可选值：`ebi_monarch_xspecies`, `impc_x_gwas` |
| `location_id` | string | 是 | Anatomical Location |
| `sortBy` | any | 是 | 可选值：`location`, `disease`, `snp`, `effect_allele`, `gene`, `study_accession`, `p_value`, `pubmed_id`, `reported_trait`, `chromosome`, `position`, `risk_allele_frequency` |
| `sortDir` | any | 是 | 可选值：`asc`, `desc` |
| `pageNum` | integer | 是 | Page number (0-based) |
| `pageSize` | integer | 是 | Number of results per page |

??? note "上游说明（英文原文）"

    Search GWAS Catalog by disease location

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gwas_by_disease_location',
    {'graph': 'ebi_monarch_xspecies',
     'location_id': '...',
     'sortBy': 'location',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gwas_by_drugs_indicated_for_disease` {#gwas_by_drugs_indicated_for_disease}

Returns GWAS Catalog associations (SNP, gene, allele, p-value) for traits involving drugs indicated for a disease.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`graph`, `disease_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **关键词：**gwas catalog, snp, drug response, chebi, pharmacogenomics

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `graph` | any | 是 | 可选值：`ebi_monarch_xspecies` |
| `disease_id` | string | 是 | Disease |
| `sortBy` | any | 是 | 可选值：`disease`, `trait`, `drug`, `snp`, `effect_allele`, `gene`, `study_accession`, `p_value`, `pubmed_id`, `reported_trait`, `chromosome`, `position` |
| `sortDir` | any | 是 | 可选值：`asc`, `desc` |
| `pageNum` | integer | 是 | Page number (0-based) |
| `pageSize` | integer | 是 | Number of results per page |

??? note "上游说明（英文原文）"

    Search GWAS Catalog by drugs or drug classes indicated for a disease

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gwas_by_drugs_indicated_for_disease',
    {'graph': 'ebi_monarch_xspecies',
     'disease_id': '...',
     'sortBy': 'disease',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gwas_by_gene_and_cell_type` {#gwas_by_gene_and_cell_type}

Returns GWAS Catalog associations mapped to a gene for traits linked to a given cell type.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`graph`, `gene_id`, `cell_type_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **关键词：**gwas catalog, gene, snp, cell type, trait

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `graph` | any | 是 | 可选值：`ebi_monarch_xspecies` |
| `gene_id` | string | 是 | Gene |
| `cell_type_id` | string | 是 | Cell Type |
| `sortBy` | any | 是 | 可选值：`cell_type`, `trait`, `snp`, `effect_allele`, `gene`, `study_accession`, `p_value`, `pubmed_id`, `reported_trait`, `chromosome`, `position`, `risk_allele_frequency` |
| `sortDir` | any | 是 | 可选值：`asc`, `desc` |
| `pageNum` | integer | 是 | Page number (0-based) |
| `pageSize` | integer | 是 | Number of results per page |

??? note "上游说明（英文原文）"

    Search GWAS Catalog by gene and cell type

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gwas_by_gene_and_cell_type',
    {'graph': 'ebi_monarch_xspecies',
     'gene_id': '...',
     'cell_type_id': '...',
     'sortBy': 'cell_type',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gwas_by_gene_and_disease` {#gwas_by_gene_and_disease}

Returns GWAS Catalog associations mapped to a gene for traits within a given disease context.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`graph`, `gene_id`, `disease_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **关键词：**gwas catalog, gene, snp, disease, trait

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `graph` | any | 是 | 可选值：`ebi_monarch_xspecies` |
| `gene_id` | string | 是 | Gene |
| `disease_id` | string | 是 | Disease |
| `sortBy` | any | 是 | 可选值：`disease`, `trait`, `snp`, `effect_allele`, `gene`, `study_accession`, `p_value`, `pubmed_id`, `reported_trait`, `chromosome`, `position`, `risk_allele_frequency` |
| `sortDir` | any | 是 | 可选值：`asc`, `desc` |
| `pageNum` | integer | 是 | Page number (0-based) |
| `pageSize` | integer | 是 | Number of results per page |

??? note "上游说明（英文原文）"

    Search GWAS Catalog by gene and disease context

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gwas_by_gene_and_disease',
    {'graph': 'ebi_monarch_xspecies',
     'gene_id': '...',
     'disease_id': '...',
     'sortBy': 'disease',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gwas_by_gene_and_location` {#gwas_by_gene_and_location}

Returns GWAS Catalog associations mapped to a gene for traits linked to an anatomical location.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`graph`, `gene_id`, `location_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **关键词：**gwas catalog, gene, snp, anatomy, trait

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `graph` | any | 是 | 可选值：`ebi_monarch_xspecies` |
| `gene_id` | string | 是 | Gene |
| `location_id` | string | 是 | Anatomical Location |
| `sortBy` | any | 是 | 可选值：`location`, `trait`, `snp`, `effect_allele`, `gene`, `study_accession`, `p_value`, `pubmed_id`, `reported_trait`, `chromosome`, `position`, `risk_allele_frequency` |
| `sortDir` | any | 是 | 可选值：`asc`, `desc` |
| `pageNum` | integer | 是 | Page number (0-based) |
| `pageSize` | integer | 是 | Number of results per page |

??? note "上游说明（英文原文）"

    Search GWAS Catalog by gene and anatomical location

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gwas_by_gene_and_location',
    {'graph': 'ebi_monarch_xspecies',
     'gene_id': '...',
     'location_id': '...',
     'sortBy': 'location',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gwas_by_location` {#gwas_by_location}

Returns GWAS Catalog associations (SNP, gene, allele, p-value, study) for traits linked to an anatomical location or its parts.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`graph`, `location_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **关键词：**gwas catalog, snp, uberon, anatomy, p-value

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `graph` | any | 是 | 可选值：`ebi_monarch_xspecies` |
| `location_id` | string | 是 | Anatomical Location |
| `sortBy` | any | 是 | 可选值：`location`, `trait`, `snp`, `effect_allele`, `gene`, `study_accession`, `p_value`, `pubmed_id`, `reported_trait`, `chromosome`, `position`, `risk_allele_frequency` |
| `sortDir` | any | 是 | 可选值：`asc`, `desc` |
| `pageNum` | integer | 是 | Page number (0-based) |
| `pageSize` | integer | 是 | Number of results per page |

??? note "上游说明（英文原文）"

    Search GWAS Catalog by anatomical location

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gwas_by_location',
    {'graph': 'ebi_monarch_xspecies',
     'location_id': '...',
     'sortBy': 'location',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gwas_by_pathway` {#gwas_by_pathway}

Returns GWAS Catalog SNP associations for genes in a Reactome pathway, with alleles, traits and p-values.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`graph`, `pathway_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **关键词：**gwas catalog, reactome, pathway genes, snp, p-value

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `graph` | any | 是 | 可选值：`ebi_monarch_xspecies` |
| `pathway_id` | string | 是 | Reactome Pathway |
| `sortBy` | any | 是 | 可选值：`pathway`, `gene`, `snp`, `effect_allele`, `trait`, `study_accession`, `p_value`, `pubmed_id`, `reported_trait`, `chromosome`, `position`, `risk_allele_frequency` |
| `sortDir` | any | 是 | 可选值：`asc`, `desc` |
| `pageNum` | integer | 是 | Page number (0-based) |
| `pageSize` | integer | 是 | Number of results per page |

??? note "上游说明（英文原文）"

    Search GWAS Catalog by Reactome pathway

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gwas_by_pathway',
    {'graph': 'ebi_monarch_xspecies',
     'pathway_id': '...',
     'sortBy': 'pathway',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gwas_studies_by_cell_type` {#gwas_studies_by_cell_type}

Returns GWAS Catalog studies annotated with traits linked to a cell type, with sample sizes, association counts and summary statistics availability.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`graph`, `cell_type_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **关键词：**gwas study, cell type, summary statistics, gcst, sample size

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `graph` | any | 是 | 可选值：`ebi_monarch_xspecies` |
| `cell_type_id` | string | 是 | Cell Type |
| `sortBy` | any | 是 | 可选值：`cell_type`, `trait`, `study`, `study_accession`, `reported_trait`, `first_author`, `journal`, `publication_date`, `pubmed_id`, `initial_sample_size`, `replication_sample_size`, `association_count` |
| `sortDir` | any | 是 | 可选值：`asc`, `desc` |
| `pageNum` | integer | 是 | Page number (0-based) |
| `pageSize` | integer | 是 | Number of results per page |

??? note "上游说明（英文原文）"

    Search GWAS Catalog studies by cell type

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gwas_studies_by_cell_type',
    {'graph': 'ebi_monarch_xspecies',
     'cell_type_id': '...',
     'sortBy': 'cell_type',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gwas_studies_by_disease_location` {#gwas_studies_by_disease_location}

Returns GWAS Catalog studies for diseases located in an anatomical location, with study metadata.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`graph`, `location_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **关键词：**gwas study, anatomy, disease location, summary statistics, gcst

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `graph` | any | 是 | 可选值：`ebi_monarch_xspecies` |
| `location_id` | string | 是 | Anatomical Location |
| `sortBy` | any | 是 | 可选值：`location`, `disease`, `study`, `study_accession`, `reported_trait`, `first_author`, `journal`, `publication_date`, `pubmed_id`, `initial_sample_size`, `replication_sample_size`, `association_count` |
| `sortDir` | any | 是 | 可选值：`asc`, `desc` |
| `pageNum` | integer | 是 | Page number (0-based) |
| `pageSize` | integer | 是 | Number of results per page |

??? note "上游说明（英文原文）"

    Search GWAS Catalog studies by disease location

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gwas_studies_by_disease_location',
    {'graph': 'ebi_monarch_xspecies',
     'location_id': '...',
     'sortBy': 'location',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gwas_studies_by_drugs_indicated_for_disease` {#gwas_studies_by_drugs_indicated_for_disease}

Returns GWAS Catalog studies whose traits involve drugs indicated for a disease, with study metadata.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`graph`, `disease_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **关键词：**gwas study, drug response, disease, gcst, summary statistics

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `graph` | any | 是 | 可选值：`ebi_monarch_xspecies` |
| `disease_id` | string | 是 | Disease |
| `sortBy` | any | 是 | 可选值：`disease`, `trait`, `drug`, `study`, `study_accession`, `reported_trait`, `first_author`, `journal`, `publication_date`, `pubmed_id`, `initial_sample_size`, `replication_sample_size` |
| `sortDir` | any | 是 | 可选值：`asc`, `desc` |
| `pageNum` | integer | 是 | Page number (0-based) |
| `pageSize` | integer | 是 | Number of results per page |

??? note "上游说明（英文原文）"

    Search GWAS Catalog studies by drugs or drug classes indicated for a disease

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gwas_studies_by_drugs_indicated_for_disease',
    {'graph': 'ebi_monarch_xspecies',
     'disease_id': '...',
     'sortBy': 'disease',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gwas_trait_to_mouse_models_via_embeddings` {#gwas_trait_to_mouse_models_via_embeddings}

Finds mouse gene models with phenotypes semantically similar to a GWAS trait whose human orthologs have GWAS associations with that trait.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`graph`, `trait_id`, `n_embedding_candidates`, `n_top_mp_hits`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **关键词：**gwas, mouse model, impc, orthologs, embedding similarity, mp phenotype

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `graph` | any | 是 | 可选值：`ebi_monarch_xspecies` |
| `trait_id` | string | 是 | GWAS Trait (EFO term) |
| `n_embedding_candidates` | string | 是 | Number of closest embedding matches to consider (across all ontologies) |
| `n_top_mp_hits` | string | 是 | Number of top MP phenotype matches to expand into mouse models |
| `sortBy` | any | 是 | 可选值：`trait`, `mp_phenotype`, `embedding_distance`, `mouse_gene`, `human_gene`, `snp`, `effect_allele`, `study_accession`, `p_value`, `or_or_beta` |
| `sortDir` | any | 是 | 可选值：`asc`, `desc` |
| `pageNum` | integer | 是 | Page number (0-based) |
| `pageSize` | integer | 是 | Number of results per page |

??? note "上游说明（英文原文）"

    GWAS trait to mouse models via embedding similarity

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gwas_trait_to_mouse_models_via_embeddings',
    {'graph': 'ebi_monarch_xspecies',
     'trait_id': '...',
     'n_embedding_candidates': '...',
     'n_top_mp_hits': '...',
     'sortBy': 'trait',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `mouse_gene_to_opentargets` {#mouse_gene_to_opentargets}

Returns Open Targets disease associations, with scores, for the human ortholog of a mouse gene.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`graph`, `mouse_gene_id`, `min_score`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **关键词：**open targets, mouse gene, human ortholog, target-disease, impc

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `graph` | any | 是 | 可选值：`ebi_monarch_xspecies`, `impc_x_gwas` |
| `mouse_gene_id` | string | 是 | Mouse Gene |
| `min_score` | string | 是 | Minimum Score |
| `sortBy` | any | 是 | 可选值：`mouse_gene`, `human_gene`, `disease`, `score` |
| `sortDir` | any | 是 | 可选值：`asc`, `desc` |
| `pageNum` | integer | 是 | Page number (0-based) |
| `pageSize` | integer | 是 | Number of results per page |

??? note "上游说明（英文原文）"

    Search OpenTargets disease associations by mouse gene

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'mouse_gene_to_opentargets',
    {'graph': 'ebi_monarch_xspecies',
     'mouse_gene_id': '...',
     'min_score': '...',
     'sortBy': 'mouse_gene',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```
