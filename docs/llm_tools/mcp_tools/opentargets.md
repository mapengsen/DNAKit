# Open Targets MCP (`opentargets`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **30** 个 DNA/基因相关工具（服务器共列出 68 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [nickzren/opentargets-mcp @ `ecb92f8`](https://github.com/nickzren/opentargets-mcp/tree/ecb92f8effdd6f8ea8f6a8c8de2fedb23334c126) |
| 许可证 | MIT |
| 启动命令 | `python -m opentargets_mcp.server --transport stdio` |
| 传输方式 | stdio |
| 安装 | Install nickzren/opentargets-mcp separately in Python &gt;=3.10; set executable_path to its interpreter. Uses the verified Python module entry point with explicit stdio. |
| 副作用 | Queries Open Targets evidence; results are research evidence and do not establish treatment suitability. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`get_credible_set_by_id`](#get_credible_set_by_id) | `association` | Fetches details of an Open Targets study-locus credible set by ID, including lead variant and fine-mapped variants. |
| [`get_credible_set_colocalisation`](#get_credible_set_colocalisation) | `association` | Retrieves colocalisation results (e.g. with eQTL or other GWAS) for a study-locus credible set. |
| [`get_credible_sets`](#get_credible_sets) | `association` | Queries Open Targets credible sets filtered by study-locus IDs, study IDs, variant IDs, study types or genomic regions. |
| [`get_disease_associated_targets`](#get_disease_associated_targets) | `association` | Lists targets (genes) associated with a disease (EFO ID) with overall and per-datatype association scores. |
| [`get_drug_pharmacogenomics`](#get_drug_pharmacogenomics) | `variant` | Retrieves pharmacogenomics evidence linking genetic variants and genotypes to response for a ChEMBL drug. |
| [`get_similar_targets`](#get_similar_targets) | `association` | Identifies targets (genes) whose disease-association profiles are similar to a seed target. |
| [`get_studies_by_disease`](#get_studies_by_disease) | `association` | Lists GWAS and molecular QTL studies linked to one or more diseases, optionally including indirect ontology children. |
| [`get_study_credible_sets`](#get_study_credible_sets) | `association` | Retrieves fine-mapped credible sets (loci) for a GWAS or QTL study. |
| [`get_study_info`](#get_study_info) | `association` | Retrieves metadata, trait, cohort and sample details for a GWAS study. |
| [`get_target_alternative_genes`](#get_target_alternative_genes) | `genome_annotation` | Lists alternate gene symbols, synonyms and database cross-references for an Ensembl gene target. |
| [`get_target_associated_diseases`](#get_target_associated_diseases) | `association` | Lists diseases linked to an Ensembl gene target with association scores. |
| [`get_target_depmap_essentiality`](#get_target_depmap_essentiality) | `genome_annotation` | Fetches DepMap CRISPR knockout gene essentiality scores across cancer cell lines for a gene. |
| [`get_target_disease_biomarkers`](#get_target_disease_biomarkers) | `association` | Inspects evidence for biomarker annotations (e.g. variants, expression) linking a target gene to a disease. |
| [`get_target_disease_evidence`](#get_target_disease_evidence) | `association` | Retrieves evidence strings linking an Ensembl target to an EFO disease, filterable by data source. |
| [`get_target_expression`](#get_target_expression) | `expression` | Returns baseline RNA and protein expression profiles for a target gene across tissues. |
| [`get_target_genetic_constraint`](#get_target_genetic_constraint) | `variant` | Fetches gnomAD genetic constraint metrics (e.g. pLI, LOEUF, missense Z) for a target gene. |
| [`get_target_hallmarks`](#get_target_hallmarks) | `genome_annotation` | Returns COSMIC cancer hallmark annotations and attributes associated with a target gene. |
| [`get_target_homologues`](#get_target_homologues) | `genome_annotation` | Retrieves cross-species orthologues and paralogues for a target gene. |
| [`get_target_info`](#get_target_info) | `genome_annotation` | Retrieves core identity details for a target gene: symbol, name, biotype and genomic location. |
| [`get_target_mouse_phenotypes`](#get_target_mouse_phenotypes) | `association` | Retrieves mouse model knockout phenotypes (MGI/IMPC) associated with a target gene. |
| [`get_target_pathways_and_go_terms`](#get_target_pathways_and_go_terms) | `genome_annotation` | Lists Reactome pathway memberships and Gene Ontology annotations for a target gene. |
| [`get_targets_batch`](#get_targets_batch) | `genome_annotation` | Retrieves core records for multiple targets by Ensembl gene IDs in one request. |
| [`get_variant_credible_sets`](#get_variant_credible_sets) | `association` | Lists fine-mapped credible sets that include a specific variant, filterable by study type. |
| [`get_variant_evidences`](#get_variant_evidences) | `association` | Retrieves evidence strings linking a variant to target genes or diseases, filterable by data source. |
| [`get_variant_info`](#get_variant_info) | `variant` | Retrieves core metadata, allele frequencies and functional consequence annotations for a variant. |
| [`get_variant_intervals`](#get_variant_intervals) | `regulation_epigenome` | Retrieves enhancer-to-gene (E2G) regulatory interval predictions overlapping a variant. |
| [`get_variant_pharmacogenomics`](#get_variant_pharmacogenomics) | `variant` | Retrieves pharmacogenomic annotations (drug response) for a variant. |
| [`get_variant_protein_coordinates`](#get_variant_protein_coordinates) | `variant` | Retrieves protein-level consequences and coordinates for a variant (amino acid change, position). |
| [`graphql_batch_query`](#graphql_batch_query) | `association` | Executes one Open Targets GraphQL query against many variable sets concurrently, keyed by a field. |
| [`graphql_query`](#graphql_query) | `association` | Executes a raw GraphQL query with variables against the Open Targets Platform API. |

## `get_credible_set_by_id` {#get_credible_set_by_id}

Fetches details of an Open Targets study-locus credible set by ID, including lead variant and fine-mapped variants.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`study_locus_id`
- **关键词：**fine-mapping, credible set, gwas, study locus, l2g
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `study_locus_id` | string | 是 |  |

??? note "上游说明（英文原文）"

    Fetch detailed information for a specific study locus credible set.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_credible_set_by_id',
    {'study_locus_id': '...'},
)
```

## `get_credible_set_colocalisation` {#get_credible_set_colocalisation}

Retrieves colocalisation results (e.g. with eQTL or other GWAS) for a study-locus credible set.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`study_locus_id`
- **关键词：**colocalisation, coloc, eqtl, credible set, gwas
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `study_locus_id` | string | 是 |  |
| `page_index` | integer | 否 | 默认 `0` |
| `page_size` | integer | 否 | 默认 `10` |

??? note "上游说明（英文原文）"

    Retrieve colocalisation results for a credible set.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_credible_set_colocalisation',
    {'study_locus_id': '...'},
)
```

## `get_credible_sets` {#get_credible_sets}

Queries Open Targets credible sets filtered by study-locus IDs, study IDs, variant IDs, study types or genomic regions.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**fine-mapping, credible set, gwas, qtl, genomic region
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `study_locus_ids` | array \| null | 否 | 默认 `None` |
| `study_ids` | array \| null | 否 | 默认 `None` |
| `variant_ids` | array \| null | 否 | 默认 `None` |
| `study_types` | array \| null | 否 | 默认 `None` |
| `regions` | array \| null | 否 | 默认 `None` |
| `page_index` | integer | 否 | 默认 `0` |
| `page_size` | integer | 否 | 默认 `10` |
| `fields` | array \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    Query credible sets with flexible filtering options.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_credible_sets',
    {},
)
```

## `get_disease_associated_targets` {#get_disease_associated_targets}

Lists targets (genes) associated with a disease (EFO ID) with overall and per-datatype association scores.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`efo_id`
- **关键词：**target-disease, efo, association score, genes, genetic evidence
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `efo_id` | string | 是 |  |
| `fields` | array \| null | 否 | 默认 `None` |
| `page_index` | integer | 否 | 默认 `0` |
| `page_size` | integer | 否 | 默认 `10` |

??? note "上游说明（英文原文）"

    List targets associated with a disease, including evidence scores.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_disease_associated_targets',
    {'efo_id': '...'},
)
```

## `get_drug_pharmacogenomics` {#get_drug_pharmacogenomics}

Retrieves pharmacogenomics evidence linking genetic variants and genotypes to response for a ChEMBL drug.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`chembl_id`
- **关键词：**pharmacogenomics, pharmgkb, genotype, drug response, chembl
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chembl_id` | string | 是 |  |
| `page_index` | integer | 否 | 默认 `0` |
| `page_size` | integer | 否 | 默认 `10` |

??? note "上游说明（英文原文）"

    Retrieve pharmacogenomics data linking genetic variants to drug response.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_drug_pharmacogenomics',
    {'chembl_id': '...'},
)
```

## `get_similar_targets` {#get_similar_targets}

Identifies targets (genes) whose disease-association profiles are similar to a seed target.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`entity_id`
- **关键词：**similar genes, association profile, target, embeddings
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `entity_id` | string | 是 |  |
| `threshold` | number \| null | 否 | 默认 `0.5` |
| `size` | integer | 否 | 默认 `10` |

??? note "上游说明（英文原文）"

    Identify targets with similar association profiles to the seed target.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_similar_targets',
    {'entity_id': '...'},
)
```

## `get_studies_by_disease` {#get_studies_by_disease}

Lists GWAS and molecular QTL studies linked to one or more diseases, optionally including indirect ontology children.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`disease_ids`
- **关键词：**gwas, study, efo, qtl, gwas catalog
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `disease_ids` | array of string | 是 |  |
| `enable_indirect` | boolean | 否 | 默认 `False` |
| `study_id` | string \| null | 否 | 默认 `None` |
| `page_index` | integer | 否 | 默认 `0` |
| `page_size` | integer | 否 | 默认 `10` |
| `fields` | array \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    List studies linked to one or more diseases.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_studies_by_disease',
    {'disease_ids': []},
)
```

## `get_study_credible_sets` {#get_study_credible_sets}

Retrieves fine-mapped credible sets (loci) for a GWAS or QTL study.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`study_id`
- **关键词：**fine-mapping, credible set, gwas, loci
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `study_id` | string | 是 |  |
| `page_index` | integer | 否 | 默认 `0` |
| `page_size` | integer | 否 | 默认 `10` |
| `fields` | array \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    Retrieve fine-mapped credible sets for a study.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_study_credible_sets',
    {'study_id': '...'},
)
```

## `get_study_info` {#get_study_info}

Retrieves metadata, trait, cohort and sample details for a GWAS study.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`study_id`
- **关键词：**gwas, cohort, sample size, trait, ancestry
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `study_id` | string | 是 |  |
| `fields` | array \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    Retrieve metadata and cohort details for a GWAS study.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_study_info',
    {'study_id': '...'},
)
```

## `get_target_alternative_genes` {#get_target_alternative_genes}

Lists alternate gene symbols, synonyms and database cross-references for an Ensembl gene target.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`ensembl_id`
- **关键词：**gene symbol, synonyms, xrefs, ensembl id
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ensembl_id` | string | 是 |  |

??? note "上游说明（英文原文）"

    List alternate gene symbols and database cross-references for a target.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_target_alternative_genes',
    {'ensembl_id': '...'},
)
```

## `get_target_associated_diseases` {#get_target_associated_diseases}

Lists diseases linked to an Ensembl gene target with association scores.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`ensembl_id`
- **关键词：**target-disease, association score, efo, genetic evidence
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ensembl_id` | string | 是 |  |
| `fields` | array \| null | 否 | 默认 `None` |
| `page_index` | integer | 否 | 默认 `0` |
| `page_size` | integer | 否 | 默认 `10` |

??? note "上游说明（英文原文）"

    List diseases linked to a target with association scores.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_target_associated_diseases',
    {'ensembl_id': '...'},
)
```

## `get_target_depmap_essentiality` {#get_target_depmap_essentiality}

Fetches DepMap CRISPR knockout gene essentiality scores across cancer cell lines for a gene.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`ensembl_id`
- **关键词：**depmap, crispr screen, gene essentiality, cell lines, dependency
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ensembl_id` | string | 是 |  |

??? note "上游说明（英文原文）"

    Fetch DepMap CRISPR essentiality scores across cell lines.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_target_depmap_essentiality',
    {'ensembl_id': '...'},
)
```

## `get_target_disease_biomarkers` {#get_target_disease_biomarkers}

Inspects evidence for biomarker annotations (e.g. variants, expression) linking a target gene to a disease.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`ensembl_id`, `efo_id`
- **关键词：**biomarker, target-disease, evidence, efo
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ensembl_id` | string | 是 |  |
| `efo_id` | string | 是 |  |
| `size` | integer | 否 | 默认 `10` |
| `cursor` | string \| null | 否 | 默认 `None` |
| `fields` | array \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    Inspect evidence for biomarker annotations linking a target and disease.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_target_disease_biomarkers',
    {'ensembl_id': '...', 'efo_id': '...'},
)
```

## `get_target_disease_evidence` {#get_target_disease_evidence}

Retrieves evidence strings linking an Ensembl target to an EFO disease, filterable by data source.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`ensembl_id`, `efo_id`
- **关键词：**evidence, target-disease, genetic association, datasource
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ensembl_id` | string | 是 |  |
| `efo_id` | string | 是 |  |
| `datasource_ids` | array \| null | 否 | 默认 `None` |
| `size` | integer | 否 | 默认 `10` |
| `cursor` | string \| null | 否 | 默认 `None` |
| `fields` | array \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    Retrieve evidence strings linking a target to a disease.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_target_disease_evidence',
    {'ensembl_id': '...', 'efo_id': '...'},
)
```

## `get_target_expression` {#get_target_expression}

Returns baseline RNA and protein expression profiles for a target gene across tissues.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`ensembl_id`
- **关键词：**tissue expression, rna, gtex, hpa, baseline
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ensembl_id` | string | 是 |  |
| `fields` | array \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    Return RNA and protein expression profiles for a target across tissues.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_target_expression',
    {'ensembl_id': '...'},
)
```

## `get_target_genetic_constraint` {#get_target_genetic_constraint}

Fetches gnomAD genetic constraint metrics (e.g. pLI, LOEUF, missense Z) for a target gene.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`ensembl_id`
- **关键词：**gnomad, pli, loeuf, loss-of-function, constraint
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ensembl_id` | string | 是 |  |

??? note "上游说明（英文原文）"

    Fetch genetic constraint metrics (gnomAD) for a target.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_target_genetic_constraint',
    {'ensembl_id': '...'},
)
```

## `get_target_hallmarks` {#get_target_hallmarks}

Returns COSMIC cancer hallmark annotations and attributes associated with a target gene.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`ensembl_id`
- **关键词：**cancer hallmarks, cosmic, gene function, oncogene
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ensembl_id` | string | 是 |  |

??? note "上游说明（英文原文）"

    Return cancer hallmark annotations associated with a target.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_target_hallmarks',
    {'ensembl_id': '...'},
)
```

## `get_target_homologues` {#get_target_homologues}

Retrieves cross-species orthologues and paralogues for a target gene.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`ensembl_id`
- **关键词：**orthologs, paralogs, homology, species
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ensembl_id` | string | 是 |  |

??? note "上游说明（英文原文）"

    Retrieve cross-species homologues for a target gene.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_target_homologues',
    {'ensembl_id': '...'},
)
```

## `get_target_info` {#get_target_info}

Retrieves core identity details for a target gene: symbol, name, biotype and genomic location.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`ensembl_id`
- **关键词：**gene symbol, biotype, genomic location, ensembl id
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ensembl_id` | string | 是 |  |
| `fields` | array \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    Retrieve core identity details for a target gene.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_target_info',
    {'ensembl_id': '...'},
)
```

## `get_target_mouse_phenotypes` {#get_target_mouse_phenotypes}

Retrieves mouse model knockout phenotypes (MGI/IMPC) associated with a target gene.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`ensembl_id`
- **关键词：**mouse knockout, mgi, impc, phenotype, model organism
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ensembl_id` | string | 是 |  |
| `page_index` | integer | 否 | 默认 `0` |
| `page_size` | integer | 否 | 默认 `10` |

??? note "上游说明（英文原文）"

    Retrieve mouse knockout phenotypes associated with a target.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_target_mouse_phenotypes',
    {'ensembl_id': '...'},
)
```

## `get_target_pathways_and_go_terms` {#get_target_pathways_and_go_terms}

Lists Reactome pathway memberships and Gene Ontology annotations for a target gene.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`ensembl_id`
- **关键词：**reactome, go annotations, pathways, gene function
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ensembl_id` | string | 是 |  |
| `page_index` | integer | 否 | 默认 `0` |
| `page_size` | integer | 否 | 默认 `10` |

??? note "上游说明（英文原文）"

    List pathway memberships and Gene Ontology annotations for a target.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_target_pathways_and_go_terms',
    {'ensembl_id': '...'},
)
```

## `get_targets_batch` {#get_targets_batch}

Retrieves core records for multiple targets by Ensembl gene IDs in one request.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`ensembl_ids`
- **关键词：**ensembl id, genes, batch, target info
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ensembl_ids` | array of string | 是 |  |

??? note "上游说明（英文原文）"

    Retrieve multiple targets in a single request.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_targets_batch',
    {'ensembl_ids': []},
)
```

## `get_variant_credible_sets` {#get_variant_credible_sets}

Lists fine-mapped credible sets that include a specific variant, filterable by study type.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`variant_id`
- **关键词：**fine-mapping, credible set, gwas, qtl, variant
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_id` | string | 是 |  |
| `study_types` | array \| null | 否 | 默认 `None` |
| `page_index` | integer | 否 | 默认 `0` |
| `page_size` | integer | 否 | 默认 `10` |

??? note "上游说明（英文原文）"

    List credible sets that include a specific variant.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_variant_credible_sets',
    {'variant_id': '...'},
)
```

## `get_variant_evidences` {#get_variant_evidences}

Retrieves evidence strings linking a variant to target genes or diseases, filterable by data source.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`variant_id`
- **关键词：**variant-to-gene, evidence, variant-disease, l2g
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_id` | string | 是 |  |
| `datasource_ids` | array \| null | 否 | 默认 `None` |
| `size` | integer | 否 | 默认 `10` |
| `cursor` | string \| null | 否 | 默认 `None` |
| `fields` | array \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    Retrieve evidence strings linking a variant to targets or diseases.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_variant_evidences',
    {'variant_id': '...'},
)
```

## `get_variant_info` {#get_variant_info}

Retrieves core metadata, allele frequencies and functional consequence annotations for a variant.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant_id`
- **关键词：**variant annotation, rsid, allele frequency, vep, consequence
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_id` | string | 是 |  |
| `fields` | array \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    Retrieve core metadata and functional annotations for a variant.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_variant_info',
    {'variant_id': '...'},
)
```

## `get_variant_intervals` {#get_variant_intervals}

Retrieves enhancer-to-gene (E2G) regulatory interval predictions overlapping a variant.

- **类别：**`regulation_epigenome`（Regulatory elements, TF binding, chromatin, methylation and 3D genome）
- **必填参数：**`variant_id`
- **关键词：**enhancer, e2g, regulatory element, variant-to-gene
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_id` | string | 是 |  |
| `page_index` | integer | 否 | 默认 `0` |
| `page_size` | integer | 否 | 默认 `10` |

??? note "上游说明（英文原文）"

    Retrieve enhancer-to-gene (E2G) predictions overlapping a variant.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_variant_intervals',
    {'variant_id': '...'},
)
```

## `get_variant_pharmacogenomics` {#get_variant_pharmacogenomics}

Retrieves pharmacogenomic annotations (drug response) for a variant.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant_id`
- **关键词：**pharmacogenomics, pharmgkb, drug response, genotype
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_id` | string | 是 |  |
| `page_index` | integer | 否 | 默认 `0` |
| `page_size` | integer | 否 | 默认 `10` |

??? note "上游说明（英文原文）"

    Retrieve pharmacogenomic annotations for a variant.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_variant_pharmacogenomics',
    {'variant_id': '...'},
)
```

## `get_variant_protein_coordinates` {#get_variant_protein_coordinates}

Retrieves protein-level consequences and coordinates for a variant (amino acid change, position).

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant_id`
- **关键词：**amino acid change, protein consequence, missense, variant effect
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_id` | string | 是 |  |
| `page_index` | integer | 否 | 默认 `0` |
| `page_size` | integer | 否 | 默认 `10` |

??? note "上游说明（英文原文）"

    Retrieve protein-level consequences for a variant.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_variant_protein_coordinates',
    {'variant_id': '...'},
)
```

## `graphql_batch_query` {#graphql_batch_query}

Executes one Open Targets GraphQL query against many variable sets concurrently, keyed by a field.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`query_string`, `variables_list`
- **关键词：**graphql, batch query, target-disease, genetic evidence
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query_string` | string | 是 |  |
| `variables_list` | array of object | 是 |  |
| `key_field` | string \| null | 否 | 默认 `None` |
| `operation_name` | string \| null | 否 | 默认 `None` |
| `max_concurrency` | integer | 否 | 默认 `4` |

??? note "上游说明（英文原文）"

    Execute one GraphQL query against many variable sets.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'graphql_batch_query',
    {'query_string': '...', 'variables_list': []},
)
```

## `graphql_query` {#graphql_query}

Executes a raw GraphQL query with variables against the Open Targets Platform API.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`query_string`
- **关键词：**graphql, raw query, target-disease, genetic evidence
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query_string` | string | 是 |  |
| `variables` | object \| null | 否 | 默认 `None` |
| `operation_name` | string \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    ADVANCED: Execute a raw GraphQL query against Open Targets.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'graphql_query',
    {'query_string': '...'},
)
```
