# ToolUniverse: expression (`tooluniverse`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **845** 个 DNA/基因相关工具（服务器共列出 2927 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

本页列出类别 `expression` 的 56 个工具，其他类别：[`sequence`](tooluniverse-sequence.md) (44), [`genome_annotation`](tooluniverse-genome_annotation.md) (253), [`variant`](tooluniverse-variant.md) (155), [`association`](tooluniverse-association.md) (138), [`regulation_epigenome`](tooluniverse-regulation_epigenome.md) (54), [`sequencing_data`](tooluniverse-sequencing_data.md) (91), [`prediction_model`](tooluniverse-prediction_model.md) (17), [`molecular_tools`](tooluniverse-molecular_tools.md) (23), [`breeding_population`](tooluniverse-breeding_population.md) (14)。

| 项目 | 内容 |
| --- | --- |
| 源码 | [mims-harvard/ToolUniverse @ `9959ebe`](https://github.com/mims-harvard/ToolUniverse/tree/9959ebee7638773314362f1746980a15f5b17e5b) |
| 许可证 | Apache-2.0 |
| 启动命令 | `tooluniverse` |
| 传输方式 | stdio |
| 安装 | Install the pinned ToolUniverse distribution and required provider dependencies separately. Use the installed tooluniverse console script; unlike the README uvx --refresh example, this bridge does not refresh or install anything. |
| 副作用 | Starts the compact ToolUniverse MCP catalog. Selected tools may contact services, load models, execute programs or write files according to their provider definitions. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`AllenBrain_get_expression_datasets`](#AllenBrain_get_expression_datasets) | `expression` | Lists Allen Mouse Brain Atlas in situ hybridization expression datasets for a given gene. |
| [`AllenBrain_get_structure_expression_values`](#AllenBrain_get_structure_expression_values) | `expression` | Returns per-brain-structure expression energy and density values for one Allen ISH experiment. |
| [`ARCHS4_get_gene_correlations`](#ARCHS4_get_gene_correlations) | `expression` | Returns the top co-expressed genes for a query gene from ARCHS4 RNA-seq correlations. |
| [`ARCHS4_get_gene_expression`](#ARCHS4_get_gene_expression) | `expression` | Returns tissue and cell-line expression distributions (log2 TPM) for a gene from ARCHS4. |
| [`BAR_get_rnaseq_expression`](#BAR_get_rnaseq_expression) | `expression` | Returns RNA-seq or single-cell cluster expression for a plant gene from BAR. |
| [`Bgee_get_gene_expression`](#Bgee_get_gene_expression) | `expression` | Retrieves Bgee expression calls for a gene across tissues and organs. |
| [`CancerPrognosis_get_gene_expression`](#CancerPrognosis_get_gene_expression) | `expression` | Fetches per-sample RNA-seq expression for a gene in a TCGA or cBioPortal study. |
| [`CellMarker_search_by_cell_type`](#CellMarker_search_by_cell_type) | `expression` | Retrieves curated marker genes for a cell type from CellMarker 2.0. |
| [`CellMarker_search_by_gene`](#CellMarker_search_by_gene) | `expression` | Finds cell types that express a given marker gene in CellMarker 2.0. |
| [`CellMarker_search_cancer_markers`](#CellMarker_search_cancer_markers) | `expression` | Searches CellMarker 2.0 for marker genes and cell types in cancer tissue contexts. |
| [`CELLxGENE_get_expression_data`](#CELLxGENE_get_expression_data) | `expression` | Queries gene expression data from CELLxGENE Census with cell and gene filters. |
| [`CELLxGENE_get_gene_metadata`](#CELLxGENE_get_gene_metadata) | `expression` | Queries gene metadata such as symbols, Ensembl IDs and summary statistics from CELLxGENE Census. |
| [`CELLxGENE_get_presence_matrix`](#CELLxGENE_get_presence_matrix) | `expression` | Returns the feature presence matrix of genes measured per CELLxGENE dataset. |
| [`Coexpression_modules`](#Coexpression_modules) | `expression` | Detects WGCNA-style co-expression modules from a gene expression matrix, returning module genes and eigengenes. |
| [`DeepSpotM_predict_gene_expression`](#DeepSpotM_predict_gene_expression) | `expression` | Predicts per-gene spatial expression for an H&amp;E tile locally with the DeepSpot-M model. |
| [`expression_anova_per_gene`](#expression_anova_per_gene) | `expression` | Runs per-gene ANOVA or log2 fold change on a gene x sample expression matrix. |
| [`ExpressionAtlas_get_baseline`](#ExpressionAtlas_get_baseline) | `expression` | Lists baseline expression experiments in EBI Expression Atlas for a species. |
| [`ExpressionAtlas_get_experiment`](#ExpressionAtlas_get_experiment) | `expression` | Gets metadata for an Expression Atlas experiment by accession. |
| [`ExpressionAtlas_search_differential`](#ExpressionAtlas_search_differential) | `expression` | Searches Expression Atlas differential expression experiments by condition and gene text. |
| [`ExpressionAtlas_search_experiments`](#ExpressionAtlas_search_experiments) | `expression` | Searches all Expression Atlas experiments by condition and gene text. |
| [`FlyBase_get_gene_expression`](#FlyBase_get_gene_expression) | `expression` | Returns Drosophila gene expression annotations by developmental stage and anatomical location from the Alliance. |
| [`GDC_get_gene_expression`](#GDC_get_gene_expression) | `expression` | Queries GDC/TCGA RNA-seq gene expression quantification file metadata by project and gene ID. |
| [`GeneNetwork_get_sample_data`](#GeneNetwork_get_sample_data) | `expression` | Gets expression or phenotype values with SEs for a gene/trait across strains of a GeneNetwork dataset. |
| [`GSVA_score`](#GSVA_score) | `expression` | Computes GSVA gene set by sample activity scores from a genes x samples expression matrix. |
| [`GTEx_get_expression_summary`](#GTEx_get_expression_summary) | `expression` | Summarizes a gene's tissue-specific expression (median TPM) across GTEx tissues. |
| [`GTEx_get_gene_expression`](#GTEx_get_gene_expression) | `expression` | Gets per-sample normalized TPM expression for a gene in GTEx tissues, optionally by sex/age subset. |
| [`GTEx_get_median_gene_expression`](#GTEx_get_median_gene_expression) | `expression` | Gets median TPM expression for genes across GTEx tissues. |
| [`GTEx_get_median_transcript_expression`](#GTEx_get_median_transcript_expression) | `expression` | Gets per-transcript median TPM across GTEx tissues for isoform usage analysis. |
| [`GTEx_get_single_nucleus_expression`](#GTEx_get_single_nucleus_expression) | `expression` | Gets GTEx snRNA-seq gene expression summarized by cell type (mean/median, cell counts). |
| [`GTEx_get_top_expressed_genes`](#GTEx_get_top_expressed_genes) | `expression` | Gets genes ranked by median expression in a GTEx tissue, optionally excluding mitochondrial genes. |
| [`GxA_get_experiment_expression`](#GxA_get_experiment_expression) | `expression` | Gets gene expression levels across tissues/conditions from an EBI Expression Atlas experiment. |
| [`GxA_get_experiment_info`](#GxA_get_experiment_info) | `expression` | Gets Expression Atlas experiment metadata: species, gene and condition counts, design. |
| [`GxA_list_experiments`](#GxA_list_experiments) | `expression` | Lists EBI Expression Atlas experiments filtered by species and baseline/differential type. |
| [`HPA_get_cancer_prognostics_by_gene`](#HPA_get_cancer_prognostics_by_gene) | `expression` | Gets whether a gene's expression correlates with survival across cancer types in HPA. |
| [`HPA_get_comparative_expression_by_gene_and_cellline`](#HPA_get_comparative_expression_by_gene_and_cellline) | `expression` | Compares a gene's expression between a cell line and healthy tissues in HPA. |
| [`HPA_get_contextual_biological_process_analysis`](#HPA_get_contextual_biological_process_analysis) | `expression` | Analyzes a gene's biological processes in a tissue or cell line context by combining annotation with expression. |
| [`HPA_get_disease_expression_by_gene_tissue_disease`](#HPA_get_disease_expression_by_gene_tissue_disease) | `expression` | Compares a gene's expression in a disease state versus healthy tissue in HPA. |
| [`HPA_get_rna_expression_by_source`](#HPA_get_rna_expression_by_source) | `expression` | Gets a gene's RNA expression in a tissue, blood, brain or single-cell source from HPA. |
| [`HPA_get_rna_expression_in_specific_tissues`](#HPA_get_rna_expression_in_specific_tissues) | `expression` | Gets a gene's RNA expression (nTPM) in specified HPA consensus tissues. |
| [`L1000FWD_sig_search`](#L1000FWD_sig_search) | `expression` | Searches the LINCS L1000 connectivity map with up and down regulated gene sets to find similar or opposite signatures. |
| [`LINCS_search_signatures`](#LINCS_search_signatures) | `expression` | Searches LINCS for drug perturbation gene expression signatures by compound and cell line. |
| [`OpenTargets_get_target_expression_by_ensemblID`](#OpenTargets_get_target_expression_by_ensemblID) | `expression` | Retrieves baseline expression statistics for a gene across tissues and cell types from GTEx, Tabula Sapiens and others. |
| [`PanglaoDB_cell_types_for_gene`](#PanglaoDB_cell_types_for_gene) | `expression` | Lists cell types for which a gene is a curated PanglaoDB marker, with sensitivity/specificity. |
| [`PanglaoDB_markers_for_cell_type`](#PanglaoDB_markers_for_cell_type) | `expression` | Returns curated PanglaoDB marker genes for a cell type with sensitivity/specificity for human and mouse. |
| [`Pharos_get_target_expression`](#Pharos_get_target_expression) | `expression` | Returns GTEx per-tissue TPM expression for a target gene from Pharos. |
| [`ReactomeAnalysis_expression_analysis`](#ReactomeAnalysis_expression_analysis) | `expression` | Runs Reactome expression analysis mapping fold-change values for identifiers onto pathways. |
| [`RNAseq_edger_limma_de`](#RNAseq_edger_limma_de) | `expression` | Runs edgeR QL-F or limma-voom differential expression on a count matrix and metadata. |
| [`run_deseq2_analysis`](#run_deseq2_analysis) | `expression` | Runs R DESeq2 differential expression on a count matrix and metadata with optional enrichment. |
| [`SCXA_get_cluster_marker_genes`](#SCXA_get_cluster_marker_genes) | `expression` | Gets marker genes per cell cluster or cell type for an EBI Single Cell Expression Atlas experiment, with ranking statistics. |
| [`SCXA_search_gene`](#SCXA_search_gene) | `expression` | Searches the EBI Single Cell Expression Atlas for scRNA-seq experiments in which a given gene symbol or Ensembl ID is expressed. |
| [`ssGSEA_score`](#ssGSEA_score) | `expression` | Computes single-sample GSEA scores for each gene set in each sample from a genes-by-samples expression matrix. |
| [`SYNERGxDB_get_biomarker_association`](#SYNERGxDB_get_biomarker_association) | `expression` | Correlates per-cell-line RNA-seq expression of a gene with drug-combination synergy scores in SYNERGxDB. |
| [`TIMER2_gene_correlation`](#TIMER2_gene_correlation) | `expression` | Computes Spearman correlation between expression of two genes across TCGA cancer samples using cBioPortal data. |
| [`TIMER2_immune_estimation`](#TIMER2_immune_estimation) | `expression` | Estimates immune cell infiltration proxies in TCGA samples from marker gene expression, optionally correlating with a gene. |
| [`WormBase_get_expression`](#WormBase_get_expression) | `expression` | Gets expression data for a C. elegans gene from WormBase: tissues, developmental stages, localization and expression clusters. |
| [`ZFIN_get_gene_expression`](#ZFIN_get_gene_expression) | `expression` | Gets zebrafish gene expression annotations (developmental stage and anatomical location) from the Alliance API. |

## `AllenBrain_get_expression_datasets` {#AllenBrain_get_expression_datasets}

Lists Allen Mouse Brain Atlas in situ hybridization expression datasets for a given gene.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`gene_acronym`
- **关键词：**in situ hybridization, ish, mouse brain, gene expression
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_acronym` | string | 是 | Gene symbol/acronym. Examples: 'Gad1', 'Pvalb', 'Sst', 'Bdnf', 'Th', 'Slc17a7'. |
| `product_id` | integer | 否 | Product ID for the dataset type. 1=Mouse Brain ISH, 2=Mouse Brain microarray. Default: 1. 默认 `1` |
| `num_rows` | integer | 否 | Max results to return. Default: 50. 默认 `50` |

??? note "上游说明（英文原文）"

    Get gene expression experiment datasets from the Allen Mouse Brain Atlas for a given gene. Returns section data sets including plane of section, section thickness, and linked gene info. These datasets contain ISH (in situ hybridization) gene expression images. Example: 'Gad1' returns 22 ISH experiment datasets.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'AllenBrain_get_expression_datasets',
     'arguments': {'gene_acronym': '...'}},
)
```

## `AllenBrain_get_structure_expression_values` {#AllenBrain_get_structure_expression_values}

Returns per-brain-structure expression energy and density values for one Allen ISH experiment.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`section_data_set_id`
- **关键词：**expression energy, brain region, ish quantification, spatial expression
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `section_data_set_id` | integer | 是 | Allen Brain Atlas SectionDataSet (ISH experiment) ID. Get one from AllenBrain_get_expression_datasets. Example: 75457539 (a Gad1 mouse-brain ISH dataset). |
| `include_structure` | boolean | 否 | Whether to include the linked brain structure record (acronym, name, hierarchy) for each row. Default: true. 默认 `True` |
| `num_rows` | integer | 否 | Max per-structure rows to return. There are ~2400 structures per dataset. Default: 50. 默认 `50` |

??? note "上游说明（英文原文）"

    Get quantified per-structure gene-expression values (StructureUnionize) for one Allen Mouse Brain ISH experiment (SectionDataSet). Returns the actual numeric spatial expression aggregated over each annotated brain structure: expression_energy (intensity x density), expression_density (fraction of expressing pixels), sum_expressing_pixels, plus the structure record (acronym, name). This is the real spatial expression signal, not just experiment metadata. Get a section_data_set_id from AllenBrain_get_expression_datasets. Example: section_data_set_id=75457539 (a Gad1 mouse-brain ISH dataset) returns 2406 per-structure rows.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'AllenBrain_get_structure_expression_values',
     'arguments': {'section_data_set_id': 0}},
)
```

## `ARCHS4_get_gene_correlations` {#ARCHS4_get_gene_correlations}

Returns the top co-expressed genes for a query gene from ARCHS4 RNA-seq correlations.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`gene`
- **关键词：**co-expression, pearson correlation, rna-seq, archs4
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 是 | Gene symbol (e.g., TP53, BRCA1, EGFR, MYC) |
| `count` | integer \| null | 否 | Number of top correlated genes to return (default: 20, max: 200) 默认 `20` |

??? note "上游说明（英文原文）"

    Get genes most co-expressed with a query gene from ARCHS4, based on Pearson correlation computed across 300K+ uniformly processed RNA-seq samples. Returns the top N most correlated genes with their Pearson correlation coefficients. Useful for identifying functionally related genes, predicting gene function, and building co-expression networks.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ARCHS4_get_gene_correlations',
     'arguments': {'gene': '...'}},
)
```

## `ARCHS4_get_gene_expression` {#ARCHS4_get_gene_expression}

Returns tissue and cell-line expression distributions (log2 TPM) for a gene from ARCHS4.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`gene`
- **关键词：**tissue expression, tpm, cell lines, rna-seq
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 是 | Gene symbol (e.g., TP53, BRCA1, EGFR, MYC) |
| `species` | string \| null | 否 | Species: human or mouse (default: human) 可选值：`human`, `mouse`, `None` 默认 `human` |
| `type` | string \| null | 否 | Expression type: tissue (organs/tissues) or cellline (cell lines). Default: tissue 可选值：`tissue`, `cellline`, `None` 默认 `tissue` |

??? note "上游说明（英文原文）"

    Get gene expression levels across human or mouse tissues and cell lines from ARCHS4 (All RNA-seq and ChIP-seq Sample and Signature Search). Returns box-plot statistics (min, Q1, median, Q3, max) in log2(TPM+1) units for each tissue/cell line, computed from 300K+ uniformly processed RNA-seq samples. Choose type='tissue' for organ/tissue expression or type='cellline' for cell line expression. Results are sorted by median expression (highest first).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ARCHS4_get_gene_expression',
     'arguments': {'gene': '...'}},
)
```

## `BAR_get_rnaseq_expression` {#BAR_get_rnaseq_expression}

Returns RNA-seq or single-cell cluster expression for a plant gene from BAR.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`gene_id`
- **关键词：**plant expression, rna-seq, single-cell, arabidopsis
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_id` | string | 是 | Gene locus identifier, e.g. 'At1g01010' (Arabidopsis). |
| `species` | string \| null | 否 | Species, e.g. 'arabidopsis' (default). |
| `database` | string \| null | 否 | RNA-seq dataset, e.g. 'single_cell' (default). |

??? note "上游说明（英文原文）"

    Get RNA-seq gene expression data from the BAR (Bio-Analytic Resource for Plant Biology). database='single_cell' (the default) returns per-cluster mean expression across replicates from single-cell atlases -- a data type no other ToolUniverse plant tool exposes. Example: gene_id='At1g01010' returns its mean expression value in each of dozens of single-cell clusters.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'BAR_get_rnaseq_expression',
     'arguments': {'gene_id': '...'}},
)
```

## `Bgee_get_gene_expression` {#Bgee_get_gene_expression}

Retrieves Bgee expression calls for a gene across tissues and organs.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`gene_id`, `species_id`
- **关键词：**tissue expression, expression calls, anatomy, bgee
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_id` | string | 是 | Ensembl gene ID. Examples: 'ENSG00000141510' (human TP53), 'ENSMUSG00000059552' (mouse Trp53), 'ENSDARG00000115148' (zebrafish tp53). |
| `species_id` | string | 是 | NCBI taxonomy ID for the species. Examples: '9606' (human), '10090' (mouse), '7955' (zebrafish), '7227' (fruit fly). |

??? note "上游说明（英文原文）"

    Get expression data for a gene across tissues and organs from the Bgee database. Returns where a gene is expressed with expression scores, confidence levels, and quality ratings. Supports RNA-Seq, Affymetrix, EST, and in situ hybridization data types. Requires both an Ensembl gene ID and NCBI taxonomy ID. Example: ENSG00000141510 (TP53) with species 9606 (human) returns expression across 145+ tissues including ventricular zone (score 95.11), ganglionic eminence, and monocytes.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Bgee_get_gene_expression',
     'arguments': {'gene_id': '...', 'species_id': '...'}},
)
```

## `CancerPrognosis_get_gene_expression` {#CancerPrognosis_get_gene_expression}

Fetches per-sample RNA-seq expression for a gene in a TCGA or cBioPortal study.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**无
- **关键词：**tcga, tumor expression, rna-seq, cbioportal
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`get_gene_expression` 默认 `get_gene_expression` |
| `cancer` | string | 否 | TCGA cancer type abbreviation (e.g., 'BRCA', 'LUAD', 'COADREAD') or cBioPortal study ID. Use cancer_type as alias. |
| `cancer_type` | string | 否 | Alias for cancer. TCGA cancer type (e.g., 'BRCA', 'LUAD', 'COADREAD') or cBioPortal study ID. |
| `study_id` | string | 否 | Alias for cancer. cBioPortal study ID (e.g., 'brca_tcga', 'luad_tcga'). |
| `gene` | string | 否 | Gene symbol (e.g., 'TP53', 'BRCA1', 'EGFR', 'CD8A') |
| `gene_symbol` | string | 否 | Alias for gene. Gene symbol (e.g., 'TP53', 'BRCA1', 'EGFR'). |
| `gene_name` | string | 否 | Alias for gene. Gene symbol or name (e.g., 'TP53', 'BRCA1'). |
| `max_samples` | integer \| null | 否 | Upper bound on sample records to return (default 500, max 2000). Actual count may be lower if some samples lack expression data for the requested gene. |

??? note "上游说明（英文原文）"

    Fetch gene expression values (RNA-seq) for a specific gene across cancer samples in a TCGA or cBioPortal study. Returns per-sample expression values with summary statistics (mean, median, min, max). Output can be combined with survival data for expression-based prognostic analysis. Supports 33 TCGA cancer types and any cBioPortal study with mRNA expression data.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CancerPrognosis_get_gene_expression', 'arguments': {}},
)
```

## `CellMarker_search_by_cell_type` {#CellMarker_search_by_cell_type}

Retrieves curated marker genes for a cell type from CellMarker 2.0.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`cell_name`
- **关键词：**marker genes, cell type, single-cell, cellmarker
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`search_by_cell_type` 默认 `search_by_cell_type` |
| `cell_name` | string | 是 | Cell type name to search for (e.g., 'T cell', 'Macrophage', 'B cell', 'Fibroblast') |
| `species` | string | 否 | Species filter. If omitted, searches both species. 可选值：`Human`, `Mouse` |
| `tissue_type` | string | 否 | Filter results by tissue type (e.g., 'Blood', 'Lung', 'Brain'). Case-insensitive partial match. |

??? note "上游说明（英文原文）"

    Search the CellMarker 2.0 database for marker genes of a specific cell type. Returns curated lists of marker genes that define a cell type, sourced from single-cell RNA-seq studies and experimental literature. Use this to find canonical markers for cell type annotation, or to build gene sets for cell scoring. Example: 'T cell' returns CD3D, CD3E, CD3G, etc. 'Macrophage' returns CD68, CD163, MARCO, etc.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CellMarker_search_by_cell_type',
     'arguments': {'cell_name': '...'}},
)
```

## `CellMarker_search_by_gene` {#CellMarker_search_by_gene}

Finds cell types that express a given marker gene in CellMarker 2.0.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`gene_symbol`
- **关键词：**marker gene, cell type annotation, scrna-seq
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`search_by_gene` 默认 `search_by_gene` |
| `gene_symbol` | string | 是 | Gene symbol to search for (e.g., 'CD3D', 'CD68', 'EPCAM', 'MS4A1') |
| `species` | string | 否 | Species filter: 'Human' or 'Mouse'. If omitted, searches both species. 可选值：`Human`, `Mouse` |
| `tissue_type` | string | 否 | Filter results by tissue type (e.g., 'Blood', 'Lung', 'Brain'). Case-insensitive partial match. |

??? note "上游说明（英文原文）"

    Search the CellMarker 2.0 database for cell types that express a given marker gene. Returns curated cell type annotations from single-cell RNA-seq studies and experimental literature. CellMarker 2.0 covers &gt;500 cell types across &gt;400 tissue types in human and mouse. Use this to identify which cell types a gene marks, enabling cell type annotation in scRNA-seq analysis. Example: CD3D marks T cells, CD68 marks macrophages, EPCAM marks epithelial cells.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CellMarker_search_by_gene',
     'arguments': {'gene_symbol': '...'}},
)
```

## `CellMarker_search_cancer_markers` {#CellMarker_search_cancer_markers}

Searches CellMarker 2.0 for marker genes and cell types in cancer tissue contexts.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**无
- **关键词：**cancer markers, tumor microenvironment, marker genes
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`search_cancer_markers` 默认 `search_cancer_markers` |
| `cancer_type` | string | 否 | Cancer tissue type (e.g., 'Breast', 'Lung', 'Brain', 'Liver'). Filters by tissue where cancer occurs. |
| `gene_symbol` | string | 否 | Marker gene to search in cancer context (e.g., 'CD274' for PD-L1, 'EPCAM') |
| `cell_type` | string | 否 | Cancer cell type to search for (e.g., 'Cancer stem cell', 'T cell', 'Macrophage') |

??? note "上游说明（英文原文）"

    Search the CellMarker 2.0 database for cancer-specific cell markers. Returns marker genes and cell types found in cancer tissue contexts. Use this to identify tumor-specific markers, immune cell markers in the tumor microenvironment, or cancer stem cell markers. Filters results to 'Cancer cell' annotations only. At least one search parameter is required: cancer_type (tissue), gene_symbol, or cell_type.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CellMarker_search_cancer_markers', 'arguments': {}},
)
```

## `CELLxGENE_get_expression_data` {#CELLxGENE_get_expression_data}

Queries gene expression data from CELLxGENE Census with cell and gene filters.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**无
- **关键词：**single-cell expression, census, expression matrix
- **需要 API 密钥：**`CELLXGENE_CENSUS_PACKAGE_INSTALLED`
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`get_anndata` 默认 `get_anndata` |
| `organism` | string | 否 | Organism name 可选值：`Homo sapiens`, `Mus musculus` 默认 `Homo sapiens` |
| `obs_value_filter` | string | 否 | REQUIRED (or use var_value_filter) - Cell filter. Common values: tissue_general (brain, lung, heart, blood, liver), cell_type (T cell, B cell, neuron, macrophage), disease (normal, COVID-19, Alzheimer disease). Example: 'tissue_general == "lung" and disease == "normal"' |
| `var_value_filter` | string | 否 | REQUIRED (or use obs_value_filter) - Gene filter by symbol or Ensembl ID. Examples: 'feature_name in \["TP53", "BRCA1", "EGFR"\]', 'feature_name == "CD4"'. Common genes: TP53, EGFR, BRCA1, CD4, CD8A, IL6, TNF, GAPDH. |
| `obs_column_names` | array of string | 否 | Cell metadata columns to include |
| `var_column_names` | array of string | 否 | Gene metadata columns to include |
| `census_version` | string | 否 | Census version to query. 'stable' (recommended, Long-Term Support release), 'latest' (newest data, may change), or specific date 'YYYY-MM-DD' (for reproducibility). Default 'stable' is best for production analyses. 默认 `stable` |

??? note "上游说明（英文原文）"

    Query gene expression data from CELLxGENE Census (50M+ cells, 60K+ genes). CRITICAL: At least one filter REQUIRED - unfiltered queries timeout. Start with: obs_value_filter='tissue_general == "lung" and disease == "normal"' AND var_value_filter='feature_name in \["TP53", "EGFR"\]'. Returns AnnData summary (dimensions, metadata). For full matrix data, use Python API with sufficient memory. Note: requires the cellxgene-census Python package; install with `pip install cellxgene-census` for live access.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CELLxGENE_get_expression_data', 'arguments': {}},
)
```

## `CELLxGENE_get_gene_metadata` {#CELLxGENE_get_gene_metadata}

Queries gene metadata such as symbols, Ensembl IDs and summary statistics from CELLxGENE Census.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**无
- **关键词：**gene metadata, var, ensembl id, feature
- **需要 API 密钥：**`CELLXGENE_CENSUS_PACKAGE_INSTALLED`
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`get_var_metadata` 默认 `get_var_metadata` |
| `organism` | string | 否 | Organism name 可选值：`Homo sapiens`, `Mus musculus` 默认 `Homo sapiens` |
| `var_value_filter` | string | 否 | Filter genes using SQL-like syntax. Format: 'field == "value"'. Common fields: feature_name (gene symbol), feature_id (Ensembl ID), feature_biotype. Examples: 'feature_name == "TP53"', 'feature_name in \["CD4", "CD8A"\]', 'feature_biotype == "protein_coding"'. |
| `column_names` | array of string | 否 | Specific columns to return (default: all columns) |
| `census_version` | string | 否 | Census version to query. 'stable' (recommended, Long-Term Support release), 'latest' (newest data, may change), or specific date 'YYYY-MM-DD' (for reproducibility). Default 'stable' is best for production analyses. 默认 `stable` |

??? note "上游说明（英文原文）"

    Query gene (variable) metadata from CELLxGENE Census. Returns gene symbols, Ensembl IDs, feature types, and summary statistics. Filter by gene name, feature type, or expression characteristics. Use for: finding genes of interest, validating gene symbols, exploring feature presence across datasets. Note: requires the cellxgene-census Python package; install with `pip install cellxgene-census` for live access.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CELLxGENE_get_gene_metadata', 'arguments': {}},
)
```

## `CELLxGENE_get_presence_matrix` {#CELLxGENE_get_presence_matrix}

Returns the feature presence matrix of genes measured per CELLxGENE dataset.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**无
- **关键词：**presence matrix, gene coverage, datasets
- **需要 API 密钥：**`CELLXGENE_CENSUS_PACKAGE_INSTALLED`
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`get_presence_matrix` 默认 `get_presence_matrix` |
| `organism` | string | 否 | Organism name 可选值：`Homo sapiens`, `Mus musculus` 默认 `Homo sapiens` |
| `census_version` | string | 否 | Census version to query. 'stable' (recommended, Long-Term Support release), 'latest' (newest data, may change), or specific date 'YYYY-MM-DD' (for reproducibility). Default 'stable' is best for production analyses. 默认 `stable` |

??? note "上游说明（英文原文）"

    Get feature presence matrix showing which genes are measured in which datasets. Returns sparse matrix dimensions and density. Useful for understanding data completeness and selecting datasets with genes of interest. Use for: checking gene coverage, dataset selection, quality assessment. Note: requires the cellxgene-census Python package; install with `pip install cellxgene-census` for live access.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CELLxGENE_get_presence_matrix', 'arguments': {}},
)
```

## `Coexpression_modules` {#Coexpression_modules}

Detects WGCNA-style co-expression modules from a gene expression matrix, returning module genes and eigengenes.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`expression`
- **关键词：**wgcna, co-expression, module eigengene, gene network, correlation
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `expression` | object | 是 | \{gene: \[value_per_sample, ...\]\} — a genes x samples expression matrix (log-normalized). |
| `power` | integer \| null | 否 | Soft-threshold exponent on \|correlation\| (default 6, the WGCNA scale-free default). |
| `correlation_threshold` | number \| null | 否 | Minimum \|correlation\| to connect two genes (default 0.5). |
| `min_module_size` | integer \| null | 否 | Minimum genes for a module to be reported (default 5). |

??? note "上游说明（英文原文）"

    Detect co-expression modules from a genes x samples expression matrix (WGCNA-style): builds a soft-thresholded correlation network and partitions it into modules of co-regulated genes by modularity, returning each module's genes and its module eigengene (first PC = per-sample summary profile). Pure-compute (NumPy + NetworkX). Simplified WGCNA (no TOM / dynamic tree cut).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Coexpression_modules', 'arguments': {'expression': {}}},
)
```

## `DeepSpotM_predict_gene_expression` {#DeepSpotM_predict_gene_expression}

Predicts per-gene spatial expression for an H&amp;E tile locally with the DeepSpot-M model.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`image_path`, `genes`
- **关键词：**spatial transcriptomics, h&amp;e, histology, virtual expression
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `image_path` | string | 是 | Path to a local 224x224 H&amp;E histology tile (.png, .tif, .tiff, .jpg, .jpeg). Must be cut at roughly 20x (~0.5 microns per pixel). |
| `genes` | any | 是 | Gene symbol(s) to predict, for example \['EPCAM', 'CD3D', 'PTPRC'\]. Required: the full ~19k-gene panel is deliberately not returned in one call. |
| `source` | string | 否 | Which frozen gene-embedding source to query the decoder with. Default 'scgpt'. 可选值：`evo2`, `orthrus`, `prott5`, `scgpt`, `apertus` |
| `device` | string | 否 | Torch device, for example 'cuda' or 'cpu'. Default 'auto', which uses CUDA when available. |
| `model_repo` | string | 否 | Hugging Face repo or local directory holding the checkpoint. Defaults to 'ratschlab/DeepSpotM'. |

??? note "上游说明（英文原文）"

    Predict spatial gene expression for a single H&amp;E histology tile LOCALLY using the DeepSpot-M multimodal foundation model (no API, no key), giving a virtual spatial transcriptomics readout without running a spatial assay. Given a local path to a 224x224 H&amp;E image and a list of HGNC gene symbols, returns expression for those genes in log1p-CPM units. DeepSpot-M represents genes as queryable embeddings drawn from DNA, RNA, protein, single-cell and text foundation models rather than as fixed outputs, so it spans the protein-coding transcriptome and covers genes unseen during training. Requesting specific genes is much cheaper than the full ~19k panel because only the requested gene queries run through the cross-attention decoder, so pass the genes you need. The tile must be exactly 224x224 px cut at roughly 20x (~0.5 microns per pixel), the magnification the model was trained on. Requires the optional 'deepspotm' package (pip install deepspotm) and access to its weights, which are gated on the Hugging Face Hub: request access at https://huggingface.co/ratschlab/DeepSpotM and run 'huggingface-cli login' first. The model is loaded lazily on first call and cached. Code is PolyForm Noncommercial 1.0.0 and weights are CC-BY-NC-SA-4.0.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DeepSpotM_predict_gene_expression',
     'arguments': {'image_path': '...', 'genes': '...'}},
)
```

## `expression_anova_per_gene` {#expression_anova_per_gene}

Runs per-gene ANOVA or log2 fold change on a gene x sample expression matrix.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`counts_file`, `meta_file`, `group_col`, `mode`
- **关键词：**anova, fold change, expression matrix, groups
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `counts_file` | string | 是 | Path to CSV count matrix (genes x samples, first column = gene IDs). |
| `meta_file` | string | 是 | Path to CSV metadata (rows=samples, contains group_col). |
| `group_col` | string | 是 | Metadata column identifying groups (e.g. 'Strain', 'cell_type'). |
| `mode` | string | 是 | 'anova': one-way ANOVA across groups. 'fold_change': median per-gene log2FC between group_a and group_b. 可选值：`anova`, `fold_change` |
| `exclude_groups` | array \| null | 否 | Groups in group_col to exclude before analysis. |
| `group_a` | string \| null | 否 | Numerator group for fold_change mode. |
| `group_b` | string \| null | 否 | Denominator group for fold_change mode. |

??? note "上游说明（英文原文）"

    Per-gene ANOVA or log2 fold-change on a gene x sample expression matrix. Each gene contributes one value per group (its mean expression across group samples), then mode='anova' runs one-way ANOVA (f_oneway) across K groups of N gene-level means; mode='fold_change' computes per-gene log2((mean_a+1)/(mean_b+1)) between two groups and returns median/mean across genes. Use for gene-level comparisons across experimental conditions where each gene is treated as one observation (e.g. comparing strains, cell types, or treatments). Accepts CSV counts matrix (rows=genes, cols=samples) and CSV metadata (rows=samples with group_col). Auto-aligns by overlapping sample IDs (transposes counts if needed).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'expression_anova_per_gene',
     'arguments': {'counts_file': '...',
                   'meta_file': '...',
                   'group_col': '...',
                   'mode': 'anova'}},
)
```

## `ExpressionAtlas_get_baseline` {#ExpressionAtlas_get_baseline}

Lists baseline expression experiments in EBI Expression Atlas for a species.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`gene`
- **关键词：**expression atlas, baseline, tissue expression
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 是 | Gene symbol (e.g., 'TP53', 'WDR7') or Ensembl ID (e.g., 'ENSG00000141510') |
| `species` | string | 否 | Species name (default: 'homo sapiens'). Also supports 'mus musculus', 'rattus norvegicus', etc. 默认 `homo sapiens` |
| `limit` | integer \| null | 否 | Maximum number of experiments to return (default 50, max 500). Compare with the total reported in the response to see how many exist. 默认 `50` |
| `offset` | integer \| null | 否 | Zero-based index of the first experiment to return, for paging through results beyond the first page (e.g. offset=50 returns 51-100). Default: 0. 默认 `0` |

??? note "上游说明（英文原文）"

    List baseline gene expression experiments in EBI Expression Atlas for a given species. The `gene` argument is used to label the request and suggest a follow-up call, but does NOT filter the returned experiment list (the underlying API has no reliable way to filter experiments by gene). Do NOT use GxA_get_experiment_expression's gene_id filter to check whether a specific experiment has data for a gene -- that endpoint only returns a small default sample of the experiment's genes (verified live: ~29 of 9,570 for E-MTAB-2836) and silently ignores gene-filter parameters, so a miss there is not evidence of absence. For a reliable per-gene baseline expression check, use GTEx_get_expression_summary (human tissues, verified working) or HPA_get_rna_expression_in_specific_tissues instead. Complements GTEx and HPA. No API key needed.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ExpressionAtlas_get_baseline',
     'arguments': {'gene': '...'}},
)
```

## `ExpressionAtlas_get_experiment` {#ExpressionAtlas_get_experiment}

Gets metadata for an Expression Atlas experiment by accession.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`accession`
- **关键词：**expression atlas, experiment, e-mtab
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `accession` | string | 是 | Expression Atlas experiment accession (e.g., 'E-MTAB-2836', 'E-GEOD-26284') |

??? note "上游说明（英文原文）"

    Get metadata for a specific Expression Atlas experiment by accession (e.g., E-MTAB-2836): accession, type, species, and description/URLs. Caveat: the underlying GXA API does not expose per-experiment design/technology/contrasts/PubMed/assay-count data through this endpoint (confirmed live: those fields are always empty) -- for that level of detail, visit the experiment's page on https://www.ebi.ac.uk/gxa directly.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ExpressionAtlas_get_experiment',
     'arguments': {'accession': '...'}},
)
```

## `ExpressionAtlas_search_differential` {#ExpressionAtlas_search_differential}

Searches Expression Atlas differential expression experiments by condition and gene text.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**无
- **关键词：**expression atlas, differential expression, condition
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 否 | Gene symbol or Ensembl ID (e.g., 'TP53', 'ENSG00000141510') |
| `condition` | string | 否 | Condition/disease to filter by (e.g., 'cancer', 'inflammation', 'breast') |
| `species` | string | 否 | Species name (default: 'homo sapiens') 默认 `homo sapiens` |
| `limit` | integer \| null | 否 | Maximum number of experiments to return (default 50, max 500). Compare with the total reported in the response to see how many exist. 默认 `50` |
| `offset` | integer \| null | 否 | Zero-based index of the first experiment to return, for paging through results beyond the first page (e.g. offset=50 returns 51-100). Default: 0. 默认 `0` |

??? note "上游说明（英文原文）"

    Search for differential expression EXPERIMENTS from EBI Expression Atlas by condition text and/or a gene mentioned in the experiment description. Supports disease-context expression analysis. Caveat: `condition` does real text filtering, but `gene` is a text-search coincidence (experiment descriptions rarely name individual genes) and does not return actual per-gene differential-expression values -- treat gene-tagged results as a starting point for manual review, not a validated gene-specific result set.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ExpressionAtlas_search_differential', 'arguments': {}},
)
```

## `ExpressionAtlas_search_experiments` {#ExpressionAtlas_search_experiments}

Searches all Expression Atlas experiments by condition and gene text.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**无
- **关键词：**expression atlas, experiments, search
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 否 | Gene symbol or Ensembl ID |
| `condition` | string | 否 | Biological condition or tissue (e.g., 'brain', 'cancer', 'kidney') |
| `species` | string | 否 | Species name (optional, e.g., 'homo sapiens') |
| `limit` | integer \| null | 否 | Maximum number of experiments to return (default 50, max 500). Compare with the total reported in the response to see how many exist. 默认 `50` |
| `offset` | integer \| null | 否 | Zero-based index of the first experiment to return, for paging through results beyond the first page (e.g. offset=50 returns 51-100). Default: 0. 默认 `0` |

??? note "上游说明（英文原文）"

    Search all EBI Expression Atlas experiments (baseline + differential) by condition text and/or a gene mentioned in the experiment description. Returns experiment accessions, types, descriptions, and assay counts. Caveat: `condition` does real text filtering, but `gene` is a text-search coincidence (experiment descriptions rarely name individual genes), so a gene-only query effectively returns the unfiltered catalog -- use to find relevant expression datasets, not as a validated per-gene result set.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ExpressionAtlas_search_experiments', 'arguments': {}},
)
```

## `FlyBase_get_gene_expression` {#FlyBase_get_gene_expression}

Returns Drosophila gene expression annotations by developmental stage and anatomical location from the Alliance.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`gene_id`
- **关键词：**drosophila, expression pattern, developmental stage, anatomy
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_id` | string | 是 | FlyBase gene ID, with or without the 'FB:' prefix. Examples: 'FB:FBgn0000490' (dpp), 'FB:FBgn0003996' (white/w), 'FB:FBgn0004647' (Notch/N). |
| `limit` | integer | 否 | Max annotations to return (default 50, max 100). 默认 `50` |

??? note "上游说明（英文原文）"

    Get expression annotations for a Drosophila gene from the Alliance of Genome Resources. Returns per-annotation records of developmental stage (whenExpressed) and anatomical location (whereExpressed) for the gene, plus the total annotation count. Example: 'FB:FBgn0000490' (dpp) has 317 expression annotations such as 'embryonic stage 14' in 'abdominal transverse muscle cell'. Useful for understanding when and where a fly gene is expressed.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'FlyBase_get_gene_expression',
     'arguments': {'gene_id': '...'}},
)
```

## `GDC_get_gene_expression` {#GDC_get_gene_expression}

Queries GDC/TCGA RNA-seq gene expression quantification file metadata by project and gene ID.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`project_id`
- **关键词：**rna-seq, tcga, expression quantification, cancer
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `project_id` | string | 是 | GDC project (e.g., 'TCGA-BRCA', 'TCGA-LUAD', 'TCGA-GBM') |
| `gene_id` | string | 否 | Optional: Ensembl gene ID (e.g., 'ENSG00000141510' for TP53) |
| `size` | integer | 否 | Number of results 默认 `10` |

??? note "上游说明（英文原文）"

    Query RNA-Seq gene expression files from GDC/TCGA projects. Returns file metadata for expression quantification data. Use to find expression data for downstream analysis.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GDC_get_gene_expression',
     'arguments': {'project_id': '...'}},
)
```

## `GeneNetwork_get_sample_data` {#GeneNetwork_get_sample_data}

Gets expression or phenotype values with SEs for a gene/trait across strains of a GeneNetwork dataset.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`dataset_name`, `trait_name`
- **关键词：**genenetwork, bxd, inbred strains, systems genetics, expression values
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `dataset_name` | string | 是 | Dataset long abbreviation (e.g., 'HC_M2_0606_P' for BXD hippocampus, 'BXDPublish' for phenotypes). Obtain from GeneNetwork_list_datasets Long_Abbreviation field. |
| `trait_name` | string | 是 | Probe/trait identifier within the dataset (e.g., '1436869_at' for a microarray probe, '10001' for a BXD phenotype trait). Must exist in the specified dataset. |

??? note "上游说明（英文原文）"

    Get expression or phenotype values for a specific gene/trait across all samples in a GeneNetwork dataset. Returns measurement values with standard errors for each inbred strain or individual in the genetic cross. Use this to find how a gene's expression varies across strains, identify strain differences, or get data for QTL analysis. Dataset abbreviation and trait/probe ID must match the dataset.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GeneNetwork_get_sample_data',
     'arguments': {'dataset_name': '...', 'trait_name': '...'}},
)
```

## `GSVA_score` {#GSVA_score}

Computes GSVA gene set by sample activity scores from a genes x samples expression matrix.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`expression`
- **关键词：**gsva, gene set variation, expression matrix, pathway activity
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `expression` | object | 是 | \{gene: \[value_per_sample, ...\]\} — a genes x samples expression matrix on a continuous (e.g. log-normalized) scale. |
| `samples` | array \| null | 否 | Sample names (same order as the per-gene value lists); defaults to sample_1.. |
| `gene_sets` | object \| null | 否 | \{set_name: \[gene, ...\]\} collection to score. Alternative to gene_set. |
| `gene_set` | array \| null | 否 | A single gene set (list of gene symbols). |
| `tau` | number \| null | 否 | Rank-weighting exponent in the random walk (default 1.0, the GSVA standard). |
| `mx_diff` | boolean \| null | 否 | True (default) = ES_pos+ES_neg difference score (signed, bimodal, GSVA default); False = classic single max-deviation KS statistic. |

??? note "上游说明（英文原文）"

    Gene Set Variation Analysis (Hänzelmann 2013): turn an expression matrix (genes x samples) into a label-free gene_set x sample pathway-activity matrix. Unlike ssGSEA, GSVA first converts each gene to a kernel-CDF quantile across samples, then runs a symmetric KS walk with a max-deviation difference score (mx_diff), giving signed scores centered near zero (positive = set up-regulated in that sample vs the cohort). Pure-compute (no R); supply log-scale expression + gene sets from MSigDB/Enrichr.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GSVA_score', 'arguments': {'expression': {}}},
)
```

## `GTEx_get_expression_summary` {#GTEx_get_expression_summary}

Summarizes a gene's tissue-specific expression (median TPM) across GTEx tissues.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**无
- **关键词：**gtex, tissue expression, tpm, baseline expression
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 否 | Gene symbol (e.g., 'TP53', 'BRCA1', 'FBN1'). Auto-resolved to versioned GENCODE ID. |
| `ensembl_gene_id` | string | 否 | Ensembl gene identifier (e.g., 'ENSG00000141510' or 'ENSG00000141510.11'). Use gene_symbol instead if you have a gene name. |

??? note "上游说明（英文原文）"

    Summarize tissue-specific expression (e.g., median TPM) for a gene across GTEx tissues. Accepts gene symbols (TP53, BRCA1) or Ensembl IDs. Use to profile baseline expression patterns for targets/biomarkers.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GTEx_get_expression_summary', 'arguments': {}},
)
```

## `GTEx_get_gene_expression` {#GTEx_get_gene_expression}

Gets per-sample normalized TPM expression for a gene in GTEx tissues, optionally by sex/age subset.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`gencode_id`
- **关键词：**gtex, tpm, sample-level, expression distribution
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`get_gene_expression` 默认 `get_gene_expression` |
| `gencode_id` | string \| array | 是 | Required: Versioned GENCODE ID(s) |
| `tissue_site_detail_id` | array of string | 否 | Optional: Filter by tissues |
| `attribute_subset` | string | 否 | Optional: Subset by donor sex or age bracket 可选值：`sex`, `ageBracket` |
| `dataset_id` | string | 否 | 可选值：`gtex_v8`, `gtex_v10`, `gtex_snrnaseq_pilot` 默认 `gtex_v8` |
| `page` | integer | 否 | 默认 `0` |
| `items_per_page` | integer | 否 | 默认 `250` |

??? note "上游说明（英文原文）"

    Get gene expression data at individual sample level (not aggregated). Returns normalized expression values (TPM) for each sample, allowing detailed distribution analysis. Can subset by donor attributes (sex, age_bracket) for demographic-specific analysis. Much more data than median expression. Use for: expression variability analysis, outlier detection, sample-level QC, demographic comparisons. Warning: large result sets for multiple genes/tissues.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GTEx_get_gene_expression',
     'arguments': {'gencode_id': '...'}},
)
```

## `GTEx_get_median_gene_expression` {#GTEx_get_median_gene_expression}

Gets median TPM expression for genes across GTEx tissues.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**无
- **关键词：**gtex, median tpm, tissue expression
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`get_median_gene_expression` 默认 `get_median_gene_expression` |
| `gencode_id` | string \| array | 否 | Gene identifier(s): gene symbol (e.g. 'TP53'), unversioned Ensembl ID (e.g. 'ENSG00000141510'), or versioned GENCODE ID (e.g. 'ENSG00000141510.18'). Auto-resolved to versioned GENCODE ID. Can be single string or array. |
| `tissue_site_detail_id` | array of string | 否 | Optional: Tissue IDs to filter (e.g. \['Liver', 'Brain_Cortex'\]). Omit for all tissues. See GTEx_get_tissue_sites for valid IDs |
| `dataset_id` | string | 否 | GTEx dataset version (default: gtex_v8, the most widely cited release). gtex_v10 is the newer release (943 eQTL donors vs 838) and is fully queryable -- gene IDs are resolved to each dataset's own GENCODE version automatically. 可选值：`gtex_v8`, `gtex_v10`, `gtex_snrnaseq_pilot` 默认 `gtex_v8` |
| `page` | integer | 否 | Page number for pagination (0-based) 默认 `0` |
| `items_per_page` | integer | 否 | Results per page 默认 `250` |
| `gene_symbol` | string | 否 | Gene symbol alias for gencode_id (e.g., "TP53", "COL5A1") |

??? note "上游说明（英文原文）"

    Get median gene expression levels across GTEx tissues. Returns median expression in TPM (Transcripts Per Million) for one or more genes across 54 tissue sites. Filter by specific tissues or get expression across all tissues. Defaults to GTEx V8 (2020) -- confirmed live no V11 dataset is queryable via this API (dataset_id accepts gtex_v7/gtex_v8/gtex_v10). Set dataset_id='gtex_v10' for the newer release; gene IDs are resolved to each dataset's own GENCODE version automatically, so V10 returns full results. Use for: tissue-specific expression profiling, identifying highly expressed genes, comparing expression patterns across tissues. Example: Get brain vs liver expression for TP53 gene.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GTEx_get_median_gene_expression', 'arguments': {}},
)
```

## `GTEx_get_median_transcript_expression` {#GTEx_get_median_transcript_expression}

Gets per-transcript median TPM across GTEx tissues for isoform usage analysis.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**无
- **关键词：**isoform, transcript expression, enst, gtex, tpm
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`get_median_transcript_expression` 默认 `get_median_transcript_expression` |
| `gencode_id` | string \| array | 否 | Gene identifier(s): gene symbol (e.g. 'BRCA1'), Ensembl ID, or versioned GENCODE ID (e.g. 'ENSG00000012048.20'). Auto-resolved. |
| `gene_symbol` | string | 否 | Gene symbol alias for gencode_id. |
| `tissue_site_detail_id` | array of string | 否 | Optional tissue ID(s) to filter (e.g. \['Whole_Blood'\]). Omit for all tissues. See GTEx_get_tissue_sites. |
| `page` | integer | 否 | Page number (0-based) 默认 `0` |
| `items_per_page` | integer | 否 | Results per page 默认 `250` |
| `dataset_id` | string | 否 | GTEx dataset version (default gtex_v8) 可选值：`gtex_v8`, `gtex_v10`, `gtex_snrnaseq_pilot` 默认 `gtex_v8` |

??? note "上游说明（英文原文）"

    Get GTEx transcript/isoform-level median expression: per-transcript (ENST) median TPM across tissues, enabling isoform-usage and transcript-switching analysis that gene-level expression cannot reveal. Based on gtex_v8. Example: BRCA1 (ENSG00000012048.20) in Whole_Blood -&gt; transcript ENST00000352993.7 median 0.37 TPM, ENST00000354071.7 median 0.0, ENST00000357654.7 median 0.07.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GTEx_get_median_transcript_expression', 'arguments': {}},
)
```

## `GTEx_get_single_nucleus_expression` {#GTEx_get_single_nucleus_expression}

Gets GTEx snRNA-seq gene expression summarized by cell type (mean/median, cell counts).

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**无
- **关键词：**snrna-seq, cell type, single nucleus, gtex
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`get_single_nucleus_expression` 默认 `get_single_nucleus_expression` |
| `result_type` | string | 否 | 'detail' for per-cell-type expression stats (default), 'summary' for tissue x cell-type cell-count summary. 可选值：`detail`, `summary` 默认 `detail` |
| `gencode_id` | string \| array | 否 | Gene identifier(s): gene symbol (e.g. 'BRCA1'), Ensembl ID, or versioned GENCODE ID (e.g. 'ENSG00000012048.20'). Auto-resolved. |
| `gene_symbol` | string | 否 | Gene symbol alias for gencode_id. |
| `tissue_site_detail_id` | array of string | 否 | Optional tissue ID(s) (e.g. \['Muscle_Skeletal'\]). Omit for all snRNA-seq tissues. |
| `dataset_id` | string | 否 | GTEx dataset version (default gtex_snrnaseq_pilot; snRNA-seq data only exists in the pilot dataset). 可选值：`gtex_v8`, `gtex_v10`, `gtex_snrnaseq_pilot` 默认 `gtex_snrnaseq_pilot` |
| `page` | integer | 否 | Page number (0-based) 默认 `0` |
| `items_per_page` | integer | 否 | Results per page 默认 `250` |

??? note "上游说明（英文原文）"

    Get GTEx single-nucleus (snRNA-seq) gene expression resolved by cell type from the snRNA-seq pilot dataset (bulk-tissue GTEx cannot resolve cell types). result_type='detail' (default) returns per-cell-type mean/median (with and without zeros), cell counts, and number of zeros for a gene. result_type='summary' returns a tissue x cell-type cell-count summary (numCells per cell type). Uses datasetId gtex_snrnaseq_pilot. Example: BRCA1 (ENSG00000012048.20) Muscle_Skeletal -&gt; Myocyte (NMJ-rich) count 6 meanWithoutZeros 1.455; Endothelial cell (vascular) count 38 meanWithoutZeros 2.264.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GTEx_get_single_nucleus_expression', 'arguments': {}},
)
```

## `GTEx_get_top_expressed_genes` {#GTEx_get_top_expressed_genes}

Gets genes ranked by median expression in a GTEx tissue, optionally excluding mitochondrial genes.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`tissue_site_detail_id`
- **关键词：**gtex, top genes, tissue signature, tpm
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`get_top_expressed_genes` 默认 `get_top_expressed_genes` |
| `tissue_site_detail_id` | string | 是 | Required: Tissue ID (e.g. 'Liver', 'Brain_Cortex'). Use GTEx_get_tissue_sites to find valid IDs |
| `filter_mt_genes` | boolean | 否 | Exclude mitochondrial genes from results 默认 `True` |
| `dataset_id` | string | 否 | 可选值：`gtex_v8`, `gtex_v10`, `gtex_snrnaseq_pilot` 默认 `gtex_v8` |
| `page` | integer | 否 | 默认 `0` |
| `items_per_page` | integer | 否 | 默认 `250` |

??? note "上游说明（英文原文）"

    Get top expressed genes in a specific tissue sorted by median expression. Returns gene list with expression levels (TPM), useful for understanding tissue-specific gene signatures. Option to filter mitochondrial genes (default: true). Use for: tissue characterization, finding housekeeping genes, identifying tissue markers, quality control. Essential for understanding what genes are active in each tissue type.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GTEx_get_top_expressed_genes',
     'arguments': {'tissue_site_detail_id': '...'}},
)
```

## `GxA_get_experiment_expression` {#GxA_get_experiment_expression}

Gets gene expression levels across tissues/conditions from an EBI Expression Atlas experiment.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`experiment_accession`
- **关键词：**expression atlas, baseline expression, differential, ebi
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `experiment_accession` | string | 是 | Expression Atlas experiment accession (e.g., 'E-MTAB-2836' for human tissues RNA-seq, 'E-MTAB-2706' for cancer cell lines). |
| `gene_id` | string | 否 | Optional Ensembl gene ID to look for within the small default SAMPLE of genes the upstream API returns for this experiment (e.g. 'ENSG00000141510' for TP53). This is NOT a true server-side or full-experiment filter -- the upstream endpoint ignores gene-filter parameters and only ever returns a default sample (commonly ~20-30 genes) of the experiment's full gene set (which can be thousands). A mis… |

??? note "上游说明（英文原文）"

    Get gene expression data from a specific Expression Atlas experiment. Returns expression levels across tissues or conditions. IMPORTANT limitation (verified live): the upstream /experiments/\{accession\} endpoint always returns only a small, arbitrary default SAMPLE of the experiment's genes (e.g. ~29 of 9,570 for E-MTAB-2836) and silently ignores every gene-filter parameter we can send it -- so the optional `gene_id` filter can only ever match within that tiny sample, never the full gene set. The response discloses `profiles_returned_by_upstream` and `profiles_available_upstream` so callers can see the sample size vs. the true total; when a `gene_id` filter matches nothing, `coverage_warning` explains that this is NOT evidence the gene is unexpressed -- it may simply have fallen outside the sample. For a reliable per-gene answer, use GTEx_get_expression_summary (human tissues) or HPA_get_rna_expression_in_specific_tissues instead. Experiment E-MTAB-2836 contains RNA-seq of 122 human tissues, E-MTAB-2706 covers 675 cancer cell lines.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GxA_get_experiment_expression',
     'arguments': {'experiment_accession': '...'}},
)
```

## `GxA_get_experiment_info` {#GxA_get_experiment_info}

Gets Expression Atlas experiment metadata: species, gene and condition counts, design.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`experiment_accession`
- **关键词：**expression atlas, experiment design, ebi
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `experiment_accession` | string | 是 | Expression Atlas experiment accession (e.g., 'E-MTAB-2836', 'E-MTAB-2706', 'E-GEOD-46817'). |

??? note "上游说明（英文原文）"

    Get metadata about a Gene Expression Atlas experiment including species, number of genes and conditions, and experimental design. Useful for understanding what data is available before querying specific gene expression. Example: E-MTAB-2836 has 29 conditions (tissues) and profiles for thousands of genes in Homo sapiens.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GxA_get_experiment_info',
     'arguments': {'experiment_accession': '...'}},
)
```

## `GxA_list_experiments` {#GxA_list_experiments}

Lists EBI Expression Atlas experiments filtered by species and baseline/differential type.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**无
- **关键词：**expression atlas, baseline, differential, experiments
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `species` | string | 否 | Species name to filter experiments (e.g., 'Homo sapiens', 'Mus musculus', 'Arabidopsis thaliana'). Leave empty for all species. |
| `experiment_type` | string | 否 | Filter by experiment type: 'baseline' (tissue/cell type expression), 'differential' (condition comparisons), 'proteomics', or 'rnaseq'. Leave empty for all types. |
| `limit` | integer | 否 | Maximum number of experiments to return (1-100). Default: 20. |

??? note "上游说明（英文原文）"

    List gene expression experiments from EBI Gene Expression Atlas. Can filter by species (e.g., 'Homo sapiens') and experiment type ('baseline' for tissue expression, 'differential' for condition comparisons). The Atlas contains 4,500+ experiments across hundreds of species. Example: filtering by 'Homo sapiens' and 'baseline' returns 123 human tissue-level expression experiments including E-MTAB-2836 (122 human tissues), E-MTAB-2706 (675 cancer cell lines), and E-GEOD-46817 (melanoma).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GxA_list_experiments', 'arguments': {}},
)
```

## `HPA_get_cancer_prognostics_by_gene` {#HPA_get_cancer_prognostics_by_gene}

Gets whether a gene's expression correlates with survival across cancer types in HPA.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`ensembl_id`
- **关键词：**prognostic gene, cancer survival, pathology atlas, expression
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ensembl_id` | string | 是 | Ensembl Gene ID of the gene to check, e.g., 'ENSG00000141510' for TP53, 'ENSG00000012048' for BRCA1. |

??? note "上游说明（英文原文）"

    Retrieve prognostic value of a gene across various cancer types, indicating if its expression level correlates with patient survival outcomes. Reports both the significant results ('prognostic_summary') and the cancer types that were assessed but found NOT prognostic ('non_prognostic_summary'), plus 'cancers_assessed_count' as the denominator, so that a negative result can be distinguished from missing data.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'HPA_get_cancer_prognostics_by_gene',
     'arguments': {'ensembl_id': '...'}},
)
```

## `HPA_get_comparative_expression_by_gene_and_cellline` {#HPA_get_comparative_expression_by_gene_and_cellline}

Compares a gene's expression between a cell line and healthy tissues in HPA.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`gene_name`, `cell_line`
- **关键词：**cell line, expression comparison, human protein atlas
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_name` | string | 是 | Gene name or gene symbol, e.g., 'TP53', 'BRCA1', 'EGFR', etc. |
| `cell_line` | string | 是 | Cell line name, supported cell lines include: ishikawa, hela, mcf7, a549, hepg2, jurkat, pc3, rh30, siha, u251. |

??? note "上游说明（英文原文）"

    Compare the expression level differences of a gene between a specific cell line and healthy tissues using gene name and cell line name.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'HPA_get_comparative_expression_by_gene_and_cellline',
     'arguments': {'gene_name': '...', 'cell_line': '...'}},
)
```

## `HPA_get_contextual_biological_process_analysis` {#HPA_get_contextual_biological_process_analysis}

Analyzes a gene's biological processes in a tissue or cell line context by combining annotation with expression.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`gene_name`, `context_name`
- **关键词：**human protein atlas, tissue context, gene function, expression
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_name` | string | 是 | Gene name or symbol, e.g., 'TP53', 'EGFR', 'BRCA1'. |
| `context_name` | string | 是 | Name of the tissue or cell line to provide context, e.g., 'brain', 'liver', 'hela', 'mcf7'. |

??? note "上游说明（英文原文）"

    Analyze a gene's biological processes in the context of a specific tissue or cell line by integrating functional annotation with expression data to determine functional relevance.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'HPA_get_contextual_biological_process_analysis',
     'arguments': {'gene_name': '...', 'context_name': '...'}},
)
```

## `HPA_get_disease_expression_by_gene_tissue_disease` {#HPA_get_disease_expression_by_gene_tissue_disease}

Compares a gene's expression in a disease state versus healthy tissue in HPA.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`gene_name`, `tissue_type`, `disease_name`
- **关键词：**disease expression, human protein atlas, tissue
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_name` | string | 是 | Gene name or gene symbol, e.g., 'TP53', 'BRCA1', 'KRAS', etc. |
| `tissue_type` | string | 是 | Tissue type, e.g., 'brain', 'breast', 'colon', 'lung', etc., optional parameter. |
| `disease_name` | string | 是 | Disease name, supported diseases include: brain_cancer, breast_cancer, colon_cancer, lung_cancer, liver_cancer, prostate_cancer, kidney_cancer, pancreatic_cancer, stomach_cancer, ovarian_cancer. |

??? note "上游说明（英文原文）"

    Compare the expression level of a gene in specific disease state versus healthy state using gene name, tissue type, and disease name.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'HPA_get_disease_expression_by_gene_tissue_disease',
     'arguments': {'gene_name': '...',
                   'tissue_type': '...',
                   'disease_name': '...'}},
)
```

## `HPA_get_rna_expression_by_source` {#HPA_get_rna_expression_by_source}

Gets a gene's RNA expression in a tissue, blood, brain or single-cell source from HPA.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`gene_name`, `source_type`, `source_name`
- **关键词：**rna expression, ntpm, blood, brain, single cell
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_name` | string | 是 | Gene name or gene symbol, e.g., 'GFAP', 'TP53', 'BRCA1', etc. |
| `source_type` | string | 是 | The type of biological source. Choose from: 'tissue' (51 consensus tissues, nTPM), 'blood' (19 immune cell subsets, nTPM), 'brain' (13 brain regions, nTPM), 'single_cell' (~160 single cell types, nCPM). |
| `source_name` | string | 是 | The specific source name, e.g. 'liver', 'heart_muscle', 'skin', 't_cell', 'hepatocyte', 'cerebellum'. Case-insensitive; spaces and hyphens are treated as underscores. Names HPA publishes only as subsets (blood 't_cell', 'b_cell', 'monocyte', 'dendritic_cell') resolve to a representative subset column and the response reports which one in a 'note'. An invalid name returns an error listing valid na… |

??? note "上游说明（英文原文）"

    Get RNA expression for a gene in a specific biological source, read from HPA's dedicated per-source download columns (t_RNA_*, blood_RNA_*, brain_RNA_*, sc_RNA_*). Supports source_type 'tissue', 'blood', 'brain' and 'single_cell'. Units differ by source: single_cell values are nCPM, all others nTPM -- the returned expression_unit is read from HPA's own column header, never assumed.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'HPA_get_rna_expression_by_source',
     'arguments': {'gene_name': '...',
                   'source_type': '...',
                   'source_name': '...'}},
)
```

## `HPA_get_rna_expression_in_specific_tissues` {#HPA_get_rna_expression_in_specific_tissues}

Gets a gene's RNA expression (nTPM) in specified HPA consensus tissues.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`ensembl_id`, `tissue_names`
- **关键词：**ntpm, tissue expression, human protein atlas
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ensembl_id` | string | 是 | Ensembl Gene ID for the gene, e.g., 'ENSG00000141510' for TP53. |
| `tissue_names` | array of string | 是 | List of tissue names to query, e.g., \['skin', 'liver', 'heart muscle', 'kidney'\]. Case-insensitive; underscores, spaces and hyphens are interchangeable. HPA has no RNA column for 'bronchus', 'nasopharynx', 'oral mucosa' or 'soft tissue' (protein/IHC only). |

??? note "上游说明（英文原文）"

    Query RNA expression levels (nTPM) for a specific gene in one or more user-specified tissues. Reads HPA's full per-tissue nTPM panel (the 'Tissue RNA - &lt;tissue&gt; \[nTPM\]' columns), so any of HPA's 51 published consensus tissues and brain regions can be queried, e.g. 'cerebral cortex', 'basal ganglia', 'midbrain', 'liver'. Common names are mapped to HPA's actual column names, including the ones HPA suffixes with '1' ('skin' -&gt; 'Skin 1', 'stomach' -&gt; 'Stomach 1', 'endometrium' -&gt; 'Endometrium 1'); when a name resolves to a differently-named HPA column the response says so in a 'note'. A tissue reported as 'Not found' means HPA publishes no RNA column for that name (it is a naming problem or a tissue HPA covers only with protein/IHC data), not that the gene lacks expression there.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'HPA_get_rna_expression_in_specific_tissues',
     'arguments': {'ensembl_id': '...', 'tissue_names': []}},
)
```

## `L1000FWD_sig_search` {#L1000FWD_sig_search}

Searches the LINCS L1000 connectivity map with up and down regulated gene sets to find similar or opposite signatures.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`up_genes`, `down_genes`
- **关键词：**connectivity map, l1000, gene signature, up genes, down genes
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`sig_search` 默认 `sig_search` |
| `up_genes` | array of string | 是 | List of up-regulated gene symbols (e.g., \['MYC', 'TP53'\]) |
| `down_genes` | array of string | 是 | List of down-regulated gene symbols (e.g., \['CDKN1A', 'RB1'\]) |
| `n_results` | integer | 否 | Number of top results to return (default 10) 默认 `10` |
| `mode` | string | 否 | Return signatures 'similar' to or 'opposite' of the query signature (default: 'similar') 可选值：`similar`, `opposite` 默认 `similar` |

??? note "上游说明（英文原文）"

    Transcriptomic signature query against the LINCS L1000 Connectivity Map via the L1000FWD Fireworks REST API (maayanlab.cloud/L1000FWD). Accepts two gene sets representing transcriptional changes in a condition of interest: genes increasing and genes decreasing. Internally issues a POST to /sig_search to register the query, then fetches ranked hits from /result/topn/\{result_id\}. Each returned hit is an L1000 perturbation record containing sig_id, combined score, p-value, and q-value. mode='similar' retrieves perturbations that activate the same transcriptional programme (e.g. disease mimics); mode='opposite' retrieves perturbations that invert it, the classical computational drug-repurposing approach. No authentication required.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'L1000FWD_sig_search',
     'arguments': {'up_genes': [], 'down_genes': []}},
)
```

## `LINCS_search_signatures` {#LINCS_search_signatures}

Searches LINCS for drug perturbation gene expression signatures by compound and cell line.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`drug_name`
- **关键词：**l1000, perturbation signature, cell line, dose
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `drug_name` | string | 是 | Drug or perturbagen name to search for. Case-sensitive, use lowercase for best results. Examples: 'vorinostat', 'imatinib', 'dexamethasone', 'trichostatin-a'. |
| `cell_line` | string | 否 | Optional cell line filter. Examples: 'MCF7', 'A549', 'HepG2', 'PC3'. Leave empty to search all cell lines. |
| `limit` | integer | 否 | Maximum number of signatures to return (1-100). Default: 20. |

??? note "上游说明（英文原文）"

    Search LINCS (Library of Integrated Network-Based Cellular Signatures) for drug perturbation gene expression signatures. Find L1000 and other assay signatures by drug/compound name, optionally filtered by cell line. Returns signature IDs, doses, time points, tissue, disease context, and PubChem IDs. LINCS contains 1.5M+ signatures from 431+ libraries covering thousands of drugs across hundreds of cell lines. Example: searching 'vorinostat' returns HDAC inhibitor signatures across multiple cell lines and doses.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'LINCS_search_signatures',
     'arguments': {'drug_name': '...'}},
)
```

## `OpenTargets_get_target_expression_by_ensemblID` {#OpenTargets_get_target_expression_by_ensemblID}

Retrieves baseline expression statistics for a gene across tissues and cell types from GTEx, Tabula Sapiens and others.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`ensemblId`
- **关键词：**gtex, tissue expression, tabula sapiens, baseline expression
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ensemblId` | string | 是 | Ensembl gene ID of the target (e.g., 'ENSG00000146648' for EGFR). |
| `size` | integer | 否 | Number of expression rows to return per page (default 250, max 3000). Targets carry ~1400 rows spread across all datasources, so use size=3000 to fetch the complete profile in one call; values above 3000 are clamped because the API rejects them. 默认 `250` |
| `index` | integer | 否 | 0-based page index (default 0). Use with `size` to walk the remaining rows when `baselineExpression.truncated` is true. 默认 `0` |

??? note "上游说明（英文原文）"

    Retrieve baseline expression summary statistics for a target across tissues and cell types (GTEx, Tabula Sapiens, and other bulk/single-cell sources) by Ensembl gene ID. Each row gives the data source, tissue (and cell type for single-cell data), value unit, and min/q1/median/q3/max distribution stats plus specificity/distribution scores. Results are PAGINATED: a target typically has ~1400 rows across sources and only `size` of them are returned per call (default 250, max 3000 per the API). `baselineExpression.count` is the total, `baselineExpression.returned` is how many came back, and `baselineExpression.truncated` is true when rows remain unfetched; a truncation note is added to metadata. Rows from a given datasource are interleaved throughout the full result, so a partial page is NOT a complete view of any one source — set size=3000 to retrieve every row before filtering client-side by datasourceId. Use to profile where a target is normally expressed.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_target_expression_by_ensemblID',
     'arguments': {'ensemblId': '...'}},
)
```

## `PanglaoDB_cell_types_for_gene` {#PanglaoDB_cell_types_for_gene}

Lists cell types for which a gene is a curated PanglaoDB marker, with sensitivity/specificity.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`gene`
- **关键词：**marker gene, cell type, scrna-seq
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 是 | Official human gene symbol, e.g. 'CD19', 'ALB', 'EPCAM', 'INS'. Case-insensitive. Required. |
| `species` | string \| null | 否 | Optional species filter: 'human' (Hs) or 'mouse' (Mm). If omitted, returns all species. |
| `limit` | integer \| null | 否 | Maximum number of cell-type records to return. Default 50. |

??? note "上游说明（英文原文）"

    Reverse cell-type lookup for a gene using PanglaoDB (panglaodb.se). Given an official human gene symbol (e.g. 'CD19', 'ALB', 'EPCAM'), return every cell type for which that gene is a curated marker, with organ, canonical-marker flag, and per-gene sensitivity/specificity for human and mouse. Use this to interpret a marker gene found in a scRNA-seq cluster (e.g. 'which cell type does CD19 indicate?'). Optionally filter by species ('human' or 'mouse'). Results are sorted by descending human specificity. Note: the marker set is the frozen 2020 PanglaoDB release (no longer updated) — reliable for canonical markers, but not for cell types/markers described since.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PanglaoDB_cell_types_for_gene',
     'arguments': {'gene': '...'}},
)
```

## `PanglaoDB_markers_for_cell_type` {#PanglaoDB_markers_for_cell_type}

Returns curated PanglaoDB marker genes for a cell type with sensitivity/specificity for human and mouse.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`cell_type`
- **关键词：**marker genes, scrna-seq, cell annotation
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `cell_type` | string | 是 | PanglaoDB cell-type label, e.g. 'B cells', 'Hepatocytes', 'T cells'. Case-insensitive exact match. Required. |
| `species` | string \| null | 否 | Optional species filter: 'human' (Hs) or 'mouse' (Mm). If omitted, returns markers for all species. |
| `organ` | string \| null | 否 | Optional organ filter, e.g. 'Liver', 'Immune system', 'Pancreas'. Case-insensitive exact match. |
| `canonical_only` | boolean \| null | 否 | If true, return only canonical marker genes. Default false. |
| `limit` | integer \| null | 否 | Maximum number of marker genes to return. Default 50. |

??? note "上游说明（英文原文）"

    Get curated marker genes for a single-cell type from PanglaoDB (panglaodb.se). Returns marker genes used to annotate scRNA-seq clusters, each with canonical-marker flag, ubiquitousness index, and per-gene sensitivity/specificity for human and mouse. Useful for cell-type annotation and validating cluster identity. cell_type must match a PanglaoDB cell-type label (use PanglaoDB_list_cell_types to discover them; examples: 'B cells', 'Hepatocytes', 'T cells', 'Macrophages'). Optionally filter by species ('human' or 'mouse'), organ, or canonical_only. Results are sorted canonical-first, then by descending human specificity. Note: the marker set is the frozen 2020 PanglaoDB release (no longer updated) — reliable for canonical markers, but not for cell types/markers described since.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PanglaoDB_markers_for_cell_type',
     'arguments': {'cell_type': '...'}},
)
```

## `Pharos_get_target_expression` {#Pharos_get_target_expression}

Returns GTEx per-tissue TPM expression for a target gene from Pharos.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**无
- **关键词：**gtex, tpm, tissue expression
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 否 | Gene symbol (e.g., 'DRD2', 'EGFR'). Use either gene or uniprot. |
| `uniprot` | string | 否 | UniProt accession (e.g., 'P14416'). Use either gene or uniprot. |

??? note "上游说明（英文原文）"

    Get GTEx baseline tissue RNA expression (per-tissue TPM) for a drug target from Pharos/TCRD (by gene symbol or UniProt ID). Returns one record per GTEx tissue with tissue name, TPM value, and gender stratification when available. Use to assess where a target is naturally expressed (e.g., for on-target tissue-toxicity reasoning).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Pharos_get_target_expression', 'arguments': {}},
)
```

## `ReactomeAnalysis_expression_analysis` {#ReactomeAnalysis_expression_analysis}

Runs Reactome expression analysis mapping fold-change values for identifiers onto pathways.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`identifiers`
- **关键词：**fold change, expression analysis, pathway
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `identifiers` | string | 是 | Tab-delimited expression rows, one identifier per line as 'GENE&lt;TAB&gt;VALUE'. Newline-separated rows. Gene symbols, UniProt, or Ensembl IDs in the first column; numeric fold-change/expression value(s) in following tab-separated column(s). Example: 'PTEN\\t2.5\\nTP53\\t-1.8\\nEGFR\\t3.1'. |
| `page_size` | integer \| null | 否 | Number of pathways to return (default 20, max 50). |
| `include_disease` | boolean \| null | 否 | Include disease pathways in results (default true). |
| `projection` | boolean \| null | 否 | Project identifiers to human pathways via orthology (default false; set true for non-human input). |

??? note "上游说明（英文原文）"

    Quantitative EXPRESSION analysis: map fold-change or numeric expression values onto Reactome pathways (analysis type=EXPRESSION). Unlike ReactomeAnalysis_pathway_enrichment (which only accepts a plain newline-separated identifier list and runs overrepresentation), this tool submits TAB-DELIMITED 'GENE\\tVALUE' rows so Reactome treats the input as an expression matrix and overlays each value per pathway. Each returned pathway carries 'entities_exp' (the matched expression values for that pathway) and the result reports the expression column names plus the min/max of the submitted values. Example: submitting 'PTEN\\t2.5\\nTP53\\t-1.8\\nEGFR\\t3.1' returns analysis_type=EXPRESSION with pathways such as 'Regulation of TP53 Expression' carrying entities_exp=\[-1.8\]. INTERPRETING RESULTS: rank pathways by 'fdr' or 'p_value' - those are the enrichment statistics. 'entities_coverage' (= entities_found / entities_total) is the fraction of a pathway your identifiers cover. 'pathway_size_fraction_of_reactome' (deprecated alias 'entities_ratio') is Reactome's raw 'ratio': pathway size over Reactome's entity universe for the species, independent of the identifiers you submit - never rank by it.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ReactomeAnalysis_expression_analysis',
     'arguments': {'identifiers': '...'}},
)
```

## `RNAseq_edger_limma_de` {#RNAseq_edger_limma_de}

Runs edgeR QL-F or limma-voom differential expression on a count matrix and metadata.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`counts_file`, `metadata_file`, `contrast`, `method`
- **关键词：**edger, limma-voom, differential expression, rna-seq
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `counts_file` | string | 是 | Path to counts CSV (genes x samples; first column = gene IDs) |
| `metadata_file` | string | 是 | Path to sample metadata CSV (one row per sample) |
| `design` | string \| null | 否 | Design formula, e.g. '~ condition' or '~ batch + condition' (default '~ condition') |
| `contrast` | string | 是 | Contrast as 'factor,level1,level2' (tests level1 vs level2) |
| `method` | string | 是 | DE framework: 'edger' (QL-F) or 'limma' (limma-voom) 可选值：`edger`, `limma` |
| `lfc_threshold` | number \| null | 否 | abs(logFC) threshold for the sig sets (default 0.5) |
| `fdr_threshold` | number \| null | 否 | FDR (adjusted p) threshold (default 0.05) |
| `expr_threshold` | number \| null | 否 | Expression threshold (logCPM/AveExpr) for the strict set (default 0) |
| `exclude_samples` | string \| null | 否 | Comma-separated sample IDs to drop before analysis |
| `report_genes` | string \| null | 否 | Comma-separated gene IDs to report individual logFC/FDR for |

??? note "上游说明（英文原文）"

    Bulk RNA-seq differential expression via edgeR (QL-F) or limma-voom on a count matrix + sample metadata — the two standard DE frameworks alongside DESeq2 (run_deseq2_analysis). method='edger' (good for small replicate counts) or 'limma' (large n / complex designs). Returns DEG counts at FDR-only, FDR+\|logFC\|, and strict (FDR+\|logFC\|+expression) filters, optional per-gene logFC/FDR, and the full ranked table on disk. Requires Rscript with Bioconductor edgeR + limma.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'RNAseq_edger_limma_de',
     'arguments': {'counts_file': '...',
                   'metadata_file': '...',
                   'contrast': '...',
                   'method': 'edger'}},
)
```

## `run_deseq2_analysis` {#run_deseq2_analysis}

Runs R DESeq2 differential expression on a count matrix and metadata with optional enrichment.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**无
- **关键词：**deseq2, differential expression, rna-seq, count matrix, log fold change
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | deseq2: differential expression, enrichgo: GO enrichment with clusterProfiler::simplify 可选值：`deseq2`, `enrichgo` 默认 `deseq2` |
| `counts_file` | string | 否 | Path to CSV count matrix (genes x samples, first column = gene IDs) |
| `metadata_file` | string | 否 | Path to CSV metadata (samples x variables, must include design variables) |
| `design` | string | 否 | R formula string (e.g., '~ condition', '~ Replicate + Media + Strain') |
| `contrast` | string \| null | 否 | Contrast specification: 'variable, level, reference' (e.g., 'Strain, 97, 1') |
| `ref_level` | string \| null | 否 | Reference level: 'variable, level' (e.g., 'condition, Control') |
| `alpha` | number | 否 | Significance threshold (default: 0.05) |
| `lfc_threshold` | number | 否 | Log2 fold change threshold for filtering (default: 0 = no filter) |
| `lfc_shrinkage` | boolean | 否 | Apply apeglm LFC shrinkage (default: false) |
| `refit_cooks` | boolean | 否 | Refit Cook's outliers (default: false) |
| `gene_list_file` | string \| null | 否 | For enrichgo: path to file with one gene ID per line |
| `background_file` | string \| null | 否 | For enrichgo: path to background gene list file |
| `ontology` | string | 否 | For enrichgo: GO ontology (BP, CC, MF). Default: BP |
| `simplify_cutoff` | number | 否 | For enrichgo: simplify similarity cutoff (default: 0.7) |
| `id_type` | string | 否 | Gene ID type (ENSEMBL, ENTREZID, SYMBOL). Default: ENSEMBL |
| `organism` | string | 否 | Organism annotation DB (org.Hs.eg.db, org.Mm.eg.db). Default: org.Hs.eg.db |

??? note "上游说明（英文原文）"

    Run R DESeq2 differential expression analysis on a count matrix + metadata. Returns gene counts, significance stats, and dispersion estimates. Uses actual R DESeq2 (not pydeseq2) for maximum accuracy. Supports custom design formulas, contrasts, reference levels, LFC shrinkage, and refit_cooks.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'run_deseq2_analysis', 'arguments': {}},
)
```

## `SCXA_get_cluster_marker_genes` {#SCXA_get_cluster_marker_genes}

Gets marker genes per cell cluster or cell type for an EBI Single Cell Expression Atlas experiment, with ranking statistics.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`experiment_accession`
- **关键词：**marker genes, single-cell rna-seq, cell cluster, cell type, expression atlas, scrna-seq
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `experiment_accession` | string | 是 | SCXA experiment accession (e.g. 'E-MTAB-5061', 'E-ENAD-13', 'E-CURD-4'). Use SCXA_list_experiments to find accessions. |
| `marker_type` | string | 否 | 'clusters' for marker genes per cluster at resolution k (default); 'cell_types' for marker genes per inferred cell type (requires organism_part). 可选值：`clusters`, `cell_types` 默认 `clusters` |
| `k` | integer | 否 | Clustering resolution / number of clusters for marker_type='clusters' (e.g. 8, 4, 6). Required for clusters mode. |
| `organism_part` | string | 否 | Organism part for marker_type='cell_types' (e.g. 'pancreas', 'lung'). Required for cell_types mode. |
| `limit` | integer | 否 | Optional maximum number of marker rows to return. |

??? note "上游说明（英文原文）"

    Get EBI Single Cell Expression Atlas computed marker genes per cell cluster (or per cell type) for an experiment - results that SCXA_list_experiments and SCXA_search_gene cannot retrieve. marker_type='clusters' (default): top marker genes per cell cluster at clustering resolution k (requires `k`); each row has geneName, the cell-group it marks, expression value, p-value, and unit (CPM). marker_type='cell_types': marker genes for inferred cell types (requires `organism_part`). Example: E-MTAB-5061 k=8 -&gt; rows like \{gene_name: ENSG00000224981, cell_group_value_where_marker: '1', p_value: 4.37e-71, expression_unit: CPM\}. E-ENAD-13 k=4 -&gt; Olfm4 p_value 3.0e-179.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SCXA_get_cluster_marker_genes',
     'arguments': {'experiment_accession': '...'}},
)
```

## `SCXA_search_gene` {#SCXA_search_gene}

Searches the EBI Single Cell Expression Atlas for scRNA-seq experiments in which a given gene symbol or Ensembl ID is expressed.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`gene`
- **关键词：**gene expression, single-cell, ensembl id, scrna-seq, expression atlas
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 是 | Gene symbol (e.g., 'TP53', 'BRCA1', 'INS') or Ensembl gene ID (e.g., 'ENSG00000141510'). |
| `species` | string | 否 | Species name to filter (e.g., 'Homo sapiens', 'Mus musculus'). Optional. |

??? note "上游说明（英文原文）"

    Search the EBI Single Cell Expression Atlas for single-cell experiments where a specific gene is expressed. Accepts gene symbols (TP53, BRCA1) or Ensembl IDs (ENSG00000141510). Returns matching scRNA-seq experiments with cell counts, technology, PubMed IDs, and DOIs. Useful for finding single-cell datasets studying a gene of interest. Example: TP53 returns 15 experiments including ovarian cancer, lung atlas datasets.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SCXA_search_gene', 'arguments': {'gene': '...'}},
)
```

## `ssGSEA_score` {#ssGSEA_score}

Computes single-sample GSEA scores for each gene set in each sample from a genes-by-samples expression matrix.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`expression`
- **关键词：**gene set enrichment, signature score, expression matrix, gsea, gene sets
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `expression` | object | 是 | \{gene: \[value_per_sample, ...\]\} — a genes x samples expression matrix (e.g. log-normalized). |
| `samples` | array \| null | 否 | Sample names (same order as the per-gene value lists); defaults to sample_1.. |
| `gene_sets` | object \| null | 否 | \{set_name: \[gene, ...\]\} collection to score. Alternative to gene_set. |
| `gene_set` | array \| null | 否 | A single gene set (list of gene symbols). |
| `alpha` | number \| null | 否 | Rank-weighting exponent (default 0.25, the ssGSEA standard). |

??? note "上游说明（英文原文）"

    Single-sample GSEA (Barbie 2009): score each sample for each gene set from an expression matrix (genes x samples), producing a gene_set x sample signature-score matrix. Positive = the set's genes are among that sample's highly-expressed genes. Basis for per-sample signature scoring (immune infiltration, pathway activity). Pure-compute (no R); supply gene sets from MSigDB/Enrichr.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ssGSEA_score', 'arguments': {'expression': {}}},
)
```

## `SYNERGxDB_get_biomarker_association` {#SYNERGxDB_get_biomarker_association}

Correlates per-cell-line RNA-seq expression of a gene with drug-combination synergy scores in SYNERGxDB.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**无
- **关键词：**pharmacogenomics, biomarker, fpkm, drug synergy, cancer cell lines
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`get_biomarker_association` 默认 `get_biomarker_association` |
| `gene` | string | 否 | HGNC gene symbol to correlate against synergy scores (e.g., 'EGFR', 'BRAF', 'TP53'). Alias: gene_name. |
| `gene_name` | string | 否 | Alias for gene. HGNC gene symbol. |

??? note "上游说明（英文原文）"

    Pharmacogenomic biomarker discovery for drug combinations in SYNERGxDB: correlates per-cell-line gene expression (rnaseq FPKM) against the measured drug-combination synergy scores (Bliss, Loewe, HSA, ZIP) for a queried gene. Answers 'does expression of gene X track with combination synergy across cancer cell lines?'. Returns one row per cell line with that cell line's FPKM for the gene plus its synergy metrics. Pass the HGNC gene symbol via 'gene' (e.g., 'EGFR', 'BRAF'). Example: gene='EGFR' returns rows like \{fpkm:0, cellName:'A2058', bliss:-1.365, loewe:-0.259, hsa:-0.628, zip:-0.498\}. An unknown gene returns an empty list with a message.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SYNERGxDB_get_biomarker_association', 'arguments': {}},
)
```

## `TIMER2_gene_correlation` {#TIMER2_gene_correlation}

Computes Spearman correlation between expression of two genes across TCGA cancer samples using cBioPortal data.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`cancer`, `gene1`, `gene2`
- **关键词：**co-expression, spearman, tcga, cbioportal, tumor
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`gene_correlation` 默认 `gene_correlation` |
| `cancer` | string | 是 | TCGA cancer type (e.g., 'BRCA', 'LUAD') |
| `gene1` | string | 是 | First gene symbol (e.g., 'CD8A') |
| `gene2` | string | 是 | Second gene symbol (e.g., 'PDCD1') |

??? note "上游说明（英文原文）"

    Analyze Spearman correlation between two genes across TCGA cancer samples using cBioPortal data. Returns correlation coefficient and p-value. Useful for identifying co-expressed genes in tumor microenvironment. Note: TIMER2.0/3.0 REST API is no longer public; uses cBioPortal as data source. Example: Correlate CD8A and PDCD1 expression in BRCA.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'TIMER2_gene_correlation',
     'arguments': {'cancer': '...', 'gene1': '...', 'gene2': '...'}},
)
```

## `TIMER2_immune_estimation` {#TIMER2_immune_estimation}

Estimates immune cell infiltration proxies in TCGA samples from marker gene expression, optionally correlating with a gene.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`cancer`
- **关键词：**immune infiltration, marker genes, tcga, tumor microenvironment, cbioportal
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`immune_estimation` 默认 `immune_estimation` |
| `cancer` | string | 是 | TCGA cancer type abbreviation (e.g., 'BRCA', 'LUAD', 'COAD', 'SKCM', 'GBM') |
| `gene` | string \| null | 否 | Optional gene symbol to correlate with immune infiltration (e.g., 'CD8A', 'PDCD1') |

??? note "上游说明（英文原文）"

    Estimate immune cell abundance in TCGA cancer samples using TIMER-equivalent marker gene expression from cBioPortal. Returns infiltration proxies for 6 immune cell types (B cells, CD4+ T cells, CD8+ T cells, neutrophils, macrophages, dendritic cells) using canonical marker genes (CD19, CD4, CD8A, FCGR3B, CD68, ITGAX). Note: TIMER2.0/3.0 REST API is no longer public; this tool uses cBioPortal TCGA data as a proxy.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'TIMER2_immune_estimation',
     'arguments': {'cancer': '...'}},
)
```

## `WormBase_get_expression` {#WormBase_get_expression}

Gets expression data for a C. elegans gene from WormBase: tissues, developmental stages, localization and expression clusters.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`gene_id`
- **关键词：**c. elegans, nematode, tissue expression, developmental stage, gfp reporter
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_id` | string | 是 | WormBase gene identifier. Examples: 'WBGene00006763' (unc-26), 'WBGene00000898' (daf-2). |

??? note "上游说明（英文原文）"

    Get expression data for a C. elegans gene from WormBase. Returns tissues where the gene is expressed, developmental stages, subcellular localization, and expression clusters. Expression data is curated from GFP reporters, antibody staining, RNA in situ, and microarray experiments. Example: WBGene00006763 (unc-26) is expressed in 38 tissues including neurons, pharynx, and body wall muscle.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'WormBase_get_expression',
     'arguments': {'gene_id': '...'}},
)
```

## `ZFIN_get_gene_expression` {#ZFIN_get_gene_expression}

Gets zebrafish gene expression annotations (developmental stage and anatomical location) from the Alliance API.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`gene_id`
- **关键词：**zebrafish, developmental stage, anatomical expression, in situ
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_id` | string | 是 | ZFIN gene ID, with or without the 'ZFIN:' prefix. Examples: 'ZFIN:ZDB-GENE-990415-8' (pax2a), 'ZFIN:ZDB-GENE-980526-166' (shha), 'ZFIN:ZDB-GENE-990415-72' (tp53). |
| `limit` | integer | 否 | Max annotations to return (default 50, max 100). 默认 `50` |

??? note "上游说明（英文原文）"

    Get expression annotations for a zebrafish gene from the Alliance of Genome Resources. Returns per-annotation records of developmental stage (whenExpressed) and anatomical location (whereExpressed) for the gene, plus the total annotation count. Example: 'ZFIN:ZDB-GENE-980526-166' (shha) has hundreds of annotations across developmental stages. Zebrafish is a premier model for vertebrate development.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ZFIN_get_gene_expression',
     'arguments': {'gene_id': '...'}},
)
```
