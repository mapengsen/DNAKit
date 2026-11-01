# GWAS MCP (muslus) (`gwas-mcp`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **33** 个 DNA/基因相关工具（服务器共列出 44 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [muslus/gwas-mcp @ `5d8160c`](https://github.com/muslus/gwas-mcp/tree/5d8160cda988e959234467cefae554ebdc1876e5) |
| 许可证 | MIT |
| 启动命令 | `python -m gwas_mcp.server` |
| 传输方式 | stdio |
| 安装 | Prepare the pinned muslus/gwas-mcp source and its dependencies in an isolated Python &gt;=3.10 environment, then point executable_path at that interpreter. Avoid assuming that the separately published gwas-mcp distribution has identical code. Pin 'mcp&lt;2': under mcp 2.x the server fails at start ('Server' has no list_tools). |
| 副作用 | Queries genomic evidence services; selected QC/association tools may execute local analyses and write files. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`annotate_snps`](#annotate_snps) | `variant` | Annotates rsIDs with gene names, functional consequences (Ensembl VEP) and gnomAD population allele frequencies. |
| [`calculate_genomic_inflation`](#calculate_genomic_inflation) | `association` | Calculates the genomic inflation factor (lambda GC) from GWAS summary statistics p-values. |
| [`calculate_heritability_ldsc`](#calculate_heritability_ldsc) | `association` | Estimates SNP-heritability from GWAS summary statistics using LD Score regression. |
| [`calculate_missingness`](#calculate_missingness) | `variant` | Calculates per-individual and per-variant genotype missingness rates from PLINK files. |
| [`calculate_pca`](#calculate_pca) | `breeding_population` | Performs principal component analysis on PLINK genotype data to assess population stratification. |
| [`calculate_prs`](#calculate_prs) | `variant` | Calculates polygenic risk scores for individuals from GWAS summary statistics and target PLINK genotypes with clumping. |
| [`clump_snps`](#clump_snps) | `association` | Performs LD-based clumping of GWAS summary statistics to identify independent significant loci. |
| [`create_manhattan_plot`](#create_manhattan_plot) | `association` | Creates a Manhattan plot of -log10(p) across chromosomes from GWAS summary statistics. |
| [`create_pca_plot`](#create_pca_plot) | `breeding_population` | Creates a PC1 vs PC2 scatter plot of genotype PCA results with optional population labels. |
| [`create_qq_plot`](#create_qq_plot) | `association` | Creates a quantile-quantile plot of observed vs expected GWAS p-values. |
| [`create_regional_plot`](#create_regional_plot) | `association` | Creates a LocusZoom-style regional association plot around a lead SNP from GWAS summary statistics. |
| [`extract_locus`](#extract_locus) | `variant` | Extracts PLINK genotype data for a genomic window around a lead SNP for fine-mapping. |
| [`gene_set_enrichment`](#gene_set_enrichment) | `genome_annotation` | Performs gene set enrichment analysis of a gene list against GO and KEGG libraries via Enrichr. |
| [`get_clinvar_variant`](#get_clinvar_variant) | `variant` | Gets detailed ClinVar information (clinical significance, conditions) for a variant by rsID or ClinVar ID. |
| [`get_disease_associations`](#get_disease_associations) | `association` | Gets diseases associated with a gene from Open Targets with evidence scores. |
| [`get_drug_gene_interactions`](#get_drug_gene_interactions) | `variant` | Gets PharmGKB drug-gene interaction and pharmacogenomic annotations for a gene. |
| [`get_eqtl_data`](#get_eqtl_data) | `association` | Gets expression quantitative trait loci (eQTL) associations from GTEx for an rsID and tissue, optionally filtered by gene. |
| [`get_functional_enrichment`](#get_functional_enrichment) | `genome_annotation` | Runs STRING functional enrichment for a list of protein/gene names, returning enriched GO terms, pathways and domains. |
| [`get_gene_diseases`](#get_gene_diseases) | `association` | Gets genetic diseases associated with a gene from OMIM and other sources. |
| [`get_gene_pathways`](#get_gene_pathways) | `genome_annotation` | Finds all KEGG pathways that a gene participates in. |
| [`get_kegg_pathway`](#get_kegg_pathway) | `genome_annotation` | Gets detailed KEGG pathway information, including its member genes and description, by pathway ID. |
| [`get_variant_info`](#get_variant_info) | `variant` | Gets Ensembl variant details for an rsID: position, alleles, clinical significance and consequences. |
| [`identify_significant_snps`](#identify_significant_snps) | `association` | Filters GWAS summary statistics to genome-wide significant and suggestive SNPs by p-value thresholds. |
| [`ld_pruning`](#ld_pruning) | `association` | Prunes SNPs in high linkage disequilibrium from PLINK data to create an independent variant set. |
| [`query_gwas_catalog`](#query_gwas_catalog) | `association` | Queries the NHGRI-EBI GWAS Catalog for reported associations by rsID, trait or gene. |
| [`run_gwas`](#run_gwas) | `association` | Performs a GWAS with linear or logistic regression on PLINK genotypes and a phenotype file, with optional covariates. |
| [`search_clinvar`](#search_clinvar) | `variant` | Searches ClinVar for variant clinical interpretations by gene, rsID or condition. |
| [`search_ensembl_gene`](#search_ensembl_gene) | `genome_annotation` | Searches Ensembl for a gene by symbol or Ensembl ID, returning location, biotype and description. |
| [`search_kegg_pathway`](#search_kegg_pathway) | `genome_annotation` | Searches KEGG for metabolic and signaling pathways by pathway name or gene symbol. |
| [`search_ncbi_gene`](#search_ncbi_gene) | `genome_annotation` | Searches the NCBI Gene database by gene symbol, name or Gene ID. |
| [`search_omim`](#search_omim) | `association` | Searches OMIM for Mendelian genetic diseases and phenotypes by gene symbol or disease name. |
| [`search_pharmgkb`](#search_pharmgkb) | `variant` | Searches PharmGKB for pharmacogenomic drug-gene-variant data by gene, drug or rsID. |
| [`vcf_quality_check`](#vcf_quality_check) | `variant` | Performs quality control on a VCF file, filtering variants by MAF, missingness and Hardy-Weinberg equilibrium. |

## `annotate_snps` {#annotate_snps}

Annotates rsIDs with gene names, functional consequences (Ensembl VEP) and gnomAD population allele frequencies.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`rsids`
- **关键词：**rsid, vep, gnomad, allele frequency, consequence

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `rsids` | array of string | 是 | List of rsIDs to annotate (e.g., \['rs12345', 'rs67890'\]) |
| `include_frequencies` | boolean | 否 | Include population allele frequencies (default: true) 默认 `True` |

??? note "上游说明（英文原文）"

    Annotate SNPs with gene names, functional consequences, and allele frequencies from Ensembl VEP and gnomAD.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'annotate_snps',
    {'rsids': []},
)
```

## `calculate_genomic_inflation` {#calculate_genomic_inflation}

Calculates the genomic inflation factor (lambda GC) from GWAS summary statistics p-values.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`sumstats_path`
- **关键词：**lambda gc, gwas, summary statistics, population stratification

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sumstats_path` | string | 是 | Path to GWAS summary statistics file |
| `pvalue_column` | string | 否 | Name of p-value column (default: P) 默认 `P` |

??? note "上游说明（英文原文）"

    Calculate genomic inflation factor (lambda GC) from GWAS summary statistics. Used to assess population stratification.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'calculate_genomic_inflation',
    {'sumstats_path': '...'},
)
```

## `calculate_heritability_ldsc` {#calculate_heritability_ldsc}

Estimates SNP-heritability from GWAS summary statistics using LD Score regression.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`sumstats_path`
- **关键词：**ldsc, heritability, ld score, gwas, summary statistics

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sumstats_path` | string | 是 | Path to GWAS summary statistics |
| `ld_scores_path` | string | 否 | Path to LD scores file |
| `sample_size` | integer | 否 | GWAS sample size (required if not in sumstats) |

??? note "上游说明（英文原文）"

    Estimate SNP-heritability using LD Score regression. Requires GWAS summary statistics and LD scores.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'calculate_heritability_ldsc',
    {'sumstats_path': '...'},
)
```

## `calculate_missingness` {#calculate_missingness}

Calculates per-individual and per-variant genotype missingness rates from PLINK files.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`plink_prefix`
- **关键词：**plink, genotype qc, call rate, missing data

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `plink_prefix` | string | 是 | Path prefix for PLINK files |
| `output_path` | string | 否 | Optional path to save missingness report |

??? note "上游说明（英文原文）"

    Calculate per-individual and per-variant missingness rates from genotype data.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'calculate_missingness',
    {'plink_prefix': '...'},
)
```

## `calculate_pca` {#calculate_pca}

Performs principal component analysis on PLINK genotype data to assess population stratification.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**`plink_prefix`
- **关键词：**plink, population structure, principal components, stratification

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `plink_prefix` | string | 是 | Path prefix for PLINK files (without .bed/.bim/.fam extension) |
| `n_components` | integer | 否 | Number of principal components to compute (default: 10) 默认 `10` |
| `output_path` | string | 否 | Optional path to save PCA results |

??? note "上游说明（英文原文）"

    Perform Principal Component Analysis for population stratification. Uses PLINK files as input.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'calculate_pca',
    {'plink_prefix': '...'},
)
```

## `calculate_prs` {#calculate_prs}

Calculates polygenic risk scores for individuals from GWAS summary statistics and target PLINK genotypes with clumping.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`sumstats_path`, `target_plink_prefix`
- **关键词：**polygenic risk score, prs, plink, clumping, p-value threshold

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sumstats_path` | string | 是 | Path to GWAS summary statistics (base data) |
| `target_plink_prefix` | string | 是 | Path prefix for target PLINK files |
| `pvalue_threshold` | number | 否 | P-value threshold for SNP inclusion (default: 1.0 - all SNPs) 默认 `1.0` |
| `clump_r2` | number | 否 | LD r² threshold for clumping (default: 0.1) 默认 `0.1` |
| `output_path` | string | 否 | Path to save PRS results |

??? note "上游说明（英文原文）"

    Calculate Polygenic Risk Scores (PRS) for individuals using GWAS summary statistics and target genotypes.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'calculate_prs',
    {'sumstats_path': '...', 'target_plink_prefix': '...'},
)
```

## `clump_snps` {#clump_snps}

Performs LD-based clumping of GWAS summary statistics to identify independent significant loci.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`sumstats_path`
- **关键词：**ld clumping, plink, index snp, gwas, independent signals

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sumstats_path` | string | 是 | Path to GWAS summary statistics |
| `plink_prefix` | string | 否 | Path prefix for PLINK files (LD reference) |
| `p1` | number | 否 | P-value threshold for index SNPs (default: 5e-8) 默认 `5e-08` |
| `p2` | number | 否 | P-value threshold for clumped SNPs (default: 1e-4) 默认 `0.0001` |
| `r2` | number | 否 | LD r² threshold (default: 0.5) 默认 `0.5` |
| `kb` | integer | 否 | Clumping window in kb (default: 250) 默认 `250` |
| `output_path` | string | 否 | Path to save clumped results |

??? note "上游说明（英文原文）"

    Perform LD-based clumping of GWAS results to identify independent significant signals.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'clump_snps',
    {'sumstats_path': '...'},
)
```

## `create_manhattan_plot` {#create_manhattan_plot}

Creates a Manhattan plot of -log10(p) across chromosomes from GWAS summary statistics.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`sumstats_path`
- **关键词：**manhattan plot, gwas, p-value, visualization

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sumstats_path` | string | 是 | Path to GWAS summary statistics file |
| `significance_threshold` | number | 否 | Genome-wide significance threshold (default: 5e-8) 默认 `5e-08` |
| `suggestive_threshold` | number | 否 | Suggestive significance threshold (default: 1e-5) 默认 `1e-05` |
| `output_path` | string | 否 | Path to save the plot (PNG format) |
| `title` | string | 否 | Plot title 默认 `Manhattan Plot` |

??? note "上游说明（英文原文）"

    Create a Manhattan plot from GWAS summary statistics. Shows -log10(p-value) across chromosomes.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'create_manhattan_plot',
    {'sumstats_path': '...'},
)
```

## `create_pca_plot` {#create_pca_plot}

Creates a PC1 vs PC2 scatter plot of genotype PCA results with optional population labels.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**`pca_results_path`
- **关键词：**pca plot, population structure, visualization, ancestry

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `pca_results_path` | string | 是 | Path to PCA results file (with PC1, PC2 columns) |
| `population_labels_path` | string | 否 | Optional path to population labels file |
| `output_path` | string | 否 | Path to save the plot |
| `title` | string | 否 | Plot title 默认 `PCA Plot` |

??? note "上游说明（英文原文）"

    Create a PCA scatter plot showing population structure. Plots PC1 vs PC2.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'create_pca_plot',
    {'pca_results_path': '...'},
)
```

## `create_qq_plot` {#create_qq_plot}

Creates a quantile-quantile plot of observed vs expected GWAS p-values.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`sumstats_path`
- **关键词：**qq plot, gwas, p-value, inflation, visualization

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sumstats_path` | string | 是 | Path to GWAS summary statistics file |
| `output_path` | string | 否 | Path to save the plot |
| `title` | string | 否 | Plot title 默认 `QQ Plot` |

??? note "上游说明（英文原文）"

    Create a Quantile-Quantile (QQ) plot from GWAS p-values. Shows observed vs expected p-values.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'create_qq_plot',
    {'sumstats_path': '...'},
)
```

## `create_regional_plot` {#create_regional_plot}

Creates a LocusZoom-style regional association plot around a lead SNP from GWAS summary statistics.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`sumstats_path`, `lead_snp`
- **关键词：**locuszoom, regional plot, lead snp, gwas

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sumstats_path` | string | 是 | Path to GWAS summary statistics |
| `lead_snp` | string | 是 | Lead SNP rsID or CHR:POS |
| `window_kb` | integer | 否 | Window size in kilobases (default: 500) 默认 `500` |
| `output_path` | string | 否 | Path to save the plot |

??? note "上游说明（英文原文）"

    Create a regional association plot (LocusZoom-style) for a specific genomic region around a lead SNP.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'create_regional_plot',
    {'sumstats_path': '...', 'lead_snp': '...'},
)
```

## `extract_locus` {#extract_locus}

Extracts PLINK genotype data for a genomic window around a lead SNP for fine-mapping.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`plink_prefix`, `lead_snp`
- **关键词：**plink, genotype, fine-mapping, lead snp, region

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `plink_prefix` | string | 是 | Path prefix for PLINK files |
| `lead_snp` | string | 是 | Lead SNP identifier (rsID or CHR:POS) |
| `window_kb` | integer | 否 | Window size in kb (default: 500) 默认 `500` |
| `output_prefix` | string | 否 | Output prefix for extracted files |

??? note "上游说明（英文原文）"

    Extract genotype data for a specific genomic region around a lead SNP for fine-mapping.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'extract_locus',
    {'plink_prefix': '...', 'lead_snp': '...'},
)
```

## `gene_set_enrichment` {#gene_set_enrichment}

Performs gene set enrichment analysis of a gene list against GO and KEGG libraries via Enrichr.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`genes`
- **关键词：**enrichr, enrichment, gene set, kegg, go

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `genes` | array of string | 是 | List of gene symbols (e.g., \['BRCA1', 'TP53', 'EGFR'\]) |
| `databases` | array of string | 否 | Databases to query (default: GO_Biological_Process, KEGG) 默认 `\['GO_Biological_Process_2021', 'KEGG_2021_Human'\]` |

??? note "上游说明（英文原文）"

    Perform gene set enrichment analysis using GO and KEGG pathways via Enrichr.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'gene_set_enrichment',
    {'genes': []},
)
```

## `get_clinvar_variant` {#get_clinvar_variant}

Gets detailed ClinVar information (clinical significance, conditions) for a variant by rsID or ClinVar ID.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant_id`
- **关键词：**clinvar, pathogenicity, rsid, clinical significance

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_id` | string | 是 | Variant identifier - rsID (e.g., 'rs80357906') or ClinVar ID |

??? note "上游说明（英文原文）"

    Get detailed ClinVar information for a specific variant by rsID or ClinVar ID.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'get_clinvar_variant',
    {'variant_id': '...'},
)
```

## `get_disease_associations` {#get_disease_associations}

Gets diseases associated with a gene from Open Targets with evidence scores.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`gene`
- **关键词：**open targets, gene-disease, association score, evidence

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 是 | Gene symbol |
| `limit` | integer | 否 | Maximum number of diseases (default: 20) 默认 `20` |

??? note "上游说明（英文原文）"

    Get disease associations for a gene from Open Targets with evidence scores.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'get_disease_associations',
    {'gene': '...'},
)
```

## `get_drug_gene_interactions` {#get_drug_gene_interactions}

Gets PharmGKB drug-gene interaction and pharmacogenomic annotations for a gene.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`gene`
- **关键词：**pharmgkb, pharmacogenomics, drug response, cyp2d6

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 是 | Gene symbol (e.g., 'CYP2D6', 'VKORC1') |

??? note "上游说明（英文原文）"

    Get drug-gene interaction annotations from PharmGKB for a specific gene.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'get_drug_gene_interactions',
    {'gene': '...'},
)
```

## `get_eqtl_data` {#get_eqtl_data}

Gets expression quantitative trait loci (eQTL) associations from GTEx for an rsID and tissue, optionally filtered by gene.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`rsid`
- **关键词：**gtex, eqtl, tissue, rsid, gene expression

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `rsid` | string | 是 | rsID to query |
| `tissue` | string | 否 | GTEx tissue type (e.g., 'Whole_Blood', 'Brain_Cortex') |
| `gene` | string | 否 | Optional gene to filter results |

??? note "上游说明（英文原文）"

    Get expression quantitative trait loci (eQTL) data from GTEx for a SNP and tissue.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'get_eqtl_data',
    {'rsid': '...'},
)
```

## `get_functional_enrichment` {#get_functional_enrichment}

Runs STRING functional enrichment for a list of protein/gene names, returning enriched GO terms, pathways and domains.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`proteins`
- **关键词：**string, enrichment, go, pathways, gene list

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `proteins` | array of string | 是 | List of protein/gene names |
| `species` | integer | 否 | NCBI species taxonomy ID (default: 9606 for human) 默认 `9606` |

??? note "上游说明（英文原文）"

    Get functional enrichment analysis for a list of proteins using STRING. Returns enriched GO terms, pathways, and domains.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'get_functional_enrichment',
    {'proteins': []},
)
```

## `get_gene_diseases` {#get_gene_diseases}

Gets genetic diseases associated with a gene from OMIM and other sources.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`gene`
- **关键词：**omim, mendelian disease, gene-disease, inheritance

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 是 | Gene symbol (e.g., 'BRCA1', 'CFTR') |

??? note "上游说明（英文原文）"

    Get genetic diseases associated with a gene from OMIM and other sources.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'get_gene_diseases',
    {'gene': '...'},
)
```

## `get_gene_pathways` {#get_gene_pathways}

Finds all KEGG pathways that a gene participates in.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`gene`
- **关键词：**kegg, pathways, gene function

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 是 | Gene symbol (e.g., 'TP53', 'BRCA1') |

??? note "上游说明（英文原文）"

    Find all KEGG pathways that a gene participates in.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'get_gene_pathways',
    {'gene': '...'},
)
```

## `get_kegg_pathway` {#get_kegg_pathway}

Gets detailed KEGG pathway information, including its member genes and description, by pathway ID.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`pathway_id`
- **关键词：**kegg, pathway genes, pathway id, description

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `pathway_id` | string | 是 | KEGG pathway ID (e.g., 'hsa04110' for cell cycle) |

??? note "上游说明（英文原文）"

    Get detailed KEGG pathway information including genes and description.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'get_kegg_pathway',
    {'pathway_id': '...'},
)
```

## `get_variant_info` {#get_variant_info}

Gets Ensembl variant details for an rsID: position, alleles, clinical significance and consequences.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`rsid`
- **关键词：**ensembl, rsid, alleles, consequence, dbsnp

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `rsid` | string | 是 | dbSNP rsID (e.g., 'rs1234567') |

??? note "上游说明（英文原文）"

    Get detailed variant/SNP information from Ensembl by rsID. Returns position, alleles, clinical significance, and consequences.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'get_variant_info',
    {'rsid': '...'},
)
```

## `identify_significant_snps` {#identify_significant_snps}

Filters GWAS summary statistics to genome-wide significant and suggestive SNPs by p-value thresholds.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`sumstats_path`
- **关键词：**genome-wide significance, gwas, p-value, 5e-8, summary statistics

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sumstats_path` | string | 是 | Path to GWAS summary statistics file |
| `pvalue_threshold` | number | 否 | P-value threshold for significance (default: 5e-8) 默认 `5e-08` |
| `suggestive_threshold` | number | 否 | Suggestive significance threshold (default: 1e-5) 默认 `1e-05` |
| `output_path` | string | 否 | Path to save significant SNPs |

??? note "上游说明（英文原文）"

    Filter GWAS results to identify genome-wide significant SNPs based on p-value threshold.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'identify_significant_snps',
    {'sumstats_path': '...'},
)
```

## `ld_pruning` {#ld_pruning}

Prunes SNPs in high linkage disequilibrium from PLINK data to create an independent variant set.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`plink_prefix`
- **关键词：**linkage disequilibrium, plink, r2, independent snps

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `plink_prefix` | string | 是 | Path prefix for PLINK files |
| `r2_threshold` | number | 否 | r² threshold for LD pruning (default: 0.2) 默认 `0.2` |
| `window_size` | integer | 否 | Window size in kb for LD calculation (default: 500) 默认 `500` |
| `step_size` | integer | 否 | Step size in variants (default: 50) 默认 `50` |

??? note "上游说明（英文原文）"

    Prune SNPs based on linkage disequilibrium (LD). Removes variants in high LD to create an independent set.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'ld_pruning',
    {'plink_prefix': '...'},
)
```

## `query_gwas_catalog` {#query_gwas_catalog}

Queries the NHGRI-EBI GWAS Catalog for reported associations by rsID, trait or gene.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**gwas catalog, nhgri-ebi, trait, rsid, reported associations

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `rsid` | string | 否 | rsID to query (e.g., 'rs12345') |
| `trait` | string | 否 | Trait/disease to search for (e.g., 'diabetes') |
| `gene` | string | 否 | Gene name to search for (e.g., 'BRCA1') |

??? note "上游说明（英文原文）"

    Query the NHGRI-EBI GWAS Catalog for previously reported associations for a SNP or trait.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'query_gwas_catalog',
    {},
)
```

## `run_gwas` {#run_gwas}

Performs a GWAS with linear or logistic regression on PLINK genotypes and a phenotype file, with optional covariates.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`plink_prefix`, `phenotype_path`
- **关键词：**gwas, plink, linear regression, logistic regression, summary statistics

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `plink_prefix` | string | 是 | Path prefix for PLINK genotype files |
| `phenotype_path` | string | 是 | Path to phenotype file (tab-separated with FID, IID, PHENO columns) |
| `covariate_path` | string | 否 | Optional path to covariate file |
| `model` | string | 否 | Regression model type (default: linear) 可选值：`linear`, `logistic` 默认 `linear` |
| `output_path` | string | 否 | Path to save GWAS results |

??? note "上游说明（英文原文）"

    Perform genome-wide association study using linear or logistic regression. Returns summary statistics including p-values, beta coefficients, and standard errors.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'run_gwas',
    {'plink_prefix': '...', 'phenotype_path': '...'},
)
```

## `search_clinvar` {#search_clinvar}

Searches ClinVar for variant clinical interpretations by gene, rsID or condition.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`query`
- **关键词：**clinvar, pathogenic, benign, rsid, condition

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Search query - can be gene name (e.g., 'BRCA1'), rsID (e.g., 'rs80357906'), or condition (e.g., 'breast cancer') |
| `limit` | integer | 否 | Maximum number of results (default: 10) 默认 `10` |

??? note "上游说明（英文原文）"

    Search ClinVar for clinical variant interpretations. Find pathogenic/benign classifications for genetic variants.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'search_clinvar',
    {'query': '...'},
)
```

## `search_ensembl_gene` {#search_ensembl_gene}

Searches Ensembl for a gene by symbol or Ensembl ID, returning location, biotype and description.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`gene`
- **关键词：**ensembl, gene location, biotype, gene symbol

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 是 | Gene symbol (e.g., 'BRCA1') or Ensembl ID (e.g., 'ENSG00000012048') |
| `species` | string | 否 | Species (default: 'homo_sapiens') 默认 `homo_sapiens` |

??? note "上游说明（英文原文）"

    Search Ensembl for gene information by gene symbol or Ensembl ID. Returns gene location, biotype, and description.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'search_ensembl_gene',
    {'gene': '...'},
)
```

## `search_kegg_pathway` {#search_kegg_pathway}

Searches KEGG for metabolic and signaling pathways by pathway name or gene symbol.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`query`
- **关键词：**kegg, pathway search, gene symbol, organism code

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Pathway name or gene symbol |
| `organism` | string | 否 | Organism code (default: 'hsa' for human) 默认 `hsa` |

??? note "上游说明（英文原文）"

    Search KEGG for metabolic and signaling pathways by name or gene.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'search_kegg_pathway',
    {'query': '...'},
)
```

## `search_ncbi_gene` {#search_ncbi_gene}

Searches the NCBI Gene database by gene symbol, name or Gene ID.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`query`
- **关键词：**ncbi gene, entrez, gene symbol, gene id

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Gene symbol (e.g., 'TP53'), name, or NCBI Gene ID |
| `organism` | string | 否 | Organism to filter results (default: 'Homo sapiens') 默认 `Homo sapiens` |

??? note "上游说明（英文原文）"

    Search NCBI Gene database for gene information by gene symbol, name, or ID.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'search_ncbi_gene',
    {'query': '...'},
)
```

## `search_omim` {#search_omim}

Searches OMIM for Mendelian genetic diseases and phenotypes by gene symbol or disease name.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`query`
- **关键词：**omim, mendelian, genetic disease, gene-phenotype

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Gene symbol or disease name |
| `limit` | integer | 否 | Maximum results (default: 10) 默认 `10` |

??? note "上游说明（英文原文）"

    Search OMIM (Online Mendelian Inheritance in Man) for genetic diseases and phenotypes.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'search_omim',
    {'query': '...'},
)
```

## `search_pharmgkb` {#search_pharmgkb}

Searches PharmGKB for pharmacogenomic drug-gene-variant data by gene, drug or rsID.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`query`
- **关键词：**pharmgkb, pharmacogenomics, rsid, drug response

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Gene name, drug name, or rsID |
| `resource_type` | string | 否 | Type of resource: 'gene', 'drug', 'variant', or 'all' (default: 'all') 默认 `all` |

??? note "上游说明（英文原文）"

    Search PharmGKB for drug-gene interactions and pharmacogenomics data.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'search_pharmgkb',
    {'query': '...'},
)
```

## `vcf_quality_check` {#vcf_quality_check}

Performs quality control on a VCF file, filtering variants by MAF, missingness and Hardy-Weinberg equilibrium.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`vcf_path`
- **关键词：**vcf, qc, maf, hardy-weinberg, missingness

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `vcf_path` | string | 是 | Path to VCF file (.vcf or .vcf.gz) |
| `maf_threshold` | number | 否 | Minor allele frequency threshold (default: 0.01) 默认 `0.01` |
| `missingness_threshold` | number | 否 | Maximum allowed missingness rate (default: 0.1) 默认 `0.1` |
| `hwe_pvalue` | number | 否 | Hardy-Weinberg equilibrium p-value threshold (default: 1e-6) 默认 `1e-06` |
| `output_path` | string | 否 | Optional path to save QC report |

??? note "上游说明（英文原文）"

    Perform quality control on a VCF file. Filters variants based on MAF, missingness, and Hardy-Weinberg equilibrium.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'vcf_quality_check',
    {'vcf_path': '...'},
)
```
