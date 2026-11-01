# GWAS MCP (muslus) (`gwas-mcp`)

DNAKit's MCP tool index keeps **33** DNA/gene tools of this server (44 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [muslus/gwas-mcp @ `5d8160c`](https://github.com/muslus/gwas-mcp/tree/5d8160cda988e959234467cefae554ebdc1876e5) |
| License | MIT |
| Command | `python -m gwas_mcp.server` |
| Transports | stdio |
| Installation | Prepare the pinned muslus/gwas-mcp source and its dependencies in an isolated Python &gt;=3.10 environment, then point executable_path at that interpreter. Avoid assuming that the separately published gwas-mcp distribution has identical code. Pin 'mcp&lt;2': under mcp 2.x the server fails at start ('Server' has no list_tools). |
| Effects | Queries genomic evidence services; selected QC/association tools may execute local analyses and write files. |

| Tool | Category | Summary |
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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `rsids`
- **Keywords:** rsid, vep, gnomad, allele frequency, consequence

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `rsids` | array of string | yes | List of rsIDs to annotate (e.g., \['rs12345', 'rs67890'\]) |
| `include_frequencies` | boolean | no | Include population allele frequencies (default: true) Default `True` |

??? note "Upstream description"

    Annotate SNPs with gene names, functional consequences, and allele frequencies from Ensembl VEP and gnomAD.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `sumstats_path`
- **Keywords:** lambda gc, gwas, summary statistics, population stratification

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sumstats_path` | string | yes | Path to GWAS summary statistics file |
| `pvalue_column` | string | no | Name of p-value column (default: P) Default `P` |

??? note "Upstream description"

    Calculate genomic inflation factor (lambda GC) from GWAS summary statistics. Used to assess population stratification.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `sumstats_path`
- **Keywords:** ldsc, heritability, ld score, gwas, summary statistics

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sumstats_path` | string | yes | Path to GWAS summary statistics |
| `ld_scores_path` | string | no | Path to LD scores file |
| `sample_size` | integer | no | GWAS sample size (required if not in sumstats) |

??? note "Upstream description"

    Estimate SNP-heritability using LD Score regression. Requires GWAS summary statistics and LD scores.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `plink_prefix`
- **Keywords:** plink, genotype qc, call rate, missing data

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `plink_prefix` | string | yes | Path prefix for PLINK files |
| `output_path` | string | no | Optional path to save missingness report |

??? note "Upstream description"

    Calculate per-individual and per-variant missingness rates from genotype data.

Call it through DNAKit once the server is installed:

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

- **Category:** `breeding_population` (Breeding genotyping, markers, population genetics and strain resources)
- **Required arguments:** `plink_prefix`
- **Keywords:** plink, population structure, principal components, stratification

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `plink_prefix` | string | yes | Path prefix for PLINK files (without .bed/.bim/.fam extension) |
| `n_components` | integer | no | Number of principal components to compute (default: 10) Default `10` |
| `output_path` | string | no | Optional path to save PCA results |

??? note "Upstream description"

    Perform Principal Component Analysis for population stratification. Uses PLINK files as input.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `sumstats_path`, `target_plink_prefix`
- **Keywords:** polygenic risk score, prs, plink, clumping, p-value threshold

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sumstats_path` | string | yes | Path to GWAS summary statistics (base data) |
| `target_plink_prefix` | string | yes | Path prefix for target PLINK files |
| `pvalue_threshold` | number | no | P-value threshold for SNP inclusion (default: 1.0 - all SNPs) Default `1.0` |
| `clump_r2` | number | no | LD r² threshold for clumping (default: 0.1) Default `0.1` |
| `output_path` | string | no | Path to save PRS results |

??? note "Upstream description"

    Calculate Polygenic Risk Scores (PRS) for individuals using GWAS summary statistics and target genotypes.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `sumstats_path`
- **Keywords:** ld clumping, plink, index snp, gwas, independent signals

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sumstats_path` | string | yes | Path to GWAS summary statistics |
| `plink_prefix` | string | no | Path prefix for PLINK files (LD reference) |
| `p1` | number | no | P-value threshold for index SNPs (default: 5e-8) Default `5e-08` |
| `p2` | number | no | P-value threshold for clumped SNPs (default: 1e-4) Default `0.0001` |
| `r2` | number | no | LD r² threshold (default: 0.5) Default `0.5` |
| `kb` | integer | no | Clumping window in kb (default: 250) Default `250` |
| `output_path` | string | no | Path to save clumped results |

??? note "Upstream description"

    Perform LD-based clumping of GWAS results to identify independent significant signals.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `sumstats_path`
- **Keywords:** manhattan plot, gwas, p-value, visualization

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sumstats_path` | string | yes | Path to GWAS summary statistics file |
| `significance_threshold` | number | no | Genome-wide significance threshold (default: 5e-8) Default `5e-08` |
| `suggestive_threshold` | number | no | Suggestive significance threshold (default: 1e-5) Default `1e-05` |
| `output_path` | string | no | Path to save the plot (PNG format) |
| `title` | string | no | Plot title Default `Manhattan Plot` |

??? note "Upstream description"

    Create a Manhattan plot from GWAS summary statistics. Shows -log10(p-value) across chromosomes.

Call it through DNAKit once the server is installed:

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

- **Category:** `breeding_population` (Breeding genotyping, markers, population genetics and strain resources)
- **Required arguments:** `pca_results_path`
- **Keywords:** pca plot, population structure, visualization, ancestry

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `pca_results_path` | string | yes | Path to PCA results file (with PC1, PC2 columns) |
| `population_labels_path` | string | no | Optional path to population labels file |
| `output_path` | string | no | Path to save the plot |
| `title` | string | no | Plot title Default `PCA Plot` |

??? note "Upstream description"

    Create a PCA scatter plot showing population structure. Plots PC1 vs PC2.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `sumstats_path`
- **Keywords:** qq plot, gwas, p-value, inflation, visualization

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sumstats_path` | string | yes | Path to GWAS summary statistics file |
| `output_path` | string | no | Path to save the plot |
| `title` | string | no | Plot title Default `QQ Plot` |

??? note "Upstream description"

    Create a Quantile-Quantile (QQ) plot from GWAS p-values. Shows observed vs expected p-values.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `sumstats_path`, `lead_snp`
- **Keywords:** locuszoom, regional plot, lead snp, gwas

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sumstats_path` | string | yes | Path to GWAS summary statistics |
| `lead_snp` | string | yes | Lead SNP rsID or CHR:POS |
| `window_kb` | integer | no | Window size in kilobases (default: 500) Default `500` |
| `output_path` | string | no | Path to save the plot |

??? note "Upstream description"

    Create a regional association plot (LocusZoom-style) for a specific genomic region around a lead SNP.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `plink_prefix`, `lead_snp`
- **Keywords:** plink, genotype, fine-mapping, lead snp, region

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `plink_prefix` | string | yes | Path prefix for PLINK files |
| `lead_snp` | string | yes | Lead SNP identifier (rsID or CHR:POS) |
| `window_kb` | integer | no | Window size in kb (default: 500) Default `500` |
| `output_prefix` | string | no | Output prefix for extracted files |

??? note "Upstream description"

    Extract genotype data for a specific genomic region around a lead SNP for fine-mapping.

Call it through DNAKit once the server is installed:

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

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `genes`
- **Keywords:** enrichr, enrichment, gene set, kegg, go

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `genes` | array of string | yes | List of gene symbols (e.g., \['BRCA1', 'TP53', 'EGFR'\]) |
| `databases` | array of string | no | Databases to query (default: GO_Biological_Process, KEGG) Default `\['GO_Biological_Process_2021', 'KEGG_2021_Human'\]` |

??? note "Upstream description"

    Perform gene set enrichment analysis using GO and KEGG pathways via Enrichr.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `variant_id`
- **Keywords:** clinvar, pathogenicity, rsid, clinical significance

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variant_id` | string | yes | Variant identifier - rsID (e.g., 'rs80357906') or ClinVar ID |

??? note "Upstream description"

    Get detailed ClinVar information for a specific variant by rsID or ClinVar ID.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `gene`
- **Keywords:** open targets, gene-disease, association score, evidence

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string | yes | Gene symbol |
| `limit` | integer | no | Maximum number of diseases (default: 20) Default `20` |

??? note "Upstream description"

    Get disease associations for a gene from Open Targets with evidence scores.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `gene`
- **Keywords:** pharmgkb, pharmacogenomics, drug response, cyp2d6

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string | yes | Gene symbol (e.g., 'CYP2D6', 'VKORC1') |

??? note "Upstream description"

    Get drug-gene interaction annotations from PharmGKB for a specific gene.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `rsid`
- **Keywords:** gtex, eqtl, tissue, rsid, gene expression

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `rsid` | string | yes | rsID to query |
| `tissue` | string | no | GTEx tissue type (e.g., 'Whole_Blood', 'Brain_Cortex') |
| `gene` | string | no | Optional gene to filter results |

??? note "Upstream description"

    Get expression quantitative trait loci (eQTL) data from GTEx for a SNP and tissue.

Call it through DNAKit once the server is installed:

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

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `proteins`
- **Keywords:** string, enrichment, go, pathways, gene list

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `proteins` | array of string | yes | List of protein/gene names |
| `species` | integer | no | NCBI species taxonomy ID (default: 9606 for human) Default `9606` |

??? note "Upstream description"

    Get functional enrichment analysis for a list of proteins using STRING. Returns enriched GO terms, pathways, and domains.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `gene`
- **Keywords:** omim, mendelian disease, gene-disease, inheritance

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string | yes | Gene symbol (e.g., 'BRCA1', 'CFTR') |

??? note "Upstream description"

    Get genetic diseases associated with a gene from OMIM and other sources.

Call it through DNAKit once the server is installed:

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

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene`
- **Keywords:** kegg, pathways, gene function

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string | yes | Gene symbol (e.g., 'TP53', 'BRCA1') |

??? note "Upstream description"

    Find all KEGG pathways that a gene participates in.

Call it through DNAKit once the server is installed:

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

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `pathway_id`
- **Keywords:** kegg, pathway genes, pathway id, description

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `pathway_id` | string | yes | KEGG pathway ID (e.g., 'hsa04110' for cell cycle) |

??? note "Upstream description"

    Get detailed KEGG pathway information including genes and description.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `rsid`
- **Keywords:** ensembl, rsid, alleles, consequence, dbsnp

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `rsid` | string | yes | dbSNP rsID (e.g., 'rs1234567') |

??? note "Upstream description"

    Get detailed variant/SNP information from Ensembl by rsID. Returns position, alleles, clinical significance, and consequences.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `sumstats_path`
- **Keywords:** genome-wide significance, gwas, p-value, 5e-8, summary statistics

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sumstats_path` | string | yes | Path to GWAS summary statistics file |
| `pvalue_threshold` | number | no | P-value threshold for significance (default: 5e-8) Default `5e-08` |
| `suggestive_threshold` | number | no | Suggestive significance threshold (default: 1e-5) Default `1e-05` |
| `output_path` | string | no | Path to save significant SNPs |

??? note "Upstream description"

    Filter GWAS results to identify genome-wide significant SNPs based on p-value threshold.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `plink_prefix`
- **Keywords:** linkage disequilibrium, plink, r2, independent snps

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `plink_prefix` | string | yes | Path prefix for PLINK files |
| `r2_threshold` | number | no | r² threshold for LD pruning (default: 0.2) Default `0.2` |
| `window_size` | integer | no | Window size in kb for LD calculation (default: 500) Default `500` |
| `step_size` | integer | no | Step size in variants (default: 50) Default `50` |

??? note "Upstream description"

    Prune SNPs based on linkage disequilibrium (LD). Removes variants in high LD to create an independent set.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** none
- **Keywords:** gwas catalog, nhgri-ebi, trait, rsid, reported associations

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `rsid` | string | no | rsID to query (e.g., 'rs12345') |
| `trait` | string | no | Trait/disease to search for (e.g., 'diabetes') |
| `gene` | string | no | Gene name to search for (e.g., 'BRCA1') |

??? note "Upstream description"

    Query the NHGRI-EBI GWAS Catalog for previously reported associations for a SNP or trait.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `plink_prefix`, `phenotype_path`
- **Keywords:** gwas, plink, linear regression, logistic regression, summary statistics

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `plink_prefix` | string | yes | Path prefix for PLINK genotype files |
| `phenotype_path` | string | yes | Path to phenotype file (tab-separated with FID, IID, PHENO columns) |
| `covariate_path` | string | no | Optional path to covariate file |
| `model` | string | no | Regression model type (default: linear) One of `linear`, `logistic` Default `linear` |
| `output_path` | string | no | Path to save GWAS results |

??? note "Upstream description"

    Perform genome-wide association study using linear or logistic regression. Returns summary statistics including p-values, beta coefficients, and standard errors.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `query`
- **Keywords:** clinvar, pathogenic, benign, rsid, condition

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search query - can be gene name (e.g., 'BRCA1'), rsID (e.g., 'rs80357906'), or condition (e.g., 'breast cancer') |
| `limit` | integer | no | Maximum number of results (default: 10) Default `10` |

??? note "Upstream description"

    Search ClinVar for clinical variant interpretations. Find pathogenic/benign classifications for genetic variants.

Call it through DNAKit once the server is installed:

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

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene`
- **Keywords:** ensembl, gene location, biotype, gene symbol

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string | yes | Gene symbol (e.g., 'BRCA1') or Ensembl ID (e.g., 'ENSG00000012048') |
| `species` | string | no | Species (default: 'homo_sapiens') Default `homo_sapiens` |

??? note "Upstream description"

    Search Ensembl for gene information by gene symbol or Ensembl ID. Returns gene location, biotype, and description.

Call it through DNAKit once the server is installed:

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

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** kegg, pathway search, gene symbol, organism code

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Pathway name or gene symbol |
| `organism` | string | no | Organism code (default: 'hsa' for human) Default `hsa` |

??? note "Upstream description"

    Search KEGG for metabolic and signaling pathways by name or gene.

Call it through DNAKit once the server is installed:

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

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** ncbi gene, entrez, gene symbol, gene id

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Gene symbol (e.g., 'TP53'), name, or NCBI Gene ID |
| `organism` | string | no | Organism to filter results (default: 'Homo sapiens') Default `Homo sapiens` |

??? note "Upstream description"

    Search NCBI Gene database for gene information by gene symbol, name, or ID.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `query`
- **Keywords:** omim, mendelian, genetic disease, gene-phenotype

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Gene symbol or disease name |
| `limit` | integer | no | Maximum results (default: 10) Default `10` |

??? note "Upstream description"

    Search OMIM (Online Mendelian Inheritance in Man) for genetic diseases and phenotypes.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `query`
- **Keywords:** pharmgkb, pharmacogenomics, rsid, drug response

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Gene name, drug name, or rsID |
| `resource_type` | string | no | Type of resource: 'gene', 'drug', 'variant', or 'all' (default: 'all') Default `all` |

??? note "Upstream description"

    Search PharmGKB for drug-gene interactions and pharmacogenomics data.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `vcf_path`
- **Keywords:** vcf, qc, maf, hardy-weinberg, missingness

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `vcf_path` | string | yes | Path to VCF file (.vcf or .vcf.gz) |
| `maf_threshold` | number | no | Minor allele frequency threshold (default: 0.01) Default `0.01` |
| `missingness_threshold` | number | no | Maximum allowed missingness rate (default: 0.1) Default `0.1` |
| `hwe_pvalue` | number | no | Hardy-Weinberg equilibrium p-value threshold (default: 1e-6) Default `1e-06` |
| `output_path` | string | no | Optional path to save QC report |

??? note "Upstream description"

    Perform quality control on a VCF file. Filters variants based on MAF, missingness, and Hardy-Weinberg equilibrium.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gwas-mcp',
    'vcf_quality_check',
    {'vcf_path': '...'},
)
```
