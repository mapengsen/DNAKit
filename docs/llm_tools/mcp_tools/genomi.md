# Genomi (`genomi`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **57** 个 DNA/基因相关工具（服务器共列出 112 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [exon-research/genomi @ `1df4f5b`](https://github.com/exon-research/genomi/tree/1df4f5bfb635e2d08f93a1102b70b840fbe80d86) |
| 许可证 | Apache-2.0 |
| 启动命令 | `genomi serve` |
| 传输方式 | stdio |
| 安装 | Install the exon-research/genomi source project (not on PyPI; Python &gt;=3.10, depends on a git-pinned gxl-paperclip) so that the genomi console script exists, then run genomi install to materialize reference libraries; optionally set GENOMI_HOME to relocate its data directory (default ~/.genomi). |
| 副作用 | Parses user-supplied genome files into a local Active Genome Index and evidence database under GENOMI_HOME, keeps a local journal, runs background jobs, and genomi.install downloads public reference libraries; returned results can contain individual genotypes, and optional lab operations may contact external evidence providers. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`active_genome_index.classify_callset_qc`](#active_genome_index.classify_callset_qc) | `variant` | Classifies a genome callset's shape, depth/quality field availability and limits on absence claims using the genome index. |
| [`active_genome_index.classify_genotype_support`](#active_genome_index.classify_genotype_support) | `variant` | Classifies whether one exact allele (chrom/pos/ref/alt) has enough sample read/genotype support for personal interpretation. |
| [`active_genome_index.classify_region_callability`](#active_genome_index.classify_region_callability) | `variant` | Classifies whether a genomic region in the personal genome is callable enough to support reference or absence claims. |
| [`ancestry.check_sample_overlap`](#ancestry.check_sample_overlap) | `variant` | Checks how many installed 1000 Genomes ancestry panel markers are present and usable in an approved personal genome index. |
| [`ancestry.estimate_population_context`](#ancestry.estimate_population_context) | `variant` | Estimates qualitative similarity of a GRCh37/GRCh38 personal genome to reference populations via local 1000 Genomes PCA projection. |
| [`ancestry.project_pca`](#ancestry.project_pca) | `variant` | Projects an approved personal genome into the installed 1000 Genomes ancestry PCA space and returns nearest reference samples. |
| [`cell_type.retrieve_markers`](#cell_type.retrieve_markers) | `expression` | Retrieves canonical marker genes for a cell type from HPA, CellMarker, PanglaoDB or ENCODE sources. |
| [`clinvar.match_variants`](#clinvar.match_variants) | `variant` | Materializes exact ClinVar matches for comparable variants in a genome index using the installed build-specific ClinVar library. |
| [`clinvar.scan_candidates`](#clinvar.scan_candidates) | `variant` | Builds a deterministic ClinVar variant evidence inventory and review groups from materialized personal genome matches. |
| [`decode.build_dashboard_evidence`](#decode.build_dashboard_evidence) | `variant` | Gathers dashboard panel evidence (variants, pharmacogenomics, risk scores, ancestry, nutrigenomics) from the approved personal genome index. |
| [`decode.render_dashboard`](#decode.render_dashboard) | `variant` | Builds and renders a self-contained HTML dashboard summarizing an approved personal genome's variants, PGx, risk and ancestry. |
| [`functional_genomics.compare_gene_perturbation`](#functional_genomics.compare_gene_perturbation) | `association` | Compares candidate genes using verified functional-genomics perturbation experiment evidence for a given context. |
| [`functional_genomics.import_perturbation_table`](#functional_genomics.import_perturbation_table) | `association` | Extracts verified gene perturbation experiment records from a local CSV/TSV result table with a column mapping. |
| [`functional_genomics.query_geo`](#functional_genomics.query_geo) | `sequencing_data` | Discovers public NCBI GEO study tables for a functional-genomics context, genes, organism, cell line or perturbation. |
| [`functional_genomics.retrieve_perturbation_records`](#functional_genomics.retrieve_perturbation_records) | `association` | Retrieves public functional-genomics perturbation records for candidate genes from BioGRID ORCS and DepMap. |
| [`genomi.parse_source`](#genomi.parse_source) | `variant` | Detects and parses a genome source such as VCF/gVCF, BAM, paired FASTQ or consumer genotype export into a local Active Genome Index. |
| [`genomilab.verify_sequence_substitution`](#genomilab.verify_sequence_substitution) | `variant` | Verifies a gene's coding change and protein substitution against a public reference protein sequence and stores digests as a research artifact. |
| [`gnomad.fetch_population_frequency`](#gnomad.fetch_population_frequency) | `variant` | Fetches gnomAD public population allele frequency for one allele (chrom/pos/ref/alt) and stores it as evidence. |
| [`gwas.compare_gene_associations`](#gwas.compare_gene_associations) | `association` | Compares candidate genes using GWAS Catalog reported- and mapped-gene trait association evidence for a phenotype. |
| [`gwas.compare_variant_associations`](#gwas.compare_variant_associations) | `association` | Compares candidate rsIDs by population-trait association evidence from the GWAS Catalog. |
| [`nutrigenomics.retrieve_domain_markers`](#nutrigenomics.retrieve_domain_markers) | `association` | Retrieves curated single-marker genetic variant records for a nutrigenomic domain, filtered by minimum evidence tier. |
| [`nutrigenomics.retrieve_variant_records`](#nutrigenomics.retrieve_variant_records) | `association` | Retrieves nutrigenomic catalogue records that reference a specific rsID. |
| [`pathway.retrieve_members`](#pathway.retrieve_members) | `genome_annotation` | Retrieves canonical member genes of a Reactome, KEGG or MSigDB Hallmark pathway/gene set. |
| [`pharmacogenomics.describe_gene_requirements`](#pharmacogenomics.describe_gene_requirements) | `variant` | Returns pharmacogene-specific sample evidence requirements for PharmCAT allele matching, outside calls, CYP2D6 SV/CNV, HLA, MT-RNR1 and G6PD. |
| [`pharmacogenomics.fetch_clinpgx`](#pharmacogenomics.fetch_clinpgx) | `variant` | Fetches ClinPGx (PharmGKB) pharmacogenomic guidelines, clinical annotations and label evidence for a drug, gene or rsID. |
| [`pharmacogenomics.fetch_pgxdb`](#pharmacogenomics.fetch_pgxdb) | `variant` | Fetches PGxDB pharmacogenomic evidence for a drug, ATC code, DrugBank ID, rsID, variant marker or gene. |
| [`pharmacogenomics.import_pharmcat_artifacts`](#pharmacogenomics.import_pharmcat_artifacts) | `variant` | Imports existing PharmCAT report, calls and match artifacts (diplotypes, phenotypes, recommendations) without running PharmCAT. |
| [`pharmacogenomics.preflight_pharmcat`](#pharmacogenomics.preflight_pharmcat) | `variant` | Inspects the selected genome index (sample columns, format) for PharmCAT PGx calling readiness without running PharmCAT. |
| [`pharmacogenomics.prepare_outside_call_tsv`](#pharmacogenomics.prepare_outside_call_tsv) | `variant` | Prepares a PharmCAT outside-call TSV from specialized caller output such as OptiType HLA or StellarPGx CYP2D6 diplotypes. |
| [`pharmacogenomics.review_medication`](#pharmacogenomics.review_medication) | `variant` | Reviews pharmacogenomic evidence for a medication, drug-gene or drug-variant target, using sample and stored evidence. |
| [`pharmacogenomics.run_pharmcat`](#pharmacogenomics.run_pharmcat) | `variant` | Runs a local PharmCAT installation on an approved genome index to call pharmacogene diplotypes and phenotypes. |
| [`pharmacogenomics.validate_outside_call_tsv`](#pharmacogenomics.validate_outside_call_tsv) | `variant` | Validates PharmCAT outside-call TSV structure and summarizes diplotype, phenotype or activity-score evidence. |
| [`phenotype.compare_disease_evidence`](#phenotype.compare_disease_evidence) | `association` | Compares supplied or gene-derived candidate diseases against phenotype/HPO evidence without selecting a diagnosis. |
| [`phenotype.compare_drug_target_evidence`](#phenotype.compare_drug_target_evidence) | `association` | Compares candidate genes using direct drug-target or mechanism-of-action evidence for a drug, class or indication. |
| [`phenotype.compare_gene_hpo_evidence`](#phenotype.compare_gene_hpo_evidence) | `association` | Compares candidate genes using phenotype, HPO and curated rare-disease annotation evidence. |
| [`phenotype.plan_risk_investigation`](#phenotype.plan_risk_investigation) | `variant` | Plans phenotype risk and ClinVar review targets (rare disease, cancer risk, carrier) from genes, conditions and optional genome evidence. |
| [`phenotype.retrieve_disease_drug_targets`](#phenotype.retrieve_disease_drug_targets) | `association` | Retrieves disease-scoped clinical drug-target genes from Open Targets drug candidate records. |
| [`phenotype.retrieve_gene_disease_associations`](#phenotype.retrieve_gene_disease_associations) | `association` | Retrieves primary gene-disease validity associations from GenCC for supplied gene symbols. |
| [`phenotype.retrieve_trait_gene_records`](#phenotype.retrieve_trait_gene_records) | `association` | Retrieves trait-to-gene records from integrated public sources such as Open Targets, optionally filtered to given genes. |
| [`prs.calculate_score`](#prs.calculate_score) | `variant` | Applies a published polygenic score to an approved genome index and returns the raw weighted score plus QC. |
| [`prs.check_score_overlap`](#prs.check_score_overlap) | `variant` | Checks how many variants from a polygenic score are present and usable in an approved genome index. |
| [`prs.fetch_score_metadata`](#prs.fetch_score_metadata) | `variant` | Fetches detailed public PGS Catalog metadata for one score ID, including scoring-file URLs and publication context. |
| [`prs.import_scoring_file`](#prs.import_scoring_file) | `variant` | Imports a PGS Catalog or local scoring file (variant weights) into the local PRS score cache for a genome build. |
| [`prs.search_scores`](#prs.search_scores) | `variant` | Searches public PGS Catalog score metadata by trait, score ID, EFO term or free text. |
| [`region.retrieve_features`](#region.retrieve_features) | `genome_annotation` | Retrieves GENCODE gene/transcript and ENCODE cCRE feature annotations overlapping a genomic region. |
| [`research.build_target_packet`](#research.build_target_packet) | `association` | Builds a target-centric evidence packet for a variant, gene, drug, condition or topic target. |
| [`sequence.analyze`](#sequence.analyze) | `sequence` | Runs a compact DNA sequence analysis bundle (summary, translation, ORFs, restriction sites, Kozak context). |
| [`sequence.check_primers`](#sequence.check_primers) | `molecular_tools` | Checks basic primer properties (e.g. length, GC, Tm) and optional amplicons on a template sequence. |
| [`sequence.classify_kozak`](#sequence.classify_kozak) | `sequence` | Classifies Kozak sequence context around ATG start codons in a supplied sequence. |
| [`sequence.find_orfs`](#sequence.find_orfs) | `sequence` | Finds ATG-to-stop open reading frames above a minimum length on one or both strands of a DNA sequence. |
| [`sequence.find_restriction_sites`](#sequence.find_restriction_sites) | `molecular_tools` | Finds common restriction enzyme or custom motif sites in a supplied DNA sequence. |
| [`sequence.match_reference`](#sequence.match_reference) | `sequence` | Matches a supplied DNA sequence against local FASTA records and returns record identifiers and annotations. |
| [`sequence.translate`](#sequence.translate) | `sequence` | Translates a DNA sequence in a selected frame and strand using the standard genetic code. |
| [`variant.find_gene_variants`](#variant.find_gene_variants) | `variant` | Finds passing variants in an approved genome index that overlap GENCODE intervals of up to 10 candidate genes. |
| [`variant.gather_allele_context`](#variant.gather_allele_context) | `variant` | Gathers existing sample, static, population-frequency and reviewed research evidence for one allele (chrom/pos/ref/alt). |
| [`variant.gather_gene_context`](#variant.gather_gene_context) | `variant` | Gathers existing sample, ClinVar and reviewed research evidence for one gene. |
| [`variant.resolve`](#variant.resolve) | `variant` | Resolves one variant target (coordinates, rsID, query or region) and returns public, local sample and stored evidence facts. |

## `active_genome_index.classify_callset_qc` {#active_genome_index.classify_callset_qc}

Classifies a genome callset's shape, depth/quality field availability and limits on absence claims using the genome index.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**callset qc, depth, genotype quality, gvcf, absence claims

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `db` | string | 否 | Local filesystem path. |
| `agi_path` | string | 否 | Local filesystem path. |
| `output` | string | 否 | Local filesystem path. |
| `genome_build` | string | 否 | 默认 `auto` |
| `scan_records` | integer | 否 | 默认 `1000` |

??? note "上游说明（英文原文）"

    Classify genome callset shape, depth/quality field availability, and absence-claim boundaries using an Active Genome Index.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'active_genome_index.classify_callset_qc', 'params': {}},
)
```

## `active_genome_index.classify_genotype_support` {#active_genome_index.classify_genotype_support}

Classifies whether one exact allele (chrom/pos/ref/alt) has enough sample read/genotype support for personal interpretation.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`chrom`, `pos`, `ref`, `alt`
- **关键词：**genotype support, allele, read depth, personal genome, zygosity

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chrom` | string | 是 |  |
| `pos` | integer | 是 |  |
| `ref` | string | 是 |  |
| `alt` | string | 是 |  |
| `db` | string | 否 | Local filesystem path. |
| `agi_path` | string | 否 | Local filesystem path. |
| `output` | string | 否 | Local filesystem path. |
| `genome_build` | string | 否 | 默认 `auto` |
| `reference_fasta` | string | 否 | Local filesystem path. |
| `min_depth` | integer | 否 | 默认 `10` |
| `min_gq` | integer | 否 | 默认 `20` |

??? note "上游说明（英文原文）"

    Classify whether one exact allele has enough sample support to be used in a personal interpretation.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'active_genome_index.classify_genotype_support',
     'params': {'chrom': '...', 'pos': 0, 'ref': '...', 'alt': '...'}},
)
```

## `active_genome_index.classify_region_callability` {#active_genome_index.classify_region_callability}

Classifies whether a genomic region in the personal genome is callable enough to support reference or absence claims.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`region`
- **关键词：**callability, coverage, reference call, absence, region

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `region` | string | 是 |  |
| `db` | string | 否 | Local filesystem path. |
| `agi_path` | string | 否 | Local filesystem path. |
| `output` | string | 否 | Local filesystem path. |
| `genome_build` | string | 否 | 默认 `auto` |
| `min_depth` | integer | 否 | 默认 `10` |
| `min_covered_fraction` | number | 否 | 默认 `0.95` |
| `limit` | integer | 否 | 默认 `5000` |

??? note "上游说明（英文原文）"

    Classify whether a region can support reference or absence claims.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'active_genome_index.classify_region_callability',
     'params': {'region': '...'}},
)
```

## `ancestry.check_sample_overlap` {#ancestry.check_sample_overlap}

Checks how many installed 1000 Genomes ancestry panel markers are present and usable in an approved personal genome index.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**ancestry markers, marker overlap, 1000 genomes, snp, personal genome

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `agi_path` | string | 否 | Local filesystem path. |
| `genome_build` | string | 否 | 默认 `GRCh38` |

??? note "上游说明（英文原文）"

    Check how many installed 1000 Genomes ancestry panel markers are usable in an approved Active Genome Index.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'ancestry.check_sample_overlap', 'params': {}},
)
```

## `ancestry.estimate_population_context` {#ancestry.estimate_population_context}

Estimates qualitative similarity of a GRCh37/GRCh38 personal genome to reference populations via local 1000 Genomes PCA projection.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**genetic ancestry, population structure, pca, 1000 genomes, superpopulation

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `agi_path` | string | 否 | Local filesystem path. |
| `genome_build` | string | 否 | 默认 `GRCh38` |
| `nearest_reference_count` | integer | 否 | 默认 `10` |

??? note "上游说明（英文原文）"

    Estimate qualitative reference-panel similarity for an approved GRCh37 or GRCh38 sample using local 1000 Genomes PCA projection.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'ancestry.estimate_population_context', 'params': {}},
)
```

## `ancestry.project_pca` {#ancestry.project_pca}

Projects an approved personal genome into the installed 1000 Genomes ancestry PCA space and returns nearest reference samples.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**pca projection, genetic ancestry, nearest neighbors, 1000 genomes, genotypes

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `agi_path` | string | 否 | Local filesystem path. |
| `genome_build` | string | 否 | 默认 `GRCh38` |
| `nearest_reference_count` | integer | 否 | 默认 `10` |

??? note "上游说明（英文原文）"

    Project an approved sample into the installed 1000 Genomes ancestry PCA space and return nearest reference neighbors.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'ancestry.project_pca', 'params': {}},
)
```

## `cell_type.retrieve_markers` {#cell_type.retrieve_markers}

Retrieves canonical marker genes for a cell type from HPA, CellMarker, PanglaoDB or ENCODE sources.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**无
- **关键词：**marker genes, cell type, cellmarker, panglaodb, human protein atlas

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `cell_type_id_or_name` | string | 否 |  |
| `cell_type_id` | string | 否 |  |
| `cell_type_name` | string | 否 |  |
| `source` | string | 否 | 可选值：`hpa`, `cellmarker`, `panglaodb`, `encode` |
| `species` | string | 否 |  |
| `marker_table` | string | 否 | Local filesystem path. |
| `hpa_api_base` | string | 否 | 默认 `https://www.proteinatlas.org/api` |
| `limit` | integer | 否 | 默认 `100` |
| `semantic_context` | object | 否 | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "上游说明（英文原文）"

    Retrieve canonical marker genes for a controlled cell-type source entity.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'cell_type.retrieve_markers', 'params': {}},
)
```

## `clinvar.match_variants` {#clinvar.match_variants}

Materializes exact ClinVar matches for comparable variants in a genome index using the installed build-specific ClinVar library.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**clinvar, pathogenicity, variant matching, clinical significance, personal genome

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `agi_path` | string | 否 | Local filesystem path. |
| `db` | string | 否 | Local filesystem path. |
| `output` | string | 否 | Local filesystem path. |
| `genome_build` | string | 否 | 默认 `auto` |
| `force` | boolean | 否 | 默认 `False` |

??? note "上游说明（英文原文）"

    Materialize exact ClinVar matches for comparable Active Genome Index variants using the installed build-specific ClinVar library.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'clinvar.match_variants', 'params': {}},
)
```

## `clinvar.scan_candidates` {#clinvar.scan_candidates}

Builds a deterministic ClinVar variant evidence inventory and review groups from materialized personal genome matches.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**clinvar, pathogenic, carrier, review groups, variant triage

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `agi_path` | string | 否 | Local filesystem path. |
| `db` | string | 否 | Local filesystem path. |
| `output` | string | 否 | Local filesystem path. |
| `genome_build` | string | 否 | 默认 `GRCh38` |
| `force` | boolean | 否 | 默认 `False` |

??? note "上游说明（英文原文）"

    Build deterministic ClinVar variant evidence inventory and review groups from materialized matches.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'clinvar.scan_candidates', 'params': {}},
)
```

## `decode.build_dashboard_evidence` {#decode.build_dashboard_evidence}

Gathers dashboard panel evidence (variants, pharmacogenomics, risk scores, ancestry, nutrigenomics) from the approved personal genome index.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**dashboard, personal genome report, pgx, polygenic risk, ancestry

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `panels` | array of string | 否 | Optional subset of dashboard panels to gather. Omitted means every dashboard panel. |
| `risk_score_ids` | array of string | 否 | Optional PGS IDs to score for the risk panel. Omitted means use imported local scores up to risk_score_limit. |
| `risk_score_limit` | integer | 否 | Maximum imported local PGS scores to calculate when risk_score_ids is omitted. 默认 `5` |
| `risk_review_types` | array of string | 否 | Phenotype review modes to add to the risk panel from selected Active Genome Index ClinVar review groups. Empty array disables phenotype review rows. 默认 `\['carrier_review', 'observed_condition_review'\]` |
| `pgx_review_targets` | array of object | 否 | Optional medication-review targets to add to the PGx panel. |
| `pgx_review_target_limit` | integer | 否 | Maximum medication-review targets to run for the PGx panel. 默认 `12` |
| `nutrigenomics_domain_ids` | array of string | 否 | Optional nutrigenomics domains to retrieve. Omitted means every declared in-scope domain. |

??? note "上游说明（英文原文）"

    Build Genomi Dashboard panel states from the approved Active Genome Index and existing capability operations. Used by decode.render_dashboard; presented output reports panel readiness, running jobs, and gaps.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'decode.build_dashboard_evidence', 'params': {}},
)
```

## `decode.render_dashboard` {#decode.render_dashboard}

Builds and renders a self-contained HTML dashboard summarizing an approved personal genome's variants, PGx, risk and ancestry.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**dashboard, html report, personal genome, pgx, polygenic risk

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `panels` | array of string | 否 | Optional subset of dashboard panels to gather. Omitted means every dashboard panel. |
| `risk_score_ids` | array of string | 否 | Optional PGS IDs to score for the risk panel. Omitted means use imported local scores up to risk_score_limit. |
| `risk_score_limit` | integer | 否 | Maximum imported local PGS scores to calculate when risk_score_ids is omitted. 默认 `5` |
| `risk_review_types` | array of string | 否 | Phenotype review modes to add to the risk panel from selected Active Genome Index ClinVar review groups. Empty array disables phenotype review rows. 默认 `\['carrier_review', 'observed_condition_review'\]` |
| `pgx_review_targets` | array of object | 否 | Optional medication-review targets to add to the PGx panel. |
| `pgx_review_target_limit` | integer | 否 | Maximum medication-review targets to run for the PGx panel. 默认 `12` |
| `nutrigenomics_domain_ids` | array of string | 否 | Optional nutrigenomics domains to retrieve. Omitted means every declared in-scope domain. |
| `output` | string | 否 | Optional override for the dashboard HTML path. Defaults to a transient path under the system temp dir: &lt;tmp&gt;/genomi-dashboards/&lt;sample&gt;/dashboard.html. |

??? note "上游说明（英文原文）"

    Build, shape, render, and expose the Genomi Dashboard self-contained HTML artifact from the approved Active Genome Index.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'decode.render_dashboard', 'params': {}},
)
```

## `functional_genomics.compare_gene_perturbation` {#functional_genomics.compare_gene_perturbation}

Compares candidate genes using verified functional-genomics perturbation experiment evidence for a given context.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`context`, `genes`
- **关键词：**crispr screen, knockout, gene essentiality, perturbation, cell line

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `source_records` | array of object | 否 |  |
| `context` | string | 是 |  |
| `genes` | array of string | 是 |  |
| `organism` | string | 否 |  |
| `cell_line` | string | 否 |  |
| `perturbation` | string | 否 |  |
| `assay` | string | 否 |  |
| `phenotype` | string | 否 |  |
| `perturbation_sources` | array of string | 否 |  |
| `biogrid_orcs_access_key` | string | 否 |  |
| `depmap_gene_effect_url` | string | 否 |  |
| `depmap_model_url` | string | 否 |  |
| `db` | string | 否 | Local filesystem path. |
| `search_stored_research` | boolean | 否 | 默认 `True` |
| `retrieve_native` | boolean | 否 | 默认 `True` |
| `limit` | integer | 否 | 默认 `25` |
| `semantic_context` | object | 否 | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "上游说明（英文原文）"

    Compare candidate genes by verified functional-genomics perturbation experiment evidence.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'functional_genomics.compare_gene_perturbation',
     'params': {'context': '...', 'genes': []}},
)
```

## `functional_genomics.import_perturbation_table` {#functional_genomics.import_perturbation_table}

Extracts verified gene perturbation experiment records from a local CSV/TSV result table with a column mapping.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`table`, `context`, `genes`
- **关键词：**perturbation screen, crispr, gene hits, csv import, phenotype

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `context` | string | 是 |  |
| `genes` | array of string | 是 |  |
| `organism` | string | 否 |  |
| `cell_line` | string | 否 |  |
| `perturbation` | string | 否 |  |
| `assay` | string | 否 |  |
| `phenotype` | string | 否 |  |
| `table` | string | 是 | Local filesystem path. |
| `column_map` | object | 否 |  |
| `delimiter` | string | 否 |  |
| `source_title` | string | 否 |  |
| `source_url` | string | 否 |  |
| `source_type` | string | 否 |  |
| `limit` | integer | 否 | 默认 `500` |

??? note "上游说明（英文原文）"

    Extract verified perturbation experiment source records from a local CSV or TSV result table.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'functional_genomics.import_perturbation_table',
     'params': {'table': '...', 'context': '...', 'genes': []}},
)
```

## `functional_genomics.query_geo` {#functional_genomics.query_geo}

Discovers public NCBI GEO study tables for a functional-genomics context, genes, organism, cell line or perturbation.

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`context`
- **关键词：**geo, gse accession, perturbation, dataset discovery, ncbi

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `context` | string | 是 |  |
| `genes` | array of string | 否 |  |
| `organism` | string | 否 |  |
| `cell_line` | string | 否 |  |
| `perturbation` | string | 否 |  |
| `assay` | string | 否 |  |
| `phenotype` | string | 否 |  |
| `accession` | string | 否 |  |
| `limit` | integer | 否 | 默认 `25` |
| `semantic_context` | object | 否 | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |
| `ncbi_api_key` | string | 否 |  |
| `ncbi_email` | string | 否 |  |
| `ncbi_tool` | string | 否 |  |

??? note "上游说明（英文原文）"

    Discover public study tables from NCBI GEO when curated perturbation sources do not cover a requested functional-genomics context.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'functional_genomics.query_geo', 'params': {'context': '...'}},
)
```

## `functional_genomics.retrieve_perturbation_records` {#functional_genomics.retrieve_perturbation_records}

Retrieves public functional-genomics perturbation records for candidate genes from BioGRID ORCS and DepMap.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`context`, `genes`
- **关键词：**biogrid orcs, depmap, crispr screen, gene effect, essentiality

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `context` | string | 是 |  |
| `genes` | array of string | 是 |  |
| `organism` | string | 否 |  |
| `cell_line` | string | 否 |  |
| `perturbation` | string | 否 |  |
| `assay` | string | 否 |  |
| `phenotype` | string | 否 |  |
| `perturbation_sources` | array of string | 否 |  |
| `biogrid_orcs_access_key` | string | 否 |  |
| `depmap_gene_effect_url` | string | 否 |  |
| `depmap_model_url` | string | 否 |  |
| `limit` | integer | 否 | 默认 `100` |
| `semantic_context` | object | 否 | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "上游说明（英文原文）"

    Retrieve native public functional-genomics perturbation records from BioGRID ORCS and DepMap for candidate genes and declared experimental context.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'functional_genomics.retrieve_perturbation_records',
     'params': {'context': '...', 'genes': []}},
)
```

## `genomi.parse_source` {#genomi.parse_source}

Detects and parses a genome source such as VCF/gVCF, BAM, paired FASTQ or consumer genotype export into a local Active Genome Index.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`source`
- **关键词：**personal genome, vcf, gvcf, bam, fastq, consumer genotype, genome build

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `db` | string | 否 | Local filesystem path. |
| `shared_db` | string | 否 | Local filesystem path. |
| `source` | string | 是 | Local filesystem path. |
| `source_evidence_db` | string | 否 | Local filesystem path. |
| `reference_fasta` | string | 否 | Local filesystem path. |
| `auto_reference_fasta` | boolean | 否 | 默认 `True` |
| `genome_build` | string | 否 | 默认 `auto` |
| `max_records` | integer | 否 |  |
| `parallel_workers` | integer | 否 |  |
| `force` | boolean | 否 | 默认 `False` |
| `user_nickname` | string | 否 |  |
| `set_default_user` | boolean | 否 | 默认 `False` |

??? note "上游说明（英文原文）"

    Detect, parse, and digitize a genome source such as VCF/gVCF, BAM, paired FASTQ, supported consumer exports, or a genome.computer .genome/1.0 bundle. After parsing, read skills/active-genome-index/SKILL.md before using the invoke-only active_genome_index.* interpretation/selection tools; the result's next_actions repeats this.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.parse_source',
    {'source': '...'},
)
```

## `genomilab.verify_sequence_substitution` {#genomilab.verify_sequence_substitution}

Verifies a gene's coding change and protein substitution against a public reference protein sequence and stores digests as a research artifact.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`investigation_id`, `round_id`, `deduplication_key`, `gene`, `transcript_accession`, `protein_accession`, `coding_change`, `protein_substitution`, `public_reference_protein_sequence`, `reference_source_label`, `reference_source_version`, `reference_source_record_id`
- **关键词：**hgvs, coding change, missense, protein substitution, transcript accession

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `investigation_id` | string | 是 |  |
| `round_id` | string | 是 |  |
| `deduplication_key` | string | 是 |  |
| `gene` | string | 是 |  |
| `transcript_accession` | string | 是 |  |
| `protein_accession` | string | 是 |  |
| `coding_change` | string | 是 |  |
| `protein_substitution` | string | 是 |  |
| `public_reference_protein_sequence` | string | 是 |  |
| `reference_source_label` | string | 是 |  |
| `reference_source_version` | string | 是 |  |
| `reference_source_record_id` | string | 是 |  |

??? note "上游说明（英文原文）"

    Verify one intended protein substitution against a supplied public reference protein sequence using deterministic local Genomi rules. Persist only sequence digests and normalized descriptors as a round-bound nonclinical research artifact; never store the sequence or treat the result as clinical evidence.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomilab.verify_sequence_substitution',
    {'investigation_id': '...',
     'round_id': '...',
     'deduplication_key': '...',
     'gene': '...',
     'transcript_accession': '...',
     'protein_accession': '...',
     'coding_change': '...',
     'protein_substitution': '...',
     'public_reference_protein_sequence': '...',
     'reference_source_label': '...',
     'reference_source_version': '...',
     'reference_source_record_id': '...'},
)
```

## `gnomad.fetch_population_frequency` {#gnomad.fetch_population_frequency}

Fetches gnomAD public population allele frequency for one allele (chrom/pos/ref/alt) and stores it as evidence.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`chrom`, `pos`, `ref`, `alt`
- **关键词：**gnomad, allele frequency, population frequency, af, chrom pos ref alt

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `db` | string | 否 | Local filesystem path. |
| `shared_db` | string | 否 | Local filesystem path. |
| `chrom` | string | 是 |  |
| `pos` | integer | 是 |  |
| `ref` | string | 是 |  |
| `alt` | string | 是 |  |
| `sync_shared` | boolean | 否 | 默认 `True` |
| `dataset` | string | 否 | 默认 `gnomad_r4` |
| `genome_build` | string | 否 | 默认 `GRCh38` |
| `api_url` | string | 否 | 默认 `https://gnomad.broadinstitute.org/api` |

??? note "上游说明（英文原文）"

    Fetch reusable gnomAD public population frequency for one allele and write it into evidence storage.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'gnomad.fetch_population_frequency',
     'params': {'chrom': '...', 'pos': 0, 'ref': '...', 'alt': '...'}},
)
```

## `gwas.compare_gene_associations` {#gwas.compare_gene_associations}

Compares candidate genes using GWAS Catalog reported- and mapped-gene trait association evidence for a phenotype.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`phenotype`, `genes`
- **关键词：**gwas catalog, mapped gene, trait association, phenotype, candidate genes

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `phenotype` | string | 是 |  |
| `genes` | array of string | 是 |  |
| `association_limit` | integer | 否 | 默认 `200` |
| `api_url` | string | 否 | 默认 `https://www.ebi.ac.uk/gwas/rest/api/v2` |
| `semantic_context` | object | 否 | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "上游说明（英文原文）"

    Compare candidate genes using GWAS Catalog reported_gene and mapped_gene trait-association evidence.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'gwas.compare_gene_associations',
     'params': {'phenotype': '...', 'genes': []}},
)
```

## `gwas.compare_variant_associations` {#gwas.compare_variant_associations}

Compares candidate rsIDs by population-trait association evidence from the GWAS Catalog.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`phenotype`, `variants`
- **关键词：**gwas catalog, rsid, snp, trait association, p-value

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `phenotype` | string | 是 |  |
| `variants` | array of string | 是 |  |
| `association_limit` | integer | 否 | 默认 `200` |
| `api_url` | string | 否 | 默认 `https://www.ebi.ac.uk/gwas/rest/api` |
| `semantic_context` | object | 否 | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "上游说明（英文原文）"

    Compare candidate rsIDs by population-trait GWAS Catalog association evidence.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'gwas.compare_variant_associations',
     'params': {'phenotype': '...', 'variants': []}},
)
```

## `nutrigenomics.retrieve_domain_markers` {#nutrigenomics.retrieve_domain_markers}

Retrieves curated single-marker genetic variant records for a nutrigenomic domain, filtered by minimum evidence tier.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`domain_id`
- **关键词：**nutrigenomics, genetic markers, diet, nutrient metabolism, snp

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `domain_id` | string | 是 |  |
| `min_evidence_tier` | string | 否 | 可选值：`established`, `probable`, `emerging` 默认 `established` |
| `semantic_context` | object | 否 | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "上游说明（英文原文）"

    Retrieve curated single-marker records for a declared nutrigenomic domain, filtered by minimum evidence tier.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'nutrigenomics.retrieve_domain_markers',
     'params': {'domain_id': '...'}},
)
```

## `nutrigenomics.retrieve_variant_records` {#nutrigenomics.retrieve_variant_records}

Retrieves nutrigenomic catalogue records that reference a specific rsID.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`rsid`
- **关键词：**rsid, nutrigenomics, diet, snp, trait

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `rsid` | string | 是 |  |

??? note "上游说明（英文原文）"

    Retrieve any nutrigenomic catalogue records referencing a specific rsID.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'nutrigenomics.retrieve_variant_records',
     'params': {'rsid': '...'}},
)
```

## `pathway.retrieve_members` {#pathway.retrieve_members}

Retrieves canonical member genes of a Reactome, KEGG or MSigDB Hallmark pathway/gene set.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**无
- **关键词：**gene set, reactome, kegg, msigdb hallmark, pathway genes

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `pathway_id_or_name` | string | 否 |  |
| `pathway_id` | string | 否 |  |
| `pathway_name` | string | 否 |  |
| `source` | string | 否 | 可选值：`reactome`, `kegg`, `msigdb_hallmark` |
| `species` | string | 否 |  |
| `reactome_api_base` | string | 否 | 默认 `https://reactome.org/ContentService` |
| `kegg_api_base` | string | 否 | 默认 `https://rest.kegg.jp` |
| `msigdb_gmt` | string | 否 | Local filesystem path. |
| `msigdb_gmt_url` | string | 否 |  |
| `msigdb_version` | string | 否 |  |
| `limit` | integer | 否 | 默认 `500` |
| `semantic_context` | object | 否 | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "上游说明（英文原文）"

    Retrieve canonical member genes for a controlled pathway or gene-set source entity.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'pathway.retrieve_members', 'params': {}},
)
```

## `pharmacogenomics.describe_gene_requirements` {#pharmacogenomics.describe_gene_requirements}

Returns pharmacogene-specific sample evidence requirements for PharmCAT allele matching, outside calls, CYP2D6 SV/CNV, HLA, MT-RNR1 and G6PD.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**pharmacogene, star allele, cyp2d6, hla typing, g6pd

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 否 |  |
| `semantic_context` | object | 否 | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "上游说明（英文原文）"

    Return pharmacogene-specific sample evidence requirements for PharmCAT named allele matching, outside calls, CYP2D6 SV/CNV handling, HLA typing, MT-RNR1, and G6PD chrX representation.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'pharmacogenomics.describe_gene_requirements', 'params': {}},
)
```

## `pharmacogenomics.fetch_clinpgx` {#pharmacogenomics.fetch_clinpgx}

Fetches ClinPGx (PharmGKB) pharmacogenomic guidelines, clinical annotations and label evidence for a drug, gene or rsID.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**clinpgx, pharmgkb, cpic, dpwg, rsid, drug response

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `drug` | string | 否 |  |
| `gene` | string | 否 |  |
| `rsid` | string | 否 |  |
| `chemical_id` | string | 否 |  |
| `gene_id` | string | 否 |  |
| `variant_id` | string | 否 |  |
| `guideline_source` | string | 否 | 可选值：`all`, `cpic`, `dpwg`, `pro` 默认 `all` |
| `include_clinical_annotations` | boolean | 否 | 默认 `True` |
| `include_labels` | boolean | 否 | 默认 `True` |
| `include_raw_records` | boolean | 否 | 默认 `False` |
| `limit` | integer | 否 | 默认 `10` |
| `api_url` | string | 否 | 默认 `https://api.pharmgkb.org/v1` |
| `semantic_context` | object | 否 | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "上游说明（英文原文）"

    Fetch traceable ClinPGx pharmacogenomic guideline, clinical annotation, and FDA label evidence for a selected drug, gene, or rsID; compact normalized records are returned by default.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'pharmacogenomics.fetch_clinpgx', 'params': {}},
)
```

## `pharmacogenomics.fetch_pgxdb` {#pharmacogenomics.fetch_pgxdb}

Fetches PGxDB pharmacogenomic evidence for a drug, ATC code, DrugBank ID, rsID, variant marker or gene.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**pgxdb, rsid, variant marker, drug response, drugbank

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `drug` | string | 否 |  |
| `atc_code` | string | 否 |  |
| `drugbank_id` | string | 否 |  |
| `rsid` | string | 否 |  |
| `gene` | string | 否 |  |
| `variant_marker` | string | 否 |  |
| `include_raw_records` | boolean | 否 | 默认 `False` |
| `limit` | integer | 否 | 默认 `25` |
| `api_url` | string | 否 | 默认 `https://pgx-db.org/rest-api` |
| `semantic_context` | object | 否 | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "上游说明（英文原文）"

    Fetch targeted PGxDB pharmacogenomic evidence for a selected drug, ATC code, DrugBank ID, rsID, variant marker, or gene.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'pharmacogenomics.fetch_pgxdb', 'params': {}},
)
```

## `pharmacogenomics.import_pharmcat_artifacts` {#pharmacogenomics.import_pharmcat_artifacts}

Imports existing PharmCAT report, calls and match artifacts (diplotypes, phenotypes, recommendations) without running PharmCAT.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**pharmcat, diplotype, star allele, metabolizer phenotype, pgx report

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `output_dir` | string | 否 | Local filesystem path. |
| `base_filename` | string | 否 |  |
| `report_json` | string | 否 | Local filesystem path. |
| `calls_only_tsv` | string | 否 | Local filesystem path. |
| `match_json` | string | 否 | Local filesystem path. |
| `phenotype_json` | string | 否 | Local filesystem path. |
| `missing_pgx_positions_vcf` | string | 否 | Local filesystem path. |

??? note "上游说明（英文原文）"

    Import existing PharmCAT artifacts without executing PharmCAT.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'pharmacogenomics.import_pharmcat_artifacts', 'params': {}},
)
```

## `pharmacogenomics.preflight_pharmcat` {#pharmacogenomics.preflight_pharmcat}

Inspects the selected genome index (sample columns, format) for PharmCAT PGx calling readiness without running PharmCAT.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**pharmcat, preflight, vcf format, pgx positions, sample columns

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `agi_path` | string | 否 | Local filesystem path. |

??? note "上游说明（英文原文）"

    Inspect selected Active Genome Index structure for broad PharmCAT PGx calling without running PharmCAT or writing artifacts.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'pharmacogenomics.preflight_pharmcat', 'params': {}},
)
```

## `pharmacogenomics.prepare_outside_call_tsv` {#pharmacogenomics.prepare_outside_call_tsv}

Prepares a PharmCAT outside-call TSV from specialized caller output such as OptiType HLA or StellarPGx CYP2D6 diplotypes.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**outside call, optitype, stellarpgx, cyp2d6, hla, diplotype

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `caller_output_file` | string | 否 | Local filesystem path. |
| `caller_format` | string | 否 | 可选值：`auto`, `pharmcat_tsv`, `generic_table`, `optitype`, `stellarpgx_summary` 默认 `auto` |
| `output_file` | string | 否 | Local filesystem path. |
| `sample` | string | 否 |  |
| `max_rows` | integer | 否 | 默认 `200` |

??? note "上游说明（英文原文）"

    Prepare a PharmCAT outside-call TSV from supported specialized caller output such as OptiType HLA calls, StellarPGx CYP2D6 summaries, or generic gene/diplotype tables.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'pharmacogenomics.prepare_outside_call_tsv', 'params': {}},
)
```

## `pharmacogenomics.review_medication` {#pharmacogenomics.review_medication}

Reviews pharmacogenomic evidence for a medication, drug-gene or drug-variant target, using sample and stored evidence.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**medication review, drug-gene, rsid, pgx, dosing

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `db` | string | 否 | Local filesystem path. |
| `shared_db` | string | 否 | Local filesystem path. |
| `drug` | string | 否 |  |
| `atc_code` | string | 否 |  |
| `drugbank_id` | string | 否 |  |
| `rsid` | string | 否 |  |
| `gene` | string | 否 |  |
| `indication` | string | 否 |  |
| `dose` | string | 否 |  |
| `current_medications` | string | 否 |  |
| `allergies_or_contraindications` | string | 否 |  |
| `known_genotype` | string | 否 |  |
| `known_diplotype` | string | 否 |  |
| `known_phenotype` | string | 否 |  |
| `known_activity_score` | string | 否 |  |
| `known_pgx_source` | string | 否 |  |
| `source_sample_pgx_row_id` | string | 否 |  |
| `semantic_context` | object | 否 | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |
| `genome_build` | string | 否 | 默认 `GRCh38` |
| `include_active_genome_index` | boolean | 否 |  |
| `include_known_active_genome_indexes` | boolean | 否 | 默认 `False` |
| `include_stored_research` | boolean | 否 | 默认 `True` |
| `include_record_research_payloads` | boolean | 否 | 默认 `False` |
| `limit` | integer | 否 | 默认 `10` |
| `clinpgx_api_url` | string | 否 |  |
| `pgxdb_api_url` | string | 否 |  |
| `fda_biomarkers_url` | string | 否 |  |
| `fda_associations_url` | string | 否 |  |

??? note "上游说明（英文原文）"

    Review pharmacogenomic evidence for a selected medication, drug-gene, or drug-variant target.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'pharmacogenomics.review_medication', 'params': {}},
)
```

## `pharmacogenomics.run_pharmcat` {#pharmacogenomics.run_pharmcat}

Runs a local PharmCAT installation on an approved genome index to call pharmacogene diplotypes and phenotypes.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**pharmcat, star allele, diplotype, pgx calling, metabolizer

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `agi_path` | string | 否 | Local filesystem path. |
| `output_dir` | string | 否 | Local filesystem path. |
| `base_filename` | string | 否 |  |
| `mode` | string | 否 | 可选值：`auto`, `pipeline`, `jar` 默认 `auto` |
| `pipeline_command` | string | 否 | Local filesystem path. |
| `pharmcat_jar` | string | 否 | Local filesystem path. |
| `java_command` | string | 否 | Local filesystem path. |
| `sample` | string | 否 |  |
| `sample_file` | string | 否 | Local filesystem path. |
| `sample_metadata` | string | 否 | Local filesystem path. |
| `outside_call_file` | string | 否 | Local filesystem path. |
| `reporter_sources` | string | 否 |  |
| `research_mode` | string | 否 |  |
| `max_memory` | string | 否 |  |
| `max_concurrent_processes` | integer | 否 |  |
| `include_reporter_json` | boolean | 否 | 默认 `True` |
| `include_calls_only_tsv` | boolean | 否 | 默认 `True` |
| `probe_version` | boolean | 否 | 默认 `True` |
| `dry_run` | boolean | 否 | 默认 `False` |
| `timeout_seconds` | integer | 否 | 默认 `7200` |

??? note "上游说明（英文原文）"

    Run a local PharmCAT installation from an approved Active Genome Index.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'pharmacogenomics.run_pharmcat', 'params': {}},
)
```

## `pharmacogenomics.validate_outside_call_tsv` {#pharmacogenomics.validate_outside_call_tsv}

Validates PharmCAT outside-call TSV structure and summarizes diplotype, phenotype or activity-score evidence.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**outside call, diplotype, activity score, validation, pharmcat

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `outside_call_file` | string | 否 | Local filesystem path. |
| `max_rows` | integer | 否 | 默认 `200` |

??? note "上游说明（英文原文）"

    Validate PharmCAT outside-call TSV structure and summarize selected diplotype, phenotype, or activity-score evidence.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'pharmacogenomics.validate_outside_call_tsv', 'params': {}},
)
```

## `phenotype.compare_disease_evidence` {#phenotype.compare_disease_evidence}

Compares supplied or gene-derived candidate diseases against phenotype/HPO evidence without selecting a diagnosis.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**hpo, rare disease, differential, phenotype matching, gene-disease

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `phenotypes` | array of string | 否 |  |
| `terms` | array of string | 否 |  |
| `hpo_ids` | array of string | 否 |  |
| `search_stored_research` | boolean | 否 | 默认 `True` |
| `use_hpo_annotations` | boolean | 否 | 默认 `True` |
| `source_records` | array of object | 否 |  |
| `db` | string | 否 | Local filesystem path. |
| `phenotype_text` | string | 否 |  |
| `text` | string | 否 |  |
| `candidate_diseases` | array of string | 否 |  |
| `diseases` | array of string | 否 |  |
| `genes` | array of string | 否 |  |
| `gene_symbols` | array of string | 否 |  |
| `hpo_disease_file` | string | 否 | Local filesystem path. |
| `use_primary_gene_disease` | boolean | 否 | 默认 `True` |
| `gencc_file` | string | 否 | Local filesystem path. |
| `limit` | integer | 否 | 默认 `25` |
| `semantic_context` | object | 否 | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "上游说明（英文原文）"

    Compare supplied or primary gene-derived diseases against phenotype/HPO evidence without selecting the diagnosis.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'phenotype.compare_disease_evidence', 'params': {}},
)
```

## `phenotype.compare_drug_target_evidence` {#phenotype.compare_drug_target_evidence}

Compares candidate genes using direct drug-target or mechanism-of-action evidence for a drug, class or indication.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`genes`
- **关键词：**drug target, mechanism of action, candidate genes, target evidence, indication

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `source_records` | array of object | 否 |  |
| `db` | string | 否 | Local filesystem path. |
| `drug` | string | 否 |  |
| `drug_class` | string | 否 |  |
| `indication` | string | 否 |  |
| `phenotype` | string | 否 |  |
| `mechanism` | string | 否 |  |
| `genes` | array of string | 是 |  |
| `search_stored_research` | boolean | 否 | 默认 `True` |
| `limit` | integer | 否 | 默认 `25` |
| `semantic_context` | object | 否 | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "上游说明（英文原文）"

    Compare candidate genes using direct drug-target or mechanism evidence only.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'phenotype.compare_drug_target_evidence',
     'params': {'genes': []}},
)
```

## `phenotype.compare_gene_hpo_evidence` {#phenotype.compare_gene_hpo_evidence}

Compares candidate genes using phenotype, HPO and curated rare-disease annotation evidence.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`genes`
- **关键词：**hpo, gene-phenotype, rare disease, candidate genes, annotation

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `phenotypes` | array of string | 否 |  |
| `terms` | array of string | 否 |  |
| `hpo_ids` | array of string | 否 |  |
| `search_stored_research` | boolean | 否 | 默认 `True` |
| `use_hpo_annotations` | boolean | 否 | 默认 `True` |
| `hpo_gene_file` | string | 否 | Local filesystem path. |
| `source_records` | array of object | 否 |  |
| `db` | string | 否 | Local filesystem path. |
| `phenotype` | string | 否 |  |
| `phenotype_text` | string | 否 |  |
| `text` | string | 否 |  |
| `task_text` | string | 否 |  |
| `question` | string | 否 |  |
| `prompt` | string | 否 |  |
| `condition` | string | 否 |  |
| `genes` | array of string | 是 |  |
| `limit` | integer | 否 | 默认 `25` |
| `semantic_context` | object | 否 | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "上游说明（英文原文）"

    Compare candidate genes using phenotype, HPO, and curated rare-disease annotation evidence only.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'phenotype.compare_gene_hpo_evidence', 'params': {'genes': []}},
)
```

## `phenotype.plan_risk_investigation` {#phenotype.plan_risk_investigation}

Plans phenotype risk and ClinVar review targets (rare disease, cancer risk, carrier) from genes, conditions and optional genome evidence.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**risk investigation, carrier screening, cancer risk, clinvar, candidate genes

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `question` | string | 否 |  |
| `investigation_type` | string | 否 | 可选值：`auto`, `rare_disease`, `cancer_risk`, `carrier_review`, `observed_condition_review` 默认 `auto` |
| `gene` | string | 否 |  |
| `genes` | array of string | 否 |  |
| `condition` | string | 否 |  |
| `topic` | string | 否 |  |
| `db` | string | 否 | Local filesystem path. |
| `matches` | string | 否 | Local filesystem path. |
| `genome_build` | string | 否 | 默认 `GRCh38` |
| `include_active_genome_index` | boolean | 否 |  |
| `search_stored_research` | boolean | 否 | 默认 `True` |
| `limit` | integer | 否 | 默认 `25` |
| `semantic_context` | object | 否 | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "上游说明（英文原文）"

    Plan phenotype risk and ClinVar review targets from public targets and optionally selected Active Genome Index evidence.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'phenotype.plan_risk_investigation', 'params': {}},
)
```

## `phenotype.retrieve_disease_drug_targets` {#phenotype.retrieve_disease_drug_targets}

Retrieves disease-scoped clinical drug-target genes from Open Targets drug candidate records.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**open targets, drug target, disease, target genes, clinical stage

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `disease` | string | 否 |  |
| `disease_id` | string | 否 |  |
| `genes` | array of string | 否 |  |
| `mode` | string | 否 | 可选值：`records`, `gene_membership` 默认 `records` |
| `minimum_clinical_stage` | string | 否 | 默认 `PHASE_2` |
| `opentargets_api_url` | string | 否 | 默认 `https://api.platform.opentargets.org/api/v4/graphql` |
| `limit` | integer | 否 | 默认 `100` |
| `semantic_context` | object | 否 | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "上游说明（英文原文）"

    Retrieve disease-scoped clinical drug-target genes from Open Targets drug candidate records.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'phenotype.retrieve_disease_drug_targets', 'params': {}},
)
```

## `phenotype.retrieve_gene_disease_associations` {#phenotype.retrieve_gene_disease_associations}

Retrieves primary gene-disease validity associations from GenCC for supplied gene symbols.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**gencc, gene-disease validity, mendelian, gene symbol, disease

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `genes` | array of string | 否 |  |
| `gene_symbols` | array of string | 否 |  |
| `semantic_context` | object | 否 | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |
| `gencc_file` | string | 否 | Local filesystem path. |
| `classifications` | array of string | 否 |  |
| `limit` | integer | 否 | 默认 `100` |

??? note "上游说明（英文原文）"

    Retrieve primary gene-disease associations from GenCC for supplied gene symbols.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'phenotype.retrieve_gene_disease_associations', 'params': {}},
)
```

## `phenotype.retrieve_trait_gene_records` {#phenotype.retrieve_trait_gene_records}

Retrieves trait-to-gene records from integrated public sources such as Open Targets, optionally filtered to given genes.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`trait`
- **关键词：**trait, gene association, open targets, genetic evidence, phenotype

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `trait` | string | 是 |  |
| `genes` | array of string | 否 |  |
| `opentargets_api_url` | string | 否 | 默认 `https://api.platform.opentargets.org/api/v4/graphql` |
| `limit` | integer | 否 | 默认 `25` |
| `semantic_context` | object | 否 | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "上游说明（英文原文）"

    Retrieve native trait-to-gene records from integrated public sources, optionally filtered to gene symbols.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'phenotype.retrieve_trait_gene_records',
     'params': {'trait': '...'}},
)
```

## `prs.calculate_score` {#prs.calculate_score}

Applies a published polygenic score to an approved genome index and returns the raw weighted score plus QC.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**polygenic risk score, pgs, weighted score, genotype dosage, personal genome

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `agi_path` | string | 否 | Local filesystem path. |
| `pgs_id` | string | 否 |  |
| `score_dir` | string | 否 | Local filesystem path. |
| `genome_build` | string | 否 | 默认 `GRCh38` |
| `skip_ambiguous_palindromic` | boolean | 否 | 默认 `True` |
| `score_mean` | number | 否 |  |
| `score_sd` | number | 否 |  |

??? note "上游说明（英文原文）"

    Apply a published polygenic score to an approved Active Genome Index and return raw weighted score plus QC.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'prs.calculate_score', 'params': {}},
)
```

## `prs.check_score_overlap` {#prs.check_score_overlap}

Checks how many variants from a polygenic score are present and usable in an approved genome index.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**pgs, variant overlap, match rate, polygenic score, coverage

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `agi_path` | string | 否 | Local filesystem path. |
| `pgs_id` | string | 否 |  |
| `score_dir` | string | 否 | Local filesystem path. |
| `genome_build` | string | 否 | 默认 `GRCh38` |
| `skip_ambiguous_palindromic` | boolean | 否 | 默认 `True` |

??? note "上游说明（英文原文）"

    Check how many variants from a polygenic score are usable in an approved Active Genome Index.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'prs.check_score_overlap', 'params': {}},
)
```

## `prs.fetch_score_metadata` {#prs.fetch_score_metadata}

Fetches detailed public PGS Catalog metadata for one score ID, including scoring-file URLs and publication context.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`pgs_id`
- **关键词：**pgs catalog, pgs id, scoring file, polygenic score, metadata

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `pgs_id` | string | 是 |  |

??? note "上游说明（英文原文）"

    Fetch detailed public PGS Catalog metadata for one score ID, including scoring-file URLs and source publication context.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'prs.fetch_score_metadata', 'params': {'pgs_id': '...'}},
)
```

## `prs.import_scoring_file` {#prs.import_scoring_file}

Imports a PGS Catalog or local scoring file (variant weights) into the local PRS score cache for a genome build.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**scoring file, variant weights, pgs catalog, import, genome build

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `pgs_id` | string | 否 |  |
| `genome_build` | string | 否 | 默认 `GRCh38` |
| `scoring_file` | string | 否 | Local filesystem path. |
| `scoring_url` | string | 否 |  |
| `force` | boolean | 否 | 默认 `False` |

??? note "上游说明（英文原文）"

    Import a PGS Catalog or local scoring file into Genomi's local PRS score cache for a declared genome build.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'prs.import_scoring_file', 'params': {}},
)
```

## `prs.search_scores` {#prs.search_scores}

Searches public PGS Catalog score metadata by trait, score ID, EFO term or free text.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**pgs catalog, polygenic score, trait, efo, search

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 否 |  |
| `trait` | string | 否 |  |
| `pgs_id` | string | 否 |  |
| `efo_id` | string | 否 |  |
| `semantic_context` | object | 否 | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |
| `limit` | integer | 否 | 默认 `20` |

??? note "上游说明（英文原文）"

    Search public PGS Catalog score metadata by trait, score ID, EFO term, or free-text query without reading Active Genome Index.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'prs.search_scores', 'params': {}},
)
```

## `region.retrieve_features` {#region.retrieve_features}

Retrieves GENCODE gene/transcript and ENCODE cCRE feature annotations overlapping a genomic region.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**无
- **关键词：**gencode, encode ccre, genomic region, gene features, grch38

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chrom` | string | 否 |  |
| `start` | integer | 否 |  |
| `end` | integer | 否 |  |
| `region` | string | 否 |  |
| `assembly` | string | 否 | 可选值：`GRCh37`, `GRCh38` |
| `gencode_gtf` | string | 否 | Local filesystem path. |
| `encode_ccre_bed` | string | 否 | Local filesystem path. |
| `limit` | integer | 否 | 默认 `100` |

??? note "上游说明（英文原文）"

    Retrieve genomic-region feature annotations from supplied or installed GENCODE and ENCODE annotation files.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'region.retrieve_features', 'params': {}},
)
```

## `research.build_target_packet` {#research.build_target_packet}

Builds a target-centric evidence packet for a variant, gene, drug, condition or topic target.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`target_type`
- **关键词：**evidence packet, variant evidence, gene evidence, chrom pos ref alt, target

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chrom` | string | 否 |  |
| `pos` | integer | 否 |  |
| `ref` | string | 否 |  |
| `alt` | string | 否 |  |
| `gene` | string | 否 |  |
| `drug` | string | 否 |  |
| `condition` | string | 否 |  |
| `topic` | string | 否 |  |
| `db` | string | 否 | Local filesystem path. |
| `target_type` | string | 是 | 可选值：`condition`, `drug`, `gene`, `topic`, `variant` |
| `genome_build` | string | 否 | 默认 `GRCh38` |
| `source_id` | string | 否 |  |
| `limit` | integer | 否 | 默认 `20` |

??? note "上游说明（英文原文）"

    Build a target-centric evidence packet after the agent identifies the user's target.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'research.build_target_packet',
    {'target_type': 'condition'},
)
```

## `sequence.analyze` {#sequence.analyze}

Runs a compact DNA sequence analysis bundle (summary, translation, ORFs, restriction sites, Kozak context).

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**gc content, orf, translation, restriction sites, kozak

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 |  |
| `mode` | string | 否 | 可选值：`summary`, `translate`, `orfs`, `restriction_sites`, `kozak` 默认 `summary` |
| `frame` | integer | 否 | 默认 `1` |
| `strand` | string | 否 | 可选值：`forward`, `reverse`, `both` |
| `min_aa` | integer | 否 | 默认 `30` |
| `enzymes` | array of string | 否 |  |
| `motifs` | array of string | 否 |  |
| `start_pos` | integer | 否 |  |
| `reference_fasta` | string | 否 | Local filesystem path. |
| `max_matches` | integer | 否 | 默认 `10` |

??? note "上游说明（英文原文）"

    Run a compact deterministic sequence analysis bundle and point to focused sequence tools when needed.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'sequence.analyze', 'params': {'sequence': '...'}},
)
```

## `sequence.check_primers` {#sequence.check_primers}

Checks basic primer properties (e.g. length, GC, Tm) and optional amplicons on a template sequence.

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**`forward_primer`
- **关键词：**pcr primers, melting temperature, amplicon, gc content, primer design

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `forward_primer` | string | 是 |  |
| `reverse_primer` | string | 否 |  |
| `template` | string | 否 |  |

??? note "上游说明（英文原文）"

    Check basic primer properties and optional template amplicons.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'sequence.check_primers', 'params': {'forward_primer': '...'}},
)
```

## `sequence.classify_kozak` {#sequence.classify_kozak}

Classifies Kozak sequence context around ATG start codons in a supplied sequence.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**kozak, start codon, atg, translation initiation

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 |  |
| `start_pos` | integer | 否 |  |

??? note "上游说明（英文原文）"

    Classify Kozak sequence context around ATG start codons.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'sequence.classify_kozak', 'params': {'sequence': '...'}},
)
```

## `sequence.find_orfs` {#sequence.find_orfs}

Finds ATG-to-stop open reading frames above a minimum length on one or both strands of a DNA sequence.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**orf, open reading frame, start codon, stop codon, strand

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 |  |
| `min_aa` | integer | 否 | 默认 `30` |
| `strand` | string | 否 | 可选值：`forward`, `reverse`, `both` 默认 `both` |

??? note "上游说明（英文原文）"

    Find ATG-to-stop open reading frames in a supplied DNA sequence.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'sequence.find_orfs', 'params': {'sequence': '...'}},
)
```

## `sequence.find_restriction_sites` {#sequence.find_restriction_sites}

Finds common restriction enzyme or custom motif sites in a supplied DNA sequence.

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**`sequence`
- **关键词：**restriction enzyme, cut sites, motif, cloning, digest

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 |  |
| `enzymes` | array of string | 否 |  |
| `motifs` | array of string | 否 |  |

??? note "上游说明（英文原文）"

    Find common restriction enzyme or custom motif sites in a supplied DNA sequence.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'sequence.find_restriction_sites',
     'params': {'sequence': '...'}},
)
```

## `sequence.match_reference` {#sequence.match_reference}

Matches a supplied DNA sequence against local FASTA records and returns record identifiers and annotations.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`, `reference_fasta`
- **关键词：**fasta, sequence matching, reference, alignment, identification

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 |  |
| `reference_fasta` | string | 是 | Local filesystem path. |
| `max_matches` | integer | 否 | 默认 `10` |

??? note "上游说明（英文原文）"

    Match a supplied DNA sequence against local FASTA records and return record identifiers plus annotations.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'sequence.match_reference',
     'params': {'sequence': '...', 'reference_fasta': '...'}},
)
```

## `sequence.translate` {#sequence.translate}

Translates a DNA sequence in a selected frame and strand using the standard genetic code.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**translation, genetic code, reading frame, codon, protein sequence

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 |  |
| `frame` | integer | 否 | 默认 `1` |
| `strand` | string | 否 | 可选值：`forward`, `reverse` 默认 `forward` |

??? note "上游说明（英文原文）"

    Translate a DNA sequence in a selected frame and strand using the standard genetic code.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'sequence.translate', 'params': {'sequence': '...'}},
)
```

## `variant.find_gene_variants` {#variant.find_gene_variants}

Finds passing variants in an approved genome index that overlap GENCODE intervals of up to 10 candidate genes.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`genes`
- **关键词：**gene variants, gencode, personal genome, candidate genes, vcf

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `genes` | array of string | 是 |  |
| `genome_build` | string | 否 | 可选值：`GRCh37`, `GRCh38` 默认 `GRCh38` |
| `agi_id` | string | 否 |  |
| `per_gene_limit` | integer | 否 | 默认 `100` |

??? note "上游说明（英文原文）"

    Find passing variants in an approved Active Genome Index that overlap managed GENCODE intervals for a bounded candidate-gene set.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'variant.find_gene_variants', 'params': {'genes': []}},
)
```

## `variant.gather_allele_context` {#variant.gather_allele_context}

Gathers existing sample, static, population-frequency and reviewed research evidence for one allele (chrom/pos/ref/alt).

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`chrom`, `pos`, `ref`, `alt`
- **关键词：**allele, population frequency, evidence, chrom pos ref alt, genotype

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chrom` | string | 是 |  |
| `pos` | integer | 是 |  |
| `ref` | string | 是 |  |
| `alt` | string | 是 |  |
| `db` | string | 否 | Local filesystem path. |
| `matches` | string | 否 | Local filesystem path. |
| `genome_build` | string | 否 | 默认 `GRCh38` |
| `population_source` | string | 否 |  |
| `population` | string | 否 |  |

??? note "上游说明（英文原文）"

    Gather existing sample, static, population, and reviewed research evidence for one allele.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'variant.gather_allele_context',
     'params': {'chrom': '...', 'pos': 0, 'ref': '...', 'alt': '...'}},
)
```

## `variant.gather_gene_context` {#variant.gather_gene_context}

Gathers existing sample, ClinVar and reviewed research evidence for one gene.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`gene`
- **关键词：**gene, clinvar, sample variants, evidence, personal genome

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `db` | string | 否 | Local filesystem path. |
| `gene` | string | 是 |  |
| `matches` | string | 否 | Local filesystem path. |
| `genome_build` | string | 否 | 默认 `GRCh38` |

??? note "上游说明（英文原文）"

    Gather existing sample, ClinVar, and reviewed research evidence for one gene.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'variant.gather_gene_context', 'params': {'gene': '...'}},
)
```

## `variant.resolve` {#variant.resolve}

Resolves one variant target (coordinates, rsID, query or region) and returns public, local sample and stored evidence facts.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**rsid, variant lookup, hgvs, genotype, clinvar

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `db` | string | 否 | Local filesystem path. |
| `shared_db` | string | 否 | Local filesystem path. |
| `chrom` | string | 否 |  |
| `pos` | integer | 否 |  |
| `ref` | string | 否 |  |
| `alt` | string | 否 |  |
| `query` | string | 否 |  |
| `rsid` | string | 否 |  |
| `region` | string | 否 |  |
| `genome_build` | string | 否 | 默认 `GRCh38` |
| `agi_id` | string | 否 |  |
| `include_active_genome_index` | boolean | 否 |  |
| `include_known_active_genome_indexes` | boolean | 否 | 默认 `False` |
| `include_fail` | boolean | 否 | 默认 `False` |
| `limit` | integer | 否 | 默认 `20` |

??? note "上游说明（英文原文）"

    Resolve one variant target and return deterministic public, local sample, and stored evidence facts.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'variant.resolve', 'params': {}},
)
```
