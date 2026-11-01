# Genomi (`genomi`)

DNAKit's MCP tool index keeps **57** DNA/gene tools of this server (112 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [exon-research/genomi @ `1df4f5b`](https://github.com/exon-research/genomi/tree/1df4f5bfb635e2d08f93a1102b70b840fbe80d86) |
| License | Apache-2.0 |
| Command | `genomi serve` |
| Transports | stdio |
| Installation | Install the exon-research/genomi source project (not on PyPI; Python &gt;=3.10, depends on a git-pinned gxl-paperclip) so that the genomi console script exists, then run genomi install to materialize reference libraries; optionally set GENOMI_HOME to relocate its data directory (default ~/.genomi). |
| Effects | Parses user-supplied genome files into a local Active Genome Index and evidence database under GENOMI_HOME, keeps a local journal, runs background jobs, and genomi.install downloads public reference libraries; returned results can contain individual genotypes, and optional lab operations may contact external evidence providers. |

| Tool | Category | Summary |
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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** callset qc, depth, genotype quality, gvcf, absence claims

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `db` | string | no | Local filesystem path. |
| `agi_path` | string | no | Local filesystem path. |
| `output` | string | no | Local filesystem path. |
| `genome_build` | string | no | Default `auto` |
| `scan_records` | integer | no | Default `1000` |

??? note "Upstream description"

    Classify genome callset shape, depth/quality field availability, and absence-claim boundaries using an Active Genome Index.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `chrom`, `pos`, `ref`, `alt`
- **Keywords:** genotype support, allele, read depth, personal genome, zygosity

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chrom` | string | yes |  |
| `pos` | integer | yes |  |
| `ref` | string | yes |  |
| `alt` | string | yes |  |
| `db` | string | no | Local filesystem path. |
| `agi_path` | string | no | Local filesystem path. |
| `output` | string | no | Local filesystem path. |
| `genome_build` | string | no | Default `auto` |
| `reference_fasta` | string | no | Local filesystem path. |
| `min_depth` | integer | no | Default `10` |
| `min_gq` | integer | no | Default `20` |

??? note "Upstream description"

    Classify whether one exact allele has enough sample support to be used in a personal interpretation.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `region`
- **Keywords:** callability, coverage, reference call, absence, region

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `region` | string | yes |  |
| `db` | string | no | Local filesystem path. |
| `agi_path` | string | no | Local filesystem path. |
| `output` | string | no | Local filesystem path. |
| `genome_build` | string | no | Default `auto` |
| `min_depth` | integer | no | Default `10` |
| `min_covered_fraction` | number | no | Default `0.95` |
| `limit` | integer | no | Default `5000` |

??? note "Upstream description"

    Classify whether a region can support reference or absence claims.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** ancestry markers, marker overlap, 1000 genomes, snp, personal genome

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `agi_path` | string | no | Local filesystem path. |
| `genome_build` | string | no | Default `GRCh38` |

??? note "Upstream description"

    Check how many installed 1000 Genomes ancestry panel markers are usable in an approved Active Genome Index.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** genetic ancestry, population structure, pca, 1000 genomes, superpopulation

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `agi_path` | string | no | Local filesystem path. |
| `genome_build` | string | no | Default `GRCh38` |
| `nearest_reference_count` | integer | no | Default `10` |

??? note "Upstream description"

    Estimate qualitative reference-panel similarity for an approved GRCh37 or GRCh38 sample using local 1000 Genomes PCA projection.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** pca projection, genetic ancestry, nearest neighbors, 1000 genomes, genotypes

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `agi_path` | string | no | Local filesystem path. |
| `genome_build` | string | no | Default `GRCh38` |
| `nearest_reference_count` | integer | no | Default `10` |

??? note "Upstream description"

    Project an approved sample into the installed 1000 Genomes ancestry PCA space and return nearest reference neighbors.

Call it through DNAKit once the server is installed:

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

- **Category:** `expression` (Gene-level expression and transcript abundance)
- **Required arguments:** none
- **Keywords:** marker genes, cell type, cellmarker, panglaodb, human protein atlas

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `cell_type_id_or_name` | string | no |  |
| `cell_type_id` | string | no |  |
| `cell_type_name` | string | no |  |
| `source` | string | no | One of `hpa`, `cellmarker`, `panglaodb`, `encode` |
| `species` | string | no |  |
| `marker_table` | string | no | Local filesystem path. |
| `hpa_api_base` | string | no | Default `https://www.proteinatlas.org/api` |
| `limit` | integer | no | Default `100` |
| `semantic_context` | object | no | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "Upstream description"

    Retrieve canonical marker genes for a controlled cell-type source entity.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** clinvar, pathogenicity, variant matching, clinical significance, personal genome

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `agi_path` | string | no | Local filesystem path. |
| `db` | string | no | Local filesystem path. |
| `output` | string | no | Local filesystem path. |
| `genome_build` | string | no | Default `auto` |
| `force` | boolean | no | Default `False` |

??? note "Upstream description"

    Materialize exact ClinVar matches for comparable Active Genome Index variants using the installed build-specific ClinVar library.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** clinvar, pathogenic, carrier, review groups, variant triage

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `agi_path` | string | no | Local filesystem path. |
| `db` | string | no | Local filesystem path. |
| `output` | string | no | Local filesystem path. |
| `genome_build` | string | no | Default `GRCh38` |
| `force` | boolean | no | Default `False` |

??? note "Upstream description"

    Build deterministic ClinVar variant evidence inventory and review groups from materialized matches.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** dashboard, personal genome report, pgx, polygenic risk, ancestry

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `panels` | array of string | no | Optional subset of dashboard panels to gather. Omitted means every dashboard panel. |
| `risk_score_ids` | array of string | no | Optional PGS IDs to score for the risk panel. Omitted means use imported local scores up to risk_score_limit. |
| `risk_score_limit` | integer | no | Maximum imported local PGS scores to calculate when risk_score_ids is omitted. Default `5` |
| `risk_review_types` | array of string | no | Phenotype review modes to add to the risk panel from selected Active Genome Index ClinVar review groups. Empty array disables phenotype review rows. Default `\['carrier_review', 'observed_condition_review'\]` |
| `pgx_review_targets` | array of object | no | Optional medication-review targets to add to the PGx panel. |
| `pgx_review_target_limit` | integer | no | Maximum medication-review targets to run for the PGx panel. Default `12` |
| `nutrigenomics_domain_ids` | array of string | no | Optional nutrigenomics domains to retrieve. Omitted means every declared in-scope domain. |

??? note "Upstream description"

    Build Genomi Dashboard panel states from the approved Active Genome Index and existing capability operations. Used by decode.render_dashboard; presented output reports panel readiness, running jobs, and gaps.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** dashboard, html report, personal genome, pgx, polygenic risk

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `panels` | array of string | no | Optional subset of dashboard panels to gather. Omitted means every dashboard panel. |
| `risk_score_ids` | array of string | no | Optional PGS IDs to score for the risk panel. Omitted means use imported local scores up to risk_score_limit. |
| `risk_score_limit` | integer | no | Maximum imported local PGS scores to calculate when risk_score_ids is omitted. Default `5` |
| `risk_review_types` | array of string | no | Phenotype review modes to add to the risk panel from selected Active Genome Index ClinVar review groups. Empty array disables phenotype review rows. Default `\['carrier_review', 'observed_condition_review'\]` |
| `pgx_review_targets` | array of object | no | Optional medication-review targets to add to the PGx panel. |
| `pgx_review_target_limit` | integer | no | Maximum medication-review targets to run for the PGx panel. Default `12` |
| `nutrigenomics_domain_ids` | array of string | no | Optional nutrigenomics domains to retrieve. Omitted means every declared in-scope domain. |
| `output` | string | no | Optional override for the dashboard HTML path. Defaults to a transient path under the system temp dir: &lt;tmp&gt;/genomi-dashboards/&lt;sample&gt;/dashboard.html. |

??? note "Upstream description"

    Build, shape, render, and expose the Genomi Dashboard self-contained HTML artifact from the approved Active Genome Index.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `context`, `genes`
- **Keywords:** crispr screen, knockout, gene essentiality, perturbation, cell line

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `source_records` | array of object | no |  |
| `context` | string | yes |  |
| `genes` | array of string | yes |  |
| `organism` | string | no |  |
| `cell_line` | string | no |  |
| `perturbation` | string | no |  |
| `assay` | string | no |  |
| `phenotype` | string | no |  |
| `perturbation_sources` | array of string | no |  |
| `biogrid_orcs_access_key` | string | no |  |
| `depmap_gene_effect_url` | string | no |  |
| `depmap_model_url` | string | no |  |
| `db` | string | no | Local filesystem path. |
| `search_stored_research` | boolean | no | Default `True` |
| `retrieve_native` | boolean | no | Default `True` |
| `limit` | integer | no | Default `25` |
| `semantic_context` | object | no | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "Upstream description"

    Compare candidate genes by verified functional-genomics perturbation experiment evidence.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `table`, `context`, `genes`
- **Keywords:** perturbation screen, crispr, gene hits, csv import, phenotype

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `context` | string | yes |  |
| `genes` | array of string | yes |  |
| `organism` | string | no |  |
| `cell_line` | string | no |  |
| `perturbation` | string | no |  |
| `assay` | string | no |  |
| `phenotype` | string | no |  |
| `table` | string | yes | Local filesystem path. |
| `column_map` | object | no |  |
| `delimiter` | string | no |  |
| `source_title` | string | no |  |
| `source_url` | string | no |  |
| `source_type` | string | no |  |
| `limit` | integer | no | Default `500` |

??? note "Upstream description"

    Extract verified perturbation experiment source records from a local CSV or TSV result table.

Call it through DNAKit once the server is installed:

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

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `context`
- **Keywords:** geo, gse accession, perturbation, dataset discovery, ncbi

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `context` | string | yes |  |
| `genes` | array of string | no |  |
| `organism` | string | no |  |
| `cell_line` | string | no |  |
| `perturbation` | string | no |  |
| `assay` | string | no |  |
| `phenotype` | string | no |  |
| `accession` | string | no |  |
| `limit` | integer | no | Default `25` |
| `semantic_context` | object | no | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |
| `ncbi_api_key` | string | no |  |
| `ncbi_email` | string | no |  |
| `ncbi_tool` | string | no |  |

??? note "Upstream description"

    Discover public study tables from NCBI GEO when curated perturbation sources do not cover a requested functional-genomics context.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `context`, `genes`
- **Keywords:** biogrid orcs, depmap, crispr screen, gene effect, essentiality

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `context` | string | yes |  |
| `genes` | array of string | yes |  |
| `organism` | string | no |  |
| `cell_line` | string | no |  |
| `perturbation` | string | no |  |
| `assay` | string | no |  |
| `phenotype` | string | no |  |
| `perturbation_sources` | array of string | no |  |
| `biogrid_orcs_access_key` | string | no |  |
| `depmap_gene_effect_url` | string | no |  |
| `depmap_model_url` | string | no |  |
| `limit` | integer | no | Default `100` |
| `semantic_context` | object | no | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "Upstream description"

    Retrieve native public functional-genomics perturbation records from BioGRID ORCS and DepMap for candidate genes and declared experimental context.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `source`
- **Keywords:** personal genome, vcf, gvcf, bam, fastq, consumer genotype, genome build

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `db` | string | no | Local filesystem path. |
| `shared_db` | string | no | Local filesystem path. |
| `source` | string | yes | Local filesystem path. |
| `source_evidence_db` | string | no | Local filesystem path. |
| `reference_fasta` | string | no | Local filesystem path. |
| `auto_reference_fasta` | boolean | no | Default `True` |
| `genome_build` | string | no | Default `auto` |
| `max_records` | integer | no |  |
| `parallel_workers` | integer | no |  |
| `force` | boolean | no | Default `False` |
| `user_nickname` | string | no |  |
| `set_default_user` | boolean | no | Default `False` |

??? note "Upstream description"

    Detect, parse, and digitize a genome source such as VCF/gVCF, BAM, paired FASTQ, supported consumer exports, or a genome.computer .genome/1.0 bundle. After parsing, read skills/active-genome-index/SKILL.md before using the invoke-only active_genome_index.* interpretation/selection tools; the result's next_actions repeats this.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `investigation_id`, `round_id`, `deduplication_key`, `gene`, `transcript_accession`, `protein_accession`, `coding_change`, `protein_substitution`, `public_reference_protein_sequence`, `reference_source_label`, `reference_source_version`, `reference_source_record_id`
- **Keywords:** hgvs, coding change, missense, protein substitution, transcript accession

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `investigation_id` | string | yes |  |
| `round_id` | string | yes |  |
| `deduplication_key` | string | yes |  |
| `gene` | string | yes |  |
| `transcript_accession` | string | yes |  |
| `protein_accession` | string | yes |  |
| `coding_change` | string | yes |  |
| `protein_substitution` | string | yes |  |
| `public_reference_protein_sequence` | string | yes |  |
| `reference_source_label` | string | yes |  |
| `reference_source_version` | string | yes |  |
| `reference_source_record_id` | string | yes |  |

??? note "Upstream description"

    Verify one intended protein substitution against a supplied public reference protein sequence using deterministic local Genomi rules. Persist only sequence digests and normalized descriptors as a round-bound nonclinical research artifact; never store the sequence or treat the result as clinical evidence.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `chrom`, `pos`, `ref`, `alt`
- **Keywords:** gnomad, allele frequency, population frequency, af, chrom pos ref alt

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `db` | string | no | Local filesystem path. |
| `shared_db` | string | no | Local filesystem path. |
| `chrom` | string | yes |  |
| `pos` | integer | yes |  |
| `ref` | string | yes |  |
| `alt` | string | yes |  |
| `sync_shared` | boolean | no | Default `True` |
| `dataset` | string | no | Default `gnomad_r4` |
| `genome_build` | string | no | Default `GRCh38` |
| `api_url` | string | no | Default `https://gnomad.broadinstitute.org/api` |

??? note "Upstream description"

    Fetch reusable gnomAD public population frequency for one allele and write it into evidence storage.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `phenotype`, `genes`
- **Keywords:** gwas catalog, mapped gene, trait association, phenotype, candidate genes

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `phenotype` | string | yes |  |
| `genes` | array of string | yes |  |
| `association_limit` | integer | no | Default `200` |
| `api_url` | string | no | Default `https://www.ebi.ac.uk/gwas/rest/api/v2` |
| `semantic_context` | object | no | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "Upstream description"

    Compare candidate genes using GWAS Catalog reported_gene and mapped_gene trait-association evidence.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `phenotype`, `variants`
- **Keywords:** gwas catalog, rsid, snp, trait association, p-value

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `phenotype` | string | yes |  |
| `variants` | array of string | yes |  |
| `association_limit` | integer | no | Default `200` |
| `api_url` | string | no | Default `https://www.ebi.ac.uk/gwas/rest/api` |
| `semantic_context` | object | no | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "Upstream description"

    Compare candidate rsIDs by population-trait GWAS Catalog association evidence.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `domain_id`
- **Keywords:** nutrigenomics, genetic markers, diet, nutrient metabolism, snp

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `domain_id` | string | yes |  |
| `min_evidence_tier` | string | no | One of `established`, `probable`, `emerging` Default `established` |
| `semantic_context` | object | no | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "Upstream description"

    Retrieve curated single-marker records for a declared nutrigenomic domain, filtered by minimum evidence tier.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `rsid`
- **Keywords:** rsid, nutrigenomics, diet, snp, trait

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `rsid` | string | yes |  |

??? note "Upstream description"

    Retrieve any nutrigenomic catalogue records referencing a specific rsID.

Call it through DNAKit once the server is installed:

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

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** gene set, reactome, kegg, msigdb hallmark, pathway genes

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `pathway_id_or_name` | string | no |  |
| `pathway_id` | string | no |  |
| `pathway_name` | string | no |  |
| `source` | string | no | One of `reactome`, `kegg`, `msigdb_hallmark` |
| `species` | string | no |  |
| `reactome_api_base` | string | no | Default `https://reactome.org/ContentService` |
| `kegg_api_base` | string | no | Default `https://rest.kegg.jp` |
| `msigdb_gmt` | string | no | Local filesystem path. |
| `msigdb_gmt_url` | string | no |  |
| `msigdb_version` | string | no |  |
| `limit` | integer | no | Default `500` |
| `semantic_context` | object | no | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "Upstream description"

    Retrieve canonical member genes for a controlled pathway or gene-set source entity.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** pharmacogene, star allele, cyp2d6, hla typing, g6pd

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string | no |  |
| `semantic_context` | object | no | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "Upstream description"

    Return pharmacogene-specific sample evidence requirements for PharmCAT named allele matching, outside calls, CYP2D6 SV/CNV handling, HLA typing, MT-RNR1, and G6PD chrX representation.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** clinpgx, pharmgkb, cpic, dpwg, rsid, drug response

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `drug` | string | no |  |
| `gene` | string | no |  |
| `rsid` | string | no |  |
| `chemical_id` | string | no |  |
| `gene_id` | string | no |  |
| `variant_id` | string | no |  |
| `guideline_source` | string | no | One of `all`, `cpic`, `dpwg`, `pro` Default `all` |
| `include_clinical_annotations` | boolean | no | Default `True` |
| `include_labels` | boolean | no | Default `True` |
| `include_raw_records` | boolean | no | Default `False` |
| `limit` | integer | no | Default `10` |
| `api_url` | string | no | Default `https://api.pharmgkb.org/v1` |
| `semantic_context` | object | no | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "Upstream description"

    Fetch traceable ClinPGx pharmacogenomic guideline, clinical annotation, and FDA label evidence for a selected drug, gene, or rsID; compact normalized records are returned by default.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** pgxdb, rsid, variant marker, drug response, drugbank

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `drug` | string | no |  |
| `atc_code` | string | no |  |
| `drugbank_id` | string | no |  |
| `rsid` | string | no |  |
| `gene` | string | no |  |
| `variant_marker` | string | no |  |
| `include_raw_records` | boolean | no | Default `False` |
| `limit` | integer | no | Default `25` |
| `api_url` | string | no | Default `https://pgx-db.org/rest-api` |
| `semantic_context` | object | no | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "Upstream description"

    Fetch targeted PGxDB pharmacogenomic evidence for a selected drug, ATC code, DrugBank ID, rsID, variant marker, or gene.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** pharmcat, diplotype, star allele, metabolizer phenotype, pgx report

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `output_dir` | string | no | Local filesystem path. |
| `base_filename` | string | no |  |
| `report_json` | string | no | Local filesystem path. |
| `calls_only_tsv` | string | no | Local filesystem path. |
| `match_json` | string | no | Local filesystem path. |
| `phenotype_json` | string | no | Local filesystem path. |
| `missing_pgx_positions_vcf` | string | no | Local filesystem path. |

??? note "Upstream description"

    Import existing PharmCAT artifacts without executing PharmCAT.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** pharmcat, preflight, vcf format, pgx positions, sample columns

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `agi_path` | string | no | Local filesystem path. |

??? note "Upstream description"

    Inspect selected Active Genome Index structure for broad PharmCAT PGx calling without running PharmCAT or writing artifacts.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** outside call, optitype, stellarpgx, cyp2d6, hla, diplotype

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `caller_output_file` | string | no | Local filesystem path. |
| `caller_format` | string | no | One of `auto`, `pharmcat_tsv`, `generic_table`, `optitype`, `stellarpgx_summary` Default `auto` |
| `output_file` | string | no | Local filesystem path. |
| `sample` | string | no |  |
| `max_rows` | integer | no | Default `200` |

??? note "Upstream description"

    Prepare a PharmCAT outside-call TSV from supported specialized caller output such as OptiType HLA calls, StellarPGx CYP2D6 summaries, or generic gene/diplotype tables.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** medication review, drug-gene, rsid, pgx, dosing

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `db` | string | no | Local filesystem path. |
| `shared_db` | string | no | Local filesystem path. |
| `drug` | string | no |  |
| `atc_code` | string | no |  |
| `drugbank_id` | string | no |  |
| `rsid` | string | no |  |
| `gene` | string | no |  |
| `indication` | string | no |  |
| `dose` | string | no |  |
| `current_medications` | string | no |  |
| `allergies_or_contraindications` | string | no |  |
| `known_genotype` | string | no |  |
| `known_diplotype` | string | no |  |
| `known_phenotype` | string | no |  |
| `known_activity_score` | string | no |  |
| `known_pgx_source` | string | no |  |
| `source_sample_pgx_row_id` | string | no |  |
| `semantic_context` | object | no | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |
| `genome_build` | string | no | Default `GRCh38` |
| `include_active_genome_index` | boolean | no |  |
| `include_known_active_genome_indexes` | boolean | no | Default `False` |
| `include_stored_research` | boolean | no | Default `True` |
| `include_record_research_payloads` | boolean | no | Default `False` |
| `limit` | integer | no | Default `10` |
| `clinpgx_api_url` | string | no |  |
| `pgxdb_api_url` | string | no |  |
| `fda_biomarkers_url` | string | no |  |
| `fda_associations_url` | string | no |  |

??? note "Upstream description"

    Review pharmacogenomic evidence for a selected medication, drug-gene, or drug-variant target.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** pharmcat, star allele, diplotype, pgx calling, metabolizer

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `agi_path` | string | no | Local filesystem path. |
| `output_dir` | string | no | Local filesystem path. |
| `base_filename` | string | no |  |
| `mode` | string | no | One of `auto`, `pipeline`, `jar` Default `auto` |
| `pipeline_command` | string | no | Local filesystem path. |
| `pharmcat_jar` | string | no | Local filesystem path. |
| `java_command` | string | no | Local filesystem path. |
| `sample` | string | no |  |
| `sample_file` | string | no | Local filesystem path. |
| `sample_metadata` | string | no | Local filesystem path. |
| `outside_call_file` | string | no | Local filesystem path. |
| `reporter_sources` | string | no |  |
| `research_mode` | string | no |  |
| `max_memory` | string | no |  |
| `max_concurrent_processes` | integer | no |  |
| `include_reporter_json` | boolean | no | Default `True` |
| `include_calls_only_tsv` | boolean | no | Default `True` |
| `probe_version` | boolean | no | Default `True` |
| `dry_run` | boolean | no | Default `False` |
| `timeout_seconds` | integer | no | Default `7200` |

??? note "Upstream description"

    Run a local PharmCAT installation from an approved Active Genome Index.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** outside call, diplotype, activity score, validation, pharmcat

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `outside_call_file` | string | no | Local filesystem path. |
| `max_rows` | integer | no | Default `200` |

??? note "Upstream description"

    Validate PharmCAT outside-call TSV structure and summarize selected diplotype, phenotype, or activity-score evidence.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** none
- **Keywords:** hpo, rare disease, differential, phenotype matching, gene-disease

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `phenotypes` | array of string | no |  |
| `terms` | array of string | no |  |
| `hpo_ids` | array of string | no |  |
| `search_stored_research` | boolean | no | Default `True` |
| `use_hpo_annotations` | boolean | no | Default `True` |
| `source_records` | array of object | no |  |
| `db` | string | no | Local filesystem path. |
| `phenotype_text` | string | no |  |
| `text` | string | no |  |
| `candidate_diseases` | array of string | no |  |
| `diseases` | array of string | no |  |
| `genes` | array of string | no |  |
| `gene_symbols` | array of string | no |  |
| `hpo_disease_file` | string | no | Local filesystem path. |
| `use_primary_gene_disease` | boolean | no | Default `True` |
| `gencc_file` | string | no | Local filesystem path. |
| `limit` | integer | no | Default `25` |
| `semantic_context` | object | no | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "Upstream description"

    Compare supplied or primary gene-derived diseases against phenotype/HPO evidence without selecting the diagnosis.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `genes`
- **Keywords:** drug target, mechanism of action, candidate genes, target evidence, indication

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `source_records` | array of object | no |  |
| `db` | string | no | Local filesystem path. |
| `drug` | string | no |  |
| `drug_class` | string | no |  |
| `indication` | string | no |  |
| `phenotype` | string | no |  |
| `mechanism` | string | no |  |
| `genes` | array of string | yes |  |
| `search_stored_research` | boolean | no | Default `True` |
| `limit` | integer | no | Default `25` |
| `semantic_context` | object | no | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "Upstream description"

    Compare candidate genes using direct drug-target or mechanism evidence only.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `genes`
- **Keywords:** hpo, gene-phenotype, rare disease, candidate genes, annotation

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `phenotypes` | array of string | no |  |
| `terms` | array of string | no |  |
| `hpo_ids` | array of string | no |  |
| `search_stored_research` | boolean | no | Default `True` |
| `use_hpo_annotations` | boolean | no | Default `True` |
| `hpo_gene_file` | string | no | Local filesystem path. |
| `source_records` | array of object | no |  |
| `db` | string | no | Local filesystem path. |
| `phenotype` | string | no |  |
| `phenotype_text` | string | no |  |
| `text` | string | no |  |
| `task_text` | string | no |  |
| `question` | string | no |  |
| `prompt` | string | no |  |
| `condition` | string | no |  |
| `genes` | array of string | yes |  |
| `limit` | integer | no | Default `25` |
| `semantic_context` | object | no | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "Upstream description"

    Compare candidate genes using phenotype, HPO, and curated rare-disease annotation evidence only.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** risk investigation, carrier screening, cancer risk, clinvar, candidate genes

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `question` | string | no |  |
| `investigation_type` | string | no | One of `auto`, `rare_disease`, `cancer_risk`, `carrier_review`, `observed_condition_review` Default `auto` |
| `gene` | string | no |  |
| `genes` | array of string | no |  |
| `condition` | string | no |  |
| `topic` | string | no |  |
| `db` | string | no | Local filesystem path. |
| `matches` | string | no | Local filesystem path. |
| `genome_build` | string | no | Default `GRCh38` |
| `include_active_genome_index` | boolean | no |  |
| `search_stored_research` | boolean | no | Default `True` |
| `limit` | integer | no | Default `25` |
| `semantic_context` | object | no | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "Upstream description"

    Plan phenotype risk and ClinVar review targets from public targets and optionally selected Active Genome Index evidence.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** none
- **Keywords:** open targets, drug target, disease, target genes, clinical stage

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `disease` | string | no |  |
| `disease_id` | string | no |  |
| `genes` | array of string | no |  |
| `mode` | string | no | One of `records`, `gene_membership` Default `records` |
| `minimum_clinical_stage` | string | no | Default `PHASE_2` |
| `opentargets_api_url` | string | no | Default `https://api.platform.opentargets.org/api/v4/graphql` |
| `limit` | integer | no | Default `100` |
| `semantic_context` | object | no | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "Upstream description"

    Retrieve disease-scoped clinical drug-target genes from Open Targets drug candidate records.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** none
- **Keywords:** gencc, gene-disease validity, mendelian, gene symbol, disease

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `genes` | array of string | no |  |
| `gene_symbols` | array of string | no |  |
| `semantic_context` | object | no | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |
| `gencc_file` | string | no | Local filesystem path. |
| `classifications` | array of string | no |  |
| `limit` | integer | no | Default `100` |

??? note "Upstream description"

    Retrieve primary gene-disease associations from GenCC for supplied gene symbols.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `trait`
- **Keywords:** trait, gene association, open targets, genetic evidence, phenotype

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `trait` | string | yes |  |
| `genes` | array of string | no |  |
| `opentargets_api_url` | string | no | Default `https://api.platform.opentargets.org/api/v4/graphql` |
| `limit` | integer | no | Default `25` |
| `semantic_context` | object | no | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |

??? note "Upstream description"

    Retrieve native trait-to-gene records from integrated public sources, optionally filtered to gene symbols.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** polygenic risk score, pgs, weighted score, genotype dosage, personal genome

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `agi_path` | string | no | Local filesystem path. |
| `pgs_id` | string | no |  |
| `score_dir` | string | no | Local filesystem path. |
| `genome_build` | string | no | Default `GRCh38` |
| `skip_ambiguous_palindromic` | boolean | no | Default `True` |
| `score_mean` | number | no |  |
| `score_sd` | number | no |  |

??? note "Upstream description"

    Apply a published polygenic score to an approved Active Genome Index and return raw weighted score plus QC.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** pgs, variant overlap, match rate, polygenic score, coverage

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `agi_path` | string | no | Local filesystem path. |
| `pgs_id` | string | no |  |
| `score_dir` | string | no | Local filesystem path. |
| `genome_build` | string | no | Default `GRCh38` |
| `skip_ambiguous_palindromic` | boolean | no | Default `True` |

??? note "Upstream description"

    Check how many variants from a polygenic score are usable in an approved Active Genome Index.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `pgs_id`
- **Keywords:** pgs catalog, pgs id, scoring file, polygenic score, metadata

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `pgs_id` | string | yes |  |

??? note "Upstream description"

    Fetch detailed public PGS Catalog metadata for one score ID, including scoring-file URLs and source publication context.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** scoring file, variant weights, pgs catalog, import, genome build

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `pgs_id` | string | no |  |
| `genome_build` | string | no | Default `GRCh38` |
| `scoring_file` | string | no | Local filesystem path. |
| `scoring_url` | string | no |  |
| `force` | boolean | no | Default `False` |

??? note "Upstream description"

    Import a PGS Catalog or local scoring file into Genomi's local PRS score cache for a declared genome build.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** pgs catalog, polygenic score, trait, efo, search

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | no |  |
| `trait` | string | no |  |
| `pgs_id` | string | no |  |
| `efo_id` | string | no |  |
| `semantic_context` | object | no | Host-supplied semantic retrieval terms. Genomi reports which terms hit local/public source records and which terms had no hit in the consulted scope. |
| `limit` | integer | no | Default `20` |

??? note "Upstream description"

    Search public PGS Catalog score metadata by trait, score ID, EFO term, or free-text query without reading Active Genome Index.

Call it through DNAKit once the server is installed:

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

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** gencode, encode ccre, genomic region, gene features, grch38

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chrom` | string | no |  |
| `start` | integer | no |  |
| `end` | integer | no |  |
| `region` | string | no |  |
| `assembly` | string | no | One of `GRCh37`, `GRCh38` |
| `gencode_gtf` | string | no | Local filesystem path. |
| `encode_ccre_bed` | string | no | Local filesystem path. |
| `limit` | integer | no | Default `100` |

??? note "Upstream description"

    Retrieve genomic-region feature annotations from supplied or installed GENCODE and ENCODE annotation files.

Call it through DNAKit once the server is installed:

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

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `target_type`
- **Keywords:** evidence packet, variant evidence, gene evidence, chrom pos ref alt, target

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chrom` | string | no |  |
| `pos` | integer | no |  |
| `ref` | string | no |  |
| `alt` | string | no |  |
| `gene` | string | no |  |
| `drug` | string | no |  |
| `condition` | string | no |  |
| `topic` | string | no |  |
| `db` | string | no | Local filesystem path. |
| `target_type` | string | yes | One of `condition`, `drug`, `gene`, `topic`, `variant` |
| `genome_build` | string | no | Default `GRCh38` |
| `source_id` | string | no |  |
| `limit` | integer | no | Default `20` |

??? note "Upstream description"

    Build a target-centric evidence packet after the agent identifies the user's target.

Call it through DNAKit once the server is installed:

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

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** gc content, orf, translation, restriction sites, kozak

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes |  |
| `mode` | string | no | One of `summary`, `translate`, `orfs`, `restriction_sites`, `kozak` Default `summary` |
| `frame` | integer | no | Default `1` |
| `strand` | string | no | One of `forward`, `reverse`, `both` |
| `min_aa` | integer | no | Default `30` |
| `enzymes` | array of string | no |  |
| `motifs` | array of string | no |  |
| `start_pos` | integer | no |  |
| `reference_fasta` | string | no | Local filesystem path. |
| `max_matches` | integer | no | Default `10` |

??? note "Upstream description"

    Run a compact deterministic sequence analysis bundle and point to focused sequence tools when needed.

Call it through DNAKit once the server is installed:

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

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** `forward_primer`
- **Keywords:** pcr primers, melting temperature, amplicon, gc content, primer design

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `forward_primer` | string | yes |  |
| `reverse_primer` | string | no |  |
| `template` | string | no |  |

??? note "Upstream description"

    Check basic primer properties and optional template amplicons.

Call it through DNAKit once the server is installed:

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

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** kozak, start codon, atg, translation initiation

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes |  |
| `start_pos` | integer | no |  |

??? note "Upstream description"

    Classify Kozak sequence context around ATG start codons.

Call it through DNAKit once the server is installed:

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

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** orf, open reading frame, start codon, stop codon, strand

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes |  |
| `min_aa` | integer | no | Default `30` |
| `strand` | string | no | One of `forward`, `reverse`, `both` Default `both` |

??? note "Upstream description"

    Find ATG-to-stop open reading frames in a supplied DNA sequence.

Call it through DNAKit once the server is installed:

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

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** `sequence`
- **Keywords:** restriction enzyme, cut sites, motif, cloning, digest

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes |  |
| `enzymes` | array of string | no |  |
| `motifs` | array of string | no |  |

??? note "Upstream description"

    Find common restriction enzyme or custom motif sites in a supplied DNA sequence.

Call it through DNAKit once the server is installed:

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

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`, `reference_fasta`
- **Keywords:** fasta, sequence matching, reference, alignment, identification

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes |  |
| `reference_fasta` | string | yes | Local filesystem path. |
| `max_matches` | integer | no | Default `10` |

??? note "Upstream description"

    Match a supplied DNA sequence against local FASTA records and return record identifiers plus annotations.

Call it through DNAKit once the server is installed:

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

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** translation, genetic code, reading frame, codon, protein sequence

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes |  |
| `frame` | integer | no | Default `1` |
| `strand` | string | no | One of `forward`, `reverse` Default `forward` |

??? note "Upstream description"

    Translate a DNA sequence in a selected frame and strand using the standard genetic code.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `genes`
- **Keywords:** gene variants, gencode, personal genome, candidate genes, vcf

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `genes` | array of string | yes |  |
| `genome_build` | string | no | One of `GRCh37`, `GRCh38` Default `GRCh38` |
| `agi_id` | string | no |  |
| `per_gene_limit` | integer | no | Default `100` |

??? note "Upstream description"

    Find passing variants in an approved Active Genome Index that overlap managed GENCODE intervals for a bounded candidate-gene set.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `chrom`, `pos`, `ref`, `alt`
- **Keywords:** allele, population frequency, evidence, chrom pos ref alt, genotype

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chrom` | string | yes |  |
| `pos` | integer | yes |  |
| `ref` | string | yes |  |
| `alt` | string | yes |  |
| `db` | string | no | Local filesystem path. |
| `matches` | string | no | Local filesystem path. |
| `genome_build` | string | no | Default `GRCh38` |
| `population_source` | string | no |  |
| `population` | string | no |  |

??? note "Upstream description"

    Gather existing sample, static, population, and reviewed research evidence for one allele.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `gene`
- **Keywords:** gene, clinvar, sample variants, evidence, personal genome

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `db` | string | no | Local filesystem path. |
| `gene` | string | yes |  |
| `matches` | string | no | Local filesystem path. |
| `genome_build` | string | no | Default `GRCh38` |

??? note "Upstream description"

    Gather existing sample, ClinVar, and reviewed research evidence for one gene.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** rsid, variant lookup, hgvs, genotype, clinvar

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `db` | string | no | Local filesystem path. |
| `shared_db` | string | no | Local filesystem path. |
| `chrom` | string | no |  |
| `pos` | integer | no |  |
| `ref` | string | no |  |
| `alt` | string | no |  |
| `query` | string | no |  |
| `rsid` | string | no |  |
| `region` | string | no |  |
| `genome_build` | string | no | Default `GRCh38` |
| `agi_id` | string | no |  |
| `include_active_genome_index` | boolean | no |  |
| `include_known_active_genome_indexes` | boolean | no | Default `False` |
| `include_fail` | boolean | no | Default `False` |
| `limit` | integer | no | Default `20` |

??? note "Upstream description"

    Resolve one variant target and return deterministic public, local sample, and stored evidence facts.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genomi',
    'genomi.invoke',
    {'tool': 'variant.resolve', 'params': {}},
)
```
