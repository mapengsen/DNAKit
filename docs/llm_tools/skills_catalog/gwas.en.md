# GWAS and causal genetics

**26** Skills. A Skill is a set of instructions that guides an Agent through an analysis; the software, reference data and APIs it uses are separate dependencies. See [Skills](../skills.md) for searching, installing and using them; the Skills tab of the web demo shows each Skill in full.

| Skill | Source | Summary |
| --- | --- | --- |
| [`causal-genomics-colocalization-analysis`](#causal-genomics-colocalization-analysis) | bioskills | Test whether two or more traits share a causal variant at a locus using Bayesian colocalization (coloc.abf, coloc.susie, HyPrColoc, moloc, eCAVIAR, SMR/HEIDI,… |
| [`causal-genomics-effector-gene-prioritization`](#causal-genomics-effector-gene-prioritization) | bioskills | Maps GWAS-implicated loci to candidate effector (causal) genes by integrating variant-to-gene (V2G) features via Open Targets L2G (Mountjoy 2021), MAGMA gene-b… |
| [`causal-genomics-fine-mapping`](#causal-genomics-fine-mapping) | bioskills | Resolves GWAS associations to candidate causal variants and credible sets via SuSiE, susie_rss, FINEMAP, CAVIAR, DAP-G, PAINTOR, PolyFun, SuSiEx, MultiSuSiE, a… |
| [`causal-genomics-genetic-correlation`](#causal-genomics-genetic-correlation) | bioskills | Estimates bivariate genetic correlation (rg) between traits from GWAS summary statistics or individual-level genotypes using cross-trait LDSC, HDL, LAVA, rho-H… |
| [`causal-genomics-genomic-sem`](#causal-genomics-genomic-sem) | bioskills | Fits structural equation models to GWAS summary statistics using GenomicSEM (Grotzinger 2019), including common-factor models, confirmatory factor models, ESEM… |
| [`causal-genomics-heritability-partitioning`](#causal-genomics-heritability-partitioning) | bioskills | Estimates SNP heritability and partitions it across functional annotations, cell types, and loci from GWAS summary statistics or individual-level genotypes. Im… |
| [`causal-genomics-mediation-analysis`](#causal-genomics-mediation-analysis) | bioskills | Decompose total effects into direct and indirect paths through mediators using mediation, CMAverse 4-way, HIMA/HIMA2 high-dimensional, BAMA, two-step / MVMR me… |
| [`causal-genomics-mendelian-randomization`](#causal-genomics-mendelian-randomization) | bioskills | Estimate causal effects of an exposure on an outcome from GWAS summary statistics using genetic instruments. Implements IVW (fixed/random), MR-Egger, weighted… |
| [`causal-genomics-pleiotropy-detection`](#causal-genomics-pleiotropy-detection) | bioskills | Detect and adjust for horizontal pleiotropy in two-sample Mendelian randomization by distinguishing uncorrelated (UHP) from correlated (CHP) pleiotropy and cho… |
| [`causal-genomics-proteome-mr-drug-target`](#causal-genomics-proteome-mr-drug-target) | bioskills | Runs cis-pQTL Mendelian randomization for drug-target validation using UKB-PPP (Olink), deCODE (SomaScan), Fenland, INTERVAL, ARIC, and FinnGen-PPP proteomes p… |
| [`causal-genomics-transcriptome-wide-association`](#causal-genomics-transcriptome-wide-association) | bioskills | Performs gene-level association from GWAS summary statistics via genetically predicted tissue expression using FUSION, PrediXcan, S-PrediXcan, S-MultiXcan, UTM… |
| [`data-visualization-manhattan-qq-locuszoom`](#data-visualization-manhattan-qq-locuszoom) | bioskills | Build Manhattan, Miami, QQ, and locuszoom-style regional plots from GWAS, TWAS, PWAS, and QTL summary statistics with correct genomic-inflation diagnostics, mu… |
| [`finemap`](#finemap) | statgen-skills | Fine-map GWAS loci with SuShiE/SuSiE or SuSiE-inf using individual genotypes or summary statistics and matched LD panels. Supports single/multiple ancestry, cr… |
| [`gtex-database`](#gtex-database) | deepmind | Use when you want to retrieve quantitative RNA expression data and variant eQTL information from the GTEx (Genotype-Tissue Expression) Project across 54 non-di… |
| [`gwas-drug-discovery`](#gwas-drug-discovery) | tooluniverse | Transform GWAS signals into drug targets and repurposing opportunities. Connects GWAS-significant loci to causal genes via fine-mapping/eQTL, then to druggable… |
| [`gwas-snp-interpretation`](#gwas-snp-interpretation) | tooluniverse | Interpret a single GWAS SNP across multiple databases — GWAS Catalog hits, LD/haplotype context, eQTL evidence, regulatory annotation, ClinVar pathogenicity, g… |
| [`gwas-trait-to-gene`](#gwas-trait-to-gene) | tooluniverse | Discover causal genes for diseases/traits from GWAS data using Open Targets L2G (locus-to-gene) scoring — integrates eQTL, chromatin interaction, and distance… |
| [`ldsc`](#ldsc) | statgen-skills | Run LD Score Regression (LDSC) for sumstats munging, SNP heritability (h2), genetic correlation (rg), and stratified/cell-type-specific partitioned heritabilit… |
| [`locuscompare-region-render`](#locuscompare-region-render) | clawbio | Render a 4-panel regional LocusCompare diagnostic for one (lead variant, exposure study, outcome study) tuple - overlays GWAS Manhattan, QTL Manhattan, GENCODE… |
| [`locuszoom`](#locuszoom) | statgen-skills | Draw LocusZoom-style regional GWAS association plots with locuszoomr, Ensembl gene tracks and LDlink LD. Select a gene, region or lead SNP and optionally add f… |
| [`manhattan-qq`](#manhattan-qq) | statgen-skills | Produce a Manhattan plot, QQ plot, and a QC summary (genomic inflation λ_GC, count of genome-wide-significant and suggestive hits, list of independent lead SNP… |
| [`pathway-disease-genetics`](#pathway-disease-genetics) | tooluniverse | Connect GWAS variants to biological pathways and druggable targets. Maps GWAS hits to causal genes (via fine-mapping/eQTL), then to pathways (Reactome, KEGG, W… |
| [`phewas`](#phewas) | tooluniverse | Cross-ancestry / cross-biobank phenome-wide association (PheWAS) and replication. Given ONE variant (rsID) or ONE gene, look up every phenotype it associates w… |
| [`prs`](#prs) | statgen-skills | Build a polygenic risk score from GWAS summary statistics using a parameterised PRS method. Currently implements SBayesRC (Zheng et al. 2024, Nat Genet) via th… |
| [`twas`](#twas) | statgen-skills | Run a transcriptome-wide association study (TWAS) on a GWAS sumstats file using TWAS-FUSION (Gusev et al. 2016, Nat Genet) with pre-computed GTEx v8 multi-tiss… |
| [`ukb-ppp-region-fetch`](#ukb-ppp-region-fetch) | user335-clawbio | Fetch a regional slice of plasma pQTL summary statistics from the UK Biobank Pharma Proteomics Project (UKB-PPP; Sun 2023 Nature) for a specific (protein, ance… |

## `causal-genomics-colocalization-analysis` {#causal-genomics-colocalization-analysis}

Test whether two or more traits share a causal variant at a locus using Bayesian colocalization (coloc.abf, coloc.susie, HyPrColoc, moloc, eCAVIAR, SMR/HEIDI, PWCoCo, SharePro). Use when integrating GWAS with eQTL/sQTL/pQTL/mQTL, distinguishing shared causal variants from LD-driven coincidence, handling allelic heterogeneity, choosing between single-causal vs multi-causal methods, picking PP.H4 thresholds, running sensitivity over p12, or harmonising summary statistics for colocalization.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/causal-genomics/colocalization-analysis)
- **License:** `MIT`
- **Also known as:** `bio-causal-genomics-colocalization-analysis`, `pgwas-coloc`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info causal-genomics-colocalization-analysis
dnakit skills install causal-genomics-colocalization-analysis --agent claude-code   # or --agent codex / cursor
```

## `causal-genomics-effector-gene-prioritization` {#causal-genomics-effector-gene-prioritization}

Maps GWAS-implicated loci to candidate effector (causal) genes by integrating variant-to-gene (V2G) features via Open Targets L2G (Mountjoy 2021), MAGMA gene-based association (de Leeuw 2015), FUMA SNP2GENE, cS2G combined SNP-to-gene scores (Gazal 2022), Polygenic Priority Scores (PoPS, Weeks 2023), FLAMES, INQUISIT, DEPICT, and enhancer-gene predictors (ABC, ENCODE-rE2G). Use when narrowing a GWAS lead locus to a candidate causal gene, picking between proximity, eQTL-based, and similarity-based prioritizers, integrating multi-evidence streams (fine-mapping, colocalization, ABC enhancer-gene, distance, chromatin), reconciling discordant L2G vs PoPS calls, prioritizing tissue-specific eQTL evidence, or triangulating across at least three independent lines of evidence for a publication-grade effector-gene nomination.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/causal-genomics/effector-gene-prioritization)
- **License:** `MIT`
- **Also known as:** `bio-causal-genomics-effector-gene-prioritization`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info causal-genomics-effector-gene-prioritization
dnakit skills install causal-genomics-effector-gene-prioritization --agent claude-code   # or --agent codex / cursor
```

## `causal-genomics-fine-mapping` {#causal-genomics-fine-mapping}

Resolves GWAS associations to candidate causal variants and credible sets via SuSiE, susie_rss, FINEMAP, CAVIAR, DAP-G, PAINTOR, PolyFun, SuSiEx, MultiSuSiE, and FOCUS. Use when narrowing a GWAS lead SNP to a 95 percent credible set, choosing between in-sample and reference LD, calibrating non-sparse loci with SuSiE-inf or FINEMAP-inf, integrating functional priors via PolyFun, fine-mapping across ancestries with SuSiEx, diagnosing LD mismatch via estimate_s_rss and kriging_rss, handling HLA or long-range LD, or feeding credible sets into coloc.susie for colocalization.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/causal-genomics/fine-mapping)
- **License:** `MIT`
- **Also known as:** `bio-causal-genomics-fine-mapping`, `fine-mapping`, `tooluniverse-gwas-finemapping`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info causal-genomics-fine-mapping
dnakit skills install causal-genomics-fine-mapping --agent claude-code   # or --agent codex / cursor
```

## `causal-genomics-genetic-correlation` {#causal-genomics-genetic-correlation}

Estimates bivariate genetic correlation (rg) between traits from GWAS summary statistics or individual-level genotypes using cross-trait LDSC, HDL, LAVA, rho-HESS, GREML-bivariate, Popcorn, and HDL-L. Use when quantifying shared genetic architecture between two traits, screening MR validity before causal inference, distinguishing global from locus-level rg, estimating trans-ancestry rg, separating partial from full causation via LCV gcp, or producing a STROBE-MR-compliant cross-trait sensitivity battery. Cross-trait LDSC intercept absorbs sample overlap and is NOT a bias; HDL is biased under sample overlap above ~5%. High rg between exposure and outcome motivates CHP-aware MR sensitivity (CAUSE, LHC-MR).

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/causal-genomics/genetic-correlation)
- **License:** `MIT`
- **Also known as:** `bio-causal-genomics-genetic-correlation`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info causal-genomics-genetic-correlation
dnakit skills install causal-genomics-genetic-correlation --agent claude-code   # or --agent codex / cursor
```

## `causal-genomics-genomic-sem` {#causal-genomics-genomic-sem}

Fits structural equation models to GWAS summary statistics using GenomicSEM (Grotzinger 2019), including common-factor models, confirmatory factor models, ESEM, common-factor GWAS with Q_SNP heterogeneity, multivariate Wald tests, and stratified GenomicSEM partitioned heritability. Reconciles results against MTAG multi-trait analysis. Handles sample overlap via the LDSC sampling-covariance matrix, identifies and resolves Heywood cases, and verifies model fit with CFI / RMSEA. Use when modeling latent genetic architecture across correlated traits, running multivariate GWAS on a shared factor, distinguishing factor-mediated from trait-specific SNP effects, or comparing GenomicSEM common-factor results against MTAG when both depend on accurate sampling covariance.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/causal-genomics/genomic-sem)
- **License:** `MIT`
- **Also known as:** `bio-causal-genomics-genomic-sem`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info causal-genomics-genomic-sem
dnakit skills install causal-genomics-genomic-sem --agent claude-code   # or --agent codex / cursor
```

## `causal-genomics-heritability-partitioning` {#causal-genomics-heritability-partitioning}

Estimates SNP heritability and partitions it across functional annotations, cell types, and loci from GWAS summary statistics or individual-level genotypes. Implements LDSC, stratified LDSC with the baseline-LD model, Finucane 2018 cell-type prioritization, LDAK SumHer, HDL, HESS local heritability, BOLT-REML, GCTA-GREML, graphREML, and Popcorn cross-population genetic correlation. Use when computing total h2_SNP from summary stats, partitioning heritability across functional categories, prioritizing trait-relevant tissues or cell types from ENCODE/Roadmap chromatin marks, reconciling LDSC vs LDAK enrichment estimates, computing local heritability with HESS, estimating genetic correlation between traits, or producing publication-grade enrichment with calibrated sensitivity to model assumptions.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/causal-genomics/heritability-partitioning)
- **License:** `MIT`
- **Also known as:** `bio-causal-genomics-heritability-partitioning`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info causal-genomics-heritability-partitioning
dnakit skills install causal-genomics-heritability-partitioning --agent claude-code   # or --agent codex / cursor
```

## `causal-genomics-mediation-analysis` {#causal-genomics-mediation-analysis}

Decompose total effects into direct and indirect paths through mediators using mediation, CMAverse 4-way, HIMA/HIMA2 high-dimensional, BAMA, two-step / MVMR mediation, or double-ML medDML. Use when testing whether a molecular phenotype (expression, methylation, protein) mediates a treatment-outcome relationship, decomposing exposure-mediator interaction via VanderWeele 4-way, screening high-dimensional EWAS mediators, or running MR-based mediation when sequential ignorability is implausible.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/causal-genomics/mediation-analysis)
- **License:** `MIT`
- **Also known as:** `bio-causal-genomics-mediation-analysis`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info causal-genomics-mediation-analysis
dnakit skills install causal-genomics-mediation-analysis --agent claude-code   # or --agent codex / cursor
```

## `causal-genomics-mendelian-randomization` {#causal-genomics-mendelian-randomization}

Estimate causal effects of an exposure on an outcome from GWAS summary statistics using genetic instruments. Implements IVW (fixed/random), MR-Egger, weighted median/mode, MR-RAPS, CAUSE, GSMR-HEIDI, MR-PRESSO, MVMR, MR-Clust, LCV, and LHC-MR via TwoSampleMR, MendelianRandomization, MR-PRESSO, cause, and lhcMR. Use when testing causal direction between traits, evaluating drug-target effects via cis-pQTL/cis-eQTL, performing multivariable mediation MR, distinguishing causation from correlated horizontal pleiotropy, or producing STROBE-MR-compliant sensitivity batteries.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/causal-genomics/mendelian-randomization)
- **License:** `MIT`
- **Also known as:** `bio-causal-genomics-mendelian-randomization`, `mendelian-randomisation`, `mendelian-randomization-twosamplemr`, `pgwas-mr`, `tooluniverse-mendelian-randomization`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info causal-genomics-mendelian-randomization
dnakit skills install causal-genomics-mendelian-randomization --agent claude-code   # or --agent codex / cursor
```

## `causal-genomics-pleiotropy-detection` {#causal-genomics-pleiotropy-detection}

Detect and adjust for horizontal pleiotropy in two-sample Mendelian randomization by distinguishing uncorrelated (UHP) from correlated (CHP) pleiotropy and choosing among Egger, MR-PRESSO, MR-RAPS, CAUSE, LHC-MR, LCV, MR-Clust, MR-Mix, and contamination-mixture methods. Use when validating an MR causal claim, running the STROBE-MR sensitivity battery, suspecting a shared heritable confounder, working under weak-instrument or polygenic-exposure regimes, or reconciling discordant estimates across robust methods.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/causal-genomics/pleiotropy-detection)
- **License:** `MIT`
- **Also known as:** `bio-causal-genomics-pleiotropy-detection`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info causal-genomics-pleiotropy-detection
dnakit skills install causal-genomics-pleiotropy-detection --agent claude-code   # or --agent codex / cursor
```

## `causal-genomics-proteome-mr-drug-target` {#causal-genomics-proteome-mr-drug-target}

Runs cis-pQTL Mendelian randomization for drug-target validation using UKB-PPP (Olink), deCODE (SomaScan), Fenland, INTERVAL, ARIC, and FinnGen-PPP proteomes plus colocalization triangulation, phenome-wide on-target adverse-effect scans, cross-platform Olink/SomaScan replication, and PAV (protein-altering variant) sensitivity. Use when nominating or de-risking a drug target from plasma-proteome GWAS, mimicking pharmacological inhibition via cis-pQTL instruments, separating shared-causal from LD-confounded signal under the Schmidt 2020 cis-MR framework, screening on-target adverse phenotypes pheWAS-style, or producing publication-grade STROBE-MR plus PP.H4 evidence for a target gene.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/causal-genomics/proteome-mr-drug-target)
- **License:** `MIT`
- **Also known as:** `bio-causal-genomics-proteome-mr-drug-target`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info causal-genomics-proteome-mr-drug-target
dnakit skills install causal-genomics-proteome-mr-drug-target --agent claude-code   # or --agent codex / cursor
```

## `causal-genomics-transcriptome-wide-association` {#causal-genomics-transcriptome-wide-association}

Performs gene-level association from GWAS summary statistics via genetically predicted tissue expression using FUSION, PrediXcan, S-PrediXcan, S-MultiXcan, UTMOST, MOSTWAS, kTWAS, EpiXcan, TIGAR-V2, and probabilistic fine-mapping with FOCUS and MA-FOCUS. Use when running TWAS from GWAS sumstats, prioritising candidate causal genes from a GWAS lead locus, picking single-tissue vs cross-tissue models, identifying LD-induced TWAS false positives, choosing ancestry-matched prediction weights, fine-mapping co-regulated TWAS hits, or triangulating TWAS with cis-eQTL Mendelian randomization and colocalization to nominate a causal gene.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/causal-genomics/transcriptome-wide-association)
- **License:** `MIT`
- **Also known as:** `bio-causal-genomics-transcriptome-wide-association`, `gwas-to-function-twas`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info causal-genomics-transcriptome-wide-association
dnakit skills install causal-genomics-transcriptome-wide-association --agent claude-code   # or --agent codex / cursor
```

## `data-visualization-manhattan-qq-locuszoom` {#data-visualization-manhattan-qq-locuszoom}

Build Manhattan, Miami, QQ, and locuszoom-style regional plots from GWAS, TWAS, PWAS, and QTL summary statistics with correct genomic-inflation diagnostics, multi-trait overlays, lead-SNP labeling, and LD-aware regional rendering. Use when visualizing association results across the genome, comparing two traits, computing genomic inflation lambda, or zooming into a locus with LD coloring.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/data-visualization/manhattan-qq-locuszoom)
- **License:** `MIT`
- **Also known as:** `bio-data-visualization-manhattan-qq-locuszoom`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info data-visualization-manhattan-qq-locuszoom
dnakit skills install data-visualization-manhattan-qq-locuszoom --agent claude-code   # or --agent codex / cursor
```

## `finemap` {#finemap}

Fine-map GWAS loci with SuShiE/SuSiE or SuSiE-inf using individual genotypes or summary statistics and matched LD panels. Supports single/multiple ancestry, credible sets, optional gnomAD/GTEx annotation and provenance manifests. Use for posterior inclusion probabilities, credible sets or fine-mapping method selection; solver/reference installation is separate.

- **Source:** [Copyright 2026 Stephen Dorn; statgen-skills contributors @ `dc5f0c6`](https://github.com/svdorn/statgen-skills/tree/dc5f0c60d2995aebfe58b3f8a07b46f27562f4e0/finemap)
- **License:** `MIT`
- **Runs with:** Python modules `numpy`, `pandas`; commands `python3`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info finemap
dnakit skills install finemap --agent claude-code   # or --agent codex / cursor
```

## `gtex-database` {#gtex-database}

Use when you want to retrieve quantitative RNA expression data and variant eQTL information from the GTEx (Genotype-Tissue Expression) Project across 54 non-diseased tissue sites.

- **Source:** [Copyright 2026 Google LLC @ `28b8482`](https://github.com/google-deepmind/science-skills/tree/28b8482603a420708c8896f6fe5e06c276d9933d/skills/gtex_database)
- **License:** `Apache-2.0 AND CC-BY-4.0`
- **Also known as:** `alterlab-gtex`, `gtex-expression`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info gtex-database
dnakit skills install gtex-database --agent claude-code   # or --agent codex / cursor
```

## `gwas-drug-discovery` {#gwas-drug-discovery}

Transform GWAS signals into drug targets and repurposing opportunities. Connects GWAS-significant loci to causal genes via fine-mapping/eQTL, then to druggable proteins via DGIdb/OpenTargets, then to existing drugs via ChEMBL. Use for GWAS-to-target hypothesis generation, druggable-fraction analysis of disease loci, and human-genetics-validated drug-repurposing prioritization.

- **Source:** [ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-gwas-drug-discovery)
- **License:** `Apache-2.0`
- **Runs with:** Python modules `dnakit`, `indigo`
- **Also known as:** `tooluniverse-gwas-drug-discovery`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info gwas-drug-discovery
dnakit skills install gwas-drug-discovery --agent claude-code   # or --agent codex / cursor
```

## `gwas-snp-interpretation` {#gwas-snp-interpretation}

Interpret a single GWAS SNP across multiple databases — GWAS Catalog hits, LD/haplotype context, eQTL evidence, regulatory annotation, ClinVar pathogenicity, gnomAD frequency. Use for 'what does this SNP do', SNP-to-mechanism tracing, and resolving lead-SNP-vs-causal-variant ambiguity. Always considers LD structure before claiming a SNP is mechanistically responsible.

- **Source:** [ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-gwas-snp-interpretation)
- **License:** `Apache-2.0`
- **Also known as:** `gwas-lookup`, `tooluniverse-gwas-snp-interpretation`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info gwas-snp-interpretation
dnakit skills install gwas-snp-interpretation --agent claude-code   # or --agent codex / cursor
```

## `gwas-trait-to-gene` {#gwas-trait-to-gene}

Discover causal genes for diseases/traits from GWAS data using Open Targets L2G (locus-to-gene) scoring — integrates eQTL, chromatin interaction, and distance evidence. Use for trait-to-gene mapping, drug-target hypothesis generation from GWAS, and replacing the 'nearest gene' heuristic with multi-evidence L2G scores.

- **Source:** [ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-gwas-trait-to-gene)
- **License:** `Apache-2.0`
- **Also known as:** `tooluniverse-gwas-trait-to-gene`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info gwas-trait-to-gene
dnakit skills install gwas-trait-to-gene --agent claude-code   # or --agent codex / cursor
```

## `ldsc` {#ldsc}

Run LD Score Regression (LDSC) for sumstats munging, SNP heritability (h2), genetic correlation (rg), and stratified/cell-type-specific partitioned heritability (partitioned, Finucane 2015/2018 S-LDSC via `--h2-cts`). Auto-installs LDSC from the CBIIT Python 3 / Mac-compatible fork and caches the canonical EUR LD-score reference files on first use. Resolves method/software/paper provenance from the OKG when $OKG_REPO is set; otherwise runs without OKG provenance. Use when a user asks to munge GWAS sumstats, estimate heritability from sumstats, compute genetic correlation between two GWAS traits, or rank tissues / cell types by heritability enrichment (S-LDSC / Finucane 2018 multi-tissue).

- **Source:** [Copyright 2026 Stephen Dorn; statgen-skills contributors @ `dc5f0c6`](https://github.com/svdorn/statgen-skills/tree/dc5f0c60d2995aebfe58b3f8a07b46f27562f4e0/ldsc)
- **License:** `MIT`
- **Runs with:** Python modules `numpy`, `pandas`, `scipy`; commands `python3`, `git`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info ldsc
dnakit skills install ldsc --agent claude-code   # or --agent codex / cursor
```

## `locuscompare-region-render` {#locuscompare-region-render}

Render a 4-panel regional LocusCompare diagnostic for one (lead variant,
exposure study, outcome study) tuple - overlays GWAS Manhattan, QTL Manhattan,
GENCODE gene track, and cross-trait scatter colored by LD r². Use when an
agent needs visual confirmation that two GWAS / QTL signals share the same
causal variant (the Liu 2019 LocusCompare convention). Inputs: lead variant +
two pre-fetched harmonised sumstats slices (or eQTL Catalogue / GWAS Catalog
identifiers for bundled fetch). Output: PNG + JSON manifest.

- **Source:** [Aviv Madar; ClawBio contributors @ `1fcecb7`](https://github.com/ClawBio/ClawBio/tree/1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d/skills/locuscompare-region-render)
- **License:** `MIT`
- **Runs with:** commands `plink`, `python3`, `tabix`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info locuscompare-region-render
dnakit skills install locuscompare-region-render --agent claude-code   # or --agent codex / cursor
```

## `locuszoom` {#locuszoom}

Draw LocusZoom-style regional GWAS association plots with locuszoomr, Ensembl gene tracks and LDlink LD. Select a gene, region or lead SNP and optionally add fine-mapping posterior probabilities and credible sets. Use for regional association/credible-set figures; requires R and the user's LDlink API token.

- **Source:** [Copyright 2026 Stephen Dorn; statgen-skills contributors @ `dc5f0c6`](https://github.com/svdorn/statgen-skills/tree/dc5f0c60d2995aebfe58b3f8a07b46f27562f4e0/locuszoom)
- **License:** `MIT`
- **Runs with:** commands `python3`, `Rscript`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info locuszoom
dnakit skills install locuszoom --agent claude-code   # or --agent codex / cursor
```

## `manhattan-qq` {#manhattan-qq}

Produce a Manhattan plot, QQ plot, and a QC summary (genomic inflation λ_GC, count of genome-wide-significant and suggestive hits, list of independent lead SNPs) from a GWAS summary statistics file. Auto-detects common column-name conventions (`SNP/rsid/variant_id`, `CHR/chromosome/hm_chrom`, `BP/position/base_pair_location`, `P/p_value`). Subsamples low-significance variants to keep plots fast (~seconds) on full HM3-scale or 1KG-scale sumstats. Writes a sidecar `.manhattan_qq.json` manifest citing the input dataset's OKG node when `$OKG_REPO` is set. Use when a user asks for a Manhattan plot, QQ plot, genomic-inflation check, or a quick QC pass on a GWAS sumstats file.

- **Source:** [Copyright 2026 Stephen Dorn; statgen-skills contributors @ `dc5f0c6`](https://github.com/svdorn/statgen-skills/tree/dc5f0c60d2995aebfe58b3f8a07b46f27562f4e0/manhattan-qq)
- **License:** `MIT`
- **Runs with:** Python modules `numpy`, `pandas`, `scipy`, `matplotlib`; commands `python3`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info manhattan-qq
dnakit skills install manhattan-qq --agent claude-code   # or --agent codex / cursor
```

## `pathway-disease-genetics` {#pathway-disease-genetics}

Connect GWAS variants to biological pathways and druggable targets. Maps GWAS hits to causal genes (via fine-mapping/eQTL), then to pathways (Reactome, KEGG, WikiPathways), then to existing drugs hitting those pathways. Use for pathway-level disease mechanisms, druggable-pathway prioritization from GWAS, SNP-to-pathway-to-target tracing, and tissue-specific eQTL evidence for drug target hypotheses.

- **Source:** [ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-pathway-disease-genetics)
- **License:** `Apache-2.0`
- **Runs with:** Python modules `dnakit`, `networkx`
- **Also known as:** `tooluniverse-pathway-disease-genetics`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info pathway-disease-genetics
dnakit skills install pathway-disease-genetics --agent claude-code   # or --agent codex / cursor
```

## `phewas` {#phewas}

Cross-ancestry / cross-biobank phenome-wide association (PheWAS) and replication. Given ONE variant (rsID) or ONE gene, look up every phenotype it associates with across European/UK (UKB-TOPMed), Finnish (FinnGen), Japanese (BioBank Japan), and Taiwanese (TPMI) biobanks, plus exome-wide gene-burden PheWAS (Genebass), then judge whether an association replicates across ancestries or is population-specific. Use whenever the user asks "what else is this variant/gene associated with", "does this association replicate in other ancestries / biobanks", "is this effect East-Asian-specific", "pleiotropy of rsXXX", "phenome scan", or wants to compare effect sizes/allele frequencies of a variant across populations. NOT for the forward direction (trait → which SNPs: use the gwas-* skills), NOT for fine-mapping a locus (use tooluniverse-gwas-finemapping), and NOT for single-SNP mechanism tracing in one population (use gwas-snp-interpretation).

- **Source:** [ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-phewas)
- **License:** `Apache-2.0`
- **Also known as:** `tooluniverse-phewas`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info phewas
dnakit skills install phewas --agent claude-code   # or --agent codex / cursor
```

## `prs` {#prs}

Build a polygenic risk score from GWAS summary statistics using a parameterised PRS method. Currently implements SBayesRC (Zheng et al. 2024, Nat Genet) via the zhilizheng/SBayesRC R package; LDpred2, PRS-CS, and PRS-CSx can be added behind the same `--method` flag. OKG-aware - resolves the ancestry-matched LD eigendecomposition reference (HapMap3 EUR/EAS/AFR) via `ld_panel:sbayesrc_hm3_(anc)` and pulls sample size from the input dataset_metadata node. Writes per-SNP PRS weights plus a sidecar manifest citing the method/software/paper/LD-panel node IDs. Use when a user asks to compute PRS weights, run SBayesRC, or build a polygenic score from sumstats.

- **Source:** [Copyright 2026 Stephen Dorn; statgen-skills contributors @ `dc5f0c6`](https://github.com/svdorn/statgen-skills/tree/dc5f0c60d2995aebfe58b3f8a07b46f27562f4e0/prs)
- **License:** `MIT`
- **Runs with:** commands `python3`, `Rscript`, `uv`; environment `OKG_REPO`, `OKG_DSN`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info prs
dnakit skills install prs --agent claude-code   # or --agent codex / cursor
```

## `twas` {#twas}

Run a transcriptome-wide association study (TWAS) on a GWAS sumstats file using TWAS-FUSION (Gusev et al. 2016, Nat Genet) with pre-computed GTEx v8 multi-tissue expression weights. Resolves the weight panel via the OKG (`--tissue Brain_Cortex --ancestry EUR` → `dataset:fusion_gtex_v8_eur:brain_cortex`), downloads + caches the panel + the FUSION 1000G LDREF on first use, runs `FUSION.assoc_test.R` per chromosome, concatenates per-gene results, and writes a `.twas.json` sidecar citing the method/software/paper/ld_panel/dataset/tissue/cohort OKG nodes. OKG-aware - if `$OKG_REPO` is set, can also auto-suggest tissues for a given trait via the `tissue → trait relevant_to` mapping (e.g., AD → Brain_Cortex / Brain_Hippocampus / Brain_Frontal_Cortex_BA9 / Whole_Blood). Use when a user asks to run TWAS, compute gene-level association from GWAS sumstats, integrate expression weights with a GWAS, or fine-map a GWAS locus to candidate causal genes.

- **Source:** [Copyright 2026 Stephen Dorn; statgen-skills contributors @ `dc5f0c6`](https://github.com/svdorn/statgen-skills/tree/dc5f0c60d2995aebfe58b3f8a07b46f27562f4e0/twas)
- **License:** `MIT`
- **Runs with:** commands `python3`, `Rscript`, `git`; environment `OKG_REPO`, `OKG_DSN`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info twas
dnakit skills install twas --agent claude-code   # or --agent codex / cursor
```

## `ukb-ppp-region-fetch` {#ukb-ppp-region-fetch}

Fetch a regional slice of plasma pQTL summary statistics from the UK
Biobank Pharma Proteomics Project (UKB-PPP; Sun 2023 Nature) for a
specific (protein, ancestry) measurement. Use when an agent needs
per-variant beta / SE / p-value around a coloc-lead variant for
downstream colocalisation, Mendelian randomisation, or regional
plotting against a pQTL exposure. The canonical use case is the
cis-window around the protein's coding gene TSS, but UKB-PPP releases
full-genome summary stats per protein so any GRCh38 window (including
trans loci) is supported when the user supplies an explicit
(chromosome, start_bp, end_bp). Input: protein_label (HGNC or
UniProt), ancestry, chromosome, start_bp, end_bp. Output: harmonised
TSV slice + manifest + human-readable report.

- **Source:** [Copyright 2026 Aviv Madar; UK Biobank Pharma Proteomics Project, Sun et al. (2023) @ `4848cc2`](https://github.com/ClawBio/ClawBio/tree/4848cc280c446093b1364a2445db3137d455cc17/skills/ukb-ppp-region-fetch)
- **License:** `MIT AND CC-BY-4.0`
- **Runs with:** Python modules `requests`, `synapseclient`; commands `python3`; environment `SYNAPSE_AUTH_TOKEN`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info ukb-ppp-region-fetch
dnakit skills install ukb-ppp-region-fetch --agent claude-code   # or --agent codex / cursor
```
