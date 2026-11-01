# Variant interpretation and disease research

**53** Skills. A Skill is a set of instructions that guides an Agent through an analysis; the software, reference data and APIs it uses are separate dependencies. See [Skills](../skills.md) for searching, installing and using them; the Skills tab of the web demo shows each Skill in full.

| Skill | Source | Summary |
| --- | --- | --- |
| [`analytical-validation`](#analytical-validation) | bioskills | Treats a ctDNA assay as a molecule-counting experiment at the Poisson edge and builds its analytical-validation case the measurement-science way. Covers the ge… |
| [`biomcp`](#biomcp) | biomcp | Search and retrieve biomedical data - genes, variants, clinical trials, diagnostic tests, articles, drugs, diseases, pathways, proteins, adverse events, pharma… |
| [`cancer-genomics-tcga`](#cancer-genomics-tcga) | tooluniverse | TCGA/GDC cancer genomics analysis — cohort construction, clinical metadata retrieval, somatic mutation frequencies, survival analysis, and multi-omics integrat… |
| [`cancer-variant-interpretation`](#cancer-variant-interpretation) | tooluniverse | Clinical interpretation of somatic cancer mutations for precision oncology. Transforms a gene + variant + cancer-type input into an actionable report: clinical… |
| [`cfdna-preprocessing`](#cfdna-preprocessing) | bioskills | Decides how to preprocess plasma cfDNA sequencing data so the recoverable signal survives - library-prep-aware fragment expectations (dsDNA vs ssDNA/adaptase p… |
| [`clinical-databases-acmg-classification`](#clinical-databases-acmg-classification) | bioskills | Applies ACMG/AMP 2015 framework with ClinGen SVI specifications, Tavtigian 2018/2020 Bayesian point system, Abou Tayoun 2018 PVS1 decision tree, Pejaver 2022 a… |
| [`clinical-databases-clinvar-lookup`](#clinical-databases-clinvar-lookup) | bioskills | Queries ClinVar for variant pathogenicity classifications, ClinGen VCEP curations, and somatic-vs-germline interpretations via REST API, weekly VCF, or bulk XM… |
| [`clinical-databases-dbsnp-queries`](#clinical-databases-dbsnp-queries) | bioskills | Resolves rsIDs, navigates RsMergeArch/SNPHistory merge chains, and converts between rsID, SPDI, HGVS, and VCF representations using the dbSNP Build 156 JSON ar… |
| [`clinical-databases-gnomad-frequencies`](#clinical-databases-gnomad-frequencies) | bioskills | Queries gnomAD v4 (807k samples), v3, v2.1.1, and constraint metrics with grpmax FAF95, bottleneck-group exclusion, LOEUF interpretation, SV/CNV/mtDNA catalogs… |
| [`clinical-databases-hla-typing`](#clinical-databases-hla-typing) | bioskills | Calls HLA class I and class II alleles at 2/4/6/8-field resolution from WGS/WES/RNA-seq/long-read data using OptiType, HLA-LA, T1K, Polysolver, HLA-HD, arcasHL… |
| [`clinical-databases-msi-detection`](#clinical-databases-msi-detection) | bioskills | Calls microsatellite instability from WES/WGS/targeted-panel with MSIsensor, MSIsensor-pro, MSIsensor-ct (panel-aware), mSINGS, and MANTIS for FDA pembrolizuma… |
| [`clinical-databases-myvariant-queries`](#clinical-databases-myvariant-queries) | bioskills | Queries myvariant.info BioThings aggregator for ClinVar, gnomAD, dbSNP, dbNSFP, COSMIC, CADD, and CIViC annotations in batched, version-tracked requests. Use w… |
| [`clinical-databases-pharmacogenomics`](#clinical-databases-pharmacogenomics) | bioskills | Queries PharmGKB / CPIC / DPWG for drug-gene interactions; calls CYP2D6/CYP2C9/CYP2C19/DPYD/TPMT/NUDT15/UGT1A1/SLCO1B1 star alleles and phenotype with PharmCAT… |
| [`clinical-databases-polygenic-risk`](#clinical-databases-polygenic-risk) | bioskills | Constructs and validates polygenic risk scores using LDpred2-auto, SBayesRC, MegaPRS, PRS-CS, PROSPER, MUSSEL, BridgePRS, JointPRS, PRSmix, or PGS Catalog Calc… |
| [`clinical-databases-somatic-signatures`](#clinical-databases-somatic-signatures) | bioskills | Extracts and assigns COSMIC v3.4 mutational signatures (86 SBS / 11 DBS / 18 ID / 21 CN / 16 SV) from somatic VCFs using SigProfilerSuite, MutationalPatterns,… |
| [`clinical-databases-tumor-mutational-burden`](#clinical-databases-tumor-mutational-burden) | bioskills | Calculates tumor mutational burden from WES/WGS/panel data with Friends of Cancer Research harmonization equations, per-assay calibration (FDA 10/Mb = 7.8 TSO5… |
| [`clinical-databases-variant-prioritization`](#clinical-databases-variant-prioritization) | bioskills | Prioritizes rare-disease variants from trio/quad WES/WGS with de novo (DeNovoGear, Triodenovo), compound-heterozygous phasing (WhatsHap), mosaic VAF tiering, p… |
| [`clinpgx`](#clinpgx) | clawbio | Query the ClinPGx API for pharmacogenomic gene-drug data, clinical annotations, CPIC guidelines, and FDA drug labels |
| [`ctdna-mutation-detection`](#ctdna-mutation-detection) | bioskills | Detects somatic mutations in circulating tumor DNA, treating low-VAF detection as a signal-versus-noise problem set by error suppression and molecules sampled,… |
| [`data-visualization-oncoprint-mutation-matrices`](#data-visualization-oncoprint-mutation-matrices) | bioskills | Build OncoPrint and co-mutation matrix plots from somatic-variant cohorts using ComplexHeatmap, maftools, and comut.py with alteration-type stacking, sample or… |
| [`draft-counseling-summary`](#draft-counseling-summary) | genetic-counseling | Draft clinician-facing and optional patient-facing rare-disease germline genetic counseling summaries from reviewed case facts and interpretation proposals. Us… |
| [`draw-pedigree`](#draw-pedigree) | genetic-counseling | Create a reviewable clinical pedigree figure or validated PED/JSON file from user-supplied family-history facts. Use when asked to draw, chart, or diagram a pe… |
| [`find-variant-literature`](#find-variant-literature) | genetic-counseling | Discover, retrieve, and extract claim-level literature evidence for a rare-disease germline variant, gene-disease pair, ACMG criterion, or CNV region. Use when… |
| [`folklore-variant-evidence`](#folklore-variant-evidence) | kdense | Retrieve ClinGen gene-disease validity assertions for a public gene or disease, and review source-linked public evidence and literature for one supported GRCh3… |
| [`fragment-analysis`](#fragment-analysis) | bioskills | Extracts cfDNA fragmentomics features (DELFI genome-wide short/long ratios, WPS nucleosome positioning, Griffin GC-corrected accessibility profiles, end-motifs… |
| [`gene-disease-association`](#gene-disease-association) | tooluniverse | Gene-disease association analysis across DisGeNET, OpenTargets, Monarch, OMIM, GenCC, Orphanet. Cross-references multiple sources for evidence-graded associati… |
| [`hla-immunogenomics`](#hla-immunogenomics) | tooluniverse | HLA gene-family analysis and MHC-peptide binding for transplant compatibility, vaccine epitope coverage, and cancer immunotherapy. Uses IMGT (HLA polymorphism)… |
| [`immune-repertoire-analysis`](#immune-repertoire-analysis) | tooluniverse | TCR/BCR repertoire analysis — V(D)J segment usage, CDR3 sequence diversity, clonality scoring, antigen specificity matching to IEDB, public-clone identificatio… |
| [`immunoinformatics-epitope-prediction`](#immunoinformatics-epitope-prediction) | bioskills | Predict B-cell and T-cell epitopes for vaccine antigen design and epitope mapping with BepiPred-3.0, DiscoTope-3.0, the IEDB tools, and EL-mode MHC presentatio… |
| [`immunoinformatics-immunogenicity-scoring`](#immunoinformatics-immunogenicity-scoring) | bioskills | Rank and prioritize neoantigen/epitope candidates by likely T-cell response using NeoFox feature annotation, PRIME2.0, BigMHC-IM, the Łuksza/Balachandran fitne… |
| [`immunoinformatics-mhc-binding-prediction`](#immunoinformatics-mhc-binding-prediction) | bioskills | Predict peptide-MHC class I binding and natural presentation with MHCflurry, NetMHCpan-4.1, and MixMHCpred to nominate candidate CD8 T-cell epitopes. Covers th… |
| [`immunoinformatics-mhc-class-ii-prediction`](#immunoinformatics-mhc-class-ii-prediction) | bioskills | Predict peptide-MHC class II (HLA-DR/DQ/DP) binding and presentation for CD4 T-cell epitopes with NetMHCIIpan-4.3 and MixMHC2pred-2.0. Covers why class II is f… |
| [`immunoinformatics-neoantigen-prediction`](#immunoinformatics-neoantigen-prediction) | bioskills | Identify tumor neoantigens from somatic variants with pVACtools (pVACseq/pVACfuse/pVACbind/pVACvector/pVACview) for personalized cancer vaccines and checkpoint… |
| [`interpret-negative-or-vus`](#interpret-negative-or-vus) | genetic-counseling | Interpret a completed negative or VUS rare-disease germline analysis by checking QC and assay blind spots, generating restrained phenotype-specific hypotheses,… |
| [`longitudinal-monitoring`](#longitudinal-monitoring) | bioskills | Tracks ctDNA across serial liquid-biopsy timepoints for molecular residual disease (MRD) and treatment-response monitoring, treating MRD as a binary integrated… |
| [`methylation-based-detection`](#methylation-based-detection) | bioskills | Detects cancer and infers tissue-of-origin from cfDNA methylation by choosing conversion chemistry (bisulfite vs EM-seq vs TAPS vs cfMeDIP), calling read-level… |
| [`multiomic-disease-characterization`](#multiomic-disease-characterization) | tooluniverse | Comprehensive disease characterization across genomics, transcriptomics, proteomics, and pathways for systems-level understanding. Identifies therapeutic oppor… |
| [`ngs-tertiary-analysis`](#ngs-tertiary-analysis) | ngs-tertiary | Run NGS tertiary analysis pipeline from BAM/VCF to ESMO-compliant clinical HTML report. Use when user has a cancer panel BAM/VCF file (TSO500) and needs somati… |
| [`prepare-rare-disease-case`](#prepare-rare-disease-case) | genetic-counseling | Structure English or Chinese rare-disease germline clinical narratives into reviewable findings, explicit negatives, suspected diagnoses, treatments, confounde… |
| [`protein-lof-mechanism`](#protein-lof-mechanism) | tooluniverse | Propose the mechanism by which a missense variant causes loss-of-function (LoF), synthesizing evidence from 5 independent layers: AlphaMissense pathogenicity,… |
| [`protein-sae-variant-interpretation`](#protein-sae-variant-interpretation) | tooluniverse | Interpret a missense variant via ESMC-6B Sparse Autoencoder (SAE) feature activations. For a given protein + variant, computes which interpretable SAE features… |
| [`rare-disease-diagnosis`](#rare-disease-diagnosis) | tooluniverse | Rare disease differential diagnosis from patient phenotype — HPO term matching to candidate diseases (Orphanet, OMIM), gene panel prioritization, ACMG variant… |
| [`rare-disease-genomics`](#rare-disease-genomics) | tooluniverse | Rare disease genomics — disease identification (Orphanet), causative gene discovery, gene-disease validity (GenCC), variant interpretation (ClinVar), and trans… |
| [`rare-high-impact-variants`](#rare-high-impact-variants) | clawbio | Count rare, high-impact loss-of-function variants carried in a VCF, annotated with molecular consequence and population allele frequency |
| [`review-variant-table`](#review-variant-table) | genetic-counseling | Inspect and normalize testing-company Excel, CSV, or TSV germline variant exports into a reviewable variant inventory. Use when a user uploads a laboratory var… |
| [`tcr-bcr-analysis-immcantation-analysis`](#tcr-bcr-analysis-immcantation-analysis) | bioskills | Reconstructs B-cell clonal families, quantifies somatic hypermutation and selection, and builds antibody lineage trees with the Immcantation R suite (alakazam,… |
| [`tcr-bcr-analysis-mixcr-analysis`](#tcr-bcr-analysis-mixcr-analysis) | bioskills | Align V(D)J reads and assemble TCR/BCR clonotypes with MiXCR, driven by a chemistry-matched preset. Use when choosing/auditing the preset for a library (5'RACE… |
| [`tcr-bcr-analysis-repertoire-visualization`](#tcr-bcr-analysis-repertoire-visualization) | bioskills | Draws TCR/BCR repertoire figures - V-J chord/circos, CDR3 spectratype, clonal-space stratification, clonal tracking across timepoints, rarefaction/extrapolatio… |
| [`tcr-bcr-analysis-scirpy-analysis`](#tcr-bcr-analysis-scirpy-analysis) | bioskills | Integrates single-cell paired TCR/BCR (10x VDJ, AIRR, dandelion, BD Rhapsody) with gene expression in an AnnData/MuData object using scirpy - chain-pairing QC,… |
| [`tcr-bcr-analysis-specificity-annotation`](#tcr-bcr-analysis-specificity-annotation) | bioskills | Maps TCR/BCR receptor sequences toward candidate antigen specificity and clusters repertoires by shared-specificity signal, while enforcing that a database mat… |
| [`tcr-bcr-analysis-vdjtools-analysis`](#tcr-bcr-analysis-vdjtools-analysis) | bioskills | Computes immune-repertoire diversity, clonal structure, overlap, and segment usage from TCR/BCR clonotype tables with VDJtools (immunarch as the modern R alter… |
| [`tumor-fraction-estimation`](#tumor-fraction-estimation) | bioskills | Estimates tumor fraction (the genome-wide proportion of cfDNA molecules that are tumor-derived, the cfDNA analogue of bulk-tumor purity) from shallow whole-gen… |
| [`variant-interpretation`](#variant-interpretation) | tooluniverse | Clinical variant interpretation from raw variant calls to ACMG-classified recommendations with structural impact analysis. Use for VUS classification, pathogen… |

## `analytical-validation` {#analytical-validation}

Treats a ctDNA assay as a molecule-counting experiment at the Poisson edge and builds its analytical-validation case the measurement-science way. Covers the genome-equivalent currency (~330 haploid copies/ng), the lambda = input_GE x VAF sampling ceiling (lambda&gt;=3 for ~95% detection), the error-suppression ladder (raw NGS ~1e-3 -&gt; single-strand UMI ~1e-4/1e-5 -&gt; duplex &lt;1e-7), the CLSI EP17 LoB/LoD/LoD95/LoQ framework, the per-locus-vs-panel-integrated LoD distinction that lets bespoke MRD reach ppm, contrived/SEQC2 reference standards, and honest LoD reporting conditioned on input mass + consensus depth + replicate detection rate. Use when stating or trusting a sensitivity claim, designing a dilution-series validation, deciding how many genome equivalents are needed at a target VAF, choosing a single-locus vs panel-integrated LoD, or auditing a "detects 0.1% VAF" claim.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/liquid-biopsy/analytical-validation)
- **License:** `MIT`
- **Runs with:** Python modules `dnakit`, `scipy`
- **Also known as:** `bio-analytical-validation`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info analytical-validation
dnakit skills install analytical-validation --agent claude-code   # or --agent codex / cursor
```

## `biomcp` {#biomcp}

Search and retrieve biomedical data - genes, variants, clinical trials, diagnostic tests, articles, drugs, diseases, pathways, proteins, adverse events, pharmacogenomics, and phenotype-disease matching. Use for gene function, variant pathogenicity, trials, diagnostics, drug safety, pathway context, disease workups, and literature evidence.

- **Source:** [GenomOncology and BioMCP contributors; see retained LICENSE/CITATION.cff @ `3b1d645`](https://github.com/genomoncology/biomcp/tree/3b1d6452ab63741562d2ef97495a97281b304729/skills)
- **License:** `MIT`
- **Runs with:** commands `biomcp`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info biomcp
dnakit skills install biomcp --agent claude-code   # or --agent codex / cursor
```

## `cancer-genomics-tcga` {#cancer-genomics-tcga}

TCGA/GDC cancer genomics analysis — cohort construction, clinical metadata retrieval, somatic mutation frequencies, survival analysis, and multi-omics integration. Use for TCGA-BRCA-style cohort studies, mutation prevalence by cancer type, survival-by-mutation analysis, and pan-cancer driver discovery. Always cancer-type-specific (don't use pan-cancer counts without cohort context).

- **Source:** [ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-cancer-genomics-tcga)
- **License:** `Apache-2.0`
- **Also known as:** `tooluniverse-cancer-genomics-tcga`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info cancer-genomics-tcga
dnakit skills install cancer-genomics-tcga --agent claude-code   # or --agent codex / cursor
```

## `cancer-variant-interpretation` {#cancer-variant-interpretation}

Clinical interpretation of somatic cancer mutations for precision oncology. Transforms a gene + variant + cancer-type input into an actionable report: clinical evidence tier (CIViC, OncoKB), therapeutic options (FDA-approved + investigational), resistance mechanisms, prognosis, and matching clinical trials. Use for tumor-board variant calls, somatic-mutation actionability assessment, and treatment selection. Always cancer-type-specific.

- **Source:** [ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-cancer-variant-interpretation)
- **License:** `Apache-2.0`
- **Runs with:** Python modules `dnakit`, `indigo`, `lxml`
- **Also known as:** `tooluniverse-cancer-variant-interpretation`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info cancer-variant-interpretation
dnakit skills install cancer-variant-interpretation --agent claude-code   # or --agent codex / cursor
```

## `cfdna-preprocessing` {#cfdna-preprocessing}

Decides how to preprocess plasma cfDNA sequencing data so the recoverable signal survives - library-prep-aware fragment expectations (dsDNA vs ssDNA/adaptase prep), UMI/duplex consensus with fgbio (ExtractUmisFromBam, GroupReadsByUmi --strategy paired for duplex, CallMolecularConsensusReads vs CallDuplexConsensusReads, FilterConsensusReads min-reads "total s1 s2"), the align-&gt;group-&gt;consensus-&gt;RE-align ordering, and the cfDNA dedup trap where naive coordinate dedup collapses nucleosome-coincident independent molecules. Covers when single-strand consensus suffices vs when duplex is mandatory, the singleton/sensitivity tax at low input, and reading the insert-size histogram as a pre-analytical QC instrument. Use when processing plasma cfDNA reads before fragmentomics, ctDNA mutation calling, or tumor-fraction estimation.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/liquid-biopsy/cfdna-preprocessing)
- **License:** `MIT`
- **Also known as:** `bio-cfdna-preprocessing`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info cfdna-preprocessing
dnakit skills install cfdna-preprocessing --agent claude-code   # or --agent codex / cursor
```

## `clinical-databases-acmg-classification` {#clinical-databases-acmg-classification}

Applies ACMG/AMP 2015 framework with ClinGen SVI specifications, Tavtigian 2018/2020 Bayesian point system, Abou Tayoun 2018 PVS1 decision tree, Pejaver 2022 and Bergquist 2025 calibrated PP3/BP4 thresholds for REVEL/BayesDel/AlphaMissense, Brnich 2020 PS3/BS3 OddsPath, Walker 2023 SpliceAI splicing framework, and AMP/ASCO/CAP 2017 tumor tiers. Use when classifying germline variants P / LP / VUS / LB / B, applying VCEP-specific CSpec rules, computing Whiffin BS1, or assigning cancer Tier I-IV per Li 2017.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/clinical-databases/acmg-classification)
- **License:** `MIT`
- **Runs with:** Python modules `dnakit`, `requests`
- **Also known as:** `bio-clinical-databases-acmg-classification`, `bio-variant-calling-clinical-interpretation`, `classify-sequence-variant`, `clinical-variant-reporter`, `genomic-variant-interpretation`, `tooluniverse-acmg-variant-classification`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info clinical-databases-acmg-classification
dnakit skills install clinical-databases-acmg-classification --agent claude-code   # or --agent codex / cursor
```

## `clinical-databases-clinvar-lookup` {#clinical-databases-clinvar-lookup}

Queries ClinVar for variant pathogenicity classifications, ClinGen VCEP curations, and somatic-vs-germline interpretations via REST API, weekly VCF, or bulk XML. Use when determining clinical significance, triangulating conflicting interpretations, or aggregating evidence against the ACMG/AMP framework with ClinGen SVI specifications.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/clinical-databases/clinvar-lookup)
- **License:** `MIT`
- **Runs with:** Python modules `dnakit`, `requests`
- **Also known as:** `alterlab-clinvar`, `bio-clinical-databases-clinvar-lookup`, `clinvar-annotation`, `clinvar-database`, `query-clinvar`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info clinical-databases-clinvar-lookup
dnakit skills install clinical-databases-clinvar-lookup --agent claude-code   # or --agent codex / cursor
```

## `clinical-databases-dbsnp-queries` {#clinical-databases-dbsnp-queries}

Resolves rsIDs, navigates RsMergeArch/SNPHistory merge chains, and converts between rsID, SPDI, HGVS, and VCF representations using the dbSNP Build 156 JSON architecture. Use when normalizing variant identifiers, joining variant databases by cluster ID, or tracking deprecated rsIDs through historical merges.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/clinical-databases/dbsnp-queries)
- **License:** `MIT`
- **Runs with:** Python modules `dnakit`, `myvariant`
- **Also known as:** `bio-clinical-databases-dbsnp-queries`, `dbsnp-database`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info clinical-databases-dbsnp-queries
dnakit skills install clinical-databases-dbsnp-queries --agent claude-code   # or --agent codex / cursor
```

## `clinical-databases-gnomad-frequencies` {#clinical-databases-gnomad-frequencies}

Queries gnomAD v4 (807k samples), v3, v2.1.1, and constraint metrics with grpmax FAF95, bottleneck-group exclusion, LOEUF interpretation, SV/CNV/mtDNA catalogs, and Whiffin max-credible-AF framework. Use when filtering rare variants, applying ACMG BS1/BA1, ranking genes by LoF intolerance, or selecting between v2 (GRCh37 + chrX/Y constraint) and v4 (GRCh38 + 807k samples).

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/clinical-databases/gnomad-frequencies)
- **License:** `MIT`
- **Runs with:** Python modules `dnakit`, `requests`
- **Also known as:** `alterlab-gnomad`, `bio-clinical-databases-gnomad-frequencies`, `gnomad-database`, `gnomad-variants`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info clinical-databases-gnomad-frequencies
dnakit skills install clinical-databases-gnomad-frequencies --agent claude-code   # or --agent codex / cursor
```

## `clinical-databases-hla-typing` {#clinical-databases-hla-typing}

Calls HLA class I and class II alleles at 2/4/6/8-field resolution from WGS/WES/RNA-seq/long-read data using OptiType, HLA-LA, T1K, Polysolver, HLA-HD, arcasHLA, StarPhase, or HIBAG imputation. Use when typing for HSCT, solid-organ transplant, neoantigen prediction, PGx screening (B*57:01, B*15:02, etc.), or disease-association studies, with reconciliation across tools and IPD-IMGT/HLA version mismatch handling.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/clinical-databases/hla-typing)
- **License:** `MIT`
- **Runs with:** commands `T1K`
- **Also known as:** `bio-clinical-databases-hla-typing`, `hla-typing`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info clinical-databases-hla-typing
dnakit skills install clinical-databases-hla-typing --agent claude-code   # or --agent codex / cursor
```

## `clinical-databases-msi-detection` {#clinical-databases-msi-detection}

Calls microsatellite instability from WES/WGS/targeted-panel with MSIsensor, MSIsensor-pro, MSIsensor-ct (panel-aware), mSINGS, and MANTIS for FDA pembrolizumab MSI-H pan-tumor / Lynch syndrome / dMMR ICI biomarker. Use when stratifying ICI eligibility (Le 2015), pairing MSI with TMB-H (Sha 2020 / Salem 2018), screening Lynch syndrome (universal IHC + MSI), or distinguishing MSI-H tumors from POLE-exo hypermutator with overlapping signatures.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/clinical-databases/msi-detection)
- **License:** `MIT`
- **Runs with:** commands `MSIsensor-pro`
- **Also known as:** `bio-clinical-databases-msi-detection`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info clinical-databases-msi-detection
dnakit skills install clinical-databases-msi-detection --agent claude-code   # or --agent codex / cursor
```

## `clinical-databases-myvariant-queries` {#clinical-databases-myvariant-queries}

Queries myvariant.info BioThings aggregator for ClinVar, gnomAD, dbSNP, dbNSFP, COSMIC, CADD, and CIViC annotations in batched, version-tracked requests. Use when annotating variant lists from multiple databases simultaneously without managing per-source APIs, and when reproducibility-grade analyses require recording source data versions via _meta.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/clinical-databases/myvariant-queries)
- **License:** `MIT`
- **Runs with:** Python modules `dnakit`, `myvariant`
- **Also known as:** `bio-clinical-databases-myvariant-queries`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info clinical-databases-myvariant-queries
dnakit skills install clinical-databases-myvariant-queries --agent claude-code   # or --agent codex / cursor
```

## `clinical-databases-pharmacogenomics` {#clinical-databases-pharmacogenomics}

Queries PharmGKB / CPIC / DPWG for drug-gene interactions; calls CYP2D6/CYP2C9/CYP2C19/DPYD/TPMT/NUDT15/UGT1A1/SLCO1B1 star alleles and phenotype with PharmCAT, Cyrius (CYP2D6 structural variants), Aldy, Stargazer; applies Caudle 2020 activity-score translation. Use when implementing pharmacogenomic-guided prescribing, applying CPIC vs DPWG guidance, screening HLA risk alleles for ICI / antiepileptics / abacavir, or interpreting compound TPMT+NUDT15 thiopurine risk.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/clinical-databases/pharmacogenomics)
- **License:** `MIT`
- **Also known as:** `alterlab-clinpgx`, `bio-clinical-databases-pharmacogenomics`, `pharmgx-reporter`, `tooluniverse-pharmacogenomics`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info clinical-databases-pharmacogenomics
dnakit skills install clinical-databases-pharmacogenomics --agent claude-code   # or --agent codex / cursor
```

## `clinical-databases-polygenic-risk` {#clinical-databases-polygenic-risk}

Constructs and validates polygenic risk scores using LDpred2-auto, SBayesRC, MegaPRS, PRS-CS, PROSPER, MUSSEL, BridgePRS, JointPRS, PRSmix, or PGS Catalog Calculator with ancestry-aware reference panels (HapMap3, UKB-LD), ancestry-conditional calibration, and PRS-RS reporting standards. Use when computing PRS for cohorts, applying absolute-risk transformation, assessing cross-ancestry portability (Martin 2017 / Ding 2023 continuous ancestry), or auditing PRS manuscripts against the 22-item PRS-RS reviewer checklist.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/clinical-databases/polygenic-risk)
- **License:** `MIT`
- **Also known as:** `bio-clinical-databases-polygenic-risk`, `gwas-prs`, `just-prs-mcp`, `polygenic-risk-score-prs-catalog`, `tooluniverse-polygenic-risk-score`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info clinical-databases-polygenic-risk
dnakit skills install clinical-databases-polygenic-risk --agent claude-code   # or --agent codex / cursor
```

## `clinical-databases-somatic-signatures` {#clinical-databases-somatic-signatures}

Extracts and assigns COSMIC v3.4 mutational signatures (86 SBS / 11 DBS / 18 ID / 21 CN / 16 SV) from somatic VCFs using SigProfilerSuite, MutationalPatterns, MuSiCal mvNMF, SigNet, or HRDetect. Use when characterizing DNA-damage etiology (BRCA1/2 HRD, MMR-D, POLE, APOBEC3A, UV, tobacco, aflatoxin, 5-FU/SBS17b, platinum, colibactin SBS88), routing PARP inhibitor decisions, or auditing de novo extraction vs refit choice for cohort size.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/clinical-databases/somatic-signatures)
- **License:** `MIT`
- **Also known as:** `bio-clinical-databases-somatic-signatures`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info clinical-databases-somatic-signatures
dnakit skills install clinical-databases-somatic-signatures --agent claude-code   # or --agent codex / cursor
```

## `clinical-databases-tumor-mutational-burden` {#clinical-databases-tumor-mutational-burden}

Calculates tumor mutational burden from WES/WGS/panel data with Friends of Cancer Research harmonization equations, per-assay calibration (FDA 10/Mb = 7.8 TSO500 = 8.4 OncomineTML), synonymous/indel/germline filtering, hypermutator tiering, blood TMB, and integration with HLA-LOH and neoantigen quality (Luksza 2017 fitness). Use when assessing ICI eligibility under tumor-specific cutoffs (McGrail 2021), comparing tissue vs bTMB, or auditing TMB-H reporting against ESMO 2024 and FDA pembrolizumab pan-tumor 2020.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/clinical-databases/tumor-mutational-burden)
- **License:** `MIT`
- **Runs with:** Python modules `cyvcf2`, `dnakit`
- **Also known as:** `bio-clinical-databases-tumor-mutational-burden`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info clinical-databases-tumor-mutational-burden
dnakit skills install clinical-databases-tumor-mutational-burden --agent claude-code   # or --agent codex / cursor
```

## `clinical-databases-variant-prioritization` {#clinical-databases-variant-prioritization}

Prioritizes rare-disease variants from trio/quad WES/WGS with de novo (DeNovoGear, Triodenovo), compound-heterozygous phasing (WhatsHap), mosaic VAF tiering, phenotype-driven ranking (Exomiser, Phen2Gene, AMELIE), ClinGen gene-disease validity gating, and ACMG SF v3.2 secondary findings reporting. Use when running diagnostic exome / genome pipelines, identifying candidate Mendelian disease genes, screening for incidental findings, or auditing VUS reclassification cycles. The ACMG/AMP classification framework (PVS1 decision tree, Pejaver PP3/BP4 calibration, Tavtigian point system) is in clinical-databases/acmg-classification.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/clinical-databases/variant-prioritization)
- **License:** `MIT`
- **Runs with:** Python modules `dnakit`, `pandas`
- **Also known as:** `bio-clinical-databases-variant-prioritization`, `clinical-variant-prioritizer`, `filter-variant-callset`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info clinical-databases-variant-prioritization
dnakit skills install clinical-databases-variant-prioritization --agent claude-code   # or --agent codex / cursor
```

## `clinpgx` {#clinpgx}

Query the ClinPGx API for pharmacogenomic gene-drug data, clinical annotations, CPIC guidelines, and FDA drug labels

- **Source:** [ClawBio contributors @ `1fcecb7`](https://github.com/ClawBio/ClawBio/tree/1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d/skills/clinpgx)
- **License:** `MIT`
- **Runs with:** commands `python3`
- **Also known as:** `clinpgx-database`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info clinpgx
dnakit skills install clinpgx --agent claude-code   # or --agent codex / cursor
```

## `ctdna-mutation-detection` {#ctdna-mutation-detection}

Detects somatic mutations in circulating tumor DNA, treating low-VAF detection as a signal-versus-noise problem set by error suppression and molecules sampled, not by the choice of caller. Distinguishes de novo CALLING (scanning a panel for unknown variants, bounded by per-locus error and multiple testing) from tumor-informed DETECTION (tracking a pre-specified variant set, where panel integration reaches single-ppm). Covers VarDict and Mutect2 for de novo calling, UMI-aware callers, and a pysam-based known-variant VAF tracker, with matched-WBC subtraction as the mandatory defense against clonal hematopoiesis (the dominant false positive). Use when calling or tracking tumor mutations from plasma cfDNA, setting a VAF threshold, or deciding whether a low-VAF call is tumor versus CHIP.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/liquid-biopsy/ctdna-mutation-detection)
- **License:** `MIT`
- **Also known as:** `bio-ctdna-mutation-detection`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info ctdna-mutation-detection
dnakit skills install ctdna-mutation-detection --agent claude-code   # or --agent codex / cursor
```

## `data-visualization-oncoprint-mutation-matrices` {#data-visualization-oncoprint-mutation-matrices}

Build OncoPrint and co-mutation matrix plots from somatic-variant cohorts using ComplexHeatmap, maftools, and comut.py with alteration-type stacking, sample ordering by mutational burden, mutual-exclusivity overlays, and clinical annotation tracks. Use when visualizing per-sample mutation patterns across recurrent driver genes, comparing alteration classes, or identifying mutually-exclusive / co-occurring driver pairs.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/data-visualization/oncoprint-mutation-matrices)
- **License:** `MIT`
- **Also known as:** `bio-data-visualization-oncoprint-mutation-matrices`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info data-visualization-oncoprint-mutation-matrices
dnakit skills install data-visualization-oncoprint-mutation-matrices --agent claude-code   # or --agent codex / cursor
```

## `draft-counseling-summary` {#draft-counseling-summary}

Draft clinician-facing and optional patient-facing rare-disease germline genetic counseling summaries from reviewed case facts and interpretation proposals. Use when asked to summarize a positive, VUS, or negative result; explain inheritance or family implications; prepare counseling notes; communicate testing limitations; or communicate reviewed follow-up recommendations for professional review.

- **Source:** [muwenbo/genetic-counseling-skills contributors @ `eda3cc5`](https://github.com/muwenbo/genetic-counseling-skills/tree/eda3cc5dace5ebed128dd14fa0f12e5bf2e16ff4/plugins/rare-disease-genetic-counseling/skills/draft-counseling-summary)
- **License:** `MIT`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info draft-counseling-summary
dnakit skills install draft-counseling-summary --agent claude-code   # or --agent codex / cursor
```

## `draw-pedigree` {#draw-pedigree}

Create a reviewable clinical pedigree figure or validated PED/JSON file from user-supplied family-history facts. Use when asked to draw, chart, or diagram a pedigree or family tree; visualize inheritance or segregation across relatives; render affection status, carriers, deceased individuals, proband, twins, consanguinity, or variant genotypes; or convert a PED/family history into a figure—without inferring missing relationships, phenotypes, or genotypes.

- **Source:** [muwenbo/genetic-counseling-skills contributors @ `eda3cc5`](https://github.com/muwenbo/genetic-counseling-skills/tree/eda3cc5dace5ebed128dd14fa0f12e5bf2e16ff4/plugins/rare-disease-genetic-counseling/skills/draw-pedigree)
- **License:** `MIT`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info draw-pedigree
dnakit skills install draw-pedigree --agent claude-code   # or --agent codex / cursor
```

## `find-variant-literature` {#find-variant-literature}

Discover, retrieve, and extract claim-level literature evidence for a rare-disease germline variant, gene-disease pair, ACMG criterion, or CNV region. Use when asked to find papers, search PubMed, fetch known PMIDs, review functional or segregation evidence, count reported cases, investigate de novo/phase evidence, or support a variant interpretation with literature.

- **Source:** [muwenbo/genetic-counseling-skills contributors @ `eda3cc5`](https://github.com/muwenbo/genetic-counseling-skills/tree/eda3cc5dace5ebed128dd14fa0f12e5bf2e16ff4/plugins/rare-disease-genetic-counseling/skills/find-variant-literature)
- **License:** `MIT`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info find-variant-literature
dnakit skills install find-variant-literature --agent claude-code   # or --agent codex / cursor
```

## `folklore-variant-evidence` {#folklore-variant-evidence}

Retrieve ClinGen gene-disease validity assertions for a public gene or disease, and review source-linked public evidence and literature for one supported GRCh38 germline nuclear SNV or simple indel through Folklore Clinical Variant Interpretation MCP. Use when a scientific agent must branch deterministically on resolved, ambiguous, not-found, invalid, unsupported, or unavailable variant outcomes; chain a resolved public variant into related literature or publication details; or preserve evidence provenance without accepting patient, phenotype, family, segregation, or private case data.

- **Source:** [Helena Bioinformatics; K-Dense Inc. and contributors @ `c1ed16d`](https://github.com/K-Dense-AI/scientific-agent-skills/tree/c1ed16d97dd61ff50a3bd46dd353e4a55fd77f34/skills/folklore-variant-evidence)
- **License:** `MIT`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info folklore-variant-evidence
dnakit skills install folklore-variant-evidence --agent claude-code   # or --agent codex / cursor
```

## `fragment-analysis` {#fragment-analysis}

Extracts cfDNA fragmentomics features (DELFI genome-wide short/long ratios, WPS nucleosome positioning, Griffin GC-corrected accessibility profiles, end-motifs/MDS, OCF) for cancer detection and tissue-of-origin from plasma WGS. Centers on the nuclease-footprint reframe (every feature re-reads one nucleosome object), the mandatory GC correction, and the cross-protocol non-comparability that breaks naive classifiers. Runs FinaleToolkit (real CLI/Python, MIT) and the Griffin Snakemake pipeline; DELFI is a method, not a package. Use when deriving fragment-based signal from cfDNA, choosing a feature family for detection vs subtyping, or diagnosing why a fragmentomic model failed validation.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/liquid-biopsy/fragment-analysis)
- **License:** `MIT`
- **Runs with:** Python modules `FinaleToolkit`, `dnakit`
- **Also known as:** `bio-fragment-analysis`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info fragment-analysis
dnakit skills install fragment-analysis --agent claude-code   # or --agent codex / cursor
```

## `gene-disease-association` {#gene-disease-association}

Gene-disease association analysis across DisGeNET, OpenTargets, Monarch, OMIM, GenCC, Orphanet. Cross-references multiple sources for evidence-graded association reports with concordance scoring (5/5 sources agree → strong, 1/5 → weak). Use for 'which diseases is gene X associated with' or 'which genes cause disease Y' queries with quantitative confidence.

- **Source:** [ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-gene-disease-association)
- **License:** `Apache-2.0`
- **Also known as:** `tooluniverse-gene-disease-association`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info gene-disease-association
dnakit skills install gene-disease-association --agent claude-code   # or --agent codex / cursor
```

## `hla-immunogenomics` {#hla-immunogenomics}

HLA gene-family analysis and MHC-peptide binding for transplant compatibility, vaccine epitope coverage, and cancer immunotherapy. Uses IMGT (HLA polymorphism), IEDB (epitope-MHC binding), UniProt (annotation), DGIdb (druggability). Use for HLA typing/imputation review, vaccine HLA coverage, and immunotherapy prediction biomarkers (HLA-LOH, neoantigen presentation).

- **Source:** [ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-hla-immunogenomics)
- **License:** `Apache-2.0`
- **Also known as:** `tooluniverse-hla-immunogenomics`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info hla-immunogenomics
dnakit skills install hla-immunogenomics --agent claude-code   # or --agent codex / cursor
```

## `immune-repertoire-analysis` {#immune-repertoire-analysis}

TCR/BCR repertoire analysis — V(D)J segment usage, CDR3 sequence diversity, clonality scoring, antigen specificity matching to IEDB, public-clone identification. Use for adaptive immune response characterization, post-treatment immune monitoring, antigen-specific clone tracking, and clonal-expansion analysis in immunotherapy or vaccination studies.

- **Source:** [ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-immune-repertoire-analysis)
- **License:** `Apache-2.0`
- **Also known as:** `tooluniverse-immune-repertoire-analysis`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info immune-repertoire-analysis
dnakit skills install immune-repertoire-analysis --agent claude-code   # or --agent codex / cursor
```

## `immunoinformatics-epitope-prediction` {#immunoinformatics-epitope-prediction}

Predict B-cell and T-cell epitopes for vaccine antigen design and epitope mapping with BepiPred-3.0, DiscoTope-3.0, the IEDB tools, and EL-mode MHC presentation. Encodes the load-bearing asymmetry that T-cell epitope prediction is mature (it reduces to MHC presentation, AUC&gt;0.9) while B-cell prediction is unreliable (linear predictors ~AUC 0.6 because ~90% of real epitopes are conformational) — so structure-based DiscoTope-3.0 on AlphaFold models is the only defensible B-cell path, propensity scales are obsolete, and NetChop is largely redundant on EL-trained models. Use when mapping epitopes or selecting vaccine antigens. MHC binding lives in mhc-binding-prediction.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/immunoinformatics/epitope-prediction)
- **License:** `MIT`
- **Runs with:** Python modules `BepiPred`, `dnakit`
- **Also known as:** `bio-immunoinformatics-epitope-prediction`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info immunoinformatics-epitope-prediction
dnakit skills install immunoinformatics-epitope-prediction --agent claude-code   # or --agent codex / cursor
```

## `immunoinformatics-immunogenicity-scoring` {#immunoinformatics-immunogenicity-scoring}

Rank and prioritize neoantigen/epitope candidates by likely T-cell response using NeoFox feature annotation, PRIME2.0, BigMHC-IM, the Łuksza/Balachandran fitness model (agretopicity + foreignness), and pVACtools tiering. Encodes the field's hard truths that immunogenicity is the least-solved layer (dedicated scores ~AUROC 0.6-0.7, modest PPV), that scores are valid only for RANKING within one patient (never absolute go/no-go or cross-patient), that DAI has anchor-inflation and WT-denominator traps, and that stacking weak correlated scores into one number is a red flag. Use when ordering a candidate list for a vaccine. Binding lives in mhc-binding-prediction; calling in neoantigen-prediction.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/immunoinformatics/immunogenicity-scoring)
- **License:** `MIT`
- **Runs with:** Python modules `NeoFox`, `dnakit`
- **Also known as:** `bio-immunoinformatics-immunogenicity-scoring`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info immunoinformatics-immunogenicity-scoring
dnakit skills install immunoinformatics-immunogenicity-scoring --agent claude-code   # or --agent codex / cursor
```

## `immunoinformatics-mhc-binding-prediction` {#immunoinformatics-mhc-binding-prediction}

Predict peptide-MHC class I binding and natural presentation with MHCflurry, NetMHCpan-4.1, and MixMHCpred to nominate candidate CD8 T-cell epitopes. Covers the binding-affinity (BA) vs eluted-ligand (EL/presentation) distinction, why %Rank beats raw nM for cross-allele work, the MS abundance bias that misranks low-expression neoantigens, allele-coverage inequity, and length bias. Use when scanning a protein or peptide set for class I epitopes, scoring neoantigen candidates, or choosing a binding predictor. For CD4/HLA class II see mhc-class-ii-prediction.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/immunoinformatics/mhc-binding-prediction)
- **License:** `MIT`
- **Runs with:** Python modules `dnakit`, `mhcflurry`
- **Also known as:** `bio-immunoinformatics-mhc-binding-prediction`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info immunoinformatics-mhc-binding-prediction
dnakit skills install immunoinformatics-mhc-binding-prediction --agent claude-code   # or --agent codex / cursor
```

## `immunoinformatics-mhc-class-ii-prediction` {#immunoinformatics-mhc-class-ii-prediction}

Predict peptide-MHC class II (HLA-DR/DQ/DP) binding and presentation for CD4 T-cell epitopes with NetMHCIIpan-4.3 and MixMHC2pred-2.0. Covers why class II is far less reliable than class I (open binding groove, 9-mer register ambiguity, sparse noisy training data, DR&gt;DP&gt;DQ accuracy asymmetry), the DQ/DP heterodimer alpha/beta pairing trap, and the looser 1%/5% %Rank thresholds. Use when predicting CD4 epitopes for vaccine help, mapping class II neoantigens, or scoring long peptides against DR/DQ/DP. For CD8/class I see mhc-binding-prediction.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/immunoinformatics/mhc-class-ii-prediction)
- **License:** `MIT`
- **Runs with:** commands `NetMHCIIpan`
- **Also known as:** `bio-immunoinformatics-mhc-class-ii-prediction`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info immunoinformatics-mhc-class-ii-prediction
dnakit skills install immunoinformatics-mhc-class-ii-prediction --agent claude-code   # or --agent codex / cursor
```

## `immunoinformatics-neoantigen-prediction` {#immunoinformatics-neoantigen-prediction}

Identify tumor neoantigens from somatic variants with pVACtools (pVACseq/pVACfuse/pVACbind/pVACvector/pVACview) for personalized cancer vaccines and checkpoint biomarkers. Encodes the field's hard truth that binding prediction is the easy, near-solved part and single-digit-percent PPV lives downstream — so it centers clonality/CCF, HLA LOH (the silent invalidator), expression, proximal-variant phasing, agretopicity/foreignness quality, and the predicted-&gt;presented-&gt;immunogenic validation tiers. Use when nominating vaccine targets, ranking neoantigens, or building a tumor-to-candidate pipeline. Binding details in mhc-binding-prediction; ranking in immunogenicity-scoring.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/immunoinformatics/neoantigen-prediction)
- **License:** `MIT`
- **Also known as:** `bio-immunoinformatics-neoantigen-prediction`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info immunoinformatics-neoantigen-prediction
dnakit skills install immunoinformatics-neoantigen-prediction --agent claude-code   # or --agent codex / cursor
```

## `interpret-negative-or-vus` {#interpret-negative-or-vus}

Interpret a completed negative or VUS rare-disease germline analysis by checking QC and assay blind spots, generating restrained phenotype-specific hypotheses, and prioritizing orthogonal confirmatory tests. Use when asked what WES/WGS/panel testing may have missed, why a case is negative, what to test next, or how to follow up an inconclusive result.

- **Source:** [muwenbo/genetic-counseling-skills contributors @ `eda3cc5`](https://github.com/muwenbo/genetic-counseling-skills/tree/eda3cc5dace5ebed128dd14fa0f12e5bf2e16ff4/plugins/rare-disease-genetic-counseling/skills/interpret-negative-or-vus)
- **License:** `MIT`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info interpret-negative-or-vus
dnakit skills install interpret-negative-or-vus --agent claude-code   # or --agent codex / cursor
```

## `longitudinal-monitoring` {#longitudinal-monitoring}

Tracks ctDNA across serial liquid-biopsy timepoints for molecular residual disease (MRD) and treatment-response monitoring, treating MRD as a binary integrated detection call across the patient's full variant set (with a defined LoD95 and per-sample specificity) rather than a per-timepoint VAF threshold, and handling undetectable samples as left-censored at the per-sample limit of detection rather than true zeros. Covers tumor-informed bespoke vs tumor-naive design, landmark vs surveillance sampling, molecular-response definitions and their non-standardization, censoring-aware clearance kinetics, and the multiple-testing structure of repeated surveillance. Use when monitoring ctDNA during therapy, calling molecular relapse before imaging, or estimating clearance half-life from serial samples.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/liquid-biopsy/longitudinal-monitoring)
- **License:** `MIT`
- **Runs with:** Python modules `dnakit`, `pandas`
- **Also known as:** `bio-longitudinal-monitoring`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info longitudinal-monitoring
dnakit skills install longitudinal-monitoring --agent claude-code   # or --agent codex / cursor
```

## `methylation-based-detection` {#methylation-based-detection}

Detects cancer and infers tissue-of-origin from cfDNA methylation by choosing conversion chemistry (bisulfite vs EM-seq vs TAPS vs cfMeDIP), calling read-level methylation haplotypes rather than averaged beta values, and deconvolving a hematopoietic-dominated cfDNA mixture against a methylation atlas via NNLS/quadratic programming. Encodes the GRAIL/CCGA thesis that thousands of tissue-specific markers make methylation outperform sparse mutations for multi-cancer early detection (MCED) and localization, and that single concordantly-methylated fragments give ppm-level sensitivity. Uses MethylDackel for extraction (mbias-then-extract), MEDIPS/QSEA for enrichment data, scipy.optimize.nnls for deconvolution. Use when building an MCED or methylation-MRD assay, picking a conversion chemistry for low-input plasma, or deconvolving tissue-of-origin from cfDNA.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/liquid-biopsy/methylation-based-detection)
- **License:** `MIT`
- **Also known as:** `bio-methylation-based-detection`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info methylation-based-detection
dnakit skills install methylation-based-detection --agent claude-code   # or --agent codex / cursor
```

## `multiomic-disease-characterization` {#multiomic-disease-characterization}

Comprehensive disease characterization across genomics, transcriptomics, proteomics, and pathways for systems-level understanding. Identifies therapeutic opportunities and biomarker candidates by integrating multi-layer molecular data. Use for full-omics disease deep-dive reports, mechanism mapping, and biomarker-and-target identification from multi-omics data.

- **Source:** [ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-multiomic-disease-characterization)
- **License:** `Apache-2.0`
- **Runs with:** Python modules `dnakit`, `networkx`
- **Also known as:** `tooluniverse-multiomic-disease-characterization`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info multiomic-disease-characterization
dnakit skills install multiomic-disease-characterization --agent claude-code   # or --agent codex / cursor
```

## `ngs-tertiary-analysis` {#ngs-tertiary-analysis}

Run NGS tertiary analysis pipeline from BAM/VCF to ESMO-compliant clinical HTML report. Use when user has a cancer panel BAM/VCF file (TSO500) and needs somatic variant calling, annotation, OncoKB + CiVIC clinical actionability, ESCAT + AMP/ASCO/CAP classification, interactive visualizations, literature evidence, and a publication-quality genomics report. Two modes: /ngs-quick (one-shot) or /ngs-interactive (step-by-step with choices).

- **Source:** [htlin222/ngs-tertiary-analysis-skills contributors @ `d483c21`](https://github.com/htlin222/ngs-tertiary-analysis-skills/tree/d483c21fddd14a034e7f2f038f237ad7bcd05b1a/.claude/skills/ngs-tertiary-analysis)
- **License:** `MIT`
- **Runs with:** commands `Rscript`
- **Also known as:** `ngs-interactive`, `ngs-quick`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info ngs-tertiary-analysis
dnakit skills install ngs-tertiary-analysis --agent claude-code   # or --agent codex / cursor
```

## `prepare-rare-disease-case` {#prepare-rare-disease-case}

Structure English or Chinese rare-disease germline clinical narratives into reviewable findings, explicit negatives, suspected diagnoses, treatments, confounders, missing information, and candidate HPO terms. Use when preparing a case, extracting phenotypes, reviewing a referral note, normalizing clinical findings, or assembling phenotype input before disease or variant analysis.

- **Source:** [muwenbo/genetic-counseling-skills contributors @ `eda3cc5`](https://github.com/muwenbo/genetic-counseling-skills/tree/eda3cc5dace5ebed128dd14fa0f12e5bf2e16ff4/plugins/rare-disease-genetic-counseling/skills/prepare-rare-disease-case)
- **License:** `MIT`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info prepare-rare-disease-case
dnakit skills install prepare-rare-disease-case --agent claude-code   # or --agent codex / cursor
```

## `protein-lof-mechanism` {#protein-lof-mechanism}

Propose the mechanism by which a missense variant causes loss-of-function (LoF), synthesizing evidence from 5 independent layers: AlphaMissense pathogenicity, AlphaFold structural context, ESMC sequence likelihood, SAE feature disruption, and DynaMut2 stability ΔΔG. Distinguishes 'structural stability LoF' (mis-folding) from 'direct functional disruption' (catalytic / binding / PTM site damage). Use for coding missense variants where you need a mechanistic causal model, not just a pathogenicity score.

- **Source:** [ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-protein-lof-mechanism)
- **License:** `Apache-2.0`
- **Also known as:** `tooluniverse-protein-lof-mechanism`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info protein-lof-mechanism
dnakit skills install protein-lof-mechanism --agent claude-code   # or --agent codex / cursor
```

## `protein-sae-variant-interpretation` {#protein-sae-variant-interpretation}

Interpret a missense variant via ESMC-6B Sparse Autoencoder (SAE) feature activations. For a given protein + variant, computes which interpretable SAE features (catalytic, ligand-binding, PTM, structural motif, domain, etc.) are lost or gained at the mutation site. Use when standard pathogenicity scores (AlphaMissense, ClinVar) say a variant is damaging but you need a MECHANISTIC explanation — e.g. 'why is this variant LoF?' Complements (does not replace) variant-interpretation and variant-to-mechanism skills, which focus on ACMG classification or regulatory mechanism.

- **Source:** [ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-protein-sae-variant-interpretation)
- **License:** `Apache-2.0`
- **Also known as:** `tooluniverse-protein-sae-variant-interpretation`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info protein-sae-variant-interpretation
dnakit skills install protein-sae-variant-interpretation --agent claude-code   # or --agent codex / cursor
```

## `rare-disease-diagnosis` {#rare-disease-diagnosis}

Rare disease differential diagnosis from patient phenotype — HPO term matching to candidate diseases (Orphanet, OMIM), gene panel prioritization, ACMG variant interpretation, and structure-based variant analysis. Use for diagnostic odyssey assistance, phenotype-to-disease ranking, and genetic-counseling differential generation.

- **Source:** [ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-rare-disease-diagnosis)
- **License:** `Apache-2.0`
- **Runs with:** Python modules `cellxgene_census`, `dnakit`
- **Also known as:** `rank-disease-candidates`, `tooluniverse-rare-disease-diagnosis`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info rare-disease-diagnosis
dnakit skills install rare-disease-diagnosis --agent claude-code   # or --agent codex / cursor
```

## `rare-disease-genomics` {#rare-disease-genomics}

Rare disease genomics — disease identification (Orphanet), causative gene discovery, gene-disease validity (GenCC), variant interpretation (ClinVar), and translational research (ClinicalTrials.gov, drug repurposing for orphans). Use for rare-disease-gene curation, novel-gene-discovery analysis, and rare-disease drug-development support.

- **Source:** [ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-rare-disease-genomics)
- **License:** `Apache-2.0`
- **Also known as:** `tooluniverse-rare-disease-genomics`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info rare-disease-genomics
dnakit skills install rare-disease-genomics --agent claude-code   # or --agent codex / cursor
```

## `rare-high-impact-variants` {#rare-high-impact-variants}

Count rare, high-impact loss-of-function variants carried in a VCF, annotated with molecular consequence and population allele frequency

- **Source:** [Manuel Corpas; ClawBio contributors @ `1fcecb7`](https://github.com/ClawBio/ClawBio/tree/1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d/skills/rare-high-impact-variants)
- **License:** `MIT`
- **Runs with:** commands `python3`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info rare-high-impact-variants
dnakit skills install rare-high-impact-variants --agent claude-code   # or --agent codex / cursor
```

## `review-variant-table` {#review-variant-table}

Inspect and normalize testing-company Excel, CSV, or TSV germline variant exports into a reviewable variant inventory. Use when a user uploads a laboratory variant table and wants to identify the data sheet or header, map Chinese or English columns, audit company annotations, verify table structure, shortlist variants, or prepare selected records for disease ranking or sequence/CNV classification.

- **Source:** [muwenbo/genetic-counseling-skills contributors @ `eda3cc5`](https://github.com/muwenbo/genetic-counseling-skills/tree/eda3cc5dace5ebed128dd14fa0f12e5bf2e16ff4/plugins/rare-disease-genetic-counseling/skills/review-variant-table)
- **License:** `MIT`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info review-variant-table
dnakit skills install review-variant-table --agent claude-code   # or --agent codex / cursor
```

## `tcr-bcr-analysis-immcantation-analysis` {#tcr-bcr-analysis-immcantation-analysis}

Reconstructs B-cell clonal families, quantifies somatic hypermutation and selection, and builds antibody lineage trees with the Immcantation R suite (alakazam, shazam, scoper, dowser, tigger) on AIRR-format BCR data. Use when deriving the clonal-clustering threshold from the distToNearest bimodal valley (never a hardcoded 0.15); choosing hierarchicalClones vs spectralClones (vj vs novj) for SHM-diverged repertoires; personalizing the germline with TIGGER before mutation counting; reconstructing D-masked germlines with createGermlines; measuring R/S mutation frequency by CDR and FWR region; testing antigen-driven selection with BASELINe; comparing Hill-number diversity at equal sampling depth; and inferring IgPhyML lineage trees for affinity maturation, class-switch, and ancestral-antibody analysis.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/tcr-bcr-analysis/immcantation-analysis)
- **License:** `MIT`
- **Also known as:** `bio-tcr-bcr-analysis-immcantation-analysis`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info tcr-bcr-analysis-immcantation-analysis
dnakit skills install tcr-bcr-analysis-immcantation-analysis --agent claude-code   # or --agent codex / cursor
```

## `tcr-bcr-analysis-mixcr-analysis` {#tcr-bcr-analysis-mixcr-analysis}

Align V(D)J reads and assemble TCR/BCR clonotypes with MiXCR, driven by a chemistry-matched preset. Use when choosing/auditing the preset for a library (5'RACE/template-switch vs multiplex-primer amplicon -&gt; rigid vs floating boundaries; RNA vs gDNA -&gt; --rna/--dna; bulk vs 10x single-cell; UMI vs no-UMI -&gt; tag pattern and barcode collapse; kit presets Takara/NEBNext/QIAseq/BD/MiLaboratory); assembling clonotypes by CDR3 vs VDJRegion; setting the reads-vs-UMI-vs-cell quantitation denominator; exporting native MiXCR fields vs AIRR rearrangement TSV for downstream Immcantation/scirpy/VDJtools; and running alignment/chain-usage QC. Keywords: MiXCR, analyze, align, refineTagsAndSort, assemblePartial, assemble, assembleCells, exportClones, exportAirr, exportQc, CDR3, V(D)J, clonotype, UMI, cell barcode, 10x VDJ, license.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/tcr-bcr-analysis/mixcr-analysis)
- **License:** `MIT`
- **Runs with:** commands `MiXCR`
- **Also known as:** `bio-tcr-bcr-analysis-mixcr-analysis`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info tcr-bcr-analysis-mixcr-analysis
dnakit skills install tcr-bcr-analysis-mixcr-analysis --agent claude-code   # or --agent codex / cursor
```

## `tcr-bcr-analysis-repertoire-visualization` {#tcr-bcr-analysis-repertoire-visualization}

Draws TCR/BCR repertoire figures - V-J chord/circos, CDR3 spectratype, clonal-space stratification, clonal tracking across timepoints, rarefaction/extrapolation curves, overlap heatmaps, and clonotype-similarity networks - and encodes how to read them. Use when choosing between a raw Shannon bar and a rarefaction curve for a diversity comparison; deciding a depth-robust overlap metric (Morisita-Horn) vs a set metric (Jaccard) for a heatmap; setting the distance threshold that defines a clonotype-similarity network; interpreting a Gaussian vs skewed spectratype as polyclonal vs clonally expanded; or laying out clonal-space and clone-tracking plots. Covers VDJtools PlotFancyVJUsage/RarefactionPlot, R circlize and iNEXT, and matplotlib/seaborn recipes.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/tcr-bcr-analysis/repertoire-visualization)
- **License:** `MIT`
- **Also known as:** `bio-tcr-bcr-analysis-repertoire-visualization`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info tcr-bcr-analysis-repertoire-visualization
dnakit skills install tcr-bcr-analysis-repertoire-visualization --agent claude-code   # or --agent codex / cursor
```

## `tcr-bcr-analysis-scirpy-analysis` {#tcr-bcr-analysis-scirpy-analysis}

Integrates single-cell paired TCR/BCR (10x VDJ, AIRR, dandelion, BD Rhapsody) with gene expression in an AnnData/MuData object using scirpy - chain-pairing QC, clonotype definition, clonal expansion, diversity, repertoire overlap, V(D)J usage, and VDJdb specificity. Operates on the awkward-array AIRR model (adata.obsm\['airr'\], accessed via get.airr after pp.index_chains), not legacy per-chain obs columns. Use when deciding clonotype definition for TCR (exact CDR3-nt identity via define_clonotypes) versus BCR (nucleotide distance clustering via define_clonotype_clusters with normalized_hamming plus same_v_gene/same_j_gene, because somatic hypermutation shatters identity clonotypes); tuning receptor_arms (all vs any), dual_ir, and within_group; filtering chain_qc categories (multichain doublets, orphan dropout, extra-VJ dual-TCR) without biasing clonal-expansion and diversity estimates; and overlaying clonality onto the transcriptomic UMAP.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/tcr-bcr-analysis/scirpy-analysis)
- **License:** `MIT`
- **Runs with:** Python modules `dnakit`, `scirpy`
- **Also known as:** `bio-tcr-bcr-analysis-scirpy-analysis`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info tcr-bcr-analysis-scirpy-analysis
dnakit skills install tcr-bcr-analysis-scirpy-analysis --agent claude-code   # or --agent codex / cursor
```

## `tcr-bcr-analysis-specificity-annotation` {#tcr-bcr-analysis-specificity-annotation}

Maps TCR/BCR receptor sequences toward candidate antigen specificity and clusters repertoires by shared-specificity signal, while enforcing that a database match or a cluster label is a HYPOTHESIS, not a specificity call. Use when deciding among database annotation (VDJdb/McPAS/IEDB+TCRMatch, requiring V-gene and HLA concordance plus a confidence score) versus sequence clustering (tcrdist3 meta-clonotypes, GLIPH2, GIANA, clusTCR, which find enrichment not per-receptor labels) versus generation-probability nulls (OLGA Pgen, IGoR, SONIA Ppost) for testing public/convergent/shared claims; and when guarding against overclaiming specificity, base-rate false positives from bare CDR3 matches, unpaired beta-only annotation, ML predictor failure on unseen epitopes, and ignored MHC restriction. TCR-focused with a BCR/antibody note (SHM, conformational epitopes, IGHV3-53/3-66 public clonotypes). Keywords CDR3, pMHC, HLA restriction, cross-reactivity, meta-clonotype, Pgen, public clonotype, convergent recombination.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/tcr-bcr-analysis/specificity-annotation)
- **License:** `MIT`
- **Also known as:** `bio-tcr-bcr-analysis-specificity-annotation`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info tcr-bcr-analysis-specificity-annotation
dnakit skills install tcr-bcr-analysis-specificity-annotation --agent claude-code   # or --agent codex / cursor
```

## `tcr-bcr-analysis-vdjtools-analysis` {#tcr-bcr-analysis-vdjtools-analysis}

Computes immune-repertoire diversity, clonal structure, overlap, and segment usage from TCR/BCR clonotype tables with VDJtools (immunarch as the modern R alternative). Use when deciding which diversity estimator answers a question (q=0 observed richness/chao1/chaoE, q=1 shannonWienerIndex, q=2 inverseSimpson as a Hill profile); normalizing sequencing depth before any cross-sample claim (DownSample or the resampled CalcDiversityStats table); choosing an overlap metric (depth-robust MorisitaHorn/F2 vs depth-biased Jaccard/public counts) and a clonotype match key (-i nt/aa, +/-V/J); summarizing clonality as 1 - normalizedShannonWienerIndex; reading spectratype and V-J usage under primer bias; interpreting public clonotypes; and choosing VDJtools (stable Java CLI) vs immunarch (active tidy R).

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/tcr-bcr-analysis/vdjtools-analysis)
- **License:** `MIT`
- **Runs with:** commands `VDJtools`
- **Also known as:** `bio-tcr-bcr-analysis-vdjtools-analysis`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info tcr-bcr-analysis-vdjtools-analysis
dnakit skills install tcr-bcr-analysis-vdjtools-analysis --agent claude-code   # or --agent codex / cursor
```

## `tumor-fraction-estimation` {#tumor-fraction-estimation}

Estimates tumor fraction (the genome-wide proportion of cfDNA molecules that are tumor-derived, the cfDNA analogue of bulk-tumor purity) from shallow whole-genome sequencing with ichorCNA, an HMM over 1 Mb bins that jointly EM-estimates tumor fraction, ploidy, and subclonal prevalence over a normal/ploidy grid. Encodes the load-bearing reframes: tumor fraction is the quantity that travels across assays and is NOT mutation VAF (clonal-het VAF approximately TF/2), CNA-based estimation has a hard ~3 percent limit-of-detection floor, and near-diploid or copy-neutral-LOH genomes return a falsely low value. Selects the estimator by data type (sWGS to ichorCNA, deep panel to max-VAF, methylation to deconvolution, sub-3 percent to fragmentomics or methylation). Use when quantifying tumor burden from a liquid biopsy, picking a tumor-fraction estimator for a given assay, or reconciling a TF estimate against a panel VAF.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/liquid-biopsy/tumor-fraction-estimation)
- **License:** `MIT`
- **Also known as:** `bio-tumor-fraction-estimation`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info tumor-fraction-estimation
dnakit skills install tumor-fraction-estimation --agent claude-code   # or --agent codex / cursor
```

## `variant-interpretation` {#variant-interpretation}

Clinical variant interpretation from raw variant calls to ACMG-classified recommendations with structural impact analysis. Use for VUS classification, pathogenicity assessment with cited criteria, structure-based variant impact (AlphaFold/PDB), non-coding/regulatory variant effect prediction with sequence deep-learning models (AlphaGenome, Enformer, Borzoi, ChromBPNet, Evo 2), and producing clinical-grade variant reports for return of results or molecular tumor boards. Use this whenever a user asks about a variant's significance, an intronic/promoter/enhancer/UTR non-coding variant's functional impact, or needs ACMG classification — even if they don't say "ACMG".

- **Source:** [ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-variant-interpretation)
- **License:** `Apache-2.0`
- **Runs with:** Python modules `cellxgene_census`, `dnakit`
- **Also known as:** `tooluniverse-variant-interpretation`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info variant-interpretation
dnakit skills install variant-interpretation --agent claude-code   # or --agent codex / cursor
```
