# Population and evolutionary genetics

**21** Skills. A Skill is a set of instructions that guides an Agent through an analysis; the software, reference data and APIs it uses are separate dependencies. See [Skills](../skills.md) for searching, installing and using them; the Skills tab of the web demo shows each Skill in full.

| Skill | Source | Summary |
| --- | --- | --- |
| [`ancestry-risk-profiler`](#ancestry-risk-profiler) | clawbio | Infers genetic super-population ancestry from a 23andMe/AncestryDNA file and computes ancestry-stratified odds ratios with an exploratory Ancestry Elevation Sc… |
| [`equity-scorer`](#equity-scorer) | clawbio | Compute HEIM diversity and equity metrics from VCF or ancestry data. Generates heterozygosity, FST, PCA plots, and a composite HEIM Equity Score with markdown… |
| [`genotype-imputation-benchmarking`](#genotype-imputation-benchmarking) | research-suite | Design, run, diagnose, and audit genotype-imputation benchmarks across low coverage, arrays, panels, families, or reference-free workflows. Use when a task nee… |
| [`gigwa-diversity-report`](#gigwa-diversity-report) | gigwa-mcp | Produce a population diversity and structure report for a Gigwa variant set. Use this when a user asks to summarise genetic diversity, describe population stru… |
| [`gigwa-region-scan`](#gigwa-region-scan) | gigwa-mcp | Characterise the variants and diversity within a single genomic region of a Gigwa variant set. Use this when a user asks about a specific chromosome, contig or… |
| [`grapevine-adna`](#grapevine-adna) | plant-skills | Ancient DNA processing and authentication for plants (grapevine-oriented). Use for aDNA FASTQ→BAM, damage, contamination, capture, or projecting ancient sample… |
| [`ld-1000g-region-compute`](#ld-1000g-region-compute) | clawbio | Compute pairwise r² between a lead variant and every variant in a window using the 1000 Genomes Phase 3 GRCh38 reference panel, ancestry-stratified. Use when a… |
| [`onekgpd-individual-queries`](#onekgpd-individual-queries) | user335-kdense | Query the 1000 Genomes Project dataset (3,202 whole-genome-sequenced individuals, GRCh38) at the level of individual participants. Use when a question is about… |
| [`phasing-imputation-genotype-imputation`](#phasing-imputation-genotype-imputation) | bioskills | Imputes untyped genotypes against a phased reference panel with Beagle, Minimac4, or IMPUTE5 (array data) or from genotype likelihoods with GLIMPSE2, QUILT2, o… |
| [`phasing-imputation-haplotype-phasing`](#phasing-imputation-haplotype-phasing) | bioskills | Estimates haplotype phase from population linkage disequilibrium with SHAPEIT5, SHAPEIT4, Eagle2, or Beagle - turning unphased genotypes (0/1) into phased hapl… |
| [`phasing-imputation-imputation-qc`](#phasing-imputation-imputation-qc) | bioskills | Assesses and filters phasing/imputation output - the quality metrics (Beagle DR2, Minimac R2 and EmpRsq, IMPUTE/GLIMPSE INFO), MAF-stratified filtering, true a… |
| [`phasing-imputation-reference-panels`](#phasing-imputation-reference-panels) | bioskills | Selects and prepares the reference panel that phasing/imputation copies haplotypes from (1000 Genomes, HRC, TOPMed, HGDP+1kGP/gnomAD, CAAPA), matching panel an… |
| [`population-genetics-association-testing`](#population-genetics-association-testing) | bioskills | Single-variant common-variant GWAS with plink2 --glm (linear/logistic, Firth) and the linear mixed models GEMMA, BOLT-LMM, SAIGE, regenie (SPA). A GWAS statist… |
| [`population-genetics-linkage-disequilibrium`](#population-genetics-linkage-disequilibrium) | bioskills | Computes linkage disequilibrium (r2, D', composite Rogers-Huff r2), prunes correlated variants, clumps GWAS summary statistics to lead SNPs, and defines haplot… |
| [`population-genetics-plink-basics`](#population-genetics-plink-basics) | bioskills | Manages PLINK genotype filesets - format conversion (VCF, BED/BIM/FAM, PED/MAP, pgen/pvar/psam) and sample/variant QC (missingness, MAF, HWE, sex check, hetero… |
| [`population-genetics-population-structure`](#population-genetics-population-structure) | bioskills | Infers and describes population structure with PCA (plink2 --pca, smartpca/EIGENSOFT, FlashPCA2), model-based clustering (ADMIXTURE, fastSTRUCTURE), FST estima… |
| [`population-genetics-rare-variant-association`](#population-genetics-rare-variant-association) | bioskills | Gene and region-based rare-variant aggregation - burden/collapsing, SKAT, SKAT-O, ACAT-V/ACAT-O, annotation-weighted STAAR - with regenie (--vc-tests), SAIGE-G… |
| [`population-genetics-scikit-allel-analysis`](#population-genetics-scikit-allel-analysis) | bioskills | In-memory Python population genetics with scikit-allel - GenotypeArray/HaplotypeArray/AlleleCountsArray, diversity (pi, theta, Tajima's D), SFS, FST (Weir-Cock… |
| [`population-genetics-selection-statistics`](#population-genetics-selection-statistics) | bioskills | Scans genomes for natural selection with SFS tests (Tajima's D, Fay &amp; Wu H, Zeng E, SweepFinder2 CLR), haplotype tests (iHS, nSL, XP-EHH, Rsb, H12), and differ… |
| [`recombinator`](#recombinator) | clawbio | Produce offspring genomes from parent pairs via meiotic recombination, mutation, and clinical evaluation |
| [`vitis-popgen`](#vitis-popgen) | plant-skills | Population-genetics workflow from a clean VCF: clones, PCA, ADMIXTURE, Fst, TreeMix, Dsuite, selection. Use for Vitis / perennial crops where clones and famili… |

## `ancestry-risk-profiler` {#ancestry-risk-profiler}

Infers genetic super-population ancestry from a 23andMe/AncestryDNA file and computes ancestry-stratified odds ratios with an exploratory Ancestry Elevation Score (AES) showing where ancestry-specific GWAS effect sizes diverge from European reference estimates.

- **Source:** [ClawBio; ClawBio contributors @ `1fcecb7`](https://github.com/ClawBio/ClawBio/tree/1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d/skills/ancestry-risk-profiler)
- **License:** `MIT`
- **Runs with:** commands `python3`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info ancestry-risk-profiler
dnakit skills install ancestry-risk-profiler --agent claude-code   # or --agent codex / cursor
```

## `equity-scorer` {#equity-scorer}

Compute HEIM diversity and equity metrics from VCF or ancestry data. Generates heterozygosity, FST, PCA plots, and a composite HEIM Equity Score with markdown reports.

- **Source:** [ClawBio contributors @ `1fcecb7`](https://github.com/ClawBio/ClawBio/tree/1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d/skills/equity-scorer)
- **License:** `MIT`
- **Runs with:** commands `python3`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info equity-scorer
dnakit skills install equity-scorer --agent claude-code   # or --agent codex / cursor
```

## `genotype-imputation-benchmarking` {#genotype-imputation-benchmarking}

Design, run, diagnose, and audit genotype-imputation benchmarks across low coverage, arrays, panels, families, or reference-free workflows. Use when a task needs data contracts, leakage control, stratification, evaluability, accuracy, calibration, downstream GWAS, or complementary-error testing. Do not use for generic genotype QC or one-tool installation.

- **Source:** [Ynzuses/research-skills-suite contributors @ `a5cb708`](https://github.com/Ynzuses/research-skills-suite/tree/a5cb708252a58bcfb771d46d39f5e38d64b8b2f6/skills/genotype-imputation-benchmarking)
- **License:** `Apache-2.0`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genotype-imputation-benchmarking
dnakit skills install genotype-imputation-benchmarking --agent claude-code   # or --agent codex / cursor
```

## `gigwa-diversity-report` {#gigwa-diversity-report}

Produce a population diversity and structure report for a Gigwa variant set. Use this when a user asks to summarise genetic diversity, describe population structure, or assess relatedness for a run identified by its variantSetDbId. Walks through dataset-wide MAF/He/Ho/PIC/Fis, principal component analysis (PCA), population-structure clustering with a suggested number of clusters K, a VanRaden kinship/relatedness matrix highlighting the closest pairs, and a UPGMA dendrogram; when a metadata TSV and group column are supplied it adds per-group diversity and pairwise Weir and Cockerham Fst. Trigger on 'run a diversity analysis in Gigwa', 'show population structure / PCA / kinship / phylogenetic tree', or 'compare diversity between populations'. Keywords: Gigwa, population genetics, diversity, PCA, population structure, kinship, Fst, UPGMA tree, He, Ho, PIC, MAF.

- **Source:** [Copyright 2026 Anestis Gkanogiannis; Gigwa-MCP contributors @ `d6ccdd1`](https://github.com/gkanogiannis/Gigwa-MCP/tree/d6ccdd1a5eb19b8a0cc224b96ad459868d9c4aaf/skills/gigwa-diversity-report)
- **License:** `Apache-2.0`
- **Runs with:** Python modules `gigwa_mcp`; commands `gigwa-mcp`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info gigwa-diversity-report
dnakit skills install gigwa-diversity-report --agent claude-code   # or --agent codex / cursor
```

## `gigwa-region-scan` {#gigwa-region-scan}

Characterise the variants and diversity within a single genomic region of a Gigwa variant set. Use this when a user asks about a specific chromosome, contig or coordinate window (e.g. 'chr1' or 'chr1:1000000-2000000') in a run identified by its variantSetDbId. Confirms the contig name, counts variants in the window and again at a common-variant MAF threshold, lists the matching variants to a CSV, and characterises just that window with a diversity summary and call-rate QC. Trigger on 'scan this region in Gigwa', 'how many variants are in chr1:...', 'is this locus variable', or 'diversity around a gene or window'. Note the interval interface differs by tool: count_variants and search_variants take reference_name/start/end, while diversity_summary and qc_call_rate take a single region string. Keywords: Gigwa, genomic region, variant density, region scan, chromosome, contig, MAF, variant search, call rate.

- **Source:** [Copyright 2026 Anestis Gkanogiannis; Gigwa-MCP contributors @ `d6ccdd1`](https://github.com/gkanogiannis/Gigwa-MCP/tree/d6ccdd1a5eb19b8a0cc224b96ad459868d9c4aaf/skills/gigwa-region-scan)
- **License:** `Apache-2.0`
- **Runs with:** Python modules `gigwa_mcp`; commands `gigwa-mcp`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info gigwa-region-scan
dnakit skills install gigwa-region-scan --agent claude-code   # or --agent codex / cursor
```

## `grapevine-adna` {#grapevine-adna}

Ancient DNA processing and authentication for plants (grapevine-oriented). Use for aDNA FASTQ→BAM, damage, contamination, capture, or projecting ancient samples onto a modern panel. Not for modern GATK WGS.

- **Source:** [Xuzhen-Li/bioinfo-agent-skills contributors @ `6dea4be`](https://github.com/Xuzhen-Li/bioinfo-agent-skills/tree/6dea4bee7aade611d22f2f1958612328bc61ec83/skills/grapevine-adna)
- **License:** `MIT`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info grapevine-adna
dnakit skills install grapevine-adna --agent claude-code   # or --agent codex / cursor
```

## `ld-1000g-region-compute` {#ld-1000g-region-compute}

Compute pairwise r² between a lead variant and every variant in a window
using the 1000 Genomes Phase 3 GRCh38 reference panel, ancestry-stratified.
Use when an agent needs LD coloring for a regional plot or LD pruning
around a candidate causal variant. Single client (on-demand region fetch
from EBI 1000G FTP); no multi-GB cold-start.

- **Source:** [Aviv Madar; ClawBio contributors @ `1fcecb7`](https://github.com/ClawBio/ClawBio/tree/1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d/skills/ld-1000g-region-compute)
- **License:** `MIT`
- **Runs with:** commands `plink`, `python3`, `tabix`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info ld-1000g-region-compute
dnakit skills install ld-1000g-region-compute --agent claude-code   # or --agent codex / cursor
```

## `onekgpd-individual-queries` {#onekgpd-individual-queries}

Query the 1000 Genomes Project dataset (3,202 whole-genome-sequenced individuals, GRCh38) at the level of individual participants. Use when a question is about individuals or variants in the 1000 Genomes Project cohort: which individuals carry variants matching specific criteria in a gene or region, which individuals are homozygous-reference at a position, which variants exist in the dataset or carried by specified individuals in a gene or region, the relatedness between two specified individuals. Variants are returned with 1000 Genomes allele frequencies (AF), gnomAD v4.1 exome and genome AF, AlphaMissense score, and HGVSp annotations.

- **Source:** [Dnaerys (OneKGPd skill author); K-Dense Inc. repository license @ `330c8e7`](https://github.com/K-Dense-AI/scientific-agent-skills/tree/330c8e764435a731eff571e3efdda70b363d0792/skills/onekgpd)
- **License:** `MIT`
- **Runs with:** Python modules `dnaerys`; commands `python3`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info onekgpd-individual-queries
dnakit skills install onekgpd-individual-queries --agent claude-code   # or --agent codex / cursor
```

## `phasing-imputation-genotype-imputation` {#phasing-imputation-genotype-imputation}

Imputes untyped genotypes against a phased reference panel with Beagle, Minimac4, or IMPUTE5 (array data) or from genotype likelihoods with GLIMPSE2, QUILT2, or STITCH (low-coverage WGS), producing per-variant dosages (DS) with a self-estimated quality (Beagle DR2, Minimac R2, IMPUTE INFO). Covers why the honest output is a dosage posterior not a hard call, why GWAS regresses on DS, why the quality metric is an ESTIMATE of r2 from posterior spread (not validation against truth), the DS/GP/HDS fields, the phasing prerequisite, chunking, chrX ploidy, the Michigan/TOPMed servers (the only access to HRC/TOPMed), and low-coverage WGS as the modern array replacement. Use when increasing variant density for GWAS, harmonizing arrays, inferring untyped variants, or imputing low-coverage sequence. Phase first with haplotype-phasing; prepare the panel with reference-panels; filter with imputation-qc; the GWAS test is population-genetics/association-testing; end-to-end orchestration is workflows/gwas-pipeline.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/phasing-imputation/genotype-imputation)
- **License:** `MIT`
- **Runs with:** commands `Beagle`
- **Also known as:** `bio-phasing-imputation-genotype-imputation`, `phasing-imputation`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info phasing-imputation-genotype-imputation
dnakit skills install phasing-imputation-genotype-imputation --agent claude-code   # or --agent codex / cursor
```

## `phasing-imputation-haplotype-phasing` {#phasing-imputation-haplotype-phasing}

Estimates haplotype phase from population linkage disequilibrium with SHAPEIT5, SHAPEIT4, Eagle2, or Beagle - turning unphased genotypes (0/1) into phased haplotypes (0\|1) for imputation input, compound-heterozygote calls, HLA typing, or population genetics. Covers why statistical phase is an INFERENCE (not a measurement) whose error concentrates at rare variants, why a genome-wide switch-error rate hides catastrophic rare-variant error and must be reported MAC-stratified, the SHAPEIT5 common-scaffold-then-rare design (phase_common, ligate, phase_rare, switch), reference-based vs within-cohort phasing, the build-matched genetic map, chrX male-haploid handling, and the switch-vs-flip-vs-Hamming distinction. Use when phasing genotypes before imputation, for compound-het/ASE/HLA, or benchmarking against trios. Read-backed / molecular phasing (long reads, Hi-C) is long-read-sequencing/haplotype-phasing; panel choice is reference-panels; imputation is genotype-imputation.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/phasing-imputation/haplotype-phasing)
- **License:** `MIT`
- **Runs with:** commands `SHAPEIT5`
- **Also known as:** `bio-phasing-imputation-haplotype-phasing`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info phasing-imputation-haplotype-phasing
dnakit skills install phasing-imputation-haplotype-phasing --agent claude-code   # or --agent codex / cursor
```

## `phasing-imputation-imputation-qc` {#phasing-imputation-imputation-qc}

Assesses and filters phasing/imputation output - the quality metrics (Beagle DR2, Minimac R2 and EmpRsq, IMPUTE/GLIMPSE INFO), MAF-stratified filtering, true accuracy by masking, the differential-imputation confound, dosage-based downstream usage, and phasing switch-error QC. Covers why every routine quality score is an ESTIMATE of r2 from the posterior spread (not validation against truth), why it is confounded with MAF so a flat INFO&gt;=0.3 cutoff is a hidden rare-variant filter, why concordance lies for rare variants while masked dosage-r2 by MAF is the gold standard, why separate case/control imputation manufactures false GWAS hits, and that the field name tells the tool (DR2=Beagle, R2=Minimac, INFO=GLIMPSE/IMPUTE). Use when filtering imputed variants before GWAS, validating accuracy, benchmarking phasing against trios, or diagnosing inflated association. Imputation is genotype-imputation; phasing is haplotype-phasing; panel ancestry is reference-panels; the test is population-genetics/association-testing.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/phasing-imputation/imputation-qc)
- **License:** `MIT`
- **Also known as:** `bio-phasing-imputation-imputation-qc`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info phasing-imputation-imputation-qc
dnakit skills install phasing-imputation-imputation-qc --agent claude-code   # or --agent codex / cursor
```

## `phasing-imputation-reference-panels` {#phasing-imputation-reference-panels}

Selects and prepares the reference panel that phasing/imputation copies haplotypes from (1000 Genomes, HRC, TOPMed, HGDP+1kGP/gnomAD, CAAPA), matching panel ancestry to the target, reconciling genome build and chromosome naming, and running the strand/allele harmonization gate. Covers why ancestry-match beats panel size (imputation can only copy haplotypes the panel contains), why palindromic A/T and C/G SNPs flip strand without erroring, why liftover is a strand-flip generator in between-build inverted regions, that HRC is SNP-only and TOPMed is never downloadable (governance can override accuracy), and panel formats (msav, bref3, imp5). Use when choosing a panel for a target ancestry, preparing or converting a panel, aligning study data, or deciding between downloadable and server-only panels. Phasing is haplotype-phasing; imputation is genotype-imputation; PCA for ancestry is population-genetics/population-structure; HLA panels are clinical-databases/hla-typing.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/phasing-imputation/reference-panels)
- **License:** `MIT`
- **Runs with:** commands `bcftools`
- **Also known as:** `bio-phasing-imputation-reference-panels`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info phasing-imputation-reference-panels
dnakit skills install phasing-imputation-reference-panels --agent claude-code   # or --agent codex / cursor
```

## `population-genetics-association-testing` {#population-genetics-association-testing}

Single-variant common-variant GWAS with plink2 --glm (linear/logistic, Firth) and the linear mixed models GEMMA, BOLT-LMM, SAIGE, regenie (SPA). A GWAS statistic is valid only when genotype is independent of unmodeled phenotype drivers after the chosen covariates and random effects, so the engine follows sample structure and case:control imbalance, not taste: PC covariates absorb continuous ancestry but cannot remove relatedness (a covariance structure needing an LMM), genomic inflation above 1 is mostly true polygenic signal not confounding (the LDSC intercept is the diagnostic), LOCO prevents proximal contamination, and SPA/Firth keep the tail calibrated at extreme imbalance and low MAC. Use when running single-variant GWAS, choosing between a GLM and a mixed model, or controlling stratification, relatedness, and case:control imbalance. For rare-variant aggregation (burden, SKAT, SKAT-O, ACAT) see rare-variant-association; fine-mapping and MR see causal-genomics; PRS see clinical-databases/polygenic-risk.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/population-genetics/association-testing)
- **License:** `MIT`
- **Runs with:** commands `plink2`
- **Also known as:** `bio-population-genetics-association-testing`, `plink2-gwas-analysis`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info population-genetics-association-testing
dnakit skills install population-genetics-association-testing --agent claude-code   # or --agent codex / cursor
```

## `population-genetics-linkage-disequilibrium` {#population-genetics-linkage-disequilibrium}

Computes linkage disequilibrium (r2, D', composite Rogers-Huff r2), prunes correlated variants, clumps GWAS summary statistics to lead SNPs, and defines haplotype blocks with PLINK 1.9/2.0 and scikit-allel. r2 and D' answer different questions - r2 (= chi2/N) is the tagging and GWAS-power currency, D' marks observed recombination and is upward-biased for rare variants. PLINK 2.0 has no bare --r2 (split into --r2-phased and --r2-unphased); pruning (--indep-pairwise, genotype-blind) and clumping (--clump, p-value-aware) are distinct operations that are constantly confused. The clumping or fine-mapping LD reference must be ancestry-matched or it fails silently into false credible sets. Use when calculating LD, pruning variants for PCA or structure, clumping GWAS hits, or selecting tag SNPs. For QC see plink-basics; for PCA see population-structure; fine-mapping is causal-genomics/fine-mapping.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/population-genetics/linkage-disequilibrium)
- **License:** `MIT`
- **Also known as:** `bio-population-genetics-linkage-disequilibrium`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info population-genetics-linkage-disequilibrium
dnakit skills install population-genetics-linkage-disequilibrium --agent claude-code   # or --agent codex / cursor
```

## `population-genetics-plink-basics` {#population-genetics-plink-basics}

Manages PLINK genotype filesets - format conversion (VCF, BED/BIM/FAM, PED/MAP, pgen/pvar/psam) and sample/variant QC (missingness, MAF, HWE, sex check, heterozygosity, KING relatedness) with PLINK 1.9 and 2.0. PLINK rewrites allele bookkeeping: PLINK 1.x A1 defaults to the minor allele and is recomputed every load, silently flipping effect-allele meaning unless --keep-allele-order, while PLINK 2.0 tracks explicit REF/ALT. QC order matters (variant before sample missingness), HWE is controls-only in 1.9 but not 2.0, add midp, and differential case/control missingness injects false hits. Use when converting between PLINK formats or running genotype QC before association, structure, or LD analysis. For LD pruning/clumping see linkage-disequilibrium; for GWAS see association-testing; VCF input from variant-calling/vcf-basics.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/population-genetics/plink-basics)
- **License:** `MIT`
- **Runs with:** commands `plink`
- **Also known as:** `bio-population-genetics-plink-basics`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info population-genetics-plink-basics
dnakit skills install population-genetics-plink-basics --agent claude-code   # or --agent codex / cursor
```

## `population-genetics-population-structure` {#population-genetics-population-structure}

Infers and describes population structure with PCA (plink2 --pca, smartpca/EIGENSOFT, FlashPCA2), model-based clustering (ADMIXTURE, fastSTRUCTURE), FST estimators (Weir-Cockerham vs Hudson), and f-statistics (f3/f4/D via AdmixTools/admixr), plus Python plotting of PCs and Q barplots. Every output is a model-conditioned description of variance, not truth: PCs conflate ancestry with LD/inversions/relatedness/batch, ADMIXTURE Q-values are panel- and K-dependent artifacts, and CV-minimum K is a guide not the true population count. FST must combine SNPs as a ratio of averages (sum numerators / sum denominators), never an average of per-SNP FST; negative per-SNP FST is normal and must not be clamped. f3/f4/D need a block jackknife or the significance is fake. Use when running PCA, ADMIXTURE, FST, or f-statistics on QC'd genotypes. For QC and KING relatedness see plink-basics; for LD pruning see linkage-disequilibrium; for array-based Python pipelines see scikit-allel-analysis.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/population-genetics/population-structure)
- **License:** `MIT`
- **Also known as:** `bio-population-genetics-population-structure`, `claw-ancestry-pca`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info population-genetics-population-structure
dnakit skills install population-genetics-population-structure --agent claude-code   # or --agent codex / cursor
```

## `population-genetics-rare-variant-association` {#population-genetics-rare-variant-association}

Gene and region-based rare-variant aggregation - burden/collapsing, SKAT, SKAT-O, ACAT-V/ACAT-O, annotation-weighted STAAR - with regenie (--vc-tests), SAIGE-GENE+, and the SKAT R package. Single-variant tests are powerless at low minor allele count, so rare variants are aggregated across a gene or region under an explicit mask (functional class plus a MAF cutoff). A burden test collapses variants into one score assuming a single effect direction (powerful when true, near-zero power when risk and protective variants cancel); SKAT is a variance-component test robust to mixed directions; SKAT-O blends the two; ACAT/STAAR are dependence-robust and annotation-weighted. The mask is the hypothesis, imbalance needs SPA or Firth, and testing burden is per-gene-per-mask. Use when aggregating rare coding or regulatory variants into gene or region tests, choosing burden vs SKAT vs SKAT-O, or building masks. For single-variant GWAS see association-testing; for mask annotations see variant-calling/variant-annotation.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/population-genetics/rare-variant-association)
- **License:** `MIT`
- **Also known as:** `bio-population-genetics-rare-variant-association`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info population-genetics-rare-variant-association
dnakit skills install population-genetics-rare-variant-association --agent claude-code   # or --agent codex / cursor
```

## `population-genetics-scikit-allel-analysis` {#population-genetics-scikit-allel-analysis}

In-memory Python population genetics with scikit-allel - GenotypeArray/HaplotypeArray/AlleleCountsArray, diversity (pi, theta, Tajima's D), SFS, FST (Weir-Cockerham, Hudson, Patterson), f3/D admixture stats, LD pruning, PCA, and selection scans (iHS, XP-EHH, nSL, Garud H). Nearly every statistic is a ratio or density with one silent denominator bug in two faces: omit is_accessible= and per-base pi/theta divide by total span not accessible bp (deflated 2-5x); average per-SNP FST instead of sum(a)/(sum(a)+sum(b)+sum(c)) and the estimate is rare-variant-biased - scikit-allel returns the (a,b,c) and (num,den) components on purpose to force ratio-of-sums. to_n_alt default fill=0 imputes missing to reference; sfs() is unfolded and wants derived not alt counts; iHS/XP-EHH need phased data and standardization. Use when computing population-genetics statistics in Python, scanning for selection, or building array pipelines. For PLINK QC see plink-basics; for VCF input see variant-calling/vcf-basics.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/population-genetics/scikit-allel-analysis)
- **License:** `MIT`
- **Also known as:** `bio-population-genetics-scikit-allel-analysis`, `dnasp`, `tooluniverse-population-genetics`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info population-genetics-scikit-allel-analysis
dnakit skills install population-genetics-scikit-allel-analysis --agent claude-code   # or --agent codex / cursor
```

## `population-genetics-selection-statistics` {#population-genetics-selection-statistics}

Scans genomes for natural selection with SFS tests (Tajima's D, Fay &amp; Wu H, Zeng E, SweepFinder2 CLR), haplotype tests (iHS, nSL, XP-EHH, Rsb, H12), and differentiation (FST, PBS) using scikit-allel, selscan, and SweepFinder2. No single statistic separates selection from demography at one locus, so the deliverable is empirical genome-wide outliers plus multiple orthogonal signals, not an absolute cutoff. iHS detects incomplete sweeps and collapses to zero at fixation while XP-EHH catches fixed sweeps; iHS/nSL standardize within derived-allele-frequency bins but XP-EHH gets a genome-wide z-score; derived-allele tests need substitution-model polarization; background selection mimics FST and CLR. Use when computing selection statistics like FST, Tajima's D, iHS, or XP-EHH, or scanning for selective sweeps. For phasing inputs see phasing-imputation/haplotype-phasing; for dN/dS see comparative-genomics/positive-selection.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/population-genetics/selection-statistics)
- **License:** `MIT`
- **Also known as:** `bio-population-genetics-selection-statistics`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info population-genetics-selection-statistics
dnakit skills install population-genetics-selection-statistics --agent claude-code   # or --agent codex / cursor
```

## `recombinator` {#recombinator}

Produce offspring genomes from parent pairs via meiotic recombination, mutation, and clinical evaluation

- **Source:** [Manuel Corpas; ClawBio contributors @ `1fcecb7`](https://github.com/ClawBio/ClawBio/tree/1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d/skills/recombinator)
- **License:** `MIT`
- **Runs with:** commands `python3`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info recombinator
dnakit skills install recombinator --agent claude-code   # or --agent codex / cursor
```

## `vitis-popgen` {#vitis-popgen}

Population-genetics workflow from a clean VCF: clones, PCA, ADMIXTURE, Fst, TreeMix, Dsuite, selection. Use for Vitis / perennial crops where clones and families wreck structure plots.

- **Source:** [Xuzhen-Li/bioinfo-agent-skills contributors @ `6dea4be`](https://github.com/Xuzhen-Li/bioinfo-agent-skills/tree/6dea4bee7aade611d22f2f1958612328bc61ec83/skills/vitis-popgen)
- **License:** `MIT`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info vitis-popgen
dnakit skills install vitis-popgen --agent claude-code   # or --agent codex / cursor
```
