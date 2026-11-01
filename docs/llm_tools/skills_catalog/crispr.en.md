# CRISPR and molecular design

**37** Skills. A Skill is a set of instructions that guides an Agent through an analysis; the software, reference data and APIs it uses are separate dependencies. See [Skills](../skills.md) for searching, installing and using them; the Skills tab of the web demo shows each Skill in full.

| Skill | Source | Summary |
| --- | --- | --- |
| [`crispr-screen-triage`](#crispr-screen-triage) | clawbio | Deterministic CRISPR screen hit ranking from local guide-level count tables |
| [`crispr-screens-bagel-essentiality`](#crispr-screens-bagel-essentiality) | bioskills | Identifies essential genes from CRISPR-Cas9 fitness screens using BAGEL2 (Kim &amp; Hart 2021 Genome Med), a Bayesian classifier scoring per-gene Bayes Factors via… |
| [`crispr-screens-base-editing-analysis`](#crispr-screens-base-editing-analysis) | bioskills | Analyzes base-editing screens for variant function. Covers library design (Hanna 2021 ClinVar-scale CBE screen benchmarked on BRCA1/2, Cuella-Martin 2021 DDR s… |
| [`crispr-screens-batch-correction`](#crispr-screens-batch-correction) | bioskills | Batch effect correction for CRISPR screens covering ComBat empirical-Bayes, RUV, SVA, control-sgRNA normalization, and the model-based alternative of including… |
| [`crispr-screens-combinatorial-screens`](#crispr-screens-combinatorial-screens) | bioskills | Designs and analyzes combinatorial CRISPR screens covering paired-Cas9 (Big Papi, Najm 2018), enhanced AsCas12a multiplex (enCas12a, DeWeirdt 2021), in4mer 4-g… |
| [`crispr-screens-copy-number-correction`](#crispr-screens-copy-number-correction) | bioskills | Corrects the gene-independent copy-number artifact in CRISPR-Cas9 screens (Aguirre 2016 / Munoz 2016 Cancer Discov) where amplified loci appear essential from… |
| [`crispr-screens-crispresso-editing`](#crispr-screens-crispresso-editing) | bioskills | Quantifies CRISPR editing outcomes with CRISPResso2 (Clement 2019 Nat Biotechnol) across Cas9-nuclease (indels, HDR), CBE and ABE base editors (target conversi… |
| [`crispr-screens-drugz-chemogenomic`](#crispr-screens-drugz-chemogenomic) | bioskills | Analyzes CRISPR drug-modifier (chemogenomic) screens with drugZ (Colic et al. 2019 Genome Med), a bidirectional Z-score method that identifies synthetic-lethal… |
| [`crispr-screens-hit-calling`](#crispr-screens-hit-calling) | bioskills | Cross-method decision tree for calling hits in pooled CRISPR screens. Catalogs statistical models (MAGeCK RRA, MAGeCK MLE, BAGEL2, drugZ, JACKS, Chronos, CERES… |
| [`crispr-screens-in-vivo-screens`](#crispr-screens-in-vivo-screens) | bioskills | Designs and analyzes in vivo CRISPR screens in animal tumor models, organoids, and immune-cell adoptive transfers. Covers bottleneck math (250x cells/sgRNA req… |
| [`crispr-screens-jacks-analysis`](#crispr-screens-jacks-analysis) | bioskills | Runs JACKS (Joint Analysis of CRISPR/Cas9 Knockout Screens; Allen et al 2019 Genome Research) which models per-sgRNA log-fold-change as the product of a treatm… |
| [`crispr-screens-library-design`](#crispr-screens-library-design) | bioskills | Designs pooled sgRNA libraries for CRISPR knockout, interference (CRISPRi), activation (CRISPRa), Cas12a multiplex, base-editor, and prime-editor screens. Cove… |
| [`crispr-screens-mageck-analysis`](#crispr-screens-mageck-analysis) | bioskills | Analyzes pooled CRISPR screens with MAGeCK (Li et al 2014), covering count generation (mageck count), the RRA two-condition workflow (mageck test using alpha-R… |
| [`crispr-screens-perturb-seq-analysis`](#crispr-screens-perturb-seq-analysis) | bioskills | Analyzes single-cell pooled CRISPR screens (Perturb-seq, CROP-seq, Perturb-CITE-seq, ECCITE-seq, multiome) where each cell carries an sgRNA and a scRNA-seq / s… |
| [`crispr-screens-prime-editing-screens`](#crispr-screens-prime-editing-screens) | bioskills | Designs and analyzes pooled prime-editor (PE) screens for installing precise genetic variants without bystander confounding. Covers pegRNA design with PRIDICT… |
| [`crispr-screens-screen-qc`](#crispr-screens-screen-qc) | bioskills | Quality control for pooled CRISPR screens covering library representation, Gini index, log-skew, replicate Pearson and Spearman concordance, essentialome preci… |
| [`functional-genomics-screens`](#functional-genomics-screens) | tooluniverse | Interpret hits from CRISPR-KO/CRISPRi/shRNA screens by integrating DepMap essentiality, gnomAD constraint scores, pathway context (Reactome, STRING), druggabil… |
| [`functional-screen-analysis`](#functional-screen-analysis) | encode-toolkit | Analyze ENCODE functional genomics screens including CRISPR screens, MPRA (Massively Parallel Reporter Assays), and STARR-seq. Find screen data in ENCODE, proc… |
| [`genome-engineering-base-editing-design`](#genome-engineering-base-editing-design) | bioskills | Designs cytosine (CBE, C-to-T) and adenine (ABE, A-to-G) base-editor guides by positioning the target base at the activity-peak of the editing window (protospa… |
| [`genome-engineering-grna-design`](#genome-engineering-grna-design) | bioskills | Designs and ranks guide RNAs (sgRNAs) for CRISPR-Cas9/Cas12a gene knockout by scanning a target for PAM sites (NGG SpCas9, NNGRRT SaCas9, TTTV Cas12a, NG SpCas… |
| [`genome-engineering-hdr-template-design`](#genome-engineering-hdr-template-design) | bioskills | Designs donor/repair templates for precise CRISPR knock-ins -- choosing the format (ssODN, long-ssDNA/Easi-CRISPR, dsDNA/plasmid, AAV6), sizing homology arms,… |
| [`genome-engineering-off-target-prediction`](#genome-engineering-off-target-prediction) | bioskills | Nominates and assesses CRISPR off-target sites genome-wide. Enumerates candidate sites by mismatch and bulge tolerance with Cas-OFFinder/CRISPRitz, ranks them… |
| [`genome-engineering-prime-editing-design`](#genome-engineering-prime-editing-design) | bioskills | Designs pegRNAs and nicking guides for prime editing (PE) -- choosing the nick/strand, tuning the primer-binding site (PBS) and reverse-transcription template… |
| [`molecular-cloning`](#molecular-cloning) | tooluniverse | Molecular cloning, in both directions. DESIGN — Gibson Assembly (overlap design for seamless multi-fragment joining) and Golden Gate Assembly (Type IIS / BsaI… |
| [`plannotate-plasmid-annotation`](#plannotate-plasmid-annotation) | sciagent | Auto-annotate plasmids with features (promoters, terminators, resistance, origins, tags, fluorescent proteins) via BLAST against curated DBs (Addgene, fpbase,… |
| [`plant-genome-guided-primer-design`](#plant-genome-guided-primer-design) | plant-qpcr | A genome-guided workflow for plant gene identity validation, homolog analysis, and experimental primer design. The workflow prevents incorrect primer design ca… |
| [`polyploid-gene-family-cloning`](#polyploid-gene-family-cloning) | polyploid-cloning | AI Skill for genome-aware polyploid gene family identification, candidate prioritization, validated full-length CDS cloning, construction-intent locking, homoe… |
| [`primer-design-primer-basics`](#primer-design-primer-basics) | bioskills | Designs and ranks PCR primer pairs for a target template with primer3-py (design_primers), returning pairs with nearest-neighbor Tm, GC, product size, and comp… |
| [`primer-design-primer-specificity`](#primer-design-primer-specificity) | bioskills | Checks whether a PCR primer PAIR amplifies only the intended target genome-wide, using pair-aware in-silico PCR (MFEprimer-3.0, UCSC isPcr, NCBI Primer-BLAST)… |
| [`primer-design-primer-validation`](#primer-design-primer-validation) | bioskills | Validates chosen PCR/qPCR oligos for intramolecular thermodynamic liabilities with primer3-py - hairpins, self-dimers, cross-dimers (calc_hairpin/homodimer/het… |
| [`primer-design-qpcr-primers`](#primer-design-qpcr-primers) | bioskills | Co-designs qPCR/RT-qPCR primers and hydrolysis (TaqMan) or molecular-beacon probes with primer3-py (PRIMER_PICK_INTERNAL_OLIGO, PRIMER_INTERNAL_* tags), for as… |
| [`restriction-enzyme-selection`](#restriction-enzyme-selection) | bioskills | Select restriction enzymes for cloning or diagnostics using Biopython Bio.Restriction. Finds enzymes by cut frequency, overhang type, recognition-site length,… |
| [`restriction-fragment-analysis`](#restriction-fragment-analysis) | bioskills | Predict restriction digest fragment sizes and gel patterns using Biopython Bio.Restriction. Computes fragment lengths and sequences for single and double diges… |
| [`restriction-golden-gate-assembly`](#restriction-golden-gate-assembly) | bioskills | Design and validate Type IIS scarless DNA assembly (Golden Gate, MoClo) using Biopython Bio.Restriction. Screens parts for internal BsaI/BsmBI/BbsI/SapI sites… |
| [`restriction-mapping`](#restriction-mapping) | bioskills | Build restriction maps showing enzyme cut positions and inter-site distances along DNA using Biopython Bio.Restriction. Produces text or graphical maps for lin… |
| [`restriction-sites`](#restriction-sites) | bioskills | Find restriction enzyme cut sites in DNA sequences using Biopython Bio.Restriction. Searches single enzymes, batches, or commercial enzyme sets and returns cut… |
| [`single-cell-lineage-tracing`](#single-cell-lineage-tracing) | bioskills | Reconstructs single-cell lineage trees and clonal relationships from CRISPR/Cas9 scars, static expressed barcodes (LARRY/CellTag), or somatic mtDNA mutations u… |

## `crispr-screen-triage` {#crispr-screen-triage}

Deterministic CRISPR screen hit ranking from local guide-level count tables

- **Source:** [ClawBio; ClawBio contributors @ `1fcecb7`](https://github.com/ClawBio/ClawBio/tree/1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d/skills/crispr-screen-triage)
- **License:** `MIT`
- **Runs with:** commands `python3`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info crispr-screen-triage
dnakit skills install crispr-screen-triage --agent claude-code   # or --agent codex / cursor
```

## `crispr-screens-bagel-essentiality` {#crispr-screens-bagel-essentiality}

Identifies essential genes from CRISPR-Cas9 fitness screens using BAGEL2 (Kim &amp; Hart 2021 Genome Med), a Bayesian classifier scoring per-gene Bayes Factors via log-likelihood ratios over per-sgRNA fold changes, calibrated against CEGv2 core-essentials (Hart 2017 G3, ~684 genes) and NEGv1 non-essentials (Hart 2014, ~927 genes). Covers the fc + bf + pr workflow, the linear-extrapolation improvement over BAGEL1 truncation, multi-target off-target correction, tumor-suppressor sensitivity (BAGEL2 detects enrichment), and BF calibration (BF &gt;6 ≈ 90% posterior per Hart 2017; ~5% FDR by BAGEL convention). Use when classifying essential vs non-essential genes, calibrating BAGEL2 thresholds against PR curves, identifying tumor suppressors alongside essentials, comparing BAGEL2 hits to MAGeCK / drugZ, or generating publication-quality essentiality calls.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/crispr-screens/bagel-essentiality)
- **License:** `MIT`
- **Runs with:** commands `BAGEL2`
- **Also known as:** `bio-crispr-screens-bagel-essentiality`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info crispr-screens-bagel-essentiality
dnakit skills install crispr-screens-bagel-essentiality --agent claude-code   # or --agent codex / cursor
```

## `crispr-screens-base-editing-analysis` {#crispr-screens-base-editing-analysis}

Analyzes base-editing screens for variant function. Covers library design (Hanna 2021 ClinVar-scale CBE screen benchmarked on BRCA1/2, Cuella-Martin 2021 DDR saturation), CBE vs ABE chemistry choice (BE3/BE4 vs ABE7.10/ABE8.20/ABE8e), editing-window math (positions 4-8 from PAM-distal end; 4-7 for ABE7.10), bystander-edit quantification and the variant-call ambiguity it creates, sgRNA-efficiency filtering before hit calling, indel byproduct interpretation, the substitution-vs-indel diagnostic, variant annotation against ClinVar / COSMIC, and the Broad be-validation-pipeline. Use when designing a BE variant screen, choosing CBE vs ABE for a specific edit, interpreting bystander-confounded hits, distinguishing functional signal from indel artifact, integrating CRISPResso2 output with screen scoring, or deciding BE vs PE for SNV installation.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/crispr-screens/base-editing-analysis)
- **License:** `MIT`
- **Also known as:** `bio-crispr-screens-base-editing-analysis`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info crispr-screens-base-editing-analysis
dnakit skills install crispr-screens-base-editing-analysis --agent claude-code   # or --agent codex / cursor
```

## `crispr-screens-batch-correction` {#crispr-screens-batch-correction}

Batch effect correction for CRISPR screens covering ComBat empirical-Bayes, RUV, SVA, control-sgRNA normalization, and the model-based alternative of including batch as a covariate in MAGeCK MLE or Chronos. Covers screen-specific batch sources (passage cohort, library lot, infection day, sequencing run, Cas9 lot, FBS lot), PCA + variance-decomposition diagnostic to decide if correction is needed, when correction harms biology by over-correcting condition into batch, limma removeBatchEffect for visualization-only correction, and relationship to multi-condition design matrices. Use when combining screens for joint analysis, when passage cohort confounds biology, when DepMap-style panels need Chronos with batch covariates, when picking ComBat vs RUV, or when correction harms biology and should be replaced with explicit covariate modeling.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/crispr-screens/batch-correction)
- **License:** `MIT`
- **Also known as:** `bio-crispr-screens-batch-correction`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info crispr-screens-batch-correction
dnakit skills install crispr-screens-batch-correction --agent claude-code   # or --agent codex / cursor
```

## `crispr-screens-combinatorial-screens` {#crispr-screens-combinatorial-screens}

Designs and analyzes combinatorial CRISPR screens covering paired-Cas9 (Big Papi, Najm 2018), enhanced AsCas12a multiplex (enCas12a, DeWeirdt 2021), in4mer 4-guide-array Cas12a (Esmaeili Anvar N et al 2024 Nat Commun 15:3577) and the Inzolia paralog-pair library, paralog-buffering detection (Dede 2020 Genome Biol; Thompson 2021 Nat Commun 12:1302), genetic-interaction (GI) scoring as observed_double_LFC minus expected_additive_double_LFC, synthetic-lethal and synthetic-rescue interaction interpretation, the half-of-essentiality buffered by paralogs phenomenon, multiplex screen statistical analysis with MAGeCK MLE interaction terms, and the relationship to single-cell combinatorial Perturb-seq. Use when designing a paralog or pathway-pair screen, choosing between paired-Cas9 (Big Papi) and Cas12a multiplex (Inzolia), interpreting genetic interaction scores, identifying synthetic-lethal targets for drug development, or scaling beyond single-gene CRISPR screens.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/crispr-screens/combinatorial-screens)
- **License:** `MIT`
- **Also known as:** `bio-crispr-screens-combinatorial-screens`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info crispr-screens-combinatorial-screens
dnakit skills install crispr-screens-combinatorial-screens --agent claude-code   # or --agent codex / cursor
```

## `crispr-screens-copy-number-correction` {#crispr-screens-copy-number-correction}

Corrects the gene-independent copy-number artifact in CRISPR-Cas9 screens (Aguirre 2016 / Munoz 2016 Cancer Discov) where amplified loci appear essential from DNA-damage burden of simultaneous cuts. Covers the gene-independent DNA-damage / G2-arrest mechanism, CRISPRcleanR (Iorio 2018) unsupervised pre-hoc correction, CERES (Meyers 2017) joint CN + gene-effect model, Chronos (Dempster 2021) DepMap-standard population-dynamics + CN model with lowest residual bias, the decision tree by data availability, the Spearman LFC-vs-CN diagnostic, focal-amplification examples (ERBB2 in HER2+, MYC in colorectal, FGFR1 in head and neck), and CRISPRi/a alternatives that bypass the artifact. Use when screening cancer cell lines, diagnosing essentiality at amplified loci, choosing CRISPRcleanR / CERES / Chronos, deciding whether CN correction is needed before MAGeCK / BAGEL2 / drugZ, or switching from Cas9 to CRISPRi.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/crispr-screens/copy-number-correction)
- **License:** `MIT`
- **Also known as:** `bio-crispr-screens-copy-number-correction`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info crispr-screens-copy-number-correction
dnakit skills install crispr-screens-copy-number-correction --agent claude-code   # or --agent codex / cursor
```

## `crispr-screens-crispresso-editing` {#crispr-screens-crispresso-editing}

Quantifies CRISPR editing outcomes with CRISPResso2 (Clement 2019 Nat Biotechnol) across Cas9-nuclease (indels, HDR), CBE and ABE base editors (target conversion + bystander), and prime editor (pegRNA-templated) modes. Covers single-amplicon (CRISPResso), multi-sample batch (CRISPRessoBatch), pooled-amplicon (CRISPRessoPooled), WGS off-target (CRISPRessoWGS), and sample-comparison (CRISPRessoCompare) workflows; quantification-window math that controls what is called edited; substitution-vs-indel diagnostic to distinguish BE from Cas9 contamination; MMEJ deletion pattern interpretation; allele-frequency tables; and failure modes from amplicon misalignment or contamination. Use when quantifying editing from amplicon sequencing, choosing CRISPResso mode by design, distinguishing intended edits from bystanders and indel byproducts, debugging low-alignment runs, or generating publication-grade editing reports.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/crispr-screens/crispresso-editing)
- **License:** `MIT`
- **Runs with:** commands `CRISPResso2`
- **Also known as:** `bio-crispr-screens-crispresso-editing`, `crispresso`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info crispr-screens-crispresso-editing
dnakit skills install crispr-screens-crispresso-editing --agent claude-code   # or --agent codex / cursor
```

## `crispr-screens-drugz-chemogenomic` {#crispr-screens-drugz-chemogenomic}

Analyzes CRISPR drug-modifier (chemogenomic) screens with drugZ (Colic et al. 2019 Genome Med), a bidirectional Z-score method that identifies synthetic-lethal sensitizing genes and resistance-conferring suppressor genes from vehicle vs drug comparisons. Covers vehicle-anchored design (not Day-0), the bidirectional Z math giving greater sensitivity to small-effect hits than MAGeCK / STARS / edgeR / RIGER on drug screens, per-gene sumZ and normZ, synth (sensitizer) vs supp (suppressor) FDR, multi-dose handling, integration with control sgRNAs, and comparison with MAGeCK MLE with dose covariate. Use when running a drug-modifier CRISPR screen, identifying sensitizing or resistance genes for a drug candidate, choosing drugZ vs MAGeCK MLE for chemogenomic analysis, troubleshooting low-effect drug screens where MAGeCK lacks sensitivity, or designing a drug-screen layout (vehicle vs drug arms).

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/crispr-screens/drugz-chemogenomic)
- **License:** `MIT`
- **Runs with:** commands `drugZ`
- **Also known as:** `bio-crispr-screens-drugz-chemogenomic`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info crispr-screens-drugz-chemogenomic
dnakit skills install crispr-screens-drugz-chemogenomic --agent claude-code   # or --agent codex / cursor
```

## `crispr-screens-hit-calling` {#crispr-screens-hit-calling}

Cross-method decision tree for calling hits in pooled CRISPR screens. Catalogs statistical models (MAGeCK RRA, MAGeCK MLE, BAGEL2, drugZ, JACKS, Chronos, CERES), experimental designs each is built for, failure modes outside design domain, reconciliation when methods disagree, multiple-testing and effect-size thresholds, the order of operations (count -&gt; QC -&gt; CN-correct -&gt; hit-call -&gt; validate), the second-best-sgRNA conservative rule, and consensus-hit strategy. Use when choosing among MAGeCK / BAGEL2 / drugZ / JACKS / Chronos for a given design, reconciling disagreement across two or three methods on the same screen, deciding whether to require consensus, gating downstream validation by hit-confidence tier, or interpreting unstable hit lists across reruns.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/crispr-screens/hit-calling)
- **License:** `MIT`
- **Also known as:** `bio-crispr-screens-hit-calling`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info crispr-screens-hit-calling
dnakit skills install crispr-screens-hit-calling --agent claude-code   # or --agent codex / cursor
```

## `crispr-screens-in-vivo-screens` {#crispr-screens-in-vivo-screens}

Designs and analyzes in vivo CRISPR screens in animal tumor models, organoids, and immune-cell adoptive transfers. Covers bottleneck math (250x cells/sgRNA requires ~25M cells implanted; impossible for most syngeneic models, forcing focused libraries), focused library design (Manguso 2017 Nature 547:413 immune screen; Chen 2015 tumor screens), CRISPR-StAR intrinsic-control screening (Uijttewaal 2025 Nat Biotechnol 43:1848), clonal-dynamics-limited detection, tumor-explant DNA recovery, syngeneic vs xenograft vs PDX considerations, and the relationship to downstream MAGeCK / drugZ analysis. Use when designing in vivo CRISPR screens for tumor / immune / metastasis biology, choosing focused vs genome-wide for animal models, addressing bottleneck-induced clonal collapse, picking the syngeneic / xenograft / PDX model, integrating in vivo with in vitro results, or applying CRISPR-StAR for animal experiments.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/crispr-screens/in-vivo-screens)
- **License:** `MIT`
- **Also known as:** `bio-crispr-screens-in-vivo-screens`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info crispr-screens-in-vivo-screens
dnakit skills install crispr-screens-in-vivo-screens --agent claude-code   # or --agent codex / cursor
```

## `crispr-screens-jacks-analysis` {#crispr-screens-jacks-analysis}

Runs JACKS (Joint Analysis of CRISPR/Cas9 Knockout Screens; Allen et al 2019 Genome Research) which models per-sgRNA log-fold-change as the product of a treatment-dependent gene-essentiality term and a treatment-independent guide-efficacy term. Covers the Bayesian decomposition math, the hierarchical efficacy prior shared across screens performed with the same library, when JACKS outperforms MAGeCK (multi-screen joint analysis, libraries with broad efficacy variance) and when it does not (single screen, novel libraries with no prior efficacy), library-reuse efficacy transfer, downstream essentiality interpretation, and the 2.5x sample-size reduction enabled by efficacy-aware testing. Use when running multiple screens with the same library, when guide-level noise is suspected to dominate per-gene signal, when reusing published essentiality reference screens for efficacy priors, or when comparing screens performed across cell lines that share library but differ biologically.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/crispr-screens/jacks-analysis)
- **License:** `MIT`
- **Runs with:** Python modules `JACKS`, `dnakit`
- **Also known as:** `bio-crispr-screens-jacks-analysis`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info crispr-screens-jacks-analysis
dnakit skills install crispr-screens-jacks-analysis --agent claude-code   # or --agent codex / cursor
```

## `crispr-screens-library-design` {#crispr-screens-library-design}

Designs pooled sgRNA libraries for CRISPR knockout, interference (CRISPRi), activation (CRISPRa), Cas12a multiplex, base-editor, and prime-editor screens. Covers on-target scoring (Rule Set 2, Azimuth, DeepSpCas9, CRISPRon), off-target scoring (CFD, MIT), TSS-relative positioning for CRISPRi/a (Horlbeck, Dolcetto, Calabrese), PAM-variant chemistries, control-guide composition, oligo cloning architecture, and library QC. Use when choosing a genome-wide library (GeCKOv2 vs Avana vs Brunello vs TKOv3 vs Inzolia), designing a focused or paralog-focused custom library, picking CRISPRi vs CRISPRa TSS windows, deciding control-guide proportions, or diagnosing library skew and dropout in a freshly cloned pool.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/crispr-screens/library-design)
- **License:** `MIT`
- **Also known as:** `bio-crispr-screens-library-design`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info crispr-screens-library-design
dnakit skills install crispr-screens-library-design --agent claude-code   # or --agent codex / cursor
```

## `crispr-screens-mageck-analysis` {#crispr-screens-mageck-analysis}

Analyzes pooled CRISPR screens with MAGeCK (Li et al 2014), covering count generation (mageck count), the RRA two-condition workflow (mageck test using alpha-RRA over per-sgRNA negative-binomial p-values), the MLE multi-condition workflow (mageck mle with explicit design matrix and beta-score output), normalization choice (median vs total vs control-sgRNA vs spike-in), sgRNA efficiency injection, paired-sample testing, time-course design, drug-screen versus dropout-screen design matrices, MAGeCKFlute and MAGeCK-VISPR downstream visualization, and decision logic for when to use MAGeCK vs JACKS / BAGEL2 / drugZ / Chronos. Use when running a fresh CRISPR screen analysis, picking RRA vs MLE for the experimental design, choosing a normalization method from QC signatures, debugging MLE convergence failure or NaN beta scores, comparing MAGeCK output across tools, or building a batch-aware multi-cell-line / multi-condition MLE design matrix.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/crispr-screens/mageck-analysis)
- **License:** `MIT`
- **Runs with:** commands `MAGeCK`
- **Also known as:** `bio-crispr-screens-mageck-analysis`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info crispr-screens-mageck-analysis
dnakit skills install crispr-screens-mageck-analysis --agent claude-code   # or --agent codex / cursor
```

## `crispr-screens-perturb-seq-analysis` {#crispr-screens-perturb-seq-analysis}

Analyzes single-cell pooled CRISPR screens (Perturb-seq, CROP-seq, Perturb-CITE-seq, ECCITE-seq, multiome) where each cell carries an sgRNA and a scRNA-seq / surface-protein / chromatin readout. Covers experimental design (direct-capture Perturb-seq Dixit 2016 vs CROP-seq 3'UTR-barcoded Datlinger 2017 vs ECCITE-seq vs Multiome), MOI for sgRNA assignment, escaper-cell filtering (Mixscape, Papalexi 2021), SCEPTRE NB GLM + permutation for low-MOI (Barry 2024 Genome Biol 25:124), the Pertpy framework, factor decomposition, genome-scale Perturb-seq (Replogle 2022 Cell, 2.5M cells), and per-perturbation single-cell DE. Use when running a single-cell CRISPR screen, choosing direct-capture vs CROP-seq architecture, filtering escaper cells, performing single-cell DE, integrating Perturb-seq with pathway analysis, scaling to GW CRISPRi via Replogle protocol, or analyzing multi-omics screens.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/crispr-screens/perturb-seq-analysis)
- **License:** `MIT`
- **Runs with:** Python modules `Pertpy`, `dnakit`
- **Also known as:** `bio-crispr-screens-perturb-seq-analysis`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info crispr-screens-perturb-seq-analysis
dnakit skills install crispr-screens-perturb-seq-analysis --agent claude-code   # or --agent codex / cursor
```

## `crispr-screens-prime-editing-screens` {#crispr-screens-prime-editing-screens}

Designs and analyzes pooled prime-editor (PE) screens for installing precise genetic variants without bystander confounding. Covers pegRNA design with PRIDICT and PRIDICT2 for predicting per-pegRNA editing efficiency, pegRNA architecture (spacer + scaffold + PBS + RTT), PE2/PE3/PE3b/PEmax variants, MOSAIC in situ saturation mutagenesis, the PRIME pooled-screen methodology (Ren 2023; ~3,699 ClinVar variant screens), chromatin context as a major locus-level determinant of PE efficiency, scaffold-incorporation and indel byproduct quantification with CRISPResso2, and the cross-modal validation strategy of PE + base-editor screens for variant function. Use when designing a pegRNA library for variant installation, choosing between BE and PE for a specific edit, predicting pegRNA efficiency before library synthesis, analyzing PE screen output, distinguishing intended-edit from scaffold-incorporation, or scaling PE screens to thousands of variants.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/crispr-screens/prime-editing-screens)
- **License:** `MIT`
- **Also known as:** `bio-crispr-screens-prime-editing-screens`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info crispr-screens-prime-editing-screens
dnakit skills install crispr-screens-prime-editing-screens --agent claude-code   # or --agent codex / cursor
```

## `crispr-screens-screen-qc` {#crispr-screens-screen-qc}

Quality control for pooled CRISPR screens covering library representation, Gini index, log-skew, replicate Pearson and Spearman concordance, essentialome precision-recall AUC against CEGv2 (Hart 2017), Cas9 cut-toxicity diagnostics, copy-number amplicon detection (Aguirre 2016 / Munoz 2016), bottleneck propagation through plasmid pool, infection, selection, and endpoint stages, MOI verification, and DepMap-style screen-quality scoring. Use when assessing screen quality before hit calling, deciding whether to repeat or rescue a screen, diagnosing low-confidence hits, choosing between MAGeCK / BAGEL2 / Chronos based on quality grade, picking a normalization strategy from QC signatures, or evaluating whether an in-vivo screen retained adequate library complexity.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/crispr-screens/screen-qc)
- **License:** `MIT`
- **Also known as:** `bio-crispr-screens-screen-qc`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info crispr-screens-screen-qc
dnakit skills install crispr-screens-screen-qc --agent claude-code   # or --agent codex / cursor
```

## `functional-genomics-screens` {#functional-genomics-screens}

Interpret hits from CRISPR-KO/CRISPRi/shRNA screens by integrating DepMap essentiality, gnomAD constraint scores, pathway context (Reactome, STRING), druggability (DGIdb), and clinical evidence (CIViC, COSMIC). Use for screen-hit prioritization, essentiality ranking, and turning a list of screen hits into a prioritized target shortlist.

- **Source:** [ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-functional-genomics-screens)
- **License:** `Apache-2.0`
- **Also known as:** `tooluniverse-functional-genomics-screens`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info functional-genomics-screens
dnakit skills install functional-genomics-screens --agent claude-code   # or --agent codex / cursor
```

## `functional-screen-analysis` {#functional-screen-analysis}

Analyze ENCODE functional genomics screens including CRISPR screens, MPRA (Massively Parallel Reporter Assays), and STARR-seq. Find screen data in ENCODE, process results, identify functional elements, and integrate with epigenomic annotations.

- **Source:** [ammawla/encode-toolkit contributors @ `b1d90bc`](https://github.com/ammawla/encode-toolkit/tree/b1d90bc9e7ead40d49b8c89f0b9b9e7b8e592a5b/skills/functional-screen-analysis)
- **License:** `AGPL-3.0-only`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info functional-screen-analysis
dnakit skills install functional-screen-analysis --agent claude-code   # or --agent codex / cursor
```

## `genome-engineering-base-editing-design` {#genome-engineering-base-editing-design}

Designs cytosine (CBE, C-to-T) and adenine (ABE, A-to-G) base-editor guides by positioning the target base at the activity-peak of the editing window (protospacer positions ~5-7, PAM-distal numbering), minimizing bystander edits for product purity, reading dinucleotide context (APOBEC1 TC favored / GC disfavored), and selecting the editor variant (BE4max, ABEmax, ABE8e, YE1/SECURE, TadCBE, CGBE, SpG/SpRY-BE). Covers knockout by premature stop (CRISPR-STOP/iSTOP) and splice-site disruption, the three off-target classes (Cas-dependent, Cas-independent DNA, RNA), outcome prediction (BE-Hive/DeepBE), and the base-vs-prime-vs-HDR decision. Use when installing a transition mutation without a double-strand break, knocking out a gene without indels, or choosing CBE vs ABE. Generic guide scoring, prime editing, and HDR donors are separate skills.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-engineering/base-editing-design)
- **License:** `MIT`
- **Runs with:** Python modules `BioPython`, `dnakit`
- **Also known as:** `bio-genome-engineering-base-editing-design`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-engineering-base-editing-design
dnakit skills install genome-engineering-base-editing-design --agent claude-code   # or --agent codex / cursor
```

## `genome-engineering-grna-design` {#genome-engineering-grna-design}

Designs and ranks guide RNAs (sgRNAs) for CRISPR-Cas9/Cas12a gene knockout by scanning a target for PAM sites (NGG SpCas9, NNGRRT SaCas9, TTTV Cas12a, NG SpCas9-NG, near-PAMless SpRY), enumerating candidate spacers, applying hard filters (Pol-III TTTT terminator, 5' G, GC), ranking on-target activity with the context-appropriate model (Rule Set 2/Azimuth for U6/lentiviral, CRISPRscan for T7/embryo, DeepHF for high-fidelity variants, DeepCpf1 for Cas12a), and predicting the indel/frameshift outcome (Bae out-of-frame score, inDelphi, FORECasT, Lindel). Use when selecting sgRNAs to knock out a gene, choosing a nuclease/PAM for a constrained locus, picking which exon to target, or shortlisting guides before an off-target check. Off-target specificity, base/prime editing, and HDR donors are separate skills.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-engineering/grna-design)
- **License:** `MIT`
- **Runs with:** Python modules `CRISPOR`, `dnakit`
- **Also known as:** `bio-genome-engineering-grna-design`, `sgrna-design-guide`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-engineering-grna-design
dnakit skills install genome-engineering-grna-design --agent claude-code   # or --agent codex / cursor
```

## `genome-engineering-hdr-template-design` {#genome-engineering-hdr-template-design}

Designs donor/repair templates for precise CRISPR knock-ins -- choosing the format (ssODN, long-ssDNA/Easi-CRISPR, dsDNA/plasmid, AAV6), sizing homology arms, placing the cut within ~10 bp of the edit, and adding a mandatory codon-checked blocking (PAM/seed) mutation so the edited allele is not re-cut. Frames the HDR-vs-NHEJ-vs-MMEJ pathway competition, the MMEJ (PITCh) and homology-independent (HITI/HMEJ) alternatives for post-mitotic cells, ssODN strand/asymmetry choice, phosphorothioate end-protection, and ranked HDR enhancers. Use when designing a donor for a point mutation, epitope/fluorophore tag, allele replacement, or knock-in, or when HDR efficiency is low. Guide design and base/prime editing are separate skills.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-engineering/hdr-template-design)
- **License:** `MIT`
- **Also known as:** `bio-genome-engineering-hdr-template-design`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-engineering-hdr-template-design
dnakit skills install genome-engineering-hdr-template-design --agent claude-code   # or --agent codex / cursor
```

## `genome-engineering-off-target-prediction` {#genome-engineering-off-target-prediction}

Nominates and assesses CRISPR off-target sites genome-wide. Enumerates candidate sites by mismatch and bulge tolerance with Cas-OFFinder/CRISPRitz, ranks them with the published CFD score (SpCas9-only, relative ranker) or MIT/CRISTA/energy models, runs variant-aware screening against gnomAD/individual genomes (CRISPRme), and frames the empirical genome-wide discovery assays (GUIDE-seq, CIRCLE-seq, CHANGE-seq, DISCOVER-seq, Digenome-seq) and high-fidelity nuclease choice (HiFi Cas9, Sniper-Cas9, eSpCas9, SpCas9-HF1). Use when assessing guide RNA specificity, choosing among candidate guides, screening a therapeutic guide against population variation, or planning empirical off-target validation. Distinguishes predicted vs detected vs validated. On-target activity scoring and deaminase (Cas-independent) base/prime-editor off-targets are separate skills.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-engineering/off-target-prediction)
- **License:** `MIT`
- **Also known as:** `bio-genome-engineering-off-target-prediction`, `score-variant-aware-crispr-offtargets`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-engineering-off-target-prediction
dnakit skills install genome-engineering-off-target-prediction --agent claude-code   # or --agent codex / cursor
```

## `genome-engineering-prime-editing-design` {#genome-engineering-prime-editing-design}

Designs pegRNAs and nicking guides for prime editing (PE) -- choosing the nick/strand, tuning the primer-binding site (PBS) and reverse-transcription template (RTT) as a per-locus panel, selecting the PE system (PE2/PE3/PE3b/PE4/PE5/PEmax/PE7), adding MMR-evading and PAM-disrupting silent edits, appending epegRNA 3' motifs (tevopreQ1/mpknot), and ranking with PRIDICT/DeepPrime. Covers twinPE/PASTE for large insertions and the prime-vs-base-editing decision. Use when designing a scarless point mutation, small insertion/deletion, or any of the 12 base conversions without a double-strand break, when efficiency is low and MMR inhibition or pegRNA stabilization is needed, or when routing a large insertion to an integrase method. Generic guide scoring and base editing are separate skills.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-engineering/prime-editing-design)
- **License:** `MIT`
- **Also known as:** `bio-genome-engineering-prime-editing-design`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-engineering-prime-editing-design
dnakit skills install genome-engineering-prime-editing-design --agent claude-code   # or --agent codex / cursor
```

## `molecular-cloning` {#molecular-cloning}

Molecular cloning, in both directions. DESIGN — Gibson Assembly (overlap design for seamless multi-fragment joining) and Golden Gate Assembly (Type IIS / BsaI / BbsI / Esp3I / BsmBI / SapI design with unique 4-bp fusion overhangs). ANALYSIS — work out what an existing reaction produces: given input plasmid sequences and an enzyme, digest them, join the fragments by their overhangs, and identify features of the product (expressed ORF, gRNA spacer and its target gene). Use when you need to plan how to join DNA fragments into a construct, design assembly overlaps/overhangs, decide between cloning methods, or determine the product of a stated Gibson/Golden Gate reaction. Covers the domestication (internal-site removal), overhang-uniqueness, and overlap-Tm rules. For PCR primers to generate the fragments, see tooluniverse-primer-design.

- **Source:** [ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-molecular-cloning)
- **License:** `Apache-2.0`
- **Also known as:** `tooluniverse-molecular-cloning`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info molecular-cloning
dnakit skills install molecular-cloning --agent claude-code   # or --agent codex / cursor
```

## `plannotate-plasmid-annotation` {#plannotate-plasmid-annotation}

Auto-annotate plasmids with features (promoters, terminators, resistance, origins, tags, fluorescent proteins) via BLAST against curated DBs (Addgene, fpbase, SnapGene). FASTA or raw sequence in; annotated GenBank, interactive HTML maps, CSV tables out. Handles circular topology. Use to verify synthetic constructs, prep Addgene submissions, share maps, or batch-annotate cloning libraries.

- **Source:** [jaechang-hits; credited K-Dense-AI and Biomni contributors @ `fe505ca`](https://github.com/jaechang-hits/SciAgent-Skills/tree/fe505cae14d20b6c33be2e49666425be98f005bb/skills/molecular-biology/plannotate-plasmid-annotation)
- **License:** `CC-BY-4.0`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info plannotate-plasmid-annotation
dnakit skills install plannotate-plasmid-annotation --agent claude-code   # or --agent codex / cursor
```

## `plant-genome-guided-primer-design` {#plant-genome-guided-primer-design}

A genome-guided workflow for plant gene identity validation, homolog analysis, and experimental primer design. The workflow prevents incorrect primer design caused by gene ID misannotation, duplicated genes, paralogs, and polyploid genomes.

- **Source:** [Li-Chia-ching/Plant-qPCR-Agent-Skill contributors @ `23acf68`](https://github.com/Li-Chia-ching/Plant-qPCR-Agent-Skill/tree/23acf68cc4cefe67f8f3b5f5057b0eabf0ac9a4e/)
- **License:** `Apache-2.0`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info plant-genome-guided-primer-design
dnakit skills install plant-genome-guided-primer-design --agent claude-code   # or --agent codex / cursor
```

## `polyploid-gene-family-cloning` {#polyploid-gene-family-cloning}

AI Skill for genome-aware polyploid gene family identification, candidate prioritization, validated full-length CDS cloning, construction-intent locking, homoeolog-aware primer design, terminal ORF integrity control, and end-to-end seamless cloning reconstruction.

- **Source:** [Li-Chia-ching/Polyploid-GeneFamily-AI-Pipeline contributors @ `ecce861`](https://github.com/Li-Chia-ching/Polyploid-GeneFamily-AI-Pipeline/tree/ecce861e5215e25cc180f1381d73a8751c5257bd/)
- **License:** `MIT`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info polyploid-gene-family-cloning
dnakit skills install polyploid-gene-family-cloning --agent claude-code   # or --agent codex / cursor
```

## `primer-design-primer-basics` {#primer-design-primer-basics}

Designs and ranks PCR primer pairs for a target template with primer3-py (design_primers), returning pairs with nearest-neighbor Tm, GC, product size, and complementarity scores. Covers why primer3 is a LOCAL weighted-penalty minimizer over the single template supplied (so PRIMER_PAIR_0 is the lowest-penalty pair under the given bounds, never a genome-specificity guarantee), why Tm is a salt/concentration-dependent SantaLucia prediction not a fixed property, why the two primers must be Tm-matched, the seq_args/global_args tag semantics (SEQUENCE_TARGET/INCLUDED/EXCLUDED/OVERLAP_JUNCTION/FORCE_*, 0-based \[start,length\]), 3'-end and GC-clamp mechanism, 5'-tail handling, masking SNPs under 3' ends, and diagnosing zero-pair runs. Use when designing standard PCR, cloning, genotyping, or sequencing primers, flanking a target, or screening pairs by Tm/size/GC. Genome off-target checking is primer-specificity; dimers/hairpins primer-validation; qPCR and probes qpcr-primers.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/primer-design/primer-basics)
- **License:** `MIT`
- **Also known as:** `bio-primer-design-primer-basics`, `pcr-primer-design`, `tooluniverse-primer-design`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info primer-design-primer-basics
dnakit skills install primer-design-primer-basics --agent claude-code   # or --agent codex / cursor
```

## `primer-design-primer-specificity` {#primer-design-primer-specificity}

Checks whether a PCR primer PAIR amplifies only the intended target genome-wide, using pair-aware in-silico PCR (MFEprimer-3.0, UCSC isPcr, NCBI Primer-BLAST) plus a primer3-py 3'-end-stability prefilter, against the correct database. Covers why plain BLAST is the wrong tool (it scores per-primer similarity, blind to 3'-terminal anchoring and to whether the two primers form a convergent amplicon in range), why a single 3'-terminal mismatch suppresses amplification while internal mismatches are tolerated, why intron-spanning RT-qPCR is defeated by processed pseudogenes that force a GENOME search not transcriptome-only, how to read a Primer-BLAST report (empty unintended-products means none passed its filter, not none exist), and that in-silico checking reduces but never replaces empirical validation. Use when confirming specificity, screening off-target amplicons, avoiding paralog/pseudogene hits, or checking SNPs under the 3' end. Design is primer-basics; dimers primer-validation; alignment read-alignment.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/primer-design/primer-specificity)
- **License:** `MIT`
- **Also known as:** `bio-primer-design-primer-specificity`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info primer-design-primer-specificity
dnakit skills install primer-design-primer-specificity --agent claude-code   # or --agent codex / cursor
```

## `primer-design-primer-validation` {#primer-design-primer-validation}

Validates chosen PCR/qPCR oligos for intramolecular thermodynamic liabilities with primer3-py - hairpins, self-dimers, cross-dimers (calc_hairpin/homodimer/heterodimer), and 3'-end stability (calc_end_stability) - returning ThermoResult dG/Tm and ASCII structures. Covers why a "dimer-free" verdict is a PREDICTION at the supplied salt/Mg/dNTP/oligo conditions and temp_c (so the same primer is fine or dimer-prone depending on conditions), why a 3'-END dimer or hairpin is the lethal class (polymerase-extendable into primer-dimer) so structures are ranked by dG at the annealing temperature and 3'-end involvement rather than global Tm, that ThermoResult dG is in cal/mol not kcal/mol, and that .structure_found must gate the numbers. Use when checking primer pairs before ordering, troubleshooting primer-dimers or smears, or screening oligos for secondary structure. Genome off-target/mispriming is primer-specificity; design is primer-basics; probe assays are qpcr-primers.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/primer-design/primer-validation)
- **License:** `MIT`
- **Also known as:** `bio-primer-design-primer-validation`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info primer-design-primer-validation
dnakit skills install primer-design-primer-validation --agent claude-code   # or --agent codex / cursor
```

## `primer-design-qpcr-primers` {#primer-design-qpcr-primers}

Co-designs qPCR/RT-qPCR primers and hydrolysis (TaqMan) or molecular-beacon probes with primer3-py (PRIMER_PICK_INTERNAL_OLIGO, PRIMER_INTERNAL_* tags), for assays whose deliverable is a quantitative measurement device. Covers why amplification efficiency (90-110%, slope -3.6 to -3.1) and single-product specificity make the 2^-ddCq / Pfaffl math valid, why the short amplicon (70-150 bp), tight Tm, and zero-dimer requirement exist, the coupled probe rules (probe Tm 8-10 C above primers so it is bound when Taq's exonuclease cleaves it; no 5' G as it quenches the reporter; C-rich strand; primer3 has NO no-5'-G tag so enforce PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME=HNNNN), gDNA exclusion by exon-junction spanning AND why pseudogenes defeat it, SYBR melt-curve QC, and reference-gene validation (geNorm/NormFinder). Use when designing TaqMan/SYBR assays, exon-spanning primers, probes, or matched-efficiency multiplex panels. Genome specificity is primer-specificity; dimers primer-validation; standard PCR primer-basics.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/primer-design/qpcr-primers)
- **License:** `MIT`
- **Also known as:** `bio-primer-design-qpcr-primers`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info primer-design-qpcr-primers
dnakit skills install primer-design-qpcr-primers --agent claude-code   # or --agent codex / cursor
```

## `restriction-enzyme-selection` {#restriction-enzyme-selection}

Select restriction enzymes for cloning or diagnostics using Biopython Bio.Restriction. Finds enzymes by cut frequency, overhang type, recognition-site length, commercial availability, compatible ends, and methylation sensitivity, and identifies isoschizomers and compatible pairs. Use when choosing which enzymes to use to linearize a vector, drop in an insert, set up a diagnostic digest, or pick a methylation-insensitive enzyme.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/restriction-analysis/enzyme-selection)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.Restriction`, `dnakit`
- **Also known as:** `bio-restriction-enzyme-selection`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info restriction-enzyme-selection
dnakit skills install restriction-enzyme-selection --agent claude-code   # or --agent codex / cursor
```

## `restriction-fragment-analysis` {#restriction-fragment-analysis}

Predict restriction digest fragment sizes and gel patterns using Biopython Bio.Restriction. Computes fragment lengths and sequences for single and double digests on linear or circular DNA, and interprets them against an agarose gel. Use when predicting the fragments from a digest, planning a diagnostic digest to verify a clone, or matching observed gel bands to an expected pattern.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/restriction-analysis/fragment-analysis)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.Restriction`, `dnakit`
- **Also known as:** `bio-restriction-fragment-analysis`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info restriction-fragment-analysis
dnakit skills install restriction-fragment-analysis --agent claude-code   # or --agent codex / cursor
```

## `restriction-golden-gate-assembly` {#restriction-golden-gate-assembly}

Design and validate Type IIS scarless DNA assembly (Golden Gate, MoClo) using Biopython Bio.Restriction. Screens parts for internal BsaI/BsmBI/BbsI/SapI sites (domestication), previews the fusion overhangs a digest exposes, and validates a fusion-overhang set for distinctness and fidelity. Use when designing a Golden Gate or MoClo assembly, domesticating a part by removing internal Type IIS sites, or choosing and checking fusion overhangs for one-pot assembly.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/restriction-analysis/golden-gate-assembly)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.Restriction`, `dnakit`
- **Also known as:** `bio-restriction-golden-gate-assembly`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info restriction-golden-gate-assembly
dnakit skills install restriction-golden-gate-assembly --agent claude-code   # or --agent codex / cursor
```

## `restriction-mapping` {#restriction-mapping}

Build restriction maps showing enzyme cut positions and inter-site distances along DNA using Biopython Bio.Restriction. Produces text or graphical maps for linear and circular molecules, orders sites from single and double digests, and overlays GenBank features. Use when creating a restriction map of a sequence, ordering cut sites along a plasmid, or relating sites to annotated features.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/restriction-analysis/restriction-mapping)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.Restriction`, `dnakit`
- **Also known as:** `bio-restriction-mapping`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info restriction-mapping
dnakit skills install restriction-mapping --agent claude-code   # or --agent codex / cursor
```

## `restriction-sites` {#restriction-sites}

Find restriction enzyme cut sites in DNA sequences using Biopython Bio.Restriction. Searches single enzymes, batches, or commercial enzyme sets and returns cut positions for linear or circular DNA. Use when locating where one or more restriction enzymes cut a sequence, screening a sequence for the presence or absence of a site, or counting how often an enzyme cuts.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/restriction-analysis/restriction-sites)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.Restriction`, `dnakit`
- **Also known as:** `bio-restriction-sites`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info restriction-sites
dnakit skills install restriction-sites --agent claude-code   # or --agent codex / cursor
```

## `single-cell-lineage-tracing` {#single-cell-lineage-tracing}

Reconstructs single-cell lineage trees and clonal relationships from CRISPR/Cas9 scars, static expressed barcodes (LARRY/CellTag), or somatic mtDNA mutations using Cassiopeia, Startle, and CoSpar. Use when building a phylogeny from barcode scars, choosing a tree-reconstruction solver, handling homoplasy and dropout, grouping clones from mtDNA, integrating clone with transcriptomic state, or judging whether a state-based fate call is trustworthy.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/single-cell/lineage-tracing)
- **License:** `MIT`
- **Runs with:** Python modules `Cassiopeia`, `dnakit`
- **Also known as:** `bio-single-cell-lineage-tracing`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info single-cell-lineage-tracing
dnakit skills install single-cell-lineage-tracing --agent claude-code   # or --agent codex / cursor
```
