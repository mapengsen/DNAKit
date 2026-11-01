# 微生物与环境基因组

本类共 **67** 个 Skills。每个 Skill 是一份让 Agent 按步骤完成分析的说明，分析所需的软件、参考数据和 API 需要另外准备。搜索、安装和使用方法见 [Skills](../skills.md)；网页版的 Skills 标签可以直接阅读每个 Skill 的完整内容。

| Skill | 来源 | 简介 |
| --- | --- | --- |
| [`antismash`](#antismash) | pku-metagenomics | Detect and annotate biosynthetic gene clusters (BGCs) on genomes and MAGs. Use for natural-product genome mining, choosing genefinding, interpreting region HTM… |
| [`anvio`](#anvio) | pku-metagenomics | Anvi'o — interactive multi-omics platform for metagenomic assemblies, MAG refinement, visualization, and SNP/profile exploration. Use when: (1) Manual or human… |
| [`basalt`](#basalt) | pku-metagenomics | PKU-EMBL BASALT and BASALT-Air — recommended ensemble MAG recovery across one or more assemblies (multi-binner + refinement). Use this skill when: (1) Running… |
| [`bigscape`](#bigscape) | pku-metagenomics | Cluster and compare biosynthetic gene clusters across genomes with BiG-SCAPE (similarity networks of antiSMASH regions). Use after batch antiSMASH for BGC fami… |
| [`binchicken`](#binchicken) | pku-metagenomics | Bin Chicken — targeted coassembly planning for recovering low-abundance / novel MAGs by grouping metagenomes that share marker-gene signal (SingleM). Use when:… |
| [`binning-qc`](#binning-qc) | omics | Bin and refine metagenomic contigs, then assess MAG quality. Use when recovering genomes with QuickBin and checking completeness, contamination, and bin consis… |
| [`biosurfdb`](#biosurfdb) | pku-metagenomics | BioSurfDB — curated portal for biosurfactant production and hydrocarbon biodegradation genes/pathways with BLAST and related metagenome analysis helpers. Use f… |
| [`cami`](#cami) | pku-metagenomics | CAMI / CAMI II metagenome benchmarking datasets and evaluation stack (CAMISIM, AMBER, OPAL, CAMI Benchmarking Portal). Use when: (1) Downloading CAMI challenge… |
| [`checkm2`](#checkm2) | pku-metagenomics | ML-based MAG completeness and contamination estimation. Use this skill when working with checkm2: installing, choosing parameters, running commands, interpreti… |
| [`checkv`](#checkv) | pku-metagenomics | Quality assessment of metagenome-assembled viral genomes. Use this skill when working with checkv: installing, choosing parameters, running commands, interpret… |
| [`dbcan`](#dbcan) | pku-metagenomics | Annotate carbohydrate-active enzymes (CAZymes) with dbCAN3 / run_dbcan on genomes or protein sets, including CGC finding and glycan substrate prediction (dbCAN… |
| [`deepurify`](#deepurify) | pku-metagenomics | Deepurify — multi-modal deep language model for MAG decontamination by matching contigs to taxonomic lineages and splitting chimeric bins into higher-quality s… |
| [`defensefinder`](#defensefinder) | pku-metagenomics | Detect anti-phage defense systems in prokaryotic genomes. Use this skill when working with defensefinder: installing, choosing parameters, running commands, in… |
| [`disinfinder`](#disinfinder) | pku-metagenomics | DisinFinder — CGE disinfectant / biocide resistance gene screen via the DisinFinder database, typically run through ResFinder. Use for isolate biocide-resistan… |
| [`dom-biosurfactant`](#dom-biosurfactant) | pku-metagenomics | Route dissolved organic matter (DOM) / glycan substrate capacity, biosurfactant production, and hydrocarbon/plastic biodegradation annotation for shotgun metag… |
| [`dram`](#dram) | pku-metagenomics | Distill microbial metabolism annotations for MAGs/genomes into curated product tables (DRAM annotate + distill). Use for ecological functional potential, pathw… |
| [`ecological-genomics-biodiversity-metrics`](#ecological-genomics-biodiversity-metrics) | bioskills | Quantifies biodiversity from species abundance/incidence tables using Hill numbers (iNEXT) with coverage-based rarefaction-extrapolation (Chao &amp; Jost 2012), as… |
| [`ecological-genomics-community-ecology`](#ecological-genomics-community-ecology) | bioskills | Analyzes species-environment relationships with constrained ordination (CCA, RDA, db-RDA), variance partitioning, indicator species (indicspecies IndVal.g grou… |
| [`ecological-genomics-conservation-genetics`](#ecological-genomics-conservation-genetics) | bioskills | Assesses genetic health of populations for conservation with Ne estimation across time horizons (LDNe NeEstimator V2 option-file API + SNeP physical-linkage co… |
| [`ecological-genomics-edna-metabarcoding`](#ecological-genomics-edna-metabarcoding) | bioskills | Processes eDNA metabarcoding from raw paired-end reads to species tables, navigating ASV (DADA2, UNOISE3) vs OTU (swarm v2) decision (Callahan 2017 vs Schloss… |
| [`ecological-genomics-landscape-genomics`](#ecological-genomics-landscape-genomics) | bioskills | Tests genotype-environment associations and identifies adaptive loci while correcting for the four-confound landscape (structure, demography, background select… |
| [`ecological-genomics-species-delimitation`](#ecological-genomics-species-delimitation) | bioskills | Delimits putative species boundaries from molecular data within the de Queiroz 2007 unified-lineage framework using ASAP (Puillandre 2021 successor to ABGD), m… |
| [`element-cycling`](#element-cycling) | pku-metagenomics | Route biogeochemical / elemental cycling annotation for shotgun metagenomics and MAGs: carbon (C), nitrogen (N), sulfur (S), iron (Fe), phosphorus (P), hydroge… |
| [`epidemiological-genomics-amr-surveillance`](#epidemiological-genomics-amr-surveillance) | bioskills | Detects acquired antimicrobial-resistance determinants and chromosomal point-mutation resistance in bacterial assemblies using AMRFinderPlus, ResFinder 4.0 (ac… |
| [`epidemiological-genomics-pathogen-typing`](#epidemiological-genomics-pathogen-typing) | bioskills | Assigns isolate identity at the right resolution for the question -- ANI/Mash species triage, 7-locus MLST historical comparability, cgMLST/wgMLST outbreak res… |
| [`epidemiological-genomics-phylodynamics`](#epidemiological-genomics-phylodynamics) | bioskills | Estimates time-scaled phylogenies, molecular-clock rates, effective reproduction number R_e, and population dynamics from dated pathogen genomes using TreeTime… |
| [`epidemiological-genomics-transmission-inference`](#epidemiological-genomics-transmission-inference) | bioskills | Infers person-to-person transmission from pathogen genomes using outbreaker2, TransPhylo, phybreak, BadTrIP, SCOTTI, BEASTLIER, and SNP-distance / cluster-pick… |
| [`epidemiological-genomics-variant-surveillance`](#epidemiological-genomics-variant-surveillance) | bioskills | Assigns pathogen lineages (SARS-CoV-2 Pangolin UShER mode; Nextclade clade + QC; pango-designation alias resolution) and tracks variant frequencies over time u… |
| [`fegenie`](#fegenie) | pku-metagenomics | FeGenie — HMM-based identification of iron-related genes and neighborhoods in genomes/metagenome assemblies (acquisition, siderophores, Fe redox, storage, magn… |
| [`genomad`](#genomad) | pku-metagenomics | Identify plasmids, viruses, and other mobile genetic elements in assemblies or genomes. Use for plasmidome/virome mining, MGE-aware AMR context, and pre-CheckV… |
| [`gunc`](#gunc) | pku-metagenomics | Detect chimerism/contamination in prokaryotic genomes/MAGs. Use this skill when working with gunc: installing, choosing parameters, running commands, interpret… |
| [`hadeg`](#hadeg) | pku-metagenomics | HADEG — curated Hydrocarbon Aerobic Degradation Enzymes and Genes database (alkanes, alkenes, aromatics, plastics) plus biosurfactant production genes with exp… |
| [`insilicoseq`](#insilicoseq) | pku-metagenomics | InSilicoSeq — simulate realistic Illumina (meta)genomic reads for benchmarking and experimental design. Use when: (1) Generating Illumina-like test sets from g… |
| [`integronfinder`](#integronfinder) | pku-metagenomics | IntegronFinder — detect integrons (intI, attC, cassettes) in bacterial genomes/assemblies, including AMR cassette context. Use for integron-focused MGE mining.… |
| [`isescan`](#isescan) | pku-metagenomics | ISEScan — automated annotation of insertion sequence (IS) elements in prokaryotic genomes via profile HMMs. Use for IS family catalogues on MAGs or contigs as… |
| [`lorbin`](#lorbin) | pku-metagenomics | LorBin — unsupervised long-read metagenome binning with multiscale adaptive clustering (DBSCAN + BIRCH) and single-copy-gene evaluation. Use when: (1) Binning… |
| [`mebs`](#mebs) | pku-metagenomics | MEBS (Multigenomic Entropy-Based Score) — score genomes/metagenomes for biogeochemical cycle machinery via informative Pfam domains; strongest for sulfur (S),… |
| [`metabat2`](#metabat2) | pku-metagenomics | Adaptive contig binning from tetranucleotide composition and multi-sample coverage. Use for primary single-binner runs feeding DAS Tool, MetaWRAP, or BASALT en… |
| [`metabolic`](#metabolic) | pku-metagenomics | METABOLIC (AnantharamanLab) — high-throughput metabolic and biogeochemical trait profiling of MAGs/SAGs/isolates with C, N, S, and other elemental cycle diagra… |
| [`metagenomics-abundance`](#metagenomics-abundance) | bioskills | Turns shotgun classifier output into a defensible abundance table with Bracken Bayesian re-estimation, then compositional treatment (CLR, zero handling), libra… |
| [`metagenomics-amr-detection`](#metagenomics-amr-detection) | bioskills | Profiles the antimicrobial-resistance gene content (resistome) of shotgun metagenomes - read-based quantification with RGI bwt, AMR++/MEGARes, ARGs-OAP/SARG, d… |
| [`metagenomics-contamination-controls`](#metagenomics-contamination-controls) | bioskills | Cleans a shotgun metagenome of everything that is not the target community before profiling - host-read depletion (Hostile, bowtie2/T2T-CHM13), reagent/kitome… |
| [`metagenomics-functional-profiling`](#metagenomics-functional-profiling) | bioskills | Profiles the functional potential of shotgun metagenomes with HUMAnN 3's tiered search (MetaPhlAn prescreen, Bowtie2 pangenome, translated DIAMOND vs UniRef),… |
| [`metagenomics-kraken`](#metagenomics-kraken) | bioskills | Classifies shotgun metagenomic reads to taxa with Kraken2's minimizer/LCA matching against a chosen reference database, then hands off to Bracken for abundance… |
| [`metagenomics-metaphlan`](#metagenomics-metaphlan) | bioskills | Profiles shotgun metagenomes to species/SGB relative abundance with MetaPhlAn 4's clade-specific marker genes (bowtie2 short reads, minimap2 long reads). Cover… |
| [`metagenomics-strain-tracking`](#metagenomics-strain-tracking) | bioskills | Resolves and compares bacterial strains below the species level from shotgun metagenomes with inStrain (popANI/conANI microdiversity), StrainPhlAn (marker-SNV… |
| [`metagenomics-visualization`](#metagenomics-visualization) | bioskills | Turns a shotgun profiler table (MetaPhlAn relative abundance, Bracken counts, HUMAnN function tables) into honest figures and defensible community statistics w… |
| [`microbial-genome-characterization`](#microbial-genome-characterization) | tooluniverse | Genome-ASSEMBLY discovery, QC, and replicon mapping for any organism (bacteria, archaea, fungi, and beyond) using NCBI Datasets. Resolves an organism name or t… |
| [`microbial-mining`](#microbial-mining) | pku-metagenomics | Route microbial / genome mining after assembly or binning: BGCs (antiSMASH, BiG-SCAPE), AMR (RGI, DeepARG, ResFinder), mobile elements (geNomad, VirSorter2, Ch… |
| [`microbiome-amplicon-processing`](#microbiome-amplicon-processing) | bioskills | Infers exact amplicon sequence variants (ASVs) from demultiplexed 16S rRNA or ITS amplicon FASTQ with DADA2 - removing primers with cutadapt (--discard-untrimm… |
| [`microbiome-differential-abundance`](#microbiome-differential-abundance) | bioskills | Tests which individual taxa differ between groups on an amplicon ASV/feature table (phyloseq) using compositionally-aware methods - ALDEx2 (Dirichlet-MC CLR, c… |
| [`microbiome-diversity-analysis`](#microbiome-diversity-analysis) | bioskills | Alpha and beta diversity of an amplicon (16S/ITS) ASV/OTU community table - observed features, Shannon, Pielou evenness, Faith PD, Bray-Curtis, Jaccard, weight… |
| [`microbiome-functional-prediction`](#microbiome-functional-prediction) | bioskills | Predicts community functional POTENTIAL from 16S/ITS amplicon ASVs with PICRUSt2 (or q2-picrust2) by phylogenetic interpolation of reference-genome gene conten… |
| [`microbiome-qiime2-workflow`](#microbiome-qiime2-workflow) | bioskills | Operates the QIIME2 framework as the glue for an amplicon analysis - the .qza/.qzv artifact model, semantic types (FeatureTable\[Frequency\], SampleData\[PairedEn… |
| [`microbiome-research`](#microbiome-research) | tooluniverse | Microbiome research using MGnify, GTDB, ENA, OLS (ENVO biomes), and EuropePMC. Covers study discovery, taxonomic profiling, host-microbe interaction analysis,… |
| [`microbiome-taxonomy-assignment`](#microbiome-taxonomy-assignment) | bioskills | Assigns taxonomy to amplicon ASVs/OTUs (16S, ITS, 18S) with a classifier conditioned on a reference database and primer region - DADA2 assignTaxonomy + addSpec… |
| [`mob-suite`](#mob-suite) | pku-metagenomics | MOB-suite — plasmid reconstruction and typing from draft assemblies (replicon, relaxase/MOB type, mate-pair, transferability, mash clusters). Use for plasmid c… |
| [`mobileelementfinder`](#mobileelementfinder) | pku-metagenomics | MobileElementFinder — CGE detection of named mobile genetic elements (MITEs, IS, Tn, ComTn, ICE, IME, CIME, …) in assembled bacterial genomes, with optional AR… |
| [`mobileog`](#mobileog) | pku-metagenomics | mobileOG-db — manually curated protein families for bacterial MGE life-cycle functions (integration/excision, replication/recombination/repair, stability/defen… |
| [`pathogen-variant-surveillance`](#pathogen-variant-surveillance) | kdense | Query live pathogen genomic surveillance data through the GenSpectrum LAPIS API to find which viral lineages are circulating now, how fast they are growing, an… |
| [`plasmaag`](#plasmaag) | pku-metagenomics | PlasMAAG — recover plasmids and cellular genomes from metagenomes using assembly–alignment graphs and contrastive learning across samples. Use when: (1) Plasmi… |
| [`qsap`](#qsap) | pku-metagenomics | QSAP — quorum-sensing (QS) related gene annotation pipeline against the QSP database (DIAMOND and/or HMMER). Use when mining autoinducer / QS circuitry (AHL, A… |
| [`singlem`](#singlem) | pku-metagenomics | SingleM — species-level metagenome profiling from conserved windows in universal marker genes, including lineages without genomic representation; Sandpiper agg… |
| [`taxvamb`](#taxvamb) | pku-metagenomics | TaxVAMB — semisupervised bimodal VAE binning that integrates TNF, coabundance, and contig taxonomy for improved MAG recovery (incl. incomplete genomes and long… |
| [`vcontact`](#vcontact) | pku-metagenomics | vConTACT2 / vConTACT3 — gene-sharing network and ML-based hierarchical virus taxonomy for prokaryotic (and in v3, broader realm) viral genomes/contigs. Use aft… |
| [`viromics`](#viromics) | omics | Detect, quality-control, and classify viral contigs. Use when identifying viruses in assemblies, checking viral completeness and contamination, or assigning vi… |
| [`virulencefinder`](#virulencefinder) | pku-metagenomics | VirulenceFinder — CGE virulence gene detection in bacterial assemblies/reads for supported species (historically E. coli, Enterococcus, S. aureus, Listeria, …)… |

## `antismash` {#antismash}

Detect and annotate biosynthetic gene clusters (BGCs) on genomes and MAGs. Use for natural-product genome mining, choosing genefinding, interpreting region HTML/JSON, and stating validation limits. Upstream: https://github.com/antismash/antismash. Route via microbial-mining.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/antismash)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info antismash
dnakit skills install antismash --agent claude-code   # 或 --agent codex / cursor
```

## `anvio` {#anvio}

Anvi'o — interactive multi-omics platform for metagenomic assemblies, MAG refinement, visualization, and SNP/profile exploration. Use when: (1) Manual or human-guided bin refinement, (2) Interactive contig profiles across samples, (3) Linking MAGs with coverage/SNVs/metatranscriptomes, (4) Producing Anvi'o interactive HTML reports. Upstream: https://github.com/merenlab/anvio. Complements automated ensemble binning (basalt/metawrap/dastool); not a replacement for CheckM2/GUNC QC.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/anvio)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info anvio
dnakit skills install anvio --agent claude-code   # 或 --agent codex / cursor
```

## `basalt` {#basalt}

PKU-EMBL BASALT and BASALT-Air — recommended ensemble MAG recovery across one or more assemblies (multi-binner + refinement). Use this skill when: (1) Running BASALT multi-binner recovery, (2) Choosing Conda BASALT vs BASALT-Air, (3) Multi-assembly refinement, (4) Resuming checkpoints, (5) Writing reproducible BASALT Methods. Upstream: https://github.com/PKU-EMBL/BASALT and https://github.com/PKU-EMBL/BASALT-Air. Ensemble peers: dastool (SCG aggregation), metawrap (binning+refinement suite). Prefer basalt as the default ensemble route; use metawrap/dastool when a classical suite or contig2bin consensus is enough. For study design of assemblies, use study-design. For QC thresholds, use checkm2 and mag-qc.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/basalt)
- **许可证：**`MIT`
- **别名：**`dastool`, `metabinner`, `metawrap`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info basalt
dnakit skills install basalt --agent claude-code   # 或 --agent codex / cursor
```

## `bigscape` {#bigscape}

Cluster and compare biosynthetic gene clusters across genomes with BiG-SCAPE (similarity networks of antiSMASH regions). Use after batch antiSMASH for BGC family discovery. Upstream: https://github.com/medema-group/BiG-SCAPE.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/bigscape)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info bigscape
dnakit skills install bigscape --agent claude-code   # 或 --agent codex / cursor
```

## `binchicken` {#binchicken}

Bin Chicken — targeted coassembly planning for recovering low-abundance / novel MAGs by grouping metagenomes that share marker-gene signal (SingleM). Use when: (1) Public or large cohorts where all-vs-all coassembly is impossible, (2) Prioritizing sample subsets for coassembly to recover rare lineages, (3) Novelty-focused recovery against reference genomes, (4) Study-design of coassembly strata. Upstream: https://github.com/AroneyS/binchicken. For general coassembly vs individual design use study-design. For binning after coassembly use basalt / metabat2 / etc.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/binchicken)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info binchicken
dnakit skills install binchicken --agent claude-code   # 或 --agent codex / cursor
```

## `binning-qc` {#binning-qc}

Bin and refine metagenomic contigs, then assess MAG quality. Use when recovering genomes with QuickBin and checking completeness, contamination, and bin consistency.

- **来源：**[Copyright (c) 2026 Frederik Schulz @ `7a4eb5b`](https://github.com/fmschulz/omics-skills/tree/7a4eb5bace3f874c3be2ec697e9d2e695c374637/skills/bio-binning-qc)
- **许可证：**`MIT`
- **别名：**`bio-binning-qc`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info binning-qc
dnakit skills install binning-qc --agent claude-code   # 或 --agent codex / cursor
```

## `biosurfdb` {#biosurfdb}

BioSurfDB — curated portal for biosurfactant production and hydrocarbon biodegradation genes/pathways with BLAST and related metagenome analysis helpers. Use for surfactant/biodegradation-focused searches and pathway browsing. Portal: https://www.biosurfdb.org/. Route via dom-biosurfactant. Offline experimental gene sets → hadeg; glycan DOM → dbcan.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/biosurfdb)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info biosurfdb
dnakit skills install biosurfdb --agent claude-code   # 或 --agent codex / cursor
```

## `cami` {#cami}

CAMI / CAMI II metagenome benchmarking datasets and evaluation stack (CAMISIM, AMBER, OPAL, CAMI Benchmarking Portal). Use when: (1) Downloading CAMI challenge data / gold standards, (2) Choosing marine, plant-associated, strain-madness, mouse-gut, or human toy sets, (3) Evaluating binners with AMBER or profilers with OPAL, (4) Simulating communities with CAMISIM, (5) Citing CAMI results in Methods without over-generalizing to new biomes. Portal: https://cami-challenge.org/ and https://data.cami-challenge.org/. For real-world public reads use sra-ena; for biome MAG catalogues use mgnify.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/cami)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info cami
dnakit skills install cami --agent claude-code   # 或 --agent codex / cursor
```

## `checkm2` {#checkm2}

ML-based MAG completeness and contamination estimation. Use this skill when working with checkm2: installing, choosing parameters, running commands, interpreting outputs, or troubleshooting. Upstream: https://github.com/chklovski/CheckM2. For routing across tools, use tool-selection. For end-to-end pipelines, use metagenomics-workflow.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/checkm2)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info checkm2
dnakit skills install checkm2 --agent claude-code   # 或 --agent codex / cursor
```

## `checkv` {#checkv}

Quality assessment of metagenome-assembled viral genomes. Use this skill when working with checkv: installing, choosing parameters, running commands, interpreting outputs, or troubleshooting. Upstream: https://bitbucket.org/berkeleylab/checkv. For routing across mining tasks use microbial-mining. For stage order use metagenomics-workflow.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/checkv)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info checkv
dnakit skills install checkv --agent claude-code   # 或 --agent codex / cursor
```

## `dbcan` {#dbcan}

Annotate carbohydrate-active enzymes (CAZymes) with dbCAN3 / run_dbcan on genomes or protein sets, including CGC finding and glycan substrate prediction (dbCAN-sub / dbCAN-PUL) relevant to dissolved organic matter (DOM) polymer processing. Use for CAZyme mining complementary to DRAM/eggNOG and to dom-biosurfactant routing. Upstream: https://github.com/linnabrown/run_dbcan.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/dbcan)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info dbcan
dnakit skills install dbcan --agent claude-code   # 或 --agent codex / cursor
```

## `deepurify` {#deepurify}

Deepurify — multi-modal deep language model for MAG decontamination by matching contigs to taxonomic lineages and splitting chimeric bins into higher-quality sub-MAGs. Use when: (1) Cleaning contaminated MAGs after binning, (2) Rescuing near-threshold contaminated bins, (3) Complementing CheckM2/GUNC (detection vs contig removal). Upstream: https://github.com/ericcombiolab/Deepurify.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/deepurify)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info deepurify
dnakit skills install deepurify --agent claude-code   # 或 --agent codex / cursor
```

## `defensefinder` {#defensefinder}

Detect anti-phage defense systems in prokaryotic genomes. Use this skill when working with defensefinder: installing, choosing parameters, running commands, interpreting outputs, or troubleshooting. Upstream: https://github.com/mdmparis/defense-finder. For routing across mining tasks use microbial-mining. For stage order use metagenomics-workflow.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/defensefinder)
- **许可证：**`MIT`
- **别名：**`padloc`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info defensefinder
dnakit skills install defensefinder --agent claude-code   # 或 --agent codex / cursor
```

## `disinfinder` {#disinfinder}

DisinFinder — CGE disinfectant / biocide resistance gene screen via the DisinFinder database, typically run through ResFinder. Use for isolate biocide-resistance markers alongside acquired AMR. DB: https://bitbucket.org/genomicepidemiology/disinfinder_db. Route via cge-finders / resfinder / microbial-mining.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/disinfinder)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info disinfinder
dnakit skills install disinfinder --agent claude-code   # 或 --agent codex / cursor
```

## `dom-biosurfactant` {#dom-biosurfactant}

Route dissolved organic matter (DOM) / glycan substrate capacity, biosurfactant production, and hydrocarbon/plastic biodegradation annotation for shotgun metagenomics. Covers “like-dissolves-like” solubility-related microbial functions: CAZyme–glycan (dbCAN), biosurfactants + HC degradation (HADEG, BioSurfDB), and related BGC context (antiSMASH). Child skills: hadeg, biosurfdb; also dbcan, antismash, dram. Route via microbial-mining.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/dom-biosurfactant)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info dom-biosurfactant
dnakit skills install dom-biosurfactant --agent claude-code   # 或 --agent codex / cursor
```

## `dram` {#dram}

Distill microbial metabolism annotations for MAGs/genomes into curated product tables (DRAM annotate + distill). Use for ecological functional potential, pathway summaries, and MAG metabolism figures. Upstream: https://github.com/WrightonLabCSU/DRAM.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/dram)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info dram
dnakit skills install dram --agent claude-code   # 或 --agent codex / cursor
```

## `ecological-genomics-biodiversity-metrics` {#ecological-genomics-biodiversity-metrics}

Quantifies biodiversity from species abundance/incidence tables using Hill numbers (iNEXT) with coverage-based rarefaction-extrapolation (Chao &amp; Jost 2012), asymptotic richness via Chao1/ACE/jackknife as a lower bound, Baselga turnover/nestedness partition with the Podani alternative as sensitivity check, mandatory Hellinger transformation before ordination (Legendre &amp; Gallagher 2001), Faith PD and SES_MPD/SES_MNTD with explicit null-model choice, and Maire 2015 functional-diversity dimensionality optimization. Use when comparing diversity across sites with unequal sampling effort, picking the right richness estimator for singleton-heavy amplicon data, partitioning beta diversity into turnover vs nestedness, reporting Hill-number effective species counts rather than raw entropies, computing SES_MPD with explicit null-model justification, or deciding whether to apply standard metrics to compositional amplicon data. Not for clinical 16S microbiome diversity (see microbiome/diversity-analysis).

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/ecological-genomics/biodiversity-metrics)
- **许可证：**`MIT`
- **别名：**`bio-ecological-genomics-biodiversity-metrics`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info ecological-genomics-biodiversity-metrics
dnakit skills install ecological-genomics-biodiversity-metrics --agent claude-code   # 或 --agent codex / cursor
```

## `ecological-genomics-community-ecology` {#ecological-genomics-community-ecology}

Analyzes species-environment relationships with constrained ordination (CCA, RDA, db-RDA), variance partitioning, indicator species (indicspecies IndVal.g group-equalized), PERMANOVA paired MANDATORILY with PERMDISP (Anderson &amp; Walsh 2013; dispersion confounds centroid tests), Joint Species Distribution Models (HMSC, sjSDM, gjam) with explicit rejection of "residual covariance equals biotic interaction", phylogenetic community ecology (SES_MPD/MNTD), trait-environment via RLQ + fourth-corner with corrected modeltype=6 (Dray 2014), bipartite network metrics (NODF, modularity) with curveball null (Strona 2014), and Mantel-test replacements (dbRDA, GDM) for spatial data. Use when testing how environmental gradients structure communities, identifying habitat indicator taxa, partitioning variance among predictors, deciding whether PERMANOVA significance is location vs dispersion, picking among HMSC/sjSDM/gjam, or replacing Mantel tests for landscape data.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/ecological-genomics/community-ecology)
- **许可证：**`MIT`
- **别名：**`bio-ecological-genomics-community-ecology`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info ecological-genomics-community-ecology
dnakit skills install ecological-genomics-community-ecology --agent claude-code   # 或 --agent codex / cursor
```

## `ecological-genomics-conservation-genetics` {#ecological-genomics-conservation-genetics}

Assesses genetic health of populations for conservation with Ne estimation across time horizons (LDNe NeEstimator V2 option-file API + SNeP physical-linkage correction; recent trajectory via GONE/GONE2; deep history via Stairway Plot 2 / dadi / fastsimcoal2 / PSMC), F-statistics, runs of homozygosity binned by length class to date inbreeding, genetic-load decomposition (Bertorelle 2022 realized vs masked), the modern 100/1000 Ne rule (Frankham 2014), Ne/Nc 2-6 orders of magnitude in marine fish (Hauser &amp; Carvalho 2008), tree-sequence forward simulations (SLiM 4 + pyslim + tskit), and the Sukumaran-Knowles caveat against MSC methods for management-unit definition. Use when estimating Ne by time horizon, detecting inbreeding via F_ROH, decomposing genetic load, justifying conservation thresholds, distinguishing ESU/MU/DPS, configuring NeEstimator V2, or correcting LDNe physical linkage.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/ecological-genomics/conservation-genetics)
- **许可证：**`MIT`
- **别名：**`bio-ecological-genomics-conservation-genetics`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info ecological-genomics-conservation-genetics
dnakit skills install ecological-genomics-conservation-genetics --agent claude-code   # 或 --agent codex / cursor
```

## `ecological-genomics-edna-metabarcoding` {#ecological-genomics-edna-metabarcoding}

Processes eDNA metabarcoding from raw paired-end reads to species tables, navigating ASV (DADA2, UNOISE3) vs OTU (swarm v2) decision (Callahan 2017 vs Schloss multi-copy-16S critique), marker/primer choice (Leray COI, MiFish 12S, 515F/806R 16S, ITS2) with primer-specific bias, OBITools3 v3 command-name break (obi stats plural; .tar.gz taxonomy), tag-jumping with dual-indexing (Schnell 2015; NovaSeq 10x MiSeq), decontam as screening-not-classifier (Davis 2018), read-counts-not-abundance critique (Lamb 2019), site-occupancy modeling (Ficetola 2015), Naive-Bayes calibration limits (Bokulich 2018), and eDNA decay (Strickler 2015). Use when going from raw eDNA FASTQ to species tables, picking marker + denoising pipeline, deciding whether read counts represent abundance, applying occupancy modeling, configuring OBITools3 v3, or interpreting decontam output. Not for clinical 16S microbiome (see microbiome/amplicon-processing).

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/ecological-genomics/edna-metabarcoding)
- **许可证：**`MIT`
- **别名：**`bio-ecological-genomics-edna-metabarcoding`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info ecological-genomics-edna-metabarcoding
dnakit skills install ecological-genomics-edna-metabarcoding --agent claude-code   # 或 --agent codex / cursor
```

## `ecological-genomics-landscape-genomics` {#ecological-genomics-landscape-genomics}

Tests genotype-environment associations and identifies adaptive loci while correcting for the four-confound landscape (structure, demography, background selection, sampling design) using LFMM2 with mandatory K via sNMF cross-entropy elbow (LEA 3), BayPass Core/AUX/C2/IS with Omega covariance matrix, RDA / pRDA for polygenic adaptation (Forester 2018; requires imputed genotypes), OutFLANK with trimmed FST null, pcadapt, gradient forests (Ellis-Smith-Pitcher 2012, NOT mis-cited Ellis-Manel), Capblancq &amp; Forester 2021 RDA Swiss-army-knife, genomic-offset prediction with Lind &amp; Lotterhos 2025 three-regime caveat, Lotterhos-Whitlock sampling optima, Wang &amp; Bradburd 2014 IBD vs IBE, and Circuitscape + ResistanceGA. Use when identifying adaptive loci across gradients, choosing K for LFMM2, deciding among GEA methods, predicting maladaptation with the novel-environment caveat, distinguishing IBD vs IBE, or optimizing sampling design.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/ecological-genomics/landscape-genomics)
- **许可证：**`MIT`
- **别名：**`bio-ecological-genomics-landscape-genomics`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info ecological-genomics-landscape-genomics
dnakit skills install ecological-genomics-landscape-genomics --agent claude-code   # 或 --agent codex / cursor
```

## `ecological-genomics-species-delimitation` {#ecological-genomics-species-delimitation}

Delimits putative species boundaries from molecular data within the de Queiroz 2007 unified-lineage framework using ASAP (Puillandre 2021 successor to ABGD), mPTP C++ (Kapli 2017 successor to bPTP; bPTP is Python NOT R), GMYC single/multi-threshold (Pons 2006; Fujisawa 2013), multilocus BPP v4 with prior calibration from data (NOT defaults; Yang 2015), SNAPP + BFD* for SNP delimitation, DELINEATE (Sukumaran 2021) speciation-process modeling to address Sukumaran &amp; Knowles 2017 PNAS critique that MSC delimits structure not species, integrative-taxonomy congruence (Padial 2010; Carstens 2013), Dsuite for introgression testing before sister claims (Malinsky 2021), and Meyer &amp; Paulay 2005 barcoding-gap-absence caveat. Use when delineating species from DNA barcoding data, resolving cryptic complexes, choosing among ASAP/mPTP/BPP/DELINEATE, calibrating BPP priors, distinguishing introgression from ILS, or applying the Sukumaran-Knowles oversplitting correction.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/ecological-genomics/species-delimitation)
- **许可证：**`MIT`
- **别名：**`bio-ecological-genomics-species-delimitation`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info ecological-genomics-species-delimitation
dnakit skills install ecological-genomics-species-delimitation --agent claude-code   # 或 --agent codex / cursor
```

## `element-cycling` {#element-cycling}

Route biogeochemical / elemental cycling annotation for shotgun metagenomics and MAGs: carbon (C), nitrogen (N), sulfur (S), iron (Fe), phosphorus (P), hydrogen (H2/H), oxygen (O), methane (CH4), and other metals/metalloids (Mn, As, Se, …) via METABOLIC, MEBS, FeGenie, plus DRAM/KEGG/eggNOG backups. Use when asking which elemental transformations a community or MAG can perform — not for taxonomy tables. Child skills: metabolic, mebs, fegenie. Related: dram, kegg, eggnog-mapper, gapseq. Route via microbial-mining.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/element-cycling)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info element-cycling
dnakit skills install element-cycling --agent claude-code   # 或 --agent codex / cursor
```

## `epidemiological-genomics-amr-surveillance` {#epidemiological-genomics-amr-surveillance}

Detects acquired antimicrobial-resistance determinants and chromosomal point-mutation resistance in bacterial assemblies using AMRFinderPlus, ResFinder 4.0 (acquired + PointFinder), CARD-RGI, abritAMR, staramr, and species-specific callers (TB-Profiler, Mykrobe). Harmonises cross-tool output via hAMRonization, contextualises determinants with mobile-genetic-element annotation (MOB-suite, PlasmidFinder, MobileElementFinder, ICEberg), predicts phenotype against EUCAST or CLSI breakpoints, and translates calls into WHO GLASS reporting categories. Use when screening clinical or surveillance isolates for AMR, distinguishing acquired vs intrinsic vs point-mutation resistance, calling rpoB / katG / pncA / gyrA / mgrB mutations, reconciling AMRFinderPlus vs RGI vs ResFinder disagreement, contextualising carbapenemases or mcr alleles on plasmids, predicting susceptibility from genotype against the WHO Mtb 2nd-edition catalogue, or building a hAMRonized multi-lab AMR surveillance pipeline.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/epidemiological-genomics/amr-surveillance)
- **许可证：**`MIT`
- **别名：**`abricate`, `bio-epidemiological-genomics-amr-surveillance`, `pointfinder`, `resfinder`, `staramr`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info epidemiological-genomics-amr-surveillance
dnakit skills install epidemiological-genomics-amr-surveillance --agent claude-code   # 或 --agent codex / cursor
```

## `epidemiological-genomics-pathogen-typing` {#epidemiological-genomics-pathogen-typing}

Assigns isolate identity at the right resolution for the question -- ANI/Mash species triage, 7-locus MLST historical comparability, cgMLST/wgMLST outbreak resolution (chewBBACA, BIGSdb, Ridom SeqSphere, EnteroBase HierCC), in-silico serotyping (SISTR/SeqSero2 Salmonella, SerotypeFinder E. coli, Kaptive Klebsiella, SeroBA pneumococcus, spa+SCCmec S. aureus), and lineage callers (TB-Profiler/Mykrobe barcode for MTBC, Pangolin + Nextclade for SARS-CoV-2, PopPUNK GPSC for S. pneumoniae). Use when typing bacterial isolates for surveillance or outbreak investigation, choosing between cgMLST allele distance and core-SNP distance for cluster definition, harmonising calls across schemas/database versions, assigning MTBC lineage with the Napier 90-SNP barcode, calling Salmonella serovar via SISTR with monophasic Typhimurium awareness, running Pangolin UShER mode with pangolin-data version pinning, or selecting a typing resolution to match the surveillance question.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/epidemiological-genomics/pathogen-typing)
- **许可证：**`MIT`
- **别名：**`bio-epidemiological-genomics-pathogen-typing`, `mlst`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info epidemiological-genomics-pathogen-typing
dnakit skills install epidemiological-genomics-pathogen-typing --agent claude-code   # 或 --agent codex / cursor
```

## `epidemiological-genomics-phylodynamics` {#epidemiological-genomics-phylodynamics}

Estimates time-scaled phylogenies, molecular-clock rates, effective reproduction number R_e, and population dynamics from dated pathogen genomes using TreeTime (maximum-likelihood) and BEAST2 (Bayesian; strict/relaxed clocks; coalescent, Bayesian-Skyline, Skygrid, Birth-Death-Skyline, and sampled-ancestor priors; structured coalescent via MASCOT). Covers root-to-tip clock QC via TempEst, date-randomisation tests, recombination masking via Gubbins/ClonalFrameML before clock inference for recombining bacteria, BDSKY origin-vs-rootHeight pitfalls, sampling-bias correction, multi-chain convergence diagnostics, and reconciling phylodynamic R_e with case-based R_t. Use when dating outbreak origins, estimating substitution rates, inferring R_e through time, building time-calibrated Nextstrain Augur trees, choosing between strict and relaxed clocks, fitting Birth-Death-Skyline models, diagnosing temporal-signal failure, running MASCOT for structured-population analyses, or using UShER for pandemic-scale placement.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/epidemiological-genomics/phylodynamics)
- **许可证：**`MIT`
- **别名：**`bio-epidemiological-genomics-phylodynamics`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info epidemiological-genomics-phylodynamics
dnakit skills install epidemiological-genomics-phylodynamics --agent claude-code   # 或 --agent codex / cursor
```

## `epidemiological-genomics-transmission-inference` {#epidemiological-genomics-transmission-inference}

Infers person-to-person transmission from pathogen genomes using outbreaker2, TransPhylo, phybreak, BadTrIP, SCOTTI, BEASTLIER, and SNP-distance / cluster-picker approaches (HIV-TRACE for HIV; transcluster). Defines outbreak clusters using pathogen-specific SNP thresholds (NOT a universal cutoff -- TB &lt;=12 SNPs; MRSA &lt;=15; C. difficile &lt;=2; Klebsiella &lt;=21), models within-host diversity and transmission bottlenecks, integrates contact-tracing data, distinguishes generation from serial interval, and attributes source via Bayesian source attribution (islandR). Use when investigating outbreaks for who-infected-whom, defining SNP-cluster outbreak definitions, accounting for unsampled intermediates, choosing between outbreaker2 (rich epi data) and TransPhylo (genomic-only after a dated phylogeny), running source attribution between host populations, calling HIV-TRACE thresholds appropriate to the local subtype, or distinguishing recent transmission from reactivation in TB or chronic HIV.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/epidemiological-genomics/transmission-inference)
- **许可证：**`MIT`
- **别名：**`bio-epidemiological-genomics-transmission-inference`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info epidemiological-genomics-transmission-inference
dnakit skills install epidemiological-genomics-transmission-inference --agent claude-code   # 或 --agent codex / cursor
```

## `epidemiological-genomics-variant-surveillance` {#epidemiological-genomics-variant-surveillance}

Assigns pathogen lineages (SARS-CoV-2 Pangolin UShER mode; Nextclade clade + QC; pango-designation alias resolution) and tracks variant frequencies over time using Nextstrain (Augur + Auspice), wastewater deconvolution (Freyja, COJAC, alcov, lineagespot), lineage-fitness modelling (multinomial logistic), and recombinant detection (3SEQ, RDP4, Bolotie). Covers Pangolin pangolin-data and Nextclade dataset version pinning (mandatory; lineage-defining mutations change with dataset), Freyja barcode forward-only date constraint, ARTIC primer scheme churn (V3/V4/V4.1/V5.3.2/Midnight) with dropout regions, and recombinant X-prefix designation lag. Use when assigning Pango lineages and Nextclade clades to viral consensus sequences, building Nextstrain Augur surveillance pipelines, deconvolving wastewater into lineage frequencies with Freyja, tracking lineage frequencies over time, handling ARTIC primer dropouts, or running surveillance for SARS-CoV-2/influenza/Mpox/RSV/H5N1/measles.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/epidemiological-genomics/variant-surveillance)
- **许可证：**`MIT`
- **别名：**`bio-epidemiological-genomics-variant-surveillance`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info epidemiological-genomics-variant-surveillance
dnakit skills install epidemiological-genomics-variant-surveillance --agent claude-code   # 或 --agent codex / cursor
```

## `fegenie` {#fegenie}

FeGenie — HMM-based identification of iron-related genes and neighborhoods in genomes/metagenome assemblies (acquisition, siderophores, Fe redox, storage, magnetosome). Use for Fe cycling potential on contigs or MAGs. Upstream: https://github.com/Arkadiy-Garber/FeGenie. Route via element-cycling / microbial-mining. Broader C/N/S → metabolic; S scores → mebs.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/fegenie)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info fegenie
dnakit skills install fegenie --agent claude-code   # 或 --agent codex / cursor
```

## `genomad` {#genomad}

Identify plasmids, viruses, and other mobile genetic elements in assemblies or genomes. Use for plasmidome/virome mining, MGE-aware AMR context, and pre-CheckV viral candidate lists. Upstream: https://github.com/apcamargo/genomad.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/genomad)
- **许可证：**`MIT`
- **别名：**`virsorter2`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info genomad
dnakit skills install genomad --agent claude-code   # 或 --agent codex / cursor
```

## `gunc` {#gunc}

Detect chimerism/contamination in prokaryotic genomes/MAGs. Use this skill when working with gunc: installing, choosing parameters, running commands, interpreting outputs, or troubleshooting. Upstream: https://github.com/grp-bork/gunc. For stage routing use tool-selection / metagenomics-workflow.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/gunc)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info gunc
dnakit skills install gunc --agent claude-code   # 或 --agent codex / cursor
```

## `hadeg` {#hadeg}

HADEG — curated Hydrocarbon Aerobic Degradation Enzymes and Genes database (alkanes, alkenes, aromatics, plastics) plus biosurfactant production genes with experimental evidence. Use to annotate MAG/contig proteins for HC biodegradation and biosurfactant potential. Upstream: https://github.com/jarojasva/HADEG. Route via dom-biosurfactant / microbial-mining. Portal-style BioSurfDB → biosurfdb; glycan DOM → dbcan.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/hadeg)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info hadeg
dnakit skills install hadeg --agent claude-code   # 或 --agent codex / cursor
```

## `insilicoseq` {#insilicoseq}

InSilicoSeq — simulate realistic Illumina (meta)genomic reads for benchmarking and experimental design. Use when: (1) Generating Illumina-like test sets from genomes, (2) Simple abundance-controlled mock communities, (3) Complementing CAMISIM for CAMI-style community simulation. Upstream: https://github.com/HadrienG/InSilicoSeq. For full CAMI gold standards prefer cami / CAMISIM.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/insilicoseq)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info insilicoseq
dnakit skills install insilicoseq --agent claude-code   # 或 --agent codex / cursor
```

## `integronfinder` {#integronfinder}

IntegronFinder — detect integrons (intI, attC, cassettes) in bacterial genomes/assemblies, including AMR cassette context. Use for integron-focused MGE mining. Upstream: https://github.com/gem-pasteur/Integron_Finder. Route via microbial-mining; pair with rgi / mobileog / genomad.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/integronfinder)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info integronfinder
dnakit skills install integronfinder --agent claude-code   # 或 --agent codex / cursor
```

## `isescan` {#isescan}

ISEScan — automated annotation of insertion sequence (IS) elements in prokaryotic genomes via profile HMMs. Use for IS family catalogues on MAGs or contigs as part of MGE refinement beyond geNomad. Upstream: https://github.com/xiezhq/ISEScan. Route via microbial-mining; pair with mobileog / integronfinder.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/isescan)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info isescan
dnakit skills install isescan --agent claude-code   # 或 --agent codex / cursor
```

## `lorbin` {#lorbin}

LorBin — unsupervised long-read metagenome binning with multiscale adaptive clustering (DBSCAN + BIRCH) and single-copy-gene evaluation. Use when: (1) Binning long-read / HiFi metagenome assemblies, (2) Recovering MAGs from species-rich or abundance-imbalanced communities, (3) Comparing LorBin to SemiBin2/VAMB long-read modes, (4) Feeding LorBin bins into ensemble refinement (BASALT `-e l` needs the PKU-EMBL LorBin fork). Upstream: https://github.com/LorMeBioAI/LorBin. For ensemble recovery use basalt. For short-read defaults use metabat2 / semibin2 / vamb.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/lorbin)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info lorbin
dnakit skills install lorbin --agent claude-code   # 或 --agent codex / cursor
```

## `mebs` {#mebs}

MEBS (Multigenomic Entropy-Based Score) — score genomes/metagenomes for biogeochemical cycle machinery via informative Pfam domains; strongest for sulfur (S), also supports N, O, CH4, and Fe cycle scores. Use for comparative cycle enrichment across large (meta)genome collections. Upstream: https://github.com/valdeanda/mebs (also eead-csic-compbio/metagenome_Pfam_score). Route via element-cycling. Gene-level Fe → fegenie; multi-element diagrams → metabolic.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/mebs)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info mebs
dnakit skills install mebs --agent claude-code   # 或 --agent codex / cursor
```

## `metabat2` {#metabat2}

Adaptive contig binning from tetranucleotide composition and multi-sample coverage. Use for primary single-binner runs feeding DAS Tool, MetaWRAP, or BASALT ensembles. Upstream: https://bitbucket.org/berkeleylab/metabat.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/metabat2)
- **许可证：**`MIT`
- **别名：**`aamb`, `comebin`, `concoct`, `maxbin2`, `semibin2`, `vamb`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info metabat2
dnakit skills install metabat2 --agent claude-code   # 或 --agent codex / cursor
```

## `metabolic` {#metabolic}

METABOLIC (AnantharamanLab) — high-throughput metabolic and biogeochemical trait profiling of MAGs/SAGs/isolates with C, N, S, and other elemental cycle diagrams (METABOLIC-G genomes-only; METABOLIC-C with reads / community networks). Use for multi-element cycling potential and MW-score networks. Upstream: https://github.com/AnantharamanLab/METABOLIC. Route via element-cycling. Fe-deep dive → fegenie; S entropy scores → mebs.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/metabolic)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info metabolic
dnakit skills install metabolic --agent claude-code   # 或 --agent codex / cursor
```

## `metagenomics-abundance` {#metagenomics-abundance}

Turns shotgun classifier output into a defensible abundance table with Bracken Bayesian re-estimation, then compositional treatment (CLR, zero handling), library-size normalization, reference-frame differential abundance, and optional absolute quantification. Covers why a relative-abundance change is not a change, why Bracken read fractions and MetaPhlAn percentages are different physical quantities, the silent -r read-length bias, the genome-size confound no library-size method fixes, and the rarefaction debate. Use when estimating species abundance from a Kraken2 report, normalizing a community count table, choosing a compositional transform, or converting relative to absolute load. For classification see kraken-classification; for diversity/ordination/DA mechanics see metagenome-visualization.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/metagenomics/abundance-estimation)
- **许可证：**`MIT`
- **别名：**`bio-metagenomics-abundance`, `coverm`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info metagenomics-abundance
dnakit skills install metagenomics-abundance --agent claude-code   # 或 --agent codex / cursor
```

## `metagenomics-amr-detection` {#metagenomics-amr-detection}

Profiles the antimicrobial-resistance gene content (resistome) of shotgun metagenomes - read-based quantification with RGI bwt, AMR++/MEGARes, ARGs-OAP/SARG, deepARG, or GROOT, and presence calling with AMRFinderPlus/ABRicate on assembled contigs or MAGs. Covers why an ARG hit is a sequence match not a phenotype, why a metagenomic ARG has no host and no genomic context until assembly (and assembly breaks at ARGs), per-gene curated thresholds vs a flat 80/80, gene-fraction false-positive control, and cross-study normalization pitfalls. Use when quantifying a community resistome, normalizing ARG abundance, or calling ARGs from metagenome contigs. For pure-culture isolate AMR, point mutations, and phenotype/MIC prediction see epidemiological-genomics/amr-surveillance.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/metagenomics/amr-detection)
- **许可证：**`MIT`
- **运行依赖：**命令 `AMRFinderPlus`
- **别名：**`bio-metagenomics-amr-detection`, `card-amr-profiler`, `deeparg`, `rgi`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info metagenomics-amr-detection
dnakit skills install metagenomics-amr-detection --agent claude-code   # 或 --agent codex / cursor
```

## `metagenomics-contamination-controls` {#metagenomics-contamination-controls}

Cleans a shotgun metagenome of everything that is not the target community before profiling - host-read depletion (Hostile, bowtie2/T2T-CHM13), reagent/kitome contamination control with blanks and decontam, mock-community validation, and depth-adequacy checks (Nonpareil). Covers why a metagenomic result is a position in a choice-chain rather than a direct observation, why extraction is the experiment, why a low-biomass community can be entirely kitome, why absence means not-detectable-by-this-chain, and why a confident classifier call can still be wrong when the reference is contaminated. Use when designing controls, removing host reads, identifying reagent contaminants, validating with mocks, or judging whether a low-biomass result is real. For adapter/quality trimming see read-qc; for MAG-level decontamination see genome-assembly/metagenome-assembly.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/metagenomics/contamination-controls)
- **许可证：**`MIT`
- **别名：**`bio-metagenomics-contamination-controls`, `kneaddata`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info metagenomics-contamination-controls
dnakit skills install metagenomics-contamination-controls --agent claude-code   # 或 --agent codex / cursor
```

## `metagenomics-functional-profiling` {#metagenomics-functional-profiling}

Profiles the functional potential of shotgun metagenomes with HUMAnN 3's tiered search (MetaPhlAn prescreen, Bowtie2 pangenome, translated DIAMOND vs UniRef), giving gene-family (RPK) and MetaCyc pathway abundances stratified by species. Covers why a metagenome measures potential not activity, why dropping UNMAPPED/UNINTEGRATED biases everything, why stratification is an estimate, coverage-vs-abundance and MinPath/gap-fill, UniRef90-vs-50 and biome database bias, and the assembly/eggNOG/dbCAN/antiSMASH alternatives. Use when obtaining pathway or gene-family abundances, regrouping to KO/EC/GO, normalizing functional tables, or choosing read-based vs assembly-based functional profiling. For AMR genes see amr-detection; for host-gene enrichment see pathway-analysis.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/metagenomics/functional-profiling)
- **许可证：**`MIT`
- **运行依赖：**命令 `HUMAnN`
- **别名：**`bio-metagenomics-functional-profiling`, `humann`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info metagenomics-functional-profiling
dnakit skills install metagenomics-functional-profiling --agent claude-code   # 或 --agent codex / cursor
```

## `metagenomics-kraken` {#metagenomics-kraken}

Classifies shotgun metagenomic reads to taxa with Kraken2's minimizer/LCA matching against a chosen reference database, then hands off to Bracken for abundance re-estimation. Covers why the database (not the algorithm) decides what can be detected, the --confidence and --minimum-hit-groups precision levers, unique-minimizer false-positive control, host-read removal, and why raw Kraken2 read counts are not abundances. Use when profiling who-is-there from shotgun reads, choosing a Kraken2 database, setting a confidence threshold, controlling false positives, or feeding reports to Bracken. For marker-gene profiling see metaphlan-profiling; for abundance mechanics see abundance-estimation; for assembly/MAG recovery see genome-assembly/metagenome-assembly.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/metagenomics/kraken-classification)
- **许可证：**`MIT`
- **运行依赖：**命令 `Kraken2`
- **别名：**`bio-metagenomics-kraken`, `bracken`, `kaiju`, `kraken2`, `metabuli`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info metagenomics-kraken
dnakit skills install metagenomics-kraken --agent claude-code   # 或 --agent codex / cursor
```

## `metagenomics-metaphlan` {#metagenomics-metaphlan}

Profiles shotgun metagenomes to species/SGB relative abundance with MetaPhlAn 4's clade-specific marker genes (bowtie2 short reads, minimap2 long reads). Covers why a MetaPhlAn percentage is a cell fraction (genome-size-normalized taxonomic abundance) and must never be merged with Kraken/Bracken read fractions, kSGB vs uSGB units for quantifying database-absent taxa, the unknown-fraction rescaling and its version-default flip, --index pinning as a batch variable, and when mOTUs3 or sourmash gather beat marker profiling. Use when profiling who-is-there with high precision, needing HMP-comparable species abundances, quantifying novel taxa, or deciding marker-gene vs k-mer profiling. For k-mer classification see kraken-classification; for strains see strain-tracking; for 16S amplicon see the microbiome category.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/metagenomics/metaphlan-profiling)
- **许可证：**`MIT`
- **运行依赖：**命令 `MetaPhlAn`
- **别名：**`bio-metagenomics-metaphlan`, `metaphlan`, `motus`, `sylph`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info metagenomics-metaphlan
dnakit skills install metagenomics-metaphlan --agent claude-code   # 或 --agent codex / cursor
```

## `metagenomics-strain-tracking` {#metagenomics-strain-tracking}

Resolves and compares bacterial strains below the species level from shotgun metagenomes with inStrain (popANI/conANI microdiversity), StrainPhlAn (marker-SNV consensus phylogeny and nGD), MIDAS2, metaSNV, and StrainGE, plus genome-vs-genome ANI (skani/fastANI/MASH) for isolate/MAG comparison. Covers why a strain is a threshold not a thing, why ANI answers same-genome while popANI/nGD answer same-population-in-situ, the 99.999% popANI and per-species nGD definitions, the coverage detection limit (absence is not absence), why sharing is not transmission direction, and mapping to the dataset's own dRep MAGs. Use when detecting shared strains, tracking transmission, resolving within-host strain dynamics, or deconvoluting co-occurring strains. For pure-culture isolate outbreak SNP trees see epidemiological-genomics; for MAG assembly see genome-assembly/metagenome-assembly.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/metagenomics/strain-tracking)
- **许可证：**`MIT`
- **别名：**`bio-metagenomics-strain-tracking`, `instrain`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info metagenomics-strain-tracking
dnakit skills install metagenomics-strain-tracking --agent claude-code   # 或 --agent codex / cursor
```

## `metagenomics-visualization` {#metagenomics-visualization}

Turns a shotgun profiler table (MetaPhlAn relative abundance, Bracken counts, HUMAnN function tables) into honest figures and defensible community statistics with phyloseq, vegan, microViz, and Python. Covers why an ordination/bar/diversity number is a modeling choice that can manufacture a result, the MetaPhlAn-percent-vs-Bracken-counts fork that decides everything, CLR/Aitchison vs Bray-Curtis, Hill numbers and why shotgun richness is a database readout, pairing PERMANOVA with betadisper, and the multi-tool differential-abundance consensus. Use when plotting taxonomic/functional profiles, computing alpha/beta diversity, running ordination/PERMANOVA, or testing differential abundance. For amplicon/QIIME2 stats see the microbiome category; for compositional theory see abundance-estimation.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/metagenomics/metagenome-visualization)
- **许可证：**`MIT`
- **别名：**`bio-metagenomics-visualization`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info metagenomics-visualization
dnakit skills install metagenomics-visualization --agent claude-code   # 或 --agent codex / cursor
```

## `microbial-genome-characterization` {#microbial-genome-characterization}

Genome-ASSEMBLY discovery, QC, and replicon mapping for any organism (bacteria, archaea, fungi, and beyond) using NCBI Datasets. Resolves an organism name or taxid to assemblies, picks the reference/representative or best-quality assembly, pulls assembly QC metrics (total length, contig/scaffold N50, contig count, GC%, assembly level, RefSeq category), enumerates chromosomes and plasmids via per-replicon sequence reports, and compares candidate assemblies on quality. Use for "what genomes are available for \[organism\]", "assembly stats / N50 / GC content for \[GCF_/GCA_ accession\]", "how many plasmids does \[strain\] have", "compare assemblies for \[species\]", "find the reference genome for \[taxon\]", "is this assembly Complete Genome or just contigs". NOT for gene-level orthology/synteny (use comparative-genomics), plant gene structure (use plant-genomics), de novo assembly from raw reads (no tool exists), or taxonomy-only name/lineage lookups.

- **来源：**[ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-microbial-genome-characterization)
- **许可证：**`Apache-2.0`
- **别名：**`tooluniverse-microbial-genome-characterization`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info microbial-genome-characterization
dnakit skills install microbial-genome-characterization --agent claude-code   # 或 --agent codex / cursor
```

## `microbial-mining` {#microbial-mining}

Route microbial / genome mining after assembly or binning: BGCs (antiSMASH, BiG-SCAPE), AMR (RGI, DeepARG, ResFinder), mobile elements (geNomad, VirSorter2, CheckV, mobileOG-db, ISEScan, IntegronFinder, MOB-suite, CGE Finders), antimicrobial peptides (Macrel), element cycling (METABOLIC, MEBS, FeGenie), DOM/biosurfactant/HC biodegradation (dbCAN, HADEG, BioSurfDB), metabolism (DRAM, gapseq), defense (DefenseFinder, PADLOC), quorum sensing (QSAP), and annotation substrates (Bakta/Prokka). Use when the question is discovery of functional elements, not community tables alone.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/microbial-mining)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info microbial-mining
dnakit skills install microbial-mining --agent claude-code   # 或 --agent codex / cursor
```

## `microbiome-amplicon-processing` {#microbiome-amplicon-processing}

Infers exact amplicon sequence variants (ASVs) from demultiplexed 16S rRNA or ITS amplicon FASTQ with DADA2 - removing primers with cutadapt (--discard-untrimmed), learning a per-run error model (filterAndTrim -&gt; learnErrors -&gt; dada -&gt; mergePairs), merging run-level tables with mergeSequenceTables, then one removeBimeraDenovo. Covers why primers come OFF before truncation, why the error model is per-run, truncLen as a merge-overlap detection budget (V4 vs V3-V4), DADA2 vs Deblur and q2-dada2 (denoise-paired/single/pyro/ccs), ASV vs OTU, NovaSeq binned-quality error-fit breakage, ITSxpress for variable-length ITS, and decontam removal of reagent/kit contaminants. Use when turning demultiplexed amplicon reads into an ASV/feature table, choosing truncation lengths, handling multi-run studies, or ITS. For shotgun reads see metagenomics/kraken-classification; for QIIME2 CLI mechanics see qiime2-workflow; for primer trimming theory see read-qc/adapter-trimming.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/microbiome/amplicon-processing)
- **许可证：**`MIT`
- **别名：**`bio-microbiome-amplicon-processing`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info microbiome-amplicon-processing
dnakit skills install microbiome-amplicon-processing --agent claude-code   # 或 --agent codex / cursor
```

## `microbiome-differential-abundance` {#microbiome-differential-abundance}

Tests which individual taxa differ between groups on an amplicon ASV/feature table (phyloseq) using compositionally-aware methods - ALDEx2 (Dirichlet-MC CLR, conservative), ANCOM-BC2/ANCOMBC (sampling-fraction bias correction, structural zeros, passed_ss, default p_adj_method=holm), MaAsLin2/MaAsLin3 (multivariable GLM, random effects, prevalence/abundance split), LinDA (CLR mixed-model regression), ZicoSeq (permutation FDR), LEfSe, and q2-composition ancombc. Covers why the hit list depends more on the DA tool than the biology (Nearing benchmark) so the deliverable is a CONSENSUS of &gt;=2 tools, why a relative change is not absolute without a load anchor, the prevalence-filter knob, BH/FDR plus an effect-size floor, and why DESeq2/edgeR misfire here. Use when finding differentially abundant taxa, handling covariates or longitudinal designs, or choosing a method. Whole-community diversity -&gt; diversity-analysis; shotgun DA -&gt; metagenomics/metagenome-visualization; CoDA theory -&gt; metagenomics/abundance-estimation

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/microbiome/differential-abundance)
- **许可证：**`MIT`
- **别名：**`ancombc`, `bio-microbiome-differential-abundance`, `lefse`, `maaslin2`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info microbiome-differential-abundance
dnakit skills install microbiome-differential-abundance --agent claude-code   # 或 --agent codex / cursor
```

## `microbiome-diversity-analysis` {#microbiome-diversity-analysis}

Alpha and beta diversity of an amplicon (16S/ITS) ASV/OTU community table - observed features, Shannon, Pielou evenness, Faith PD, Bray-Curtis, Jaccard, weighted/unweighted/generalized UniFrac, Aitchison/RPCA - via QIIME2 core-metrics-phylogenetic, phyloseq/vegan, and scikit-bio. Covers the three knobs that set the answer before it is seen (rarefaction sampling depth, the tree, the metric), why core-metrics silently deletes samples below the sampling depth, why de novo trees lose to SEPP fragment-insertion and Greengenes2, why unweighted and weighted UniFrac can flip the story, why observed features is an ASV count not a species count, the QIIME2-log2 vs R-ln Shannon mismatch, and pairing PERMANOVA (adonis2) with betadisper. Use when summarizing whole-community richness/evenness or testing group differences in community structure. Per-taxon testing -&gt; differential-abundance. Shotgun tables -&gt; metagenomics/metagenome-visualization. Shared CoDA/rarefaction theory -&gt; metagenomics/abundance-estimation.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/microbiome/diversity-analysis)
- **许可证：**`MIT`
- **别名：**`bio-microbiome-diversity-analysis`, `microeco`, `phyloseq`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info microbiome-diversity-analysis
dnakit skills install microbiome-diversity-analysis --agent claude-code   # 或 --agent codex / cursor
```

## `microbiome-functional-prediction` {#microbiome-functional-prediction}

Predicts community functional POTENTIAL from 16S/ITS amplicon ASVs with PICRUSt2 (or q2-picrust2) by phylogenetic interpolation of reference-genome gene content - EPA-ng placement, gappa, castor hidden-state prediction of KO/EC/Pfam copy number, 16S copy-number normalization, and MinPath MetaCyc/KEGG pathways - gated by the NSTI quality index. Covers why predicted function is taxonomy re-encoded (never measured gene content and never activity), the mandatory NSTI report (--max_nsti 2 silently drops novel ASVs), why accuracy IS reference coverage (gut Spearman ~0.8, soil/marine collapse), the circularity trap, and Tax4Fun2/FAPROTAX/BugBase alternatives. Use when inferring KO/EC/MetaCyc potential from an ASV table, gating on NSTI, or choosing a prediction method. For MEASURED shotgun function see metagenomics/functional-profiling; for enrichment of KO lists see pathway-analysis/go-enrichment; for DA of predicted tables see differential-abundance.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/microbiome/functional-prediction)
- **许可证：**`MIT`
- **运行依赖：**命令 `PICRUSt2`
- **别名：**`bio-microbiome-functional-prediction`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info microbiome-functional-prediction
dnakit skills install microbiome-functional-prediction --agent claude-code   # 或 --agent codex / cursor
```

## `microbiome-qiime2-workflow` {#microbiome-qiime2-workflow}

Operates the QIIME2 framework as the glue for an amplicon analysis - the .qza/.qzv artifact model, semantic types (FeatureTable\[Frequency\], SampleData\[PairedEndSequencesWithQuality\], Phylogeny\[Rooted\], FeatureData\[Taxonomy\]), embedded provenance plus provenance replay, import (Casava/manifest/EMP/BIOM), export, the Metadata object, and the q2cli vs Artifact API interfaces. Covers why a .qza is data-plus-executable-history not a file, why export drops provenance, why a .qzv is terminal, why classifier .qza are version-pinned, and the 2026 distribution/rachis rename. Use when importing reads, choosing a manifest/Casava/EMP/BIOM path, reading or replaying provenance, exporting to BIOM/phyloseq, fixing semantic-type or Phred or sklearn-version errors, or orchestrating the pipeline. Denoising -&gt; amplicon-processing; classifier/DB -&gt; taxonomy-assignment; diversity metric/depth -&gt; diversity-analysis; DA tool -&gt; differential-abundance; PICRUSt2 -&gt; functional-prediction; shotgun moshpit -&gt; metagenomics.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/microbiome/qiime2-workflow)
- **许可证：**`MIT`
- **运行依赖：**命令 `QIIME2`
- **别名：**`alterlab-qiime2-amplicon`, `bio-microbiome-qiime2-workflow`, `qiime2-16s-microbiome`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info microbiome-qiime2-workflow
dnakit skills install microbiome-qiime2-workflow --agent claude-code   # 或 --agent codex / cursor
```

## `microbiome-research` {#microbiome-research}

Microbiome research using MGnify, GTDB, ENA, OLS (ENVO biomes), and EuropePMC. Covers study discovery, taxonomic profiling, host-microbe interaction analysis, and biome-by-condition queries. Use for microbiome study selection, organism-environment associations, and clinical-microbiome literature review. Distinct from analytical workflow (use tooluniverse-metagenomics-analysis for that).

- **来源：**[ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-microbiome-research)
- **许可证：**`Apache-2.0`
- **运行依赖：**Python 模块 `dnakit`, `pandas`
- **别名：**`tooluniverse-microbiome-research`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info microbiome-research
dnakit skills install microbiome-research --agent claude-code   # 或 --agent codex / cursor
```

## `microbiome-taxonomy-assignment` {#microbiome-taxonomy-assignment}

Assigns taxonomy to amplicon ASVs/OTUs (16S, ITS, 18S) with a classifier conditioned on a reference database and primer region - DADA2 assignTaxonomy + addSpecies (RDP naive Bayes), DECIPHER IDTAXA, and QIIME2 q2-feature-classifier (classify-sklearn naive Bayes, classify-consensus-vsearch alignment-consensus). Covers region-specific training (extract-reads, fit-classifier-naive-bayes), why a full-length classifier fabricates calls on a V4 read, the scikit-learn version-pinning trap on pre-trained .qza classifiers, confidence thresholds (classify-sklearn 0.7, assignTaxonomy minBoot 50), and choosing SILVA/GTDB/Greengenes2/UNITE/PR2/RDP. Use when classifying ASVs after DADA2, picking a reference database, training a region-matched classifier, setting a confidence threshold, or deciding whether a 16S species call is defensible (usually not - genus at best). For shotgun read classification see metagenomics/kraken-classification and metagenomics/metaphlan-profiling.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/microbiome/taxonomy-assignment)
- **许可证：**`MIT`
- **别名：**`bio-microbiome-taxonomy-assignment`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info microbiome-taxonomy-assignment
dnakit skills install microbiome-taxonomy-assignment --agent claude-code   # 或 --agent codex / cursor
```

## `mob-suite` {#mob-suite}

MOB-suite — plasmid reconstruction and typing from draft assemblies (replicon, relaxase/MOB type, mate-pair, transferability, mash clusters). Use for plasmid compatibility / MOB typing and plasmid contig binning beyond geNomad scores. Upstream: https://github.com/phac-nml/mob-suite. Route via microbial-mining; pair with genomad, plasmaag, mobileog, rgi.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/mob-suite)
- **许可证：**`MIT`
- **别名：**`plasmidfinder`, `pmlst`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info mob-suite
dnakit skills install mob-suite --agent claude-code   # 或 --agent codex / cursor
```

## `mobileelementfinder` {#mobileelementfinder}

MobileElementFinder — CGE detection of named mobile genetic elements (MITEs, IS, Tn, ComTn, ICE, IME, CIME, …) in assembled bacterial genomes, with optional ARG/virulence context. Use for isolate mobilome screens. Upstream: PyPI MobileElementFinder / CGE web service. Route via cge-finders / microbial-mining; complement isescan, mobileog, integronfinder, genomad.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/mobileelementfinder)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info mobileelementfinder
dnakit skills install mobileelementfinder --agent claude-code   # 或 --agent codex / cursor
```

## `mobileog` {#mobileog}

mobileOG-db — manually curated protein families for bacterial MGE life-cycle functions (integration/excision, replication/recombination/repair, stability/defense, transfer, phage-related) with diamond + mobileOG-pl parsing. Use to annotate MGE machinery on MAG/contig proteins beyond geNomad contig scores. Upstream: https://github.com/clb21565/mobileOG-db. Portal: https://mobileogdb.flsi.cloud.vt.edu/. Route via microbial-mining.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/mobileog)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info mobileog
dnakit skills install mobileog --agent claude-code   # 或 --agent codex / cursor
```

## `pathogen-variant-surveillance` {#pathogen-variant-surveillance}

Query live pathogen genomic surveillance data through the GenSpectrum LAPIS API to find which viral lineages are circulating now, how fast they are growing, and what mutations they carry. Use whenever a question depends on the current state of a pathogen population rather than on remembered facts - which SARS-CoV-2 variant is dominant, whether a Pango lineage is still designated or has been withdrawn, what clade or genotype of H5N1 is in a host or region, whether a PCR primer or assay target still matches circulating sequence, or how a lineage's prevalence has moved week to week. Triggers include "variant surveillance", "genomic surveillance", "what variant is circulating", "dominant variant", "Pango lineage", "lineage prevalence", "growth advantage", "SARS-CoV-2 variant", "XFG", "clade 2.3.4.4b", "H5N1 genotype", "influenza clade", "RSV/mpox/measles/dengue lineage", "CoV-Spectrum", "LAPIS", "Nextclade", "pango-designation", and any request to report what a pathogen population looks like today.

- **来源：**[K-Dense Inc.; K-Dense Inc. and contributors @ `c1ed16d`](https://github.com/K-Dense-AI/scientific-agent-skills/tree/c1ed16d97dd61ff50a3bd46dd353e4a55fd77f34/skills/pathogen-variant-surveillance)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info pathogen-variant-surveillance
dnakit skills install pathogen-variant-surveillance --agent claude-code   # 或 --agent codex / cursor
```

## `plasmaag` {#plasmaag}

PlasMAAG — recover plasmids and cellular genomes from metagenomes using assembly–alignment graphs and contrastive learning across samples. Use when: (1) Plasmid MAG reconstruction beyond circular-path assemblers, (2) Improving plasmid contig classification vs geNomad alone, (3) Host–plasmid association across multi-sample datasets. Upstream: https://github.com/RasmussenLab/PlasMAAG. Pair with genomad for MGE screening and checkv for viral (not plasmid) QC.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/plasmaag)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info plasmaag
dnakit skills install plasmaag --agent claude-code   # 或 --agent codex / cursor
```

## `qsap` {#qsap}

QSAP — quorum-sensing (QS) related gene annotation pipeline against the QSP database (DIAMOND and/or HMMER). Use when mining autoinducer / QS circuitry (AHL, AIP, AI-2-related families as covered by QSP) from proteins or assemblies. Upstream: https://github.com/chunxiao-dcx/QSAP (QSP DB: https://github.com/chunxiao-dcx/QSP). Route via microbial-mining. Related secondary metabolites → antismash; Fe siderophores → fegenie.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/qsap)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info qsap
dnakit skills install qsap --agent claude-code   # 或 --agent codex / cursor
```

## `singlem` {#singlem}

SingleM — species-level metagenome profiling from conserved windows in universal marker genes, including lineages without genomic representation; Sandpiper aggregates public SingleM profiles. Use when: (1) Profiling environmental communities rich in unknowns, (2) Feeding Bin Chicken coassembly targeting, (3) GTDB-linked marker OTUs from short reads, (4) Comparing to MetaPhlAn/mOTUs/sylph. Upstream: https://github.com/wwood/singlem. Sandpiper: https://sandpiper.qut.edu.au.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/singlem)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info singlem
dnakit skills install singlem --agent claude-code   # 或 --agent codex / cursor
```

## `taxvamb` {#taxvamb}

TaxVAMB — semisupervised bimodal VAE binning that integrates TNF, coabundance, and contig taxonomy for improved MAG recovery (incl. incomplete genomes and long-read sets). Use when: (1) Taxonomy-informed binning, (2) Comparing to VAMB/AAMB/SemiBin2, (3) Single-sample setups needing more HQ bins. Upstream: https://github.com/RasmussenLab/vamb. Taxonomy labels are an input feature — garbage taxonomy in, garbage bins out.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/taxvamb)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info taxvamb
dnakit skills install taxvamb --agent claude-code   # 或 --agent codex / cursor
```

## `vcontact` {#vcontact}

vConTACT2 / vConTACT3 — gene-sharing network and ML-based hierarchical virus taxonomy for prokaryotic (and in v3, broader realm) viral genomes/contigs. Use after viral discovery (geNomad / VirSorter2) and quality (CheckV): (1) Genus-to-order style clustering, (2) Placing vMAGs vs RefSeq/ICTV, (3) Choosing v2 vs v3. Upstream: https://bitbucket.org/MAVERICLab/vcontact2 and https://bitbucket.org/MAVERICLab/vcontact3.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/vcontact)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info vcontact
dnakit skills install vcontact --agent claude-code   # 或 --agent codex / cursor
```

## `viromics` {#viromics}

Detect, quality-control, and classify viral contigs. Use when identifying viruses in assemblies, checking viral completeness and contamination, or assigning viral taxonomy.

- **来源：**[Copyright (c) 2026 Frederik Schulz @ `7a4eb5b`](https://github.com/fmschulz/omics-skills/tree/7a4eb5bace3f874c3be2ec697e9d2e695c374637/skills/bio-viromics)
- **许可证：**`MIT`
- **别名：**`bio-viromics`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info viromics
dnakit skills install viromics --agent claude-code   # 或 --agent codex / cursor
```

## `virulencefinder` {#virulencefinder}

VirulenceFinder — CGE virulence gene detection in bacterial assemblies/reads for supported species (historically E. coli, Enterococcus, S. aureus, Listeria, …). Upstream: https://bitbucket.org/genomicepidemiology/virulencefinder / PyPI Virulencefinder. Route via cge-finders / microbial-mining; not a general metagenome pathogenicity score.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/virulencefinder)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info virulencefinder
dnakit skills install virulencefinder --agent claude-code   # 或 --agent codex / cursor
```
