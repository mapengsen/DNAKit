# 比较基因组与系统发育

本类共 **30** 个 Skills。每个 Skill 是一份让 Agent 按步骤完成分析的说明，分析所需的软件、参考数据和 API 需要另外准备。搜索、安装和使用方法见 [Skills](../skills.md)；网页版的 Skills 标签可以直接阅读每个 Skill 的完整内容。

| Skill | 来源 | 简介 |
| --- | --- | --- |
| [`comparative-genomics`](#comparative-genomics) | tooluniverse | Cross-species gene comparison and ortholog analysis. Integrates Ensembl Compara orthologs, NCBI Gene, UniProt, OLS, Monarch, and OpenTargets to identify orthol… |
| [`comparative-genomics-ancestral-reconstruction`](#comparative-genomics-ancestral-reconstruction) | bioskills | Reconstruct ancestral states at internal phylogenetic nodes for sequences (PAML codeml, IQ-TREE --ancestral, GRASP, FastML), discrete traits (corHMM hidden-rat… |
| [`comparative-genomics-comparative-annotation-projection`](#comparative-genomics-comparative-annotation-projection) | bioskills | Project gene annotations across genomes using TOGA (Kirilenko 2023 whole-genome-alignment chain-based projection with intactness classification), CESAR 2.0 (Sh… |
| [`comparative-genomics-gene-family-evolution`](#comparative-genomics-gene-family-evolution) | bioskills | Model gene-family birth-death dynamics across a species tree using CAFE5 (Mendes et al 2020 Bioinformatics 36:5516 gamma-distributed rate categories), CAFE5-er… |
| [`comparative-genomics-gene-tree-species-tree-reconciliation`](#comparative-genomics-gene-tree-species-tree-reconciliation) | bioskills | Reconcile gene trees against a species tree under probabilistic models of duplication, transfer, and loss (DTL) using ALE (Szöllősi 2013 amalgamated likelihood… |
| [`comparative-genomics-genome-distance-and-species-delineation`](#comparative-genomics-genome-distance-and-species-delineation) | bioskills | Compute genome-to-genome distances (ANI, AAI, dDDH, k-mer Mash) and assign taxonomic classifications using skani (Shaw 2023), FastANI (Jain 2018), pyani / pyAN… |
| [`comparative-genomics-hgt-detection`](#comparative-genomics-hgt-detection) | bioskills | Detect horizontal gene transfer (HGT / LGT) using compositional methods (GC%, codon usage, tetranucleotide z-scores via SIGI-HMM, AlienHunter, IslandViewer 4,… |
| [`comparative-genomics-introgression-detection`](#comparative-genomics-introgression-detection) | bioskills | Detect introgression and admixture between species or populations using Dsuite (Malinsky 2021 fast D-statistics), Patterson's D / ABBA-BABA test (Green 2010; D… |
| [`comparative-genomics-ortholog-inference`](#comparative-genomics-ortholog-inference) | bioskills | Infer orthologous genes and gene families across species using OrthoFinder3 (HOG-based phylogenetic orthology), SonicParanoid2, Broccoli, ProteinOrtho, OMA / F… |
| [`comparative-genomics-pangenome-analysis`](#comparative-genomics-pangenome-analysis) | bioskills | Build and analyze pangenomes for prokaryotes (Panaroo, PPanGGOLiN, PEPPAN, GET_HOMOLOGUES, anvi'o pangenomics) and eukaryotes (Minigraph-Cactus, PGGB, vg pange… |
| [`comparative-genomics-positive-selection`](#comparative-genomics-positive-selection) | bioskills | Detect positive (diversifying / episodic / pervasive) selection using codon dN/dS frameworks. Implements PAML codeml site models (M0/M1a/M2a/M7/M8/M8a), branch… |
| [`comparative-genomics-synteny-analysis`](#comparative-genomics-synteny-analysis) | bioskills | Detect syntenic blocks and structural rearrangements between genomes using MCScanX (Wang 2012), JCVI/MCScan (Tang 2008 Python), GENESPACE (Lovell 2022) for ort… |
| [`comparative-genomics-whole-genome-alignment`](#comparative-genomics-whole-genome-alignment) | bioskills | Build whole-genome alignments using Progressive Cactus (Armstrong 2020 reference-free clade-level WGA), Minigraph-Cactus (Hickey 2024 pangenome-aware), LASTZ c… |
| [`comparative-genomics-whole-genome-duplication`](#comparative-genomics-whole-genome-duplication) | bioskills | Detect, date, and contextualize whole-genome duplication (WGD / paleopolyploidy) events using wgd v2 (Chen et al 2024), KsRates (Sensalari 2022 substitution-ra… |
| [`data-visualization-circos-plots`](#data-visualization-circos-plots) | bioskills | Build circular genome visualizations using circlize (R), pyCirclize (Python), or Circos (Perl CLI) with ideogram tracks, multi-data tracks (scatter, histogram,… |
| [`interdomain-hgt`](#interdomain-hgt) | omics | Detect and polarize interdomain horizontal gene transfer with homology, context, and phylogenetic checks. Use when studying lateral gene transfer, virus-host g… |
| [`phylo-bayesian-inference`](#phylo-bayesian-inference) | bioskills | Frames Bayesian phylogenetics as approximating a posterior distribution over trees conditioned on data AND priors via an MCMC that must be proven to have conve… |
| [`phylo-distance-calculations`](#phylo-distance-calculations) | bioskills | Build model-corrected evolutionary distance matrices and distance trees (NJ, BIONJ, FastME, UPGMA) with Biopython Bio.Phylo plus R ape/phangorn/FastME. Covers… |
| [`phylo-divergence-dating`](#phylo-divergence-dating) | bioskills | Estimate divergence times under molecular-clock models with BEAST2, MCMCTree/PAML, TreePL, and LSD2, framing a date as a product of the calibration prior and t… |
| [`phylo-modern-tree-inference`](#phylo-modern-tree-inference) | bioskills | Infers maximum-likelihood phylogenetic trees with IQ-TREE2 and RAxML-NG -- model selection (ModelFinder), branch support (UFBoot2, SH-aLRT), concordance factor… |
| [`phylo-species-trees`](#phylo-species-trees) | bioskills | Estimates species trees under the multispecies coalescent from per-locus gene trees with the modern ASTER astral binary (ASTRAL-III/wASTRAL/ASTRAL-Pro), plus S… |
| [`phylo-tree-io`](#phylo-tree-io) | bioskills | Read, write, and convert phylogenetic tree files with Biopython Bio.Phylo, and choose an annotation-preserving parser (treeio, DendroPy) when metadata matters.… |
| [`phylo-tree-manipulation`](#phylo-tree-manipulation) | bioskills | Edit phylogenetic tree structure with Biopython Bio.Phylo, and treat rooting as a separate statistical inference rather than a display choice. Covers why most… |
| [`phylo-tree-visualization`](#phylo-tree-visualization) | bioskills | Draw and export phylogenetic trees with Bio.Phylo plus matplotlib, and route rich figures to ggtree, ETE4, or iTOL. Covers why a tree figure is an argument not… |
| [`phylogenomics`](#phylogenomics) | omics | Build and validate marker-gene alignments and phylogenetic trees. Use when inferring evolutionary relationships, choosing models, or checking tree support and… |
| [`plant-genomics`](#plant-genomics) | tooluniverse | Plant genomics and biology research — PlantReactome pathways, Ensembl Plants gene structure, POWO species taxonomy, UniProt annotation, KEGG plant pathways. Ha… |
| [`spacedust`](#spacedust) | pku-metagenomics | Spacedust — de novo discovery of conserved gene clusters in microbial genomes using Foldseek structure comparison plus clustering/order-conservation statistics… |
| [`tracking-taxonomy-updates`](#tracking-taxonomy-updates) | omics | Track taxonomy changes and triage sequence assignments across NCBI, GTDB, ICTV, and eukaryotic frameworks. Use when comparing releases, resolving renamed taxa,… |
| [`unicore`](#unicore) | pku-metagenomics | Unicore — universal structure-based core-gene phylogeny using Foldseek and ProstT5 (Steinegger lab). Use when: (1) Building phylogenies from proteomes via stru… |
| [`vitis-pangenome`](#vitis-pangenome) | plant-skills | Vitis / crop pangenome decisions: when to build a graph, PAV, reference choice for mixed cohorts and chips. Use for minigraph-cactus vs PGGB, PanGenie, or "sho… |

## `comparative-genomics` {#comparative-genomics}

Cross-species gene comparison and ortholog analysis. Integrates Ensembl Compara orthologs, NCBI Gene, UniProt, OLS, Monarch, and OpenTargets to identify orthologs, paralogs, sequence conservation, functional conservation across species, and lineage-specific gene gains/losses. Use for phylogenetic gene tracing, model-organism mapping, and evolutionary-genomics queries.

- **来源：**[ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-comparative-genomics)
- **许可证：**`Apache-2.0`
- **别名：**`tooluniverse-comparative-genomics`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info comparative-genomics
dnakit skills install comparative-genomics --agent claude-code   # 或 --agent codex / cursor
```

## `comparative-genomics-ancestral-reconstruction` {#comparative-genomics-ancestral-reconstruction}

Reconstruct ancestral states at internal phylogenetic nodes for sequences (PAML codeml, IQ-TREE --ancestral, GRASP, FastML), discrete traits (corHMM hidden-rate Markov, ape::ace, phytools::make.simmap stochastic mapping, BayesTraits), and continuous traits (phytools::fastAnc, geiger Brownian/OU, RPANDA). Use when designing constructs for ancestral protein resurrection, tracing trait evolution along a tree, performing stochastic character mapping, testing models of trait evolution (BM vs OU vs EB), inferring ancestral genome content via Dollo or DTL reconciliation, or quantifying ancestral-state uncertainty for downstream comparative analyses.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/comparative-genomics/ancestral-reconstruction)
- **许可证：**`MIT`
- **别名：**`bio-comparative-genomics-ancestral-reconstruction`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info comparative-genomics-ancestral-reconstruction
dnakit skills install comparative-genomics-ancestral-reconstruction --agent claude-code   # 或 --agent codex / cursor
```

## `comparative-genomics-comparative-annotation-projection` {#comparative-genomics-comparative-annotation-projection}

Project gene annotations across genomes using TOGA (Kirilenko 2023 whole-genome-alignment chain-based projection with intactness classification), CESAR 2.0 (Sharma, Schwede &amp; Hiller 2017 codon-aware exon projection), LiftOff (Shumate &amp; Salzberg 2021 reference-based annotation transfer), Liftover (UCSC), GeMoMa (Keilwagen 2019 evidence-based projection), and Comparative Annotation Toolkit (CAT). Use when transferring annotations from a well-annotated reference to query genome(s), classifying gene-loss vs gene-intact across many genomes at scale, building Zoonomia-style comparative annotations across hundreds of mammals or birds (Kirilenko 2023), detecting pseudogenization, projecting alternative isoforms, or selecting between WGA-anchored (TOGA) vs ortholog-based (LiftOff) annotation transfer strategies.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/comparative-genomics/comparative-annotation-projection)
- **许可证：**`MIT`
- **运行依赖：**命令 `TOGA`
- **别名：**`bio-comparative-genomics-comparative-annotation-projection`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info comparative-genomics-comparative-annotation-projection
dnakit skills install comparative-genomics-comparative-annotation-projection --agent claude-code   # 或 --agent codex / cursor
```

## `comparative-genomics-gene-family-evolution` {#comparative-genomics-gene-family-evolution}

Model gene-family birth-death dynamics across a species tree using CAFE5 (Mendes et al 2020 Bioinformatics 36:5516 gamma-distributed rate categories), CAFE5-error (annotation-error-aware), Count (Csurös 2010 ancestral state reconstruction), BadiRate (Librado 2012 likelihood + parsimony), DupliPHY-Family, and ALE/AleRax (for per-family DTL; see \[\[gene-tree-species-tree-reconciliation\]\]). Test lineage-specific gene-family expansions and contractions, distinguish biological dynamics from annotation artifacts, account for assembly fragmentation, identify functional enrichment in expanded / contracted families. Use when correlating gene-family changes with phenotype evolution, ranking lineages by adaptive gene-family-rate shifts, post-WGD dosage-balance analysis, or building Birth-death models from OrthoFinder presence/absence matrices.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/comparative-genomics/gene-family-evolution)
- **许可证：**`MIT`
- **运行依赖：**命令 `CAFE5`
- **别名：**`bio-comparative-genomics-gene-family-evolution`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info comparative-genomics-gene-family-evolution
dnakit skills install comparative-genomics-gene-family-evolution --agent claude-code   # 或 --agent codex / cursor
```

## `comparative-genomics-gene-tree-species-tree-reconciliation` {#comparative-genomics-gene-tree-species-tree-reconciliation}

Reconcile gene trees against a species tree under probabilistic models of duplication, transfer, and loss (DTL) using ALE (Szöllősi 2013 amalgamated likelihood), GeneRax (Morel 2020 ML reconciliation), AleRax (Morel 2024 co-estimation), Whale.jl (Bayesian DL+WGD), RANGER-DTL 2 parsimony, NOTUNG, ecceTERA, and Treerecs. Use when inferring ancestral gene-family content, distinguishing duplication from horizontal transfer from differential loss, rooting deep species trees from gene-content signals (STRIDE / Williams 2017 ALE-rooting), counting DTL events per branch, refining noisy gene trees against a species tree, modeling WGD events jointly with DTL, or producing publication-grade gene-family histories for phylogenomic / comparative analyses.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/comparative-genomics/gene-tree-species-tree-reconciliation)
- **许可证：**`MIT`
- **运行依赖：**命令 `ALE`
- **别名：**`bio-comparative-genomics-gene-tree-species-tree-reconciliation`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info comparative-genomics-gene-tree-species-tree-reconciliation
dnakit skills install comparative-genomics-gene-tree-species-tree-reconciliation --agent claude-code   # 或 --agent codex / cursor
```

## `comparative-genomics-genome-distance-and-species-delineation` {#comparative-genomics-genome-distance-and-species-delineation}

Compute genome-to-genome distances (ANI, AAI, dDDH, k-mer Mash) and assign taxonomic classifications using skani (Shaw 2023), FastANI (Jain 2018), pyani / pyANI ANIb / ANIm, OrthoANI (Lee 2016), AAI (amino-acid identity), dDDH via TYGS / GGDC, GTDB-Tk (Chaumeil 2020 standard prokaryote taxonomy), and Mash MinHash (Ondov 2016). Use when delineating prokaryote species (95% ANI threshold; Jain 2018 Nat Commun 9:5114), assigning genomes to GTDB taxonomy with ANI radius, computing genome similarity matrices for clustering, classifying archaea, evaluating MAG (metagenome-assembled genome) species assignment, applying skani for fast metagenomic ANI screening, or reconciling 16S rRNA-based taxonomy with whole-genome ANI.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/comparative-genomics/genome-distance-and-species-delineation)
- **许可证：**`MIT`
- **运行依赖：**命令 `skani`
- **别名：**`bio-comparative-genomics-genome-distance-and-species-delineation`, `drep`, `fastani`, `gtdbtk`, `skani`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info comparative-genomics-genome-distance-and-species-delineation
dnakit skills install comparative-genomics-genome-distance-and-species-delineation --agent claude-code   # 或 --agent codex / cursor
```

## `comparative-genomics-hgt-detection` {#comparative-genomics-hgt-detection}

Detect horizontal gene transfer (HGT / LGT) using compositional methods (GC%, codon usage, tetranucleotide z-scores via SIGI-HMM, AlienHunter, IslandViewer 4, IslandPath-DIMOB), phylogenetic-incongruence methods (AvP, HGTphyloDetect, ALE / GeneRax / AleRax reconciliation, RANGER-DTL), and BLAST-distribution methods (HGTector v2, DarkHorse, Alien Index). Use when screening prokaryote genomes for genomic islands and HGT events, distinguishing HGT from incomplete lineage sorting / differential gene loss / hybridization, mapping donor lineages via phylogenetic placement, separating eukaryotic HGT from contamination, ruling out gBGC as a false signal, or quantifying DTL rates with ALE/GeneRax on bacterial trees.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/comparative-genomics/hgt-detection)
- **许可证：**`MIT`
- **别名：**`bio-comparative-genomics-hgt-detection`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info comparative-genomics-hgt-detection
dnakit skills install comparative-genomics-hgt-detection --agent claude-code   # 或 --agent codex / cursor
```

## `comparative-genomics-introgression-detection` {#comparative-genomics-introgression-detection}

Detect introgression and admixture between species or populations using Dsuite (Malinsky 2021 fast D-statistics), Patterson's D / ABBA-BABA test (Green 2010; Durand 2011), f4-ratio and f-branch statistic (Malinsky 2018), TreeMix (Pickrell &amp; Pritchard 2012), HyDe (Blischak 2018), QuIBL (Edelman 2019), sprime (Browning 2018), Twisst (Martin 2017), PhyloNet (Than 2008) for explicit phylogenetic networks, and qpAdm / qpGraph (Patterson 2012). Distinguish introgression from incomplete lineage sorting (ILS), ancestral structure, ghost-lineage admixture, and rate variation. Use when testing inter-species gene flow, dating admixture events, identifying introgressed segments, building phylogenetic networks for reticulate evolution, or applying the ABBAclustering (Koppetsch-Malinsky-Matschiner 2024) framework for divergent-species gene flow.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/comparative-genomics/introgression-detection)
- **许可证：**`MIT`
- **运行依赖：**命令 `Dsuite`
- **别名：**`archaic-introgression`, `bio-comparative-genomics-introgression-detection`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info comparative-genomics-introgression-detection
dnakit skills install comparative-genomics-introgression-detection --agent claude-code   # 或 --agent codex / cursor
```

## `comparative-genomics-ortholog-inference` {#comparative-genomics-ortholog-inference}

Infer orthologous genes and gene families across species using OrthoFinder3 (HOG-based phylogenetic orthology), SonicParanoid2, Broccoli, ProteinOrtho, OMA / FastOMA hierarchical orthologous groups, eggNOG-mapper, JustOrthologs, and TOGA whole-genome-alignment orthology. Use when building single-copy ortholog sets for phylogenomics, classifying co-orthologs and in/out-paralogs after gene duplication, propagating functional annotation via orthology with awareness of the ortholog conjecture, distinguishing speciation from duplication via gene-tree species-tree reconciliation, computing Quest-for-Orthologs benchmark performance, or running synteny-aware ortholog detection in WGD-affected lineages.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/comparative-genomics/ortholog-inference)
- **许可证：**`MIT`
- **别名：**`bio-comparative-genomics-ortholog-inference`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info comparative-genomics-ortholog-inference
dnakit skills install comparative-genomics-ortholog-inference --agent claude-code   # 或 --agent codex / cursor
```

## `comparative-genomics-pangenome-analysis` {#comparative-genomics-pangenome-analysis}

Build and analyze pangenomes for prokaryotes (Panaroo, PPanGGOLiN, PEPPAN, GET_HOMOLOGUES, anvi'o pangenomics) and eukaryotes (Minigraph-Cactus, PGGB, vg pangenome graphs). Implement Tettelin core/accessory/cloud genome decomposition (Tettelin 2005), Heap's law open/closed pangenome modeling, gene presence/absence GWAS (Scoary, pyseer), pangenome graph variant calling (vg, PanGenie), and structural-variation graph indexing. Use when assembling species- or genus-level pan-gene catalogs, separating core from accessory/shell/cloud genes, testing gene-content associations with phenotypes, building pangenome graphs from haplotype-resolved assemblies, calling SVs from pangenome graphs, or selecting between bacterial-pangenome and eukaryotic-pangenome workflows.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/comparative-genomics/pangenome-analysis)
- **许可证：**`MIT`
- **运行依赖：**命令 `Panaroo`
- **别名：**`bio-comparative-genomics-pangenome-analysis`, `bio-protein-clustering-pangenome`, `roary-pangenome`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info comparative-genomics-pangenome-analysis
dnakit skills install comparative-genomics-pangenome-analysis --agent claude-code   # 或 --agent codex / cursor
```

## `comparative-genomics-positive-selection` {#comparative-genomics-positive-selection}

Detect positive (diversifying / episodic / pervasive) selection using codon dN/dS frameworks. Implements PAML codeml site models (M0/M1a/M2a/M7/M8/M8a), branch models, branch-site model A (Zhang 2005), and HyPhy methods (BUSTED, BUSTED-S, BUSTED-MH, BUSTED-PH, MEME, FEL, FUBAR, aBSREL, SLAC, RELAX, GARD, FUBAR-MH). Includes McDonald-Kreitman framework (asymptotic alpha, impMKT, polyDFE, DFE-alpha, GRAPES) for within-species + divergence inference, RERconverge for trait-correlated rate shifts, CSUBST for convergent substitution, and PhyloAcc for accelerated noncoding evolution. Use when testing adaptive evolution at codons, branches, or full gene; running GARD recombination pre-screen; controlling alignment-error and gBGC false positives; reconciling PAML vs HyPhy results; or performing genome-scale selection scans.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/comparative-genomics/positive-selection)
- **许可证：**`MIT`
- **别名：**`bio-comparative-genomics-positive-selection`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info comparative-genomics-positive-selection
dnakit skills install comparative-genomics-positive-selection --agent claude-code   # 或 --agent codex / cursor
```

## `comparative-genomics-synteny-analysis` {#comparative-genomics-synteny-analysis}

Detect syntenic blocks and structural rearrangements between genomes using MCScanX (Wang 2012), JCVI/MCScan (Tang 2008 Python), GENESPACE (Lovell 2022) for orthology-anchored riparian visualization, SyRI for structural variation, AnchorWave for sequence-level synteny, i-ADHoRe 3.0 for highly diverged species, SynNet for synteny networks, and ntSynt for multi-genome macrosynteny. Use when identifying collinear gene blocks across species, distinguishing macrosynteny from microsynteny, detecting inversions/translocations/duplications, anchoring orthology in WGD lineages, producing publication riparian plots, computing synteny block age via Ks (cross-references whole-genome-duplication), or running synteny-aware ortholog inference in polyploids.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/comparative-genomics/synteny-analysis)
- **许可证：**`MIT`
- **别名：**`bio-comparative-genomics-synteny-analysis`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info comparative-genomics-synteny-analysis
dnakit skills install comparative-genomics-synteny-analysis --agent claude-code   # 或 --agent codex / cursor
```

## `comparative-genomics-whole-genome-alignment` {#comparative-genomics-whole-genome-alignment}

Build whole-genome alignments using Progressive Cactus (Armstrong 2020 reference-free clade-level WGA), Minigraph-Cactus (Hickey 2024 pangenome-aware), LASTZ chain/net (UCSC pipeline), MUMmer4 (Marçais 2018 pairwise), minimap2 -x asm5/10/20 (Li 2018 fast pairwise), AnchorWave (Song 2022 WGD-aware), and Mauve / progressiveMauve (bacterial). Operates the HAL toolkit (Hickey 2013) for downstream extraction including halSynteny, halLiftover, halBranchMutations, and hal2maf. Use when constructing multi-species alignments for comparative-annotation projection (TOGA), synteny detection, conservation analyses (phyloP / PhastCons), or pangenome graph construction; selecting between reference-free (Cactus) and reference-anchored (LASTZ chains/nets) approaches; tuning sensitivity for closely vs distantly related genomes; or producing HAL files for genome-wide downstream tools.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/comparative-genomics/whole-genome-alignment)
- **许可证：**`MIT`
- **运行依赖：**命令 `Cactus`
- **别名：**`bio-comparative-genomics-whole-genome-alignment`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info comparative-genomics-whole-genome-alignment
dnakit skills install comparative-genomics-whole-genome-alignment --agent claude-code   # 或 --agent codex / cursor
```

## `comparative-genomics-whole-genome-duplication` {#comparative-genomics-whole-genome-duplication}

Detect, date, and contextualize whole-genome duplication (WGD / paleopolyploidy) events using wgd v2 (Chen et al 2024), KsRates (Sensalari 2022 substitution-rate-corrected Ks dating), DupGen_finder (Qiao 2019), MAPS (Li 2018 phylogenomic), POInT (Conant 2008 ordered-block), SLEDGe (2024 ML-based), Whale.jl (Bayesian DL+WGD), and synteny-anchored paranome construction. Use when identifying ancient polyploidy from Ks distributions and synteny block analysis, positioning WGD events relative to speciation, distinguishing tandem from segmental from WGD duplications, dating the 2R/3R vertebrate / fish / salmonid WGDs, building paranome and Ks-age mixture models, applying KsRates substitution-rate correction across lineages, or testing alternative biased-fractionation / dosage-balance models post-WGD.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/comparative-genomics/whole-genome-duplication)
- **许可证：**`MIT`
- **别名：**`bio-comparative-genomics-whole-genome-duplication`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info comparative-genomics-whole-genome-duplication
dnakit skills install comparative-genomics-whole-genome-duplication --agent claude-code   # 或 --agent codex / cursor
```

## `data-visualization-circos-plots` {#data-visualization-circos-plots}

Build circular genome visualizations using circlize (R), pyCirclize (Python), or Circos (Perl CLI) with ideogram tracks, multi-data tracks (scatter, histogram, heatmap), chord/link arcs for interactions, and explicit circos.clear() between plots. Covers when circular is appropriate vs when Cartesian wins (Cleveland-McGill 1984), karyograms, and chromosome adjacency in chord diagrams. Use when adjacency on the circle conveys meaning — chromosome-level overview, structural variants, Hi-C interactions, cross-genome comparisons.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/data-visualization/circos-plots)
- **许可证：**`MIT`
- **别名：**`bio-data-visualization-circos-plots`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info data-visualization-circos-plots
dnakit skills install data-visualization-circos-plots --agent claude-code   # 或 --agent codex / cursor
```

## `interdomain-hgt` {#interdomain-hgt}

Detect and polarize interdomain horizontal gene transfer with homology, context, and phylogenetic checks. Use when studying lateral gene transfer, virus-host gene exchange, endogenous viral elements, or donor direction.

- **来源：**[Copyright (c) 2026 Frederik Schulz @ `7a4eb5b`](https://github.com/fmschulz/omics-skills/tree/7a4eb5bace3f874c3be2ec697e9d2e695c374637/skills/bio-interdomain-hgt)
- **许可证：**`MIT`
- **别名：**`bio-interdomain-hgt`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info interdomain-hgt
dnakit skills install interdomain-hgt --agent claude-code   # 或 --agent codex / cursor
```

## `phylo-bayesian-inference` {#phylo-bayesian-inference}

Frames Bayesian phylogenetics as approximating a posterior distribution over trees conditioned on data AND priors via an MCMC that must be proven to have converged, using MrBayes, BEAST2, RevBayes, and PhyloBayes-MPI. Covers why convergence (ESS, PSRF, ASDSF, topology vs scalar) is the load-bearing claim, why posterior probabilities are systematically higher than bootstrap and overconfident under model misspecification, why the default branch-length prior inflates tree length, why the harmonic-mean estimator must never select models (use stepping-stone), and when site-heterogeneous CAT-GTR is required at depth. Use when needing posterior clade support, model averaging, marginal-likelihood model comparison, or CAT models for deep phylogeny. Routes topology-only ML to modern-tree-inference, divergence times to divergence-dating, and tree summarization to tree-io.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/phylogenetics/bayesian-inference)
- **许可证：**`MIT`
- **别名：**`bio-phylo-bayesian-inference`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info phylo-bayesian-inference
dnakit skills install phylo-bayesian-inference --agent claude-code   # 或 --agent codex / cursor
```

## `phylo-distance-calculations` {#phylo-distance-calculations}

Build model-corrected evolutionary distance matrices and distance trees (NJ, BIONJ, FastME, UPGMA) with Biopython Bio.Phylo plus R ape/phangorn/FastME. Covers why a distance is a model-corrected estimate of substitutions per site that undercounts raw because of multiple/back/parallel hits (saturation); why the matrix discards the per-site information ML keeps; the LogDet/paralinear fix for compositional heterogeneity; the UPGMA molecular-clock trap; and the Bio.Phylo landmine that DistanceCalculator offers only identity/matrix distances, not JC/K80/TN93. Use when computing a distance matrix, building a fast NJ/FastME tree, seeding an ML search, barcoding, or testing substitution saturation before a deep tree. Routes ML and starting-tree work to modern-tree-inference, alignment quality to alignment/alignment-io, and tree I/O to tree-io.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/phylogenetics/distance-calculations)
- **许可证：**`MIT`
- **别名：**`bio-phylo-distance-calculations`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info phylo-distance-calculations
dnakit skills install phylo-distance-calculations --agent claude-code   # 或 --agent codex / cursor
```

## `phylo-divergence-dating` {#phylo-divergence-dating}

Estimate divergence times under molecular-clock models with BEAST2, MCMCTree/PAML, TreePL, and LSD2, framing a date as a product of the calibration prior and the clock model far more than of the sequence data. Covers why branch length = rate x time is nonidentifiable so only calibrations convert relative rate-time into absolute age; why the effective (marginal) prior on a calibrated node differs from the density specified, mandating a sample-from-prior run; the fossil-as-minimum rule, soft bounds, tip-dating, and the fossilized birth-death process; the temporal-signal check (TempEst root-to-tip regression + date-randomization) required before dating viruses or ancient DNA; and clock-model choice via the coefficient of variation. Use when dating nodes, calibrating with fossils or sampling dates, choosing a clock or dating engine, or routing topology to modern-tree-inference, posteriors to bayesian-inference, and rooting to tree-manipulation.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/phylogenetics/divergence-dating)
- **许可证：**`MIT`
- **别名：**`bio-phylo-divergence-dating`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info phylo-divergence-dating
dnakit skills install phylo-divergence-dating --agent claude-code   # 或 --agent codex / cursor
```

## `phylo-modern-tree-inference` {#phylo-modern-tree-inference}

Infers maximum-likelihood phylogenetic trees with IQ-TREE2 and RAxML-NG -- model selection (ModelFinder), branch support (UFBoot2, SH-aLRT), concordance factors (gCF/sCF), partitioning, topology tests, and long-branch-attraction control. Covers why an ML tree inherits every flaw of the assumed model and the fixed alignment, why reported support measures repeatability under resampling and not correctness, why UFBoot uses a &gt;=95 cutoff and not the bootstrap-70 rule, and why a node with UFBoot 100 but gCF ~35 is essentially unresolved ILS rather than a clade. Use when inferring an ML tree, selecting a substitution or partition model, choosing or interpreting support measures, testing an a-priori topology, or diagnosing LBA. Routes model-free distance trees to distance-calculations, posteriors to bayesian-inference, and species trees under ILS to species-trees.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/phylogenetics/modern-tree-inference)
- **许可证：**`MIT`
- **运行依赖：**命令 `IQ-TREE2`
- **别名：**`bio-phylo-modern-tree-inference`, `fastreer`, `fasttree`, `iqtree`, `phylogenetics-builder`, `tooluniverse-phylogenetics`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info phylo-modern-tree-inference
dnakit skills install phylo-modern-tree-inference --agent claude-code   # 或 --agent codex / cursor
```

## `phylo-species-trees` {#phylo-species-trees}

Estimates species trees under the multispecies coalescent from per-locus gene trees with the modern ASTER astral binary (ASTRAL-III/wASTRAL/ASTRAL-Pro), plus SVDQuartets, BPP, and StarBEAST2. Covers why a species tree is not a gene tree, why each locus has its own genealogy that disagrees by incomplete lineage sorting (ILS) even with zero error, why concatenation is statistically inconsistent and positively misleading in the anomaly zone where more loci converge on the wrong tree with full support, why gene-tree estimation error biases summary methods, that localPP is not bootstrap and ASTRAL branch lengths are coalescent units, and how minority-quartet symmetry separates ILS from introgression. Use when multi-locus discordance, rapid radiations, anomaly-zone risk, concordance-factor interpretation, or concatenation-vs-coalescent choice arise. Routes per-locus gene-tree inference and gCF/sCF to modern-tree-inference, dating to divergence-dating, orthology to comparative-genomics/ortholog-inference.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/phylogenetics/species-trees)
- **许可证：**`MIT`
- **别名：**`bio-phylo-species-trees`, `phylophlan`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info phylo-species-trees
dnakit skills install phylo-species-trees --agent claude-code   # 或 --agent codex / cursor
```

## `phylo-tree-io` {#phylo-tree-io}

Read, write, and convert phylogenetic tree files with Biopython Bio.Phylo, and choose an annotation-preserving parser (treeio, DendroPy) when metadata matters. Covers why a tree file is a lossy serialization, why format conversion silently drops BEAST/MrBayes node annotations (posteriors, HPD intervals, rates), the Newick support-vs-label ambiguity that mislabels bootstrap values, and the Nexus TRANSLATE and rooted/unrooted traps. Use when parsing Newick, Nexus, NHX, phyloXML, or NeXML, converting between formats, handling posterior tree sets, or moving annotated BEAST trees without losing the credible intervals. Routes annotation-critical reads to DendroPy or treeio and orthology/alignment context to sibling skills.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/phylogenetics/tree-io)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `Bio.Phylo`, `dnakit`
- **别名：**`bio-phylo-tree-io`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info phylo-tree-io
dnakit skills install phylo-tree-io --agent claude-code   # 或 --agent codex / cursor
```

## `phylo-tree-manipulation` {#phylo-tree-manipulation}

Edit phylogenetic tree structure with Biopython Bio.Phylo, and treat rooting as a separate statistical inference rather than a display choice. Covers why most inference returns an unrooted tree so placing the root creates every ancestor/descendant and basal claim; why a distant or lonely outgroup misroots inside the ingroup via long-branch attraction; the outgroup/midpoint/MAD/MinVar/non-reversible-likelihood rooting tradeoffs; why pruning must suppress degree-2 nodes and sum their branch lengths or all patristic distances silently corrupt; and why collapsing by support makes SOFT (uncertainty) polytomies, not HARD (radiation) ones. Use when rooting, re-rooting, pruning or subsetting taxa, extracting a clade or induced subtree, collapsing low-support branches, resolving polytomies, or ladderizing. Routes clock-based rooting to divergence-dating, inference to modern-tree-inference, and reading/plotting to tree-io and tree-visualization.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/phylogenetics/tree-manipulation)
- **许可证：**`MIT`
- **别名：**`alterlab-etetoolkit`, `bio-phylo-tree-manipulation`, `etetoolkit`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info phylo-tree-manipulation
dnakit skills install phylo-tree-manipulation --agent claude-code   # 或 --agent codex / cursor
```

## `phylo-tree-visualization` {#phylo-tree-visualization}

Draw and export phylogenetic trees with Bio.Phylo plus matplotlib, and route rich figures to ggtree, ETE4, or iTOL. Covers why a tree figure is an argument not a neutral picture, the cladogram-vs-phylogram-vs-chronogram choice that hides or reveals rate and time, how ladderization manufactures a false arrow of progress, why an unlabeled support number always flatters the result (bootstrap vs posterior vs SH-aLRT vs UFBoot are different scales), why Bio.Phylo silently drops BEAST HPD bars so annotated Bayesian trees must go through treeio plus ggtree, and the tip-count and raster-vs-vector thresholds for legible publication figures. Use when drawing a tree, choosing a layout, coloring branches, showing or labeling support, exporting vector figures, or deciding a drawing tool. Routes annotation-preserving reads to tree-io, and rooting and ladderizing to tree-manipulation.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/phylogenetics/tree-visualization)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `Bio.Phylo`, `dnakit`
- **别名：**`bio-phylo-tree-visualization`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info phylo-tree-visualization
dnakit skills install phylo-tree-visualization --agent claude-code   # 或 --agent codex / cursor
```

## `phylogenomics` {#phylogenomics}

Build and validate marker-gene alignments and phylogenetic trees. Use when inferring evolutionary relationships, choosing models, or checking tree support and contamination.

- **来源：**[Copyright (c) 2026 Frederik Schulz @ `7a4eb5b`](https://github.com/fmschulz/omics-skills/tree/7a4eb5bace3f874c3be2ec697e9d2e695c374637/skills/bio-phylogenomics)
- **许可证：**`MIT`
- **别名：**`bio-phylogenomics`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info phylogenomics
dnakit skills install phylogenomics --agent claude-code   # 或 --agent codex / cursor
```

## `plant-genomics` {#plant-genomics}

Plant genomics and biology research — PlantReactome pathways, Ensembl Plants gene structure, POWO species taxonomy, UniProt annotation, KEGG plant pathways. Handles polyploidy (wheat hexaploidy etc.) and homeologous gene copies. Use for crop-gene annotation, plant secondary metabolism queries, and plant-disease/stress-response biology.

- **来源：**[ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-plant-genomics)
- **许可证：**`Apache-2.0`
- **别名：**`tooluniverse-plant-genomics`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info plant-genomics
dnakit skills install plant-genomics --agent claude-code   # 或 --agent codex / cursor
```

## `spacedust` {#spacedust}

Spacedust — de novo discovery of conserved gene clusters in microbial genomes using Foldseek structure comparison plus clustering/order-conservation statistics. Use when: (1) Finding conserved neighborhoods / modules in MAG collections, (2) Annotating dark genes via genomic context, (3) Recovering defense or other systems as clusters. Upstream: https://github.com/soedinglab/Spacedust. Route via microbial-mining / metagenomics-llm; complement defensefinder / antismash for known systems.

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/spacedust)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info spacedust
dnakit skills install spacedust --agent claude-code   # 或 --agent codex / cursor
```

## `tracking-taxonomy-updates` {#tracking-taxonomy-updates}

Track taxonomy changes and triage sequence assignments across NCBI, GTDB, ICTV, and eukaryotic frameworks. Use when comparing releases, resolving renamed taxa, or routing genomes, bins, and contigs by domain.

- **来源：**[Copyright (c) 2026 Frederik Schulz @ `7a4eb5b`](https://github.com/fmschulz/omics-skills/tree/7a4eb5bace3f874c3be2ec697e9d2e695c374637/skills/tracking-taxonomy-updates)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info tracking-taxonomy-updates
dnakit skills install tracking-taxonomy-updates --agent claude-code   # 或 --agent codex / cursor
```

## `unicore` {#unicore}

Unicore — universal structure-based core-gene phylogeny using Foldseek and ProstT5 (Steinegger lab). Use when: (1) Building phylogenies from proteomes via structural core genes, (2) 3Di MSAs for tree inference, (3) Scaling beyond sequence-only marker sets. Upstream: https://github.com/steineggerlab/unicore. Route via phylogenomics; depends on foldseek (± GPU ProstT5).

- **来源：**[Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/unicore)
- **许可证：**`MIT`
- **别名：**`lucaphylo`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info unicore
dnakit skills install unicore --agent claude-code   # 或 --agent codex / cursor
```

## `vitis-pangenome` {#vitis-pangenome}

Vitis / crop pangenome decisions: when to build a graph, PAV, reference choice for mixed cohorts and chips. Use for minigraph-cactus vs PGGB, PanGenie, or "should we build another genus graph".

- **来源：**[Xuzhen-Li/bioinfo-agent-skills contributors @ `6dea4be`](https://github.com/Xuzhen-Li/bioinfo-agent-skills/tree/6dea4bee7aade611d22f2f1958612328bc61ec83/skills/vitis-pangenome)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info vitis-pangenome
dnakit skills install vitis-pangenome --agent claude-code   # 或 --agent codex / cursor
```
