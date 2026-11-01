# Genome assembly

**11** Skills. A Skill is a set of instructions that guides an Agent through an analysis; the software, reference data and APIs it uses are separate dependencies. See [Skills](../skills.md) for searching, installing and using them; the Skills tab of the web demo shows each Skill in full.

| Skill | Source | Summary |
| --- | --- | --- |
| [`busco-assessor`](#busco-assessor) | clawbio | Genome, transcriptome, and protein completeness assessment via BUSCO v6. Agentic lineage routing from organism description, all three BUSCO modes, auto-lineage… |
| [`genome-assembly-assembly-polishing`](#genome-assembly-assembly-polishing) | bioskills | Decides whether and how to polish a draft genome assembly to raise consensus accuracy (QV) with read-type-matched tools - Racon and medaka (ONT consensus), dor… |
| [`genome-assembly-assembly-qc`](#genome-assembly-assembly-qc) | bioskills | Evaluates genome assembly quality across the three orthogonal axes - contiguity (QUAST auN/NG50/NGx, not bare N50), completeness (BUSCO/compleasm gene-space pl… |
| [`genome-assembly-contamination-detection`](#genome-assembly-contamination-detection) | bioskills | Detects and removes contamination in genome assemblies via two disjoint workflows - foreign-sequence screening of a single-organism (eukaryote/isolate) assembl… |
| [`genome-assembly-genome-profiling`](#genome-assembly-genome-profiling) | bioskills | Profiles a genome from raw reads BEFORE assembly with a k-mer spectrum (KMC or Jellyfish histogram), then models it with GenomeScope2 to estimate genome size,… |
| [`genome-assembly-hifi-assembly`](#genome-assembly-hifi-assembly) | bioskills | Assembles haplotype-resolved diploid and telomere-to-telomere (T2T) genomes from PacBio HiFi reads with hifiasm (HiFi-only, Hi-C, or trio phasing) and verkko (… |
| [`genome-assembly-long-read-assembly`](#genome-assembly-long-read-assembly) | bioskills | Assembles genomes de novo from noisy long reads (Oxford Nanopore R9/R10/Dorado, PacBio CLR) with Flye (repeat graph), Canu (correct-trim-assemble OLC), NextDen… |
| [`genome-assembly-metagenome-assembly`](#genome-assembly-metagenome-assembly) | bioskills | Assembles microbial-community sequencing into metagenome-assembled genomes (MAGs) with metaFlye (ONT), metaSPAdes/MEGAHIT (Illumina), and hifiasm-meta/metaMDBG… |
| [`genome-assembly-scaffolding`](#genome-assembly-scaffolding) | bioskills | Orders and orients assembled contigs into chromosome-scale scaffolds from long-range linking data, inserting N-gap spacers (adds no sequence). Covers Hi-C/Omni… |
| [`genome-assembly-short-read-assembly`](#genome-assembly-short-read-assembly) | bioskills | Assembles a genome de novo from Illumina short reads with SPAdes (isolate/careful/sc/meta/plasmid/rna modes), MEGAHIT (low-memory, huge datasets), Unicycler (b… |
| [`optimize-hifi-diploid-assembly`](#optimize-hifi-diploid-assembly) | myriad | Activate when an agent must generate, compare, quality-control, or select a diploid PacBio HiFi de novo genome assembly using hifiasm; when requests mention Hi… |

## `busco-assessor` {#busco-assessor}

Genome, transcriptome, and protein completeness assessment via BUSCO v6. Agentic lineage routing from organism description, all three BUSCO modes, auto-lineage support, and full demo mode without the BUSCO binary.

- **Source:** [ClawBio Contributors; ClawBio contributors @ `1fcecb7`](https://github.com/ClawBio/ClawBio/tree/1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d/skills/busco-assessor)
- **License:** `MIT`
- **Runs with:** commands `python3`
- **Also known as:** `busco`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info busco-assessor
dnakit skills install busco-assessor --agent claude-code   # or --agent codex / cursor
```

## `genome-assembly-assembly-polishing` {#genome-assembly-assembly-polishing}

Decides whether and how to polish a draft genome assembly to raise consensus accuracy (QV) with read-type-matched tools - Racon and medaka (ONT consensus), dorado polish, Polypolish and pypolca (Illumina, repeat-aware), Pilon (legacy short-read), NextPolish/NextPolish2, Hapo-G (haplotype-aware), ntEdit, and DeepPolisher/PEPPER-Margin-DeepVariant for human. Covers the do-not-polish-HiFi rule, the medaka basecaller-model footgun, held-out Merqury QV as the only honest stop signal, and the haplotype-collapse trap. Use when correcting homopolymer indels or residual SNPs in a long-read assembly, deciding if a HiFi assembly needs polishing, or choosing an ONT vs hybrid vs short-read polishing chain.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-assembly/assembly-polishing)
- **License:** `MIT`
- **Runs with:** commands `Pilon`
- **Also known as:** `bio-genome-assembly-assembly-polishing`, `pilon`, `racon`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-assembly-assembly-polishing
dnakit skills install genome-assembly-assembly-polishing --agent claude-code   # or --agent codex / cursor
```

## `genome-assembly-assembly-qc` {#genome-assembly-assembly-qc}

Evaluates genome assembly quality across the three orthogonal axes - contiguity (QUAST auN/NG50/NGx, not bare N50), completeness (BUSCO/compleasm gene-space plus Merqury k-mer completeness), and correctness (reference-free Merqury QV, Inspector/CRAQ structural errors, asmgene false-duplication/collapse). Covers why N50 is the most-gamed metric, why QV measured on the polishing reads is circular, distinguishing uncollapsed haplotigs from real WGD, and the EBP/VGP 6.C.Q40 standard. Use when judging whether an assembly is good enough to annotate or publish, comparing assemblers, diagnosing a fragmented or duplicated assembly, or assessing a phased diploid assembly.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-assembly/assembly-qc)
- **License:** `MIT`
- **Runs with:** commands `QUAST`
- **Also known as:** `bio-genome-assembly-assembly-qc`, `quast`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-assembly-assembly-qc
dnakit skills install genome-assembly-assembly-qc --agent claude-code   # or --agent codex / cursor
```

## `genome-assembly-contamination-detection` {#genome-assembly-contamination-detection}

Detects and removes contamination in genome assemblies via two disjoint workflows - foreign-sequence screening of a single-organism (eukaryote/isolate) assembly with NCBI FCS-GX (GenBank-submission-mandatory), FCS-adaptor, and BlobToolKit blob plots; and MAG/bin quality assessment with CheckM2 plus GUNC (chimerism) plus GTDB-Tk taxonomy, judged against MIMAG. Covers why CheckM2 alone is blind to disjoint-marker chimeras, the FCS-GX RAM wall, organelle/NUMT triage, strain heterogeneity, and the HGT-vs-contamination (tardigrade) trap. Use when screening an assembly for foreign contamination before GenBank submission, assessing MAG completeness/contamination/chimerism, deciding which contigs to remove, or distinguishing real HGT from contaminant contigs.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-assembly/contamination-detection)
- **License:** `MIT`
- **Runs with:** commands `CheckM2`
- **Also known as:** `bio-genome-assembly-contamination-detection`, `mag-qc`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-assembly-contamination-detection
dnakit skills install genome-assembly-contamination-detection --agent claude-code   # or --agent codex / cursor
```

## `genome-assembly-genome-profiling` {#genome-assembly-genome-profiling}

Profiles a genome from raw reads BEFORE assembly with a k-mer spectrum (KMC or Jellyfish histogram), then models it with GenomeScope2 to estimate genome size, heterozygosity, repeat content, and ploidy, and Smudgeplot to infer ploidy from heterozygous k-mer pairs (diploid AB vs triploid AAB vs tetraploid AABB). Covers choosing k via Merqury best_k.sh, the k-mer-coverage vs sequencing-coverage confusion, reading het/repeat/contamination/organelle peaks, why noisy ONT must not be used for counting, and how the estimate becomes the NG50 denominator, the Flye -g value, the hifiasm --hom-cov/purge setting, and the 1.5-2x-too-big haplotig sanity check. Use when starting any de novo assembly, deciding whether short reads can work, estimating genome size for an unknown organism, diagnosing ploidy, or sanity-checking an assembly's size against expectation.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-assembly/genome-profiling)
- **License:** `MIT`
- **Runs with:** commands `GenomeScope2`
- **Also known as:** `bio-genome-assembly-genome-profiling`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-assembly-genome-profiling
dnakit skills install genome-assembly-genome-profiling --agent claude-code   # or --agent codex / cursor
```

## `genome-assembly-hifi-assembly` {#genome-assembly-hifi-assembly}

Assembles haplotype-resolved diploid and telomere-to-telomere (T2T) genomes from PacBio HiFi reads with hifiasm (HiFi-only, Hi-C, or trio phasing) and verkko (HiFi + ultralong ONT for T2T), extracting contigs from GFA and routing phasing QC to k-mer/trio metrics. Covers why a primary assembly is a haplotype mosaic that exists in no cell, partial-vs-full phasing (the .bp. vs .dip. filename convention), the purge-default trap on inbred samples, the --hom-cov coverage-estimate alarm, and verkko-vs-hifiasm for T2T. Use when assembling a diploid eukaryote from HiFi, phasing haplotypes with parents (trio) or Hi-C, deciding whether to chase T2T, or diagnosing switch errors invisible to N50/BUSCO/QV.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-assembly/hifi-assembly)
- **License:** `MIT`
- **Runs with:** commands `hifiasm`
- **Also known as:** `bio-genome-assembly-hifi-assembly`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-assembly-hifi-assembly
dnakit skills install genome-assembly-hifi-assembly --agent claude-code   # or --agent codex / cursor
```

## `genome-assembly-long-read-assembly` {#genome-assembly-long-read-assembly}

Assembles genomes de novo from noisy long reads (Oxford Nanopore R9/R10/Dorado, PacBio CLR) with Flye (repeat graph), Canu (correct-trim-assemble OLC), NextDenovo, Shasta, Raven, wtdbg2, or miniasm, and reconciles bacterial assemblies into a consensus with Trycycler/Autocycler. Covers matching the input flag to the basecaller era (--nano-hq vs --nano-raw), why a raw long-read assembly is contiguous but low-QV and not finished until polished, haplotig false-duplication and purge_dups, coverage and read-N50 as non-substitutable inputs, and mid-read adapter de-chimerization. Use when assembling a bacterial or eukaryotic genome from ONT or PacBio noisy reads, choosing a long-read assembler, or diagnosing an over-collapsed or duplicated assembly. For PacBio HiFi use hifi-assembly instead.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-assembly/long-read-assembly)
- **License:** `MIT`
- **Runs with:** commands `Flye`
- **Also known as:** `bio-genome-assembly-long-read-assembly`, `metaflye`, `metamdbg`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-assembly-long-read-assembly
dnakit skills install genome-assembly-long-read-assembly --agent claude-code   # or --agent codex / cursor
```

## `genome-assembly-metagenome-assembly` {#genome-assembly-metagenome-assembly}

Assembles microbial-community sequencing into metagenome-assembled genomes (MAGs) with metaFlye (ONT), metaSPAdes/MEGAHIT (Illumina), and hifiasm-meta/metaMDBG (PacBio HiFi), then recovers genomes via multi-binner consolidation (MetaBAT2, MaxBin2, CONCOCT, SemiBin2, VAMB -&gt; DAS_Tool) and QCs them against MIMAG with CheckM2, GUNC, and GTDB-Tk. Covers why a metagenome is not a genome (uneven coverage, micro-diversity, strain collapse to consensus), differential-coverage binning, co-assembly vs per-sample, the rRNA-operon collapse that fails short-read MAGs, and strain resolution with inStrain. Use when reconstructing genomes from a microbiome, soil, ocean, or gut community, recovering MAGs, or resolving strain-level variation.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-assembly/metagenome-assembly)
- **License:** `MIT`
- **Runs with:** commands `metaFlye`
- **Also known as:** `bio-genome-assembly-metagenome-assembly`, `megahit`, `metaspades`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-assembly-metagenome-assembly
dnakit skills install genome-assembly-metagenome-assembly --agent claude-code   # or --agent codex / cursor
```

## `genome-assembly-scaffolding` {#genome-assembly-scaffolding}

Orders and orients assembled contigs into chromosome-scale scaffolds from long-range linking data, inserting N-gap spacers (adds no sequence). Covers Hi-C/Omni-C scaffolding (YaHS, SALSA2, 3D-DNA/Juicer), Hi-C read-mapping prerequisites (map each end separately, no mate rescue, dedup, enzyme-aware), reading the contact map for misjoins/inversions/false-duplications, manual curation in Juicebox/PretextView (the VGP/DToL standard), reference-guided scaffolding (RagTag) and its karyotype-erasure hazard, genetic-map (ALLMAPS) and Bionano optical-map integration, chimera-breaking before scaffolding, gap-filling, and telomere/contig-vs-scaffold-N50 QC (tidk). Use when turning contigs into chromosomes with Hi-C, integrating a linkage map or optical map, choosing a scaffolder by available linking data, or judging whether a chromosome-scale assembly is trustworthy.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-assembly/scaffolding)
- **License:** `MIT`
- **Runs with:** commands `YaHS`
- **Also known as:** `bio-genome-assembly-scaffolding`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-assembly-scaffolding
dnakit skills install genome-assembly-scaffolding --agent claude-code   # or --agent codex / cursor
```

## `genome-assembly-short-read-assembly` {#genome-assembly-short-read-assembly}

Assembles a genome de novo from Illumina short reads with SPAdes (isolate/careful/sc/meta/plasmid/rna modes), MEGAHIT (low-memory, huge datasets), Unicycler (bacterial finishing/hybrid), MaSuRCA (large hybrid), ABySS (Bloom-filter), and Platanus (heterozygous diploids), using multi-k de Bruijn graphs. Covers the repeat-resolution limit, why N50 plateaus at the genome not the depth, GenomeScope2 k-mer profiling first, the heterozygosity/haplotig trap, error-correction erasing rare alleles, GC dropout, and NG50/auN/BUSCO reporting. Use when assembling a bacterial isolate, fungal, small-eukaryotic, single-cell, or metagenome genome from Illumina reads, or when deciding whether short reads can even produce the assembly being asked for.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-assembly/short-read-assembly)
- **License:** `MIT`
- **Runs with:** commands `SPAdes`
- **Also known as:** `bio-genome-assembly-short-read-assembly`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-assembly-short-read-assembly
dnakit skills install genome-assembly-short-read-assembly --agent claude-code   # or --agent codex / cursor
```

## `optimize-hifi-diploid-assembly` {#optimize-hifi-diploid-assembly}

Activate when an agent must generate, compare, quality-control, or select a diploid PacBio HiFi de novo genome assembly using hifiasm; when requests mention HiFi read assembly, haplotype-resolved assembly, hifiasm parameter optimization, primary and alternate contigs, trio binning, Hi-C phasing, Merqury QV, k-mer completeness, BUSCO, QUAST, contamination screening, assembly Pareto ranking, or reproducible assembly QC.

- **Source:** [MYRIAD by Srihari Muralikrishnan and Kalaris Labs. https://github.com/KalarisLabs/myriad @ `5efbc5e`](https://github.com/KalarisLabs/Myriad/tree/5efbc5e9940fd37608f1d75cccac2befebad1215/skills/optimize-hifi-diploid-assembly)
- **License:** `Apache-2.0 AND CC-BY-4.0`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info optimize-hifi-diploid-assembly
dnakit skills install optimize-hifi-diploid-assembly --agent claude-code   # or --agent codex / cursor
```
