# Genome annotation

**15** Skills. A Skill is a set of instructions that guides an Agent through an analysis; the software, reference data and APIs it uses are separate dependencies. See [Skills](../skills.md) for searching, installing and using them; the Skills tab of the web demo shows each Skill in full.

| Skill | Source | Summary |
| --- | --- | --- |
| [`alphagem`](#alphagem) | pku-metagenomics | AlphaGEM — genome-scale metabolic modelling aided by protein structure alignment and deep-learning dark-metabolism mining. Use when: (1) Building or refining G… |
| [`dnakit-prokaryotic-genes`](#dnakit-prokaryotic-genes) | dnakit | Predict prokaryotic and metagenomic CDS with optional Pyrodigal through DNAKit, returning coordinates, coding sequences and proteins. Use for bacterial or arch… |
| [`expression-matrix-gene-id-mapping`](#expression-matrix-gene-id-mapping) | bioskills | Maps between gene identifier systems (Ensembl, Entrez, HGNC symbol, UniProt, RefSeq, MANE) using AnnotationDbi, biomaRt, mygene, pyensembl, and Ensembl REST. E… |
| [`genome-annotation-annotation-qc`](#genome-annotation-annotation-qc) | bioskills | Assesses the quality and completeness of a genome annotation with BUSCO (conserved single-copy ortholog recovery), OMArk (proteome completeness, consistency, a… |
| [`genome-annotation-annotation-transfer`](#genome-annotation-annotation-transfer) | bioskills | Transfers gene annotations between genome assemblies via coordinate liftover (UCSC liftOver, CrossMap for same-species version updates) or feature/sequence pro… |
| [`genome-annotation-eukaryotic-gene-prediction`](#genome-annotation-eukaryotic-gene-prediction) | bioskills | Predicts protein-coding gene structures (exons, introns, UTRs) in eukaryotic genomes with BRAKER3 (RNA-seq + protein evidence), BRAKER1/BRAKER2, GALBA (protein… |
| [`genome-annotation-functional-annotation`](#genome-annotation-functional-annotation) | bioskills | Assigns GO terms, Pfam/InterPro domains, KEGG orthologs, EC numbers, and product names to predicted proteins using eggNOG-mapper (orthology), InterProScan (dom… |
| [`genome-annotation-ncrna-annotation`](#genome-annotation-ncrna-annotation) | bioskills | Identifies non-coding RNAs (tRNA, rRNA, snoRNA, snRNA, riboswitches, sRNAs) using Infernal covariance-model search against Rfam, tRNAscan-SE 2.0 for tRNA, barr… |
| [`genome-annotation-prokaryotic-annotation`](#genome-annotation-prokaryotic-annotation) | bioskills | Annotates bacterial and archaeal genomes (isolates, MAGs, plasmids) with Bakta (active versioned databases, NCBI-compliant output) or Prokka (legacy), producin… |
| [`genome-annotation-repeat-annotation`](#genome-annotation-repeat-annotation) | bioskills | Discovers, classifies, and masks repetitive elements and transposable elements with RepeatModeler2 (de novo family library), RepeatMasker (masking against a li… |
| [`interpro-database`](#interpro-database) | user335-deepmind | Identify domains, families, and sites in proteins; find all proteins in a family or sharing a domain; explore species distribution for a domain; annotate genom… |
| [`qtl-candidate-gene`](#qtl-candidate-gene) | haploqtl | Interpret a fine-mapped QTL / introgression interval for a stated trait in tomato. Given an interval (SL4.0 coordinates) and the trait it was mapped for, lists… |
| [`systems-biology-gene-essentiality`](#systems-biology-gene-essentiality) | bioskills | Performs in-silico single and double gene deletions, condition-dependent essentiality, and synthetic-lethality screens on genome-scale metabolic models with CO… |
| [`systems-biology-metabolic-reconstruction`](#systems-biology-metabolic-reconstruction) | bioskills | Builds draft genome-scale metabolic models from an annotated genome using CarveMe (top-down carving of a BiGG universal model) or gapseq (bottom-up pathway-evi… |
| [`vitis-te`](#vitis-te) | plant-skills | Transposable-element annotation and curation for plant genomes (Vitis). Use for EDTA, TEtrimmer, TEsorter, curated libraries, RepeatMasker, LAI. Do not hook a… |

## `alphagem` {#alphagem}

AlphaGEM — genome-scale metabolic modelling aided by protein structure alignment and deep-learning dark-metabolism mining. Use when: (1) Building or refining GEMs beyond sequence-only gap filling, (2) Integrating structure signals into metabolic networks, (3) Comparing to gapseq / VMH workflows. Upstream: https://github.com/hongzhonglu/AlphaGEM. Route via metagenomics-llm; pair with gapseq / vmh for classical GEM context.

- **Source:** [Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/alphagem)
- **License:** `MIT`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info alphagem
dnakit skills install alphagem --agent claude-code   # or --agent codex / cursor
```

## `dnakit-prokaryotic-genes` {#dnakit-prokaryotic-genes}

Predict prokaryotic and metagenomic CDS with optional Pyrodigal through DNAKit, returning coordinates, coding sequences and proteins. Use for bacterial or archaeal sequence annotation.

- **Source:** [DNAKit contributors; method authors cited in SKILL.md and DNAKIT.md @ `dnakit-`](https://github.com/mapengsen/DNAKit/tree/main/skills/dnakit-prokaryotic-genes)
- **License:** `MIT`
- **Runs with:** Python modules `dnakit`, `pyrodigal`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info dnakit-prokaryotic-genes
dnakit skills install dnakit-prokaryotic-genes --agent claude-code   # or --agent codex / cursor
```

## `expression-matrix-gene-id-mapping` {#expression-matrix-gene-id-mapping}

Maps between gene identifier systems (Ensembl, Entrez, HGNC symbol, UniProt, RefSeq, MANE) using AnnotationDbi, biomaRt, mygene, pyensembl, and Ensembl REST. Encodes Ensembl version stripping with GENCODE _PAR_Y preservation, the Ziemann 2016 Excel autocorrect debacle and Bruford 2020 HGNC renames (SEPT*-&gt;SEPTIN*, MARCH*-&gt;MARCHF*, MARC*-&gt;MTARC*, DEC1-&gt;DELEC1), OCT4/POU5F1 alias resolution, biomaRt archive endpoints for release pinning, the `filters` (plural) gotcha, MANE Select for clinical reporting, cross-species orthology via Ensembl Compara / OMA / OrthoDB, and tx2gene construction for tximport. Use when converting gene IDs across systems, handling renamed symbols, building tx2gene, pinning to a specific Ensembl release for reproducibility, or mapping cross-species orthologs.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/expression-matrix/gene-id-mapping)
- **License:** `MIT`
- **Also known as:** `bio-expression-matrix-gene-id-mapping`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info expression-matrix-gene-id-mapping
dnakit skills install expression-matrix-gene-id-mapping --agent claude-code   # or --agent codex / cursor
```

## `genome-annotation-annotation-qc` {#genome-annotation-annotation-qc}

Assesses the quality and completeness of a genome annotation with BUSCO (conserved single-copy ortholog recovery), OMArk (proteome completeness, consistency, and contamination), CheckM2 (prokaryotic completeness/contamination), and a gene-set sanity panel (gene count, mono-exonic fraction, protein-length distribution, mRNA:gene ratio, coding density). Covers the assembly-BUSCO-vs-proteome-BUSCO diagnostic, what BUSCO-Duplicated really means, why gene count is a vanity metric, and the QC of transferred annotations. Use when judging whether an annotation is good enough to publish or submit, diagnosing a suspect annotation, or comparing annotation completeness across pipelines.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-annotation/annotation-qc)
- **License:** `MIT`
- **Runs with:** commands `BUSCO`
- **Also known as:** `bio-genome-annotation-annotation-qc`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-annotation-annotation-qc
dnakit skills install genome-annotation-annotation-qc --agent claude-code   # or --agent codex / cursor
```

## `genome-annotation-annotation-transfer` {#genome-annotation-annotation-transfer}

Transfers gene annotations between genome assemblies via coordinate liftover (UCSC liftOver, CrossMap for same-species version updates) or feature/sequence projection (Liftoff for same/close species, miniprot for protein-level cross-species, TOGA/GeMoMa/CAT for distant clades). Covers the coordinate-vs-projection decision by divergence, why a successful lift is not biological confirmation, reference bias, the silent-dropping of unmapped features, build/PAR/MHC/inversion hazards, and transfer-vs-de-novo validation. Use when annotating a new assembly of a species with an existing reference, harmonizing coordinates across builds, or mapping annotations across related species.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-annotation/annotation-transfer)
- **License:** `MIT`
- **Runs with:** commands `Liftoff`
- **Also known as:** `bio-genome-annotation-annotation-transfer`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-annotation-annotation-transfer
dnakit skills install genome-annotation-annotation-transfer --agent claude-code   # or --agent codex / cursor
```

## `genome-annotation-eukaryotic-gene-prediction` {#genome-annotation-eukaryotic-gene-prediction}

Predicts protein-coding gene structures (exons, introns, UTRs) in eukaryotic genomes with BRAKER3 (RNA-seq + protein evidence), BRAKER1/BRAKER2, GALBA (protein-only), Funannotate (fungi), GeMoMa (homology projection), or Helixer/Tiberius (deep-learning ab initio). Covers the evidence-first tool decision, mandatory soft-masking, the training-set-quality-dominates principle, OrthoDB clade-partition selection, the one-isoform-per-locus and missing-UTR traps, merge/split errors, and reference bias against orphan genes. Use when annotating a newly assembled eukaryotic genome, choosing a gene-prediction pipeline based on available evidence, or diagnosing a poor annotation.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-annotation/eukaryotic-gene-prediction)
- **License:** `MIT`
- **Runs with:** commands `BRAKER3`
- **Also known as:** `bio-genome-annotation-eukaryotic-gene-prediction`, `metaeuk`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-annotation-eukaryotic-gene-prediction
dnakit skills install genome-annotation-eukaryotic-gene-prediction --agent claude-code   # or --agent codex / cursor
```

## `genome-annotation-functional-annotation` {#genome-annotation-functional-annotation}

Assigns GO terms, Pfam/InterPro domains, KEGG orthologs, EC numbers, and product names to predicted proteins using eggNOG-mapper (orthology), InterProScan (domain signatures), and KofamScan (KEGG), routing specialized functions to dbCAN/antiSMASH/AMRFinderPlus/SignalP. Covers the orthology-vs-domain-vs-homology paradigms, the annotation-error percolation cascade, domain-presence-is-not-function, GO IEA circularity in enrichment, evidence tiering, and bit-score/coverage thresholds. Use when adding functional annotation to predicted genes, choosing between eggNOG-mapper and InterProScan, or judging how much to trust a functional label.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-annotation/functional-annotation)
- **License:** `MIT`
- **Runs with:** commands `eggNOG-mapper`
- **Also known as:** `bio-annotation`, `bio-genome-annotation-functional-annotation`, `eggnog-mapper`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-annotation-functional-annotation
dnakit skills install genome-annotation-functional-annotation --agent claude-code   # or --agent codex / cursor
```

## `genome-annotation-ncrna-annotation` {#genome-annotation-ncrna-annotation}

Identifies non-coding RNAs (tRNA, rRNA, snoRNA, snRNA, riboswitches, sRNAs) using Infernal covariance-model search against Rfam, tRNAscan-SE 2.0 for tRNA, barrnap for rRNA, and ARAGORN for tmRNA, plus the small-RNA-seq boundary for miRNA and the transcript-assembly boundary for lncRNA. Covers the structure-conserved-not-sequence-conserved principle (why BLAST fails), GA-threshold and clan-competition correctness, tRNAscan-SE domain modes and pseudogene flags, rDNA copy-number collapse, and why homology annotation is a recall floor. Use when performing genome-wide ncRNA annotation, choosing the right tool for an RNA class, or interpreting ncRNA counts.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-annotation/ncrna-annotation)
- **License:** `MIT`
- **Runs with:** commands `Infernal`
- **Also known as:** `barrnap`, `bio-genome-annotation-ncrna-annotation`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-annotation-ncrna-annotation
dnakit skills install genome-annotation-ncrna-annotation --agent claude-code   # or --agent codex / cursor
```

## `genome-annotation-prokaryotic-annotation` {#genome-annotation-prokaryotic-annotation}

Annotates bacterial and archaeal genomes (isolates, MAGs, plasmids) with Bakta (active versioned databases, NCBI-compliant output) or Prokka (legacy), producing GFF3/GenBank/EMBL/FASTA with INSDC locus tags. Covers Bakta-vs-Prokka-vs-PGAP-vs-DFAST choice, light-vs-full database tiers, translation-table selection (11/4/25), archaeal and leaderless-gene caveats, the small-ORF blind spot, pseudogene-vs-phase-variation, the pangenome re-annotation trap, and submission compliance. Use when annotating a newly assembled prokaryotic genome, choosing an annotation tool, re-annotating a collection for pangenomics, or preparing annotations for NCBI/DDBJ submission.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-annotation/prokaryotic-annotation)
- **License:** `MIT`
- **Runs with:** commands `Bakta`
- **Also known as:** `bakta`, `bakta-genome-annotation`, `bio-genome-annotation-prokaryotic-annotation`, `prodigal`, `prokka`, `prokka-genome-annotation`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-annotation-prokaryotic-annotation
dnakit skills install genome-annotation-prokaryotic-annotation --agent claude-code   # or --agent codex / cursor
```

## `genome-annotation-repeat-annotation` {#genome-annotation-repeat-annotation}

Discovers, classifies, and masks repetitive elements and transposable elements with RepeatModeler2 (de novo family library), RepeatMasker (masking against a library), EDTA (plant/structural TEs), or EarlGrey (auto-curating wrapper), and quantifies TE expression from RNA-seq with TEtranscripts/SQuIRE. Covers de-novo-library-as-curation-project, soft-vs-hard masking, the domesticated-gene over-masking massacre, Dfam-vs-RepBase, TE classification (Class I/II, family-vs-copy), Kimura repeat landscapes, LAI, and the RNA-seq multimapping problem. Use when masking repeats before gene prediction, building a TE library for a non-model genome, or analyzing transposable-element content or expression.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-annotation/repeat-annotation)
- **License:** `MIT`
- **Runs with:** commands `RepeatMasker`
- **Also known as:** `bio-genome-annotation-repeat-annotation`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-annotation-repeat-annotation
dnakit skills install genome-annotation-repeat-annotation --agent claude-code   # or --agent codex / cursor
```

## `interpro-database` {#interpro-database}

Identify domains, families, and sites in proteins; find all proteins in a family or sharing a domain; explore species distribution for a domain; annotate genomes with protein families and GO terms. InterPro combines 14 databases (e.g., Pfam, CDD) into one searchable resource. InterPro-N significantly expands annotation and sequence coverage with deep learning. Includes domain architecture (IDA) search.

- **Source:** [Copyright 2026 Google LLC @ `6883275`](https://github.com/google-deepmind/science-skills/tree/68832757cbbf941c620b71df5756cf6e5cc287b0/skills/interpro_database)
- **License:** `Apache-2.0 AND CC-BY-4.0`
- **Runs with:** Python modules `polite_http`; commands `uv`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info interpro-database
dnakit skills install interpro-database --agent claude-code   # or --agent codex / cursor
```

## `qtl-candidate-gene` {#qtl-candidate-gene}

Interpret a fine-mapped QTL / introgression interval for a stated trait in tomato. Given an interval (SL4.0 coordinates) and the trait it was mapped for, lists the genes in the interval from SGN ITAG4.1, enriches them with protein function from UniProt, ranks candidate genes by plausibility for that trait, grounds every literature citation with a live PubMed search, optionally flags marker-assisted-selection (MAS) markers from a two-group phenotype contrast, then self-verifies the draft — dropping any gene not in the interval and any citation that doesn't resolve — and delivers a candidate-gene report with a verification stamp. Works for any tomato trait — disease resistance, fruit size/quality, color, plant architecture, abiotic-stress tolerance. Use when a user has a QTL interval plus a trait and wants candidate or causal genes, mechanistic hypotheses, or MAS/KASP markers, or wants to interpret haploqtl clustering output.

- **Source:** [Copyright 2026 Taylor A. Anderson; haploqtl contributors @ `210eace`](https://github.com/tayloranthonyanderson/haploqtl/tree/210eace2912c7bf11d119cb503c89f8f6cad479f/skills/qtl-candidate-gene)
- **License:** `MIT`
- **Runs with:** Python modules `numpy`, `allel`; commands `python3`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info qtl-candidate-gene
dnakit skills install qtl-candidate-gene --agent claude-code   # or --agent codex / cursor
```

## `systems-biology-gene-essentiality` {#systems-biology-gene-essentiality}

Performs in-silico single and double gene deletions, condition-dependent essentiality, and synthetic-lethality screens on genome-scale metabolic models with COBRApy, evaluating gene-protein-reaction rules and comparing FBA re-optimization against MOMA/ROOM minimal-adjustment. Use when predicting essential genes, finding synthetic-lethal pairs for drug targets, choosing a growth cutoff, deciding FBA vs MOMA vs ROOM for a knockout, making essentiality medium-specific to match an experiment, or validating predictions against Keio/Tn-seq/CRISPR screens with MCC.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/systems-biology/gene-essentiality)
- **License:** `MIT`
- **Runs with:** Python modules `cobrapy`, `dnakit`
- **Also known as:** `bio-systems-biology-gene-essentiality`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info systems-biology-gene-essentiality
dnakit skills install systems-biology-gene-essentiality --agent claude-code   # or --agent codex / cursor
```

## `systems-biology-metabolic-reconstruction` {#systems-biology-metabolic-reconstruction}

Builds draft genome-scale metabolic models from an annotated genome using CarveMe (top-down carving of a BiGG universal model) or gapseq (bottom-up pathway-evidence reconstruction), then loads and sanity-checks the draft in COBRApy. Use when creating a model for an organism without one, choosing between CarveMe and gapseq, gap-filling to a target medium, understanding why a draft that grows is still only a hypothesis, handling BiGG-vs-ModelSEED namespace mismatch, or preparing a draft for curation and community modeling.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/systems-biology/metabolic-reconstruction)
- **License:** `MIT`
- **Runs with:** commands `CarveMe`
- **Also known as:** `bio-systems-biology-metabolic-reconstruction`, `gapseq`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info systems-biology-metabolic-reconstruction
dnakit skills install systems-biology-metabolic-reconstruction --agent claude-code   # or --agent codex / cursor
```

## `vitis-te` {#vitis-te}

Transposable-element annotation and curation for plant genomes (Vitis). Use for EDTA, TEtrimmer, TEsorter, curated libraries, RepeatMasker, LAI. Do not hook a classifier-only tool as a genome-wide scanner.

- **Source:** [Xuzhen-Li/bioinfo-agent-skills contributors @ `6dea4be`](https://github.com/Xuzhen-Li/bioinfo-agent-skills/tree/6dea4bee7aade611d22f2f1958612328bc61ec83/skills/vitis-te)
- **License:** `MIT`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info vitis-te
dnakit skills install vitis-te --agent claude-code   # or --agent codex / cursor
```
