# DNAKit function module index

> Updated: 2026-09-13. Checked against current workspace public APIs. See the [feature tree overview](function_tree.md) for the complete hierarchy. This index links feature guides and public API extensions; cross-category links do not count as additional capabilities.

## 1) Data preparation

Added APIs: [genes, sequences and downloads](../native/sequences.md) (41), [ENCODE](../native/encode.md) (20), [organisms and microbiomes](../native/organisms.md) (49), and [literature and studies](../native/literature.md) (17). Each has individual parameters and returns.

- [Core Data Object](01_core_objects.md) (8 items)
- [File reading and writing](02_io_data.md) (Category 4; large file capability as optional mode for IO-001)
- [Index and metadata APIs](../index.md) (FASTA/FASTQ random access, BED/GFF3/AGP, metadata and result tables)
- [Download](15_download.md) (common species names, reference genome FASTA and public database data)
- [Data query](22_database_query.md) (public database query capability)
- [Download and external index APIs](../index.md) (general file downloads, BLAST/BWA/Bowtie2/minimap2/samtools indexes)
- [Sequencing records and regions](27_bioinformatics.md) (FASTX/gzip, pairing, four-line FASTQ recovery, sorting, indexing, BAM consensus)
- [Legality check](03_validation.md) (unified check entrance)

## 2) Data processing

Added sequence statistics, VCF processing, digestion and assembly APIs are in [sequence analysis, VCF, cloning and inheritance](../native/analysis.md). Cross-links from analysis do not count functions again.

- [DNA sequence standardization](03_standardization.md) (STD-001 character standardization)
- [Gap and AGP APIs](../index.md#standardization) (uppercase N runs to Gaps, AGP assembly, provenance)
- **Intra-sequence operations**
  - [OPS-001 sequence completion](04_ops_001_sequence_direction.md)
  - [OPS-002 Transcription and Translation](04_ops_002_transcription_translation.md)
  - [OPS-003 subsequence extraction](04_ops_003_subsequence_extraction.md)
  - [OPS-004 Sequence Editing](04_ops_004_sequence_editing.md)
  - [OPS-005 Sequence Generation](04_ops_005_mutation_generation.md)
  - [OPS-006 sequence splicing](04_ops_006_sequence_concatenation.md)
  - **OPS-007**
    - [Trimming](04_ops_007_trimming.md)
    - [Masking](04_ops_007_masking.md)
  - [OPS-008 Circular Sequence Operation](04_ops_008_circular_sequence_operations.md)
  - [OPS-010 Sequence Segmentation](17_sequence_chunking.md)
- [Sequence Deduplication](10_deduplication.md) (DATA-001–005)
- [Sequence Clustering](10_clustering.md) (DATA-007–011, DATA-027 Neural Network Clustering)
- [Data set organization and division](10_datasets.md) (DATA-012–018, DATA-023, including custom label division)
- [Dataset API extensions](../index.md#datasets) (temporal/joint splits, metadata exclusion, and split quality)
- [Molecular biology simulation and design](../index.md#molecular-biology) (digestion, ligation, PCR, primer design, assembly, CRISPR, codon/rule optimization, and mutation libraries)
- [Constraint optimization and duplex cloning](27_bioinformatics.md) (native rules/weighted objectives, circular DNA, cloning steps, homologous assembly, Sanger consensus)
- [External NGS pipelines](27_bioinformatics.md) (FastQC, mapping/BAM, variant calling, SPAdes, DNA-seq/RNA-seq/ChIP-seq)

## 3) Data analysis

Added APIs: [variants and GWAS](../native/variants.md) (80), [diseases and phenotypes](../native/phenotypes.md) (45), [expression and regulation](../native/regulation.md) (65), [pathways, interactions and structures](../native/pathways.md) (64), [targets and drugs](../native/targets.md) (62), [models](../native/models.md) (14), and [sequence analysis and inheritance](../native/analysis.md) (15).

- **DNA Descriptors, Characterizations and Fingerprints**
  - [DNA descriptor](05_all_descriptors.md)
  - [Sequence representations](08_fingerprints.md)
    - [Neural Network Representation](08_fingerprints.md#neural-representations) (DATA-027: 11 DNA foundation-model representations)
  - [DNA fingerprint](08_feature_engineering.md)
  - [Advanced fingerprints and preprocessing](../index.md#fingerprints) (motifs, restriction sites, spatial GC, repeats, coding, thermodynamics, multiscale, hybrid, training-set preprocessing)
  - [Traditional DNA encodings and feature engineering](27_bioinformatics.md) (44 encodings, selection, reduction, clustering)
- **Sequence Search**
  - [Universal Search](09_search.md) (SIM-001, SIM-002, SIM-004, SIM-005, SIM-014, SIM-015)
  - [Sequence Function Search](06_patterns.md) (PAT-001, PAT-003–007, PAT-009–012)
- [Physical and chemical properties](07_physicochemical.md)
- [Double-chain thermodynamics extension](19_duplex_thermodynamics.md)
- [Secondary structure properties](20_secondary_structure.md)
- [Three-dimensional structure and mechanical properties](21_structure3d.md)
- [Conversion](18_optics_concentration.md)
- [Property prediction: deep learning and machine learning](23_deep_learning_property_prediction.md) (custom classification/regression, 37 machine-learning algorithms, reloaded inference, 54 official direct-inference tasks)
- [Custom training API](../index.md#deep-learning-property-prediction) (full fine-tuning, linear probing, machine learning, model catalogs, and optional LLM field analysis)
- [Optional bioinformatics functions](24_optional_bioinformatics.md) (variant annotation, dN/dS, and Golden Gate)
- [Gene names and intervals](26_gene_annotation.md) (genomic DNA reference alignment, existing annotation mapping, overlapping intervals, and cache reuse)
- [MSA, phylogeny, and probabilistic motifs](27_bioinformatics.md) (multiple alignment, BLAST/HMMER, NJ/UPGMA, tree comparison/consensus, PWM/PSSM)
- [Statistics and count inference](27_bioinformatics.md) (permutation tests, ordination, composition, ANCOM/ANCOM-BC, NB2 differential expression, BH correction)

## 4) Data evaluation

- [Commonly used evaluation indicators](12_evaluation.md) (EVAL-001, EVAL-002, EVAL-005–008, EVAL-016–018)
- [Similarity calculation](09_similarity_alignment.md) (SIM-010–013 and EVAL-011 reference similarity)
- [Sequence Distance and Alignment](09_alignment.md) (SIM-006–008)
- [Evaluation API extensions](../index.md#evaluation) (quality, complexity, reference libraries, memorization, distribution similarity, and scorecards)
- [Generated-set and multi-objective metrics](27_bioinformatics.md) (MMD/KID, PRDC/clipped, authenticity, property KL/hit rates, conformity, Pareto, hypervolume)
- [Community diversity and distances](27_bioinformatics.md) (alpha/beta diversity, weighted/unweighted UniFrac; distinct from generated DNA-set diversity)

## 5) Visualization

- [Visualization](13_visualization.md) (VIZ-001–008)
- [Visualization API extensions](../index.md#visualization) (similarity matrices, HTML reports, and image export)
- [Bioinformatics plots](27_bioinformatics.md) (trees, dotplots, Sanger chromatograms, qualitative gels)

## 6) Engineering and Expansion

Added [discovery, workflow and report APIs](../native/execution.md) cover 12 functions. See the [complete index of 484 added functions](../native/index.md) and [calling examples](../native/usage.md).

- [Backends, performance, and reproducibility](14_engineering.md) (backends, batches, bounded threads, caches, workflows, provenance)
- [CLI commands](../index.md#cli) (data processing, analysis, workflows, custom property training/prediction, and gene annotation)
- [LLM tools → MCP](../../llm_tools/mcp.md) (ENG-013; server setup, eight MCP entry points, schemas, tool execution, and authorization)
- [LLM tools → Python API](../../llm_tools/python_api.md) (no MCP installation required; natural-language tool execution, general-purpose LLM calls, sync/async, custom prompts, text/JSON, local files, provider presets, and extension parameters)
- [Python DAGs and interval queries](27_bioinformatics.md) (dependency checks, failure propagation, interval intersection/union/difference)
- Galaxy wrappers (`galaxy/dnakit`): convert, deduplicate, split, report; requires Galaxy.

## Bioinformatics Python APIs

- [Complete tree of 13 submodules](function_tree.md#bioinformatics-extension): alignment, phylogeny, reads, ngs, encoding, features, metrics, statistics, counts, optimization, motifs, cloning, workflows.
- [Implementation scope, dependencies, and examples](27_bioinformatics.md): development-version Python APIs, outside the default MCP catalog; native algorithms, optional Python libraries, and external programs are identified separately.
