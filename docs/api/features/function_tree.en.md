# DNAKit function tree

> Updated: 2026-09-09. Checked against public exports in the current workspace's `src/dnakit`, model/encoding catalogs, CLI, Agent catalog, and feature guides, including uncommitted implementations. The six functional areas are followed by an expanded bioinformatics Python API tree. IDs follow feature guides or requirement tracking and are not consecutive feature counts. This overview groups user capabilities rather than enumerating configuration classes, result types, compatibility aliases, or low-level helpers. Listed capabilities have interfaces and implementations; this does not imply all optional dependencies are installed or all external backends have been tested locally.

Details: [Feature module index](index.md) · [Agent/MCP tools and general-purpose LLM calls](../../agent_tools.md) · [Custom property training and prediction](23_deep_learning_property_prediction.md) · [Gene names and intervals](26_gene_annotation.md) · [Bioinformatics API](27_bioinformatics.md) · [Complete API reference](../index.md).

```
DNAKit feature overview (quick start tutorial + standalone tools and API extensions)
├── Function Tree Overview
├── Feature Module Index
├── 1. Data preparation
│   ├── (1) Core data object
│   │   ├── Unified entrance for ordinary users: dnakit.DNA(...)
│   │   │   ├── A sequence + optional ID / topology / metadata / feature
│   │   │   ├── Multiple sequences or record mapping one by one
│   │   │   └── Subscripting and slicing still return DNA
│   │   ├── 1. CORE-001 DNA sequence object
│   │   ├── 2. CORE-002 DNA recording object
│   │   ├── 3. CORE-003 DNA Dataset Object
│   │   ├── 4. CORE-004 Characteristic Object
│   │   ├── 5. CORE-005 Gap object
│   │   ├── 6. CORE-006 sequence type declaration
│   │   └── 7. CORE-007 coordinate system
│   ├── (2) File reading and writing
│   │   ├── 1. IO-001 sequence format
│   │   │   ├── read(..., mode="dna"): ordinary single/multiple reads
│   │   │   ├── read(..., mode="stream"): Streaming reading of large files
│   │   │   ├── write(...): Ordinary data or streaming data use the same write entry and exit
│   │   │   ├── Batch iteration: iter_chunks (advanced modes retain tracking ID IO-005)
│   │   │   └── FASTA/FASTQ index construction, loading, and random access: build_fasta_index / build_fastq_index / load_fasta_index / load_fastq_index / fetch
│   │   ├── 2. IO-002 Annotation formats: BED / GFF3 / AGP reading/writing and GenBank features
│   │   ├── 3. IO-003 table format
│   │   ├── 4. IO-004 compressed file
│   │   └── Metadata and result tables: filter, merge, select, validate, CSV/TSV/JSON/Parquet import/export
│   ├── (3)Download
│   │   ├── Common species names (download input assistance)
│   │   ├── 1. DBD-001 complete genome
│   │   ├── 2. DBD-002 Partial Genome Data Package
│   │   ├── 3. DBD-003 region sequence
│   │   ├── 4. DBD-004 genome annotation
│   │   ├── 5. DBD-005 gene sequence
│   │   ├── 6. DBD-006 protein sequence
│   │   ├── 7. DBD-007 assembly information
│   │   ├── 8. DBD-008 classification information
│   │   ├── 9. DBD-009 Gene Metadata
│   │   ├── 10. DBD-010 virus packet
│   │   ├── 11. DBD-011 raw sequencing data
│   │   ├── 12. DBD-012 comparison data
│   │   ├── 13. DBD-013 sequencing metadata
│   │   ├── 14. DBD-014 mutation data
│   │   ├── 15. DBD-015 ClinVar data
│   │   ├── 16. DBD-016 Expression Data
│   │   ├── 17. DBD-017 control data
│   │   ├── 18. DBD-018 Repeated Sequence Data
│   │   ├── 19. DBD-019 Conservative Data
│   │   ├── 20. DBD-020 Multi-species comparison
│   │   ├── 21. DBD-021 coordinate conversion file
│   │   ├── General file download, verification, and provenance manifests: download_file
│   │   └── External index construction: build_index (makeblastdb / BWA / Bowtie2 / minimap2 / samtools faidx, dict)
│   ├── (4) Data query
│   │   ├── 1. Basic queries: genes, sequences, accessions, assemblies, taxonomy, projects, samples, literature, database versions
│   │   ├── 2. Coordinates, transcripts, and annotations: nearby genes, ID/coordinate conversion, Ensembl lookup
│   │   ├── 3. Reads, expression, regulation, variants, viruses, homology/orthology, comparative genomics
│   │   ├── Dedicated ENA / ENCODE / UCSC queries, tracks, and regional sequences
│   │   └── Remote BLAST: submit_blast / blast_status / wait_for_blast / blast_results; identify / novelty queries
│   └── (5) Legality check
│       └── 1. Unified sequence, record and data set legality check
├── 2. Data processing
│   ├── (1)DNA sequence standardization
│   │   ├── 1. STD-001 character standardization
│   │   ├── Uppercase N runs to known-length Gaps: normalize_gaps
│   │   └── AGP sequence assembly with Gaps and provenance: sequence_from_agp
│   ├── (2) Operation within sequence
│   │   ├── Unified suffix-free entry: insert/delete/substitute/mask/trim/reverse_complement/rotate/canonical_origin
│   │   ├── OPS-001 sequence completion
│   │   │   ├── 1. OPS-001.1 reverse order
│   │   │   ├── 2. OPS-001.2 Complementary
│   │   │   └── 3. OPS-001.3 Reverse complementation
│   │   ├── OPS-002 Transcription and Translation
│   │   │   ├── 1. OPS-002.1 Transcription
│   │   │   └── 2. OPS-002.2 Translation
│   │   ├── OPS-003 subsequence extraction
│   │   ├── OPS-004 Sequence Editing
│   │   ├── OPS-005 sequence generation
│   │   │   ├── 1. OPS-005.1 Mutation Generation
│   │   │   ├── 2. OPS-005.2 Insert generation
│   │   │   ├── 3. OPS-005.3 Delete generation
│   │   │   ├── 4. OPS-005.4 Evolutionary algorithm-like operation
│   │   │   ├── 5. OPS-005.5 Fragment rearrangement
│   │   │   ├── 6. OPS-005.6 k-mer randomly scrambled
│   │   │   └── 7. OPS-005.7 Multiple sequence recombination/crossover
│   │   ├── OPS-006 Sequence splicing
│   │   │   ├── 1. Ordinary sequence splicing
│   │   │   └── 2. Remove overlapping splicing
│   │   ├── OPS-007
│   │   │   ├── 1. Trimming
│   │   │   └── 2. Masking
│   │   ├── OPS-008 Circular Sequence Operation
│   │   └── OPS-010 Sequence segmentation
│   │       ├── 1. Fixed length, no overlap
│   │       ├── 2. Fixed length, overlapping sliding windows
│   │       ├── 3. Random interval
│   │       ├── 4. Multi-scale
│   │       └── 5. Short first and then long
│   ├── (3) Sequence deduplication
│   │   ├── 1. DATA-001 · Standard deduplication
│   │   ├── 2. DATA-002 · Reverse complementary deduplication
│   │   ├── 3. DATA-003 · Ring equivalent deduplication
│   │   ├── 4. DATA-004 · IUPAC-aware deduplication
│   │   └── 5. DATA-005 · Approximate deduplication
│   ├── (4) Sequence clustering
│   │   ├── 1. DATA-007 · Identity clustering
│   │   ├── 2. DATA-008 · k-mer clustering
│   │   ├── 3. DATA-009 · Fingerprint clustering
│   │   ├── 4. DATA-010 · Hierarchical clustering
│   │   ├── 5. DATA-011 · Representative sequence selection
│   │   └── 6. DATA-027 · Neural network clustering
│   ├── (5) Data set organization and division
│   │   ├── 1. DATA-012 · Random and stable hash partitioning
│   │   ├── 2. DATA-013 · Stratified random division
│   │   ├── 3. DATA-014 · Similarity division
│   │   ├── 4. DATA-015 · Cluster split
│   │   ├── 5. DATA-016 · Species classification
│   │   ├── 6. DATA-017 · Chromosome division
│   │   ├── 7. DATA-018 · Individual Division
│   │   ├── 8. Divide by custom label
│   │   ├── 9. DATA-023 · Leak detection
│   │   ├── Temporal and joint constraint splitting: temporal_split / joint_split
│   │   ├── Exclude records by species, chromosome, or custom metadata
│   │   └── Split quality evaluation: evaluate_split_quality
│   └── (6) Molecular biology simulation and design (API extension)
│       ├── Restriction digestion, end classification, compatibility checks, and ligation
│       ├── Primer properties, primer matching, and PCR simulation
│       ├── Primer design requests and the Primer3 CLI design adapter
│       ├── Sequence-level assembly simulation: simulate_assembly
│       ├── CRISPR candidate and off-target scanning: scan_crispr_candidates
│       ├── Codon optimization, rule-constrained optimization, and mutation libraries
│       └── Extended constraint optimization, duplex cloning, and Sanger validation: bioinformatics tree below
├── 3. Data analysis
│   ├──  (1) DNA descriptor + characterization + fingerprint
│   │   ├── DNA descriptor
│   │   │   ├── 1. STD-005 Fuzzy base statistics
│   │   │   ├── 2. DESC-001 length characteristics
│   │   │   ├── 3. DESC-002 base composition
│   │   │   ├── 4. DESC-003 GC/AT features
│   │   │   ├── 5. DESC-004 base skew
│   │   │   ├── 6. DESC-005 CpG Characteristics
│   │   │   ├── 7. DESC-006 k-mer statistics
│   │   │   ├── 8. DESC-007 Sequence Entropy
│   │   │   ├── 9. DESC-008 sequence complexity
│   │   │   ├── 10. DESC-009 Homopolymer
│   │   │   ├── 11. DESC-010 Repeat Ratio
│   │   │   ├── 12. DESC-011 window descriptor
│   │   │   ├── 13. DESC-012 Codon Statistics
│   │   │   ├── 14. All descriptors: all_descriptors (240 output fields in a fixed schema)
│   │   │   └── Descriptor schema discovery and provenance-aware dinucleotide property table loading
│   │   ├── DNA characterization
│   │   │   ├── 1. FP-001 integer encoding
│   │   │   ├── 2. FP-002 One-hot encoding
│   │   │   ├── 3. FP-003 k-mer features
│   │   │   ├── 4. FP-005 k-mer Sketch (MinHash/FracMinHash)
│   │   │   └── 5. Neural network representation
│   │   │       ├── 1. AlphaGenome
│   │   │       ├── 2. Caduceus
│   │   │       ├── 3.DNABERT-2
│   │   │       ├── 4. Enformer
│   │   │       ├── 5. Evo 2
│   │   │       ├── 6. GENERator
│   │   │       ├── 7.GROVER
│   │   │       ├── 8. HyenaDNA
│   │   │       ├── 9. JanusDNA
│   │   │       ├── 10. LucaOne (default)
│   │   │       └── 11. Nucleotide Transformer v2
│   │   └── DNA fingerprint
│   │       ├── 1. Hashed k-mer bit fingerprint
│   │       ├── 2. Panel existence fingerprint
│   │       ├── Motifs / restriction sites: motif_fingerprint / restriction_fingerprint
│   │       ├── Spatial GC / repeats / coding: gc_spatial_fingerprint / repeat_fingerprint / coding_fingerprint
│   │       ├── Thermodynamic / multiscale / hybrid: thermodynamic_fingerprint / multiscale_fingerprint / hybrid_fingerprint
│   │       ├── Feature preprocessing: fit_preprocessor / transform (missing values, scaling, normalization, low-variance filtering)
│   │       └── 44 traditional DNA encodings, selection, reduction, and clustering: encoding / features below
│   ├── (2) Sequence search
│   │   ├── Universal search
│   │   │   ├── 1. SIM-001 · Exact search
│   │   │   ├── 2. SIM-002 · Subsequence search
│   │   │   ├── 3. SIM-004 · Approximate matching
│   │   │   ├── 4. SIM-005 · Reverse complementary search
│   │   │   ├── 5. SIM-014 · Nearest neighbor search
│   │   │   ├── 6. SIM-015 · Database Index
│   │   │   └── Sketch index construction and persistence: build_sketch_index / save_sketch_index / load_sketch_index
│   │   └── Sequence function search
│   │       ├── 1. PAT-001 / SIM-003 motif search
│   │       ├── 2. PAT-003 start and stop codons
│   │       ├── 3. PAT-004 promoter motif
│   │       ├── 4. PAT-005 TF motif
│   │       ├── 5. PAT-006 restriction enzyme site
│   │       ├── 6. PAT-007 CRISPR PAM
│   │       ├── 7. PAT-009 palindrome sequence
│   │       ├── 8. PAT-010 Repeat upside down
│   │       ├── 9. PAT-011 Tandem Repeat
│   │       ├── 10. PAT-012 Microsatellite
│   │       ├── Open reading frames: scan_orfs
│   │       └── CpG islands and low-complexity regions: find_cpg_islands / find_low_complexity_regions
│   ├── (3)Physical and chemical properties
│   │   ├── 1. THERMO-001 molecular weight
│   │   ├── 2. THERMO-014 260 nm extinction coefficient
│   │   ├── 3. THERMO-002 melting temperature Tm
│   │   ├── 4. THERMO-003 salt concentration correction
│   │   ├── 5. THERMO-012 local melting characteristics
│   │   ├── 6. EVAL-013 Synthesis Risk
│   │   └── Thermodynamic properties
│   │       ├── 1. THERMO-004 thermodynamic parameters
│   │       ├── 2. THERMO-005 Nearest-neighbor
│   │       ├── 3. THERMO-006 Duplex stability
│   │       ├── 4. THERMO-007 base stacking
│   │       ├── 5. THERMO-008 Hairpin
│   │       ├── 6. THERMO-009 Self-dimer
│   │       └── 7. THERMO-010 Heterodimer
│   ├── (4) Double-chain thermodynamic expansion
│   │   ├── 1. ΔH, ΔS, ΔG and Tm of perfectly complementary duplexes
│   │   ├── 2. Unify double-strand stability results
│   │   ├── 3. Contribution of adjacent base pairs to steps
│   │   ├── 4. Conditions and Na⁺/K⁺ unit price salt
│   │   ├── 5. Ka, Kd and duplex ratio
│   │   ├── 6. Theoretical melting curve
│   │   ├── 7. 5′/3′ end stability
│   │   ├── 8. DMSO and formamide experience correction
│   │   └── 9. Mg²⁺, dNTP, mismatch and dangling end for Primer3 CLI
│   ├── (5) Secondary structure properties
│   │   ├── 1. Dot-bracket structure analysis
│   │   ├── 2. Pairing probability and window accessibility
│   │   ├── 3. Collection defects of target structure
│   │   ├── 4. Thermodynamic probability of target structure
│   │   ├── 5. NUPACK Passive Availability Check
│   │   ├── 6. NUPACK single-complex ensemble analysis
│   │   └── 7. NUPACK tube multi-complex equilibrium
│   ├── (6) Three-dimensional structure and mechanical properties
│   │   ├── 1. Read a single PDB model
│   │   ├── 2. Read PDB multi-model collection
│   │   ├── 3. Explicit coordinate geometry analysis
│   │   ├── 4. NMR/Multi-model RMSF Flexibility
│   │   ├── 5. 3DNA bp_step.par standard parameter analysis
│   │   └── 6. DSSR JSON summary parsing
│   ├── (7) Conversion
│   │   ├── 1. Single chain molar extinction coefficient at 260 nm
│   │   ├── 2. Single-chain/double-chain theoretical optical properties
│   │   ├── 3. 1 nmol and mass corresponding to OD260
│   │   ├── 4. A260 to molar concentration, mass concentration and total amount
│   │   ├── 5. Concentration, quantity and mass interchange of substances
│   │   └── 6. Explicit correction of dyes and modifying groups
│   ├── (8) Property prediction: deep learning and machine learning
│   │   ├── User-defined classification and regression
│   │   │   ├── Multi-format inspection, automatic/explicit field mapping, and DataReader extensions
│   │   │   ├── Binary, multiclass, multilabel, multi-output classification; single/multi-target regression
│   │   │   ├── Full fine-tuning (full) and frozen-representation linear probing (linear_probe)
│   │   │   ├── Machine learning (machine_learning): 37 algorithms and sequence/embedding/numeric features
│   │   │   │   ├── Classification and regression: ridge / sgd / svm / linear_svm / nu_svm / knn / decision_tree / extra_tree
│   │   │   │   ├── Classification and regression: random_forest / extra_trees / gradient_boosting / hist_gradient_boosting / adaboost / bagging
│   │   │   │   ├── Classification and regression: gaussian_process / mlp / dummy / xgboost / lightgbm / catboost
│   │   │   │   ├── Classification only: logistic_regression / perceptron / gaussian_nb / multinomial_nb / bernoulli_nb / complement_nb / lda / qda
│   │   │   │   └── Regression only: linear_regression / lasso / elastic_net / bayesian_ridge / ard / huber / ransac / theil_sen / kernel_ridge
│   │   │   ├── Data splits, training parameters, classification/regression metrics, and artifact export
│   │   │   ├── Training and reloaded inference: train_custom_property_model / predict_custom_properties
│   │   │   └── Optional LLM field analysis: provider presets, structure/bounded samples, declarative mapping
│   │   ├── Official tasks 1–11. Genomic output modalities and SNV track differences
│   │   │   ├── RNA-seq, CAGE, PRO-cap, ATAC-seq, and DNase-seq
│   │   │   ├── Histone modifications and transcription factor binding
│   │   │   └── Splice sites, splice usage, splice junction counts, and 3D DNA contact maps
│   │   ├── Official tasks 12–13. Human/mouse regulatory tracks and SNV differences
│   │   ├── Official task 14. Fourteen-class single-nucleotide segmentation
│   │   │   ├── 14.1) Protein-coding gene
│   │   │   ├── 14.2) Long noncoding RNA (lncRNA)
│   │   │   ├── 14.3) Exon
│   │   │   ├── 14.4) Intron
│   │   │   ├── 14.5) Splice donor
│   │   │   ├── 14.6) Splice acceptor
│   │   │   ├── 14.7) 5′ UTR
│   │   │   ├── 14.8) 3′ UTR
│   │   │   ├── 14.9) CTCF-bound region
│   │   │   ├── 14.10) PolyA signal
│   │   │   ├── 14.11) Tissue-specific enhancer
│   │   │   ├── 14.12) Tissue-invariant enhancer
│   │   │   ├── 14.13) Tissue-specific promoter
│   │   │   └── 14.14) Tissue-invariant promoter
│   │   ├── Official tasks 15–16. Long-context zero-shot variant score and exon probability
│   │   ├── Official task 17. Allele-conditional-probability variant score
│   │   ├── Official tasks 18–27. Ten trained nucleic-acid/protein downstream tasks
│   │   │   ├── Nucleic-acid/protein central-dogma relation, SupKTax, GenusTax, and species classification
│   │   │   ├── Prokaryotic protein localization, protein stability, and ncRNA families
│   │   │   └── Influenza A antigenic relationships, protein/protein and ncRNA/protein interactions
│   │   ├── Official tasks 28–45. Eighteen fully fine-tuned NT Revised classifiers
│   │   │   ├── H2A.Z, H3K27ac, H3K27me3, H3K36me3, H3K4me1/2/3, H3K9ac/me3, and H4K20me1
│   │   │   ├── Enhancer detection/types; promoters, non-TATA promoters, and TATA promoters
│   │   │   └── Splice acceptors, splice-site types, and splice donors
│   │   └── Official tasks 46–54. Nine fully fine-tuned Genomic Benchmarks classifiers
│   │       ├── Coding/intergenic regions, human/worm, Drosophila enhancers, and mouse enhancers
│   │       ├── Human Cohn/Ensembl enhancers and human regulatory-element types
│   │       └── Human non-TATA promoters and human open-chromatin regions
│   ├── (9) Deep-learning task checkpoint download
│   │   ├── Eighteen NT Revised and nine Genomic Benchmarks task weights
│   │   └── Default/custom directories, exact paths, loading validation, and prediction output
│   ├── (10) Optional bioinformatics functions
│   │   ├── 1. Ensembl VEP and variant identifier recoding
│   │   ├── 2. ClinVar, dbSNP, and gnomAD annotation
│   │   ├── 3. Nei–Gojobori dN/dS
│   │   ├── 4. Golden Gate design and reaction assembly
│   │   └── MSA, phylogeny, NGS, statistics, and molecular design: bioinformatics tree below
│   └── (11) Gene names and intervals
│       ├── Genomic DNA + species → annotate_genes / GeneAnnotationConfig
│       ├── Default RefSeq assembly resolution, explicit assembly, or local FASTA + GFF3
│       ├── minimap2 reference alignment and existing annotation mapping (not de novo prediction)
│       ├── Gene names, IDs, strands, overlapping genes, and gene/intergenic/unknown segments
│       ├── Python: 0-based half-open intervals; JSON/CLI: 1-based closed intervals
│       └── Pinned assembly/annotation cache, checksums, offline reuse, and provenance
├── 4. Data evaluation
│   ├── (1) Commonly used evaluation indicators
│   │   ├── 1. EVAL-001 Validity
│   │   ├── 2. EVAL-005 Uniqueness
│   │   ├── 3. EVAL-006 Diversity (similarity / mean pairwise Levenshtein)
│   │   ├── 4. EVAL-008 Novelty (similarity / mean nearest-reference Levenshtein)
│   │   ├── 5. EVAL-016 Fréchet DNA distance
│   │   ├── 6. EVAL-017 Frag
│   │   ├── 7. EVAL-018 SNN
│   │   ├── 8. EVAL-002 Ambiguity
│   │   ├── 9. EVAL-007 Redundancy
│   │   ├── Complexity and quality evaluation: evaluate_complexity / evaluate_quality
│   │   ├── Reference libraries, nearest-reference lookup, and memorization evaluation
│   │   ├── Composite scorecard: evaluate_scorecard
│   │   └── MMD/KID, PRDC, authenticity, property distributions, and multi-objective evaluation: metrics below
│   ├── (2) Similarity calculation
│   │   ├── 1. SIM-010 · k-mer similarity
│   │   ├── 2. SIM-011 · Fingerprint similarity
│   │   ├── 3. SIM-012 · Sketch similarity
│   │   ├── 4. SIM-013 · Dashing sketch similarity
│   │   ├── 5. EVAL-011 Reference similarity
│   │   ├── Descriptor distribution similarity: evaluate_distribution_similarity
│   │   └── Vectors and matrices: kmer_vector_similarity / similarity_matrix / pairwise_matrix
│   └── (3) Sequence distance and alignment
│       ├── 1. Sequence distance
│       │   ├── 1.1 SIM-006 · Hamming distance
│       │   └── 1.2 SIM-007 · Edit distance
│       └── 2. Pairwise sequence alignment (SIM-008)
│           ├── 2.1 Global comparison
│           ├── 2.2 Local comparison
│           └── 2.3 Semi-global comparison
├── 5. Visualization
│   ├── 1. VIZ-001 sequence text diagram
│   ├── 2. VIZ-002 position highlighting
│   ├── 3. VIZ-003 style control
│   ├── 4. VIZ-004 Gap display
│   ├── 5. VIZ-005 linear sequence diagram
│   ├── 6. VIZ-006 circular DNA diagram
│   ├── 7. VIZ-007 Custom character representation
│   ├── 8. VIZ-008 Alignment diagram
│   ├── Similarity matrix heatmap: plot_similarity_matrix
│   ├── Self-contained, filterable, expandable HTML reports: build_html_report / save_html_report
│   ├── Phylogenetic trees, dotplots, Sanger chromatograms, and qualitative gels: phylogeny / motifs / cloning below
│   └── Image export: save_svg / save_image (PNG, SVG, JPG; TIFF/PDF also supported)
└── 6. Engineering and Expansion
    ├── 1. ENG-001 unified backend interface
    ├── 2. ENG-002 Native and external tags
    ├── 3. ENG-003 Python API
    ├── 4. ENG-004 CLI
    │   ├── Environment and backends: info / backends
    │   ├── Data preparation: download-genome / convert / normalize / validate
    │   ├── Analysis and comparison: describe / fingerprint / search / orfs / compare
    │   ├── Datasets and reports: deduplicate / split / report / workflow
    │   ├── Custom properties: property-models / inspect-property-data / train-property / predict-property
    │   └── Gene annotation: annotate-genes
    ├── 5. ENG-005 Configuration workflows: load_workflow / run_workflow; artifact manifest I/O, verification, reuse
    ├── 6. ENG-006 Batch computation: iter_batch / run_batch; progress, cancellation, per-record errors
    ├── 7. ENG-007 Parallel computation: bounded threads and stable output order
    ├── 8. ENG-008 Chunking and Streaming
    ├── 9. ENG-009 cache
    ├── 10. ENG-010 Random Seed
    ├── 11. ENG-011 version tracking
    ├── 12. ENG-012 Errors and Warnings
    ├── 13. ENG-013 Agent and MCP
    │   ├── Local stdio server: dnakit-mcp; optional dependency dnakit[agent]
    │   ├── MCP client integration (Codex, Claude, Cursor, VS Code, etc.)
    │   ├── Catalog statistics: dnakit_catalog_stats
    │   ├── Categories: list_dnakit_categories
    │   ├── Paginated catalog: list_dnakit_tools
    │   ├── Capability search: search_dnakit_tools
    │   ├── Parameter schema and execution attributes: describe_dnakit_tool
    │   ├── JSON-argument execution: call_dnakit_tool
    │   ├── Python registry and server construction: default_tool_registry / create_server
    │   └── Bounded JSON output, provenance, Python-only markers, and allow_side_effects authorization
    ├── 14. ENG-014 consistency verification
    ├── 16. ENG-016 Documentation and Tutorials
    ├── 17. ENG-017 Galaxy wrappers: convert / deduplicate / split / report (requires Galaxy)
    ├── Python DAG and interval-set workflows: workflows below
    ├── General-purpose LLM calls (standalone Python API: dnakit.llm)
    │   ├── Configuration and results: LLMConfig / LLMResult
    │   ├── Sync/async: call_llm / acall_llm; custom prompts and current-step data
    │   ├── Text/JSON input, text/JSON-object output, and token usage
    │   ├── Local text files: file_path; local preview with read_llm_file
    │   ├── Encoding and file/request size limits; explicit file input sends the complete text
    │   ├── Provider presets: DeepSeek / Gemini / Qwen / Kimi / OpenRouter / OpenAI / custom
    │   ├── Preset discovery: available_llm_providers / get_llm_provider
    │   └── Custom compatible endpoint, model, credentials, and extra_body parameters
    └── Development validation: unit tests and consistency checks
```

## Bioinformatics Python API extension {#bioinformatics-extension}

[Implementation boundaries and examples](27_bioinformatics.md). These 13 submodules are development-version Python APIs. They are not registered in the default MCP catalog and have no dedicated DNAKit CLI commands. Import submodules explicitly. Native domain algorithms may require NumPy/SciPy; codecs, general machine learning, and specialist assembly use the optional dependencies identified below. External programs require separate installation.

```text
dnakit.bioinformatics
├── 1. alignment: multiple alignment, BLAST, and profile HMMs
│   ├── Native MSA, width/ID validation, and p-distance: MultipleAlignment.from_sequences / distances
│   ├── MSA codecs: read_alignment / write_alignment (Biopython)
│   ├── Multiple alignment: align_sequences (external MAFFT / MUSCLE 5 / Clustal Omega)
│   ├── Local nucleotide BLAST: build_blast_database / local_blast (external BLAST+)
│   ├── BLAST XML/tab and protein HMMER2/3 parsing: parse_search_results (Biopython)
│   └── Profile construction/search and nucleotide table parsing: build_profile / search_profile / parse_nhmmer_table (hmmbuild / nhmmer / hmmsearch)
├── 2. phylogeny: phylogenetic analysis
│   ├── Native trees and construction: TreeNode / build_tree (NJ / cluster-size-weighted UPGMA)
│   ├── Tree codecs: read_tree / write_tree (Newick / Nexus / phyloXML / NeXML; Biopython)
│   ├── Native tree comparison, consensus, and pruning: robinson_foulds / consensus_tree / prune_tree
│   ├── Native community phylogenetic distances: unifrac (weighted / unweighted; rooted trees)
│   └── Tree plots: plot_tree (Matplotlib)
├── 3. reads: sequencing records, indexes, and regions
│   ├── FASTA / multiline FASTQ / gzip and Phred quality: Read / read_fastx / write_fastx
│   ├── Pairing and damaged four-line FASTQ recovery: pair_reads / repair_fastq (recovery/discard statistics)
│   ├── Sorting, shuffling, and common records: arrange_reads / common_reads
│   ├── Regex renaming and order-sensitive, ID-independent sequence digests: rename_reads / sequence_digest
│   ├── Intervals and BED/GTF parsing: Region / read_regions
│   ├── FASTA indexing and regional access: index_reference / IndexedReference.fetch (pysam)
│   └── BAM summaries, coverage, and reference-coordinate consensus: alignment_summary / coverage_consensus (pysam; no indels in consensus)
├── 4. ngs: external program adapters and pipelines
│   ├── Raw read quality control: quality_control (FastQC)
│   ├── Single/paired-end mapping and private indexes: map_reads (BWA / Bowtie2)
│   ├── BAM sorting, mate repair, optional duplicate marking, indexing, and statistics: prepare_bam (samtools)
│   ├── Variant calling, reference normalization, and VCF indexing: call_variants (bcftools)
│   ├── Short-read/hybrid de novo assembly: assemble_reads (SPAdes; optional Nanopore)
│   ├── DNA-seq: dna_seq_workflow (QC → mapping → BAM → variants)
│   ├── Single-sample RNA-seq: rna_seq_workflow (HISAT2 / StringTie / featureCounts; differential analysis under counts)
│   └── ChIP-seq from prepared BAMs: chip_seq_workflow (MACS3 peaks + samtools coverage)
├── 5. encoding: 44 native DNA encodings (encode_dna / DNA_ENCODINGS)
│   ├── Word frequencies and approximate/noncontiguous words: Kmer / RCKmer / Mismatch / Subsequence
│   ├── Bases and positions: NAC / ANF / ENAC / binary / PS2 / PS3 / PS4
│   ├── Dinucleotides: CKSNAP / ASDC / DBE / LPDF
│   ├── Chemical/EIIP/mutual information: NCP / EIIP / PseEIIP / MMI
│   ├── Z-curve: Z_curve_9bit / Z_curve_12bit / Z_curve_36bit / Z_curve_48bit / Z_curve_144bit
│   ├── Positional physicochemical properties: DPCP / TPCP
│   ├── Auto/cross covariance: DAC / DCC / DACC / TAC / TCC / TACC
│   ├── Autocorrelation: NMBroto / Moran / Geary
│   ├── Pseudo nucleotide composition: PseDNC / PseKNC / PCPseDNC / PCPseTNC / SCPseDNC / SCPseTNC
│   └── Supervised encodings: PSTNPss / PSTNPds / KNN (explicit training_ids; physicochemical encodings require explicit property tables)
├── 6. features: feature processing (scikit-learn; native MCL)
│   ├── Fitting and held-out reuse: analyze_features / FeatureAnalysis.transform
│   ├── Scaling: standard / minmax / l1 / l2 / none
│   ├── Reduction: pca / tsne / lda / latent_dirichlet (t-SNE has no held-out transform)
│   ├── Selection: chi2 / f_classif / f_regression / mutual_info_classif / mutual_info_regression
│   ├── Selection: information_gain / pearson / variance; none for preprocessing only
│   └── Feature clustering: cluster_features
│       ├── kmeans / minibatch_kmeans / gmm / spectral / affinity_propagation
│       └── mean_shift / dbscan / agglomerative / birch / mcl
├── 7. metrics: generated sets, property distributions, and multiple objectives
│   ├── Custom scalar property matrices: sequence_properties
│   ├── Gaussian MMD² / polynomial KID: kernel_distance
│   ├── Precision / Recall / Density / Coverage: manifold_metrics
│   ├── Calibrated clipped density / coverage: clipped_manifold_metrics
│   ├── Authenticity relative to generator training data: authenticity
│   ├── Pareto layers and exact hypervolume: pareto_ranks / hypervolume
│   ├── Property interval hit rates and histogram KL: property_hit_rate / property_kl
│   └── KDE distribution conformity and standard error: conformity_score
├── 8. statistics: distance tests, ordination, and composition
│   ├── Distance validation and permutation tests: distance_matrix / distance_test (PERMANOVA / ANOSIM / centroid PERMDISP)
│   ├── Matrix correlations and environmental subset search: mantel / bioenv
│   ├── Principal coordinates and constrained ordination: pcoa / ordination (CA / CCA / RDA)
│   ├── Compositional transforms and inverses: composition_transform (CLR / ILR / ALR)
│   ├── Native compositional differential analysis: ancom
│   ├── Alpha diversity: alpha_diversity (observed / Shannon / Simpson / inverse Simpson)
│   └── Beta distances: beta_diversity (Bray-Curtis / Jaccard / Euclidean / cityblock / cosine)
├── 9. counts: count-data inference
│   ├── Original ANCOM-BC: ancom_bc (OLS / sandwich variance / EM sampling-bias correction)
│   ├── Differential expression: differential_expression (unshrunk NB2 GLM / size factors / explicit contrast / Wald)
│   └── Multiple testing: adjust_pvalues (Benjamini–Hochberg)
├── 10. optimization: native constraint optimization
│   ├── Rules and custom evaluators: specification / Specification / SpecEvaluation
│   ├── GC, pattern avoidance/counts, sequence/translation/length, locked positions: gc / avoid_pattern / pattern_count / sequence / translation / length / avoid_changes
│   ├── Unique k-mers, stop/forbidden/rare codons: unique_kmers / avoid_stop / avoid_codons / avoid_rare_codons
│   ├── Exact-complement structure/match rules: hairpin / heterodimer / avoid_matches / regions_compatibility
│   ├── Tm and codon objectives: tm / cai / codon_usage / harmonize
│   ├── External homology constraint: blast_match_specification (BLAST)
│   ├── Hard constraints first, weighted objectives, synonymous edits, and circular sequences: optimize_sequence
│   ├── GenBank qualifier rule import: specifications_from_record
│   └── Edits, failure locations, search budgets, and JSON reports: OptimizationResult.write_report
├── 11. motifs: probabilistic motifs and dotplots
│   ├── Native counts / PWM / PSSM, backgrounds, and degenerate consensus: MotifModel / build_motif
│   ├── Conservative FPR thresholds and both-strand/circular scans: fpr_threshold / scan_motif_model
│   ├── JASPAR / MEME / TRANSFAC reading: read_motifs (Biopython)
│   └── Forward/reverse exact k-mer dotplots: dotplot (Matplotlib)
├── 12. cloning: explicit duplexes and cloning
│   ├── Duplexes, offsets, sticky/blunt ends, and topology: Duplex.from_sequence / ends / cut / ligate / reverse_complement
│   ├── Executable cloning steps and intermediates: CloningPlan.add / apply / steps
│   ├── Gibson / homologous overlap assembly: assemble_duplexes (pydna)
│   ├── ABI reading, Mott trimming, and chromatograms: read_sanger / plot_chromatogram (Biopython / Matplotlib)
│   ├── Both-strand local alignment and reference-coordinate consensus: sanger_consensus
│   ├── SnapGene / GenBank constructs: read_constructs (Biopython; optional permissive pydna parsing)
│   └── Qualitative gel simulation: simulate_gel (Matplotlib)
└── 13. workflows: native Python workflows
    ├── DAG definitions, dependency/cycle checks, and failure propagation: Step / run_workflow
    └── Chromosome/strand-aware interval sets: query_regions (intersection / union / difference)
```

## Audit scope and counts

This audit covered public module exports, public functions/class methods in each `bioinformatics` submodule, `tools/catalog.py`, CLI commands, and Galaxy wrappers. The original six-area organization remains; the extension tree expands related capabilities without counting them twice.

| Code catalog | Verified count | Meaning |
| --- | --- | --- |
| `descriptors.DESCRIPTOR_NAMES_V1` | 240 | Fixed-schema output fields, including statistics and property-table fields; not 240 independent APIs |
| `bioinformatics.encoding.DNA_ENCODINGS` | 44 | Encoding names; output dimensions depend on parameters/sequences |
| `representations.available_embedding_models()` | 11 | Representation models; full fine-tuning is not supported for every model |
| `predictions.available_machine_learning_models()` | 37 | Algorithm entries with separate classification/regression support |
| `predictions.available_prediction_tasks()` | 54 | Registered direct prediction tasks, excluding user-trained properties |

Formulas, input requirements, and limits are defined in the [bioinformatics API guide](27_bioinformatics.md) and code docstrings. For example, gene naming uses reference alignment; NB2 differential expression omits the full DESeq2/edgeR shrinkage pipelines; constraint optimization uses a bounded search. Entries reflect these implementation boundaries.
