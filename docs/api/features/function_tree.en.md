# DNAKit function tree

## Core feature tree

```
DNAKit feature overview
├── Function Tree Overview
├── 1. Data preparation
│   ├── 1) Core data object
│   │   ├── Unified entrance for ordinary users: dnakit.DNA(...)
│   │   │   ├── A sequence + optional ID / topology / metadata / feature
│   │   │   ├── Multiple sequences or record mapping one by one
│   │   │   └── Subscripting and slicing still return DNA
│   │   ├── 1. DNA sequence object
│   │   ├── 2. DNA recording object
│   │   ├── 3. DNA Dataset Object
│   │   ├── 4. Characteristic Object
│   │   ├── 5. Gap object
│   │   ├── 6. sequence type declaration
│   │   └── 7. coordinate system
│   ├── 2) File reading and writing
│   │   ├── 1. sequence format
│   │   │   ├── read(..., mode="dna"): ordinary single/multiple reads
│   │   │   ├── read(..., mode="stream"): Streaming reading of large files
│   │   │   ├── write(...): Ordinary data or streaming data use the same write entry and exit
│   │   │   ├── Batch iteration: iter_chunks
│   │   │   └── FASTA/FASTQ index construction, loading, and random access: build_fasta_index / build_fastq_index / load_fasta_index / load_fastq_index / fetch
│   │   ├── 2. Annotation formats: BED / GFF3 / AGP reading/writing and GenBank features
│   │   ├── 3. table format
│   │   ├── 4. compressed file
│   │   ├── Metadata and result tables: filter, merge, select, validate, CSV/TSV/JSON/Parquet import/export
│   │   ├── MSA codecs: read_alignment / write_alignment (Biopython)
│   │   ├── Tree codecs: read_tree / write_tree (Newick / Nexus / phyloXML / NeXML; Biopython)
│   │   ├── FASTA / multiline FASTQ / gzip and Phred quality: Read / read_fastx / write_fastx
│   │   ├── JASPAR / MEME / TRANSFAC reading: read_motifs (Biopython)
│   │   ├── SnapGene / GenBank constructs: read_constructs (Biopython; optional permissive pydna parsing)
│   │   └── ABI file reading and Mott trimming: read_sanger
│   ├── 3) Download
│   │   ├── Common species names (download input assistance)
│   │   ├── 1. complete genome
│   │   ├── 2. Partial Genome Data Package
│   │   ├── 3. region sequence
│   │   ├── 4. genome annotation
│   │   ├── 5. gene sequence
│   │   ├── 6. protein sequence
│   │   ├── 7. assembly information
│   │   ├── 8. classification information
│   │   ├── 9. Gene Metadata
│   │   ├── 10. virus packet
│   │   ├── 11. raw sequencing data
│   │   ├── 12. comparison data
│   │   ├── 13. sequencing metadata
│   │   ├── 14. mutation data
│   │   ├── 15. ClinVar data
│   │   ├── 16. Expression Data
│   │   ├── 17. control data
│   │   ├── 18. Repeated Sequence Data
│   │   ├── 19. Conservative Data
│   │   ├── 20. Multi-species comparison
│   │   ├── 21. coordinate conversion file
│   │   ├── NCBI/ENA accession sequence downloads and ENCODE batch downloads
│   │   ├── General URL/file downloads, checksums, and manifests: download_url / download_file / dataset
│   │   └── External index construction: build_index (makeblastdb / BWA / Bowtie2 / minimap2 / samtools faidx, dict)
│   ├── 4) Data query
│   │   ├── 1. Discover and call database APIs
│   │   ├── 2. Sequence, assembly and sequencing queries
│   │   ├── 3. Disease, phenotype and ontology queries
│   │   ├── 4. Expression, regulation and cancer queries
│   │   ├── 5. Organism, taxonomy and microbiome queries
│   │   ├── 6. Pathway, interaction, protein and structure queries
│   │   ├── 7. Target, drug and association queries
│   │   ├── 8. Literature, preprint and clinical study queries
│   │   ├── 9. ENCODE experiment and file queries
│   │   └── 10. Genetic Variation and Population Associations
│   │       ├── 1) Variant Records and Population Frequencies
│   │       ├── 2) Structural Variant Queries
│   │       ├── 3) ClinVar Clinical Significance
│   │       ├── 4) Gene–Disease Evidence, Dosage and Actionability
│   │       ├── 5) Cancer Variants and Treatment Evidence
│   │       ├── 6) GWAS Studies and Associations
│   │       ├── 7) Cohort Phenotype Associations and Fine-mapping
│   │       └── 8) Populations, Samples and Linkage Disequilibrium
│   └── 5) Legality check
│       └── 1. Unified sequence, record and data set legality check
├── 2. Data processing
│   ├── 1) DNA sequence standardization
│   │   ├── 1. character standardization
│   │   ├── Uppercase N runs to known-length Gaps: normalize_gaps
│   │   └── AGP sequence assembly with Gaps and provenance: sequence_from_agp
│   ├── 2) Operation within sequence
│   │   ├── Unified suffix-free entry: insert/delete/substitute/mask/trim/reverse_complement/rotate/canonical_origin
│   │   ├── sequence completion
│   │   │   ├── 1) reverse order
│   │   │   ├── 2) Complementary
│   │   │   ├── 3) Reverse complementation
│   │   │   └── Explicit duplex reverse complementation: Duplex.reverse_complement
│   │   ├── Transcription and Translation
│   │   │   ├── 1) Transcription
│   │   │   └── 2) Translation
│   │   ├── subsequence extraction
│   │   │   ├── Intervals and BED/GTF parsing: Region / read_regions
│   │   │   ├── FASTA indexing and regional access: index_reference / IndexedReference.fetch (pysam)
│   │   │   └── Chromosome/strand-aware interval sets: query_regions (intersection / union / difference)
│   │   ├── Sequence Editing
│   │   │   ├── Explicit duplexes, offsets, topology, sticky/blunt ends and digestion: Duplex.from_sequence / ends / cut
│   │   │   ├── Trimming: trim
│   │   │   └── Masking: mask
│   │   ├── sequence generation
│   │   │   ├── 1) Mutation Generation
│   │   │   ├── 2) Insert generation
│   │   │   ├── 3) Delete generation
│   │   │   ├── 4) Evolutionary algorithm-like operation
│   │   │   ├── 5) Fragment rearrangement
│   │   │   ├── 6) k-mer randomly scrambled
│   │   │   ├── 7) Multiple sequence recombination/crossover
│   │   │   └── Mutation library generation
│   │   ├── Sequence Concatenation and Assembly
│   │   │   ├── 1) Ordinary sequence splicing
│   │   │   ├── 2) Remove overlapping splicing
│   │   │   ├── Short-read/hybrid de novo assembly: assemble_reads (SPAdes; optional Nanopore)
│   │   │   ├── Executable cloning steps and intermediates: CloningPlan.add / apply / steps
│   │   │   ├── Gibson / homologous overlap assembly: assemble_duplexes (pydna)
│   │   │   ├── Golden Gate part design and reaction assembly: design_golden_gate / assemble_golden_gate
│   │   │   ├── Sequence-level assembly simulation: simulate_assembly
│   │   │   └── Duplex end compatibility and ligation: Duplex.ligate
│   │   ├── Circular Sequence Operation
│   │   └── Sequence segmentation
│   │       ├── 1) Fixed length, no overlap
│   │       ├── 2) Fixed length, overlapping sliding windows
│   │       ├── 3) Random interval
│   │       ├── 4) Multi-scale
│   │       └── 5) Short first and then long
│   ├── 3) Sequence deduplication
│   │   ├── 1. Standard deduplication
│   │   ├── 2. Reverse complementary deduplication
│   │   ├── 3. Ring equivalent deduplication
│   │   ├── 4. IUPAC-aware deduplication
│   │   └── 5. Approximate deduplication
│   ├── 4) Sequence clustering
│   │   ├── 1. Identity clustering
│   │   ├── 2. k-mer clustering
│   │   ├── 3. Fingerprint clustering
│   │   ├── 4. Hierarchical clustering
│   │   ├── 5. Representative sequence selection
│   │   ├── 6. Neural network clustering
│   │   ├── Label-stratified diversity sampling: diversity_sample
│   │   └── Feature clustering: cluster_features
│   │       ├── kmeans / minibatch_kmeans / gmm / spectral / affinity_propagation
│   │       └── mean_shift / dbscan / agglomerative / birch / mcl
│   ├── 5) Data set organization and division
│   │   ├── ENCODE tracking, collection summaries, citation export and provenance
│   │   ├── 1. Random and stable hash partitioning
│   │   ├── 2. Stratified random division
│   │   ├── 3. Similarity division
│   │   ├── 4. Cluster split
│   │   ├── 5. Species classification
│   │   ├── 6. Chromosome division
│   │   ├── 7. Individual Division
│   │   ├── 8. Divide by custom label
│   │   ├── 9. Leak detection
│   │   ├── Temporal and joint constraint splitting: temporal_split / joint_split
│   │   ├── Exclude records by species, chromosome, or custom metadata
│   │   └── Split quality evaluation: evaluate_split_quality
│   ├── 6) Sequencing Data Processing
│   │   ├── Pairing and damaged four-line FASTQ recovery: pair_reads / repair_fastq (recovery/discard statistics)
│   │   ├── Sorting, shuffling, and common records: arrange_reads / common_reads
│   │   ├── Regex renaming and order-sensitive, ID-independent sequence digests: rename_reads / sequence_digest
│   │   └── BAM sorting, mate repair, optional duplicate marking, indexing, and statistics: prepare_bam (samtools)
│   └── 7) Sequence Optimization
│       ├── Rules and custom evaluators: specification / Specification / SpecEvaluation
│       ├── GC, pattern avoidance/counts, sequence/translation/length, locked positions: gc / avoid_pattern / pattern_count / sequence / translation / length / avoid_changes
│       ├── Unique k-mers, stop/forbidden/rare codons: unique_kmers / avoid_stop / avoid_codons / avoid_rare_codons
│       ├── Exact-complement structure/match rules: hairpin / heterodimer / avoid_matches / regions_compatibility
│       ├── Tm and codon objectives: tm / cai / codon_usage / harmonize
│       ├── External homology constraint: blast_match_specification (BLAST)
│       ├── Hard constraints first, weighted objectives, synonymous edits, and circular sequences: optimize_sequence
│       ├── GenBank qualifier rule import: specifications_from_record
│       └── Edits, failure locations, search budgets, and JSON reports: OptimizationResult.write_report
├── 3. Data analysis
│   ├── 1) DNA descriptors + representations + fingerprints
│   │   ├── DNA descriptor
│   │   │   ├── 1) Fuzzy base statistics
│   │   │   ├── 2) length characteristics
│   │   │   ├── 3) base composition
│   │   │   ├── 4) GC/AT features
│   │   │   ├── 5) base skew
│   │   │   ├── 6) CpG Characteristics
│   │   │   ├── 7) k-mer statistics
│   │   │   ├── 8) Sequence Entropy
│   │   │   ├── 9) sequence complexity
│   │   │   ├── 10) Homopolymer
│   │   │   ├── 11) Repeat Ratio
│   │   │   ├── 12) window descriptor
│   │   │   ├── 13) Codon Statistics
│   │   │   ├── 14) All descriptors: all_descriptors (240 output fields in a fixed schema)
│   │   │   └── Descriptor schema discovery and provenance-aware dinucleotide property table loading
│   │   ├── Sequence Representations and Feature Processing
│   │   │   ├── 1) integer encoding
│   │   │   ├── 2) One-hot encoding
│   │   │   ├── 3) k-mer features
│   │   │   ├── 4) k-mer Sketch (MinHash/FracMinHash)
│   │   │   ├── 5) Neural network representation
│   │   │   │   ├── 1. AlphaGenome
│   │   │   │   ├── 2. Caduceus
│   │   │   │   ├── 3. DNABERT-2
│   │   │   │   ├── 4. Enformer
│   │   │   │   ├── 5. Evo 2
│   │   │   │   ├── 6. GENERator
│   │   │   │   ├── 7. GROVER
│   │   │   │   ├── 8. HyenaDNA
│   │   │   │   ├── 9. JanusDNA
│   │   │   │   ├── 10. LucaOne (default)
│   │   │   │   └── 11. Nucleotide Transformer v2
│   │   │   ├── 44 DNA encodings: encode_dna / DNA_ENCODINGS
│   │   │   │   ├── Word frequencies and approximate/noncontiguous words: Kmer / RCKmer / Mismatch / Subsequence
│   │   │   │   ├── Bases and positions: NAC / ANF / ENAC / binary / PS2 / PS3 / PS4
│   │   │   │   ├── Dinucleotides: CKSNAP / ASDC / DBE / LPDF
│   │   │   │   ├── Chemical/EIIP/mutual information: NCP / EIIP / PseEIIP / MMI
│   │   │   │   ├── Z-curve: Z_curve_9bit / Z_curve_12bit / Z_curve_36bit / Z_curve_48bit / Z_curve_144bit
│   │   │   │   ├── Positional physicochemical properties: DPCP / TPCP
│   │   │   │   ├── Auto/cross covariance: DAC / DCC / DACC / TAC / TCC / TACC
│   │   │   │   ├── Autocorrelation: NMBroto / Moran / Geary
│   │   │   │   ├── Pseudo nucleotide composition: PseDNC / PseKNC / PCPseDNC / PCPseTNC / SCPseDNC / SCPseTNC
│   │   │   │   └── Supervised encodings: PSTNPss / PSTNPds / KNN (explicit training_ids; physicochemical encodings require explicit property tables)
│   │   │   ├── Fitting and held-out reuse: analyze_features / FeatureAnalysis.transform
│   │   │   ├── Scaling: standard / minmax / l1 / l2 / none
│   │   │   ├── Reduction: pca / tsne / lda / latent_dirichlet (t-SNE has no held-out transform)
│   │   │   ├── Selection: chi2 / f_classif / f_regression / mutual_info_classif / mutual_info_regression
│   │   │   └── Selection: information_gain / pearson / variance; none for preprocessing only
│   │   └── DNA fingerprint
│   │       ├── 1) Hashed k-mer bit fingerprint
│   │       ├── 2) Panel existence fingerprint
│   │       ├── Motifs / restriction sites: motif_fingerprint / restriction_fingerprint
│   │       ├── Spatial GC / repeats / coding: gc_spatial_fingerprint / repeat_fingerprint / coding_fingerprint
│   │       ├── Thermodynamic / multiscale / hybrid: thermodynamic_fingerprint / multiscale_fingerprint / hybrid_fingerprint
│   │       └── Feature preprocessing: fit_preprocessor / transform (missing values, scaling, normalization, low-variance filtering)
│   ├── 2) Sequence search
│   │   ├── Universal search
│   │   │   ├── 1) Exact search
│   │   │   ├── 2) Subsequence search
│   │   │   ├── 3) Approximate matching
│   │   │   ├── 4) Reverse complementary search
│   │   │   ├── 5) Nearest neighbor search
│   │   │   ├── 6) Database Index
│   │   │   ├── Sketch index construction and persistence: build_sketch_index / save_sketch_index / load_sketch_index
│   │   │   ├── Local nucleotide BLAST: build_blast_database / local_blast (external BLAST+)
│   │   │   ├── BLAST XML/tab and protein HMMER2/3 parsing: parse_search_results (Biopython)
│   │   │   └── Profile construction/search and nucleotide table parsing: build_profile / search_profile / parse_nhmmer_table (hmmbuild / nhmmer / hmmsearch)
│   │   └── Sequence function search
│   │       ├── 1) motif search
│   │       ├── 2) start and stop codons
│   │       ├── 3) promoter motif
│   │       ├── 4) TF motif
│   │       ├── 5) restriction enzyme site
│   │       ├── 6) CRISPR PAM
│   │       ├── 7) palindrome sequence
│   │       ├── 8) Repeat upside down
│   │       ├── 9) Tandem Repeat
│   │       ├── 10) Microsatellite
│   │       ├── Open reading frames: scan_orfs
│   │       ├── CpG islands and low-complexity regions: find_cpg_islands / find_low_complexity_regions
│   │       ├── Native counts / PWM / PSSM, backgrounds, and degenerate consensus: MotifModel / build_motif
│   │       ├── Conservative FPR thresholds and both-strand/circular scans: fpr_threshold / scan_motif_model
│   │       ├── Primer properties, primer matching, and PCR simulation
│   │       ├── Primer design requests and the Primer3 CLI design adapter
│   │       └── CRISPR candidate and off-target scanning: scan_crispr_candidates
│   ├── 3) Physical and chemical properties
│   │   ├── 1. molecular weight
│   │   ├── 2. 260 nm extinction coefficient
│   │   ├── 3. melting temperature Tm
│   │   ├── 4. salt concentration correction
│   │   ├── 5. local melting characteristics
│   │   ├── 6. Synthesis Risk
│   │   └── Thermodynamic properties
│   │       ├── 1) thermodynamic parameters
│   │       ├── 2) Nearest-neighbor
│   │       ├── 3) Duplex stability
│   │       ├── 4) base stacking
│   │       ├── 5) Hairpin
│   │       ├── 6) Self-dimer
│   │       └── 7) Heterodimer
│   ├── 4) Double-chain thermodynamic expansion
│   │   ├── 1. ΔH, ΔS, ΔG and Tm of perfectly complementary duplexes
│   │   ├── 2. Unify double-strand stability results
│   │   ├── 3. Contribution of adjacent base pairs to steps
│   │   ├── 4. Conditions and Na⁺/K⁺ monovalent salt
│   │   ├── 5. Ka, Kd and duplex ratio
│   │   ├── 6. Theoretical melting curve
│   │   ├── 7. 5′/3′ end stability
│   │   ├── 8. DMSO and formamide experience correction
│   │   └── 9. Mg²⁺, dNTP, mismatch and dangling end for Primer3 CLI
│   ├── 5) Secondary structure properties
│   │   ├── 1. Dot-bracket structure analysis
│   │   ├── 2. Pairing probability and window accessibility
│   │   ├── 3. Collection defects of target structure
│   │   ├── 4. Thermodynamic probability of target structure
│   │   ├── 5. NUPACK Passive Availability Check
│   │   ├── 6. NUPACK single-complex ensemble analysis
│   │   └── 7. NUPACK tube multi-complex equilibrium
│   ├── 6) Three-dimensional structure and mechanical properties
│   │   ├── 1. Read a single PDB model
│   │   ├── 2. Read PDB multi-model collection
│   │   ├── 3. Explicit coordinate geometry analysis
│   │   ├── 4. NMR/Multi-model RMSF Flexibility
│   │   ├── 5. 3DNA bp_step.par standard parameter analysis
│   │   └── 6. DSSR JSON summary parsing
│   ├── 7) Conversion
│   │   ├── 1. Single chain molar extinction coefficient at 260 nm
│   │   ├── 2. Single-chain/double-chain theoretical optical properties
│   │   ├── 3. 1 nmol and mass corresponding to OD260
│   │   ├── 4. A260 to molar concentration, mass concentration and total amount
│   │   ├── 5. Concentration, quantity and mass interchange of substances
│   │   └── 6. Explicit correction of dyes and modifying groups
│   ├── 8) Property Prediction: Deep Learning and Machine Learning
│   │   ├── 1. Custom Model Training and Prediction
│   │   │   ├── User-defined classification and regression
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
│   │   ├── 2. Expression, Regulatory Tracks and Contact Maps
│   │   │   ├── AlphaGenome genomic tracks: RNA-seq, CAGE, PRO-cap, ATAC-seq, and DNase-seq
│   │   │   ├── Histone modifications, transcription factor binding, and 3D DNA contact maps
│   │   │   └── Enformer human/mouse regulatory tracks
│   │   ├── 3. Splicing and Genomic Regions
│   │   │   ├── Splice sites, splice usage, and splice junction counts
│   │   │   ├── Fourteen-class single-nucleotide segmentation
│   │   │   │   ├── 1. Protein-coding gene
│   │   │   │   ├── 2. Long noncoding RNA (lncRNA)
│   │   │   │   ├── 3. Exon
│   │   │   │   ├── 4. Intron
│   │   │   │   ├── 5. Splice donor
│   │   │   │   ├── 6. Splice acceptor
│   │   │   │   ├── 7. 5′ UTR
│   │   │   │   ├── 8. 3′ UTR
│   │   │   │   ├── 9. CTCF-bound region
│   │   │   │   ├── 10. PolyA signal
│   │   │   │   ├── 11. Tissue-specific enhancer
│   │   │   │   ├── 12. Tissue-invariant enhancer
│   │   │   │   ├── 13. Tissue-specific promoter
│   │   │   │   └── 14. Tissue-invariant promoter
│   │   │   ├── Long-context exon probability
│   │   │   ├── Fully fine-tuned NT Revised tasks: splice acceptors, splice-site types, and splice donors
│   │   │   └── Fully fine-tuned Genomic Benchmarks task: coding/intergenic regions
│   │   ├── 4. Epigenetic Marks and Open Chromatin
│   │   │   ├── Fully fine-tuned NT Revised tasks: H2A.Z, H3K27ac, H3K27me3, H3K36me3, H3K4me1/2/3, H3K9ac/me3, and H4K20me1
│   │   │   └── Fully fine-tuned Genomic Benchmarks task: human open-chromatin regions
│   │   ├── 5. Promoters and Enhancers
│   │   │   ├── Fully fine-tuned NT Revised tasks: enhancer detection/types; promoters, non-TATA promoters, and TATA promoters
│   │   │   ├── Fully fine-tuned Genomic Benchmarks tasks: Drosophila enhancers and mouse enhancers
│   │   │   └── Human Cohn/Ensembl enhancers, human regulatory-element types, and human non-TATA promoters
│   │   ├── 6. Variant Effects and Protein Stability
│   │   │   ├── SNV differences in genomic modalities and human/mouse regulatory tracks
│   │   │   ├── Evo 2 long-context zero-shot variant score
│   │   │   ├── Allele-conditional-probability variant score
│   │   │   ├── Protein stability and ESM protein variant effects
│   │   │   └── Splicing variant effects: SpliceAI / Pangolin / maximum delta score
│   │   ├── 7. Taxonomy and Molecular Function Classification
│   │   │   ├── Trained nucleic-acid/protein downstream tasks: SupKTax, GenusTax, and species classification
│   │   │   ├── Prokaryotic protein localization and ncRNA families
│   │   │   └── Fully fine-tuned Genomic Benchmarks task: human/worm
│   │   ├── 8. Molecular Relationships and Interactions
│   │   │   ├── Trained nucleic-acid/protein downstream task: nucleic-acid/protein central-dogma relation
│   │   │   └── Influenza A antigenic relationships, protein/protein and ncRNA/protein interactions
│   │   └── 9. Protein Features and Structure
│   │       ├── ESM protein features and sequence scores
│   │       └── NVIDIA NIM: AlphaFold2 / ESMFold protein structure prediction
│   ├── 9) Gene names and intervals
│   │   ├── Gene identifiers, transcripts, homology and database annotations
│   │   ├── Genomic DNA + species → annotate_genes / GeneAnnotationConfig
│   │   ├── Default RefSeq assembly resolution, explicit assembly, or local FASTA + GFF3
│   │   ├── minimap2 reference alignment and existing annotation mapping (not de novo prediction)
│   │   ├── Gene names, IDs, strands, overlapping genes, and gene/intergenic/unknown segments
│   │   ├── Python: 0-based half-open intervals; JSON/CLI: 1-based closed intervals
│   │   └── Pinned assembly/annotation cache, checksums, offline reuse, and provenance
│   ├── 10) Phylogenetics and Molecular Evolution
│   │   ├── Native trees and construction: TreeNode / build_tree (NJ / cluster-size-weighted UPGMA)
│   │   ├── Native tree comparison, consensus, and pruning: robinson_foulds / consensus_tree / prune_tree
│   │   ├── Native community phylogenetic distances: unifrac (weighted / unweighted; rooted trees)
│   │   └── Nei–Gojobori dN/dS: selection-pressure APIs in comparative and database namespaces
│   ├── 11) Statistics, Differential Analysis and Enrichment
│   │   ├── 1. Statistical and Differential Analysis
│   │   │   ├── Distance validation and permutation tests: distance_matrix / distance_test (PERMANOVA / ANOSIM / centroid PERMDISP)
│   │   │   ├── Matrix correlations and environmental subset search: mantel / bioenv
│   │   │   ├── Principal coordinates and constrained ordination: pcoa / ordination (CA / CCA / RDA)
│   │   │   ├── Compositional transforms and inverses: composition_transform (CLR / ILR / ALR)
│   │   │   ├── Native compositional differential analysis: ancom
│   │   │   ├── Alpha diversity: alpha_diversity (observed / Shannon / Simpson / inverse Simpson)
│   │   │   ├── Beta distances: beta_diversity (Bray-Curtis / Jaccard / Euclidean / cityblock / cosine)
│   │   │   ├── Original ANCOM-BC: ancom_bc (OLS / sandwich variance / EM sampling-bias correction)
│   │   │   ├── Differential expression: differential_expression (unshrunk NB2 GLM / size factors / explicit contrast / Wald)
│   │   │   └── Multiple testing: adjust_pvalues (Benjamini–Hochberg)
│   │   └── 2. Gene Set and Pathway Enrichment
│   │       └── Enrichr / PANTHER / Reactome / STRING enrichment analysis
│   ├── 12) Sequencing Data Analysis
│   │   ├── BAM summaries, coverage, and reference-coordinate consensus: alignment_summary / coverage_consensus (pysam; no indels in consensus)
│   │   ├── DNA-seq: dna_seq_workflow (QC → mapping → BAM → variants)
│   │   ├── Single-sample RNA-seq: rna_seq_workflow (HISAT2 / StringTie / featureCounts)
│   │   ├── ChIP-seq from prepared BAMs: chip_seq_workflow (MACS3 peaks + samtools coverage)
│   │   └── Both-strand local alignment and reference-coordinate consensus: sanger_consensus
│   ├── 13) Variant Calling, Annotation and Scoring
│   │   ├── 1. Variant Calling and Statistics
│   │   │   ├── Variant calling, reference normalization, and VCF indexing: call_variants (bcftools)
│   │   │   ├── Count VCF variants with filters: vcf_count_variants
│   │   │   ├── Normalize VCF records and multiallelic sites: vcf_normalize
│   │   │   ├── Summarize variant types and sample genotypes: vcf_summary_stats
│   │   │   └── Coding-variant fraction: coding_variant_fraction
│   │   ├── 2. Variant Annotation and Identifier Conversion
│   │   │   └── Ensembl VEP / MyVariant / FAVOR / VariantValidator annotations and identifier conversion
│   │   └── 3. Variant Effect Scores and Experimental Data
│   │       └── AlphaMissense / CADD / EVE / MaveDB scores and experimental effect data
│   ├── 14) Inheritance Simulation
│   │   └── Simulate allele inheritance at independent loci: simulate_inheritance
│   └── 15) Genomics and DNA Engineering
│       ├── Genomic relationships and breeding prediction: genomic_relationship_matrix / genomic_prediction
│       ├── Prokaryotic CDS prediction: predict_prokaryotic_genes (Pyrodigal)
│       ├── Methylation summaries, DNA storage and nanostructure design
│       ├── Genotype PCA: genotype_pca (NumPy; dosage-matrix input)
│       └── Genome ANI: average_nucleotide_identity (PyFastANI; no FastANI executable)
├── 4. Data evaluation
│   ├── 1) Commonly used evaluation indicators
│   │   ├── 1. Validity
│   │   │   ├── 1.1 Syntax validity: evaluate_validity
│   │   │   └── 1.2 Reference-genome consistency: evaluate_reference_consistency
│   │   ├── 2. Uniqueness
│   │   ├── 3. Diversity (mean raw pairwise Levenshtein distance)
│   │   ├── 4. Novelty (similarity / mean nearest-reference Levenshtein)
│   │   ├── 5. Fréchet DNA distance
│   │   ├── 6. Frag
│   │   ├── 7. SNN
│   │   ├── 8. Ambiguity
│   │   ├── 9. Redundancy
│   │   ├── Complexity and quality evaluation: evaluate_complexity / evaluate_quality
│   │   ├── Reference libraries, nearest-reference lookup, and memorization evaluation
│   │   ├── Composite scorecard: evaluate_scorecard
│   │   ├── Raw read quality control: quality_control (FastQC)
│   │   ├── Custom scalar property matrices: sequence_properties
│   │   ├── Gaussian MMD² / polynomial KID: kernel_distance
│   │   ├── Precision / Recall / Density / Coverage: manifold_metrics
│   │   ├── Calibrated clipped density / coverage: clipped_manifold_metrics
│   │   ├── Authenticity relative to generator training data: authenticity
│   │   ├── Pareto layers and exact hypervolume: pareto_ranks / hypervolume
│   │   ├── Property interval hit rates and histogram KL: property_hit_rate / property_kl
│   │   └── KDE distribution conformity and standard error: conformity_score
│   ├── 2) Similarity calculation
│   │   ├── 1. k-mer similarity
│   │   ├── 2. Fingerprint similarity
│   │   ├── 3. Sketch similarity
│   │   ├── 4. Dashing sketch similarity
│   │   ├── 5. Reference similarity
│   │   ├── Descriptor distribution similarity: evaluate_distribution_similarity
│   │   └── Vectors and matrices: kmer_vector_similarity / similarity_matrix / pairwise_matrix
│   └── 3) Sequence distance and alignment
│       ├── 1. Sequence distance
│       │   ├── 1) Hamming distance
│       │   └── 2) Edit distance
│       ├── 2. Pairwise sequence alignment
│       │   ├── 2) 1 Global comparison
│       │   ├── 2) 2 Local comparison
│       │   └── 2) 3 Semi-global comparison
│       ├── Native MSA, width/ID validation, and p-distance: MultipleAlignment.from_sequences / distances
│       ├── Multiple alignment: align_sequences (external MAFFT / MUSCLE 5 / Clustal Omega)
│       └── Single/paired-end mapping and private indexes: map_reads (BWA / Bowtie2)
├── 5. Visualization
│   ├── 1) sequence text diagram
│   ├── 2) position highlighting
│   ├── 3) style control
│   ├── 4) Gap display
│   ├── 5) linear sequence diagram
│   ├── 6) circular DNA diagram
│   ├── 7) Custom character representation
│   ├── 8) Alignment diagram
│   ├── Similarity matrix heatmap: plot_similarity_matrix
│   ├── Self-contained, filterable, expandable HTML reports: build_html_report / save_html_report
│   ├── Image export: save_svg / save_image (PNG, SVG, JPG; TIFF/PDF also supported)
│   ├── Tree plots: plot_tree (Matplotlib)
│   ├── Forward/reverse exact k-mer dotplots: dotplot (Matplotlib)
│   ├── Sanger chromatograms: plot_chromatogram
│   └── Qualitative gel simulation: simulate_gel (Matplotlib)
├── 6. Large Language Model Tools
│   ├── Python API: general-purpose LLM calls (dnakit.llm)
│   │   ├── Configuration and results: LLMConfig / LLMResult
│   │   ├── Sync/async: call_llm / acall_llm; custom prompts and current-step data
│   │   ├── Text/JSON input, text/JSON-object output, and token usage
│   │   ├── Local text files: file_path; local preview with read_llm_file
│   │   ├── Encoding and file/request size limits; explicit file input sends the complete text
│   │   ├── Provider presets: DeepSeek / Gemini / Qwen / Kimi / OpenRouter / OpenAI / custom
│   │   ├── Preset discovery: available_llm_providers / get_llm_provider
│   │   └── Custom compatible endpoint, model, credentials, and extra_body parameters
│   ├── MCP: Agent and MCP
│   │   ├── Local stdio server: dnakit-mcp; optional dependency dnakit[agent]
│   │   ├── MCP client integration (Codex, Claude, Cursor, VS Code, etc.)
│   │   ├── Catalog statistics: dnakit_catalog_stats
│   │   ├── Categories: list_dnakit_categories
│   │   ├── Paginated catalog: list_dnakit_tools
│   │   ├── Capability search: search_dnakit_tools
│   │   ├── Parameter schema and execution attributes: describe_dnakit_tool
│   │   ├── JSON-argument execution: call_dnakit_tool
│   │   ├── Python registry and server construction: default_tool_registry / create_server
│   │   └── Bounded JSON output, provenance, Python-only markers, and allow_side_effects authorization
│   ├── Skills: install the DNAKit skill and find analysis workflows
│   └── Agent: built-in natural-language task execution (run_agent)
│       ├── Tool discovery, parameter inspection, execution, and feedback
│       ├── Optional Python execution, failure recovery, and run logs
│       └── Keyword/semantic/hybrid search, permissions, and step limits
└── FAQ
    └── Engineering and Extensions
        ├── Bundled scripts, resource lookup, ENCODE credentials, notebooks and PDF reports
        ├── 1. unified backend interface
        ├── 2. Native and external tags
        ├── 3. Python API
        ├── 4. CLI
        │   ├── Environment and backends: info / backends
        │   ├── Data preparation: download-genome / convert / normalize / validate
        │   ├── Analysis and comparison: describe / fingerprint / search / orfs / compare
        │   ├── Datasets and reports: deduplicate / split / report / workflow
        │   ├── Custom properties: property-models / inspect-property-data / train-property / predict-property
        │   └── Gene annotation: annotate-genes
        ├── 5. Configuration workflows: load_workflow / run_workflow; artifact manifest I/O, verification, reuse
        ├── 6. Batch computation: iter_batch / run_batch; progress, cancellation, per-record errors
        ├── 7. Parallel computation: bounded threads and stable output order
        ├── 8. Chunking and Streaming
        ├── 9. cache
        ├── 10. Random Seed
        ├── 11. version tracking
        ├── 12. Errors and Warnings
        ├── 14. consistency verification
        ├── 16. Documentation and Tutorials
        ├── 17. Galaxy wrappers: convert / deduplicate / split / report (requires Galaxy)
        ├── Development validation: unit tests and consistency checks
        └── DAG definitions, dependency/cycle checks, and failure propagation: Step / run_workflow
```

See the [FAQ](../../faq/tutorial_overview.md) for catalog scope and counts.
