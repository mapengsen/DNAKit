# DNAKit function tree

> Update date: 2026-08-31. This function tree is synchronized with the current web navigation and function pages; the numbering follows the corresponding function page or demand tracking number, and web extension functions that are not assigned numbers retain their page names.

```
Quick Start Tutorial
├── Function Tree Overview
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
│   │   │   └── Advanced chunking/indexing: IO-005 is reserved for requirement tracking numbers only
│   │   ├── 2. IO-002 comment format
│   │   ├── 3. IO-003 table format
│   │   └── 4. IO-004 compressed file
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
│   │   └── 21. DBD-021 coordinate conversion file
│   ├── (4) Data query
│   │   ├── 1. Basic query
│   │   ├── 2. Coordinates, transcripts and regional annotations
│   │   └── 3. Sequencing, expression, regulation and comparison of genomes
│   └── (5) Legality check
│       └── 1. Unified sequence, record and data set legality check
├── 2. Data processing
│   ├── (1)DNA sequence standardization
│   │   └── 1. STD-001 character standardization
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
│   └── (5) Data set organization and division
│       ├── 1. DATA-012 · Random and stable hash partitioning
│       ├── 2. DATA-013 · Stratified random division
│       ├── 3. DATA-014 · Similarity division
│       ├── 4. DATA-015 · Cluster split
│       ├── 5. DATA-016 · Species classification
│       ├── 6. DATA-017 · Chromosome division
│       ├── 7. DATA-018 · Individual Division
│       ├── 8. Divide by custom label
│       └── 9. DATA-023 · Leak detection
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
│   │   │   └── 14. All descriptor calculations (240 items)
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
│   │       └── 2. Panel existence fingerprint
│   ├── (2) Sequence search
│   │   ├── Universal search
│   │   │   ├── 1. SIM-001 · Exact search
│   │   │   ├── 2. SIM-002 · Subsequence search
│   │   │   ├── 3. SIM-004 · Approximate matching
│   │   │   ├── 4. SIM-005 · Reverse complementary search
│   │   │   ├── 5. SIM-014 · Nearest neighbor search
│   │   │   └── 6. SIM-015 · Database Index
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
│   │       └── 10. PAT-012 Microsatellite
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
│   ├── (8) Deep-learning property prediction
│   │   ├── 1. AlphaGenome 11 output modalities and SNV track differences
│   │   ├── 2. Enformer human/mouse regulatory tracks and SNV differences
│   │   ├── 3. SegmentNT 14-class single-nucleotide segmentation
│   │   ├── 4. Evo 2 zero-shot variant score and exon probability
│   │   ├── 5. GENERator zero-shot variant score
│   │   └── 6. Ten pretrained LucaOneTasks heads
│   └── (9) Optional bioinformatics functions
│       ├── 1. Ensembl VEP and variant identifier recoding
│       ├── 2. ClinVar, dbSNP, and gnomAD annotation
│       ├── 3. Nei–Gojobori dN/dS
│       └── 4. Golden Gate design and reaction assembly
├── 4. Data evaluation
│   ├── (1) Commonly used evaluation indicators
│   │   ├── 1. EVAL-001 Validity
│   │   ├── 2. EVAL-005 Uniqueness
│   │   ├── 3. EVAL-006 Diversity
│   │   ├── 4. EVAL-008 Novelty
│   │   ├── 5. EVAL-016 Fréchet DNA distance
│   │   ├── 6. EVAL-017 Frag
│   │   ├── 7. EVAL-018 SNN
│   │   ├── 8. EVAL-002 Ambiguity
│   │   └── 9. EVAL-007 Redundancy
│   ├── (2) Similarity calculation
│   │   ├── 1. SIM-010 · k-mer similarity
│   │   ├── 2. SIM-011 · Fingerprint similarity
│   │   ├── 3. SIM-012 · Sketch similarity
│   │   ├── 4. SIM-013 · Dashing sketch similarity
│   │   └── 5. EVAL-011 Reference similarity
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
│   └── 8. VIZ-008 Alignment diagram
└── 6. Engineering and Expansion
    ├── 1. ENG-001 unified backend interface
    ├── 2. ENG-002 Native and external tags
    ├── 3. ENG-003 Python API
    ├── 4. ENG-004 CLI
    ├── 5. ENG-005 Configuration Workflow
    ├── 6. ENG-006 Batch calculation
    ├── 7. ENG-007 Parallel Computing
    ├── 8. ENG-008 Chunking and Streaming
    ├── 9. ENG-009 cache
    ├── 10. ENG-010 Random Seed
    ├── 11. ENG-011 version tracking
    ├── 12. ENG-012 Errors and Warnings
    ├── 13. ENG-013 unit testing
    ├── 14. ENG-014 consistency verification
    ├── 15. ENG-015 performance benchmark
    ├── 16. ENG-016 Documentation and Tutorials
    └── 17. ENG-017 optional graphics entry
```
