# DNAKit function module index

## 1) Data preparation

- [Core Data Object](01_core_objects.md) (8 items)
- [File reading and writing](02_io_data.md) (Category 4; large file capability as optional mode for IO-001)
- [Download](15_download.md) (common species names, reference genome FASTA and public database data)
- [Data query](22_database_query.md) (public database query capability)
- [Legality check](03_validation.md) (unified check entrance)

## 2) Data processing

- [DNA sequence standardization](03_standardization.md) (STD-001 character standardization)
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
- [Sequence Deduplication](10_deduplication.md) (DATA-001–006)
- [Sequence Clustering](10_clustering.md) (DATA-007–011, DATA-027 Neural Network Clustering)
- [Data set organization and division](10_datasets.md) (DATA-012–018, DATA-023, including custom label division)

## 3) Data analysis

- **DNA Descriptors, Characterizations and Fingerprints**
  - [DNA descriptor](05_all_descriptors.md)
  - [Sequence characterization](08_fingerprints.md)
    - [Neural Network Representation](08_fingerprints.md#neural-representations) (DATA-027: 11 DNA basic models rep)
  - [DNA fingerprint](08_feature_engineering.md)
- **Sequence Search**
  - [Universal Search](09_search.md) (SIM-001, SIM-002, SIM-004, SIM-005, SIM-014, SIM-015)
  - [Sequence Function Search](06_patterns.md) (PAT-001, PAT-003–007, PAT-009–012)
- [Physical and chemical properties](07_physicochemical.md)
- [Double-chain thermodynamics extension](19_duplex_thermodynamics.md)
- [Secondary structure properties](20_secondary_structure.md)
- [Three-dimensional structure and mechanical properties](21_structure3d.md)
- [Conversion](18_optics_concentration.md)

## 4) Data evaluation

- [Commonly used evaluation indicators](12_evaluation.md) (EVAL-001, EVAL-002, EVAL-005–008, EVAL-016–018)
- [Similarity calculation](09_similarity_alignment.md) (SIM-010–013, SIM-016 and reference/distribution similarity)
- [Sequence Distance and Alignment](09_alignment.md) (SIM-006–008)

## 5) Visualization

- [Visualization](13_visualization.md) (VIZ-001–009)

## 6) Engineering and Expansion

- [Backend, Performance and Reproducibility](14_engineering.md) (ENG-001–017)
