# Update log

This project follows semantic versioning. It is still in the early development stage and the public interface may continue to be adjusted.

## [Unreleased]

## [0.1.2] - 2026-08-30

### Changed

- Added JPG to the unified visualization export interface. `image_type` now selects PNG, SVG, or JPG, while extensionless targets default to PNG. `SequencePlotConfig` adds `column_spacing`, `line_spacing`, and `symbol_map`; all plots now use square canvases without decorative in-figure titles.
- Changed the default GitHub repository and PyPI project description to English while retaining a Simplified Chinese switch.
- Added complete Chinese and English documentation, navigation, locale-preserving links, and language-isolated search results.
- Added Bioconda, GNU Guix, and Galaxy Tool Shed packaging files and removed obsolete planning, demo, and delivery-report pages.

## [0.1.1] - 2026-08-29

### Changed

- Release of the first DNAKit version without the development version suffix.
- Merged third-party declarations into `DISCLAIMER.md` and moved project reference information to `README.md`.
- Streamline root directory files and simultaneously update installation, packaging and release configurations.

## [0.1.0.dev0] - 2026-08-28

### Added

- Unified `dnakit.DNA(...)` facade for ordinary users: the same entrance receives single/multiple sequences, record mapping, ID, topology, metadata, feature and old core objects; subscripts and slices still return `DNA`.
- Newly added `dnakit.representations`: 11 basic DNA models can be selected to extract fixed-length reps, and checkpoints are downloaded to the current directory `ckpt/` by default and reused; newly added `neural_cluster_sequences()` to complete L2, optional PCA, seed-fixed k-means, central nearest representative sequence and complete result audit.
- New `evaluate_frechet_distance()`: By default, the LucaOne representation and L2 normalization of DATA-027 are multiplexed to calculate the Fréchet representation distribution distance of two DNA collections with the sample space equivalent algorithm.
- Added `evaluate_fragment_similarity()` and `evaluate_snn()`: implement DNA adaptation of MOSES Frag/SNN using exact k-mer occurrence count cosine similarity and hashed k-mer fingerprint nearest neighbor Tanimoto mean respectively.
- `dnakit.datasets.exclude_species()` and `exclude_chromosomes()`: Exactly exclude one or more species/chromosomes based on explicit metadata; structured error reporting when target fields are missing.
- The MkDocs document site enables the search box, search suggestions and result keyword highlighting in the Chinese site.
- `ValidationConfig(sequence_length=...)`, `dnakit validate --sequence-length` and workflow `validate` steps now determine whether a sequence is legal with standardized accuracy `symbol_length`.
- Add a new MIT license file to the root directory, mark the copyright owner as Pengsen Ma, and include the license in the release archive.
- Added `DISCLAIMER.md` (including third-party statements) and the "Acknowledgments and Citations" merge page of the documentation site to clarify the boundaries of external tools, user data, costs, and non-clinical use.
- Bounded GenBank subsets, GFF3, BED3–6, AGP 2.1, plain FASTA/strict four-row FASTQ persistent indexing, coordinate extraction, chunked iteration and metadata management.
- `read_table()`, `export_table()` and immutable audit results for explicit `TableSchema` bounded reading and writing of CSV/TSV/JSON/Parquet.
- Complexity/repeat descriptors, motif/PWM, ORF, restriction enzyme, PAM, CpG island, palindrome, inverted repeat, tandem repeat, STR and low complexity scans.
- Molecular weight, theoretical 260 nm extinction coefficient of unmodified ssDNA based on published nearest-neighbor absorbance parameters, and internally reproduced SantaLucia 1998 Tm, salt correction, nearest-neighbor, ΔG/ΔH/ΔS, stacking, full complementary duplex and local Tm; `duplex_stability()` Canonical mismatch/dangling heterodimer can be evaluated explicitly using the user-installed Primer3 CLI adapter.
- Single/double strand optical properties, OD260/A260/molarity/mass concentration/amount/mass conversion of species, explicit dye/modification corrections, and empirical corrections for Ka/Kd, double strand ratio, theoretical melting curves, end stability, Na⁺+K⁺ and DMSO/formamide.
- dot-bracket secondary structure and pairing probability derived metrics, explicit NUPACK 4 adapter performed after separate license installation, and PDB coordinate geometry, SASA, volume, shape, backbone dihedral angles, approximate helix parameters, NMR ensemble RMSF and 3DNA/DSSR analysis.
- RCSB 1BNA, 1AC7, 139D local check samples, SHA-256 manifest, structural analysis script with progress bar and machine-readable results.
- `Primer3CLIAdapter` for pure external CLI, requires explicit `oligotm`/`ntthal` paths, supports Tm, hairpin, self-dimer and heterodimer, and records permissions, paths and provenance; downstream unified results verify sequence, conditions and structure option bindings.
- MinHash/FracMinHash, motif/restriction enzyme/GC/repeat/coding/fixed 16-dimensional thermodynamics/hybrid/multi-scale fingerprint and preprocessor for training set fitting; structural features use explicit adapter or fixed missing strategy.
- Bounded approximate matching, global/local/semi-global linear/affine gap alignment, sketch similarity, persistent indexing and Top-k nearest neighbors.
- Ring/IUPAC/approximate deduplication, identity/edit/k-mer/fingerprint clustering, hierarchical clustering, representative sequences, temporal partitioning, multi-constraint heuristic partitioning, leakage and partitioning quality evaluation.
- `SplitConfig(method="hash")` Order-independent stable hash partitioning: using SHA-256, seed, and unique `record.id`, and maintaining exact ratio quotas.
- validity, ambiguity, quality, complexity, uniqueness, diversity, redundancy, reference-scoped novelty/memorization, distribution similarity, synthesis-risk and regular scorecard.
- Digestion, end/ligation, PCR, primer matching, optional structural primer attributes, stringent design requests and `Primer3CLIDesignAdapter` for explicit paths, Gibson/LCR/Golden Gate/BioBrick sequence-level assembly, CRISPR candidates, rule optimization, codon optimization and mutation libraries; Primer3 design candidates back-check template coordinates, left and right primer sequences, and product lengths.
- Linear/circular feature plots, alignment plots, similarity matrices, self-contained HTML reports, and optional PNG/TIFF/PDF 600 dpi export.
- Backend registry, BLAST/MMseqs2/sourmash passive metadata/explicit version handles, strict Dashing Jaccard/Top-k adapter that requires the caller to provide an executable, content-addressed JSON caching, threaded batching, stable resume, JSON/YAML multi-step workflows, and auditable run manifests.
- Unified CLI capabilities such as `describe`, `fingerprint`, `search`, `orfs`, `compare`, `report` and `workflow`, and compatible development portal `python -m dnakit.cli.workflow run`.
- Local correctness validator, microbenchmark with seed/environment/parameters/sample-by-sample logging, DNAKit/Biopython peer-to-peer task comparison, selected public callable source code line calibers, advanced notebook, full workflow example and stage 4/5 delivery report.
- Manual workflow configuration for TestPyPI and official PyPI; the build review job has no release credentials, only the upload job has `id-token: write`, both with license, version and explicit confirmation gates, this stage is not triggered.

### Changed

- `evolution_generate()` retains the `mut_frac` parameter name, and adds base-by-base `insert_frac`, `delete_frac` probabilities; insertion is available `insert_min/insert_max` to select single base or random length fragments, and the algorithm audit version is upgraded to `dnakit-evoaug-v3`. Single-shot contiguous segment insertion/deletion continues to be provided by `indel_generate()`.
- The document merges FP-004 Canonical k-mer into the Canonical mode of FP-003 k-mer features, and renames FP-005 k-mer Sketch; the public Python API and stable tracking number remain compatible.
- Ordinary file reading is unified as `read(..., mode="dna"|"stream")`, verification is unified as `validate()`; suffix-free editing, reverse complementation and ring operations input synchronization comments on `DNA` and return `DNA`, and the old `read_one/read_set/validate_set/*_record` entry continues to be compatible.
- The default checkpoint for DNA basic model characterization and neural network clustering has been changed from GROVER to `LucaGroup/LucaOne-gene-step36.8M`; the default cache directory is `ckpt/lucaone-gene-step36-8m/`, and custom checkpoint code still requires explicit authorization.
- Removed the built-in DiProDB 240 values; the fixed 240 field schema remains unchanged, the default last 60 items are `None`, and users can load strict 15×16 JSON tables that they have permission to use.
- Removed `primer3` extra, `primer3-py` dependencies and Python binding adapter; Primer3 is unified and changed to be installed separately by users and provides an explicit CLI path.
- `normalize()` Added `keep_ambiguous`, `keep_u` and `keep_other`; retain IUPAC ambiguous bases by default and audit to delete `U` and other non-DNA characters.
- The requirements tracking matrix is ​​changed to check the current source code, testing and local backend status item by item, and expired `planned` statements are no longer retained.
- Merge independent identity, query coverage and target coverage functional items into `SIM-008` sequence alignment results and remove duplicate numbers.
- Batch processing supports bounded thread execution, stable raw index/derivative seed and completed ID skipping.
- The visualization scope is reduced to sequence, feature, alignment and similarity matrices, while general image export and HTML reporting remain unchanged.
- Visual export expanded to atomic multi-format output with hard caps on number of pixels, records, and bytes of embedded results.
- README, Installation, Quick Start, API, FAQ, References/License and Documentation site navigation synced to Phase 4/5 real state.

### Validation

- 11 model registry, checkpoint reuse/download contract, IUPAC/Gap/matrix boundary and rep → k-means/PCA are covered by unit tests; GROVER uses real local checkpoint to complete rep + k-means smoke of 4 sequences, and the remaining models remain in conditional state.
- Boundary validation for artificial minisequences, IUPAC, empty sequences, 200,000 nt input, and cross-origin circular restriction enzyme sites.
- Restriction enzymes, molecular weight terminal conventions, global alignments, literal searches, and single/complete/average linkage passed with Biopython.
- Primer3 CLI adapter validates command whitelisting, Boulder-IO, field/unit parsing, failures, timeouts, and output caps with temporary controlled aliases; no longer writes old Python bindings against current scientific validation.
- Dashing adapter's fixed commands, matrix parsing, failure/timeout/output/input mutation boundaries are covered by controlled surrogates, local `v1.0.2-4-g0635` two-sequence exact document example smoke passes; scientific difference has not been done yet.
- Primer3 is not a Python dependency and is not packaged with DNAKit; the CLI adapter records a `GPL-2.0-or-later` prompt, and users must still check the actual installed version.
- The NUPACK adapter's field mapping, bounds, and error handling have been tested with controlled doubles; the current environment does not have NUPACK, so true NUPACK numerical differential conclusions still do not exist.

### Remaining boundaries

- DNAKit's own license has been determined to be MIT; this version is released as a development preview, and dependency compatibility still needs to be continuously reviewed.
- GitHub Pages has not yet been deployed, nor has the paper replication experiment been performed.
- NUPACK, Primer3, DSSR/3DNA and Dashing are not installed/downloaded with the project; DiProDB values ​​are no longer built-in; the license, possible costs and real scientific differences of each external tool/data are still reviewed by the user based on the actual version.
- DNAKit Parquet write out /`read_table()` round trip passed under PyArrow 25.0.1; other engine/version cross matrices have not yet been executed.
