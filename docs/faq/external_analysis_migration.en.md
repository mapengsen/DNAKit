# Making external analyses available through pip

Review date: 2026-09-21. This review covers all **43 analysis IDs** in `src/dnakit/integrations/analyses.json`, checking their inputs, actual commands, output contracts, and related Python implementations. Recommendations in this page are migration proposals, not completed features or a claim that an extra has been published to PyPI.

The goal is to let users compute after installing DNAKit or a selected extra. Prefer existing Python computing libraries and reusable APIs; retain specialist implementations for complex algorithms. A program installed by pip can still run in a subprocess without a separate manual executable installation. This differs from an in-process Python binding.

## Implemented in this checkout

This checkout provides `dnakit.genomics.genotype_pca()` and `average_nucleotide_identity()`, with a `genomics` extra installing NumPy, Pyrodigal, pysam, and PyFastANI. The former defines a bounded genotype-matrix PCA method; the latter uses the PyFastANI computing library. The existing `plink2-pca` and `fastani` external command interfaces remain available. See the [genomics tutorial](../api/features/38_genomics_extensions.md#pip-genomics). These additions are not yet published to PyPI.

The new APIs passed 80 unit tests and synthetic calls with real NumPy, PyFastANI, Pyrodigal and pysam in an isolated Python 3.11 environment. Run `python others/validation/pip_genomics.py` from the checkout to verify PCA, ANI, CDS prediction and BAM statistics; the script blocks `subprocess.Popen`. This validates interfaces and synthetic data, not scientific accuracy on real samples or every supported platform.

Existing gene prediction, VanRaden relationships, and genomic prediction can reuse these dependencies. Users still supply genotype data, assembly FASTAs, and reference resources.

## Categories and priorities

Each analysis appears once: **A: 3, B: 8, C: 6, D: 26; total: 43**.

| Category | Meaning | Implementation rule |
| --- | --- | --- |
| A | Candidate for a pip computing library or program package | Verify wheels, Python versions, results, and installed entry-point discovery |
| B | A related DNAKit API already exists | Reuse components while documenting differences in models, inputs, and outputs |
| C | Suitable for a bounded independent feature | Use explicit method names and limit algorithm and data scope |
| D | Retain the specialist backend for now | Address packaging and resources before attempting algorithm replacement |

P0 means first priority; P1 means a subsequent bounded implementation or integration; P2 requires a dedicated algorithm, packaging, or resource project. Priorities are engineering recommendations, not release commitments.

## Assessment of all 43 analyses

| Analysis ID | Category / priority | Source-based assessment and proposed scope |
| --- | --- | --- |
| `plink2-pca` | B / P0 | The new independent `genotype_pca()` specifies frequency scaling, missingness and invariant-marker filtering. It accepts dosage matrices and uses sample-scaled covariance eigenvalues, leaving the external BED/eigenvec/eigenval contract intact. `analyze_features(method="pca")` continues to provide generic feature scaling. |
| `plink2-king` | C / P1 | The recipe runs `--make-king-table`. A bounded KING-robust implementation could support diploid biallelic autosomal hard calls, with pair-specific missingness, heterozygosity, IBS0, and zero-denominator handling. Compare against PLINK; VanRaden relationships are not KING estimates. |
| `plink2-ld-prune` | C / P1 | The recipe uses variant-count windows and steps for `--indep-pairwise`. An independent windowed correlation-pruning method is feasible; selection order, ties, and missingness can change retained markers. Do not label it PLINK-compatible without comparisons. |
| `plink2-gwas` | C / P1 | The current `--glm` recipe selects regression from phenotype coding. Start with quantitative-trait per-marker linear tests with covariates. Binary traits, separation/Firth correction, sample filtering, and multiple phenotypes require separate work. |
| `gcta-grm` | B / P1 | `genomic_relationship_matrix()` implements VanRaden method 1, dividing by the sum of marker `2p(1-p)` terms. GCTA standardization, pair-specific marker counts, and binary files require a separate implementation. The existing matrix is not a GCTA GRM replacement. |
| `gcta-reml` | D / P2 | The recipe runs GREML on a GCTA GRM, quantitative phenotypes, and covariates. Fixed-penalty `genomic_prediction()` does not estimate REML or heritability. Retain the backend pending convergence, boundary-solution, and standard-error validation. |
| `minigraph-pangenome` | D / P2 | The recipe builds a graph from a reference and one assembly using `-cxggs`. Retain graph construction; FASTA concatenation or ordinary alignment cannot replace GFA graph generation. |
| `hifiasm-assembly` | D / P2 | The recipe uses default single-sample HiFi assembly and returns GFA. Retain the assembler and evaluate platform wheels separately; simple sequence joining is not an assembly substitute. |
| `flye-assembly` | D / P2 | The recipe selects `--nano-raw`. A Python installation entry point exists upstream, but source installation involves compiled backends. Verify complete distributions before claiming installation without system build dependencies. [Installation](https://github.com/mikolmogorov/Flye/blob/flye/docs/INSTALL.md) |
| `sniffles-sv` | A / P1 | Upstream supports pip installation. An extra can install Sniffles and dependencies while retaining the same executable and single-sample long-read SV semantics. This automates installation rather than reimplementing SV calling. Verify dependency compatibility and virtual-environment entry points. [Upstream](https://github.com/fritzsedlazeck/Sniffles) |
| `modkit-pileup` | B / P2 | `read_bedmethyl()` and `summarize_methylation()` consume already generated statistics. The recipe generates pileups from MM/ML-tagged BAM, including probability thresholds and modification-specific coverage. Retain modkit; downstream summaries do not replace pileup generation. |
| `mosdepth-coverage` | B / P1 | Existing pysam-based `coverage_consensus()` counts filtered A/C/G/T bases in a bounded region. A dedicated interval-coverage API is feasible. Pair overlap, deletion/skipped CIGAR handling, genome-wide distributions, and mosdepth output files require separate comparison. |
| `fastani` | A / P0 | The new `average_nucleotide_identity()` uses FastANI internals through PyFastANI without a `fastANI` executable. It records query/reference direction, matched and total fragments, and null ANI when no hit is reported. [Library](https://pyfastani.readthedocs.io/en/stable/) / [Installation coverage](https://pyfastani.readthedocs.io/en/stable/guide/install.html) |
| `kraken2-classify` | D / P2 | The recipe classifies reads using an explicit local Kraken2 database. Retain classification and database-version semantics; generic k-mer similarity is not equivalent classification or abundance estimation. |
| `amrfinder-nucleotide` | D / P2 | The recipe runs AMRFinder on nucleotide assemblies with an explicit database. Retain curated data and interpretation rules; ordinary similarity searching does not reproduce the annotation pipeline. |
| `prokka-annotation` | B / P2 | Existing Pyrodigal CDS prediction serves a gene-prediction-only use case. Prokka returns a complete annotation set with additional steps; retain its full functional annotation backend. |
| `barrnap-rrna` | D / P2 | rRNA calling needs models and hit-filtering rules. Pin the upstream version and models first: current upstream includes an Infernal/Rfam route, so PyHMMER alone cannot reproduce all versions. [Upstream](https://github.com/tseemann/barrnap) |
| `fimo-motifs` | A / P1 | Evaluate the pip-installable LightMotif library. Existing `scan_motif_model()` provides basic PSSM scanning, while the FIMO contract filters site p-values. Compare backgrounds, pseudocounts, strands, p/q-values, and formats. PyMEMEsuite is marked experimental upstream and recommends LightMotif; its bundled MEME component has separate terms. [LightMotif](https://lightmotif.readthedocs.io/en/stable/) / [PyMEMEsuite](https://github.com/althonos/pymemesuite) |
| `cas-offinder` | C / P1 | A bounded CPU-only PAM/mismatch candidate enumerator without bulges is feasible. Define strands, coordinates, and input limits. Do not claim full Cas-OFFinder device or specification-file compatibility. |
| `bakta-annotation` | B / P2 | Pyrodigal only covers CDS prediction. The Bakta recipe uses a dedicated database and produces comprehensive annotations; retain functional transfer, RNA annotation, and database matching. |
| `busco-proteins` | D / P2 | The recipe performs offline assessment against an explicit lineage. PyHMMER can supply a search component, but BUSCO classification, duplication handling, and lineage thresholds need separate validation. A fraction of HMM hits is not BUSCO completeness. [PyHMMER](https://pyhmmer.readthedocs.io/en/stable/) |
| `busco-genome` | D / P2 | Genome mode also includes gene prediction. Retain BUSCO and lineage resources; prokaryotic CDS prediction does not cover eukaryotic genome mode. |
| `augustus-genes` | D / P2 | The recipe requires a species model and configuration directory. Retain eukaryotic gene-structure prediction; ORF scans or Pyrodigal do not replace intron-containing gene models. |
| `diamond-blastp` | D / P2 | The recipe builds a DIAMOND database from user proteins and returns 12-column hits. Retain the search engine. Small-scale dynamic programming could be a separate feature without equivalent performance or E-value claims. |
| `diamond-blastx` | D / P2 | The recipe searches translated nucleotides against proteins. Preserve reading-frame, translation, and significance semantics; ordinary nucleotide alignment is insufficient. |
| `eggnog-proteins` | D / P2 | Protein-mode emapper requires DIAMOND and version-matched eggNOG data. Retain the annotation pipeline and record compatible database/software versions together. |
| `miniprot-align` | D / P2 | Protein-to-genome alignment produces GFF with splicing and frameshift semantics. Retain the specialist aligner; ordinary protein alignment does not provide this coordinate model. |
| `liftoff-annotation` | D / P2 | The recipe transfers reference annotations to a target assembly. Installing Python code alone does not establish complete alignment dependencies. Retain transfer logic and unmapped/low-identity flags. |
| `checkm2-predict` | D / P2 | The recipe uses bins and a version-compatible DIAMOND database for model-based completeness/contamination. Manage search dependencies, models, and reference data; length or GC summaries are not replacements. |
| `quast-assembly` | D / P2 | The recipe performs reference-supported assembly assessment. Pure Python length, N50/L50, GC, and N-content summaries could be separate functions, but do not cover reference alignment, misassemblies, or the QUAST report. |
| `vsearch-dereplicate` | B / P0 | `deduplicate(equivalence="exact")` already groups identical sequences. Reuse it with FASTA/count handling while defining case, input abundance, minimum length, and ordering. Do not assume UC or VSEARCH label compatibility. |
| `vsearch-cluster` | B / P1 | Existing `cluster_sequences()` computes connected components of an exhaustive threshold graph. VSEARCH uses length-sorted greedy clustering and can produce different clusters. Preserve distinct method names and add a separate greedy mode if needed. |
| `mmseqs-protein-search` | D / P2 | The recipe runs `easy-search --search-type 1`. Retain protein search, prefiltering, and significance calculations; DNA similarity APIs do not cover it. |
| `mmseqs-protein-cluster` | D / P2 | The recipe sets identity, coverage, and `--cov-mode 1`. Retain protein coverage and clustering semantics; existing DNA connected-components clustering is not equivalent. |
| `abricate-screen` | D / P2 | The recipe searches an explicit local resistance/virulence database. Preserve database versions and identity/coverage semantics. Result parsing may be reusable, but parsing is not screening computation. |
| `mlst-typing` | D / P2 | The recipe matches installed schemes and a BLAST database to report ST. Exact known-allele matching could be a bounded future feature; full scheme selection, partial hits, and novel alleles retain the specialist backend. |
| `genomad-mobile-elements` | D / P2 | The recipe runs end-to-end with explicit models/database. Preserve classification, marker evidence, and thresholds; installing Python code does not provision those resources. |
| `repeatmasker-library` | D / P2 | The recipe uses a supplied repeat library and configured search engine. Retain repeat annotation; low-complexity masking alone can only be a separate subfeature. |
| `mageck-test` | D / P2 | The recipe performs unpaired two-condition RRA with median normalization. Retain the MAGeCK test; simple count differences or mean ratios are not equivalent gene-significance tests. |
| `mageck-count` | C / P1 | A single-sample, fixed-length guide exact-counter with explicit extraction position is feasible. Define duplicate guides, reverse complements, low-quality reads, and unmatched reads. Automatic trimming and multiple guide lengths need separate validation. |
| `crispresso2-amplicon` | D / P2 | The recipe quantifies one amplicon/guide. Preserve alignment, cut-window, and indel-assignment rules; mismatch counts cannot directly be called editing efficiency. |
| `drugz-unpaired` | C / P1 | Upstream has an inspectable Python implementation. Evaluate a pinned vendored adapter or an independently implemented unpaired method, preserving normalization, variance windows, guide aggregation, and FDR. Copied code must retain its source version and notices. [Upstream source](https://github.com/hart-lab/drugz/blob/master/drugz.py) |
| `gemma-lmm` | D / P2 | `-lmm 4` produces Wald, likelihood-ratio, and score tests. Fixed-lambda prediction or ordinary regression does not provide kinship-adjusted association; dedicated mixed-model validation is needed. |

## Existing components and subsequent delivery

These implementations existed when the review started:

| Existing Python API | Dependency extra | Verified implementation scope |
| --- | --- | --- |
| `predict_prokaryotic_genes()` | `genomics-genes` | In-process Pyrodigal CDS prediction, not complete functional annotation |
| `genomic_relationship_matrix()` / `genomic_prediction()` | `bio-analysis` | VanRaden method 1 / explicitly fixed-lambda ridge prediction |
| `analyze_features(method="pca")` | `bio-analysis` | Generic matrix PCA and feature scaling |
| `build_motif()` / `scan_motif_model()` | `bio-analysis` | PSSM construction and score-threshold scans; FPR is not genome-wide q-value |
| `read_bedmethyl()` / `summarize_methylation()` | Base installation | bedMethyl parsing, validation, and count summaries |
| `coverage_consensus()` | `bio-reads`, `bio-analysis` | Regional A/C/G/T coverage and majority-base consensus |
| `deduplicate()` | Base installation | Normalized DNA-record deduplication and audit, including exact/reverse-complement modes |
| `cluster_sequences()` | Depends on the similarity method | Bounded exhaustive threshold-graph connected components |

Source files: `src/dnakit/genomics/{genes,breeding,methylation}.py`, `src/dnakit/bioinformatics/{features,motifs,reads}.py`, and `src/dnakit/datasets/{deduplicate,clustering}.py`. See the [genomics tutorial](../api/features/38_genomics_extensions.md) and [bioinformatics reference](bioinformatics.md) for usage.

New implementations should use dedicated APIs or explicit backend selection and return method/version information, input-handling rules, and citations. Keep reproducible external-ID commands. Do not silently change statistical methods based on which executable is installed. Mark a migration complete only after installation, execution, and validation work is finished.

## Requirements for an install-and-run claim

1. Install and run in clean environments for declared Python/OS/architecture combinations. An sdist does not remove compiler requirements; PyFastANI explicitly documents source builds on some platforms.
2. Compare upstream results, coordinates, units, missing data, edge cases, no-hit behavior, and failures. Compare PCA subspaces or sign-invariant quantities rather than requiring identical eigenvector signs.
3. Report synthetic interface validation separately from scientific accuracy and real-data performance. Successful existing commands do not validate a new algorithm.
4. Record database/model/reference versions, checksums, sources, sizes, and cache locations separately. Downloads require an explicit operation; software installation does not establish resource readiness.
5. When including upstream code, retain pinned commits, modifications, licensing, and component provenance, and define an update process. Prefer established installable computing libraries over maintaining copied sources.

This implementation adds the two APIs and combined installation extra described above. Other rows retain their proposed scope and limitations; they are not all implemented in-process.
