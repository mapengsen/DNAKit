# Genomics and DNA engineering: parameters and notes

See the [tutorial](../../api/features/38_genomics_extensions.md) for purpose, methods and examples. These APIs show progress by default; use `verbose=False` to disable it.

## Genomic relationships and breeding prediction {#breeding}

Install numeric dependencies:

```bash
python -m pip install -e ".[bio-analysis]"
```

Rows are individuals; columns must have identical marker order and allele orientation. Only diploid numeric dosages in [0, 2] are supported. `None` denotes missing dosage; `missing="error"` rejects missingness. Training means control imputation and centering, including prediction rows. All-missing and allele-monomorphic training markers are excluded and reported.

The relationship matrix is VanRaden method 1's additive relationship matrix; conventional kinship is G/2. Prediction is fixed-lambda ridge regression/RR-BLUP: `ridge_lambda` is required, and 1.0 is only an example. Select it using training-only validation and account for relatedness when splitting samples. This interface does not fit REML, polyploid or multi-trait models, environmental covariates, or mating allocation.

## Prokaryotic gene prediction {#genes}

Install the real Pyrodigal backend:

```bash
python -m pip install -e ".[genomics-genes]"
```

The tutorial's artificial example demonstrates the interface; it does not guarantee a gene call. Supply one bacterial/archaeal linear contig per call. `meta=True` uses upstream pretrained models. `meta=False` requires an explicit suitable `training_sequence`; failures do not silently switch modes. `translation_table` is restricted to training mode.

Outputs include CDS, proteins, model scores, translation tables and backend version. Coordinates are **0-based half-open** on the input; negative-strand CDSs are reverse-complemented. Ambiguous IUPAC symbols become N without shifting coordinates. The function does not identify eukaryotic intron/exon structures or assign biological gene functions. Circular input requires deliberate linearization.

## Methylation counts {#methylation}

The reader accepts modkit's 18-column extended bedMethyl, including gzip/BGZF. Column 11 is a percentage in [0, 100]. Summaries preserve chromosome, strand, modification and motif groups; the fraction is total modified count divided by total valid coverage. Zero coverage produces no fraction, and 5mC/5hmC remain separate. Defaults are `max_records=1_000_000`, `max_bytes=268_435_456` decompressed bytes and `max_line_bytes=65_536`. Exceeding a limit raises an error without returning partial data; partition larger inputs or select suitable limits.

This performs descriptive summaries, not differential methylation tests or disease prediction. The external `modkit-pileup` adapter produces bedMethyl from modified BAM files.

## DNA storage {#storage}

Install the codec:

```bash
python -m pip install -e ".[dna-storage]"
```

This adapts Chamaeleo's fixed-Huffman Goldman codec. Preserve the entire JSON packet: indices, lengths and SHA-256 are JSON metadata, not DNA-encoded. Intact reordered segments can be restored; missing/corrupted segments or digest mismatches fail. The digest detects corruption but provides neither correction nor authentication. This does not implement the original paper's complete overlapping archive protocol or any physical synthesis/sequencing experiment.

## DNA nanostructures and oxDNA export {#nanotechnology}

Install the geometry dependency:

```bash
python -m pip install -e ".[dna-nanotechnology]"
```

The design is scadnano JSON suitable for a `.sc` file. Export returns `configuration` and `topology` strings, without writing files or running simulation. Caller-supplied designs are also accepted, with complete canonical DNA assignments on every strand.

Only linear strands are supported: the verified scadnano 0.21.1 export loses circular closure bonds, so circular designs are rejected. Non-finite geometry or unsupported upstream offsets return explicit errors. Geometry construction/export does not establish folding, relaxation or molecular dynamics results.

## Genotype PCA {#pca}

`genotype_pca()` uses NumPy, installed through `.[genomics]` or `.[bio-analysis]`. Input rows are samples and columns are aligned diploid biallelic dosages in [0, 2], including fractional dosages. `None` and NaN denote missing values. Optional `sample_ids` must be unique, nonempty and match the row count.

Estimate each observed marker mean μ and frequency p=μ/2, remove all-missing or observation-invariant markers, impute other missing entries with μ, and calculate an SVD of `(dosage−μ)/sqrt(2p(1−p))`. Scores are `U*S`, retained-marker loadings are `V`, and covariance eigenvalues are `S²/(n_samples−1)`. Variance ratios use the full spectrum. `components` must be in 1..`min(n_samples−1, n_markers_used)`.

`marker_indices` maps loading rows to original columns. Means and frequencies retain the full input width, with null for all-missing markers. Imputation counts cover retained markers only. The largest absolute loading is positive; degenerate eigenspaces may still rotate between numerical libraries. Input arrays are not mutated and output is JSON-compatible.

This dense in-memory PCA does not read BED/VCF, prune LD or exclude related samples. Its eigenvalue scale differs from PLINK's relationship-matrix scale. Harmonize alleles and perform QC first; large cohorts should use specialized backends.

## Genome ANI {#ani}

`average_nucleotide_identity(query, reference)` uses PyFastANI from `.[genomics]`, without a separate FastANI executable. Supply DNA strings or local FASTA paths, including gzip. Strings containing only IUPAC letters are DNA; use `./filename` or a `Path` object for such filenames. All FASTA records describe contigs of ONE genome, with boundaries preserved, rather than separate samples.

The interface retains FastANI defaults, including 3,000 bp query fragments, for related prokaryotic genomes. Short contigs and incomplete terminal fragments cannot contribute complete query fragments. Inputs are loaded into memory; ambiguous bases become N and their count is recorded.

Only query→reference is computed; swap inputs for the reverse direction. `ani_percent` measures identity in homologous fragments. `query_aligned_fraction` is matching fragments divided by complete query fragments, not base-level coverage. A missing report returns `status="no_match"` with null identity and aligned fraction, never a zero-percent identity claim.

Results record versions, parameters, input digests, fragment counts and method citations. Progress shows loading, indexing and querying stages; `verbose=False` disables progress but does not suppress upstream warnings such as short-reference warnings. Precompiled wheel availability depends on the dependency and platform.

## Full parameters

::: dnakit.genomics
    options:
      members:
        - genotype_pca
        - average_nucleotide_identity
        - genomic_relationship_matrix
        - genomic_prediction
        - predict_prokaryotic_genes
        - read_bedmethyl
        - summarize_methylation
        - encode_dna_storage
        - decode_dna_storage
        - design_dna_duplex
        - export_oxdna_design
      show_source: false
