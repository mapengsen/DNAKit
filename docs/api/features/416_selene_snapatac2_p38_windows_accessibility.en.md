# Selene genome windows and SnapATAC2 accessibility

These DNAKit interfaces call their documented upstream implementations. Optional dependencies are imported only when needed. Coordinate conventions, defaults and parameter choices are specified per function.

Validation scope: Reference FASTA windows and fragment import are tested offline. Training and network/regression inference with user models are not validated by these small examples.

## `dnakit.predictions.selene_sdk_p38_functions.sample_selene_dna_windows`

**Purpose:** Draw labelled DNA windows with a selene sampler and optionally export the sampled positions.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def sample_selene_dna_windows(sampler: str, *, reference_fasta: str | None=None, features_bed: str | None=None, feature_names: Sequence[str] | None=None, intervals_path: str | None=None, bed_path: str | None=None, bed_has_targets: bool=False, mat_path: str | None=None, sequence_key: str | None=None, targets_key: str | None=None, sequence_length: int=1000, center_bin_to_predict: int=200, feature_threshold: float=0.5, batch_size: int=8, mode: str='train', validation_holdout: Sequence[str]=('chr6', 'chr7'), test_holdout: Sequence[str]=('chr8', 'chr9'), seed: int=436, export_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Draw labelled DNA windows with a selene sampler and optionally export the sampled positions.

``sampler`` selects the selene class: ``random_positions`` (``RandomPositionsSampler``: random positions
of the genome), ``intervals`` (``IntervalsSampler``: positions inside the BED ``intervals_path``),
``bed_file`` (``BedFileSampler``: windows listed in ``bed_path``; with ``bed_has_targets`` the last column
holds ``;``-separated feature indices) or ``mat_file`` (``MatFileSampler``: sequences and targets from a
MAT or HDF5 file, keys ``sequence_key`` and ``targets_key``).

The first two need ``reference_fasta``, a tabix-indexed ``features_bed`` (bgzip) and ``feature_names``;
chromosomes in ``validation_holdout`` and ``test_holdout`` are kept out of ``train`` mode (``mode`` is
``train``, ``validate`` or ``test``). ``feature_threshold`` is the overlap fraction that makes a label
positive. With ``export_dir`` the sampled positions are written to ``<export_dir>/<mode>_data.bed``
(``OnlineSampler.save_dataset_to_file``; ``random_positions`` and ``intervals`` only). Returns the shapes,
the decoded windows (at most 50) and their labels.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.selene_sdk_p38_functions.lookup_selene_genomic_features`

**Purpose:** Look up which genomic features overlap an interval, and whether any feature is positive there.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def lookup_selene_genomic_features(features_bed: str, feature_names: Sequence[str], chrom: str, start: int, end: int, *, strand: str | None=None, feature_threshold: float | None=0.5, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Look up which genomic features overlap an interval, and whether any feature is positive there.

Uses ``selene_sdk.targets.GenomicFeatures.get_feature_data`` (a 0/1 label per feature name for the interval
``chrom:start-end``, 0-based half-open) and ``is_positive``. ``features_bed`` is a bgzip-compressed,
tabix-indexed BED with the feature name in column 4. A feature is labelled 1 when it overlaps at least
``feature_threshold`` of the interval (a fraction in [0, 1]; null requires any overlap), and the interval is
positive when at least one feature is. ``strand`` is ``+``, ``-`` or null.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.selene_sdk_p38_functions.get_selene_genome_window`

**Purpose:** Get the DNA sequence or the one-hot encoding of a genomic interval from a FASTA, using

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def get_selene_genome_window(reference_fasta: str, chrom: str, start: int, end: int, *, output: str='sequence', strand: str='+', pad: bool=False, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Get the DNA sequence or the one-hot encoding of a genomic interval from a FASTA, using
``selene_sdk.sequences.Genome``.

``output`` is ``sequence`` (``get_sequence_from_coords``: the upper-case sequence; the reverse complement when
``strand`` is ``-``) or ``encoding`` (``get_encoding_from_coords_check_unk``: a length-by-4 one-hot matrix in
the order A, C, G, T, with an all-zero row per N, plus a flag telling whether unknown bases occurred).
Coordinates are 0-based, half-open. With ``pad`` out-of-range bases become N; otherwise an out-of-range
interval returns an empty result.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.selene_sdk_p38_functions.read_selene_ism_scores`

**Purpose:** Load an in silico mutagenesis result table and extract a per-base score matrix for one feature.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def read_selene_ism_scores(ism_path: str, *, feature: str | None=None, reference_mask: float | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Load an in silico mutagenesis result table and extract a per-base score matrix for one feature.

Uses ``selene_sdk.interpret.ISMResult.from_file`` (a tab-separated table with columns ``pos``, ``ref``,
``alt`` and one column per feature; an optional first row with ``ref`` and ``alt`` equal to ``NA`` holds the
reference prediction) and ``ISMResult.get_score_matrix_for``. Returns the reference sequence and the feature
names; with ``feature`` also the length-by-4 matrix (A, C, G, T) where alternative-base cells hold the
mutant score and reference-base cells hold the reference prediction (or ``reference_mask``).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.selene_sdk_p38_functions.predict_selene_sequences_strand_symmetric`

**Purpose:** Strand-symmetric prediction: combine a model's forward and reverse-complement outputs, using

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def predict_selene_sequences_strand_symmetric(model_path: str, sequences: Sequence[str], *, mode: str='mean', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Strand-symmetric prediction: combine a model's forward and reverse-complement outputs, using
``selene_sdk.utils.NonStrandSpecific``.

``model_path`` holds a whole ``torch.nn.Module`` saved with ``torch.save(model)`` (it is unpickled, so only
use trusted files). Each sequence is one-hot encoded, the model is run on it and on its reverse complement,
and the two outputs are combined with ``mode`` ``mean`` or ``max``. Sequences share one length. Returns the
combined predictions and the plain forward predictions for comparison.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.regulation.snapatac2_p38_functions.import_atac_fragments`

**Purpose:** Import an ATAC fragment file into a snapATAC2 AnnData file with per-barcode QC, using

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def import_atac_fragments(fragment_file: str, chrom_sizes: Mapping[str, int], output_h5ad: str, *, min_num_fragments: int=200, is_paired: bool=True, sorted_by_barcode: bool=False, whitelist: Sequence[str] | None=None, mitochondrial_chromosomes: Sequence[str]=('chrM', 'M'), tile_bin_size: int | None=None, limit: int=100, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Import an ATAC fragment file into a snapATAC2 AnnData file with per-barcode QC, using
``snapatac2.pp.import_fragments`` (and ``add_tile_matrix`` when ``tile_bin_size`` is given).

``fragment_file`` is a (optionally gzip-compressed) tab-separated file with chromosome, start, end,
barcode and count columns, 0-based start. ``chrom_sizes`` maps chromosome to length. Barcodes with fewer than
``min_num_fragments`` unique fragments are dropped; ``whitelist`` restricts barcodes. ``sorted_by_barcode``
tells snapATAC2 the file is already barcode-sorted (the default sorts internally). The result file has the
QC columns ``n_fragment``, ``frac_dup`` and ``frac_mito``. With ``tile_bin_size`` (bp) a cell-by-genomic-bin
count matrix is added (insertion counting, as snapATAC2's default), which the later steps need.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.regulation.snapatac2_p38_functions.filter_atac_cells`

**Purpose:** Filter ATAC cells by fragment count and TSS enrichment, using ``snapatac2.pp.filter_cells``.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def filter_atac_cells(input_h5ad: str, output_h5ad: str, *, min_counts: int | None=1000, min_tsse: float | None=5.0, max_counts: int | None=None, max_tsse: float | None=None, gene_gtf: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Filter ATAC cells by fragment count and TSS enrichment, using ``snapatac2.pp.filter_cells``.

The input file is copied to ``output_h5ad`` and the copy is filtered. The defaults are snapATAC2's (at least
1000 fragments and TSS enrichment 5). TSS enrichment is taken from ``obs["tsse"]``; when absent and a
thresholds on it is set, give ``gene_gtf`` and it is computed with ``snapatac2.metrics.tsse`` first, or set
``min_tsse`` and ``max_tsse`` to null to filter by counts only. Returns the cell counts before and after.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.regulation.snapatac2_p38_functions.select_atac_features`

**Purpose:** Select informative genomic bins or peaks, with blacklist and whitelist regions, using

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def select_atac_features(input_h5ad: str, output_h5ad: str, *, n_features: int=500000, filter_lower_quantile: float=0.005, filter_upper_quantile: float=0.005, whitelist: str | None=None, blacklist: str | None=None, max_iter: int=1, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Select informative genomic bins or peaks, with blacklist and whitelist regions, using
``snapatac2.pp.select_features``.

The input needs a count matrix (for example from ``import_atac_fragments`` with ``tile_bin_size``). Features
are ranked by accessibility after dropping the lower and upper quantiles; ``blacklist`` and ``whitelist`` are
BED files of regions to exclude or to keep. The mask is stored in ``var["selected"]`` of the copy written
to ``output_h5ad``. Returns the number and names (first 50) of the selected features.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.regulation.snapatac2_p38_functions.score_atac_doublets`

**Purpose:** Doublet scoring of scATAC cells (Scrublet on accessibility features), using ``snapatac2.pp.scrublet``.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def score_atac_doublets(input_h5ad: str, output_h5ad: str, *, use_selected_features: bool=True, n_comps: int=15, sim_doublet_ratio: float=2.0, expected_doublet_rate: float=0.1, random_state: int=0, limit: int=100, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Doublet scoring of scATAC cells (Scrublet on accessibility features), using ``snapatac2.pp.scrublet``.

Adds ``doublet_score`` and ``doublet_probability`` to ``obs`` of the copy written to ``output_h5ad``. With
``use_selected_features`` the features marked by ``select_atac_features`` are used (they must exist).
``sim_doublet_ratio`` is the number of simulated doublets per cell, ``expected_doublet_rate`` the prior.
Returns the scores for the first ``limit`` cells and the number of cells with probability above 0.5.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.regulation.snapatac2_p38_functions.compute_atac_spectral_embedding`

**Purpose:** Spectral (Laplacian eigenmap) embedding of an accessibility matrix, or one joint embedding of several

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def compute_atac_spectral_embedding(input_h5ads: Sequence[str], output_tsv: str, *, n_comps: int=30, use_selected_features: bool=True, distance_metric: str='cosine', weights: Sequence[float] | None=None, weighted_by_sd: bool=True, random_state: int=0, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Spectral (Laplacian eigenmap) embedding of an accessibility matrix, or one joint embedding of several
modalities, using ``snapatac2.tl.spectral`` (one file) or ``snapatac2.tl.multi_spectral`` (several files).

One file in ``input_h5ads`` gives ``tl.spectral`` with ``distance_metric`` ``cosine`` or ``jaccard``; two or
more files (same cells in the same order, different feature spaces) give ``tl.multi_spectral``, optionally
with per-modality ``weights``. ``use_selected_features`` uses the ``var["selected"]`` mask of
``select_atac_features``. The embedding (cells by components, with barcodes) is written to ``output_tsv``;
the eigenvalues are returned.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.regulation.snapatac2_p38_functions.find_atac_marker_regions`

**Purpose:** Select marker genomic regions per group by z-score enrichment, using ``snapatac2.tl.marker_regions``.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def find_atac_marker_regions(input_h5ad: str, groupby: str, *, pvalue: float=0.01, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Select marker genomic regions per group by z-score enrichment, using ``snapatac2.tl.marker_regions``.

``input_h5ad`` holds an aggregated or cell-level matrix with regions as variables and the grouping column
``groupby`` in ``obs``. ``pvalue`` is the z-score p-value cut-off. Returns the marker regions per group.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.regulation.snapatac2_p38_functions.test_atac_differential_accessibility`

**Purpose:** Differential accessibility between two cell groups over regions, using ``snapatac2.tl.diff_test``

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def test_atac_differential_accessibility(input_h5ad: str, *, group1: Sequence[str], group2: Sequence[str], direction: str='both', min_log_fc: float=0.25, min_pct: float=0.05, output_tsv: str | None=None, limit: int=100, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Differential accessibility between two cell groups over regions, using ``snapatac2.tl.diff_test``
(logistic-regression likelihood-ratio tests).

``group1`` and ``group2`` are lists of cell barcodes. Regions failing ``min_pct`` (detection fraction) or
``min_log_fc`` are not reported. Regions detected in every cell, or in none, of the two groups are dropped
first because the logistic regression is undefined for them (snapATAC2 would stop with an error).
``direction`` is ``positive``, ``negative`` or ``both``. The table (feature name, log2 fold change, p-value,
adjusted p-value) is returned (first ``limit`` rows) and written to ``output_tsv`` when given.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.regulation.snapatac2_p38_functions.score_regulatory_network_edges`

**Purpose:** Score region-gene edges of a regulatory network by correlation or regression, using

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def score_regulatory_network_edges(peak_h5ad: str, gene_h5ad: str, gene_gtf: str, output_tsv: str, *, method: str='correlation', upstream: int=250000, downstream: int=250000, coding_gene_only: bool=True, alpha: float=1.0, l1_ratio: float=0.5, select: Sequence[str] | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Score region-gene edges of a regulatory network by correlation or regression, using
``snapatac2.tl.init_network_from_annotation`` with ``add_cor_scores`` or ``add_regr_scores``.

The network links every region of ``peak_h5ad`` (cell-by-region matrix) to genes whose TSS lies within
``upstream``/``downstream`` bp, using the GTF ``gene_gtf``. ``gene_h5ad`` is the matching cell-by-gene
matrix (same cells in the same order). ``method`` is ``correlation`` (Spearman correlation per edge:
``add_cor_scores``) or a regression importance score, ``elastic_net`` (with ``alpha`` and ``l1_ratio``) or
``gb_tree`` (``add_regr_scores``). ``select`` restricts the target genes. The edge table (source, target,
distance, scores) is written to ``output_tsv``.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.regulation.snapatac2_p38_functions.link_tf_to_gene_network`

**Purpose:** Build a TF-to-gene network from region-gene and TF-region (motif) links, using

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def link_tf_to_gene_network(regions: Sequence[str], gene_gtf: str, motif_file: str, genome_fasta: str, *, pvalue: float=1e-05, upstream: int=250000, downstream: int=250000, coding_gene_only: bool=True, limit: int=1000, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Build a TF-to-gene network from region-gene and TF-region (motif) links, using
``snapatac2.tl.init_network_from_annotation``, ``add_tf_binding`` and ``link_tf_to_gene``.

``regions`` are ``chrom:start-end`` strings. A TF (motif) is linked to a gene when the motif is found in a
region (``pvalue`` cut-off, genome scanned in ``genome_fasta``, motifs from the MEME file) that is within
``upstream``/``downstream`` bp of the gene's TSS in ``gene_gtf``. Returns the TF-to-gene edges.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.regulation.snapatac2_p38_functions.aggregate_atac_profiles`

**Purpose:** Pseudobulk accessibility per group (or all cells), optionally normalised, using ``snapatac2.tl.aggregate_X``.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def aggregate_atac_profiles(input_h5ad: str, output_tsv: str, *, groupby: str | None=None, normalize: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Pseudobulk accessibility per group (or all cells), optionally normalised, using ``snapatac2.tl.aggregate_X``.

``groupby`` is an ``obs`` column (null aggregates all cells into one profile). ``normalize`` is null (summed
counts), ``RPM`` or ``RPKM`` (needs regions named ``chrom:start-end`` as variables). The group-by-region
matrix is written to ``output_tsv`` (rows are groups, columns are regions).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.regulation.snapatac2_p38_functions.read_meme_motifs`

**Purpose:** Read DNA motifs from a MEME-format file, using ``snapatac2.read_motifs``.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def read_meme_motifs(motif_file: str, *, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Read DNA motifs from a MEME-format file, using ``snapatac2.read_motifs``.

Returns, for each motif, its ID, name, family and information content (bits).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.regulation.snapatac2_p38_functions.scan_dna_with_motif`

**Purpose:** Scan DNA sequences with a MEME motif, using ``snapatac2.PyDNAMotifScanner``.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def scan_dna_with_motif(motif_file: str, sequences: Sequence[str], *, motif_id: str | None=None, mode: str='find', pvalue: float=1e-05, reverse_complement: bool=True, report_pvalue: bool=False, base_probabilities: Sequence[float] | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Scan DNA sequences with a MEME motif, using ``snapatac2.PyDNAMotifScanner``.

``mode`` is ``find`` (positions and log-likelihood scores of forward-strand hits above the ``pvalue`` cut-off:
``find``; ``report_pvalue`` adds log10 p-values) or ``exist`` (boolean presence per sequence: ``exist`` /
``exists``; ``reverse_complement`` also tests the opposite strand). ``motif_id`` picks a motif by ID or name
(the first motif by default). ``base_probabilities`` is the background A, C, G, T (default 0.25 each).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.regulation.snapatac2_p38_functions.test_dna_motif_enrichment`

**Purpose:** Motif enrichment (log2 fold change and p-value) of a sequence set against a background set, using

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def test_dna_motif_enrichment(motif_file: str, foreground: Sequence[str], background: Sequence[str], *, motif_id: str | None=None, pvalue: float=1e-05, base_probabilities: Sequence[float] | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Motif enrichment (log2 fold change and p-value) of a sequence set against a background set, using
``snapatac2.PyDNAMotifScanner.with_background`` and ``PyDNAMotifTest.test``.

The motif occurrence rate in ``background`` (sequences containing the motif above the ``pvalue`` cut-off,
either strand) is compared with ``foreground``. Do not build the background from the foreground sequences.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## Executed example

**Example code:**

```python
from pathlib import Path
from dnakit.predictions.selene_sdk_p38_functions import get_selene_genome_window
Path('genome.fa').write_text('>chr1\nACGTACGT\n')
r = get_selene_genome_window('genome.fa', 'chr1', 2, 6, strand='-', verbose=False)
print(r['sequence'], r['length'])
```

**Example output:**

```text
GTAC 4
```
