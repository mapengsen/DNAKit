"""Single-cell and spatial transcriptomics functions: count-matrix I/O, QC, normalisation, markers, clustering and more.

Count matrices are genes x cells (columns are cells or spots). Input forms and options are described in each function's docstring;
the shared parser lives in ``_impl_matrix``. Operations that need a trained model, a reference database or an external program
return ``status="not_implemented"``.
"""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import candidate_not_implemented
from dnakit.singlecell_spatial import _impl_embed as _embed, _impl_sc as _sc, _impl_seq as _seq
from dnakit.singlecell_spatial._impl_matrix import run as _run

def analyze_crispr_perturbations(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Summarise a CRISPR perturbation screen: cells per guide and target-gene knockdown against control cells.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv).
    options: perturbations (required; {cell: guide label} or a list in column order), control ('control'; label of control cells),
    targets ({guide: target gene}; default the guide label is the gene), target_sum (10000), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    Result: per guide n_cells, mean library-size-normalised target expression in guide and control cells, log2fc and
    percent_knockdown = 100 * (1 - mean_guide / mean_control).
    """
    return _run('singlecell_spatial', 'analyze_crispr_perturbations', input, options, _sc.op_analyze_crispr)

def annotate_cells(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Annotate cells with cell-type labels by scoring user-supplied marker genes (no trained model).

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv).
    options: markers (required; {cell type: [marker genes]}), min_score (0.0; cells whose best score is not above it become 'Unassigned'),
    target_sum (10000), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    Score of a type in a cell = mean over its markers of the per-gene z-score (across cells) of ln(1 + CP10k); the highest score wins.
    """
    return _run('singlecell_spatial', 'annotate_cells', input, options, _sc.op_annotate_markers)

def annotate_cells_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Annotate spatial spots with cell-type labels by scoring user-supplied marker genes (no trained model).

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv). Spots are treated as cells; the columns of the matrix are spots.
    options: markers (required; {cell type: [marker genes]}), min_score (0.0), target_sum (10000), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    Same definition as annotate_cells: mean z-score of ln(1 + CP10k) over the markers; the highest score wins.
    """
    return _run('singlecell_spatial', 'annotate_cells_spatial', input, options, _sc.op_annotate_markers)

def build_cell_gene_matrix(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build a gene-by-cell count matrix from per-read or per-UMI (cell, gene) assignment records.

    input: records as TSV/CSV text or a file path (one record per line: cell, gene[, umi]; an optional header may name the columns
    cell/barcode, gene/feature, umi and count), a list of [cell, gene[, umi]] lists, or a list of mappings with those keys.
    options: max_matrix_elements (5000000), return_matrix (default True), output_path (optional TSV file to write; never overwrites).
    When a UMI is present the count is the number of distinct UMIs per (cell, gene) (run correct_umis first to merge sequencing errors);
    otherwise it is the number of records (or the sum of the 'count' column).
    """
    return _run('singlecell_spatial', 'build_cell_gene_matrix', input, options, _seq.op_build_cell_gene_matrix)

def build_cell_peak_matrix(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build a peak-by-cell count matrix from scATAC fragments and a list of peaks.

    input: fragments as TSV text or a file path (chrom, start, end, barcode[, count]; BED-like, 0-based half-open; '#' lines ignored)
    or a list of such lists.
    options: peaks (required; 'chr1:100-600' strings, [chrom, start, end] lists, or BED text), count ('fragments' default = fragments
    overlapping a peak by >= 1 bp, each counted once per peak, or 'cut_sites' = the first and last base of each fragment inside a peak),
    max_matrix_elements (5000000), return_matrix (default True), output_path (optional TSV file to write; never overwrites).
    Result: definition, n_fragments, mode, n_peaks (matrix rows), n_cells, matrix (its 'genes' list holds the peak names) and output_path.
    """
    return _run('singlecell_spatial', 'build_cell_peak_matrix', input, options, _seq.op_build_cell_peak_matrix)

def calculate_cell_cycle_score(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Score S and G2/M cell-cycle phases per cell from user-supplied gene sets and assign G1/S/G2M.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv).
    options: s_genes (required), g2m_genes (required), n_bins (25), ctrl_size (50), target_sum (10000), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    Score = mean ln(1 + CP10k) of the gene set minus the mean of expression-matched control genes (genes in the same mean-expression
    rank bin, nearest rank first, up to ctrl_size per set gene; deterministic, so values differ slightly from scanpy's random controls).
    Phase = G1 if both scores < 0, otherwise S or G2M, whichever score is larger.
    """
    return _run('singlecell_spatial', 'calculate_cell_cycle_score', input, options, _sc.op_cell_cycle)

def calculate_cell_interactions(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Score ligand-receptor interactions between cell types from mean expression (no significance test).

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv).
    options: cell_types (required; {cell: type} or a list in column order), pairs (required; [[ligand, receptor], ...]),
    min_pct (0.0; minimum fraction of sender/receiver cells expressing the gene), top_n (100), target_sum (10000), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    score = mean ln(1 + CP10k) of the ligand in the sender type x mean of the receptor in the receiver type, for every ordered type pair.
    A ligand-receptor database is not bundled; supply the pairs.
    """
    return _run('singlecell_spatial', 'calculate_cell_interactions', input, options, _sc.op_cell_interactions)

def calculate_cell_interactions_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Score ligand-receptor interactions between neighbouring spots/cells within a distance radius.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv). Spots are treated as cells; the columns of the matrix are spots. At most 5000 spots.
    options: cell_types (required), pairs (required; [[ligand, receptor], ...]), radius (required; maximum distance between sender and
    receiver in coordinate units), coordinates (required): {spot: [x, y]}, a list of [x, y] in column order, or a CSV/TSV/Visium tissue_positions file, min_pct (0.0), top_n (100), target_sum (10000), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    score = mean over ordered neighbour pairs (a, b) of different cells within radius of ln(1+CP10k)(ligand in a) x ln(1+CP10k)(receptor in b),
    grouped by (sender type, receiver type). No significance test.
    """
    return _run('singlecell_spatial', 'calculate_cell_interactions_spatial', input, options, lambda v, o: _sc.op_cell_interactions(v, o, spatial=True))

def calculate_markers(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Find marker genes of each group (cluster or cell type) by testing it against all other cells.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv).
    options: groups (required; {cell: label} or a list in column order), only_positive (True; keep genes up in the group), method ('welch' default, or 'wilcoxon'), normalization ('cpm' default = library-size scaling to target_sum, or 'none'), target_sum (10000), pseudocount (1.0), top_n (100), min_pct (0.0), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    Per gene: normalised value = count / library size * target_sum; the test runs on log2(value + pseudocount); log2fc = log2(mean(value, group) + pseudocount) - log2(mean(value, other) + pseudocount) (group = the tested group, other = all remaining cells); pct = fraction of cells with count > 0; q_value = Benjamini-Hochberg over all genes. Welch = unequal-variance t-test; wilcoxon = Mann-Whitney U with normal approximation, tie and continuity correction. Groups too small for Welch's test (fewer than 2 cells in the group or among the other cells) are skipped with a warning and the call only raises GROUPS_TOO_SMALL when no group could be tested; find_differential_features and, in transcriptomics, calculate_differential_expression instead raise InputFormatError (GROUP_TOO_SMALL) for a too-small group.
    """
    return _run('singlecell_spatial', 'calculate_markers', input, options, _sc.op_calculate_markers)

def calculate_markers_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Find marker genes of each spatial domain or spot cluster by testing it against all other spots.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv). Spots are treated as cells; the columns of the matrix are spots.
    options: groups (required; {spot: domain label} or a list in column order), only_positive (True), method ('welch' default, or 'wilcoxon'), normalization ('cpm' default = library-size scaling to target_sum, or 'none'), target_sum (10000), pseudocount (1.0), top_n (100), min_pct (0.0), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    Per gene: normalised value = count / library size * target_sum; the test runs on log2(value + pseudocount); log2fc = log2(mean(value, group) + pseudocount) - log2(mean(value, other) + pseudocount) (group = the tested domain, other = all remaining spots); pct = fraction of cells with count > 0; q_value = Benjamini-Hochberg over all genes. Welch = unequal-variance t-test; wilcoxon = Mann-Whitney U with normal approximation, tie and continuity correction. Domains too small for Welch's test are skipped with a warning (GROUPS_TOO_SMALL is raised only when none could be tested), unlike find_differential_features_spatial which raises GROUP_TOO_SMALL.
    """
    return _run('singlecell_spatial', 'calculate_markers_spatial', input, options, _sc.op_calculate_markers)

def calculate_perturbation_effect(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Test the transcriptome-wide effect of each perturbation against control cells.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv).
    options: perturbations (required; {cell: perturbation label} or a list), control ('control'), method ('welch' default, or 'wilcoxon'), normalization ('cpm' default = library-size scaling to target_sum, or 'none'), target_sum (10000), pseudocount (1.0), top_n (100), min_pct (0.0), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    Per gene: normalised value = count / library size * target_sum; the test runs on log2(value + pseudocount); log2fc = log2(mean(value, perturbation) + pseudocount) - log2(mean(value, control) + pseudocount) (each perturbation is group1 and the control cells are group2: *_group1 fields describe the perturbation, *_group2 the control; positive log2fc = up in the perturbed cells); pct = fraction of cells with count > 0; q_value = Benjamini-Hochberg over all genes. Perturbations with fewer than 2 cells (or fewer than 2 control cells) are skipped with a warning under Welch, and an error is raised only when none can be tested. Welch = unequal-variance t-test; wilcoxon = Mann-Whitney U with normal approximation, tie and continuity correction.
    """
    return _run('singlecell_spatial', 'calculate_perturbation_effect', input, options, _sc.op_perturbation_effect)

def cluster(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Cluster cells with k-means on the top principal components of log-normalised expression.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv). Requires numpy; at most 20000 cells (5000 with n_clusters='auto').
    options: n_clusters ('auto' default = k in 2..max_clusters with the best mean silhouette, or an integer), max_clusters (10), n_pcs (10; 0 = no PCA),
    n_top_genes (None = all; otherwise the highest-variance genes), scale (False), target_sum (10000), return_embedding (False), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    Values are ln(1 + CP10k). k-means uses a deterministic farthest-point start and Lloyd iterations; labels are numbered by first appearance.
    """
    return _run('singlecell_spatial', 'cluster', input, options, _embed.op_cluster)

def cluster_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Cluster spots into spatial domains with k-means on spatially smoothed principal components.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv). Spots are treated as cells; the columns of the matrix are spots. Requires numpy; at most 5000 spots.
    options: coordinates (required): {spot: [x, y]}, a list of [x, y] in column order, or a CSV/TSV/Visium tissue_positions file, n_neighbors (6), alpha (0.5), n_clusters ('auto' or integer), max_clusters (10), n_pcs (10), n_top_genes, scale, target_sum (10000), return_embedding, orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    Heuristic: each spot's PC vector becomes (1 - alpha) * own + alpha * mean of its n_neighbors nearest spatial neighbours, then k-means as in cluster().
    """
    return _run('singlecell_spatial', 'cluster_spatial', input, options, _embed.op_cluster_spatial)

def correct_umis(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Correct sequencing errors in UMIs by collapsing abundant-neighbour groups (UMI-tools directional/adjacency/cluster).

    input: a list of UMI strings (each occurrence counts as one read), {UMI: read count}, or text/file with one 'UMI' or 'UMI<TAB>count' per line.
    All UMIs must have the same length and use A, C, G, T, N.
    options: method ('directional' default, 'adjacency', 'cluster' or 'unique'), threshold (1; maximum Hamming distance, 1-4).
    directional: from the most abundant UMI absorb neighbours b with count_a >= 2 * count_b - 1, transitively; adjacency: absorb all unassigned
    neighbours one step; cluster: connected components. Returns the UMI -> representative mapping and unique-UMI counts before/after.
    """
    return _run('singlecell_spatial', 'correct_umis', input, options, _seq.op_correct_umis)

def correct_umis_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Correct sequencing errors in spatial-barcoded UMIs by collapsing abundant-neighbour groups (UMI-tools algorithms).

    input: a list of UMI strings, {UMI: read count}, or text/file with one 'UMI' or 'UMI<TAB>count' per line (UMIs of one spot and gene).
    options: method ('directional' default, 'adjacency', 'cluster' or 'unique'), threshold (1; maximum Hamming distance, 1-4).
    Same algorithms as correct_umis.
    """
    return _run('singlecell_spatial', 'correct_umis_spatial', input, options, _seq.op_correct_umis)

def deconvolute_spots(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Estimate cell-type fractions of each spot by non-negative least squares on a signature matrix.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv). Requires numpy and scipy.
    options: signature (genes x cell types matrix in the same input forms) OR reference (count matrix of annotated cells) + reference_labels
    ({cell: type} or a list); the signature is then the mean library-size-normalised profile per type. normalize (True: scale every column to sum 1), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    Fractions = NNLS coefficients / their sum over the genes shared with the signature; rmse is the fit residual per spot.
    """
    return _run('singlecell_spatial', 'deconvolute_spots', input, options, _sc.op_deconvolve)

def deconvolute_spots_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Estimate cell-type fractions of each spatial spot by non-negative least squares on a signature matrix.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv). Spots are treated as cells; the columns of the matrix are spots. Requires numpy and scipy.
    options: signature OR reference + reference_labels, normalize (True), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000). Same method as deconvolute_spots.
    """
    return _run('singlecell_spatial', 'deconvolute_spots_spatial', input, options, _sc.op_deconvolve)

def demultiplex(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Assign cells to samples from hashtag-oligo (HTO) counts with a simple count-threshold rule.

    input: a hashtag-by-cell count matrix (a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv)); rows are hashtags.
    options: min_counts (5.0), ratio (3.0), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    Negative if the top hashtag count < min_counts; Doublet if the second count >= min_counts and top < ratio * second; otherwise Singlet of the top hashtag.
    A threshold heuristic, not the mixture-model fit of HTODemux.
    """
    return _run('singlecell_spatial', 'demultiplex', input, options, _sc.op_demultiplex)

def demultiplex_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Assign spots to samples from hashtag-oligo counts with a simple count-threshold rule.

    input: a hashtag-by-spot count matrix (a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv)); rows are hashtags. Spots are treated as cells; the columns of the matrix are spots.
    options: min_counts (5.0), ratio (3.0), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000). Same rule as demultiplex (heuristic).
    """
    return _run('singlecell_spatial', 'demultiplex_spatial', input, options, _sc.op_demultiplex)

def extract_barcodes(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Extract the cell barcode (a fixed read window) from reads and count barcodes, optionally against a whitelist.

    input: FASTQ/FASTA text or a file path, a list of sequences, or {read name: sequence} (read 1 of a droplet protocol).
    options: barcode_start (0), barcode_length (16), whitelist (list/text/path of valid barcodes; exact match or unique 1-mismatch correction),
    top_n (100).
    Barcode = read[start:start+length]; reads shorter than the window are counted as too short.
    Result counts: n_reads (all reads) = n_extracted (barcodes counted; with a whitelist only exact or corrected matches)
    + n_unmatched (whitelist given but no exact or unique 1-mismatch match; 0 without a whitelist) + n_too_short.
    """
    return _run('singlecell_spatial', 'extract_barcodes', input, options, _seq.op_extract_barcodes)

def extract_barcodes_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Extract the spatial (spot) barcode from reads and count barcodes, optionally against a whitelist.

    input: FASTQ/FASTA text or a file path, a list of sequences, or {read name: sequence} (e.g. Visium read 1).
    options: barcode_start (0), barcode_length (16), whitelist (list/text/path; exact or unique 1-mismatch match), top_n (100).
    Result counts: n_reads = n_extracted (spot barcodes counted, whitelist matches only) + n_unmatched + n_too_short.
    """
    return _run('singlecell_spatial', 'extract_barcodes_spatial', input, options, _seq.op_extract_barcodes)

def extract_umis(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Extract the UMI (a fixed read window) from reads and count UMIs.

    input: FASTQ/FASTA text or a file path, a list of sequences, or {read name: sequence} (read 1 of a droplet protocol).
    options: umi_start (16), umi_length (12), top_n (100).
    UMI = read[start:start+length]; the defaults are the 10x Chromium v3 layout (16 bp barcode + 12 bp UMI).
    """
    return _run('singlecell_spatial', 'extract_umis', input, options, _seq.op_extract_umis)

def extract_umis_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Extract the UMI (a fixed read window) from spatial reads and count UMIs.

    input: FASTQ/FASTA text or a file path, a list of sequences, or {read name: sequence} (e.g. Visium read 1).
    options: umi_start (16), umi_length (12), top_n (100). UMI = read[start:start+length].
    """
    return _run('singlecell_spatial', 'extract_umis_spatial', input, options, _seq.op_extract_umis)

def find_differential_features(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Test differential expression between two groups of cells (Welch t-test or Wilcoxon rank-sum).

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv).
    options: groups (required; {cell: label} or a list), group1 and group2 (labels to compare; optional when there are exactly two labels), method ('welch' default, or 'wilcoxon'), normalization ('cpm' default = library-size scaling to target_sum, or 'none'), target_sum (10000), pseudocount (1.0), top_n (100), min_pct (0.0), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    Per gene: normalised value = count / library size * target_sum; the test runs on log2(value + pseudocount); log2fc = log2(mean(value, group1) + pseudocount) - log2(mean(value, group2) + pseudocount) (group1 versus group2, not group versus all other cells; the result fields mean_*/pct_* carry the _group1/_group2 suffix); pct = fraction of cells with count > 0; q_value = Benjamini-Hochberg over all genes. Welch = unequal-variance t-test; wilcoxon = Mann-Whitney U with normal approximation, tie and continuity correction. With method='welch' a group with fewer than 2 cells raises InputFormatError (GROUP_TOO_SMALL); calculate_markers skips such groups with a warning instead.
    """
    return _run('singlecell_spatial', 'find_differential_features', input, options, _sc.op_find_differential)

def find_differential_features_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Test differential expression between two groups of spots or spatial domains.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv). Spots are treated as cells; the columns of the matrix are spots.
    options: groups (required), group1 and group2 (optional when there are exactly two labels), method ('welch' default, or 'wilcoxon'), normalization ('cpm' default = library-size scaling to target_sum, or 'none'), target_sum (10000), pseudocount (1.0), top_n (100), min_pct (0.0), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    Per gene: normalised value = count / library size * target_sum; the test runs on log2(value + pseudocount); log2fc = log2(mean(value, group1) + pseudocount) - log2(mean(value, group2) + pseudocount) (group1 versus group2, not group versus all other spots; the result fields mean_*/pct_* carry the _group1/_group2 suffix); pct = fraction of cells with count > 0; q_value = Benjamini-Hochberg over all genes. Welch = unequal-variance t-test; wilcoxon = Mann-Whitney U with normal approximation, tie and continuity correction. With method='welch' a group with fewer than 2 spots raises InputFormatError (GROUP_TOO_SMALL); calculate_markers_spatial skips such domains with a warning instead.
    """
    return _run('singlecell_spatial', 'find_differential_features_spatial', input, options, _sc.op_find_differential)

def infer_cell_type_from_markers(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Infer the cell type of each cell from user-supplied marker genes by marker-score ranking.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv).
    options: markers (required; {cell type: [marker genes]}), min_score (0.0; below it a cell is 'Unassigned'), target_sum (10000), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    Score of a type in a cell = mean over its markers of the per-gene z-score (across cells) of ln(1 + CP10k); the highest score wins.
    """
    return _run('singlecell_spatial', 'infer_cell_type_from_markers', input, options, _sc.op_annotate_markers)

def infer_cnv(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Infer copy-number-like expression profiles along chromosomes relative to reference cells (inferCNV-style, no HMM).

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv).
    options: gene_positions (required; {gene: [chromosome, start]}), groups (required; {cell: label}), reference_group (required; label of normal cells),
    window (101 genes), min_mean (0.1; minimum mean raw count of a gene), target_sum (1000000), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000), return_matrix (default True), output_path (optional TSV file to write; never overwrites).
    modified expression = log2(CPM + 1) minus the reference-cell mean per gene, smoothed by a moving average over `window` genes within each
    chromosome, then median-centred per cell. cnv_score = mean absolute modified expression. It does not call CNV segments.
    """
    return _run('singlecell_spatial', 'infer_cnv', input, options, _sc.op_infer_cnv)

def infer_cnv_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Infer copy-number-like expression profiles along chromosomes for spots relative to reference spots (inferCNV-style).

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv). Spots are treated as cells; the columns of the matrix are spots.
    options: gene_positions (required), groups (required; {spot: label}), reference_group (required), window (101), min_mean (0.1), target_sum (1000000), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000), return_matrix (default True), output_path (optional TSV file to write; never overwrites).
    Same definition as infer_cnv.
    """
    return _run('singlecell_spatial', 'infer_cnv_spatial', input, options, _sc.op_infer_cnv)

def infer_trajectory(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compute geodesic pseudotime from a root cell through a k-nearest-neighbour graph of expression space.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv). Requires numpy; at most 5000 cells.
    options: root (required; cell name), n_neighbors (15), n_pcs (0 = use the centred ln(1+CP10k) vectors; >0 = PCA first), n_top_genes, scale, target_sum (10000), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    pseudotime = shortest-path distance from the root in the symmetric kNN graph (edge weight = Euclidean distance), divided by the largest finite distance.
    Disconnected cells get null.
    """
    return _run('singlecell_spatial', 'infer_trajectory', input, options, _embed.op_trajectory)

def infer_trajectory_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compute pseudotime for spots along the tissue: shortest path from a root spot over a spatial neighbour graph.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv). Spots are treated as cells; the columns of the matrix are spots. Requires numpy; at most 5000 spots.
    options: root (required; spot name), coordinates (required): {spot: [x, y]}, a list of [x, y] in column order, or a CSV/TSV/Visium tissue_positions file, n_neighbors (6; spatial neighbours), n_pcs (0), n_top_genes, scale, target_sum (10000), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    Graph = symmetric kNN on the coordinates; edge weight = Euclidean distance between the spots' ln(1+CP10k) expression vectors;
    pseudotime = shortest-path length from the root divided by the maximum.
    """
    return _run('singlecell_spatial', 'infer_trajectory_spatial', input, options, _embed.op_trajectory_spatial)

def infer_velocity(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Estimate RNA velocity per gene and cell with the steady-state model from spliced and unspliced counts.

    input: spliced count matrix (a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv)).
    options: unspliced (required; matrix with the same genes and cells in the same order), quantile (0.05), normalize (False; scale each layer to 10000 per cell),
    orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000), return_matrix (default True), output_path (optional TSV file to write; never overwrites).
    gamma_g = sum(u * s) / sum(s^2) over the cells in the top `quantile` of s or u for that gene (least squares through the origin, La Manno et al. 2018);
    velocity = u - gamma * s. Pass normalised layers (or normalize=True); no embedding projection is done.
    """
    return _run('singlecell_spatial', 'infer_velocity', input, options, _sc.op_velocity)

def infer_velocity_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Estimate RNA velocity per gene and spot with the steady-state model from spliced and unspliced counts.

    input: spliced count matrix (a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv)). Spots are treated as cells; the columns of the matrix are spots.
    options: unspliced (required), quantile (0.05), normalize (False), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000), return_matrix (default True), output_path (optional TSV file to write; never overwrites). Same model as infer_velocity.
    """
    return _run('singlecell_spatial', 'infer_velocity_spatial', input, options, _sc.op_velocity)

def integrate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Integrate batches of cells by per-gene batch mean-centring of log-normalised expression.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv).
    options: batches (required; {cell: batch label} or a list in column order), log_normalize (True: work on ln(1 + CP10k)), target_sum (10000),
    orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000), return_matrix (default True), output_path (optional TSV file to write; never overwrites).
    Per-gene batch mean-centring: x' = x - mean(gene in batch) + mean(gene in all cells). It aligns batch means only; it is not Harmony, scVI or MNN.
    """
    return _run('singlecell_spatial', 'integrate', input, options, _sc.op_integrate)

def integrate_single_cell_batches(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Integrate single-cell batches by per-gene batch mean-centring of log-normalised expression.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv).
    options: batches (required; {cell: batch label} or a list in column order), log_normalize (True: work on ln(1 + CP10k)), target_sum (10000),
    orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000), return_matrix (default True), output_path (optional TSV file to write; never overwrites).
    Per-gene batch mean-centring: x' = x - mean(gene in batch) + mean(gene in all cells). It aligns batch means only; it is not Harmony, scVI or MNN.
    """
    return _run('singlecell_spatial', 'integrate_single_cell_batches', input, options, _sc.op_integrate)

def integrate_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Integrate spatial slices by per-gene batch mean-centring of log-normalised expression.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv).
    options: batches (required; {spot: batch label} or a list in column order), log_normalize (True: work on ln(1 + CP10k)), target_sum (10000),
    orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000), return_matrix (default True), output_path (optional TSV file to write; never overwrites).
    Per-gene batch mean-centring: x' = x - mean(gene in batch) + mean(gene in all cells). It aligns batch means only; it is not Harmony, scVI or MNN. Spots are treated as cells; the columns of the matrix are spots.
    """
    return _run('singlecell_spatial', 'integrate_spatial', input, options, _sc.op_integrate)

def map_cells(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Map query cells onto an annotated reference by k-nearest-neighbour label transfer.

    input: query a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv). Requires numpy; query and reference have at most 5000 cells each.
    options: reference (required; annotated reference count matrix), reference_labels (required; {cell: label} or a list), n_neighbors (5),
    target_sum (10000), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    k-nearest-neighbour label transfer: Euclidean distance on ln(1 + CP10k) of the genes shared with the reference; majority vote of the k nearest reference cells
    (ties go to the label of the nearest tied neighbour); confidence = winning votes / k.
    """
    return _run('singlecell_spatial', 'map_cells', input, options, _embed.op_map_cells)

def map_cells_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Map spatial spots onto an annotated single-cell reference by k-nearest-neighbour label transfer.

    input: query a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv). Requires numpy; query and reference have at most 5000 cells each.
    options: reference (required; annotated reference count matrix), reference_labels (required; {cell: label} or a list), n_neighbors (5),
    target_sum (10000), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    k-nearest-neighbour label transfer: Euclidean distance on ln(1 + CP10k) of the genes shared with the reference; majority vote of the k nearest reference cells
    (ties go to the label of the nearest tied neighbour); confidence = winning votes / k. Spots are treated as cells; the columns of the matrix are spots.
    """
    return _run('singlecell_spatial', 'map_cells_spatial', input, options, _embed.op_map_cells)

def map_single_cell_to_reference(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Map single cells onto an annotated reference atlas by k-nearest-neighbour label transfer.

    input: query a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv). Requires numpy; query and reference have at most 5000 cells each.
    options: reference (required; annotated reference count matrix), reference_labels (required; {cell: label} or a list), n_neighbors (5),
    target_sum (10000), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    k-nearest-neighbour label transfer: Euclidean distance on ln(1 + CP10k) of the genes shared with the reference; majority vote of the k nearest reference cells
    (ties go to the label of the nearest tied neighbour); confidence = winning votes / k.
    """
    return _run('singlecell_spatial', 'map_single_cell_to_reference', input, options, _embed.op_map_cells)

def normalize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Normalise a count matrix by library size, optionally followed by ln(1 + x).

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv).
    options: method ('log1p' default = scale each cell to target_sum then ln(1 + x); 'total' = scale only; 'cpm' = scale to 1000000), target_sum (10000), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000), return_matrix (default True), output_path (optional TSV file to write; never overwrites).
    Cells with zero total counts stay all-zero.
    """
    return _run('singlecell_spatial', 'normalize', input, options, _sc.op_normalize)

def normalize_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Normalise a spot-by-gene count matrix by library size, optionally followed by ln(1 + x).

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv). Spots are treated as cells; the columns of the matrix are spots.
    options: method ('log1p' default, 'total' or 'cpm'), target_sum (10000), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000), return_matrix (default True), output_path (optional TSV file to write; never overwrites). Same definition as normalize.
    """
    return _run('singlecell_spatial', 'normalize_spatial', input, options, _sc.op_normalize)

def plot_embedding(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Draw a 2-D embedding (UMAP, t-SNE, PCA coordinates) as an SVG scatter plot.

    input: coordinates with one row per cell: {cell: [x, y]}, a list of [x, y], CSV/TSV text or file (header = dimension names, first column = cell), a numpy array.
    options: labels ({cell: label} or a list; categorical colours + legend), values ({cell: number} or a list; continuous colour ramp), output_path (optional .svg file to write; never overwrites), title.
    Returns the SVG text under result['svg'].
    """
    return _run('singlecell_spatial', 'plot_embedding', input, options, _sc.op_plot_embedding)

def plot_embedding_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Draw spot coordinates of a spatial dataset as an SVG scatter plot coloured by label or value.

    input: spot coordinates with one row per spot: {spot: [x, y]}, a list of [x, y], CSV/TSV text or file, a numpy array (y axis points down like an image).
    options: labels, values, output_path (optional .svg file to write; never overwrites), title. Returns the SVG text under result['svg'].
    """
    return _run('singlecell_spatial', 'plot_embedding_spatial', input, options, lambda v, o: _sc.op_plot_embedding(v, o, spatial=True))

def plot_heatmap(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Draw a gene-by-cell expression heatmap as an SVG.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv).
    options: max_genes (50; genes with the highest variance are shown), max_cells (200; first cells shown), log1p (True), row_zscore (False), output_path (optional .svg file to write; never overwrites), title, orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    """
    return _run('singlecell_spatial', 'plot_heatmap', input, options, _sc.op_plot_heatmap)

def plot_heatmap_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Draw the spatial expression map of one gene (spots coloured by ln(1 + CP10k)) as an SVG.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv). Spots are treated as cells; the columns of the matrix are spots.
    options: gene (required), coordinates (required): {spot: [x, y]}, a list of [x, y] in column order, or a CSV/TSV/Visium tissue_positions file, target_sum (10000), output_path (optional .svg file to write; never overwrites), title, orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    """
    return _run('singlecell_spatial', 'plot_heatmap_spatial', input, options, _sc.op_plot_spatial_gene)

def qc(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compute per-cell and per-gene quality-control metrics of a count matrix and flag cells failing thresholds.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv).
    options: min_genes, min_counts, max_counts, max_pct_mito, min_cells (all default None = report only), mito_prefix ('MT-' and 'mt-'), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    Per cell: total counts, genes detected (count > 0) and percent mitochondrial counts; lists cells failing the given thresholds.
    """
    return _run('singlecell_spatial', 'qc', input, options, _sc.op_qc)

def qc_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compute per-spot and per-gene quality-control metrics of a spatial count matrix and flag spots failing thresholds.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv). Spots are treated as cells; the columns of the matrix are spots.
    options: min_genes, min_counts, max_counts, max_pct_mito, min_cells (all default None = report only), mito_prefix ('MT-' and 'mt-'), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    Per cell: total counts, genes detected (count > 0) and percent mitochondrial counts; lists cells failing the given thresholds.
    Optional coordinates (same forms as other spatial functions) add the spot count and coordinate ranges.
    """
    return _run('singlecell_spatial', 'qc_spatial', input, options, lambda v, o: _sc.op_qc(v, o, spatial=True))

def read(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Read a count matrix (mapping, rows, CSV/TSV, array or 10x directory) and summarise its shape and sparsity.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv).
    options: max_return_elements (10000; the matrix values are returned when it has at most this many elements), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    """
    return _run('singlecell_spatial', 'read', input, options, _sc.op_read)

def read_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Read a spatial count matrix together with spot coordinates and summarise it.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv). Spots are treated as cells; the columns of the matrix are spots.
    options: coordinates (optional; {spot: [x, y]}, a list of [x, y], or a CSV/TSV/Visium tissue_positions file), x_column and y_column (header names), max_return_elements (10000), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    """
    return _run('singlecell_spatial', 'read_spatial', input, options, lambda v, o: _sc.op_read(v, o, spatial=True))

def remove_ambient(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Subtract ambient-RNA counts from a droplet count matrix using an ambient profile and a contamination fraction.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv).
    options: contamination (required; fraction in [0, 1) of each cell's counts that is ambient), ambient_profile ({gene: weight}) OR empty_cells (list of empty-droplet
    columns) OR empty_max_counts (columns with total <= this are empty droplets and are removed from the output), round (False), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000), return_matrix (default True), output_path (optional TSV file to write; never overwrites).
    corrected = max(0, count - contamination * cell_total * ambient_fraction(gene)); a simplified subtractive, SoupX-style correction with a user-given contamination
    fraction (no contamination estimation, no integer redistribution).
    """
    return _run('singlecell_spatial', 'remove_ambient', input, options, _sc.op_remove_ambient)

def remove_ambient_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Subtract ambient-RNA counts from a spatial count matrix using an ambient profile and a contamination fraction.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv). Spots are treated as cells; the columns of the matrix are spots.
    options: contamination (required; fraction in [0, 1) of each cell's counts that is ambient), ambient_profile ({gene: weight}) OR empty_cells (list of empty-droplet
    columns) OR empty_max_counts (columns with total <= this are empty droplets and are removed from the output), round (False), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000), return_matrix (default True), output_path (optional TSV file to write; never overwrites).
    corrected = max(0, count - contamination * cell_total * ambient_fraction(gene)); a simplified subtractive, SoupX-style correction with a user-given contamination
    fraction (no contamination estimation, no integer redistribution).
    """
    return _run('singlecell_spatial', 'remove_ambient_spatial', input, options, _sc.op_remove_ambient)

def remove_doublets(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Remove doublets from a count matrix (not implemented).

    Doublet detection needs a simulation-and-classifier implementation (for example Scrublet or DoubletFinder): artificial doublets are simulated and a trained classifier or neighbour model scores each cell. DNAKit does not call any such package, so installing one changes nothing; no deterministic reference definition is bundled.
    Returns status "not_implemented" for every input.
    """
    return candidate_not_implemented('singlecell_spatial', 'remove_doublets', input, options=options, requires="a doublet-detection implementation based on simulation and a classifier (for example Scrublet or DoubletFinder); DNAKit does not call any such package", message="Doublet detection needs a simulation-and-classifier implementation (for example Scrublet or DoubletFinder) which DNAKit does not call; no deterministic reference definition is bundled.")

def remove_doublets_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Remove doublets from a count matrix (not implemented).

    Spot-level doublet detection needs a simulation-and-classifier implementation (for example Scrublet or DoubletFinder) which DNAKit does not call, so installing one changes nothing; none is bundled.
    Returns status "not_implemented" for every input.
    """
    return candidate_not_implemented('singlecell_spatial', 'remove_doublets_spatial', input, options=options, requires="a doublet-detection implementation based on simulation and a classifier (for example Scrublet or DoubletFinder); DNAKit does not call any such package", message="Spot-level doublet detection needs a simulation-and-classifier implementation (for example Scrublet or DoubletFinder) which DNAKit does not call; none is bundled.")

def validate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Validate a count matrix: shape, numeric values, non-negativity, integer counts, unique names, empty cells/genes.

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv).
    options: require_integer (True; non-integer values are errors), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    Returns valid, errors and warnings; a matrix that cannot be parsed is reported as invalid.
    """
    return _run('singlecell_spatial', 'validate', input, options, _sc.op_validate)

def validate_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Validate a spatial count matrix and its spot coordinates (every spot needs finite coordinates).

    input: a gene-by-cell count matrix: {gene: [counts per cell]}, {gene: {cell: count}}, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array / DataFrame, or a 10x directory (matrix.mtx, features.tsv, barcodes.tsv). Spots are treated as cells; the columns of the matrix are spots.
    options: coordinates (required): {spot: [x, y]}, a list of [x, y] in column order, or a CSV/TSV/Visium tissue_positions file, require_integer (True), orientation ('genes_by_cells' default, or 'cells_by_genes' when rows are cells), genes and cells (name lists for unnamed forms), max_matrix_elements (default 5000000).
    Returns valid, errors and warnings.
    """
    return _run('singlecell_spatial', 'validate_spatial', input, options, lambda v, o: _sc.op_validate(v, o, spatial=True))

__all__ = ['analyze_crispr_perturbations', 'annotate_cells', 'annotate_cells_spatial', 'build_cell_gene_matrix', 'build_cell_peak_matrix', 'calculate_cell_cycle_score', 'calculate_cell_interactions', 'calculate_cell_interactions_spatial', 'calculate_markers', 'calculate_markers_spatial', 'calculate_perturbation_effect', 'cluster', 'cluster_spatial', 'correct_umis', 'correct_umis_spatial', 'deconvolute_spots', 'deconvolute_spots_spatial', 'demultiplex', 'demultiplex_spatial', 'extract_barcodes', 'extract_barcodes_spatial', 'extract_umis', 'extract_umis_spatial', 'find_differential_features', 'find_differential_features_spatial', 'infer_cell_type_from_markers', 'infer_cnv', 'infer_cnv_spatial', 'infer_trajectory', 'infer_trajectory_spatial', 'infer_velocity', 'infer_velocity_spatial', 'integrate', 'integrate_single_cell_batches', 'integrate_spatial', 'map_cells', 'map_cells_spatial', 'map_single_cell_to_reference', 'normalize', 'normalize_spatial', 'plot_embedding', 'plot_embedding_spatial', 'plot_heatmap', 'plot_heatmap_spatial', 'qc', 'qc_spatial', 'read', 'read_spatial', 'remove_ambient', 'remove_ambient_spatial', 'remove_doublets', 'remove_doublets_spatial', 'validate', 'validate_spatial']

__candidate_adapter_operations__ = set(['analyze_crispr_perturbations', 'annotate_cells', 'annotate_cells_spatial', 'build_cell_gene_matrix', 'build_cell_peak_matrix', 'calculate_cell_cycle_score', 'calculate_cell_interactions', 'calculate_cell_interactions_spatial', 'calculate_markers', 'calculate_markers_spatial', 'calculate_perturbation_effect', 'cluster', 'cluster_spatial', 'correct_umis', 'correct_umis_spatial', 'deconvolute_spots', 'deconvolute_spots_spatial', 'demultiplex', 'demultiplex_spatial', 'extract_barcodes', 'extract_barcodes_spatial', 'extract_umis', 'extract_umis_spatial', 'find_differential_features', 'find_differential_features_spatial', 'infer_cell_type_from_markers', 'infer_cnv', 'infer_cnv_spatial', 'infer_trajectory', 'infer_trajectory_spatial', 'infer_velocity', 'infer_velocity_spatial', 'integrate', 'integrate_single_cell_batches', 'integrate_spatial', 'map_cells', 'map_cells_spatial', 'map_single_cell_to_reference', 'normalize', 'normalize_spatial', 'plot_embedding', 'plot_embedding_spatial', 'plot_heatmap', 'plot_heatmap_spatial', 'qc', 'qc_spatial', 'read', 'read_spatial', 'remove_ambient', 'remove_ambient_spatial', 'remove_doublets', 'remove_doublets_spatial', 'validate', 'validate_spatial'])


from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
