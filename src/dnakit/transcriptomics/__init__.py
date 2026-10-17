"""Transcriptomics functions: expression-matrix normalisation, QC, differential expression/splicing, annotation and RNA-seq tool runners.

Count matrices are genes x samples. Operations that need a trained model, a reference database or an external program
return ``status="not_implemented"`` (programs are only run when found on PATH).
"""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import candidate_not_implemented
from dnakit.singlecell_spatial import _impl_sc as _sc
from dnakit.singlecell_spatial._impl_matrix import run as _run
from dnakit.transcriptomics import _impl_tx as _tx

def align(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Align RNA-seq reads to a genome by running STAR or HISAT2 (an external program).

    input: a FASTQ path (single-end) or a list of two FASTQ paths (paired-end); .gz is accepted for STAR.
    options: index (required; STAR genome directory or HISAT2 index prefix), output_prefix (required; must not exist), aligner ('auto' = STAR if on PATH else hisat2, 'star' or 'hisat2'), threads (1), timeout_seconds (3600).
    Returns status "not_implemented" when neither program is on PATH. STAR writes <prefix>Aligned.sortedByCoord.out.bam, HISAT2 writes <prefix>.sam.
    """
    return _run('transcriptomics', 'align', input, options, _tx.op_align)

def annotate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Annotate gene identifiers (Ensembl IDs or symbols) with coordinates, strand and biotype from a GTF file.

    input: a list of gene IDs or symbols, text with one per line, a file path, or a mapping/count table keyed by gene.
    options: gtf (required; path to a .gtf or .gtf.gz that contains 'gene' features).
    An identifier matches gene_id (with or without the .version suffix) or gene_name; coordinates are 1-based inclusive as in GTF.
    """
    return _run('transcriptomics', 'annotate', input, options, _tx.op_annotate)

def calculate_differential_expression(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Test differential gene expression between two sample groups (Welch t-test or Wilcoxon rank-sum on log2 CPM).

    input: a gene-by-sample count matrix: {gene: [counts per sample]}, {gene: {sample: count}}, a list of rows, CSV/TSV text or a file path (header = sample names, first column = gene names), a numpy array / DataFrame, or a 10x directory.
    options: groups (required; {sample: condition} or a list in column order), group1 and group2 (labels; optional when there are exactly two), method ('welch' default, or 'wilcoxon'), normalization ('cpm' default = scale each sample to target_sum, or 'none' = test the supplied values), target_sum (1000000), pseudocount (1.0), top_n (100), min_pct (0.0; minimum fraction of samples with count > 0 in either group), orientation ('genes_by_cells' default = rows are genes, or 'cells_by_genes' when rows are samples), genes and cells (name lists for unnamed forms; 'cells' names the samples), max_matrix_elements (default 5000000).
    value = count / library size * target_sum; test on log2(value + pseudocount); log2fc = log2(mean(group1) + pseudocount) - log2(mean(group2) + pseudocount) (group1 versus group2);
    q_value = Benjamini-Hochberg over all genes. This is not a negative-binomial model (DESeq2/edgeR) and needs at least 2 samples per group for Welch's test: a smaller group raises InputFormatError (GROUP_TOO_SMALL) instead of being skipped with a warning (singlecell_spatial.calculate_markers skips too-small groups, this function and find_differential_features raise).
    """
    return _run('transcriptomics', 'calculate_differential_expression', input, options, _tx.op_differential_expression)

def calculate_differential_splicing(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Test differential exon inclusion (percent spliced in, PSI) between two groups from junction read counts.

    input: inclusion-junction count matrix of splicing events x samples (a gene-by-sample count matrix: {gene: [counts per sample]}, {gene: {sample: count}}, a list of rows, CSV/TSV text or a file path (header = sample names, first column = gene names), a numpy array / DataFrame, or a 10x directory).
    options: exclusion (required; exclusion-junction counts with identical events and samples), groups (required), group1 and group2 (optional when there are exactly two labels),
    inclusion_length (1.0) and exclusion_length (1.0; effective lengths used to normalise the read counts), min_reads (10; pooled reads required in each group), top_n (100), orientation ('genes_by_cells' default = rows are genes, or 'cells_by_genes' when rows are samples), genes and cells (name lists for unnamed forms; 'cells' names the samples), max_matrix_elements (default 5000000).
    PSI = (inc/Li) / (inc/Li + exc/Le) per sample; psi_group = mean over samples; delta_psi = psi_group1 - psi_group2; p_value = two-sided Fisher exact test on the pooled
    inclusion/exclusion counts of the groups (ignores replicate variance, so it is anti-conservative); q_value = Benjamini-Hochberg.
    """
    return _run('transcriptomics', 'calculate_differential_splicing', input, options, _tx.op_differential_splicing)

def calculate_translation_efficiency(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Calculate translation efficiency (ribosome-profiling abundance divided by RNA-seq abundance) per gene.

    input: ribosome-footprint count matrix, genes x samples (a gene-by-sample count matrix: {gene: [counts per sample]}, {gene: {sample: count}}, a list of rows, CSV/TSV text or a file path (header = sample names, first column = gene names), a numpy array / DataFrame, or a 10x directory).
    options: rna (required; RNA-seq count matrix of the same genes; samples are paired by identical names, otherwise by position and reported under the ribo-seq names), gene_lengths ({gene: bp}; only relabels the abundances RPKM instead of CPM: the same length divides the ribo and the RNA abundance, so it cancels and does not change TE unless a pseudocount > 0 is added), log2 (False; report log2(TE)), pseudocount (0.0), orientation ('genes_by_cells' default = rows are genes, or 'cells_by_genes' when rows are samples), genes and cells (name lists for unnamed forms; 'cells' names the samples), max_matrix_elements (default 5000000).
    TE = ribo abundance / RNA abundance with each library scaled to a million; null when the RNA abundance is zero. Also returns the mean TE across samples.
    Because gene length cancels in the ratio, gene_lengths has no effect on the TE values of a gene (with the default pseudocount 0); it only changes the unit named in the result's definition.
    """
    return _run('transcriptomics', 'calculate_translation_efficiency', input, options, _tx.op_translation_efficiency)

def compare(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compare expression profiles by Pearson or Spearman correlation across shared genes.

    input: a gene-by-sample count matrix: {gene: [counts per sample]}, {gene: {sample: count}}, a list of rows, CSV/TSV text or a file path (header = sample names, first column = gene names), a numpy array / DataFrame, or a 10x directory.
    options: other (optional second matrix; without it all sample pairs of the input are compared), method ('pearson' default or 'spearman'), log2 (True: correlate log2(value + 1)),
    normalization ('none' default, or 'cpm'), orientation ('genes_by_cells' default = rows are genes, or 'cells_by_genes' when rows are samples), genes and cells (name lists for unnamed forms; 'cells' names the samples), max_matrix_elements (default 5000000).
    Returns a correlation table {input sample: {other sample: r}} over the genes present in both matrices (null when a profile is constant).
    """
    return _run('transcriptomics', 'compare', input, options, _tx.op_compare)

def deconvolute(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Estimate cell-type proportions of bulk RNA-seq samples by non-negative least squares on a signature matrix.

    input: a gene-by-sample count matrix: {gene: [counts per sample]}, {gene: {sample: count}}, a list of rows, CSV/TSV text or a file path (header = sample names, first column = gene names), a numpy array / DataFrame, or a 10x directory. Requires numpy and scipy.
    options: signature (genes x cell types matrix in the same input forms) OR reference (count matrix of annotated cells) + reference_labels ({cell: type} or a list);
    the signature is then the mean library-size-normalised profile per type. normalize (True: scale every column to sum 1), orientation ('genes_by_cells' default = rows are genes, or 'cells_by_genes' when rows are samples), genes and cells (name lists for unnamed forms; 'cells' names the samples), max_matrix_elements (default 5000000).
    Fractions = NNLS coefficients / their sum over the genes shared with the signature; rmse is the fit residual per sample.
    """
    return _run('transcriptomics', 'deconvolute', input, options, _sc.op_deconvolve)

def detect(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Detect expressed genes in RNA-seq samples with a read-count (or CPM) threshold.

    input: a gene-by-sample count matrix: {gene: [counts per sample]}, {gene: {sample: count}}, a list of rows, CSV/TSV text or a file path (header = sample names, first column = gene names), a numpy array / DataFrame, or a 10x directory.
    options: min_count (1; a gene is detected in a sample when count >= min_count), min_cpm (None; if set, CPM >= min_cpm is used instead), min_samples (1; samples that must detect the gene), orientation ('genes_by_cells' default = rows are genes, or 'cells_by_genes' when rows are samples), genes and cells (name lists for unnamed forms; 'cells' names the samples), max_matrix_elements (default 5000000).
    """
    return _run('transcriptomics', 'detect', input, options, _tx.op_detect)

def normalize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Normalise an RNA-seq count matrix to CPM, TPM, FPKM, median-of-ratios or upper-quartile scale.

    input: a gene-by-sample count matrix: {gene: [counts per sample]}, {gene: {sample: count}}, a list of rows, CSV/TSV text or a file path (header = sample names, first column = gene names), a numpy array / DataFrame, or a 10x directory.
    options: method ('cpm' default; 'tpm' and 'fpkm'/'rpkm' need gene_lengths; 'median_of_ratios' = DESeq2 size factors; 'upper_quartile'), gene_lengths ({gene: bp}),
    log2 (False; then log2(value + pseudocount)), pseudocount (1.0), orientation ('genes_by_cells' default = rows are genes, or 'cells_by_genes' when rows are samples), genes and cells (name lists for unnamed forms; 'cells' names the samples), max_matrix_elements (default 5000000), return_matrix (default True), output_path (optional TSV file to write; never overwrites).
    CPM = count / library size * 1e6; FPKM = count / (length kb * library size in millions); TPM = (count / length kb) / sum(count / length kb) * 1e6;
    median-of-ratios: size factor = median over genes with all counts > 0 of count / geometric mean across samples; upper quartile: size factor = 75th percentile of the sample's non-zero counts / mean of those percentiles.
    """
    return _run('transcriptomics', 'normalize', input, options, _tx.op_normalize)

def plot(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Draw an RNA-seq expression heatmap or a sample-versus-sample scatter plot as SVG.

    input: a gene-by-sample count matrix: {gene: [counts per sample]}, {gene: {sample: count}}, a list of rows, CSV/TSV text or a file path (header = sample names, first column = gene names), a numpy array / DataFrame, or a 10x directory.
    options: kind ('heatmap' default, or 'scatter'), x and y (sample names, required for 'scatter'; points are genes), log2 (True: plot log2(value + 1)), max_genes (50; the most variable genes are shown in the heatmap),
    output_path (optional .svg file to write; never overwrites), title, orientation ('genes_by_cells' default = rows are genes, or 'cells_by_genes' when rows are samples), genes and cells (name lists for unnamed forms; 'cells' names the samples), max_matrix_elements (default 5000000). Returns the SVG text under result['svg'].
    """
    return _run('transcriptomics', 'plot', input, options, _tx.op_plot)

def predict_protein_function(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Predict protein function (not implemented).

    Protein function prediction needs a trained model or functional annotation database (for example InterPro/DeepGOPlus weights) that is not bundled.
    Returns status "not_implemented" for every input.
    """
    return candidate_not_implemented('transcriptomics', 'predict_protein_function', input, options=options, requires="trained protein-function model or annotation database", message="Protein function prediction needs a trained model or functional annotation database (for example InterPro/DeepGOPlus weights) that is not bundled.")

def predict_protein_localization(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Predict protein localization (not implemented).

    Protein localisation prediction needs a trained predictor (for example DeepLoc or SignalP weights) that is not bundled.
    Returns status "not_implemented" for every input.
    """
    return candidate_not_implemented('transcriptomics', 'predict_protein_localization', input, options=options, requires="trained subcellular-localisation model", message="Protein localisation prediction needs a trained predictor (for example DeepLoc or SignalP weights) that is not bundled.")

def qc(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compute per-sample RNA-seq quality metrics: library size, detected genes, top-gene fractions and mitochondrial percentage.

    input: a gene-by-sample count matrix: {gene: [counts per sample]}, {gene: {sample: count}}, a list of rows, CSV/TSV text or a file path (header = sample names, first column = gene names), a numpy array / DataFrame, or a 10x directory.
    options: min_count (1; count threshold for 'detected'), min_library_size, min_genes, max_pct_mito (all default None = report only), mito_prefix ('MT-' and 'mt-'), orientation ('genes_by_cells' default = rows are genes, or 'cells_by_genes' when rows are samples), genes and cells (name lists for unnamed forms; 'cells' names the samples), max_matrix_elements (default 5000000).
    Reports the share of reads in the top gene and top 10 genes and lists samples that fail the given thresholds.
    """
    return _run('transcriptomics', 'qc', input, options, _tx.op_qc)

def quantify(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Quantify transcript or gene expression by running salmon, kallisto (FASTQ input) or featureCounts (BAM/SAM input).

    input: FASTQ path(s) (one single-end, two paired-end) for salmon/kallisto, or BAM/SAM path(s) for featureCounts.
    options: method ('auto' = first of salmon, kallisto on PATH for FASTQ, featureCounts for alignments; or name one), plus the options of run_salmon_quantification,
    run_kallisto_quantification or run_featurecounts. Returns status "not_implemented" when no suitable program is on PATH.
    """
    return _run('transcriptomics', 'quantify', input, options, _tx.op_quantify)

def read(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Read an expression count matrix (mapping, rows, CSV/TSV, array or 10x directory) and summarise genes, samples and sparsity.

    input: a gene-by-sample count matrix: {gene: [counts per sample]}, {gene: {sample: count}}, a list of rows, CSV/TSV text or a file path (header = sample names, first column = gene names), a numpy array / DataFrame, or a 10x directory.
    options: max_return_elements (10000; the matrix values are returned when it has at most this many elements), orientation ('genes_by_cells' default = rows are genes, or 'cells_by_genes' when rows are samples), genes and cells (name lists for unnamed forms; 'cells' names the samples), max_matrix_elements (default 5000000).
    """
    return _run('transcriptomics', 'read', input, options, _sc.op_read)

def run_featurecounts(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Count reads per gene by running Subread featureCounts on alignment files (an external program).

    input: a BAM/SAM path or a list of paths.
    options: annotation (required; GTF/GFF path), output (required; counts file to create, must not exist), strand (0 unstranded, 1 stranded, 2 reversely stranded), feature_type ('exon'),
    attribute ('gene_id'), paired (False; count read pairs), threads (1), timeout_seconds (3600).
    Returns status "not_implemented" when featureCounts is not on PATH; otherwise the per-sample assigned counts parsed from the output table.
    """
    return _run('transcriptomics', 'run_featurecounts', input, options, _tx.op_run_featurecounts)

def run_kallisto_quantification(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Quantify transcript abundance by running kallisto quant on FASTQ reads (an external program).

    input: a FASTQ path (single-end) or a list of two FASTQ paths (paired-end).
    options: index (required; kallisto index file), output_dir (required; must not exist or be empty), fragment_length and fragment_sd (required for single-end), threads (1), timeout_seconds (3600).
    Returns status "not_implemented" when kallisto is not on PATH; otherwise a summary of abundance.tsv.
    """
    return _run('transcriptomics', 'run_kallisto_quantification', input, options, _tx.op_run_kallisto)

def run_salmon_quantification(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Quantify transcript abundance by running salmon quant on FASTQ reads (an external program).

    input: a FASTQ path (single-end) or a list of two FASTQ paths (paired-end).
    options: index (required; salmon index directory), output_dir (required; must not exist or be empty), library_type ('A' = automatic), threads (1), timeout_seconds (3600).
    Returns status "not_implemented" when salmon is not on PATH; otherwise a summary of quant.sf.
    """
    return _run('transcriptomics', 'run_salmon_quantification', input, options, _tx.op_run_salmon)

__all__ = ['align', 'annotate', 'calculate_differential_expression', 'calculate_differential_splicing', 'calculate_translation_efficiency', 'compare', 'deconvolute', 'detect', 'normalize', 'plot', 'predict_protein_function', 'predict_protein_localization', 'qc', 'quantify', 'read', 'run_featurecounts', 'run_kallisto_quantification', 'run_salmon_quantification']

__candidate_adapter_operations__ = set(['align', 'annotate', 'calculate_differential_expression', 'calculate_differential_splicing', 'calculate_translation_efficiency', 'compare', 'deconvolute', 'detect', 'normalize', 'plot', 'predict_protein_function', 'predict_protein_localization', 'qc', 'quantify', 'read', 'run_featurecounts', 'run_kallisto_quantification', 'run_salmon_quantification'])


from dnakit._verbosity import register_public as _register_public_verbose
# Optional-backend exports (see docs/api/features) 
from .backends import (
    aggregate_pseudobulk,
    calculate_cell_qc_metrics,
    calculate_gene_qc_metrics,
    filter_genes_by_min_cells,
    rank_marker_genes,
    read_anndata_summary,
    select_highly_variable_genes,
    summarize_mudata,
    train_scvi_model,
)

__all__ = [*__all__, "aggregate_pseudobulk", "calculate_cell_qc_metrics", "calculate_gene_qc_metrics", "filter_genes_by_min_cells", "rank_marker_genes", "read_anndata_summary", "select_highly_variable_genes", "summarize_mudata", "train_scvi_model"]

_register_public_verbose(globals())
