"""snapATAC2 single-cell ATAC operations that DNAKit did not yet expose (integration pack P38).

Fragment import with per-barcode QC, cell filtering, feature selection, doublet scoring, spectral embedding
(single and multi-modal), marker regions, differential accessibility, regulatory-network edge scores (Spearman
correlation or regression), TF-to-gene networks, pseudobulk aggregation and the DNA motif scanner/test.

Install with ``python -m pip install "dnakit[bio-snapatac2]"``. The single-cell functions work on snapATAC2
``.h5ad`` files (as written by ``import_atac_fragments``) and always write results to NEW files: functions that
change an AnnData file first copy it to the output path. Tables are written as tab-separated files. snapATAC2
failures (including Rust panics) are reported as ``BackendExecutionError``. Tested with snapatac2 2.10.0.
"""

from __future__ import annotations

import shutil
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p38_helpers as h
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_EXTRA = "bio-snapatac2"
_MAX_CELLS = 200_000
_MAX_FEATURES = 5_000_000
_MAX_SEQS = 10_000
_MAX_SEQ_LEN = 1_000_000
_MAX_ROWS = 1_000_000
_REGRESSION_METHODS = ("correlation", "elastic_net", "gb_tree")


def _snap() -> Any:
    return require("snapatac2", _EXTRA)


def _run(call: Any) -> Any:
    """Run a snapATAC2 call; pyo3 panics derive from BaseException and are converted too."""

    try:
        return call()
    except (KeyboardInterrupt, SystemExit):
        raise
    except ConfigurationError:
        raise
    except BaseException as exc:
        raise BackendExecutionError(f"snapatac2 failed: {exc}") from exc


def _chrom_sizes(value: Mapping[str, int]) -> dict[str, int]:
    if not isinstance(value, Mapping) or not value:
        raise ConfigurationError("chrom_sizes must be a non-empty mapping of chromosome to length.")
    return {h.text(k, "chromosome"): h.integer(v, "chromosome length", 1, 10**12) for k, v in value.items()}


def _open(snap: Any, path: Path, *, write: bool) -> Any:
    return _run(lambda: snap.read(str(path), backed="r+" if write else "r"))


def _copy_for_edit(source: object, target: object, source_name: str, target_name: str) -> tuple[Path, Path]:
    src = h.existing(source, source_name)
    dst = h.new_path(target, target_name)
    if dst.suffix != ".h5ad":
        raise ConfigurationError(f"{target_name} must end with .h5ad.")
    shutil.copyfile(src, dst)
    return src, dst


def _obs_records(adata: Any, limit: int, columns: Sequence[str] | None = None) -> list[dict[str, Any]]:
    frame = adata.obs[:]
    names = list(adata.obs_names)
    use = [c for c in frame.columns if columns is None or c in columns]
    out = []
    for i in range(min(len(names), limit)):
        row = {"barcode": names[i]}
        for column in use:
            row[str(column)] = h.plain(frame[column][i])
        out.append(row)
    return out


def import_atac_fragments(
    fragment_file: str,
    chrom_sizes: Mapping[str, int],
    output_h5ad: str,
    *,
    min_num_fragments: int = 200,
    is_paired: bool = True,
    sorted_by_barcode: bool = False,
    whitelist: Sequence[str] | None = None,
    mitochondrial_chromosomes: Sequence[str] = ("chrM", "M"),
    tile_bin_size: int | None = None,
    limit: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Import an ATAC fragment file into a snapATAC2 AnnData file with per-barcode QC, using
    ``snapatac2.pp.import_fragments`` (and ``add_tile_matrix`` when ``tile_bin_size`` is given).

    ``fragment_file`` is a (optionally gzip-compressed) tab-separated file with chromosome, start, end,
    barcode and count columns, 0-based start. ``chrom_sizes`` maps chromosome to length. Barcodes with fewer than
    ``min_num_fragments`` unique fragments are dropped; ``whitelist`` restricts barcodes. ``sorted_by_barcode``
    tells snapATAC2 the file is already barcode-sorted (the default sorts internally). The result file has the
    QC columns ``n_fragment``, ``frac_dup`` and ``frac_mito``. With ``tile_bin_size`` (bp) a cell-by-genomic-bin
    count matrix is added (insertion counting, as snapATAC2's default), which the later steps need.
    """

    source = h.existing(fragment_file, "fragment_file")
    target = h.new_path(output_h5ad, "output_h5ad")
    if target.suffix != ".h5ad":
        raise ConfigurationError("output_h5ad must end with .h5ad.")
    sizes = _chrom_sizes(chrom_sizes)
    minimum = h.integer(min_num_fragments, "min_num_fragments", 1, 10**9)
    paired = h.flag(is_paired, "is_paired")
    barcode_sorted = h.flag(sorted_by_barcode, "sorted_by_barcode")
    keep = None if whitelist is None else [h.text(b, "barcode") for b in whitelist]
    mito = [h.text(c, "mitochondrial chromosome") for c in mitochondrial_chromosomes]
    bins = None if tile_bin_size is None else h.integer(tile_bin_size, "tile_bin_size", 1, 100_000_000)
    cap = h.integer(limit, "limit", 0, 100_000)
    snap = _snap()
    import tempfile

    with tempfile.TemporaryDirectory(prefix="dnakit-snapatac2-") as scratch:
        adata = _run(lambda: snap.pp.import_fragments(
            source, sizes, is_paired=paired, file=target, min_num_fragments=minimum,
            sorted_by_barcode=barcode_sorted, whitelist=keep, chrM=mito, tempdir=Path(scratch), n_jobs=1,
        ))
        try:
            if bins is not None:
                _run(lambda: snap.pp.add_tile_matrix(adata, bin_size=bins, n_jobs=1))
            n_obs, n_vars = int(adata.n_obs), int(adata.n_vars)
            records = _obs_records(adata, cap)
        finally:
            adata.close()
    return {
        "status": "completed",
        "operation": "import_atac_fragments",
        "cell_count": n_obs,
        "feature_count": n_vars,
        "output_h5ad": str(target),
        "qc": records,
        "provenance": h.provenance("snapatac2", "snapatac2.pp.import_fragments"),
    }


def filter_atac_cells(
    input_h5ad: str,
    output_h5ad: str,
    *,
    min_counts: int | None = 1000,
    min_tsse: float | None = 5.0,
    max_counts: int | None = None,
    max_tsse: float | None = None,
    gene_gtf: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter ATAC cells by fragment count and TSS enrichment, using ``snapatac2.pp.filter_cells``.

    The input file is copied to ``output_h5ad`` and the copy is filtered. The defaults are snapATAC2's (at least
    1000 fragments and TSS enrichment 5). TSS enrichment is taken from ``obs["tsse"]``; when absent and a
    thresholds on it is set, give ``gene_gtf`` and it is computed with ``snapatac2.metrics.tsse`` first, or set
    ``min_tsse`` and ``max_tsse`` to null to filter by counts only. Returns the cell counts before and after.
    """

    _, dst = _copy_for_edit(input_h5ad, output_h5ad, "input_h5ad", "output_h5ad")
    lo = None if min_counts is None else h.integer(min_counts, "min_counts", 0, 10**9)
    hi = None if max_counts is None else h.integer(max_counts, "max_counts", 0, 10**10)
    lo_t = None if min_tsse is None else h.number(min_tsse, "min_tsse", 0, 1e6)
    hi_t = None if max_tsse is None else h.number(max_tsse, "max_tsse", 0, 1e6)
    gtf = None if gene_gtf is None else h.existing(gene_gtf, "gene_gtf")
    snap = _snap()
    adata = _open(snap, dst, write=True)
    try:
        before = int(adata.n_obs)
        if (lo_t is not None or hi_t is not None) and "tsse" not in adata.obs[:].columns:
            if gtf is None:
                raise ConfigurationError("obs has no tsse column; give gene_gtf to compute it or set min_tsse/max_tsse to null.")
            _run(lambda: snap.metrics.tsse(adata, gtf, n_jobs=1))
        _run(lambda: snap.pp.filter_cells(adata, min_counts=lo, min_tsse=lo_t, max_counts=hi, max_tsse=hi_t, n_jobs=1))
        after = int(adata.n_obs)
    finally:
        adata.close()
    return {
        "status": "completed",
        "operation": "filter_atac_cells",
        "cells_before": before,
        "cells_after": after,
        "cells_removed": before - after,
        "output_h5ad": str(dst),
        "provenance": h.provenance("snapatac2", "snapatac2.pp.filter_cells"),
    }


def select_atac_features(
    input_h5ad: str,
    output_h5ad: str,
    *,
    n_features: int = 500_000,
    filter_lower_quantile: float = 0.005,
    filter_upper_quantile: float = 0.005,
    whitelist: str | None = None,
    blacklist: str | None = None,
    max_iter: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Select informative genomic bins or peaks, with blacklist and whitelist regions, using
    ``snapatac2.pp.select_features``.

    The input needs a count matrix (for example from ``import_atac_fragments`` with ``tile_bin_size``). Features
    are ranked by accessibility after dropping the lower and upper quantiles; ``blacklist`` and ``whitelist`` are
    BED files of regions to exclude or to keep. The mask is stored in ``var["selected"]`` of the copy written
    to ``output_h5ad``. Returns the number and names (first 50) of the selected features.
    """

    _, dst = _copy_for_edit(input_h5ad, output_h5ad, "input_h5ad", "output_h5ad")
    keep = h.integer(n_features, "n_features", 1, _MAX_FEATURES)
    low = h.number(filter_lower_quantile, "filter_lower_quantile", 0.0, 0.5)
    high = h.number(filter_upper_quantile, "filter_upper_quantile", 0.0, 0.5)
    rounds = h.integer(max_iter, "max_iter", 1, 100)
    black = None if blacklist is None else h.existing(blacklist, "blacklist")
    white = None if whitelist is None else h.existing(whitelist, "whitelist")
    snap = _snap()
    adata = _open(snap, dst, write=True)
    try:
        _run(lambda: snap.pp.select_features(
            adata, n_features=keep, filter_lower_quantile=low, filter_upper_quantile=high,
            whitelist=white, blacklist=black, max_iter=rounds, n_jobs=1, verbose=False,
        ))
        mask = list(adata.var["selected"][:])
        names = list(adata.var_names)
        chosen = [names[i] for i, flag in enumerate(mask) if flag]
        total = len(names)
    finally:
        adata.close()
    return {
        "status": "completed",
        "operation": "select_atac_features",
        "feature_count": total,
        "selected_count": len(chosen),
        "selected_features_preview": chosen[:50],
        "output_h5ad": str(dst),
        "provenance": h.provenance("snapatac2", "snapatac2.pp.select_features"),
    }


def score_atac_doublets(
    input_h5ad: str,
    output_h5ad: str,
    *,
    use_selected_features: bool = True,
    n_comps: int = 15,
    sim_doublet_ratio: float = 2.0,
    expected_doublet_rate: float = 0.1,
    random_state: int = 0,
    limit: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Doublet scoring of scATAC cells (Scrublet on accessibility features), using ``snapatac2.pp.scrublet``.

    Adds ``doublet_score`` and ``doublet_probability`` to ``obs`` of the copy written to ``output_h5ad``. With
    ``use_selected_features`` the features marked by ``select_atac_features`` are used (they must exist).
    ``sim_doublet_ratio`` is the number of simulated doublets per cell, ``expected_doublet_rate`` the prior.
    Returns the scores for the first ``limit`` cells and the number of cells with probability above 0.5.
    """

    _, dst = _copy_for_edit(input_h5ad, output_h5ad, "input_h5ad", "output_h5ad")
    comps = h.integer(n_comps, "n_comps", 2, 200)
    ratio = h.number(sim_doublet_ratio, "sim_doublet_ratio", 0.1, 100)
    rate = h.number(expected_doublet_rate, "expected_doublet_rate", 0.001, 0.99)
    seed = h.integer(random_state, "random_state", 0, 2**31 - 1)
    cap = h.integer(limit, "limit", 0, 100_000)
    use_mask = h.flag(use_selected_features, "use_selected_features")
    snap = _snap()
    adata = _open(snap, dst, write=True)
    try:
        if use_mask and "selected" not in adata.var[:].columns:
            raise ConfigurationError("var has no 'selected' mask; run select_atac_features first or set use_selected_features to false.")
        _run(lambda: snap.pp.scrublet(
            adata, features="selected" if use_mask else None, n_comps=comps, sim_doublet_ratio=ratio,
            expected_doublet_rate=rate, random_state=seed, n_jobs=1, verbose=False,
        ))
        records = _obs_records(adata, cap, ("doublet_score", "doublet_probability"))
        probabilities = [float(v) for v in adata.obs["doublet_probability"][:]]
    finally:
        adata.close()
    return {
        "status": "completed",
        "operation": "score_atac_doublets",
        "cell_count": len(probabilities),
        "likely_doublets": int(sum(p > 0.5 for p in probabilities)),
        "scores": records,
        "output_h5ad": str(h.existing(dst, "output_h5ad")),
        "provenance": h.provenance("snapatac2", "snapatac2.pp.scrublet"),
    }


def compute_atac_spectral_embedding(
    input_h5ads: Sequence[str],
    output_tsv: str,
    *,
    n_comps: int = 30,
    use_selected_features: bool = True,
    distance_metric: str = "cosine",
    weights: Sequence[float] | None = None,
    weighted_by_sd: bool = True,
    random_state: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Spectral (Laplacian eigenmap) embedding of an accessibility matrix, or one joint embedding of several
    modalities, using ``snapatac2.tl.spectral`` (one file) or ``snapatac2.tl.multi_spectral`` (several files).

    One file in ``input_h5ads`` gives ``tl.spectral`` with ``distance_metric`` ``cosine`` or ``jaccard``; two or
    more files (same cells in the same order, different feature spaces) give ``tl.multi_spectral``, optionally
    with per-modality ``weights``. ``use_selected_features`` uses the ``var["selected"]`` mask of
    ``select_atac_features``. The embedding (cells by components, with barcodes) is written to ``output_tsv``;
    the eigenvalues are returned.
    """

    if isinstance(input_h5ads, (str, bytes)) or not input_h5ads:
        raise ConfigurationError("input_h5ads must be a non-empty list of .h5ad paths.")
    paths = [h.existing(p, "input_h5ad") for p in input_h5ads]
    target = h.new_path(output_tsv, "output_tsv")
    comps = h.integer(n_comps, "n_comps", 2, 200)
    metric = h.choice(distance_metric, "distance_metric", ("cosine", "jaccard"))
    seed = h.integer(random_state, "random_state", 0, 2**31 - 1)
    use_mask = h.flag(use_selected_features, "use_selected_features")
    sd = h.flag(weighted_by_sd, "weighted_by_sd")
    wts = None if weights is None else [h.number(w, "weight", 0, 1e6) for w in weights]
    if wts is not None and len(wts) != len(paths):
        raise ConfigurationError("weights must have one value per input file.")
    np = require("numpy", _EXTRA)
    snap = _snap()
    adatas = [_open(snap, p, write=False) for p in paths]
    try:
        names = list(adatas[0].obs_names)
        if any(list(a.obs_names) != names for a in adatas):
            raise ConfigurationError("All input files must hold the same cells in the same order.")
        mask = "selected" if use_mask else None
        if len(adatas) == 1:
            result = _run(lambda: snap.tl.spectral(
                adatas[0], n_comps=comps, features=mask, random_state=seed, distance_metric=metric,
                weighted_by_sd=sd, inplace=False, num_threads=1,
            ))
            kind = "spectral"
        else:
            result = _run(lambda: snap.tl.multi_spectral(
                adatas, n_comps=comps, features=mask, weights=wts, random_state=seed, weighted_by_sd=sd,
            ))
            kind = "multi_spectral"
    finally:
        for a in adatas:
            a.close()
    eigenvalues, vectors = np.asarray(result[0]), np.asarray(result[1])
    with target.open("x", encoding="utf-8") as handle:
        handle.write("barcode\t" + "\t".join(f"comp{i + 1}" for i in range(vectors.shape[1])) + "\n")
        for name, row in zip(names, vectors):
            handle.write(name + "\t" + "\t".join(repr(float(v)) for v in row) + "\n")
    return {
        "status": "completed",
        "operation": "compute_atac_spectral_embedding",
        "method": kind,
        "cell_count": len(names),
        "component_count": int(vectors.shape[1]),
        "eigenvalues": [h.plain(float(v)) for v in eigenvalues],
        "output_tsv": str(target),
        "provenance": h.provenance("snapatac2", f"snapatac2.tl.{kind}"),
    }


def find_atac_marker_regions(
    input_h5ad: str,
    groupby: str,
    *,
    pvalue: float = 0.01,
    verbose: bool = True,
) -> dict[str, Any]:
    """Select marker genomic regions per group by z-score enrichment, using ``snapatac2.tl.marker_regions``.

    ``input_h5ad`` holds an aggregated or cell-level matrix with regions as variables and the grouping column
    ``groupby`` in ``obs``. ``pvalue`` is the z-score p-value cut-off. Returns the marker regions per group.
    """

    path = h.existing(input_h5ad, "input_h5ad")
    key = h.text(groupby, "groupby")
    cutoff = h.number(pvalue, "pvalue", 1e-12, 1.0)
    snap = _snap()
    adata = _open(snap, path, write=False)
    try:
        if key not in adata.obs[:].columns:
            raise ConfigurationError(f"obs has no column {key!r}; columns are {list(adata.obs[:].columns)}.")
        markers = _run(lambda: snap.tl.marker_regions(adata, groupby=key, pvalue=cutoff))
    finally:
        adata.close()
    clean = {str(k): [str(v) for v in values] for k, values in markers.items()}
    return {
        "status": "completed",
        "operation": "find_atac_marker_regions",
        "groupby": key,
        "pvalue": cutoff,
        "marker_counts": {k: len(v) for k, v in clean.items()},
        "markers": clean,
        "provenance": h.provenance("snapatac2", "snapatac2.tl.marker_regions"),
    }


def test_atac_differential_accessibility(
    input_h5ad: str,
    *,
    group1: Sequence[str],
    group2: Sequence[str],
    direction: str = "both",
    min_log_fc: float = 0.25,
    min_pct: float = 0.05,
    output_tsv: str | None = None,
    limit: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Differential accessibility between two cell groups over regions, using ``snapatac2.tl.diff_test``
    (logistic-regression likelihood-ratio tests).

    ``group1`` and ``group2`` are lists of cell barcodes. Regions failing ``min_pct`` (detection fraction) or
    ``min_log_fc`` are not reported. Regions detected in every cell, or in none, of the two groups are dropped
    first because the logistic regression is undefined for them (snapATAC2 would stop with an error).
    ``direction`` is ``positive``, ``negative`` or ``both``. The table (feature name, log2 fold change, p-value,
    adjusted p-value) is returned (first ``limit`` rows) and written to ``output_tsv`` when given.
    """

    path = h.existing(input_h5ad, "input_h5ad")
    target = None if output_tsv is None else h.new_path(output_tsv, "output_tsv")
    first = [h.text(b, "barcode") for b in group1]
    second = [h.text(b, "barcode") for b in group2]
    if not first or not second or set(first) & set(second):
        raise ConfigurationError("group1 and group2 must be non-empty and must not share barcodes.")
    way = h.choice(direction, "direction", ("positive", "negative", "both"))
    fc = h.number(min_log_fc, "min_log_fc", 0.0, 100)
    pct = h.number(min_pct, "min_pct", 0.0, 1.0)
    cap = h.integer(limit, "limit", 0, 100_000)
    np = require("numpy", _EXTRA)
    snap = _snap()
    adata = _open(snap, path, write=False)
    try:
        names = list(adata.obs_names)
        index = {n: i for i, n in enumerate(names)}
        missing = [b for b in first + second if b not in index]
        if missing:
            raise ConfigurationError(f"Unknown barcodes: {missing[:5]}")
        rows = sorted(index[b] for b in first + second)
        matrix = adata.X[rows, :]
        detected = np.asarray((matrix > 0).sum(axis=0)).ravel()
        testable = [str(n) for n, c in zip(adata.var_names, detected) if 0 < c < len(rows)]
        if not testable:
            raise ConfigurationError("No region varies between the cells of the two groups.")
        table = _run(lambda: snap.tl.diff_test(
            adata, first, second, features=testable, direction=way, min_log_fc=fc, min_pct=pct,
        ))
    finally:
        adata.close()
    frame = table.to_pandas() if hasattr(table, "to_pandas") else table
    if target is not None:
        frame.to_csv(target, sep="\t", index=False)
    return {
        "status": "completed",
        "operation": "test_atac_differential_accessibility",
        "group_sizes": [len(first), len(second)],
        "tested_regions": len(testable),
        "reported_regions": int(len(frame)),
        "columns": [str(c) for c in frame.columns],
        "results": h.records(frame, cap),
        "output_tsv": None if target is None else str(target),
        "provenance": h.provenance("snapatac2", "snapatac2.tl.diff_test"),
    }


def _edge_rows(network: Any) -> list[dict[str, Any]]:
    rows = []
    for src, dst, data in network.edge_index_map().values():
        rows.append({
            "source": str(network[src].id), "source_type": str(network[src].type),
            "target": str(network[dst].id), "target_type": str(network[dst].type),
            "distance": h.plain(getattr(data, "distance", None)),
            "cor_score": h.plain(getattr(data, "cor_score", None)),
            "regr_score": h.plain(getattr(data, "regr_score", None)),
        })
    rows.sort(key=lambda r: (r["target"], r["source"]))
    return rows


def score_regulatory_network_edges(
    peak_h5ad: str,
    gene_h5ad: str,
    gene_gtf: str,
    output_tsv: str,
    *,
    method: str = "correlation",
    upstream: int = 250_000,
    downstream: int = 250_000,
    coding_gene_only: bool = True,
    alpha: float = 1.0,
    l1_ratio: float = 0.5,
    select: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Score region-gene edges of a regulatory network by correlation or regression, using
    ``snapatac2.tl.init_network_from_annotation`` with ``add_cor_scores`` or ``add_regr_scores``.

    The network links every region of ``peak_h5ad`` (cell-by-region matrix) to genes whose TSS lies within
    ``upstream``/``downstream`` bp, using the GTF ``gene_gtf``. ``gene_h5ad`` is the matching cell-by-gene
    matrix (same cells in the same order). ``method`` is ``correlation`` (Spearman correlation per edge:
    ``add_cor_scores``) or a regression importance score, ``elastic_net`` (with ``alpha`` and ``l1_ratio``) or
    ``gb_tree`` (``add_regr_scores``). ``select`` restricts the target genes. The edge table (source, target,
    distance, scores) is written to ``output_tsv``.
    """

    peak_path = h.existing(peak_h5ad, "peak_h5ad")
    gene_path = h.existing(gene_h5ad, "gene_h5ad")
    gtf = h.existing(gene_gtf, "gene_gtf")
    target = h.new_path(output_tsv, "output_tsv")
    how = h.choice(method, "method", _REGRESSION_METHODS)
    up = h.integer(upstream, "upstream", 0, 1_000_000)
    down = h.integer(downstream, "downstream", 0, 1_000_000)
    coding = h.flag(coding_gene_only, "coding_gene_only")
    a = h.number(alpha, "alpha", 0.0, 1e6)
    l1 = h.number(l1_ratio, "l1_ratio", 0.0, 1.0)
    chosen = None if select is None else [h.text(g, "gene") for g in select]
    pandas = require("pandas", _EXTRA)
    snap = _snap()
    peaks = _open(snap, peak_path, write=False)
    genes = _open(snap, gene_path, write=False)
    try:
        if list(peaks.obs_names) != list(genes.obs_names):
            raise ConfigurationError("peak_h5ad and gene_h5ad must hold the same cells in the same order.")
        regions = [str(v) for v in peaks.var_names]
        network = _run(lambda: snap.tl.init_network_from_annotation(
            regions, str(gtf), upstream=up, downstream=down, id_type="gene_name", coding_gene_only=coding,
        ))
        if how == "correlation":
            _run(lambda: snap.tl.add_cor_scores(network, gene_mat=genes, peak_mat=peaks, select=chosen))
        else:
            _run(lambda: snap.tl.add_regr_scores(
                network, peak_mat=peaks, gene_mat=genes, select=chosen, method=how, alpha=a, l1_ratio=l1,
            ))
        rows = _edge_rows(network)
    finally:
        peaks.close()
        genes.close()
    pandas.DataFrame(rows).to_csv(target, sep="\t", index=False)
    scored = [r for r in rows if r["cor_score"] is not None or r["regr_score"] is not None]
    return {
        "status": "completed",
        "operation": "score_regulatory_network_edges",
        "method": how,
        "edge_count": len(rows),
        "scored_edge_count": len(scored),
        "edges_preview": rows[:10],
        "output_tsv": str(target),
        "provenance": h.provenance("snapatac2", "snapatac2.tl.add_cor_scores / add_regr_scores"),
    }


def link_tf_to_gene_network(
    regions: Sequence[str],
    gene_gtf: str,
    motif_file: str,
    genome_fasta: str,
    *,
    pvalue: float = 1e-5,
    upstream: int = 250_000,
    downstream: int = 250_000,
    coding_gene_only: bool = True,
    limit: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a TF-to-gene network from region-gene and TF-region (motif) links, using
    ``snapatac2.tl.init_network_from_annotation``, ``add_tf_binding`` and ``link_tf_to_gene``.

    ``regions`` are ``chrom:start-end`` strings. A TF (motif) is linked to a gene when the motif is found in a
    region (``pvalue`` cut-off, genome scanned in ``genome_fasta``, motifs from the MEME file) that is within
    ``upstream``/``downstream`` bp of the gene's TSS in ``gene_gtf``. Returns the TF-to-gene edges.
    """

    if isinstance(regions, (str, bytes)) or not regions:
        raise ConfigurationError("regions must be a non-empty list of chrom:start-end strings.")
    region_list = [h.text(r, "region") for r in regions]
    gtf = h.existing(gene_gtf, "gene_gtf")
    motifs_path = h.existing(motif_file, "motif_file")
    fasta = h.existing(genome_fasta, "genome_fasta")
    cutoff = h.number(pvalue, "pvalue", 1e-15, 1.0)
    up = h.integer(upstream, "upstream", 0, 1_000_000)
    down = h.integer(downstream, "downstream", 0, 1_000_000)
    coding = h.flag(coding_gene_only, "coding_gene_only")
    cap = h.integer(limit, "limit", 0, 100_000)
    snap = _snap()
    motifs = _run(lambda: snap.read_motifs(str(motifs_path)))
    if not motifs:
        raise ConfigurationError("The motif file holds no motifs.")
    network = _run(lambda: snap.tl.init_network_from_annotation(
        region_list, str(gtf), upstream=up, downstream=down, id_type="gene_name", coding_gene_only=coding,
    ))
    _run(lambda: snap.tl.add_tf_binding(network, motifs=motifs, genome_fasta=str(fasta), pvalue=cutoff))
    tf_network = _run(lambda: snap.tl.link_tf_to_gene(network))
    rows = [r for r in _edge_rows(tf_network)]
    return {
        "status": "completed",
        "operation": "link_tf_to_gene_network",
        "motif_count": len(motifs),
        "edge_count": len(rows),
        "edges": [{"tf": r["source"], "gene": r["target"]} for r in rows[:cap]],
        "provenance": h.provenance("snapatac2", "snapatac2.tl.link_tf_to_gene"),
    }


def aggregate_atac_profiles(
    input_h5ad: str,
    output_tsv: str,
    *,
    groupby: str | None = None,
    normalize: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pseudobulk accessibility per group (or all cells), optionally normalised, using ``snapatac2.tl.aggregate_X``.

    ``groupby`` is an ``obs`` column (null aggregates all cells into one profile). ``normalize`` is null (summed
    counts), ``RPM`` or ``RPKM`` (needs regions named ``chrom:start-end`` as variables). The group-by-region
    matrix is written to ``output_tsv`` (rows are groups, columns are regions).
    """

    path = h.existing(input_h5ad, "input_h5ad")
    target = h.new_path(output_tsv, "output_tsv")
    key = None if groupby is None else h.text(groupby, "groupby")
    norm = None if normalize is None else h.choice(normalize, "normalize", ("RPM", "RPKM"))
    np = require("numpy", _EXTRA)
    snap = _snap()
    adata = _open(snap, path, write=False)
    try:
        if key is not None and key not in adata.obs[:].columns:
            raise ConfigurationError(f"obs has no column {key!r}.")
        result = _run(lambda: snap.tl.aggregate_X(adata, groupby=key, normalize=norm))
        matrix = result.X[:]
        matrix = matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)
        groups = [str(v) for v in result.obs_names]
        regions = [str(v) for v in result.var_names]
    finally:
        adata.close()
    with target.open("x", encoding="utf-8") as handle:
        handle.write("group\t" + "\t".join(regions) + "\n")
        for name, row in zip(groups, matrix):
            handle.write(name + "\t" + "\t".join(repr(float(v)) for v in row) + "\n")
    return {
        "status": "completed",
        "operation": "aggregate_atac_profiles",
        "groups": groups,
        "region_count": len(regions),
        "normalize": norm,
        "column_sums": [h.plain(float(v)) for v in matrix.sum(axis=1)],
        "output_tsv": str(target),
        "provenance": h.provenance("snapatac2", "snapatac2.tl.aggregate_X"),
    }


def _motif_scanner(snap: Any, path: Path, motif_id: str | None, base_probabilities: Sequence[float] | None) -> tuple[Any, Any]:
    motifs = _run(lambda: snap.read_motifs(str(path)))
    if not motifs:
        raise ConfigurationError("The motif file holds no motifs.")
    if motif_id is None:
        chosen = motifs[0]
    else:
        wanted = h.text(motif_id, "motif_id")
        matches = [m for m in motifs if wanted in (m.id, m.name)]
        if not matches:
            raise ConfigurationError(f"No motif with id or name {wanted!r}; first ids: {[m.id for m in motifs[:5]]}.")
        chosen = matches[0]
    if base_probabilities is None:
        return chosen, chosen.with_nucl_prob()
    probs = [h.number(p, "base probability", 0.0, 1.0) for p in base_probabilities]
    if len(probs) != 4 or abs(sum(probs) - 1.0) > 1e-6:
        raise ConfigurationError("base_probabilities must be four numbers (A, C, G, T) that sum to 1.")
    return chosen, chosen.with_nucl_prob(*probs)


def _dna_list(values: Sequence[str], name: str) -> list[str]:
    if isinstance(values, (str, bytes)) or not values or len(values) > _MAX_SEQS:
        raise ConfigurationError(f"{name} must be a list of 1 to {_MAX_SEQS} DNA strings.")
    out = []
    for item in values:
        text = h.text(item, f"{name} entry").upper()
        if set(text) - set("ACGT") or len(text) > _MAX_SEQ_LEN:
            raise ConfigurationError(f"{name} entries may contain only A, C, G and T.")
        out.append(text)
    return out


def read_meme_motifs(motif_file: str, *, verbose: bool = True) -> dict[str, Any]:
    """Read DNA motifs from a MEME-format file, using ``snapatac2.read_motifs``.

    Returns, for each motif, its ID, name, family and information content (bits).
    """

    path = h.existing(motif_file, "motif_file")
    snap = _snap()
    motifs = _run(lambda: snap.read_motifs(str(path)))
    return {
        "status": "completed",
        "operation": "read_meme_motifs",
        "motif_count": len(motifs),
        "motifs": [
            {"id": str(m.id), "name": str(m.name), "family": None if m.family is None else str(m.family),
             "information_content": h.plain(float(m.info_content))}
            for m in motifs
        ],
        "provenance": h.provenance("snapatac2", "snapatac2.read_motifs"),
    }


def scan_dna_with_motif(
    motif_file: str,
    sequences: Sequence[str],
    *,
    motif_id: str | None = None,
    mode: str = "find",
    pvalue: float = 1e-5,
    reverse_complement: bool = True,
    report_pvalue: bool = False,
    base_probabilities: Sequence[float] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Scan DNA sequences with a MEME motif, using ``snapatac2.PyDNAMotifScanner``.

    ``mode`` is ``find`` (positions and log-likelihood scores of forward-strand hits above the ``pvalue`` cut-off:
    ``find``; ``report_pvalue`` adds log10 p-values) or ``exist`` (boolean presence per sequence: ``exist`` /
    ``exists``; ``reverse_complement`` also tests the opposite strand). ``motif_id`` picks a motif by ID or name
    (the first motif by default). ``base_probabilities`` is the background A, C, G, T (default 0.25 each).
    """

    path = h.existing(motif_file, "motif_file")
    kind = h.choice(mode, "mode", ("find", "exist"))
    cutoff = h.number(pvalue, "pvalue", 1e-15, 1.0)
    rc = h.flag(reverse_complement, "reverse_complement")
    with_p = h.flag(report_pvalue, "report_pvalue")
    seqs = _dna_list(sequences, "sequences")
    snap = _snap()
    motif, scanner = _motif_scanner(snap, path, motif_id, base_probabilities)
    results: list[Any]
    if kind == "find":
        results = [
            [{"position": int(pos), "score": h.plain(float(score)),
              "log10_pvalue": None if lp is None else h.plain(float(lp))}
             for pos, score, lp in _run(lambda s=s: scanner.find(s, pvalue=cutoff, report_pvalue=with_p))]
            for s in seqs
        ]
    else:
        results = [bool(v) for v in _run(lambda: scanner.exists(seqs, pvalue=cutoff, rc=rc))]
    return {
        "status": "completed",
        "operation": "scan_dna_with_motif",
        "motif": {"id": str(motif.id), "name": str(motif.name)},
        "mode": kind,
        "sequence_count": len(seqs),
        "results": results,
        "provenance": h.provenance("snapatac2", "snapatac2.PyDNAMotifScanner"),
    }


def test_dna_motif_enrichment(
    motif_file: str,
    foreground: Sequence[str],
    background: Sequence[str],
    *,
    motif_id: str | None = None,
    pvalue: float = 1e-5,
    base_probabilities: Sequence[float] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Motif enrichment (log2 fold change and p-value) of a sequence set against a background set, using
    ``snapatac2.PyDNAMotifScanner.with_background`` and ``PyDNAMotifTest.test``.

    The motif occurrence rate in ``background`` (sequences containing the motif above the ``pvalue`` cut-off,
    either strand) is compared with ``foreground``. Do not build the background from the foreground sequences.
    """

    path = h.existing(motif_file, "motif_file")
    cutoff = h.number(pvalue, "pvalue", 1e-15, 1.0)
    fore = _dna_list(foreground, "foreground")
    back = _dna_list(background, "background")
    snap = _snap()
    motif, scanner = _motif_scanner(snap, path, motif_id, base_probabilities)
    tester = _run(lambda: scanner.with_background(back, pvalue=cutoff))
    log_fc, p = _run(lambda: tester.test(fore))
    return {
        "status": "completed",
        "operation": "test_dna_motif_enrichment",
        "motif": {"id": str(motif.id), "name": str(motif.name)},
        "foreground_count": len(fore),
        "background_count": len(back),
        "log2_fold_change": h.plain(float(log_fc)),
        "pvalue": h.plain(float(p)),
        "provenance": h.provenance("snapatac2", "snapatac2.PyDNAMotifTest.test"),
    }


__all__ = [
    "aggregate_atac_profiles",
    "compute_atac_spectral_embedding",
    "filter_atac_cells",
    "find_atac_marker_regions",
    "import_atac_fragments",
    "link_tf_to_gene_network",
    "read_meme_motifs",
    "scan_dna_with_motif",
    "score_atac_doublets",
    "score_regulatory_network_edges",
    "select_atac_features",
    "test_atac_differential_accessibility",
    "test_dna_motif_enrichment",
]

_register_public_verbose(globals(), tuple(__all__))
