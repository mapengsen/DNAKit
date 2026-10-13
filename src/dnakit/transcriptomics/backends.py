"""Optional single-cell and gene-level backends: AnnData, scanpy, muon and scvi-tools.

The functions in this module are explicit adapters over optional packages:

- ``read_anndata_summary()``, ``filter_genes_by_min_cells()`` and
  ``aggregate_pseudobulk()`` read and write ``.h5ad`` files with ``anndata``.
  The pseudobulk sums come from ``scanpy.get.aggregate``; cell counts per group
  are counted with pandas.
- ``calculate_gene_qc_metrics()``, ``calculate_cell_qc_metrics()``,
  ``select_highly_variable_genes()`` and ``rank_marker_genes()`` call ``scanpy``.
- ``summarize_mudata()`` reads ``.h5mu`` multi-modal files with ``muon``.
- ``train_scvi_model()`` sets up and trains a scVI model with ``scvi-tools``.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[...]"`` hint. Every
input is checked against explicit ``max_file_bytes`` and ``max_matrix_elements``
limits before the matrix is loaded.
"""

from __future__ import annotations

import importlib.metadata
import math
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_SCANPY_EXTRA = "bio-scanpy"
_SCVI_EXTRA = "bio-scvi"
_MUON_EXTRA = "bio-muon"
_READ_ERRORS = (OSError, KeyError, RuntimeError, TypeError, ValueError)
_DEFAULT_MAX_FILE_BYTES = 1024**3
_MAX_FILE_BYTES_LIMIT = 16 * 1024**3
_DEFAULT_MAX_ELEMENTS = 50_000_000
_MAX_ELEMENTS_LIMIT = 500_000_000
_HVG_FLAVORS = ("seurat", "cell_ranger")
_MARKER_METHODS = ("wilcoxon", "t-test")
_ACCELERATORS = ("cpu", "gpu", "auto")
_PREVIEW_GENES = 5


def _bounded_int(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or "\x00" in value:
        raise ConfigurationError(f"{name} must be a non-empty string without NUL characters.")
    return value


def _optional_text(value: object, name: str) -> str | None:
    return None if value is None else _text(value, name)


def _choose(value: object, allowed: tuple[str, ...], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {', '.join(allowed)}; got {value!r}.")
    return value


def _positive_float(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a positive finite number.")
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise ConfigurationError(f"{name} must be a positive finite number.")
    return number


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _file_limit(max_file_bytes: object) -> int:
    return _bounded_int(max_file_bytes, "max_file_bytes", 1, _MAX_FILE_BYTES_LIMIT)


def _limits(max_file_bytes: object, max_matrix_elements: object) -> tuple[int, int]:
    return (
        _bounded_int(max_file_bytes, "max_file_bytes", 1, _MAX_FILE_BYTES_LIMIT),
        _bounded_int(max_matrix_elements, "max_matrix_elements", 1, _MAX_ELEMENTS_LIMIT),
    )


def _input_file(source: object, suffix: str, max_file_bytes: int) -> Path:
    if isinstance(source, Path):
        source = str(source)
    _text(source, "source")
    path = input_path(source)
    if path.suffix.lower() != suffix:
        raise ConfigurationError(f"source must be a {suffix} file; got {path.name}.")
    size = path.stat().st_size
    if size > max_file_bytes:
        raise ConfigurationError(
            f"{path.name} has {size} bytes, above max_file_bytes={max_file_bytes}."
        )
    return path


def _new_output(value: object, suffix: str) -> Path | None:
    """Validate an output path that must not exist yet, so results are never overwritten."""

    if value is None:
        return None
    if isinstance(value, Path):
        value = str(value)
    path = Path(_text(value, "output_path")).expanduser().resolve()
    if path.suffix.lower() != suffix:
        raise ConfigurationError(f"output_path must end with {suffix}.")
    if path.exists():
        raise ConfigurationError(f"Output file already exists: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {path.parent}")
    return path


def _json(value: object) -> Any:
    """Convert numpy and pandas scalars to JSON-compatible Python values."""

    if hasattr(value, "item") and not isinstance(value, (str, bytes)):
        value = value.item()
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if value is None or isinstance(value, (bool, int, str)):
        return value
    return str(value)


def _provenance(distribution: str, convention: str) -> dict[str, str]:
    try:
        version = importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        version = "version unknown"
    return {"package": f"{distribution} {version}", "convention": convention}


def _load_anndata(path: Path, extra: str, max_matrix_elements: int) -> Any:
    """Check dimensions in backed mode before loading, so oversize matrices are never read."""

    anndata = require("anndata", extra)
    try:
        probe = anndata.read_h5ad(path, backed="r")
    except _READ_ERRORS as exc:
        raise ConfigurationError(f"anndata could not read {path.name}: {exc}") from exc
    try:
        n_obs, n_vars = int(probe.n_obs), int(probe.n_vars)
    finally:
        probe.file.close()
    if n_obs == 0 or n_vars == 0:
        raise ConfigurationError(f"{path.name} must contain at least one cell and one gene.")
    if n_obs * n_vars > max_matrix_elements:
        raise ConfigurationError(
            f"{path.name} has {n_obs * n_vars} matrix elements; "
            f"max_matrix_elements={max_matrix_elements}."
        )
    try:
        return anndata.read_h5ad(path)
    except _READ_ERRORS as exc:
        raise ConfigurationError(f"anndata could not load {path.name}: {exc}") from exc


def _mapping_keys(mapping: Any) -> list[str]:
    # anndata 0.13 reports an empty layers or obsm mapping as a single None key.
    return sorted(str(key) for key in mapping if key is not None)


def _log_normalize(scanpy: Any, adata: Any, target_sum: float) -> None:
    scanpy.pp.normalize_total(adata, target_sum=target_sum)
    scanpy.pp.log1p(adata)


def _require_counts(adata: Any, layer: str | None) -> None:
    np = require("numpy", _SCANPY_EXTRA)
    sparse = require("scipy.sparse", _SCANPY_EXTRA)
    matrix = adata.layers[layer] if layer is not None else adata.X
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if (
        values.size == 0
        or not np.isfinite(values).all()
        or (values < 0).any()
        or (values != np.floor(values)).any()
    ):
        raise ConfigurationError(
            "train_scvi_model needs finite, non-negative integer counts; "
            "pass layer= to select a counts layer."
        )


def read_anndata_summary(
    source: str | Path,
    *,
    max_file_bytes: int = _DEFAULT_MAX_FILE_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarize an .h5ad file and its gene metadata without loading the count matrix."""

    anndata = require("anndata", _SCANPY_EXTRA)
    path = _input_file(source, ".h5ad", _file_limit(max_file_bytes))
    try:
        adata = anndata.read_h5ad(path, backed="r")
    except _READ_ERRORS as exc:
        raise ConfigurationError(f"anndata could not read {path.name}: {exc}") from exc
    try:
        preview = [
            {"gene": str(name), **{str(key): _json(value) for key, value in row.items()}}
            for name, row in zip(
                adata.var_names[:_PREVIEW_GENES],
                adata.var.iloc[:_PREVIEW_GENES].to_dict(orient="records"),
                strict=True,
            )
        ]
        return {
            "status": "completed",
            "operation": "read_anndata_summary",
            "source": str(path),
            "n_obs": int(adata.n_obs),
            "n_vars": int(adata.n_vars),
            "x_dtype": None if adata.X is None else str(adata.X.dtype),
            "obs_columns": [str(column) for column in adata.obs.columns],
            "var_columns": [str(column) for column in adata.var.columns],
            "var_preview": preview,
            "layers": _mapping_keys(adata.layers),
            "obsm": _mapping_keys(adata.obsm),
            "varm": _mapping_keys(adata.varm),
            "uns": _mapping_keys(adata.uns),
            "provenance": _provenance(
                "anndata", "X is not loaded; metadata is read in backed mode"
            ),
        }
    finally:
        adata.file.close()


def calculate_gene_qc_metrics(
    source: str | Path,
    *,
    top_n: int = 10,
    max_file_bytes: int = _DEFAULT_MAX_FILE_BYTES,
    max_matrix_elements: int = _DEFAULT_MAX_ELEMENTS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Per-gene QC metrics from scanpy: cells detected, total counts and dropout."""

    scanpy = require("scanpy", _SCANPY_EXTRA)
    np = require("numpy", _SCANPY_EXTRA)
    top = _bounded_int(top_n, "top_n", 1, 1000)
    file_limit, element_limit = _limits(max_file_bytes, max_matrix_elements)
    path = _input_file(source, ".h5ad", file_limit)
    adata = _load_anndata(path, _SCANPY_EXTRA, element_limit)
    _cells, genes = scanpy.pp.calculate_qc_metrics(
        adata, percent_top=None, log1p=False, inplace=False
    )
    names = [str(name) for name in genes.index]
    totals = genes["total_counts"].to_numpy(dtype=float)
    cells = genes["n_cells_by_counts"].to_numpy(dtype=float)
    means = genes["mean_counts"].to_numpy(dtype=float)
    dropout = genes["pct_dropout_by_counts"].to_numpy(dtype=float)
    order = sorted(range(len(names)), key=lambda index: (-totals[index], names[index]))
    return {
        "status": "completed",
        "operation": "calculate_gene_qc_metrics",
        "source": str(path),
        "n_genes": len(names),
        "genes_detected": int(np.count_nonzero(cells > 0)),
        "median_cells_per_gene": _json(np.median(cells)),
        "top_genes": [
            {
                "gene": names[index],
                "total_counts": _json(totals[index]),
                "n_cells_by_counts": int(cells[index]),
                "mean_counts": _json(means[index]),
                "pct_dropout_by_counts": _json(dropout[index]),
            }
            for index in order[:top]
        ],
        "provenance": _provenance("scanpy", "gene metrics from scanpy.pp.calculate_qc_metrics"),
    }


def calculate_cell_qc_metrics(
    source: str | Path,
    *,
    mito_prefix: str = "MT-",
    mito_threshold_percent: float = 20.0,
    max_file_bytes: int = _DEFAULT_MAX_FILE_BYTES,
    max_matrix_elements: int = _DEFAULT_MAX_ELEMENTS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Per-cell QC summary with the mitochondrial count fraction, using scanpy."""

    scanpy = require("scanpy", _SCANPY_EXTRA)
    np = require("numpy", _SCANPY_EXTRA)
    prefix = _text(mito_prefix, "mito_prefix")
    threshold = _positive_float(mito_threshold_percent, "mito_threshold_percent")
    if threshold > 100:
        raise ConfigurationError("mito_threshold_percent must be at most 100.")
    file_limit, element_limit = _limits(max_file_bytes, max_matrix_elements)
    path = _input_file(source, ".h5ad", file_limit)
    adata = _load_anndata(path, _SCANPY_EXTRA, element_limit)
    mito = np.asarray(adata.var_names.str.startswith(prefix), dtype=bool)
    adata.var["mt"] = mito
    cells, _genes = scanpy.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=False
    )
    genes_per_cell = cells["n_genes_by_counts"].to_numpy(dtype=float)
    totals = cells["total_counts"].to_numpy(dtype=float)
    pct_mito = cells["pct_counts_mt"].to_numpy(dtype=float)
    return {
        "status": "completed",
        "operation": "calculate_cell_qc_metrics",
        "source": str(path),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "mito_prefix": prefix,
        "mito_genes": int(mito.sum()),
        "median_genes_per_cell": _json(np.median(genes_per_cell)),
        "median_total_counts": _json(np.median(totals)),
        "median_pct_counts_mt": _json(np.median(pct_mito)),
        "max_pct_counts_mt": _json(np.max(pct_mito)),
        "mito_threshold_percent": threshold,
        "cells_above_mito_threshold": int(np.count_nonzero(pct_mito > threshold)),
        "provenance": _provenance(
            "scanpy", "percentages are of total counts per cell; mito genes match var_names prefix"
        ),
    }


def filter_genes_by_min_cells(
    source: str | Path,
    *,
    min_cells: int = 3,
    output_path: str | Path | None = None,
    max_file_bytes: int = _DEFAULT_MAX_FILE_BYTES,
    max_matrix_elements: int = _DEFAULT_MAX_ELEMENTS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep genes detected (count > 0) in at least ``min_cells`` cells, using scanpy."""

    scanpy = require("scanpy", _SCANPY_EXTRA)
    minimum = _bounded_int(min_cells, "min_cells", 1, 1_000_000_000)
    file_limit, element_limit = _limits(max_file_bytes, max_matrix_elements)
    path = _input_file(source, ".h5ad", file_limit)
    output = _new_output(output_path, ".h5ad")
    adata = _load_anndata(path, _SCANPY_EXTRA, element_limit)
    mask, _cells_per_gene = scanpy.pp.filter_genes(adata, min_cells=minimum, inplace=False)
    kept = adata[:, mask].copy()
    if output is not None:
        kept.write_h5ad(output)
    return {
        "status": "completed",
        "operation": "filter_genes_by_min_cells",
        "source": str(path),
        "min_cells": minimum,
        "n_cells": int(adata.n_obs),
        "n_genes_before": int(adata.n_vars),
        "n_genes_after": int(kept.n_vars),
        "n_genes_removed": int(adata.n_vars - kept.n_vars),
        "output_path": None if output is None else str(output),
        "provenance": _provenance(
            "scanpy", "a gene is detected in a cell when its count is above 0"
        ),
    }


def select_highly_variable_genes(
    source: str | Path,
    *,
    n_top_genes: int = 2000,
    flavor: str = "seurat",
    log_normalize: bool = True,
    target_sum: float = 10_000.0,
    max_file_bytes: int = _DEFAULT_MAX_FILE_BYTES,
    max_matrix_elements: int = _DEFAULT_MAX_ELEMENTS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Select highly variable genes with scanpy, after optional normalization and log1p."""

    scanpy = require("scanpy", _SCANPY_EXTRA)
    top = _bounded_int(n_top_genes, "n_top_genes", 1, 1_000_000)
    method = _choose(flavor, _HVG_FLAVORS, "flavor")
    normalize = _flag(log_normalize, "log_normalize")
    target = _positive_float(target_sum, "target_sum")
    file_limit, element_limit = _limits(max_file_bytes, max_matrix_elements)
    path = _input_file(source, ".h5ad", file_limit)
    adata = _load_anndata(path, _SCANPY_EXTRA, element_limit)
    if top > adata.n_vars:
        raise ConfigurationError(
            f"n_top_genes={top} exceeds the {adata.n_vars} genes in the input."
        )
    if normalize:
        _log_normalize(scanpy, adata, target)
    table = scanpy.pp.highly_variable_genes(adata, n_top_genes=top, flavor=method, inplace=False)
    selected = table[table["highly_variable"].astype(bool)]
    ordered = selected.assign(gene=selected.index.astype(str)).sort_values(
        ["dispersions_norm", "gene"], ascending=[False, True], kind="mergesort"
    )
    names = ordered["gene"].tolist()
    return {
        "status": "completed",
        "operation": "select_highly_variable_genes",
        "source": str(path),
        "flavor": method,
        "n_top_genes": top,
        "n_genes_input": int(adata.n_vars),
        "n_highly_variable": len(names),
        "log_normalized": normalize,
        "target_sum": target if normalize else None,
        "highly_variable_genes": names,
        "provenance": _provenance(
            "scanpy",
            "ranked by normalized dispersion; input normalized with normalize_total"
            if normalize
            else "ranked by normalized dispersion; input used as given",
        ),
    }


def rank_marker_genes(
    source: str | Path,
    *,
    groupby: str,
    method: str = "wilcoxon",
    n_genes: int = 10,
    log_normalize: bool = True,
    target_sum: float = 10_000.0,
    max_file_bytes: int = _DEFAULT_MAX_FILE_BYTES,
    max_matrix_elements: int = _DEFAULT_MAX_ELEMENTS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Rank marker genes for each group against the rest, using scanpy.tl.rank_genes_groups."""

    scanpy = require("scanpy", _SCANPY_EXTRA)
    key = _text(groupby, "groupby")
    test = _choose(method, _MARKER_METHODS, "method")
    top = _bounded_int(n_genes, "n_genes", 1, 1000)
    normalize = _flag(log_normalize, "log_normalize")
    target = _positive_float(target_sum, "target_sum")
    file_limit, element_limit = _limits(max_file_bytes, max_matrix_elements)
    path = _input_file(source, ".h5ad", file_limit)
    adata = _load_anndata(path, _SCANPY_EXTRA, element_limit)
    if key not in adata.obs.columns:
        raise ConfigurationError(f"groupby {key!r} is not an obs column.")
    labels = adata.obs[key]
    if labels.isna().any():
        raise ConfigurationError("groupby must not contain missing values.")
    if labels.nunique() < 2:
        raise ConfigurationError("groupby needs at least two groups.")
    if normalize:
        _log_normalize(scanpy, adata, target)
    scanpy.tl.rank_genes_groups(
        adata,
        groupby=key,
        method=test,
        n_genes=top,
        use_raw=False,
        corr_method="benjamini-hochberg",
        key_added="rank_genes_groups",
    )
    table = scanpy.get.rank_genes_groups_df(adata, group=None)
    groups: dict[str, list[dict[str, Any]]] = {}
    for row in table.to_dict(orient="records"):
        groups.setdefault(str(row["group"]), []).append(
            {
                "gene": str(row["names"]),
                "score": _json(row["scores"]),
                "logfoldchange": _json(row["logfoldchanges"]),
                "pval_adj": _json(row["pvals_adj"]),
            }
        )
    return {
        "status": "completed",
        "operation": "rank_marker_genes",
        "source": str(path),
        "groupby": key,
        "method": test,
        "n_genes": top,
        "log_normalized": normalize,
        "n_groups": len(groups),
        "groups": {name: groups[name] for name in sorted(groups)},
        "provenance": _provenance(
            "scanpy",
            "each group against the rest; logfoldchanges use log2 of expm1 means of X",
        ),
    }


def aggregate_pseudobulk(
    source: str | Path,
    *,
    group_key: str,
    sample_key: str,
    min_cells_per_sample: int = 1,
    output_path: str | Path | None = None,
    max_file_bytes: int = _DEFAULT_MAX_FILE_BYTES,
    max_matrix_elements: int = _DEFAULT_MAX_ELEMENTS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sum counts per (group, sample) with scanpy.get.aggregate; cell counts use pandas."""

    scanpy = require("scanpy", _SCANPY_EXTRA)
    np = require("numpy", _SCANPY_EXTRA)
    group = _text(group_key, "group_key")
    sample = _text(sample_key, "sample_key")
    if group == sample:
        raise ConfigurationError("group_key and sample_key must be different obs columns.")
    minimum = _bounded_int(min_cells_per_sample, "min_cells_per_sample", 1, 1_000_000_000)
    file_limit, element_limit = _limits(max_file_bytes, max_matrix_elements)
    path = _input_file(source, ".h5ad", file_limit)
    output = _new_output(output_path, ".h5ad")
    adata = _load_anndata(path, _SCANPY_EXTRA, element_limit)
    for key in (group, sample):
        if key not in adata.obs.columns:
            raise ConfigurationError(f"{key!r} is not an obs column.")
    if adata.obs[[group, sample]].isna().any().any():
        raise ConfigurationError("pseudobulk keys must not contain missing values.")
    cell_counts = {
        (str(key_group), str(key_sample)): int(count)
        for (key_group, key_sample), count in adata.obs.groupby(
            [group, sample], observed=True, sort=True
        )
        .size()
        .items()
    }
    pseudo = scanpy.get.aggregate(adata, by=[group, sample], func="sum")
    sums = pseudo.layers["sum"]
    totals = np.asarray(sums.sum(axis=1)).ravel()
    rows = []
    for index, (key_group, key_sample) in enumerate(
        pseudo.obs[[group, sample]].itertuples(index=False, name=None)
    ):
        pair = (str(key_group), str(key_sample))
        rows.append((pair, index, cell_counts[pair]))
    rows.sort(key=lambda row: row[0])
    kept = [row for row in rows if row[2] >= minimum]
    if output is not None:
        subset = pseudo[[index for _pair, index, _count in kept], :].copy()
        subset.obs["n_cells"] = [count for _pair, _index, count in kept]
        subset.write_h5ad(output)
    return {
        "status": "completed",
        "operation": "aggregate_pseudobulk",
        "source": str(path),
        "group_key": group,
        "sample_key": sample,
        "min_cells_per_sample": minimum,
        "n_genes": int(adata.n_vars),
        "n_pseudobulk_samples": len(kept),
        "dropped_samples": len(rows) - len(kept),
        "samples": [
            {
                "group": pair[0],
                "sample": pair[1],
                "n_cells": count,
                "total_counts": _json(totals[index]),
            }
            for pair, index, count in kept
        ],
        "output_path": None if output is None else str(output),
        "provenance": _provenance(
            "scanpy", "sum of raw X per group and sample; counts are not normalized"
        ),
    }


def summarize_mudata(
    source: str | Path,
    *,
    max_file_bytes: int = _DEFAULT_MAX_FILE_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarize a multi-modal .h5mu file: modalities, cell counts and cells shared by all."""

    muon = require("muon", _MUON_EXTRA)
    path = _input_file(source, ".h5mu", _file_limit(max_file_bytes))
    try:
        mdata = muon.read_h5mu(path)
    except _READ_ERRORS as exc:
        raise ConfigurationError(f"muon could not read {path.name}: {exc}") from exc
    names = sorted(str(name) for name in mdata.mod)
    if not names:
        raise ConfigurationError(f"{path.name} contains no modalities.")
    modalities: dict[str, dict[str, int]] = {}
    obs_sets = []
    for name in names:
        modality = mdata.mod[name]
        modalities[name] = {"n_obs": int(modality.n_obs), "n_vars": int(modality.n_vars)}
        obs_sets.append({str(item) for item in modality.obs_names})
    shared = set.intersection(*obs_sets)
    return {
        "status": "completed",
        "operation": "summarize_mudata",
        "source": str(path),
        "n_obs": int(mdata.n_obs),
        "n_obs_shared": len(shared),
        "modalities": modalities,
        "obs_columns": [str(column) for column in mdata.obs.columns],
        "provenance": _provenance(
            "muon",
            "n_obs is the union of cells; n_obs_shared counts cells present in every modality",
        ),
    }


def train_scvi_model(
    source: str | Path,
    *,
    layer: str | None = None,
    batch_key: str | None = None,
    n_latent: int = 10,
    n_hidden: int = 128,
    max_epochs: int = 1,
    seed: int = 0,
    accelerator: str = "cpu",
    max_file_bytes: int = _DEFAULT_MAX_FILE_BYTES,
    max_matrix_elements: int = _DEFAULT_MAX_ELEMENTS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Set up and train an scVI model on raw counts, then return the latent-space summary."""

    scvi = require("scvi", _SCVI_EXTRA)
    np = require("numpy", _SCVI_EXTRA)
    counts_layer = _optional_text(layer, "layer")
    batch = _optional_text(batch_key, "batch_key")
    latent_dims = _bounded_int(n_latent, "n_latent", 1, 256)
    hidden = _bounded_int(n_hidden, "n_hidden", 1, 4096)
    epochs = _bounded_int(max_epochs, "max_epochs", 1, 1000)
    random_seed = _bounded_int(seed, "seed", 0, 2**32 - 1)
    device = _choose(accelerator, _ACCELERATORS, "accelerator")
    file_limit, element_limit = _limits(max_file_bytes, max_matrix_elements)
    path = _input_file(source, ".h5ad", file_limit)
    adata = _load_anndata(path, _SCVI_EXTRA, element_limit)
    if counts_layer is not None and counts_layer not in adata.layers:
        raise ConfigurationError(f"layer {counts_layer!r} is not in the input layers.")
    if batch is not None and batch not in adata.obs.columns:
        raise ConfigurationError(f"batch_key {batch!r} is not an obs column.")
    _require_counts(adata, counts_layer)
    # scvi.settings.seed is process-global; it fixes torch, numpy and random for this call.
    scvi.settings.seed = random_seed
    model_class = scvi.model.SCVI
    model_class.setup_anndata(adata, layer=counts_layer, batch_key=batch)
    model = model_class(adata, n_latent=latent_dims, n_hidden=hidden)
    model.train(
        max_epochs=epochs,
        accelerator=device,
        devices=1,
        enable_progress_bar=verbose,
        enable_model_summary=False,
        enable_checkpointing=False,
    )
    latent = np.asarray(model.get_latent_representation(), dtype=float)
    losses = model.history["elbo_train"]
    return {
        "status": "completed",
        "operation": "train_scvi_model",
        "source": str(path),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "layer": counts_layer,
        "batch_key": batch,
        "n_latent": latent_dims,
        "n_hidden": hidden,
        "max_epochs": epochs,
        "epochs_trained": len(losses),
        "seed": random_seed,
        "accelerator": device,
        "latent_shape": [int(size) for size in latent.shape],
        "latent_finite": bool(np.isfinite(latent).all()),
        "latent_dimension_means": [_json(value) for value in latent.mean(axis=0)],
        "final_train_elbo": _json(losses.iloc[-1, 0]),
        "provenance": _provenance(
            "scvi-tools", "SCVI on raw integer counts; latent is the posterior mean"
        ),
    }


__all__ = [
    "aggregate_pseudobulk",
    "calculate_cell_qc_metrics",
    "calculate_gene_qc_metrics",
    "filter_genes_by_min_cells",
    "rank_marker_genes",
    "read_anndata_summary",
    "select_highly_variable_genes",
    "summarize_mudata",
    "train_scvi_model",
]

_register_public_verbose(globals(), tuple(__all__))
