"""P41 Kipoi attribution, bounded dataloader batches and genomic writers.

Version note: adds six real upstream operations with JSON results and explicit
output paths. Model loading can access the configured source; use source='dir'
for a local model. No pretrained weights are bundled.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError


def _array(value: Any) -> Any:
    np = require("numpy", "bio-kipoi")
    if isinstance(value, Mapping):
        return {str(k): _array(v) for k, v in value.items()}
    return np.asarray(value)


def _plain(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if hasattr(value, "tolist"):
        return value.tolist()
    return value


def _model_inputs(value: Any, dtype: str | None) -> Any:
    if dtype not in {None, "float32", "float64", "int32", "int64"}:
        raise ConfigurationError("input_dtype must be float32, float64, int32, int64 or None.")
    if isinstance(value, Mapping):
        return {str(k): _model_inputs(v, dtype) for k, v in value.items()}
    np = require("numpy", "bio-kipoi")
    return np.asarray(value, dtype=dtype)


def _path(value: str) -> Path:
    path = Path(value).expanduser().resolve()
    if path.exists() or not path.parent.is_dir():
        raise ConfigurationError("output_path must be new and its parent must exist.")
    return path


def _result(kipoi: Any, **data: Any) -> dict[str, Any]:
    return {
        "status": "completed",
        **data,
        "provenance": {"backend": "kipoi", "version": kipoi.__version__},
    }


def kipoi_layer_activations(
    model: str,
    inputs: Any,
    layer: str | int,
    *,
    source: str = "dir",
    pre_nonlinearity: bool = False,
    input_dtype: str | None = "float32",
    verbose: bool = True,
) -> dict[str, Any]:
    """Evaluate an intermediate layer; input_dtype=None preserves supplied dtype."""
    kipoi = require("kipoi", "bio-kipoi")
    obj = kipoi.get_model(model, source=source, with_dataloader=False)
    values = obj.predict_activation_on_batch(
        _model_inputs(inputs, input_dtype), layer, pre_nonlinearity=pre_nonlinearity
    )
    return _result(kipoi, activations=_plain(values))


def kipoi_input_gradients(
    model: str,
    inputs: Any,
    *,
    source: str = "dir",
    filter_idx: int | None = None,
    avg_func: str | None = None,
    layer: str | int | None = None,
    final_layer: bool = True,
    input_dtype: str | None = "float32",
    pre_nonlinearity: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Calculate upstream input_grad attribution; select a filter or averaging rule."""
    if avg_func not in {None, "sum", "max", "min", "absmax"}:
        raise ConfigurationError("Unsupported gradient averaging function.")
    kipoi = require("kipoi", "bio-kipoi")
    obj = kipoi.get_model(model, source=source, with_dataloader=False)
    values = obj.input_grad(
        _model_inputs(inputs, input_dtype),
        filter_idx=filter_idx,
        avg_func=avg_func,
        layer=layer,
        final_layer=final_layer,
        pre_nonlinearity=pre_nonlinearity,
    )
    return _result(kipoi, gradients=_plain(values))


def kipoi_dataloader_batches(
    dataloader: str,
    *,
    source: str = "dir",
    dataloader_args: Mapping[str, Any] | None = None,
    mode: str = "predict",
    batch_size: int = 32,
    max_batches: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Materialise at most max_batches from raw, train or predict iteration.

    Training uses cycle=False so the wrapper never creates an infinite iterator.
    truncated indicates that at least one additional batch exists.
    """
    if mode not in {"raw", "train", "predict"}:
        raise ConfigurationError("mode must be raw, train or predict.")
    for name, val, cap in [("batch_size", batch_size, 10000), ("max_batches", max_batches, 100)]:
        if isinstance(val, bool) or not isinstance(val, int) or not 1 <= val <= cap:
            raise ConfigurationError(f"{name} must be between 1 and {cap}.")
    kipoi = require("kipoi", "bio-kipoi")
    obj = kipoi.get_dataloader_factory(dataloader, source=source)(**dict(dataloader_args or {}))
    method = {
        "raw": obj.batch_iter,
        "train": obj.batch_train_iter,
        "predict": obj.batch_predict_iter,
    }[mode]
    kwargs = {"batch_size": batch_size}
    if mode == "train":
        kwargs["cycle"] = False
    rows = []
    truncated = False
    for batch in method(**kwargs):
        if len(rows) == max_batches:
            truncated = True
            break
        rows.append(_plain(batch))
    return _result(kipoi, batches=rows, batch_count=len(rows), truncated=truncated)


def kipoi_write_bed_predictions(
    ranges: Mapping[str, Any],
    predictions: Any,
    output_path: str,
    *,
    header: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write BED6 plus predictions using Kipoi BedBatchWriter.

    ranges contains parallel chr/start/end/id/strand arrays. Scores retain the
    upstream BED writer default '.' when absent.
    """
    kipoi = require("kipoi", "bio-kipoi")
    writers = require("kipoi.writers", "bio-kipoi")
    specs = require("kipoi.specs", "bio-kipoi")
    path = _path(output_path)
    schema = {
        "ranges": specs.MetadataStruct(
            doc="Genomic intervals", type=specs.MetadataType.GENOMIC_RANGES
        )
    }
    writer = writers.BedBatchWriter(str(path), schema, header=header)
    try:
        writer.batch_write({"metadata": {"ranges": _array(ranges)}, "preds": _array(predictions)})
    finally:
        writer.close()
    return _result(kipoi, outputs=[str(path)])


def kipoi_write_region_track(
    regions: Sequence[Mapping[str, Any]],
    values: Sequence[Any],
    output_path: str,
    *,
    output_format: str = "bedgraph",
    chrom_sizes: Mapping[str, int] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write per-base bedGraph or bigWig with upstream region_write.

    Regions use chr/start/end; each value array spans end-start bases. bigWig
    requires chromosome sizes and regions ordered by chromosome and start.
    """
    if output_format not in {"bedgraph", "bigwig"} or len(regions) != len(values):
        raise ConfigurationError("Choose bedgraph/bigwig and one value array per region.")
    np = require("numpy", "bio-kipoi")
    arrays = []
    for region, data in zip(regions, values, strict=False):
        arr = np.asarray(data)
        if arr.ndim not in ({1, 2} if output_format == "bedgraph" else {1}):
            raise ConfigurationError("Track values have unsupported dimensionality.")
        start, end = region.get("start"), region.get("end")
        if (
            isinstance(start, bool)
            or isinstance(end, bool)
            or not isinstance(start, int)
            or not isinstance(end, int)
            or not 0 <= start < end
            or len(arr) != end - start
            or not isinstance(region.get("chr"), str)
            or not np.isfinite(arr).all()
        ):
            raise ConfigurationError("Invalid interval, value length or nonfinite value.")
        arrays.append(arr)
    if output_format == "bigwig" and not chrom_sizes:
        raise ConfigurationError("bigWig requires chrom_sizes.")
    kipoi = require("kipoi", "bio-kipoi")
    writers = require("kipoi.writers", "bio-kipoi")
    path = _path(output_path)
    writer = (
        writers.BedGraphWriter(str(path))
        if output_format == "bedgraph"
        else writers.BigWigWriter(str(path), chrom_sizes=list(chrom_sizes.items()), is_sorted=True)
    )
    try:
        for region, array in zip(regions, arrays, strict=False):
            writer.region_write(dict(region), array)
    finally:
        writer.close()
    return _result(kipoi, outputs=[str(path)], region_count=len(regions))


__all__ = [
    "kipoi_dataloader_batches",
    "kipoi_input_gradients",
    "kipoi_layer_activations",
    "kipoi_write_bed_predictions",
    "kipoi_write_region_track",
]
_register_public_verbose(globals(), tuple(__all__))
