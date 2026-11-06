"""Kipoi model-zoo operations as typed DNAKit functions (kipoi 0.8.6).

The user installs kipoi with ``python -m pip install kipoi``. Each function calls kipoi's own
loaders and predictors: ``kipoi.get_model``, ``kipoi.get_dataloader_factory``, the dataloader and
component descriptions, the steps of ``kipoi predict``, ``kipoi test`` and ``kipoi get-example``,
and kipoi's prediction comparison.

Local model directories are read with ``source="dir"`` (a path relative to the current directory, or
an absolute path) or with a local source built from a directory. Those paths never touch the
network. The default source ``"kipoi"`` reads kipoi's git-backed model zoo and is used only when the
caller asks for it.

Prediction output is written as TSV. kipoi 0.8.6 cannot write HDF5 files with NumPy 2 (its writer
refers to the removed ``np.string_``), so HDF5 output is not offered; HDF5 files are read where they
are inputs."""

from __future__ import annotations

import math
import os
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require as _require
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError


def _np():  # numpy is imported lazily: it is not a core dnakit dependency
    return _require("numpy", "bio-kipoi")


_HINT = "Install kipoi with: python -m pip install kipoi"
_TSV = ".tsv"
_DEFAULT_SOURCE = "kipoi"
_MODEL_COLUMNS = ("model", "version", "type", "doc", "inputs", "targets", "license")
_DATALOADER_COLUMNS = ("dataloader", "version", "type", "doc", "inputs", "targets", "license")


def _kipoi() -> Any:
    try:
        import kipoi
        import kipoi.cli.main
    except ImportError as exc:
        raise BackendUnavailableError(f"kipoi could not be imported: {exc}", hint=_HINT) from exc
    return kipoi


def _call(label: str, function: Any, *args: Any, **kwargs: Any) -> Any:
    try:
        return function(*args, **kwargs)
    except (ValueError, KeyError, OSError) as exc:
        raise BackendExecutionError(f"kipoi {label} failed: {exc}") from exc


def _name(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty string.")
    return value.strip()


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be true or false.")
    return value


def _count(value: object, name: str, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ConfigurationError(f"{name} must be an integer of at least {minimum}.")
    return value


def _choice(value: object, name: str, options: tuple[str, ...]) -> str:
    if value not in options:
        raise ConfigurationError(f"{name} must be one of: {', '.join(options)}.")
    return str(value)


def _existing_dir(value: object, name: str) -> Path:
    path = Path(_name(value, name)).expanduser().resolve()
    if not path.is_dir():
        raise ConfigurationError(f"{name} must be an existing directory: {path}")
    return path


def _existing_file(value: object, name: str) -> Path:
    path = Path(_name(value, name)).expanduser().resolve()
    if not path.is_file():
        raise ConfigurationError(f"{name} must be an existing file: {path}")
    return path


def _new_tsv(value: object, name: str) -> Path:
    path = Path(_name(value, name)).expanduser().resolve()
    if path.suffix.lower() != _TSV:
        raise ConfigurationError(
            f"{name} must end with {_TSV}: kipoi 0.8.6 cannot write HDF5 with NumPy 2, so only TSV output is offered."
        )
    if path.exists():
        raise ConfigurationError(f"{name} must not exist yet: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {path.parent}")
    return path


def _mapping(value: object, name: str) -> dict[str, Any]:
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"{name} must be a mapping of argument names to values.")
    return {str(key): item for key, item in value.items()}


def _plain(value: Any) -> Any:
    """Convert kipoi description values to JSON-serialisable Python values."""
    from kipoi_utils.external.related.fields import UNSPECIFIED

    if value is None or isinstance(value, (bool, int, str)):
        return value
    if isinstance(value, float):
        return None if math.isnan(value) else value
    if isinstance(value, UNSPECIFIED):
        return None
    if hasattr(value, "url") and hasattr(value, "md5"):
        return {
            "url": _plain(value.url),
            "md5": _plain(value.md5),
            "name": _plain(getattr(value, "name", None)),
        }
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_plain(item) for item in value]
    if hasattr(value, "tolist"):
        return _plain(value.tolist())
    return str(value)


def _schema(value: Any) -> Any:
    """Describe a kipoi input or target schema (a single array, a list or a mapping of arrays)."""
    if value is None:
        return None
    if isinstance(value, Mapping):
        return {str(key): _schema(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_schema(item) for item in value]
    if hasattr(value, "name") and hasattr(value, "shape"):
        return {
            "name": _plain(value.name),
            "shape": _plain(value.shape),
            "doc": _plain(getattr(value, "doc", None)),
        }
    return _plain(value)


def _table_rows(table: Any, columns: Sequence[str]) -> list[dict[str, Any]]:
    present = [column for column in columns if column in table.columns]
    return [
        {key: _plain(item) for key, item in row.items()}
        for row in table[present].to_dict(orient="records")
    ]


def _leading(values: Any) -> int:
    if isinstance(values, Mapping):
        return _leading(next(iter(values.values()))) if values else 0
    if isinstance(values, (list, tuple)):
        return len(values)
    shape = _np().shape(values)
    return int(shape[0]) if shape else 0


def _shape(values: Any) -> Any:
    if isinstance(values, Mapping):
        return {str(key): _shape(item) for key, item in values.items()}
    return [int(size) for size in _np().shape(values)]


def _provenance(kipoi: Any) -> dict[str, Any]:
    return {"backend": "kipoi", "version": _plain(getattr(kipoi, "__version__", None))}


def load_kipoi_model(
    model: str,
    *,
    source: str = _DEFAULT_SOURCE,
    with_dataloader: bool = True,
    verbose: bool = True,
) -> Any:
    """Load a kipoi model with its default dataloader and return the kipoi model object (kipoi.get_model).

    With ``source="dir"`` the model is a directory relative to the current directory (or an absolute path) and
    nothing is downloaded. With the default source kipoi clones or updates its git-backed model zoo first. The
    returned object has ``predict_on_batch``, ``schema``, ``info``, ``default_dataloader`` and the other attributes
    that kipoi attaches.
    """
    kipoi = _kipoi()
    return _call(
        "get_model",
        kipoi.get_model,
        _name(model, "model"),
        _name(source, "source"),
        with_dataloader=_flag(with_dataloader, "with_dataloader"),
    )


def load_kipoi_dataloader(
    dataloader: str, *, source: str = _DEFAULT_SOURCE, verbose: bool = True
) -> Any:
    """Return the dataloader class of a kipoi dataloader (kipoi.get_dataloader_factory).

    Instantiate the class with the dataloader arguments and iterate it with ``batch_iter``. With ``source="dir"``
    the dataloader is a directory relative to the current directory, for example ``"my_model/."``.
    """
    kipoi = _kipoi()
    return _call(
        "get_dataloader_factory",
        kipoi.get_dataloader_factory,
        _name(dataloader, "dataloader"),
        source=_name(source, "source"),
    )


def describe_kipoi_dataloader(
    dataloader: str,
    *,
    source: str = _DEFAULT_SOURCE,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the description of a kipoi dataloader: type, defined class, output schema and arguments.

    Each argument reports its documentation, type, whether it is optional and its example value. Remote example
    files are reported with their URL and md5 hash; nothing is downloaded.
    """
    kipoi = _kipoi()
    name = _name(dataloader, "dataloader")
    source_name = _name(source, "source")
    descr = _call("get_dataloader_descr", kipoi.get_dataloader_descr, name, source=source_name)
    arguments = {
        str(arg_name): {
            "doc": _plain(getattr(arg, "doc", None)),
            "type": _plain(getattr(arg, "type", None)),
            "optional": _plain(getattr(arg, "optional", None)),
            "example": _plain(getattr(arg, "example", None)),
        }
        for arg_name, arg in descr.args.items()
    }
    info = descr.info
    return {
        "status": "completed",
        "operation": "describe_kipoi_dataloader",
        "dataloader": name,
        "source": source_name,
        "description": {
            "type": _plain(descr.type),
            "defined_as": _plain(descr.defined_as),
            "version": _plain(info.version),
            "doc": _plain(info.doc),
            "license": _plain(info.license),
            "inputs": _schema(descr.output_schema.inputs),
            "targets": _schema(descr.output_schema.targets),
            "arguments": arguments,
        },
        "provenance": _provenance(kipoi),
    }


def describe_kipoi_component(
    component_dir: str, *, which: str = "model", verbose: bool = True
) -> dict[str, Any]:
    """Read the ``model.yaml`` (``which="model"``) or ``dataloader.yaml`` (``which="dataloader"``) of a local directory.

    The directory is read as it is; nothing is downloaded. The result has the type, documentation, license and the
    input and target schemas, plus the defined class for a dataloader.
    """
    kipoi = _kipoi()
    directory = _existing_dir(component_dir, "component_dir")
    kind = _choice(which, "which", ("model", "dataloader"))
    descr = _call(
        "load_component_descr", kipoi.sources.load_component_descr, str(directory), which=kind
    )
    info = descr.info
    result: dict[str, Any] = {
        "status": "completed",
        "operation": "describe_kipoi_component",
        "component_dir": str(directory),
        "which": kind,
        "type": _plain(descr.type),
        "version": _plain(info.version),
        "doc": _plain(info.doc),
        "license": _plain(info.license),
    }
    if kind == "model":
        result["inputs"] = _schema(descr.schema.inputs)
        result["targets"] = _schema(descr.schema.targets)
    else:
        result["defined_as"] = _plain(descr.defined_as)
        result["inputs"] = _schema(descr.output_schema.inputs)
        result["targets"] = _schema(descr.output_schema.targets)
    result["provenance"] = _provenance(kipoi)
    return result


def list_kipoi_directory_models(directory: str, *, verbose: bool = True) -> dict[str, Any]:
    """List the models of a local kipoi model zoo, that is a directory whose subfolders hold ``model.yaml`` files.

    The listing reads only the local directory (kipoi.config.list_models on a local source).
    """
    kipoi = _kipoi()
    root = _existing_dir(directory, "directory")
    table = _call(
        "list_models",
        kipoi.config.list_models,
        sources={"directory": kipoi.sources.LocalSource(str(root))},
    )
    models = _table_rows(table, _MODEL_COLUMNS)
    return {
        "status": "completed",
        "operation": "list_kipoi_directory_models",
        "directory": str(root),
        "total_models": len(models),
        "models": models,
        "provenance": _provenance(kipoi),
    }


def list_kipoi_directory_dataloaders(directory: str, *, verbose: bool = True) -> dict[str, Any]:
    """List the dataloaders of a local kipoi model zoo, that is a directory whose subfolders hold ``dataloader.yaml``.

    The listing reads only the local directory (kipoi.config.list_dataloaders on a local source), so no configured
    network source is touched.
    """
    kipoi = _kipoi()
    root = _existing_dir(directory, "directory")
    table = _call(
        "list_dataloaders",
        kipoi.config.list_dataloaders,
        sources={"directory": kipoi.sources.LocalSource(str(root))},
    )
    dataloaders = _table_rows(table, _DATALOADER_COLUMNS)
    return {
        "status": "completed",
        "operation": "list_kipoi_directory_dataloaders",
        "directory": str(root),
        "total_dataloaders": len(dataloaders),
        "dataloaders": dataloaders,
        "provenance": _provenance(kipoi),
    }


def predict_kipoi_model(
    model: str,
    output_file: str,
    *,
    source: str = _DEFAULT_SOURCE,
    dataloader_args: Mapping[str, Any] | None = None,
    batch_size: int = 32,
    num_workers: int = 0,
    keep_inputs: bool = False,
    keep_metadata: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run a kipoi model over the data of its default dataloader and write the predictions to a TSV file.

    This performs the steps of ``kipoi predict``: the model and its default dataloader are loaded, the dataloader
    arguments are validated, the dataloader is iterated in batches, the model predicts each batch and kipoi's TSV
    writer stores the results. ``dataloader_args`` maps argument names to values, for example
    ``{"fasta_file": "reads.fa"}``; relative paths are relative to the current directory. ``output_file`` must not
    exist yet and must end with ``.tsv``.
    """
    kipoi = _kipoi()
    name = _name(model, "model")
    source_name = _name(source, "source")
    out = _new_tsv(output_file, "output_file")
    size = _count(batch_size, "batch_size", 1)
    workers = _count(num_workers, "num_workers", 0)
    keep_in = _flag(keep_inputs, "keep_inputs")
    keep_meta = _flag(keep_metadata, "keep_metadata")
    arguments = _mapping(dataloader_args, "dataloader_args")

    model_obj = _call("get_model", kipoi.get_model, name, source_name)
    dataloader_cls = model_obj.default_dataloader
    if dataloader_cls is None:
        raise ConfigurationError(f"model {name!r} does not declare a default dataloader.")
    dataloader = dataloader_cls(**kipoi.pipeline.validate_kwargs(dataloader_cls, arguments))
    writer = kipoi.writers.get_writer(str(out), dataloader.get_output_schema().metadata)
    batches = 0
    samples = 0
    try:
        for chunk in dataloader.batch_iter(batch_size=size, num_workers=workers):
            predictions = model_obj.predict_on_batch(chunk["inputs"])
            samples += _leading(predictions)
            batches += 1
            writer.batch_write(
                kipoi.cli.main.prepare_batch(
                    chunk, predictions, keep_inputs=keep_in, keep_metadata=keep_meta
                )
            )
    finally:
        writer.close()
    written = out.is_file()
    return {
        "status": "completed" if written else "failed",
        "operation": "predict_kipoi_model",
        "model": name,
        "source": source_name,
        "dataloader": dataloader_cls.__name__,
        "batches": batches,
        "samples": samples,
        "outputs": [str(out)] if written else [],
        "provenance": _provenance(kipoi),
    }


def validate_kipoi_model_example(
    model_dir: str,
    *,
    batch_size: int = 32,
    output_file: str | None = None,
    expect_file: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run a local model directory on its example input and optionally compare with expected predictions (kipoi test).

    The example is the dataloader's example arguments. ``output_file`` (optional, ``.tsv``) receives the inputs and
    predictions as kipoi writes them. ``expect_file`` is an HDF5 file with ``inputs`` and ``preds`` datasets, as
    written by ``kipoi predict -o file.h5 --keep_inputs``; when it is omitted, a local ``test.expect`` path in
    model.yaml (relative to the model directory) is used. A remote ``test.expect`` is reported and not downloaded.
    The comparison uses the precision in model.yaml (``test.precision_decimal``) and kipoi's array comparison.
    """
    kipoi = _kipoi()
    directory = _existing_dir(model_dir, "model_dir")
    size = _count(batch_size, "batch_size", 1)
    out = _new_tsv(output_file, "output_file") if output_file is not None else None
    expect_path = _existing_file(expect_file, "expect_file") if expect_file is not None else None

    model_obj = _call("get_model", kipoi.get_model, str(directory), "dir")
    precision = model_obj.test.precision_decimal
    declared = model_obj.test.expect
    note = None
    if expect_path is None and declared is not None:
        if hasattr(declared, "url"):
            note = "test.expect is a remote file in model.yaml; it was not downloaded. Pass expect_file to compare."
        else:
            candidate = Path(str(declared)).expanduser()
            if not candidate.is_absolute():
                candidate = directory / candidate
            if candidate.is_file():
                expect_path = candidate.resolve()
            else:
                note = f"test.expect file was not found: {candidate}"

    predictions = model_obj.pipeline.predict_example(
        batch_size=size,
        output_file=None if out is None else str(out),
    )

    matched: bool | None = None
    compared = 0
    if expect_path is not None:
        matched = True
        reader = kipoi.readers.HDF5Reader(str(expect_path))
        reader.open()
        try:
            for chunk in reader.batch_iter(batch_size=size):
                if "inputs" not in chunk or "preds" not in chunk:
                    raise ConfigurationError(
                        "expect_file must contain the 'inputs' and 'preds' datasets."
                    )
                if isinstance(chunk["preds"], Mapping):
                    raise ConfigurationError(
                        "expect_file holds several prediction arrays; only a single prediction array is compared."
                    )
                observed = model_obj.predict_on_batch(chunk["inputs"])
                compared += _leading(chunk["preds"])
                try:
                    kipoi.cli.main.compare_numpy_dict(
                        observed, chunk["preds"], exact=False, decimal=precision
                    )
                except AssertionError:
                    matched = False
                    break
        finally:
            reader.close()

    written = out is not None and out.is_file()
    return {
        "status": "failed" if matched is False else "completed",
        "operation": "validate_kipoi_model_example",
        "model_dir": str(directory),
        "samples": _leading(predictions),
        "prediction_shape": _shape(predictions),
        "outputs": [str(out)] if written else [],
        "expect_file": None if expect_path is None else str(expect_path),
        "compared_samples": compared,
        "matched": matched,
        "precision_decimal": _plain(precision),
        "note": note,
        "provenance": _provenance(kipoi),
    }


def get_kipoi_example_kwargs(
    model: str,
    *,
    source: str = _DEFAULT_SOURCE,
    output_dir: str = "example",
    dry_run: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the example keyword arguments of a model's default dataloader (kipoi get-example).

    Remote example files are downloaded into ``output_dir`` (an absolute path is recorded in the result). With
    ``dry_run=True`` only their target paths are returned and nothing is downloaded. Plain example values, such as
    ``"example.fa"``, are returned as they are written in dataloader.yaml, relative to the model directory.
    """
    from kipoi_utils.utils import cd

    kipoi = _kipoi()
    name = _name(model, "model")
    source_name = _name(source, "source")
    out = Path(_name(output_dir, "output_dir")).expanduser().absolute()
    dry = _flag(dry_run, "dry_run")

    md = _call("get_model_descr", kipoi.get_model_descr, name, source=source_name)
    if isinstance(md.default_dataloader, kipoi.specs.DataLoaderImport):
        source_obj = _call("get_source", kipoi.get_source, source_name)
        with cd(source_obj.get_model_dir(name)):
            dl_descr = md.default_dataloader.get()
    else:
        dl_descr = _call(
            "get_dataloader_descr",
            kipoi.get_dataloader_descr,
            os.path.join(name, md.default_dataloader),
            source=source_name,
        )
    kwargs = dl_descr.download_example(output_dir=str(out), absolute_path=True, dry_run=dry)
    return {
        "status": "completed",
        "operation": "get_kipoi_example_kwargs",
        "model": name,
        "source": source_name,
        "output_dir": str(out),
        "dry_run": dry,
        "example_kwargs": _plain(kwargs),
        "provenance": _provenance(kipoi),
    }


def compare_kipoi_predictions(
    observed: Any,
    expected: Any,
    *,
    decimal: int = 7,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare two prediction arrays to ``decimal`` decimal places with kipoi's comparison (kipoi test).

    Two arrays match when ``abs(expected - observed) < 1.5 * 10**(-decimal)`` holds everywhere and the shapes agree.
    Nested dictionaries of arrays are not accepted here: kipoi's dictionary comparison reports equal arrays as
    different.
    """
    kipoi = _kipoi()
    places = _count(decimal, "decimal", 0)
    left = _np().asarray(observed)
    right = _np().asarray(expected)
    if left.dtype == object or right.dtype == object:
        raise ConfigurationError(
            "observed and expected must be numeric arrays, not nested objects."
        )
    matched = True
    mismatch = None
    try:
        kipoi.cli.main.compare_numpy_dict(left, right, exact=False, decimal=places)
    except AssertionError as exc:
        matched = False
        mismatch = str(exc).strip().splitlines()[0] if str(exc).strip() else "arrays differ"
    return {
        "status": "completed",
        "operation": "compare_kipoi_predictions",
        "matched": matched,
        "decimal": places,
        "shape": [int(size) for size in left.shape],
        "expected_shape": [int(size) for size in right.shape],
        "mismatch": mismatch,
        "provenance": _provenance(kipoi),
    }


__all__ = [
    "compare_kipoi_predictions",
    "describe_kipoi_component",
    "describe_kipoi_dataloader",
    "get_kipoi_example_kwargs",
    "list_kipoi_directory_dataloaders",
    "list_kipoi_directory_models",
    "load_kipoi_dataloader",
    "load_kipoi_model",
    "predict_kipoi_model",
    "validate_kipoi_model_example",
]

_register_public_verbose(globals())
