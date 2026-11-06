"""Kipoi operations on dataloader batches, dataloader example inputs and dataloader default files (kipoi 0.8.6).

The user installs kipoi with ``python -m pip install kipoi``. Each function calls kipoi's own code:

- ``prepare_kipoi_batch`` attaches predictions to a dataloader batch as ``kipoi predict`` does
  (``kipoi.cli.main.prepare_batch``).
- ``concatenate_kipoi_batches`` joins dataloader batches field by field (``kipoi.cli.main.numpy_collate_concat``).
- ``kipoi_dataloader_example_kwargs`` resolves the example inputs of a local dataloader directory
  (``kipoi.specs.example_kwargs``, the computation behind ``kipoi get-example``).
- ``kipoi_dataloader_default_args`` resolves the default arguments of a local dataloader directory and fetches
  remote default files (``kipoi.specs.download_default_args``, the step kipoi runs when it loads a model).
- ``test_kipoi_source_models`` runs ``kipoi test`` on the models of a local source directory in the current
  environment. It is the prediction part of ``kipoi test-source``; conda environments and git ranges are not used.

Relative paths follow kipoi: plain example values are returned as written in dataloader.yaml, and default-file
checks are made against the current working directory. Nothing touches the network unless a dataloader names a
remote example or default file and the caller asks for it to be downloaded.
"""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require as _require
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError


def _np():  # numpy is imported lazily: it is not a core dnakit dependency
    return _require("numpy", "bio-kipoi")


_HINT = "Install kipoi with: python -m pip install kipoi"

__all__ = [
    "concatenate_kipoi_batches",
    "kipoi_dataloader_default_args",
    "kipoi_dataloader_example_kwargs",
    "prepare_kipoi_batch",
    "test_kipoi_source_models",
]


def _kipoi() -> Any:
    try:
        import kipoi
        import kipoi.cli.main
        import kipoi.sources
        import kipoi.specs
    except ImportError as exc:
        raise BackendUnavailableError(f"kipoi could not be imported: {exc}", hint=_HINT) from exc
    return kipoi


def _call(label: str, function: Any, *args: Any, **kwargs: Any) -> Any:
    try:
        return function(*args, **kwargs)
    except (ValueError, KeyError, OSError, TypeError, AssertionError) as exc:
        raise BackendExecutionError(f"kipoi {label} failed: {exc}") from exc


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty string.")
    return value.strip()


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _count(value: object, name: str, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ConfigurationError(f"{name} must be an integer of at least {minimum}.")
    return value


def _existing_dir(value: object, name: str) -> Path:
    path = Path(_text(value, name)).expanduser().absolute()
    if not path.is_dir():
        raise ConfigurationError(f"{name} must be an existing directory: {path}")
    return path


def _plain(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, _np().ndarray):
        return value.tolist()
    if isinstance(value, _np().generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    return str(value)


def _provenance(kipoi: Any) -> dict[str, Any]:
    return {"backend": "kipoi", "version": _plain(getattr(kipoi, "__version__", None))}


def prepare_kipoi_batch(
    dataloader_batch: dict[str, Any],
    predictions: Any,
    *,
    keep_inputs: bool = False,
    keep_metadata: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Attach model predictions to one dataloader batch the way ``kipoi predict`` does (kipoi.cli.main.prepare_batch).

    The predictions are stored under ``"preds"``. Unless ``keep_inputs`` is true, the ``"inputs"`` and ``"targets"``
    entries are removed; unless ``keep_metadata`` is true, ``"metadata"`` is removed. The batch is copied first, so the
    caller's dictionary is not changed.
    """
    kipoi = _kipoi()
    if not isinstance(dataloader_batch, Mapping):
        raise ConfigurationError(
            "dataloader_batch must be a mapping, such as one batch yielded by a kipoi dataloader."
        )
    batch = dict(dataloader_batch)
    return _call(
        "prepare_batch",
        kipoi.cli.main.prepare_batch,
        batch,
        predictions,
        keep_inputs=_flag(keep_inputs, "keep_inputs"),
        keep_metadata=_flag(keep_metadata, "keep_metadata"),
    )


def concatenate_kipoi_batches(batches: Sequence[Any], *, verbose: bool = True) -> Any:
    """Join dataloader batches along the first axis, field by field (kipoi.cli.main.numpy_collate_concat).

    Every batch must have the same structure, as the batches yielded by a kipoi dataloader do: arrays are joined with
    ``numpy.concatenate``, and nested dictionaries or lists are joined entry by entry. The result has the structure of
    one batch with all the samples of all batches.
    """
    kipoi = _kipoi()
    if isinstance(batches, (str, bytes, Mapping)):
        raise ConfigurationError("batches must be a sequence of dataloader batches.")
    items = list(batches)
    if not items:
        raise ConfigurationError("batches must contain at least one dataloader batch.")
    return _call("numpy_collate_concat", kipoi.cli.main.numpy_collate_concat, items)


def _load_dataloader(kipoi: Any, root: Path) -> Any:
    return _call(
        "load_component_descr",
        kipoi.sources.load_component_descr,
        str(root),
        which="dataloader",
    )


def kipoi_dataloader_example_kwargs(
    dataloader_dir: str,
    *,
    cache_dir: str | None = None,
    absolute_path: bool = True,
    dry_run: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Resolve the example arguments of a local kipoi dataloader (kipoi.specs.example_kwargs).

    ``dataloader_dir`` is a directory holding ``dataloader.yaml``. Every argument that has an ``example`` is returned:
    plain values as they are written in the YAML, and remote example files as their local path under ``cache_dir``.
    With ``cache_dir``, remote examples that are missing, or whose md5 does not match, are downloaded there, unless
    ``dry_run`` is true, in which case only the target paths are returned. Without ``cache_dir``, a remote example is
    reported as its url and md5 and nothing is downloaded.
    """
    kipoi = _kipoi()
    root = _existing_dir(dataloader_dir, "dataloader_dir")
    cache = (
        None if cache_dir is None else Path(_text(cache_dir, "cache_dir")).expanduser().absolute()
    )
    absolute = _flag(absolute_path, "absolute_path")
    dry = _flag(dry_run, "dry_run")
    descr = _load_dataloader(kipoi, root)
    kwargs = _call(
        "example_kwargs",
        kipoi.specs.example_kwargs,
        descr.args,
        cache_path=None if cache is None else str(cache),
        absolute_path=absolute,
        dry_run=dry,
    )
    examples: dict[str, Any] = {}
    for name, value in kwargs.items():
        if isinstance(value, kipoi.specs.RemoteFile):
            examples[str(name)] = {"url": value.url, "md5": value.md5, "name": value.name}
        else:
            examples[str(name)] = _plain(value)
    return {
        "status": "completed",
        "operation": "kipoi_dataloader_example_kwargs",
        "dataloader_dir": str(root),
        "cache_dir": None if cache is None else str(cache),
        "dry_run": dry,
        "example_kwargs": examples,
        "provenance": _provenance(kipoi),
    }


def kipoi_dataloader_default_args(
    dataloader_dir: str, output_dir: str, *, verbose: bool = True
) -> dict[str, Any]:
    """Resolve the default arguments of a local kipoi dataloader (kipoi.specs.download_default_args).

    Remote default files are downloaded into ``output_dir`` (created when needed) and their local paths replace the
    remote entries. Other defaults that are not ``None`` are returned as they are, with existing file paths made
    absolute and booleans kept as flags. A dataloader without remote defaults downloads nothing.
    """
    kipoi = _kipoi()
    root = _existing_dir(dataloader_dir, "dataloader_dir")
    out = Path(_text(output_dir, "output_dir")).expanduser().absolute()
    descr = _load_dataloader(kipoi, root)
    # kipoi 0.8.6 passes every non-bool default to os.path.exists, which reads an int as a file descriptor and
    # then fails in os.path.abspath when that descriptor happens to be open. Numeric and container defaults
    # are plain values, so they bypass the upstream call and are returned unchanged.
    upstream_args: OrderedDict[str, Any] = OrderedDict()
    plain: dict[str, Any] = {}
    for key, argument in descr.args.items():
        default = argument.default
        if (isinstance(default, (int, float)) and not isinstance(default, bool)) or isinstance(
            default, (list, tuple, dict)
        ):
            plain[key] = default
        else:
            upstream_args[key] = argument
    overrides = dict(
        _call("download_default_args", kipoi.specs.download_default_args, upstream_args, str(out))
    )
    overrides.update(plain)
    return {
        "status": "completed",
        "operation": "kipoi_dataloader_default_args",
        "dataloader_dir": str(root),
        "output_dir": str(out),
        "overrides": _plain(overrides),
        "provenance": _provenance(kipoi),
    }


def _model_names(kipoi: Any, root: Path, models: Sequence[str] | None) -> list[str]:
    if models is None:
        table = _call(
            "list_models",
            kipoi.config.list_models,
            sources={"directory": kipoi.sources.LocalSource(str(root))},
        )
        columns = list(getattr(table, "columns", []))
        names = sorted({str(name) for name in table["model"]}) if "model" in columns else []
    else:
        if isinstance(models, (str, bytes)):
            raise ConfigurationError("models must be a sequence of model directory names.")
        names = []
        for raw in models:
            name = _text(raw, "model")
            relative = Path(name)
            if relative.is_absolute() or ".." in relative.parts or not (root / relative).is_dir():
                raise ConfigurationError(f"model {name!r} is not a directory inside {root}")
            names.append(relative.as_posix())
    if not names:
        raise ConfigurationError(f"no models were found in {root}")
    return names


def test_kipoi_source_models(
    source_dir: str,
    *,
    models: Sequence[str] | None = None,
    batch_size: int = 4,
    stop_on_failure: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run ``kipoi test`` on the models of a local model source in the current Python environment.

    ``source_dir`` is a directory whose subfolders hold ``model.yaml`` files. ``models`` names the subfolders to test;
    when it is omitted every model found in the directory is tested. Each model is tested with
    ``test_kipoi_model_in_current_env`` (kipoi test with ``--source dir``, WARNING records count as failures), and
    failures are collected in the result instead of stopping the run, unless ``stop_on_failure`` is true. Nothing is
    installed and no conda environment is created: a missing kipoi raises ``BackendUnavailableError``, and the models'
    dependencies must already be importable.
    """
    from dnakit.predictions.kipoi_precondition_b_functions import test_kipoi_model_in_current_env

    kipoi = _kipoi()
    root = _existing_dir(source_dir, "source_dir")
    size = _count(batch_size, "batch_size", 1)
    stop = _flag(stop_on_failure, "stop_on_failure")
    names = _model_names(kipoi, root, models)
    results: list[dict[str, Any]] = []
    failed: list[str] = []
    for name in names:
        try:
            outcome = test_kipoi_model_in_current_env(
                str(root / name), batch_size=size, verbose=False
            )
        except BackendExecutionError as exc:
            results.append({"model": name, "status": "failed", "warnings": None, "error": str(exc)})
            failed.append(name)
            if stop:
                break
            continue
        results.append(
            {"model": name, "status": "passed", "warnings": outcome["warnings"], "error": None}
        )
    return {
        "status": "failed" if failed else "passed",
        "operation": "test_kipoi_source_models",
        "source_dir": str(root),
        "batch_size": size,
        "tested": len(results),
        "passed": len(results) - len(failed),
        "failed": failed,
        "results": results,
        "provenance": _provenance(kipoi),
    }


_register_public_verbose(globals())
