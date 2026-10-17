"""Model-facing operations of ``ml_models``: thin, checked delegations to ``dnakit.predictions`` etc."""

from __future__ import annotations

import dataclasses
import json
import os
from collections.abc import Mapping
from pathlib import Path

from ._impl_common import Unavailable, bad_input, bad_option, get_int, get_str

_ENVELOPE_KEYS = {"backend", "verbose"}
_EXCLUDED_CONFIG_KEYS = {"hf_token", "llm", "allow_remote_code", "show_progress"}


def _source_path(value: object, what: str = "data file") -> Path:
    if not isinstance(value, (str, os.PathLike)) or not str(value).strip() or "\0" in str(value):
        raise bad_input(f"input must be the path of a {what}.", code="INVALID_ML_INPUT")
    path = Path(value).expanduser()
    if not path.is_file():
        raise bad_input(f"{what} does not exist: {path}", code="ML_INPUT_NOT_FOUND")
    return path


def _mapping(value: object) -> object:
    from dnakit.predictions import DataMapping

    if value is None or isinstance(value, DataMapping):
        return value
    if not isinstance(value, Mapping):
        raise bad_option("option 'mapping' must be a DataMapping or a mapping of its fields.")
    kwargs = {k: (tuple(v) if isinstance(v, (list, tuple)) else v) for k, v in value.items()}
    return DataMapping(**kwargs)


def build_config(options: Mapping[str, object], *, method: str, require_output: bool):
    """Build a ``CustomPropertyTrainingConfig`` from whitelisted option keys."""

    from dnakit.predictions import CustomPropertyTrainingConfig

    allowed = {
        f.name for f in dataclasses.fields(CustomPropertyTrainingConfig)
    } - _EXCLUDED_CONFIG_KEYS
    unknown = sorted(set(options) - allowed - _ENVELOPE_KEYS - {"trust_model"})
    if unknown:
        raise bad_option(
            f"unknown option(s) {unknown}.",
            code="UNKNOWN_ML_OPTION",
            hint=f"Supported keys: {sorted(allowed)}.",
        )
    kwargs = {k: v for k, v in options.items() if k in allowed}
    if "task_type" not in kwargs:
        raise bad_option("option 'task_type' ('classification' or 'regression') is required.")
    if require_output and "output_dir" not in kwargs:
        raise bad_option("option 'output_dir' is required (the model directory is written there).")
    kwargs["finetune_method"] = method
    kwargs["mapping"] = _mapping(kwargs.get("mapping"))
    if method == "machine_learning":
        kwargs.setdefault("model", "random_forest")
        kwargs.setdefault("workers", 1)
    return CustomPropertyTrainingConfig(**kwargs, show_progress=False)


def train(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Train a classical machine-learning property model via ``train_custom_property_model``."""

    method = str(options.get("finetune_method", "machine_learning"))
    if method != "machine_learning":
        raise Unavailable(
            "pretrained foundation-model weights and the optional neural runtime",
            "ml_models.train trains classical ML estimators; neural fine-tuning belongs to ml_models.finetune.",
        )
    from dnakit.predictions import train_custom_property_model

    config = build_config(
        {k: v for k, v in options.items() if k != "finetune_method"},
        method=method,
        require_output=True,
    )
    result = train_custom_property_model(_source_path(value), config=config, verbose=False)
    return {
        "delegate": "dnakit.predictions.train_custom_property_model",
        "training": _jsonable(result),
    }


def _jsonable(result: object) -> object:
    to_dict = getattr(result, "to_dict", None)
    payload = to_dict() if callable(to_dict) else dataclasses.asdict(result)  # type: ignore[call-overload]
    return json.loads(json.dumps(payload, default=str))


def predict(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Predict with a saved DNAKit property model via ``predict_custom_properties``."""

    from dnakit.predictions import predict_custom_properties

    model_dir = options.get("model_dir")
    if model_dir is None:
        raise bad_option("option 'model_dir' (a directory written by ml_models.train) is required.")
    allowed = {
        "model_dir",
        "trust_model",
        "mapping",
        "threshold",
        "batch_size",
        "output_path",
        "overwrite",
        "device",
    }
    unknown = sorted(set(options) - allowed - _ENVELOPE_KEYS)
    if unknown:
        raise bad_option(
            f"unknown option(s) {unknown}.",
            code="UNKNOWN_ML_OPTION",
            hint=f"Supported keys: {sorted(allowed)}.",
        )
    kwargs = {
        k: options[k]
        for k in ("threshold", "batch_size", "output_path", "overwrite", "device")
        if k in options
    }
    result = predict_custom_properties(
        _source_path(value),
        model_dir=model_dir,
        mapping=_mapping(options.get("mapping")),  # type: ignore[arg-type]
        trust_model=options.get("trust_model", False) is True,
        show_progress=False,
        verbose=False,
        **kwargs,  # type: ignore[arg-type]
    )
    return {
        "delegate": "dnakit.predictions.predict_custom_properties",
        "prediction": _jsonable(result),
    }


def build(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Describe the untrained scikit-learn style estimator selected by ``model`` and ``model_params``."""

    from dnakit.predictions.custom_training.ml_models import create_estimator

    spec = dict(options)
    if isinstance(value, Mapping):
        spec = {**value, **spec}
    spec.setdefault("model", "random_forest")
    config = build_config(spec, method="machine_learning", require_output=False)
    estimator = create_estimator(config)
    params = estimator.get_params() if hasattr(estimator, "get_params") else {}
    return {
        "model": config.model,
        "task_type": config.task_type,
        "estimator_class": f"{type(estimator).__module__}.{type(estimator).__qualname__}",
        "params": json.loads(json.dumps(params, default=repr)),
        "trained": False,
    }


def load(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Read the manifest of a saved DNAKit property model directory (no weights are deserialised)."""

    if not isinstance(value, (str, os.PathLike)) or not str(value).strip():
        raise bad_input("input must be a model directory written by ml_models.train.")
    root = Path(value).expanduser()
    manifest_path = root / "model.json"
    if not manifest_path.is_file():
        raise bad_input(f"no model.json in {root}", code="ML_MODEL_NOT_FOUND")
    try:
        if manifest_path.stat().st_size > 10_000_000:
            raise bad_input("model.json exceeds 10 MB.", code="ML_INPUT_TOO_LARGE")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise bad_input(
            f"could not parse model.json: {exc}", code="INVALID_CUSTOM_MODEL_ARTIFACT"
        ) from exc
    if not isinstance(manifest, dict):
        raise bad_input(
            "model.json must contain a JSON object.", code="INVALID_CUSTOM_MODEL_ARTIFACT"
        )
    is_ml = manifest.get("finetune_method") == "machine_learning"
    artifact = "estimator.joblib" if is_ml else "head.npz"
    return {
        "model_dir": str(root),
        "manifest": manifest,
        "artifact": artifact,
        "artifact_present": (root / artifact).is_file(),
        "files": sorted(p.name for p in root.iterdir())[:200],
        "weights_loaded": False,
        "note": "Only the JSON manifest was read; use ml_models.predict with trust_model=True to run the model.",
    }


def search(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Search the catalogue of models known to DNAKit (ML estimators, fine-tunable and embedding DNA models, predictors)."""

    from dnakit.predictions import (
        available_custom_training_models,
        available_machine_learning_models,
        available_prediction_models,
    )
    from dnakit.representations import available_embedding_models, get_embedding_model

    query = value if isinstance(value, str) else options.get("query", "")
    if not isinstance(query, str):
        raise bad_input("input (or options['query']) must be a search string.")
    kind = get_str(options, "kind", None)
    task = get_str(options, "task", None)
    limit = get_int(options, "limit", 50, low=1, high=1000)
    entries: list[dict[str, object]] = []
    for m in available_machine_learning_models(verbose=False):
        entries.append(
            {
                "name": m.name,
                "kind": "machine_learning_estimator",
                "tasks": list(m.tasks),
                "description": f"{m.dependency} estimator",
            }
        )
    for m in available_custom_training_models(verbose=False):
        entries.append(
            {
                "name": m.name,
                "kind": "dna_model_finetuning",
                "tasks": ["classification", "regression"],
                "description": f"{m.display_name}. {m.notes}",
            }
        )
    for name in available_embedding_models(verbose=False):
        spec = get_embedding_model(name, verbose=False)
        entries.append(
            {
                "name": name,
                "kind": "dna_embedding_model",
                "tasks": ["embedding"],
                "description": f"{spec.display_name}. {spec.notes or ''}".strip(),
            }
        )
    for name in available_prediction_models(verbose=False):
        entries.append(
            {
                "name": str(name),
                "kind": "pretrained_prediction_model",
                "tasks": ["prediction"],
                "description": "pretrained property/variant-effect prediction backend",
            }
        )
    kinds = sorted({str(e["kind"]) for e in entries})
    if kind is not None and kind not in kinds:
        raise bad_option(f"unknown kind '{kind}'.", hint=f"Use one of {kinds}.")
    terms = query.lower().split()
    matches = []
    for e in entries:
        haystack = f"{e['name']} {e['kind']} {e['description']}".lower()
        if (
            (kind is None or e["kind"] == kind)
            and (task is None or task in e["tasks"])
            and all(t in haystack for t in terms)
        ):  # type: ignore[operator]
            matches.append(e)
    return {"query": query, "kinds": kinds, "n_matches": len(matches), "matches": matches[:limit]}


def embeddings(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Embed sequences with a DNAKit representation model that has a local checkpoint (or a supplied backend)."""

    from dnakit.core import DNARecord, DNASequence
    from dnakit.representations import RepresentationConfig, extract_representations

    from ._impl_common import sequences_from

    backend = options.get("representation_backend")
    if backend is None and options.get("checkpoint_path") is None:
        raise Unavailable(
            "pretrained DNA foundation-model weights (options['checkpoint_path']) or options['representation_backend']",
            "Sequence embeddings come from a DNA language model; no local checkpoint or backend was supplied.",
        )
    pairs = sequences_from(
        value, options, limit=get_int(options, "max_sequences", 10_000, low=1, high=100_000)
    )
    keys = (
        "model",
        "checkpoint_path",
        "model_source_path",
        "pooling",
        "ambiguity_policy",
        "device",
        "dtype",
        "batch_size",
        "max_length",
    )
    kwargs = {k: options[k] for k in keys if k in options}
    config = RepresentationConfig(show_progress=False, **kwargs)  # type: ignore[arg-type]
    records = [
        DNARecord(id=name, sequence=DNASequence(seq.upper(), alphabet="iupac"))
        for name, seq in pairs
    ]
    result = extract_representations(records, config=config, backend=backend, verbose=False)  # type: ignore[arg-type]
    matrix = result.representations
    return {
        "delegate": "dnakit.representations.extract_representations",
        "model": result.model_name,
        "pooling": result.pooling,
        "dimension": result.embedding_dimension,
        "record_ids": list(result.record_ids),
        "embeddings": [[float(x) for x in row] for row in matrix],
    }


def finetune(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Fine-tune (linear probe or full) a DNA model on a labelled file via ``train_custom_property_model``."""

    from dnakit.predictions import train_custom_property_model

    backend = options.get("representation_backend")
    encoder = options.get("training_encoder")
    if backend is None and encoder is None and options.get("checkpoint_path") is None:
        raise Unavailable(
            "pretrained DNA foundation-model weights (options['checkpoint_path']) or a representation_backend/training_encoder",
            "Fine-tuning needs a pretrained model; no local checkpoint or encoder was supplied (DNAKit does not download during this call).",
        )
    method = str(options.get("finetune_method", "linear_probe"))
    if method not in {"linear_probe", "full"}:
        raise bad_option("finetune_method must be 'linear_probe' or 'full' for ml_models.finetune.")
    config_options = {
        k: v
        for k, v in options.items()
        if k not in {"finetune_method", "representation_backend", "training_encoder"}
    }
    config = build_config(config_options, method=method, require_output=True)
    result = train_custom_property_model(
        _source_path(value),
        config=config,
        representation_backend=backend,
        training_encoder=encoder,  # type: ignore[arg-type]
        verbose=False,
    )
    return {
        "delegate": "dnakit.predictions.train_custom_property_model",
        "training": _jsonable(result),
    }


def export_torchscript(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Trace a torch module with an example input and save it as a TorchScript archive."""

    try:
        import torch
    except ImportError:
        raise Unavailable("pip install torch", "TorchScript export needs PyTorch.") from None
    module = options.get("model", value)
    if not isinstance(module, torch.nn.Module):
        raise bad_input(
            "input (or options['model']) must be a torch.nn.Module.", code="ML_NOT_A_TORCH_MODULE"
        )
    output = options.get("output")
    if output is None:
        raise bad_option("option 'output' (the .pt path to write) is required.")
    example = options.get("example_input")
    if example is None:
        raise bad_option("option 'example_input' (nested list of numbers) is required for tracing.")
    try:
        tensor = torch.tensor(example, dtype=torch.float32)
    except (TypeError, ValueError, RuntimeError) as exc:
        raise bad_option(
            f"example_input must be a rectangular nested list of numbers: {exc}"
        ) from exc
    path = Path(str(output)).expanduser()
    if path.exists() and not options.get("overwrite", False):
        raise bad_option(
            f"output file already exists: {path}",
            code="ML_OUTPUT_EXISTS",
            hint="Set options['overwrite']=True or choose another path.",
        )
    was_training = module.training
    module.eval()
    try:
        with torch.no_grad():
            traced = torch.jit.trace(module, tensor)
            path.parent.mkdir(parents=True, exist_ok=True)
            traced.save(str(path))
            reloaded = torch.jit.load(str(path))
            difference = float((module(tensor) - reloaded(tensor)).abs().max())
    finally:
        module.train(was_training)
    return {
        "output_path": str(path),
        "bytes": path.stat().st_size,
        "input_shape": list(tensor.shape),
        "roundtrip_max_abs_difference": difference,
        "torch_version": str(torch.__version__),
    }
