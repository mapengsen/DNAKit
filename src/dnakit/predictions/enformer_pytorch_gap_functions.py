"""enformer-pytorch operations: pretrained weight download and embedding extraction from local weights.

``download_pretrained_enformer`` downloads published Enformer weights from the Hugging Face hub (network) and
saves the state dict and config to local files. ``enformer_embeddings`` rebuilds an Enformer from a local
state dict and a config, then returns the embeddings of the one-hot-encoded DNA sequences as a ``.npy`` file.
Training-mode options follow enformer_pytorch.finetune.get_enformer_embeddings.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Mapping, Sequence

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "enformer"
_DNA = re.compile(r"^[ACGTNacgtn]+$")


def _existing(path: str | Path, name: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise ConfigurationError(f"{name} must be an existing file: {resolved}")
    return resolved


def _new(path: str | Path, name: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if resolved.exists():
        raise ConfigurationError(f"{name} must not exist yet: {resolved}")
    if not resolved.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {resolved.parent}")
    return resolved


def _parameter_count(model: Any) -> int:
    return int(sum(parameter.numel() for parameter in model.parameters()))


def download_pretrained_enformer(
    name: str,
    output_state_dict: str,
    *,
    use_tf_gamma: bool | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download pretrained Enformer weights from the Hugging Face hub and save them locally (enformer_pytorch.from_pretrained).

    This needs network access the first time (the Hugging Face cache is used afterwards). The state dict goes to
    ``output_state_dict`` and the model config to ``<output_state_dict>.config.json``; both paths must be new.
    ``use_tf_gamma`` defaults to True for ``EleutherAI/enformer-official-rough``, as in enformer_pytorch.
    """
    torch = require("torch", _EXTRA)
    modeling = require("enformer_pytorch.modeling_enformer", _EXTRA)
    if not isinstance(name, str) or not name.strip():
        raise ConfigurationError("name must be a Hugging Face model id such as EleutherAI/enformer-official-rough.")
    state_path = _new(output_state_dict, "output_state_dict")
    config_path = state_path.with_name(state_path.name + ".config.json")
    if config_path.exists():
        raise ConfigurationError(f"The config file already exists: {config_path}")
    kwargs = {} if use_tf_gamma is None else {"use_tf_gamma": use_tf_gamma}
    model = modeling.from_pretrained(name, **kwargs)
    torch.save(model.state_dict(), state_path)
    config = model.config.to_dict()
    config_path.write_text(json.dumps(config, indent=2, sort_keys=True), encoding="utf-8")
    return {
        "status": "completed",
        "backend": "enformer_pytorch",
        "model": name,
        "state_dict": str(state_path),
        "config": str(config_path),
        "parameter_count": _parameter_count(model),
        "output_heads": sorted(config.get("output_heads", {})),
    }


def _sequences(sequences: Sequence[str]) -> list[str]:
    if isinstance(sequences, str) or not sequences:
        raise ConfigurationError("sequences must be a non-empty list of DNA strings.")
    values = list(sequences)
    if not all(isinstance(value, str) and _DNA.match(value) for value in values):
        raise ConfigurationError("sequences may contain only A, C, G, T and N.")
    if len({len(value) for value in values}) != 1:
        raise ConfigurationError("all sequences must have the same length.")
    return values


def enformer_embeddings(
    state_dict: str,
    sequences: Sequence[str],
    output_npy: str,
    *,
    config: Mapping[str, Any],
    freeze: bool = False,
    train_layernorms_only: bool = False,
    train_last_n_layers_only: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract Enformer embeddings for DNA sequences from a local state dict (enformer_pytorch.finetune.get_enformer_embeddings).

    ``config`` holds the EnformerConfig arguments the weights were trained with (the saved
    ``*.config.json`` from ``download_pretrained_enformer`` can be read and passed here). Sequences must have one
    common length; bases are encoded as in enformer_pytorch (A, C, G, T, N). ``freeze`` disables gradients,
    ``train_layernorms_only`` and ``train_last_n_layers_only`` select the trainable layers. The embeddings are
    written to ``output_npy`` (a new ``.npy`` file) and the shape is returned.
    """
    torch = require("torch", _EXTRA)
    numpy = require("numpy", "bio-analysis")
    enformer = require("enformer_pytorch", _EXTRA)
    finetune = require("enformer_pytorch.finetune", _EXTRA)
    weights = _existing(state_dict, "state_dict")
    values = _sequences(sequences)
    out = _new(output_npy, "output_npy")
    if not isinstance(config, Mapping) or not config:
        raise ConfigurationError("config must be a non-empty mapping of EnformerConfig arguments.")
    if freeze and (train_layernorms_only or train_last_n_layers_only is not None):
        raise ConfigurationError("freeze cannot be combined with train_layernorms_only or train_last_n_layers_only.")
    if train_last_n_layers_only is not None and train_last_n_layers_only < 0:
        raise ConfigurationError("train_last_n_layers_only must not be negative.")
    try:
        model = enformer.Enformer(enformer.EnformerConfig(**dict(config)))
    except TypeError as exc:
        raise ConfigurationError(f"config does not match EnformerConfig: {exc}") from exc
    state = torch.load(weights, map_location="cpu", weights_only=True)
    try:
        model.load_state_dict(state)
    except RuntimeError as exc:
        raise ConfigurationError(f"state_dict does not match the config: {exc}") from exc
    model.eval()
    seq = enformer.data.str_to_seq_indices(values)
    embeddings = finetune.get_enformer_embeddings(
        model,
        seq,
        freeze=freeze,
        train_layernorms_only=train_layernorms_only,
        train_last_n_layers_only=train_last_n_layers_only,
    )
    array = embeddings.detach().cpu().numpy()
    numpy.save(out, array)
    return {
        "status": "completed",
        "backend": "enformer_pytorch",
        "output": str(out),
        "shape": [int(dim) for dim in array.shape],
        "sequence_count": len(values),
        "sequence_length": len(values[0]),
        "freeze": freeze,
        "train_layernorms_only": train_layernorms_only,
        "train_last_n_layers_only": train_last_n_layers_only,
    }


__all__ = [
    "download_pretrained_enformer",
    "enformer_embeddings",
]
_register_public_verbose(globals())
