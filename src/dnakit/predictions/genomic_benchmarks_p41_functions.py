"""Local genomic-benchmarks CNN training, evaluation, prediction and TF export.

Version note: adds real upstream CNN operations with weights-only checkpoints.
Token IDs are explicit DNA vocabulary indices; no implicit sequence recoding.
Upstream CNN training and evaluation defaults (including binary activation and
loss) are retained, even where they may be surprising.
"""

from __future__ import annotations

import contextlib
import io
import json
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError


def _count(value: Any, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _config(config: Mapping[str, Any]) -> dict[str, int]:
    limits = {
        "number_of_classes": (2, 1000),
        "vocab_size": (2, 256),
        "embedding_dim": (1, 1024),
        "input_len": (64, 100000),
    }
    if set(config) != set(limits):
        raise ConfigurationError(
            "config requires number_of_classes, vocab_size, embedding_dim, input_len."
        )
    return {k: _count(config[k], k, *limit) for k, limit in limits.items()}


def _output(path: str) -> Path:
    result = Path(path).expanduser().resolve()
    if result.exists() or not result.parent.is_dir():
        raise ConfigurationError("Output must be new and its parent must exist.")
    return result


def _model(config: Mapping[str, Any], checkpoint: str | None = None) -> tuple[Any, Any]:
    torch = require("torch", "bio-genomic-benchmarks")
    module = require("genomic_benchmarks.models.torch", "bio-genomic-benchmarks")
    obj = module.CNN(**config, device="cpu")
    if checkpoint is not None:
        obj.load_state_dict(
            torch.load(str(input_path(checkpoint)), map_location="cpu", weights_only=True)
        )
    return torch, obj


def _data(
    torch: Any,
    token_ids: Sequence[Sequence[int]],
    labels: Sequence[int] | None,
    config: Mapping[str, int],
    batch_size: int,
) -> Any:
    _count(batch_size, "batch_size", 1, 10000)
    np = require("numpy", "bio-genomic-benchmarks")
    try:
        arr = np.asarray(token_ids)
    except ValueError as exc:
        raise ConfigurationError("token_ids must be a rectangular integer matrix.") from exc
    if (
        arr.ndim != 2
        or not 1 <= arr.shape[0] <= 100000
        or arr.shape[1] != config["input_len"]
        or arr.dtype.kind not in "iu"
        or arr.min() < 0
        or arr.max() >= config["vocab_size"]
    ):
        raise ConfigurationError("token_ids must match input_len and the vocabulary range.")
    x = torch.tensor(arr, dtype=torch.long)
    if labels is None:
        return x
    yarr = np.asarray(labels)
    if (
        yarr.shape != (len(arr),)
        or yarr.dtype.kind not in "iu"
        or yarr.min() < 0
        or yarr.max() >= config["number_of_classes"]
    ):
        raise ConfigurationError("One integer class label per sequence is required.")
    y = torch.tensor(yarr, dtype=torch.float32).reshape(-1, 1)
    return torch.utils.data.DataLoader(
        torch.utils.data.TensorDataset(x, y), batch_size=batch_size, shuffle=False
    )


def train_genomic_benchmark_cnn(
    token_ids: Sequence[Sequence[int]],
    labels: Sequence[int],
    config: Mapping[str, Any],
    output_checkpoint: str,
    *,
    epochs: int = 1,
    batch_size: int = 32,
    seed: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call upstream CNN.fit on local labelled DNA token batches and save state_dict.

    Training runs on CPU with Adam's upstream defaults. The local RNG context
    makes initialisation reproducible without changing the caller's RNG state.
    Upstream epoch progress is returned in training_log.
    """
    cfg = _config(config)
    _count(epochs, "epochs", 1, 10000)
    _count(seed, "seed", 0, 2**32 - 1)
    path = _output(output_checkpoint)
    torch = require("torch", "bio-genomic-benchmarks")
    log = io.StringIO()
    with torch.random.fork_rng(devices=[]):
        torch.manual_seed(seed)
        torch, obj = _model(cfg)
        loader = _data(torch, token_ids, labels, cfg, batch_size)
        with contextlib.redirect_stdout(log):
            obj.fit(loader, epochs)
        torch.save(obj.state_dict(), path)
    if verbose:
        print(log.getvalue(), end="")
    return {
        "status": "completed",
        "outputs": [str(path)],
        "config": cfg,
        "epochs": epochs,
        "training_log": log.getvalue()[-16000:],
        "backend": "genomic_benchmarks.models.torch.CNN.fit",
    }


def evaluate_genomic_benchmark_cnn(
    token_ids: Sequence[Sequence[int]],
    labels: Sequence[int],
    config: Mapping[str, Any],
    checkpoint: str,
    *,
    batch_size: int = 32,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return upstream CNN.test accuracy, F1 and its printed average loss.

    Loss is returned at the upstream printed six-decimal precision.
    """
    cfg = _config(config)
    torch, obj = _model(cfg, checkpoint)
    loader = _data(torch, token_ids, labels, cfg, batch_size)
    log = io.StringIO()
    with contextlib.redirect_stdout(log):
        accuracy, f1 = obj.test(loader)
    match = re.search(r"Avg loss:\s*([-+0-9.eE]+)", log.getvalue())
    if verbose:
        print(log.getvalue(), end="")
    return {
        "status": "completed",
        "accuracy": float(accuracy),
        "f1": float(f1),
        "loss": float(match[1]) if match else None,
        "evaluation_log": log.getvalue(),
        "backend": "genomic_benchmarks.models.torch.CNN.test",
    }


def predict_genomic_benchmark_cnn(
    token_ids: Sequence[Sequence[int]],
    config: Mapping[str, Any],
    checkpoint: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run upstream CNN.forward in evaluation mode and return class scores.

    Binary scores are sigmoid outputs; multiclass scores retain upstream logits.
    """
    cfg = _config(config)
    torch, obj = _model(cfg, checkpoint)
    x = _data(torch, token_ids, None, cfg, 32)
    obj.eval()
    with torch.inference_mode():
        values = obj(x).cpu().tolist()
    return {
        "status": "completed",
        "scores": values,
        "backend": "genomic_benchmarks.models.torch.CNN.forward",
    }


def build_genomic_benchmark_tf_cnn(
    number_of_classes: int,
    vocab_size: int,
    output_path: str,
    *,
    input_len: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build, compile and save upstream get_basic_cnn_model_v0 as a Keras file.

    Requires the upstream compatible TensorFlow and tensorflow-addons versions.
    This exports an initialised architecture; it does not train the model.
    """
    _count(number_of_classes, "number_of_classes", 2, 1000)
    _count(vocab_size, "vocab_size", 2, 256)
    _count(input_len, "input_len", 64, 100000)
    path = _output(output_path)
    if path.suffix != ".keras":
        raise ConfigurationError("output_path must end with .keras.")
    module = require("genomic_benchmarks.models.tf", "bio-genomic-benchmarks-tf")
    obj = module.get_basic_cnn_model_v0(number_of_classes, vocab_size)
    obj.build((None, input_len))
    obj.save(str(path))
    return {
        "status": "completed",
        "outputs": [str(path)],
        "parameter_count": int(obj.count_params()),
        "model_config": json.loads(obj.to_json()),
        "trained": False,
    }


__all__ = [
    "build_genomic_benchmark_tf_cnn",
    "evaluate_genomic_benchmark_cnn",
    "predict_genomic_benchmark_cnn",
    "train_genomic_benchmark_cnn",
]
_register_public_verbose(globals(), tuple(__all__))
