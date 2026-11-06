"""Local Enformer prediction and upstream track metrics (integration P40).

Weights are supplied explicitly; no pretrained model or data are downloaded. Adapter
checkpoints must contain the entire adapter state, including its Enformer backbone.
"""

from __future__ import annotations

import importlib.metadata
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p40_helpers as _h
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError
from dnakit.predictions.enformer_pytorch_gap_functions import _sequences


def _tensor(value: Any) -> Any:
    torch = require("torch", "enformer")
    tensor = torch.as_tensor(value, dtype=torch.float32)
    if tensor.numel() == 0 or not torch.isfinite(tensor).all():
        raise ConfigurationError("Array values must be nonempty and finite.")
    return tensor


def _result(value: Any) -> dict[str, Any]:
    def plain(item: Any) -> Any:
        return [plain(v) for v in item] if isinstance(item, list) else _h.plain(item)

    return {
        "status": "completed",
        "backend": "enformer_pytorch",
        "version": importlib.metadata.version("enformer-pytorch"),
        "values": plain(value.detach().cpu().tolist()),
        "shape": list(value.shape),
    }


def enformer_poisson_loss(
    predictions: list[Any], targets: list[Any], *, verbose: bool = True
) -> dict[str, Any]:
    """Compute upstream poisson_loss, including its epsilon and mean reduction."""
    module = require("enformer_pytorch.modeling_enformer", "enformer")
    pred, target = _tensor(predictions), _tensor(targets)
    if pred.shape != target.shape or (pred < 0).any() or (target < 0).any():
        raise ConfigurationError(
            "Predictions and targets must have equal shapes and nonnegative values."
        )
    return _result(module.poisson_loss(pred, target))


def enformer_pearson_correlation(
    predictions: list[Any], targets: list[Any], *, verbose: bool = True
) -> dict[str, Any]:
    """Compute upstream per-sequence Pearson coefficient for batch x positions x tracks."""
    module = require("enformer_pytorch.modeling_enformer", "enformer")
    pred, target = _tensor(predictions), _tensor(targets)
    if pred.shape != target.shape or pred.ndim != 3:
        raise ConfigurationError("Arrays must have the same batch x positions x tracks shape.")
    return _result(module.pearson_corr_coef(pred, target))


def enformer_channel_correlation(
    prediction_batches: list[Any], target_batches: list[Any], *, verbose: bool = True
) -> dict[str, Any]:
    """Accumulate upstream MeanPearsonCorrCoefPerChannel over successive track batches."""
    module = require("enformer_pytorch.metrics", "enformer")
    if not prediction_batches or len(prediction_batches) != len(target_batches):
        raise ConfigurationError("Provide equally many nonempty prediction and target batches.")
    first = _tensor(prediction_batches[0])
    if first.ndim != 3:
        raise ConfigurationError("Batches must have batch x positions x tracks shapes.")
    metric = module.MeanPearsonCorrCoefPerChannel(first.shape[-1])
    for pred, target in zip(prediction_batches, target_batches, strict=True):
        p, t = _tensor(pred), _tensor(target)
        if p.shape != t.shape or p.ndim != 3 or p.shape[-1] != first.shape[-1]:
            raise ConfigurationError(
                "Each batch must have equal batch x positions x tracks shapes."
            )
        metric.update(p, t)
    return _result(metric.compute())


def _predict(
    state_dict: str,
    sequences: Sequence[str],
    output_npy: str,
    config: Mapping[str, Any],
    adapter: str | None,
    adapter_config: Mapping[str, Any] | None,
    context: list[Any] | None,
    head: str | None,
) -> dict[str, Any]:
    torch, np = require("torch", "enformer"), require("numpy", "bio-analysis")
    enformer = require("enformer_pytorch", "enformer")
    values = _sequences(sequences)
    weights = _h.existing(state_dict, "state_dict")
    output = _h.new_file(output_npy, "output_npy")
    if output.suffix != ".npy":
        raise ConfigurationError("output_npy must end in .npy.")
    model = enformer.Enformer(enformer.EnformerConfig(**dict(config)))
    if adapter is not None:
        finetune = require("enformer_pytorch.finetune", "enformer")
        allowed = {
            "HeadAdapterWrapper": {"num_tracks", "post_transformer_embed"},
            "ContextAdapterWrapper": {"context_dim"},
            "ContextAttentionAdapterWrapper": {"context_dim", "heads", "dim_head"},
        }
        options = dict(adapter_config or {})
        if set(options) - allowed[adapter]:
            raise ConfigurationError(
                f"Unsupported adapter options: {sorted(set(options) - allowed[adapter])}"
            )
        model = getattr(finetune, adapter)(enformer=model, **options)
    model.load_state_dict(torch.load(weights, weights_only=True, map_location="cpu"))
    model.eval()
    kwargs = {}
    if context is not None:
        kwargs["context"] = _tensor(context)
    if head is not None:
        kwargs["head"] = head
    with torch.no_grad():
        predictions = model(enformer.data.str_to_seq_indices(values), **kwargs)
    if isinstance(predictions, dict):
        raise ConfigurationError("Select one output head, for example human or mouse.")
    np.save(output, predictions.detach().cpu().numpy())
    return {
        "status": "completed",
        "backend": "enformer_pytorch",
        "outputs": [str(output)],
        "shape": list(predictions.shape),
        "version": importlib.metadata.version("enformer-pytorch"),
    }


def enformer_predict_tracks(
    state_dict: str,
    sequences: Sequence[str],
    output_npy: str,
    *,
    config: Mapping[str, Any],
    head: str = "human",
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict the human or mouse tracks using a complete local Enformer checkpoint."""
    return _predict(state_dict, sequences, output_npy, config, None, None, None, head)


def enformer_predict_custom_tracks(
    state_dict: str,
    sequences: Sequence[str],
    output_npy: str,
    *,
    config: Mapping[str, Any],
    num_tracks: int,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict user-defined tracks with a trained HeadAdapterWrapper checkpoint."""
    _h.integer(num_tracks, "num_tracks", minimum=1)
    return _predict(
        state_dict,
        sequences,
        output_npy,
        config,
        "HeadAdapterWrapper",
        {"num_tracks": num_tracks},
        None,
        None,
    )


def enformer_predict_context_tracks(
    state_dict: str,
    sequences: Sequence[str],
    output_npy: str,
    *,
    config: Mapping[str, Any],
    context: list[Any],
    context_dim: int,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict context-specific tracks using trained hypernetwork adapter weights."""
    _h.integer(context_dim, "context_dim", minimum=1)
    return _predict(
        state_dict,
        sequences,
        output_npy,
        config,
        "ContextAdapterWrapper",
        {"context_dim": context_dim},
        context,
        None,
    )


def enformer_predict_context_attention(
    state_dict: str,
    sequences: Sequence[str],
    output_npy: str,
    *,
    config: Mapping[str, Any],
    context: list[Any],
    context_dim: int,
    heads: int = 8,
    dim_head: int = 64,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict tracks with trained cross-attention context-token adapter weights."""
    for name, value in (("context_dim", context_dim), ("heads", heads), ("dim_head", dim_head)):
        _h.integer(value, name, minimum=1)
    return _predict(
        state_dict,
        sequences,
        output_npy,
        config,
        "ContextAttentionAdapterWrapper",
        {"context_dim": context_dim, "heads": heads, "dim_head": dim_head},
        context,
        None,
    )


__all__ = [
    "enformer_channel_correlation",
    "enformer_pearson_correlation",
    "enformer_poisson_loss",
    "enformer_predict_context_attention",
    "enformer_predict_context_tracks",
    "enformer_predict_custom_tracks",
    "enformer_predict_tracks",
]

_register_public_verbose(globals(), tuple(__all__))
