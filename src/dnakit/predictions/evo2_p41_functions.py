"""Local Evo 2 scoring, forward embeddings and per-position uncertainty.

Version note: adds real upstream entry points, using an explicit local checkpoint
and bounded inputs. No weights are downloaded. CUDA and a compatible Evo2/Vortex
installation are required for real model execution.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from importlib.metadata import version
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError
from dnakit.predictions.dna_model_backends import _EVO2_MODELS


def _sequences(sequences: Sequence[str]) -> list[str]:
    if isinstance(sequences, (str, bytes)) or not 1 <= len(sequences) <= 1000:
        raise ConfigurationError("sequences requires 1-1000 DNA strings.")
    values = list(sequences)
    if any(
        not isinstance(s, str) or not 2 <= len(s) <= 100000 or not re.fullmatch("[ACGTNacgtn]+", s)
        for s in values
    ):
        raise ConfigurationError("Each DNA sequence must contain 2-100000 A/C/G/T/N characters.")
    if sum(map(len, values)) > 1000000:
        raise ConfigurationError("The total input is limited to 1000000 bases.")
    return values


def _load(checkpoint: str, model_name: str) -> Any:
    path = input_path(checkpoint)
    if model_name not in _EVO2_MODELS:
        raise ConfigurationError("Unsupported Evo 2 model name.")
    evo = require("evo2", "neural-evo2")
    return evo.Evo2(model_name=model_name, local_path=str(path))


def _plain(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if hasattr(value, "detach"):
        value = value.detach().float().cpu()
    if hasattr(value, "tolist"):
        return value.tolist()
    return value


def _result(**values: Any) -> dict[str, Any]:
    return {
        "status": "completed",
        **values,
        "provenance": {"backend": "evo2", "version": version("evo2")},
    }


def score_evo2_sequences(
    sequences: Sequence[str],
    checkpoint: str,
    *,
    model_name: str = "evo2_7b",
    batch_size: int = 1,
    prepend_bos: bool = False,
    reduce_method: str = "mean",
    average_reverse_complement: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Score local DNA sequences with upstream log-likelihood, optionally both strands.

    Upstream defaults (mean reduction and no BOS) are preserved. Lowercase
    tokens are preserved; callers choose their own model preprocessing.
    """
    seqs = _sequences(sequences)
    if (
        isinstance(batch_size, bool)
        or not isinstance(batch_size, int)
        or not 1 <= batch_size <= 1000
        or reduce_method not in {"mean", "sum"}
    ):
        raise ConfigurationError("batch_size must be 1-1000 and reduce_method mean or sum.")
    obj = _load(checkpoint, model_name)
    values = obj.score_sequences(
        seqs,
        batch_size=batch_size,
        prepend_bos=prepend_bos,
        reduce_method=reduce_method,
        average_reverse_complement=average_reverse_complement,
    )
    return _result(scores=_plain(values), average_reverse_complement=average_reverse_complement)


def forward_evo2_sequences(
    sequences: Sequence[str],
    checkpoint: str,
    *,
    model_name: str = "evo2_7b",
    layer_names: Sequence[str] | None = None,
    prepend_bos: bool = False,
    device: str = "cuda:0",
    max_output_values: int = 1000000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return upstream forward logits and requested layer embeddings as JSON arrays.

    Inputs are padded with the upstream tokenizer. sequence_lengths identifies
    real bases; returned logits include padding/BOS. Large results are rejected
    explicitly, rather than silently truncated.
    """
    seqs = _sequences(sequences)
    if (
        isinstance(max_output_values, bool)
        or not isinstance(max_output_values, int)
        or not 1 <= max_output_values <= 10000000
    ):
        raise ConfigurationError("max_output_values must be 1-10000000.")
    if layer_names is not None and (
        isinstance(layer_names, (str, bytes))
        or any(not isinstance(x, str) or not x for x in layer_names)
    ):
        raise ConfigurationError("layer_names must be a list of module names.")
    if len(seqs) * (max(map(len, seqs)) + int(prepend_bos)) * 512 > max_output_values:
        raise ConfigurationError("Logit output exceeds max_output_values; use shorter inputs.")
    obj = _load(checkpoint, model_name)
    scoring = require("evo2.scoring", "neural-evo2")
    ids, lengths = scoring.prepare_batch(
        seqs, obj.tokenizer, prepend_bos=prepend_bos, device=device
    )
    logits, embeddings = obj.forward(
        ids,
        return_embeddings=bool(layer_names),
        layer_names=None if layer_names is None else list(layer_names),
    )
    # Vortex forward may return (logits, inference-cache); the cache is not a prediction.
    if isinstance(logits, tuple):
        logits = logits[0]
    total = int(logits.numel()) + sum(int(v.numel()) for v in (embeddings or {}).values())
    if total > max_output_values:
        raise ConfigurationError("Logits plus embeddings exceed max_output_values.")
    return _result(logits=_plain(logits), embeddings=_plain(embeddings), sequence_lengths=lengths)


def evo2_positional_uncertainty(
    sequence: str,
    checkpoint: str,
    *,
    model_name: str = "evo2_7b",
    metric: str = "entropy",
    reverse_complement: bool = False,
    prepend_bos: bool = False,
    device: str = "cuda:0",
    verbose: bool = True,
) -> dict[str, Any]:
    """Return per-position predictive entropy or perplexity using upstream entropy.

    Both-strand mode reverses the reverse-strand positions, averages entropy,
    then exponentiates when perplexity is requested, as upstream
    score_perplexity_along_sequence does. device and prepend_bos are explicit.
    """
    seq = _sequences([sequence])[0]
    if metric not in {"entropy", "perplexity"}:
        raise ConfigurationError("metric must be entropy or perplexity.")
    obj = _load(checkpoint, model_name)
    scoring = require("evo2.scoring", "neural-evo2")
    np = require("numpy", "neural-evo2")
    seqs = [seq]
    if reverse_complement:
        seqs.append(seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1])
    rows = scoring.positional_entropies(
        seqs, obj.model, obj.tokenizer, prepend_bos=prepend_bos, device=device
    )
    values = (rows[0] + rows[1][::-1]) * 0.5 if reverse_complement else rows[0]
    if metric == "perplexity":
        values = np.exp(values)
    return _result(values=_plain(values), metric=metric, reverse_complement=reverse_complement)


__all__ = ["evo2_positional_uncertainty", "forward_evo2_sequences", "score_evo2_sequences"]
_register_public_verbose(globals(), tuple(__all__))
