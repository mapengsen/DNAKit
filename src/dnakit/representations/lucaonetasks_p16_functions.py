"""DNABERT-S / DNABERT-2 sequence embeddings with the LucaOneTasks protocol.

LucaOneTasks (``src/llm/dnaberts/inference_embedding.py`` and ``src/llm/dnabert2/inference_embedding.py``) embeds one
DNA sequence at a time: the sequence is upper-cased and reduced to the letters A-Z, U is read as T, it is truncated to
4,094 nucleotides (right side by default) and encoded by the checkpoint's tokenizer; the model's last hidden state is
returned either as a matrix (one row per BPE token, with or without the ``[CLS]`` and ``[SEP]`` rows) or as the
``[CLS]`` vector. :func:`dnakit.representations.extract_representations` already pools DNABERT-2 into one vector with the
special tokens excluded; this function adds the per-token matrix and ``[CLS]`` outputs and the DNABERT-S checkpoint.
"""

from __future__ import annotations

import sys
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

_MODELS = {"dnabert-s": "zhihan1996/DNABERT-S", "dnabert2": "zhihan1996/DNABERT-2-117M"}
_EMBEDDING_TYPES = ("matrix", "vector", "mean")
_MAX_SEQUENCES = 10_000
_MAX_INLINE_VALUES = 200_000


def _clean(sequence: str) -> str:
    """LucaOneTasks ``clean_seq_luca`` after ``U -> T``: upper case, keep the letters A-Z."""
    return "".join(ch for ch in sequence.upper().replace("U", "T") if "A" <= ch <= "Z")


def _records(sequences: Mapping[str, str] | Sequence[str]) -> list[tuple[str, str]]:
    if isinstance(sequences, str):
        raise ConfigurationError(
            "sequences must be a list of DNA strings or an id-to-sequence mapping."
        )
    items = (
        list(sequences.items())
        if isinstance(sequences, Mapping)
        else [(f"sequence_{i + 1}", s) for i, s in enumerate(sequences)]
    )
    if not items or len(items) > _MAX_SEQUENCES:
        raise ConfigurationError(f"Provide between 1 and {_MAX_SEQUENCES} sequences.")
    seen: set[str] = set()
    result: list[tuple[str, str]] = []
    for identifier, sequence in items:
        if not isinstance(identifier, str) or not identifier or identifier in seen:
            raise ConfigurationError("Sequence ids must be unique non-empty strings.")
        seen.add(identifier)
        if not isinstance(sequence, str) or not _clean(sequence):
            raise ConfigurationError(f"Sequence {identifier!r} has no usable letters.")
        result.append((identifier, _clean(sequence)))
    return result


def _disable_flash_attention(model: Any, checkpoint: Path) -> None:
    """The checkpoint's Triton flash-attention path needs CUDA; use the PyTorch attention instead."""
    names = {type(layer).__module__ for layer in model.modules()}
    for module in tuple(sys.modules.values()):
        module_file = getattr(module, "__file__", "") or ""
        in_checkpoint = getattr(module, "__name__", "") in names or checkpoint.name in module_file
        if in_checkpoint and hasattr(module, "flash_attn_qkvpacked_func"):
            module.__dict__["flash_attn_qkvpacked_func"] = None


def embed_dna_with_dnabert(
    sequences: Mapping[str, str] | Sequence[str],
    *,
    model: str = "dnabert-s",
    checkpoint_path: str | None = None,
    embedding_type: str = "matrix",
    include_special_tokens: bool = False,
    truncation_length: int = 4094,
    truncate_side: str = "right",
    output_path: str | None = None,
    device: str = "auto",
    allow_remote_code: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Embed DNA sequences with DNABERT-S or DNABERT-2 (LucaOneTasks ``dnaberts`` / ``dnabert2`` inference scripts).

    ``model`` is ``"dnabert-s"`` (checkpoint ``zhihan1996/DNABERT-S``) or ``"dnabert2"`` (``zhihan1996/DNABERT-2-117M``);
    ``checkpoint_path`` is a local checkpoint directory (a Hugging Face cache snapshot or a download); without it the
    model is looked up in the local Hugging Face cache only, nothing is downloaded. ``embedding_type`` is ``"matrix"``
    (token matrix, ``(tokens, 768)``; ``include_special_tokens`` keeps the ``[CLS]`` and ``[SEP]`` rows), ``"vector"``
    (the ``[CLS]`` vector of the last layer) or ``"mean"`` (mean of every token row including ``[CLS]`` and ``[SEP]``,
    the DNABERT-S paper protocol). ``truncation_length`` nucleotides are kept from the ``truncate_side`` (right keeps
    the start). The checkpoint ships Python code, so ``allow_remote_code=True`` is required. With ``output_path`` (a new
    ``.npz`` file) the arrays are saved under the sequence ids; ``matrix`` results are otherwise only summarised.
    Needs torch and transformers (``dnakit[neural]``).
    """
    records = _records(sequences)
    key = model.strip().lower().replace("_", "-") if isinstance(model, str) else ""
    key = "dnabert2" if key in {"dnabert-2", "dnabert2"} else key
    if key not in _MODELS:
        raise ConfigurationError(f"model must be one of {sorted(_MODELS)}; got {model!r}.")
    if embedding_type not in _EMBEDDING_TYPES:
        raise ConfigurationError(f"embedding_type must be one of {list(_EMBEDDING_TYPES)}.")
    if truncate_side not in {"right", "left"}:
        raise ConfigurationError("truncate_side must be 'right' or 'left'.")
    if (
        isinstance(truncation_length, bool)
        or not isinstance(truncation_length, int)
        or truncation_length < 1
    ):
        raise ConfigurationError("truncation_length must be a positive integer.")
    if not isinstance(include_special_tokens, bool) or not isinstance(allow_remote_code, bool):
        raise ConfigurationError("include_special_tokens and allow_remote_code must be booleans.")
    if not allow_remote_code:
        raise ConfigurationError(
            "The DNABERT checkpoints load Python code that ships with the checkpoint.",
            hint="Review the repository zhihan1096/"
            + _MODELS[key].split("/")[1]
            + ", then pass allow_remote_code=True.",
        )
    target = None
    if output_path is not None:
        target = Path(output_path).expanduser().resolve()
        if target.exists() or not target.parent.is_dir() or target.suffix != ".npz":
            raise ConfigurationError("output_path must be a new .npz file in an existing folder.")
    np = require("numpy", "bio-analysis")
    torch = require("torch", "neural")
    transformers = require("transformers", "neural")

    if checkpoint_path is not None:
        checkpoint = Path(checkpoint_path).expanduser().resolve()
        if not checkpoint.is_dir():
            raise ConfigurationError(f"checkpoint_path must be an existing directory: {checkpoint}")
        source: str = str(checkpoint)
        local_only = True
    else:
        source = _MODELS[key]
        local_only = True
        checkpoint = Path(source)
    try:
        tokenizer = transformers.AutoTokenizer.from_pretrained(
            source, trust_remote_code=True, local_files_only=local_only
        )
        network = transformers.AutoModel.from_pretrained(
            source, trust_remote_code=True, local_files_only=local_only
        )
    except OSError as exc:
        raise BackendUnavailableError(
            f"The {key} checkpoint {source!r} is not available locally.",
            hint=f"Download {_MODELS[key]} (for example with huggingface_hub.snapshot_download) and pass checkpoint_path.",
        ) from exc
    resolved = torch.device(
        ("cuda" if torch.cuda.is_available() else "cpu") if device == "auto" else device
    )
    network = network.to(resolved).eval()
    _disable_flash_attention(network, checkpoint)

    arrays: dict[str, Any] = {}
    summaries: list[dict[str, Any]] = []
    for identifier, sequence in records:
        kept = (
            sequence[-truncation_length:]
            if truncate_side == "left"
            else sequence[:truncation_length]
        )
        token_ids = tokenizer(kept, return_tensors="pt")["input_ids"].to(resolved)
        try:
            with torch.no_grad():
                hidden = network(token_ids)[0].to("cpu")
        except RuntimeError as exc:
            raise BackendExecutionError(
                f"{key} failed on {identifier!r}: {exc}", context={"sequence_length": len(kept)}
            ) from exc
        tokens = int(token_ids.shape[1])
        usable = min(truncation_length, tokens - 2)
        if embedding_type == "matrix":
            rows = (
                hidden[0, 0 : usable + 2] if include_special_tokens else hidden[0, 1 : usable + 1]
            )
            value = rows.clone().numpy()
        elif embedding_type == "vector":
            value = hidden[0, 0].clone().numpy()
        else:
            value = hidden[0].mean(dim=0).numpy()
        arrays[identifier] = value.astype(np.float32)
        summaries.append(
            {
                "id": identifier,
                "nucleotides_used": len(kept),
                "tokens": tokens,
                "shape": list(value.shape),
            }
        )
    if target is not None:
        np.savez_compressed(target, **arrays)
    inline = None
    if embedding_type != "matrix" and sum(a.size for a in arrays.values()) <= _MAX_INLINE_VALUES:
        inline = {name: [float(x) for x in value] for name, value in arrays.items()}
    return {
        "status": "completed",
        "backend": "transformers",
        "operation": "embed_dna_with_dnabert",
        "model": key,
        "checkpoint": source,
        "embedding_type": embedding_type,
        "include_special_tokens": include_special_tokens,
        "embedding_dimension": int(next(iter(arrays.values())).shape[-1]),
        "device": str(resolved),
        "sequences": summaries,
        "embeddings": inline,
        "outputs": [str(target)] if target is not None else [],
        "provenance": {
            "tool": "LucaOneTasks src/llm/"
            + ("dnaberts" if key == "dnabert-s" else "dnabert2")
            + "/inference_embedding.py",
            "transformers": getattr(transformers, "__version__", None),
            "torch": getattr(torch, "__version__", None),
            "citation": "https://doi.org/10.48550/arXiv.2402.08777",
        },
    }


__all__ = ["embed_dna_with_dnabert"]

_register_public_verbose(globals(), tuple(__all__))
