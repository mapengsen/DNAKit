"""Lazy training-encoder dispatch shared by training and saved-model inference."""

from __future__ import annotations

from pathlib import Path
from typing import TypedDict

from dnakit.exceptions import ConfigurationError
from dnakit.representations import DNAEmbeddingModel

from .encoder import TrainableSequenceEncoder, TransformersSequenceEncoder

TRAINABLE_LOADERS = frozenset({"transformers", "enformer", "janusdna", "evo2", "alphagenome"})


class _EncoderOptions(TypedDict):
    device: str
    dtype: str
    pooling: str
    max_length: int | None
    allow_remote_code: bool


def create_trainable_encoder(
    spec: DNAEmbeddingModel,
    checkpoint: str | Path,
    *,
    device: str,
    dtype: str,
    pooling: str,
    max_length: int | None,
    allow_remote_code: bool,
    model_source_path: str | Path | None = None,
    seed: int = 42,
) -> TrainableSequenceEncoder:
    """Create or reload a differentiable encoder without importing unused runtimes."""
    options: _EncoderOptions = {
        "device": device,
        "dtype": dtype,
        "pooling": pooling,
        "max_length": max_length,
        "allow_remote_code": allow_remote_code,
    }
    if spec.loader == "transformers":
        return TransformersSequenceEncoder(spec, checkpoint, **options)
    if spec.loader == "enformer":
        from .encoder_enformer import EnformerSequenceEncoder

        return EnformerSequenceEncoder(spec, checkpoint, **options)
    if spec.loader in {"janusdna", "evo2"}:
        from .encoder_native import Evo2SequenceEncoder, JanusDNASequenceEncoder

        encoder = JanusDNASequenceEncoder if spec.loader == "janusdna" else Evo2SequenceEncoder
        return encoder(spec, checkpoint, model_source_path=model_source_path, **options)
    if spec.loader == "alphagenome":
        from .encoder_alphagenome import AlphaGenomeSequenceEncoder

        return AlphaGenomeSequenceEncoder(
            spec, checkpoint, model_source_path=model_source_path, seed=seed, **options
        )
    raise ConfigurationError(
        "The selected model has no full fine-tuning encoder.",
        code="CUSTOM_FULL_FINETUNE_UNSUPPORTED",
        context={"model": spec.name, "loader": spec.loader},
    )
