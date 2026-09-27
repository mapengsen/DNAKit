"""Differentiable Enformer sequence encoder and portable training checkpoints."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.representations import DNAEmbeddingModel
from dnakit.representations.backends import (
    _pool_hidden_states,
    _require_module,
    _resolve_torch_device,
    _torch_dtype,
)


class EnformerSequenceEncoder:
    """Train the complete trunk, replacing genomic-track heads with a task head.

    Sequence chunks use an uncropped trunk and length-aware pooling. The native
    gamma basis supports these shorter windows; the fixed TensorFlow lookup only
    supports the original 196,608 bp context. Saved artifacts retain this choice.
    """

    def __init__(
        self,
        spec: DNAEmbeddingModel,
        checkpoint: str | Path,
        *,
        device: str,
        dtype: str,
        pooling: str,
        max_length: int | None,
        allow_remote_code: bool,
        model_source_path: str | Path | None = None,
    ) -> None:
        del allow_remote_code, model_source_path
        self.spec = spec
        self.torch = _require_module("torch", spec)
        library = _require_module("enformer_pytorch", spec)
        self.device = _resolve_torch_device(self.torch, device)
        self.checkpoint = Path(checkpoint).expanduser().resolve()
        self.pooling = pooling
        metadata_path = self.checkpoint / "enformer_encoder.json"
        try:
            if metadata_path.is_file():
                metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
                if metadata.get("format") != "dnakit.enformer-encoder.v1":
                    raise ValueError("Unsupported Enformer encoder format.")
                from enformer_pytorch.config_enformer import EnformerConfig

                self.model = library.Enformer(EnformerConfig(**metadata["config"]))
                state = self.torch.load(
                    self.checkpoint / "encoder.pt", map_location="cpu", weights_only=True
                )
                self.model.load_state_dict(state, strict=True)
                self.max_length = int(metadata["max_length"])
                self.pooling = str(metadata["pooling"])
            else:
                self.model = library.from_pretrained(str(self.checkpoint))
                self.max_length = int(max_length or min(spec.chunk_length, 4_096))
                self.model.add_heads()
                self.model.config.output_heads = {}
                self.model.config.target_length = -1
                self.model.target_length = -1
                self.model.crop_final.target_length = -1
                self.model.config.use_tf_gamma = False
                for module in self.model.modules():
                    if hasattr(module, "use_tf_gamma"):
                        module.use_tf_gamma = False
        except Exception as exc:
            raise BackendExecutionError(
                "Could not load the Enformer encoder for full fine-tuning.",
                code="MODEL_LOAD_FAILED",
                context={"model": spec.name, "checkpoint": str(self.checkpoint)},
            ) from exc
        if not 0 < self.max_length <= spec.chunk_length:
            raise ConfigurationError(
                "Enformer max_length must be positive and within the registered context limit.",
                code="INVALID_CUSTOM_PROPERTY_CONFIG",
                context={"max_length": self.max_length, "limit": spec.chunk_length},
            )
        if self.device.type == "cpu" and dtype in {"float16", "bfloat16"}:
            raise ConfigurationError(
                "Use float32 or auto for Enformer full fine-tuning on CPU.",
                code="INVALID_CUSTOM_PROPERTY_CONFIG",
            )
        self.downsample = 2 ** int(self.model.config.num_downsamples)
        # At least two bins also make a one-sample training batch valid for BatchNorm.
        self.padded_length = max(
            2 * self.downsample,
            math.ceil(self.max_length / self.downsample) * self.downsample,
        )
        requested_dtype = _torch_dtype(self.torch, dtype)
        if requested_dtype is not None:
            self.model = self.model.to(dtype=requested_dtype)
        self.model = self.model.to(self.device)
        for parameter in self.model.parameters():
            parameter.requires_grad_(True)

    def parameters(self) -> Any:
        return self.model.parameters()

    def train(self) -> None:
        self.model.train()

    def eval(self) -> None:
        self.model.eval()

    def encode(self, sequences: list[str]) -> Any:
        torch = self.torch
        chunks: list[str] = []
        counts: list[int] = []
        for sequence in sequences:
            if not sequence:
                raise ConfigurationError(
                    "Enformer training requires non-empty sequences.",
                    code="INVALID_CUSTOM_ENCODER_INPUT",
                )
            selected = [
                sequence[start : start + self.max_length]
                for start in range(0, len(sequence), self.max_length)
            ]
            chunks.extend(selected)
            counts.append(len(selected))
        if not chunks:
            raise ConfigurationError(
                "Enformer training requires at least one sequence.",
                code="INVALID_CUSTOM_ENCODER_INPUT",
            )
        parameter = next(self.model.parameters())
        encoded = torch.zeros(
            (len(chunks), self.padded_length, 4),
            device=self.device,
            dtype=parameter.dtype,
        )
        for row, chunk in enumerate(chunks):
            for column, base in enumerate("ACGT"):
                encoded[row, : len(chunk), column] = torch.tensor(
                    [character == base for character in chunk],
                    device=self.device,
                    dtype=parameter.dtype,
                )
        hidden = self.model(encoded, return_only_embeddings=True)
        valid_bins = torch.tensor(
            [math.ceil(len(chunk) / self.downsample) for chunk in chunks],
            device=self.device,
        )
        mask = torch.arange(hidden.shape[1], device=self.device)[None, :] < valid_bins[:, None]
        pooled = _pool_hidden_states(hidden, torch=torch, pooling=self.pooling, attention_mask=mask)
        return torch.stack([part.mean(dim=0) for part in pooled.split(counts)])

    def save_pretrained(self, directory: Path) -> None:
        directory.mkdir(parents=True, exist_ok=True)
        metadata = {
            "format": "dnakit.enformer-encoder.v1",
            "config": self.model.config.to_dict(),
            "max_length": self.max_length,
            "pooling": self.pooling,
        }
        (directory / "enformer_encoder.json").write_text(
            json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
        )
        # Enformer registers trunk modules under multiple names; torch.save
        # preserves these shared tensors without safetensors alias loss.
        self.torch.save(self.model.state_dict(), directory / "encoder.pt")

    def load_pretrained(self, directory: Path) -> None:
        """Restore the selected epoch without allocating a second GPU backbone."""
        state = self.torch.load(directory / "encoder.pt", map_location="cpu", weights_only=True)
        self.model.load_state_dict(state, strict=True)
