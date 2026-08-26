"""Registry and configuration for DNA foundation-model representations."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Literal, TypeAlias

from dnakit.exceptions import ConfigurationError

CheckpointKind: TypeAlias = Literal["huggingface", "dataverse"]
LoaderKind: TypeAlias = Literal[
    "transformers",
    "enformer",
    "alphagenome",
    "janusdna",
    "evo2",
]
PoolingMethod: TypeAlias = Literal["mean", "cls", "max", "last"]
AmbiguityPolicy: TypeAlias = Literal["replace_with_n", "error"]
TorchDType: TypeAlias = Literal["auto", "float32", "float16", "bfloat16"]


@dataclass(frozen=True, slots=True)
class DNAEmbeddingModel:
    """One selectable DNA model and its official checkpoint metadata."""

    name: str
    display_name: str
    loader: LoaderKind
    checkpoint_kind: CheckpointKind
    checkpoint_id: str
    cache_name: str
    source_repository: str
    checkpoint_url: str
    chunk_length: int
    chunk_unit: Literal["bp", "token"]
    default_pooling: PoolingMethod = "mean"
    trust_remote_code: bool = False
    required_package: str | None = None
    output_dimension: int | None = None
    allow_patterns: tuple[str, ...] | None = None
    archive_name: str | None = None
    archive_md5: str | None = None
    evo2_layer: str | None = None
    gated: bool = False
    notes: str = ""


_HF_MODEL_PATTERNS = (
    "*.bin",
    "*.json",
    "*.model",
    "*.pt",
    "*.py",
    "*.safetensors",
    "*.txt",
)


MODEL_REGISTRY: dict[str, DNAEmbeddingModel] = {
    "dnabert2": DNAEmbeddingModel(
        name="dnabert2",
        display_name="DNABERT-2-117M",
        loader="transformers",
        checkpoint_kind="huggingface",
        checkpoint_id="zhihan1996/DNABERT-2-117M",
        cache_name="dnabert2-117m",
        source_repository="https://github.com/MAGICS-LAB/DNABERT_2",
        checkpoint_url="https://huggingface.co/zhihan1996/DNABERT-2-117M",
        chunk_length=512,
        chunk_unit="token",
        trust_remote_code=True,
        allow_patterns=_HF_MODEL_PATTERNS,
    ),
    "ntv2": DNAEmbeddingModel(
        name="ntv2",
        display_name="Nucleotide Transformer v2 500M multi-species",
        loader="transformers",
        checkpoint_kind="huggingface",
        checkpoint_id="InstaDeepAI/nucleotide-transformer-v2-500m-multi-species",
        cache_name="ntv2-500m-multi-species",
        source_repository="https://github.com/instadeepai/nucleotide-transformer",
        checkpoint_url=(
            "https://huggingface.co/InstaDeepAI/nucleotide-transformer-v2-500m-multi-species"
        ),
        chunk_length=2_048,
        chunk_unit="token",
        trust_remote_code=True,
        allow_patterns=_HF_MODEL_PATTERNS,
    ),
    "hyenadna": DNAEmbeddingModel(
        name="hyenadna",
        display_name="HyenaDNA medium 450k",
        loader="transformers",
        checkpoint_kind="huggingface",
        checkpoint_id="LongSafari/hyenadna-medium-450k-seqlen-hf",
        cache_name="hyenadna-medium-450k",
        source_repository="https://github.com/HazyResearch/hyena-dna",
        checkpoint_url="https://huggingface.co/LongSafari/hyenadna-medium-450k-seqlen-hf",
        chunk_length=450_000,
        chunk_unit="bp",
        trust_remote_code=True,
        allow_patterns=_HF_MODEL_PATTERNS,
    ),
    "caduceus": DNAEmbeddingModel(
        name="caduceus",
        display_name="Caduceus-Ph 131k",
        loader="transformers",
        checkpoint_kind="huggingface",
        checkpoint_id="kuleshov-group/caduceus-ph_seqlen-131k_d_model-256_n_layer-16",
        cache_name="caduceus-ph-131k",
        source_repository="https://github.com/kuleshov-group/caduceus",
        checkpoint_url=(
            "https://huggingface.co/kuleshov-group/caduceus-ph_seqlen-131k_d_model-256_n_layer-16"
        ),
        chunk_length=131_072,
        chunk_unit="bp",
        trust_remote_code=True,
        allow_patterns=_HF_MODEL_PATTERNS,
    ),
    "grover": DNAEmbeddingModel(
        name="grover",
        display_name="GROVER",
        loader="transformers",
        checkpoint_kind="huggingface",
        checkpoint_id="PoetschLab/GROVER",
        cache_name="grover",
        source_repository="https://huggingface.co/PoetschLab/GROVER",
        checkpoint_url="https://huggingface.co/PoetschLab/GROVER",
        chunk_length=512,
        chunk_unit="token",
        allow_patterns=_HF_MODEL_PATTERNS,
        notes="Only model and tokenizer files are downloaded; tokenized chromosomes are excluded.",
    ),
    "lucaone": DNAEmbeddingModel(
        name="lucaone",
        display_name="LucaOne gene step 36.8M",
        loader="transformers",
        checkpoint_kind="huggingface",
        checkpoint_id="LucaGroup/LucaOne-gene-step36.8M",
        cache_name="lucaone-gene-step36-8m",
        source_repository="https://github.com/LucaOne/LucaOne",
        checkpoint_url="https://huggingface.co/LucaGroup/LucaOne-gene-step36.8M",
        chunk_length=4_096,
        chunk_unit="token",
        trust_remote_code=True,
        allow_patterns=_HF_MODEL_PATTERNS,
    ),
    "generator": DNAEmbeddingModel(
        name="generator",
        display_name="GENERator v2 eukaryote 1.2B",
        loader="transformers",
        checkpoint_kind="huggingface",
        checkpoint_id="GenerTeam/GENERator-v2-eukaryote-1.2b-base",
        cache_name="generator-v2-eukaryote-1-2b",
        source_repository="https://github.com/GenerTeam/GENERator",
        checkpoint_url="https://huggingface.co/GenerTeam/GENERator-v2-eukaryote-1.2b-base",
        chunk_length=16_384,
        chunk_unit="token",
        trust_remote_code=True,
        allow_patterns=_HF_MODEL_PATTERNS,
    ),
    "enformer": DNAEmbeddingModel(
        name="enformer",
        display_name="Enformer PyTorch",
        loader="enformer",
        checkpoint_kind="huggingface",
        checkpoint_id="EleutherAI/enformer-official-rough",
        cache_name="enformer-official-rough",
        source_repository="https://github.com/lucidrains/enformer-pytorch",
        checkpoint_url="https://huggingface.co/EleutherAI/enformer-official-rough",
        chunk_length=196_608,
        chunk_unit="bp",
        default_pooling="mean",
        required_package="enformer-pytorch>=0.8.11",
        output_dimension=3_072,
        allow_patterns=_HF_MODEL_PATTERNS,
    ),
    "alphagenome": DNAEmbeddingModel(
        name="alphagenome",
        display_name="AlphaGenome all folds",
        loader="alphagenome",
        checkpoint_kind="huggingface",
        checkpoint_id="google/alphagenome-all-folds",
        cache_name="alphagenome-all-folds",
        source_repository="https://github.com/google-deepmind/alphagenome_research",
        checkpoint_url="https://huggingface.co/google/alphagenome-all-folds",
        chunk_length=1_048_576,
        chunk_unit="bp",
        default_pooling="mean",
        required_package=("git+https://github.com/google-deepmind/alphagenome_research.git"),
        output_dimension=1_536,
        gated=True,
        allow_patterns=None,
        notes="Checkpoint access requires accepting Google's non-commercial model terms.",
    ),
    "janusdna": DNAEmbeddingModel(
        name="janusdna",
        display_name="JanusDNA 131k",
        loader="janusdna",
        checkpoint_kind="dataverse",
        checkpoint_id="doi:10.7910/DVN/HDT0RN",
        cache_name="janusdna-131k",
        source_repository="https://github.com/Qihao-Duan/JanusDNA",
        checkpoint_url="https://dataverse.harvard.edu/api/access/datafile/12016924",
        chunk_length=131_072,
        chunk_unit="bp",
        trust_remote_code=False,
        required_package="JanusDNA source checkout plus its official environment",
        output_dimension=144,
        archive_name=(
            "janusdna_len-131k_d_model-144_inter_dim-576_n_layer-8_lr-8e-3_"
            "step-50K_moeloss-true_1head_onlymoe_finalmlp.tar.bz2"
        ),
        archive_md5="1fce7f29c728f312d72262d6820e2ba0",
    ),
    "evo2": DNAEmbeddingModel(
        name="evo2",
        display_name="Evo 2 7B",
        loader="evo2",
        checkpoint_kind="huggingface",
        checkpoint_id="arcinstitute/evo2_7b",
        cache_name="evo2-7b",
        source_repository="https://github.com/ArcInstitute/evo2",
        checkpoint_url="https://huggingface.co/arcinstitute/evo2_7b",
        chunk_length=70_000,
        chunk_unit="bp",
        default_pooling="mean",
        required_package="evo2>=0.5.3,<0.6",
        output_dimension=4_096,
        allow_patterns=_HF_MODEL_PATTERNS,
        evo2_layer="blocks.28.mlp.l3",
    ),
}


_MODEL_ALIASES = {
    "alpha-genome": "alphagenome",
    "caduceus-ph": "caduceus",
    "dna-bert-2": "dnabert2",
    "dnabert-2": "dnabert2",
    "enformer-pytorch": "enformer",
    "evo-2": "evo2",
    "generator-v2": "generator",
    "hyena-dna": "hyenadna",
    "janus-dna": "janusdna",
    "luca-one": "lucaone",
    "nt-v2": "ntv2",
    "nucleotide-transformer": "ntv2",
}


def available_embedding_models() -> tuple[str, ...]:
    """Return the stable model names accepted by :class:`RepresentationConfig`."""

    return tuple(sorted(MODEL_REGISTRY))


def get_embedding_model(name: str) -> DNAEmbeddingModel:
    """Resolve one canonical model name or documented alias."""

    if not isinstance(name, str) or not name.strip():
        raise ConfigurationError(
            "model must be a non-empty string.",
            code="INVALID_EMBEDDING_MODEL",
        )
    key = name.strip().lower().replace("_", "-")
    canonical = _MODEL_ALIASES.get(key, key.replace("-", ""))
    if canonical not in MODEL_REGISTRY:
        raise ConfigurationError(
            "Unknown DNA embedding model.",
            code="INVALID_EMBEDDING_MODEL",
            context={"model": name, "available": available_embedding_models()},
        )
    return MODEL_REGISTRY[canonical]


@dataclass(frozen=True, slots=True)
class RepresentationConfig:
    """Configure checkpoint resolution and sequence-level representation extraction."""

    model: str = "lucaone"
    checkpoint_dir: str | os.PathLike[str] | None = None
    checkpoint_path: str | os.PathLike[str] | None = None
    model_source_path: str | os.PathLike[str] | None = None
    pooling: PoolingMethod = "mean"
    ambiguity_policy: AmbiguityPolicy = "replace_with_n"
    device: str = "auto"
    dtype: TorchDType = "auto"
    batch_size: int = 4
    max_length: int | None = None
    max_records: int = 10_000
    show_progress: bool = True
    allow_remote_code: bool = False
    hf_token: str | None = field(default=None, repr=False)

    def __post_init__(self) -> None:
        model_spec = get_embedding_model(self.model)
        object.__setattr__(self, "model", model_spec.name)
        for field_name in ("checkpoint_dir", "checkpoint_path", "model_source_path"):
            value = getattr(self, field_name)
            if value is None:
                continue
            try:
                resolved = os.fspath(value)
            except TypeError as exc:
                raise ConfigurationError(
                    f"{field_name} must be path-like or None.",
                    code="INVALID_REPRESENTATION_CONFIG",
                ) from exc
            if not resolved.strip():
                raise ConfigurationError(
                    f"{field_name} must not be empty.",
                    code="INVALID_REPRESENTATION_CONFIG",
                )
            object.__setattr__(self, field_name, resolved)
        if self.checkpoint_dir is not None and self.checkpoint_path is not None:
            raise ConfigurationError(
                "checkpoint_dir and checkpoint_path are mutually exclusive.",
                code="INVALID_REPRESENTATION_CONFIG",
            )
        if self.pooling not in {"mean", "cls", "max", "last"}:
            raise ConfigurationError(
                "Unknown representation pooling method.",
                code="INVALID_REPRESENTATION_POOLING",
            )
        if self.ambiguity_policy not in {"replace_with_n", "error"}:
            raise ConfigurationError(
                "Unknown ambiguity_policy.",
                code="INVALID_REPRESENTATION_CONFIG",
            )
        if not isinstance(self.device, str) or not self.device.strip():
            raise ConfigurationError(
                "device must be non-empty text.",
                code="INVALID_REPRESENTATION_CONFIG",
            )
        if self.dtype not in {"auto", "float32", "float16", "bfloat16"}:
            raise ConfigurationError(
                "Unknown representation dtype.",
                code="INVALID_REPRESENTATION_CONFIG",
            )
        for name, value, maximum in (
            ("batch_size", self.batch_size, 1_024),
            ("max_records", self.max_records, 1_000_000),
        ):
            if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= maximum:
                raise ConfigurationError(
                    f"{name} must be in [1, {maximum}].",
                    code="INVALID_REPRESENTATION_CONFIG",
                )
        if self.max_length is not None and (
            isinstance(self.max_length, bool)
            or not isinstance(self.max_length, int)
            or not 1 <= self.max_length <= 10_000_000
        ):
            raise ConfigurationError(
                "max_length must be None or an integer in [1, 10000000].",
                code="INVALID_REPRESENTATION_CONFIG",
            )
        if not isinstance(self.show_progress, bool):
            raise ConfigurationError(
                "show_progress must be boolean.",
                code="INVALID_REPRESENTATION_CONFIG",
            )
        if not isinstance(self.allow_remote_code, bool):
            raise ConfigurationError(
                "allow_remote_code must be boolean.",
                code="INVALID_REPRESENTATION_CONFIG",
            )
        if self.hf_token is not None and (
            not isinstance(self.hf_token, str) or not self.hf_token.strip()
        ):
            raise ConfigurationError(
                "hf_token must be None or non-empty text.",
                code="INVALID_REPRESENTATION_CONFIG",
            )

    @property
    def model_spec(self) -> DNAEmbeddingModel:
        """Return the resolved registry entry."""

        return get_embedding_model(self.model)


__all__ = [
    "MODEL_REGISTRY",
    "AmbiguityPolicy",
    "DNAEmbeddingModel",
    "PoolingMethod",
    "RepresentationConfig",
    "TorchDType",
    "available_embedding_models",
    "get_embedding_model",
]
