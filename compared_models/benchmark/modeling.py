from __future__ import annotations

import importlib.util
import json
import math
import sys
import types
from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import torch
import torch.nn as nn
import torch.nn.functional as F
from transformers import AutoConfig, AutoModel, AutoModelForMaskedLM, AutoTokenizer


PROJECT_ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class ModelSpec:
    name: str
    display_name: str
    local_checkpoint: str | None
    hf_model_id: str | None
    loader_kind: str
    chunk_length: int
    chunk_unit: str
    default_token_readout: str
    trust_remote_code: bool = False
    tokenizer_loader: str | None = None
    output_dim: int | None = None
    unsupported_reason: str | None = None
    allowed_finetune_methods: tuple[str, ...] = ("frozen_linear_probe",)
    ia3_target_modules: tuple[str, ...] = ()
    ia3_feedforward_modules: tuple[str, ...] = ()


MODEL_SPECS: dict[str, ModelSpec] = {
    "dnabert2": ModelSpec(
        name="dnabert2",
        display_name="DNABERT-2-117M",
        local_checkpoint="compared_models/models/DNABERT_2/weights/DNABERT-2-117M",
        hf_model_id="zhihan1996/DNABERT-2-117M",
        loader_kind="hf_token",
        chunk_length=512,
        chunk_unit="token",
        default_token_readout="mean",
        trust_remote_code=True,
        allowed_finetune_methods=("frozen_linear_probe", "ia3"),
        ia3_target_modules=("key", "value", "output.dense"),
        ia3_feedforward_modules=("output.dense",),
    ),
    "ntv2": ModelSpec(
        name="ntv2",
        display_name="NT-v2-500M",
        local_checkpoint="compared_models/models/nucleotide-transformer/weights/nucleotide-transformer-v2-500m-multi-species",
        hf_model_id="InstaDeepAI/nucleotide-transformer-v2-500m-multi-species",
        loader_kind="hf_token",
        chunk_length=2048,
        chunk_unit="token",
        default_token_readout="mean",
        trust_remote_code=True,
        allowed_finetune_methods=("frozen_linear_probe", "ia3"),
        ia3_target_modules=("key", "value", "output.dense"),
        ia3_feedforward_modules=("output.dense",),
    ),
    "hyenadna": ModelSpec(
        name="hyenadna",
        display_name="HyenaDNA-medium-450k",
        local_checkpoint="compared_models/models/hyena-dna/weights/hyenadna-medium-450k-seqlen-hf",
        hf_model_id="LongSafari/hyenadna-medium-450k-seqlen-hf",
        loader_kind="hf_token",
        chunk_length=450000,
        chunk_unit="bp",
        default_token_readout="mean",
        trust_remote_code=True,
        allowed_finetune_methods=("frozen_linear_probe", "full"),
    ),
    "caduceus": ModelSpec(
        name="caduceus",
        display_name="Caduceus-Ph-131k",
        local_checkpoint="compared_models/models/caduceus/weights/caduceus-ph_seqlen-131k_d_model-256_n_layer-16",
        hf_model_id="kuleshov-group/caduceus-ph_seqlen-131k_d_model-256_n_layer-16",
        loader_kind="hf_token",
        chunk_length=131072,
        chunk_unit="bp",
        default_token_readout="mean",
        trust_remote_code=True,
        tokenizer_loader="caduceus",
        allowed_finetune_methods=("frozen_linear_probe", "full"),
    ),
    "grover": ModelSpec(
        name="grover",
        display_name="GROVER",
        local_checkpoint="compared_models/models/GROVER/model",
        hf_model_id="PoetschLab/GROVER",
        loader_kind="hf_token",
        chunk_length=512,
        chunk_unit="token",
        default_token_readout="mean",
        allowed_finetune_methods=("frozen_linear_probe", "ia3"),
        ia3_target_modules=("key", "value", "output.dense"),
        ia3_feedforward_modules=("output.dense",),
    ),
    "lucaone": ModelSpec(
        name="lucaone",
        display_name="LucaOne-gene-step36.8M",
        local_checkpoint="compared_models/models/LucaOne/weights/LucaOne-gene-step36.8M",
        hf_model_id="LucaGroup/LucaOne-gene-step36.8M",
        loader_kind="hf_token",
        chunk_length=4096,
        chunk_unit="token",
        default_token_readout="mean",
        trust_remote_code=True,
        allowed_finetune_methods=("frozen_linear_probe", "ia3"),
        ia3_target_modules=("k_proj", "v_proj", "fc2"),
        ia3_feedforward_modules=("fc2",),
    ),
    "generator": ModelSpec(
        name="generator",
        display_name="GENERator-v2-eukaryote-1.2b-base",
        local_checkpoint="compared_models/models/GENERator/weights/GENERator-v2-eukaryote-1.2b-base",
        hf_model_id="GenerTeam/GENERator-v2-eukaryote-1.2b-base",
        loader_kind="hf_token",
        chunk_length=16384,
        chunk_unit="token",
        default_token_readout="mean",
        trust_remote_code=True,
        allowed_finetune_methods=("frozen_linear_probe", "ia3"),
        ia3_target_modules=("k_proj", "v_proj", "down_proj"),
        ia3_feedforward_modules=("down_proj",),
    ),
    "enformer": ModelSpec(
        name="enformer",
        display_name="Enformer-PyTorch",
        local_checkpoint="compared_models/models/enformer-pytorch/weights/enformer-official-rough",
        hf_model_id="EleutherAI/enformer-official-rough",
        loader_kind="enformer",
        chunk_length=196608,
        chunk_unit="bp",
        default_token_readout="embedding",
        output_dim=3072,
        allowed_finetune_methods=("frozen_linear_probe", "ia3"),
        ia3_target_modules=("to_k", "to_v", "3"),
        ia3_feedforward_modules=("3",),
    ),
    "alphagenome": ModelSpec(
        name="alphagenome",
        display_name="AlphaGenome-all-folds",
        local_checkpoint="compared_models/models/alphagenome_research/weights/alphagenome-all-folds",
        hf_model_id="google/alphagenome-all-folds",
        loader_kind="alphagenome",
        chunk_length=1048576,
        chunk_unit="bp",
        default_token_readout="embedding",
        output_dim=1536,
        allowed_finetune_methods=("frozen_linear_probe",),
    ),
    "janusdna": ModelSpec(
        name="janusdna",
        display_name="JanusDNA-131k",
        local_checkpoint="compared_models/models/JanusDNA/weights/janusdna_len-131k_d_model-144_inter_dim-576_n_layer-8_lr-8e-3_step-50K_moeloss-true_1head_onlymoe_finalmlp",
        hf_model_id=None,
        loader_kind="janusdna",
        chunk_length=1024,
        chunk_unit="bp",
        default_token_readout="mean",
        output_dim=144,
        allowed_finetune_methods=("frozen_linear_probe", "full"),
    ),
    "evo2": ModelSpec(
        name="evo2",
        display_name="Evo2-7B",
        local_checkpoint="compared_models/models/evo2/weights/evo2_7b/evo2_7b.pt",
        hf_model_id="arcinstitute/evo2_7b",
        loader_kind="evo2",
        chunk_length=70000,
        chunk_unit="bp",
        default_token_readout="mean",
        output_dim=4096,
        allowed_finetune_methods=("frozen_linear_probe",),
    ),
}


FINETUNE_METHOD_ALIASES = {
    "frozen": "frozen_linear_probe",
    "linear_probe": "frozen_linear_probe",
    "frozen_linear_probe": "frozen_linear_probe",
    "frozen_embeddings_linear_probe": "frozen_linear_probe",
    "frozen_embeddings_probe": "frozen_linear_probe",
    "full": "full",
    "full_finetuning": "full",
    "full_fine_tuning": "full",
    "ia3": "ia3",
    "ia3_like": "ia3",
}


def available_compared_models(*, include_unsupported: bool = True) -> list[str]:
    if include_unsupported:
        return sorted(MODEL_SPECS)
    return sorted(name for name, spec in MODEL_SPECS.items() if spec.loader_kind != "unsupported")


def get_model_spec(model_name: str) -> ModelSpec:
    key = str(model_name).strip().lower()
    if key not in MODEL_SPECS:
        available = ", ".join(available_compared_models())
        raise KeyError(f"未知对比模型 {model_name!r}，可选项: {available}")
    return MODEL_SPECS[key]


def normalize_finetune_method(method: str | None) -> str:
    value = "frozen_linear_probe" if method is None else str(method).strip().lower()
    value = value.replace("-", "_").replace("+", "_").replace(" ", "_")
    value = value.replace("__", "_")
    if value not in FINETUNE_METHOD_ALIASES:
        available = ", ".join(sorted(set(FINETUNE_METHOD_ALIASES.values())))
        raise ValueError(f"未知微调方式 {method!r}，可选项: {available}")
    return FINETUNE_METHOD_ALIASES[value]


def validate_finetune_method(spec: ModelSpec, method: str | None) -> str:
    resolved = normalize_finetune_method(method)
    if resolved not in spec.allowed_finetune_methods:
        allowed = ", ".join(spec.allowed_finetune_methods)
        raise ValueError(
            f"{spec.display_name} 不支持微调方式 {resolved!r}；"
            f"该模型允许的方式: {allowed}"
        )
    return resolved


def _parse_module_list(value: str | list[str] | tuple[str, ...] | None) -> tuple[str, ...] | None:
    if value is None:
        return None
    if isinstance(value, str):
        parts = [item.strip() for item in value.split(",")]
    else:
        parts = [str(item).strip() for item in value]
    return tuple(item for item in parts if item)


def _module_from_file(module_name: str, file_path: Path):
    spec = importlib.util.spec_from_file_location(module_name, str(file_path))
    if spec is None or spec.loader is None:
        raise ImportError(f"无法从 {file_path} 构造模块 {module_name}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _prepend_repo_to_sys_path(repo_root: Path) -> None:
    for candidate in (repo_root, repo_root / "src"):
        if not candidate.exists():
            continue
        candidate_str = str(candidate.resolve())
        if candidate_str not in sys.path:
            sys.path.insert(0, candidate_str)


def _disable_dnabert2_flash_attention() -> None:
    for module in list(sys.modules.values()):
        module_file = getattr(module, "__file__", "") or ""
        if "DNABERT-2-117M" not in module_file or not hasattr(module, "flash_attn_qkvpacked_func"):
            continue
        try:
            setattr(module, "flash_attn_qkvpacked_func", None)
        except Exception:
            continue


def _patch_mamba_rms_norm_cpu_fallback() -> None:
    try:
        import mamba_ssm.ops.triton.layer_norm as layer_norm_module
    except Exception:
        return

    if getattr(layer_norm_module.RMSNorm, "_visualdna_cpu_fallback", False):
        return

    original_forward = layer_norm_module.RMSNorm.forward

    def _cpu_safe_forward(self, x, residual=None, prenorm=False, residual_in_fp32=False):
        if x.is_cuda:
            return original_forward(
                self,
                x,
                residual=residual,
                prenorm=prenorm,
                residual_in_fp32=residual_in_fp32,
            )

        if residual is not None:
            x = x + residual
        residual_out = x.to(torch.float32) if residual_in_fp32 else x
        x_norm = x.to(self.weight.dtype)
        variance = x_norm.pow(2).mean(dim=-1, keepdim=True)
        y = x_norm * torch.rsqrt(variance + self.eps)
        y = y * self.weight
        if prenorm:
            return y, residual_out
        return y

    layer_norm_module.RMSNorm.forward = _cpu_safe_forward
    layer_norm_module.RMSNorm._visualdna_cpu_fallback = True

    try:
        import mamba_ssm.modules.mamba_simple as mamba_simple_module
        import mamba_ssm.ops.selective_scan_interface as selective_scan_interface
    except Exception:
        return

    if getattr(mamba_simple_module, "_visualdna_cpu_conv_fallback", False):
        return

    mamba_simple_module.causal_conv1d_fn = None
    mamba_simple_module.causal_conv1d_update = None
    fast_selective_scan_fn = getattr(mamba_simple_module, "selective_scan_fn", None)
    ref_selective_scan_fn = selective_scan_interface.selective_scan_ref
    if fast_selective_scan_fn is None or fast_selective_scan_fn is ref_selective_scan_fn:
        mamba_simple_module.selective_scan_fn = ref_selective_scan_fn
    else:
        def _cuda_safe_selective_scan_fn(
            u,
            delta,
            A,
            B,
            C,
            D=None,
            z=None,
            delta_bias=None,
            delta_softplus=False,
            return_last_state=False,
        ):
            if isinstance(u, torch.Tensor) and u.is_cuda:
                return fast_selective_scan_fn(
                    u,
                    delta,
                    A,
                    B,
                    C,
                    D=D,
                    z=z,
                    delta_bias=delta_bias,
                    delta_softplus=delta_softplus,
                    return_last_state=return_last_state,
                )
            return ref_selective_scan_fn(
                u,
                delta,
                A,
                B,
                C,
                D=D,
                z=z,
                delta_bias=delta_bias,
                delta_softplus=delta_softplus,
                return_last_state=return_last_state,
            )

        mamba_simple_module.selective_scan_fn = _cuda_safe_selective_scan_fn
    mamba_simple_module._visualdna_cpu_conv_fallback = True


def _disable_mamba_fast_path(model: nn.Module) -> None:
    for module in model.modules():
        if hasattr(module, "use_fast_path"):
            try:
                module.use_fast_path = False
            except Exception:
                continue


def _load_caduceus_tokenizer(checkpoint_path: Path, model_max_length: int):
    module_name = f"_visualdna_caduceus_tokenizer_{checkpoint_path.name}"
    if module_name in sys.modules:
        module = sys.modules[module_name]
    else:
        module = _module_from_file(module_name, checkpoint_path / "tokenization_caduceus.py")
    tokenizer_class = getattr(module, "CaduceusTokenizer")
    return tokenizer_class(model_max_length=model_max_length, padding_side="right")


def _load_tokenizer(
    spec: ModelSpec,
    checkpoint: str,
    *,
    cache_dir: str | None,
) -> Any:
    checkpoint_path = Path(checkpoint)
    model_max_length = None
    if checkpoint_path.exists():
        try:
            config = AutoConfig.from_pretrained(
                str(checkpoint_path),
                trust_remote_code=spec.trust_remote_code,
                local_files_only=True,
            )
            model_max_length = int(
                getattr(
                    config,
                    "max_position_embeddings",
                    getattr(config, "max_seq_len", spec.chunk_length),
                )
            )
        except Exception:
            model_max_length = int(spec.chunk_length)

    if spec.tokenizer_loader == "caduceus":
        if not checkpoint_path.exists():
            raise FileNotFoundError(
                f"Caduceus 自定义 tokenizer 需要本地 checkpoint 目录，收到: {checkpoint}"
            )
        return _load_caduceus_tokenizer(checkpoint_path, model_max_length or spec.chunk_length)

    local_only = checkpoint_path.exists()
    tokenizer = AutoTokenizer.from_pretrained(
        checkpoint,
        trust_remote_code=spec.trust_remote_code,
        local_files_only=local_only,
        cache_dir=cache_dir,
    )
    try:
        tokenizer.padding_side = "right"
    except Exception:
        pass
    return tokenizer


def _unwrap_base_model(model: nn.Module) -> nn.Module:
    for attr_name in ("esm", "bert", "roberta", "transformer", "base_model", "model"):
        candidate = getattr(model, attr_name, None)
        if isinstance(candidate, nn.Module):
            return candidate
    return model


def _load_hf_backbone(
    spec: ModelSpec,
    *,
    checkpoint: str,
    cache_dir: str | None,
) -> tuple[Any, nn.Module]:
    tokenizer = _load_tokenizer(spec, checkpoint, cache_dir=cache_dir)
    local_only = Path(checkpoint).exists()
    errors: list[str] = []
    model_config = None

    if spec.name == "caduceus":
        try:
            model_config = AutoConfig.from_pretrained(
                checkpoint,
                trust_remote_code=spec.trust_remote_code,
                local_files_only=local_only,
                cache_dir=cache_dir,
            )
            model_config.fused_add_norm = False
        except Exception as exc:
            errors.append(f"AutoConfig: {type(exc).__name__}: {exc}")

    try:
        model = AutoModel.from_pretrained(
            checkpoint,
            trust_remote_code=spec.trust_remote_code,
            local_files_only=local_only,
            cache_dir=cache_dir,
            config=model_config,
        )
        if spec.name == "dnabert2":
            _disable_dnabert2_flash_attention()
        if spec.name == "caduceus":
            _patch_mamba_rms_norm_cpu_fallback()
            _disable_mamba_fast_path(model)
        return tokenizer, model
    except Exception as exc:
        errors.append(f"AutoModel: {type(exc).__name__}: {exc}")

    try:
        mlm_model = AutoModelForMaskedLM.from_pretrained(
            checkpoint,
            trust_remote_code=spec.trust_remote_code,
            local_files_only=local_only,
            cache_dir=cache_dir,
            config=model_config,
        )
        if spec.name == "dnabert2":
            _disable_dnabert2_flash_attention()
        if spec.name == "caduceus":
            _patch_mamba_rms_norm_cpu_fallback()
            _disable_mamba_fast_path(mlm_model)
        return tokenizer, _unwrap_base_model(mlm_model)
    except Exception as exc:
        errors.append(f"AutoModelForMaskedLM: {type(exc).__name__}: {exc}")

    raise RuntimeError("; ".join(errors))


def _infer_output_dim(model: nn.Module, spec: ModelSpec) -> int:
    if spec.output_dim is not None:
        return int(spec.output_dim)
    config = getattr(model, "config", None)
    if config is None and hasattr(model, "get_base_model"):
        try:
            config = getattr(model.get_base_model(), "config", None)
        except Exception:
            config = None
    if config is not None:
        for attr_name in ("hidden_size", "d_model", "n_embd", "dim", "embed_dim"):
            value = getattr(config, attr_name, None)
            if value is not None:
                return int(value)
    for attr_name in ("hidden_size", "d_model", "n_embd", "dim", "embed_dim"):
        value = getattr(model, attr_name, None)
        if value is not None:
            return int(value)
    raise ValueError(f"无法推断模型 {spec.name} 的输出维度。")


def _safe_model_max_length(tokenizer: Any, spec: ModelSpec) -> int:
    value = getattr(tokenizer, "model_max_length", None)
    if value is None or int(value) <= 0 or int(value) > 10_000_000:
        return int(spec.chunk_length)
    return int(value)


def _token_window_length(tokenizer: Any, spec: ModelSpec) -> int:
    special_count = 0
    if hasattr(tokenizer, "num_special_tokens_to_add"):
        try:
            special_count = int(tokenizer.num_special_tokens_to_add(pair=False))
        except Exception:
            special_count = 0
    max_length = _safe_model_max_length(tokenizer, spec)
    usable = max(1, int(max_length) - int(special_count))
    return int(min(spec.chunk_length, usable))


def _resolve_hidden_states(outputs: Any) -> torch.Tensor:
    if hasattr(outputs, "last_hidden_state") and outputs.last_hidden_state is not None:
        return outputs.last_hidden_state
    if hasattr(outputs, "hidden_states") and outputs.hidden_states:
        return outputs.hidden_states[-1]
    if isinstance(outputs, tuple) and outputs:
        candidate = outputs[0]
        if isinstance(candidate, torch.Tensor):
            return candidate
    raise ValueError("模型输出中未找到 last_hidden_state / hidden_states。")


def _build_mean_mask(
    *,
    attention_mask: torch.Tensor | None,
    special_tokens_mask: torch.Tensor | None,
    hidden_states: torch.Tensor,
) -> torch.Tensor:
    if attention_mask is None:
        mask = torch.ones(hidden_states.shape[:2], device=hidden_states.device, dtype=torch.bool)
    else:
        mask = attention_mask.to(hidden_states.device).bool()
    if special_tokens_mask is not None:
        candidate = mask & (~special_tokens_mask.to(hidden_states.device).bool())
        if candidate.any(dim=1).all():
            mask = candidate
    return mask


def _apply_token_readout(
    hidden_states: torch.Tensor,
    *,
    attention_mask: torch.Tensor | None,
    special_tokens_mask: torch.Tensor | None,
    token_readout: str,
) -> torch.Tensor:
    if token_readout == "mean":
        mask = _build_mean_mask(
            attention_mask=attention_mask,
            special_tokens_mask=special_tokens_mask,
            hidden_states=hidden_states,
        )
        weights = mask.to(hidden_states.dtype).unsqueeze(-1)
        denom = weights.sum(dim=1).clamp_min(1.0)
        return (hidden_states * weights).sum(dim=1) / denom

    if token_readout == "last":
        if attention_mask is None:
            indices = torch.full(
                (hidden_states.shape[0],),
                hidden_states.shape[1] - 1,
                device=hidden_states.device,
                dtype=torch.long,
            )
        else:
            positions = torch.arange(
                hidden_states.shape[1],
                device=hidden_states.device,
                dtype=torch.long,
            )[None, :]
            masked_positions = positions.masked_fill(
                ~attention_mask.to(hidden_states.device).bool(),
                -1,
            )
            indices = masked_positions.max(dim=1).values.clamp_min(0)
        gather_index = indices[:, None, None].expand(-1, 1, hidden_states.shape[-1])
        return hidden_states.gather(dim=1, index=gather_index).squeeze(1)

    raise ValueError(f"不支持的 token_readout: {token_readout}")


def _fill_special_tokens_mask(tokenizer: Any, batch: dict[str, torch.Tensor]) -> dict[str, torch.Tensor]:
    if "special_tokens_mask" in batch:
        return batch
    if "input_ids" not in batch:
        return batch
    special_ids = set(getattr(tokenizer, "all_special_ids", []) or [])
    if not special_ids:
        batch["special_tokens_mask"] = torch.zeros_like(batch["input_ids"], dtype=torch.long)
        return batch
    mask = torch.zeros_like(batch["input_ids"], dtype=torch.long)
    for special_id in special_ids:
        mask = mask | (batch["input_ids"] == int(special_id))
    batch["special_tokens_mask"] = mask.long()
    return batch


def _set_requires_grad(module: nn.Module, requires_grad: bool) -> None:
    for parameter in module.parameters():
        parameter.requires_grad = bool(requires_grad)


def _module_name_matches(module_name: str, patterns: tuple[str, ...]) -> bool:
    if not patterns:
        return False
    for pattern in patterns:
        if pattern == "all-linear":
            return True
        if module_name == pattern or module_name.endswith(f".{pattern}"):
            return True
    return False


def _set_child_module(parent: nn.Module, child_name: str, child: nn.Module) -> None:
    if isinstance(parent, (nn.Sequential, nn.ModuleList)) and child_name.isdigit():
        parent[int(child_name)] = child
    else:
        setattr(parent, child_name, child)


def _resolve_parent_module(root: nn.Module, module_name: str) -> tuple[nn.Module, str]:
    if "." not in module_name:
        return root, module_name
    parent_name, child_name = module_name.rsplit(".", 1)
    return root.get_submodule(parent_name), child_name


class IA3LinearLike(nn.Module):
    """Minimal IA3-style scaling wrapper for non-PEFT modules."""

    def __init__(self, linear: nn.Linear, *, feedforward: bool) -> None:
        super().__init__()
        self.linear = linear
        self.feedforward = bool(feedforward)
        scale_dim = linear.in_features if self.feedforward else linear.out_features
        self.ia3_l = nn.Parameter(torch.ones(int(scale_dim)))

    def _scale_view(self, value: torch.Tensor) -> torch.Tensor:
        shape = [1 for _ in range(value.ndim)]
        shape[-1] = -1
        return self.ia3_l.to(dtype=value.dtype, device=value.device).view(*shape)

    def forward(self, inputs: torch.Tensor) -> torch.Tensor:
        if self.feedforward:
            return self.linear(inputs * self._scale_view(inputs))
        outputs = self.linear(inputs)
        return outputs * self._scale_view(outputs)


def _apply_ia3_like_layers(
    model: nn.Module,
    *,
    target_modules: tuple[str, ...],
    feedforward_modules: tuple[str, ...],
) -> int:
    replacements: list[tuple[str, nn.Linear, bool]] = []
    for module_name, module in model.named_modules():
        if not module_name or not isinstance(module, nn.Linear):
            continue
        if not _module_name_matches(module_name, target_modules):
            continue
        replacements.append(
            (
                module_name,
                module,
                _module_name_matches(module_name, feedforward_modules),
            )
        )

    for module_name, linear, is_feedforward in replacements:
        parent, child_name = _resolve_parent_module(model, module_name)
        _set_child_module(
            parent,
            child_name,
            IA3LinearLike(linear, feedforward=is_feedforward),
        )
    return len(replacements)


def _apply_peft_ia3(
    model: nn.Module,
    *,
    target_modules: tuple[str, ...],
    feedforward_modules: tuple[str, ...],
) -> nn.Module:
    try:
        from peft import IA3Config, get_peft_model
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "使用 --finetune-method ia3 需要安装 peft；"
            "请在当前环境中运行: conda run -n visualdna pip install peft"
        ) from exc

    config = IA3Config(
        target_modules=list(target_modules),
        feedforward_modules=list(feedforward_modules),
    )
    return get_peft_model(model, config)


def _configure_backbone_finetuning(
    model: nn.Module,
    *,
    spec: ModelSpec,
    finetune_method: str,
    ia3_target_modules: tuple[str, ...] | None,
    ia3_feedforward_modules: tuple[str, ...] | None,
    prefer_peft: bool = True,
) -> nn.Module:
    if finetune_method == "frozen_linear_probe":
        _set_requires_grad(model, False)
        return model.eval()

    if finetune_method == "full":
        _set_requires_grad(model, True)
        return model.train()

    if finetune_method != "ia3":
        raise ValueError(f"不支持的微调方式: {finetune_method}")

    target_modules = ia3_target_modules or spec.ia3_target_modules
    feedforward_modules = ia3_feedforward_modules or spec.ia3_feedforward_modules
    if not target_modules:
        raise ValueError(f"{spec.display_name} 未配置 IA3 target_modules。")
    if not feedforward_modules:
        raise ValueError(f"{spec.display_name} 未配置 IA3 feedforward_modules。")

    _set_requires_grad(model, False)
    peft_error: Exception | None = None
    if prefer_peft:
        try:
            return _apply_peft_ia3(
                model,
                target_modules=target_modules,
                feedforward_modules=feedforward_modules,
            ).train()
        except Exception as exc:
            peft_error = exc

    injected = _apply_ia3_like_layers(
        model,
        target_modules=target_modules,
        feedforward_modules=feedforward_modules,
    )
    if injected == 0:
        suffix = f" PEFT 错误: {type(peft_error).__name__}: {peft_error}" if peft_error else ""
        raise RuntimeError(f"{spec.display_name} 未匹配到可注入 IA3-like 的 Linear 层。{suffix}")
    return model.train()


class BaseBackboneAdapter(nn.Module):
    def __init__(self, spec: ModelSpec, *, finetune_method: str = "frozen_linear_probe") -> None:
        super().__init__()
        self.spec = spec
        self.finetune_method = finetune_method
        self.output_dim = int(spec.output_dim or 0)

    @property
    def needs_backbone_grad(self) -> bool:
        return self.finetune_method != "frozen_linear_probe"

    def forward_context(self):
        return nullcontext() if self.needs_backbone_grad else torch.no_grad()

    def encode(self, sequences: list[str], *, token_readout: str) -> torch.Tensor:
        raise NotImplementedError

    def encode_serial(self, sequences: list[str], *, token_readout: str) -> torch.Tensor:
        raise NotImplementedError

    def encode_batched(self, sequences: list[str], *, token_readout: str) -> torch.Tensor:
        raise NotImplementedError


class HFTokenBackboneAdapter(BaseBackboneAdapter):
    def __init__(
        self,
        spec: ModelSpec,
        *,
        checkpoint: str,
        token_readout: str,
        cache_dir: str | None,
        allow_remote_fallback: bool,
        forward_batch_size: int,
        finetune_method: str,
        ia3_target_modules: tuple[str, ...] | None,
        ia3_feedforward_modules: tuple[str, ...] | None,
    ) -> None:
        self.token_readout = token_readout
        self.forward_batch_size = max(int(forward_batch_size), 1)
        self.cache_dir = cache_dir
        self.allow_remote_fallback = bool(allow_remote_fallback)
        self.checkpoint = checkpoint
        tokenizer, model = self._load_tokenizer_and_model(spec)
        super().__init__(spec, finetune_method=finetune_method)
        self.tokenizer = tokenizer
        self.backbone = _configure_backbone_finetuning(
            model,
            spec=spec,
            finetune_method=finetune_method,
            ia3_target_modules=ia3_target_modules,
            ia3_feedforward_modules=ia3_feedforward_modules,
            prefer_peft=True,
        )
        self.output_dim = _infer_output_dim(self.backbone, spec)
        self.window_length = _token_window_length(self.tokenizer, spec)

    def _device(self) -> torch.device:
        parameter = next(self.backbone.parameters(), None)
        if parameter is not None:
            return parameter.device
        buffer = next(self.backbone.buffers(), None)
        if buffer is not None:
            return buffer.device
        return torch.device("cpu")

    def _load_tokenizer_and_model(self, spec: ModelSpec) -> tuple[Any, nn.Module]:
        load_errors: list[str] = []
        candidates = [self.checkpoint]
        if (
            self.allow_remote_fallback
            and spec.hf_model_id is not None
            and spec.hf_model_id not in candidates
        ):
            candidates.append(spec.hf_model_id)

        for checkpoint in candidates:
            try:
                return _load_hf_backbone(spec, checkpoint=checkpoint, cache_dir=self.cache_dir)
            except Exception as exc:
                load_errors.append(f"{checkpoint}: {type(exc).__name__}: {exc}")

        raise RuntimeError(" | ".join(load_errors))

    def _chunk_sequences_by_bp(self, sequence: str) -> list[str]:
        raw_window = int(self.window_length)
        chunks = [
            sequence[start:start + raw_window]
            for start in range(0, len(sequence), raw_window)
        ]
        return chunks or [sequence]

    def _chunk_sequences_by_token(self, sequence: str) -> list[list[int]]:
        encoded = self.tokenizer(
            sequence,
            add_special_tokens=False,
            return_attention_mask=False,
            return_token_type_ids=False,
        )
        token_ids = list(encoded["input_ids"])
        window = int(self.window_length)
        chunks = [
            token_ids[start:start + window]
            for start in range(0, len(token_ids), window)
        ]
        return chunks or [token_ids]

    def _tokenize_bp_chunks(self, chunks: list[str]) -> dict[str, torch.Tensor]:
        batch = self.tokenizer(
            chunks,
            add_special_tokens=True,
            padding=True,
            truncation=True,
            max_length=self._max_length_for_tokenizer(),
            return_attention_mask=True,
            return_special_tokens_mask=True,
            return_tensors="pt",
        )
        batch.pop("token_type_ids", None)
        return _fill_special_tokens_mask(self.tokenizer, batch)

    def _tokenize_token_chunks(self, chunks: list[list[int]]) -> dict[str, torch.Tensor]:
        prepared = [
            self.tokenizer.prepare_for_model(
                token_ids,
                add_special_tokens=True,
                return_attention_mask=True,
                return_token_type_ids=False,
                truncation=True,
                max_length=self._max_length_for_tokenizer(),
            )
            for token_ids in chunks
        ]
        batch = self.tokenizer.pad(prepared, padding=True, return_tensors="pt")
        batch.pop("token_type_ids", None)
        return _fill_special_tokens_mask(self.tokenizer, batch)

    def _max_length_for_tokenizer(self) -> int:
        return int(_safe_model_max_length(self.tokenizer, self.spec))

    def _forward_backbone(self, inputs: dict[str, torch.Tensor]) -> Any:
        try:
            return self.backbone(**inputs)
        except TypeError as exc:
            if "unexpected keyword argument 'attention_mask'" not in str(exc):
                raise
        reduced_inputs = {"input_ids": inputs["input_ids"]}
        return self.backbone(**reduced_inputs)

    def _chunks_for_sequence(self, sequence: str) -> tuple[list[Any], Any]:
        sequence = str(sequence).upper()
        if self.spec.chunk_unit == "bp":
            chunks: list[Any] = self._chunk_sequences_by_bp(sequence)
            tokenize_fn = self._tokenize_bp_chunks
        else:
            chunks = self._chunk_sequences_by_token(sequence)
            tokenize_fn = self._tokenize_token_chunks
        return chunks, tokenize_fn

    def _encode_chunk_batch(
        self,
        chunks: list[Any],
        *,
        tokenize_fn: Any,
        token_readout: str,
    ) -> torch.Tensor:
        inputs = tokenize_fn(chunks)
        device = self._device()
        inputs = {
            key: value.to(device)
            for key, value in inputs.items()
            if isinstance(value, torch.Tensor)
        }
        special_tokens_mask = inputs.pop("special_tokens_mask", None)
        with self.forward_context():
            outputs = self._forward_backbone(inputs)
            hidden_states = _resolve_hidden_states(outputs)
        return _apply_token_readout(
            hidden_states,
            attention_mask=inputs.get("attention_mask"),
            special_tokens_mask=special_tokens_mask,
            token_readout=token_readout,
        )

    def _chunk_bucket_key(self, chunk: Any) -> int:
        return len(chunk)

    def _encode_single_sequence(self, sequence: str, *, token_readout: str) -> torch.Tensor:
        chunks, tokenize_fn = self._chunks_for_sequence(sequence)
        chunk_representations: list[torch.Tensor] = []
        for start in range(0, len(chunks), self.forward_batch_size):
            batch_chunks = chunks[start:start + self.forward_batch_size]
            chunk_vector = self._encode_chunk_batch(
                batch_chunks,
                tokenize_fn=tokenize_fn,
                token_readout=token_readout,
            )
            chunk_representations.append(chunk_vector)

        return torch.cat(chunk_representations, dim=0).mean(dim=0)

    def encode_serial(self, sequences: list[str], *, token_readout: str) -> torch.Tensor:
        return torch.stack(
            [
                self._encode_single_sequence(sequence, token_readout=token_readout)
                for sequence in sequences
            ],
            dim=0,
        )

    def _encode_batched_sequences(self, sequences: list[str], *, token_readout: str) -> torch.Tensor:
        if not sequences:
            return torch.empty((0, self.output_dim), device=self._device())

        tokenize_fn = (
            self._tokenize_bp_chunks
            if self.spec.chunk_unit == "bp"
            else self._tokenize_token_chunks
        )
        chunk_buckets: dict[int, list[tuple[int, Any]]] = {}
        for sequence_index, sequence in enumerate(sequences):
            chunks, _ = self._chunks_for_sequence(sequence)
            for chunk in chunks:
                chunk_buckets.setdefault(self._chunk_bucket_key(chunk), []).append((sequence_index, chunk))

        sequence_chunks: list[list[torch.Tensor]] = [[] for _ in sequences]
        for bucket_items in chunk_buckets.values():
            for start in range(0, len(bucket_items), self.forward_batch_size):
                batch_items = bucket_items[start:start + self.forward_batch_size]
                batch_chunks = [chunk for _, chunk in batch_items]
                chunk_vectors = self._encode_chunk_batch(
                    batch_chunks,
                    tokenize_fn=tokenize_fn,
                    token_readout=token_readout,
                )
                for (sequence_index, _), chunk_vector in zip(batch_items, chunk_vectors):
                    sequence_chunks[sequence_index].append(chunk_vector)

        features = []
        for chunks in sequence_chunks:
            if not chunks:
                raise RuntimeError("内部错误：序列没有生成任何 chunk 表示。")
            features.append(torch.stack(chunks, dim=0).mean(dim=0))
        return torch.stack(features, dim=0)

    def encode_batched(self, sequences: list[str], *, token_readout: str) -> torch.Tensor:
        return self._encode_batched_sequences(sequences, token_readout=token_readout)

    def encode(self, sequences: list[str], *, token_readout: str) -> torch.Tensor:
        if self.spec.name == "generator":
            return self.encode_serial(sequences, token_readout=token_readout)
        return self.encode_batched(sequences, token_readout=token_readout)


def _ensure_enformer_modules(repo_root: Path) -> tuple[Any, Any]:
    package_name = "enformer_pytorch"
    config_name = "enformer_pytorch.config_enformer"
    data_name = "enformer_pytorch.data"
    modeling_name = "enformer_pytorch.modeling_enformer"

    if modeling_name in sys.modules and config_name in sys.modules:
        return sys.modules[config_name], sys.modules[modeling_name]

    package = sys.modules.get(package_name)
    if package is None:
        package = types.ModuleType(package_name)
        package.__path__ = [str(repo_root / "enformer_pytorch")]
        sys.modules[package_name] = package

    if data_name not in sys.modules:
        data_module = types.ModuleType(data_name)

        seq_indices_embed = torch.zeros(256).long()
        for base, index in {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4, ".": -1}.items():
            seq_indices_embed[ord(base)] = index
            seq_indices_embed[ord(base.lower())] = index

        one_hot_embed = torch.zeros(256, 4)
        for base, vector in {
            "A": torch.tensor([1.0, 0.0, 0.0, 0.0]),
            "C": torch.tensor([0.0, 1.0, 0.0, 0.0]),
            "G": torch.tensor([0.0, 0.0, 1.0, 0.0]),
            "T": torch.tensor([0.0, 0.0, 0.0, 1.0]),
            "N": torch.tensor([0.0, 0.0, 0.0, 0.0]),
            ".": torch.tensor([0.25, 0.25, 0.25, 0.25]),
        }.items():
            one_hot_embed[ord(base)] = vector
            one_hot_embed[ord(base.lower())] = vector

        def _torch_fromstring(seq_strs):
            if isinstance(seq_strs, str):
                seq_strs = [seq_strs]
                batched = False
            else:
                batched = True
            arrays = [
                torch.tensor(list(seq.encode("ascii")), dtype=torch.uint8)
                for seq in seq_strs
            ]
            stacked = torch.stack(arrays, dim=0)
            return stacked if batched else stacked[0]

        def str_to_one_hot(seq_strs):
            seq_chrs = _torch_fromstring(seq_strs)
            return one_hot_embed[seq_chrs.long()]

        def seq_indices_to_one_hot(t, padding=-1):
            is_padding = t == padding
            t = t.clamp(min=0)
            one_hot = F.one_hot(t, num_classes=5)
            out = one_hot[..., :4].float()
            return out.masked_fill(is_padding[..., None], 0.25)

        data_module.str_to_one_hot = str_to_one_hot
        data_module.seq_indices_to_one_hot = seq_indices_to_one_hot
        sys.modules[data_name] = data_module

    if config_name not in sys.modules:
        _module_from_file(
            config_name,
            repo_root / "enformer_pytorch" / "config_enformer.py",
        )
    if modeling_name not in sys.modules:
        _module_from_file(
            modeling_name,
            repo_root / "enformer_pytorch" / "modeling_enformer.py",
        )

    return sys.modules[config_name], sys.modules[modeling_name]


class EnformerBackboneAdapter(BaseBackboneAdapter):
    def __init__(
        self,
        spec: ModelSpec,
        *,
        checkpoint: str,
        forward_batch_size: int,
        finetune_method: str,
        ia3_target_modules: tuple[str, ...] | None,
        ia3_feedforward_modules: tuple[str, ...] | None,
    ) -> None:
        super().__init__(spec, finetune_method=finetune_method)
        self.forward_batch_size = max(int(forward_batch_size), 1)
        repo_root = PROJECT_ROOT / "compared_models" / "models" / "enformer-pytorch"
        _, modeling_module = _ensure_enformer_modules(repo_root)
        model = modeling_module.from_pretrained(checkpoint, use_tf_gamma=True)
        self.backbone = _configure_backbone_finetuning(
            model,
            spec=spec,
            finetune_method=finetune_method,
            ia3_target_modules=ia3_target_modules,
            ia3_feedforward_modules=ia3_feedforward_modules,
            prefer_peft=False,
        )
        self.output_dim = int(spec.output_dim or self.backbone.dim * 2)
        one_hot_embed = torch.zeros(256, 4)
        for base, vector in {
            "A": torch.tensor([1.0, 0.0, 0.0, 0.0]),
            "C": torch.tensor([0.0, 1.0, 0.0, 0.0]),
            "G": torch.tensor([0.0, 0.0, 1.0, 0.0]),
            "T": torch.tensor([0.0, 0.0, 0.0, 1.0]),
            "N": torch.tensor([0.0, 0.0, 0.0, 0.0]),
            ".": torch.tensor([0.25, 0.25, 0.25, 0.25]),
        }.items():
            one_hot_embed[ord(base)] = vector
            one_hot_embed[ord(base.lower())] = vector
        self.register_buffer("_one_hot_embed", one_hot_embed, persistent=False)

    def _chunk_sequence(self, sequence: str) -> list[str]:
        sequence = str(sequence).upper()
        chunks = [
            sequence[start:start + self.spec.chunk_length]
            for start in range(0, len(sequence), self.spec.chunk_length)
        ]
        return chunks or [sequence]

    def _pad_chunk(self, chunk: str) -> str:
        if len(chunk) >= self.spec.chunk_length:
            return chunk[:self.spec.chunk_length]
        return chunk + ("N" * (self.spec.chunk_length - len(chunk)))

    def _one_hot_batch(self, chunks: list[str]) -> torch.Tensor:
        encoded = [
            torch.tensor(list(chunk.encode("ascii")), dtype=torch.long, device=self._one_hot_embed.device)
            for chunk in chunks
        ]
        return self._one_hot_embed[torch.stack(encoded, dim=0)]

    def _encode_chunk_batch(self, chunks: list[str]) -> torch.Tensor:
        batch_inputs = self._one_hot_batch(chunks)
        with self.forward_context():
            embeddings = self.backbone(batch_inputs, return_only_embeddings=True)
        if embeddings.ndim == 2:
            embeddings = embeddings.unsqueeze(0)
        return embeddings.mean(dim=1)

    def _encode_single_sequence(self, sequence: str) -> torch.Tensor:
        chunks = [self._pad_chunk(chunk) for chunk in self._chunk_sequence(sequence)]
        chunk_vectors: list[torch.Tensor] = []
        for start in range(0, len(chunks), self.forward_batch_size):
            batch_chunks = chunks[start:start + self.forward_batch_size]
            chunk_vectors.append(self._encode_chunk_batch(batch_chunks))
        return torch.cat(chunk_vectors, dim=0).mean(dim=0)

    def encode_serial(self, sequences: list[str], *, token_readout: str) -> torch.Tensor:
        del token_readout
        return torch.stack(
            [self._encode_single_sequence(sequence) for sequence in sequences],
            dim=0,
        )

    def encode_batched(self, sequences: list[str], *, token_readout: str) -> torch.Tensor:
        del token_readout
        if not sequences:
            return torch.empty((0, self.output_dim), device=self._one_hot_embed.device)

        chunk_buckets: dict[int, list[tuple[int, str]]] = {}
        for sequence_index, sequence in enumerate(sequences):
            for chunk in self._chunk_sequence(sequence):
                padded_chunk = self._pad_chunk(chunk)
                chunk_buckets.setdefault(len(padded_chunk), []).append((sequence_index, padded_chunk))

        sequence_chunks: list[list[torch.Tensor]] = [[] for _ in sequences]
        for bucket_items in chunk_buckets.values():
            for start in range(0, len(bucket_items), self.forward_batch_size):
                batch_items = bucket_items[start:start + self.forward_batch_size]
                batch_chunks = [chunk for _, chunk in batch_items]
                chunk_vectors = self._encode_chunk_batch(batch_chunks)
                for (sequence_index, _), chunk_vector in zip(batch_items, chunk_vectors):
                    sequence_chunks[sequence_index].append(chunk_vector)

        features = []
        for chunks in sequence_chunks:
            if not chunks:
                raise RuntimeError("内部错误：序列没有生成任何 Enformer chunk 表示。")
            features.append(torch.stack(chunks, dim=0).mean(dim=0))
        return torch.stack(features, dim=0)

    def encode(self, sequences: list[str], *, token_readout: str) -> torch.Tensor:
        return self.encode_batched(sequences, token_readout=token_readout)


def _clean_janusdna_config(config_dict: dict[str, Any]) -> dict[str, Any]:
    cleaned = dict(config_dict)
    cleaned.pop("_target_", None)
    for bool_like_key in (
        "bidirectional",
        "bidirectional_weight_tie",
        "layer_fusion",
        "final_attention",
        "mid_single_direction_attention",
        "bidirectional_attn_tie",
        "gradient_checkpointing",
        "output_router_logits",
    ):
        value = cleaned.get(bool_like_key)
        if isinstance(value, str):
            cleaned[bool_like_key] = value.lower().strip(", ") == "true"
    cleaned["gradient_checkpointing"] = False
    cleaned["output_router_logits"] = False
    cleaned["use_cache"] = False
    cleaned["use_mamba_kernels"] = False
    cleaned["attn_implementation"] = "eager"
    cleaned["final_attention_class"] = "eager"
    return cleaned


class JanusDNABackboneAdapter(HFTokenBackboneAdapter):
    def _chunk_bucket_key(self, chunk: Any) -> int:
        del chunk
        return 0

    def _load_tokenizer_and_model(self, spec: ModelSpec) -> tuple[Any, nn.Module]:
        checkpoint_path = Path(self.checkpoint)
        checkpoint_dir = checkpoint_path if checkpoint_path.is_dir() else checkpoint_path.parent
        repo_root = PROJECT_ROOT / "compared_models" / "models" / "JanusDNA"
        _prepend_repo_to_sys_path(repo_root)

        from janusdna.configuration_janusdna import JanusDNAConfig
        from janusdna.modeling_janusdna import JanusDNAForCausalLM
        from src.dataloaders.datasets.hg38_char_tokenizer import CharacterTokenizer

        config_path = checkpoint_dir / "model_config.json"
        if not config_path.exists():
            raise FileNotFoundError(f"缺少 JanusDNA model_config.json: {config_path}")
        config_payload = json.loads(config_path.read_text(encoding="utf-8"))
        config_data = _clean_janusdna_config(config_payload["config"])

        model_config = JanusDNAConfig(**config_data)
        model_config._attn_implementation = "eager"
        model_config.final_attention_class = "eager"

        model = JanusDNAForCausalLM(model_config)
        state_path = checkpoint_dir / "checkpoints" / "last.ckpt"
        if not state_path.exists():
            state_path = checkpoint_path
        state = torch.load(state_path, map_location="cpu")
        model_state_dict = dict(state["state_dict"])
        torch.nn.modules.utils.consume_prefix_in_state_dict_if_present(
            model_state_dict,
            "model.",
        )
        for key in list(model_state_dict.keys()):
            if ".self_attn.o_projs.0." in key:
                mapped_key = key.replace(".self_attn.o_projs.0.", ".self_attn.o_proj.")
                model_state_dict[mapped_key] = model_state_dict.pop(key)
        for key in list(model_state_dict.keys()):
            if "torchmetrics" in key:
                model_state_dict.pop(key)
        missing_keys, unexpected_keys = model.load_state_dict(model_state_dict, strict=False)
        ignored_missing_prefixes = ("lm_head.",)
        missing_keys = [
            key for key in missing_keys
            if not key.startswith(ignored_missing_prefixes)
        ]
        if missing_keys or unexpected_keys:
            raise RuntimeError(
                "JanusDNA 权重加载不完整: "
                f"missing={missing_keys[:10]}, unexpected={unexpected_keys[:10]}"
            )

        tokenizer = CharacterTokenizer(
            characters=["A", "C", "G", "T", "N"],
            model_max_length=int(spec.chunk_length) + 2,
            padding_side="right",
        )
        return tokenizer, model.model.eval()


def _build_alphagenome_apply_fn(
    dna_model_module: Any,
    metadata: Any,
    *,
    force_float32: bool,
):
    if not force_float32:
        _, apply_fn, _ = dna_model_module.create_model(metadata)
        return apply_fn

    hk = dna_model_module.hk
    jmp = dna_model_module.jmp
    research_model = dna_model_module.model

    jmp_policy = jmp.get_policy("params=float32,compute=float32,output=float32")

    @hk.transform_with_state
    def _forward(dna_sequence, organism_index):
        with hk.mixed_precision.push_policy(research_model.AlphaGenome, jmp_policy):
            return research_model.AlphaGenome(
                metadata,
                num_splice_sites=research_model.DEFAULT_NUM_SPLICE_SITES,
                splice_site_threshold=research_model.DEFAULT_SPLICE_SITE_THRESHOLD,
            )(dna_sequence, organism_index)

    def _apply_fn(params, state, dna_sequence, organism_index):
        (predictions, _), _ = _forward.apply(
            params,
            state,
            None,
            dna_sequence,
            organism_index,
        )
        return predictions

    return _apply_fn


def _cast_jax_tree_to_float32(tree: Any, *, jax_module: Any, jnp_module: Any) -> Any:
    def _cast_leaf(x: Any) -> Any:
        dtype = getattr(x, "dtype", None)
        if dtype is None:
            return x
        try:
            if jnp_module.issubdtype(dtype, jnp_module.floating):
                return x.astype(jnp_module.float32)
        except Exception:
            pass
        if str(dtype) == "bfloat16":
            return x.astype(jnp_module.float32)
        return x

    return jax_module.tree_util.tree_map(_cast_leaf, tree)


class AlphaGenomeBackboneAdapter(BaseBackboneAdapter):
    def __init__(
        self,
        spec: ModelSpec,
        *,
        checkpoint: str,
        finetune_method: str,
        forward_batch_size: int,
    ) -> None:
        super().__init__(spec, finetune_method=finetune_method)
        self.forward_batch_size = max(int(forward_batch_size), 1)
        checkpoint_path = Path(checkpoint)
        api_repo_root = PROJECT_ROOT / "compared_models" / "models" / "alphagenome"
        research_repo_root = PROJECT_ROOT / "compared_models" / "models" / "alphagenome_research"
        _prepend_repo_to_sys_path(api_repo_root)
        _prepend_repo_to_sys_path(research_repo_root)

        import jax
        import jax.numpy as jnp
        import orbax.checkpoint as ocp
        from alphagenome.models import dna_model as api_dna_model
        from alphagenome_research.model import dna_model
        from alphagenome_research.model.metadata import metadata as metadata_lib

        self._jax = jax
        self._jnp = jnp
        self._api_dna_model = api_dna_model
        self.output_dim = int(spec.output_dim or 1536)
        self.chunk_length = int(spec.chunk_length)
        self.min_sequence_length = 2048
        self.sequence_alignment = 2048
        if any(device.platform == "gpu" for device in jax.devices()):
            self.device = next(device for device in jax.devices() if device.platform == "gpu")
        else:
            self.device = jax.devices("cpu")[0]

        metadata = {
            organism: metadata_lib.load(organism)
            for organism in api_dna_model.Organism
        }
        checkpointer = ocp.StandardCheckpointer()
        params, state = checkpointer.restore(
            str(checkpoint_path),
            target=None,
            strict=False,
        )
        params = _cast_jax_tree_to_float32(params, jax_module=jax, jnp_module=jnp)
        state = _cast_jax_tree_to_float32(state, jax_module=jax, jnp_module=jnp)
        self.params = jax.device_put(params, self.device)
        self.state = jax.device_put(state, self.device)
        apply_fn = _build_alphagenome_apply_fn(
            dna_model,
            metadata,
            force_float32=(self.device.platform != "gpu"),
        )
        self.apply_fn = jax.jit(apply_fn, device=self.device)

    def _chunk_sequence(self, sequence: str) -> list[str]:
        sequence = str(sequence).upper()
        chunks = [
            sequence[start:start + self.chunk_length]
            for start in range(0, len(sequence), self.chunk_length)
        ]
        return chunks or [sequence]

    def _pad_chunk(self, chunk: str) -> str:
        effective_length = max(len(chunk), self.min_sequence_length)
        effective_length = int(math.ceil(effective_length / self.sequence_alignment) * self.sequence_alignment)
        effective_length = min(effective_length, self.chunk_length)
        if len(chunk) >= effective_length:
            return chunk[:effective_length]
        return chunk + ("N" * (effective_length - len(chunk)))

    def _one_hot(self, sequence: str):
        import numpy as np

        base_to_vec = {
            "A": [1.0, 0.0, 0.0, 0.0],
            "C": [0.0, 1.0, 0.0, 0.0],
            "G": [0.0, 0.0, 1.0, 0.0],
            "T": [0.0, 0.0, 0.0, 1.0],
            "N": [0.0, 0.0, 0.0, 0.0],
        }
        encoded = np.asarray(
            [base_to_vec.get(base, base_to_vec["N"]) for base in str(sequence).upper()],
            dtype=np.float32,
        )
        return encoded[None, ...]

    def _one_hot_batch(self, sequences: list[str]):
        import numpy as np

        base_to_vec = {
            "A": [1.0, 0.0, 0.0, 0.0],
            "C": [0.0, 1.0, 0.0, 0.0],
            "G": [0.0, 0.0, 1.0, 0.0],
            "T": [0.0, 0.0, 0.0, 1.0],
            "N": [0.0, 0.0, 0.0, 0.0],
        }
        encoded = [
            np.asarray(
                [base_to_vec.get(base, base_to_vec["N"]) for base in str(sequence).upper()],
                dtype=np.float32,
            )
            for sequence in sequences
        ]
        return np.stack(encoded, axis=0)

    def _encode_chunk_batch(self, chunks: list[str]) -> list[torch.Tensor]:
        import numpy as np

        dna_sequence = self._jax.device_put(self._one_hot_batch(chunks), self.device)
        organism_index = self._jax.device_put(
            self._jnp.asarray(
                [self._api_dna_model.Organism.HOMO_SAPIENS.value] * len(chunks),
                dtype=self._jnp.int32,
            ),
            self.device,
        )
        predictions = self.apply_fn(self.params, self.state, dna_sequence, organism_index)
        embeddings = np.asarray(predictions["embeddings_1bp"], dtype=np.float32)
        return [
            torch.from_numpy(chunk_embeddings.mean(axis=0))
            for chunk_embeddings in embeddings
        ]

    def _encode_single_sequence(self, sequence: str) -> torch.Tensor:
        import numpy as np

        chunk_vectors: list[torch.Tensor] = []
        for chunk in self._chunk_sequence(sequence):
            dna_sequence = self._jax.device_put(self._one_hot(self._pad_chunk(chunk)), self.device)
            organism_index = self._jax.device_put(
                self._jnp.asarray(
                    [self._api_dna_model.Organism.HOMO_SAPIENS.value],
                    dtype=self._jnp.int32,
                ),
                self.device,
            )
            predictions = self.apply_fn(self.params, self.state, dna_sequence, organism_index)
            embeddings = np.asarray(predictions["embeddings_1bp"], dtype=np.float32)
            chunk_vectors.append(torch.from_numpy(embeddings.mean(axis=1).squeeze(0)))
        return torch.stack(chunk_vectors, dim=0).mean(dim=0)

    def encode_serial(self, sequences: list[str], *, token_readout: str) -> torch.Tensor:
        del token_readout
        return torch.stack(
            [self._encode_single_sequence(sequence) for sequence in sequences],
            dim=0,
        )

    def encode_batched(self, sequences: list[str], *, token_readout: str) -> torch.Tensor:
        del token_readout
        if not sequences:
            return torch.empty((0, self.output_dim))

        chunk_buckets: dict[int, list[tuple[int, str]]] = {}
        for sequence_index, sequence in enumerate(sequences):
            for chunk in self._chunk_sequence(sequence):
                padded_chunk = self._pad_chunk(chunk)
                chunk_buckets.setdefault(len(padded_chunk), []).append((sequence_index, padded_chunk))

        sequence_chunks: list[list[torch.Tensor]] = [[] for _ in sequences]
        for bucket_items in chunk_buckets.values():
            for start in range(0, len(bucket_items), self.forward_batch_size):
                batch_items = bucket_items[start:start + self.forward_batch_size]
                batch_chunks = [chunk for _, chunk in batch_items]
                chunk_vectors = self._encode_chunk_batch(batch_chunks)
                for (sequence_index, _), chunk_vector in zip(batch_items, chunk_vectors):
                    sequence_chunks[sequence_index].append(chunk_vector)

        features = []
        for chunks in sequence_chunks:
            if not chunks:
                raise RuntimeError("内部错误：序列没有生成任何 AlphaGenome chunk 表示。")
            features.append(torch.stack(chunks, dim=0).mean(dim=0))
        return torch.stack(features, dim=0)

    def encode(self, sequences: list[str], *, token_readout: str) -> torch.Tensor:
        return self.encode_serial(sequences, token_readout=token_readout)


class Evo2BackboneAdapter(BaseBackboneAdapter):
    def __init__(
        self,
        spec: ModelSpec,
        *,
        checkpoint: str,
        finetune_method: str,
    ) -> None:
        super().__init__(spec, finetune_method=finetune_method)
        checkpoint_path = Path(checkpoint)
        repo_root = PROJECT_ROOT / "compared_models" / "models" / "evo2"
        _prepend_repo_to_sys_path(repo_root)

        import yaml
        from evo2 import Evo2
        from evo2.utils import CONFIG_MAP

        config_name = CONFIG_MAP.get("evo2_7b")
        if config_name is None:
            raise RuntimeError("Evo2 配置映射中缺少 evo2_7b。")
        config = yaml.safe_load((repo_root / "evo2" / config_name).read_text(encoding="utf-8"))
        self.layer_name = f"blocks.{int(config['num_layers']) - 1}.mlp.l3"
        evo2_model = Evo2(model_name="evo2_7b", local_path=str(checkpoint_path))
        self.backbone = evo2_model.model.eval()
        self.tokenizer = evo2_model.tokenizer
        for parameter in self.backbone.parameters():
            parameter.requires_grad = False
        self.output_dim = int(spec.output_dim or config["hidden_size"])
        self.chunk_length = int(spec.chunk_length)

    def _device(self) -> torch.device:
        parameter = next(self.backbone.parameters(), None)
        if parameter is not None:
            return parameter.device
        buffer = next(self.backbone.buffers(), None)
        if buffer is not None:
            return buffer.device
        return torch.device("cpu")

    def _chunk_sequence(self, sequence: str) -> list[str]:
        sequence = str(sequence).upper()
        chunks = [
            sequence[start:start + self.chunk_length]
            for start in range(0, len(sequence), self.chunk_length)
        ]
        return chunks or [sequence]

    def _readout(self, hidden_states: torch.Tensor, token_readout: str) -> torch.Tensor:
        if token_readout == "mean":
            return hidden_states.mean(dim=1)
        if token_readout == "last":
            return hidden_states[:, -1, :]
        raise ValueError(f"不支持的 Evo2 token_readout: {token_readout}")

    def _encode_single_sequence(self, sequence: str, *, token_readout: str) -> torch.Tensor:
        chunk_vectors: list[torch.Tensor] = []
        device = self._device()
        for chunk in self._chunk_sequence(sequence):
            input_ids = torch.tensor(
                self.tokenizer.tokenize(chunk),
                dtype=torch.long,
                device=device,
            ).unsqueeze(0)
            embeddings: dict[str, torch.Tensor] = {}

            def hook(_, __, output):
                if isinstance(output, tuple):
                    output = output[0]
                embeddings[self.layer_name] = output.detach()

            handle = self.backbone.get_submodule(self.layer_name).register_forward_hook(hook)
            try:
                with torch.no_grad():
                    _ = self.backbone.forward(input_ids)
            finally:
                handle.remove()
            hidden_states = embeddings[self.layer_name]
            chunk_vectors.append(self._readout(hidden_states, token_readout).squeeze(0).cpu())
        return torch.stack(chunk_vectors, dim=0).mean(dim=0)

    def encode(self, sequences: list[str], *, token_readout: str) -> torch.Tensor:
        return torch.stack(
            [self._encode_single_sequence(sequence, token_readout=token_readout) for sequence in sequences],
            dim=0,
        )


class FrozenSequenceLinearProbe(nn.Module):
    def __init__(
        self,
        backbone: BaseBackboneAdapter,
        *,
        num_labels: int,
        token_readout: str,
        finetune_method: str,
        pair_mode: bool = False,
    ) -> None:
        super().__init__()
        self.backbone = backbone
        self.token_readout = token_readout
        self.finetune_method = finetune_method
        self.pair_mode = bool(pair_mode)
        head_input_dim = int(backbone.output_dim) * (3 if self.pair_mode else 1)
        self.head = nn.Linear(head_input_dim, int(num_labels))

    @property
    def output_dim(self) -> int:
        return int(self.backbone.output_dim)

    def train(self, mode: bool = True):
        super().train(mode)
        if self.finetune_method == "frozen_linear_probe":
            self.backbone.eval()
        return self

    def forward(
        self,
        sequences: list[str],
        sequences_alt: list[str] | None = None,
    ) -> torch.Tensor:
        if self.pair_mode:
            if sequences_alt is None:
                raise ValueError("pair_mode=True 时必须提供 sequences_alt。")
            features_ref = self.backbone.encode(sequences, token_readout=self.token_readout)
            features_alt = self.backbone.encode(sequences_alt, token_readout=self.token_readout)
            features = torch.cat(
                [features_ref, features_alt, features_alt - features_ref],
                dim=-1,
            )
        else:
            features = self.backbone.encode(sequences, token_readout=self.token_readout)
        if (
            features.device != self.head.weight.device
            or features.dtype != self.head.weight.dtype
        ):
            features = features.to(
                device=self.head.weight.device,
                dtype=self.head.weight.dtype,
            )
        return self.head(features)


def build_linear_probe_model(
    *,
    model_name: str,
    num_labels: int,
    token_readout: str = "auto",
    finetune_method: str = "frozen_linear_probe",
    cache_dir: str | None = None,
    allow_remote_fallback: bool = True,
    chunk_forward_batch_size: int = 8,
    checkpoint_override: str | None = None,
    pair_mode: bool = False,
    ia3_target_modules: str | list[str] | tuple[str, ...] | None = None,
    ia3_feedforward_modules: str | list[str] | tuple[str, ...] | None = None,
) -> FrozenSequenceLinearProbe:
    spec = get_model_spec(model_name)
    if spec.loader_kind == "unsupported":
        raise RuntimeError(spec.unsupported_reason or f"模型 {spec.name} 当前未接入 benchmark 脚本。")
    resolved_finetune_method = validate_finetune_method(spec, finetune_method)
    resolved_ia3_target_modules = _parse_module_list(ia3_target_modules)
    resolved_ia3_feedforward_modules = _parse_module_list(ia3_feedforward_modules)

    checkpoint = checkpoint_override
    if checkpoint is None:
        if spec.local_checkpoint is not None:
            checkpoint = str(PROJECT_ROOT / spec.local_checkpoint)
        elif spec.hf_model_id is not None:
            checkpoint = spec.hf_model_id
        else:
            raise RuntimeError(f"模型 {spec.name} 未配置 checkpoint。")

    resolved_readout = spec.default_token_readout if token_readout == "auto" else token_readout
    if spec.loader_kind == "hf_token":
        backbone = HFTokenBackboneAdapter(
            spec,
            checkpoint=checkpoint,
            token_readout=resolved_readout,
            cache_dir=cache_dir,
            allow_remote_fallback=allow_remote_fallback,
            forward_batch_size=chunk_forward_batch_size,
            finetune_method=resolved_finetune_method,
            ia3_target_modules=resolved_ia3_target_modules,
            ia3_feedforward_modules=resolved_ia3_feedforward_modules,
        )
    elif spec.loader_kind == "enformer":
        backbone = EnformerBackboneAdapter(
            spec,
            checkpoint=checkpoint,
            forward_batch_size=chunk_forward_batch_size,
            finetune_method=resolved_finetune_method,
            ia3_target_modules=resolved_ia3_target_modules,
            ia3_feedforward_modules=resolved_ia3_feedforward_modules,
        )
        resolved_readout = "embedding"
    elif spec.loader_kind == "janusdna":
        backbone = JanusDNABackboneAdapter(
            spec,
            checkpoint=checkpoint,
            token_readout=resolved_readout,
            cache_dir=cache_dir,
            allow_remote_fallback=False,
            forward_batch_size=chunk_forward_batch_size,
            finetune_method=resolved_finetune_method,
            ia3_target_modules=resolved_ia3_target_modules,
            ia3_feedforward_modules=resolved_ia3_feedforward_modules,
        )
    elif spec.loader_kind == "alphagenome":
        backbone = AlphaGenomeBackboneAdapter(
            spec,
            checkpoint=checkpoint,
            finetune_method=resolved_finetune_method,
            forward_batch_size=chunk_forward_batch_size,
        )
        resolved_readout = "embedding"
    elif spec.loader_kind == "evo2":
        backbone = Evo2BackboneAdapter(
            spec,
            checkpoint=checkpoint,
            finetune_method=resolved_finetune_method,
        )
    else:
        raise RuntimeError(f"未实现的 loader_kind: {spec.loader_kind}")

    return FrozenSequenceLinearProbe(
        backbone,
        num_labels=num_labels,
        token_readout=resolved_readout,
        finetune_method=resolved_finetune_method,
        pair_mode=pair_mode,
    )
