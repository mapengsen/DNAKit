"""Checkpoint resolution for official direct-prediction task heads."""

from __future__ import annotations

import importlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from dnakit.download import DownloadConfig, RemoteFile, download_file
from dnakit.exceptions import BackendUnavailableError, ConfigurationError, DownloadError
from dnakit.representations import default_checkpoint_root, ensure_model_checkpoint

from .models import PropertyPredictionConfig

_HF_PATTERNS = (
    "*.bin",
    "*.json",
    "*.model",
    "*.pt",
    "*.py",
    "*.safetensors",
    "*.txt",
)
_ENFORMER_TARGET_URLS = {
    "human": (
        "https://raw.githubusercontent.com/calico/basenji/0.5/"
        "manuscripts/cross2020/targets_human.txt"
    ),
    "mouse": (
        "https://raw.githubusercontent.com/calico/basenji/0.5/"
        "manuscripts/cross2020/targets_mouse.txt"
    ),
}


@dataclass(frozen=True, slots=True)
class PredictionCheckpointInfo:
    """Local paths used by one prediction backend."""

    paths: tuple[str, ...]
    downloaded: bool
    sources: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        """Return JSON-compatible checkpoint provenance."""

        return {
            "paths": self.paths,
            "downloaded": self.downloaded,
            "sources": self.sources,
        }


def default_prediction_checkpoint_root() -> Path:
    """Return the shared user-visible ``./ckpt`` directory."""

    return default_checkpoint_root()


def _root(config: PropertyPredictionConfig) -> Path:
    return (
        default_prediction_checkpoint_root()
        if config.checkpoint_dir is None
        else Path(config.checkpoint_dir).expanduser().resolve()
    )


def _require_huggingface_hub() -> Any:
    try:
        return importlib.import_module("huggingface_hub")
    except ImportError as exc:
        raise BackendUnavailableError(
            "Prediction checkpoint download requires huggingface-hub.",
            code="MISSING_NEURAL_DEPENDENCY",
            hint='Install the neural extra with: python -m pip install "dnakit[neural]"',
        ) from exc


def _snapshot(
    repo_id: str,
    target: Path,
    *,
    token: str | None,
    show_progress: bool,
    allow_patterns: tuple[str, ...] = _HF_PATTERNS,
) -> None:
    hub = _require_huggingface_hub()
    target.mkdir(parents=True, exist_ok=True)
    progress_was_disabled = False
    progress_utils: Any | None = None
    if not show_progress:
        try:
            progress_utils = importlib.import_module("huggingface_hub.utils")
            progress_was_disabled = bool(progress_utils.are_progress_bars_disabled())
            progress_utils.disable_progress_bars()
        except (ImportError, AttributeError):
            progress_utils = None
    try:
        hub.snapshot_download(
            repo_id=repo_id,
            repo_type="model",
            local_dir=str(target),
            allow_patterns=list(allow_patterns),
            token=token,
            max_workers=4,
        )
    except Exception as exc:
        raise DownloadError(
            "Could not download a direct-prediction checkpoint from Hugging Face.",
            code="MODEL_CHECKPOINT_DOWNLOAD_FAILED",
            context={"checkpoint": repo_id},
            hint="Check network access, model terms, and hf_token.",
        ) from exc
    finally:
        if not show_progress and not progress_was_disabled and progress_utils is not None:
            progress_utils.enable_progress_bars()


def _has_weights(path: Path) -> bool:
    return any(
        any(path.glob(pattern))
        for pattern in ("*.bin", "*.pt", "*.safetensors", "*.safetensors.index.json")
    )


def _segmentnt_checkpoint(config: PropertyPredictionConfig) -> PredictionCheckpointInfo:
    target = _root(config) / "segment-nt"
    ready = (
        (target / "config.json").is_file()
        and (target / "modeling_segment_nt.py").is_file()
        and _has_weights(target)
    )
    if not ready:
        _snapshot(
            "InstaDeepAI/segment_nt",
            target,
            token=config.hf_token,
            show_progress=config.show_progress,
        )
    if not (
        (target / "config.json").is_file()
        and (target / "modeling_segment_nt.py").is_file()
        and _has_weights(target)
    ):
        raise DownloadError(
            "SegmentNT checkpoint download did not produce the required files.",
            code="INCOMPLETE_MODEL_CHECKPOINT",
        )
    return PredictionCheckpointInfo(
        (str(target),),
        not ready,
        ("https://huggingface.co/InstaDeepAI/segment_nt",),
    )


def _evo2_exon_checkpoint(config: PropertyPredictionConfig) -> PredictionCheckpointInfo:
    target = _root(config) / "evo2-exon-classifier"
    base = target / "base"
    head = target / "head"
    base_ready = (base / "evo2_7b_base.pt").is_file()
    head_ready = (head / "config.json").is_file() and _has_weights(head)
    if not base_ready:
        _snapshot(
            "arcinstitute/evo2_7b_base",
            base,
            token=config.hf_token,
            show_progress=config.show_progress,
            allow_patterns=("config.json", "evo2_7b_base.pt"),
        )
    if not head_ready:
        _snapshot(
            "schmojo/evo2-exon-classifier",
            head,
            token=config.hf_token,
            show_progress=config.show_progress,
        )
    if not (base / "evo2_7b_base.pt").is_file() or not (
        (head / "config.json").is_file() and _has_weights(head)
    ):
        raise DownloadError(
            "Evo 2 exon-classifier checkpoint download is incomplete.",
            code="INCOMPLETE_MODEL_CHECKPOINT",
        )
    return PredictionCheckpointInfo(
        (str(base), str(head)),
        not (base_ready and head_ready),
        (
            "https://huggingface.co/arcinstitute/evo2_7b_base",
            "https://huggingface.co/schmojo/evo2-exon-classifier",
        ),
    )


def _enformer_metadata(checkpoint: Path) -> None:
    for organism, url in _ENFORMER_TARGET_URLS.items():
        target = checkpoint / f"targets_{organism}.txt"
        if target.is_file():
            continue
        try:
            download_file(
                RemoteFile(url, filename=target.name),
                target,
                config=DownloadConfig(
                    timeout=60.0,
                    max_file_bytes=10_000_000,
                    max_total_bytes=10_000_000,
                ),
            )
        except (DownloadError, FileExistsError, OSError):
            # Track predictions remain valid without names; the backend then
            # records target indices and the official metadata URL.
            continue


def _representation_checkpoint(
    config: PropertyPredictionConfig,
) -> PredictionCheckpointInfo:
    checkpoint = ensure_model_checkpoint(
        config.model,
        checkpoint_dir=config.checkpoint_dir,
        hf_token=config.hf_token,
        show_progress=config.show_progress,
    )
    if config.model == "enformer":
        _enformer_metadata(Path(checkpoint.path))
    return PredictionCheckpointInfo(
        (checkpoint.path,),
        checkpoint.downloaded,
        (checkpoint.source,),
    )


def _explicit_checkpoint(config: PropertyPredictionConfig) -> PredictionCheckpointInfo:
    assert config.checkpoint_path is not None
    path = Path(config.checkpoint_path).expanduser().resolve()
    if not path.is_dir():
        raise ConfigurationError(
            "checkpoint_path must be an existing checkpoint directory.",
            code="MODEL_CHECKPOINT_NOT_FOUND",
            context={"path": str(path)},
        )
    if config.model == "evo2" and config.task == "exon_probability":
        base = path / "base"
        head = path / "head"
        if not (base / "evo2_7b_base.pt").is_file() or not (
            (head / "config.json").is_file() and _has_weights(head)
        ):
            raise ConfigurationError(
                "The Evo 2 exon checkpoint root must contain base/ and head/.",
                code="MODEL_CHECKPOINT_NOT_FOUND",
                context={"path": str(path)},
            )
        return PredictionCheckpointInfo(
            (str(base), str(head)),
            False,
            ("explicit", "explicit"),
        )
    if config.model == "enformer":
        _enformer_metadata(path)
    return PredictionCheckpointInfo((str(path),), False, ("explicit",))


def ensure_prediction_checkpoint(
    config: PropertyPredictionConfig,
) -> PredictionCheckpointInfo:
    """Resolve and, when needed, download checkpoints for one direct task."""

    if not isinstance(config, PropertyPredictionConfig):
        raise ConfigurationError(
            "config must be PropertyPredictionConfig.", code="INVALID_PREDICTION_CONFIG"
        )
    if config.model == "lucaone":
        return PredictionCheckpointInfo((), False, ("official-source-download",))
    if config.checkpoint_path is not None:
        return _explicit_checkpoint(config)
    if config.model == "segmentnt":
        return _segmentnt_checkpoint(config)
    if config.model == "evo2" and config.task == "exon_probability":
        return _evo2_exon_checkpoint(config)
    if config.model in {"alphagenome", "enformer", "evo2", "generator"}:
        return _representation_checkpoint(config)
    raise AssertionError(f"Unsupported prediction checkpoint strategy: {config.model}")


__all__ = [
    "PredictionCheckpointInfo",
    "default_prediction_checkpoint_root",
    "ensure_prediction_checkpoint",
]
