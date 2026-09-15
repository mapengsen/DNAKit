"""Checkpoint resolution for official direct-prediction task heads."""

from __future__ import annotations

import importlib
import os
import tempfile
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from dnakit._checkpoint_status import checkpoint_status
from dnakit._dependency_hints import dependency_hint
from dnakit._verbosity import progress_enabled
from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.download import DownloadConfig, RemoteFile, download_file
from dnakit.exceptions import BackendUnavailableError, ConfigurationError, DownloadError
from dnakit.representations import default_checkpoint_root, ensure_model_checkpoint
from dnakit.representations.checkpoints import _indexed_weights_are_complete

from .enformer_benchmarks import (
    EnformerBenchmarkTask,
    get_enformer_benchmark_task,
    is_enformer_benchmark_task,
)
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


def default_prediction_checkpoint_root(*, verbose: bool = True) -> Path:
    """Return the shared user-visible ``./ckpt`` directory."""

    return default_checkpoint_root()


def enformer_benchmark_checkpoint_path(
    task: str, *, checkpoint_dir: str | Path | None = None, verbose: bool = True
) -> Path:
    """Return the canonical local path for one of the 27 Enformer checkpoints."""

    spec = get_enformer_benchmark_task(task)
    if checkpoint_dir is None:
        root = default_prediction_checkpoint_root() / "enformer-benchmarks"
    else:
        root = Path(checkpoint_dir).expanduser().resolve()
    return (root / spec.checkpoint_filename).resolve()


def _root(config: PropertyPredictionConfig) -> Path:
    return (
        default_prediction_checkpoint_root()
        if config.checkpoint_dir is None
        else Path(config.checkpoint_dir).expanduser().resolve()
    )


def _benchmark_candidates(root: Path, spec: EnformerBenchmarkTask) -> tuple[Path, ...]:
    family_dir = "nt-revised" if spec.family == "nt" else "genomic-benchmarks"
    return (
        root / spec.checkpoint_filename,
        root / "enformer-benchmarks" / spec.checkpoint_filename,
        root / family_dir / spec.checkpoint_filename,
    )


def _benchmark_file_ready(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def _download_benchmark_checkpoint(path: Path, spec: EnformerBenchmarkTask) -> None:
    try:
        gdown = importlib.import_module("gdown")
    except ImportError as exc:
        raise BackendUnavailableError(
            "Enformer task checkpoint download requires gdown.",
            code="MISSING_NEURAL_DEPENDENCY",
            hint=dependency_hint("neural,neural-enformer"),
        ) from exc
    checkpoint_status("Starting checkpoint download", path, spec.checkpoint_url)
    temporary: Path | None = None
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        # Supply a stream so gdown cannot leave its own untracked partial files.
        with tempfile.NamedTemporaryFile(
            mode="w+b",
            prefix=f".{path.name}.",
            suffix=".part",
            dir=path.parent,
            delete=False,
        ) as output:
            temporary = Path(output.name)
            result = gdown.download(
                id=spec.checkpoint_file_id,
                output=output,
                quiet=not progress_enabled(),
                use_cookies=False,
            )
        # These checkpoints use torch's ZIP serialization (required by mmap).
        # Checking the archive directory rejects HTML, empty and truncated files
        # without deserializing tensors or reading the entire checkpoint.
        if result is None or not zipfile.is_zipfile(temporary):
            raise ValueError("Download is not a complete PyTorch checkpoint archive")
        os.replace(temporary, path)
    except Exception as exc:
        checkpoint_status("Checkpoint download failed", path, spec.checkpoint_url)
        raise DownloadError(
            "Could not download the Enformer task checkpoint.",
            code="MODEL_CHECKPOINT_DOWNLOAD_FAILED",
            context={
                "task": spec.name,
                "download_url": spec.checkpoint_url,
                "expected_path": str(path),
            },
            hint=(
                f"Check network access and Google Drive sharing/quota, or download "
                f"{spec.checkpoint_url} manually to {path}."
            ),
        ) from exc
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    checkpoint_status("Checkpoint download complete", path, spec.checkpoint_url)


def _resolve_benchmark_file(
    path: Path, spec: EnformerBenchmarkTask, *, explicit: bool = False
) -> PredictionCheckpointInfo:
    ready = _benchmark_file_ready(path)
    if ready:
        checkpoint_status("Using local checkpoint", path, None if explicit else spec.checkpoint_url)
    else:
        _download_benchmark_checkpoint(path, spec)
    source = "explicit" if ready and explicit else spec.checkpoint_url
    return PredictionCheckpointInfo((str(path),), not ready, (source,))


def _enformer_benchmark_checkpoint(
    config: PropertyPredictionConfig,
) -> PredictionCheckpointInfo:
    spec = get_enformer_benchmark_task(config.task)
    candidates = _benchmark_candidates(_root(config), spec)
    if config.checkpoint_dir is None:
        candidates = (candidates[1], candidates[0], candidates[2])
    path = next(
        (candidate.resolve() for candidate in candidates if _benchmark_file_ready(candidate)),
        candidates[0].resolve(),
    )
    return _resolve_benchmark_file(path, spec)


def _require_huggingface_hub() -> Any:
    try:
        return importlib.import_module("huggingface_hub")
    except ImportError as exc:
        raise BackendUnavailableError(
            "Prediction checkpoint download requires huggingface-hub.",
            code="MISSING_NEURAL_DEPENDENCY",
            hint=dependency_hint("neural"),
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
    if any(
        (path / name).is_file()
        for name in ("model.safetensors.index.json", "pytorch_model.bin.index.json")
    ):
        return _indexed_weights_are_complete(path)
    return any(any(path.glob(pattern)) for pattern in ("*.bin", "*.pt", "*.safetensors"))


def _segmentnt_checkpoint(config: PropertyPredictionConfig) -> PredictionCheckpointInfo:
    target = _root(config) / "segment-nt"
    ready = (
        (target / "config.json").is_file()
        and (target / "modeling_segment_nt.py").is_file()
        and _has_weights(target)
    )
    if not ready:
        checkpoint_status(
            "Starting checkpoint download", target, "https://huggingface.co/InstaDeepAI/segment_nt"
        )
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
    checkpoint_status(
        "Using local checkpoint" if ready else "Checkpoint download complete",
        target,
        "https://huggingface.co/InstaDeepAI/segment_nt",
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
        checkpoint_status(
            "Starting checkpoint download", base, "https://huggingface.co/arcinstitute/evo2_7b_base"
        )
        _snapshot(
            "arcinstitute/evo2_7b_base",
            base,
            token=config.hf_token,
            show_progress=config.show_progress,
            allow_patterns=("config.json", "evo2_7b_base.pt"),
        )
    if not head_ready:
        checkpoint_status(
            "Starting checkpoint download",
            head,
            "https://huggingface.co/schmojo/evo2-exon-classifier",
        )
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
    checkpoint_status(
        "Using local checkpoint" if base_ready else "Checkpoint download complete",
        base,
        "https://huggingface.co/arcinstitute/evo2_7b_base",
    )
    checkpoint_status(
        "Using local checkpoint" if head_ready else "Checkpoint download complete",
        head,
        "https://huggingface.co/schmojo/evo2-exon-classifier",
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
    if config.model == "enformer" and is_enformer_benchmark_task(config.task):
        spec = get_enformer_benchmark_task(config.task)
        if path.is_dir():
            path = next(
                (
                    candidate.resolve()
                    for candidate in _benchmark_candidates(path, spec)
                    if _benchmark_file_ready(candidate)
                ),
                path / spec.checkpoint_filename,
            )
        return _resolve_benchmark_file(path, spec, explicit=True)
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
        checkpoint_status("Using local checkpoint", base)
        checkpoint_status("Using local checkpoint", head)
        return PredictionCheckpointInfo(
            (str(base), str(head)),
            False,
            ("explicit", "explicit"),
        )
    if config.model == "enformer":
        _enformer_metadata(path)
    checkpoint_status("Using local checkpoint", path)
    return PredictionCheckpointInfo((str(path),), False, ("explicit",))


def ensure_prediction_checkpoint(
    config: PropertyPredictionConfig, *, verbose: bool = True
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
    if config.model == "enformer" and is_enformer_benchmark_task(config.task):
        return _enformer_benchmark_checkpoint(config)
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
    "enformer_benchmark_checkpoint_path",
    "ensure_prediction_checkpoint",
]

_register_public_verbose(globals())
