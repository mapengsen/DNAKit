"""On-demand checkpoint downloads into a user-visible ``ckpt`` directory."""

from __future__ import annotations

import hashlib
import importlib
import json
import os
import shutil
import tarfile
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from dnakit.download import DownloadConfig, DownloadProgress, RemoteFile, download_file
from dnakit.exceptions import BackendUnavailableError, DownloadError

from .models import DNAEmbeddingModel, get_embedding_model

_MANIFEST_NAME = ".dnakit-checkpoint.json"


@dataclass(frozen=True, slots=True)
class CheckpointInfo:
    """Resolved checkpoint path and whether this call downloaded it."""

    model_name: str
    path: str
    source: str
    source_repository: str
    downloaded: bool
    gated: bool

    def to_dict(self) -> dict[str, object]:
        """Return JSON-compatible checkpoint provenance."""

        return {
            "model_name": self.model_name,
            "path": self.path,
            "source": self.source,
            "source_repository": self.source_repository,
            "downloaded": self.downloaded,
            "gated": self.gated,
        }


def default_checkpoint_root() -> Path:
    """Return ``./ckpt`` using the working directory at call time."""

    return (Path.cwd() / "ckpt").resolve()


def _checkpoint_root(value: str | os.PathLike[str] | None) -> Path:
    return default_checkpoint_root() if value is None else Path(value).expanduser().resolve()


def _read_json(path: Path) -> dict[str, Any] | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError):
        return None
    return payload if isinstance(payload, dict) else None


def _indexed_weights_are_complete(target: Path) -> bool:
    for index_path in (
        target / "model.safetensors.index.json",
        target / "pytorch_model.bin.index.json",
    ):
        if not index_path.is_file():
            continue
        payload = _read_json(index_path)
        weight_map = None if payload is None else payload.get("weight_map")
        if not isinstance(weight_map, dict) or not weight_map:
            return False
        filenames = {value for value in weight_map.values() if isinstance(value, str)}
        return bool(filenames) and all((target / filename).is_file() for filename in filenames)
    return False


def _checkpoint_ready(spec: DNAEmbeddingModel, target: Path) -> bool:
    if not target.is_dir():
        return False
    if spec.loader == "janusdna":
        return (target / "model_config.json").is_file() and (
            target / "checkpoints" / "last.ckpt"
        ).is_file()
    if spec.loader == "alphagenome":
        return (target / "_CHECKPOINT_METADATA").is_file() and (target / "manifest.ocdbt").is_file()
    if spec.loader == "evo2":
        return (target / "evo2_7b.pt").is_file()
    if not (target / "config.json").is_file():
        return False
    if _indexed_weights_are_complete(target):
        return True
    weight_patterns = ("*.safetensors", "*.bin", "*.pt", "*.ckpt")
    return any(any(target.glob(pattern)) for pattern in weight_patterns)


def _write_manifest(target: Path, spec: DNAEmbeddingModel) -> None:
    payload = {
        "schema_version": 1,
        "model_name": spec.name,
        "checkpoint_kind": spec.checkpoint_kind,
        "checkpoint_id": spec.checkpoint_id,
        "checkpoint_url": spec.checkpoint_url,
        "source_repository": spec.source_repository,
        "archive_md5": spec.archive_md5,
        "gated": spec.gated,
    }
    target.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        prefix=f".{_MANIFEST_NAME}.",
        suffix=".part",
        dir=target,
        delete=False,
    ) as output:
        json.dump(payload, output, ensure_ascii=False, indent=2, sort_keys=True)
        output.write("\n")
        temp_path = Path(output.name)
    os.replace(temp_path, target / _MANIFEST_NAME)


def _missing_huggingface_dependency() -> BackendUnavailableError:
    return BackendUnavailableError(
        "Checkpoint download requires huggingface-hub.",
        code="MISSING_NEURAL_DEPENDENCY",
        hint='Install the neural extra with: python -m pip install "dnakit[neural]"',
    )


def _download_huggingface(
    spec: DNAEmbeddingModel,
    target: Path,
    *,
    token: str | None,
    show_progress: bool,
) -> None:
    try:
        huggingface_hub = importlib.import_module("huggingface_hub")
    except ImportError as exc:
        raise _missing_huggingface_dependency() from exc

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
        huggingface_hub.snapshot_download(
            repo_id=spec.checkpoint_id,
            repo_type="model",
            local_dir=str(target),
            allow_patterns=(None if spec.allow_patterns is None else list(spec.allow_patterns)),
            token=token,
            max_workers=4,
        )
    except Exception as exc:
        hint = (
            "Accept the model terms on Hugging Face and provide hf_token."
            if spec.gated
            else "Check network access and the official checkpoint page."
        )
        raise DownloadError(
            "Could not download the model checkpoint from Hugging Face.",
            code="MODEL_CHECKPOINT_DOWNLOAD_FAILED",
            context={"model": spec.name, "checkpoint": spec.checkpoint_id},
            hint=hint,
        ) from exc
    finally:
        if not show_progress and not progress_was_disabled and progress_utils is not None:
            progress_utils.enable_progress_bars()


def _md5(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as source:
        while chunk := source.read(1_048_576):
            digest.update(chunk)
    return digest.hexdigest()


def _janusdna_members(spec: DNAEmbeddingModel) -> dict[str, str]:
    if spec.archive_name is None:
        raise AssertionError("JanusDNA registry entry is missing archive_name.")
    root_name = spec.archive_name.removesuffix(".tar.bz2")
    return {
        f"{root_name}/model_config.json": "model_config.json",
        f"{root_name}/config.json": "config.json",
        f"{root_name}/checkpoints/last.ckpt": "checkpoints/last.ckpt",
    }


def _extract_janusdna_archive(archive: Path, target: Path, spec: DNAEmbeddingModel) -> None:
    selected = _janusdna_members(spec)
    extracted: set[str] = set()
    with tarfile.open(archive, mode="r:bz2") as bundle:
        members = {member.name: member for member in bundle.getmembers()}
        for archive_name, relative_name in selected.items():
            member = members.get(archive_name)
            if member is None or not member.isfile():
                raise DownloadError(
                    "JanusDNA archive is missing a required checkpoint member.",
                    code="INVALID_MODEL_CHECKPOINT_ARCHIVE",
                    context={"member": archive_name},
                )
            source = bundle.extractfile(member)
            if source is None:
                raise DownloadError(
                    "Could not read a required JanusDNA archive member.",
                    code="INVALID_MODEL_CHECKPOINT_ARCHIVE",
                    context={"member": archive_name},
                )
            destination = target / relative_name
            destination.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(
                prefix=f".{destination.name}.",
                suffix=".part",
                dir=destination.parent,
                delete=False,
            ) as output:
                shutil.copyfileobj(source, output)
                temp_path = Path(output.name)
            os.replace(temp_path, destination)
            extracted.add(relative_name)
    if extracted != set(selected.values()):
        raise DownloadError(
            "JanusDNA checkpoint extraction did not complete.",
            code="INVALID_MODEL_CHECKPOINT_ARCHIVE",
        )


def _download_janusdna(
    spec: DNAEmbeddingModel,
    target: Path,
    *,
    show_progress: bool,
) -> None:
    if spec.archive_name is None or spec.archive_md5 is None:
        raise AssertionError("JanusDNA registry entry is incomplete.")
    target.mkdir(parents=True, exist_ok=True)
    archive = target / spec.archive_name
    if archive.exists():
        if not archive.is_file() or _md5(archive) != spec.archive_md5:
            raise DownloadError(
                "Existing JanusDNA archive failed its official MD5 checksum.",
                code="CHECKSUM_MISMATCH",
                context={"path": str(archive)},
                hint="Move the invalid archive aside, then retry the download.",
            )
    else:
        progress_ui: Any | None = None
        callback = None
        if show_progress:
            from rich.progress import Progress

            progress_ui = Progress()
            progress_ui.start()
            task_id = progress_ui.add_task(f"Downloading {spec.display_name}", total=None)

            def update(event: DownloadProgress) -> None:
                progress_ui.update(
                    task_id,
                    completed=event.bytes_completed,
                    total=event.total_bytes,
                )

            callback = update
        try:
            download_file(
                RemoteFile(
                    spec.checkpoint_url,
                    filename=spec.archive_name,
                    expected_md5=spec.archive_md5,
                ),
                archive,
                config=DownloadConfig(
                    timeout=300.0,
                    max_file_bytes=2_000_000_000,
                    max_total_bytes=2_000_000_000,
                ),
                progress=callback,
            )
        finally:
            if progress_ui is not None:
                progress_ui.stop()
    _extract_janusdna_archive(archive, target, spec)


def ensure_model_checkpoint(
    model: str,
    *,
    checkpoint_dir: str | os.PathLike[str] | None = None,
    hf_token: str | None = None,
    show_progress: bool = True,
) -> CheckpointInfo:
    """Return a local checkpoint, downloading it to ``./ckpt`` only when absent.

    The model is stored under ``<checkpoint_dir>/<model-cache-name>``.  A complete
    existing checkpoint is never downloaded again.  Incomplete directories are
    resumed by the provider downloader instead of being reported as complete.
    """

    spec = get_embedding_model(model)
    root = _checkpoint_root(checkpoint_dir)
    target = root / spec.cache_name
    if _checkpoint_ready(spec, target):
        if not (target / _MANIFEST_NAME).is_file():
            _write_manifest(target, spec)
        return CheckpointInfo(
            spec.name,
            str(target),
            spec.checkpoint_url,
            spec.source_repository,
            False,
            spec.gated,
        )

    if spec.checkpoint_kind == "huggingface":
        _download_huggingface(
            spec,
            target,
            token=hf_token,
            show_progress=show_progress,
        )
    elif spec.checkpoint_kind == "dataverse":
        _download_janusdna(spec, target, show_progress=show_progress)
    else:  # pragma: no cover - registry type and tests guard this branch.
        raise AssertionError(f"Unsupported checkpoint kind: {spec.checkpoint_kind}")

    if not _checkpoint_ready(spec, target):
        raise DownloadError(
            "Downloaded checkpoint did not contain the required model files.",
            code="INCOMPLETE_MODEL_CHECKPOINT",
            context={"model": spec.name, "path": str(target)},
        )
    _write_manifest(target, spec)
    return CheckpointInfo(
        spec.name,
        str(target),
        spec.checkpoint_url,
        spec.source_repository,
        True,
        spec.gated,
    )


__all__ = ["CheckpointInfo", "default_checkpoint_root", "ensure_model_checkpoint"]
