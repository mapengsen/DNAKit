"""Atomic persistence for standalone SVG visualization artifacts."""

from __future__ import annotations

import hashlib
import importlib
import os
import tempfile
from contextlib import suppress
from datetime import datetime, timezone
from io import BytesIO
from pathlib import Path
from types import TracebackType
from typing import BinaryIO, Protocol, TypeAlias, cast

from dnakit.core.provenance import ArtifactRef
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

from .config import ImageExportConfig, SaveConfig
from .results import HTMLReportArtifact, SVGArtifact, VisualizationSaveResult

PathLike: TypeAlias = str | os.PathLike[str]


class _CairoSVGModule(Protocol):
    def svg2pdf(
        self,
        *,
        bytestring: bytes,
        dpi: int,
        output_width: int,
        output_height: int,
    ) -> bytes: ...

    def svg2png(
        self,
        *,
        bytestring: bytes,
        dpi: int,
        output_width: int,
        output_height: int,
    ) -> bytes: ...


class _PILImage(Protocol):
    def __enter__(self) -> _PILImage: ...

    def __exit__(
        self,
        exception_type: type[BaseException] | None,
        exception: BaseException | None,
        traceback: TracebackType | None,
    ) -> bool | None: ...

    def save(self, stream: BinaryIO, *, format: str, **parameters: object) -> None: ...


class _PILImageModule(Protocol):
    def open(self, stream: BinaryIO) -> _PILImage: ...


def _artifact_ref(
    path: Path,
    *,
    sha256: str,
    byte_size: int,
    media_type: str,
    schema_version: str,
) -> ArtifactRef:
    stat = path.stat()
    return ArtifactRef(
        relative_path=os.path.relpath(path.resolve(), Path.cwd()),
        media_type=media_type,
        schema_version=schema_version,
        sha256=sha256,
        byte_size=byte_size,
        created_at=datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat(),
    )


def _atomic_write(path: Path, payload: bytes, *, overwrite: bool, create_parents: bool) -> bool:
    parent = path.parent
    if not parent.exists():
        if create_parents:
            parent.mkdir(parents=True, exist_ok=True)
        else:
            raise FileNotFoundError(f"Output parent directory does not exist: {parent}")
    if not parent.is_dir():
        raise NotADirectoryError(f"Output parent is not a directory: {parent}")
    existed = path.exists()
    if existed and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing output: {path}")
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=parent
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        if overwrite:
            os.replace(temporary_path, path)
        else:
            os.link(temporary_path, path)
            temporary_path.unlink()
    except BaseException:
        with suppress(OSError):
            os.close(descriptor)
        temporary_path.unlink(missing_ok=True)
        raise
    return existed


def save_svg(
    artifact: SVGArtifact,
    target: PathLike,
    *,
    config: SaveConfig | None = None,
) -> VisualizationSaveResult:
    """Atomically save exact SVG bytes, refusing overwrite by default."""

    if not isinstance(artifact, SVGArtifact):
        raise TypeError("artifact must be SVGArtifact.")
    if not isinstance(target, (str, os.PathLike)):
        raise TypeError("target must be a string or path-like object.")
    resolved = SaveConfig() if config is None else config
    if not isinstance(resolved, SaveConfig):
        raise TypeError("config must be SaveConfig or None.")
    path = Path(target)
    if path.suffix.lower() != ".svg":
        raise ConfigurationError(
            "SVG output path must end with .svg.",
            code="INVALID_VISUALIZATION_FORMAT",
            context={"target": str(path)},
        )
    payload = artifact.svg.encode("utf-8")
    existed = _atomic_write(
        path,
        payload,
        overwrite=resolved.overwrite,
        create_parents=resolved.create_parents,
    )
    target_artifact = _artifact_ref(
        path,
        sha256=artifact.sha256,
        byte_size=len(payload),
        media_type="image/svg+xml",
        schema_version="dnakit-visualization-v1",
    )
    return VisualizationSaveResult(
        artifact.sha256,
        target_artifact,
        overwritten=existed,
        provenance=artifact.provenance,
        issues=artifact.issues,
    )


def _target_dimensions(artifact: SVGArtifact, config: ImageExportConfig) -> tuple[int, int]:
    if config.width is None and config.height is None:
        scale = config.dpi / 96.0
        width = max(1, round(artifact.width * scale))
        height = max(1, round(artifact.height * scale))
    elif config.width is None:
        assert config.height is not None
        height = config.height
        width = max(1, round(artifact.width * height / artifact.height))
    elif config.height is None:
        width = config.width
        height = max(1, round(artifact.height * width / artifact.width))
    else:
        width, height = config.width, config.height
    if width * height > config.max_output_pixels:
        raise ConfigurationError(
            "Rendered image exceeds max_output_pixels.",
            code="VISUALIZATION_PIXEL_LIMIT",
            context={
                "width": width,
                "height": height,
                "max_output_pixels": config.max_output_pixels,
            },
        )
    return width, height


def _svg_payload(artifact: SVGArtifact, *, transparent: bool) -> bytes:
    if not transparent:
        return artifact.svg.encode("utf-8")
    import xml.etree.ElementTree as ET

    root = ET.fromstring(artifact.svg)
    for child in tuple(root):
        if child.tag.endswith("rect") and child.attrib.get("class") == "background":
            root.remove(child)
            break
    return cast(bytes, ET.tostring(root, encoding="utf-8", xml_declaration=False))


def _convert_image(
    artifact: SVGArtifact, output_format: str, config: ImageExportConfig
) -> tuple[bytes, int, int]:
    width, height = _target_dimensions(artifact, config)
    try:
        cairosvg = cast(_CairoSVGModule, importlib.import_module("cairosvg"))
    except ImportError as exc:
        raise BackendUnavailableError(
            "PNG, TIFF, and PDF export requires the optional visualization backend.",
            code="VISUALIZATION_BACKEND_UNAVAILABLE",
            hint="Install DNAKit with the 'viz' extra.",
        ) from exc
    source = _svg_payload(artifact, transparent=config.transparent)
    if output_format == "pdf":
        payload = cairosvg.svg2pdf(
            bytestring=source,
            dpi=config.dpi,
            output_width=width,
            output_height=height,
        )
        return payload, width, height
    png = cairosvg.svg2png(
        bytestring=source,
        dpi=config.dpi,
        output_width=width,
        output_height=height,
    )
    try:
        image_module = cast(_PILImageModule, importlib.import_module("PIL.Image"))
    except ImportError as exc:
        raise BackendUnavailableError(
            "PNG and TIFF export require Pillow from the optional visualization extra.",
            code="VISUALIZATION_BACKEND_UNAVAILABLE",
            hint="Install DNAKit with the 'viz' extra.",
        ) from exc
    output = BytesIO()
    with image_module.open(BytesIO(png)) as image:
        if output_format == "png":
            image.save(output, format="PNG", dpi=(config.dpi, config.dpi))
        else:
            image.save(output, format="TIFF", compression="tiff_lzw", dpi=(config.dpi, config.dpi))
    return output.getvalue(), width, height


def save_image(
    artifact: SVGArtifact,
    target: PathLike,
    *,
    config: ImageExportConfig | None = None,
) -> VisualizationSaveResult:
    """Atomically export an SVG artifact as SVG, PNG, TIFF, or PDF."""

    if not isinstance(artifact, SVGArtifact):
        raise TypeError("artifact must be SVGArtifact.")
    if not isinstance(target, (str, os.PathLike)):
        raise TypeError("target must be a string or path-like object.")
    resolved = ImageExportConfig() if config is None else config
    if not isinstance(resolved, ImageExportConfig):
        raise TypeError("config must be ImageExportConfig or None.")
    path = Path(target)
    suffix = path.suffix.lower()
    format_by_suffix = {
        ".svg": "svg",
        ".png": "png",
        ".tif": "tiff",
        ".tiff": "tiff",
        ".pdf": "pdf",
    }
    output_format = format_by_suffix.get(suffix)
    if output_format is None:
        raise ConfigurationError(
            "Visualization output must end with .svg, .png, .tif, .tiff, or .pdf.",
            code="INVALID_VISUALIZATION_FORMAT",
        )
    if output_format == "svg":
        return save_svg(
            artifact,
            path,
            config=SaveConfig(overwrite=resolved.overwrite, create_parents=resolved.create_parents),
        )
    payload, _width, _height = _convert_image(artifact, output_format, resolved)
    existed = _atomic_write(
        path,
        payload,
        overwrite=resolved.overwrite,
        create_parents=resolved.create_parents,
    )
    media_type = {"png": "image/png", "tiff": "image/tiff", "pdf": "application/pdf"}[output_format]
    digest = hashlib.sha256(payload).hexdigest()
    target_artifact = _artifact_ref(
        path,
        sha256=digest,
        byte_size=len(payload),
        media_type=media_type,
        schema_version="dnakit-visualization-export-v1",
    )
    return VisualizationSaveResult(
        artifact.sha256,
        target_artifact,
        format=output_format,
        overwritten=existed,
        provenance=artifact.provenance,
        issues=artifact.issues,
    )


def save_html_report(
    artifact: HTMLReportArtifact,
    target: PathLike,
    *,
    config: SaveConfig | None = None,
) -> VisualizationSaveResult:
    """Atomically persist a self-contained DNAKit HTML report."""

    if not isinstance(artifact, HTMLReportArtifact):
        raise TypeError("artifact must be HTMLReportArtifact.")
    resolved = SaveConfig() if config is None else config
    if not isinstance(resolved, SaveConfig):
        raise TypeError("config must be SaveConfig or None.")
    path = Path(target)
    if path.suffix.lower() not in {".html", ".htm"}:
        raise ConfigurationError(
            "HTML report path must end with .html or .htm.",
            code="INVALID_VISUALIZATION_FORMAT",
        )
    payload = artifact.html.encode("utf-8")
    existed = _atomic_write(
        path,
        payload,
        overwrite=resolved.overwrite,
        create_parents=resolved.create_parents,
    )
    target_artifact = _artifact_ref(
        path,
        sha256=artifact.sha256,
        byte_size=len(payload),
        media_type="text/html",
        schema_version="dnakit-html-report-v1",
    )
    return VisualizationSaveResult(
        artifact.sha256,
        target_artifact,
        format="html",
        overwritten=existed,
        provenance=artifact.provenance,
        issues=artifact.issues,
    )


__all__ = ["PathLike", "save_html_report", "save_image", "save_svg"]
