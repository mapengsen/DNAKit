"""Immutable SVG and persistence result objects."""

from __future__ import annotations

import hashlib
import xml.etree.ElementTree as ET
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any, cast

from dnakit.core._json import FrozenDict, freeze_mapping, to_json_compatible
from dnakit.core.issues import Issue
from dnakit.core.provenance import ArtifactRef, Provenance
from dnakit.exceptions import ConfigurationError

MAX_SVG_INPUT_BYTES = 50_000_000
MAX_SVG_ELEMENTS = 1_000_000


@dataclass(frozen=True, init=False)
class SVGArtifact:
    """A deterministic, standalone SVG document and its audit metadata."""

    kind: str
    svg: str
    width: int
    height: int
    element_count: int
    metadata: FrozenDict
    provenance: Provenance
    issues: tuple[Issue, ...]

    def __init__(
        self,
        kind: str,
        svg: str,
        width: int,
        height: int,
        element_count: int,
        *,
        metadata: Mapping[str, object] | None = None,
        provenance: Provenance | None = None,
        issues: Iterable[Issue] = (),
    ) -> None:
        for name, text_value in (("kind", kind), ("svg", svg)):
            if not isinstance(text_value, str) or not text_value.strip():
                raise ConfigurationError(f"SVGArtifact {name} must be non-empty.")
        byte_count = len(svg.encode("utf-8"))
        if byte_count > MAX_SVG_INPUT_BYTES:
            raise ConfigurationError(
                "SVGArtifact exceeds the SVG byte limit.",
                code="SVG_INPUT_SIZE_LIMIT",
                context={"max_svg_input_bytes": MAX_SVG_INPUT_BYTES},
            )
        if "<!DOCTYPE" in svg.upper() or "<!ENTITY" in svg.upper():
            raise ConfigurationError(
                "SVGArtifact must not contain document type or entity declarations.",
                code="UNSAFE_SVG_INPUT",
            )
        try:
            root = ET.fromstring(svg)
        except (ET.ParseError, RecursionError) as exc:
            raise ConfigurationError(
                "SVGArtifact svg must be valid XML with one SVG root.",
                code="INVALID_SVG_INPUT",
            ) from exc
        if root.tag.rsplit("}", 1)[-1] != "svg":
            raise ConfigurationError(
                "SVGArtifact root element must be SVG.", code="INVALID_SVG_INPUT"
            )
        actual_element_count = 0
        for actual_element_count, _element in enumerate(root.iter(), start=1):
            if actual_element_count > MAX_SVG_ELEMENTS:
                raise ConfigurationError(
                    "SVGArtifact exceeds the SVG element limit.",
                    code="SVG_ELEMENT_LIMIT",
                    context={"max_svg_elements": MAX_SVG_ELEMENTS},
                )
        for name, integer_value in (
            ("width", width),
            ("height", height),
            ("element_count", element_count),
        ):
            if (
                isinstance(integer_value, bool)
                or not isinstance(integer_value, int)
                or integer_value <= 0
            ):
                raise ConfigurationError(f"SVGArtifact {name} must be a positive integer.")
        if element_count != actual_element_count:
            raise ConfigurationError(
                "SVGArtifact element_count does not match the parsed SVG.",
                code="SVG_ELEMENT_COUNT_MISMATCH",
                context={"declared": element_count, "actual": actual_element_count},
            )
        if provenance is not None and not isinstance(provenance, Provenance):
            raise ConfigurationError("SVGArtifact provenance must be Provenance or None.")
        issue_tuple = tuple(issues)
        if any(not isinstance(issue, Issue) for issue in issue_tuple):
            raise ConfigurationError("SVGArtifact issues must contain Issue objects.")
        object.__setattr__(self, "kind", kind)
        object.__setattr__(self, "svg", svg)
        object.__setattr__(self, "width", width)
        object.__setattr__(self, "height", height)
        object.__setattr__(self, "element_count", element_count)
        object.__setattr__(self, "metadata", freeze_mapping(metadata))
        object.__setattr__(
            self, "provenance", provenance if provenance is not None else Provenance()
        )
        object.__setattr__(self, "issues", issue_tuple)

    @property
    def sha256(self) -> str:
        """SHA-256 of the exact UTF-8 SVG payload."""

        return hashlib.sha256(self.svg.encode("utf-8")).hexdigest()

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True, init=False)
class VisualizationSaveResult:
    """Audit record for a successfully persisted visualization artifact."""

    format: str
    source_sha256: str
    target_artifact: ArtifactRef
    overwritten: bool
    provenance: Provenance
    issues: tuple[Issue, ...]

    def __init__(
        self,
        source_sha256: str,
        target_artifact: ArtifactRef,
        *,
        format: str = "svg",
        overwritten: bool,
        provenance: Provenance | None = None,
        issues: Iterable[Issue] = (),
    ) -> None:
        if (
            not isinstance(source_sha256, str)
            or len(source_sha256) != 64
            or any(character not in "0123456789abcdefABCDEF" for character in source_sha256)
        ):
            raise ConfigurationError("source_sha256 must contain 64 hexadecimal digits.")
        if not isinstance(target_artifact, ArtifactRef):
            raise ConfigurationError("target_artifact must be ArtifactRef.")
        if format not in {"svg", "png", "jpg", "tiff", "pdf", "html"}:
            raise ConfigurationError("Unknown visualization output format.")
        if not isinstance(overwritten, bool):
            raise ConfigurationError("overwritten must be a boolean.")
        if provenance is not None and not isinstance(provenance, Provenance):
            raise ConfigurationError("provenance must be Provenance or None.")
        issue_tuple = tuple(issues)
        if any(not isinstance(issue, Issue) for issue in issue_tuple):
            raise ConfigurationError("issues must contain Issue objects.")
        object.__setattr__(self, "format", format)
        object.__setattr__(self, "source_sha256", source_sha256.lower())
        object.__setattr__(self, "target_artifact", target_artifact)
        object.__setattr__(self, "overwritten", overwritten)
        object.__setattr__(
            self, "provenance", provenance if provenance is not None else Provenance()
        )
        object.__setattr__(self, "issues", issue_tuple)

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True, init=False)
class HTMLReportArtifact:
    """A bounded self-contained HTML report with no external resources."""

    html: str
    record_count: int
    result_names: tuple[str, ...]
    provenance: Provenance
    issues: tuple[Issue, ...]

    def __init__(
        self,
        html: str,
        record_count: int,
        result_names: Iterable[str],
        *,
        provenance: Provenance | None = None,
        issues: Iterable[Issue] = (),
    ) -> None:
        if not isinstance(html, str) or not html.startswith("<!doctype html>"):
            raise ConfigurationError("HTMLReportArtifact html must be a complete HTML document.")
        if isinstance(record_count, bool) or not isinstance(record_count, int) or record_count < 0:
            raise ConfigurationError("record_count must be a non-negative integer.")
        names = tuple(result_names)
        if any(not isinstance(name, str) or not name.strip() for name in names):
            raise ConfigurationError("result_names must contain non-empty strings.")
        if len(set(names)) != len(names):
            raise ConfigurationError("result_names must be unique.")
        resolved_provenance = Provenance() if provenance is None else provenance
        if not isinstance(resolved_provenance, Provenance):
            raise ConfigurationError("provenance must be Provenance or None.")
        issue_tuple = tuple(issues)
        if any(not isinstance(issue, Issue) for issue in issue_tuple):
            raise ConfigurationError("issues must contain Issue objects.")
        object.__setattr__(self, "html", html)
        object.__setattr__(self, "record_count", record_count)
        object.__setattr__(self, "result_names", names)
        object.__setattr__(self, "provenance", resolved_provenance)
        object.__setattr__(self, "issues", issue_tuple)

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.html.encode("utf-8")).hexdigest()

    def to_dict(self, *, include_html: bool = False) -> dict[str, Any]:
        payload: dict[str, object] = {
            "record_count": self.record_count,
            "result_names": self.result_names,
            "sha256": self.sha256,
            "provenance": self.provenance,
            "issues": self.issues,
        }
        if include_html:
            payload["html"] = self.html
        return cast(dict[str, Any], to_json_compatible(payload))


__all__ = [
    "MAX_SVG_ELEMENTS",
    "MAX_SVG_INPUT_BYTES",
    "HTMLReportArtifact",
    "SVGArtifact",
    "VisualizationSaveResult",
]
