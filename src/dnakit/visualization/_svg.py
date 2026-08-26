"""Small deterministic SVG construction helpers."""

from __future__ import annotations

import json
import math
import xml.etree.ElementTree as ET
from collections.abc import Iterable, Mapping

from dnakit.core.issues import Issue
from dnakit.exceptions import ConfigurationError

from .config import SVGTheme
from .results import SVGArtifact

SVG_NAMESPACE = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG_NAMESPACE)


def number(value: int | float) -> str:
    """Stable compact numeric formatting for SVG attributes and labels."""

    if isinstance(value, int):
        return str(value)
    if not math.isfinite(value):
        raise ValueError("SVG numbers must be finite.")
    if math.isclose(value, 0.0, abs_tol=1e-12):
        return "0"
    return format(value, ".8g")


def _attributes(values: Mapping[str, object]) -> dict[str, str]:
    resolved: dict[str, str] = {}
    for name, value in values.items():
        if value is None:
            continue
        attribute_name = name[:-1] if name.endswith("_") else name
        attribute_name = attribute_name.replace("_", "-")
        text = number(value) if isinstance(value, (int, float)) else str(value)
        validate_xml_text(text, name=attribute_name)
        resolved[attribute_name] = text
    return resolved


def validate_xml_text(value: str, *, name: str = "text") -> None:
    """Reject characters forbidden by XML 1.0 before serialization."""

    if not isinstance(value, str):
        raise TypeError(f"{name} must be a string.")
    for character in value:
        codepoint = ord(character)
        if not (
            codepoint in (0x09, 0x0A, 0x0D)
            or 0x20 <= codepoint <= 0xD7FF
            or 0xE000 <= codepoint <= 0xFFFD
            or 0x10000 <= codepoint <= 0x10FFFF
        ):
            raise ValueError(f"{name} contains a character forbidden by XML 1.0.")


class SVGBuilder:
    """Minimal ElementTree wrapper that keeps SVG output deterministic."""

    def __init__(
        self,
        width: int,
        height: int,
        *,
        title: str,
        description: str,
        kind: str,
        theme: SVGTheme,
        metadata: Mapping[str, object],
    ) -> None:
        self.width = width
        self.height = height
        self.kind = kind
        self.metadata = dict(metadata)
        try:
            validate_xml_text(title, name="title")
            validate_xml_text(description, name="description")
        except ValueError as exc:
            raise ConfigurationError(
                "Visualization text contains a character forbidden by XML.",
                code="INVALID_VISUALIZATION_TEXT",
            ) from exc
        self.root = ET.Element(
            f"{{{SVG_NAMESPACE}}}svg",
            {
                "width": str(width),
                "height": str(height),
                "viewBox": f"0 0 {width} {height}",
                "role": "img",
                "aria-labelledby": "dnakit-title dnakit-description",
                "data-kind": kind,
            },
        )
        title_node = ET.SubElement(self.root, f"{{{SVG_NAMESPACE}}}title", {"id": "dnakit-title"})
        title_node.text = title
        description_node = ET.SubElement(
            self.root, f"{{{SVG_NAMESPACE}}}desc", {"id": "dnakit-description"}
        )
        description_node.text = description
        metadata_node = ET.SubElement(self.root, f"{{{SVG_NAMESPACE}}}metadata")
        metadata_text = json.dumps(
            metadata, ensure_ascii=False, sort_keys=True, separators=(",", ":")
        )
        try:
            validate_xml_text(metadata_text, name="metadata")
        except ValueError as exc:
            raise ConfigurationError(
                "Visualization metadata contains a character forbidden by XML.",
                code="INVALID_VISUALIZATION_TEXT",
            ) from exc
        metadata_node.text = metadata_text
        self.rect(0, 0, width, height, fill=theme.background, class_="background")

    def element(self, tag: str, **attributes: object) -> ET.Element:
        return ET.SubElement(
            self.root,
            f"{{{SVG_NAMESPACE}}}{tag}",
            _attributes(attributes),
        )

    def rect(
        self, x: float, y: float, width: float, height: float, **attributes: object
    ) -> ET.Element:
        return self.element("rect", x=x, y=y, width=width, height=height, **attributes)

    def line(self, x1: float, y1: float, x2: float, y2: float, **attributes: object) -> ET.Element:
        return self.element("line", x1=x1, y1=y1, x2=x2, y2=y2, **attributes)

    def text(self, x: float, y: float, value: object, **attributes: object) -> ET.Element:
        node = self.element("text", x=x, y=y, **attributes)
        text = str(value)
        validate_xml_text(text)
        node.text = text
        return node

    def polyline(self, points: Iterable[tuple[float, float]], **attributes: object) -> ET.Element:
        encoded = " ".join(f"{number(x)},{number(y)}" for x, y in points)
        return self.element("polyline", points=encoded, **attributes)

    def artifact(
        self,
        *,
        issues: Iterable[Issue] = (),
        metadata: Mapping[str, object] | None = None,
    ) -> SVGArtifact:
        svg = ET.tostring(self.root, encoding="unicode", short_empty_elements=True)
        element_count = sum(1 for _ in self.root.iter())
        return SVGArtifact(
            self.kind,
            svg,
            self.width,
            self.height,
            element_count,
            metadata=self.metadata if metadata is None else metadata,
            issues=issues,
        )


def parse_hex_color(value: str) -> tuple[int, int, int]:
    """Return RGB for validated 3- or 6-digit hexadecimal colors."""

    if len(value) == 4:
        red, green, blue = (int(symbol * 2, 16) for symbol in value[1:])
        return red, green, blue
    if len(value) == 7:
        red, green, blue = (int(value[index : index + 2], 16) for index in (1, 3, 5))
        return red, green, blue
    raise ValueError("Heatmap interpolation requires #RGB or #RRGGBB colors.")


def interpolate_color(low: str, high: str, fraction: float) -> str:
    low_rgb = parse_hex_color(low)
    high_rgb = parse_hex_color(high)
    bounded = min(1.0, max(0.0, fraction))
    channels = tuple(
        round(start + (end - start) * bounded) for start, end in zip(low_rgb, high_rgb, strict=True)
    )
    return "#" + "".join(f"{channel:02x}" for channel in channels)


__all__ = ["SVG_NAMESPACE", "SVGBuilder", "interpolate_color", "number", "validate_xml_text"]
