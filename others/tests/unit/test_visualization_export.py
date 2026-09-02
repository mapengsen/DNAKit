"""Tests for safe SVG persistence and its structured audit result."""

from __future__ import annotations

import importlib
import os
from collections.abc import Iterator
from pathlib import Path
from typing import cast

import pytest

from dnakit.core import DNARecord, DNASequence
from dnakit.exceptions import ConfigurationError
from dnakit.visualization import (
    ImageExportConfig,
    ImageExportFormat,
    SaveConfig,
    build_html_report,
    plot_sequence,
    save_html_report,
    save_image,
    save_svg,
)
from dnakit.visualization.results import SVGArtifact


def test_save_svg_writes_exact_payload_and_structured_artifact(tmp_path: Path) -> None:
    artifact = plot_sequence(DNASequence("ACGT"))
    target = tmp_path / "sequence.svg"

    result = save_svg(artifact, target)

    assert target.read_bytes() == artifact.svg.encode("utf-8")
    assert result.source_sha256 == artifact.sha256
    assert result.target_artifact.sha256 == artifact.sha256
    assert result.target_artifact.byte_size == len(artifact.svg.encode("utf-8"))
    assert result.target_artifact.media_type == "image/svg+xml"
    assert result.target_artifact.schema_version == "dnakit-visualization-v1"
    assert result.target_artifact.relative_path == os.path.relpath(target.resolve(), Path.cwd())
    assert result.overwritten is False


def test_save_svg_refuses_overwrite_by_default_without_changing_file(tmp_path: Path) -> None:
    target = tmp_path / "sequence.svg"
    target.write_text("original", encoding="utf-8")

    with pytest.raises(FileExistsError):
        save_svg(plot_sequence(DNASequence("AC")), target)

    assert target.read_text(encoding="utf-8") == "original"
    assert not list(tmp_path.glob(f".{target.name}.*.tmp"))


def test_save_svg_can_explicitly_overwrite(tmp_path: Path) -> None:
    target = tmp_path / "sequence.svg"
    target.write_text("old", encoding="utf-8")
    artifact = plot_sequence(DNASequence("TGCA"))

    result = save_svg(artifact, target, config=SaveConfig(overwrite=True))

    assert target.read_text(encoding="utf-8") == artifact.svg
    assert result.overwritten is True


def test_save_svg_requires_explicit_parent_creation_and_svg_extension(tmp_path: Path) -> None:
    artifact = plot_sequence(DNASequence("A"))
    nested = tmp_path / "nested" / "sequence.svg"

    with pytest.raises(FileNotFoundError):
        save_svg(artifact, nested)
    created = save_svg(artifact, nested, config=SaveConfig(create_parents=True))
    assert nested.exists()
    assert created.overwritten is False

    with pytest.raises(ConfigurationError) as extension_error:
        save_svg(artifact, tmp_path / "sequence.png")
    assert extension_error.value.code == "INVALID_VISUALIZATION_FORMAT"


def test_save_svg_atomic_link_failure_leaves_no_target_or_temporary_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    export_module = importlib.import_module("dnakit.visualization.export")
    target = tmp_path / "sequence.svg"

    def fail_link(_source: Path, _target: Path) -> None:
        raise OSError("simulated link failure")

    monkeypatch.setattr(export_module.os, "link", fail_link)
    with pytest.raises(OSError, match="simulated link failure"):
        save_svg(plot_sequence(DNASequence("AC")), target)

    assert not target.exists()
    assert not list(tmp_path.glob(f".{target.name}.*.tmp"))


@pytest.mark.parametrize(
    ("suffix", "magic", "expected_format"),
    [
        ("png", b"\x89PNG", "png"),
        ("jpg", b"\xff\xd8\xff", "jpg"),
        ("tiff", b"II*\x00", "tiff"),
        ("pdf", b"%PDF", "pdf"),
    ],
)
def test_save_image_optional_backend_exports_600_dpi_formats(
    tmp_path: Path, suffix: str, magic: bytes, expected_format: str
) -> None:
    pytest.importorskip("cairosvg")
    if suffix in {"png", "jpg", "tiff"}:
        pytest.importorskip("PIL")
    target = tmp_path / f"sequence.{suffix}"
    result = save_image(
        plot_sequence(DNASequence("ACGT")),
        target,
        config=ImageExportConfig(dpi=600, width=320),
    )

    assert target.read_bytes().startswith(magic)
    assert result.format == expected_format
    assert (
        result.target_artifact.media_type
        == {
            "png": "image/png",
            "jpg": "image/jpeg",
            "tiff": "image/tiff",
            "pdf": "application/pdf",
        }[expected_format]
    )
    assert result.target_artifact.byte_size == target.stat().st_size


def test_save_image_defaults_extensionless_target_to_png(tmp_path: Path) -> None:
    pytest.importorskip("cairosvg")
    pytest.importorskip("PIL")
    target = tmp_path / "sequence"

    result = save_image(plot_sequence(DNASequence("ACGT")), target)

    resolved_target = target.with_suffix(".png")
    assert resolved_target.read_bytes().startswith(b"\x89PNG")
    assert result.format == "png"
    assert result.target_artifact.media_type == "image/png"
    assert result.target_artifact.relative_path == os.path.relpath(
        resolved_target.resolve(), Path.cwd()
    )


def test_save_image_applies_target_dpi_once(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cairosvg = pytest.importorskip("cairosvg")
    image_module = pytest.importorskip("PIL.Image")
    original_svg2png = cairosvg.svg2png
    rendered_dpi: list[int] = []

    def capture_svg2png(**parameters: object) -> bytes:
        rendered_dpi.append(parameters["dpi"])  # type: ignore[arg-type]
        return cast(bytes, original_svg2png(**parameters))

    monkeypatch.setattr(cairosvg, "svg2png", capture_svg2png)
    target = tmp_path / "sequence.png"
    save_image(
        plot_sequence(DNASequence("ACGT")),
        target,
        config=ImageExportConfig(dpi=600, width=320),
    )

    assert rendered_dpi == [96]
    with image_module.open(target) as image:
        assert image.size == (320, 320)
        assert image.info["dpi"] == pytest.approx((600, 600), abs=0.01)


def test_image_export_rejects_non_square_explicit_dimensions() -> None:
    with pytest.raises(ConfigurationError) as error:
        ImageExportConfig(width=320, height=240)

    assert error.value.code == "INVALID_VISUALIZATION_CONFIG"


@pytest.mark.parametrize(
    ("image_type", "expected_suffix", "expected_format", "magic"),
    [
        ("png", ".png", "png", b"\x89PNG"),
        ("svg", ".svg", "svg", b"<svg"),
        ("jpg", ".jpg", "jpg", b"\xff\xd8\xff"),
    ],
)
def test_save_image_selects_png_svg_or_jpg_and_adds_extension(
    tmp_path: Path,
    image_type: ImageExportFormat,
    expected_suffix: str,
    expected_format: str,
    magic: bytes,
) -> None:
    if image_type != "svg":
        pytest.importorskip("cairosvg")
        pytest.importorskip("PIL")
    target = tmp_path / f"sequence-{image_type}"

    result = save_image(
        plot_sequence(DNASequence("ACGT")),
        target,
        image_type=image_type,
    )

    resolved_target = target.with_suffix(expected_suffix)
    assert resolved_target.read_bytes().startswith(magic)
    assert result.format == expected_format


def test_save_image_rejects_type_extension_mismatch_and_invalid_type(tmp_path: Path) -> None:
    artifact = plot_sequence(DNASequence("ACGT"))

    with pytest.raises(ConfigurationError) as mismatch_error:
        save_image(artifact, tmp_path / "sequence.svg", image_type="jpg")
    assert mismatch_error.value.code == "INVALID_VISUALIZATION_FORMAT"

    with pytest.raises(ConfigurationError) as invalid_error:
        save_image(artifact, tmp_path / "sequence", image_type="gif")  # type: ignore[arg-type]
    assert invalid_error.value.code == "INVALID_VISUALIZATION_FORMAT"


def test_save_image_rejects_transparent_jpg(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError) as error:
        save_image(
            plot_sequence(DNASequence("ACGT")),
            tmp_path / "sequence",
            image_type="jpg",
            config=ImageExportConfig(transparent=True),
        )
    assert error.value.code == "INVALID_VISUALIZATION_CONFIG"


def test_save_image_rejects_pixel_limit_before_backend_conversion(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError) as error:
        save_image(
            plot_sequence(DNASequence("A")),
            tmp_path / "large.png",
            config=ImageExportConfig(width=20_000, height=20_000, max_output_pixels=100),
        )
    assert error.value.code == "VISUALIZATION_PIXEL_LIMIT"


def test_svg_artifact_validates_xml_and_actual_element_count() -> None:
    with pytest.raises(ConfigurationError) as count_error:
        SVGArtifact("test", '<svg xmlns="http://www.w3.org/2000/svg"><rect /></svg>', 10, 10, 1)
    assert count_error.value.code == "SVG_ELEMENT_COUNT_MISMATCH"

    with pytest.raises(ConfigurationError) as xml_error:
        SVGArtifact("test", "<svg><rect></svg>", 10, 10, 2)
    assert xml_error.value.code == "INVALID_SVG_INPUT"

    with pytest.raises(ConfigurationError) as entity_error:
        SVGArtifact("test", '<!DOCTYPE svg><svg xmlns="http://www.w3.org/2000/svg" />', 10, 10, 1)
    assert entity_error.value.code == "UNSAFE_SVG_INPUT"


def test_html_report_is_self_contained_searchable_escaped_and_atomic(tmp_path: Path) -> None:
    report = build_html_report(
        [DNARecord(DNASequence("ACGT"), "<unsafe>", description="A & B")],
        results={"metric": {"value": 0.5}},
    )
    assert "https://" not in report.html
    assert "<unsafe>" not in report.html
    assert "&lt;unsafe&gt;" in report.html
    assert 'id="filter"' in report.html
    assert "<details>" in report.html
    target = tmp_path / "report.html"

    result = save_html_report(report, target)

    assert target.read_text(encoding="utf-8") == report.html
    assert result.format == "html"
    assert result.target_artifact.sha256 == report.sha256


def test_html_report_bounds_input_and_result_payload() -> None:
    def records() -> Iterator[DNARecord]:
        for index in range(3):
            yield DNARecord(DNASequence("A"), str(index))

    with pytest.raises(ConfigurationError) as records_error:
        build_html_report(records(), max_records=2)
    assert records_error.value.code == "HTML_REPORT_RECORD_LIMIT"
    with pytest.raises(ConfigurationError) as result_error:
        build_html_report([], results={"large": "A" * 20}, max_result_bytes=10)
    assert result_error.value.code == "HTML_REPORT_RESULT_LIMIT"


@pytest.mark.parametrize(
    ("argument", "value"),
    [
        ("max_total_sequence_symbols", 0),
        ("max_total_sequence_symbols", True),
        ("max_total_record_text_characters", 1.5),
        ("max_output_bytes", -1),
    ],
)
def test_html_report_rejects_invalid_resource_limits(argument: str, value: object) -> None:
    with pytest.raises(ConfigurationError, match=f"{argument} must be a positive integer"):
        build_html_report([], **{argument: value})  # type: ignore[arg-type]


def test_html_report_rejects_one_oversized_sequence_before_rendering() -> None:
    record = DNARecord(DNASequence("A" * 11), "long")

    with pytest.raises(ConfigurationError) as error:
        build_html_report([record], max_total_sequence_symbols=10)

    assert error.value.code == "HTML_REPORT_SEQUENCE_LIMIT"
    assert error.value.context["total_sequence_symbols"] == 11


def test_html_report_limits_cumulative_id_and_description_characters() -> None:
    records = [
        DNARecord(DNASequence("A"), "abc", description="de"),
        DNARecord(DNASequence("C"), "f", description="gh"),
    ]

    with pytest.raises(ConfigurationError) as error:
        build_html_report(records, max_total_record_text_characters=7)

    assert error.value.code == "HTML_REPORT_RECORD_TEXT_LIMIT"
    assert error.value.context["total_record_text_characters"] == 8


def test_html_report_max_output_bytes_is_exact_utf8_boundary() -> None:
    records = [DNARecord(DNASequence("ACGT"), "id<&", description="描述 & value")]
    unbounded = build_html_report(records, results={"metric": {"value": "<&"}})
    exact_size = len(unbounded.html.encode("utf-8"))

    at_limit = build_html_report(
        records,
        results={"metric": {"value": "<&"}},
        max_output_bytes=exact_size,
    )
    assert at_limit.html == unbounded.html
    with pytest.raises(ConfigurationError) as error:
        build_html_report(
            records,
            results={"metric": {"value": "<&"}},
            max_output_bytes=exact_size - 1,
        )
    assert error.value.code == "HTML_REPORT_OUTPUT_LIMIT"


def test_html_report_rejects_deep_result_before_json_encoding() -> None:
    value: object = "leaf"
    for _ in range(100):
        value = [value]

    with pytest.raises(ConfigurationError) as error:
        build_html_report([], results={"deep": value}, max_result_depth=32)
    assert error.value.code == "HTML_REPORT_RESULT_STRUCTURE_LIMIT"


def test_html_report_accepts_shared_dag_and_rejects_cycles() -> None:
    shared = {"x": 1}
    report = build_html_report([], results={"ok": {"left": shared, "right": shared}})
    assert "&quot;left&quot;" in report.html

    cycle: list[object] = []
    cycle.append(cycle)
    with pytest.raises(ConfigurationError) as error:
        build_html_report([], results={"cycle": cycle})
    assert error.value.code == "HTML_REPORT_RESULT_STRUCTURE_LIMIT"
