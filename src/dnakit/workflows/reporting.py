"""Render user-supplied WES report content without generating medical conclusions."""

from __future__ import annotations

import hashlib
import json
import re
import tempfile
from pathlib import Path
from typing import Any
from xml.sax.saxutils import escape

from dnakit.exceptions import BackendUnavailableError, ConfigurationError


def write_wes_report(
    input_report: str | Path,
    output_dir: str | Path,
    *,
    language: str = "en",
    title: str | None = None,
    logo_left: str | Path | None = None,
    logo_right: str | Path | None = None,
    font_path: str | Path | None = None,
    overwrite: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Render supplied Markdown sections and tables as a PDF, with a source hash.

    language selects English or Spanish report labels; input content is preserved
    without translation or new clinical interpretations. Optional logos and a TTF
    font are supported. Requires the reportlab rendering backend (dnakit[reports]).
    This implementation was written for DNAKit; no proprietary report code is used.
    """
    if language not in {"en", "es"}:
        raise ConfigurationError("language must be en or es.")
    source = Path(input_report).expanduser().resolve(strict=True)
    if source.stat().st_size > 5_000_000:
        raise ConfigurationError("Report Markdown exceeds 5 MB.")
    content = source.read_bytes()
    markdown = content.decode("utf-8")
    directory = Path(output_dir).expanduser().resolve()
    pdf_path, manifest_path = directory / "wes_report.pdf", directory / "wes_report.manifest.json"
    if not overwrite and (pdf_path.exists() or manifest_path.exists()):
        raise ConfigurationError(
            "Report output already exists; select a new output directory or overwrite=True."
        )
    try:
        from reportlab.lib import colors
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.styles import getSampleStyleSheet
        from reportlab.pdfbase import pdfmetrics
        from reportlab.pdfbase.ttfonts import TTFont
        from reportlab.platypus import (
            Image,
            Paragraph,
            SimpleDocTemplate,
            Spacer,
            Table,
            TableStyle,
        )
    except ImportError as exc:
        raise BackendUnavailableError(
            'Install the PDF renderer with python -m pip install "dnakit[reports]".'
        ) from exc
    styles = getSampleStyleSheet()
    if font_path is not None:
        font = Path(font_path).expanduser().resolve(strict=True)
        font_name = "DNAKitReport" + hashlib.sha256(str(font).encode()).hexdigest()[:12]
        pdfmetrics.registerFont(TTFont(font_name, str(font)))
        for style in styles.byName.values():
            style.fontName = font_name
    story = []
    logos = []
    for logo in (logo_left, logo_right):
        if logo is not None:
            picture = Image(str(Path(logo).expanduser().resolve(strict=True)))
            ratio = min(150 / picture.imageWidth, 65 / picture.imageHeight)
            picture.drawWidth, picture.drawHeight = (
                picture.imageWidth * ratio,
                picture.imageHeight * ratio,
            )
            logos.append(picture)
    if logos:
        story.append(Table([logos]))
    report_title = title or ("WES report" if language == "en" else "Informe WES")
    story.extend([Paragraph(escape(report_title), styles["Title"]), Spacer(1, 14)])
    lines = markdown.splitlines()
    index = 0
    while index < len(lines):
        line = lines[index]
        if line.strip().startswith("|"):
            rows = []
            while index < len(lines) and lines[index].strip().startswith("|"):
                raw = lines[index].strip().strip("|").split("|")
                if not all(re.fullmatch(r"\s*:?-+:?\s*", value) for value in raw):
                    rows.append(
                        [Paragraph(escape(value.strip()), styles["BodyText"]) for value in raw]
                    )
                index += 1
            if rows:
                columns = max(map(len, rows))
                rows = [row + [""] * (columns - len(row)) for row in rows]
                table = Table(rows, colWidths=[(A4[0] - 88) / columns] * columns, repeatRows=1)
                table.setStyle(
                    TableStyle(
                        [
                            ("GRID", (0, 0), (-1, -1), 0.3, colors.grey),
                            ("VALIGN", (0, 0), (-1, -1), "TOP"),
                        ]
                    )
                )
                story.extend([table, Spacer(1, 8)])
            continue
        heading = re.match(r"^(#{1,6})\s+(.*)", line)
        if heading:
            story.append(
                Paragraph(escape(heading[2]), styles["Heading" + str(min(len(heading[1]), 3))])
            )
        elif line.strip():
            story.append(Paragraph(escape(line), styles["BodyText"]))
        else:
            story.append(Spacer(1, 6))
        index += 1
    directory.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".dnakit-report-", dir=directory) as temporary_dir:
        temporary = Path(temporary_dir) / "report.pdf"
        SimpleDocTemplate(
            str(temporary), pagesize=A4, leftMargin=44, rightMargin=44, title=report_title
        ).build(story)
        manifest = {
            "implementation": "dnakit.workflows.write_wes_report",
            "language": language,
            "source": str(source),
            "source_sha256": hashlib.sha256(content).hexdigest(),
            "output": str(pdf_path),
            "output_sha256": hashlib.sha256(temporary.read_bytes()).hexdigest(),
        }
        temporary.replace(pdf_path)
        staged_manifest = Path(temporary_dir) / "manifest.json"
        staged_manifest.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
        staged_manifest.replace(manifest_path)
    return {**manifest, "manifest": str(manifest_path)}
