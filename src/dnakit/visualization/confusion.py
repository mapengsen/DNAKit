"""Confusion matrices from ordinary prediction CSV files."""

from __future__ import annotations

import csv
from collections.abc import Sequence
from pathlib import Path

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import ConfigurationError, InputFormatError

from ._svg import SVGBuilder, interpolate_color, number, validate_xml_text
from .config import ImageExportConfig, SVGTheme
from .export import save_image
from .results import SVGArtifact


def _read_labels(
    source: str | Path, true_column: str, predicted_column: str
) -> tuple[list[str], list[str]]:
    path = Path(source)
    actual: list[str] = []
    predicted: list[str] = []
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, strict=True)
        headers = reader.fieldnames
        if (
            not headers
            or len(set(headers)) != len(headers)
            or not {true_column, predicted_column} <= set(headers)
        ):
            raise InputFormatError(
                f"Prediction CSV needs unique headers including {true_column!r} "
                f"and {predicted_column!r}."
            )
        try:
            for index, row in enumerate(reader, start=1):
                if index > 1_000_000:
                    raise ConfigurationError("Confusion matrix CSV exceeds one million rows.")
                if None in row or any(value is None for value in row.values()):
                    raise InputFormatError(f"Malformed prediction CSV row {index + 1}.")
                truth, guess = row[true_column].strip(), row[predicted_column].strip()
                if not truth or not guess or len(truth) > 128 or len(guess) > 128:
                    raise InputFormatError(
                        f"Invalid class label in prediction CSV row {index + 1}."
                    )
                try:
                    validate_xml_text(truth, name="true label")
                    validate_xml_text(guess, name="predicted label")
                except ValueError as exc:
                    raise InputFormatError(
                        f"Invalid XML label in prediction CSV row {index + 1}."
                    ) from exc
                actual.append(truth)
                predicted.append(guess)
        except csv.Error as exc:
            raise InputFormatError(f"Malformed prediction CSV: {exc}.") from exc
    if not actual:
        raise InputFormatError("Prediction CSV has no samples.")
    return actual, predicted


def plot_confusion_matrix_csv(
    source: str | Path,
    *,
    true_column: str = "label",
    predicted_column: str = "prediction",
    labels: Sequence[str] | None = None,
    title: str | None = None,
    theme: SVGTheme | None = None,
    output: str | Path | None = None,
    overwrite: bool = False,
    verbose: bool = True,
) -> SVGArtifact:
    """Plot counts and true-class row proportions from one prediction CSV.

    The CSV needs true and predicted class columns; extra columns are ignored.
    Rows are true classes, columns are predicted classes. ``labels`` fixes their
    order and may include unobserved classes. Empty true-class rows have zero
    proportions. ``output`` saves SVG or another format selected by its suffix.
    """
    if (
        not isinstance(true_column, str)
        or not true_column.strip()
        or not isinstance(predicted_column, str)
        or not predicted_column.strip()
        or true_column == predicted_column
    ):
        raise ConfigurationError("true_column and predicted_column must be distinct names.")
    if title is not None and (not isinstance(title, str) or not title.strip()):
        raise ConfigurationError("title must be a non-empty string or None.")
    if not isinstance(overwrite, bool):
        raise ConfigurationError("overwrite must be a boolean.")
    resolved_theme = SVGTheme(font_family="sans-serif") if theme is None else theme
    if not isinstance(resolved_theme, SVGTheme):
        raise TypeError("theme must be SVGTheme or None.")
    actual, predicted = _read_labels(source, true_column, predicted_column)

    if labels is None:
        categories = list(dict.fromkeys((*actual, *predicted)))
    else:
        if isinstance(labels, str):
            raise ConfigurationError("labels must be a sequence of class names.")
        categories = list(labels)
        if not categories or any(
            not isinstance(label, str) or not label.strip() or len(label) > 128
            for label in categories
        ):
            raise ConfigurationError("labels must contain non-empty class names.")
        if len(set(categories)) != len(categories):
            raise ConfigurationError("labels must be unique.")
        for label in categories:
            try:
                validate_xml_text(label, name="class label")
            except ValueError as exc:
                raise ConfigurationError("labels contain invalid XML characters.") from exc
        if (set(actual) | set(predicted)) - set(categories):
            raise InputFormatError("labels must include every observed class.")
    if len(categories) > 25:
        raise ConfigurationError("At most 25 classes may be plotted.")

    indices = {label: index for index, label in enumerate(categories)}
    size = len(categories)
    counts = [[0] * size for _ in range(size)]
    for truth, guess in zip(actual, predicted, strict=True):
        counts[indices[truth]][indices[guess]] += 1
    row_totals = [sum(row) for row in counts]
    normalized = [
        [count / row_total if row_total else 0.0 for count in row]
        for row, row_total in zip(counts, row_totals, strict=True)
    ]

    cell = max(42, min(100, 480 // size))
    matrix_width = size * cell
    displayed = [label if len(label) <= 18 else label[:15] + "..." for label in categories]
    label_gutter = max(80, min(230, max(map(len, displayed)) * 9 + 20))
    panel_lefts = (label_gutter, label_gutter + matrix_width + label_gutter + 30)
    top = 108
    bottom = top + matrix_width
    rotated = size > 4 or any(len(label) > 9 for label in displayed)
    width = panel_lefts[1] + matrix_width + 40
    height = bottom + (175 if rotated else 110)
    metadata: dict[str, object] = {
        "labels": categories,
        "counts": counts,
        "row_normalized": normalized,
        "row_totals": row_totals,
        "sample_count": len(actual),
        "true_column": true_column,
        "predicted_column": predicted_column,
    }
    heading = title or "Confusion matrix"
    builder = SVGBuilder(
        width,
        height,
        title=heading,
        description=(
            f"Confusion matrix for {len(actual)} samples, with counts and "
            "proportions normalized within each true class."
        ),
        kind="confusion_matrix",
        theme=resolved_theme,
        metadata=metadata,
    )
    text_style = {"font_family": resolved_theme.font_family, "fill": resolved_theme.foreground}
    builder.text(
        width / 2,
        36,
        heading,
        text_anchor="middle",
        font_size=22,
        class_="plot-title",
        **text_style,
    )
    maximum = max(max(row) for row in counts)
    for panel, left, values, panel_title, ceiling in (
        ("counts", panel_lefts[0], counts, "Counts", maximum),
        ("row-normalized", panel_lefts[1], normalized, "Row proportions", 1.0),
    ):
        builder.text(
            left + matrix_width / 2,
            top - 24,
            panel_title,
            text_anchor="middle",
            font_size=20,
            class_="panel-title",
            **text_style,
        )
        for row_index, row in enumerate(values):
            y = top + row_index * cell
            for column_index, value in enumerate(row):
                x = left + column_index * cell
                intensity = value / ceiling
                builder.rect(
                    x,
                    y,
                    cell,
                    cell,
                    fill=interpolate_color("#f7fbff", "#2171b5", intensity),
                    stroke=resolved_theme.background,
                    class_="confusion-cell",
                    data_panel=panel,
                    data_true=categories[row_index],
                    data_predicted=categories[column_index],
                    data_value=value,
                )
                display_value = str(value) if panel == "counts" else f"{value:.1%}"
                builder.text(
                    x + cell / 2,
                    y + cell / 2 + 6,
                    display_value,
                    text_anchor="middle",
                    font_size=min(18, max(12, cell * 0.25)),
                    font_family=resolved_theme.font_family,
                    fill="#ffffff" if intensity >= 0.55 else resolved_theme.foreground,
                    class_="cell-value",
                )
            builder.text(
                left - 12,
                y + cell / 2 + 6,
                displayed[row_index],
                text_anchor="end",
                font_size=16,
                class_="true-tick",
                **text_style,
            )
        for column_index, label in enumerate(displayed):
            x = left + (column_index + 0.5) * cell
            y = bottom + 25
            builder.text(
                x,
                y,
                label,
                text_anchor="end" if rotated else "middle",
                transform=f"rotate(-45 {number(x)} {number(y)})" if rotated else None,
                font_size=16,
                class_="predicted-tick",
                **text_style,
            )
        builder.text(
            left + matrix_width / 2,
            height - 26,
            "Predicted",
            text_anchor="middle",
            font_size=18,
            class_="x-label",
            **text_style,
        )
        true_x = left - label_gutter + 25
        true_y = top + matrix_width / 2
        builder.text(
            true_x,
            true_y,
            "True",
            text_anchor="middle",
            font_size=18,
            transform=f"rotate(-90 {number(true_x)} {number(true_y)})",
            class_="y-label",
            **text_style,
        )

    artifact = builder.artifact()
    if output is not None:
        save_image(
            artifact,
            output,
            config=ImageExportConfig(overwrite=overwrite, create_parents=True),
            verbose=False,
        )
    return artifact


__all__ = ["plot_confusion_matrix_csv"]

_register_public_verbose(globals())
