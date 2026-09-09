"""Plot CSV sample proportions from the DNAKit command line."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, cast

import typer

from dnakit.exceptions import DNAKitError
from dnakit.visualization import (
    DistributionPlotConfig,
    DistributionTask,
    ImageExportConfig,
    SVGTheme,
    plot_csv_distributions,
    save_image,
)


def plot_distribution_command(
    input_path: Annotated[Path, typer.Option("--input", exists=True, dir_okay=False)],
    column: Annotated[list[str], typer.Option("--column", help="Repeat for multiple columns.")],
    task: Annotated[str, typer.Option("--task", help="classification or regression.")],
    output: Annotated[
        list[Path], typer.Option("--output", help="One filename per column, in the same order.")
    ],
    bins: Annotated[int, typer.Option("--bins", min=1, max=1000)] = 30,
    bin_edges: Annotated[
        str | None,
        typer.Option("--bin-edges", help="Comma-separated numeric boundaries; overrides --bins."),
    ] = None,
    category: Annotated[
        list[str] | None,
        typer.Option(
            "--category", help="Repeat to set category order, including absent categories."
        ),
    ] = None,
    size: Annotated[
        int, typer.Option("--size", min=320, help="Square SVG canvas size in pixels.")
    ] = 720,
    font_family: Annotated[str, typer.Option("--font-family")] = "sans-serif",
    font_size: Annotated[int, typer.Option("--font-size", min=1)] = 16,
    tick_font_size: Annotated[int, typer.Option("--tick-font-size", min=1)] = 12,
    value_font_size: Annotated[int, typer.Option("--value-font-size", min=1)] = 11,
    color: Annotated[str, typer.Option("--color")] = "#2563eb",
    background: Annotated[str, typer.Option("--background")] = "#ffffff",
    bar_width: Annotated[float, typer.Option("--bar-width", min=0.01, max=1.0)] = 0.9,
    label_rotation: Annotated[
        int, typer.Option("--label-rotation", help="0, 45, or 90 degrees.")
    ] = 0,
    max_tick_labels: Annotated[int, typer.Option("--max-tick-labels", min=1)] = 12,
    show_values: Annotated[bool, typer.Option("--show-values/--hide-values")] = False,
    show_grid: Annotated[bool, typer.Option("--grid/--no-grid")] = True,
    x_label: Annotated[str | None, typer.Option("--x-label")] = None,
    y_label: Annotated[str, typer.Option("--y-label")] = "Sample proportion (%)",
    dpi: Annotated[int, typer.Option("--dpi", min=1)] = 600,
    encoding: Annotated[str, typer.Option("--encoding")] = "utf-8-sig",
    delimiter: Annotated[str, typer.Option("--delimiter")] = ",",
    max_values: Annotated[int, typer.Option("--max-values", min=1)] = 1_000_000,
    overwrite: Annotated[bool, typer.Option("--overwrite")] = False,
) -> None:
    """Plot category proportions or regression histograms from CSV columns."""
    try:
        if len(column) != len(output):
            raise ValueError("Supply one --output for each --column, in the same order.")
        # Extensionless filenames follow save_image's PNG default.
        targets = [path if path.suffix else path.with_suffix(".png") for path in output]
        paths = [path.resolve() for path in targets]
        if len(set(paths)) != len(paths):
            raise ValueError("Output filenames must be distinct.")
        for path in paths:
            if path == input_path.resolve():
                raise ValueError("Output must not overwrite the input CSV.")
            if path.suffix.lower() not in (
                ".png",
                ".svg",
                ".jpg",
                ".jpeg",
                ".tif",
                ".tiff",
                ".pdf",
            ):
                raise ValueError(f"Unsupported image extension: {path.suffix}.")
            if path.exists() and (not overwrite or path.is_dir()):
                raise FileExistsError(
                    f"Refusing to overwrite {path}; choose a new file or use --overwrite."
                )
        edges = None if bin_edges is None else tuple(float(value) for value in bin_edges.split(","))
        config = DistributionPlotConfig(
            bins=bins,
            bin_edges=edges,
            categories=None if category is None else tuple(category),
            size=size,
            font_size=font_size,
            tick_font_size=tick_font_size,
            value_font_size=value_font_size,
            bar_width=bar_width,
            label_rotation=label_rotation,
            max_tick_labels=max_tick_labels,
            show_values=show_values,
            show_grid=show_grid,
            x_label=x_label,
            y_label=y_label,
            max_values=max_values,
            theme=SVGTheme(font_family=font_family, accent=color, background=background),
        )
        plots = plot_csv_distributions(
            input_path,
            columns=column,
            task=cast(DistributionTask, task),
            config=config,
            encoding=encoding,
            delimiter=delimiter,
        )
        for name, target in zip(column, targets, strict=True):
            artifact = plots[name]
            save_image(
                artifact,
                target,
                config=ImageExportConfig(dpi=dpi, overwrite=overwrite, create_parents=True),
            )
            typer.echo(
                f"{name}: {artifact.metadata['valid_count']} valid, "
                f"{artifact.metadata['missing_count']} blank; saved {target}"
            )
    except (DNAKitError, OSError, TypeError, ValueError, LookupError) as exc:
        typer.echo(str(exc), err=True)
        raise typer.Exit(code=2) from exc
