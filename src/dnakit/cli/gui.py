"""Open the DNAKit graphical interface from the command line."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from dnakit.exceptions import DNAKitError


def gui_command(
    browser: Annotated[
        bool,
        typer.Option("--browser", help="Open in the default web browser instead of a window."),
    ] = False,
    workdir: Annotated[
        Path | None,
        typer.Option(
            "--workdir",
            file_okay=False,
            help="Folder for uploaded inputs and results (default: ./dnakit-gui).",
        ),
    ] = None,
    host: Annotated[
        str, typer.Option("--host", help="Address to serve on with --browser or --no-open.")
    ] = "127.0.0.1",
    port: Annotated[
        int | None,
        typer.Option("--port", min=1, max=65535, help="Port (default: first free port from 7860)."),
    ] = None,
    no_open: Annotated[
        bool,
        typer.Option("--no-open", help="Only serve the page; do not open a window or browser."),
    ] = False,
    public_demo: Annotated[
        bool,
        typer.Option(
            "--public-demo",
            help="Apply the hosted web-demo restrictions (time limits, no file-writing tools).",
        ),
    ] = False,
    self_test: Annotated[
        float | None,
        typer.Option("--self-test", hidden=True, help="Open, check rendering, close (seconds)."),
    ] = None,
    screenshot: Annotated[
        Path | None, typer.Option("--screenshot", hidden=True, help="Save a screen capture.")
    ] = None,
) -> None:
    """Open the DNAKit graphical interface in a desktop window."""

    from dnakit.gui import GuiSettings, launch

    settings = GuiSettings.public_demo(workdir) if public_demo else GuiSettings.local(workdir)
    try:
        ok = launch(
            settings,
            window=not (browser or no_open),
            open_browser=not no_open,
            host=host,
            port=port,
            self_test=self_test,
            screenshot=screenshot,
        )
    except DNAKitError as exc:
        typer.echo(str(exc), err=True)
        raise typer.Exit(1) from exc
    if self_test is not None:
        typer.echo(f"DNAKIT_GUI_SELF_TEST {'ok' if ok else 'failed'}")
        if not ok:
            raise typer.Exit(1)
