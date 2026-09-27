"""Persistence for caller-created Matplotlib figures."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

from dnakit.io.artifacts import _staged_path


def save_figure(
    figure: Any,
    target: str | os.PathLike[str],
    *,
    overwrite: bool = False,
    tight_layout: bool = True,
    close: bool = True,
    verbose: bool = True,
    **kwargs: Any,
) -> Path:
    """Save a Matplotlib figure as SVG/PNG/PDF/etc. using its target suffix.

    Creates parents, stages output, optionally applies tight_layout and closes
    the figure. SVG/PDF dates are omitted by default for reproducible output.
    Other options (for example dpi) are forwarded to Figure.savefig.
    """
    path = Path(target)
    format = kwargs.pop("format", path.suffix.lstrip("."))
    if format == "svg":
        kwargs.setdefault("metadata", {"Date": None})
    elif format == "pdf":
        kwargs.setdefault("metadata", {"CreationDate": None, "ModDate": None})
    try:
        if tight_layout:
            figure.tight_layout()
        with _staged_path(path, overwrite=overwrite) as temporary:
            figure.savefig(temporary, format=format, **kwargs)
    finally:
        if close:
            import matplotlib.pyplot as plt

            plt.close(figure)
    return path
