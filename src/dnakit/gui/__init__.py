"""Graphical interface for DNAKit: a desktop window on your computer, or a hosted web demo.

``dnakit gui`` (or :func:`launch`) opens the interface locally; the same page, built with
:func:`build_app` and ``GuiSettings.public_demo()``, serves the public web demo. Requires the
optional dependencies: ``pip install "dnakit[gui]"``.
"""

from dnakit.gui._launch import build_app, launch
from dnakit.gui._settings import GuiMode, GuiSettings

__all__ = ["GuiMode", "GuiSettings", "build_app", "launch"]
