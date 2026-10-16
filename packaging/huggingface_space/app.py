"""Hosted DNAKit web demo: the ``dnakit.gui`` interface with the public-demo restrictions.

The page, its handlers and the desktop app (``dnakit gui``) share the code in
``src/dnakit/gui``; this file only starts it as a web server. Used by Hugging Face Spaces
(``app_file``) and by the self-hosted systemd service (``dnakit-web.service``).
"""

from __future__ import annotations

import os

from dnakit.gui import GuiSettings, launch

if __name__ == "__main__":
    launch(
        GuiSettings.public_demo(),
        window=False,
        open_browser=False,
        host=os.environ.get("GRADIO_SERVER_NAME", "0.0.0.0"),
        port=int(os.environ.get("GRADIO_SERVER_PORT", "7860")),
    )
