"""Publish the Agent installation guide as raw Markdown with the documentation."""

from pathlib import Path
from shutil import copyfile

from mkdocs.config.defaults import MkDocsConfig


def on_post_build(config: MkDocsConfig) -> None:
    source = Path(__file__).resolve().parents[2] / "setup.md"
    copyfile(source, Path(config.site_dir) / "setup.md")
