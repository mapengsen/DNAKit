"""Settings that separate the local desktop interface from the hosted web demo."""

from __future__ import annotations

import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

GuiMode = Literal["local", "public"]

# Single-sequence tools pass the sequence on the ``dnakit`` command line; Windows limits a
# command line to 32,767 characters, so longer sequences go through the Dataset tab (files).
COMMAND_LINE_SEQUENCE_LIMIT = 30_000


@dataclass(frozen=True)
class GuiSettings:
    """Limits and permissions of one DNAKit graphical interface.

    ``local`` is the desktop application started by ``dnakit gui``: one user on their own
    computer, results saved in ``workdir``, every function allowed after confirmation and no
    time limit. ``public`` is the hosted web demo: per-visitor temporary folders, functions that
    write files, download models or run external programs are not executed, path arguments must
    stay inside the visitor's folder, and every call has a time limit.
    """

    mode: GuiMode
    workdir: Path
    call_timeout_s: float | None
    cli_timeout_s: float | None
    max_upload_mb: int | None
    max_sequence_length: int
    max_record_length: int | None
    max_table_records: int | None
    allow_side_effects: bool
    confine_paths: bool
    session_max_age_s: int | None

    @property
    def is_public(self) -> bool:
        """Whether the public-demo restrictions apply."""

        return self.mode == "public"

    @classmethod
    def local(cls, workdir: str | Path | None = None) -> GuiSettings:
        """Return settings for the desktop interface on the user's own computer."""

        folder = Path.cwd() / "dnakit-gui" if workdir is None else Path(workdir)
        return cls(
            mode="local",
            workdir=folder.expanduser().resolve(),
            call_timeout_s=None,
            cli_timeout_s=None,
            max_upload_mb=None,
            max_sequence_length=COMMAND_LINE_SEQUENCE_LIMIT,
            max_record_length=None,
            max_table_records=None,
            allow_side_effects=True,
            confine_paths=False,
            session_max_age_s=None,
        )

    @classmethod
    def public_demo(cls, workdir: str | Path | None = None) -> GuiSettings:
        """Return settings for a shared, internet-facing web demo."""

        folder = Path(tempfile.gettempdir()) / "dnakit_space" if workdir is None else Path(workdir)
        return cls(
            mode="public",
            workdir=folder.expanduser().resolve(),
            call_timeout_s=60,
            cli_timeout_s=120,
            max_upload_mb=20,
            max_sequence_length=COMMAND_LINE_SEQUENCE_LIMIT,
            max_record_length=50_000,
            max_table_records=1_000,
            allow_side_effects=False,
            confine_paths=True,
            session_max_age_s=6 * 3600,
        )
