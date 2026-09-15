"""Checkpoint feedback that remains visible inside prediction operations."""

from __future__ import annotations

import sys
from pathlib import Path

from dnakit._verbosity import progress_enabled


def checkpoint_status(state: str, path: Path, source: str | None = None) -> None:
    """Report a checkpoint event, respecting the enclosing verbosity setting."""
    if progress_enabled():
        origin = f" (source: {source})" if source else ""
        print(
            f"[DNAKit] {state}: {path.expanduser().resolve()}{origin}",
            file=sys.stderr,
            flush=True,
        )
