# Adapted from ClawBio/ClawBio @ 1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d.
# Original file: clawbio/common/checksums.py. See LICENSE and UPSTREAM.json at the package root.
# DNAKit changes: package namespace, isolated cache, lazy imports and native execution.
"""SHA-256 checksums and reproducibility helpers for ClawBio."""

from __future__ import annotations

import hashlib
from pathlib import Path


def sha256_file(filepath: str | Path) -> str:
    """Compute SHA-256 checksum of a file (full 64-char hex digest)."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_hex(filepath: str | Path, length: int = 16) -> str:
    """Compute truncated SHA-256 hex digest (default: 16 chars).

    This is the convention used across ClawBio skills for
    report reproducibility sections.
    """
    return sha256_file(filepath)[:length]
