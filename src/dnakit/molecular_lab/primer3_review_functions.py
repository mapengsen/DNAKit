"""primer3-py strict DNA validation for primer design input, offline and without network.

``validate_acgt_sequence()`` wraps primer3's ``ensure_acgt_uppercase``: it uppercases a sequence
and rejects anything other than A, C, G and T. Unlike the normalised check in
``dnakit.sequence``, it does not convert U to T or strip whitespace, so input with those symbols
is refused rather than silently changed.

primer3-py is GPL-2.0. It is imported only when a function runs, so dnakit does not depend on it
at import time, and no primer3 code is copied here.
"""

from __future__ import annotations

import importlib.metadata
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_MAX_SEQUENCE = 100_000


def _version() -> str | None:
    try:
        return importlib.metadata.version("primer3-py")
    except importlib.metadata.PackageNotFoundError:
        return None


def validate_acgt_sequence(sequence: str, *, verbose: bool = True) -> dict[str, Any]:
    """Uppercase a DNA sequence and require it to contain only A, C, G and T.

    Calls ``primer3.p3helpers.ensure_acgt_uppercase``. The input must be a non-empty string of at
    most 100000 characters. Lowercase a, c, g and t are accepted, and the returned ``sequence`` is
    uppercase. N, IUPAC codes, U, gaps and whitespace are rejected with the first offending symbol.
    """

    primer3 = require("primer3", "bio-primer3")
    if not isinstance(sequence, str) or not sequence or len(sequence) > _MAX_SEQUENCE:
        raise ConfigurationError(
            f"sequence must be a non-empty string of at most {_MAX_SEQUENCE} characters."
        )
    try:
        checked = str(primer3.p3helpers.ensure_acgt_uppercase(sequence))
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"sequence must contain only A, C, G and T: {exc}") from exc
    return {
        "status": "completed",
        "operation": "validate_acgt_sequence",
        "length": len(checked),
        "sequence": checked,
        "provenance": {
            "package": "primer3-py",
            "version": _version(),
            "license": "GPL-2.0",
            "method": "p3helpers.ensure_acgt_uppercase",
            "coordinate_system": "sequence positions; no genomic coordinates",
            "convention": "strict: rejects U, N, IUPAC codes and whitespace instead of "
            "converting them",
        },
    }


__all__ = [
    "validate_acgt_sequence",
]

_register_public_verbose(globals(), tuple(__all__))
