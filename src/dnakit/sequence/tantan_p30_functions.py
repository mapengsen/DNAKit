"""Tandem-repeat masking through Biotite 2's fixed tantan application."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as cli
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError


def mask_tandem_repeats(sequences: Sequence[str], *, executable_path: str | None = None,
                        timeout_seconds: float = 300, verbose: bool = True) -> dict[str, Any]:
    """Run tantan on DNA strings and return repeat masks and lowercase-masked sequences.

    Bases are A/C/G/T/N (case insensitive); the original length and order are preserved.
    True mask entries denote tandem repeat positions, following TantanApp.run.
    The executable is external and must be installed by the caller.
    """
    from dnakit.secondary_structure.viennarna_p30_functions import _number, _strings

    texts = _strings(sequences)
    timeout = _number(timeout_seconds, "timeout_seconds", 0.01, 86400)
    executable = cli.program("tantan", explicit=executable_path,
                             hint="Install tantan (conda install -c bioconda tantan) and dnakit[bio-biotite].")
    seq = require("biotite.sequence", "bio-biotite")
    apps = require("biotite.application.tantan", "bio-biotite")
    future = apps.TantanApp(path=executable).run([seq.NucleotideSequence(s) for s in texts])
    try:
        result = future.result(timeout=timeout)
    except TimeoutError:
        future.cancel()
        raise
    masks = [[bool(flag) for flag in mask] for mask in result]
    if len(masks) != len(texts) or any(len(mask) != len(text) for mask, text in zip(masks, texts)):
        raise ConfigurationError("tantan returned masks with unexpected dimensions.")
    return {"status": "completed", "backend": "tantan", "masks": masks,
            "masked_sequences": ["".join(base.lower() if flag else base for base, flag in zip(text, mask))
                                 for text, mask in zip(texts, masks)],
            "provenance": {"executable": executable, "version": cli._version_of(executable)}}


__all__ = ["mask_tandem_repeats"]
_register_public_verbose(globals(), tuple(__all__))
