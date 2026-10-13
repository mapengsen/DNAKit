"""Pyhmmer alignment of DNA sequences to a DNA profile HMM (hmmalign).

Only the profile alphabet decides the search space: the profile must be a DNA HMM in HMMER3
format, and every sequence is digitised with that alphabet. The other pyhmmer entries in the
inventory are protein searches (phmmer, jackhmmer, hmmscan), a press step (hmmpress), or aliases of
functions already wrapped in backends.py (nhmmer and hmmsearch are covered by search_dna_hmm).
pyhmmer is imported only when a function runs.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.alignment.backends import (
    _IUPAC_LETTERS,
    _PYHMMER_EXTRA,
    _SEQUENCE_LIMIT,
    _flag,
    _integer,
    _named_sequences,
    _package,
    _path_argument,
    _single_profile,
)
from dnakit.bioinformatics._common import input_path, require
from dnakit.core import DNA
from dnakit.exceptions import BackendExecutionError, ConfigurationError


def pyhmmer_align_to_hmm(
    profile: str | Path,
    sequences: DNA | Mapping[str, str],
    *,
    trim: bool = False,
    all_consensus_cols: bool = True,
    max_sequences: int = 10_000,
    max_length: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align named DNA sequences to one DNA profile HMM with pyhmmer's hmmalign."""

    trim_columns = _flag(trim, "trim")
    consensus_columns = _flag(all_consensus_cols, "all_consensus_cols")
    count_limit = _integer(max_sequences, "max_sequences", minimum=1, maximum=1_000_000)
    length_limit = _integer(max_length, "max_length", minimum=1, maximum=_SEQUENCE_LIMIT)
    profile_path = input_path(_path_argument(profile, "profile"))
    named = _named_sequences(
        sequences,
        allowed=_IUPAC_LETTERS,
        max_count=count_limit,
        limit_name="max_sequences",
        max_length=length_limit,
    )
    pyhmmer = require("pyhmmer", _PYHMMER_EXTRA)
    hmm = _single_profile(pyhmmer, profile_path)
    if not hmm.alphabet.is_dna():
        raise ConfigurationError(
            f"pyhmmer_align_to_hmm needs a DNA profile; {profile_path.name} uses "
            f"{hmm.alphabet.type}."
        )
    alphabet = hmm.alphabet
    digital = [
        pyhmmer.easel.TextSequence(name=name, sequence=text).digitize(alphabet)
        for name, text in named
    ]
    try:
        msa = pyhmmer.hmmalign(
            hmm, digital, trim=trim_columns, all_consensus_cols=consensus_columns
        )
    except (ValueError, pyhmmer.errors.EaselError) as exc:
        raise BackendExecutionError(f"pyhmmer hmmalign failed: {exc}") from exc
    aligned = {str(name): str(row) for name, row in zip(msa.names, msa.alignment, strict=True)}
    rows = [{"name": name, "aligned": aligned[name]} for name, _ in named]
    columns = len(rows[0]["aligned"]) if rows else 0
    return {
        "status": "completed",
        "operation": "pyhmmer_align_to_hmm",
        "profile": str(hmm.name),
        "profile_length": int(hmm.M),
        "consensus": str(hmm.consensus),
        "sequence_count": len(rows),
        "alignment_columns": columns,
        "trim": trim_columns,
        "all_consensus_cols": consensus_columns,
        "rows": rows,
        "provenance": {
            "package": _package("pyhmmer"),
            "algorithm": "pyhmmer hmmalign (profile-based alignment, default Plan 7 settings)",
            "alignment": "upper-case letters are match columns; lower-case letters are insert "
            "residues outside the model; '.' marks a gap in a row; trim=True keeps only "
            "model columns",
            "coordinate_system": "model positions are 1-based match columns of the profile",
            "alphabet": "the profile alphabet (DNA); IUPAC codes are accepted",
        },
    }


__all__ = ["pyhmmer_align_to_hmm"]

_register_public_verbose(globals(), tuple(__all__))
