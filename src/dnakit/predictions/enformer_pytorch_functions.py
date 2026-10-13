"""Offline enformer-pytorch encodings for DNA: index, one-hot and reverse-complement forms.

Each public function imports ``enformer_pytorch`` and ``torch`` only when it runs and
takes local inputs. Pretrained weights (``from_pretrained``), model inference, fine-tuning
helpers and loss functions are not wrapped here.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError
from dnakit.predictions.dna_model_backends import (
    _clean,
    _dna,
    _integer,
    _provenance,
    _rejected,
)

_ENFORMER = "enformer-pytorch"
_EXTRAS = "neural-enformer"
_CHANNELS = ["A", "C", "G", "T"]
_INDEX_CONVENTION = "A=0, C=1, G=2, T=3, N=4"
_MAX_ONE_HOT = 100_000
_MAX_INDEX = 4


def encode_enformer_indices(
    sequence: str,
    *,
    max_length: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Integer-encode DNA with enformer-pytorch's ``str_to_seq_indices`` (A=0 ... N=4)."""

    limit = _integer(max_length, "max_length", minimum=1, maximum=_MAX_ONE_HOT)
    value = _dna(sequence, "sequence", maximum=limit)
    data = require("enformer_pytorch.data", _EXTRAS)
    try:
        indices = data.str_to_seq_indices(value)
    except (ValueError, TypeError) as exc:
        raise _rejected(_ENFORMER, exc) from exc
    return {
        "status": "completed",
        "operation": "encode_enformer_indices",
        "length": len(value),
        "indices": [int(item) for item in indices.tolist()],
        "index_convention": _INDEX_CONVENTION,
        "provenance": _provenance(_ENFORMER, encoding="str_to_seq_indices"),
    }


def indices_to_one_hot_enformer(
    indices: Sequence[int],
    *,
    padding: int = -1,
    verbose: bool = True,
) -> dict[str, Any]:
    """One-hot encode integer DNA indices with enformer-pytorch's ``seq_indices_to_one_hot``.

    Indices 0 to 3 give A, C, G, T; index 4 (N) gives all zeros. Positions equal to
    ``padding`` (a negative integer, default -1) give 0.25 in every channel. The result
    has one row per position, shape (L, 4).
    """

    if isinstance(indices, (str, bytes)) or not isinstance(indices, Sequence) or not indices:
        raise ConfigurationError("indices must be a non-empty list of integers.")
    if len(indices) > _MAX_ONE_HOT:
        raise ConfigurationError(f"indices may contain at most {_MAX_ONE_HOT} entries.")
    pad = _integer(padding, "padding", minimum=-(10**6), maximum=-1)
    values: list[int] = []
    for index, item in enumerate(indices):
        if isinstance(item, bool) or not isinstance(item, int):
            raise ConfigurationError(f"indices[{index}] must be an integer.")
        if item != pad and not 0 <= item <= _MAX_INDEX:
            raise ConfigurationError(
                f"indices[{index}] must be 0 to {_MAX_INDEX} or the padding value {pad}."
            )
        values.append(item)
    torch = require("torch", _EXTRAS)
    data = require("enformer_pytorch.data", _EXTRAS)
    try:
        tensor = data.seq_indices_to_one_hot(torch.tensor(values, dtype=torch.long), padding=pad)
    except (ValueError, TypeError, RuntimeError) as exc:
        raise _rejected(_ENFORMER, exc) from exc
    return {
        "status": "completed",
        "operation": "indices_to_one_hot_enformer",
        "length": len(values),
        "channels": list(_CHANNELS),
        "one_hot": [[_clean(float(item)) for item in row] for row in tensor.tolist()],
        "padding_value": pad,
        "provenance": _provenance(
            _ENFORMER,
            unknown_base_encoding="all zeros for N (index 4)",
            padding_encoding="0.25 in every channel",
        ),
    }


def reverse_complement_enformer(
    sequence: str,
    *,
    representation: str = "indices",
    verbose: bool = True,
) -> dict[str, Any]:
    """Reverse-complement a DNA sequence in enformer-pytorch's integer or one-hot form.

    ``representation="indices"`` returns the reverse complement as integer indices (N stays
    4). ``representation="one_hot"`` returns (L, 4) one-hot rows, with N as all zeros.
    """

    value = _dna(sequence, "sequence", maximum=100_000)
    if representation not in ("indices", "one_hot"):
        raise ConfigurationError("representation must be 'indices' or 'one_hot'.")
    data = require("enformer_pytorch.data", _EXTRAS)
    try:
        if representation == "indices":
            encoded = data.seq_indices_reverse_complement(data.str_to_seq_indices(value))
            output: Any = [int(item) for item in encoded.tolist()]
        else:
            encoded = data.one_hot_reverse_complement(data.str_to_one_hot(value))
            output = [[float(item) for item in row] for row in encoded.tolist()]
    except (ValueError, TypeError, RuntimeError) as exc:
        raise _rejected(_ENFORMER, exc) from exc
    return {
        "status": "completed",
        "operation": "reverse_complement_enformer",
        "representation": representation,
        "length": len(value),
        "reverse_complement": output,
        "provenance": _provenance(_ENFORMER, index_convention=_INDEX_CONVENTION),
    }


__all__ = [
    "encode_enformer_indices",
    "indices_to_one_hot_enformer",
    "reverse_complement_enformer",
]

_register_public_verbose(globals(), tuple(__all__))
