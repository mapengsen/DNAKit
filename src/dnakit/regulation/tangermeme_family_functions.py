"""Tangermeme wrappers for variants that the existing tangermeme wrappers do not reach.

``design_beam_substitution`` wraps the beam-search designer, which keeps several candidate
sequences per round where ``design_greedy_substitution`` keeps one. ``insert_motif`` returns the
sequences after a motif insertion without trimming them back to the input length. The caller
supplies a PyTorch model in eval mode, as for the other tangermeme functions. Packages are imported
only when a function runs. A missing package raises ``BackendUnavailableError`` with the
``pip install "dnakit[bio-tangermeme]"`` hint.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError
from dnakit.regulation.backends import _dna, _flag, _integer, _tensor_provenance, _torch_model
from dnakit.regulation.tangermeme_functions import (
    _MAX_DESIGN_ITERATIONS,
    _MAX_MOTIF,
    _TANGERMEME_EXTRA,
    _call,
    _encode,
    _fraction_like,
    _sequence_batch,
)

_MAX_BEAM = 64


def design_beam_substitution(
    model: Any,
    sequence: str,
    target: Sequence[float],
    *,
    motifs: Sequence[str] | None = None,
    reverse_complement: bool = True,
    beam_size: int = 4,
    n_best: int = 1,
    max_iter: int = 10,
    tol: float = 0.001,
    verbose: bool = True,
) -> dict[str, Any]:
    """Beam search over motif substitutions that move the model output towards a target.

    Calls ``tangermeme.design.beam_substitution``. Each round tiles every motif (and its reverse
    complement when ``reverse_complement=True``) at every position of every sequence in the beam,
    scores the full sequences by their mean squared error to ``target``, and keeps the
    ``beam_size`` best distinct ones. The search stops once the best loss improves by no more than
    ``tol`` in a round or after ``max_iter`` rounds. ``beam_size=1`` equals
    ``design_greedy_substitution``. Without ``motifs`` single bases are the candidates. The
    ``n_best`` best sequences are returned, lowest loss first.
    """

    torch = _torch_model(model)
    design = require("tangermeme.design", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    if not isinstance(sequence, str):
        raise ConfigurationError("sequence must be one DNA string.")
    _names, texts = _sequence_batch({"sequence": sequence})
    if isinstance(target, (str, bytes)) or not isinstance(target, Sequence) or not target:
        raise ConfigurationError("target must be a non-empty list of model output values.")
    values = [float(_fraction_like(value)) for value in target]
    if motifs is not None:
        if isinstance(motifs, (str, bytes)) or not isinstance(motifs, Sequence) or not motifs:
            raise ConfigurationError("motifs must be a non-empty list of DNA strings or None.")
        motif_list = [
            _dna(value, f"motif {index}", _MAX_MOTIF) for index, value in enumerate(motifs)
        ]
    else:
        motif_list = None
    rc = _flag(reverse_complement, "reverse_complement")
    width = _integer(beam_size, "beam_size", 1, _MAX_BEAM)
    best = _integer(n_best, "n_best", 1, width)
    iterations = _integer(max_iter, "max_iter", 1, _MAX_DESIGN_ITERATIONS)
    tolerance = _fraction_like(tol)
    if tolerance < 0:
        raise ConfigurationError("tol must be a non-negative number.")
    inputs = _encode(torch, utils, texts)
    target_tensor = torch.tensor([values], dtype=torch.float32)
    designed = _call(
        "tangermeme",
        design.beam_substitution,
        model,
        inputs,
        y=target_tensor,
        motifs=motif_list,
        reverse_complement=rc,
        beam_size=width,
        n_best=best,
        max_iter=iterations,
        tol=tolerance,
        batch_size=32,
    )
    texts_out = [utils.characters(designed[index].detach().cpu()) for index in range(len(designed))]
    return {
        "status": "completed",
        "operation": "design_beam_substitution",
        "input_sequence": texts[0],
        "designed_sequences": texts_out,
        "designed_sequence": texts_out[0],
        "length": len(texts_out[0]),
        "target": values,
        "beam_size": width,
        "n_best": len(texts_out),
        "max_iter": iterations,
        "provenance": _tensor_provenance(
            torch,
            function="tangermeme.design.beam_substitution",
            candidates="motifs and their reverse complements" if motif_list else "single bases",
            order="lowest loss first",
        ),
    }


def insert_motif(
    sequences: object,
    motif: str,
    *,
    start: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Insert a motif into each sequence and keep the longer result.

    Calls ``tangermeme.ersatz.insert``. The motif is inserted before base ``start`` (0-based; the
    centre of each sequence when ``start`` is None), so every output is ``len(motif)`` bases longer
    than its input and no base is lost. ``place_motif(mode="insert")`` trims the result back to the
    input length from the right, which drops the last bases.
    """

    ersatz = require("tangermeme.ersatz", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    torch = require("torch", _TANGERMEME_EXTRA)
    names, texts = _sequence_batch(sequences)
    text = _dna(motif, "motif", _MAX_MOTIF)
    length = len(texts[0])
    first = None if start is None else _integer(start, "start", 0, length)
    inputs = _encode(torch, utils, texts)
    edited = _call("tangermeme", ersatz.insert, inputs, text, start=first)
    result = [utils.characters(edited[index]) for index in range(edited.shape[0])]
    return {
        "status": "completed",
        "operation": "insert_motif",
        "motif": text,
        "start": first,
        "input_length": length,
        "output_length": len(result[0]),
        "sequence_count": len(names),
        "sequences": [
            {"name": name, "sequence": value} for name, value in zip(names, result, strict=True)
        ],
        "provenance": _tensor_provenance(
            torch,
            function="tangermeme.ersatz.insert",
            convention="the motif is inserted before base start; sequences grow by len(motif)",
        ),
    }


__all__ = ["design_beam_substitution", "insert_motif"]

_register_public_verbose(globals(), tuple(__all__))
