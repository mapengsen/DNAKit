"""scikit-bio-backed DNA functions for ``dnakit.sequence``.

Each public function imports scikit-bio only when it is called, so the module
imports without the optional ``bio-skibio`` extra. Inputs are explicit
strings, mappings or local files. Ecology, count-table statistics and
protein-only utilities are not wrapped here.

Distances are reported as plain floats. Undefined values, such as a distance
that scikit-bio returns as NaN, are reported as ``None`` with an explicit
status instead of being replaced by a number.
"""

from __future__ import annotations

import importlib.metadata
import math
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-skibio"
_PACKAGE = "scikit-bio"
_IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")
_GAPPED_DNA = _IUPAC_DNA | frozenset("-")
_MAX_SEQUENCE = 1_000_000
_MAX_ROWS = 500
_DISTANCE_MODELS = ("pdist", "jc69", "k2p", "f81", "f84", "tn93", "logdet", "paralin")
_MSA_METRICS = ("hamming", *_DISTANCE_MODELS)
_MATRIX_MODELS = ("jc69", "k2p", "f81", "hky85", "tn93")


def _provenance(**notes: str) -> dict[str, Any]:
    try:
        installed: str | None = importlib.metadata.version(_PACKAGE)
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": _PACKAGE, "package_version": installed, **notes}


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _number(value: object, name: str, *, minimum: float | None = None) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(f"{name} must be a finite number.")
    number = float(value)
    if minimum is not None and number < minimum:
        raise ConfigurationError(f"{name} must be at least {minimum}.")
    return number


def _choice(value: object, allowed: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}; got {value!r}.")
    return value


def _sequence(
    value: object,
    name: str,
    *,
    alphabet: frozenset[str],
    maximum: int = _MAX_SEQUENCE,
) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty string.")
    text = "".join(value.split()).upper()
    if len(text) > maximum:
        raise ConfigurationError(f"{name} exceeds the limit of {maximum} symbols.")
    invalid = sorted(set(text) - alphabet)
    if invalid:
        raise ConfigurationError(f"{name} contains unsupported symbols {invalid}.")
    return text


def _named_rows(
    value: object,
    name: str,
    *,
    alphabet: frozenset[str],
    minimum_rows: int,
    max_rows: int,
    maximum: int = _MAX_SEQUENCE,
    equal_length: bool = True,
) -> dict[str, str]:
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"{name} must be a mapping from identifier to sequence.")
    if not minimum_rows <= len(value) <= max_rows:
        raise ConfigurationError(
            f"{name} must contain between {minimum_rows} and {max_rows} sequences."
        )
    rows: dict[str, str] = {}
    for key in sorted(value, key=str):
        identifier = str(key).strip()
        if not identifier or identifier in rows:
            raise ConfigurationError(f"{name} identifiers must be unique and non-empty.")
        rows[identifier] = _sequence(
            value[key], f"{name}[{identifier}]", alphabet=alphabet, maximum=maximum
        )
    if equal_length and len({len(text) for text in rows.values()}) != 1:
        raise ConfigurationError(f"all rows in {name} must have equal length.")
    return rows


def _dna_objects(rows: Mapping[str, str]) -> Any:
    skbio = require("skbio", _EXTRA)
    return skbio.TabularMSA([skbio.DNA(text) for text in rows.values()], index=list(rows.keys()))


def nucleotide_distance(
    seq1: str,
    seq2: str,
    *,
    model: str = "pdist",
    verbose: bool = True,
) -> dict[str, Any]:
    """Return one evolutionary distance between two aligned, gap-free DNA sequences."""

    selected = _choice(model, _DISTANCE_MODELS, "model")
    left = _sequence(seq1, "seq1", alphabet=_IUPAC_DNA)
    right = _sequence(seq2, "seq2", alphabet=_IUPAC_DNA)
    if len(left) != len(right):
        raise ConfigurationError("seq1 and seq2 must be aligned sequences of equal length.")
    distance = require("skbio.sequence.distance", _EXTRA)
    skbio = require("skbio", _EXTRA)
    function = getattr(distance, selected)
    try:
        value = float(function(skbio.DNA(left), skbio.DNA(right)))
    except ValueError as exc:
        raise ConfigurationError(f"scikit-bio could not compute {selected}: {exc}") from exc
    undefined = math.isnan(value)
    return {
        "status": "undefined" if undefined else "completed",
        "operation": "nucleotide_distance",
        "model": selected,
        "sequence_length": len(left),
        "distance": None if undefined else value,
        "provenance": _provenance(
            parser=f"skbio.sequence.distance.{selected}",
            gaps="gap symbols are rejected; use gap-free aligned sequences",
            undefined="a NaN from scikit-bio is reported as distance null",
        ),
    }


def nucleotide_transition_matrix(
    distance: float,
    *,
    model: str = "jc69",
    kappa: float | None = None,
    kappa_purine: float | None = None,
    kappa_pyrimidine: float | None = None,
    base_frequencies: Sequence[float] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the 4 x 4 nucleotide transition probability matrix for one distance."""

    selected = _choice(model, _MATRIX_MODELS, "model")
    time = _number(distance, "distance", minimum=0.0)
    freqs = _frequencies(base_frequencies)
    transition = require("skbio.sequence.transition", _EXTRA)
    try:
        if selected == "jc69":
            result = transition.jc69(time, dtype="float64")
        elif selected == "k2p":
            result = transition.k2p(time, _kappa(kappa, "kappa"), dtype="float64")
        elif selected == "f81":
            result = transition.f81(time, freqs, dtype="float64")
        elif selected == "hky85":
            result = transition.hky85(time, freqs, _kappa(kappa, "kappa"), dtype="float64")
        else:
            result = transition.tn93(
                time,
                freqs,
                _kappa(kappa_purine, "kappa_purine"),
                _kappa(kappa_pyrimidine, "kappa_pyrimidine"),
                dtype="float64",
            )
    except ValueError as exc:
        raise ConfigurationError(f"scikit-bio could not build the matrix: {exc}") from exc
    matrix = [[float(value) for value in row] for row in result.data.tolist()]
    return {
        "status": "completed",
        "operation": "nucleotide_transition_matrix",
        "model": selected,
        "distance": time,
        "alphabet": "ACGT",
        "base_frequencies": list(freqs),
        "matrix": matrix,
        "provenance": _provenance(
            parser=f"skbio.sequence.transition.{selected}",
            order="rows are the source base and columns the target base, in A, C, G, T order",
        ),
    }


def _frequencies(values: Sequence[float] | None) -> tuple[float, ...]:
    if values is None:
        return (0.25, 0.25, 0.25, 0.25)
    if isinstance(values, str) or not isinstance(values, Sequence) or len(values) != 4:
        raise ConfigurationError("base_frequencies must be four numbers in A, C, G, T order.")
    numbers = tuple(_number(value, "base_frequencies", minimum=0.0) for value in values)
    if not math.isclose(sum(numbers), 1.0, abs_tol=1e-6) or min(numbers) <= 0.0:
        raise ConfigurationError("base_frequencies must be positive and sum to 1.")
    return numbers


def _kappa(value: object, name: str) -> float:
    """Transition/transversion parameters are relative and must lie in (0, 1] for scikit-bio."""

    if value is None:
        raise ConfigurationError(f"{name} is required for this model.")
    number = _number(value, name, minimum=0.0)
    if not 0.0 < number <= 1.0:
        raise ConfigurationError(f"{name} must be greater than 0 and at most 1.")
    return number


def alignment_distance_matrix(
    aligned_sequences: Mapping[str, str],
    *,
    metric: str = "hamming",
    verbose: bool = True,
) -> dict[str, Any]:
    """Return a symmetric distance matrix between the rows of a DNA multiple alignment."""

    selected = _choice(metric, _MSA_METRICS, "metric")
    rows = _named_rows(
        aligned_sequences,
        "aligned_sequences",
        alphabet=_GAPPED_DNA,
        minimum_rows=2,
        max_rows=_MAX_ROWS,
    )
    align = require("skbio.alignment", _EXTRA)
    try:
        matrix = align.align_dists(_dna_objects(rows), selected)
    except ValueError as exc:
        raise ConfigurationError(f"scikit-bio could not compute {selected}: {exc}") from exc
    ids = [str(identifier) for identifier in matrix.ids]
    return {
        "status": "completed",
        "operation": "alignment_distance_matrix",
        "metric": selected,
        "ids": ids,
        "matrix": [[float(value) for value in row] for row in matrix.data.tolist()],
        "provenance": _provenance(
            parser="skbio.alignment.align_dists",
            gaps="gap symbols are kept; scikit-bio applies its shared_by_all default",
        ),
    }


def multiple_nucleotide_alignment(
    sequences: Mapping[str, str],
    *,
    match_score: float = 1.0,
    mismatch_score: float = -1.0,
    gap_cost: float = 2.0,
    free_ends: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Progressively align several DNA sequences with scikit-bio, without an external binary."""

    rows = _named_rows(
        sequences,
        "sequences",
        alphabet=_IUPAC_DNA,
        minimum_rows=2,
        max_rows=100,
        maximum=100_000,
        equal_length=False,
    )
    match = _number(match_score, "match_score")
    mismatch = _number(mismatch_score, "mismatch_score")
    gap = _number(gap_cost, "gap_cost", minimum=0.0)
    if not isinstance(free_ends, bool):
        raise ConfigurationError("free_ends must be True or False.")
    align = require("skbio.alignment", _EXTRA)
    skbio = require("skbio", _EXTRA)
    objects = [skbio.DNA(text) for text in rows.values()]
    try:
        result = align.multi_align_nucl(
            objects, sub_score=(match, mismatch), gap_cost=gap, free_ends=free_ends
        )
        gapped = result.path.to_aligned(objects)
    except ValueError as exc:
        raise ConfigurationError(f"scikit-bio could not align the sequences: {exc}") from exc
    alignment = {identifier: str(text) for identifier, text in zip(rows, gapped, strict=True)}
    return {
        "status": "completed",
        "operation": "multiple_nucleotide_alignment",
        "ids": list(rows),
        "alignment": alignment,
        "alignment_length": len(next(iter(alignment.values()))),
        "match_score": match,
        "mismatch_score": mismatch,
        "gap_cost": gap,
        "free_ends": free_ends,
        "provenance": _provenance(
            parser="skbio.alignment.multi_align_nucl",
            method="progressive alignment with a guide tree computed by scikit-bio",
            gaps="gap characters are inserted only where the alignment needs them",
        ),
    }


def alignment_score(
    aligned_sequences: Mapping[str, str],
    *,
    match_score: float = 1.0,
    mismatch_score: float = -1.0,
    gap_cost: float = 2.0,
    free_ends: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Score an aligned set of DNA sequences with scikit-bio's alignment scoring."""

    rows = _named_rows(
        aligned_sequences,
        "aligned_sequences",
        alphabet=_GAPPED_DNA,
        minimum_rows=2,
        max_rows=_MAX_ROWS,
    )
    match = _number(match_score, "match_score")
    mismatch = _number(mismatch_score, "mismatch_score")
    gap = _number(gap_cost, "gap_cost", minimum=0.0)
    if not isinstance(free_ends, bool):
        raise ConfigurationError("free_ends must be True or False.")
    align = require("skbio.alignment", _EXTRA)
    try:
        score = float(
            align.align_score(
                _dna_objects(rows), sub_score=(match, mismatch), gap_cost=gap, free_ends=free_ends
            )
        )
    except ValueError as exc:
        raise ConfigurationError(f"scikit-bio could not score the alignment: {exc}") from exc
    return {
        "status": "completed",
        "operation": "alignment_score",
        "ids": list(rows),
        "alignment_length": len(next(iter(rows.values()))),
        "score": score,
        "match_score": match,
        "mismatch_score": mismatch,
        "gap_cost": gap,
        "free_ends": free_ends,
        "provenance": _provenance(
            parser="skbio.alignment.align_score",
            gaps="gap characters in the input are scored with gap_cost",
        ),
    }


def sniff_sequence_format(path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Detect the format of a local sequence file and the reader options scikit-bio suggests."""

    source = input_path(path)
    io = require("skbio.io", _EXTRA)
    try:
        detected, options = io.sniff(str(source))
    except (io.UnrecognizedFormatError, ValueError) as exc:
        raise ConfigurationError(
            f"scikit-bio could not identify the format of {source.name}: {exc}"
        ) from exc
    return {
        "status": "completed",
        "operation": "sniff_sequence_format",
        "source": str(source),
        "format": str(detected),
        "reader_options": {
            str(key): value
            if isinstance(value, (str, int, float, bool, type(None)))
            else str(value)
            for key, value in sorted(options.items(), key=lambda item: str(item[0]))
        },
        "provenance": _provenance(
            parser="skbio.io.sniff",
            scope="reads the start of the file only; the file is not converted",
        ),
    }


__all__ = [
    "alignment_distance_matrix",
    "alignment_score",
    "multiple_nucleotide_alignment",
    "nucleotide_distance",
    "nucleotide_transition_matrix",
    "sniff_sequence_format",
]

_register_public_verbose(globals(), tuple(__all__))
