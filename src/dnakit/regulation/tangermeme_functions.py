"""Model-based sequence functions over tangermeme: mutagenesis, variant effects, seqlets and loci.

The caller supplies a PyTorch model that is already in eval mode and takes one-hot DNA of one
length. Models are never trained or loaded here. Sequences are DNA strings (A, C, G, T) or a
named mapping of such strings; all sequences in one call must have the same length.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.regulation.backends import (
    _choose,
    _dna,
    _flag,
    _integer,
    _label,
    _named_dna,
    _nested,
    _tensor_provenance,
    _torch_model,
)

_TANGERMEME_EXTRA = "bio-tangermeme"
_BASES = "ACGT"
_MAX_SEQUENCES = 256
_MAX_LENGTH = 10_000
_MAX_VARIANTS = 10_000
_MAX_MOTIF = 100
_MAX_SHUFFLES = 1_000
_MAX_KMER = 8
_MAX_SPACINGS = 10_000
_MAX_DESIGN_ITERATIONS = 1_000


def _call(package: str, function: Any, *args: Any, **kwargs: Any) -> Any:
    """Run an upstream call and convert its failures into DNAKit errors."""

    try:
        return function(*args, **kwargs)
    except (KeyboardInterrupt, SystemExit):
        raise
    except BaseException as exc:  # tangermeme raises ValueError, TypeError and numba errors
        raise BackendExecutionError(f"{package} failed: {exc}") from exc


def _sequence_batch(
    sequences: object, max_sequences: int = _MAX_SEQUENCES
) -> tuple[list[str], list[str]]:
    """Validate named DNA sequences of one length; return their IDs and texts."""

    named = _named_dna(sequences, max_sequences, _MAX_LENGTH)
    lengths = {len(text) for _name, text in named}
    if len(lengths) != 1:
        raise ConfigurationError("all sequences in one call must have the same length.")
    return [name for name, _text in named], [text for _name, text in named]


def _encode(torch: Any, utils: Any, texts: Sequence[str]) -> Any:
    return torch.stack([utils.one_hot_encode(text) for text in texts]).float()


def _rows(
    variants: object,
    names: Sequence[str],
    *,
    has_base: bool,
    alphabet: str = _BASES,
) -> list[list[int]]:
    """Map variant dictionaries onto the COO rows tangermeme expects."""

    if isinstance(variants, (str, bytes)) or not isinstance(variants, Sequence) or not variants:
        raise ConfigurationError("variants must be a non-empty list of dictionaries.")
    if len(variants) > _MAX_VARIANTS:
        raise ConfigurationError(f"variant count exceeds {_MAX_VARIANTS}; split the request.")
    index = {name: position for position, name in enumerate(names)}
    rows: list[list[int]] = []
    for number, item in enumerate(variants):
        if not isinstance(item, Mapping):
            raise ConfigurationError(f"variant {number} must be a mapping.")
        name = item.get("sequence")
        if name not in index:
            raise ConfigurationError(f"variant {number} names an unknown sequence {name!r}.")
        position = _integer(item.get("position"), f"variant {number} position", 0, _MAX_LENGTH)
        row = [index[name], position]
        if has_base:
            base = item.get("base")
            if not isinstance(base, str) or len(base) != 1 or base not in alphabet:
                raise ConfigurationError(f"variant {number} base must be one of A, C, G, T.")
            row.append(alphabet.index(base))
        rows.append(row)
    return rows


def _check_length(texts: Sequence[str], rows: Sequence[Sequence[int]], *, insertion: bool) -> None:
    length = len(texts[0])
    for row in rows:
        if insertion and row[1] > length:
            raise ConfigurationError("an insertion position lies beyond the sequence end.")
        if not insertion and row[1] >= length:
            raise ConfigurationError("a variant position lies outside the sequence.")


def _effect_result(
    operation: str,
    names: Sequence[str],
    variant_count: int,
    result: Any,
    function: str,
    *,
    extra: Mapping[str, Any] | None = None,
    with_delta: bool = True,
) -> dict[str, Any]:
    torch = require("torch", _TANGERMEME_EXTRA)
    before = result.y_before.detach().cpu()
    after = result.y_after.detach().cpu()
    payload = {
        "status": "completed",
        "operation": operation,
        "sequence_count": len(names),
        "variant_count": variant_count,
        "sequences": list(names),
        "output_shape": list(before.shape),
        "y_before": _nested(before.tolist()),
        "y_after": _nested(after.tolist()),
        "provenance": _tensor_provenance(torch, function=function),
    }
    if with_delta:
        payload["delta"] = _nested((after - before).tolist())
    payload.update(extra or {})
    return payload


def saturation_mutagenesis(
    model: Any,
    sequences: object,
    *,
    start: int = 0,
    end: int | None = None,
    target: int | None = None,
    hypothetical: bool = False,
    batch_size: int = 32,
    verbose: bool = True,
) -> dict[str, Any]:
    """In-silico saturation mutagenesis: attribution of every base in a window, with tangermeme.

    Each single-base substitution in ``[start, end)`` is scored. The change in prediction
    (mutated minus original) is centred over the four bases at each position. ``target`` picks
    one model output; ``None`` averages over outputs. With ``hypothetical=False`` the values are
    kept only for the observed base at each position. The result is indexed as sequence, base
    (A, C, G, T), position, and covers the window ``[start, end)`` only.
    """

    torch = _torch_model(model)
    saturation = require("tangermeme.saturation_mutagenesis", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    names, texts = _sequence_batch(sequences)
    length = len(texts[0])
    first = _integer(start, "start", 0, length - 1)
    last = length if end is None else _integer(end, "end", first + 1, length)
    chosen = None if target is None else _integer(target, "target", 0, 100_000)
    hypo = _flag(hypothetical, "hypothetical")
    batch = _integer(batch_size, "batch_size", 1, 4096)
    inputs = _encode(torch, utils, texts)
    values = _call(
        "tangermeme",
        saturation.saturation_mutagenesis,
        model,
        inputs,
        start=first,
        end=last,
        batch_size=batch,
        target=chosen,
        hypothetical=hypo,
    )
    attributions = values.detach().cpu()
    return {
        "status": "completed",
        "operation": "saturation_mutagenesis",
        "sequence_count": len(names),
        "sequences": names,
        "start": first,
        "end": last,
        "hypothetical": hypo,
        "axis_order": "sequence, base (A, C, G, T), position",
        "shape": list(attributions.shape),
        "attributions": _nested(attributions.tolist()),
        "provenance": _tensor_provenance(
            torch,
            function="tangermeme.saturation_mutagenesis.saturation_mutagenesis",
            aggregation="mutant minus reference, centred over bases, projected onto observed base",
        ),
    }


def substitution_effect(
    model: Any,
    sequences: object,
    substitutions: object,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Model output before and after substituting bases (tangermeme ``substitution_effect``).

    ``substitutions`` is a list of ``{"sequence": name, "position": i, "base": "A"}`` items. Several
    items for one sequence are applied together. ``delta`` is the output after minus before.
    """

    torch = _torch_model(model)
    variant = require("tangermeme.variant_effect", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    names, texts = _sequence_batch(sequences)
    rows = _rows(substitutions, names, has_base=True)
    _check_length(texts, rows, insertion=False)
    inputs = _encode(torch, utils, texts)
    tensor = torch.tensor(rows, dtype=torch.long)
    result = _call("tangermeme", variant.substitution_effect, model, inputs, tensor)
    return _effect_result(
        "substitution_effect",
        names,
        len(rows),
        result,
        "tangermeme.variant_effect.substitution_effect",
    )


def insertion_effect(
    model: Any,
    sequences: object,
    insertions: object,
    *,
    left: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Model output before and after inserting bases, with tangermeme ``insertion_effect``.

    ``insertions`` is a list of ``{"sequence": name, "position": i, "base": "A"}`` items; position
    ``i`` is the gap before original base ``i``. Inserted sequences are trimmed back to the input
    length from the right, or from the left when ``left=True``, so the model sees the same length.
    """

    torch = _torch_model(model)
    variant = require("tangermeme.variant_effect", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    names, texts = _sequence_batch(sequences)
    rows = _rows(insertions, names, has_base=True)
    _check_length(texts, rows, insertion=True)
    trim_left = _flag(left, "left")
    inputs = _encode(torch, utils, texts)
    tensor = torch.tensor(rows, dtype=torch.long)
    result = _call("tangermeme", variant.insertion_effect, model, inputs, tensor, left=trim_left)
    return _effect_result(
        "insertion_effect",
        names,
        len(rows),
        result,
        "tangermeme.variant_effect.insertion_effect",
        extra={"left": trim_left},
    )


def deletion_effect(
    model: Any,
    sequences: object,
    deletions: object,
    *,
    left: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Model output before and after deleting bases, with tangermeme ``deletion_effect``.

    ``deletions`` is a list of ``{"sequence": name, "position": i}`` items. The sequences must be
    longer than the model input by at least the largest number of deletions in one sequence;
    tangermeme keeps the first (or, with ``left=True``, the last) model-length window after each
    deletion.
    """

    torch = _torch_model(model)
    variant = require("tangermeme.variant_effect", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    names, texts = _sequence_batch(sequences)
    rows = _rows(deletions, names, has_base=False)
    _check_length(texts, rows, insertion=False)
    trim_left = _flag(left, "left")
    inputs = _encode(torch, utils, texts)
    tensor = torch.tensor(rows, dtype=torch.long)
    result = _call("tangermeme", variant.deletion_effect, model, inputs, tensor, left=trim_left)
    return _effect_result(
        "deletion_effect",
        names,
        len(rows),
        result,
        "tangermeme.variant_effect.deletion_effect",
        extra={"left": trim_left},
    )


def ablate_region(
    model: Any,
    sequences: object,
    *,
    start: int,
    end: int,
    n: int = 20,
    shuffle: str = "dinucleotide",
    random_state: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Shuffle a region repeatedly and measure the model output, with tangermeme ``ablate``.

    ``shuffle`` is ``"dinucleotide"`` (preserves dinucleotide counts) or ``"shuffle"`` (random
    permutation). ``y_after`` is indexed as sequence, shuffle, output; ``y_before`` is the output
    of the unshuffled sequence. ``random_state`` makes the shuffles reproducible.
    """

    torch = _torch_model(model)
    ablate = require("tangermeme.ablate", _TANGERMEME_EXTRA)
    ersatz = require("tangermeme.ersatz", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    names, texts = _sequence_batch(sequences)
    length = len(texts[0])
    first = _integer(start, "start", 0, length - 1)
    last = _integer(end, "end", first + 1, length)
    count = _integer(n, "n", 1, _MAX_SHUFFLES)
    kind = _choose(shuffle, ("dinucleotide", "shuffle"), "shuffle")
    seed = _integer(random_state, "random_state", 0, 2**32 - 1)
    shuffle_fn = ersatz.dinucleotide_shuffle if kind == "dinucleotide" else ersatz.shuffle
    inputs = _encode(torch, utils, texts)
    result = _call(
        "tangermeme",
        ablate.ablate,
        model,
        inputs,
        start=first,
        end=last,
        n=count,
        shuffle_fn=shuffle_fn,
        random_state=seed,
    )
    return _effect_result(
        "ablate_region",
        names,
        count,
        result,
        "tangermeme.ablate.ablate",
        extra={"start": first, "end": last, "shuffle": kind, "random_state": seed},
        with_delta=False,
    )


def motif_pair_spacing(
    model: Any,
    sequences: object,
    motifs: Sequence[str],
    spacings: Sequence[Sequence[int]],
    *,
    start: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Model output for motifs placed in several spacings, with tangermeme ``space``.

    ``motifs`` are ACGT strings placed in order. ``spacings`` has one row per arrangement with
    ``len(motifs) - 1`` non-negative gaps. ``start`` is the first motif position; ``None`` centres
    the arrangement in each sequence.
    """

    torch = _torch_model(model)
    spacing_module = require("tangermeme.space", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    names, texts = _sequence_batch(sequences)
    if (
        isinstance(motifs, (str, bytes))
        or not isinstance(motifs, Sequence)
        or not 2 <= len(motifs) <= 10
    ):
        raise ConfigurationError("motifs must be a list of 2 to 10 DNA strings.")
    motif_texts = [_dna(value, f"motif {index}", _MAX_MOTIF) for index, value in enumerate(motifs)]
    if isinstance(spacings, (str, bytes)) or not isinstance(spacings, Sequence) or not spacings:
        raise ConfigurationError("spacings must be a non-empty list of gap lists.")
    if len(spacings) > _MAX_SPACINGS:
        raise ConfigurationError(f"spacing count exceeds {_MAX_SPACINGS}.")
    gaps = []
    for row in spacings:
        if (
            isinstance(row, (str, bytes))
            or not isinstance(row, Sequence)
            or len(row) != len(motifs) - 1
        ):
            raise ConfigurationError("each spacing row needs one gap per adjacent motif pair.")
        gaps.append([_integer(gap, "spacing", 0, len(texts[0])) for gap in row])
    total = sum(len(motif) for motif in motif_texts) + max(sum(row) for row in gaps)
    if total > len(texts[0]):
        raise ConfigurationError("motifs and spacing do not fit in the sequences.")
    first = None if start is None else _integer(start, "start", 0, len(texts[0]))
    inputs = _encode(torch, utils, texts)
    tensor = torch.tensor(gaps, dtype=torch.long)
    result = _call(
        "tangermeme",
        spacing_module.space,
        model,
        inputs,
        motif_texts,
        tensor,
        start=first,
    )
    before = result.y_before.detach().cpu()
    raw_afters = result.y_afters
    if isinstance(raw_afters, (list, tuple)):
        afters = [item.detach().cpu() for item in raw_afters]
    else:
        afters = [raw_afters[index].detach().cpu() for index in range(raw_afters.shape[0])]
    return {
        "status": "completed",
        "operation": "motif_pair_spacing",
        "sequence_count": len(names),
        "sequences": names,
        "motifs": motif_texts,
        "spacing_count": len(gaps),
        "spacings": gaps,
        "y_before": _nested(before.tolist()),
        "y_after": [_nested(item.tolist()) for item in afters],
        "provenance": _tensor_provenance(
            torch,
            function="tangermeme.space.space",
            layout="y_after is indexed by spacing row, then sequence, then output",
        ),
    }


def design_greedy_substitution(
    model: Any,
    sequence: str,
    target: Sequence[float],
    *,
    motifs: Sequence[str] | None = None,
    reverse_complement: bool = True,
    max_iter: int = 10,
    tol: float = 0.001,
    verbose: bool = True,
) -> dict[str, Any]:
    """Greedily insert motifs into one sequence to move the model output towards a target.

    Each iteration tries every candidate motif at every position (and its reverse complement when
    ``reverse_complement=True``), keeps the substitution with the lowest mean squared error to
    ``target``, and stops after ``max_iter`` insertions or once the error is below ``tol``.
    Without ``motifs``, single bases are the candidates.
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
    iterations = _integer(max_iter, "max_iter", 1, _MAX_DESIGN_ITERATIONS)
    tolerance = _fraction_like(tol)
    if tolerance < 0:
        raise ConfigurationError("tol must be a non-negative number.")
    inputs = _encode(torch, utils, texts)
    target_tensor = torch.tensor([values], dtype=torch.float32)
    designed = _call(
        "tangermeme",
        design.greedy_substitution,
        model,
        inputs,
        y=target_tensor,
        motifs=motif_list,
        reverse_complement=rc,
        max_iter=iterations,
        tol=tolerance,
        batch_size=32,
    )
    designed_text = utils.characters(designed[0].detach().cpu())
    return {
        "status": "completed",
        "operation": "design_greedy_substitution",
        "input_sequence": texts[0],
        "designed_sequence": designed_text,
        "length": len(designed_text),
        "target": values,
        "max_iter": iterations,
        "provenance": _tensor_provenance(
            torch,
            function="tangermeme.design.greedy_substitution",
            candidates="motifs and their reverse complements" if motif_list else "single bases",
        ),
    }


def _fraction_like(value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError("expected a finite number.")
    number = float(value)
    if not math.isfinite(number):
        raise ConfigurationError("expected a finite number.")
    return number


def call_seqlets(
    attributions: Sequence[Sequence[float]],
    *,
    threshold: float = 0.01,
    min_seqlet_len: int = 4,
    max_seqlet_len: int = 25,
    additional_flanks: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call seqlets (contiguous high-attribution regions) with tangermeme ``recursive_seqlets``.

    ``attributions`` is a list of rows, one per example, with one value per position. Seqlets
    are sorted by example, start and end; ``p_value`` is the empirical value from tangermeme.
    """

    seqlet_module = require("tangermeme.seqlet", _TANGERMEME_EXTRA)
    torch = require("torch", _TANGERMEME_EXTRA)
    if isinstance(attributions, (str, bytes)) or not isinstance(attributions, Sequence):
        raise ConfigurationError("attributions must be a list of per-position value rows.")
    if not attributions or len(attributions) > _MAX_SEQUENCES:
        raise ConfigurationError(f"attributions need 1 to {_MAX_SEQUENCES} rows.")
    width = len(attributions[0])
    rows = []
    for index, row in enumerate(attributions):
        if isinstance(row, (str, bytes)) or not isinstance(row, Sequence) or len(row) != width:
            raise ConfigurationError(f"attribution row {index} must have {width} values.")
        if width > _MAX_LENGTH:
            raise ConfigurationError(f"attribution rows are longer than {_MAX_LENGTH} positions.")
        rows.append([_fraction_like(value) for value in row])
    cut = _fraction_like(threshold)
    if not 0.0 < cut < 1.0:
        raise ConfigurationError("threshold must be in (0, 1).")
    shortest = _integer(min_seqlet_len, "min_seqlet_len", 1, 1_000)
    longest = _integer(max_seqlet_len, "max_seqlet_len", shortest, 1_000)
    flanks = _integer(additional_flanks, "additional_flanks", 0, 1_000)
    tensor = torch.tensor(rows, dtype=torch.float32)
    frame = _call(
        "tangermeme",
        seqlet_module.recursive_seqlets,
        tensor,
        threshold=cut,
        min_seqlet_len=shortest,
        max_seqlet_len=longest,
        additional_flanks=flanks,
    )
    seqlets = []
    for row in frame.to_dict(orient="records"):
        seqlets.append(
            {
                "example_index": int(row["example_idx"]),
                "start": int(row["start"]),
                "end": int(row["end"]),
                "attribution": float(row["attribution"]),
                "p_value": float(row["p-value"]),
            }
        )
    seqlets.sort(key=lambda item: (item["example_index"], item["start"], item["end"]))
    return {
        "status": "completed",
        "operation": "call_seqlets",
        "example_count": len(rows),
        "seqlet_count": len(seqlets),
        "seqlets": seqlets,
        "provenance": _tensor_provenance(
            torch,
            function="tangermeme.seqlet.recursive_seqlets",
            threshold=cut,
        ),
    }


def kmer_counts(
    sequences: object,
    k: int,
    *,
    scores: Mapping[str, Sequence[float]] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count contiguous k-mers per sequence, or sum per-position scores (tangermeme ``kmers``).

    Only k-mers that occur are returned, sorted by sequence and k-mer. With ``scores`` (a mapping
    of sequence name to one value per position), each k-mer instance adds the sum of the scores it
    spans, as in tangermeme ``kmers``. Without scores each instance adds 1.
    """

    torch = require("torch", _TANGERMEME_EXTRA)
    kmer_module = require("tangermeme.kmers", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    size = _integer(k, "k", 1, _MAX_KMER)
    names, texts = _sequence_batch(sequences)
    length = len(texts[0])
    if size > length:
        raise ConfigurationError("k must not exceed the sequence length.")
    inputs = _encode(torch, utils, texts)
    weights = None
    if scores is not None:
        if not isinstance(scores, Mapping) or set(scores) != set(names):
            raise ConfigurationError("scores must map every sequence name to per-position values.")
        rows = []
        for name in names:
            values = scores[name]
            if (
                isinstance(values, (str, bytes))
                or not isinstance(values, Sequence)
                or len(values) != length
            ):
                raise ConfigurationError(f"scores for {name!r} must have one value per position.")
            rows.append([_fraction_like(value) for value in values])
        weights = torch.tensor(rows, dtype=torch.float32)
    counts = _call("tangermeme", kmer_module.kmers, inputs, size, weights)
    matrix = counts.detach().cpu().numpy()
    entries = []
    for row_index, name in enumerate(names):
        for column in range(matrix.shape[1]):
            value = float(matrix[row_index, column])
            if value != 0.0:
                entries.append(
                    {"sequence": name, "kmer": _kmer_from_index(column, size), "value": value}
                )
    entries.sort(key=lambda item: (names.index(item["sequence"]), item["kmer"]))
    return {
        "status": "completed",
        "operation": "kmer_counts",
        "k": size,
        "sequence_count": len(names),
        "weighted": weights is not None,
        "entry_count": len(entries),
        "entries": entries,
        "provenance": _tensor_provenance(torch, function="tangermeme.kmers.kmers"),
    }


def _kmer_from_index(index: int, size: int) -> str:
    """Decode tangermeme's column index: the first base is the least significant base-4 digit."""

    letters = []
    for _ in range(size):
        letters.append(_BASES[index % 4])
        index //= 4
    return "".join(letters)


def extract_loci_sequences(
    loci: str | Path,
    genome: str | Path,
    *,
    in_window: int = 200,
    out_window: int = 1000,
    max_jitter: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cut fixed-width sequences around BED loci from a genome (tangermeme ``extract_loci``).

    Each locus is centred on its midpoint; ``in_window`` is the sequence width. Loci that touch
    ambiguous bases (N) are returned as N. ``out_window`` is passed through to tangermeme and only
    affects signal extraction; with no signal files it has no effect on the returned sequences.
    """

    io = require("tangermeme.io", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    torch = require("torch", _TANGERMEME_EXTRA)
    loci_path = input_path(loci)
    fasta = input_path(genome)
    width = _integer(in_window, "in_window", 2, 100_000)
    outer = _integer(out_window, "out_window", 1, 100_000)
    jitter = _integer(max_jitter, "max_jitter", 0, 10_000)
    rows = _bed_loci(loci_path)
    tensor = _call(
        "tangermeme",
        io.extract_loci,
        str(loci_path),
        str(fasta),
        in_window=width,
        out_window=outer,
        max_jitter=jitter,
    )
    texts = [utils.characters(tensor[index]) for index in range(tensor.shape[0])]
    if len(texts) != len(rows):
        raise BackendExecutionError(
            "tangermeme returned a different number of sequences than loci."
        )
    extracted = []
    for row, text in zip(rows, texts, strict=True):
        extracted.append({"chrom": row[0], "start": row[1], "end": row[2], "sequence": text})
    return {
        "status": "completed",
        "operation": "extract_loci_sequences",
        "locus_count": len(extracted),
        "in_window": width,
        "max_jitter": jitter,
        "shape": list(tensor.shape),
        "loci": extracted,
        "provenance": _tensor_provenance(torch, function="tangermeme.io.extract_loci"),
    }


def _bed_loci(path: Path) -> list[tuple[str, int, int]]:
    rows: list[tuple[str, int, int]] = []
    with path.open("rt", encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith(("track", "browser")):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                raise ConfigurationError(f"{path.name} line {number} needs chrom, start and end.")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ConfigurationError(
                    f"{path.name} line {number} has non-integer coordinates."
                ) from exc
            if not 0 <= start < end:
                raise ConfigurationError(
                    f"{path.name} line {number} must satisfy 0 <= start < end."
                )
            rows.append((_label(fields[0], "chromosome"), start, end))
            if len(rows) > 100_000:
                raise ConfigurationError(f"{path.name} has more than 100000 loci.")
    if not rows:
        raise ConfigurationError(f"{path.name} contains no loci.")
    return rows


def matched_background_loci(
    loci: str | Path,
    genome: str | Path,
    *,
    in_window: int = 200,
    out_window: int = 1000,
    max_n_perc: float = 0.1,
    gc_bin_width: float = 0.02,
    random_state: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Background loci with GC content matched to the input (tangermeme ``extract_matching_loci``).

    Candidate windows are sampled from the genome FASTA, excluding windows with more than
    ``max_n_perc`` N, and kept in GC bins of width ``gc_bin_width`` that match the input loci.
    ``random_state`` makes the sample reproducible.
    """

    match = require("tangermeme.match", _TANGERMEME_EXTRA)
    pandas = require("pandas", _TANGERMEME_EXTRA)
    torch = require("torch", _TANGERMEME_EXTRA)
    loci_path = input_path(loci)
    fasta = input_path(genome)
    width = _integer(in_window, "in_window", 2, 100_000)
    outer = _integer(out_window, "out_window", 1, 100_000)
    if (
        isinstance(max_n_perc, bool)
        or not isinstance(max_n_perc, (int, float))
        or not 0 <= max_n_perc <= 1
    ):
        raise ConfigurationError("max_n_perc must be a number in [0, 1].")
    if (
        isinstance(gc_bin_width, bool)
        or not isinstance(gc_bin_width, (int, float))
        or not 0 < gc_bin_width <= 1
    ):
        raise ConfigurationError("gc_bin_width must be a number in (0, 1].")
    seed = _integer(random_state, "random_state", 0, 2**32 - 1)
    rows = _bed_loci(loci_path)
    frame = pandas.DataFrame(rows, columns=["chrom", "start", "end"])
    matched = _call(
        "tangermeme",
        match.extract_matching_loci,
        frame,
        str(fasta),
        in_window=width,
        out_window=outer,
        max_n_perc=float(max_n_perc),
        gc_bin_width=float(gc_bin_width),
        random_state=seed,
        n_jobs=1,
    )
    intervals = []
    for row in matched.itertuples(index=False):
        intervals.append({"chrom": str(row.chrom), "start": int(row.start), "end": int(row.end)})
    intervals.sort(key=lambda item: (item["chrom"], item["start"], item["end"]))
    return {
        "status": "completed",
        "operation": "matched_background_loci",
        "input_locus_count": len(rows),
        "background_count": len(intervals),
        "random_state": seed,
        "loci": intervals,
        "provenance": _tensor_provenance(
            torch,
            function="tangermeme.match.extract_matching_loci",
            matching="GC bins with width gc_bin_width; windows with too many N are excluded",
        ),
    }


def pwm_information_content(
    pwm: Sequence[Sequence[float]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Information content (bits) of each PWM position, with tangermeme ``information_content``.

    ``pwm`` has one row per position and four probabilities (A, C, G, T); each row should sum to 1.
    """

    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    torch = require("torch", _TANGERMEME_EXTRA)
    if (
        isinstance(pwm, (str, bytes))
        or not isinstance(pwm, Sequence)
        or not 1 <= len(pwm) <= _MAX_MOTIF * 10
    ):
        raise ConfigurationError("pwm must be a list of probability rows.")
    rows = []
    for index, row in enumerate(pwm):
        if isinstance(row, (str, bytes)) or not isinstance(row, Sequence) or len(row) != 4:
            raise ConfigurationError(f"pwm row {index} must have four values (A, C, G, T).")
        values = [_fraction_like(value) for value in row]
        if any(value < 0 for value in values):
            raise ConfigurationError(f"pwm row {index} must not contain negative values.")
        rows.append(values)
    bits = _call("tangermeme", utils.information_content, torch.tensor(rows, dtype=torch.float32))
    values = [float(item) for item in bits.tolist()]
    return {
        "status": "completed",
        "operation": "pwm_information_content",
        "length": len(values),
        "information_content": values,
        "total": sum(values),
        "provenance": _tensor_provenance(torch, function="tangermeme.utils.information_content"),
    }


def place_motif(
    sequences: object,
    motif: str,
    *,
    start: int | None = None,
    mode: str = "substitute",
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a motif into each sequence, replacing bases or inserting them, with tangermeme.

    ``mode="substitute"`` overwrites the bases at ``start``; ``mode="insert"`` inserts the motif
    and keeps the original length by trimming from the right. ``start=None`` centres the motif.
    """

    ersatz = require("tangermeme.ersatz", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    torch = require("torch", _TANGERMEME_EXTRA)
    names, texts = _sequence_batch(sequences)
    text = _dna(motif, "motif", _MAX_MOTIF)
    chosen = _choose(mode, ("substitute", "insert"), "mode")
    length = len(texts[0])
    if len(text) > length:
        raise ConfigurationError("the motif is longer than the sequences.")
    first = None if start is None else _integer(start, "start", 0, length - len(text))
    inputs = _encode(torch, utils, texts)
    function = ersatz.substitute if chosen == "substitute" else ersatz.insert
    edited = _call("tangermeme", function, inputs, text, start=first)
    result = [utils.characters(edited[index]) for index in range(edited.shape[0])]
    if chosen == "insert":
        result = [item[:length] for item in result]
    return {
        "status": "completed",
        "operation": "place_motif",
        "mode": chosen,
        "motif": text,
        "sequence_count": len(names),
        "sequences": [
            {"name": name, "sequence": value} for name, value in zip(names, result, strict=True)
        ],
        "provenance": _tensor_provenance(
            torch,
            function=f"tangermeme.ersatz.{'substitute' if chosen == 'substitute' else 'insert'}",
        ),
    }


__all__ = [
    "ablate_region",
    "call_seqlets",
    "deletion_effect",
    "design_greedy_substitution",
    "extract_loci_sequences",
    "insertion_effect",
    "kmer_counts",
    "matched_background_loci",
    "motif_pair_spacing",
    "place_motif",
    "pwm_information_content",
    "saturation_mutagenesis",
    "substitution_effect",
]

_register_public_verbose(globals(), tuple(__all__))
