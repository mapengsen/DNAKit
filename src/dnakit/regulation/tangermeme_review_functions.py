"""Sequence-editing, random-sequence, PWM, coordinate, FASTA, k-mer, annotation and design functions
over tangermeme.

Each function runs one tangermeme routine on DNA strings, explicit matrices or small tables that the
caller supplies. Random functions take an explicit ``random_state``; no global seed is changed.
Functions that take a torch model require it in eval mode and run it on CPU.
"""

from __future__ import annotations

import math
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.regulation.backends import (
    _dna,
    _flag,
    _fraction,
    _integer,
    _nested,
    _torch_model,
    _version,
)

_TANGERMEME_EXTRA = "bio-tangermeme"
_BASES = "ACGT"
_MAX_LENGTH = 100_000
_MAX_SEQUENCES = 10_000
_MAX_COPIES = 1_000
_MAX_KMER_LENGTH = 8
_MAX_INDEX = 100_000
_MAX_TABLE_ROWS = 10_000
_MAX_CELLS = 1_000_000
_MAX_DESIGN_LENGTH = 1_000
_MAX_TARGET = 1_000
_MAX_MOTIFS = 5_000
_MAX_PISA_BASES = 100_000


def _call(function: Any, *args: Any, **kwargs: Any) -> Any:
    """Run an upstream call and convert its failures into DNAKit errors."""

    try:
        return function(*args, **kwargs)
    except KeyboardInterrupt:
        raise
    except Exception as exc:  # tangermeme and torch raise several exception types
        raise BackendExecutionError(f"tangermeme failed: {exc}") from exc


def _provenance(**notes: object) -> dict[str, Any]:
    return {"package": "tangermeme", "version": _version("tangermeme"), **notes}


def _sequence(value: object, name: str) -> str:
    return _dna(value, name, _MAX_LENGTH)


def _sequences(values: object, name: str) -> list[str]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list of DNA strings.")
    if not 1 <= len(values) <= _MAX_SEQUENCES:
        raise ConfigurationError(f"{name} must hold 1 to {_MAX_SEQUENCES} sequences.")
    return [_sequence(value, f"{name}[{index}]") for index, value in enumerate(values)]


def _span(start: object, end: object, length: int) -> tuple[int, int]:
    low = _integer(start, "start", 0, length - 1)
    high = _integer(end, "end", 1, length)
    if high <= low:
        raise ConfigurationError("end must be greater than start.")
    return low, high


def _probabilities(values: object) -> list[float]:
    if values is None:
        return [0.25, 0.25, 0.25, 0.25]
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or len(values) != 4:
        raise ConfigurationError("probs must hold four values for A, C, G and T.")
    numbers: list[float] = []
    for value in values:
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ConfigurationError("probs must contain finite non-negative numbers.")
        number = float(value)
        if not math.isfinite(number) or number < 0:
            raise ConfigurationError("probs must contain finite non-negative numbers.")
        numbers.append(number)
    total = sum(numbers)
    if total <= 0:
        raise ConfigurationError("probs must have a positive sum.")
    return [number / total for number in numbers]


def _encode(utils: Any, torch: Any, text: str) -> Any:
    return utils.one_hot_encode(text).float()


def _decode(tensor: Any) -> str:
    """Read a (4, L) one-hot tensor whose columns each hold exactly one base."""

    columns = tensor.sum(dim=0).tolist()
    if any(abs(value - 1.0) > 1e-6 for value in columns):
        raise BackendExecutionError("tangermeme returned a column that is not one-hot.")
    indices = tensor.argmax(dim=0).tolist()
    return "".join(_BASES[index] for index in indices)


def shuffle_region(
    sequence: str,
    start: int,
    end: int,
    *,
    n: int = 1,
    random_state: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Shuffle the bases between ``start`` and ``end`` of one sequence ``n`` times with tangermeme.

    Positions outside ``[start, end)`` are unchanged; the base composition inside the region is
    kept. Coordinates are 0-based and half-open.
    """

    ersatz = require("tangermeme.ersatz", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    torch = require("torch", _TANGERMEME_EXTRA)
    text = _sequence(sequence, "sequence")
    low, high = _span(start, end, len(text))
    copies = _integer(n, "n", 1, _MAX_COPIES)
    seed = None if random_state is None else _integer(random_state, "random_state", 0, 2**31 - 1)
    x = _encode(utils, torch, text)[None]
    out = _call(ersatz.shuffle, x, low, high, n=copies, random_state=seed)
    result = [_decode(out[0, index]) for index in range(copies)]
    return {
        "status": "completed",
        "operation": "shuffle_region",
        "start": low,
        "end": high,
        "n": copies,
        "sequences": result,
        "provenance": _provenance(function="tangermeme.ersatz.shuffle"),
    }


def randomize_region(
    sequence: str,
    start: int,
    end: int,
    *,
    probs: Sequence[float] | None = None,
    n: int = 1,
    random_state: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Replace the bases between ``start`` and ``end`` with bases drawn from ``probs``.

    ``probs`` gives the chance of A, C, G and T at every drawn position (uniform by default). Other
    positions are unchanged. The result holds ``n`` replacement sequences.
    """

    ersatz = require("tangermeme.ersatz", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    torch = require("torch", _TANGERMEME_EXTRA)
    text = _sequence(sequence, "sequence")
    low, high = _span(start, end, len(text))
    weights = _probabilities(probs)
    copies = _integer(n, "n", 1, _MAX_COPIES)
    seed = None if random_state is None else _integer(random_state, "random_state", 0, 2**31 - 1)
    x = _encode(utils, torch, text)[None]
    out = _call(ersatz.randomize, x, low, high, probs=[weights], n=copies, random_state=seed)
    result = [_decode(out[0, index]) for index in range(copies)]
    return {
        "status": "completed",
        "operation": "randomize_region",
        "start": low,
        "end": high,
        "probs": weights,
        "n": copies,
        "sequences": result,
        "provenance": _provenance(function="tangermeme.ersatz.randomize"),
    }


def delete_region(sequence: str, start: int, end: int, *, verbose: bool = True) -> dict[str, Any]:
    """Delete the bases between ``start`` and ``end`` with tangermeme ``delete``."""

    ersatz = require("tangermeme.ersatz", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    torch = require("torch", _TANGERMEME_EXTRA)
    text = _sequence(sequence, "sequence")
    low, high = _span(start, end, len(text))
    x = _encode(utils, torch, text)[None]
    out = _call(ersatz.delete, x, low, high)
    result = _decode(out[0])
    return {
        "status": "completed",
        "operation": "delete_region",
        "start": low,
        "end": high,
        "deleted": high - low,
        "sequence": result,
        "length": len(result),
        "provenance": _provenance(function="tangermeme.ersatz.delete"),
    }


def random_one_hot_sequences(
    count: int,
    length: int,
    *,
    probs: Sequence[float] | None = None,
    random_state: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw ``count`` random DNA sequences of ``length`` bases with tangermeme ``random_one_hot``.

    ``probs`` gives the chance of A, C, G and T at each position (uniform by default).
    """

    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    copies = _integer(count, "count", 1, _MAX_COPIES)
    size = _integer(length, "length", 1, _MAX_LENGTH)
    weights = _probabilities(probs)
    seed = None if random_state is None else _integer(random_state, "random_state", 0, 2**31 - 1)
    out = _call(utils.random_one_hot, (copies, 4, size), probs=[weights], random_state=seed)
    result = [_decode(out[index].float()) for index in range(copies)]
    return {
        "status": "completed",
        "operation": "random_one_hot_sequences",
        "count": copies,
        "length": size,
        "probs": weights,
        "sequences": result,
        "provenance": _provenance(function="tangermeme.utils.random_one_hot"),
    }


def place_motif_set(
    sequence: str,
    motifs: Sequence[str],
    spacing: Sequence[int],
    *,
    start: int,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write several motifs into one sequence at given spacings (tangermeme ``multisubstitute``).

    The first motif starts at ``start``; each later motif starts ``spacing[i]`` bases after the end
    of the previous one. Bases are substituted in place, so the length does not change.
    """

    space = require("tangermeme.ersatz", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    torch = require("torch", _TANGERMEME_EXTRA)
    text = _sequence(sequence, "sequence")
    if isinstance(motifs, (str, bytes)) or not isinstance(motifs, Sequence):
        raise ConfigurationError("motifs must be a list of DNA strings.")
    if not 2 <= len(motifs) <= 10:
        raise ConfigurationError("motifs must hold 2 to 10 DNA strings.")
    words = [_sequence(value, f"motifs[{index}]") for index, value in enumerate(motifs)]
    if any(len(word) > 50 for word in words):
        raise ConfigurationError("each motif must be at most 50 bases.")
    if isinstance(spacing, (str, bytes)) or not isinstance(spacing, Sequence):
        raise ConfigurationError("spacing must be a list of integers.")
    if len(spacing) != len(words) - 1:
        raise ConfigurationError("spacing needs one value between each pair of motifs.")
    gaps = [_integer(value, "spacing", 0, 1_000) for value in spacing]
    position = _integer(start, "start", 0, len(text) - 1)
    x = _encode(utils, torch, text)[None]
    out = _call(space.multisubstitute, x, words, gaps, start=position)
    result = _decode(out[0])
    return {
        "status": "completed",
        "operation": "place_motif_set",
        "start": position,
        "spacing": gaps,
        "sequence": result,
        "provenance": _provenance(function="tangermeme.ersatz.multisubstitute"),
    }


def local_dinucleotide_shuffle(
    sequences: Sequence[str],
    *,
    n: int = 1,
    bin_size: int = 1024,
    min_bin_size: int = 512,
    random_state: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Dinucleotide-shuffle equal-length sequences inside local bins with tangermeme.

    Each sequence is cut into bins of ``bin_size`` bases (a short last bin is merged into the
    previous
    one if it is below ``min_bin_size``) and shuffled within each bin, keeping dinucleotide counts
    per bin. Results list ``n`` shuffles for every input sequence.
    """

    ersatz = require("tangermeme.ersatz", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    torch = require("torch", _TANGERMEME_EXTRA)
    texts = _sequences(sequences, "sequences")
    lengths = {len(text) for text in texts}
    if len(lengths) != 1:
        raise ConfigurationError("all sequences must have the same length.")
    copies = _integer(n, "n", 1, _MAX_COPIES)
    size = _integer(bin_size, "bin_size", 2, _MAX_LENGTH)
    minimum = _integer(min_bin_size, "min_bin_size", 1, size)
    seed = None if random_state is None else _integer(random_state, "random_state", 0, 2**31 - 1)
    x = torch.stack([_encode(utils, torch, text) for text in texts])
    out = _call(
        ersatz.local_dinucleotide_shuffle,
        x,
        n=copies,
        bin_size=size,
        min_bin_size=minimum,
        random_state=seed,
        verbose=False,
    )
    grouped = [
        [_decode(out[row, index].float()) for index in range(copies)] for row in range(len(texts))
    ]
    return {
        "status": "completed",
        "operation": "local_dinucleotide_shuffle",
        "input_count": len(texts),
        "n": copies,
        "bin_size": size,
        "min_bin_size": minimum,
        "sequences": grouped,
        "provenance": _provenance(function="tangermeme.ersatz.local_dinucleotide_shuffle"),
    }


def pwm_consensus_sequence(
    pwm: Sequence[Sequence[float]], *, verbose: bool = True
) -> dict[str, Any]:
    """Consensus DNA string of a position weight matrix with tangermeme ``pwm_consensus``.

    ``pwm`` has one row per position and four columns (A, C, G, T); the consensus base is the column
    with the largest value. Rows that are all zero are rejected.
    """

    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    torch = require("torch", _TANGERMEME_EXTRA)
    if (
        isinstance(pwm, (str, bytes))
        or not isinstance(pwm, Sequence)
        or not 1 <= len(pwm) <= 10_000
    ):
        raise ConfigurationError("pwm must be a list of 1 to 10000 rows of four numbers.")
    rows: list[list[float]] = []
    for index, row in enumerate(pwm):
        if isinstance(row, (str, bytes)) or not isinstance(row, Sequence) or len(row) != 4:
            raise ConfigurationError(f"pwm row {index} must hold four numbers.")
        values = []
        for value in row:
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or value < 0
            ):
                raise ConfigurationError("pwm values must be finite non-negative numbers.")
            values.append(float(value))
        if sum(values) <= 0:
            raise ConfigurationError(f"pwm row {index} must not be all zero.")
        rows.append(values)
    tensor = torch.tensor(rows, dtype=torch.float32)
    consensus = _call(utils.pwm_consensus, tensor)
    letters = "".join(_BASES[index] for index in consensus.argmax(dim=1).tolist())
    return {
        "status": "completed",
        "operation": "pwm_consensus_sequence",
        "length": len(rows),
        "consensus": letters,
        "provenance": _provenance(function="tangermeme.utils.pwm_consensus"),
    }


def example_to_sequence_coords(
    examples: Sequence[Mapping[str, Any]],
    loci: Sequence[Mapping[str, Any]],
    *,
    window: int | None = None,
    one_indexed: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert coordinates inside each example to genome coordinates with tangermeme.

    ``examples`` rows hold ``example_idx``, ``start`` and ``end`` relative to the example. ``loci``
    rows hold ``chrom``, ``start`` and ``end`` of the locus that each example index refers to.
    Output
    rows hold ``chrom``, ``start``, ``end`` and ``example_idx``, sorted by genome position.
    """

    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    pandas = require("pandas", _TANGERMEME_EXTRA)
    if isinstance(examples, (str, bytes)) or not isinstance(examples, Sequence) or not examples:
        raise ConfigurationError("examples must be a non-empty list of rows.")
    if isinstance(loci, (str, bytes)) or not isinstance(loci, Sequence) or not loci:
        raise ConfigurationError("loci must be a non-empty list of rows.")
    if len(examples) > _MAX_SEQUENCES or len(loci) > _MAX_SEQUENCES:
        raise ConfigurationError(f"tables may hold at most {_MAX_SEQUENCES} rows.")
    example_rows = []
    for row in examples:
        example_rows.append(
            {
                "example_idx": _integer(row["example_idx"], "example_idx", 0, _MAX_SEQUENCES),
                "start": _integer(row["start"], "start", 0, _MAX_LENGTH),
                "end": _integer(row["end"], "end", 1, _MAX_LENGTH),
            }
        )
    locus_rows = []
    for row in loci:
        locus_rows.append(
            {
                "chrom": str(row["chrom"]),
                "start": _integer(row["start"], "start", 0, 10**12),
                "end": _integer(row["end"], "end", 1, 10**12),
            }
        )
    width = None if window is None else _integer(window, "window", 1, _MAX_LENGTH)
    if not isinstance(one_indexed, bool):
        raise ConfigurationError("one_indexed must be True or False.")
    frame = _call(
        utils.example_to_fasta_coords,
        pandas.DataFrame(example_rows),
        pandas.DataFrame(locus_rows),
        window=width,
        one_indexed=one_indexed,
    )
    result = [
        {
            "chrom": str(record["chrom"]),
            "start": int(record["start"]),
            "end": int(record["end"]),
        }
        for record in frame.to_dict("records")
    ]
    result.sort(key=lambda item: (item["chrom"], item["start"], item["end"]))
    return {
        "status": "completed",
        "operation": "example_to_sequence_coords",
        "window": width,
        "one_indexed": one_indexed,
        "row_count": len(result),
        "intervals": result,
        "provenance": _provenance(function="tangermeme.utils.example_to_fasta_coords"),
    }


def one_hot_to_fasta_file(
    sequences: Sequence[str],
    output: str | Path,
    *,
    headers: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write DNA strings to a FASTA file through tangermeme ``one_hot_to_fasta``.

    Headers default to ``seq0``, ``seq1``, and so on. The output file must not exist yet.
    """

    io_module = require("tangermeme.io", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    torch = require("torch", _TANGERMEME_EXTRA)
    texts = _sequences(sequences, "sequences")
    lengths = {len(text) for text in texts}
    if len(lengths) != 1:
        raise ConfigurationError("all sequences must have the same length.")
    path = Path(output).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"Output file already exists: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {path.parent}")
    names = [f"seq{index}" for index in range(len(texts))] if headers is None else list(headers)
    if len(names) != len(texts) or any(not str(name).strip() for name in names):
        raise ConfigurationError("headers must hold one non-empty name per sequence.")
    x = torch.stack([_encode(utils, torch, text) for text in texts])
    _call(io_module.one_hot_to_fasta, x, str(path), headers=[str(name) for name in names])
    if not path.is_file():
        raise BackendExecutionError("tangermeme did not write the FASTA file.")
    return {
        "status": "completed",
        "operation": "one_hot_to_fasta_file",
        "output": path.name,
        "record_count": len(texts),
        "headers": [str(name) for name in names],
        "provenance": _provenance(function="tangermeme.io.one_hot_to_fasta"),
    }


def gapped_kmer_features(
    sequences: Sequence[str],
    *,
    min_k: int = 2,
    max_k: int = 2,
    max_gap: int = 1,
    max_len: int = 3,
    max_gkmers: int = 5,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count gapped k-mers per sequence with tangermeme ``gapped_kmers``.

    Sequences must be longer than the window of the largest k-mer (tangermeme fails on short
    inputs). Returns the non-zero counts as ``[sequence_index, feature_index, count]`` triples.
    Feature
    indices
    refer to the columns that tangermeme creates for the given k, gap and length limits.
    """

    kmers = require("tangermeme.kmers", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    torch = require("torch", _TANGERMEME_EXTRA)
    texts = _sequences(sequences, "sequences")
    lengths = {len(text) for text in texts}
    if len(lengths) != 1:
        raise ConfigurationError("all sequences must have the same length.")
    low = _integer(min_k, "min_k", 1, _MAX_KMER_LENGTH)
    high = _integer(max_k, "max_k", low, _MAX_KMER_LENGTH)
    gap = _integer(max_gap, "max_gap", 0, 10)
    span = _integer(max_len, "max_len", high, 20)
    limit = _integer(max_gkmers, "max_gkmers", 1, 1_000)
    x = torch.stack([_encode(utils, torch, text) for text in texts])
    matrix = _call(
        kmers.gapped_kmers,
        x,
        min_k=low,
        max_k=high,
        max_gap=gap,
        max_len=span,
        max_gkmers=limit,
    ).tocoo()
    triples = sorted(
        (
            [int(r), int(c), int(v)]
            for r, c, v in zip(matrix.row, matrix.col, matrix.data, strict=True)
        ),
    )
    if any(value <= 0 for _, _, value in triples):
        raise BackendExecutionError("tangermeme returned a non-positive k-mer count.")
    return {
        "status": "completed",
        "operation": "gapped_kmer_features",
        "sequence_count": len(texts),
        "feature_count": int(matrix.shape[1]),
        "min_k": low,
        "max_k": high,
        "max_gap": gap,
        "max_len": span,
        "counts": triples,
        "provenance": _provenance(function="tangermeme.kmers.gapped_kmers"),
    }


_SPAN_KEYS = ("example_idx", "start", "end")


def _seed(value: object) -> int | None:
    return None if value is None else _integer(value, "random_state", 0, 2**31 - 1)


def _numbers(value: object, name: str, low: int, high: int) -> list[float]:
    """Finite numbers in a list whose length lies in ``[low, high]``."""

    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must hold {low} to {high} finite numbers.")
    if not low <= len(value) <= high:
        raise ConfigurationError(f"{name} must hold {low} to {high} finite numbers.")
    numbers: list[float] = []
    for item in value:
        if isinstance(item, bool) or not isinstance(item, (int, float)) or not math.isfinite(item):
            raise ConfigurationError(f"{name} must hold finite numbers only.")
        numbers.append(float(item))
    return numbers


def _matrix(value: object, name: str) -> list[list[float]]:
    """Rows of finite numbers that all have the same length."""

    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list of rows of numbers.")
    if not 1 <= len(value) <= _MAX_TABLE_ROWS:
        raise ConfigurationError(f"{name} must hold 1 to {_MAX_TABLE_ROWS} rows.")
    rows = [_numbers(row, f"{name}[{index}]", 1, _MAX_LENGTH) for index, row in enumerate(value)]
    if len({len(row) for row in rows}) != 1:
        raise ConfigurationError(f"all {name} rows must have the same length.")
    _check_cells(len(rows) * len(rows[0]))
    return rows


def _check_cells(count: int) -> None:
    if count > _MAX_CELLS:
        raise ConfigurationError(f"the result would hold {count} cells; the limit is {_MAX_CELLS}.")


def _table(rows: object, name: str, keys: tuple[str, ...]) -> list[list[int]]:
    """Rows of mappings with ``keys`` as integer lists, in key order."""

    if isinstance(rows, (str, bytes)) or not isinstance(rows, Sequence):
        raise ConfigurationError(f"{name} must be a list of rows.")
    if not 1 <= len(rows) <= _MAX_TABLE_ROWS:
        raise ConfigurationError(f"{name} must hold 1 to {_MAX_TABLE_ROWS} rows.")
    table: list[list[int]] = []
    for index, row in enumerate(rows):
        if isinstance(row, Mapping) and all(key in row for key in keys):
            values: list[object] = [row[key] for key in keys]
        elif (
            isinstance(row, Sequence)
            and not isinstance(row, (str, bytes))
            and len(row) == len(keys)
        ):
            values = list(row)
        else:
            raise ConfigurationError(f"{name}[{index}] must hold {', '.join(keys)}.")
        table.append(
            [
                _integer(value, f"{name}[{index}].{key}", 0, _MAX_INDEX)
                for key, value in zip(keys, values, strict=True)
            ]
        )
    return table


def _spans(table: list[list[int]], texts: Sequence[str], name: str) -> list[list[int]]:
    """Check ``example_idx``, ``start`` and ``end`` rows against the sequences they name."""

    for index, (example, start, end) in enumerate(table):
        if example >= len(texts):
            raise ConfigurationError(f"{name}[{index}].example_idx must be below {len(texts)}.")
        _span(start, end, len(texts[example]))
    return table


def _common_length(texts: Sequence[str], name: str) -> int:
    lengths = {len(text) for text in texts}
    if len(lengths) != 1:
        raise ConfigurationError(f"all {name} must have the same length.")
    return lengths.pop()


def _stack(utils: Any, torch: Any, texts: Sequence[str]) -> Any:
    return torch.stack([_encode(utils, torch, text) for text in texts])


def _plain(tensor: Any) -> Any:
    """Nested lists for one tensor; non-finite values become None."""

    return _nested(tensor.detach().cpu().tolist())


def _tensor_output(value: Any, torch: Any, name: str) -> Any:
    if not isinstance(value, torch.Tensor):
        raise ConfigurationError(f"{name}: models with several output tensors are not supported.")
    return value


def _target(value: object, torch: Any) -> Any:
    """A target vector as a (1, n_outputs) float tensor."""

    return torch.tensor([_numbers(value, "target", 1, _MAX_TARGET)], dtype=torch.float32)


def _motif_words(motifs: object) -> list[str]:
    if isinstance(motifs, (str, bytes)) or not isinstance(motifs, Sequence):
        raise ConfigurationError("motifs must be a list of DNA strings.")
    if not 1 <= len(motifs) <= 20:
        raise ConfigurationError("motifs must hold 1 to 20 DNA strings.")
    words = [_dna(value, f"motifs[{index}]", 50) for index, value in enumerate(motifs)]
    if any(not word for word in words):
        raise ConfigurationError("motifs must not be empty strings.")
    return words


def _pwm_rows(value: object, name: str) -> list[list[float]]:
    """PWM rows of four non-negative numbers, each row scaled to sum to one."""

    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list of rows of four numbers.")
    if not 1 <= len(value) <= 1000:
        raise ConfigurationError(f"{name} must hold 1 to 1000 rows.")
    rows: list[list[float]] = []
    for index, row in enumerate(value):
        values = _numbers(row, f"{name} row {index}", 4, 4)
        if any(item < 0 for item in values):
            raise ConfigurationError(f"{name} values must be non-negative.")
        total = sum(values)
        if total <= 0:
            raise ConfigurationError(f"{name} row {index} must not be all zero.")
        rows.append([item / total for item in values])
    return rows


def _motif_tensors(motifs: object, torch: Any) -> tuple[list[str], list[Any]]:
    if not isinstance(motifs, Mapping) or not 1 <= len(motifs) <= _MAX_MOTIFS:
        raise ConfigurationError(f"motifs must map 1 to {_MAX_MOTIFS} names to PWM rows.")
    names: list[str] = []
    pwms: list[Any] = []
    for name, pwm in motifs.items():
        if not isinstance(name, str) or not name.strip():
            raise ConfigurationError("motif names must be non-empty strings.")
        rows = _pwm_rows(pwm, f"motifs[{name!r}]")
        names.append(name)
        # tangermeme expects (alphabet, length); the caller gives one row per position.
        pwms.append(torch.tensor(rows, dtype=torch.float32).T.contiguous())
    return names, pwms


def _axis(value: object) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or value not in (0, 1):
        raise ConfigurationError("dim must be None, 0 or 1.")
    return value


def _gapped(tensor: Any) -> str:
    """Read a (4, length) one-hot tensor whose all-zero columns are gaps, written as N."""

    columns = tensor.sum(dim=0).tolist()
    if any(abs(value - 1.0) > 1e-6 and abs(value) > 1e-6 for value in columns):
        raise BackendExecutionError(
            "tangermeme returned a column that is neither one-hot nor empty."
        )
    indices = tensor.argmax(dim=0).tolist()
    return "".join(
        _BASES[index] if value > 0.5 else "N" for index, value in zip(indices, columns, strict=True)
    )


@contextmanager
def _numpy_state_kept() -> Iterator[None]:
    """Restore NumPy's global random state after an upstream call that may draw from it."""

    numpy = require("numpy", _TANGERMEME_EXTRA)
    saved = numpy.random.get_state()
    try:
        yield
    finally:
        numpy.random.set_state(saved)


def ablate_annotations(
    model: Any,
    sequences: Sequence[str],
    annotations: Sequence[Mapping[str, Any]],
    *,
    n: int = 1,
    random_state: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Shuffle each annotated span on its own and record the model before and after with tangermeme.

    ``sequences`` are equal-length DNA strings. Each annotation holds ``example_idx``, ``start``
    and ``end`` (0-based, end exclusive). ``n`` shuffles are drawn per annotation and
    ``random_state`` makes them repeatable. Shuffling keeps the base counts inside the span, so a
    model that only counts bases returns the same value after ablation.
    """

    ablate = require("tangermeme.ablate", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    torch = _torch_model(model)
    texts = _sequences(sequences, "sequences")
    _common_length(texts, "sequences")
    table = _spans(_table(annotations, "annotations", _SPAN_KEYS), texts, "annotations")
    copies = _integer(n, "n", 1, _MAX_COPIES)
    seed = _seed(random_state)
    result = _call(
        ablate.ablate_annotations,
        model,
        _stack(utils, torch, texts),
        torch.tensor(table, dtype=torch.int64),
        n=copies,
        random_state=seed,
    )
    before = _tensor_output(result.y_befores, torch, "y_befores")
    after = _tensor_output(result.y_afters, torch, "y_afters")
    return {
        "status": "completed",
        "operation": "ablate_annotations",
        "annotation_count": len(table),
        "n": copies,
        "before": _plain(before[:, 0]),
        "after": _plain(after[:, 0]),
        "provenance": _provenance(function="tangermeme.ablate.ablate_annotations"),
    }


def annotate_seqlets(
    sequences: Sequence[str],
    seqlets: Sequence[Mapping[str, Any]],
    motifs: Mapping[str, Sequence[Sequence[float]]],
    *,
    n_nearest: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Match each seqlet to its best motifs with TOMTOM in tangermeme.

    ``sequences`` are equal-length DNA strings and each seqlet holds ``example_idx``, ``start`` and
    ``end``. ``motifs`` maps names to PWMs with one row per position and four columns (A, C, G,
    T); rows are scaled to sum to one. The result lists the ``n_nearest`` motifs for each seqlet,
    most significant first, with TOMTOM p-values. Use at least 25 motifs: tangermeme warns below
    that size because the p-values depend on the set of targets.
    """

    annotate = require("tangermeme.annotate", _TANGERMEME_EXTRA)
    pandas = require("pandas", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    torch = require("torch", _TANGERMEME_EXTRA)
    texts = _sequences(sequences, "sequences")
    _common_length(texts, "sequences")
    table = _spans(_table(seqlets, "seqlets", _SPAN_KEYS), texts, "seqlets")
    names, pwms = _motif_tensors(motifs, torch)
    nearest = _integer(n_nearest, "n_nearest", 1, len(names))
    frame = pandas.DataFrame(table, columns=list(_SPAN_KEYS))
    idxs, pvals = _call(
        annotate.annotate_seqlets,
        _stack(utils, torch, texts),
        frame,
        dict(zip(names, pwms, strict=True)),
        n_nearest=nearest,
        n_jobs=1,
    )
    matches = [
        [
            {"motif": names[hit] if hit >= 0 else None, "p_value": float(score)}
            for hit, score in zip(hits, scores, strict=True)
        ]
        for hits, scores in zip(idxs.tolist(), pvals.tolist(), strict=True)
    ]
    return {
        "status": "completed",
        "operation": "annotate_seqlets",
        "seqlet_count": len(table),
        "motif_count": len(names),
        "n_nearest": nearest,
        "matches": matches,
        "provenance": _provenance(function="tangermeme.annotate.annotate_seqlets"),
    }


def count_seqlet_annotations(
    pairs: Sequence[Sequence[int]],
    *,
    dim: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count how often each annotation occurs in each example with tangermeme.

    ``pairs`` holds ``[example_idx, annotation_idx]`` rows. ``dim`` None returns the matrix of
    examples by annotations; 0 sums over examples and gives one count per annotation; 1 sums over
    annotations and gives one count per example. Counts are int64, so they cannot overflow the
    uint8 default of tangermeme.
    """

    annotate = require("tangermeme.annotate", _TANGERMEME_EXTRA)
    torch = require("torch", _TANGERMEME_EXTRA)
    table = _table(pairs, "pairs", ("example_idx", "annotation_idx"))
    axis = _axis(dim)
    _check_cells((max(row[0] for row in table) + 1) * (max(row[1] for row in table) + 1))
    counts = _call(
        annotate.count_annotations,
        torch.tensor(table, dtype=torch.int64),
        dtype=torch.int64,
        dim=axis,
    )
    return {
        "status": "completed",
        "operation": "count_seqlet_annotations",
        "dim": axis,
        "counts": _plain(counts),
        "provenance": _provenance(function="tangermeme.annotate.count_annotations"),
    }


def pairwise_annotation_counts(
    annotations: Sequence[Sequence[int] | Mapping[str, Any]],
    *,
    max_distance: int | None = None,
    symmetric: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count annotation pairs that share an example, or the spacing between them, with tangermeme.

    With ``max_distance`` None, ``annotations`` holds ``[example_idx, annotation_idx]`` rows and the
    result is an annotation by annotation matrix; each distinct pair counts once per example
    (tangermeme ``pairwise_annotations``). With an integer ``max_distance``, rows are mappings with
    ``example_idx``, ``annotation``, ``start`` and ``end``, and the result is an annotation by
    annotation by distance cube whose last axis counts gaps from 0 to ``max_distance`` minus 1
    (tangermeme ``pairwise_annotations_spacing``). Spans within one example must not overlap.
    ``symmetric`` False keeps the row order of the input (pairs) or the left-to-right order
    (spacing) instead of mirroring each pair.
    """

    annotate = require("tangermeme.annotate", _TANGERMEME_EXTRA)
    torch = require("torch", _TANGERMEME_EXTRA)
    both = _flag(symmetric, "symmetric")
    if max_distance is None:
        table = _table(annotations, "annotations", ("example_idx", "annotation_idx"))
        width = max(row[1] for row in table) + 1
        _check_cells(width * width)
        counts = _call(
            annotate.pairwise_annotations,
            torch.tensor(table, dtype=torch.int64),
            dtype=torch.int64,
            symmetric=both,
        )
        mode, spacing, function = "pairs", None, "tangermeme.annotate.pairwise_annotations"
    else:
        span = _integer(max_distance, "max_distance", 1, _MAX_LENGTH)
        table = _table(annotations, "annotations", ("example_idx", "annotation", "start", "end"))
        if any(row[3] <= row[2] for row in table):
            raise ConfigurationError("each annotation end must be greater than its start.")
        width = max(row[1] for row in table) + 1
        _check_cells(width * width * span)
        counts = _call(
            annotate.pairwise_annotations_spacing,
            torch.tensor(table, dtype=torch.int64),
            max_distance=span,
            dtype=torch.int64,
            symmetric=both,
        )
        mode, spacing = "spacing", span
        function = "tangermeme.annotate.pairwise_annotations_spacing"
    return {
        "status": "completed",
        "operation": "pairwise_annotation_counts",
        "mode": mode,
        "max_distance": spacing,
        "symmetric": both,
        "row_count": len(table),
        "counts": _plain(counts),
        "provenance": _provenance(function=function),
    }


def screen_sequence_design(
    model: Any,
    length: int,
    target: Sequence[float],
    *,
    n_best: int = 1,
    tol: float = 1e-3,
    max_iter: int = 100,
    batch_size: int = 32,
    random_state: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find random DNA sequences whose model output is closest to ``target`` with tangermeme.

    Random one-hot sequences of ``length`` bases are scored in batches of ``batch_size`` and the
    ``n_best`` with the smallest squared error are kept. The search stops once the worst kept
    candidate is within ``tol`` of the target, or after ``max_iter`` batches, so it always ends.
    ``target`` holds one number per model output. With ``n_best`` above one, tangermeme can fail
    when two candidates tie exactly on loss, because it compares tensors inside a heap; keep
    ``n_best`` at one for models with discrete outputs.
    """

    design = require("tangermeme.design", _TANGERMEME_EXTRA)
    torch = _torch_model(model)
    size = _integer(length, "length", 1, _MAX_DESIGN_LENGTH)
    goal = _target(target, torch)
    best = _integer(n_best, "n_best", 1, 16)
    if (
        isinstance(tol, bool)
        or not isinstance(tol, (int, float))
        or not math.isfinite(tol)
        or tol < 0
    ):
        raise ConfigurationError("tol must be a finite non-negative number.")
    rounds = _integer(max_iter, "max_iter", 1, 10_000)
    batch = _integer(batch_size, "batch_size", 1, 1024)
    seed = _seed(random_state)
    found = _call(
        design.screen,
        model,
        (4, size),
        y=goal,
        n_best=best,
        tol=float(tol),
        max_iter=rounds,
        batch_size=batch,
        random_state=seed,
        device="cpu",
    )
    sequences = [_decode(found[index]) for index in range(best)]
    with torch.no_grad():
        predicted = model(found)
    if not isinstance(predicted, torch.Tensor) or predicted.ndim != 2:
        raise ConfigurationError("model must return one (batch, n_outputs) tensor.")
    if predicted.shape[1] != goal.shape[1]:
        raise ConfigurationError("target must hold one number per model output.")
    squared = ((predicted - goal) ** 2).sum(dim=1)
    return {
        "status": "completed",
        "operation": "screen_sequence_design",
        "length": size,
        "n_best": best,
        "sequences": sequences,
        "predictions": _plain(predicted),
        "squared_error": _plain(squared),
        "provenance": _provenance(function="tangermeme.design.screen"),
    }


def greedy_marginalize_design(
    model: Any,
    sequence: str,
    target: Sequence[float],
    motifs: Sequence[str],
    *,
    max_spacing: int = 12,
    reverse_complement: bool = True,
    max_iter: int = 10,
    batch_size: int = 32,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a construct by greedily placing motifs that move the model output toward ``target``.

    The first motif starts at position ``len(sequence) // 2 - 1``. Each later motif is searched
    within ``max_spacing`` bases of the motifs already placed. A motif is kept only while it lowers
    the mean squared error by more than 0.001, up to ``max_iter`` motifs. Reverse complements of
    the motifs are also tried when ``reverse_complement`` is True. Bases that no motif covers are
    written as ``N`` and the ends are trimmed. ``sequence`` must be at least
    ``2 * (max_spacing + longest motif)`` bases long.
    """

    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    design = require("tangermeme.design", _TANGERMEME_EXTRA)
    torch = _torch_model(model)
    text = _sequence(sequence, "sequence")
    goal = _target(target, torch)
    words = _motif_words(motifs)
    spacing = _integer(max_spacing, "max_spacing", 0, 100)
    flip = _flag(reverse_complement, "reverse_complement")
    rounds = _integer(max_iter, "max_iter", 1, 100)
    batch = _integer(batch_size, "batch_size", 1, 1024)
    longest = max(len(word) for word in words)
    if len(text) < 2 * (spacing + longest):
        raise ConfigurationError(
            "sequence must be at least 2 * (max_spacing + longest motif) bases long."
        )
    built = _call(
        design.greedy_marginalize,
        model,
        _encode(utils, torch, text)[None],
        goal,
        words,
        max_spacing=spacing,
        reverse_complement=flip,
        max_iter=rounds,
        batch_size=batch,
        device="cpu",
    )
    construct = _gapped(built)
    return {
        "status": "completed",
        "operation": "greedy_marginalize_design",
        "max_spacing": spacing,
        "reverse_complement": flip,
        "sequence": construct,
        "length": len(construct),
        "provenance": _provenance(function="tangermeme.design.greedy_marginalize"),
    }


def marginalize_annotations(
    model: Any,
    sequences: Sequence[str],
    background: Sequence[str],
    annotations: Sequence[Mapping[str, Any]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Insert each annotated span into every background sequence and record the model output.

    ``sequences`` and ``background`` are lists of equal-length DNA strings. Each annotation holds
    ``example_idx``, ``start`` and ``end`` and names a span of ``sequences``. The span replaces
    the bases of each background sequence at ``insert_start``, which is the centre of the
    background (``len // 2 - span // 2``), as in tangermeme. The results are indexed by
    annotation and then by background sequence.
    """

    marginal = require("tangermeme.marginalize", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    torch = _torch_model(model)
    texts = _sequences(sequences, "sequences")
    _common_length(texts, "sequences")
    grounds = _sequences(background, "background")
    width = _common_length(grounds, "background")
    table = _spans(_table(annotations, "annotations", _SPAN_KEYS), texts, "annotations")
    for index, (_example, start, end) in enumerate(table):
        if end - start > width:
            raise ConfigurationError(
                f"annotations[{index}] is longer than the background sequences."
            )
    result = _call(
        marginal.marginalize_annotations,
        model,
        _stack(utils, torch, texts),
        _stack(utils, torch, grounds),
        torch.tensor(table, dtype=torch.int64),
    )
    return {
        "status": "completed",
        "operation": "marginalize_annotations",
        "annotation_count": len(table),
        "background_count": len(grounds),
        "insert_start": [width // 2 - (end - start) // 2 for _example, start, end in table],
        "before": _plain(_tensor_output(result.y_befores, torch, "y_befores")),
        "after": _plain(_tensor_output(result.y_afters, torch, "y_afters")),
        "provenance": _provenance(function="tangermeme.marginalize.marginalize_annotations"),
    }


def pisa_attribution(
    model: Any,
    sequences: Sequence[str],
    *,
    n_shuffles: int = 20,
    hypothetical: bool = False,
    batch_size: int = 32,
    random_state: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Attribute each model output to the bases of each sequence with PISA in tangermeme.

    PISA averages gradient-times-input attributions against ``n_shuffles`` dinucleotide shuffles
    of each sequence. The result has the shape sequences by outputs by 4 by length. With
    ``hypothetical`` False the attributions are multiplied by the observed one-hot bases, so only
    the base present at each position is non-zero. The model runs on CPU, and NumPy's global
    random state is restored after the call.
    """

    pisa = require("tangermeme.pisa", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    torch = _torch_model(model)
    texts = _sequences(sequences, "sequences")
    width = _common_length(texts, "sequences")
    if len(texts) * width > _MAX_PISA_BASES:
        raise ConfigurationError(
            f"sequences may hold at most {_MAX_PISA_BASES} bases in total for PISA."
        )
    shuffles = _integer(n_shuffles, "n_shuffles", 1, 200)
    flag = _flag(hypothetical, "hypothetical")
    batch = _integer(batch_size, "batch_size", 1, 1024)
    seed = _seed(random_state)
    with _numpy_state_kept():
        attributions = _call(
            pisa.pisa,
            model,
            _stack(utils, torch, texts),
            batch_size=batch,
            n_shuffles=shuffles,
            hypothetical=flag,
            random_state=seed,
            device="cpu",
            verbose=False,
        )
    return {
        "status": "completed",
        "operation": "pisa_attribution",
        "sequence_count": len(texts),
        "n_shuffles": shuffles,
        "hypothetical": flag,
        "shape": list(attributions.shape),
        "attributions": _plain(attributions),
        "provenance": _provenance(function="tangermeme.pisa.pisa"),
    }


def tfmodisco_seqlets(
    attributions: Sequence[Sequence[float]],
    *,
    window_size: int = 21,
    flank: int = 10,
    target_fdr: float = 0.2,
    min_passing_frac: float = 0.03,
    max_passing_frac: float = 0.2,
    weak_threshold_for_counting_sign: float = 0.8,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract seqlets from attribution scores with the TF-MoDISco procedure in tangermeme.

    ``attributions`` has one row per example and one column per position. Windows of
    ``window_size`` positions are summed, a threshold is chosen from a null model to reach about
    ``target_fdr``, and each seqlet keeps the window plus ``flank`` positions on both sides. The
    attribution of a seqlet is the sum over the seqlet after trimming
    ``(length - window_size) // 2`` positions from each end. Coordinates are 0-based with the end
    exclusive, and the procedure is deterministic.
    """

    seqlet = require("tangermeme.seqlet", _TANGERMEME_EXTRA)
    torch = require("torch", _TANGERMEME_EXTRA)
    rows = _matrix(attributions, "attributions")
    window = _integer(window_size, "window_size", 1, 1000)
    if len(rows[0]) < window:
        raise ConfigurationError("each attribution row must be at least window_size long.")
    flanking = _integer(flank, "flank", 0, 1000)
    frame = _call(
        seqlet.tfmodisco_seqlets,
        torch.tensor(rows, dtype=torch.float32),
        window_size=window,
        flank=flanking,
        target_fdr=_fraction(target_fdr, "target_fdr"),
        min_passing_frac=_fraction(min_passing_frac, "min_passing_frac"),
        max_passing_frac=_fraction(max_passing_frac, "max_passing_frac"),
        weak_threshold_for_counting_sign=_fraction(
            weak_threshold_for_counting_sign, "weak_threshold_for_counting_sign"
        ),
    )
    found = [
        {
            "example_idx": int(row["example_idx"]),
            "start": int(row["start"]),
            "end": int(row["end"]),
            "attribution": float(row["attribution"]),
        }
        for row in frame.to_dict("records")
    ]
    return {
        "status": "completed",
        "operation": "tfmodisco_seqlets",
        "seqlet_count": len(found),
        "seqlets": found,
        "provenance": _provenance(function="tangermeme.seqlet.tfmodisco_seqlets"),
    }


__all__ = [
    "ablate_annotations",
    "annotate_seqlets",
    "count_seqlet_annotations",
    "delete_region",
    "example_to_sequence_coords",
    "gapped_kmer_features",
    "greedy_marginalize_design",
    "local_dinucleotide_shuffle",
    "marginalize_annotations",
    "one_hot_to_fasta_file",
    "pairwise_annotation_counts",
    "pisa_attribution",
    "place_motif_set",
    "pwm_consensus_sequence",
    "random_one_hot_sequences",
    "randomize_region",
    "screen_sequence_design",
    "shuffle_region",
    "tfmodisco_seqlets",
]

_register_public_verbose(globals(), tuple(__all__))
