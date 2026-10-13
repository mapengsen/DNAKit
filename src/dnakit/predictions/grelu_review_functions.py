"""Review group G18d for gReLU: offline DNA intervals, motifs, attributions and model experiments.

Each public function imports ``grelu`` (and torch, numpy or pandas where needed) only when it
runs. Model-based functions take a ``grelu.lightning.LightningModel`` object, which can be
loaded from a local checkpoint with ``LightningModel.load_from_checkpoint``. Genome-based
functions take a local FASTA path that has a ``.fai`` index and a ``.sizes`` file. Nothing
here downloads model weights or genomes.

Upstream calls that reseed NumPy, torch or Python random state run inside a context that
restores those states afterwards, so the wrappers leave global state unchanged.
"""

from __future__ import annotations

import contextlib
import random
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError
from dnakit.predictions.dna_model_backends import (
    _choice,
    _clean,
    _dna,
    _dna_batch,
    _integer,
    _provenance,
    _real,
    _rejected,
    _text,
)
from dnakit.predictions.grelu_functions import (
    _interval_frame,
    _motif_set,
    _quiet,
    _rows,
    _sequence_items,
)

_GRELU = "grelu"
_EXTRAS = "bio-grelu"
_UPSTREAM_ERRORS = (
    AssertionError,
    AttributeError,
    IndexError,
    KeyError,
    OSError,
    RuntimeError,
    TypeError,
    ValueError,
)
_MAX_ROWS = 100_000
_MAX_HITS = 100_000
_MAX_SEQ = 100_000
_MAX_COORD = 10**9
_MAX_NAME = 256
_MAX_SEEDS = 2**32 - 1
_TASK_AGG = ("mean", "sum", "max", "min")
_ATTRIBUTION_METHODS = ("deepshap", "saliency", "inputxgradient", "integratedgradients")
_CHROM_SHORTCUTS = ("autosomes", "autosomesX", "autosomesXY")


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


@contextlib.contextmanager
def _rng_preserved() -> Iterator[None]:
    """Restore NumPy, torch and Python random state after an upstream call."""

    np = require("numpy", _EXTRAS)
    torch = require("torch", _EXTRAS)
    numpy_state = np.random.get_state()
    torch_state = torch.random.get_rng_state()
    python_state = random.getstate()
    try:
        yield
    finally:
        np.random.set_state(numpy_state)
        torch.random.set_rng_state(torch_state)
        random.setstate(python_state)


def _model(model: object) -> Any:
    lightning = require("grelu.lightning", _EXTRAS)
    if not isinstance(model, lightning.LightningModel):
        raise ConfigurationError(
            "model must be a grelu LightningModel; load one with "
            "LightningModel.load_from_checkpoint."
        )
    return model


def _task_transform(model: Any, task_index: object) -> Any:
    """``Aggregate`` transform that keeps one task, as the attribution helpers require."""

    n_tasks = int(model.model_params["n_tasks"])
    index = _integer(task_index, "task_index", minimum=0, maximum=max(n_tasks - 1, 0))
    transforms = require("grelu.transforms.prediction_transforms", _EXTRAS)
    return transforms.Aggregate(tasks=[index])


def _seed(value: object, name: str = "seed") -> int:
    return _integer(value, name, minimum=0, maximum=_MAX_SEEDS)


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _coverage_matrix(values: object, columns: int, name: str) -> Any:
    """Finite non-negative float matrix (tasks x intervals) from nested lists."""

    np = require("numpy", _EXTRAS)
    rows = _rows(values, name, limit=10_000)
    matrix = []
    for index, row in enumerate(rows):
        if isinstance(row, (str, bytes)) or not isinstance(row, Sequence):
            raise ConfigurationError(f"{name}[{index}] must be a list of numbers.")
        if len(row) != columns:
            raise ConfigurationError(
                f"{name}[{index}] must have one value per interval ({columns})."
            )
        for item in row:
            if isinstance(item, bool) or not isinstance(item, (int, float)):
                raise ConfigurationError(f"{name}[{index}] must contain numbers only.")
            if item != item or item in (float("inf"), float("-inf")) or item < 0:
                raise ConfigurationError(f"{name}[{index}] must contain finite values >= 0.")
        matrix.append([float(item) for item in row])
    return np.asarray(matrix, dtype=float)


def _task_names(values: object, count: int) -> list[str]:
    if values is None:
        return [f"task{index}" for index in range(count)]
    names = _rows(values, "tasks", limit=10_000)
    if len(names) != count:
        raise ConfigurationError(f"tasks must list {count} names, one per coverage row.")
    cleaned = [
        _text(item, f"tasks[{index}]", maximum=_MAX_NAME) for index, item in enumerate(names)
    ]
    if len(set(cleaned)) != len(cleaned):
        raise ConfigurationError("tasks must be unique.")
    return cleaned


def _attribution_matrix(value: object, name: str, length: int) -> Any:
    np = require("numpy", _EXTRAS)
    if hasattr(value, "tolist") and not isinstance(value, (str, bytes)):
        value = value.tolist()
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != 4:
        raise ConfigurationError(f"{name} must be 4 rows (A, C, G, T) of {length} numbers.")
    for row in value:
        if isinstance(row, (str, bytes)) or not isinstance(row, Sequence) or len(row) != length:
            raise ConfigurationError(f"{name} must be 4 rows (A, C, G, T) of {length} numbers.")
        for item in row:
            if isinstance(item, bool) or not isinstance(item, (int, float)):
                raise ConfigurationError(f"{name} must contain numbers only.")
            if item != item or item in (float("inf"), float("-inf")):
                raise ConfigurationError(f"{name} must contain finite numbers.")
    return np.asarray(value, dtype=float)


def _attribution_stack(values: object, ids: list[str], seqs: list[str]) -> Any:
    """Attributions aligned to the sequences: a list in order, or a mapping with the same keys."""

    np = require("numpy", _EXTRAS)
    if len({len(seq) for seq in seqs}) != 1:
        raise ConfigurationError("Attribution scoring needs sequences of equal length.")
    length = len(seqs[0])
    if isinstance(values, Mapping):
        if set(values) != set(ids):
            raise ConfigurationError("attributions must have the same keys as sequences.")
        items = [values[key] for key in ids]
    elif isinstance(values, Sequence) and not isinstance(values, (str, bytes)):
        if len(values) != len(ids):
            raise ConfigurationError("attributions must list one matrix per sequence.")
        items = list(values)
    else:
        raise ConfigurationError("attributions must be a list or a mapping keyed like sequences.")
    matrices = [
        _attribution_matrix(item, f"attributions[{index}]", length)
        for index, item in enumerate(items)
    ]
    return np.stack(matrices)


def _hit_rows(hits: object, ids: list[str], seqs: list[str]) -> list[dict[str, Any]]:
    """Validated motif hits (as returned by scan_motifs_with_pwms), sorted for determinism."""

    if isinstance(hits, (str, bytes, Mapping)) or not isinstance(hits, Sequence):
        raise ConfigurationError("hits must be a list of motif hits.")
    if len(hits) > _MAX_HITS:
        raise ConfigurationError(f"hits may contain at most {_MAX_HITS} entries.")
    index_of = {identifier: index for index, identifier in enumerate(ids)}
    lengths = {identifier: len(seq) for identifier, seq in zip(ids, seqs, strict=True)}
    rows = []
    for position, hit in enumerate(hits):
        label = f"hits[{position}]"
        if not isinstance(hit, Mapping):
            raise ConfigurationError(f"{label} must be a mapping.")
        motif = _text(hit.get("motif"), f"{label}.motif", maximum=_MAX_NAME)
        sequence = _text(hit.get("sequence"), f"{label}.sequence", maximum=_MAX_NAME)
        if sequence not in index_of:
            raise ConfigurationError(f"{label}.sequence {sequence!r} is not an input sequence.")
        start = _integer(hit.get("start"), f"{label}.start", minimum=0, maximum=_MAX_COORD)
        end = _integer(hit.get("end"), f"{label}.end", minimum=1, maximum=_MAX_COORD)
        if end <= start or end > lengths[sequence]:
            raise ConfigurationError(f"{label} must be an interval inside its sequence.")
        strand = _choice(hit.get("strand"), f"{label}.strand", ("+", "-"))
        rows.append(
            {
                "motif": motif,
                "sequence": sequence,
                "seq_idx": index_of[sequence],
                "start": start,
                "end": end,
                "strand": strand,
            }
        )
    rows.sort(key=lambda row: (row["seq_idx"], row["start"], row["strand"], row["motif"]))
    return rows


def _model_cells(array: Any) -> Any:
    return _clean(array.tolist())


# ---------------------------------------------------------------------------
# Intervals, chromosomes and coverage tables (no model)
# ---------------------------------------------------------------------------


def read_narrowpeak_intervals(
    path: str | Path,
    *,
    skiprows: int = 0,
    has_header: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a narrowPeak file (BED6+4) with gReLU's reader, keeping the file order.

    Each peak keeps its chrom, start, end, name, score, strand, signal, p-value, q-value and
    summit offset. Files with fewer than ten columns are rejected. gReLU's reader treats the
    first line as a header unless told otherwise, which silently drops the first peak of a
    headerless file, so ``has_header`` defaults to False and is passed to pandas as ``header``.
    """

    source = input_path(path)
    skip = _integer(skiprows, "skiprows", minimum=0, maximum=1_000_000)
    header = _flag(has_header, "has_header")
    bed = require("grelu.io.bed", _EXTRAS)
    try:
        table, _ = _quiet(
            bed.read_narrowpeak,
            str(source),
            skiprows=skip,
            header=0 if header else None,
            verbose=verbose,
        )
    except _UPSTREAM_ERRORS as exc:
        raise _rejected(_GRELU, exc) from exc
    peaks = [
        {
            "chrom": str(row["chrom"]),
            "start": int(row["start"]),
            "end": int(row["end"]),
            "name": _clean(row["name"]),
            "score": _clean(row["score"]),
            "strand": _clean(row["strand"]),
            "signal": _clean(row["signal"]),
            "pvalue": _clean(row["pvalue"]),
            "qvalue": _clean(row["qvalue"]),
            "summit": int(row["summit"]),
        }
        for row in table.to_dict(orient="records")
    ]
    return {
        "status": "completed",
        "operation": "read_narrowpeak_intervals",
        "source": str(source),
        "count": len(peaks),
        "peaks": peaks,
        "provenance": _provenance(_GRELU, format="narrowPeak BED6+4"),
    }


def expand_chromosome_names(
    chroms: str | Sequence[str],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Expand a chromosome shortcut into GRCh38-style names, or check an explicit list.

    The shortcuts are ``autosomes`` (chr1 to chr22), ``autosomesX`` (adds chrX) and
    ``autosomesXY`` (adds chrX and chrY). An explicit list is returned in the given order.
    """

    data_utils = require("grelu.data.utils", _EXTRAS)
    if isinstance(chroms, str):
        shortcut = _choice(chroms, "chroms", _CHROM_SHORTCUTS)
        names = [str(item) for item in data_utils.get_chromosomes(shortcut)]
        source = "shortcut"
    else:
        values = _rows(chroms, "chroms", limit=10_000)
        names = [_text(item, f"chroms[{index}]", maximum=64) for index, item in enumerate(values)]
        shortcut = None
        source = "list"
    return {
        "status": "completed",
        "operation": "expand_chromosome_names",
        "source": source,
        "shortcut": shortcut,
        "count": len(names),
        "chromosomes": names,
        "provenance": _provenance(_GRELU, naming="chr-prefixed names"),
    }


def check_genomic_interval_table(
    intervals: Sequence[Mapping[str, Any]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check whether a table of mappings has gReLU's genomic-interval layout.

    The first three columns must be ``chrom``, ``start`` and ``end``, with a string ``chrom``
    column and integer ``start`` and ``end`` columns. The rows are built as given, so a
    column in the wrong order or a float coordinate makes ``valid`` False.
    """

    rows = _rows(intervals, "intervals", limit=_MAX_ROWS)
    records = []
    for index, row in enumerate(rows):
        if not isinstance(row, Mapping):
            raise ConfigurationError(f"intervals[{index}] must be a mapping.")
        records.append(dict(row))
    pandas = require("pandas", _EXTRAS)
    frame = pandas.DataFrame.from_records(records)
    sequence_format = require("grelu.sequence.format", _EXTRAS)
    try:
        valid = bool(sequence_format.check_intervals(frame))
    except _UPSTREAM_ERRORS as exc:
        raise _rejected(_GRELU, exc) from exc
    return {
        "status": "completed",
        "operation": "check_genomic_interval_table",
        "row_count": len(frame),
        "columns": [str(column) for column in frame.columns],
        "dtypes": {str(column): str(frame[column].dtype) for column in frame.columns},
        "valid": valid,
        "provenance": _provenance(_GRELU, check="grelu.sequence.format.check_intervals"),
    }


def filter_interval_coverage(
    coverage: Sequence[Sequence[float]],
    intervals: Sequence[Mapping[str, Any]],
    *,
    tasks: Sequence[str] | None = None,
    aggfunc: str = "mean",
    method: str = "cutoff",
    cutoff: float = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep genomic intervals whose coverage across tasks passes a cutoff or top-n rule.

    ``coverage`` has one row per task and one column per interval. The per-interval value is
    the aggregate over tasks (``mean``, ``sum``, ``max`` or ``min``). ``method="cutoff"``
    keeps intervals at or above ``cutoff``; ``method="top"`` keeps the ``cutoff`` intervals
    with the highest value. gReLU's percentile mode and its random negative sampling are not
    exposed: the first is broken upstream and the second needs a seed.
    """

    frame = _interval_frame(intervals, "intervals")
    count = len(frame)
    matrix = _coverage_matrix(coverage, count, "coverage")
    names = _task_names(tasks, matrix.shape[0])
    aggregate = _choice(aggfunc, "aggfunc", _TASK_AGG)
    mode = _choice(method, "method", ("cutoff", "top"))
    if mode == "cutoff":
        threshold: float = _real(cutoff, "cutoff", low=0.0, high=1e12)
    else:
        threshold = _integer(cutoff, "cutoff", minimum=1, maximum=count)
    anndata = require("anndata", _EXTRAS)
    pandas = require("pandas", _EXTRAS)
    preprocess = require("grelu.data.preprocess", _EXTRAS)
    var = frame.copy()
    var.index = [f"interval{index}" for index in range(count)]
    adata = anndata.AnnData(
        X=matrix,
        obs=pandas.DataFrame(index=names),
        var=var,
    )
    try:
        kept, messages = _quiet(
            preprocess.filter_coverage,
            adata,
            aggfunc=aggregate,
            method=mode,
            cutoff=threshold,
            negative_frac=0.0,
            inplace=False,
            verbose=verbose,
        )
    except _UPSTREAM_ERRORS as exc:
        raise _rejected(_GRELU, exc) from exc
    kept_indices = sorted(int(str(name)[len("interval") :]) for name in kept.var_names)
    records = [
        {
            "index": index,
            "chrom": str(frame.loc[index, "chrom"]),
            "start": int(frame.loc[index, "start"]),
            "end": int(frame.loc[index, "end"]),
        }
        for index in kept_indices
    ]
    return {
        "status": "completed",
        "operation": "filter_interval_coverage",
        "tasks": names,
        "method": mode,
        "cutoff": threshold,
        "aggfunc": aggregate,
        "input_count": count,
        "count": len(records),
        "intervals": records,
        "messages": messages,
        "provenance": _provenance(_GRELU, filter="coverage across tasks; negatives not sampled"),
    }


def add_negative_intervals(
    coverage: Sequence[Sequence[float]],
    intervals: Sequence[Mapping[str, Any]],
    negatives: Sequence[Mapping[str, Any]],
    *,
    tasks: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Append negative control intervals after the positive ones, with zero coverage.

    Positive intervals are labelled ``positive`` and negatives ``negative`` in the output
    ``class`` field. Negative columns of the coverage matrix are all zero, as gReLU writes
    them when ``negative_labels`` is 0.
    """

    positives = _interval_frame(intervals, "intervals")
    negative_frame = _interval_frame(negatives, "negatives")
    count = len(positives)
    matrix = _coverage_matrix(coverage, count, "coverage")
    names = _task_names(tasks, matrix.shape[0])
    anndata = require("anndata", _EXTRAS)
    pandas = require("pandas", _EXTRAS)
    preprocess = require("grelu.data.preprocess", _EXTRAS)
    positive_var = positives.copy()
    positive_var.index = [f"interval{index}" for index in range(count)]
    negative_var = negative_frame.copy()
    negative_var.index = [f"negative{index}" for index in range(len(negative_frame))]
    adata = anndata.AnnData(
        X=matrix,
        obs=pandas.DataFrame(index=names),
        var=positive_var,
    )
    try:
        combined, messages = _quiet(
            preprocess.add_negatives,
            adata,
            negative_var,
            negative_labels=0,
            inplace=False,
            verbose=verbose,
        )
    except _UPSTREAM_ERRORS as exc:
        raise _rejected(_GRELU, exc) from exc
    records = [
        {
            "name": str(name),
            "chrom": str(row["chrom"]),
            "start": int(row["start"]),
            "end": int(row["end"]),
            "class": str(row["class"]),
        }
        for name, row in zip(
            combined.var_names, combined.var.to_dict(orient="records"), strict=True
        )
    ]
    return {
        "status": "completed",
        "operation": "add_negative_intervals",
        "tasks": names,
        "positive_count": count,
        "negative_count": len(negative_frame),
        "count": len(records),
        "intervals": records,
        "coverage": _model_cells(combined.X),
        "messages": messages,
        "provenance": _provenance(_GRELU, labels="negative_labels=0"),
    }


def score_motif_sites(
    hits: Sequence[Mapping[str, Any]],
    sequences: Mapping[str, str] | Sequence[str],
    attributions: Mapping[str, Any] | Sequence[Any],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Average each motif hit's attribution over the matched bases (site attribution score).

    ``hits`` are the ``hits`` of ``scan_motifs_with_pwms``. ``attributions`` holds one 4 x L
    matrix per sequence, aligned to ``sequences`` by list order or by key. The score is the
    mean attribution over the hit interval, so it reflects importance, not PWM shape.
    """

    ids, seqs = _sequence_items(sequences)
    stack = _attribution_stack(attributions, ids, seqs)
    rows = _hit_rows(hits, ids, seqs)
    if not rows:
        return {
            "status": "completed",
            "operation": "score_motif_sites",
            "count": 0,
            "hits": [],
            "provenance": _provenance(_GRELU, score="mean attribution over the site"),
        }
    pandas = require("pandas", _EXTRAS)
    motif_module = require("grelu.interpret.motifs", _EXTRAS)
    frame = pandas.DataFrame(rows)
    try:
        scored, messages = _quiet(motif_module.score_sites, frame, stack, seqs, verbose=verbose)
    except _UPSTREAM_ERRORS as exc:
        raise _rejected(_GRELU, exc) from exc
    output = []
    for row, value in zip(rows, scored["site_attr_score"].tolist(), strict=True):
        output.append(
            {
                "motif": row["motif"],
                "sequence": row["sequence"],
                "start": row["start"],
                "end": row["end"],
                "strand": row["strand"],
                "site_attr_score": _clean(float(value)),
            }
        )
    return {
        "status": "completed",
        "operation": "score_motif_sites",
        "count": len(output),
        "hits": output,
        "messages": messages,
        "provenance": _provenance(_GRELU, score="mean attribution over the site"),
    }


def score_motif_pwm_products(
    hits: Sequence[Mapping[str, Any]],
    sequences: Mapping[str, str] | Sequence[str],
    attributions: Mapping[str, Any] | Sequence[Any],
    motifs: Mapping[str, Sequence[Sequence[float]]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Score each hit by the sum of the element-wise product of its PWM and the attributions.

    The PWM of each motif is a 4 x L matrix (rows A, C, G, T) whose L must equal the hit
    length. Minus-strand hits use the reverse-complement PWM. The score ranks motifs that
    match the same region by how well their shape fits the attribution profile.
    """

    ids, seqs = _sequence_items(sequences)
    stack = _attribution_stack(attributions, ids, seqs)
    pwms = _motif_set(motifs, "motifs")
    rows = _hit_rows(hits, ids, seqs)
    for row in rows:
        if row["motif"] not in pwms:
            raise ConfigurationError(f"hit motif {row['motif']!r} is not in motifs.")
        if pwms[row["motif"]].shape[1] != row["end"] - row["start"]:
            raise ConfigurationError(
                f"hit motif {row['motif']!r} must span the PWM length "
                f"{pwms[row['motif']].shape[1]}."
            )
    if not rows:
        return {
            "status": "completed",
            "operation": "score_motif_pwm_products",
            "count": 0,
            "hits": [],
            "provenance": _provenance(_GRELU, score="sum of PWM times attribution"),
        }
    pandas = require("pandas", _EXTRAS)
    motif_module = require("grelu.interpret.motifs", _EXTRAS)
    frame = pandas.DataFrame(rows)
    try:
        scored, messages = _quiet(motif_module.score_motifs, frame, stack, pwms, verbose=verbose)
    except _UPSTREAM_ERRORS as exc:
        raise _rejected(_GRELU, exc) from exc
    output = []
    for row, value in zip(rows, scored["motif_attr_score"].tolist(), strict=True):
        output.append(
            {
                "motif": row["motif"],
                "sequence": row["sequence"],
                "start": row["start"],
                "end": row["end"],
                "strand": row["strand"],
                "motif_attr_score": _clean(float(value)),
            }
        )
    return {
        "status": "completed",
        "operation": "score_motif_pwm_products",
        "count": len(output),
        "hits": output,
        "messages": messages,
        "provenance": _provenance(_GRELU, score="sum of PWM times attribution"),
    }


def read_tomtom_matches(
    tomtom_dir: str | Path,
    *,
    qthresh: float = 0.05,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read every ``tomtom.tsv`` file in a TOMTOM output folder, keeping q-values up to ``qthresh``.

    Rows with any missing value are dropped, as gReLU does. Matches are sorted by query,
    target, orientation and offset. A folder without a ``tomtom.tsv`` file is rejected.
    """

    folder = Path(_text(str(tomtom_dir), "tomtom_dir", maximum=4096)).expanduser().resolve()
    if not folder.is_dir():
        raise ConfigurationError(f"tomtom_dir is not a directory: {folder}")
    if not any(name.endswith("tomtom.tsv") for name in _listdir(folder)):
        raise ConfigurationError(f"no tomtom.tsv file found in {folder}")
    threshold = _real(qthresh, "qthresh", low=0.0, high=1.0, low_open=True)
    grelu_io = require("grelu.io", _EXTRAS)
    try:
        table, _ = _quiet(grelu_io.read_tomtom, str(folder), qthresh=threshold, verbose=verbose)
    except _UPSTREAM_ERRORS as exc:
        raise _rejected(_GRELU, exc) from exc
    matches = [
        {str(key): _clean(value) for key, value in row.items()}
        for row in table.to_dict(orient="records")
    ]
    matches.sort(
        key=lambda item: (
            str(item.get("Query_ID")),
            str(item.get("Target_ID")),
            str(item.get("Orientation")),
            float(item.get("Optimal_offset") or 0),
        )
    )
    return {
        "status": "completed",
        "operation": "read_tomtom_matches",
        "source": str(folder),
        "qthresh": threshold,
        "count": len(matches),
        "matches": matches,
        "provenance": _provenance(_GRELU, filter="q-value <= qthresh"),
    }


def _listdir(folder: Path) -> list[str]:
    import os

    return sorted(os.listdir(folder))


def read_modisco_motifs(
    h5_path: str | Path,
    *,
    group: str | None = None,
    names: Sequence[str] | None = None,
    trim_threshold: float = 0.3,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read TF-MoDISco motifs from an H5 report as trimmed position probability matrices.

    ``group`` is ``pos``, ``neg`` or None for both. ``names`` picks motifs by name and
    overrides ``group``; names must start with ``pos`` or ``neg``. Each PPM is trimmed from
    both ends until a base reaches ``trim_threshold``.
    """

    path = input_path(h5_path)
    chosen = None if group is None else _choice(group, "group", ("pos", "neg"))
    selected: list[str] | None = None
    if names is not None:
        values = _rows(names, "names", limit=_MAX_ROWS)
        selected = [
            _text(item, f"names[{index}]", maximum=_MAX_NAME) for index, item in enumerate(values)
        ]
        for name in selected:
            if not name.startswith(("pos", "neg")):
                raise ConfigurationError("every name must start with pos or neg.")
    threshold = _real(trim_threshold, "trim_threshold", low=0.0, high=1.0)
    io_motifs = require("grelu.io.motifs", _EXTRAS)
    try:
        motifs, _ = _quiet(
            io_motifs.read_modisco_report,
            str(path),
            group=chosen,
            names=selected,
            trim_threshold=threshold,
            verbose=verbose,
        )
    except _UPSTREAM_ERRORS as exc:
        raise _rejected(_GRELU, exc) from exc
    records = [
        {"name": name, "length": int(motifs[name].shape[1]), "ppm": _model_cells(motifs[name])}
        for name in sorted(motifs)
    ]
    return {
        "status": "completed",
        "operation": "read_modisco_motifs",
        "source": str(path),
        "group": chosen or "all",
        "trim_threshold": threshold,
        "count": len(records),
        "motifs": records,
        "provenance": _provenance(_GRELU, trimming="CWM information threshold"),
    }


# ---------------------------------------------------------------------------
# Model-based experiments (LightningModel objects, offline on CPU)
# ---------------------------------------------------------------------------


def get_sequence_attributions(
    model: Any,
    sequences: str | Sequence[str],
    *,
    task_index: int = 0,
    method: str = "deepshap",
    n_shuffles: int = 20,
    seed: int = 0,
    correct_grad: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Per-base attribution scores (B, 4, L) of one task for equal-length DNA sequences.

    ``method`` is ``deepshap``, ``saliency``, ``inputxgradient`` or ``integratedgradients``.
    ``n_shuffles`` and ``seed`` apply to DeepSHAP's dinucleotide shuffles. ``correct_grad``
    is only supported with ``saliency``. Global NumPy and torch random state is restored.
    """

    checked = _model(model)
    seqs = _dna_batch(sequences, "sequences")
    if len({len(seq) for seq in seqs}) != 1:
        raise ConfigurationError("sequences must all have the same length.")
    transform = _task_transform(checked, task_index)
    mode = _choice(method, "method", _ATTRIBUTION_METHODS)
    shuffles = _integer(n_shuffles, "n_shuffles", minimum=1, maximum=1000)
    seed_value = _seed(seed)
    flag = _flag(correct_grad, "correct_grad")
    if flag and mode != "saliency":
        raise ConfigurationError("correct_grad is only supported with method='saliency'.")
    score = require("grelu.interpret.score", _EXTRAS)
    with _rng_preserved():
        try:
            values, messages = _quiet(
                score.get_attributions,
                checked,
                seqs,
                prediction_transform=transform,
                method=mode,
                correct_grad=flag,
                n_shuffles=shuffles,
                seed=seed_value,
                verbose=verbose,
            )
        except _UPSTREAM_ERRORS as exc:
            raise _rejected(_GRELU, exc) from exc
    np = require("numpy", _EXTRAS)
    array = np.asarray(values)
    return {
        "status": "completed",
        "operation": "get_sequence_attributions",
        "method": mode,
        "task_index": int(task_index),
        "count": int(array.shape[0]),
        "length": int(array.shape[-1]),
        "attributions": _model_cells(array),
        "messages": messages,
        "provenance": _provenance(
            _GRELU, method=mode, shuffles=str(shuffles), seed=str(seed_value)
        ),
    }


def in_silico_mutagenesis_scores(
    model: Any,
    sequence: str,
    *,
    task_index: int = 0,
    compare_func: str = "subtract",
    start: int = 0,
    end: int | None = None,
    batch_size: int = 64,
    verbose: bool = True,
) -> dict[str, Any]:
    """In silico mutagenesis: effect of each alternative base at each position of one task.

    ``compare_func`` is ``subtract``, ``divide`` or ``log2FC`` (mutant against reference). The
    output has one row per position from ``start`` to ``end`` (exclusive) and one column per
    base A, C, G, T; the reference base scores zero under ``subtract``.
    """

    checked = _model(model)
    seq = _dna(sequence, "sequence", maximum=_MAX_SEQ)
    transform = _task_transform(checked, task_index)
    start_pos = _integer(start, "start", minimum=0, maximum=len(seq) - 1)
    end_pos = None if end is None else _integer(end, "end", minimum=start_pos + 1, maximum=len(seq))
    compare = _choice(compare_func, "compare_func", ("subtract", "divide", "log2FC"))
    size = _integer(batch_size, "batch_size", minimum=1, maximum=4096)
    score = require("grelu.interpret.score", _EXTRAS)
    try:
        values, messages = _quiet(
            score.ISM_predict,
            seq,
            checked,
            prediction_transform=transform,
            start_pos=start_pos,
            end_pos=end_pos,
            compare_func=compare,
            batch_size=size,
            return_df=False,
            verbose=verbose,
        )
    except _UPSTREAM_ERRORS as exc:
        raise _rejected(_GRELU, exc) from exc
    np = require("numpy", _EXTRAS)
    array = np.asarray(values).reshape(-1, 4)
    return {
        "status": "completed",
        "operation": "in_silico_mutagenesis_scores",
        "task_index": int(task_index),
        "compare_func": compare,
        "start": start_pos,
        "positions": int(array.shape[0]),
        "bases": ["A", "C", "G", "T"],
        "scores": _model_cells(array),
        "messages": messages,
        "provenance": _provenance(_GRELU, method="ISM_predict"),
    }


def get_attention_scores(
    model: Any,
    sequence: str,
    *,
    block_idx: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Attention maps of a transformer model for one sequence.

    Returns (layers, heads, L, L) when ``block_idx`` is None, or (heads, L, L) for one block.
    The model needs ``n_transformers`` of at least 1 in its parameters.
    """

    checked = _model(model)
    seq = _dna(sequence, "sequence", maximum=_MAX_SEQ)
    blocks = int(checked.model_params.get("n_transformers") or 0)
    if blocks < 1:
        raise ConfigurationError("model has no transformer blocks (n_transformers < 1).")
    index = (
        None
        if block_idx is None
        else _integer(block_idx, "block_idx", minimum=0, maximum=blocks - 1)
    )
    score = require("grelu.interpret.score", _EXTRAS)
    try:
        attention, messages = _quiet(
            score.get_attention_scores, checked, seq, block_idx=index, verbose=verbose
        )
    except _UPSTREAM_ERRORS as exc:
        raise _rejected(_GRELU, exc) from exc
    np = require("numpy", _EXTRAS)
    array = np.asarray(attention)
    return {
        "status": "completed",
        "operation": "get_attention_scores",
        "block_idx": index,
        "shape": [int(item) for item in array.shape],
        "attention": _model_cells(array),
        "messages": messages,
        "provenance": _provenance(_GRELU, source="grelu.interpret.score.get_attention_scores"),
    }


def marginalize_pattern_insertions(
    model: Any,
    sequences: str | Sequence[str],
    patterns: str | Sequence[str],
    *,
    task_index: int = 0,
    compare_func: str = "subtract",
    n_shuffles: int = 1,
    seed: int = 0,
    rc: bool = False,
    batch_size: int = 64,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predictions for background sequences before and after inserting each pattern.

    ``compare_func`` is ``subtract`` (after minus before) or ``divide`` (after over before).
    The raw before and after predictions are not offered, and ``n_shuffles`` must be at least
    1: gReLU's raw mode and its zero-shuffle mode both fail with a concat error on the installed
    torch. Shuffling needs backgrounds that are not dinucleotide-periodic.
    """

    checked = _model(model)
    background = _dna_batch(sequences, "sequences")
    inserts = _dna_batch(patterns, "patterns")
    if max(len(item) for item in inserts) >= min(len(item) for item in background):
        raise ConfigurationError("every pattern must be shorter than every background sequence.")
    transform = _task_transform(checked, task_index)
    shuffles = _integer(n_shuffles, "n_shuffles", minimum=1, maximum=1000)
    seed_value = _seed(seed)
    reverse = _flag(rc, "rc")
    compare = _choice(compare_func, "compare_func", ("subtract", "divide"))
    size = _integer(batch_size, "batch_size", minimum=1, maximum=4096)
    simulate = require("grelu.interpret.simulate", _EXTRAS)
    with _rng_preserved():
        try:
            result, messages = _quiet(
                simulate.marginalize_patterns,
                checked,
                inserts,
                background,
                batch_size=size,
                n_shuffles=shuffles,
                seed=seed_value,
                prediction_transform=transform,
                rc=reverse,
                compare_func=compare,
                verbose=verbose,
            )
        except _UPSTREAM_ERRORS as exc:
            raise _rejected(_GRELU, exc) from exc
    base = {
        "status": "completed",
        "operation": "marginalize_pattern_insertions",
        "task_index": int(task_index),
        "patterns": inserts,
        "compare_func": compare,
        "messages": messages,
        "provenance": _provenance(_GRELU, method="marginalize_patterns"),
    }
    base["effects"] = _model_cells(result)
    base["shape"] = [int(item) for item in result.shape]
    return base


def marginalize_pattern_spacing(
    model: Any,
    sequence: str,
    fixed_pattern: str,
    moving_pattern: str,
    *,
    task_index: int = 0,
    stride: int = 1,
    n_shuffles: int = 1,
    seed: int = 0,
    rc: bool = False,
    compare_func: str | None = None,
    batch_size: int = 64,
    verbose: bool = True,
) -> dict[str, Any]:
    """Effect of the spacing between two patterns on a shuffled background sequence.

    The fixed pattern is inserted at the centre and the moving pattern at every ``stride``
    spacing. Spacing positions are returned with the predictions, before and after insertion
    when ``compare_func`` is None, otherwise as the comparison only.
    """

    checked = _model(model)
    seq = _dna(sequence, "sequence", maximum=_MAX_SEQ)
    fixed = _dna(fixed_pattern, "fixed_pattern", maximum=_MAX_SEQ)
    moving = _dna(moving_pattern, "moving_pattern", maximum=_MAX_SEQ)
    transform = _task_transform(checked, task_index)
    step = _integer(stride, "stride", minimum=1, maximum=_MAX_SEQ)
    shuffles = _integer(n_shuffles, "n_shuffles", minimum=1, maximum=1000)
    seed_value = _seed(seed)
    reverse = _flag(rc, "rc")
    compare = (
        None
        if compare_func is None
        else _choice(compare_func, "compare_func", ("subtract", "divide"))
    )
    size = _integer(batch_size, "batch_size", minimum=1, maximum=4096)
    simulate = require("grelu.interpret.simulate", _EXTRAS)
    with _rng_preserved():
        try:
            result, messages = _quiet(
                simulate.marginalize_pattern_spacing,
                checked,
                seq,
                fixed,
                moving,
                stride=step,
                n_shuffles=shuffles,
                seed=seed_value,
                rc=reverse,
                compare_func=compare,
                batch_size=size,
                prediction_transform=transform,
                verbose=verbose,
            )
        except _UPSTREAM_ERRORS as exc:
            raise _rejected(_GRELU, exc) from exc
    if compare is None:
        before, after, spacings = result
        body: dict[str, Any] = {
            "before": _model_cells(before),
            "after": _model_cells(after),
            "shape": [int(item) for item in after.shape],
        }
    else:
        effects, spacings = result
        body = {"effects": _model_cells(effects), "shape": [int(item) for item in effects.shape]}
    return {
        "status": "completed",
        "operation": "marginalize_pattern_spacing",
        "task_index": int(task_index),
        "compare_func": compare,
        "spacings": [_clean(item) for item in spacings],
        **body,
        "messages": messages,
        "provenance": _provenance(_GRELU, method="marginalize_pattern_spacing"),
    }


def shuffle_tile_effects(
    model: Any,
    sequences: str | Sequence[str],
    tile_len: int,
    *,
    task_index: int = 0,
    stride: int | None = None,
    protect_center: int | None = None,
    n_shuffles: int = 1,
    seed: int = 0,
    compare_func: str | None = None,
    batch_size: int = 64,
    verbose: bool = True,
) -> dict[str, Any]:
    """Shuffle each tile of equal-length sequences and measure the prediction change.

    Tiles start every ``stride`` bases (default ``tile_len``). ``protect_center`` leaves a
    central window unshuffled. Returns the tile table (start and end) with the predictions
    before and after shuffling, or the comparison when ``compare_func`` is given.
    """

    checked = _model(model)
    seqs = _dna_batch(sequences, "sequences")
    if len({len(seq) for seq in seqs}) != 1:
        raise ConfigurationError("sequences must all have the same length.")
    length = len(seqs[0])
    tile = _integer(tile_len, "tile_len", minimum=1, maximum=length)
    step = None if stride is None else _integer(stride, "stride", minimum=1, maximum=length)
    protect = (
        None
        if protect_center is None
        else _integer(protect_center, "protect_center", minimum=1, maximum=length)
    )
    transform = _task_transform(checked, task_index)
    shuffles = _integer(n_shuffles, "n_shuffles", minimum=1, maximum=1000)
    seed_value = _seed(seed)
    compare = (
        None
        if compare_func is None
        else _choice(compare_func, "compare_func", ("subtract", "divide"))
    )
    size = _integer(batch_size, "batch_size", minimum=1, maximum=4096)
    simulate = require("grelu.interpret.simulate", _EXTRAS)
    with _rng_preserved():
        try:
            result, messages = _quiet(
                simulate.shuffle_tiles,
                checked,
                seqs,
                tile_len=tile,
                stride=step,
                protect_center=protect,
                n_shuffles=shuffles,
                seed=seed_value,
                compare_func=compare,
                batch_size=size,
                prediction_transform=transform,
                verbose=verbose,
            )
        except _UPSTREAM_ERRORS as exc:
            raise _rejected(_GRELU, exc) from exc
    if compare is None:
        before, after, tiles = result
        body: dict[str, Any] = {
            "before": _model_cells(before),
            "after": _model_cells(after),
            "shape": [int(item) for item in after.shape],
        }
    else:
        effects, tiles = result
        body = {"effects": _model_cells(effects), "shape": [int(item) for item in effects.shape]}
    tile_rows = [
        {"start": int(row["start"]), "end": int(row["end"])}
        for row in tiles.to_dict(orient="records")
    ]
    return {
        "status": "completed",
        "operation": "shuffle_tile_effects",
        "task_index": int(task_index),
        "compare_func": compare,
        "tile_count": len(tile_rows),
        "tiles": tile_rows,
        **body,
        "messages": messages,
        "provenance": _provenance(_GRELU, method="shuffle_tiles"),
    }


def _check_reference_alleles(records: list[dict[str, Any]], fasta: Path) -> None:
    """Confirm each reference allele matches the local genome (1-based ``pos``).

    gReLU's own check is shadowed upstream by a flag of the same name, so it cannot run.
    """

    genome_module = require("grelu.io.genome", _EXTRAS)
    try:
        genome = genome_module.get_genome(str(fasta))
    except _UPSTREAM_ERRORS as exc:
        raise _rejected(_GRELU, exc) from exc
    for index, row in enumerate(records):
        try:
            # pyfaidx coordinates are 1-based and inclusive.
            observed = str(
                genome.get_seq(row["chrom"], row["pos"], row["pos"] + len(row["ref"]) - 1)
            )
        except _UPSTREAM_ERRORS as exc:
            raise ConfigurationError(
                f"variants[{index}] cannot be read from the genome: {exc}"
            ) from exc
        if observed.upper() != row["ref"]:
            raise ConfigurationError(
                f"variants[{index}] reference {row['ref']} does not match the genome ({observed})."
            )


def grelu_predict_variant_effects(
    model: Any,
    variants: Sequence[Mapping[str, Any]],
    genome_fasta: str | Path,
    *,
    task_index: int = 0,
    seq_len: int | None = None,
    compare_func: str = "divide",
    rc: bool = False,
    max_seq_shift: int = 0,
    check_reference: bool = True,
    batch_size: int = 64,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predicted effect of each SNV or indel on one task, from a local reference FASTA.

    Each variant is a mapping with ``chrom``, 1-based ``pos``, ``ref`` and ``alt``. The genome
    is a local FASTA with a ``.fai`` index and a ``.sizes`` file. With ``check_reference`` the
    reference allele must match the genome. ``seq_len`` defaults to the model's training
    length. ``compare_func`` is ``divide`` (alt over ref) or ``subtract``.
    """

    checked = _model(model)
    fasta = input_path(genome_fasta)
    if not Path(str(fasta) + ".fai").is_file() or not Path(str(fasta) + ".sizes").is_file():
        raise ConfigurationError("genome_fasta needs a .fai index and a .sizes file next to it.")
    rows = _rows(variants, "variants", limit=_MAX_ROWS)
    records = []
    for index, row in enumerate(rows):
        if not isinstance(row, Mapping):
            raise ConfigurationError(f"variants[{index}] must be a mapping.")
        records.append(
            {
                "chrom": _text(row.get("chrom"), f"variants[{index}].chrom", maximum=64),
                "pos": _integer(
                    row.get("pos"), f"variants[{index}].pos", minimum=1, maximum=_MAX_COORD
                ),
                "ref": _dna(row.get("ref"), f"variants[{index}].ref", maximum=1000),
                "alt": _dna(row.get("alt"), f"variants[{index}].alt", maximum=1000),
            }
        )
    if seq_len is None:
        training = checked.data_params.get("train") or {}
        if "seq_len" not in training:
            raise ConfigurationError("seq_len is required for models without training metadata.")
        window = int(training["seq_len"])
    else:
        window = _integer(seq_len, "seq_len", minimum=10, maximum=100_000)
    transform = _task_transform(checked, task_index)
    compare = _choice(compare_func, "compare_func", ("divide", "subtract"))
    shift = _integer(max_seq_shift, "max_seq_shift", minimum=0, maximum=10_000)
    size = _integer(batch_size, "batch_size", minimum=1, maximum=4096)
    reverse = _flag(rc, "rc")
    check = _flag(check_reference, "check_reference")
    pandas = require("pandas", _EXTRAS)
    np = require("numpy", _EXTRAS)
    if check:
        _check_reference_alleles(records, fasta)
    frame = pandas.DataFrame(records)
    grelu_variant = require("grelu.variant", _EXTRAS)
    try:
        values, messages = _quiet(
            grelu_variant.predict_variant_effects,
            frame,
            checked,
            seq_len=window,
            genome=str(fasta),
            rc=reverse,
            max_seq_shift=shift,
            compare_func=compare,
            batch_size=size,
            return_ad=False,
            check_reference=False,
            prediction_transform=transform,
            verbose=verbose,
        )
    except _UPSTREAM_ERRORS as exc:
        raise _rejected(_GRELU, exc) from exc
    array = np.asarray(values).reshape(len(records), -1)
    return {
        "status": "completed",
        "operation": "grelu_predict_variant_effects",
        "task_index": int(task_index),
        "compare_func": compare,
        "seq_len": window,
        "count": len(records),
        "variants": records,
        "effects": _model_cells(array),
        "messages": messages,
        "provenance": _provenance(_GRELU, genome="local FASTA", reference_check=str(check)),
    }


def grelu_marginalize_variants(
    model: Any,
    variants: Sequence[Mapping[str, Any]],
    genome_fasta: str | Path,
    *,
    task_index: int = 0,
    seq_len: int | None = None,
    n_shuffles: int = 20,
    seed: int = 0,
    compare_func: str = "log2FC",
    rc: bool = False,
    max_seq_shift: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Marginalization test of each variant against dinucleotide-shuffled backgrounds.

    Each variant's reference and alternate windows are predicted for ``n_shuffles``
    shuffles, then compared with a two-sided test. The genome and variant rules match
    ``grelu_predict_variant_effects``. The output is gReLU's result table.
    """

    checked = _model(model)
    fasta = input_path(genome_fasta)
    if not Path(str(fasta) + ".fai").is_file() or not Path(str(fasta) + ".sizes").is_file():
        raise ConfigurationError("genome_fasta needs a .fai index and a .sizes file next to it.")
    rows = _rows(variants, "variants", limit=_MAX_ROWS)
    records = []
    for index, row in enumerate(rows):
        if not isinstance(row, Mapping):
            raise ConfigurationError(f"variants[{index}] must be a mapping.")
        records.append(
            {
                "chrom": _text(row.get("chrom"), f"variants[{index}].chrom", maximum=64),
                "pos": _integer(
                    row.get("pos"), f"variants[{index}].pos", minimum=1, maximum=_MAX_COORD
                ),
                "ref": _dna(row.get("ref"), f"variants[{index}].ref", maximum=1000),
                "alt": _dna(row.get("alt"), f"variants[{index}].alt", maximum=1000),
            }
        )
    if seq_len is None:
        training = checked.data_params.get("train") or {}
        if "seq_len" not in training:
            raise ConfigurationError("seq_len is required for models without training metadata.")
        window = int(training["seq_len"])
    else:
        window = _integer(seq_len, "seq_len", minimum=10, maximum=100_000)
    transform = _task_transform(checked, task_index)
    shuffles = _integer(n_shuffles, "n_shuffles", minimum=1, maximum=1000)
    seed_value = _seed(seed)
    compare = _choice(compare_func, "compare_func", ("divide", "subtract", "log2FC"))
    shift = _integer(max_seq_shift, "max_seq_shift", minimum=0, maximum=10_000)
    reverse = _flag(rc, "rc")
    pandas = require("pandas", _EXTRAS)
    frame = pandas.DataFrame(records)
    grelu_variant = require("grelu.variant", _EXTRAS)
    with _rng_preserved():
        try:
            result, messages = _quiet(
                grelu_variant.marginalize_variants,
                checked,
                frame,
                genome=str(fasta),
                seq_len=window,
                n_shuffles=shuffles,
                seed=seed_value,
                prediction_transform=transform,
                compare_func=compare,
                rc=reverse,
                max_seq_shift=shift,
                verbose=verbose,
            )
        except _UPSTREAM_ERRORS as exc:
            raise _rejected(_GRELU, exc) from exc
    table = {str(key): _clean(_table_cells(value)) for key, value in dict(result).items()}
    return {
        "status": "completed",
        "operation": "grelu_marginalize_variants",
        "task_index": int(task_index),
        "compare_func": compare,
        "seq_len": window,
        "n_shuffles": shuffles,
        "count": len(records),
        "variants": records,
        "result": table,
        "messages": messages,
        "provenance": _provenance(_GRELU, test="two-sided marginalization test"),
    }


def _table_cells(value: Any) -> Any:
    if hasattr(value, "to_dict"):
        return value.to_dict(orient="list")
    if hasattr(value, "tolist"):
        return value.tolist()
    return value


def design_sequences_evolve(
    sequences: str | Sequence[str],
    model: Any,
    *,
    task_index: int = 0,
    method: str = "ism",
    patterns: Sequence[str] | None = None,
    max_iter: int = 10,
    positions: Sequence[int] | None = None,
    batch_size: int = 64,
    return_seqs: str = "all",
    verbose: bool = True,
) -> dict[str, Any]:
    """Greedy directed evolution of each sequence toward a higher prediction of one task.

    ``method="ism"`` tries every single-base change; ``method="pattern"`` tries inserting
    each of ``patterns``. ``return_seqs`` is ``all``, ``best`` or ``none``. The model must
    carry task names in ``data_params["tasks"]["name"]``, as gReLU sets when a model is
    trained or loaded. The output is one record per iteration and candidate.
    """

    checked = _model(model)
    seqs = _dna_batch(sequences, "sequences")
    tasks_meta = checked.data_params.get("tasks") or {}
    if "name" not in tasks_meta:
        raise ConfigurationError(
            "model.data_params['tasks']['name'] must list the task names; "
            "gReLU sets it when a model is trained or loaded."
        )
    transform = _task_transform(checked, task_index)
    mode = _choice(method, "method", ("ism", "pattern"))
    inserts = None
    if mode == "pattern":
        if not patterns:
            raise ConfigurationError("patterns are required with method='pattern'.")
        inserts = _dna_batch(list(patterns), "patterns")
    rounds = _integer(max_iter, "max_iter", minimum=1, maximum=100)
    locations = None
    if positions is not None:
        values = _rows(positions, "positions", limit=_MAX_SEQ)
        locations = [
            _integer(item, f"positions[{index}]", minimum=0, maximum=min(len(s) for s in seqs) - 1)
            for index, item in enumerate(values)
        ]
    size = _integer(batch_size, "batch_size", minimum=1, maximum=4096)
    keep = _choice(return_seqs, "return_seqs", ("all", "best", "none"))
    design = require("grelu.design", _EXTRAS)
    try:
        frame, messages = _quiet(
            design.evolve,
            seqs,
            checked,
            method=mode,
            patterns=inserts,
            prediction_transform=transform,
            max_iter=rounds,
            positions=locations,
            batch_size=size,
            return_seqs=keep,
            verbose=verbose,
        )
    except _UPSTREAM_ERRORS as exc:
        raise _rejected(_GRELU, exc) from exc
    records = [
        {str(key): _clean(value) for key, value in row.items()}
        for row in frame.to_dict(orient="records")
    ]
    return {
        "status": "completed",
        "operation": "design_sequences_evolve",
        "method": mode,
        "task_index": int(task_index),
        "max_iter": rounds,
        "count": len(records),
        "records": records,
        "messages": messages,
        "provenance": _provenance(_GRELU, method="greedy directed evolution"),
    }


def design_sequence_ledidi(
    sequence: str,
    model: Any,
    *,
    task_index: int = 0,
    max_iter: int = 20,
    batch_size: int = 16,
    positions: Sequence[int] | None = None,
    lr: float = 1.0,
    edit_penalty: float = 0.1,
    seed: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Design edits of one sequence with Ledidi, maximizing one task's prediction.

    The loss is the negative mean of the prediction. Only ``positions`` may change; all
    positions may change when it is None. The output lists the designed sequences and the
    number of bases edited in each. ``grelu.design.ledidi`` passes a 0-dimensional target
    that ledidi 2.2.0 rejects, so this wrapper calls ``ledidi.Ledidi`` directly with a
    one-element target and the same loss. Global torch random state is restored.
    """

    checked = _model(model)
    seq = _dna(sequence, "sequence", maximum=_MAX_SEQ)
    transform = _task_transform(checked, task_index)
    rounds = _integer(max_iter, "max_iter", minimum=1, maximum=100_000)
    size = _integer(batch_size, "batch_size", minimum=1, maximum=4096)
    rate = _real(lr, "lr", low=0.0, high=1e6, low_open=True)
    penalty = _real(edit_penalty, "edit_penalty", low=0.0, high=1e6)
    seed_value = _seed(seed)
    locations = None
    if positions is not None:
        values = _rows(positions, "positions", limit=len(seq))
        locations = [
            _integer(item, f"positions[{index}]", minimum=0, maximum=len(seq) - 1)
            for index, item in enumerate(values)
        ]
    ledidi = require("ledidi", _EXTRAS)
    torch = require("torch", _EXTRAS)
    grelu_format = require("grelu.sequence.format", _EXTRAS)
    x = grelu_format.convert_input_type(seq, "one_hot").unsqueeze(0)
    mask = torch.ones(x.shape[-1], dtype=torch.bool)
    if locations is not None:
        mask[locations] = False

    def loss_func(output: Any, target: Any) -> Any:
        return -output.mean()

    checked.add_transform(transform)
    checked.eval()
    with _rng_preserved():
        try:
            designs, messages = _quiet(
                _run_ledidi,
                ledidi,
                checked,
                x,
                loss_func,
                rounds,
                mask,
                size,
                rate,
                penalty,
                seed_value,
                verbose=verbose,
            )
        except _UPSTREAM_ERRORS as exc:
            raise _rejected(_GRELU, exc) from exc
        finally:
            checked.reset_transform()
    designed = grelu_format.convert_input_type(designs, "strings")
    designed = [str(item) for item in designed]
    edits = [sum(1 for a, b in zip(seq, item, strict=True) if a != b) for item in designed]
    return {
        "status": "completed",
        "operation": "design_sequence_ledidi",
        "task_index": int(task_index),
        "input": seq,
        "count": len(designed),
        "designs": designed,
        "edit_counts": edits,
        "messages": messages,
        "provenance": _provenance(_GRELU, designer="ledidi.Ledidi"),
    }


def _run_ledidi(
    ledidi: Any,
    model: Any,
    x: Any,
    loss_func: Any,
    rounds: int,
    mask: Any,
    size: int,
    rate: float,
    penalty: float,
    seed: int,
) -> Any:
    import torch

    torch.manual_seed(seed)
    designer = ledidi.Ledidi(
        model,
        x[0].shape,
        target=None,
        output_loss=loss_func,
        max_iter=rounds,
        input_mask=mask,
        batch_size=size,
        lr=rate,
        l=penalty,
        random_state=seed,
        verbose=False,
    )
    return designer.fit_transform(x, torch.zeros(1)).cpu()


def run_modisco_motifs(
    model: Any,
    sequences: str | Sequence[str],
    *,
    out_dir: str | Path,
    task_index: int = 0,
    window: int | None = None,
    method: str = "deepshap",
    n_shuffles: int = 10,
    seed: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Attribution-based motif discovery with TF-MoDISco over equal-length sequences.

    The attributions of one task are computed first, then TF-MoDISco finds seqlet patterns.
    When patterns are found, the report ``modisco_report.h5`` and its logos are written to
    ``out_dir``; its motifs are listed here. Random sequences usually give no patterns, in
    which case ``patterns_found`` is False and nothing is written. ``out_dir`` must be absent
    or empty.
    """

    checked = _model(model)
    seqs = _dna_batch(sequences, "sequences")
    if len({len(seq) for seq in seqs}) != 1:
        raise ConfigurationError("sequences must all have the same length.")
    target = Path(_text(str(out_dir), "out_dir", maximum=4096)).expanduser().resolve()
    if target.exists() and (not target.is_dir() or any(target.iterdir())):
        raise ConfigurationError("out_dir must be absent or an empty directory.")
    if not target.parent.is_dir():
        raise ConfigurationError(f"parent of out_dir does not exist: {target.parent}")
    transform = _task_transform(checked, task_index)
    span = None if window is None else _integer(window, "window", minimum=1, maximum=len(seqs[0]))
    mode = _choice(method, "method", _ATTRIBUTION_METHODS)
    shuffles = _integer(n_shuffles, "n_shuffles", minimum=1, maximum=1000)
    seed_value = _seed(seed)
    modisco = require("grelu.interpret.modisco", _EXTRAS)
    with _rng_preserved():
        try:
            _, messages = _quiet(
                modisco.run_modisco,
                checked,
                seqs,
                prediction_transform=transform,
                window=span,
                out_dir=str(target),
                method=mode,
                n_shuffles=shuffles,
                seed=seed_value,
                verbose=verbose,
            )
        except _UPSTREAM_ERRORS as exc:
            raise _rejected(_GRELU, exc) from exc
    report = target / "modisco_report.h5"
    found = report.is_file()
    motifs: list[str] = []
    if found:
        io_motifs = require("grelu.io.motifs", _EXTRAS)
        try:
            found_motifs, _ = _quiet(io_motifs.read_modisco_report, str(report), verbose=False)
        except _UPSTREAM_ERRORS as exc:
            raise _rejected(_GRELU, exc) from exc
        motifs = sorted(found_motifs)
    files = (
        sorted(str(path.relative_to(target)) for path in target.rglob("*"))
        if target.is_dir()
        else []
    )
    return {
        "status": "completed",
        "operation": "run_modisco_motifs",
        "out_dir": str(target),
        "task_index": int(task_index),
        "patterns_found": found,
        "motifs": motifs,
        "files": files,
        "messages": messages,
        "provenance": _provenance(
            _GRELU, method="TF-MoDISco via grelu.interpret.modisco.run_modisco"
        ),
    }


__all__ = [
    "add_negative_intervals",
    "check_genomic_interval_table",
    "design_sequence_ledidi",
    "design_sequences_evolve",
    "expand_chromosome_names",
    "filter_interval_coverage",
    "get_attention_scores",
    "get_sequence_attributions",
    "grelu_marginalize_variants",
    "grelu_predict_variant_effects",
    "in_silico_mutagenesis_scores",
    "marginalize_pattern_insertions",
    "marginalize_pattern_spacing",
    "read_modisco_motifs",
    "read_narrowpeak_intervals",
    "read_tomtom_matches",
    "run_modisco_motifs",
    "score_motif_pwm_products",
    "score_motif_sites",
    "shuffle_tile_effects",
]

_register_public_verbose(globals(), tuple(__all__))
