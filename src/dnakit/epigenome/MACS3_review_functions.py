"""Optional MACS3 wrappers for HMMRATAC: digested fragment signals, HMM models and outputs.

HMMRATAC splits ATAC-seq fragments into short, mono-, di- and tri-nucleosomal signal tracks, trains
a three-state hidden Markov model (open, nucleosomal, background) on bins of those tracks, decodes
candidate regions and writes the open regions. The functions here call the pure-Python and Cython
functions behind that pipeline, not the ``macs3 hmmratac`` command line. Objects such as the
fragment track or the fitted model cannot pass between calls, so each function takes plain lists or
files and returns plain values: fragments as ``(chrom, left, right)`` triples, training data as
rows of four numbers, the model as a JSON file written by MACS3, and the predicted probabilities
as the text file that the MACS3 pipeline keeps in a temporary file. Packages are imported only
when a function runs. A missing package raises ``BackendUnavailableError`` with the
``pip install "dnakit[bio-macs3]"`` hint.
"""

from __future__ import annotations

import contextlib
import io
import json
import math
import re
from collections.abc import Mapping, Sequence
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_MACS3_EXTRA = "bio-macs3"
_MAX_COORDINATE = 2_147_483_647
_MAX_FRAGMENTS = 5_000_000
_MAX_FRAGMENT_LENGTH = 10_000
_MAX_REGIONS = 1_000_000
_MAX_ROWS = 5_000_000
_MAX_PROBA_LINES = 20_000_000
_MAX_STATES = 5_000_000
_MAX_BINSIZE = 1_000_000
_MAX_FLANK = 1_000_000
_MAX_MODEL_BYTES = 10_000_000
_CLASSES = ("short", "mono", "di", "tri")
_STATE_LABELS = ("open", "nuc", "bg")
_HMM_TYPES = ("gaussian", "poisson")
_COVARIANCE_TYPES = ("full", "diag")
_FEATURES = len(_CLASSES)
_STATES = len(_STATE_LABELS)
_CHROM_PATTERN = re.compile(r"[!-+\--~]{1,255}")


def _integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _number(value: object, name: str, minimum: float, maximum: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a number in [{minimum}, {maximum}].")
    number = float(value)
    if not math.isfinite(number) or not minimum <= number <= maximum:
        raise ConfigurationError(f"{name} must be a number in [{minimum}, {maximum}].")
    return number


def _choice(value: object, name: str, choices: tuple[str, ...]) -> str:
    if not isinstance(value, str) or value not in choices:
        options = ", ".join(repr(item) for item in choices)
        raise ConfigurationError(f"{name} must be one of {options}.")
    return value


def _list(value: object, name: str, maximum: int) -> Sequence[Any]:
    if isinstance(value, (str, bytes, Mapping)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list.")
    if len(value) > maximum:
        raise ConfigurationError(f"{name} may contain at most {maximum} entries.")
    return value


def _chrom(value: object, name: str) -> str:
    if not isinstance(value, str) or not _CHROM_PATTERN.fullmatch(value):
        raise ConfigurationError(
            f"{name} must be a chromosome name of 1 to 255 printable ASCII characters "
            "without spaces or commas."
        )
    return value


def _path(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a file path.")
    return input_path(value)


def _new_path(value: object, name: str) -> Path:
    """Output path that must not exist yet; its parent directory must exist."""

    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a file path.")
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} already exists and will not be overwritten: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The directory for {name} does not exist: {path.parent}")
    return path


def _new_paths(prefix: object, name: str, suffixes: Sequence[str]) -> list[Path]:
    """Output paths made from a prefix; none of them may exist yet."""

    if not isinstance(prefix, (str, Path)) or not str(prefix).strip():
        raise ConfigurationError(f"{name} must be a non-empty path prefix.")
    base = Path(prefix).expanduser().resolve()
    if base.is_dir():
        raise ConfigurationError(f"{name} must be a file name prefix, not a directory: {base}")
    if not base.parent.is_dir():
        raise ConfigurationError(f"The directory for {name} does not exist: {base.parent}")
    paths = [base.with_name(base.name + suffix) for suffix in suffixes]
    for path in paths:
        if path.exists():
            raise ConfigurationError(f"{path.name} already exists and will not be overwritten.")
    return paths


def _provenance(statistic: str, *, hmm: bool = False, **notes: str) -> dict[str, Any]:
    packages = ["MACS3"] + (["hmmlearn"] if hmm else [])
    labels = []
    for distribution in packages:
        try:
            installed = version(distribution)
        except PackageNotFoundError:
            installed = "unknown"
        labels.append(f"{distribution} {installed}")
    return {"package": ", ".join(labels), "method": statistic, **notes}


def _upstream_call(function: Any, *args: Any, **kwargs: Any) -> Any:
    """Run a MACS3 function that logs, asserts its inputs and may call sys.exit."""

    captured = io.StringIO()
    try:
        with contextlib.redirect_stdout(captured), contextlib.redirect_stderr(captured):
            return function(*args, **kwargs)
    except (SystemExit, AssertionError, ValueError) as exc:
        lines = captured.getvalue().strip().splitlines()
        detail = str(exc.code if isinstance(exc, SystemExit) else exc)
        if not detail or detail == "None":
            detail = lines[-1] if lines else type(exc).__name__
        raise ConfigurationError(f"MACS3 or hmmlearn rejected the input: {detail}") from exc


def _fragments(value: object) -> list[tuple[bytes, int, int]]:
    items = _list(value, "fragments", _MAX_FRAGMENTS)
    if not items:
        raise ConfigurationError("fragments must contain at least one (chrom, left, right).")
    parsed = []
    for item in items:
        if isinstance(item, (str, bytes, Mapping)) or not isinstance(item, Sequence):
            raise ConfigurationError("Each fragment must be a (chrom, left, right) triple.")
        if len(item) != 3:
            raise ConfigurationError("Each fragment must be a (chrom, left, right) triple.")
        chrom = _chrom(item[0], "fragment chromosome")
        left = _integer(item[1], "fragment left", 0, _MAX_COORDINATE)
        right = _integer(item[2], "fragment right", 1, _MAX_COORDINATE)
        if right <= left:
            raise ConfigurationError(
                "Every fragment must have its right coordinate after its left."
            )
        if right - left > _MAX_FRAGMENT_LENGTH:
            raise ConfigurationError(
                f"Fragments may be at most {_MAX_FRAGMENT_LENGTH} bases long; drop longer ones."
            )
        parsed.append((chrom.encode("ascii"), left, right))
    return parsed


def _four(value: object, name: str, minimum: float, maximum: float) -> list[float]:
    items = _list(value, name, _FEATURES)
    if len(items) != _FEATURES:
        raise ConfigurationError(f"{name} must have four values: short, mono, di, tri.")
    return [_number(item, f"{name} entries", minimum, maximum) for item in items]


def _petrack(fragments: list[tuple[bytes, int, int]]) -> Any:
    """Build and finalize the paired-end track MACS3 digests."""

    paired = require("MACS3.Signal.PairedEndTrack", _MACS3_EXTRA)
    chromosomes = len({chrom for chrom, _, _ in fragments})
    # One array of this size is allocated per chromosome, so keep it near the mean fragment count.
    buffer_size = min(100_000, max(1_024, -(-len(fragments) // chromosomes)))
    track = paired.PETrackI(buffer_size=buffer_size)
    for chrom, left, right in fragments:
        track.add_loc(chrom, left, right)
    track.finalize()
    return track


def _digested(
    fragments: list[tuple[bytes, int, int]],
    means: list[float],
    stddevs: list[float],
    floor: float,
) -> list[Any]:
    """Four bedGraphTrackI objects (short, mono, di, tri) from fragments."""

    from dnakit.epigenome.MACS3_functions import hmmratac_fragment_weights

    processing = require("MACS3.Signal.HMMR_Signal_Processing", _MACS3_EXTRA)
    track = _petrack(fragments)
    lengths = sorted(track.count_fraglengths())
    rows = hmmratac_fragment_weights(
        lengths, means, stddevs, min_fragment_probability=floor, verbose=False
    )["weights"]
    mapping = [{row["fragment_length"]: row[label] for row in rows} for label in _CLASSES]
    return list(_upstream_call(processing.generate_digested_signals, track, mapping))


def _track_arrays(track: Any) -> dict[str, dict[str, list[Any]]]:
    result: dict[str, dict[str, list[Any]]] = {}
    for chrom in sorted(track.get_chr_names()):
        positions, values = track.get_data_by_chr(chrom)
        result[chrom.decode()] = {
            "positions": [int(item) for item in positions],
            "values": [float(item) for item in values],
        }
    return result


def digest_hmmratac_fragments(
    fragments: Sequence[Sequence[Any]],
    means: Sequence[float],
    stddevs: Sequence[float],
    *,
    min_fragment_probability: float = 0.001,
    output_prefix: str | Path | None = None,
    include_tracks: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split paired-end fragments into the short, mono, di and tri signal tracks of HMMRATAC."""

    parsed = _fragments(fragments)
    mu = _four(means, "means", 0.0, 1_000_000.0)
    sigma = _four(stddevs, "stddevs", 1e-9, 1_000_000.0)
    floor = _number(min_fragment_probability, "min_fragment_probability", 0.0, 1.0)
    if not isinstance(include_tracks, bool):
        raise ConfigurationError("include_tracks must be a boolean.")
    targets = (
        None
        if output_prefix is None
        else _new_paths(output_prefix, "output_prefix", [f"_{label}.bdg" for label in _CLASSES])
    )
    tracks = _digested(parsed, mu, sigma, floor)
    if len(tracks) != _FEATURES:
        raise ConfigurationError("MACS3 did not return four digested signal tracks.")
    arrays = [_track_arrays(track) for track in tracks]
    result: dict[str, Any] = {
        "status": "completed",
        "operation": "digest_hmmratac_fragments",
        "fragments": len(parsed),
        "chromosomes": sorted({chrom.decode() for chrom, _, _ in parsed}),
        "classes": list(_CLASSES),
        "intervals": {
            label: sum(len(chrom["positions"]) for chrom in array.values())
            for label, array in zip(_CLASSES, arrays, strict=True)
        },
    }
    if include_tracks:
        result["tracks"] = dict(zip(_CLASSES, arrays, strict=True))
    if targets is not None:
        io_module = require("MACS3.IO.BedGraphIO", _MACS3_EXTRA)
        for label, path, track in zip(_CLASSES, targets, tracks, strict=True):
            writer = io_module.bedGraphIO(str(path), data=track)
            _upstream_call(writer.write_bedGraph, label, label)
        result["files"] = {label: str(path) for label, path in zip(_CLASSES, targets, strict=True)}
    result["provenance"] = _provenance(
        "Gaussian fragment-length weights split the paired-end pileup into four tracks",
        note="in each track a position is the end of a segment; the first segment starts at 0",
    )
    return result


def _regions(value: object) -> list[tuple[bytes, int, int]]:
    items = _list(value, "regions", _MAX_REGIONS)
    if not items:
        raise ConfigurationError("regions must contain at least one (chrom, start, end).")
    parsed = []
    for item in items:
        if isinstance(item, (str, bytes, Mapping)) or not isinstance(item, Sequence):
            raise ConfigurationError("Each region must be a (chrom, start, end) triple.")
        if len(item) != 3:
            raise ConfigurationError("Each region must be a (chrom, start, end) triple.")
        chrom = _chrom(item[0], "region chromosome")
        start = _integer(item[1], "region start", 0, _MAX_COORDINATE)
        end = _integer(item[2], "region end", 1, _MAX_COORDINATE)
        if end <= start:
            raise ConfigurationError("Every region must end after it starts.")
        parsed.append((chrom.encode("ascii"), start, end))
    return parsed


def _rows(value: object, name: str, *, integers: bool) -> list[list[float]]:
    """Rows of four non-negative numbers, whole numbers when ``integers`` is set."""

    items = _list(value, name, _MAX_ROWS)
    if not items:
        raise ConfigurationError(f"{name} must contain at least one row.")
    rows: list[list[Any]] = []
    for item in items:
        if isinstance(item, (str, bytes, Mapping)) or not isinstance(item, Sequence):
            raise ConfigurationError(f"Each row of {name} must be a list of four numbers.")
        if len(item) != _FEATURES:
            raise ConfigurationError(
                f"Each row of {name} must have four values: short, mono, di, tri."
            )
        if integers:
            rows.append([_integer(entry, f"{name} entries", 0, _MAX_COORDINATE) for entry in item])
        else:
            rows.append([_number(entry, f"{name} entries", 0.0, 1e12) for entry in item])
    return rows


def _lengths(value: object, name: str, rows: int) -> list[int]:
    items = _list(value, name, _MAX_ROWS)
    lengths = [_integer(item, f"{name} entries", 1, _MAX_ROWS) for item in items]
    if not lengths or sum(lengths) != rows:
        raise ConfigurationError(
            f"{name} must list positive sequence lengths that add up to the {rows} rows."
        )
    return lengths


def extract_hmmratac_training_data(
    fragments: Sequence[Sequence[Any]],
    regions: Sequence[Sequence[Any]],
    means: Sequence[float],
    stddevs: Sequence[float],
    *,
    binsize: int = 10,
    flank: int = 0,
    hmm_type: str = "gaussian",
    min_fragment_probability: float = 0.001,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cut the four digested signal tracks into bins inside regions, as HMMRATAC does."""

    parsed = _fragments(fragments)
    wanted = _regions(regions)
    mu = _four(means, "means", 0.0, 1_000_000.0)
    sigma = _four(stddevs, "stddevs", 1e-9, 1_000_000.0)
    size = _integer(binsize, "binsize", 1, _MAX_BINSIZE)
    expand = _integer(flank, "flank", 0, _MAX_FLANK)
    kind = _choice(hmm_type, "hmm_type", _HMM_TYPES)
    floor = _number(min_fragment_probability, "min_fragment_probability", 0.0, 1.0)
    region_module = require("MACS3.Signal.Region", _MACS3_EXTRA)
    processing = require("MACS3.Signal.HMMR_Signal_Processing", _MACS3_EXTRA)
    tracks = _digested(parsed, mu, sigma, floor)
    merged = region_module.Regions()
    for chrom, start, end in wanted:
        merged.add_loc(chrom, start, end)
    if expand:
        merged.expand(expand)
    merged.merge_overlap()
    try:
        bins, data, lengths = _upstream_call(
            processing.extract_signals_from_regions, tracks, merged, binsize=size, hmm_type=kind
        )
    except ConfigurationError as exc:
        raise ConfigurationError(
            "No signal bin falls inside the regions; check the chromosome names and that the "
            f"regions cover at least one complete bin of binsize={size} within the fragments."
        ) from exc
    return {
        "status": "completed",
        "operation": "extract_hmmratac_training_data",
        "binsize": size,
        "hmm_type": kind,
        "regions": int(merged.total),
        "bins": [[chrom.decode(), int(end)] for chrom, end in bins],
        "data": [
            [float(entry) if kind == "gaussian" else int(entry) for entry in row] for row in data
        ],
        "lengths": [int(item) for item in lengths],
        "provenance": _provenance(
            "bins of the short, mono, di and tri tracks inside the merged regions",
            note=(
                "a bin is listed by its end position and covers [end - binsize, end); "
                "values below 0.0001 are raised to 0.0001; chromosomes come in reverse order"
            ),
        ),
    }


def _fit(
    data: Sequence[Sequence[float]],
    lengths: Sequence[int],
    hmm_type: str,
    covariance_type: str,
    random_seed: int,
) -> tuple[Any, list[int], bool]:
    """Fit the three-state model and label its states as the HMMRATAC pipeline does."""

    kind = _choice(hmm_type, "hmm_type", _HMM_TYPES)
    covariance = _choice(covariance_type, "covariance_type", _COVARIANCE_TYPES)
    seed = _integer(random_seed, "random_seed", 0, _MAX_COORDINATE)
    rows = _rows(data, "data", integers=kind == "poisson")
    groups = _lengths(lengths, "lengths", len(rows))
    numpy = require("numpy", _MACS3_EXTRA)
    model_module = require("MACS3.Signal.HMMR_HMM", _MACS3_EXTRA)
    model = _upstream_call(
        model_module.hmm_training,
        rows,
        groups,
        random_seed=seed,
        covar=covariance,
        hmm_type=kind,
    )
    emissions = model.means_ if kind == "gaussian" else model.lambdas_
    if not numpy.isfinite(emissions).all():
        raise ConfigurationError(
            "The fit gave non-finite emissions; use more rows or less repetitive values."
        )
    if kind == "gaussian" and any(
        float(numpy.linalg.eigvalsh(matrix).min()) <= 0.0 for matrix in model.covars_
    ):
        raise ConfigurationError(
            "A fitted covariance matrix is not positive definite, so the model cannot be "
            "saved or reloaded. A signal track is constant or nearly so; use more varied data."
        )
    totals = numpy.sum(emissions, axis=1)
    open_state = int(numpy.where(totals == max(totals))[0][0])
    background_state = int(numpy.where(totals == min(totals))[0][0])
    remaining = set(range(_STATES)) - {open_state, background_state}
    if open_state == background_state or len(remaining) != 1:
        raise ConfigurationError(
            "The fitted states have equal emissions, so open, nucleosomal and background "
            "states cannot be told apart; use more varied data."
        )
    nucleosomal_state = remaining.pop()
    return model, [open_state, nucleosomal_state, background_state], bool(model.monitor_.converged)


def _model_values(model: Any, hmm_type: str) -> dict[str, Any]:
    values: dict[str, Any] = {
        "start_probabilities": [float(item) for item in model.startprob_],
        "transition_matrix": [[float(item) for item in row] for row in model.transmat_],
    }
    if hmm_type == "gaussian":
        values["covariance_type"] = str(model.covariance_type)
        values["means"] = [[float(item) for item in row] for row in model.means_]
        values["covariances"] = [
            [[float(item) for item in row] for row in matrix] for matrix in model.covars_
        ]
    else:
        values["lambdas"] = [[float(item) for item in row] for row in model.lambdas_]
    return values


def train_hmmratac_model(
    data: Sequence[Sequence[float]],
    lengths: Sequence[int],
    *,
    hmm_type: str = "gaussian",
    covariance_type: str = "full",
    random_seed: int = 12345,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fit the three-state HMMRATAC model and label its open, nucleosomal and background states."""

    model, indexes, converged = _fit(data, lengths, hmm_type, covariance_type, random_seed)
    kind = str(hmm_type)
    history = list(model.monitor_.history)
    return {
        "status": "completed",
        "operation": "train_hmmratac_model",
        "hmm_type": kind,
        "converged": converged,
        "iterations": int(model.monitor_.iter),
        "log_likelihood": float(history[-1]) if history else None,
        "open_state": indexes[0],
        "nucleosomal_state": indexes[1],
        "background_state": indexes[2],
        **_model_values(model, kind),
        "provenance": _provenance(
            "Baum-Welch fit of a three-state model; the state with the largest summed emission "
            "is open, the smallest is background, and the other is nucleosomal",
            hmm=True,
            note="training is seeded by random_seed and does not touch global random state",
        ),
    }


def save_hmmratac_model(
    data: Sequence[Sequence[float]],
    lengths: Sequence[int],
    output_path: str | Path,
    *,
    binsize: int = 10,
    hmm_type: str = "gaussian",
    covariance_type: str = "full",
    random_seed: int = 12345,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fit the HMMRATAC model and write it, with its state labels and bin size, to a JSON file."""

    size = _integer(binsize, "binsize", 1, _MAX_BINSIZE)
    target = _new_path(output_path, "output_path")
    model, indexes, converged = _fit(data, lengths, hmm_type, covariance_type, random_seed)
    model_module = require("MACS3.Signal.HMMR_HMM", _MACS3_EXTRA)
    _upstream_call(
        model_module.hmm_model_save,
        str(target),
        model,
        size,
        indexes[0],
        indexes[1],
        indexes[2],
        hmm_type,
    )
    if not target.is_file():
        raise ConfigurationError(f"MACS3 did not write the model file {target.name}.")
    try:
        _load_model(target)
    except ConfigurationError:
        target.unlink()
        raise
    return {
        "status": "completed",
        "operation": "save_hmmratac_model",
        "output": str(target),
        "hmm_type": hmm_type,
        "binsize": size,
        "converged": converged,
        "open_state": indexes[0],
        "nucleosomal_state": indexes[1],
        "background_state": indexes[2],
        "provenance": _provenance(
            "MACS3 JSON model: start and transition probabilities, emissions, state indexes and "
            "bin size",
            hmm=True,
            note="read it back with load_hmmratac_model or use it in predict_hmmratac_states",
        ),
    }


def _shaped(value: object, shape: tuple[int, ...]) -> bool:
    """Whether ``value`` is nested lists of the given shape holding finite numbers."""

    if not shape:
        return (
            not isinstance(value, bool) and isinstance(value, (int, float)) and math.isfinite(value)
        )
    return (
        isinstance(value, list)
        and len(value) == shape[0]
        and all(_shaped(item, shape[1:]) for item in value)
    )


def _model_file(path: Path) -> dict[str, Any]:
    """Read and check a MACS3 HMMRATAC model file before MACS3 builds the model from it."""

    if path.stat().st_size > _MAX_MODEL_BYTES:
        raise ConfigurationError(f"{path.name} is larger than {_MAX_MODEL_BYTES} bytes.")
    try:
        content = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ConfigurationError(f"{path.name} is not a JSON model file: {exc}") from exc
    if not isinstance(content, dict):
        raise ConfigurationError(f"{path.name} must hold a JSON object.")
    kind = content.get("hmm_type", "gaussian")
    if kind not in _HMM_TYPES:
        raise ConfigurationError(f"hmm_type in {path.name} must be 'gaussian' or 'poisson'.")
    shapes: dict[str, tuple[int, ...]] = {
        "startprob": (_STATES,),
        "transmat": (_STATES, _STATES),
    }
    if kind == "gaussian":
        covariance = content.get("covariance_type")
        if covariance not in _COVARIANCE_TYPES:
            raise ConfigurationError(f"covariance_type in {path.name} must be 'full' or 'diag'.")
        shapes["means"] = (_STATES, _FEATURES)
        shapes["covars"] = (
            (_STATES, _FEATURES) if covariance == "diag" else (_STATES, _FEATURES, _FEATURES)
        )
    else:
        shapes["lambdas"] = (_STATES, _FEATURES)
    for key, shape in shapes.items():
        if not _shaped(content.get(key), shape):
            dimensions = " x ".join(str(size) for size in shape)
            raise ConfigurationError(
                f"{key} in {path.name} must be a {dimensions} array of numbers."
            )
    if content.get("n_features") != _FEATURES:
        raise ConfigurationError(f"n_features in {path.name} must be {_FEATURES}.")
    indexes = [
        content.get(key) for key in ("i_open_region", "i_nucleosomal_region", "i_background_region")
    ]
    if any(isinstance(index, bool) or not isinstance(index, int) for index in indexes) or sorted(
        indexes
    ) != [0, 1, 2]:
        raise ConfigurationError(
            f"{path.name} must give distinct open, nucleosomal and background state indexes 0 to 2."
        )
    _integer(content.get("hmm_binsize"), "hmm_binsize in the model file", 1, _MAX_BINSIZE)
    return content


def _load_model(path: Path) -> tuple[Any, int, int, int, int, str]:
    model_module = require("MACS3.Signal.HMMR_HMM", _MACS3_EXTRA)
    _model_file(path)
    model, open_state, background_state, nucleosomal_state, binsize, kind = _upstream_call(
        model_module.hmm_model_init, str(path)
    )
    return model, int(open_state), int(nucleosomal_state), int(background_state), int(binsize), kind


def load_hmmratac_model(model_path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Read a MACS3 HMMRATAC model file: probabilities, emissions, state indexes and bin size."""

    source = _path(model_path, "model_path")
    model, open_state, nucleosomal_state, background_state, binsize, kind = _load_model(source)
    labels = [""] * _STATES
    labels[open_state] = "open"
    labels[nucleosomal_state] = "nuc"
    labels[background_state] = "bg"
    return {
        "status": "completed",
        "operation": "load_hmmratac_model",
        "model": str(source),
        "hmm_type": kind,
        "binsize": binsize,
        "open_state": open_state,
        "nucleosomal_state": nucleosomal_state,
        "background_state": background_state,
        "state_labels": labels,
        "features": list(_CLASSES),
        **_model_values(model, kind),
        "provenance": _provenance(
            "parameters of the model rebuilt from the JSON file",
            hmm=True,
            note="a file without hmm_type is read as a Gaussian model",
        ),
    }


def _bins(value: object, rows: int) -> list[tuple[bytes, int]]:
    items = _list(value, "bins", _MAX_ROWS)
    if len(items) != rows:
        raise ConfigurationError("bins must have one (chrom, end) pair per row of data.")
    parsed = []
    for item in items:
        if isinstance(item, (str, bytes, Mapping)) or not isinstance(item, Sequence):
            raise ConfigurationError("Each bin must be a (chrom, end) pair.")
        if len(item) != 2:
            raise ConfigurationError("Each bin must be a (chrom, end) pair.")
        chrom = _chrom(item[0], "bin chromosome")
        parsed.append((chrom.encode("ascii"), _integer(item[1], "bin end", 1, _MAX_COORDINATE)))
    return parsed


def predict_hmmratac_states(
    model_path: str | Path,
    data: Sequence[Sequence[float]],
    lengths: Sequence[int],
    *,
    bins: Sequence[Sequence[Any]] | None = None,
    output_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Give each bin its probability of being an open, nucleosomal or background state."""

    source = _path(model_path, "model_path")
    if (bins is None) != (output_path is None):
        raise ConfigurationError(
            "bins and output_path go together: give both to write the probability file."
        )
    target = None if output_path is None else _new_path(output_path, "output_path")
    model, open_state, nucleosomal_state, background_state, binsize, kind = _load_model(source)
    rows = _rows(data, "data", integers=kind == "poisson")
    groups = _lengths(lengths, "lengths", len(rows))
    named = None if bins is None else _bins(bins, len(rows))
    model_module = require("MACS3.Signal.HMMR_HMM", _MACS3_EXTRA)
    raw = _upstream_call(model_module.hmm_predict, rows, groups, model)
    table = [[float(item) for item in row] for row in raw]
    if any(not math.isfinite(item) for row in table for item in row):
        raise ConfigurationError("The model gave non-finite probabilities for these rows.")
    probabilities = [
        [row[open_state], row[nucleosomal_state], row[background_state]] for row in table
    ]
    labels = [
        _STATE_LABELS[max(range(_STATES), key=lambda index, row=row: row[index])]
        for row in probabilities
    ]
    result: dict[str, Any] = {
        "status": "completed",
        "operation": "predict_hmmratac_states",
        "model": str(source),
        "hmm_type": kind,
        "binsize": binsize,
        "open_state": open_state,
        "nucleosomal_state": nucleosomal_state,
        "background_state": background_state,
        "rows": len(rows),
        "columns": ["open", "nucleosomal", "background"],
        "probabilities": probabilities,
        "labels": labels,
        "label_counts": {label: labels.count(label) for label in _STATE_LABELS},
    }
    if target is not None and named is not None:
        with target.open("wb") as handle:
            for (chrom, end), row in zip(named, table, strict=True):
                handle.write(b"%s,%d" % (chrom, end))
                handle.write(b",%f,%f,%f\n" % tuple(row))
        result["output"] = str(target)
    result["provenance"] = _provenance(
        "posterior probability of each model state for every bin (forward-backward algorithm)",
        hmm=True,
        note=(
            "probabilities are ordered open, nucleosomal, background; the file written to "
            "output_path keeps the model's own state order, six decimals per value"
        ),
    )
    return result


def _scan_proba(path: Path, binsize: int) -> tuple[int, list[str]]:
    """Check a predicted-probability file (chrom,end,p0,p1,p2) and list its chromosomes."""

    chromosomes: list[str] = []
    previous_end = 0
    count = 0
    with path.open("rb") as handle:
        for number, line in enumerate(handle, start=1):
            fields = line.strip().split(b",")
            if len(fields) != 2 + _STATES:
                raise ConfigurationError(
                    f"{path.name} line {number} must have chrom, end and three probabilities."
                )
            try:
                chrom = fields[0].decode("ascii")
                end = int(fields[1])
                values = [float(item) for item in fields[2:]]
            except ValueError as exc:
                raise ConfigurationError(
                    f"{path.name} line {number} has a malformed value."
                ) from exc
            if not _CHROM_PATTERN.fullmatch(chrom) or end < binsize or end > _MAX_COORDINATE:
                raise ConfigurationError(
                    f"{path.name} line {number} has an invalid chromosome or bin end."
                )
            if any(not math.isfinite(item) or item < 0.0 for item in values):
                raise ConfigurationError(f"{path.name} line {number} has an invalid probability.")
            if not chromosomes or chromosomes[-1] != chrom:
                if chrom in chromosomes:
                    raise ConfigurationError(
                        f"{path.name} line {number}: rows of {chrom} must stay together."
                    )
                chromosomes.append(chrom)
            elif end - binsize < previous_end:
                raise ConfigurationError(f"{path.name} line {number}: bins must not overlap.")
            previous_end = end
            count += 1
            if count > _MAX_PROBA_LINES:
                raise ConfigurationError(f"{path.name} has more than {_MAX_PROBA_LINES} lines.")
    if not count:
        raise ConfigurationError(f"{path.name} is empty.")
    return count, chromosomes


def _state_indexes(
    open_state: object, nucleosomal_state: object, background_state: object
) -> list[int]:
    indexes = [
        _integer(open_state, "open_state", 0, _STATES - 1),
        _integer(nucleosomal_state, "nucleosomal_state", 0, _STATES - 1),
        _integer(background_state, "background_state", 0, _STATES - 1),
    ]
    if sorted(indexes) != [0, 1, 2]:
        raise ConfigurationError(
            "open_state, nucleosomal_state and background_state must be 0, 1 and 2 in some order."
        )
    return indexes


def hmmratac_proba_to_bedgraph(
    proba_path: str | Path,
    output_prefix: str | Path,
    *,
    open_state: int,
    nucleosomal_state: int,
    background_state: int,
    binsize: int = 10,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the state probabilities of a prediction file as one bedGraph per state."""

    source = _path(proba_path, "proba_path")
    size = _integer(binsize, "binsize", 1, _MAX_BINSIZE)
    indexes = _state_indexes(open_state, nucleosomal_state, background_state)
    targets = _new_paths(output_prefix, "output_prefix", ["_open.bdg", "_nuc.bdg", "_bg.bdg"])
    bins, chromosomes = _scan_proba(source, size)
    command = require("MACS3.Commands.hmmratac_cmd", _MACS3_EXTRA)
    with source.open("rb") as handle:
        _upstream_call(
            command.save_proba_to_bedGraph,
            handle,
            size,
            str(targets[0]),
            str(targets[1]),
            str(targets[2]),
            indexes[0],
            indexes[1],
            indexes[2],
        )
    files: dict[str, str] = {}
    intervals: dict[str, int] = {}
    for label, path in zip(("open", "nucleosomal", "background"), targets, strict=True):
        if not path.is_file():
            raise ConfigurationError(f"MACS3 did not write {path.name}.")
        with path.open("r", encoding="utf-8") as handle:
            intervals[label] = sum(1 for _ in handle)
        files[label] = str(path)
    return {
        "status": "completed",
        "operation": "hmmratac_proba_to_bedgraph",
        "proba": str(source),
        "binsize": size,
        "bin_count": bins,
        "chromosomes": chromosomes,
        "files": files,
        "intervals": intervals,
        "provenance": _provenance(
            "one bedGraph per state; gaps between bins and before the first bin are background",
            note="a bin ending at E covers [E - binsize, E); background is 1.0 in the gaps",
        ),
    }


def _states_result(path: list[tuple[bytes, int, int, str]]) -> list[dict[str, Any]]:
    return [
        {"chrom": chrom.decode(), "start": int(start), "end": int(end), "state": str(state)}
        for chrom, start, end, state in path
    ]


def hmmratac_states_path(
    proba_path: str | Path,
    *,
    open_state: int,
    nucleosomal_state: int,
    background_state: int,
    binsize: int = 10,
    verbose: bool = True,
) -> dict[str, Any]:
    """Label each bin of a prediction file with its most probable state and merge equal runs."""

    source = _path(proba_path, "proba_path")
    size = _integer(binsize, "binsize", 1, _MAX_BINSIZE)
    indexes = _state_indexes(open_state, nucleosomal_state, background_state)
    bins, chromosomes = _scan_proba(source, size)
    command = require("MACS3.Commands.hmmratac_cmd", _MACS3_EXTRA)
    with source.open("rb") as handle:
        path = _upstream_call(
            command.generate_states_path, handle, size, indexes[0], indexes[1], indexes[2]
        )
    states = _states_result(path)
    return {
        "status": "completed",
        "operation": "hmmratac_states_path",
        "proba": str(source),
        "binsize": size,
        "bin_count": bins,
        "chromosomes": chromosomes,
        "states": states,
        "state_counts": {
            label: sum(1 for item in states if item["state"] == label) for label in _STATE_LABELS
        },
        "state_bases": {
            label: sum(item["end"] - item["start"] for item in states if item["state"] == label)
            for label in _STATE_LABELS
        },
        "provenance": _provenance(
            "the most probable state per bin (ties go to open, then nucleosomal); gaps between "
            "bins and before a chromosome's first bin become background",
            note="neighbouring background entries are not merged across a gap",
        ),
    }


def _state_entries(value: object) -> list[tuple[bytes, int, int, str]]:
    items = _list(value, "states", _MAX_STATES)
    entries: list[tuple[bytes, int, int, str]] = []
    for item in items:
        if isinstance(item, Mapping):
            fields = (item.get("chrom"), item.get("start"), item.get("end"), item.get("state"))
        elif isinstance(item, Sequence) and not isinstance(item, (str, bytes)) and len(item) == 4:
            fields = tuple(item)
        else:
            raise ConfigurationError(
                "Each state must be a mapping with chrom, start, end and state, or a 4-item list."
            )
        chrom = _chrom(fields[0], "state chromosome")
        start = _integer(fields[1], "state start", 0, _MAX_COORDINATE)
        end = _integer(fields[2], "state end", 1, _MAX_COORDINATE)
        label = _choice(fields[3], "state label", _STATE_LABELS)
        if end <= start:
            raise ConfigurationError("Every state must end after it starts.")
        if entries and entries[-1][0] == chrom.encode("ascii") and start < entries[-1][2]:
            raise ConfigurationError("States must be sorted along each chromosome without overlap.")
        entries.append((chrom.encode("ascii"), start, end, label))
    return entries


def write_hmmratac_states_bed(
    states: Sequence[Mapping[str, Any] | Sequence[Any]],
    output_path: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the open and nucleosomal runs of a states path as a four-column BED file."""

    entries = _state_entries(states)
    target = _new_path(output_path, "output_path")
    command = require("MACS3.Commands.hmmratac_cmd", _MACS3_EXTRA)
    with target.open("w", encoding="utf-8") as handle:
        _upstream_call(command.save_states_bed, entries, handle)
    written = {label: sum(1 for entry in entries if entry[3] == label) for label in _STATE_LABELS}
    return {
        "status": "completed",
        "operation": "write_hmmratac_states_bed",
        "output": str(target),
        "records": written["open"] + written["nuc"],
        "state_counts": {"open": written["open"], "nuc": written["nuc"]},
        "skipped_background": written["bg"],
        "provenance": _provenance(
            "BED rows chrom, start, end, state for every open and nuc entry",
            note="background entries are not written",
        ),
    }


def _bedgraph_file(path: Path) -> None:
    """Check that a bedGraph file has sorted, contiguous intervals per chromosome."""

    previous: dict[str, int] = {}
    last_chrom = ""
    intervals = 0
    with path.open("r", encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            if line.startswith(("track", "browser", "#")) or not line.strip():
                continue
            fields = line.split()
            if len(fields) < 4:
                raise ConfigurationError(f"{path.name} line {number} must have four fields.")
            try:
                start, end, value = int(fields[1]), int(fields[2]), float(fields[3])
            except ValueError as exc:
                raise ConfigurationError(
                    f"{path.name} line {number} has a malformed value."
                ) from exc
            chrom = fields[0]
            if end <= start or start < 0 or end > _MAX_COORDINATE or not math.isfinite(value):
                raise ConfigurationError(f"{path.name} line {number} has an invalid interval.")
            if chrom in previous and (chrom != last_chrom or start != previous[chrom]):
                raise ConfigurationError(
                    f"{path.name} line {number}: intervals must be sorted and contiguous."
                )
            previous[chrom] = end
            last_chrom = chrom
            intervals += 1
    if not intervals:
        raise ConfigurationError(f"{path.name} has no intervals.")


def write_hmmratac_accessible_regions(
    states: Sequence[Mapping[str, Any] | Sequence[Any]],
    score_bedgraph: str | Path,
    output_path: str | Path,
    *,
    min_length: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write nucleosomal-open-nucleosomal runs of a states path as a scored narrowPeak file."""

    entries = _state_entries(states)
    scores = _path(score_bedgraph, "score_bedgraph")
    minimum = _integer(min_length, "min_length", 0, _MAX_COORDINATE)
    target = _new_path(output_path, "output_path")
    _bedgraph_file(scores)
    command = require("MACS3.Commands.hmmratac_cmd", _MACS3_EXTRA)
    bedgraph = require("MACS3.IO.BedGraphIO", _MACS3_EXTRA)
    track = bedgraph.bedGraphIO(str(scores)).read_bedGraph(0)
    with target.open("w", encoding="utf-8") as handle:
        _upstream_call(command.save_accessible_regions, entries, handle, minimum, track)
    regions = []
    with target.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            regions.append(
                {
                    "chrom": fields[0],
                    "start": int(fields[1]),
                    "end": int(fields[2]),
                    "name": fields[3],
                    "score": int(fields[4]),
                    "summit": int(fields[1]) + int(fields[9]),
                }
            )
    return {
        "status": "completed",
        "operation": "write_hmmratac_accessible_regions",
        "output": str(target),
        "score_bedgraph": str(scores),
        "min_length": minimum,
        "regions": regions,
        "count": len(regions),
        "provenance": _provenance(
            "open runs flanked by nucleosomal runs, spanning more than min_length bases; the "
            "score is ten times the largest bedGraph value in the region",
            note=(
                "MACS3 3.0.5 leaves out the last accessible region of the states list "
                "(accessible_regions[:-1] in save_accessible_regions)"
            ),
        ),
    }


__all__ = [
    "digest_hmmratac_fragments",
    "extract_hmmratac_training_data",
    "hmmratac_proba_to_bedgraph",
    "hmmratac_states_path",
    "load_hmmratac_model",
    "predict_hmmratac_states",
    "save_hmmratac_model",
    "train_hmmratac_model",
    "write_hmmratac_accessible_regions",
    "write_hmmratac_states_bed",
]

_register_public_verbose(globals(), tuple(__all__))
