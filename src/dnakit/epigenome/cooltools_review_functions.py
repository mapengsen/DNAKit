"""Dot calling, track alignment, pairwise block sums, decoys and eigen functions from cooltools.

Each public function imports cooltools (and cooler, numpy, pandas) only when it is called.
Missing dependencies raise ``BackendUnavailableError`` with the ``pip install "dnakit[...]"``
hint. Inputs are checked before upstream calls; upstream ``ValueError``, ``KeyError``,
``TypeError`` and ``OSError`` are reported as ``ConfigurationError``.

Eigen-decomposition and decoy functions use a private random state for their call. The
numpy global generator and the C library ``random()`` state (glibc, BSD) are saved before the
call and restored afterwards, and the ARPACK start vector is fixed while scipy's ``eigsh`` and
``eigs`` are patched. The patch is restored in ``finally``. Patching is not thread-safe, and
on platforms without ``initstate``/``setstate`` the C generator is not isolated.
"""

from __future__ import annotations

import contextlib
import ctypes
import math
import warnings
from collections.abc import Iterator, Mapping, Sequence
from importlib.metadata import PackageNotFoundError, version
from itertools import pairwise
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_COOLER_EXTRA = "bio-cooler"
_COOLTOOLS_EXTRA = "bio-cooltools"
_NUMPY_EXTRA = "bio-analysis"
_MAX_BINS = 50_000
_MAX_TRACK_ROWS = 1_000_000
_MAX_DENSE = 2_000
_MAX_EIGS = 10
_MAX_WINDOWS = 10_000
_MAX_PARTITIONS = 10_000
_MAX_POSITION = 2_000_000_000
_MAX_SEED = 2**31 - 1
_C_STATE_BYTES = 256


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


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or any(c in value for c in "\x00\r\n\t"):
        raise ConfigurationError(f"{name} must be a non-empty string without tabs or newlines.")
    return value


def _existing(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a file path.")
    return input_path(value)


def _plain(value: Any) -> Any:
    """Convert numpy scalars to Python values and non-finite floats to None."""

    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _records(frame: Any) -> list[dict[str, Any]]:
    return [
        {str(key): _plain(value) for key, value in row.items()}
        for row in frame.to_dict(orient="records")
    ]


def _matrix(values: Any) -> list[list[Any]]:
    return [[_plain(float(value)) for value in row] for row in values]


def _vector(values: Any) -> list[Any]:
    return [_plain(float(value)) for value in values]


def _provenance(distribution: str, method: str, **notes: str) -> dict[str, Any]:
    try:
        installed = version(distribution)
    except PackageNotFoundError:
        installed = "unknown"
    return {"package": f"{distribution} {installed}", "method": method, **notes}


def _open(cool_path: Path) -> tuple[Any, Any]:
    cooler = require("cooler", _COOLER_EXTRA)
    try:
        clr = cooler.Cooler(str(cool_path))
        nbins = int(clr.info["nbins"])
    except (OSError, KeyError, ValueError) as exc:
        raise ConfigurationError(
            f"Input is not a readable cooler file: {cool_path} ({exc})"
        ) from exc
    if nbins > _MAX_BINS:
        raise ConfigurationError(f"The cooler has {nbins} bins, above the limit of {_MAX_BINS}.")
    return cooler, clr


def _weight(clr: Any, value: object) -> str:
    name = _text(value, "weight")
    if name not in clr.bins()[:1].columns:
        raise ConfigurationError(
            f"This cooler has no {name!r} column; run balance_contact_matrix() first."
        )
    return name


def _square(values: Any, name: str, *, allow_nan: bool = False) -> Any:
    numpy = require("numpy", _NUMPY_EXTRA)
    try:
        rows = [
            [float("nan") if (item is None and allow_nan) else item for item in row]
            for row in values
        ]
        array = numpy.asarray(rows, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"{name} must be a square numeric matrix.") from exc
    if array.ndim != 2 or array.shape[0] != array.shape[1]:
        raise ConfigurationError(f"{name} must be a square matrix.")
    if not 2 <= array.shape[0] <= _MAX_DENSE:
        raise ConfigurationError(f"{name} must be 2 to {_MAX_DENSE} rows wide.")
    if not allow_nan and not numpy.isfinite(array).all():
        raise ConfigurationError(f"{name} must contain only finite values.")
    return array


def _symmetric(array: Any, name: str) -> Any:
    numpy = require("numpy", _NUMPY_EXTRA)
    if not numpy.allclose(array, array.T, rtol=0, atol=1e-9):
        raise ConfigurationError(f"{name} must be symmetric.")
    return array


def _canonical_rows(rows: Any) -> Any:
    """Flip each eigenvector so its largest finite entry is positive."""

    numpy = require("numpy", _NUMPY_EXTRA)
    array = numpy.array(rows, dtype=float, copy=True)
    for index, row in enumerate(array):
        finite = numpy.isfinite(row)
        if finite.any():
            pivot = int(numpy.argmax(numpy.where(finite, numpy.abs(row), -1.0)))
            if row[pivot] < 0:
                array[index] = -row
    return array


def _c_random_api() -> tuple[Any, Any] | None:
    """glibc/BSD ``initstate`` and ``setstate``, or None when the C library lacks them."""

    try:
        libc = ctypes.CDLL(None)
        initstate = libc.initstate
        setstate = libc.setstate
    except (OSError, AttributeError):
        return None
    initstate.argtypes = [ctypes.c_uint, ctypes.c_void_p, ctypes.c_size_t]
    initstate.restype = ctypes.c_void_p
    setstate.argtypes = [ctypes.c_void_p]
    setstate.restype = ctypes.c_void_p
    return initstate, setstate


def _fixed_start_vector(function: Any) -> Any:
    """Wrap scipy eigsh/eigs so that ARPACK gets a fixed start vector when none is given."""

    numpy = require("numpy", _NUMPY_EXTRA)

    def seeded(matrix: Any, *args: Any, **kwargs: Any) -> Any:
        # positional order of eigsh/eigs: k, M, sigma, which, v0
        if "v0" not in kwargs and len(args) <= 4:
            kwargs["v0"] = numpy.random.default_rng(0).random(matrix.shape[0]) + 0.1
        return function(matrix, *args, **kwargs)

    return seeded


@contextlib.contextmanager
def _isolated_random(seed: int, *, eigen: bool) -> Iterator[None]:
    """Run one upstream call with a private random state and restore the caller's state.

    The numpy global generator is saved with ``get_state`` and restored with ``set_state``.
    The C library ``random()`` generator is switched to a private buffer seeded with ``seed``,
    and the caller's buffer, including its position, is restored with ``setstate``. With
    ``eigen=True`` scipy's ``eigsh`` and ``eigs`` get a fixed ARPACK start vector for the call.
    """

    numpy = require("numpy", _NUMPY_EXTRA)
    numpy_state = numpy.random.get_state()
    api = _c_random_api()
    previous: int | None = None
    buffer: Any = None
    originals: dict[str, Any] = {}
    sparse_linalg: Any = None
    if api is not None:
        buffer = (ctypes.c_char * _C_STATE_BYTES)()
        previous = api[0](seed, buffer, _C_STATE_BYTES)
    try:
        if eigen:
            sparse_linalg = require("scipy.sparse.linalg", _NUMPY_EXTRA)
            originals = {name: getattr(sparse_linalg, name) for name in ("eigsh", "eigs")}
            for name, function in originals.items():
                setattr(sparse_linalg, name, _fixed_start_vector(function))
        yield
    finally:
        for name, function in originals.items():
            setattr(sparse_linalg, name, function)
        if api is not None and previous is not None:
            api[1](previous)
        numpy.random.set_state(numpy_state)


def _track(path: Path) -> Any:
    pandas = require("pandas", _COOLER_EXTRA)
    try:
        table = pandas.read_csv(
            path,
            sep="\t",
            comment="#",
            header=None,
            names=["chrom", "start", "end", "value"],
            usecols=[0, 1, 2, 3],
            dtype={"chrom": str, "start": "int64", "end": "int64", "value": "float64"},
        )
    except (ValueError, OSError) as exc:
        raise ConfigurationError(
            f"Track must be a tab-separated file with chrom, start, end and value: {exc}"
        ) from exc
    if table.empty or len(table) > _MAX_TRACK_ROWS:
        raise ConfigurationError(f"The track must have 1 to {_MAX_TRACK_ROWS} rows.")
    if (table["end"] <= table["start"]).any() or (table["start"] < 0).any():
        raise ConfigurationError("Every track interval must have 0 <= start < end.")
    return table


def _bed(path: Path, name: str) -> Any:
    """BED intervals with chrom, start, end and an optional name column (coordinates if absent)."""

    pandas = require("pandas", _COOLER_EXTRA)
    try:
        table = pandas.read_csv(path, sep="\t", comment="#", header=None, dtype=str)
    except (ValueError, OSError) as exc:
        raise ConfigurationError(f"{name} must be a tab-separated BED file: {exc}") from exc
    if table.empty or len(table) > _MAX_WINDOWS or table.shape[1] < 3:
        raise ConfigurationError(
            f"{name} must have 1 to {_MAX_WINDOWS} rows of at least 3 columns."
        )
    frame = pandas.DataFrame(
        {
            "chrom": table[0].astype(str),
            "start": pandas.to_numeric(table[1]).astype("int64"),
            "end": pandas.to_numeric(table[2]).astype("int64"),
        }
    )
    if (frame["end"] <= frame["start"]).any() or (frame["start"] < 0).any():
        raise ConfigurationError(f"Every interval in {name} must have 0 <= start < end.")
    if table.shape[1] >= 4:
        frame["name"] = table[3].astype(str)
    else:
        frame["name"] = [
            f"{c}:{s}-{e}"
            for c, s, e in zip(frame["chrom"], frame["start"], frame["end"], strict=True)
        ]
    return frame


def _bin_regions(bins: Any, view: Any) -> Any:
    """View region index of each bin (NaN when no region contains the bin)."""

    numpy = require("numpy", _NUMPY_EXTRA)
    labels = numpy.full(len(bins), numpy.nan)
    chroms = bins["chrom"].astype(str).to_numpy()
    starts = bins["start"].to_numpy()
    ends = bins["end"].to_numpy()
    for index, row in enumerate(view[["chrom", "start", "end"]].itertuples(index=False)):
        inside = (chroms == str(row[0])) & (starts >= row[1]) & (ends <= row[2])
        labels[inside & numpy.isnan(labels)] = index
    return labels


def call_chromatin_dots(
    cool_path: str | Path,
    *,
    weight: str = "weight",
    kernel_half_width: int = 3,
    donut_inner_half_width: int = 1,
    max_loci_separation: int = 10_000_000,
    n_lambda_bins: int = 40,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call chromatin dots with HiCCUPS-style kernels of given width, using cooltools.dots."""

    source = _existing(cool_path, "cool_path")
    half = _integer(kernel_half_width, "kernel_half_width", 2, 10)
    inner = _integer(donut_inner_half_width, "donut_inner_half_width", 1, half - 1)
    separation = _integer(max_loci_separation, "max_loci_separation", 1, 100_000_000)
    lambda_bins = _integer(n_lambda_bins, "n_lambda_bins", 40, 50)
    cooler, clr = _open(source)
    name = _weight(clr, weight)
    cooltools = require("cooltools", _COOLTOOLS_EXTRA)
    numutils = require("cooltools.lib.numutils", _COOLTOOLS_EXTRA)
    common = require("cooltools.lib.common", _COOLTOOLS_EXTRA)
    expected_api = require("cooltools.api.expected", _COOLTOOLS_EXTRA)
    dotfinder = require("cooltools.api.dotfinder", _COOLTOOLS_EXTRA)
    del cooler
    view = common.make_cooler_view(clr)
    try:
        expected = expected_api.expected_cis(
            clr, view_df=view, intra_only=True, smooth=False, clr_weight_name=name
        )
        kernels = {
            kind: numutils.get_kernel(half, inner, kind)
            for kind in ("donut", "vertical", "horizontal", "lowleft")
        }
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            frame = dotfinder.dots(
                clr,
                expected,
                expected_value_col="balanced.avg",
                clr_weight_name=name,
                view_df=view,
                kernels=kernels,
                max_loci_separation=separation,
                n_lambda_bins=lambda_bins,
                nproc=1,
            )
    except (ValueError, KeyError, TypeError) as exc:
        raise ConfigurationError(f"cooltools could not call dots: {exc}") from exc
    del cooltools
    columns = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "count"]
    columns += sorted(str(c) for c in frame.columns if str(c).startswith("la_exp."))
    table = (
        frame[columns].sort_values(["chrom1", "start1", "chrom2", "start2"]).reset_index(drop=True)
    )
    return {
        "status": "completed",
        "operation": "call_chromatin_dots",
        "cool_path": str(source),
        "weight": name,
        "kernel_half_width": half,
        "donut_inner_half_width": inner,
        "dot_count": len(table),
        "dots": _records(table),
        "provenance": _provenance(
            "cooltools",
            "donut, vertical, horizontal and lowleft kernels (numutils.get_kernel) with cis "
            "expected values (api.expected.expected_cis, balanced.avg) passed to dots",
            kernels="custom widths apply at any resolution; no HiCCUPS recommendation is used",
        ),
    }


def align_track_to_cooler_bins(
    track_path: str | Path,
    cool_path: str | Path,
    *,
    weight: str = "weight",
    mask_bad_bins: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align a bedGraph-like track to the bins of a cooler, using cooltools.lib.common."""

    track_source = _existing(track_path, "track_path")
    source = _existing(cool_path, "cool_path")
    mask = _flag(mask_bad_bins, "mask_bad_bins")
    _, clr = _open(source)
    name = _weight(clr, weight)
    track = _track(track_source)
    common = require("cooltools.lib.common", _COOLTOOLS_EXTRA)
    try:
        aligned = common.align_track_with_cooler(
            track, clr, view_df=None, clr_weight_name=name, mask_clr_bad_bins=mask
        )
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"The track does not match the cooler bins: {exc}") from exc
    table = aligned.sort_values(["chrom", "start", "end"]).reset_index(drop=True)
    return {
        "status": "completed",
        "operation": "align_track_to_cooler_bins",
        "track_path": str(track_source),
        "cool_path": str(source),
        "weight": name,
        "bin_count": len(table),
        "valued_bin_count": int(table["value"].notna().sum()),
        "bins": _records(table),
        "provenance": _provenance(
            "cooltools",
            "track values placed on cooler bins; bins absent from the track are None",
            masking="bins with a non-finite weight are set to None when mask_bad_bins is true",
        ),
    }


def make_track_view(track_path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Build a viewframe with one region per chromosome covered by a track, using cooltools."""

    source = _existing(track_path, "track_path")
    track = _track(source)
    common = require("cooltools.lib.common", _COOLTOOLS_EXTRA)
    try:
        view = common.view_from_track(track[["chrom", "start", "end", "value"]])
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"Could not build a viewframe from the track: {exc}") from exc
    table = view[["chrom", "start", "end", "name"]].sort_values("chrom").reset_index(drop=True)
    return {
        "status": "completed",
        "operation": "make_track_view",
        "track_path": str(source),
        "region_count": len(table),
        "regions": _records(table),
        "provenance": _provenance(
            "cooltools",
            "one region per chromosome from the smallest start and largest end "
            "(lib.common.view_from_track)",
            order="sorted by chromosome name",
        ),
    }


def mask_track_bad_bins(
    track_path: str | Path,
    cool_path: str | Path,
    *,
    weight: str = "weight",
    verbose: bool = True,
) -> dict[str, Any]:
    """Set track values to None for bins whose cooler weight is not finite, using cooltools."""

    track_source = _existing(track_path, "track_path")
    source = _existing(cool_path, "cool_path")
    _, clr = _open(source)
    name = _weight(clr, weight)
    track = _track(track_source)
    bintable = clr.bins()[:]
    common = require("cooltools.lib.common", _COOLTOOLS_EXTRA)
    try:
        masked = common.mask_cooler_bad_bins((track, "value"), (bintable, name))
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"The track could not be masked: {exc}") from exc
    table = masked.sort_values(["chrom", "start", "end"]).reset_index(drop=True)
    table["chrom"] = table["chrom"].astype(str)
    return {
        "status": "completed",
        "operation": "mask_track_bad_bins",
        "track_path": str(track_source),
        "cool_path": str(source),
        "weight": name,
        "row_count": len(table),
        "masked_count": int(table["value"].isna().sum()),
        "rows": _records(table),
        "provenance": _provenance(
            "cooltools",
            "track rows that match a cooler bin, with value set to NaN for bad bins "
            "(lib.common.mask_cooler_bad_bins)",
            matching="rows are kept only when chrom, start and end match a bin",
        ),
    }


def assign_paired_features_to_view(
    features_path: str | Path,
    view_path: str | Path,
    *,
    drop_unassigned: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Name the view region of each side of paired features, using cooltools.lib.common."""

    features_source = _existing(features_path, "features_path")
    view_source = _existing(view_path, "view_path")
    drop = _flag(drop_unassigned, "drop_unassigned")
    pandas = require("pandas", _COOLER_EXTRA)
    try:
        features = pandas.read_csv(
            features_source,
            sep="\t",
            comment="#",
            header=None,
            usecols=[0, 1, 2, 3, 4, 5],
            names=["chrom1", "start1", "end1", "chrom2", "start2", "end2"],
            dtype={
                "chrom1": str,
                "chrom2": str,
                "start1": "int64",
                "end1": "int64",
                "start2": "int64",
                "end2": "int64",
            },
        )
    except (ValueError, OSError) as exc:
        raise ConfigurationError(f"features_path must have six BEDPE columns: {exc}") from exc
    if features.empty or len(features) > _MAX_WINDOWS:
        raise ConfigurationError(f"features_path must have 1 to {_MAX_WINDOWS} rows.")
    view = _bed(view_source, "view_path")
    common = require("cooltools.lib.common", _COOLTOOLS_EXTRA)
    try:
        assigned = common.assign_view_paired(features, view, drop_unassigned=drop)
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"Could not assign the features to the view: {exc}") from exc
    table = assigned.reset_index(drop=True)
    return {
        "status": "completed",
        "operation": "assign_paired_features_to_view",
        "features_path": str(features_source),
        "view_path": str(view_source),
        "feature_count": len(table),
        "assigned_count": int(table["region1"].notna().sum()),
        "features": _records(table),
        "provenance": _provenance(
            "cooltools",
            "region names of both sides of each paired feature (lib.common.assign_view_paired)",
            naming="region names are the view names; sides outside every region are None",
        ),
    }


def assign_bins_to_regions(
    bin_ids: Sequence[int],
    regions: Sequence[Mapping[str, Any]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Name the region that contains each bin index, using cooltools.lib.common."""

    numpy = require("numpy", _COOLTOOLS_EXTRA)
    pandas = require("pandas", _COOLER_EXTRA)
    common = require("cooltools.lib.common", _COOLTOOLS_EXTRA)
    if isinstance(bin_ids, (str, bytes)) or not isinstance(bin_ids, Sequence):
        raise ConfigurationError("bin_ids must be a list of bin indices.")
    if not 1 <= len(bin_ids) <= _MAX_WINDOWS:
        raise ConfigurationError(f"bin_ids must contain 1 to {_MAX_WINDOWS} indices.")
    ids = [_integer(item, "bin id", 0, _MAX_POSITION) for item in bin_ids]
    if isinstance(regions, (str, bytes)) or not isinstance(regions, Sequence):
        raise ConfigurationError(
            "regions must be a list of mappings with name, bin_start and bin_end."
        )
    if not 1 <= len(regions) <= _MAX_WINDOWS:
        raise ConfigurationError(f"regions must contain 1 to {_MAX_WINDOWS} entries.")
    rows = []
    for item in regions:
        if not isinstance(item, Mapping):
            raise ConfigurationError(
                "Each region must be a mapping with name, bin_start and bin_end."
            )
        start = _integer(item.get("bin_start"), "bin_start", 0, _MAX_POSITION)
        end = _integer(item.get("bin_end"), "bin_end", start + 1, _MAX_POSITION)
        rows.append({"name": _text(item.get("name"), "name"), "bin_start": start, "bin_end": end})
    rows.sort(key=lambda row: (row["bin_start"], row["bin_end"]))
    for previous, current in pairwise(rows):
        if current["bin_start"] < previous["bin_end"]:
            raise ConfigurationError("Regions must not overlap.")
    frame = pandas.DataFrame(rows, columns=["name", "bin_start", "bin_end"])
    try:
        found = common.assign_regions_to_bins(numpy.asarray(ids, dtype=numpy.int64), frame)
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"Could not assign the bins to regions: {exc}") from exc
    assignments = [
        {"bin_id": bin_id, "region": None if pandas.isna(region) else str(region)}
        for bin_id, region in zip(ids, list(found), strict=True)
    ]
    return {
        "status": "completed",
        "operation": "assign_bins_to_regions",
        "bin_count": len(assignments),
        "assigned_count": sum(1 for row in assignments if row["region"] is not None),
        "assignments": assignments,
        "provenance": _provenance(
            "cooltools",
            "region whose bin span contains each bin (lib.common.assign_regions_to_bins)",
            order="regions are sorted by bin_start; overlapping regions are rejected",
        ),
    }


def dense_trans_eigenvectors(
    matrix: Sequence[Sequence[float]],
    partition: Sequence[int],
    *,
    n_eigs: int = 3,
    perc_top: float = 99.95,
    perc_bottom: float = 1.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Trans compartment eigenvectors of a dense genome-wide matrix (cooltools trans_eig)."""

    numpy = require("numpy", _COOLTOOLS_EXTRA)
    array = _symmetric(_square(matrix, "matrix"), "matrix")
    size = int(array.shape[0])
    if isinstance(partition, (str, bytes)) or not isinstance(partition, Sequence):
        raise ConfigurationError("partition must be a list of bin offsets.")
    if not 2 <= len(partition) <= _MAX_PARTITIONS:
        raise ConfigurationError("partition must contain 2 or more offsets.")
    offsets = [_integer(item, "partition offset", 0, size) for item in partition]
    if offsets[0] != 0 or offsets[-1] != size or any(b <= a for a, b in pairwise(offsets)):
        raise ConfigurationError(f"partition must increase strictly from 0 to {size}.")
    eigs = _integer(n_eigs, "n_eigs", 1, min(_MAX_EIGS, size - 1))
    top = _number(perc_top, "perc_top", 50.0, 100.0)
    bottom = _number(perc_bottom, "perc_bottom", 0.0, 50.0)
    eigdecomp = require("cooltools.api.eigdecomp", _COOLTOOLS_EXTRA)
    try:
        with _isolated_random(0, eigen=True):
            values, vectors = eigdecomp.trans_eig(
                numpy.array(array, copy=True),
                offsets,
                n_eigs=eigs,
                perc_top=top,
                perc_bottom=bottom,
                phasing_track=None,
                sort_metric=False,
            )
    except (ValueError, TypeError, KeyError) as exc:
        raise ConfigurationError(f"cooltools could not decompose the matrix: {exc}") from exc
    return {
        "status": "completed",
        "operation": "dense_trans_eigenvectors",
        "size": size,
        "partition": offsets,
        "n_eigs": eigs,
        "perc_top": top,
        "perc_bottom": bottom,
        "eigenvalues": _vector(values),
        "eigenvectors": _matrix(_canonical_rows(vectors)),
        "provenance": _provenance(
            "cooltools",
            "trans eigenvectors of the O/E matrix after trans clipping "
            "(cooltools.api.eigdecomp.trans_eig)",
            sign="each eigenvector's largest absolute value is positive",
            random_state="numpy global state and the C random() state are restored after the call",
            limitation="ARPACK uses a fixed start vector, so its Fortran seed is not advanced",
        ),
    }


def pairwise_diagonal_sums(
    cool_path: str | Path,
    *,
    weight: str = "weight",
    ignore_diags: int = 2,
    verbose: bool = True,
) -> dict[str, Any]:
    """Diagonal sums of intra-chromosomal blocks for pairs of chromosome arms (cooltools)."""

    source = _existing(cool_path, "cool_path")
    diags = _integer(ignore_diags, "ignore_diags", 0, 10)
    _, clr = _open(source)
    name = _weight(clr, weight)
    common = require("cooltools.lib.common", _COOLTOOLS_EXTRA)
    expected_api = require("cooltools.api.expected", _COOLTOOLS_EXTRA)
    view = common.make_cooler_view(clr)
    try:
        frame = expected_api.diagsum_pairwise(clr, view, clr_weight_name=name, ignore_diags=diags)
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"cooltools could not sum the diagonals: {exc}") from exc
    table = frame.sort_values(["region1", "region2", "dist"]).reset_index(drop=True)
    return {
        "status": "completed",
        "operation": "pairwise_diagonal_sums",
        "cool_path": str(source),
        "weight": name,
        "ignore_diags": diags,
        "row_count": len(table),
        "rows": _records(table),
        "provenance": _provenance(
            "cooltools",
            "per-diagonal pixel counts of each intra-chromosomal block "
            "(api.expected.diagsum_pairwise)",
            order="sorted by region1, region2 and distance",
        ),
    }


def pairwise_block_sums(
    cool_path: str | Path,
    *,
    weight: str = "weight",
    verbose: bool = True,
) -> dict[str, Any]:
    """Sums over rectangular blocks of every pair of chromosomes, using cooltools."""

    source = _existing(cool_path, "cool_path")
    _, clr = _open(source)
    name = _weight(clr, weight)
    common = require("cooltools.lib.common", _COOLTOOLS_EXTRA)
    expected_api = require("cooltools.api.expected", _COOLTOOLS_EXTRA)
    view = common.make_cooler_view(clr)
    try:
        frame = expected_api.blocksum_pairwise(clr, view, clr_weight_name=name)
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"cooltools could not sum the blocks: {exc}") from exc
    table = frame.sort_values(["region1", "region2"]).reset_index(drop=True)
    return {
        "status": "completed",
        "operation": "pairwise_block_sums",
        "cool_path": str(source),
        "weight": name,
        "row_count": len(table),
        "rows": _records(table),
        "provenance": _provenance(
            "cooltools",
            "pixels summed over each block of a region pair (api.expected.blocksum_pairwise)",
            order="sorted by region1 and region2",
        ),
    }


def saddle_compartment_strength(
    saddle_sums: Sequence[Sequence[float | None]],
    saddle_counts: Sequence[Sequence[float | None]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compartment strength from a saddle sum and count matrix, using cooltools.saddle_strength."""

    saddle = _square(saddle_sums, "saddle_sums", allow_nan=True)
    counts = _square(saddle_counts, "saddle_counts", allow_nan=True)
    if saddle.shape != counts.shape:
        raise ConfigurationError("saddle_sums and saddle_counts must have the same shape.")
    saddle_api = require("cooltools.api.saddle", _COOLTOOLS_EXTRA)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ratios = saddle_api.saddle_strength(saddle, counts)
    rows = [{"k": int(k), "strength": _plain(float(ratios[k]))} for k in range(1, len(ratios))]
    return {
        "status": "completed",
        "operation": "saddle_compartment_strength",
        "size": int(saddle.shape[0]),
        "strengths": rows,
        "provenance": _provenance(
            "cooltools",
            "ratio of intra-compartment to inter-compartment mean as the corners grow "
            "(api.saddle.saddle_strength)",
            k="k is the number of bins taken from each corner of the saddle",
        ),
    }


def expand_features_to_bins(
    features_path: str | Path,
    *,
    flank: int,
    resolution: int,
    format: str = "bed",
    verbose: bool = True,
) -> dict[str, Any]:
    """Expand BED or BEDPE features by a flank and snap them to bin edges, using cooltools."""

    source = _existing(features_path, "features_path")
    flank_bp = _integer(flank, "flank", 0, _MAX_POSITION)
    size = _integer(resolution, "resolution", 1, 10_000_000)
    if format not in ("bed", "bedpe"):
        raise ConfigurationError("format must be 'bed' or 'bedpe'.")
    pandas = require("pandas", _COOLER_EXTRA)
    snipping = require("cooltools.api.snipping", _COOLTOOLS_EXTRA)
    if format == "bed":
        columns = ["chrom", "start", "end"]
        dtypes = {"chrom": str, "start": "int64", "end": "int64"}
    else:
        columns = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
        dtypes = {
            "chrom1": str,
            "chrom2": str,
            "start1": "int64",
            "end1": "int64",
            "start2": "int64",
            "end2": "int64",
        }
    try:
        features = pandas.read_csv(
            source,
            sep="\t",
            comment="#",
            header=None,
            usecols=list(range(len(columns))),
            names=columns,
            dtype=dtypes,
        )
    except (ValueError, OSError) as exc:
        raise ConfigurationError(f"features_path must have {len(columns)} columns: {exc}") from exc
    if features.empty or len(features) > _MAX_WINDOWS:
        raise ConfigurationError(f"features_path must have 1 to {_MAX_WINDOWS} rows.")
    try:
        expanded = snipping.expand_align_features(features, flank_bp, size, format=format)
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"Could not expand the features: {exc}") from exc
    table = expanded.reset_index(drop=True)
    return {
        "status": "completed",
        "operation": "expand_features_to_bins",
        "features_path": str(source),
        "format": format,
        "flank": flank_bp,
        "resolution": size,
        "feature_count": len(table),
        "features": _records(table),
        "provenance": _provenance(
            "cooltools",
            "feature centres expanded by the flank and snapped to multiples of the resolution "
            "(api.snipping.expand_align_features)",
            order="input order is kept",
        ),
    }


def adaptive_coarsegrain_matrix(
    matrix: Sequence[Sequence[float]],
    counts: Sequence[Sequence[float]],
    *,
    cutoff: float = 5.0,
    max_levels: int = 8,
    min_shape: int = 8,
    verbose: bool = True,
) -> dict[str, Any]:
    """Adaptively coarse-grain a contact matrix by pooling sparse neighbourhoods (cooltools)."""

    array = _square(matrix, "matrix")
    raw = _square(counts, "counts")
    if array.shape != raw.shape:
        raise ConfigurationError("matrix and counts must have the same shape.")
    limit = _number(cutoff, "cutoff", 0.0, 1_000_000.0)
    levels = _integer(max_levels, "max_levels", 1, 8)
    shape = _integer(min_shape, "min_shape", 2, 64)
    numutils = require("cooltools.lib.numutils", _COOLTOOLS_EXTRA)
    try:
        coarse = numutils.adaptive_coarsegrain(
            array, raw, cutoff=limit, max_levels=levels, min_shape=shape
        )
    except (ValueError, TypeError, IndexError) as exc:
        raise ConfigurationError(f"cooltools could not coarse-grain the matrix: {exc}") from exc
    return {
        "status": "completed",
        "operation": "adaptive_coarsegrain_matrix",
        "size": int(array.shape[0]),
        "cutoff": limit,
        "max_levels": levels,
        "min_shape": shape,
        "matrix": _matrix(coarse),
        "provenance": _provenance(
            "cooltools",
            "2x2 pooling of neighbourhoods until each pixel has at least cutoff raw counts "
            "(lib.numutils.adaptive_coarsegrain)",
        ),
    }


def smooth_contact_decay(
    cool_path: str | Path,
    *,
    mode: str = "genomewide",
    weight: str = "weight",
    sigma_log10: float = 0.1,
    window_sigma: int = 5,
    points_per_sigma: int = 10,
    verbose: bool = True,
) -> dict[str, Any]:
    """Smooth the balanced contact-versus-distance curve in log space, using cooltools."""

    source = _existing(cool_path, "cool_path")
    if mode not in ("genomewide", "per_region"):
        raise ConfigurationError("mode must be 'genomewide' or 'per_region'.")
    sigma = _number(sigma_log10, "sigma_log10", 0.001, 2.0)
    window = _integer(window_sigma, "window_sigma", 1, 50)
    per_sigma = _integer(points_per_sigma, "points_per_sigma", 1, 100)
    _, clr = _open(source)
    name = _weight(clr, weight)
    common = require("cooltools.lib.common", _COOLTOOLS_EXTRA)
    expected_api = require("cooltools.api.expected", _COOLTOOLS_EXTRA)
    view = common.make_cooler_view(clr)
    try:
        cvd = expected_api.expected_cis(
            clr, view_df=view, intra_only=True, smooth=False, clr_weight_name=name
        )
        cols = expected_api.DEFAULT_BALANCED_CVD_COLS
        smoother = (
            expected_api.genomewide_smooth_cvd
            if mode == "genomewide"
            else expected_api.per_region_smooth_cvd
        )
        smoothed = smoother(
            cvd,
            sigma_log10=sigma,
            window_sigma=window,
            points_per_sigma=per_sigma,
            cols=cols,
            suffix=".smoothed",
        )
    except (ValueError, KeyError, TypeError) as exc:
        raise ConfigurationError(f"cooltools could not smooth the decay curve: {exc}") from exc
    keep = ["region1", "region2", "dist", "n_valid", "balanced.sum", "balanced.avg"]
    keep += sorted(str(c) for c in smoothed.columns if "smoothed" in str(c))
    table = smoothed[keep].sort_values(["region1", "region2", "dist"]).reset_index(drop=True)
    return {
        "status": "completed",
        "operation": "smooth_contact_decay",
        "cool_path": str(source),
        "mode": mode,
        "weight": name,
        "sigma_log10": sigma,
        "window_sigma": window,
        "points_per_sigma": per_sigma,
        "row_count": len(table),
        "rows": _records(table),
        "provenance": _provenance(
            "cooltools",
            "Gaussian smoothing of balanced.avg over log10(distance) "
            "(api.expected.genomewide_smooth_cvd or per_region_smooth_cvd)",
            regions="genomewide pools all regions; per_region smooths each region separately",
        ),
    }


def replace_cis_with_decoys(
    matrix: Sequence[Sequence[float]],
    mask: Sequence[Sequence[int]],
    *,
    seed: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Replace masked cis pixels of a contact matrix with decoy contacts (cooltools fake_cis)."""

    numpy = require("numpy", _COOLTOOLS_EXTRA)
    array = _square(matrix, "matrix")
    size = int(array.shape[0])
    flags = _square(mask, "mask")
    if flags.shape != array.shape or not numpy.isin(flags, (0.0, 1.0)).all():
        raise ConfigurationError("mask must be a 0/1 matrix with the same shape as matrix.")
    value = _integer(seed, "seed", 0, _MAX_SEED)
    numutils = require("cooltools.lib.numutils", _COOLTOOLS_EXTRA)
    work = numpy.array(array, dtype=float, copy=True)
    cis = numpy.ascontiguousarray(flags, dtype=numpy.uint8)
    try:
        with _isolated_random(value, eigen=False):
            numutils.fake_cis(work, cis)
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"cooltools could not write the decoys: {exc}") from exc
    replaced = int(numpy.count_nonzero(cis))
    return {
        "status": "completed",
        "operation": "replace_cis_with_decoys",
        "size": size,
        "seed": value,
        "replaced_pixels": replaced,
        "matrix": _matrix(work),
        "provenance": _provenance(
            "cooltools",
            "masked pixels are filled with decoy values taken from the C random() stream "
            "(lib.numutils.fake_cis)",
            random_state="C random() is seeded for this call and restored afterwards",
        ),
    }


def expected_full_table(
    cool_path: str | Path,
    *,
    weight: str = "weight",
    fast: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Expected contacts for every pair of chromosome regions, using the cooltools sandbox."""

    source = _existing(cool_path, "cool_path")
    flag = _flag(fast, "fast")
    _, clr = _open(source)
    name = _weight(clr, weight)
    common = require("cooltools.lib.common", _COOLTOOLS_EXTRA)
    sandbox = require("cooltools.sandbox.obs_over_exp_cooler", _COOLTOOLS_EXTRA)
    view = common.make_cooler_view(clr)
    function = sandbox.expected_full_fast if flag else sandbox.expected_full
    try:
        frame = function(clr, view_df=view, clr_weight_name=name, nproc=1)
    except (ValueError, KeyError, TypeError) as exc:
        raise ConfigurationError(f"cooltools could not compute the expected table: {exc}") from exc
    keep = ["region1", "region2", "dist", "n_valid", "count.sum", "expected"]
    table = frame[keep].sort_values(["region1", "region2", "dist"]).reset_index(drop=True)
    return {
        "status": "completed",
        "operation": "expected_full_table",
        "cool_path": str(source),
        "weight": name,
        "fast": flag,
        "row_count": len(table),
        "rows": _records(table),
        "provenance": _provenance(
            "cooltools",
            "distance-decay expected for cis regions and block means for trans regions "
            "(sandbox.obs_over_exp_cooler.expected_full)"
            if not flag
            else "faster sandbox variant of the same expected table "
            "(sandbox.obs_over_exp_cooler.expected_full_fast)",
            order="sorted by region1, region2 and distance",
        ),
    }


def observed_over_expected_pixels(
    cool_path: str | Path,
    *,
    weight: str = "weight",
    verbose: bool = True,
) -> dict[str, Any]:
    """Observed over expected for every pixel of a cooler, using the cooltools sandbox."""

    source = _existing(cool_path, "cool_path")
    _, clr = _open(source)
    name = _weight(clr, weight)
    common = require("cooltools.lib.common", _COOLTOOLS_EXTRA)
    sandbox = require("cooltools.sandbox.obs_over_exp_cooler", _COOLTOOLS_EXTRA)
    pandas = require("pandas", _COOLER_EXTRA)
    view = common.make_cooler_view(clr)
    try:
        expected = sandbox.expected_full(clr, view_df=view, clr_weight_name=name, nproc=1)
        bins = clr.bins()[:]
        bins["r"] = _bin_regions(bins, view)
        pixels = clr.pixels()[:]
        frame = sandbox.obs_over_exp(
            pixels,
            bins,
            expected,
            view_column_name="r",
            clr_weight_name=name,
            oe_column_name="oe",
        )
    except (ValueError, KeyError, TypeError) as exc:
        raise ConfigurationError(
            f"cooltools could not compute observed over expected: {exc}"
        ) from exc
    keep = ["bin1_id", "bin2_id", "count", "balanced", "r1", "r2", "dist", "expected", "oe"]
    table = pandas.DataFrame(frame)[keep].sort_values(["bin1_id", "bin2_id"]).reset_index(drop=True)
    return {
        "status": "completed",
        "operation": "observed_over_expected_pixels",
        "cool_path": str(source),
        "weight": name,
        "pixel_count": len(table),
        "pixels": _records(table),
        "provenance": _provenance(
            "cooltools",
            "balanced pixel counts divided by the expected value of their region pair and "
            "distance (sandbox.obs_over_exp_cooler.obs_over_exp)",
            order="sorted by bin1_id and bin2_id",
        ),
    }


__all__ = [
    "adaptive_coarsegrain_matrix",
    "align_track_to_cooler_bins",
    "assign_bins_to_regions",
    "assign_paired_features_to_view",
    "call_chromatin_dots",
    "dense_trans_eigenvectors",
    "expand_features_to_bins",
    "expected_full_table",
    "make_track_view",
    "mask_track_bad_bins",
    "observed_over_expected_pixels",
    "pairwise_block_sums",
    "pairwise_diagonal_sums",
    "replace_cis_with_decoys",
    "saddle_compartment_strength",
    "smooth_contact_decay",
]

_register_public_verbose(globals(), tuple(__all__))
