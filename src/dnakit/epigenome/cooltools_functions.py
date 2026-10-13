"""Contact-matrix, compartment, view and dense-matrix functions from cooltools.

Each public function imports cooltools and cooler only when it is called. Missing
dependencies raise ``BackendUnavailableError`` with the ``pip install "dnakit[...]"``
hint. Weighted functions need a balancing column (default ``weight``); run
``balance_contact_matrix()`` first. Inputs are checked before upstream calls, and
upstream ``ValueError``, ``KeyError`` and ``OSError`` are reported as
``ConfigurationError``.

Functions that write a cooler file create a new path and never overwrite one.
Results are plain dictionaries; NaN values are returned as ``None``.
"""

from __future__ import annotations

import contextlib
import ctypes
import math
import warnings
from collections.abc import Sequence
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_COOLER_EXTRA = "bio-cooler"
_COOLTOOLS_EXTRA = "bio-cooltools"
_NUMPY_EXTRA = "bio-analysis"
_MAX_BINS = 50_000
_MAX_EIGS = 10
_MAX_DENSE = 2_000
_MAX_TRACK_ROWS = 1_000_000
_MAX_ITER = 100_000
_MAX_WINDOWS = 10_000
_MAX_VIEW_ROWS = 10_000
_MAX_PIXEL_SAMPLE = 1_000_000_000_000
_MAX_RATIO = 2.0


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


def _new(value: object, name: str) -> Path:
    """Output path that must not exist yet; its parent directory must exist."""

    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a file path.")
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} already exists and will not be overwritten: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The directory for {name} does not exist: {path.parent}")
    return path


def _plain(value: Any) -> Any:
    """Convert numpy and pandas scalars to Python values and non-finite floats to None."""

    if value is None:
        return None
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


def _seed_c_random(value: int) -> None:
    """Seed the C library random() that cooltools' compiled fake-cis step uses (glibc, BSD)."""

    try:
        ctypes.CDLL(None).srandom(value)
    except (OSError, AttributeError):
        return


@contextlib.contextmanager
def _reproducible_eigen() -> Any:
    """Fix the random inputs of cooltools' eigen-solvers for one call.

    Two sources change between calls: the ARPACK start vector used by scipy's eigsh and
    eigs, and the C library random() used by the compiled fake-cis step of the trans
    compartments. The start vector is replaced by a fixed one and the C generator is
    seeded with 0. The scipy functions are restored on exit; the C generator state is
    not restored.
    """

    numpy = require("numpy", _NUMPY_EXTRA)
    sparse_linalg = require("scipy.sparse.linalg", _NUMPY_EXTRA)
    originals = {name: getattr(sparse_linalg, name) for name in ("eigsh", "eigs")}

    def wrap(function: Any) -> Any:
        def seeded(matrix: Any, *args: Any, **kwargs: Any) -> Any:
            if "v0" not in kwargs and len(args) <= 3:
                size = matrix.shape[0]
                kwargs["v0"] = numpy.random.default_rng(0).random(size) + 0.1
            return function(matrix, *args, **kwargs)

        return seeded

    for name, function in originals.items():
        setattr(sparse_linalg, name, wrap(function))
    _seed_c_random(0)
    try:
        yield
    finally:
        for name, function in originals.items():
            setattr(sparse_linalg, name, function)


def _canonical_rows(rows: Any) -> Any:
    """Flip each eigenvector so its largest finite entry is positive.

    ARPACK (scipy.sparse.linalg.eigsh) keeps state between calls, so the sign of an
    eigenvector can differ from one process to another. This convention makes it stable.
    """

    numpy = require("numpy", _NUMPY_EXTRA)
    array = numpy.array(rows, dtype=float, copy=True)
    for index, row in enumerate(array):
        finite = numpy.isfinite(row)
        if finite.any():
            pivot = int(numpy.argmax(numpy.where(finite, numpy.abs(row), -1.0)))
            if row[pivot] < 0:
                array[index] = -row
    return array


def _canonical_frame(frame: Any) -> Any:
    """Apply the eigenvector sign convention to the E1, E2, ... columns of a table."""

    columns = [
        column
        for column in frame.columns
        if str(column).startswith("E") and str(column)[1:].isdigit()
    ]
    if columns:
        frame = frame.copy()
        frame[columns] = _canonical_rows(frame[columns].to_numpy().T).T
    return frame


def _provenance(distribution: str, method: str, **notes: str) -> dict[str, Any]:
    try:
        installed = version(distribution)
    except PackageNotFoundError:
        installed = "unknown"
    return {"package": f"{distribution} {installed}", "method": method, **notes}


def _modules(*names: str) -> tuple[Any, ...]:
    return tuple(
        require(name, _COOLTOOLS_EXTRA if name.startswith("cooltools") else _COOLER_EXTRA)
        for name in names
    )


def _open(cool_path: Path) -> tuple[Any, Any]:
    (cooler,) = _modules("cooler")
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


def _binsize(clr: Any) -> int:
    return int(clr.binsize)


def _summary(clr: Any) -> dict[str, Any]:
    info = clr.info
    return {
        "binsize": _binsize(clr),
        "bins": int(info["nbins"]),
        "nonzero_pixels": int(info["nnz"]),
        "total_contacts": int(info["sum"]),
        "chromosomes": [str(name) for name in clr.chromnames],
    }


def _track(path: Path) -> Any:
    pandas = require("pandas", _COOLER_EXTRA)
    try:
        table = pandas.read_csv(
            path,
            sep="\t",
            comment="#",
            header=None,
            names=["chrom", "start", "end", "value"],
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


def _square(values: Any, name: str) -> Any:
    numpy = require("numpy", _NUMPY_EXTRA)
    try:
        array = numpy.asarray(values, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"{name} must be a square numeric matrix.") from exc
    if array.ndim != 2 or array.shape[0] != array.shape[1]:
        raise ConfigurationError(f"{name} must be a square matrix.")
    if not 2 <= array.shape[0] <= _MAX_DENSE:
        raise ConfigurationError(f"{name} must be 2 to {_MAX_DENSE} rows wide.")
    if not numpy.isfinite(array).all():
        raise ConfigurationError(f"{name} must contain only finite values.")
    return array


def _symmetric(array: Any, name: str) -> Any:
    numpy = require("numpy", _NUMPY_EXTRA)
    if not numpy.allclose(array, array.T, rtol=0, atol=1e-9):
        raise ConfigurationError(f"{name} must be symmetric.")
    return array


def calculate_expected_trans(
    cool_path: str | Path,
    *,
    weight: str = "weight",
    verbose: bool = True,
) -> dict[str, Any]:
    """Average inter-chromosomal contact frequency for each pair of chromosomes, using cooltools."""

    source = _existing(cool_path, "cool_path")
    (cooltools,) = _modules("cooltools")
    _, clr = _open(source)
    name = _weight(clr, weight)
    try:
        table = cooltools.expected_trans(clr, clr_weight_name=name, nproc=1)
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"cooltools rejected the trans expected request: {exc}") from exc
    rows = _records(table.sort_values(["region1", "region2"]))
    return {
        "status": "completed",
        "operation": "calculate_expected_trans",
        "cool_path": str(source),
        "weight": name,
        "pair_count": len(rows),
        "pairs": rows,
        "provenance": _provenance(
            "cooltools",
            "mean balanced and raw contacts for each chromosome pair (cooltools.expected_trans)",
            order="rows are sorted by region1, region2",
        ),
    }


def calculate_compartments(
    cool_path: str | Path,
    *,
    n_eigs: int = 3,
    weight: str = "weight",
    verbose: bool = True,
) -> dict[str, Any]:
    """Cis compartment eigenvalues and eigenvectors for each chromosome, using cooltools."""

    source = _existing(cool_path, "cool_path")
    eigs = _integer(n_eigs, "n_eigs", 1, _MAX_EIGS)
    cooltools, common = _modules("cooltools", "cooltools.lib.common")
    _, clr = _open(source)
    name = _weight(clr, weight)
    view = common.make_cooler_view(clr)
    try:
        with _reproducible_eigen():
            values, vectors = cooltools.eigs_cis(
                clr, view_df=view, n_eigs=eigs, clr_weight_name=name
            )
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"cooltools could not compute compartments: {exc}") from exc
    return {
        "status": "completed",
        "operation": "calculate_compartments",
        "cool_path": str(source),
        "weight": name,
        "n_eigs": eigs,
        "eigenvalues": _records(values.sort_values(["chrom", "start"])),
        "eigenvectors": _records(_canonical_frame(vectors.sort_values(["chrom", "start"]))),
        "provenance": _provenance(
            "cooltools",
            "O/E correlation matrix eigenvectors per chromosome (cooltools.eigs_cis)",
            sign="each eigenvector's largest absolute value is positive",
        ),
    }


def calculate_trans_compartments(
    cool_path: str | Path,
    *,
    n_eigs: int = 3,
    weight: str = "weight",
    verbose: bool = True,
) -> dict[str, Any]:
    """Genome-wide trans compartment eigenvectors from inter-chromosomal contacts (cooltools)."""

    source = _existing(cool_path, "cool_path")
    eigs = _integer(n_eigs, "n_eigs", 1, _MAX_EIGS)
    (cooltools,) = _modules("cooltools")
    _, clr = _open(source)
    name = _weight(clr, weight)
    try:
        with _reproducible_eigen():
            values, vectors = cooltools.eigs_trans(clr, n_eigs=eigs, clr_weight_name=name)
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"cooltools could not compute trans compartments: {exc}") from exc
    return {
        "status": "completed",
        "operation": "calculate_trans_compartments",
        "cool_path": str(source),
        "weight": name,
        "n_eigs": eigs,
        "eigenvalues": _records(values),
        "eigenvectors": _records(_canonical_frame(vectors)),
        "provenance": _provenance(
            "cooltools",
            "eigenvectors of the inter-chromosomal correlation matrix (cooltools.eigs_trans)",
            sign="each eigenvector's largest absolute value is positive",
            start="ARPACK start vector and C random seed are fixed, so repeated calls agree",
        ),
    }


def calculate_saddle(
    cool_path: str | Path,
    track_path: str | Path,
    *,
    contact_type: str = "cis",
    n_bins: int = 5,
    qrange: Sequence[float] = (0.02, 0.98),
    min_diag: int = 3,
    weight: str = "weight",
    verbose: bool = True,
) -> dict[str, Any]:
    """Saddle plot: average contacts between genomic bins grouped by a track, using cooltools."""

    source = _existing(cool_path, "cool_path")
    track_file = _existing(track_path, "track_path")
    if contact_type not in ("cis", "trans"):
        raise ConfigurationError("contact_type must be 'cis' or 'trans'.")
    bins = _integer(n_bins, "n_bins", 1, 100)
    if isinstance(qrange, (str, bytes)) or not isinstance(qrange, Sequence) or len(qrange) != 2:
        raise ConfigurationError("qrange must be a pair of quantiles (low, high).")
    low = _number(qrange[0], "qrange low", 0.0, 1.0)
    high = _number(qrange[1], "qrange high", 0.0, 1.0)
    if not low < high:
        raise ConfigurationError("qrange low must be smaller than qrange high.")
    diag = _integer(min_diag, "min_diag", 0, 1000)
    cooltools, common = _modules("cooltools", "cooltools.lib.common")
    _, clr = _open(source)
    name = _weight(clr, weight)
    track = _track(track_file)
    view = common.make_cooler_view(clr)
    try:
        if contact_type == "cis":
            expected = cooltools.expected_cis(clr, view_df=view, clr_weight_name=name, nproc=1)
        else:
            expected = cooltools.expected_trans(clr, view_df=view, clr_weight_name=name)
        interaction_sum, interaction_count = cooltools.saddle(
            clr,
            expected,
            track,
            contact_type,
            bins,
            qrange=(low, high),
            view_df=view,
            clr_weight_name=name,
            expected_value_col="balanced.avg",
            min_diag=diag,
        )
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"cooltools could not build the saddle plot: {exc}") from exc
    return {
        "status": "completed",
        "operation": "calculate_saddle",
        "cool_path": str(source),
        "track_path": str(track_file),
        "contact_type": contact_type,
        "n_bins": bins,
        "qrange": [low, high],
        "matrix_size": int(interaction_sum.shape[0]),
        "interaction_sum": _matrix(interaction_sum),
        "interaction_count": _matrix(interaction_count),
        "provenance": _provenance(
            "cooltools",
            "mean balanced contacts per pair of track quantile bins (cooltools.saddle)",
            layout="the outer rows and columns hold the bins below and above qrange",
        ),
    }


def sample_contacts(
    cool_path: str | Path,
    output_path: str | Path,
    *,
    count: int | None = None,
    fraction: float | None = None,
    seed: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a cooler that keeps a random subset of contacts, using cooltools (seeded)."""

    source = _existing(cool_path, "cool_path")
    target = _new(output_path, "output_path")
    if (count is None) == (fraction is None):
        raise ConfigurationError("Give exactly one of count (exact) or fraction (binomial).")
    seed_value = _integer(seed, "seed", 0, 2**31 - 1)
    numpy = require("numpy", _NUMPY_EXTRA)
    (cooltools,) = _modules("cooltools")
    _, clr = _open(source)
    total = int(clr.info["sum"])
    if count is not None:
        contacts = _integer(count, "count", 1, _MAX_PIXEL_SAMPLE)
        if contacts > total:
            raise ConfigurationError(f"count {contacts} is above the {total} contacts available.")
        kwargs: dict[str, Any] = {"count": contacts, "exact": True}
        method = "exact"
    else:
        share = _number(fraction, "fraction", 0.0, 1.0)
        if share <= 0.0:
            raise ConfigurationError("fraction must be greater than 0 and at most 1.")
        kwargs = {"frac": share, "exact": False}
        method = "binomial"
    state = numpy.random.get_state()
    numpy.random.seed(seed_value)
    try:
        cooltools.sample(clr, str(target), nproc=1, **kwargs)
    except (ValueError, KeyError, OSError) as exc:
        target.unlink(missing_ok=True)
        raise ConfigurationError(f"cooltools could not sample this cooler: {exc}") from exc
    finally:
        numpy.random.set_state(state)
    (cooler,) = _modules("cooler")
    after = int(cooler.Cooler(str(target)).info["sum"])
    return {
        "status": "completed",
        "operation": "sample_contacts",
        "input": str(source),
        "output": str(target),
        "method": method,
        "seed": seed_value,
        "total_contacts_before": total,
        "total_contacts_after": after,
        "provenance": _provenance(
            "cooltools",
            "random subset of pixel counts (cooltools.sample); seeded with numpy global RNG "
            "and restored afterwards",
            note="exact keeps the requested count; binomial keeps about the fraction",
        ),
    }


def calculate_virtual_4c(
    cool_path: str | Path,
    viewpoint: str,
    *,
    weight: str = "weight",
    verbose: bool = True,
) -> dict[str, Any]:
    """Contact profile of one viewpoint across the genome (virtual 4C), using cooltools."""

    source = _existing(cool_path, "cool_path")
    point = _text(viewpoint, "viewpoint")
    cooltools, cooler_util = _modules("cooltools", "cooler.util")
    _, clr = _open(source)
    name = _weight(clr, weight)
    try:
        region = cooler_util.parse_region(point, clr.chromsizes)
        table = cooltools.virtual4c(clr, tuple(region), clr_weight_name=name, nproc=1)
    except (ValueError, KeyError, TypeError) as exc:
        raise ConfigurationError(
            f"Viewpoint is not valid for this cooler: {point} ({exc})"
        ) from exc
    rows = _records(table)
    return {
        "status": "completed",
        "operation": "calculate_virtual_4c",
        "cool_path": str(source),
        "viewpoint": {"chrom": str(region[0]), "start": int(region[1]), "end": int(region[2])},
        "weight": name,
        "bin_count": len(rows),
        "profile": rows,
        "provenance": _provenance(
            "cooltools",
            "balanced contacts of the viewpoint bins with every bin (cooltools.virtual4c)",
            status="experimental upstream; interface may change",
        ),
    }


def calculate_directionality_index(
    cool_path: str | Path,
    *,
    window_bp: int = 100_000,
    weight: str = "weight",
    verbose: bool = True,
) -> dict[str, Any]:
    """Directionality index and insulation-based boundaries per bin, using cooltools."""

    source = _existing(cool_path, "cool_path")
    window = _integer(window_bp, "window_bp", 1, 10_000_000)
    (cooltools,) = _modules("cooltools")
    _, clr = _open(source)
    name = _weight(clr, weight)
    if window < _binsize(clr):
        raise ConfigurationError("window_bp must be at least the cooler bin size.")
    try:
        table = cooltools.directionality(clr, window_bp=window, balance=name)
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"cooltools could not compute directionality: {exc}") from exc
    rows = _records(table)
    return {
        "status": "completed",
        "operation": "calculate_directionality_index",
        "cool_path": str(source),
        "window_bp": window,
        "weight": name,
        "bin_count": len(rows),
        "bins": rows,
        "provenance": _provenance(
            "cooltools",
            "diamond window directionality index on balanced contacts (cooltools.directionality)",
            order="rows follow the bin table order",
        ),
    }


def pileup_features(
    cool_path: str | Path,
    features_path: str | Path,
    *,
    flank_bp: int = 100_000,
    weight: str = "weight",
    verbose: bool = True,
) -> dict[str, Any]:
    """Average contact snippet around BED features (pileup), using cooltools."""

    source = _existing(cool_path, "cool_path")
    features_file = _existing(features_path, "features_path")
    flank = _integer(flank_bp, "flank_bp", 1, 10_000_000)
    numpy = require("numpy", _NUMPY_EXTRA)
    pandas = require("pandas", _COOLER_EXTRA)
    cooltools, common = _modules("cooltools", "cooltools.lib.common")
    _, clr = _open(source)
    name = _weight(clr, weight)
    if flank < _binsize(clr):
        raise ConfigurationError("flank_bp must be at least the cooler bin size.")
    try:
        features = pandas.read_csv(
            features_file,
            sep="\t",
            comment="#",
            header=None,
            usecols=[0, 1, 2],
            names=["chrom", "start", "end"],
            dtype={"chrom": str, "start": "int64", "end": "int64"},
        )
    except (ValueError, OSError) as exc:
        raise ConfigurationError(f"Features need chrom, start and end columns: {exc}") from exc
    if features.empty:
        raise ConfigurationError("The features file has no intervals.")
    view = common.make_cooler_view(clr)
    try:
        stack = cooltools.pileup(
            clr,
            features,
            view_df=view,
            expected_df=None,
            flank=flank,
            clr_weight_name=name,
            nproc=1,
        )
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"cooltools could not build the pileup: {exc}") from exc
    with warnings.catch_warnings(), numpy.errstate(all="ignore"):
        warnings.simplefilter("ignore", RuntimeWarning)
        mean = numpy.nanmean(stack, axis=0)
    with_data = int(numpy.isfinite(stack).any(axis=(1, 2)).sum())
    return {
        "status": "completed",
        "operation": "pileup_features",
        "cool_path": str(source),
        "features_path": str(features_file),
        "flank_bp": flank,
        "weight": name,
        "feature_count": int(stack.shape[0]),
        "features_with_data": with_data,
        "snippet_size": int(stack.shape[1]),
        "mean_snippet": _matrix(mean),
        "provenance": _provenance(
            "cooltools",
            "nan-mean of per-feature balanced snippets (cooltools.pileup)",
            note="expected normalisation is not applied; features are not filtered by score",
        ),
    }


def calculate_logbin_expected(
    cool_path: str | Path,
    *,
    bins_per_order_magnitude: int = 10,
    min_nvalid: int = 200,
    min_count: int = 50,
    weight: str = "weight",
    verbose: bool = True,
) -> dict[str, Any]:
    """Expected contact frequency vs distance, log-binned with local slopes, using cooltools."""

    source = _existing(cool_path, "cool_path")
    per_order = _integer(bins_per_order_magnitude, "bins_per_order_magnitude", 1, 50)
    valid_min = _integer(min_nvalid, "min_nvalid", 0, 1_000_000_000)
    count_min = _integer(min_count, "min_count", 0, 1_000_000_000)
    cooltools, common, expected_api = _modules(
        "cooltools", "cooltools.lib.common", "cooltools.api.expected"
    )
    _, clr = _open(source)
    name = _weight(clr, weight)
    view = common.make_cooler_view(clr)
    try:
        expected = cooltools.expected_cis(clr, view_df=view, clr_weight_name=name, nproc=1)
        binned, slopes, edges = expected_api.logbin_expected(
            expected,
            summary_name="balanced.sum",
            bins_per_order_magnitude=per_order,
            min_nvalid=valid_min,
            min_count=count_min,
        )
    except (ValueError, KeyError, AttributeError) as exc:
        raise ConfigurationError(f"cooltools could not log-bin expected: {exc}") from exc
    return {
        "status": "completed",
        "operation": "calculate_logbin_expected",
        "cool_path": str(source),
        "weight": name,
        "bins_per_order_magnitude": per_order,
        "binned_rows": len(binned),
        "binned": _records(binned.sort_values(["region1", "region2", "dist_bin_id"])),
        "slopes": _records(slopes.sort_values(["region1", "region2", "dist_bin_id"])),
        "bin_edges": [int(edge) for edge in edges],
        "provenance": _provenance(
            "cooltools",
            "log-spaced distance bins of expected contacts with local slopes "
            "(cooltools.api.expected.logbin_expected)",
        ),
    }


def rearrange_contact_matrix(
    cool_path: str | Path,
    view_path: str | Path,
    output_path: str | Path,
    *,
    new_chrom_column: str = "new_chrom",
    orientation_column: str = "strand",
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a reordered cooler that follows a viewframe (regions, names, strands)."""

    source = _existing(cool_path, "cool_path")
    view_file = _existing(view_path, "view_path")
    target = _new(output_path, "output_path")
    new_chrom = _text(new_chrom_column, "new_chrom_column")
    orientation = _text(orientation_column, "orientation_column")
    pandas = require("pandas", _COOLER_EXTRA)
    (rearrange,) = _modules("cooltools.api.rearrange")
    cooler, clr = _open(source)
    try:
        view = pandas.read_csv(view_file, sep="\t", comment="#")
    except (ValueError, OSError) as exc:
        raise ConfigurationError(f"The viewframe must be a tab-separated table: {exc}") from exc
    required = {"chrom", "start", "end", "name", new_chrom, orientation}
    missing = sorted(required - set(view.columns))
    if missing:
        raise ConfigurationError(f"The viewframe needs these columns: {missing}")
    if len(view) > _MAX_VIEW_ROWS:
        raise ConfigurationError(f"The viewframe has more than {_MAX_VIEW_ROWS} regions.")
    before = list(clr.chromnames)
    try:
        rearrange.rearrange_cooler(
            clr,
            view,
            str(target),
            new_chrom_col=new_chrom,
            orientation_col=orientation,
            mode="w",
        )
        after_clr = cooler.Cooler(str(target))
        after = _summary(after_clr)
    except (ValueError, KeyError, OSError) as exc:
        target.unlink(missing_ok=True)
        raise ConfigurationError(f"cooltools could not rearrange this cooler: {exc}") from exc
    return {
        "status": "completed",
        "operation": "rearrange_contact_matrix",
        "input": str(source),
        "output": str(target),
        "chromosomes_before": before,
        "chromosomes_after": after["chromosomes"],
        "output_bins": after["bins"],
        "output_total_contacts": after["total_contacts"],
        "provenance": _provenance(
            "cooltools",
            "bins and pixels re-indexed by the new order (cooltools rearrange_cooler)",
            note="regions with '-' strand are reversed",
        ),
    }


def align_windows_to_bins(
    binsize: int,
    chroms: Sequence[str],
    centers_bp: Sequence[int],
    *,
    flank_bp: int = 0,
    region_start_bp: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Snap genomic loci to bin-aligned windows of a fixed bin size, using cooltools."""

    size = _integer(binsize, "binsize", 1, 10_000_000)
    if isinstance(chroms, (str, bytes)) or not isinstance(chroms, Sequence):
        raise ConfigurationError("chroms must be a list of chromosome names.")
    if isinstance(centers_bp, (str, bytes)) or not isinstance(centers_bp, Sequence):
        raise ConfigurationError("centers_bp must be a list of positions.")
    if len(chroms) != len(centers_bp) or not 1 <= len(chroms) <= _MAX_WINDOWS:
        raise ConfigurationError(
            f"chroms and centers_bp must have the same length, 1 to {_MAX_WINDOWS}."
        )
    names = [_text(name, "chroms") for name in chroms]
    centers = [_integer(value, "centers_bp", 0, 4_000_000_000) for value in centers_bp]
    flank = _integer(flank_bp, "flank_bp", 0, 100_000_000)
    start = _integer(region_start_bp, "region_start_bp", 0, 4_000_000_000)
    (snipping,) = _modules("cooltools.api.snipping")
    try:
        windows = snipping.make_bin_aligned_windows(
            size,
            names,
            centers,
            flank_bp=flank,
            region_start_bp=start,
        )
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"cooltools could not align the windows: {exc}") from exc
    rows = _records(windows)
    return {
        "status": "completed",
        "operation": "align_windows_to_bins",
        "binsize": size,
        "flank_bp": flank,
        "region_start_bp": start,
        "window_count": len(rows),
        "windows": rows,
        "provenance": _provenance(
            "cooltools",
            "window limits moved to bin edges around each locus (cooltools.api.snipping)",
            indices="lo and hi are bin indices counted from region_start_bp",
        ),
    }


def assign_features_to_view(
    features_path: str | Path,
    view_path: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Name the viewframe region that contains each feature (BED), using cooltools."""

    features_file = _existing(features_path, "features_path")
    view_file = _existing(view_path, "view_path")
    pandas = require("pandas", _COOLER_EXTRA)
    (common,) = _modules("cooltools.lib.common")
    try:
        features = pandas.read_csv(
            features_file,
            sep="\t",
            comment="#",
            header=None,
            usecols=[0, 1, 2],
            names=["chrom", "start", "end"],
            dtype={"chrom": str, "start": "int64", "end": "int64"},
        )
        view = pandas.read_csv(
            view_file,
            sep="\t",
            comment="#",
            header=None,
            usecols=[0, 1, 2, 3],
            names=["chrom", "start", "end", "name"],
            dtype={"chrom": str, "start": "int64", "end": "int64", "name": str},
        )
    except (ValueError, OSError) as exc:
        raise ConfigurationError(f"Features and view need BED columns: {exc}") from exc
    if features.empty or view.empty:
        raise ConfigurationError("Features and view must each have at least one interval.")
    try:
        assigned = common.assign_view_auto(features, view)
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"cooltools could not assign features to the view: {exc}") from exc
    rows = _records(assigned[["chrom", "start", "end", "region"]])
    unassigned = sum(1 for row in rows if row["region"] is None)
    return {
        "status": "completed",
        "operation": "assign_features_to_view",
        "features_path": str(features_file),
        "view_path": str(view_file),
        "feature_count": len(rows),
        "unassigned_count": unassigned,
        "features": rows,
        "provenance": _provenance(
            "cooltools",
            "interval overlap of each feature with viewframe regions (cooltools.lib.common)",
            order="feature order is kept; unassigned features have region None",
        ),
    }


def make_chromosome_view(
    cool_path: str | Path,
    *,
    ucsc_names: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Full-chromosome viewframe (one region per chromosome) from a cooler, using cooltools."""

    source = _existing(cool_path, "cool_path")
    flag = _flag(ucsc_names, "ucsc_names")
    (common,) = _modules("cooltools.lib.common")
    _, clr = _open(source)
    try:
        view = common.make_cooler_view(clr, ucsc_names=flag)
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"cooltools could not build a viewframe: {exc}") from exc
    rows = _records(view)
    return {
        "status": "completed",
        "operation": "make_chromosome_view",
        "cool_path": str(source),
        "ucsc_names": flag,
        "region_count": len(rows),
        "regions": rows,
        "provenance": _provenance(
            "cooltools",
            "one region per chromosome from the cooler chromsizes (cooltools.lib.common)",
        ),
    }


def read_viewframe_bed(
    view_path: str | Path,
    *,
    cool_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read and validate a viewframe BED file (unique names, non-overlapping), using cooltools."""

    view_file = _existing(view_path, "view_path")
    (io_module,) = _modules("cooltools.lib.io")
    cooler_clr = None
    if cool_path is not None:
        _, cooler_clr = _open(_existing(cool_path, "cool_path"))
    try:
        view = io_module.read_viewframe_from_file(
            str(view_file),
            verify_cooler=cooler_clr,
            check_sorting=False,
        )
    except (ValueError, KeyError, OSError) as exc:
        raise ConfigurationError(f"The viewframe is not valid: {exc}") from exc
    rows = _records(view)
    return {
        "status": "completed",
        "operation": "read_viewframe_bed",
        "view_path": str(view_file),
        "verified_against": None if cool_path is None else str(cool_path),
        "region_count": len(rows),
        "regions": rows,
        "provenance": _provenance(
            "cooltools",
            "viewframe BED read and checked by cooltools.lib.io.read_viewframe_from_file",
        ),
    }


def read_expected_table(
    expected_path: str | Path,
    *,
    contact_type: str = "cis",
    limit: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read an expected-contacts table (cooltools v1.0 format) and check its schema."""

    source = _existing(expected_path, "expected_path")
    if contact_type not in ("cis", "trans"):
        raise ConfigurationError("contact_type must be 'cis' or 'trans'.")
    rows_limit = _integer(limit, "limit", 1, _MAX_VIEW_ROWS)
    (io_module,) = _modules("cooltools.lib.io")
    try:
        table = io_module.read_expected_from_file(
            str(source),
            contact_type=contact_type,
            expected_value_cols=["count.avg", "balanced.avg"],
            raise_errors=True,
        )
    except (ValueError, KeyError, OSError) as exc:
        raise ConfigurationError(f"The expected table is not valid: {exc}") from exc
    head = table.head(rows_limit)
    return {
        "status": "completed",
        "operation": "read_expected_table",
        "expected_path": str(source),
        "contact_type": contact_type,
        "row_count": len(table),
        "columns": [str(column) for column in table.columns],
        "rows": _records(head),
        "provenance": _provenance(
            "cooltools",
            "expected table checked against the v1.0 schema (cooltools.lib.io)",
            note=f"at most {rows_limit} rows are returned, in file order",
        ),
    }


def observed_over_expected_matrix(
    matrix: Sequence[Sequence[float]],
    *,
    dist_bin_edge_ratio: float = 1.03,
    verbose: bool = True,
) -> dict[str, Any]:
    """Divide a dense symmetric contact matrix by its distance-dependent mean (O/E)."""

    array = _symmetric(_square(matrix, "matrix"), "matrix")
    ratio = _number(dist_bin_edge_ratio, "dist_bin_edge_ratio", 1.0001, _MAX_RATIO)
    numpy = require("numpy", _NUMPY_EXTRA)
    (numutils,) = _modules("cooltools.lib.numutils")
    try:
        oe, dist_bins, sum_pixels, n_pixels = numutils.observed_over_expected(
            array,
            mask=numpy.empty(shape=(0), dtype=numpy.bool_),
            dist_bin_edge_ratio=ratio,
        )
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"cooltools could not compute O/E: {exc}") from exc
    return {
        "status": "completed",
        "operation": "observed_over_expected_matrix",
        "size": int(array.shape[0]),
        "dist_bin_edge_ratio": ratio,
        "observed_over_expected": _matrix(oe),
        "distance_bin_edges": [int(edge) for edge in dist_bins],
        "sum_pixels": _vector(sum_pixels),
        "n_pixels": [int(value) for value in n_pixels],
        "provenance": _provenance(
            "cooltools",
            "diagonal means in log-spaced distance bins divide each diagonal "
            "(cooltools.lib.numutils.observed_over_expected)",
        ),
    }


def iterative_correction_dense(
    matrix: Sequence[Sequence[float]],
    *,
    symmetric: bool = True,
    max_iter: int = 1000,
    ignore_diags: int = 0,
    tol: float = 1e-5,
    verbose: bool = True,
) -> dict[str, Any]:
    """Iterative correction (matrix balancing) of a dense contact matrix, using cooltools."""

    array = _square(matrix, "matrix")
    use_symmetric = _flag(symmetric, "symmetric")
    if use_symmetric:
        _symmetric(array, "matrix")
    iterations_max = _integer(max_iter, "max_iter", 1, _MAX_ITER)
    diags = _integer(ignore_diags, "ignore_diags", 0, array.shape[0] - 1)
    if not use_symmetric and diags != 0:
        raise ConfigurationError("ignore_diags applies to symmetric matrices only; use 0.")
    tolerance = _number(tol, "tol", 1e-12, 1.0)
    (numutils,) = _modules("cooltools.lib.numutils")
    try:
        if use_symmetric:
            corrected, bias, report = numutils.iterative_correction_symmetric(
                array,
                max_iter=iterations_max,
                ignore_diags=diags,
                tol=tolerance,
                verbose=False,
            )
            bias_out: Any = _vector(bias)
        else:
            corrected, bias_row, bias_col, report = numutils.iterative_correction_asymmetric(
                array,
                max_iter=iterations_max,
                tol=tolerance,
                verbose=False,
            )
            bias_out = {"row": _vector(bias_row), "column": _vector(bias_col)}
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"cooltools could not balance the matrix: {exc}") from exc
    converged, iterations = report[0], int(report[1])
    return {
        "status": "completed",
        "operation": "iterative_correction_dense",
        "size": int(array.shape[0]),
        "symmetric": use_symmetric,
        "max_iter": iterations_max,
        "tol": tolerance,
        "ignore_diags": diags,
        "converged": bool(converged),
        "iterations": iterations,
        "corrected": _matrix(corrected),
        "bias": bias_out,
        "provenance": _provenance(
            "cooltools",
            "Knight-Ruiz style iterative correction on dense input "
            "(cooltools.lib.numutils.iterative_correction_symmetric or _asymmetric)",
        ),
    }


def dense_cis_eigenvectors(
    matrix: Sequence[Sequence[float]],
    *,
    n_eigs: int = 3,
    ignore_diags: int = 2,
    clip_percentile: float = 0.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compartment eigenvectors of one dense cis contact matrix, using cooltools."""

    array = _symmetric(_square(matrix, "matrix"), "matrix")
    size = int(array.shape[0])
    eigs = _integer(n_eigs, "n_eigs", 1, min(_MAX_EIGS, size))
    diags = _integer(ignore_diags, "ignore_diags", 0, size - 1)
    clip = _number(clip_percentile, "clip_percentile", 0.0, 100.0)
    (eigdecomp,) = _modules("cooltools.api.eigdecomp")
    try:
        with _reproducible_eigen():
            values, vectors = eigdecomp.cis_eig(
                array,
                n_eigs=eigs,
                phasing_track=None,
                ignore_diags=diags,
                clip_percentile=clip,
                sort_metric=None,
            )
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"cooltools could not decompose the matrix: {exc}") from exc
    return {
        "status": "completed",
        "operation": "dense_cis_eigenvectors",
        "size": size,
        "n_eigs": eigs,
        "ignore_diags": diags,
        "clip_percentile": clip,
        "eigenvalues": _vector(values),
        "eigenvectors": _matrix(_canonical_rows(vectors)),
        "provenance": _provenance(
            "cooltools",
            "eigen-decomposition of the O/E correlation matrix (cooltools.api.eigdecomp.cis_eig)",
            sign="each eigenvector's largest absolute value is positive",
        ),
    }


def expected_from_dense_matrix(
    matrix: Sequence[Sequence[float]],
    *,
    ignore_diags: int = 2,
    region_name: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sum and count of contacts on each diagonal of a dense matrix, using cooltools."""

    array = _square(matrix, "matrix")
    size = int(array.shape[0])
    diags = _integer(ignore_diags, "ignore_diags", 0, size - 1)
    label = None if region_name is None else _text(region_name, "region_name")
    (expected,) = _modules("cooltools.api.expected")
    try:
        table = expected.diagsum_from_array(
            array,
            counts=None,
            offset=0,
            ignore_diags=diags,
            filter_counts=False,
            region_name=label,
        )
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"cooltools could not sum the diagonals: {exc}") from exc
    rows = _records(table)
    return {
        "status": "completed",
        "operation": "expected_from_dense_matrix",
        "size": size,
        "ignore_diags": diags,
        "region_name": label,
        "diagonals": rows,
        "provenance": _provenance(
            "cooltools",
            "per-diagonal valid pixel counts and sums (cooltools.api.expected.diagsum_from_array)",
            order="rows are sorted by distance",
        ),
    }


def digitize_track(
    track_path: str | Path,
    *,
    n_bins: int = 5,
    qrange: Sequence[float] = (0.02, 0.98),
    verbose: bool = True,
) -> dict[str, Any]:
    """Quantile-bin a genomic signal track (bedGraph-like) into integer classes, using cooltools."""

    source = _existing(track_path, "track_path")
    bins = _integer(n_bins, "n_bins", 1, 100)
    if isinstance(qrange, (str, bytes)) or not isinstance(qrange, Sequence) or len(qrange) != 2:
        raise ConfigurationError("qrange must be a pair of quantiles (low, high).")
    low = _number(qrange[0], "qrange low", 0.0, 1.0)
    high = _number(qrange[1], "qrange high", 0.0, 1.0)
    if not low < high:
        raise ConfigurationError("qrange low must be smaller than qrange high.")
    track = _track(source)
    (saddle,) = _modules("cooltools.api.saddle")
    try:
        digitized, edges = saddle.digitize(track, bins, qrange=(low, high))
    except (ValueError, KeyError, TypeError) as exc:
        raise ConfigurationError(f"cooltools could not digitize the track: {exc}") from exc
    rows = [
        {
            "chrom": str(row["chrom"]),
            "start": int(row["start"]),
            "end": int(row["end"]),
            "class": _plain(row["value.d"]),
        }
        for row in digitized.to_dict(orient="records")
    ]
    return {
        "status": "completed",
        "operation": "digitize_track",
        "track_path": str(source),
        "n_bins": bins,
        "qrange": [low, high],
        "bin_edges": _vector(edges),
        "intervals": rows,
        "provenance": _provenance(
            "cooltools",
            "quantile-based classes (cooltools.api.saddle.digitize)",
            codes="1 to n_bins inside qrange; 0 below and n_bins + 1 above; -1 for missing values",
            order="rows follow the input track order",
        ),
    }


__all__ = [
    "align_windows_to_bins",
    "assign_features_to_view",
    "calculate_compartments",
    "calculate_directionality_index",
    "calculate_expected_trans",
    "calculate_logbin_expected",
    "calculate_saddle",
    "calculate_trans_compartments",
    "calculate_virtual_4c",
    "dense_cis_eigenvectors",
    "digitize_track",
    "expected_from_dense_matrix",
    "iterative_correction_dense",
    "make_chromosome_view",
    "observed_over_expected_matrix",
    "pileup_features",
    "read_expected_table",
    "read_viewframe_bed",
    "rearrange_contact_matrix",
    "sample_contacts",
]

_register_public_verbose(globals(), tuple(__all__))
