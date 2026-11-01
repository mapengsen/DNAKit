"""bed-reader wrappers for cloud genotype reads, cloud location checks and bed-reader sample files.

``read_cloud_genotypes`` opens a PLINK 1 .bed file at a URL through bed-reader's ``open_bed`` API,
which
dispatches to the native ``read_cloud_i8``, ``read_cloud_f32`` or ``read_cloud_f64`` reader for the
requested
dtype. Missing genotype calls are ``None`` in the result. ``check_cloud_location`` asks the native
``check_file_cloud`` whether a URL can be read. ``sample_bed_file`` and ``sample_bed_url`` are
network
wrappers: they download bed-reader's sample data with pooch when it is not cached. ``url_to_bytes``
is not
wrapped because it returns raw bytes.

bed-reader, NumPy and pooch are imported only when a function runs. Only the two sample wrappers use
the
network; the cloud functions also accept ``file://`` URLs, which read local files.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.variants.genotype_formats import _check_cells, _limit, _provenance

_BED_EXTRA = "bio-bed-reader"
_DTYPES = ("int8", "float32", "float64")
_INT8_MISSING = -127
_DEFAULT_MAX_CELLS = 5_000_000
_MAX_CELLS_LIMIT = 100_000_000
_MAX_URL_LENGTH = 4096


def _url(value: object, name: str) -> str:
    if not isinstance(value, str) or not value or len(value) > _MAX_URL_LENGTH:
        raise ConfigurationError(f"{name} must be a URL string.")
    if any(char in value for char in "\x00\r\n\t"):
        raise ConfigurationError(f"{name} must not contain control characters.")
    scheme, separator, rest = value.partition("://")
    if not separator or not scheme.isalnum() or not rest:
        raise ConfigurationError(
            f"{name} must be a URL with a scheme, such as https:// or file://; "
            "a plain local path is not a cloud location."
        )
    return value


def _dtype(value: object) -> str:
    if not isinstance(value, str) or value not in _DTYPES:
        raise ConfigurationError("dtype must be 'int8', 'float32' or 'float64'.")
    return value


def _options(value: object) -> dict[str, str]:
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise ConfigurationError("cloud_options must be a dict of str to str.")
    options: dict[str, str] = {}
    for key, item in value.items():
        if not isinstance(key, str) or not isinstance(item, str):
            raise ConfigurationError("cloud_options must map str keys to str values.")
        options[key] = item
    return options


def _indices(value: object, name: str, count: int) -> list[int] | None:
    if value is None:
        return None
    if isinstance(value, (str, bytes, Mapping)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(f"{name} must be a non-empty list of indices, or None for all.")
    indices: list[int] = []
    for item in value:
        if isinstance(item, bool) or not isinstance(item, int) or not 0 <= item < count:
            raise ConfigurationError(f"{name} entries must be integers in [0, {count}).")
        indices.append(item)
    return indices


def _sample_name(value: object) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > 255:
        raise ConfigurationError(
            "filepath must be the name of a bed-reader sample file, such as small.bed."
        )
    if any(char in value for char in "\x00\r\n\t"):
        raise ConfigurationError("filepath must not contain control characters.")
    return value


def _pooch_call(function_name: str, name: str) -> Any:
    bed_module = require("bed_reader", _BED_EXTRA)
    function = getattr(bed_module, function_name)
    try:
        return function(name)
    except ImportError as exc:
        raise BackendUnavailableError(
            f"bed-reader sample data needs pooch: {exc}",
            hint=(
                'Install the sample-data support with: python -m pip install "bed-reader[samples]"'
            ),
        ) from exc


def read_cloud_genotypes(
    location: str,
    *,
    dtype: str = "float32",
    iid_index: list[int] | None = None,
    sid_index: list[int] | None = None,
    cloud_options: dict[str, str] | None = None,
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Genotypes of a PLINK 1 .bed file at a URL, read through bed-reader's open_bed API.

    ``location`` is a URL with a scheme (https://, s3://, file://, ...). The dtype picks
    bed-reader's
    cloud reader: int8 (read_cloud_i8), float32 or float64. ``iid_index`` and ``sid_index`` select
    samples and
    variants by position; None keeps all. Missing calls are ``None``: -127 in int8 and NaN in
    floats.
    ``cloud_options`` are passed to the cloud storage client, and ``max_cells`` bounds the returned
    matrix.
    """

    url = _url(location, "location")
    kind = _dtype(dtype)
    options = _options(cloud_options)
    limit = _limit(max_cells, "max_cells", _MAX_CELLS_LIMIT)
    bed_module = require("bed_reader", _BED_EXTRA)
    numpy = require("numpy", _BED_EXTRA)
    try:
        handle = bed_module.open_bed(url, cloud_options=options)
    except (OSError, ValueError) as exc:
        raise ConfigurationError(f"{url} could not be opened as a PLINK .bed file: {exc}") from exc
    with handle as bed:
        try:
            rows = int(bed.iid_count)
            columns = int(bed.sid_count)
            iid = _indices(iid_index, "iid_index", rows)
            sid = _indices(sid_index, "sid_index", columns)
            selected_rows = rows if iid is None else len(iid)
            selected_columns = columns if sid is None else len(sid)
            _check_cells(selected_rows, selected_columns, limit)
            raw = numpy.asarray(bed.read(index=(iid, sid), dtype=kind))
        except (OSError, ValueError) as exc:
            raise ConfigurationError(f"{url} could not be read: {exc}") from exc
    missing = raw == _INT8_MISSING if kind == "int8" else numpy.isnan(raw)
    values = [
        [None if absent else value for value, absent in zip(row, absent_row, strict=True)]
        for row, absent_row in zip(raw.tolist(), missing.tolist(), strict=True)
    ]
    missing_calls = int(numpy.count_nonzero(missing))
    return {
        "status": "completed",
        "operation": "read_cloud_genotypes",
        "location": url,
        "dtype": kind,
        "shape": [selected_rows, selected_columns],
        "missing_calls": missing_calls,
        "values": values,
        "provenance": _provenance(
            "bed-reader",
            method=(
                "open_bed on a URL, which dispatches to read_cloud_i8, read_cloud_f32 "
                "or read_cloud_f64"
            ),
            missing="None for -127 in int8 and NaN in floats",
        ),
    }


def check_cloud_location(
    location: str,
    *,
    cloud_options: dict[str, str] | None = None,
    verbose: bool = True,
) -> bool:
    """True when bed-reader can open the URL's object for reading, False when it reports it cannot.

    Uses the native ``check_file_cloud``. The native check raises ``ValueError`` for a missing or
    unreadable
    object, which is reported as False; other errors are raised.
    """

    url = _url(location, "location")
    options = _options(cloud_options)
    native = require("bed_reader.bed_reader", _BED_EXTRA)
    try:
        native.check_file_cloud(url, options)
    except ValueError:
        return False
    return True


def sample_bed_file(filepath: str, *, verbose: bool = True) -> dict[str, Any]:
    """Local path of a bed-reader sample .bed file (network: pooch downloads it when it is not
    cached).

    For a .bed name the matching .fam and .bim files are fetched too; their paths are returned when
    present.
    Requires the pooch package, which the bed-reader ``samples`` extra installs.
    """

    name = _sample_name(filepath)
    local = Path(str(_pooch_call("sample_file", name))).resolve()
    fam = bim = None
    if name.lower().endswith(".bed"):
        fam_path = local.with_suffix(".fam")
        bim_path = local.with_suffix(".bim")
        fam = str(fam_path) if fam_path.is_file() else None
        bim = str(bim_path) if bim_path.is_file() else None
    return {
        "service": "bed_reader.sample_file",
        "source": name,
        "result": {"bed": str(local), "fam": fam, "bim": bim},
    }


def sample_bed_url(filepath: str, *, verbose: bool = True) -> dict[str, Any]:
    """URL of a bed-reader sample .bed file (network: pooch fetches the file first when it is not
    cached)."""

    name = _sample_name(filepath)
    url = _pooch_call("sample_url", name)
    return {
        "service": "bed_reader.sample_url",
        "source": name,
        "result": {"url": str(url)},
    }


__all__ = ["check_cloud_location", "read_cloud_genotypes", "sample_bed_file", "sample_bed_url"]

_register_public_verbose(globals())
