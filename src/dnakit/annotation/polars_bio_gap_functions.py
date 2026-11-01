"""polars-bio wrappers for BGEN, PLINK 2 PGEN, VCF Zarr and CRAM registration, and an interval plot.

Each function reads an explicit local path with polars-bio, converts the result to JSON-compatible
values and never downloads. Remote URLs are not accepted. Rows keep the order polars-bio returns
them in. Coordinates are the ones polars-bio reports; for BGEN, describe_bgen reports 1-based
coordinates. The register functions add the file as a table in the polars-bio session only for the
duration of the call, then remove it. The scan functions keep the lazy frame and materialize at most
``max_rows`` rows. polars-bio, polars and matplotlib are imported only when a function runs.

Functions are importable as ``dnakit.annotation.polars_bio_gap_functions.<name>``.
"""

from __future__ import annotations

import contextlib
import importlib
import math
import re
from collections.abc import Iterator
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

_EXTRA = "bio-polarsbio"
_PLOT_EXTRA = "bio-plots"
_DEFAULT_ROWS = 100_000
_MAX_ROWS = 1_000_000
_MAX_CELLS = 5_000_000
_MAX_PLOT_PAIRS = 50
_PLOT_DPI = 150
_TABLE_NAME = re.compile(r"^[a-z_][a-z0-9_]{0,62}$")
_GENOTYPE_OUTPUTS = ("probability", "dosage")
_PGEN_FIELDS = ("ALT_COUNT", "DS")
_PSAM_ID_MODES = ("iid", "fid_iid", "fid_iid_sid")


# Validation and conversion helpers.


def _integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _text(value: object, name: str, *, maximum: int = 256) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > maximum:
        raise ConfigurationError(
            f"{name} must be a non-empty string of at most {maximum} characters."
        )
    if "\x00" in value:
        raise ConfigurationError(f"{name} cannot contain NUL characters.")
    return value


def _choice(value: object, name: str, options: tuple[str, ...]) -> str:
    if not isinstance(value, str) or value not in options:
        raise ConfigurationError(f"{name} must be one of {list(options)}; got {value!r}.")
    return value


def _names(value: list[str] | None, name: str) -> list[str] | None:
    if value is None:
        return None
    if isinstance(value, (str, bytes)) or not isinstance(value, list) or not value:
        raise ConfigurationError(f"{name} must be a non-empty list of names when given.")
    return [_text(item, name, maximum=128) for item in value]


def _rows_limit(value: object) -> int:
    return _integer(value, "max_rows", 1, _MAX_ROWS)


def _optional_file(value: str | None, name: str) -> str | None:
    if value is None:
        return None
    return str(input_path(_text(value, name, maximum=4096)))


def _local_file(value: object, name: str, suffix: str | None = None) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a local file path.")
    if "://" in value:
        raise ConfigurationError(f"{name} must be a local path; remote URLs are not read.")
    path = input_path(value)
    if suffix is not None and path.suffix.lower() != suffix:
        raise ConfigurationError(f"{name} must be a file ending in {suffix}; got {path.name!r}.")
    return path


def _local_directory(value: object, name: str) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a local directory path.")
    if "://" in value:
        raise ConfigurationError(f"{name} must be a local path; remote URLs are not read.")
    path = Path(value).expanduser().resolve()
    if not path.is_dir():
        raise ConfigurationError(f"{name} is not a directory: {path}")
    return path


def _new_png(value: object) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError("output_path must be a path string.")
    target = Path(value).expanduser().resolve()
    if target.suffix.lower() != ".png":
        raise ConfigurationError("output_path must end with .png.")
    if target.exists():
        raise ConfigurationError(
            f"Output file already exists and will not be overwritten: {target}"
        )
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def _plain(value: Any) -> Any:
    """Convert numpy, polars and missing values to JSON-compatible Python objects."""

    if value is None or type(value).__name__ in {"NAType", "NaTType"}:
        return None
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if hasattr(value, "tolist") and not isinstance(value, (str, bytes)):
        return _plain(value.tolist())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, (bool, int, float, str)):
        return value
    return str(value)


def _materialize(frame: Any, limit: int, what: str) -> Any:
    if hasattr(frame, "collect"):
        frame = frame.collect()
    if frame.height > limit:
        raise ConfigurationError(f"{what} has {frame.height} rows; max_rows={limit} is the limit.")
    return frame


def _records(frame: Any, limit: int, what: str) -> list[dict[str, Any]]:
    frame = _materialize(frame, limit, what)
    return [{str(key): _plain(value) for key, value in row.items()} for row in frame.to_dicts()]


def _schema(frame: Any) -> list[dict[str, str]]:
    schema = frame.collect_schema() if hasattr(frame, "collect_schema") else frame.schema
    return [{"name": str(name), "dtype": str(dtype)} for name, dtype in schema.items()]


def _package(distribution: str) -> str:
    try:
        return f"{distribution} {version(distribution)}"
    except PackageNotFoundError:
        return distribution


def _provenance(convention: str) -> dict[str, str]:
    return {"package": _package("polars-bio"), "convention": convention}


def _polars_bio() -> Any:
    return require("polars_bio", _EXTRA)


def _polars() -> Any:
    return require("polars", _EXTRA)


def _bioframe() -> None:
    """Check for bioframe, which polars-bio's viz extra installs and visualize_intervals needs."""

    try:
        importlib.import_module("bioframe")
    except ImportError as exc:
        raise BackendUnavailableError(
            f"bioframe could not be imported: {exc}",
            hint=(
                "Install the visualization extra of polars-bio: "
                'python -m pip install "polars-bio[viz]".'
            ),
        ) from exc


def _pyplot() -> Any:
    matplotlib = require("matplotlib", _PLOT_EXTRA)
    matplotlib.use("Agg")
    return require("matplotlib.pyplot", _PLOT_EXTRA)


def _session_has_table(pb: Any, name: str) -> bool:
    try:
        pb.ctx.table(name)
    except KeyError:
        return False
    return True


@contextlib.contextmanager
def _registered(pb: Any, table: str, register: Any) -> Iterator[None]:
    """Register one table for the duration of the block and always remove it afterwards."""

    if _session_has_table(pb, table):
        raise ConfigurationError(f"A table named {table!r} is already registered in this session.")
    register()
    try:
        yield
    finally:
        pb.ctx.deregister_table(table)


def _table_name(value: object) -> str:
    name = _text(value, "table_name", maximum=63)
    if not _TABLE_NAME.match(name):
        raise ConfigurationError(
            f"table_name {name!r} must be lowercase letters, digits, or underscores."
        )
    return name


def _bgen_source(path: object) -> Path:
    return _local_file(path, "path", ".bgen")


def _pgen_source(path: object) -> Path:
    return _local_file(path, "path", ".pgen")


def _bgen_options(
    genotype_output: str, samples: list[str] | None, sample_path: str | None, bgi_path: str | None
) -> dict[str, Any]:
    return {
        "genotype_output": _choice(genotype_output, "genotype_output", _GENOTYPE_OUTPUTS),
        "samples": _names(samples, "samples"),
        "sample_path": _optional_file(sample_path, "sample_path"),
        "bgi_path": _optional_file(bgi_path, "bgi_path"),
    }


def _pgen_options(
    genotype_fields: list[str] | None,
    samples: list[str] | None,
    psam_id_mode: str,
    pvar_path: str | None,
    psam_path: str | None,
    pgi_path: str | None,
) -> dict[str, Any]:
    fields = _names(genotype_fields, "genotype_fields") or ["GT"]
    return {
        "genotype_fields": fields,
        "samples": _names(samples, "samples"),
        "psam_id_mode": _choice(psam_id_mode, "psam_id_mode", _PSAM_ID_MODES),
        "pvar_path": _optional_file(pvar_path, "pvar_path"),
        "psam_path": _optional_file(psam_path, "psam_path"),
        "pgi_path": _optional_file(pgi_path, "pgi_path"),
    }


def _vcz_options(
    info_fields: list[str] | None, format_fields: list[str] | None, samples: list[str] | None
) -> dict[str, Any]:
    return {
        "info_fields": _names(info_fields, "info_fields"),
        "format_fields": _names(format_fields, "format_fields"),
        "samples": _names(samples, "samples"),
    }


def _matrix_result(matrix: Any, operation: str, convention: str) -> dict[str, Any]:
    values = matrix.values
    cells = int(values.size)
    if cells > _MAX_CELLS:
        raise ConfigurationError(f"The matrix has {cells} cells; the limit is {_MAX_CELLS}.")
    positions = [int(item) for item in matrix.positions.tolist()]
    sample_names = [str(item) for item in matrix.sample_names]
    return {
        "status": "completed",
        "operation": operation,
        "variant_count": len(positions),
        "sample_count": len(sample_names),
        "positions": positions,
        "sample_names": sample_names,
        "values": _plain(values.tolist()),
        "provenance": _provenance(convention),
    }


# BGEN files.


def polars_bio_read_bgen(
    path: str,
    *,
    genotype_output: str = "probability",
    samples: list[str] | None = None,
    sample_path: str | None = None,
    bgi_path: str | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a local BGEN file into one record per variant with its probabilities or dosages."""

    source = _bgen_source(path)
    options = _bgen_options(genotype_output, samples, sample_path, bgi_path)
    limit = _rows_limit(max_rows)
    pb = _polars_bio()
    variants = _records(pb.read_bgen(str(source), **options), limit, "BGEN file")
    return {
        "status": "completed",
        "operation": "read_bgen",
        "genotype_output": options["genotype_output"],
        "variant_count": len(variants),
        "variants": variants,
        "provenance": _provenance("BGEN read by polars-bio; coordinates are 1-based"),
    }


def polars_bio_scan_bgen(
    path: str,
    *,
    genotype_output: str = "probability",
    samples: list[str] | None = None,
    sample_path: str | None = None,
    bgi_path: str | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Scan a local BGEN file lazily: count all variants but materialize at most ``max_rows``."""

    source = _bgen_source(path)
    options = _bgen_options(genotype_output, samples, sample_path, bgi_path)
    limit = _rows_limit(max_rows)
    pb = _polars_bio()
    pl = _polars()
    lazy = pb.scan_bgen(str(source), **options)
    total = int(lazy.select(pl.len()).collect().item())
    head = _records(lazy.head(limit), limit, "BGEN head")
    return {
        "status": "completed",
        "operation": "scan_bgen",
        "columns": _schema(lazy),
        "variant_count": total,
        "returned_count": len(head),
        "truncated": total > len(head),
        "variants": head,
        "provenance": _provenance(
            "BGEN scanned lazily by polars-bio; only the first rows are materialized"
        ),
    }


def polars_bio_describe_bgen(
    path: str,
    *,
    sample_path: str | None = None,
    bgi_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Describe the schema a local BGEN file produces, one entry per output column."""

    source = _bgen_source(path)
    pb = _polars_bio()
    entries = _records(
        pb.describe_bgen(
            str(source),
            sample_path=_optional_file(sample_path, "sample_path"),
            bgi_path=_optional_file(bgi_path, "bgi_path"),
        ),
        _DEFAULT_ROWS,
        "BGEN description",
    )
    return {
        "status": "completed",
        "operation": "describe_bgen",
        "field_count": len(entries),
        "fields": entries,
        "provenance": _provenance("schema entries reported by polars-bio for the BGEN file"),
    }


def polars_bio_register_bgen(
    path: str,
    table_name: str,
    *,
    genotype_output: str = "probability",
    samples: list[str] | None = None,
    sample_path: str | None = None,
    bgi_path: str | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Register a local BGEN file as a table, query all of it with SQL, then remove the table."""

    source = _bgen_source(path)
    name = _table_name(table_name)
    options = _bgen_options(genotype_output, samples, sample_path, bgi_path)
    limit = _rows_limit(max_rows)
    pb = _polars_bio()
    with _registered(pb, name, lambda: pb.register_bgen(str(source), name, **options)):
        frame = _materialize(pb.sql(f"SELECT * FROM {name} LIMIT {limit + 1}"), limit, "BGEN table")
        rows = [{str(k): _plain(v) for k, v in row.items()} for row in frame.to_dicts()]
        columns = list(frame.columns)
    return {
        "status": "completed",
        "operation": "register_bgen",
        "table_name": name,
        "columns": columns,
        "row_count": len(rows),
        "rows": rows,
        "provenance": _provenance("BGEN registered in the polars-bio session for this call only"),
    }


def polars_bio_read_bgen_matrix(
    path: str,
    *,
    samples: list[str] | None = None,
    missing: float | None = None,
    sample_path: str | None = None,
    bgi_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read the ALT dosages of a local BGEN file into a dense variants-by-samples matrix."""

    source = _bgen_source(path)
    names = _names(samples, "samples")
    if missing is not None and (isinstance(missing, bool) or not isinstance(missing, (int, float))):
        raise ConfigurationError("missing must be a number or None.")
    pb = _polars_bio()
    matrix = pb.read_bgen_matrix(
        str(source),
        samples=names,
        missing=None if missing is None else float(missing),
        sample_path=_optional_file(sample_path, "sample_path"),
        bgi_path=_optional_file(bgi_path, "bgi_path"),
    )
    return _matrix_result(
        matrix,
        "read_bgen_matrix",
        "ALT dosages per variant (rows) and sample (columns); coordinates are 1-based",
    )


# PLINK 2 PGEN filesets.


def polars_bio_read_pgen(
    path: str,
    *,
    genotype_fields: list[str] | None = None,
    samples: list[str] | None = None,
    psam_id_mode: str = "iid",
    pvar_path: str | None = None,
    psam_path: str | None = None,
    pgi_path: str | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a local PLINK 2 PGEN fileset (the .pgen path, with .pvar and .psam beside it)."""

    source = _pgen_source(path)
    options = _pgen_options(genotype_fields, samples, psam_id_mode, pvar_path, psam_path, pgi_path)
    limit = _rows_limit(max_rows)
    pb = _polars_bio()
    variants = _records(pb.read_pgen(str(source), **options), limit, "PGEN fileset")
    return {
        "status": "completed",
        "operation": "read_pgen",
        "variant_count": len(variants),
        "variants": variants,
        "provenance": _provenance("PLINK 2 PGEN read by polars-bio; coordinates are 1-based"),
    }


def polars_bio_scan_pgen(
    path: str,
    *,
    genotype_fields: list[str] | None = None,
    samples: list[str] | None = None,
    psam_id_mode: str = "iid",
    pvar_path: str | None = None,
    psam_path: str | None = None,
    pgi_path: str | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Scan a local PLINK 2 PGEN fileset lazily: count all variants, materialize at most
    ``max_rows``."""

    source = _pgen_source(path)
    options = _pgen_options(genotype_fields, samples, psam_id_mode, pvar_path, psam_path, pgi_path)
    limit = _rows_limit(max_rows)
    pb = _polars_bio()
    pl = _polars()
    lazy = pb.scan_pgen(str(source), **options)
    total = int(lazy.select(pl.len()).collect().item())
    head = _records(lazy.head(limit), limit, "PGEN head")
    return {
        "status": "completed",
        "operation": "scan_pgen",
        "columns": _schema(lazy),
        "variant_count": total,
        "returned_count": len(head),
        "truncated": total > len(head),
        "variants": head,
        "provenance": _provenance("PLINK 2 PGEN scanned lazily by polars-bio"),
    }


def polars_bio_describe_pgen(
    path: str,
    *,
    pvar_path: str | None = None,
    psam_path: str | None = None,
    pgi_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Describe the schema a local PLINK 2 PGEN fileset produces, one entry per output column."""

    source = _pgen_source(path)
    pb = _polars_bio()
    entries = _records(
        pb.describe_pgen(
            str(source),
            pvar_path=_optional_file(pvar_path, "pvar_path"),
            psam_path=_optional_file(psam_path, "psam_path"),
            pgi_path=_optional_file(pgi_path, "pgi_path"),
        ),
        _DEFAULT_ROWS,
        "PGEN description",
    )
    return {
        "status": "completed",
        "operation": "describe_pgen",
        "field_count": len(entries),
        "fields": entries,
        "provenance": _provenance("schema entries reported by polars-bio for the PGEN fileset"),
    }


def polars_bio_register_pgen(
    path: str,
    table_name: str,
    *,
    genotype_fields: list[str] | None = None,
    samples: list[str] | None = None,
    psam_id_mode: str = "iid",
    pvar_path: str | None = None,
    psam_path: str | None = None,
    pgi_path: str | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Register a local PLINK 2 PGEN fileset as a table, query all of it with SQL, then remove
    it."""

    source = _pgen_source(path)
    name = _table_name(table_name)
    options = _pgen_options(genotype_fields, samples, psam_id_mode, pvar_path, psam_path, pgi_path)
    limit = _rows_limit(max_rows)
    pb = _polars_bio()
    with _registered(pb, name, lambda: pb.register_pgen(str(source), name, **options)):
        frame = _materialize(pb.sql(f"SELECT * FROM {name} LIMIT {limit + 1}"), limit, "PGEN table")
        rows = [{str(k): _plain(v) for k, v in row.items()} for row in frame.to_dicts()]
        columns = list(frame.columns)
    return {
        "status": "completed",
        "operation": "register_pgen",
        "table_name": name,
        "columns": columns,
        "row_count": len(rows),
        "rows": rows,
        "provenance": _provenance("PGEN registered in the polars-bio session for this call only"),
    }


def polars_bio_read_pgen_matrix(
    path: str,
    *,
    field: str = "ALT_COUNT",
    samples: list[str] | None = None,
    missing: float | None = None,
    psam_id_mode: str = "iid",
    pvar_path: str | None = None,
    psam_path: str | None = None,
    pgi_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read one genotype field (ALT_COUNT or DS) of a local PGEN fileset into a dense matrix."""

    source = _pgen_source(path)
    choice = _choice(field, "field", _PGEN_FIELDS)
    names = _names(samples, "samples")
    if missing is not None and (isinstance(missing, bool) or not isinstance(missing, (int, float))):
        raise ConfigurationError("missing must be a number or None.")
    pb = _polars_bio()
    matrix = pb.read_pgen_matrix(
        str(source),
        field=choice,
        samples=names,
        missing=missing,
        psam_id_mode=_choice(psam_id_mode, "psam_id_mode", _PSAM_ID_MODES),
        pvar_path=_optional_file(pvar_path, "pvar_path"),
        psam_path=_optional_file(psam_path, "psam_path"),
        pgi_path=_optional_file(pgi_path, "pgi_path"),
    )
    return _matrix_result(
        matrix,
        "read_pgen_matrix",
        f"{choice} per variant (rows) and sample (columns); coordinates are 1-based",
    )


# VCF Zarr stores.


def polars_bio_read_vcf_zarr(
    path: str,
    *,
    info_fields: list[str] | None = None,
    format_fields: list[str] | None = None,
    samples: list[str] | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a local VCF Zarr store (a directory) into one record per variant."""

    source = _local_directory(path, "path")
    options = _vcz_options(info_fields, format_fields, samples)
    limit = _rows_limit(max_rows)
    pb = _polars_bio()
    variants = _records(pb.read_vcf_zarr(str(source), **options), limit, "VCF Zarr store")
    return {
        "status": "completed",
        "operation": "read_vcf_zarr",
        "variant_count": len(variants),
        "variants": variants,
        "provenance": _provenance("VCF Zarr store read by polars-bio; coordinates are 1-based"),
    }


def polars_bio_scan_vcf_zarr(
    path: str,
    *,
    info_fields: list[str] | None = None,
    format_fields: list[str] | None = None,
    samples: list[str] | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Scan a local VCF Zarr store lazily: count all variants, materialize at most ``max_rows``."""

    source = _local_directory(path, "path")
    options = _vcz_options(info_fields, format_fields, samples)
    limit = _rows_limit(max_rows)
    pb = _polars_bio()
    pl = _polars()
    lazy = pb.scan_vcf_zarr(str(source), **options)
    total = int(lazy.select(pl.len()).collect().item())
    head = _records(lazy.head(limit), limit, "VCF Zarr head")
    return {
        "status": "completed",
        "operation": "scan_vcf_zarr",
        "columns": _schema(lazy),
        "variant_count": total,
        "returned_count": len(head),
        "truncated": total > len(head),
        "variants": head,
        "provenance": _provenance("VCF Zarr store scanned lazily by polars-bio"),
    }


def polars_bio_describe_vcf_zarr(path: str, *, verbose: bool = True) -> dict[str, Any]:
    """Describe the INFO and FORMAT schema of a local VCF Zarr store."""

    source = _local_directory(path, "path")
    pb = _polars_bio()
    entries = _records(pb.describe_vcf_zarr(str(source)), _DEFAULT_ROWS, "VCF Zarr description")
    return {
        "status": "completed",
        "operation": "describe_vcf_zarr",
        "field_count": len(entries),
        "fields": entries,
        "provenance": _provenance("INFO and FORMAT schema entries reported by polars-bio"),
    }


def polars_bio_register_vcf_zarr(
    path: str,
    table_name: str,
    *,
    info_fields: list[str] | None = None,
    format_fields: list[str] | None = None,
    samples: list[str] | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Register a local VCF Zarr store as a table, query all of it with SQL, then remove the
    table."""

    source = _local_directory(path, "path")
    name = _table_name(table_name)
    options = _vcz_options(info_fields, format_fields, samples)
    limit = _rows_limit(max_rows)
    pb = _polars_bio()
    with _registered(pb, name, lambda: pb.register_vcf_zarr(str(source), name, **options)):
        frame = _materialize(
            pb.sql(f"SELECT * FROM {name} LIMIT {limit + 1}"), limit, "VCF Zarr table"
        )
        rows = [{str(k): _plain(v) for k, v in row.items()} for row in frame.to_dicts()]
        columns = list(frame.columns)
    return {
        "status": "completed",
        "operation": "register_vcf_zarr",
        "table_name": name,
        "columns": columns,
        "row_count": len(rows),
        "rows": rows,
        "provenance": _provenance(
            "VCF Zarr store registered in the polars-bio session for this call only"
        ),
    }


# CRAM files.


def polars_bio_register_cram(
    path: str,
    table_name: str,
    *,
    tag_fields: list[str] | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Register a local CRAM file as a table, query all of it with SQL, then remove the table.

    polars-bio's CRAM registration takes no reference FASTA, so it decodes only CRAM files whose
    sequences are embedded. For CRAM files that need an external reference, use read_cram_records,
    which takes ``reference_path``.
    """

    source = _local_file(path, "path")
    name = _table_name(table_name)
    tags = _names(tag_fields, "tag_fields")
    limit = _rows_limit(max_rows)
    pb = _polars_bio()
    with _registered(pb, name, lambda: pb.register_cram(str(source), name, tag_fields=tags)):
        frame = _materialize(pb.sql(f"SELECT * FROM {name} LIMIT {limit + 1}"), limit, "CRAM table")
        rows = [{str(k): _plain(v) for k, v in row.items()} for row in frame.to_dicts()]
        columns = list(frame.columns)
    return {
        "status": "completed",
        "operation": "register_cram",
        "table_name": name,
        "columns": columns,
        "row_count": len(rows),
        "rows": rows,
        "provenance": _provenance(
            "CRAM registered in the polars-bio session; embedded sequences only"
        ),
    }


# Interval plot.


def polars_bio_visualize_intervals(
    left_bed: str,
    right_bed: str,
    output_path: str,
    *,
    label: str = "overlapping pair",
    max_pairs: int = 20,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw the overlapping pairs of two local BED files, one panel per pair, into one PNG.

    The pairs come from polars-bio's overlap of the two files. polars-bio's visualize_intervals
    draws each pair on its own figure (bioframe opens a new figure per call); those figures are
    captured and stacked, in pair order, into the single output image. The first ``max_pairs``
    pairs are drawn, and ``pair_count`` reports all of them.
    """

    left = _local_file(left_bed, "left_bed")
    right = _local_file(right_bed, "right_bed")
    target = _new_png(output_path)
    text = _text(label, "label", maximum=128)
    limit = _integer(max_pairs, "max_pairs", 1, _MAX_PLOT_PAIRS)
    pb = _polars_bio()
    _bioframe()
    overlaps = pb.overlap(
        pb.read_bed(str(left), use_zero_based=True),
        pb.read_bed(str(right), use_zero_based=True),
        suffixes=("_1", "_2"),
        output_type="polars.DataFrame",
    )
    if overlaps.height == 0:
        raise ConfigurationError("The two BED files have no overlapping intervals to draw.")
    shown = overlaps.head(limit)
    panels = shown.to_pandas()
    plt = _pyplot()
    np = require("numpy", _PLOT_EXTRA)
    tiles: list[Any] = []
    created: list[int] = []
    try:
        for index in range(len(panels)):
            before = set(plt.get_fignums())
            pb.visualize_intervals(panels.iloc[[index]], label=text)
            fresh = [number for number in plt.get_fignums() if number not in before]
            created.extend(fresh)
            pair_figure = plt.figure(fresh[-1]) if fresh else plt.gcf()
            pair_figure.canvas.draw()
            tiles.append(np.asarray(pair_figure.canvas.buffer_rgba()).copy())
        figure, axes = plt.subplots(len(tiles), 1, figsize=(10.0, 2.2 * len(tiles)), squeeze=False)
        for axis, tile in zip(axes[:, 0], tiles, strict=True):
            axis.imshow(tile)
            axis.axis("off")
        figure.subplots_adjust(left=0, right=1, bottom=0, top=1, hspace=0.05)
        figure.savefig(target, format="png", dpi=_PLOT_DPI)
        plt.close(figure)
    finally:
        for number in created:
            plt.close(number)
    return {
        "path": str(target),
        "data": {
            "label": text,
            "pair_count": overlaps.height,
            "plotted_count": shown.height,
            "pairs": _records(shown, limit, "overlap pairs"),
            "provenance": _provenance("BED intervals, 0-based half-open; overlaps from polars-bio"),
        },
    }


__all__ = [
    "polars_bio_describe_bgen",
    "polars_bio_describe_pgen",
    "polars_bio_describe_vcf_zarr",
    "polars_bio_read_bgen",
    "polars_bio_read_bgen_matrix",
    "polars_bio_read_pgen",
    "polars_bio_read_pgen_matrix",
    "polars_bio_read_vcf_zarr",
    "polars_bio_register_bgen",
    "polars_bio_register_cram",
    "polars_bio_register_pgen",
    "polars_bio_register_vcf_zarr",
    "polars_bio_scan_bgen",
    "polars_bio_scan_pgen",
    "polars_bio_scan_vcf_zarr",
    "polars_bio_visualize_intervals",
]

_register_public_verbose(globals())
