"""varcode table readers that the existing varcode wrappers do not reach.

``load_maf_table`` returns every column of a MAF file and ``read_vcf_table`` returns the VCF as the
raw text table that varcode builds, including full per-sample FORMAT strings. Packages are imported
only when a function runs. A missing ``varcode`` raises ``BackendUnavailableError`` with the
``pip install "dnakit[bio-varcode]"`` hint.
"""

from __future__ import annotations

import gzip
import math
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError
from dnakit.variants.varcode_functions import _MAX_ROWS, _VARCODE_EXTRA, _call, _int, _provenance


def _cell(value: Any) -> Any:
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _rows(frame: Any) -> list[dict[str, Any]]:
    return [
        {str(key): _cell(value) for key, value in record.items()}
        for record in frame.to_dict(orient="records")
    ]


def load_maf_table(
    path: str | Path,
    *,
    max_rows: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a TCGA-style MAF file as a table with all of its columns.

    Calls ``varcode.maf.load_maf_dataframe``. The file needs the 19 standard MAF columns in the
    standard order (names are matched ignoring case and renamed to the standard spelling), and may
    have more. ``load_maf_variants`` keeps only contig, start, ref and alt; this function also
    returns gene symbol, variant class, sample barcodes and every other column. Values are text or
    numbers as pandas reads them; empty cells are ``null``. ``max_rows`` limits the rows read.
    """

    source = input_path(path)
    limit = None if max_rows is None else _int(max_rows, "max_rows", minimum=1, maximum=_MAX_ROWS)
    maf = require("varcode.maf", _VARCODE_EXTRA)
    frame = _call(maf.load_maf_dataframe, str(source), nrows=limit)
    if len(frame) > _MAX_ROWS:
        raise ConfigurationError(f"the MAF file has more than {_MAX_ROWS} rows; set max_rows.")
    return {
        "status": "completed",
        "operation": "load_maf_table",
        "path": str(source),
        "row_count": len(frame),
        "columns": [str(name) for name in frame.columns],
        "rows": _rows(frame),
        "provenance": _provenance(),
    }


def _vcf_samples(source: Path) -> list[str]:
    opener = gzip.open if source.name.endswith(".gz") else open
    try:
        with opener(source, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#CHROM"):
                    return line.rstrip("\r\n").split("\t")[9:]
                if not line.startswith("#"):
                    break
    except (OSError, UnicodeDecodeError) as exc:
        raise ConfigurationError(f"{source.name} could not be read as a VCF: {exc}") from exc
    raise ConfigurationError(f"{source.name} has no #CHROM header line.")


def read_vcf_table(
    path: str | Path,
    *,
    include_info: bool = False,
    sample_names: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a VCF as a raw text table with varcode.

    Calls ``varcode.vcf.read_vcf_into_dataframe``. Every column is the text of the VCF field:
    ``ALT`` stays comma separated, ``FILTER`` and ``ID`` are as written, and ``INFO`` (only with
    ``include_info=True``) is not parsed. When ``sample_names`` is given (it needs
    ``include_info=True``), the ``FORMAT`` column and those sample columns are added. varcode only
    labels the first sample columns in file order, so the names must be the first samples of the
    VCF header. Each sample cell is the whole colon-separated string, so genotype quality and depth
    are kept, unlike ``read_vcf``. Header lines are ignored.
    """

    source = input_path(path)
    if not isinstance(include_info, bool):
        raise ConfigurationError("include_info must be a bool.")
    names: list[str] | None = None
    if sample_names is not None:
        if isinstance(sample_names, (str, bytes)) or not isinstance(sample_names, Sequence):
            raise ConfigurationError("sample_names must be a list of sample names.")
        names = [str(item) for item in sample_names]
        if not names or len(set(names)) != len(names):
            raise ConfigurationError("sample_names must be non-empty and unique.")
        if not include_info:
            raise ConfigurationError("sample_names needs include_info=True.")
        header = _vcf_samples(source)
        if header[: len(names)] != names:
            raise ConfigurationError(
                "varcode labels the first sample columns in file order, so sample_names must be "
                f"the first {len(names)} sample(s) of the VCF header: {header[: len(names)]}."
            )
    vcf = require("varcode.vcf", _VARCODE_EXTRA)
    frame = _call(
        vcf.read_vcf_into_dataframe, str(source), include_info=include_info, sample_names=names
    )
    if len(frame) > _MAX_ROWS:
        raise ConfigurationError(f"the VCF file has more than {_MAX_ROWS} rows.")
    return {
        "status": "completed",
        "operation": "read_vcf_table",
        "path": str(source),
        "include_info": include_info,
        "sample_names": names,
        "row_count": len(frame),
        "columns": [str(name) for name in frame.columns],
        "rows": _rows(frame),
        "provenance": _provenance(),
    }


__all__ = ["load_maf_table", "read_vcf_table"]

_register_public_verbose(globals(), tuple(__all__))
