"""Convert a VCF file to a Zarr store with scikit-allel's reader and zarr 3's writer.

scikit-allel 1.3.13 ``vcf_to_zarr`` writes its arrays with the zarr 2 API (``Group.create_dataset``),
which zarr 3.4 does not have, so the call fails after the first chunk is read. ``allel_vcf_to_zarr``
keeps the same layout (``variants/<FIELD>``, ``calldata/<FIELD>``, ``samples``) and reads the file
with ``allel.read_vcf``; the arrays are written with ``create_array``. Strings are stored as zarr
variable-length strings. The store is a new directory; an existing path is never overwritten.
"""

from __future__ import annotations

import importlib.metadata
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_ALLEL_EXTRA = "bio-scikit-allel"
_MAX_FIELDS = 200


def _field_names(value: object) -> list[str] | None:
    if value is None:
        return None
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(
            "fields must be a non-empty list of field names, or None for the defaults."
        )
    if len(value) > _MAX_FIELDS:
        raise ConfigurationError(f"fields lists more than {_MAX_FIELDS} names.")
    names: list[str] = []
    for item in value:
        if not isinstance(item, str) or not item.strip():
            raise ConfigurationError("each field name must be a non-empty string.")
        names.append(item.strip())
    return names


def _version(package: str) -> str:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:  # pragma: no cover
        return "unknown"


def allel_vcf_to_zarr(
    vcf_path: str | Path,
    zarr_path: str | Path,
    *,
    fields: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a VCF file and write its variant and genotype arrays to a new Zarr store directory.

    ``fields`` takes names such as ``variants/CHROM`` or ``calldata/GT`` (scikit-allel's own names);
    ``None`` reads scikit-allel's default fields. The result lists each array with its shape and
    dtype. When the VCF has no variants no store is written and ``output_path`` is ``None``, as
    scikit-allel does.
    """
    allel = require("allel", _ALLEL_EXTRA)
    zarr = require("zarr", _ALLEL_EXTRA)
    np = require("numpy", _ALLEL_EXTRA)
    source = Path(vcf_path).expanduser()
    if not source.is_file():
        raise ConfigurationError(f"vcf_path must be an existing file: {source}")
    target = Path(zarr_path).expanduser().resolve()
    if target.exists():
        raise ConfigurationError(f"Output store already exists: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {target.parent}")
    field_list = _field_names(fields)
    arrays = allel.read_vcf(str(source), fields=field_list)
    if not arrays:
        return {
            "status": "completed",
            "operation": "allel_vcf_to_zarr",
            "output_path": None,
            "variants": 0,
            "arrays": [],
            "provenance": {
                "package": f"scikit-allel {_version('scikit-allel')}",
                "convention": "no variants",
            },
        }
    root = zarr.open_group(str(target), mode="w")
    summary = []
    for name in sorted(arrays):
        values = np.asarray(arrays[name])
        if values.dtype.kind == "O":
            values = values.astype(str)
        if "/" in name:
            group_name, leaf = name.split("/", 1)
            group = root.require_group(group_name)
        else:
            leaf, group = name, root
        group.create_array(leaf, data=values)
        summary.append({"name": name, "shape": list(values.shape), "dtype": str(values.dtype)})
    variants = (
        int(np.asarray(arrays["variants/CHROM"]).shape[0]) if "variants/CHROM" in arrays else 0
    )
    return {
        "status": "completed",
        "operation": "allel_vcf_to_zarr",
        "output_path": str(target),
        "variants": variants,
        "samples": int(np.asarray(arrays["samples"]).shape[0]) if "samples" in arrays else 0,
        "arrays": summary,
        "provenance": {
            "package": f"scikit-allel {_version('scikit-allel')}; zarr {_version('zarr')}",
            "convention": "vcf_to_zarr layout (variants/, calldata/, samples); arrays written with zarr 3 create_array",
        },
    }


__all__ = ["allel_vcf_to_zarr"]

_register_public_verbose(globals())
