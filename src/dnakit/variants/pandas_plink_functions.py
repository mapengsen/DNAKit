"""Relationship-matrix readers for GCTA GRM and PLINK REL files through pandas-plink.

- ``read_gcta_relationship_matrix`` reads a GCTA GRM file set (text or binary ``.grm``/``.grm.bin``,
  with the ``.id`` file and, for binary files, the ``.N.bin`` SNP-count file).
- ``read_plink_relationship_matrix`` reads a PLINK 2 relationship matrix (``.rel``, ``.rel.bin``).

Both return dense square matrices in sample order. The number of samples is limited by
``max_cells`` (samples squared) before any matrix is read. Compressed ``.zst`` REL files need the
``zstandard`` package. Packages are imported only when a function runs.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError
from dnakit.variants.genotype_formats import _limit, _provenance

_PANDAS_PLINK_EXTRA = "bio-pandas-plink"
_DEFAULT_MAX_CELLS = 5_000_000
_MAX_CELLS_LIMIT = 100_000_000


def _file(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip() or "\x00" in str(value):
        raise ConfigurationError(f"{name} must be a non-empty file path.")
    return input_path(value)


def _optional_file(value: object, name: str) -> Path | None:
    return None if value is None else _file(value, name)


def _grm_id_name(name: str) -> str:
    # Same rule as pandas-plink: binary "x.bin" maps to "x.id"; text "x.gz" maps to "x.id".
    if name.endswith(".bin"):
        return name.rsplit(".bin", 1)[0] + ".id"
    if name.endswith(".gz"):
        return name[:-3] + ".id"
    return name + ".id"


def _rel_id_name(name: str) -> str:
    # Same rule as pandas-plink: strip ".bin", ".zst" or ".gz" before adding ".id".
    if name.endswith((".bin", ".zst")):
        return name[:-4] + ".id"
    if name.endswith(".gz"):
        return name[:-3] + ".id"
    return name + ".id"


def _sample_ids(id_path: Path, max_samples: int) -> list[dict[str, str]]:
    """Read family and individual IDs (tab separated), refusing more than max_samples rows."""

    samples: list[dict[str, str]] = []
    with id_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            text = line.rstrip("\r\n")
            if not text.strip():
                continue
            fields = text.split("\t")
            if len(fields) < 2:
                raise ConfigurationError(f"{id_path.name} needs family and individual columns.")
            if len(samples) >= max_samples:
                raise ConfigurationError(
                    f"{id_path.name} lists more samples than max_cells allows "
                    f"(max {max_samples} samples); raise max_cells or subset the file."
                )
            samples.append({"family_id": fields[0], "sample_id": fields[1]})
    if not samples:
        raise ConfigurationError(f"{id_path.name} contains no samples.")
    return samples


def _max_samples(max_cells: int) -> int:
    """Largest square matrix side allowed by max_cells."""

    return math.isqrt(_limit(max_cells, "max_cells", _MAX_CELLS_LIMIT))


def read_gcta_relationship_matrix(
    path: object,
    *,
    id_path: object = None,
    n_snps_path: object = None,
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a GCTA genetic relationship matrix (GRM) with its sample and SNP-count files."""

    max_samples = _max_samples(max_cells)
    matrix_path = _file(path, "path")
    id_file = (
        _file(id_path, "id_path")
        if id_path is not None
        else matrix_path.with_name(_grm_id_name(matrix_path.name))
    )
    if id_path is None and not id_file.is_file():
        raise ConfigurationError(f"GRM sample file does not exist: {id_file}")
    n_snps_file = _optional_file(n_snps_path, "n_snps_path")
    samples = _sample_ids(id_file, max_samples)
    module = require("pandas_plink", _PANDAS_PLINK_EXTRA)
    kinship, n_snps = module.read_grm(
        str(matrix_path),
        id_filepath=str(id_file),
        n_snps_filepath=None if n_snps_file is None else str(n_snps_file),
    )
    matrix = kinship.values.tolist()
    if len(matrix) != len(samples):
        raise ConfigurationError("GRM sample file does not match the matrix size.")
    return {
        "status": "completed",
        "operation": "read_gcta_relationship_matrix",
        "samples": samples,
        "matrix": matrix,
        "n_snps": [int(value) for value in n_snps],
        "provenance": _provenance(
            "pandas-plink",
            source="read_grm",
            convention="n_snps holds one count per stored pair, in file order",
        ),
    }


def read_plink_relationship_matrix(
    path: object,
    *,
    id_path: object = None,
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a PLINK 2 realized relationship matrix (.rel, .rel.bin) with its sample file."""

    max_samples = _max_samples(max_cells)
    matrix_path = _file(path, "path")
    id_file = (
        _file(id_path, "id_path")
        if id_path is not None
        else matrix_path.with_name(_rel_id_name(matrix_path.name))
    )
    if id_path is None and not id_file.is_file():
        raise ConfigurationError(f"REL sample file does not exist: {id_file}")
    samples = _sample_ids(id_file, max_samples)
    module = require("pandas_plink", _PANDAS_PLINK_EXTRA)
    kinship = module.read_rel(str(matrix_path), id_filepath=str(id_file))
    matrix = kinship.values.tolist()
    if len(matrix) != len(samples):
        raise ConfigurationError("REL sample file does not match the matrix size.")
    return {
        "status": "completed",
        "operation": "read_plink_relationship_matrix",
        "samples": samples,
        "matrix": matrix,
        "provenance": _provenance(
            "pandas-plink",
            source="read_rel",
            convention="text and binary REL files are both returned as full square matrices",
        ),
    }


__all__ = [
    "read_gcta_relationship_matrix",
    "read_plink_relationship_matrix",
]

_register_public_verbose(globals(), tuple(__all__))
