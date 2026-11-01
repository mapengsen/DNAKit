"""deepTools computeMatrixOperations subcommands as typed DNAKit facades (operate on computeMatrix output).

File-writing subcommands (subset, filterStrand, filterValues, sort, relabel, rbind, cbind) return the output
path; text subcommands (info, dataRange) return the printed text. Same contract as the other deepTools facades.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal, Sequence

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.backends.external import execute_bounded_command
from dnakit.epigenome.deeptools_cli_facades import _existing, _new, _program, _run
from dnakit.exceptions import ConfigurationError

_SUB = "computeMatrixOperations"


def _text(args: list[str], cwd: Path, deeptools_path: str | None) -> str:
    program = _program(_SUB, deeptools_path)
    result = execute_bounded_command(program, [*args], backend_id=f"deeptools:{_SUB}", cwd=cwd, verbose=False)
    if result.return_code != 0:
        raise ConfigurationError(f"computeMatrixOperations failed: {result.output[-500:]}")
    return result.output


def matrix_info(matrix: str, *, deeptools_path: str | None = None, verbose: bool = True) -> str:
    """Return the group and sample layout of a computeMatrix output (computeMatrixOperations info)."""
    path = _existing(matrix, "matrix")
    return _text(["info", "-m", str(path)], path.parent, deeptools_path)


def data_range(matrix: str, *, deeptools_path: str | None = None, verbose: bool = True) -> str:
    """Return the minimum and maximum values of a computeMatrix output (computeMatrixOperations dataRange)."""
    path = _existing(matrix, "matrix")
    return _text(["dataRange", "-m", str(path)], path.parent, deeptools_path)


def subset_matrix(
    matrix: str,
    output: str,
    *,
    groups: Sequence[str] | None = None,
    samples: Sequence[str] | None = None,
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep only selected groups and/or samples of a computeMatrix output (computeMatrixOperations subset)."""
    path = _existing(matrix, "matrix")
    out = _new(output, "output")
    args = ["subset", "-m", str(path), "-o", str(out)]
    if groups:
        args += ["--groups", *list(groups)]
    if samples:
        args += ["--samples", *list(samples)]
    if not groups and not samples:
        raise ConfigurationError("subset_matrix needs groups or samples to keep.")
    return _run(_SUB, args, [out], deeptools_path)


def filter_strand(
    matrix: str,
    output: str,
    *,
    strand: Literal["+", "-", "."],
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep regions on one strand (computeMatrixOperations filterStrand)."""
    path = _existing(matrix, "matrix")
    out = _new(output, "output")
    return _run(_SUB, ["filterStrand", "-m", str(path), "-o", str(out), "-s", strand], [out], deeptools_path)


def filter_values(
    matrix: str,
    output: str,
    *,
    minimum: float | None = None,
    maximum: float | None = None,
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep rows whose values lie within a range (computeMatrixOperations filterValues)."""
    path = _existing(matrix, "matrix")
    out = _new(output, "output")
    args = ["filterValues", "-m", str(path), "-o", str(out)]
    if minimum is not None:
        args += ["--min", str(minimum)]
    if maximum is not None:
        args += ["--max", str(maximum)]
    if minimum is None and maximum is None:
        raise ConfigurationError("filter_values needs a minimum or a maximum.")
    return _run(_SUB, args, [out], deeptools_path)


def sort_matrix(
    matrix: str,
    output: str,
    regions: Sequence[str],
    *,
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sort rows of a computeMatrix output by region order (computeMatrixOperations sort)."""
    path = _existing(matrix, "matrix")
    out = _new(output, "output")
    if not regions:
        raise ConfigurationError("sort_matrix needs the region file(s) that define the row order.")
    args = ["sort", "-m", str(path), "-o", str(out), "-R", *[str(_existing(r, "regions")) for r in regions]]
    return _run(_SUB, args, [out], deeptools_path)


def relabel_matrix(
    matrix: str,
    output: str,
    *,
    group_labels: Sequence[str] | None = None,
    sample_labels: Sequence[str] | None = None,
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Rename groups and/or samples of a computeMatrix output (computeMatrixOperations relabel)."""
    path = _existing(matrix, "matrix")
    out = _new(output, "output")
    args = ["relabel", "-m", str(path), "-o", str(out)]
    if group_labels:
        args += ["--groupLabels", *list(group_labels)]
    if sample_labels:
        args += ["--sampleLabels", *list(sample_labels)]
    if not group_labels and not sample_labels:
        raise ConfigurationError("relabel_matrix needs group_labels or sample_labels.")
    return _run(_SUB, args, [out], deeptools_path)


def bind_rows(matrices: Sequence[str], output: str, *, deeptools_path: str | None = None, verbose: bool = True) -> dict[str, Any]:
    """Stack computeMatrix outputs row-wise (computeMatrixOperations rbind)."""
    if len(matrices) < 2:
        raise ConfigurationError("bind_rows needs at least two matrices.")
    paths = [str(_existing(m, "matrix")) for m in matrices]
    out = _new(output, "output")
    return _run(_SUB, ["rbind", "-m", *paths, "-o", str(out)], [out], deeptools_path)


def bind_columns(matrices: Sequence[str], output: str, *, deeptools_path: str | None = None, verbose: bool = True) -> dict[str, Any]:
    """Join computeMatrix outputs column-wise (computeMatrixOperations cbind)."""
    if len(matrices) < 2:
        raise ConfigurationError("bind_columns needs at least two matrices.")
    paths = [str(_existing(m, "matrix")) for m in matrices]
    out = _new(output, "output")
    return _run(_SUB, ["cbind", "-m", *paths, "-o", str(out)], [out], deeptools_path)


__all__ = [
    "bind_columns",
    "bind_rows",
    "data_range",
    "filter_strand",
    "filter_values",
    "matrix_info",
    "relabel_matrix",
    "sort_matrix",
    "subset_matrix",
]
_register_public_verbose(globals())
