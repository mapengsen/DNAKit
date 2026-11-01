"""methylpy wrapper for the ``filter-allc`` operation: coverage, context and chromosome filtering of allc tables.

``methylpy_filter_allc_files`` runs methylpy's own ``filter_allc_files`` on local allc files, sequentially
(``num_procs`` is fixed to 1). A row is kept when its chromosome is in ``chroms`` (if given), its
mc class matches one of the ``mc_type`` contexts (if given; a list of 3-letter codes, as on the command
line) and ``min_cov <= total reads <= max_cov``. A bare string is refused: methylpy's expansion tests
``"C" in mc_type`` by substring when given a string, which keeps every cytosine context (reproduced:
``'CGN'`` keeps the CHH row of a three-row table, ``['CGN']`` keeps only the CG rows).

The SNP mismatch filters (``max_mismatch``, ``max_mismatch_frac``) are deliberately not offered.
Reproduced upstream defect: methylpy 1.4.7 ``filter_allc_file_worker`` builds ``map`` objects
(utilities.py:345-346) and then subscripts them; on Python 3 this raises ``TypeError``, which the
bare ``except: pass`` clauses swallow, so every row passes the mismatch filter unchanged.

Each call checks its inputs, refuses to overwrite existing outputs, removes partial outputs on failure
and turns methylpy's exits and low-level errors into ``ConfigurationError``. Packages are imported only
when the function runs.
"""

from __future__ import annotations

import contextlib
import gzip
import io
import re
from collections.abc import Sequence
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_METHYLPY_EXTRA = "bio-methylpy"
_MAX_FILES = 32
_MAX_CHROMS = 10_000
_MAX_COVERAGE = 1_000_000_000
_BUFFER_LINES = 100_000
_MC_CONTEXT = re.compile(r"[ACGTNHDBVRYKMSW]{3}")
_MAX_CONTEXTS = 16
_CHROM = re.compile(r"[^\t\r\n\x00]{1,255}")
_UPSTREAM_ERRORS = (
    ArithmeticError,
    EOFError,
    LookupError,
    NameError,
    OSError,
    ValueError,
    TypeError,
)

__all__ = ["methylpy_filter_allc_files"]


def _integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigurationError(f"{name} must be an integer.")
    if not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be between {minimum} and {maximum}.")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be true or false.")
    return value


def _input_files(values: object) -> list[Path]:
    if isinstance(values, (str, bytes, Path)) or not isinstance(values, Sequence):
        raise ConfigurationError("allc_files must be a list of file paths.")
    if not 1 <= len(values) <= _MAX_FILES:
        raise ConfigurationError(f"allc_files must hold 1 to {_MAX_FILES} paths.")
    paths: list[Path] = []
    for value in values:
        if not isinstance(value, (str, Path)) or not str(value).strip():
            raise ConfigurationError("each allc file must be a non-empty path.")
        path = Path(value).expanduser()
        if not path.is_file():
            raise ConfigurationError(f"allc file not found: {path}")
        paths.append(path)
    return paths


def _output_files(values: object, compress: bool, inputs: Sequence[Path]) -> tuple[list[Path], list[Path]]:
    """Return the names given to methylpy and the final names (methylpy appends .gz when compressing)."""

    if isinstance(values, (str, bytes, Path)) or not isinstance(values, Sequence):
        raise ConfigurationError("output_files must be a list of file paths.")
    if len(values) != len(inputs):
        raise ConfigurationError("output_files must have one path per allc file.")
    given: list[Path] = []
    finals: list[Path] = []
    taken = {str(path.resolve()) for path in inputs}
    for value in values:
        if not isinstance(value, (str, Path)) or not str(value).strip():
            raise ConfigurationError("each output file must be a non-empty path.")
        path = Path(value).expanduser()
        if not path.parent.is_dir():
            raise ConfigurationError(f"output folder does not exist: {path.parent}")
        name = str(path)
        if not compress and name.endswith(".gz"):
            raise ConfigurationError("set compress_output=True for an output name ending in .gz.")
        final = name + ".gz" if compress and not name.endswith(".gz") else name
        final_path = Path(final)
        if final_path.exists():
            raise ConfigurationError(f"{final_path.name} already exists and will not be overwritten.")
        if str(final_path.resolve()) in taken:
            raise ConfigurationError(f"output {final_path.name} would overwrite an input or another output.")
        taken.add(str(final_path.resolve()))
        given.append(path)
        finals.append(final_path)
    return given, finals


def _contexts(values: object) -> list[str] | None:
    """Cytosine contexts in the form the methylpy command line passes (a list of 3-letter codes)."""

    if values is None:
        return None
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ConfigurationError(
            "mc_type must be a list of 3-letter contexts such as ['CGN'] or ['CHH', 'CHN']; "
            "a bare string makes methylpy match substrings and keep every cytosine context."
        )
    if not 1 <= len(values) <= _MAX_CONTEXTS:
        raise ConfigurationError(f"mc_type must hold 1 to {_MAX_CONTEXTS} contexts.")
    contexts: list[str] = []
    for value in values:
        if not isinstance(value, str) or not _MC_CONTEXT.fullmatch(value):
            raise ConfigurationError(
                f"invalid cytosine context {value!r}: use 3 IUPAC letters, for example CGN, CHN or CHH."
            )
        contexts.append(value)
    return contexts


def _chroms(values: object) -> list[str] | None:
    if values is None:
        return None
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ConfigurationError("chroms must be a list of chromosome names.")
    if not 1 <= len(values) <= _MAX_CHROMS:
        raise ConfigurationError(f"chroms must hold 1 to {_MAX_CHROMS} names.")
    names: list[str] = []
    for value in values:
        if not isinstance(value, str) or not _CHROM.fullmatch(value):
            raise ConfigurationError(f"invalid chromosome name: {value!r}")
        names.append(value)
    return names


def _row_count(path: Path) -> int:
    # Same test methylpy uses in open_allc_file: a name ending in "gz" is read as gzip.
    opener = gzip.open if path.name.endswith("gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip())


def _discard(paths: Sequence[Path]) -> None:
    for path in paths:
        with contextlib.suppress(OSError):
            path.unlink()


def _upstream_call(function: Any, *args: Any, **kwargs: Any) -> Any:
    """Run a methylpy function that prints progress, exits on bad input and fails on bad tables."""

    captured = io.StringIO()
    try:
        with contextlib.redirect_stdout(captured), contextlib.redirect_stderr(captured):
            return function(*args, **kwargs)
    except SystemExit as exc:
        lines = captured.getvalue().strip().splitlines()
        detail = lines[-1] if lines else str(exc.code)
        raise ConfigurationError(f"methylpy rejected the input: {detail}") from exc
    except _UPSTREAM_ERRORS as exc:
        raise ConfigurationError(
            f"methylpy could not process the input ({type(exc).__name__}: {exc})."
        ) from exc


def _provenance() -> dict[str, Any]:
    try:
        installed = version("methylpy")
    except PackageNotFoundError:
        installed = "unknown"
    return {
        "package": f"methylpy {installed}",
        "method": (
            "a row is kept when its chromosome is listed in chroms (if given), its mc class matches one "
            "of the mc_type contexts (if given; H and N are wildcards), and min_cov <= total reads "
            "(column 6) <= max_cov"
        ),
        "omitted": (
            "SNP mismatch filters (max_mismatch, max_mismatch_frac) are not offered: methylpy 1.4.7 "
            "applies them to map objects that raise TypeError on Python 3 and the bare except hides "
            "it, so every row passes them (reproduced)"
        ),
    }


def methylpy_filter_allc_files(
    allc_files: Sequence[str | Path],
    output_files: Sequence[str | Path],
    *,
    mc_type: Sequence[str] | None = None,
    chroms: Sequence[str] | None = None,
    min_cov: int = 0,
    max_cov: int | None = None,
    compress_output: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter allc tables by mc context, chromosome and read coverage (methylpy filter-allc).

    ``allc_files`` are read in order and each filtered table is written to the matching entry of
    ``output_files``. ``mc_type`` is a list of 3-letter cytosine contexts with the H and N wildcards,
    as on the methylpy command line (for example ``["CGN"]`` keeps CG sites, ``["CHH"]`` keeps CHH
    sites); ``chroms`` keeps only the named chromosomes; ``min_cov`` and ``max_cov`` bound the total
    read count (inclusive, ``max_cov=None`` means no upper bound). Defaults match methylpy: no context
    or chromosome filter, ``min_cov=0``, ``compress_output=False``. The SNP mismatch filters of methylpy
    are not offered (see the module docstring). ``verbose`` is accepted for the public-function
    contract; methylpy's progress lines are captured and not printed.
    """

    inputs = _input_files(allc_files)
    compress = _flag(compress_output, "compress_output")
    given, finals = _output_files(output_files, compress, inputs)
    contexts = _contexts(mc_type)
    chrom_list = _chroms(chroms)
    low = _integer(min_cov, "min_cov", 0, _MAX_COVERAGE)
    high = None if max_cov is None else _integer(max_cov, "max_cov", low, _MAX_COVERAGE)
    counts_in = [_row_count(path) for path in inputs]
    utilities = require("methylpy.utilities", _METHYLPY_EXTRA)
    try:
        _upstream_call(
            utilities.filter_allc_files,
            [str(path) for path in inputs],
            [str(path) for path in given],
            num_procs=1,
            mc_type=contexts,
            chroms=chrom_list,
            compress_output=compress,
            min_cov=low,
            max_cov=high,
            max_mismatch=None,
            max_mismatch_frac=None,
            buffer_line_number=_BUFFER_LINES,
        )
    except BaseException:
        _discard(finals)
        raise
    counts_out = [_row_count(path) for path in finals]
    return {
        "status": "completed",
        "operation": "methylpy_filter_allc_files",
        "inputs": [str(path) for path in inputs],
        "outputs": [
            {"input": str(src), "output": str(dst), "records_in": before, "records_out": after}
            for src, dst, before, after in zip(inputs, finals, counts_in, counts_out)
        ],
        "filters": {"mc_type": contexts, "chroms": chrom_list, "min_cov": low, "max_cov": high},
        "compressed": compress,
        "provenance": _provenance(),
    }


_register_public_verbose(globals())
