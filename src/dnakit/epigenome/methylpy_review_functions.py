"""methylpy wrappers for DMR finding steps, allc table splitting, merging and binomial tests.

These functions complete the methylpy coverage of ``methylpy_functions``. They call methylpy's
pure-Python table helpers on local files; the bundled ``run_rms_tests.out`` permutation binary and
the samtools, bowtie2, tabix and bgzip pipelines are not used. Every call is sequential
(``num_procs`` is fixed to 1) so that results do not depend on process scheduling.

methylpy shells out to ``sort``, ``wc`` and ``rm``, sorts some input files in place, leaves
temporary files behind on errors and fails with bare Python errors on malformed tables. Each
wrapper therefore checks its inputs first, refuses to overwrite existing outputs, keeps methylpy
away from the caller's input files, and turns methylpy's exits and low-level errors into
``ConfigurationError``. Packages are imported only when a function runs.
"""

from __future__ import annotations

import contextlib
import gzip
import io
import math
import os
import re
import shutil
import subprocess
import tempfile
from collections.abc import Iterator, Mapping, Sequence
from importlib.metadata import PackageNotFoundError, version
from itertools import groupby
from pathlib import Path
from typing import IO, Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.epigenome.methylpy_functions import index_allc_table
from dnakit.exceptions import ConfigurationError

_METHYLPY_EXTRA = "bio-methylpy"
_MAX_FILES = 32
_MAX_MERGE_FILES = 1_000
_MAX_CONTEXT_LENGTH = 6
_MAX_BUFFER_LINES = 10_000_000
_MAX_CHUNKS = 10_000
_MAX_SIMULATIONS = 1_000_000
_MAX_RESULT_ROWS = 1_000
_SAMPLE_NAME = re.compile(r"[A-Za-z0-9_.-]{1,64}")
_CHROM_NAME = re.compile(r"[A-Za-z0-9_.+=@:-]{1,255}")
_SHELL_SAFE_PATH = re.compile(r"[A-Za-z0-9_./+=:,@%-]+")
_SORT_MEMORY = re.compile(r"[1-9][0-9]{0,5}[KMG]?")
_CONTEXT_PATTERN = re.compile(rf"[ACGTNHDBRYKMSW]{{1,{_MAX_CONTEXT_LENGTH}}}")
_CONTROL_REGION = re.compile(r"([A-Za-z0-9_.+=@]+)(?::([0-9]+)(?:-([0-9]+))?)?")
_M0_LINE = re.compile(r"m0 estimate for iteration \d+: (\S+)")
_FDR_LINE = re.compile(r"corresponds to an FDR of (\S+)")
_RMS_FIXED = ("chr", "pos", "strand", "mc_class", "pvalue")
_RMS_GROUPS = ("mc_", "h_", "frac_", "mc_residual_", "uc_residual_")
_RMS_TAIL = ("num_simulations_sig", "num_simulations_run")


def _integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _number(value: object, name: str, minimum: float, maximum: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a number in [{minimum}, {maximum}].")
    number = float(value)
    if not minimum <= number <= maximum:
        raise ConfigurationError(f"{name} must be a number in [{minimum}, {maximum}].")
    return number


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
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


def _directory(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a directory path.")
    path = Path(value).expanduser().resolve()
    if not path.is_dir():
        raise ConfigurationError(f"{name} is not an existing directory: {path}")
    return path


def _shell_safe(path: Path, name: str) -> Path:
    """methylpy passes some paths to ``sort`` and ``rm`` after splitting them on blanks."""

    if not _SHELL_SAFE_PATH.fullmatch(str(path)):
        raise ConfigurationError(
            f"{name} may only contain letters, digits and _ . / + = : , @ % - because methylpy "
            f"hands it to shell commands: {path}"
        )
    return path


def _sample(value: object, name: str) -> str:
    if not isinstance(value, str) or not _SAMPLE_NAME.fullmatch(value):
        raise ConfigurationError(f"{name} must be 1 to 64 letters, digits, '_', '.' or '-'.")
    return value


def _sample_list(values: object, name: str, count: int | None = None) -> list[str]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of sample names.")
    names = [_sample(item, f"{name} entry") for item in values]
    if len(set(names)) != len(names):
        raise ConfigurationError(f"{name} must not list the same sample twice.")
    if count is not None and len(names) != count:
        raise ConfigurationError(f"{name} must list {count} names; it lists {len(names)}.")
    return names


def _path_list(values: object, name: str, minimum: int, maximum: int) -> list[Path]:
    if isinstance(values, (str, Path)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list of file paths.")
    if not minimum <= len(values) <= maximum:
        raise ConfigurationError(f"{name} must list {minimum} to {maximum} files.")
    paths = [_path(item, name) for item in values]
    if len(set(paths)) != len(paths):
        raise ConfigurationError(f"{name} must not list the same file twice.")
    return paths


def _contexts(value: object, name: str) -> list[str]:
    items = [value] if isinstance(value, str) else value
    if isinstance(items, (str, bytes)) or not isinstance(items, Sequence) or not items:
        raise ConfigurationError(f"{name} must be a context string or a non-empty list of them.")
    patterns = []
    for item in items:
        if not isinstance(item, str) or not _CONTEXT_PATTERN.fullmatch(item):
            raise ConfigurationError(
                f"{name} entries must be 1 to {_MAX_CONTEXT_LENGTH} IUPAC nucleotide codes, "
                "for example CGN."
            )
        patterns.append(item)
    return patterns


def _sort_memory(value: object) -> str:
    if not isinstance(value, str) or not _SORT_MEMORY.fullmatch(value):
        raise ConfigurationError("sort_mem must look like 500M, 2G or 64K.")
    return value


def _plain(value: Any) -> Any:
    if hasattr(value, "item"):
        value = value.item()
    return value


def _provenance(statistic: str, **notes: str) -> dict[str, Any]:
    try:
        installed = version("methylpy")
    except PackageNotFoundError:
        installed = "unknown"
    return {"package": f"methylpy {installed}", "method": statistic, **notes}


_UPSTREAM_ERRORS = (
    ArithmeticError,
    EOFError,
    LookupError,
    NameError,
    OSError,
    ValueError,
    subprocess.CalledProcessError,
)


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


def _upstream_text(function: Any, *args: Any, **kwargs: Any) -> tuple[Any, str]:
    """Like ``_upstream_call`` but also return what the function printed."""

    buffer = io.StringIO()

    def run() -> Any:
        with contextlib.redirect_stdout(buffer):
            return function(*args, **kwargs)

    value = _upstream_call(run)
    return value, buffer.getvalue()


@contextlib.contextmanager
def _c_locale() -> Iterator[None]:
    """Make the ``sort`` commands run by methylpy byte-ordered, whatever the user's locale."""

    saved = os.environ.get("LC_ALL")
    os.environ["LC_ALL"] = "C"
    try:
        yield
    finally:
        if saved is None:
            os.environ.pop("LC_ALL", None)
        else:
            os.environ["LC_ALL"] = saved


def _is_float(text: str) -> bool:
    try:
        float(text)
    except ValueError:
        return False
    return True


def _with_gz(path: Path, compress: bool) -> Path:
    """methylpy appends .gz to an output name that lacks it when it compresses."""

    if not compress and path.name.endswith("gz"):
        raise ConfigurationError("output_path must not end with gz when compress_output is False.")
    return path.with_name(path.name + ".gz") if compress and path.name[-3:] != ".gz" else path


def _discard(paths: Sequence[Path]) -> None:
    for path in paths:
        with contextlib.suppress(OSError):
            path.unlink(missing_ok=True)


def _is_gzip(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def _open_allc_text(path: Path) -> IO[str]:
    """Open an allc table the way methylpy does: a name ending in ``gz`` means gzip."""

    if path.name.endswith("gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _scan_allc(
    path: Path,
    *,
    contiguous: bool,
    numeric: bool,
    safe_names: bool = False,
    min_fields: int = 6,
) -> tuple[list[str], dict[str, int]]:
    """Check an allc table and return its chromosomes in order of appearance with line counts."""

    order: list[str] = []
    counts: dict[str, int] = {}
    current: str | None = None
    try:
        with _open_allc_text(path) as handle:
            for number, line in enumerate(handle, start=1):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < min_fields or not fields[0]:
                    raise ConfigurationError(
                        f"{path.name} line {number} has fewer than {min_fields} tab-separated "
                        "fields."
                    )
                chrom = fields[0]
                if safe_names and not _CHROM_NAME.fullmatch(chrom):
                    raise ConfigurationError(
                        f"{path.name} line {number}: chromosome name {chrom!r} cannot be used in "
                        "a file name."
                    )
                if numeric:
                    try:
                        position, methylated, total = (
                            int(fields[1]),
                            int(fields[4]),
                            int(fields[5]),
                        )
                    except ValueError as exc:
                        raise ConfigurationError(
                            f"{path.name} line {number} needs integer position, methylated and "
                            "total counts (a header line is not allowed)."
                        ) from exc
                    if position < 1 or methylated < 0 or total < methylated:
                        raise ConfigurationError(
                            f"{path.name} line {number} needs position >= 1 and "
                            "0 <= methylated <= total."
                        )
                if chrom != current:
                    if contiguous and chrom in counts:
                        raise ConfigurationError(
                            f"{path.name} line {number}: chromosome {chrom} appears in more than "
                            "one block; sort the allc table by chromosome first."
                        )
                    if chrom not in counts:
                        order.append(chrom)
                    current = chrom
                counts[chrom] = counts.get(chrom, 0) + 1
    except (OSError, EOFError, UnicodeDecodeError) as exc:
        raise ConfigurationError(f"{path.name} could not be read as an allc table: {exc}") from exc
    if not counts:
        raise ConfigurationError(f"{path.name} is empty.")
    return order, counts


def _tsv_rows(path: Path, *, skip: int = 1) -> Iterator[list[str]]:
    with path.open("r", encoding="utf-8") as handle:
        for number, line in enumerate(handle):
            if number >= skip:
                yield line.rstrip("\n").split("\t")


def _optional_float(text: str) -> float | None:
    return None if text == "NA" else float(text)


def _optional_int(text: str) -> int | None:
    return None if text == "NA" else int(text)


def _line_count(path: Path) -> int:
    with path.open("r", encoding="utf-8") as handle:
        return sum(1 for _ in handle)


# RMS results tables -----------------------------------------------------------------------


def _rms_samples(path: Path) -> list[str]:
    """Sample names from the header of a gzip-compressed RMS results table."""

    if not _is_gzip(path):
        raise ConfigurationError(f"{path.name} must be a gzip-compressed RMS results table.")
    try:
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            header = handle.readline().rstrip("\n").split("\t")
    except (OSError, EOFError, UnicodeDecodeError) as exc:
        raise ConfigurationError(f"{path.name} could not be read as gzip text: {exc}") from exc
    count, remainder = divmod(len(header) - len(_RMS_FIXED) - len(_RMS_TAIL), len(_RMS_GROUPS))
    if count < 1 or remainder != 0:
        raise ConfigurationError(f"{path.name} does not have the header of an RMS results table.")
    offset = len(_RMS_FIXED)
    samples = [name[3:] for name in header[offset : offset + count]]
    expected = list(_RMS_FIXED)
    for prefix in _RMS_GROUPS:
        expected.extend(prefix + sample for sample in samples)
    expected.extend(_RMS_TAIL)
    if header != expected:
        raise ConfigurationError(
            f"{path.name} does not have the header of an RMS results table: expected columns "
            f"{', '.join(_RMS_FIXED)}, then mc_/h_/frac_/mc_residual_/uc_residual_ per sample, "
            f"then {', '.join(_RMS_TAIL)}."
        )
    return samples


def _scan_rms(
    path: Path,
    *,
    less_than: bool = False,
    residuals: str = "number",
    simulations: tuple[int, int] | None = None,
    null_above: float | None = None,
) -> tuple[list[str], int, int]:
    """Check every row of an RMS table; return the samples, the row count and the null-site count.

    ``residuals`` is "number" (a number is required), "missing" (a number or NA) or "skip".
    ``simulations`` is ``(num_sims, num_sig_tests)`` for tables whose last two columns are checked.
    """

    samples = _rms_samples(path)
    count = len(samples)
    width = len(_RMS_FIXED) + len(_RMS_GROUPS) * count + len(_RMS_TAIL)
    first = len(_RMS_FIXED) + 3 * count
    rows = 0
    null_sites = 0
    try:
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            handle.readline()
            for number, line in enumerate(handle, start=2):
                fields = line.rstrip("\n").split("\t")
                where = f"{path.name} line {number}"
                if len(fields) != width:
                    raise ConfigurationError(
                        f"{where} has {len(fields)} fields; the header lists {width}."
                    )
                try:
                    int(fields[1])
                except ValueError as exc:
                    raise ConfigurationError(f"{where} has a non-integer position.") from exc
                raw = fields[4]
                if raw.startswith("<"):
                    if not less_than:
                        raise ConfigurationError(f"{where} has a P-value starting with '<'.")
                    p_value = 0.0
                else:
                    try:
                        p_value = float(raw)
                    except ValueError as exc:
                        raise ConfigurationError(f"{where} has a non-numeric P-value.") from exc
                    if not 0.0 <= p_value <= 1.0:
                        raise ConfigurationError(f"{where} has a P-value outside [0, 1].")
                if residuals != "skip":
                    for item in fields[first : first + 2 * count]:
                        if residuals == "missing" and item == "NA":
                            continue
                        try:
                            value = float(item)
                        except ValueError as exc:
                            raise ConfigurationError(
                                f"{where} has a non-numeric residual."
                            ) from exc
                        if not math.isfinite(value):
                            raise ConfigurationError(f"{where} has a non-finite residual.")
                if simulations is not None:
                    num_sims, num_sig_tests = simulations
                    try:
                        sig, run = int(fields[-2]), int(fields[-1])
                    except ValueError as exc:
                        raise ConfigurationError(
                            f"{where} needs integer simulation counts in its last two columns."
                        ) from exc
                    if not (
                        (run == num_sims and 1 <= sig <= num_sig_tests)
                        or (sig == num_sig_tests and num_sig_tests <= run < num_sims)
                    ):
                        raise ConfigurationError(
                            f"{where} reports {sig} of {run} simulations, which is impossible for "
                            f"num_sims={num_sims} and num_sig_tests={num_sig_tests}."
                        )
                if null_above is not None and p_value > null_above:
                    null_sites += 1
                rows += 1
    except (OSError, EOFError, UnicodeDecodeError) as exc:
        raise ConfigurationError(f"{path.name} could not be read: {exc}") from exc
    return samples, rows, null_sites


def _simulation_settings(num_sims: object, num_sig_tests: object) -> tuple[int, int]:
    sims = _integer(num_sims, "num_sims", 2, _MAX_SIMULATIONS)
    significant = _integer(num_sig_tests, "num_sig_tests", 1, sims)
    return sims, significant


def _fdr_cutoff(
    rms: Path,
    num_sims: int,
    num_sig_tests: int,
    target_fdr: float,
    max_iterations: int,
    convergence_diff: float,
) -> dict[str, Any]:
    _, rows, _ = _scan_rms(
        rms, residuals="skip", simulations=(num_sims, num_sig_tests), less_than=False
    )
    if rows == 0:
        raise ConfigurationError(f"{rms.name} contains no tests.")
    dmr = require("methylpy.DMRfind", _METHYLPY_EXTRA)
    cutoff, text = _upstream_text(
        dmr.histogram_correction_DMRfind,
        str(rms),
        num_sims,
        num_sig_tests,
        target_fdr=target_fdr,
        max_iterations=max_iterations,
        convergence_diff=convergence_diff,
    )
    if "did not converge" in text:
        outcome = "not_converged"
    elif "No test result found" in text:
        outcome = "no_significant_result"
    else:
        outcome = "converged"
    m0_values = _M0_LINE.findall(text)
    fdr_values = _FDR_LINE.findall(text)
    return {
        "total_tests": rows,
        "p_value_cutoff": float(cutoff),
        "outcome": outcome,
        "estimated_fdr": float(fdr_values[-1]) if fdr_values else None,
        "null_tests_estimate": float(m0_values[-1]) if m0_values else None,
    }


def _residual_cutoff(rms: Path, pvalue_cutoff: float, fraction: float) -> dict[str, Any]:
    _, _, null_sites = _scan_rms(rms, null_above=pvalue_cutoff)
    dmr = require("methylpy.DMRfind", _METHYLPY_EXTRA)
    value = _upstream_call(dmr.get_resid_cutoff, fraction, pvalue_cutoff, str(rms))
    return {
        "residual_cutoff": float(_plain(value)),
        "null_sites": null_sites,
        "default_used": null_sites == 0,
    }


def _check_collapse_samples(
    rms_samples: list[str],
    collapse_samples: object,
    sample_category: object,
    min_cluster: object,
) -> tuple[list[str] | None, list[int] | None, int]:
    clusters = _integer(min_cluster, "min_cluster", 0, 1_000_000)
    chosen: list[str] | None = None
    if collapse_samples is not None:
        chosen = _sample_list(collapse_samples, "collapse_samples")
        unknown = [name for name in chosen if name not in rms_samples]
        if unknown:
            raise ConfigurationError(
                f"collapse_samples {unknown} are not in the table, which has {rms_samples}."
            )
    categories: list[int] | None = None
    if sample_category is not None:
        if chosen is None:
            raise ConfigurationError("sample_category needs collapse_samples.")
        if isinstance(sample_category, (str, bytes)) or not isinstance(sample_category, Sequence):
            raise ConfigurationError("sample_category must be a list of integers.")
        if len(sample_category) != len(chosen):
            raise ConfigurationError("sample_category must have one entry per collapse sample.")
        categories = [
            _integer(item, "sample_category entry", 0, 1_000_000) for item in sample_category
        ]
    return chosen, categories, clusters


def _read_regions(table: Path) -> list[dict[str, Any]]:
    regions = []
    for fields in _tsv_rows(table):
        regions.append(
            {
                "chrom": fields[0],
                "start": int(fields[1]),
                "end": int(fields[2]),
                "dms": int(fields[3]),
                "hypermethylated_samples": [name for name in fields[4].split(",") if name],
                "hypomethylated_samples": [name for name in fields[5].split(",") if name],
            }
        )
    return regions


def _collapse_outputs(table: Path) -> tuple[Path, Path]:
    return Path(str(table) + ".DMS.bed"), Path(str(table) + ".DMR.bed")


def _collapse_result(table: Path, count: int) -> dict[str, Any]:
    dms_bed, dmr_bed = _collapse_outputs(table)
    regions = _read_regions(table)
    return {
        "dmr_count": count,
        "dms_count": _line_count(dms_bed) if dms_bed.is_file() else 0,
        "regions": regions[:_MAX_RESULT_ROWS],
        "regions_truncated": len(regions) > _MAX_RESULT_ROWS,
        "files": {
            "table": str(table),
            "dms_bed": str(dms_bed) if dms_bed.is_file() else None,
            "dmr_bed": str(dmr_bed) if dmr_bed.is_file() else None,
        },
    }


def estimate_dmr_fdr_cutoff(
    rms_path: str | Path,
    *,
    num_sims: int = 3000,
    num_sig_tests: int = 100,
    target_fdr: float = 0.01,
    max_iterations: int = 1000,
    convergence_diff: float = 0.01,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find the permutation P-value cutoff of an RMS results table at a target FDR."""

    rms = _path(rms_path, "rms_path")
    sims, significant = _simulation_settings(num_sims, num_sig_tests)
    fdr = _number(target_fdr, "target_fdr", 0.0, 1.0)
    if fdr <= 0.0:
        raise ConfigurationError("target_fdr must be greater than zero.")
    iterations = _integer(max_iterations, "max_iterations", 1, 1_000_000)
    tolerance = _number(convergence_diff, "convergence_diff", 0.0, 1e12)
    found = _fdr_cutoff(rms, sims, significant, fdr, iterations, tolerance)
    return {
        "status": "completed",
        "operation": "estimate_dmr_fdr_cutoff",
        "rms_table": str(rms),
        "num_sims": sims,
        "num_sig_tests": significant,
        "target_fdr": fdr,
        **found,
        "provenance": _provenance(
            "histogram method (Iowa State preprint 2009-04): the null-test count is re-estimated "
            "until it changes by no more than convergence_diff, then the largest cutoff with "
            "estimated FDR at or below target_fdr is chosen",
            note="outcome is not_converged or no_significant_result when methylpy falls back to "
            "1/num_sims; the fallback value is returned as the cutoff",
        ),
    }


def dmr_residual_cutoff(
    rms_path: str | Path,
    pvalue_cutoff: float,
    *,
    resid_cutoff: float = 0.01,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find the residual cutoff of non-significant sites in an RMS results table."""

    rms = _path(rms_path, "rms_path")
    level = _number(pvalue_cutoff, "pvalue_cutoff", 0.0, 1.0)
    fraction = _number(resid_cutoff, "resid_cutoff", 0.0, 1.0)
    found = _residual_cutoff(rms, level, fraction)
    return {
        "status": "completed",
        "operation": "dmr_residual_cutoff",
        "rms_table": str(rms),
        "pvalue_cutoff": level,
        "resid_cutoff": fraction,
        **found,
        "provenance": _provenance(
            "(1 - resid_cutoff) percentile of the mc and uc residuals of sites whose P-value is "
            "above pvalue_cutoff",
            note="methylpy includes the -1 placeholder of samples without data in the percentile; "
            "it returns 2 when no site is above the P-value cutoff",
        ),
    }


def collapse_dmr_windows(
    rms_path: str | Path,
    output_path: str | Path,
    *,
    max_dist: int = 100,
    resid_cutoff: float = 0.0,
    sig_cutoff: float = 0.01,
    min_num_dms: int = 0,
    collapse_samples: Sequence[str] | None = None,
    sample_category: Sequence[int] | None = None,
    min_cluster: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Collapse significant differentially methylated sites of an RMS table into regions."""

    rms = _path(rms_path, "rms_path")
    target = _new_path(output_path, "output_path")
    dms_bed, dmr_bed = _collapse_outputs(target)
    for extra in (dms_bed, dmr_bed):
        if extra.exists():
            raise ConfigurationError(f"{extra.name} already exists and will not be overwritten.")
    distance = _integer(max_dist, "max_dist", 0, 1_000_000_000)
    residual = _number(resid_cutoff, "resid_cutoff", 0.0, 1e6)
    level = _number(sig_cutoff, "sig_cutoff", 0.0, 1.0)
    minimum_sites = _integer(min_num_dms, "min_num_dms", 0, 1_000_000)
    samples, _, _ = _scan_rms(rms, less_than=True, residuals="missing")
    chosen, categories, clusters = _check_collapse_samples(
        samples, collapse_samples, sample_category, min_cluster
    )
    dmr = require("methylpy.DMRfind", _METHYLPY_EXTRA)
    try:
        count = _upstream_call(
            dmr.collapse_dmr_windows,
            str(rms),
            str(target),
            column=4,
            max_dist=distance,
            resid_cutoff=residual,
            sig_cutoff=level,
            min_num_dms=minimum_sites,
            collapse_samples=chosen if chosen is not None else False,
            sample_category=categories if categories is not None else False,
            min_cluster=clusters,
        )
    except BaseException:
        _discard([target, dms_bed, dmr_bed])
        raise
    return {
        "status": "completed",
        "operation": "collapse_dmr_windows",
        "rms_table": str(rms),
        "samples": samples,
        "max_dist": distance,
        "resid_cutoff": residual,
        "sig_cutoff": level,
        **_collapse_result(target, int(count or 0)),
        "provenance": _provenance(
            "consecutive significant sites within max_dist whose residual direction agrees form "
            "one region; the chromosome gets a 'chr' prefix in the BED files",
            note="the first significant site of the table never contributes sample lists of its "
            "own; the BED files are not written when no site is significant",
        ),
    }


def merge_dms_to_dmrs(
    rms_path: str | Path,
    output_path: str | Path,
    *,
    num_sims: int = 3000,
    num_sig_tests: int = 100,
    sig_cutoff: float = 0.01,
    dmr_max_dist: int = 250,
    min_num_dms: int = 0,
    resid_cutoff: float = 0.01,
    collapse_samples: Sequence[str] | None = None,
    sample_category: Sequence[int] | None = None,
    min_cluster: int = 0,
    max_iterations: int = 1000,
    convergence_diff: float = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Turn the sites of an RMS results table into DMRs after FDR and residual cutoffs."""

    rms = _path(rms_path, "rms_path")
    target = _new_path(output_path, "output_path")
    dms_bed, dmr_bed = _collapse_outputs(target)
    for extra in (dms_bed, dmr_bed):
        if extra.exists():
            raise ConfigurationError(f"{extra.name} already exists and will not be overwritten.")
    sims, significant = _simulation_settings(num_sims, num_sig_tests)
    fdr = _number(sig_cutoff, "sig_cutoff", 0.0, 1.0)
    if fdr <= 0.0:
        raise ConfigurationError("sig_cutoff must be greater than zero.")
    distance = _integer(dmr_max_dist, "dmr_max_dist", 0, 1_000_000_000)
    minimum_sites = _integer(min_num_dms, "min_num_dms", 0, 1_000_000)
    fraction = _number(resid_cutoff, "resid_cutoff", 0.0, 1.0)
    iterations = _integer(max_iterations, "max_iterations", 1, 1_000_000)
    tolerance = _number(convergence_diff, "convergence_diff", 0.0, 1e12)
    samples, rows, _ = _scan_rms(
        rms, residuals="number", simulations=(sims, significant), less_than=False
    )
    if rows == 0:
        raise ConfigurationError(f"{rms.name} contains no tests.")
    chosen, categories, clusters = _check_collapse_samples(
        samples, collapse_samples, sample_category, min_cluster
    )
    dmr = require("methylpy.DMRfind", _METHYLPY_EXTRA)
    try:
        _upstream_call(
            dmr.merge_DMS_to_DMR,
            str(rms),
            str(target),
            collapse_samples=chosen if chosen is not None else False,
            sample_category=categories if categories is not None else False,
            min_cluster=clusters,
            sig_cutoff=fdr,
            dmr_max_dist=distance,
            min_num_dms=minimum_sites,
            resid_cutoff=fraction,
            num_sims=sims,
            num_sig_tests=significant,
            max_iterations=iterations,
            convergence_diff=tolerance,
        )
        histogram = _fdr_cutoff(rms, sims, significant, fdr, iterations, tolerance)
        residual = _residual_cutoff(rms, histogram["p_value_cutoff"], fraction)
    except BaseException:
        _discard([target, dms_bed, dmr_bed])
        raise
    count = len(_read_regions(target))
    return {
        "status": "completed",
        "operation": "merge_dms_to_dmrs",
        "rms_table": str(rms),
        "samples": samples,
        "sig_cutoff": fdr,
        "dmr_max_dist": distance,
        "p_value_cutoff": histogram["p_value_cutoff"],
        "fdr_outcome": histogram["outcome"],
        "residual_cutoff": residual["residual_cutoff"],
        **_collapse_result(target, count),
        "provenance": _provenance(
            "FDR cutoff by the histogram method, residual cutoff from non-significant sites, then "
            "window collapsing",
            note="the two cutoffs are recomputed from the same table with the same options after "
            "the merge so that they can be reported",
        ),
    }


# allc tables joined to DMR tables ---------------------------------------------------------


def _dmr_input(dmr_path: str | Path, input_no_header: bool) -> tuple[Path, set[str]]:
    source = _path(dmr_path, "dmr_path")
    skip = 0 if _flag(input_no_header, "input_no_header") else 1
    chroms: set[str] = set()
    previous: tuple[str, int] | None = None
    try:
        for number, fields in enumerate(_tsv_rows(source, skip=skip), start=skip + 1):
            if len(fields) < 3:
                raise ConfigurationError(f"{source.name} line {number} has fewer than 3 fields.")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ConfigurationError(
                    f"{source.name} line {number} needs integer start and end positions."
                ) from exc
            if start < 0 or end < start:
                raise ConfigurationError(f"{source.name} line {number} needs 0 <= start <= end.")
            if previous is not None and previous[0] == fields[0] and start < previous[1]:
                raise ConfigurationError(
                    f"{source.name} line {number}: regions must be sorted by start within a "
                    "chromosome."
                )
            previous = (fields[0], start)
            chroms.add(fields[0])
    except (OSError, UnicodeDecodeError) as exc:
        raise ConfigurationError(f"{source.name} could not be read: {exc}") from exc
    if not chroms:
        raise ConfigurationError(f"{source.name} contains no regions.")
    return source, chroms


def _allc_inputs(allc_paths: object, samples: object) -> tuple[list[Path], list[str]]:
    paths = _path_list(allc_paths, "allc_paths", 1, _MAX_FILES)
    return paths, _sample_list(samples, "samples", count=len(paths))


def _require_chromosomes(paths: Sequence[Path], names: Sequence[str], wanted: set[str]) -> None:
    """Every region's chromosome must be in each allc index (this writes missing indexes)."""

    for path, name in zip(paths, names, strict=True):
        known = {item["chrom"] for item in index_allc_table(path, verbose=False)["chromosomes"]}
        missing = sorted(wanted - known)
        if missing:
            raise ConfigurationError(
                f"Regions on {missing} have no records in the allc table of sample {name}."
            )


def _joined_output(output_path: object, suffixes: Sequence[str], names: Sequence[str]) -> Path:
    target = _new_path(output_path, "output_path")
    if not target.name.endswith(".tsv") or str(target).count(".tsv") != 1:
        raise ConfigurationError(
            "output_path must end with .tsv and contain '.tsv' nowhere else, because methylpy "
            "derives its temporary file names from it."
        )
    _shell_safe(target, "output_path")
    related = [Path(str(target) + extra) for extra in suffixes]
    temporary = [
        Path(str(target).replace(".tsv", "") + f"_{name}_temp_{kind}.tsv")
        for name in names
        for kind in ("methylation_levels", "c_info")
    ]
    for item in [*related, *temporary]:
        if item.exists():
            raise ConfigurationError(f"{item.name} already exists and will not be overwritten.")
    return target


def _temporary_files(target: Path, names: Sequence[str], kind: str) -> list[Path]:
    stem = str(target).replace(".tsv", "")
    return [Path(f"{stem}_{name}_temp_{kind}.tsv") for name in names]


def add_methylation_levels_to_dmrs(
    dmr_path: str | Path,
    allc_paths: Sequence[str | Path],
    samples: Sequence[str],
    output_path: str | Path,
    *,
    mc_type: str | Sequence[str] = "CGN",
    min_cov: int = 0,
    max_cov: int | None = None,
    buffer_line_number: int = 100_000,
    input_no_header: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add the methylation level of every sample to each region of a DMR table."""

    source, chroms = _dmr_input(dmr_path, input_no_header)
    paths, names = _allc_inputs(allc_paths, samples)
    patterns = _contexts(mc_type, "mc_type")
    lowest = _integer(min_cov, "min_cov", 0, 1_000_000_000)
    highest = None if max_cov is None else _integer(max_cov, "max_cov", lowest, 1_000_000_000)
    buffer_lines = _integer(buffer_line_number, "buffer_line_number", 1, _MAX_BUFFER_LINES)
    target = _joined_output(output_path, [], names)
    _require_chromosomes(paths, names, chroms)
    dmr = require("methylpy.DMRfind", _METHYLPY_EXTRA)
    temporary = _temporary_files(target, names, "methylation_levels")
    try:
        _upstream_call(
            dmr.get_methylation_levels_DMRfind,
            str(source),
            str(target),
            [str(path) for path in paths],
            samples=names,
            mc_type=patterns,
            num_procs=1,
            min_cov=lowest,
            max_cov=highest,
            buffer_line_number=buffer_lines,
            input_no_header=input_no_header,
        )
    except BaseException:
        _discard([target, *temporary])
        raise
    _discard(temporary)
    rows = []
    for fields in _tsv_rows(target):
        values = fields[-len(names) :]
        rows.append(
            {
                "chrom": fields[0],
                "start": int(fields[1]),
                "end": int(fields[2]),
                "methylation_level": {
                    name: _optional_float(value) for name, value in zip(names, values, strict=True)
                },
            }
        )
    return {
        "status": "completed",
        "operation": "add_methylation_levels_to_dmrs",
        "dmr_table": str(source),
        "output": str(target),
        "samples": names,
        "mc_type": patterns,
        "dmr_count": len(rows),
        "rows": rows[:_MAX_RESULT_ROWS],
        "rows_truncated": len(rows) > _MAX_RESULT_ROWS,
        "provenance": _provenance(
            "methylated reads divided by total reads of the sites in each region whose context "
            "matches mc_type and whose coverage is within [min_cov, max_cov]; NA without sites",
            note="allc indexes (.idx) are created next to the allc files; a region on a "
            "chromosome missing from an allc table is rejected because methylpy would repeat "
            "the previous region's value",
        ),
    }


def add_cytosine_info_to_dmrs(
    dmr_path: str | Path,
    allc_paths: Sequence[str | Path],
    samples: Sequence[str],
    output_path: str | Path,
    *,
    mc_type: str | Sequence[str] = "CGN",
    min_cov: int = 0,
    max_cov: int | None = None,
    buffer_line_number: int = 100_000,
    input_no_header: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add methylation level, total base calls and covered sites per sample to each region."""

    source, chroms = _dmr_input(dmr_path, input_no_header)
    paths, names = _allc_inputs(allc_paths, samples)
    patterns = _contexts(mc_type, "mc_type")
    lowest = _integer(min_cov, "min_cov", 0, 1_000_000_000)
    highest = None if max_cov is None else _integer(max_cov, "max_cov", lowest, 1_000_000_000)
    buffer_lines = _integer(buffer_line_number, "buffer_line_number", 1, _MAX_BUFFER_LINES)
    target = _joined_output(output_path, [".total_basecalls.tsv", ".covered_sites.tsv"], names)
    _require_chromosomes(paths, names, chroms)
    extra = [Path(str(target) + ".total_basecalls.tsv"), Path(str(target) + ".covered_sites.tsv")]
    dmr = require("methylpy.DMRfind", _METHYLPY_EXTRA)
    temporary = _temporary_files(target, names, "c_info")
    try:
        _upstream_call(
            dmr.get_c_info_DMRfind,
            str(source),
            [str(path) for path in paths],
            names,
            str(target),
            mc_type=patterns,
            num_procs=1,
            min_cov=lowest,
            max_cov=highest,
            buffer_line_number=buffer_lines,
            input_no_header=input_no_header,
        )
    except BaseException:
        _discard([target, *extra, *temporary])
        raise
    _discard(temporary)
    levels = list(_tsv_rows(target))
    totals = list(_tsv_rows(extra[0]))
    covered = list(_tsv_rows(extra[1]))
    rows = []
    for level_row, total_row, covered_row in zip(levels, totals, covered, strict=True):
        count = len(names)
        rows.append(
            {
                "chrom": level_row[0],
                "start": int(level_row[1]),
                "end": int(level_row[2]),
                "samples": {
                    name: {
                        "methylation_level": _optional_float(level),
                        "total_basecalls": _optional_int(total),
                        "covered_sites": _optional_int(sites),
                    }
                    for name, level, total, sites in zip(
                        names,
                        level_row[-count:],
                        total_row[-count:],
                        covered_row[-count:],
                        strict=True,
                    )
                },
            }
        )
    return {
        "status": "completed",
        "operation": "add_cytosine_info_to_dmrs",
        "dmr_table": str(source),
        "output": str(target),
        "total_basecalls_output": str(extra[0]),
        "covered_sites_output": str(extra[1]),
        "samples": names,
        "mc_type": patterns,
        "dmr_count": len(rows),
        "rows": rows[:_MAX_RESULT_ROWS],
        "rows_truncated": len(rows) > _MAX_RESULT_ROWS,
        "provenance": _provenance(
            "per region: methylated reads over total reads, the total reads (base calls) and the "
            "number of sites with a matching context and coverage within [min_cov, max_cov]",
            note="allc indexes (.idx) are created next to the allc files; a region on a "
            "chromosome missing from an allc table is rejected",
        ),
    }


# allc tables: splitting, testing, filtering and merging -----------------------------------


def _split_names(directory: Path, sample: str, labels: Sequence[str], compress: bool) -> list[Path]:
    suffix = ".tsv.gz" if compress else ".tsv"
    return [directory / f"allc_{sample}_{label}{suffix}" for label in labels]


def _read_chunk(path: Path) -> tuple[int, list[str]]:
    chroms: set[str] = set()
    records = 0
    with _open_allc_text(path) as handle:
        for line in handle:
            records += 1
            chroms.add(line.split("\t", 1)[0])
    return records, sorted(chroms)


def split_allc_by_chromosome(
    allc_path: str | Path,
    sample: str,
    output_dir: str | Path,
    *,
    compress_output: bool = True,
    buffer_line_number: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write one allc file per chromosome, named ``allc_<sample>_<chrom>.tsv[.gz]``."""

    source = _path(allc_path, "allc_path")
    name = _sample(sample, "sample")
    directory = _directory(output_dir, "output_dir")
    compress = _flag(compress_output, "compress_output")
    buffer_lines = _integer(buffer_line_number, "buffer_line_number", 1, _MAX_BUFFER_LINES)
    order, counts = _scan_allc(source, contiguous=True, numeric=False, safe_names=True)
    outputs = _split_names(directory, name, order, compress)
    for item in outputs:
        if item.exists():
            raise ConfigurationError(f"{item.name} already exists and will not be overwritten.")
    call_mc = require("methylpy.call_mc_se", _METHYLPY_EXTRA)
    try:
        _upstream_call(
            call_mc.do_split_allc_file,
            str(source),
            name,
            str(directory),
            compress_output=compress,
            buffer_line_number=buffer_lines,
        )
    except BaseException:
        _discard(outputs)
        raise
    missing = [item.name for item in outputs if not item.is_file()]
    if missing:
        raise ConfigurationError(f"methylpy did not create {missing}.")
    return {
        "status": "completed",
        "operation": "split_allc_by_chromosome",
        "input": str(source),
        "sample": name,
        "compressed": compress,
        "records": sum(counts.values()),
        "files": [
            {"chrom": chrom, "path": str(item), "records": counts[chrom]}
            for chrom, item in sorted(zip(order, outputs, strict=True))
        ],
        "provenance": _provenance(
            "one output file per chromosome block of the input table",
            note="methylpy lists the names without .gz even when it writes .tsv.gz; the paths "
            "returned here are the files that exist. Chromosomes must be contiguous because "
            "methylpy overwrites the file of a chromosome that reappears",
        ),
    }


def split_allc_into_chunks(
    allc_path: str | Path,
    sample: str,
    num_chunks: int,
    output_dir: str | Path,
    *,
    compress_output: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Distribute the chromosomes of an allc table over at most ``num_chunks`` files."""

    source = _path(allc_path, "allc_path")
    name = _sample(sample, "sample")
    chunks = _integer(num_chunks, "num_chunks", 1, _MAX_CHUNKS)
    directory = _directory(output_dir, "output_dir")
    compress = _flag(compress_output, "compress_output")
    order, counts = _scan_allc(source, contiguous=False, numeric=False)
    labels = [str(index) for index in range(min(chunks, len(order)))]
    outputs = _split_names(directory, name, labels, compress)
    for item in outputs:
        if item.exists():
            raise ConfigurationError(f"{item.name} already exists and will not be overwritten.")
    call_mc = require("methylpy.call_mc_se", _METHYLPY_EXTRA)
    try:
        _upstream_call(
            call_mc.do_split_allc_file_chunk,
            str(source),
            name,
            chunks,
            str(directory),
            compress_output=compress,
        )
    except BaseException:
        _discard(outputs)
        raise
    written = [item for item in outputs if item.is_file()]
    if not written:
        raise ConfigurationError("methylpy wrote no chunk files.")
    details = []
    for index, item in enumerate(written):
        records, chroms = _read_chunk(item)
        details.append(
            {"index": index, "path": str(item), "records": records, "chromosomes": chroms}
        )
    return {
        "status": "completed",
        "operation": "split_allc_into_chunks",
        "input": str(source),
        "sample": name,
        "compressed": compress,
        "requested_chunks": chunks,
        "records": sum(counts.values()),
        "chunks": details,
        "provenance": _provenance(
            "chromosomes, smallest first, are packed into a chunk until it holds "
            "ceil(records / num_chunks) records; a chromosome is never split",
            note="fewer files than num_chunks are written when chromosomes are few or uneven; "
            "methylpy lists the names without .gz even when it writes .tsv.gz",
        ),
    }


def _context_counts(path: Path) -> list[dict[str, Any]]:
    counts: dict[str, list[int]] = {}
    with _open_allc_text(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            entry = counts.setdefault(fields[3], [0, 0])
            entry[0] += 1
            entry[1] += 1 if fields[6] == "1" else 0
    return [
        {"context": context, "sites": counts[context][0], "significant_sites": counts[context][1]}
        for context in sorted(counts)
    ]


def run_allc_binomial_test(
    allc_path: str | Path,
    sample: str,
    output_dir: str | Path,
    unmethylated_control: float | str,
    *,
    min_cov: int = 2,
    sig_cutoff: float = 0.01,
    sort_mem: str = "500M",
    compress_output: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call methylated sites with a binomial test and a per-context Benjamini-Hochberg cutoff."""

    source = _path(allc_path, "allc_path")
    name = _sample(sample, "sample")
    directory = _shell_safe(_directory(output_dir, "output_dir"), "output_dir")
    lowest = _integer(min_cov, "min_cov", 1, 1_000_000)
    level = _number(sig_cutoff, "sig_cutoff", 0.0, 1.0)
    if level <= 0.0:
        raise ConfigurationError("sig_cutoff must be greater than zero.")
    memory = _sort_memory(sort_mem)
    compress = _flag(compress_output, "compress_output")
    if isinstance(unmethylated_control, str):
        if not _CONTROL_REGION.fullmatch(unmethylated_control) or _is_float(unmethylated_control):
            raise ConfigurationError(
                "unmethylated_control must be a rate in [0, 1] or a region such as chrL, "
                "chrL:100 or chrL:100-2000."
            )
        control: float | str = unmethylated_control
    else:
        control = _number(unmethylated_control, "unmethylated_control", 0.0, 1.0)
    _scan_allc(source, contiguous=True, numeric=True)
    final = directory / (f"allc_{name}.tsv.gz" if compress else f"allc_{name}.tsv")
    scratch = directory / f"allc_{name}.tsv_binom_results.tsv"
    for item in (final, scratch):
        if item.exists():
            raise ConfigurationError(f"{item.name} already exists and will not be overwritten.")
    call_mc = require("methylpy.call_mc_se", _METHYLPY_EXTRA)
    # methylpy always removes a "chr" prefix from a control region; look the region up verbatim.
    rate = (
        _upstream_call(call_mc.calculate_non_conversion_rate, control, str(source), None, False)
        if isinstance(control, str)
        else control
    )
    try:
        with _c_locale():
            _upstream_call(
                call_mc.perform_binomial_test,
                str(source),
                name,
                str(directory),
                float(rate),
                min_cov=lowest,
                sig_cutoff=level,
                num_procs=1,
                sort_mem=memory,
                compress_output=compress,
            )
    except BaseException:
        _discard([final, scratch])
        raise
    contexts = _context_counts(final)
    return {
        "status": "completed",
        "operation": "run_allc_binomial_test",
        "input": str(source),
        "index": str(source) + ".idx",
        "output": str(final),
        "sample": name,
        "non_conversion_rate": float(rate),
        "min_cov": lowest,
        "sig_cutoff": level,
        "sites": sum(item["sites"] for item in contexts),
        "significant_sites": sum(item["significant_sites"] for item in contexts),
        "contexts": contexts,
        "provenance": _provenance(
            "one-sided binomial test of the methylated reads against the non-conversion rate, "
            "Benjamini-Hochberg cutoff per context; column 7 of the output marks significant "
            "sites with 1",
            note="sites without reads are dropped; the allc index is created next to the input; "
            "sort runs with LC_ALL=C so the order does not depend on the locale",
        ),
    }


def _context_key(value: object) -> str:
    if not isinstance(value, str) or not value or any(c in value for c in "\x00\r\n\t "):
        raise ConfigurationError("cutoffs keys must be context strings such as CGN.")
    return value


def filter_allc_by_pvalue_cutoffs(
    binomial_paths: Sequence[str | Path],
    output_path: str | Path,
    cutoffs: Mapping[str, float],
    *,
    compress_output: bool = True,
    sort_mem: str = "500M",
    verbose: bool = True,
) -> dict[str, Any]:
    """Flag the sites of binomial-test files whose P-value is within their context's cutoff."""

    paths = _path_list(binomial_paths, "binomial_paths", 1, _MAX_FILES)
    compress = _flag(compress_output, "compress_output")
    memory = _sort_memory(sort_mem)
    if not isinstance(cutoffs, Mapping) or not 1 <= len(cutoffs) <= 64:
        raise ConfigurationError("cutoffs must map 1 to 64 contexts to P-value cutoffs.")
    limits = {
        _context_key(key): _number(value, f"cutoff for {key}", 0.0, 1.0)
        for key, value in cutoffs.items()
    }
    given = _new_path(output_path, "output_path")
    final = _with_gz(given, compress)
    if final.exists():
        raise ConfigurationError(f"output_path already exists and will not be overwritten: {final}")
    for path in paths:
        _check_binomial_file(path, limits)
    call_mc = require("methylpy.call_mc_se", _METHYLPY_EXTRA)
    try:
        with tempfile.TemporaryDirectory(dir=final.parent, prefix="methylpy_filter_") as folder:
            work = _shell_safe(Path(folder), "the output directory")
            copies = []
            for index, path in enumerate(paths):
                copy = work / f"input_{index}.tsv"
                shutil.copyfile(path, copy)
                copies.append(str(copy))
            with _c_locale():
                _upstream_call(
                    call_mc.filter_files_by_pvalue_combined,
                    copies,
                    str(final),
                    limits,
                    1,
                    compress_output=compress,
                    sort_mem=memory,
                )
    except BaseException:
        _discard([final])
        raise
    contexts = _context_counts(final)
    return {
        "status": "completed",
        "operation": "filter_allc_by_pvalue_cutoffs",
        "inputs": [str(path) for path in paths],
        "output": str(final),
        "records": sum(item["sites"] for item in contexts),
        "significant_sites": sum(item["significant_sites"] for item in contexts),
        "contexts": contexts,
        "provenance": _provenance(
            "each input is sorted by chromosome and position, then written with column 7 set to 1 "
            "when its P-value is at most the cutoff of its context and 0 otherwise",
            note="methylpy sorts its inputs in place; copies are filtered so that the input files "
            "stay untouched. A P-value of 2.0 marks an untested site and is never significant",
        ),
    }


def _check_binomial_file(path: Path, limits: Mapping[str, float]) -> None:
    try:
        with path.open("r", encoding="utf-8") as handle:
            for number, line in enumerate(handle, start=1):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 7 or not fields[0]:
                    raise ConfigurationError(f"{path.name} line {number} has fewer than 7 fields.")
                try:
                    int(fields[1])
                    p_value = float(fields[6])
                except ValueError as exc:
                    raise ConfigurationError(
                        f"{path.name} line {number} needs an integer position and a numeric "
                        "P-value in column 7."
                    ) from exc
                if not p_value >= 0.0:
                    raise ConfigurationError(f"{path.name} line {number} has a negative P-value.")
                if fields[6] != "2.0" and fields[3] not in limits:
                    raise ConfigurationError(
                        f"{path.name} line {number}: no cutoff was given for context {fields[3]}."
                    )
    except (OSError, UnicodeDecodeError) as exc:
        raise ConfigurationError(f"{path.name} could not be read: {exc}") from exc


def _chromosome_of(line: str) -> str:
    return line.split("\t", 1)[0]


def _sort_allc_chromosomes(path: Path, utilities: Any) -> bool:
    """Put the chromosome blocks of a merged allc table in sorted order; True when it changed."""

    with _open_allc_text(path) as handle:
        blocks = [chrom for chrom, _ in groupby(handle, key=_chromosome_of)]
    if blocks == sorted(blocks):
        return False
    rewritten = path.with_name(path.name + ".sorted.tmp")
    try:
        with tempfile.TemporaryDirectory(dir=path.parent, prefix="methylpy_merge_") as folder:
            parts: dict[str, Path] = {}
            with _open_allc_text(path) as handle:
                for index, (chrom, lines) in enumerate(groupby(handle, key=_chromosome_of)):
                    parts[chrom] = Path(folder) / f"part_{index}"
                    with parts[chrom].open("w", encoding="utf-8") as writer:
                        writer.writelines(lines)
            opener: Any = gzip.open if path.name.endswith("gz") else open
            with opener(rewritten, "wt", encoding="utf-8") as output:
                for chrom in sorted(parts):
                    with parts[chrom].open("r", encoding="utf-8") as part:
                        shutil.copyfileobj(part, output)
        os.replace(rewritten, path)
    finally:
        rewritten.unlink(missing_ok=True)
    _upstream_call(utilities.index_allc_file, str(path), reindex=True)
    return True


def merge_allc_tables(
    allc_paths: Sequence[str | Path],
    output_path: str | Path,
    *,
    mini_batch: int = 100,
    compress_output: bool = True,
    skip_snp_info: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add up the methylated and total read counts of several allc tables site by site."""

    paths = _path_list(allc_paths, "allc_paths", 2, _MAX_MERGE_FILES)
    batch = _integer(mini_batch, "mini_batch", 2, _MAX_MERGE_FILES)
    compress = _flag(compress_output, "compress_output")
    skip_snp = _flag(skip_snp_info, "skip_snp_info")
    given = _new_path(output_path, "output_path")
    final = _with_gz(given, compress)
    index = Path(str(final) + ".idx")
    leftover = Path(str(final) + ".tmp")
    for item in (final, index, leftover):
        if item.exists():
            raise ConfigurationError(f"{item.name} already exists and will not be overwritten.")
    for path in paths:
        _scan_allc(path, contiguous=True, numeric=True, min_fields=6 if skip_snp else 9)
    utilities = require("methylpy.utilities", _METHYLPY_EXTRA)
    try:
        _upstream_call(
            utilities.merge_allc_files,
            [str(path) for path in paths],
            str(given),
            num_procs=1,
            mini_batch=batch,
            compress_output=compress,
            skip_snp_info=skip_snp,
        )
        _sort_allc_chromosomes(final, utilities)
    except BaseException:
        _discard([final, index, leftover])
        raise
    order, counts = _scan_allc(final, contiguous=True, numeric=True)
    return {
        "status": "completed",
        "operation": "merge_allc_tables",
        "inputs": [str(path) for path in paths],
        "output": str(final),
        "index": str(index),
        "compressed": compress,
        "records": sum(counts.values()),
        "chromosomes": [{"chrom": chrom, "records": counts[chrom]} for chrom in sorted(order)],
        "provenance": _provenance(
            "methylated and total counts are summed over the tables at each position and the "
            "merged record gets methylated flag 1",
            note="methylpy appends .gz to the output name when compress_output is true; "
            "chromosome blocks are written in sorted order because methylpy's own order varies "
            "between runs; allc indexes are created next to the inputs",
        ),
    }


def split_mpileup_chunks(
    mpileup_path: str | Path,
    output_prefix: str | Path,
    num_chunks: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split a plain-text mpileup file into chunks whose neighbours share two lines."""

    source = _path(mpileup_path, "mpileup_path")
    chunks = _integer(num_chunks, "num_chunks", 1, _MAX_CHUNKS)
    if not isinstance(output_prefix, (str, Path)):
        raise ConfigurationError("output_prefix must be a path prefix.")
    prefix = str(Path(output_prefix).expanduser().resolve())
    if str(output_prefix).endswith(("/", os.sep)):
        raise ConfigurationError("output_prefix must end with a file-name prefix, not a slash.")
    if not Path(prefix).parent.is_dir():
        raise ConfigurationError(f"The directory for output_prefix does not exist: {prefix}")
    if _is_gzip(source):
        raise ConfigurationError(
            "mpileup_path must be plain text: methylpy sizes the chunks with wc -l, which counts "
            "the bytes of a gzip file instead of its lines."
        )
    newlines = 0
    last = b"\n"
    with source.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            newlines += block.count(b"\n")
            last = block[-1:]
    records = newlines + (0 if last == b"\n" else 1)
    chunk_size = math.ceil(newlines / chunks)
    if chunk_size > 0:
        planned = [Path(f"{prefix}{index}") for index in range(math.ceil(records / chunk_size))]
        for item in planned:
            if item.exists():
                raise ConfigurationError(f"{item.name} already exists and will not be overwritten.")
    else:
        planned = []
    utilities = require("methylpy.utilities", _METHYLPY_EXTRA)
    try:
        _upstream_call(utilities.split_mpileup_file, chunks, str(source), prefix)
    except BaseException:
        _discard(planned)
        raise
    written = [item for item in planned if item.is_file()]
    details = []
    for index, item in enumerate(written):
        lines = _line_count(item)
        new = min(chunk_size, records - index * chunk_size)
        details.append(
            {"index": index, "path": str(item), "records": lines, "overlap_records": lines - new}
        )
    return {
        "status": "completed",
        "operation": "split_mpileup_chunks",
        "input": str(source),
        "records": records,
        "chunk_size": chunk_size,
        "chunks": details,
        "provenance": _provenance(
            "consecutive chunks of ceil(lines / num_chunks) lines; every chunk after the first "
            "starts with the last two lines of its predecessor so that C contexts at chunk "
            "borders stay complete",
            note="the chunk count can be below num_chunks; the line count comes from wc -l, "
            "which ignores a last line without a newline",
        ),
    }


__all__ = [
    "add_cytosine_info_to_dmrs",
    "add_methylation_levels_to_dmrs",
    "collapse_dmr_windows",
    "dmr_residual_cutoff",
    "estimate_dmr_fdr_cutoff",
    "filter_allc_by_pvalue_cutoffs",
    "merge_allc_tables",
    "merge_dms_to_dmrs",
    "run_allc_binomial_test",
    "split_allc_by_chromosome",
    "split_allc_into_chunks",
    "split_mpileup_chunks",
]

_register_public_verbose(globals(), tuple(__all__))
