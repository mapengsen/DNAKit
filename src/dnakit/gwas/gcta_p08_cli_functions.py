"""GCTA (Genome-wide Complex Trait Analysis) command-line operations as typed DNAKit facades
(integration pack P08).

The user installs GCTA (a release binary from https://yanglab.westlake.edu.cn/software/gcta/ or
``conda install -c bioconda gcta``). Each function runs one fixed GCTA analysis with an explicit,
typed argument list (no shell). The program is resolved from ``gcta_path``, then ``gcta64`` or
``gcta`` on PATH; a missing program raises ``BackendUnavailableError`` with an install hint.
Genotypes are PLINK binary filesets (``<prefix>.bed/.bim/.fam``). Outputs are new ``out_prefix``
files: GCTA's result files are written next to the prefix and an existing result is never
overwritten. Results are JSON-serialisable dictionaries with ``status`` (``completed`` only when the
exit code is 0 and every expected output exists), ``backend``, ``version``, ``command``,
``return_code``, ``outputs`` (result key to file path), ``log_tail`` and a small parsed ``summary``
of the main result file.

Optional arguments left as ``None`` are not passed, so GCTA's own defaults apply. Tested with GCTA
1.95.4; the fastGWA-GE ``envir_file`` option needs 1.95 or later (1.94.1 has no ``--envir``). GCTA
result files get a ``.cojo`` suffix in 1.95 (``test.jma.cojo``) that 1.94 did not have
(``test.jma``); both are recognised.
"""

from __future__ import annotations

import re
import tempfile
from collections.abc import Mapping, Sequence
from functools import lru_cache
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.backends.external import execute_bounded_command
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.exceptions import BackendUnavailableError, ConfigurationError, DNAKitError

_HINT = (
    "Install GCTA 1.95 or later (release binary from "
    "https://yanglab.westlake.edu.cn/software/gcta/, or "
    "`conda install -c bioconda gcta` for 1.94.1), then put gcta64 on PATH or pass gcta_path."
)
_TESTED_VERSION = "1.95.4"
_VERSION = re.compile(r"version v?(\d+(?:\.\d+)+)")
_TIMEOUT = 3600.0
_GRM_TYPES = ("autosome", "xchr", "dominance", "inbred")
_FASTGWA_MODELS = ("mlm", "mlm_exact", "mlm_binary", "lr")
_COJO_MODES = ("slct", "cond", "joint")
_GSMR_DIRECTIONS = (0, 1, 2)
_Suffixes = Sequence[str]


# ---------------------------------------------------------------------------
# Running GCTA
# ---------------------------------------------------------------------------


def _executable(gcta_path: str | None) -> str:
    if gcta_path is not None:
        return _cli.program("gcta64", explicit=gcta_path, hint=_HINT)
    for name in ("gcta64", "gcta"):
        try:
            return _cli.program(name, explicit=None, hint=_HINT)
        except BackendUnavailableError:
            continue
    raise BackendUnavailableError("gcta64 was not found on PATH.", hint=_HINT)


@lru_cache(maxsize=8)
def _version_of(executable: str) -> str | None:
    try:
        result = execute_bounded_command(
            executable,
            ["--version"],
            backend_id="gcta",
            cwd=Path(executable).parent,
            timeout_seconds=60.0,
            max_output_bytes=100_000,
            verbose=False,
        )
    except DNAKitError:
        return None
    match = _VERSION.search(result.output or "")
    return match.group(1) if match else None


def _prefix(value: str | Path, name: str, expected: Mapping[str, _Suffixes]) -> Path:
    """An output prefix: its folder exists and none of its result files (nor the log) exist yet."""
    if not isinstance(value, (str, Path)) or not str(value):
        raise ConfigurationError(f"{name} must be a path prefix.")
    base = Path(value).expanduser().resolve()
    if not base.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {base.parent}")
    candidates = [".log"] + [suffix for group in expected.values() for suffix in group]
    for suffix in candidates:
        if Path(f"{base}{suffix}").exists():
            raise ConfigurationError(f"{name} would overwrite an existing file: {base}{suffix}")
    return base


def _found(base: Path, group: _Suffixes) -> Path | None:
    for suffix in group:
        candidate = Path(f"{base}{suffix}")
        if candidate.is_file():
            return candidate
    return None


def _run(
    arguments: Sequence[str],
    base: Path,
    *,
    required: Mapping[str, _Suffixes],
    optional: Mapping[str, _Suffixes] | None = None,
    gcta_path: str | None,
    timeout_seconds: float,
    cwd: Path | None = None,
    out_name: str | None = None,
) -> dict[str, Any]:
    """Run GCTA with ``--out base`` (or ``--out out_name``, relative to ``cwd``).

    ``status`` is completed only when the exit code is 0 and every required output exists.
    """
    executable = _executable(gcta_path)
    if (
        isinstance(timeout_seconds, bool)
        or not isinstance(timeout_seconds, (int, float))
        or timeout_seconds <= 0
    ):
        raise ConfigurationError("timeout_seconds must be a positive number.")
    working = base.parent if cwd is None else cwd
    result = execute_bounded_command(
        executable,
        [*arguments, "--out", out_name or str(base)],
        backend_id="gcta",
        cwd=working,
        timeout_seconds=float(timeout_seconds),
        max_output_bytes=_cli.MAX_OUTPUT_BYTES,
        monitored_output_paths=[Path(f"{base}.log")]
        if Path(f"{base}.log").parent == working
        else [],
        max_monitored_output_bytes=_cli.MAX_ARTIFACT_BYTES,
        verbose=False,
    )
    outputs: dict[str, str] = {}
    complete = result.return_code == 0
    for key, group in required.items():
        path = _found(base, group)
        if path is None:
            complete = False
        else:
            outputs[key] = str(path)
    for key, group in (optional or {}).items():
        path = _found(base, group)
        if path is not None:
            outputs[key] = str(path)
    log = Path(f"{base}.log")
    if log.is_file():
        outputs["log"] = str(log)
    return {
        "status": "completed" if complete else "failed",
        "backend": "gcta",
        "version": _version_of(executable),
        "tested_version": _TESTED_VERSION,
        "command": ["gcta64", *arguments, "--out", out_name or str(base)],
        "return_code": result.return_code,
        "outputs": outputs,
        "log_tail": (result.output or "")[-_cli.LOG_TAIL_CHARS :],
    }


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------


def _bfile(value: str | Path, name: str = "bfile") -> Path:
    """A PLINK binary fileset prefix whose .bed, .bim and .fam files exist."""
    base = Path(value).expanduser().resolve()
    for suffix in (".bed", ".bim", ".fam"):
        if not Path(f"{base}{suffix}").is_file():
            raise ConfigurationError(
                f"{name} needs {base}{suffix} (a PLINK binary fileset prefix)."
            )
    return base


def _grm(value: str | Path, name: str = "grm_prefix") -> Path:
    base = Path(value).expanduser().resolve()
    if not Path(f"{base}.grm.id").is_file() or not any(
        Path(f"{base}.grm{ext}").is_file() for ext in (".bin", ".gz")
    ):
        raise ConfigurationError(f"{name} needs {base}.grm.id and {base}.grm.bin (or .grm.gz).")
    return base


def _sparse_grm(value: str | Path, name: str = "sparse_grm_prefix") -> Path:
    base = Path(value).expanduser().resolve()
    if not Path(f"{base}.grm.id").is_file() or not Path(f"{base}.grm.sp").is_file():
        raise ConfigurationError(f"{name} needs {base}.grm.id and {base}.grm.sp (a sparse GRM).")
    return base


def _file(value: str | Path | None, name: str) -> Path | None:
    return None if value is None else _cli.existing_file(value, name)


def _genotypes(bfile: str | None, mbfile: str | None) -> list[str]:
    """``--bfile`` or ``--mbfile`` (a text file listing PLINK prefixes, one per chromosome)."""
    if (bfile is None) == (mbfile is None):
        raise ConfigurationError("Give exactly one of bfile and mbfile.")
    if bfile is not None:
        return ["--bfile", str(_bfile(bfile))]
    listing = _cli.existing_file(str(mbfile), "mbfile")
    for line in listing.read_text(encoding="utf-8").splitlines():
        if line.strip():
            _bfile(
                line.strip() if Path(line.strip()).is_absolute() else listing.parent / line.strip(),
                "mbfile entry",
            )
    return ["--mbfile", str(listing)]


def _flag(value: bool, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _valued(flag: str, value: object | None) -> list[str]:
    return [] if value is None else [flag, str(value)]


def _optional_int(
    flag: str, value: int | None, name: str, *, minimum: int | None = None
) -> list[str]:
    return [] if value is None else [flag, str(_cli.integer(value, name, minimum=minimum))]


def _optional_number(
    flag: str,
    value: float | None,
    name: str,
    *,
    minimum: float | None = None,
    maximum: float | None = None,
) -> list[str]:
    return (
        []
        if value is None
        else [flag, format(_cli.number(value, name, minimum=minimum, maximum=maximum), ".12g")]
    )


def _optional_file(flag: str, value: str | Path | None, name: str) -> list[str]:
    return [] if value is None else [flag, str(_cli.existing_file(value, name))]


def _common(
    *,
    keep: str | None = None,
    extract: str | None = None,
    chromosome: int | None = None,
    maf: float | None = None,
    max_maf: float | None = None,
    threads: int | None = None,
) -> list[str]:
    """The sample, SNP and thread selection options shared by many analyses."""
    arguments = _optional_file("--keep", keep, "keep") + _optional_file(
        "--extract", extract, "extract"
    )
    arguments += _optional_int("--chr", chromosome, "chromosome", minimum=1)
    arguments += _optional_number("--maf", maf, "maf", minimum=0, maximum=0.5)
    arguments += _optional_number("--max-maf", max_maf, "max_maf", minimum=0, maximum=0.5)
    arguments += _optional_int("--thread-num", threads, "threads", minimum=1)
    return arguments


def _covariates(covariate_file: str | None, quantitative_covariate_file: str | None) -> list[str]:
    return _optional_file("--covar", covariate_file, "covariate_file") + _optional_file(
        "--qcovar", quantitative_covariate_file, "quantitative_covariate_file"
    )


def _list_file(directory: Path, name: str, lines: Sequence[str]) -> Path:
    path = directory / name
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Parsing result files
# ---------------------------------------------------------------------------


def _number(value: str) -> object:
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


def _table(
    path: str | Path, *, p_column: str | None = None, top: int = 5, header: bool = True
) -> dict[str, Any]:
    """Row count, columns and the most significant rows of a whitespace-delimited result table."""
    lines = [
        line.split()
        for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
        if line.strip()
    ]
    if not lines:
        return {"rows": 0, "columns": []}
    columns = lines[0] if header else []
    body = lines[1:] if header else lines
    rows = (
        [{name: _number(value) for name, value in zip(columns, row, strict=False)} for row in body]
        if header
        else []
    )
    summary: dict[str, Any] = {"rows": len(body), "columns": columns}
    if p_column is not None and p_column in columns:
        scored = [row for row in rows if isinstance(row.get(p_column), (int, float))]
        scored.sort(key=lambda row: float(str(row[p_column])))
        summary["top_rows"] = scored[:top]
    elif rows:
        summary["first_rows"] = rows[:top]
    return summary


def _hsq(path: str | Path) -> dict[str, Any]:
    """Parse a GCTA ``.hsq`` file into variance components and named statistics."""
    components: list[dict[str, Any]] = []
    statistics: dict[str, Any] = {}
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    for line in lines[1:]:
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 3:
            components.append(
                {"source": parts[0], "variance": _number(parts[1]), "se": _number(parts[2])}
            )
        elif len(parts) == 2 and parts[0]:
            statistics[parts[0]] = _number(parts[1])
    return {"variance_components": components, "statistics": statistics}


def _hereg(path: str | Path) -> dict[str, Any]:
    """Parse a GCTA ``.HEreg`` file into the HE-CP and HE-SD tables."""
    blocks: dict[str, list[dict[str, Any]]] = {}
    current: str | None = None
    columns: list[str] = []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] in ("HE-CP", "HE-SD"):
            current = parts[0]
            blocks[current] = []
            columns = []
        elif parts[0] == "Coefficient":
            columns = [name.lower() for name in parts[1:]]
        elif current is not None:
            values = [_number(value) for value in parts[1:]]
            blocks[current].append(
                {
                    "coefficient": parts[0],
                    **dict(zip(columns or ["estimate"], values, strict=False)),
                }
            )
    return blocks


def _log_counts(log_tail: str, pattern: str) -> tuple[int, ...] | None:
    match = re.search(pattern, log_tail)
    return tuple(int(group) for group in match.groups()) if match else None


def _summarise(result: dict[str, Any], key: str, **table_options: Any) -> dict[str, Any]:
    if key in result["outputs"]:
        result["summary"] = _table(result["outputs"][key], **table_options)
    return result


# ---------------------------------------------------------------------------
# Genotype export and relationship matrices
# ---------------------------------------------------------------------------


def gcta_recode_standardized_genotypes(
    bfile: str,
    out_prefix: str,
    *,
    keep: str | None = None,
    extract: str | None = None,
    chromosome: int | None = None,
    maf: float | None = None,
    max_maf: float | None = None,
    threads: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the standardised SNP genotype matrix of a PLINK fileset as compressed text (``gcta64
    --recode-std``).

    Each genotype is standardised as (x - 2p) / sqrt(2p(1-p)) with p the allele frequency estimated
    from the data. The matrix is written to ``<out_prefix>.xmat.gz``; ``n_individuals`` and
    ``n_snps`` come from GCTA's log.
    """
    base = _prefix(out_prefix, "out_prefix", {"xmat": (".xmat.gz",)})
    selection = _common(
        keep=keep, extract=extract, chromosome=chromosome, maf=maf, max_maf=max_maf, threads=threads
    )
    arguments = ["--bfile", str(_bfile(bfile)), "--recode-std", *selection]
    result = _run(
        arguments,
        base,
        required={"xmat": (".xmat.gz",)},
        gcta_path=gcta_path,
        timeout_seconds=timeout_seconds,
    )
    counts = _log_counts(result["log_tail"], r"Genotype data for (\d+) individuals and (\d+) SNPs")
    if counts is not None:
        result["n_individuals"], result["n_snps"] = counts
    return result


def gcta_make_grm(
    out_prefix: str,
    *,
    bfile: str | None = None,
    mbfile: str | None = None,
    grm_type: str = "autosome",
    keep: str | None = None,
    extract: str | None = None,
    chromosome: int | None = None,
    maf: float | None = None,
    max_maf: float | None = None,
    threads: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute a dense genetic relationship matrix (GRM) from SNP genotypes and save it in binary
    format (``gcta64 --make-grm*``).

    ``grm_type`` selects the matrix: ``autosome`` (``--make-grm``, the additive GRM of the autosomal
    SNPs), ``xchr`` (``--make-grm-xchr``, the GRM of the X chromosome), ``dominance``
    (``--make-grm-d``, the dominance GRM, written as ``<out_prefix>.d.grm.*``) or ``inbred``
    (``--autosome --make-grm-inbred``, the GRM for an inbred population). Give genotypes as one
    PLINK prefix ``bfile`` or a list file ``mbfile`` of per-chromosome prefixes. The files
    ``<out_prefix>.grm.bin``, ``.grm.N.bin`` and ``.grm.id`` are written; ``grm_prefix`` is the
    prefix to pass on.
    """
    kind = _cli.choice(grm_type, _GRM_TYPES, "grm_type")
    stem = ".d.grm" if kind == "dominance" else ".grm"
    expected = {
        "grm_bin": (f"{stem}.bin",),
        "grm_id": (f"{stem}.id",),
        "grm_n_bin": (f"{stem}.N.bin",),
    }
    base = _prefix(out_prefix, "out_prefix", expected)
    flag = {
        "autosome": ["--make-grm"],
        "xchr": ["--make-grm-xchr"],
        "dominance": ["--make-grm-d"],
        "inbred": ["--autosome", "--make-grm-inbred"],
    }[kind]
    arguments = (
        _genotypes(bfile, mbfile)
        + flag
        + _common(
            keep=keep,
            extract=extract,
            chromosome=chromosome,
            maf=maf,
            max_maf=max_maf,
            threads=threads,
        )
    )
    result = _run(
        arguments, base, required=expected, gcta_path=gcta_path, timeout_seconds=timeout_seconds
    )
    result["grm_prefix"] = f"{base}.d" if kind == "dominance" else str(base)
    if "grm_id" in result["outputs"]:
        result["n_individuals"] = _count_lines(result["outputs"]["grm_id"])
    return result


def _count_lines(path: str) -> int:
    return sum(1 for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip())


def gcta_make_sparse_grm(
    grm_prefix: str,
    out_prefix: str,
    *,
    cutoff: float = 0.05,
    output_format: str = "sparse",
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Threshold a dense GRM into a sparse GRM for fastGWA (``gcta64 --grm P --make-bK-sparse`` or
    ``--make-bK``).

    Off-diagonal elements below ``cutoff`` (GCTA's recommended 0.05) are set to zero.
    ``output_format="sparse"`` (``--make-bK-sparse``) writes ``<out_prefix>.grm.sp`` and
    ``.grm.id``, the format ``--grm-sparse`` reads; ``output_format="dense"`` (``--make-bK``) writes
    the same thresholded matrix as a binary GRM (``.grm.bin``, ``.grm.N.bin``, ``.grm.id``).
    ``n_pairs`` counts the retained pairs of the sparse file.
    """
    source = _grm(grm_prefix)
    threshold = format(_cli.number(cutoff, "cutoff", minimum=0, maximum=1), ".12g")
    kind = _cli.choice(output_format, ("sparse", "dense"), "output_format")
    expected: dict[str, _Suffixes] = (
        {"grm_sp": (".grm.sp",), "grm_id": (".grm.id",)}
        if kind == "sparse"
        else {"grm_bin": (".grm.bin",), "grm_id": (".grm.id",), "grm_n_bin": (".grm.N.bin",)}
    )
    base = _prefix(out_prefix, "out_prefix", expected)
    arguments = [
        "--grm",
        str(source),
        "--make-bK-sparse" if kind == "sparse" else "--make-bK",
        threshold,
    ]
    result = _run(
        arguments, base, required=expected, gcta_path=gcta_path, timeout_seconds=timeout_seconds
    )
    result["grm_prefix"] = str(base)
    if "grm_sp" in result["outputs"]:
        result["n_pairs"] = _count_lines(result["outputs"]["grm_sp"])
    return result


def gcta_inbreeding_coefficient(
    bfile: str,
    out_prefix: str,
    *,
    keep: str | None = None,
    extract: str | None = None,
    maf: float | None = None,
    max_maf: float | None = None,
    threads: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate the inbreeding coefficient of every individual from autosomal SNPs (``gcta64
    --autosome --ibc``).

    Writes ``<out_prefix>.ibc`` with three estimates per individual (``Fhat1``
    variance-standardised, ``Fhat2`` excess homozygosity, ``Fhat3`` allele-sharing); ``summary``
    holds their means.
    """
    base = _prefix(out_prefix, "out_prefix", {"ibc": (".ibc",)})
    selection = _common(keep=keep, extract=extract, maf=maf, max_maf=max_maf, threads=threads)
    arguments = ["--bfile", str(_bfile(bfile)), "--autosome", "--ibc", *selection]
    result = _run(
        arguments,
        base,
        required={"ibc": (".ibc",)},
        gcta_path=gcta_path,
        timeout_seconds=timeout_seconds,
    )
    if "ibc" in result["outputs"]:
        table = _table(result["outputs"]["ibc"], top=0)
        rows = [
            line.split()
            for line in Path(result["outputs"]["ibc"]).read_text(encoding="utf-8").splitlines()[1:]
            if line.strip()
        ]
        means = {}
        for index, name in enumerate(table["columns"]):
            if name.startswith("Fhat") and rows:
                means[name] = sum(float(row[index]) for row in rows) / len(rows)
        result["summary"] = {"individuals": table["rows"], "mean": means}
    return result


def gcta_unify_grm(
    grm_prefixes: Sequence[str],
    out_prefix: str,
    *,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Restrict several GRMs to the individuals they have in common, in the same order (``gcta64
    --mgrm list --unify-grm``).

    GCTA names and places the unified matrices by itself (``<input>_<out name>`` next to the input),
    so this function runs it on links to the inputs in a scratch folder and moves the results to
    stable names: ``<out_prefix>.1.grm.bin`` / ``.grm.id`` for the first input,
    ``<out_prefix>.2...`` for the second and so on. Those prefixes are listed in ``unified`` and can
    be passed to ``gcta_reml``. ``common_individuals`` is the number of individuals kept.
    """
    if isinstance(grm_prefixes, (str, bytes)) or len(grm_prefixes) < 2:
        raise ConfigurationError("grm_prefixes must list at least two GRM prefixes.")
    sources = [_grm(prefix, "grm_prefixes") for prefix in grm_prefixes]
    count = len(sources)
    expected = {
        f"grm_{k}_{kind}": (f".{k}.grm.{kind}",)
        for k in range(1, count + 1)
        for kind in ("id", "bin")
    }
    base = _prefix(out_prefix, "out_prefix", expected)
    with tempfile.TemporaryDirectory(prefix="dnakit-gcta-") as folder:
        scratch = Path(folder)
        names = []
        for k, source in enumerate(sources, start=1):
            names.append(f"grm{k}")
            (scratch / f"grm{k}.grm.id").symlink_to(f"{source}.grm.id")
            binary = (
                Path(f"{source}.grm.bin")
                if Path(f"{source}.grm.bin").is_file()
                else Path(f"{source}.grm.gz")
            )
            (scratch / f"grm{k}{binary.name[len(source.name) :]}").symlink_to(binary)
        listing = _list_file(scratch, "grms.txt", names)
        result = _run(
            ["--mgrm", str(listing), "--unify-grm"],
            scratch / "unify",
            required={},
            gcta_path=gcta_path,
            timeout_seconds=timeout_seconds, cwd=scratch, out_name="unify",
        )  # fmt: skip
        unified: list[str] = []
        if result["return_code"] == 0:
            for k in range(1, count + 1):
                ids = scratch / f"grm{k}_unify.grm.id"
                binary = scratch / f"grm{k}_unify.grm.bin"
                if not ids.is_file() or not binary.is_file():
                    break
                ids.replace(f"{base}.{k}.grm.id")
                binary.replace(f"{base}.{k}.grm.bin")
                unified.append(f"{base}.{k}")
        if (scratch / "unify.log").is_file():
            (scratch / "unify.log").replace(f"{base}.log")
    result["outputs"] = {"log": f"{base}.log"} if Path(f"{base}.log").is_file() else {}
    result["command"] = [
        "gcta64",
        "--mgrm",
        "<list of the input GRMs>",
        "--unify-grm",
        "--out",
        "unify",
    ]
    result["unified"] = unified
    for k, prefix in enumerate(unified, start=1):
        result["outputs"][f"grm_{k}_id"] = f"{prefix}.grm.id"
        result["outputs"][f"grm_{k}_bin"] = f"{prefix}.grm.bin"
    if unified:
        result["common_individuals"] = _count_lines(f"{unified[0]}.grm.id")
    if result["status"] == "completed" and len(unified) != count:
        result["status"] = "failed"
    return result


# ---------------------------------------------------------------------------
# REML, prediction and Haseman-Elston regression
# ---------------------------------------------------------------------------


def _grm_args(
    grm_prefix: str | None, grm_prefixes: Sequence[str] | None, directory: Path
) -> tuple[list[str], Path | None]:
    """``--grm`` for one GRM or ``--mgrm`` for several (the list file goes into ``directory``)."""
    if (grm_prefix is None) == (not grm_prefixes):
        raise ConfigurationError("Give exactly one of grm_prefix and grm_prefixes.")
    if grm_prefix is not None:
        return ["--grm", str(_grm(grm_prefix))], None
    if isinstance(grm_prefixes, (str, bytes)):
        raise ConfigurationError("grm_prefixes must be a list of GRM prefixes.")
    listing = _list_file(
        directory, "grms.txt", [str(_grm(prefix, "grm_prefixes")) for prefix in grm_prefixes or []]
    )
    return ["--mgrm", str(listing)], listing


def gcta_reml(
    phenotype_file: str,
    out_prefix: str,
    *,
    grm_prefix: str | None = None,
    grm_prefixes: Sequence[str] | None = None,
    phenotype_column: int | None = None,
    covariate_file: str | None = None,
    quantitative_covariate_file: str | None = None,
    gxe_file: str | None = None,
    predict: str = "none",
    no_constrain: bool = False,
    max_iterations: int | None = None,
    algorithm: int | None = None,
    prevalence: float | None = None,
    estimate_fixed_effects: bool = False,
    keep: str | None = None,
    threads: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate variance components and SNP-based heritability by GREML, optionally with GxE or
    genetic-value prediction (``gcta64 --reml``).

    ``grm_prefix`` is one GRM; ``grm_prefixes`` fits several GRMs (``--mgrm``). ``phenotype_file``
    has columns FID, IID and the phenotype(s); ``phenotype_column`` (``--mpheno``) picks the n-th
    phenotype. ``gxe_file`` (``--gxe``) adds a genotype-by-environment variance component from a
    discrete environment factor. ``predict`` is ``none``, ``blup`` (``--reml-pred-rand``: BLUP of
    every individual's total genetic value, written to ``<out_prefix>.indi.blp``) or ``cvblup``
    (``--cvblup``: leave-one-out BLUP, ``<out_prefix>.indi.cvblp``). ``no_constrain`` lets
    components leave the parameter space (``--reml-no-constrain``); ``prevalence`` transforms the
    heritability of a binary trait to the liability scale. The ``.hsq`` file is parsed into
    ``variance_components`` (source, variance, SE) and ``statistics`` (logL, logL0, LRT, df, Pval,
    n).
    """
    kind = _cli.choice(predict, ("none", "blup", "cvblup"), "predict")
    expected: dict[str, _Suffixes] = {"hsq": (".hsq",)}
    if kind == "blup":
        expected["indi_blp"] = (".indi.blp",)
    if kind == "cvblup":
        expected["indi_cvblp"] = (".indi.cvblp",)
    base = _prefix(out_prefix, "out_prefix", expected)
    phenotype = _cli.existing_file(phenotype_file, "phenotype_file")
    with tempfile.TemporaryDirectory(prefix="dnakit-gcta-") as folder:
        grm_args, _ = _grm_args(grm_prefix, grm_prefixes, Path(folder))
        arguments = [*grm_args, "--pheno", str(phenotype), "--reml"]
        arguments += _optional_int("--mpheno", phenotype_column, "phenotype_column", minimum=1)
        arguments += _covariates(covariate_file, quantitative_covariate_file)
        arguments += _optional_file("--gxe", gxe_file, "gxe_file")
        arguments += {"none": [], "blup": ["--reml-pred-rand"], "cvblup": ["--cvblup"]}[kind]
        arguments += ["--reml-no-constrain"] if _flag(no_constrain, "no_constrain") else []
        arguments += _optional_int("--reml-maxit", max_iterations, "max_iterations", minimum=1)
        arguments += _optional_int("--reml-alg", algorithm, "algorithm", minimum=0)
        arguments += _optional_number(
            "--prevalence", prevalence, "prevalence", minimum=0, maximum=1
        )
        arguments += (
            ["--reml-est-fix"] if _flag(estimate_fixed_effects, "estimate_fixed_effects") else []
        )
        arguments += _common(keep=keep, threads=threads)
        result = _run(
            arguments, base, required=expected, gcta_path=gcta_path, timeout_seconds=timeout_seconds
        )
    if "hsq" in result["outputs"]:
        result["summary"] = _hsq(result["outputs"]["hsq"])
    for key in ("indi_blp", "indi_cvblp"):
        if key in result["outputs"]:
            result["n_predictions"] = _count_lines(result["outputs"][key])
    return result


def gcta_blup_snp(
    bfile: str,
    indi_blp_file: str,
    out_prefix: str,
    *,
    autosome: bool = False,
    keep: str | None = None,
    extract: str | None = None,
    threads: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert the BLUP of individuals' genetic values into BLUP effects of every SNP (``gcta64
    --blup-snp``).

    ``indi_blp_file`` is the ``.indi.blp`` file written by ``gcta_reml(predict='blup')``. The SNP
    effects (``<out_prefix>.snp.blp``: SNP, reference allele, one effect per fitted GRM) can be used
    to score new samples. ``autosome=True`` restricts the SNPs to chromosomes 1 to 22, which GCTA
    needs when the fileset has sex chromosomes.
    """
    base = _prefix(out_prefix, "out_prefix", {"snp_blp": (".snp.blp",)})
    arguments = [
        "--bfile",
        str(_bfile(bfile)),
        "--blup-snp",
        str(_cli.existing_file(indi_blp_file, "indi_blp_file")),
    ]
    arguments += ["--autosome"] if _flag(autosome, "autosome") else []
    arguments += _common(keep=keep, extract=extract, threads=threads)
    result = _run(
        arguments,
        base,
        required={"snp_blp": (".snp.blp",)},
        gcta_path=gcta_path,
        timeout_seconds=timeout_seconds,
    )
    if "snp_blp" in result["outputs"]:
        result["n_snps"] = _count_lines(result["outputs"]["snp_blp"])
    return result


def gcta_reml_bivar(
    phenotype_file: str,
    out_prefix: str,
    *,
    grm_prefix: str | None = None,
    grm_prefixes: Sequence[str] | None = None,
    trait_columns: Sequence[int] | None = None,
    covariate_file: str | None = None,
    quantitative_covariate_file: str | None = None,
    no_constrain: bool = False,
    max_iterations: int | None = None,
    keep: str | None = None,
    threads: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Bivariate GREML: the genetic correlation (rG) between two traits (``gcta64 --reml-bivar``).

    ``phenotype_file`` has columns FID, IID and the two traits (missing values ``NA``);
    ``trait_columns`` is a pair such as ``[1, 2]`` to pick traits from a wider file. ``summary``
    holds ``variance_components`` (V(G)_tr1, V(G)_tr2, C(G)_tr12, V(e)_..., rG with SE) and
    ``statistics``.
    """
    expected: dict[str, _Suffixes] = {"hsq": (".hsq",)}
    base = _prefix(out_prefix, "out_prefix", expected)
    phenotype = _cli.existing_file(phenotype_file, "phenotype_file")
    pair: list[str] = []
    if trait_columns is not None:
        if isinstance(trait_columns, (str, bytes)) or len(trait_columns) != 2:
            raise ConfigurationError("trait_columns must be two trait column numbers.")
        pair = [str(_cli.integer(value, "trait_columns", minimum=1)) for value in trait_columns]
    with tempfile.TemporaryDirectory(prefix="dnakit-gcta-") as folder:
        grm_args, _ = _grm_args(grm_prefix, grm_prefixes, Path(folder))
        arguments = [*grm_args, "--pheno", str(phenotype), "--reml-bivar", *pair]
        arguments += _covariates(covariate_file, quantitative_covariate_file)
        arguments += ["--reml-bivar-no-constrain"] if _flag(no_constrain, "no_constrain") else []
        arguments += _optional_int("--reml-maxit", max_iterations, "max_iterations", minimum=1)
        arguments += _common(keep=keep, threads=threads)
        result = _run(
            arguments, base, required=expected, gcta_path=gcta_path, timeout_seconds=timeout_seconds
        )
    if "hsq" in result["outputs"]:
        result["summary"] = _hsq(result["outputs"]["hsq"])
    return result


def gcta_he_regression(
    phenotype_file: str,
    out_prefix: str,
    *,
    grm_prefix: str | None = None,
    grm_prefixes: Sequence[str] | None = None,
    bivariate: bool = False,
    trait_columns: Sequence[int] = (1, 2),
    phenotype_column: int | None = None,
    keep: str | None = None,
    threads: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Haseman-Elston regression for SNP-based heritability or, with ``bivariate=True``, genetic
    correlation (``gcta64 --HEreg`` / ``--HEreg-bivar``).

    ``bivariate=False`` regresses the cross products (HE-CP) and squared differences (HE-SD) of one
    phenotype (``phenotype_column`` = ``--mpheno``) on the GRM. ``bivariate=True`` uses the two
    traits named by ``trait_columns`` (default ``[1, 2]``) and reports their genetic correlation
    (rG). ``summary`` parses ``<out_prefix>.HEreg`` into the ``HE-CP`` and ``HE-SD`` tables of
    coefficient, estimate, OLS and jackknife standard errors and p-values.
    """
    expected: dict[str, _Suffixes] = {"hereg": (".HEreg",)}
    base = _prefix(out_prefix, "out_prefix", expected)
    phenotype = _cli.existing_file(phenotype_file, "phenotype_file")
    two = bool(_flag(bivariate, "bivariate"))
    pair: list[str] = []
    if two:
        if isinstance(trait_columns, (str, bytes)) or len(trait_columns) != 2:
            raise ConfigurationError("trait_columns must be two trait column numbers.")
        pair = [str(_cli.integer(value, "trait_columns", minimum=1)) for value in trait_columns]
    with tempfile.TemporaryDirectory(prefix="dnakit-gcta-") as folder:
        grm_args, _ = _grm_args(grm_prefix, grm_prefixes, Path(folder))
        arguments = [*grm_args, "--pheno", str(phenotype)]
        arguments += ["--HEreg-bivar", *pair] if two else ["--HEreg"]
        arguments += _optional_int("--mpheno", phenotype_column, "phenotype_column", minimum=1)
        arguments += _common(keep=keep, threads=threads)
        result = _run(
            arguments, base, required=expected, gcta_path=gcta_path, timeout_seconds=timeout_seconds
        )
    if "hereg" in result["outputs"]:
        result["summary"] = _hereg(result["outputs"]["hereg"])
    return result


# ---------------------------------------------------------------------------
# Association analyses
# ---------------------------------------------------------------------------


def gcta_fastgwa(
    phenotype_file: str,
    out_prefix: str,
    *,
    bfile: str | None = None,
    mbfile: str | None = None,
    model: str = "mlm",
    sparse_grm_prefix: str | None = None,
    covariate_file: str | None = None,
    quantitative_covariate_file: str | None = None,
    envir_file: str | None = None,
    model_only: bool = False,
    keep: str | None = None,
    extract: str | None = None,
    chromosome: int | None = None,
    maf: float | None = None,
    geno: float | None = None,
    threads: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """fastGWA association tests of every SNP (``gcta64 --fastGWA-mlm`` / ``--fastGWA-mlm-exact`` /
    ``--fastGWA-mlm-binary`` / ``--fastGWA-lr``).

    ``model`` selects the analysis: ``mlm`` (mixed linear model with a sparse GRM and the
    GRAMMAR-Gamma approximation), ``mlm_exact`` (the exact test without GRAMMAR-Gamma),
    ``mlm_binary`` (fastGWA-GLMM for a binary trait coded 0/1; needs at least 2000 autosomal SNPs
    for tuning and is slow without covariates) or ``lr`` (linear regression, no sparse GRM). The
    three mixed models need ``sparse_grm_prefix`` (see ``gcta_make_sparse_grm``). ``envir_file``
    (FID, IID, continuous environment value; GCTA 1.95 or later) turns ``mlm`` or ``lr`` into the
    fastGWA-GE genotype-by-environment interaction test. ``model_only=True`` stops after
    variance-component estimation and saves the model
    (``<out_prefix>.fastGWA.mdl.bin2``/``.mdl.id``) for ``gcta_fastgwa_burden``. The statistics
    table ``<out_prefix>.fastGWA`` is summarised in ``summary`` (row count, columns, smallest-p
    rows).
    """
    kind = _cli.choice(model, _FASTGWA_MODELS, "model")
    if kind == "lr" and sparse_grm_prefix is not None:
        raise ConfigurationError("model='lr' cannot use a sparse GRM.")
    if kind != "lr" and sparse_grm_prefix is None:
        raise ConfigurationError(f"model='{kind}' needs sparse_grm_prefix.")
    only = _flag(model_only, "model_only")
    expected: dict[str, _Suffixes] = (
        {"model": (".fastGWA.mdl.bin2", ".fastGWA.mdl.bin"), "model_ids": (".fastGWA.mdl.id",)}
        if only
        else {"fastgwa": (".fastGWA",)}
    )
    base = _prefix(out_prefix, "out_prefix", expected)
    flag = {
        "mlm": "--fastGWA-mlm",
        "mlm_exact": "--fastGWA-mlm-exact",
        "mlm_binary": "--fastGWA-mlm-binary",
        "lr": "--fastGWA-lr",
    }[kind]
    arguments = _genotypes(bfile, mbfile)
    arguments += (
        ["--grm-sparse", str(_sparse_grm(sparse_grm_prefix))]
        if sparse_grm_prefix is not None
        else []
    )
    arguments += [flag, "--pheno", str(_cli.existing_file(phenotype_file, "phenotype_file"))]
    arguments += _covariates(covariate_file, quantitative_covariate_file)
    arguments += _optional_file("--envir", envir_file, "envir_file")
    arguments += ["--model-only"] if only else []
    arguments += _common(
        keep=keep, extract=extract, chromosome=chromosome, maf=maf, threads=threads
    )
    arguments += _optional_number("--geno", geno, "geno", minimum=0, maximum=1)
    result = _run(
        arguments, base, required=expected, gcta_path=gcta_path, timeout_seconds=timeout_seconds
    )
    result["model_prefix"] = f"{base}.fastGWA" if only else None
    return _summarise(result, "fastgwa", p_column="P" if envir_file is None else "P_G")


def gcta_mlma(
    bfile: str,
    phenotype_file: str,
    out_prefix: str,
    *,
    grm_prefix: str | None = None,
    grm_prefixes: Sequence[str] | None = None,
    loco: bool = False,
    subtract_grm_prefix: str | None = None,
    covariate_file: str | None = None,
    quantitative_covariate_file: str | None = None,
    no_preadj_covar: bool = False,
    keep: str | None = None,
    extract: str | None = None,
    chromosome: int | None = None,
    maf: float | None = None,
    threads: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mixed linear model association of every SNP, optionally leaving one chromosome out (``gcta64
    --mlma`` / ``--mlma-loco`` / ``--mlma-subtract-grm``).

    Without ``loco`` the polygenic effect comes from ``grm_prefix`` (or the several
    ``grm_prefixes``); if neither is given GCTA builds the GRM from ``bfile``, which includes the
    candidate SNP (MLMi). ``loco=True`` runs MLMA-LOCO: the GRM leaves out the chromosome of the
    candidate SNP and is re-estimated per chromosome (``--mlma-loco``, results in
    ``<out_prefix>.loco.mlma``). ``subtract_grm_prefix`` (``--mlma-subtract-grm``) subtracts a
    chromosome's GRM from the full ``grm_prefix`` GRM to run LOCO one chromosome at a time (use with
    ``chromosome``). ``no_preadj_covar`` fits covariates together with each SNP instead of
    pre-adjusting the phenotype. The table ``<out_prefix>.mlma`` has the columns Chr, SNP, bp, A1,
    A2, Freq, b, se, p.
    """
    use_loco = _flag(loco, "loco")
    if use_loco and (grm_prefix is not None or grm_prefixes or subtract_grm_prefix is not None):
        raise ConfigurationError(
            "loco=True builds its own GRMs; do not pass grm_prefix, grm_prefixes or "
            "subtract_grm_prefix."
        )
    if subtract_grm_prefix is not None and grm_prefix is None:
        raise ConfigurationError("subtract_grm_prefix needs the full GRM in grm_prefix.")
    if grm_prefix is not None and grm_prefixes:
        raise ConfigurationError("Give grm_prefix or grm_prefixes, not both.")
    suffix = ".loco.mlma" if use_loco else ".mlma"
    base = _prefix(out_prefix, "out_prefix", {"mlma": (suffix,)})
    arguments = ["--bfile", str(_bfile(bfile)), "--mlma-loco" if use_loco else "--mlma"]
    arguments += ["--pheno", str(_cli.existing_file(phenotype_file, "phenotype_file"))]
    with tempfile.TemporaryDirectory(prefix="dnakit-gcta-") as folder:
        if grm_prefix is not None or grm_prefixes:
            grm_args, _ = _grm_args(grm_prefix, grm_prefixes, Path(folder))
            arguments += grm_args
        if subtract_grm_prefix is not None:
            arguments += [
                "--mlma-subtract-grm",
                str(_grm(subtract_grm_prefix, "subtract_grm_prefix")),
            ]
        arguments += _covariates(covariate_file, quantitative_covariate_file)
        arguments += ["--mlma-no-preadj-covar"] if _flag(no_preadj_covar, "no_preadj_covar") else []
        arguments += _common(
            keep=keep, extract=extract, chromosome=chromosome, maf=maf, threads=threads
        )
        result = _run(
            arguments,
            base,
            required={"mlma": (suffix,)},
            gcta_path=gcta_path,
            timeout_seconds=timeout_seconds,
        )
    return _summarise(result, "mlma", p_column="p")


def gcta_fastgwa_burden(
    bfile: str,
    model_prefix: str,
    set_list_file: str,
    out_prefix: str,
    *,
    max_maf: float | None = None,
    maf: float | None = None,
    extract: str | None = None,
    geno: float | None = None,
    threads: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """fastGWA-BB burden test of rare variants per gene for a binary trait (``gcta64 --burden``).

    Uses the model saved by ``gcta_fastgwa(model='mlm_binary', model_only=True)``: pass its
    ``model_prefix`` (the ``<out_prefix>.fastGWA`` returned there). ``set_list_file`` lists genes
    with columns ``gene chr start end strand`` (header line first). ``max_maf`` removes common
    variants (the manual uses 0.01). The table ``<out_prefix>.fastGWA`` gives per gene the number of
    variants, mean MAF, score statistic, effect and SPA-corrected p.
    """
    base = _prefix(out_prefix, "out_prefix", {"fastgwa": (".fastGWA",)})
    model = Path(model_prefix).expanduser().resolve()
    if not Path(f"{model}.mdl.id").is_file():
        raise ConfigurationError(
            f"model_prefix needs {model}.mdl.id from gcta_fastgwa(model_only=True)."
        )
    arguments = ["--burden", "--bfile", str(_bfile(bfile)), "--load-model", str(model)]
    arguments += ["--set-list", str(_cli.existing_file(set_list_file, "set_list_file"))]
    arguments += _optional_number("--max-maf", max_maf, "max_maf", minimum=0, maximum=0.5)
    arguments += _optional_number("--maf", maf, "maf", minimum=0, maximum=0.5)
    arguments += _optional_file("--extract", extract, "extract")
    arguments += _optional_number("--geno", geno, "geno", minimum=0, maximum=1)
    arguments += _optional_int("--thread-num", threads, "threads", minimum=1)
    result = _run(
        arguments,
        base,
        required={"fastgwa": (".fastGWA",)},
        gcta_path=gcta_path,
        timeout_seconds=timeout_seconds,
    )
    return _summarise(result, "fastgwa", p_column="P")


def gcta_acat(
    snp_list_file: str,
    gene_list_file: str,
    output_path: str,
    *,
    max_maf: float | None = None,
    min_mac: int | None = None,
    window_kb: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """ACAT-V gene-based rare-variant test from single-variant p-values (``gcta64 --acat``).

    ``snp_list_file`` is a fastGWA(-GLMM) association table (``gcta_fastgwa`` output);
    ``gene_list_file`` has the columns chromosome, start, end and gene id. Variants with MAF above
    ``max_maf`` or minor allele count below ``min_mac`` are skipped; ``window_kb`` widens each gene
    (default 0). The result table (CHR, GENE, START, END, SNP_NUM, P_ACAT) is written to
    ``output_path`` itself and summarised in ``summary``.
    """
    target = Path(output_path).expanduser().resolve()
    if target.exists():
        raise ConfigurationError(f"output_path must not exist yet: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"The folder for output_path does not exist: {target.parent}")
    arguments = ["--acat", "--snp-list", str(_cli.existing_file(snp_list_file, "snp_list_file"))]
    arguments += ["--gene-list", str(_cli.existing_file(gene_list_file, "gene_list_file"))]
    arguments += _optional_number("--max-maf", max_maf, "max_maf", minimum=0, maximum=0.5)
    arguments += _optional_int("--min-mac", min_mac, "min_mac", minimum=0)
    arguments += _optional_int("--wind", window_kb, "window_kb", minimum=0)
    result = _run(
        arguments,
        target,
        required={"acat": ("",)},
        gcta_path=gcta_path,
        timeout_seconds=timeout_seconds,
    )
    return _summarise(result, "acat", p_column="P_ACAT")


# ---------------------------------------------------------------------------
# Summary-statistics analyses
# ---------------------------------------------------------------------------


def gcta_cojo(
    bfile: str,
    summary_file: str,
    out_prefix: str,
    *,
    mode: str = "slct",
    snp_list_file: str | None = None,
    top_snps: int | None = None,
    p_threshold: float | None = None,
    window_kb: int | None = None,
    collinearity: float | None = None,
    diff_freq: float | None = None,
    genomic_control: bool | float = False,
    keep: str | None = None,
    extract: str | None = None,
    chromosome: int | None = None,
    maf: float | None = None,
    threads: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Conditional and joint multiple-SNP analysis of GWAS summary statistics, using a genotype
    fileset as the LD reference (``gcta64 --cojo-slct`` / ``--cojo-cond`` / ``--cojo-joint``).

    ``summary_file`` is a COJO ``.ma`` file with the columns SNP, A1, A2, freq, b, se, p, N (give
    statistics for all SNPs; restrict the analysis with ``extract``/``chromosome``). ``mode='slct'``
    selects independent signals by stepwise regression (``p_threshold``, default 5e-8; ``top_snps``
    instead selects a fixed number without a threshold); ``mode='joint'`` fits all included SNPs
    jointly; ``mode='cond'`` tests every SNP conditional on the SNPs in ``snp_list_file``.
    ``window_kb`` (default 10000) is the distance beyond which SNPs are in linkage equilibrium,
    ``collinearity`` the R-squared cut-off (default 0.9), ``diff_freq`` the largest allowed
    frequency difference to the reference (default 0.2); ``genomic_control=True`` (or a value such
    as 1.05) adjusts the p-values. Results: ``jma`` (selected or joint effects) and ``ldr`` (their
    LD matrix) for ``slct``/``joint``, ``cma`` (conditional effects) for ``cond`` and, with
    ``slct``, for the remaining SNPs.
    """
    kind = _cli.choice(mode, _COJO_MODES, "mode")
    if top_snps is not None and kind != "slct":
        raise ConfigurationError("top_snps only applies to mode='slct'.")
    if kind == "cond" and snp_list_file is None:
        raise ConfigurationError("mode='cond' needs snp_list_file.")
    if kind != "cond" and snp_list_file is not None:
        raise ConfigurationError("snp_list_file only applies to mode='cond'.")
    if p_threshold is not None and (kind != "slct" or top_snps is not None):
        raise ConfigurationError("p_threshold only applies to mode='slct' without top_snps.")
    jma: _Suffixes = (".jma.cojo", ".jma")
    ldr: _Suffixes = (".ldr.cojo", ".jma.ldr")
    cma: _Suffixes = (".cma.cojo", ".cma")
    required: dict[str, _Suffixes] = {"cma": cma} if kind == "cond" else {"jma": jma, "ldr": ldr}
    optional: dict[str, _Suffixes] = (
        {"given": (".given.cojo", ".given")} if kind == "cond" else {"cma": cma}
    )
    base = _prefix(out_prefix, "out_prefix", {**required, **optional})
    arguments = [
        "--bfile",
        str(_bfile(bfile)),
        "--cojo-file",
        str(_cli.existing_file(summary_file, "summary_file")),
    ]
    if kind == "slct":
        arguments += (
            ["--cojo-top-SNPs", str(_cli.integer(top_snps, "top_snps", minimum=1))]
            if top_snps is not None
            else ["--cojo-slct"]
        )
    elif kind == "joint":
        arguments += ["--cojo-joint"]
    else:
        arguments += ["--cojo-cond", str(_cli.existing_file(snp_list_file or "", "snp_list_file"))]
    arguments += _optional_number("--cojo-p", p_threshold, "p_threshold", minimum=0, maximum=1)
    arguments += _optional_int("--cojo-wind", window_kb, "window_kb", minimum=1)
    arguments += _optional_number(
        "--cojo-collinear", collinearity, "collinearity", minimum=0, maximum=1
    )
    arguments += _optional_number("--diff-freq", diff_freq, "diff_freq", minimum=0, maximum=1)
    if isinstance(genomic_control, bool):
        arguments += ["--cojo-gc"] if genomic_control else []
    else:
        arguments += [
            "--cojo-gc",
            format(_cli.number(genomic_control, "genomic_control", minimum=1), ".12g"),
        ]
    arguments += _common(
        keep=keep, extract=extract, chromosome=chromosome, maf=maf, threads=threads
    )
    result = _run(
        arguments,
        base,
        required=required,
        optional=optional,
        gcta_path=gcta_path,
        timeout_seconds=timeout_seconds,
    )
    if "jma" in result["outputs"]:
        result["summary"] = _table(result["outputs"]["jma"], p_column="pJ")
    elif "cma" in result["outputs"]:
        result["summary"] = _table(result["outputs"]["cma"], p_column="pC")
    return result


def gcta_sblup(
    bfile: str,
    summary_file: str,
    out_prefix: str,
    shrinkage: float,
    *,
    window_kb: int = 200,
    keep: str | None = None,
    extract: str | None = None,
    chromosome: int | None = None,
    maf: float | None = None,
    threads: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summary-data-based BLUP risk-prediction SNP effects (``gcta64 --cojo-sblup``).

    ``shrinkage`` is lambda = m (1/h2 - 1), with m the number of SNPs and h2 the SNP-based
    heritability; ``window_kb`` (``--cojo-wind``, 200 here as in the GCTA example) is the distance
    beyond which SNPs are in linkage equilibrium. The SNP effects (SNP, A1, original effect, SBLUP
    effect) are written to ``<out_prefix>.sblup.cojo`` and can be used to score samples.
    """
    expected: dict[str, _Suffixes] = {"sblup": (".sblup.cojo", ".sblup")}
    base = _prefix(out_prefix, "out_prefix", expected)
    arguments = [
        "--bfile",
        str(_bfile(bfile)),
        "--cojo-file",
        str(_cli.existing_file(summary_file, "summary_file")),
    ]
    arguments += ["--cojo-sblup", format(_cli.number(shrinkage, "shrinkage", minimum=0), ".12g")]
    arguments += ["--cojo-wind", str(_cli.integer(window_kb, "window_kb", minimum=1))]
    arguments += _common(
        keep=keep, extract=extract, chromosome=chromosome, maf=maf, threads=threads
    )
    result = _run(
        arguments, base, required=expected, gcta_path=gcta_path, timeout_seconds=timeout_seconds
    )
    if "sblup" in result["outputs"]:
        result["n_snps"] = _count_lines(result["outputs"]["sblup"])
    return result


def gcta_mtcojo(
    summary_list_file: str,
    out_prefix: str,
    ref_ld_dir: str,
    w_ld_dir: str,
    *,
    bfile: str | None = None,
    mbfile: str | None = None,
    bxy_file: str | None = None,
    diff_freq: float | None = None,
    gwas_thresh: float | None = None,
    clump_r2: float | None = None,
    heidi_thresh: float | None = None,
    snp_min: int | None = None,
    ld_fdr: float | None = None,
    threads: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Multi-trait conditional and joint analysis: a GWAS of a target trait conditioned on other
    traits (``gcta64 --mtcojo-file``).

    ``summary_list_file`` lists the traits, one per row: name, summary file (COJO format) and
    optionally sample and population prevalence; the first row is the target trait, the others are
    covariate traits. ``ref_ld_dir`` and ``w_ld_dir`` are directories of LDSC-format LD scores
    (``1.l2.ldscore.gz`` ... ``22.l2.M_5_50``) for the LD-score regression and its weights. The
    effect of each covariate trait on the target (bxy) is estimated by GSMR, or taken from
    ``bxy_file`` (rows ``trait bxy``). The conditioned statistics are written to
    ``<out_prefix>.mtcojo.cma`` (summarised in ``summary``); ``diff_freq``, ``gwas_thresh``,
    ``clump_r2``, ``heidi_thresh``, ``snp_min`` and ``ld_fdr`` are the GSMR options of the same
    names.
    """
    expected: dict[str, _Suffixes] = {"mtcojo": (".mtcojo.cma",)}
    base = _prefix(out_prefix, "out_prefix", expected)
    arguments = _genotypes(bfile, mbfile)
    arguments += ["--mtcojo-file", str(_cli.existing_file(summary_list_file, "summary_list_file"))]
    arguments += [
        "--ref-ld-chr",
        str(_ld_dir(ref_ld_dir, "ref_ld_dir")),
        "--w-ld-chr",
        str(_ld_dir(w_ld_dir, "w_ld_dir")),
    ]
    arguments += _optional_file("--mtcojo-bxy", bxy_file, "bxy_file")
    arguments += _optional_number("--diff-freq", diff_freq, "diff_freq", minimum=0, maximum=1)
    arguments += _optional_number("--gwas-thresh", gwas_thresh, "gwas_thresh", minimum=0, maximum=1)
    arguments += _optional_number("--clump-r2", clump_r2, "clump_r2", minimum=0, maximum=1)
    arguments += _optional_number(
        "--heidi-thresh", heidi_thresh, "heidi_thresh", minimum=0, maximum=1
    )
    arguments += _optional_int("--gsmr-snp-min", snp_min, "snp_min", minimum=1)
    arguments += _optional_number("--gsmr-ld-fdr", ld_fdr, "ld_fdr", minimum=0, maximum=1)
    arguments += _optional_int("--thread-num", threads, "threads", minimum=1)
    result = _run(
        arguments, base, required=expected, gcta_path=gcta_path, timeout_seconds=timeout_seconds
    )
    return _summarise(result, "mtcojo", p_column="bC_pval")


def _ld_dir(value: str | Path, name: str) -> str:
    """An LD score directory as GCTA wants it: with a trailing slash."""
    path = Path(value).expanduser().resolve()
    if not path.is_dir():
        raise ConfigurationError(f"{name} must be an existing directory: {path}")
    return f"{path}/"


def gcta_fastbat(
    bfile: str,
    summary_file: str,
    out_prefix: str,
    *,
    gene_list_file: str | None = None,
    set_list_file: str | None = None,
    segment_kb: int | None = None,
    gene_window_kb: int | None = None,
    ld_cutoff: float | None = None,
    diff_freq: float | None = None,
    keep: str | None = None,
    extract: str | None = None,
    chromosome: int | None = None,
    maf: float | None = None,
    threads: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Set-based association test from GWAS summary statistics and an LD reference (``gcta64
    --fastBAT``).

    Exactly one set definition is needed: ``gene_list_file`` (tests genes; columns chromosome,
    start, end, gene id; ``gene_window_kb`` widens each gene, GCTA's default is 50 kb),
    ``set_list_file`` (named SNP sets, each ending with a line ``END``) or ``segment_kb`` (fixed
    windows of that many kb). ``ld_cutoff`` is the R-squared above which SNPs are pruned (default
    0.9; 1 turns pruning off). The result is ``<out_prefix>.gene.fastbat``, ``.fastbat`` or
    ``.seg.fastbat`` respectively, summarised with its smallest-p rows.
    """
    return _fastbat(
        bfile, summary_file, out_prefix, gene_list_file, set_list_file, segment_kb,
        gene_window_kb, ld_cutoff, diff_freq,
        False, keep, extract, chromosome, maf, threads, gcta_path, timeout_seconds,
    )  # fmt: skip


def _fastbat(
    bfile: str,
    summary_file: str,
    out_prefix: str,
    gene_list_file: str | None,
    set_list_file: str | None,
    segment_kb: int | None,
    gene_window_kb: int | None,
    ld_cutoff: float | None,
    diff_freq: float | None,
    write_snpset: bool,
    keep: str | None,
    extract: str | None,
    chromosome: int | None,
    maf: float | None,
    threads: int | None,
    gcta_path: str | None,
    timeout_seconds: float,
) -> dict[str, Any]:
    if sum(value is not None for value in (gene_list_file, set_list_file, segment_kb)) != 1:
        raise ConfigurationError(
            "Give exactly one of gene_list_file, set_list_file and segment_kb."
        )
    if gene_window_kb is not None and gene_list_file is None:
        raise ConfigurationError("gene_window_kb only applies with gene_list_file.")
    if gene_list_file is not None:
        kind, suffix = "gene", ".gene.fastbat"
    elif set_list_file is not None:
        kind, suffix = "set", ".fastbat"
    else:
        kind, suffix = "seg", ".seg.fastbat"
    expected: dict[str, _Suffixes] = {"fastbat": (suffix,)}
    snpset_group: _Suffixes = (
        (f".{kind}.snpset", ".snpset") if kind != "set" else (".snpset", ".fastbat.snpset")
    )
    optional: dict[str, _Suffixes] = {"snpset": snpset_group}
    base = _prefix(out_prefix, "out_prefix", {**expected, **optional})
    arguments = [
        "--bfile",
        str(_bfile(bfile)),
        "--fastBAT",
        str(_cli.existing_file(summary_file, "summary_file")),
    ]
    if gene_list_file is not None:
        arguments += [
            "--fastBAT-gene-list",
            str(_cli.existing_file(gene_list_file, "gene_list_file")),
        ]
        arguments += _optional_int("--fastBAT-wind", gene_window_kb, "gene_window_kb", minimum=0)
    elif set_list_file is not None:
        arguments += ["--fastBAT-set-list", str(_cli.existing_file(set_list_file, "set_list_file"))]
    else:
        arguments += ["--fastBAT-seg", str(_cli.integer(segment_kb, "segment_kb", minimum=1))]
    arguments += _optional_number(
        "--fastBAT-ld-cutoff", ld_cutoff, "ld_cutoff", minimum=0, maximum=1
    )
    arguments += _optional_number("--diff-freq", diff_freq, "diff_freq", minimum=0, maximum=1)
    arguments += ["--fastBAT-write-snpset"] if write_snpset else []
    arguments += _common(
        keep=keep, extract=extract, chromosome=chromosome, maf=maf, threads=threads
    )
    required = expected if not write_snpset else {**expected, **optional}
    result = _run(
        arguments,
        base,
        required=required,
        optional=optional,
        gcta_path=gcta_path,
        timeout_seconds=timeout_seconds,
    )
    return _summarise(result, "fastbat", p_column="Pvalue")


def gcta_mbat_combo(
    bfile: str,
    summary_file: str,
    gene_list_file: str,
    out_prefix: str,
    *,
    gene_window_kb: int | None = None,
    svd_gamma: float | None = None,
    ld_cutoff: float | None = None,
    print_all_p: bool = False,
    diff_freq: float | None = None,
    keep: str | None = None,
    extract: str | None = None,
    chromosome: int | None = None,
    maf: float | None = None,
    threads: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """mBAT-combo gene-based test that stays powerful when SNP effects mask each other (``gcta64
    --mBAT-combo``).

    Combines the mBAT test (principal components of the LD matrix explaining at least ``svd_gamma``
    of the variance, default 0.9) with the fastBAT test (``ld_cutoff`` pruning, default 0.9).
    ``gene_list_file`` has the columns chromosome, start, end and gene id; ``gene_window_kb`` widens
    each gene (default 50 kb). With ``print_all_p`` the table ``<out_prefix>.gene.assoc.mbat`` also
    holds the separate mBAT and fastBAT p-values.
    """
    expected: dict[str, _Suffixes] = {"mbat": (".gene.assoc.mbat",)}
    base = _prefix(out_prefix, "out_prefix", {**expected, "snpset": (".gene.snpset.mbat",)})
    arguments = [
        "--bfile",
        str(_bfile(bfile)),
        "--mBAT-combo",
        str(_cli.existing_file(summary_file, "summary_file")),
    ]
    arguments += ["--mBAT-gene-list", str(_cli.existing_file(gene_list_file, "gene_list_file"))]
    arguments += _optional_int("--mBAT-wind", gene_window_kb, "gene_window_kb", minimum=0)
    arguments += _optional_number("--mBAT-svd-gamma", svd_gamma, "svd_gamma", minimum=0, maximum=1)
    arguments += _optional_number(
        "--fastBAT-ld-cutoff", ld_cutoff, "ld_cutoff", minimum=0, maximum=1
    )
    arguments += ["--mBAT-print-all-p"] if _flag(print_all_p, "print_all_p") else []
    arguments += _optional_number("--diff-freq", diff_freq, "diff_freq", minimum=0, maximum=1)
    arguments += _common(
        keep=keep, extract=extract, chromosome=chromosome, maf=maf, threads=threads
    )
    result = _run(
        arguments, base, required=expected, gcta_path=gcta_path, timeout_seconds=timeout_seconds
    )
    return _summarise(result, "mbat", p_column="P_mBATcombo")


def gcta_write_gene_snpsets(
    bfile: str,
    summary_file: str,
    gene_list_file: str,
    out_prefix: str,
    *,
    method: str = "fastbat",
    gene_window_kb: int | None = None,
    diff_freq: float | None = None,
    keep: str | None = None,
    extract: str | None = None,
    chromosome: int | None = None,
    maf: float | None = None,
    threads: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build and write the SNP set of every gene used by a gene-based test (``gcta64
    --fastBAT-write-snpset`` / ``--mBAT-write-snpset``).

    GCTA writes the sets as a by-product of the gene test, so the test itself runs as well and its
    table is returned in ``outputs``. ``method='fastbat'`` (``--fastBAT --fastBAT-write-snpset``)
    writes ``<out_prefix>.gene.snpset`` (the format ``--fastBAT-set-list`` reads); ``method='mbat'``
    (``--mBAT-combo --mBAT-write-snpset``) writes ``<out_prefix>.gene.snpset.mbat`` (the same
    layout: a gene name, its SNPs, ``END``). ``snpset`` is the path of the set file and ``genes``
    the number of genes with a set.
    """
    kind = _cli.choice(method, ("fastbat", "mbat"), "method")
    if kind == "fastbat":
        result = _fastbat(
            bfile, summary_file, out_prefix, gene_list_file, None, None, gene_window_kb, None,
            diff_freq,
            True, keep, extract, chromosome, maf, threads, gcta_path, timeout_seconds,
        )  # fmt: skip
        key = "snpset"
    else:
        expected: dict[str, _Suffixes] = {"snpset": (".gene.snpset.mbat",)}
        base = _prefix(out_prefix, "out_prefix", {**expected, "mbat": (".gene.assoc.mbat",)})
        arguments = [
            "--bfile",
            str(_bfile(bfile)),
            "--mBAT-combo",
            str(_cli.existing_file(summary_file, "summary_file")),
        ]
        arguments += [
            "--mBAT-gene-list",
            str(_cli.existing_file(gene_list_file, "gene_list_file")),
            "--mBAT-write-snpset",
        ]
        arguments += _optional_int("--mBAT-wind", gene_window_kb, "gene_window_kb", minimum=0)
        arguments += _optional_number("--diff-freq", diff_freq, "diff_freq", minimum=0, maximum=1)
        arguments += _common(
            keep=keep, extract=extract, chromosome=chromosome, maf=maf, threads=threads
        )
        result = _run(
            arguments, base, required=expected, optional={"mbat": (".gene.assoc.mbat",)},
            gcta_path=gcta_path,
            timeout_seconds=timeout_seconds,
        )  # fmt: skip
        key = "snpset"
    if key in result["outputs"]:
        result["snpset"] = result["outputs"][key]
        text = Path(result["outputs"][key]).read_text(encoding="utf-8", errors="replace")
        result["genes"] = sum(1 for line in text.splitlines() if line.strip() == "END")
    return result


def gcta_gsmr(
    exposure_list_file: str,
    outcome_list_file: str,
    out_prefix: str,
    *,
    bfile: str | None = None,
    mbfile: str | None = None,
    direction: int = 0,
    gsmr2: bool = False,
    diff_freq: float | None = None,
    gwas_thresh: float | None = None,
    clump_r2: float | None = None,
    heidi_thresh: float | None = None,
    snp_min: int | None = None,
    ld_fdr: float | None = None,
    effect_plot: bool = False,
    threads: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Generalised summary-data-based Mendelian randomisation between exposures and outcomes
    (``gcta64 --gsmr-file``).

    ``exposure_list_file`` and ``outcome_list_file`` list traits (name and summary-statistics file
    in COJO format, one per row). ``direction`` is 0 (forward, default), 1 (reverse) or 2
    (bidirectional). The SNP instruments are clumped at ``gwas_thresh`` (default 5e-8) and
    ``clump_r2`` (default 0.05); ``heidi_thresh`` (default 0.01) removes pleiotropic outliers and
    ``gsmr2=True`` uses the new global HEIDI-outlier method of GSMR2; ``snp_min`` (default 10) is
    the minimum number of instruments and ``ld_fdr`` (default 0.05) shrinks chance LD between
    instruments. ``<out_prefix>.gsmr`` holds exposure, outcome, bxy, se, p and the number of SNPs
    (``summary``); ``effect_plot`` also saves the data for the SNP-effect plot (``.eff_plot.gz``).
    """
    expected: dict[str, _Suffixes] = {"gsmr": (".gsmr",)}
    optional: dict[str, _Suffixes] = {"effect_plot": (".eff_plot.gz",)}
    base = _prefix(out_prefix, "out_prefix", {**expected, **optional})
    way = _cli.integer(direction, "direction", minimum=0, maximum=2)
    arguments = _genotypes(bfile, mbfile)
    arguments += [
        "--gsmr-file",
        str(_cli.existing_file(exposure_list_file, "exposure_list_file")),
        str(_cli.existing_file(outcome_list_file, "outcome_list_file")),
    ]
    arguments += ["--gsmr-direction", str(way)]
    arguments += ["--gsmr2-beta"] if _flag(gsmr2, "gsmr2") else []
    arguments += _optional_number("--diff-freq", diff_freq, "diff_freq", minimum=0, maximum=1)
    arguments += _optional_number("--gwas-thresh", gwas_thresh, "gwas_thresh", minimum=0, maximum=1)
    arguments += _optional_number("--clump-r2", clump_r2, "clump_r2", minimum=0, maximum=1)
    arguments += _optional_number(
        "--heidi-thresh", heidi_thresh, "heidi_thresh", minimum=0, maximum=1
    )
    arguments += _optional_int("--gsmr-snp-min", snp_min, "snp_min", minimum=1)
    arguments += _optional_number("--gsmr-ld-fdr", ld_fdr, "ld_fdr", minimum=0, maximum=1)
    arguments += ["--effect-plot"] if _flag(effect_plot, "effect_plot") else []
    arguments += _optional_int("--thread-num", threads, "threads", minimum=1)
    result = _run(
        arguments,
        base,
        required=expected,
        optional=optional,
        gcta_path=gcta_path,
        timeout_seconds=timeout_seconds,
    )
    return _summarise(result, "gsmr", p_column="p")


# ---------------------------------------------------------------------------
# Linkage disequilibrium, PCA loadings and population genetics
# ---------------------------------------------------------------------------


def gcta_ld_friends(
    bfile: str,
    snp_list_file: str,
    out_prefix: str,
    *,
    window_kb: int | None = None,
    significance: float | None = None,
    keep: str | None = None,
    extract: str | None = None,
    threads: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find the SNPs in significant LD with each target SNP by simple regression (``gcta64 --ld``).

    ``snp_list_file`` names the target SNPs; ``window_kb`` is the search distance in either
    direction and ``significance`` the p-value threshold of the regression test.
    ``<out_prefix>.rsq.ld`` summarises the LD block of each target (flanking SNPs, number of SNPs,
    mean, median and maximum r2, best tagging SNP), ``.r.ld`` holds the correlations and ``.snp.ld``
    the names of the SNPs in the block.
    """
    expected: dict[str, _Suffixes] = {"rsq": (".rsq.ld",), "r": (".r.ld",), "snp": (".snp.ld",)}
    base = _prefix(out_prefix, "out_prefix", expected)
    arguments = [
        "--bfile",
        str(_bfile(bfile)),
        "--ld",
        str(_cli.existing_file(snp_list_file, "snp_list_file")),
    ]
    arguments += _optional_int("--ld-wind", window_kb, "window_kb", minimum=1)
    arguments += _optional_number("--ld-sig", significance, "significance", minimum=0, maximum=1)
    arguments += _common(keep=keep, extract=extract, threads=threads)
    result = _run(
        arguments, base, required=expected, gcta_path=gcta_path, timeout_seconds=timeout_seconds
    )
    return _summarise(result, "rsq")


def gcta_segment_ld_score(
    bfile: str,
    out_prefix: str,
    *,
    segment_kb: int = 200,
    window_kb: int | None = None,
    rsq_cutoff: float | None = None,
    adjust: bool = False,
    keep: str | None = None,
    extract: str | None = None,
    maf: float | None = None,
    threads: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Segment-based LD scores of autosomal SNPs for stratifying SNPs in GREML-LDMS (``gcta64
    --autosome --ld-score-region``).

    Segments are ``segment_kb`` long (default 200) and overlap by half. ``window_kb``
    (``--ld-wind``) is the window used for each SNP's LD score and ``rsq_cutoff``
    (``--ld-rsq-cutoff``) the r2 below which pairs are ignored (not recommended for LD score
    regression); ``adjust=True`` applies the unbiased r2 adjustment (``--ld-score-adj``). The table
    ``<out_prefix>.score.ld`` has SNP, chr, bp, freq, mean_rsq, snp_num, max_rsq, ldscore_SNP and
    ldscore_region.
    """
    expected: dict[str, _Suffixes] = {"score_ld": (".score.ld",)}
    base = _prefix(out_prefix, "out_prefix", expected)
    arguments = [
        "--bfile",
        str(_bfile(bfile)),
        "--autosome",
        "--ld-score-region",
        str(_cli.integer(segment_kb, "segment_kb", minimum=1)),
    ]
    arguments += _optional_int("--ld-wind", window_kb, "window_kb", minimum=1)
    arguments += _optional_number("--ld-rsq-cutoff", rsq_cutoff, "rsq_cutoff", minimum=0, maximum=1)
    arguments += ["--ld-score-adj"] if _flag(adjust, "adjust") else []
    arguments += _common(keep=keep, extract=extract, maf=maf, threads=threads)
    result = _run(
        arguments, base, required=expected, gcta_path=gcta_path, timeout_seconds=timeout_seconds
    )
    return _summarise(result, "score_ld")


def gcta_multi_set_ld_score(
    bfile: str,
    snp_set_files: Sequence[str],
    out_prefix: str,
    *,
    window_kb: int | None = None,
    rsq_cutoff: float | None = None,
    adjust: bool = False,
    keep: str | None = None,
    extract: str | None = None,
    maf: float | None = None,
    threads: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """LD score of every autosomal SNP against each of several SNP sets, for multi-component LD
    score regression (``gcta64 --autosome --ld-score-multi``).

    ``snp_set_files`` are files of SNP names, one per set (for example the SNPs in exons and in
    introns). The table ``<out_prefix>.multi_score.ld`` (``.score.ld`` in GCTA 1.94) has, per set,
    mean_rsq, snp_num, max_rsq and the LD score against that set. ``window_kb`` (``--ld-wind``),
    ``rsq_cutoff`` (``--ld-rsq-cutoff``) and ``adjust`` (``--ld-score-adj``) are as for
    ``gcta_segment_ld_score``.
    """
    if isinstance(snp_set_files, (str, bytes)) or not snp_set_files:
        raise ConfigurationError("snp_set_files must be a non-empty list of SNP list files.")
    sets = [str(_cli.existing_file(path, "snp_set_files")) for path in snp_set_files]
    expected: dict[str, _Suffixes] = {"score_ld": (".multi_score.ld", ".score.ld")}
    base = _prefix(out_prefix, "out_prefix", expected)
    with tempfile.TemporaryDirectory(prefix="dnakit-gcta-") as folder:
        listing = _list_file(Path(folder), "sets.txt", sets)
        arguments = ["--bfile", str(_bfile(bfile)), "--autosome", "--ld-score-multi", str(listing)]
        arguments += _optional_int("--ld-wind", window_kb, "window_kb", minimum=1)
        arguments += _optional_number(
            "--ld-rsq-cutoff", rsq_cutoff, "rsq_cutoff", minimum=0, maximum=1
        )
        arguments += ["--ld-score-adj"] if _flag(adjust, "adjust") else []
        arguments += _common(keep=keep, extract=extract, maf=maf, threads=threads)
        result = _run(
            arguments, base, required=expected, gcta_path=gcta_path, timeout_seconds=timeout_seconds
        )
    result["n_sets"] = len(sets)
    return _summarise(result, "score_ld")


def gcta_pc_loading(
    bfile: str,
    pca_prefix: str,
    out_prefix: str,
    *,
    autosome: bool = True,
    keep: str | None = None,
    extract: str | None = None,
    maf: float | None = None,
    threads: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """SNP loadings on the principal components computed in a reference sample (``gcta64
    --pc-loading``).

    ``pca_prefix`` names the ``<prefix>.eigenvec`` and ``<prefix>.eigenval`` files of a GCTA
    ``--pca`` run on the GRM of the same individuals. The loadings are written to
    ``<out_prefix>.pcl`` and are projected onto a target sample with ``gcta_project_loading``.
    ``autosome`` (default true) restricts the SNPs to chromosomes 1 to 22 (the analysis fails for
    sex chromosomes).
    """
    pca = Path(pca_prefix).expanduser().resolve()
    for suffix in (".eigenvec", ".eigenval"):
        if not Path(f"{pca}{suffix}").is_file():
            raise ConfigurationError(f"pca_prefix needs {pca}{suffix} (from gcta64 --pca).")
    expected: dict[str, _Suffixes] = {"pcl": (".pcl",)}
    base = _prefix(out_prefix, "out_prefix", expected)
    arguments = (
        ["--bfile", str(_bfile(bfile))]
        + (["--autosome"] if _flag(autosome, "autosome") else [])
        + ["--pc-loading", str(pca)]
    )
    arguments += _common(keep=keep, extract=extract, maf=maf, threads=threads)
    result = _run(
        arguments, base, required=expected, gcta_path=gcta_path, timeout_seconds=timeout_seconds
    )
    if "pcl" in result["outputs"]:
        result["n_snps"] = _count_lines(result["outputs"]["pcl"]) - 1
    return result


def gcta_project_loading(
    bfile: str,
    loading_prefix: str,
    out_prefix: str,
    *,
    n_components: int = 20,
    autosome: bool = True,
    update_freq_file: str | None = None,
    keep: str | None = None,
    extract: str | None = None,
    threads: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Project the samples of a fileset onto reference SNP loadings to obtain their principal
    components (``gcta64 --project-loading``).

    ``loading_prefix`` is the prefix of the ``<prefix>.pcl`` file made by ``gcta_pc_loading``; only
    the SNPs shared with the reference are used. ``n_components`` PCs are written (default 20, at
    most the number computed for the reference) to ``<out_prefix>.proj.eigenvec`` (FID, IID, PC1
    ...). By default loadings use the reference allele frequencies; ``update_freq_file``
    (``--update-freq``) supplies other frequencies.
    """
    loadings = Path(loading_prefix).expanduser().resolve()
    if not Path(f"{loadings}.pcl").is_file():
        raise ConfigurationError(f"loading_prefix needs {loadings}.pcl (from gcta_pc_loading).")
    expected: dict[str, _Suffixes] = {"eigenvec": (".proj.eigenvec",)}
    base = _prefix(out_prefix, "out_prefix", expected)
    count = _cli.integer(n_components, "n_components", minimum=1)
    arguments = ["--bfile", str(_bfile(bfile))] + (
        ["--autosome"] if _flag(autosome, "autosome") else []
    )
    arguments += ["--project-loading", str(loadings), str(count)]
    arguments += _optional_file("--update-freq", update_freq_file, "update_freq_file")
    arguments += _common(keep=keep, extract=extract, threads=threads)
    result = _run(
        arguments, base, required=expected, gcta_path=gcta_path, timeout_seconds=timeout_seconds
    )
    if "eigenvec" in result["outputs"]:
        result["n_individuals"] = _count_lines(result["outputs"]["eigenvec"])
    return result


def gcta_fst(
    bfile: str,
    subpopulation_file: str,
    out_prefix: str,
    *,
    keep: str | None = None,
    extract: str | None = None,
    chromosome: int | None = None,
    maf: float | None = None,
    threads: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fst of every SNP between sub-populations, following Weir (1996) (``gcta64 --fst``).

    ``subpopulation_file`` has the columns FID, IID and a population label. The table
    ``<out_prefix>.fst`` gives per SNP the chromosome, position, reference allele, the allele
    frequency in each population and Fst; ``summary`` holds its columns and first rows.
    """
    expected: dict[str, _Suffixes] = {"fst": (".fst",)}
    base = _prefix(out_prefix, "out_prefix", expected)
    arguments = [
        "--bfile",
        str(_bfile(bfile)),
        "--fst",
        "--sub-popu",
        str(_cli.existing_file(subpopulation_file, "subpopulation_file")),
    ]
    arguments += _common(
        keep=keep, extract=extract, chromosome=chromosome, maf=maf, threads=threads
    )
    result = _run(
        arguments, base, required=expected, gcta_path=gcta_path, timeout_seconds=timeout_seconds
    )
    if "fst" in result["outputs"]:
        table = _table(result["outputs"]["fst"])
        fst_values = [
            float(line.split()[-1])
            for line in Path(result["outputs"]["fst"]).read_text(encoding="utf-8").splitlines()[1:]
            if line.strip()
        ]
        table["mean_fst"] = sum(fst_values) / len(fst_values) if fst_values else None
        result["summary"] = table
    return result


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------


def gcta_simulate_phenotypes(
    bfile: str,
    out_prefix: str,
    *,
    trait: str = "quantitative",
    n_cases: int | None = None,
    n_controls: int | None = None,
    causal_loci_file: str | None = None,
    heritability: float | None = None,
    prevalence: float | None = None,
    replicates: int | None = None,
    keep: str | None = None,
    extract: str | None = None,
    threads: int | None = None,
    gcta_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Simulate phenotypes from real genotypes with chosen causal variants (``gcta64 --simu-qt`` /
    ``--simu-cc``).

    ``trait='quantitative'`` simulates a quantitative trait; ``trait='case_control'`` simulates
    ``n_cases`` cases and ``n_controls`` controls from the liability model (at most n*K cases and
    n*(1-K) controls, K = ``prevalence``, default 0.1). ``causal_loci_file`` lists causal SNPs with
    optional effect sizes (columns SNP, effect; effects are drawn from a standard normal when
    absent, and so are the loci when the file is omitted). ``heritability`` (default 0.1; of the
    liability for case-control) and ``replicates`` (default 1) shape the simulation. GCTA draws new
    random numbers on every run (its ``--simu-seed`` does not make runs repeat), so results differ
    between calls. ``<out_prefix>.phen`` holds FID, IID and one column per replicate (missing -9 for
    individuals that were not sampled in a case-control study); ``<out_prefix>.par`` holds the
    causal variants (QTL, RefAllele, Frequency, Effect).
    """
    kind = _cli.choice(trait, ("quantitative", "case_control"), "trait")
    if kind == "case_control":
        if n_cases is None or n_controls is None:
            raise ConfigurationError("trait='case_control' needs n_cases and n_controls.")
    elif n_cases is not None or n_controls is not None or prevalence is not None:
        raise ConfigurationError(
            "n_cases, n_controls and prevalence only apply to trait='case_control'."
        )
    expected: dict[str, _Suffixes] = {"phen": (".phen",), "par": (".par",)}
    base = _prefix(out_prefix, "out_prefix", expected)
    arguments = ["--bfile", str(_bfile(bfile))]
    if kind == "quantitative":
        arguments += ["--simu-qt"]
    else:
        arguments += [
            "--simu-cc",
            str(_cli.integer(n_cases, "n_cases", minimum=1)),
            str(_cli.integer(n_controls, "n_controls", minimum=1)),
        ]
    arguments += _optional_file("--simu-causal-loci", causal_loci_file, "causal_loci_file")
    arguments += _optional_number("--simu-hsq", heritability, "heritability", minimum=0, maximum=1)
    arguments += _optional_number("--simu-k", prevalence, "prevalence", minimum=0, maximum=1)
    arguments += _optional_int("--simu-rep", replicates, "replicates", minimum=1)
    arguments += _common(keep=keep, extract=extract, threads=threads)
    result = _run(
        arguments, base, required=expected, gcta_path=gcta_path, timeout_seconds=timeout_seconds
    )
    if "par" in result["outputs"]:
        result["causal_variants"] = _table(result["outputs"]["par"], top=1000).get("first_rows", [])
    if "phen" in result["outputs"]:
        result["n_individuals"] = _count_lines(result["outputs"]["phen"])
    return result


__all__ = [
    "gcta_acat",
    "gcta_blup_snp",
    "gcta_cojo",
    "gcta_fastbat",
    "gcta_fastgwa",
    "gcta_fastgwa_burden",
    "gcta_fst",
    "gcta_gsmr",
    "gcta_he_regression",
    "gcta_inbreeding_coefficient",
    "gcta_ld_friends",
    "gcta_make_grm",
    "gcta_make_sparse_grm",
    "gcta_mbat_combo",
    "gcta_mlma",
    "gcta_mtcojo",
    "gcta_multi_set_ld_score",
    "gcta_pc_loading",
    "gcta_project_loading",
    "gcta_recode_standardized_genotypes",
    "gcta_reml",
    "gcta_reml_bivar",
    "gcta_sblup",
    "gcta_segment_ld_score",
    "gcta_simulate_phenotypes",
    "gcta_unify_grm",
    "gcta_write_gene_snpsets",
]
_register_public_verbose(globals())
