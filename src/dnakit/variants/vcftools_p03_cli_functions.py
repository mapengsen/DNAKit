"""VCFtools command-line facades: the ``vcftools`` program and the VCFtools Perl scripts (``vcf-*``, ``fill-*``).

The user installs VCFtools (``conda install -c bioconda vcftools``, ``apt install vcftools`` or a build from
https://github.com/vcftools/vcftools; 0.1.16 is the tested version) and puts it on PATH or passes ``vcftools_path`` or
``script_path``. The PCA functions need a ``vcftools`` that was configured with ``--enable-pca`` (LAPACK). The ``vcf-*``
Perl scripts also need their ``Vcf.pm`` library; it is found next to the installation (``../share/perl*``) when it is
not in ``PERL5LIB`` already. Scripts that read a bgzip-compressed file call ``tabix`` or ``samtools`` from PATH.

Each public function runs one fixed program with an explicit argument list built from typed arguments: no shell, no
generic runner. Statistics are written next to ``output_prefix`` as ``<output_prefix>.<suffix>`` (the file names
VCFtools itself uses) and an existing output is never overwritten. Perl scripts that write to standard output are
redirected to a new ``output_path``. Results are JSON-serialisable dictionaries: status, backend, version, command,
return code, existing outputs, a bounded log tail and, for tabular outputs, a preview of the table.

Defaults are the defaults of the VCFtools options they map to: an option is only passed when it is given.
"""

from __future__ import annotations

import contextlib
import math
import os
import re
import shutil
import signal
import subprocess
import tempfile
import time
from collections.abc import Mapping, Sequence
from functools import lru_cache
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.backends.external import execute_bounded_command
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.exceptions import (
    BackendExecutionError,
    BackendTimeoutError,
    BackendUnavailableError,
    ConfigurationError,
    DNAKitError,
)

_TESTED_VERSION = "0.1.16"
_VCFTOOLS_HINT = (
    "Install VCFtools 0.1.16 or later (conda install -c bioconda vcftools, apt install vcftools, or build "
    "https://github.com/vcftools/vcftools; configure with --enable-pca for the PCA functions), then put it on PATH "
    "or pass vcftools_path."
)
_SCRIPT_HINT = (
    "Install VCFtools 0.1.16 or later with its Perl scripts (conda install -c bioconda vcftools, apt install "
    "vcftools, or build https://github.com/vcftools/vcftools), then put the script on PATH or pass script_path. "
    "Scripts that read compressed VCF need tabix (htslib); fill-aa, fill-fs and fill-ref-md5 need samtools."
)
_PCA_HINT = (
    "This vcftools was built without PCA support. Rebuild VCFtools with ./configure --enable-pca (needs LAPACK), "
    "or point vcftools_path to such a build."
)
_NAME_TOKEN = re.compile(r"[A-Za-z0-9_][A-Za-z0-9_.]*")
_SAMPLE = re.compile(r"[^\s,/\x00]+")
_KEPT = re.compile(r"After filtering, kept (\d+) out of (?:a possible )?(\d+) (Individuals|Sites)")
_FILTER_SPEC = re.compile(r"[A-Za-z0-9_+=./,-]+")
_REMOVE_TAG = re.compile(r"[A-Za-z0-9_.]+(?:/[A-Za-z0-9_.]+)?")
_TIMEOUT = 3600.0
_MAX_ARTIFACT_BYTES = 2_000_000_000
_LOG_TAIL_CHARS = 4000
_MAX_ROWS = 1000

# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------


def _sequence_of_text(value: object, name: str, *, pattern: re.Pattern[str] | None = None) -> list[str]:
    """One text or a list of texts, each non-empty, not starting with '-' and optionally matching ``pattern``."""
    items = [value] if isinstance(value, str) else value
    if isinstance(items, (str, bytes)) or not isinstance(items, Sequence) or not items:
        raise ConfigurationError(f"{name} must be a non-empty string or list of strings.")
    cleaned = [_cli.text(item, name) for item in items]
    if pattern is not None:
        for item in cleaned:
            if pattern.fullmatch(item) is None:
                raise ConfigurationError(f"{name} contains an unsupported value: {item!r}")
    return cleaned


def _files(value: object, name: str) -> list[str]:
    items = [value] if isinstance(value, (str, Path)) else value
    if isinstance(items, (str, bytes)) or not isinstance(items, Sequence) or not items:
        raise ConfigurationError(f"{name} must be a path or a non-empty list of paths.")
    return [str(_cli.existing_file(item, name)) for item in items]


def _switch(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be true or false.")
    return value


def _name_token(value: object, name: str) -> str:
    text = _cli.text(value, name)
    if _NAME_TOKEN.fullmatch(text) is None:
        raise ConfigurationError(f"{name} may contain only letters, digits, '_' and '.': {text!r}")
    return text


def _input_flag(path: str | Path, *, second: bool = False) -> list[str]:
    """``--vcf``/``--gzvcf``/``--bcf`` (or ``--diff``/``--gzdiff``/``--diff-bcf``) chosen from the file name."""
    resolved = _cli.existing_file(path, "second_vcf_path" if second else "vcf_path")
    name = resolved.name.lower()
    if name.endswith(".bcf"):
        flag = "--diff-bcf" if second else "--bcf"
    elif name.endswith((".gz", ".bgz")):
        flag = "--gzdiff" if second else "--gzvcf"
    else:
        flag = "--diff" if second else "--vcf"
    return [flag, str(resolved)]


def _prefix(output_prefix: str | Path, suffixes: Sequence[str]) -> Path:
    """The ``--out`` prefix: its folder exists and neither the outputs nor the ``.log`` exist yet."""
    if isinstance(output_prefix, Path):
        output_prefix = str(output_prefix)
    text = _cli.text(output_prefix, "output_prefix")
    prefix = Path(text).expanduser().resolve()
    if not prefix.parent.is_dir():
        raise ConfigurationError(f"The folder of output_prefix does not exist: {prefix.parent}")
    if prefix.is_dir():
        raise ConfigurationError(f"output_prefix must be a file-name prefix, not a folder: {prefix}")
    for suffix in (*suffixes, ".log"):
        target = Path(str(prefix) + suffix)
        if target.exists():
            raise ConfigurationError(f"output_prefix would overwrite an existing file: {target}")
    return prefix


def _output_path(value: str | Path, name: str = "output_path") -> Path:
    return _cli.new_file(str(value) if isinstance(value, Path) else _cli.text(value, name), name)


def _jsonable(value: float) -> float | None:
    return None if isinstance(value, float) and not math.isfinite(value) else value


def _cell(text: str) -> Any:
    for convert in (int, float):
        try:
            return _jsonable(convert(text))
        except ValueError:
            continue
    return text


def _read_table(path: Path, max_rows: int) -> dict[str, Any]:
    """Preview of a tab-separated table with a header line: columns, the first ``max_rows`` rows and the row count."""
    columns: list[str] = []
    rows: list[list[Any]] = []
    count = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for index, line in enumerate(handle):
            fields = line.rstrip("\n").split("\t")
            if index == 0:
                columns = fields
                continue
            count += 1
            if len(rows) < max_rows:
                rows.append([_cell(field) for field in fields])
    return {"columns": columns, "row_count": count, "rows": rows, "truncated": count > len(rows)}


def _count_records(path: Path) -> int:
    count = 0
    with path.open("rb") as handle:
        for line in handle:
            if not line.startswith(b"#") and line.strip():
                count += 1
    return count


def _preview_rows(max_rows: object) -> int:
    return _cli.integer(max_rows, "max_rows", minimum=1, maximum=1_000_000)


# ---------------------------------------------------------------------------
# Site, individual and genotype filters shared by the vcftools functions
# ---------------------------------------------------------------------------

# key -> (vcftools flag, kind). The keys are the keyword names of ``vcftools_filter_vcf``.
_FILTER_SPECS: Mapping[str, tuple[str, str]] = {
    "chromosomes": ("--chr", "texts"),
    "exclude_chromosomes": ("--not-chr", "texts"),
    "from_bp": ("--from-bp", "int"),
    "to_bp": ("--to-bp", "int"),
    "positions_file": ("--positions", "file"),
    "exclude_positions_file": ("--exclude-positions", "file"),
    "positions_overlap_file": ("--positions-overlap", "file"),
    "exclude_positions_overlap_file": ("--exclude-positions-overlap", "file"),
    "bed_file": ("--bed", "file"),
    "exclude_bed_file": ("--exclude-bed", "file"),
    "thin": ("--thin", "int"),
    "mask_file": ("--mask", "file"),
    "invert_mask_file": ("--invert-mask", "file"),
    "mask_min": ("--mask-min", "int"),
    "snp_ids": ("--snp", "texts"),
    "snp_ids_file": ("--snps", "file"),
    "exclude_ids_file": ("--exclude", "file"),
    "only_indels": ("--keep-only-indels", "switch"),
    "remove_indels": ("--remove-indels", "switch"),
    "remove_filtered_all": ("--remove-filtered-all", "switch"),
    "keep_filtered": ("--keep-filtered", "texts"),
    "remove_filtered": ("--remove-filtered", "texts"),
    "keep_info": ("--keep-INFO", "texts"),
    "remove_info": ("--remove-INFO", "texts"),
    "min_maf": ("--maf", "float"),
    "max_maf": ("--max-maf", "float"),
    "min_mac": ("--mac", "int"),
    "max_mac": ("--max-mac", "int"),
    "min_alleles": ("--min-alleles", "int"),
    "max_alleles": ("--max-alleles", "int"),
    "min_mean_depth": ("--min-meanDP", "float"),
    "max_mean_depth": ("--max-meanDP", "float"),
    "hwe_pvalue": ("--hwe", "float"),
    "phased_only": ("--phased", "switch"),
    "min_quality": ("--minQ", "float"),
    "max_missing": ("--max-missing", "float"),
    "keep_individuals": ("--indv", "texts"),
    "keep_individuals_files": ("--keep", "files"),
    "remove_individuals": ("--remove-indv", "texts"),
    "remove_individuals_files": ("--remove", "files"),
    "max_individuals": ("--max-indv", "int"),
    "remove_filtered_genotypes": ("--remove-filtered-geno-all", "switch"),
    "remove_filtered_genotype_flags": ("--remove-filtered-geno", "texts"),
    "min_genotype_quality": ("--minGQ", "float"),
    "min_genotype_depth": ("--minDP", "int"),
    "max_genotype_depth": ("--maxDP", "int"),
}
_EXCLUSIVE_FILTERS = (("bed_file", "exclude_bed_file"), ("mask_file", "invert_mask_file"))


def _filter_arguments(filters: Mapping[str, object] | None) -> list[str]:
    """Arguments for a ``filters`` mapping (keys of ``vcftools_filter_vcf``); ``None``, False and empty values are skipped."""
    if filters is None:
        return []
    if not isinstance(filters, Mapping):
        raise ConfigurationError("filters must be a mapping of filter names to values.")
    unknown = sorted(set(filters) - set(_FILTER_SPECS))
    if unknown:
        raise ConfigurationError(f"Unknown filters {unknown}; allowed: {sorted(_FILTER_SPECS)}.")
    active = {key: value for key, value in filters.items() if value is not None and value is not False}
    for first, second in _EXCLUSIVE_FILTERS:
        if first in active and second in active:
            raise ConfigurationError(f"{first} and {second} cannot be used together.")
    arguments: list[str] = []
    for key, (flag, kind) in _FILTER_SPECS.items():
        if key not in active:
            continue
        value = active[key]
        if kind == "switch":
            if _switch(value, key):
                arguments.append(flag)
        elif kind == "texts":
            for item in _sequence_of_text(value, key):
                arguments += [flag, item]
        elif kind == "file":
            arguments += [flag, str(_cli.existing_file(str(value), key))]  # type: ignore[arg-type]
        elif kind == "files":
            for item in _files(value, key):
                arguments += [flag, item]
        elif kind == "int":
            arguments += [flag, str(_cli.integer(value, key, minimum=0))]
        else:
            arguments += [flag, format(_cli.number(value, key, minimum=0.0), ".12g")]
    return arguments


# ---------------------------------------------------------------------------
# Running the vcftools program
# ---------------------------------------------------------------------------


def _vcftools(
    arguments: Sequence[str],
    prefix: Path,
    suffixes: Sequence[str],
    *,
    vcftools_path: str | None,
    timeout_seconds: float,
    table_suffix: str | None = None,
    max_rows: int = _MAX_ROWS,
) -> dict[str, Any]:
    """Run ``vcftools <arguments> --out <prefix>``; ``completed`` only when exit code is 0 and every output exists."""
    expected = [Path(str(prefix) + suffix) for suffix in suffixes]
    result = _cli.run(
        "vcftools",
        [*arguments, "--out", str(prefix)],
        cwd=prefix.parent,
        expected=expected,
        explicit=vcftools_path,
        hint=_VCFTOOLS_HINT,
        timeout_seconds=_cli.number(timeout_seconds, "timeout_seconds", minimum=1.0),
    )
    result["output_prefix"] = str(prefix)
    log_file = Path(str(prefix) + ".log")
    if log_file.is_file():
        result["log_file"] = str(log_file)
    for kept, total, what in _KEPT.findall(result["log_tail"]):
        result[f"{what.lower()}_kept"] = int(kept)
        result[f"{what.lower()}_total"] = int(total)
    if "compiled without PCA" in result["log_tail"]:
        raise BackendUnavailableError("vcftools was compiled without PCA support.", hint=_PCA_HINT)
    if result["status"] == "completed" and table_suffix is not None:
        result["table"] = _read_table(Path(str(prefix) + table_suffix), max_rows)
    return result


def _report(
    vcf_path: str | Path,
    output_prefix: str | Path,
    analysis: Sequence[str],
    suffixes: Sequence[str],
    *,
    filters: Mapping[str, object] | None,
    vcftools_path: str | None,
    timeout_seconds: float,
    max_rows: int = _MAX_ROWS,
    table: bool = True,
    extra: Sequence[str] = (),
) -> dict[str, Any]:
    """Shared body of the statistics functions: input flag, filters, one analysis, optional table preview."""
    rows = _preview_rows(max_rows)
    input_arguments = _input_flag(vcf_path)
    filter_arguments = _filter_arguments(filters)
    prefix = _prefix(output_prefix, suffixes)
    return _vcftools(
        [*input_arguments, *filter_arguments, *extra, *analysis],
        prefix,
        suffixes,
        vcftools_path=vcftools_path,
        timeout_seconds=timeout_seconds,
        table_suffix=suffixes[0] if table else None,
        max_rows=rows,
    )


# ---------------------------------------------------------------------------
# Filtering and writing a VCF/BCF
# ---------------------------------------------------------------------------


def vcftools_filter_vcf(
    vcf_path: str,
    output_prefix: str,
    *,
    output_format: str = "vcf",
    recode_info: bool | Sequence[str] = False,
    chromosomes: Sequence[str] | None = None,
    exclude_chromosomes: Sequence[str] | None = None,
    from_bp: int | None = None,
    to_bp: int | None = None,
    positions_file: str | None = None,
    exclude_positions_file: str | None = None,
    positions_overlap_file: str | None = None,
    exclude_positions_overlap_file: str | None = None,
    bed_file: str | None = None,
    exclude_bed_file: str | None = None,
    thin: int | None = None,
    mask_file: str | None = None,
    invert_mask_file: str | None = None,
    mask_min: int | None = None,
    snp_ids: Sequence[str] | None = None,
    snp_ids_file: str | None = None,
    exclude_ids_file: str | None = None,
    only_indels: bool = False,
    remove_indels: bool = False,
    remove_filtered_all: bool = False,
    keep_filtered: Sequence[str] | None = None,
    remove_filtered: Sequence[str] | None = None,
    keep_info: Sequence[str] | None = None,
    remove_info: Sequence[str] | None = None,
    min_maf: float | None = None,
    max_maf: float | None = None,
    min_mac: int | None = None,
    max_mac: int | None = None,
    min_alleles: int | None = None,
    max_alleles: int | None = None,
    min_mean_depth: float | None = None,
    max_mean_depth: float | None = None,
    hwe_pvalue: float | None = None,
    phased_only: bool = False,
    min_quality: float | None = None,
    max_missing: float | None = None,
    keep_individuals: Sequence[str] | None = None,
    keep_individuals_files: Sequence[str] | None = None,
    remove_individuals: Sequence[str] | None = None,
    remove_individuals_files: Sequence[str] | None = None,
    max_individuals: int | None = None,
    remove_filtered_genotypes: bool = False,
    remove_filtered_genotype_flags: Sequence[str] | None = None,
    min_genotype_quality: float | None = None,
    min_genotype_depth: int | None = None,
    max_genotype_depth: int | None = None,
    vcftools_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Filter sites, individuals and genotypes of a VCF/BCF with VCFtools and write the result as VCF or BCF.

    One function for the VCFtools filter options (all that are given are applied together, as in VCFtools):
    ``chromosomes``/``exclude_chromosomes`` (--chr, --not-chr), ``from_bp``/``to_bp`` (--from-bp, --to-bp),
    ``positions_file``/``exclude_positions_file`` (--positions, --exclude-positions),
    ``positions_overlap_file``/``exclude_positions_overlap_file`` (--positions-overlap, --exclude-positions-overlap),
    ``bed_file``/``exclude_bed_file`` (--bed, --exclude-bed), ``thin`` (--thin: minimum spacing in bp),
    ``mask_file``/``invert_mask_file``/``mask_min`` (--mask, --invert-mask, --mask-min), ``snp_ids``/``snp_ids_file``/
    ``exclude_ids_file`` (--snp, --snps, --exclude), ``only_indels``/``remove_indels`` (--keep-only-indels,
    --remove-indels), ``remove_filtered_all``/``keep_filtered``/``remove_filtered`` (--remove-filtered-all,
    --keep-filtered, --remove-filtered), ``keep_info``/``remove_info`` (--keep-INFO, --remove-INFO), ``min_maf``/
    ``max_maf`` (--maf, --max-maf), ``min_mac``/``max_mac`` (--mac, --max-mac), ``min_alleles``/``max_alleles``
    (--min-alleles, --max-alleles), ``min_mean_depth``/``max_mean_depth`` (--min-meanDP, --max-meanDP),
    ``hwe_pvalue`` (--hwe), ``phased_only`` (--phased), ``min_quality`` (--minQ), ``max_missing`` (--max-missing),
    individuals ``keep_individuals``/``keep_individuals_files``/``remove_individuals``/``remove_individuals_files``
    (--indv, --keep, --remove-indv, --remove), ``max_individuals`` (--max-indv: random thinning of individuals),
    and genotype filters ``remove_filtered_genotypes``/``remove_filtered_genotype_flags`` (--remove-filtered-geno-all,
    --remove-filtered-geno), ``min_genotype_quality`` (--minGQ), ``min_genotype_depth``/``max_genotype_depth``
    (--minDP, --maxDP). ``output_format`` is ``"vcf"`` (--recode) or ``"bcf"`` (--recode-bcf). Like VCFtools, INFO is
    dropped from the output unless ``recode_info`` is True (--recode-INFO-all) or a list of INFO keys (--recode-INFO).
    The output is ``<output_prefix>.recode.vcf`` or ``<output_prefix>.recode.bcf``; ``sites_kept`` and
    ``individuals_kept`` come from the VCFtools log.
    """
    fmt = _cli.choice(output_format, ("vcf", "bcf"), "output_format")
    filters = {
        "chromosomes": chromosomes,
        "exclude_chromosomes": exclude_chromosomes,
        "from_bp": from_bp,
        "to_bp": to_bp,
        "positions_file": positions_file,
        "exclude_positions_file": exclude_positions_file,
        "positions_overlap_file": positions_overlap_file,
        "exclude_positions_overlap_file": exclude_positions_overlap_file,
        "bed_file": bed_file,
        "exclude_bed_file": exclude_bed_file,
        "thin": thin,
        "mask_file": mask_file,
        "invert_mask_file": invert_mask_file,
        "mask_min": mask_min,
        "snp_ids": snp_ids,
        "snp_ids_file": snp_ids_file,
        "exclude_ids_file": exclude_ids_file,
        "only_indels": only_indels,
        "remove_indels": remove_indels,
        "remove_filtered_all": remove_filtered_all,
        "keep_filtered": keep_filtered,
        "remove_filtered": remove_filtered,
        "keep_info": keep_info,
        "remove_info": remove_info,
        "min_maf": min_maf,
        "max_maf": max_maf,
        "min_mac": min_mac,
        "max_mac": max_mac,
        "min_alleles": min_alleles,
        "max_alleles": max_alleles,
        "min_mean_depth": min_mean_depth,
        "max_mean_depth": max_mean_depth,
        "hwe_pvalue": hwe_pvalue,
        "phased_only": phased_only,
        "min_quality": min_quality,
        "max_missing": max_missing,
        "keep_individuals": keep_individuals,
        "keep_individuals_files": keep_individuals_files,
        "remove_individuals": remove_individuals,
        "remove_individuals_files": remove_individuals_files,
        "max_individuals": max_individuals,
        "remove_filtered_genotypes": remove_filtered_genotypes,
        "remove_filtered_genotype_flags": remove_filtered_genotype_flags,
        "min_genotype_quality": min_genotype_quality,
        "min_genotype_depth": min_genotype_depth,
        "max_genotype_depth": max_genotype_depth,
    }
    info: list[str] = []
    if recode_info is True:
        info = ["--recode-INFO-all"]
    elif recode_info is not False:
        for key in _sequence_of_text(recode_info, "recode_info", pattern=_NAME_TOKEN):
            info += ["--recode-INFO", key]
    suffix = ".recode.vcf" if fmt == "vcf" else ".recode.bcf"
    return _report(
        vcf_path,
        output_prefix,
        ["--recode" if fmt == "vcf" else "--recode-bcf", *info],
        [suffix],
        filters=filters,
        vcftools_path=vcftools_path,
        timeout_seconds=timeout_seconds,
        table=False,
    )


def vcftools_list_sites(
    vcf_path: str,
    output_prefix: str,
    *,
    which: str = "kept",
    filters: Mapping[str, object] | None = None,
    max_rows: int = _MAX_ROWS,
    vcftools_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """List the sites that pass (``which="kept"``, --kept-sites) or fail (``which="removed"``, --removed-sites) a filter set.

    ``filters`` uses the keyword names of ``vcftools_filter_vcf`` (for example ``{"min_maf": 0.05}``). The list is
    written to ``<output_prefix>.kept.sites`` or ``<output_prefix>.removed.sites`` and returned as a table preview
    (columns CHROM and POS).
    """
    mode = _cli.choice(which, ("kept", "removed"), "which")
    flag = "--kept-sites" if mode == "kept" else "--removed-sites"
    return _report(
        vcf_path,
        output_prefix,
        [flag],
        [f".{mode}.sites"],
        filters=filters,
        vcftools_path=vcftools_path,
        timeout_seconds=timeout_seconds,
        max_rows=max_rows,
    )


# ---------------------------------------------------------------------------
# Depth, quality and missingness reports
# ---------------------------------------------------------------------------

_DEPTH_LEVELS: Mapping[str, tuple[str, str]] = {
    "individual_mean": ("--depth", ".idepth"),
    "site_sum": ("--site-depth", ".ldepth"),
    "site_mean": ("--site-mean-depth", ".ldepth.mean"),
    "genotype": ("--geno-depth", ".gdepth"),
}


def vcftools_depth(
    vcf_path: str,
    output_prefix: str,
    *,
    level: str = "individual_mean",
    filters: Mapping[str, object] | None = None,
    max_rows: int = _MAX_ROWS,
    vcftools_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Read-depth reports from the DP genotype field of a VCF/BCF.

    ``level`` is ``"individual_mean"`` (mean depth per individual, --depth, ``.idepth``), ``"site_sum"`` (depth summed
    over individuals per site, --site-depth, ``.ldepth``), ``"site_mean"`` (mean depth per site, --site-mean-depth,
    ``.ldepth.mean``) or ``"genotype"`` (depth of every genotype, --geno-depth, ``.gdepth``). ``filters`` uses the
    keyword names of ``vcftools_filter_vcf``. The result holds a table preview.
    """
    flag, suffix = _DEPTH_LEVELS[_cli.choice(level, tuple(_DEPTH_LEVELS), "level")]
    return _report(
        vcf_path,
        output_prefix,
        [flag],
        [suffix],
        filters=filters,
        vcftools_path=vcftools_path,
        timeout_seconds=timeout_seconds,
        max_rows=max_rows,
    )


def vcftools_site_quality(
    vcf_path: str,
    output_prefix: str,
    *,
    filters: Mapping[str, object] | None = None,
    max_rows: int = _MAX_ROWS,
    vcftools_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Extract the QUAL value of every site (--site-quality) into ``<output_prefix>.lqual`` (columns CHROM, POS, QUAL)."""
    return _report(
        vcf_path,
        output_prefix,
        ["--site-quality"],
        [".lqual"],
        filters=filters,
        vcftools_path=vcftools_path,
        timeout_seconds=timeout_seconds,
        max_rows=max_rows,
    )


def vcftools_missingness(
    vcf_path: str,
    output_prefix: str,
    *,
    by: str = "individual",
    filters: Mapping[str, object] | None = None,
    max_rows: int = _MAX_ROWS,
    vcftools_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Fraction of missing genotypes per individual (``by="individual"``, --missing-indv, ``.imiss``) or per site
    (``by="site"``, --missing-site, ``.lmiss``)."""
    flag, suffix = {"individual": ("--missing-indv", ".imiss"), "site": ("--missing-site", ".lmiss")}[
        _cli.choice(by, ("individual", "site"), "by")
    ]
    return _report(
        vcf_path,
        output_prefix,
        [flag],
        [suffix],
        filters=filters,
        vcftools_path=vcftools_path,
        timeout_seconds=timeout_seconds,
        max_rows=max_rows,
    )


# ---------------------------------------------------------------------------
# Linkage disequilibrium
# ---------------------------------------------------------------------------


def vcftools_linkage_disequilibrium(
    vcf_path: str,
    output_prefix: str,
    *,
    statistic: str = "geno_r2",
    positions_file: str | None = None,
    interchromosomal: bool = False,
    ld_window_bp: int | None = None,
    ld_window_snps: int | None = None,
    ld_window_bp_min: int | None = None,
    ld_window_snps_min: int | None = None,
    min_r2: float | None = None,
    filters: Mapping[str, object] | None = None,
    max_rows: int = _MAX_ROWS,
    vcftools_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Pairwise linkage disequilibrium of a VCF/BCF with VCFtools.

    ``statistic`` is ``"geno_r2"`` (r2 from the correlation of 0/1/2 genotypes, --geno-r2, ``.geno.ld``),
    ``"geno_chisq"`` (chi-squared test of genotype independence between sites, --geno-chisq, ``.geno.chisq``) or
    ``"hap_r2"`` (r2 from phased haplotypes, --hap-r2, ``.hap.ld``). With ``positions_file`` (tab-separated CHROM, POS)
    only the listed sites are compared with all other sites (--geno-r2-positions ``.list.geno.ld``,
    --hap-r2-positions ``.list.hap.ld``). With ``interchromosomal=True`` sites on different chromosomes are compared
    (--interchrom-geno-r2 ``.interchrom.geno.ld``, --interchrom-hap-r2 ``.interchrom.hap.ld``). ``ld_window_bp``,
    ``ld_window_snps``, ``ld_window_bp_min``, ``ld_window_snps_min`` and ``min_r2`` are --ld-window-bp, --ld-window,
    --ld-window-bp-min, --ld-window-min and --min-r2. The chi-squared test has no list or inter-chromosomal form.
    """
    kind = _cli.choice(statistic, ("geno_r2", "geno_chisq", "hap_r2"), "statistic")
    interchrom = _switch(interchromosomal, "interchromosomal")
    if positions_file is not None and interchrom:
        raise ConfigurationError("positions_file and interchromosomal cannot be used together.")
    if kind == "geno_chisq" and (positions_file is not None or interchrom):
        raise ConfigurationError("geno_chisq has no positions_file or interchromosomal form.")
    haplotype = kind == "hap_r2"
    if positions_file is not None:
        listed = str(_cli.existing_file(positions_file, "positions_file"))
        analysis = ["--hap-r2-positions" if haplotype else "--geno-r2-positions", listed]
        suffix = ".list.hap.ld" if haplotype else ".list.geno.ld"
    elif interchrom:
        analysis = ["--interchrom-hap-r2" if haplotype else "--interchrom-geno-r2"]
        suffix = ".interchrom.hap.ld" if haplotype else ".interchrom.geno.ld"
    else:
        analysis = [{"geno_r2": "--geno-r2", "geno_chisq": "--geno-chisq", "hap_r2": "--hap-r2"}[kind]]
        suffix = {"geno_r2": ".geno.ld", "geno_chisq": ".geno.chisq", "hap_r2": ".hap.ld"}[kind]
    window: list[str] = []
    for flag, value, name in (
        ("--ld-window-bp", ld_window_bp, "ld_window_bp"),
        ("--ld-window", ld_window_snps, "ld_window_snps"),
        ("--ld-window-bp-min", ld_window_bp_min, "ld_window_bp_min"),
        ("--ld-window-min", ld_window_snps_min, "ld_window_snps_min"),
    ):
        window += _cli.option(flag, None if value is None else _cli.integer(value, name, minimum=0))
    if min_r2 is not None:
        window += ["--min-r2", format(_cli.number(min_r2, "min_r2", minimum=0.0, maximum=1.0), ".12g")]
    return _report(
        vcf_path,
        output_prefix,
        analysis,
        [suffix],
        filters=filters,
        vcftools_path=vcftools_path,
        timeout_seconds=timeout_seconds,
        max_rows=max_rows,
        extra=window,
    )


# ---------------------------------------------------------------------------
# Ts/Tv, diversity, heterozygosity, burden, relatedness, density
# ---------------------------------------------------------------------------

_TSTV_MODES: Mapping[str, tuple[str, str]] = {
    "distance_bins": ("--TsTv", ".TsTv"),
    "alt_allele_count": ("--TsTv-by-count", ".TsTv.count"),
    "qual_threshold": ("--TsTv-by-qual", ".TsTv.qual"),
    "filter_category": ("--FILTER-summary", ".FILTER.summary"),
}


def vcftools_transition_transversion(
    vcf_path: str,
    output_prefix: str,
    *,
    mode: str = "distance_bins",
    bin_size: int | None = None,
    filters: Mapping[str, object] | None = None,
    max_rows: int = _MAX_ROWS,
    vcftools_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Transition/transversion (Ts/Tv) counts of bi-allelic SNPs in four groupings.

    ``mode`` is ``"distance_bins"`` (Ts/Tv in bins of ``bin_size`` bp, --TsTv, ``.TsTv``), ``"alt_allele_count"``
    (by alternative allele count, --TsTv-by-count, ``.TsTv.count``), ``"qual_threshold"`` (by QUAL threshold,
    --TsTv-by-qual, ``.TsTv.qual``) or ``"filter_category"`` (SNP count and Ts/Tv per FILTER value, --FILTER-summary,
    ``.FILTER.summary``). ``bin_size`` is required for ``distance_bins`` and ignored otherwise.
    """
    kind = _cli.choice(mode, tuple(_TSTV_MODES), "mode")
    flag, suffix = _TSTV_MODES[kind]
    analysis = [flag]
    if kind == "distance_bins":
        if bin_size is None:
            raise ConfigurationError("bin_size is required for mode='distance_bins'.")
        analysis.append(str(_cli.integer(bin_size, "bin_size", minimum=1)))
    return _report(
        vcf_path,
        output_prefix,
        analysis,
        [suffix],
        filters=filters,
        vcftools_path=vcftools_path,
        timeout_seconds=timeout_seconds,
        max_rows=max_rows,
    )


def vcftools_site_nucleotide_diversity(
    vcf_path: str,
    output_prefix: str,
    *,
    filters: Mapping[str, object] | None = None,
    max_rows: int = _MAX_ROWS,
    vcftools_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Nucleotide diversity (pi) of every site (--site-pi) into ``<output_prefix>.sites.pi`` (columns CHROM, POS, PI)."""
    return _report(
        vcf_path,
        output_prefix,
        ["--site-pi"],
        [".sites.pi"],
        filters=filters,
        vcftools_path=vcftools_path,
        timeout_seconds=timeout_seconds,
        max_rows=max_rows,
    )


def vcftools_heterozygosity(
    vcf_path: str,
    output_prefix: str,
    *,
    filters: Mapping[str, object] | None = None,
    max_rows: int = _MAX_ROWS,
    vcftools_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Per-individual homozygosity and inbreeding coefficient F (--het) into ``<output_prefix>.het``.

    Columns: INDV, O(HOM), E(HOM), N_SITES, F (method-of-moments F, as in PLINK --het).
    """
    return _report(
        vcf_path,
        output_prefix,
        ["--het"],
        [".het"],
        filters=filters,
        vcftools_path=vcftools_path,
        timeout_seconds=timeout_seconds,
        max_rows=max_rows,
    )


_BURDEN_MODES: Mapping[str, tuple[str, str]] = {
    "genotype_classes": ("--indv-burden", ".iburden"),
    "frequency_classes": ("--indv-freq-burden", ".ifreqburden"),
    "frequency_classes_hom_once": ("--indv-freq-burden2", ".ifreqburden"),
}


def vcftools_individual_burden(
    vcf_path: str,
    output_prefix: str,
    *,
    burden: str = "genotype_classes",
    filters: Mapping[str, object] | None = None,
    max_rows: int = _MAX_ROWS,
    vcftools_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Variant burden of each individual.

    ``burden="genotype_classes"`` counts homozygous-reference, heterozygous, homozygous-alternative and missing
    genotypes (--indv-burden, ``.iburden``). ``"frequency_classes"`` counts non-reference alleles by the allele count
    of their site, counting both alleles of a homozygous genotype (--indv-freq-burden, ``.ifreqburden``).
    ``"frequency_classes_hom_once"`` is the same but counts a homozygous genotype once (--indv-freq-burden2).
    """
    flag, suffix = _BURDEN_MODES[_cli.choice(burden, tuple(_BURDEN_MODES), "burden")]
    return _report(
        vcf_path,
        output_prefix,
        [flag],
        [suffix],
        filters=filters,
        vcftools_path=vcftools_path,
        timeout_seconds=timeout_seconds,
        max_rows=max_rows,
    )


def vcftools_long_runs_of_homozygosity(
    vcf_path: str,
    output_prefix: str,
    chromosome: str,
    *,
    filters: Mapping[str, object] | None = None,
    max_rows: int = _MAX_ROWS,
    vcftools_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Long runs of homozygosity of one chromosome (--LROH, ``.LROH``).

    VCFtools requires the chromosome (--chr); it is the third argument. ``filters`` must not set ``chromosomes``.
    """
    contig = _cli.text(chromosome, "chromosome")
    if filters is not None and filters.get("chromosomes"):
        raise ConfigurationError("Pass the chromosome as 'chromosome', not in filters['chromosomes'].")
    return _report(
        vcf_path,
        output_prefix,
        ["--LROH"],
        [".LROH"],
        filters=filters,
        vcftools_path=vcftools_path,
        timeout_seconds=timeout_seconds,
        max_rows=max_rows,
        extra=["--chr", contig],
    )


def vcftools_relatedness(
    vcf_path: str,
    output_prefix: str,
    *,
    method: str = "yang",
    filters: Mapping[str, object] | None = None,
    max_rows: int = _MAX_ROWS,
    vcftools_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Pairwise relatedness between individuals.

    ``method="yang"`` is the statistic of Yang et al. 2010 (--relatedness, ``.relatedness``, column RELATEDNESS_AJK).
    ``method="manichaikul"`` is the KING-robust-style estimator of Manichaikul et al. 2010 (--relatedness2,
    ``.relatedness2``, column RELATEDNESS_PHI).
    """
    flag, suffix = {"yang": ("--relatedness", ".relatedness"), "manichaikul": ("--relatedness2", ".relatedness2")}[
        _cli.choice(method, ("yang", "manichaikul"), "method")
    ]
    return _report(
        vcf_path,
        output_prefix,
        [flag],
        [suffix],
        filters=filters,
        vcftools_path=vcftools_path,
        timeout_seconds=timeout_seconds,
        max_rows=max_rows,
    )


def vcftools_snp_density(
    vcf_path: str,
    output_prefix: str,
    bin_size: int,
    *,
    filters: Mapping[str, object] | None = None,
    max_rows: int = _MAX_ROWS,
    vcftools_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Number of SNPs in non-overlapping bins of ``bin_size`` bp (--SNPdensity) into ``<output_prefix>.snpden``."""
    size = str(_cli.integer(bin_size, "bin_size", minimum=1))
    return _report(
        vcf_path,
        output_prefix,
        ["--SNPdensity", size],
        [".snpden"],
        filters=filters,
        vcftools_path=vcftools_path,
        timeout_seconds=timeout_seconds,
        max_rows=max_rows,
    )


def vcftools_singletons(
    vcf_path: str,
    output_prefix: str,
    *,
    filters: Mapping[str, object] | None = None,
    max_rows: int = _MAX_ROWS,
    vcftools_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Locations of singleton and doubleton alleles (--singletons) in ``<output_prefix>.singletons``.

    Columns: CHROM, POS, SINGLETON/DOUBLETON, ALLELE, INDV.
    """
    return _report(
        vcf_path,
        output_prefix,
        ["--singletons"],
        [".singletons"],
        filters=filters,
        vcftools_path=vcftools_path,
        timeout_seconds=timeout_seconds,
        max_rows=max_rows,
    )


def vcftools_indel_length_histogram(
    vcf_path: str,
    output_prefix: str,
    *,
    filters: Mapping[str, object] | None = None,
    max_rows: int = _MAX_ROWS,
    vcftools_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Histogram of indel lengths (--hist-indel-len) in ``<output_prefix>.indel.hist`` (columns LENGTH, COUNT, PRCT).

    Negative lengths are deletions, positive lengths insertions, 0 are SNPs.
    """
    return _report(
        vcf_path,
        output_prefix,
        ["--hist-indel-len"],
        [".indel.hist"],
        filters=filters,
        vcftools_path=vcftools_path,
        timeout_seconds=timeout_seconds,
        max_rows=max_rows,
    )


def vcftools_haplotype_count(
    vcf_path: str,
    output_prefix: str,
    regions_bed: str,
    *,
    filters: Mapping[str, object] | None = None,
    max_rows: int = _MAX_ROWS,
    vcftools_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Number of unique haplotypes in each BED region of a phased VCF (--hapcount) into ``<output_prefix>.hapcount``.

    Columns: CHROM, BIN_START, BIN_END, N_SNP, N_UNIQ_HAPS, N_GROUPS and the haplotype multiplicities. The regions
    must not overlap. VCFtools 0.1.16 treats the first line of the BED file as a header and skips it (the result
    lists a ``warnings`` entry when the first line does not look like a header), and it reads freed memory when the
    last record of the VCF lies inside the last region of the last chromosome, which can zero the last row; keep the
    last record outside the last region (the result warns about the possibility).
    """
    bed = _cli.existing_file(regions_bed, "regions_bed")
    result = _report(
        vcf_path,
        output_prefix,
        ["--hapcount", str(bed)],
        [".hapcount"],
        filters=filters,
        vcftools_path=vcftools_path,
        timeout_seconds=timeout_seconds,
        max_rows=max_rows,
    )
    warnings = ["VCFtools 0.1.16 --hapcount can return zeros for the last region when the last VCF record lies in it."]
    with bed.open(encoding="utf-8", errors="replace") as handle:
        first = handle.readline()
    if not first.startswith(("#", "track", "browser")):
        warnings.append("VCFtools skips the first line of the BED file; it is not a header line, so that region was dropped.")
    result["warnings"] = warnings
    return result


# ---------------------------------------------------------------------------
# Field extraction, PCA, export, comparison of two VCFs
# ---------------------------------------------------------------------------


def vcftools_extract_fields(
    vcf_path: str,
    output_prefix: str,
    fields: Sequence[str],
    *,
    source: str = "info",
    filters: Mapping[str, object] | None = None,
    max_rows: int = _MAX_ROWS,
    vcftools_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Extract INFO or FORMAT fields into a table.

    ``source="info"`` writes the listed INFO keys for every site (--get-INFO, ``<output_prefix>.INFO``).
    ``source="format"`` writes one FORMAT key for every genotype of every site (--extract-FORMAT-info,
    ``<output_prefix>.<field>.FORMAT``); VCFtools accepts one FORMAT field per run, so ``fields`` has one entry.
    """
    kind = _cli.choice(source, ("info", "format"), "source")
    keys = _sequence_of_text(fields, "fields", pattern=_NAME_TOKEN)
    if kind == "format":
        if len(keys) != 1:
            raise ConfigurationError("source='format' extracts exactly one field per run.")
        analysis = ["--extract-FORMAT-info", keys[0]]
        suffix = f".{keys[0]}.FORMAT"
    else:
        analysis = [part for key in keys for part in ("--get-INFO", key)]
        suffix = ".INFO"
    return _report(
        vcf_path,
        output_prefix,
        analysis,
        [suffix],
        filters=filters,
        vcftools_path=vcftools_path,
        timeout_seconds=timeout_seconds,
        max_rows=max_rows,
    )


def vcftools_pca(
    vcf_path: str,
    output_prefix: str,
    *,
    output: str = "components",
    normalise: bool = True,
    snp_loadings: int | None = None,
    filters: Mapping[str, object] | None = None,
    max_rows: int = _MAX_ROWS,
    vcftools_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Principal component analysis of the genotypes (Patterson, Price and Reich 2006) with VCFtools.

    ``output="components"`` writes the principal components of the individuals (--pca, or --pca-no-norm when
    ``normalise=False``; ``<output_prefix>.pca``). ``output="snp_loadings"`` writes the loadings of the first
    ``snp_loadings`` components for every SNP (--pca-snp-loadings, ``<output_prefix>.pca.loadings``) and is only
    available with normalisation. Needs a ``vcftools`` built with ``--enable-pca``, bi-allelic diploid sites and more
    sites than individuals.
    """
    kind = _cli.choice(output, ("components", "snp_loadings"), "output")
    scaled = _switch(normalise, "normalise")
    if kind == "components":
        return _report(
            vcf_path,
            output_prefix,
            ["--pca" if scaled else "--pca-no-norm"],
            [".pca"],
            filters=filters,
            vcftools_path=vcftools_path,
            timeout_seconds=timeout_seconds,
            max_rows=max_rows,
        )
    if snp_loadings is None:
        raise ConfigurationError("snp_loadings (the number of components) is required for output='snp_loadings'.")
    if not scaled:
        raise ConfigurationError("SNP loadings are only written together with normalisation (normalise=True).")
    return _report(
        vcf_path,
        output_prefix,
        ["--pca-snp-loadings", str(_cli.integer(snp_loadings, "snp_loadings", minimum=1))],
        [".pca.loadings"],
        filters=filters,
        vcftools_path=vcftools_path,
        timeout_seconds=timeout_seconds,
        max_rows=max_rows,
    )


# format -> (flag, output suffixes, needs a single chromosome)
_EXPORT_FORMATS: Mapping[str, tuple[str, tuple[str, ...], bool]] = {
    "012": ("--012", (".012", ".012.indv", ".012.pos"), False),
    "impute": ("--IMPUTE", (".impute.hap", ".impute.hap.indv", ".impute.legend"), False),
    "ldhat": ("--ldhat", (".ldhat.sites", ".ldhat.locs"), True),
    "ldhelmet": ("--ldhelmet", (".ldhelmet.snps", ".ldhelmet.pos"), True),
    "ldhat_geno": ("--ldhat-geno", (".ldhat.sites", ".ldhat.locs"), True),
    "beagle_gl": ("--BEAGLE-GL", (".BEAGLE.GL",), True),
    "beagle_pl": ("--BEAGLE-PL", (".BEAGLE.PL",), True),
    "plink": ("--plink", (".ped", ".map"), False),
    "plink_tped": ("--plink-tped", (".tped", ".tfam"), False),
}


def vcftools_export_genotypes(
    vcf_path: str,
    output_prefix: str,
    *,
    export_format: str = "012",
    chromosome: str | None = None,
    filters: Mapping[str, object] | None = None,
    vcftools_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Export genotypes of a VCF/BCF to another file format with VCFtools.

    ``export_format``: ``"012"`` (matrix of 0/1/2 alternative-allele counts, --012; ``.012``, ``.012.indv``,
    ``.012.pos``), ``"impute"`` (phased haplotypes in IMPUTE format, --IMPUTE; ``.impute.hap``, ``.impute.hap.indv``,
    ``.impute.legend``), ``"ldhat"`` (phased haplotypes in LDhat format, --ldhat; ``.ldhat.sites``, ``.ldhat.locs``),
    ``"ldhelmet"`` (LDhelmet format, --ldhelmet; ``.ldhelmet.snps``, ``.ldhelmet.pos``), ``"ldhat_geno"`` (unphased
    genotypes in LDhat format, --ldhat-geno), ``"beagle_gl"`` / ``"beagle_pl"`` (BEAGLE genotype likelihoods from the
    GL or PL field, --BEAGLE-GL, --BEAGLE-PL; ``.BEAGLE.GL``, ``.BEAGLE.PL``), ``"plink"`` (PED/MAP, --plink) and
    ``"plink_tped"`` (TPED/TFAM, --plink-tped). The LDhat, LDhelmet and BEAGLE formats require one ``chromosome``
    (--chr); the result lists every written file in ``outputs``.
    """
    kind = _cli.choice(export_format, tuple(_EXPORT_FORMATS), "export_format")
    flag, suffixes, needs_chromosome = _EXPORT_FORMATS[kind]
    extra: list[str] = []
    if chromosome is not None:
        if filters is not None and filters.get("chromosomes"):
            raise ConfigurationError("Pass the chromosome as 'chromosome', not in filters['chromosomes'].")
        extra = ["--chr", _cli.text(chromosome, "chromosome")]
    elif needs_chromosome and not (filters and filters.get("chromosomes")):
        raise ConfigurationError(f"export_format={kind!r} requires a chromosome (VCFtools --chr).")
    return _report(
        vcf_path,
        output_prefix,
        [flag],
        list(suffixes),
        filters=filters,
        vcftools_path=vcftools_path,
        timeout_seconds=timeout_seconds,
        table=False,
        extra=extra,
    )


_DIFF_MODES: Mapping[str, tuple[str, tuple[str, ...]]] = {
    "individuals": ("--diff-indv", (".diff.indv_in_files",)),
    "site_discordance": ("--diff-site-discordance", (".diff.sites",)),
    "switch_error": ("--diff-switch-error", (".diff.switch", ".diff.indv.switch")),
}


def vcftools_compare_vcfs(
    vcf_path: str,
    second_vcf_path: str,
    output_prefix: str,
    *,
    mode: str = "individuals",
    individual_map_file: str | None = None,
    filters: Mapping[str, object] | None = None,
    max_rows: int = _MAX_ROWS,
    vcftools_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Compare two VCF/BCF files with the VCFtools ``--diff`` family.

    ``mode="individuals"`` lists which individuals are in the first file, the second file or both (--diff-indv,
    ``.diff.indv_in_files``; FILES is 1, 2 or B). ``"site_discordance"`` reports genotype discordance at every shared
    site (--diff-site-discordance, ``.diff.sites``). ``"switch_error"`` reports phase switches between phased
    genotypes (--diff-switch-error, ``.diff.switch`` and ``.diff.indv.switch``). ``individual_map_file`` renames
    individuals of the second file (--diff-indv-map). ``filters`` apply to the first file.
    """
    kind = _cli.choice(mode, tuple(_DIFF_MODES), "mode")
    flag, suffixes = _DIFF_MODES[kind]
    extra = _input_flag(second_vcf_path, second=True)
    if individual_map_file is not None:
        extra += ["--diff-indv-map", str(_cli.existing_file(individual_map_file, "individual_map_file"))]
    return _report(
        vcf_path,
        output_prefix,
        [flag],
        list(suffixes),
        filters=filters,
        vcftools_path=vcftools_path,
        timeout_seconds=timeout_seconds,
        max_rows=max_rows,
        extra=extra,
    )


# ---------------------------------------------------------------------------
# Running the VCFtools Perl scripts
# ---------------------------------------------------------------------------


def _perl_environment(executable: str) -> dict[str, str]:
    """``PERL5LIB`` with the folder holding ``Vcf.pm`` next to the installed script, in front of an existing value."""
    root = Path(executable).parent.parent
    candidates = [root / "share" / "perl5", root / "lib" / "perl5", root / "lib" / "perl", root / "perl"]
    candidates += sorted((root / "share" / "perl").glob("*")) if (root / "share" / "perl").is_dir() else []
    found = [str(path) for path in candidates if (path / "Vcf.pm").is_file()]
    if not found:
        return {}
    existing = os.environ.get("PERL5LIB")
    return {"PERL5LIB": os.pathsep.join([*found, existing] if existing else found)}


@lru_cache(maxsize=16)
def _vcf_pm_version(perl_lib: str) -> str | None:
    perl = shutil.which("perl")
    if perl is None:
        return None
    try:
        result = execute_bounded_command(
            perl,
            ["-MVcf", "-e", "print $Vcf::VERSION"],
            backend_id="perl",
            cwd=Path(tempfile.gettempdir()),
            timeout_seconds=60.0,
            max_output_bytes=10_000,
            environment_overrides={"PERL5LIB": perl_lib} if perl_lib else None,
            verbose=False,
        )
    except DNAKitError:
        return None
    text = (result.output or "").strip()
    return text if result.return_code == 0 and text else None


def _run_script(
    script: str,
    arguments: Sequence[str],
    *,
    output_path: Path | None,
    stdin_path: Path | None = None,
    extra_outputs: Sequence[Path] = (),
    script_path: str | None,
    timeout_seconds: float,
) -> dict[str, Any]:
    """Run one VCFtools Perl script without a shell; standard output goes to the new ``output_path`` file.

    ``completed`` only when the exit code is 0 and every output exists. A failed run leaves no partial output file.
    """
    executable = _cli.program(script, explicit=script_path, hint=_SCRIPT_HINT)
    limit = _cli.number(timeout_seconds, "timeout_seconds", minimum=1.0)
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    environment.update(_perl_environment(executable))
    command = [executable, *arguments]
    cwd = (output_path or Path(tempfile.gettempdir())).parent
    with contextlib.ExitStack() as stack:
        stdin = stack.enter_context(stdin_path.open("rb")) if stdin_path is not None else subprocess.DEVNULL
        if output_path is not None:
            stdout = stack.enter_context(output_path.open("xb"))
        else:
            stdout = stack.enter_context(tempfile.TemporaryFile())
        stderr = stack.enter_context(tempfile.TemporaryFile())
        process = subprocess.Popen(
            command,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            cwd=cwd,
            env=environment,
            shell=False,
            start_new_session=os.name == "posix",
        )
        deadline = time.monotonic() + limit
        try:
            while process.poll() is None:
                if time.monotonic() > deadline:
                    raise BackendTimeoutError(f"{script} exceeded its time limit.")
                if os.fstat(stdout.fileno()).st_size > _MAX_ARTIFACT_BYTES:
                    raise BackendExecutionError(f"{script} exceeded its output size limit.")
                time.sleep(0.02)
        except DNAKitError:
            if output_path is not None:
                output_path.unlink(missing_ok=True)
            raise
        finally:
            if process.poll() is None:
                if os.name == "posix":
                    with contextlib.suppress(ProcessLookupError):
                        os.killpg(process.pid, signal.SIGKILL)
                else:
                    process.kill()
            process.wait()
        stderr.seek(0)
        log = stderr.read(2_000_000).decode("utf-8", errors="replace")
        if output_path is None:
            stdout.seek(0)
            log = stdout.read(2_000_000).decode("utf-8", errors="replace") + log
    wanted = ([output_path] if output_path is not None else []) + list(extra_outputs)
    complete = process.returncode == 0 and all(path.is_file() for path in wanted)
    if not complete and output_path is not None:
        output_path.unlink(missing_ok=True)
    perl_lib = environment.get("PERL5LIB", "")
    return {
        "status": "completed" if complete else "failed",
        "backend": script,
        "version": _vcf_pm_version(perl_lib),
        "command": [script, *arguments],
        "return_code": process.returncode,
        "outputs": [str(path) for path in wanted if path.is_file()] if complete else [],
        "log_tail": log[-_LOG_TAIL_CHARS:],
    }


def _script_result(
    script: str,
    arguments: Sequence[str],
    output: Path | None,
    *,
    stdin: Path | None = None,
    script_path: str | None,
    timeout_seconds: float,
    extra_outputs: Sequence[Path] = (),
    records: bool = False,
) -> dict[str, Any]:
    result = _run_script(
        script,
        arguments,
        output_path=output,
        stdin_path=stdin,
        extra_outputs=extra_outputs,
        script_path=script_path,
        timeout_seconds=timeout_seconds,
    )
    if stdin is not None:
        result["stdin"] = str(stdin)
    if output is not None:
        result["output_path"] = str(output)
        if records and result["status"] == "completed":
            result["records"] = _count_records(output)
    return result


def _sample_names(value: object, name: str) -> list[str]:
    names = _sequence_of_text(value, name)
    for item in names:
        if _SAMPLE.fullmatch(item) is None:
            raise ConfigurationError(f"{name} contains an unsupported sample name: {item!r}")
    return names


# ---------------------------------------------------------------------------
# vcf-annotate, vcf-subset, vcf-convert, vcf-to-tab, vcf-fix-ploidy, vcf-shuffle-cols, vcf-contrast
# ---------------------------------------------------------------------------


def vcf_annotate(
    vcf_path: str,
    output_path: str,
    *,
    filters: str | None = None,
    hard_filter: bool = False,
    normalize_alleles: bool = False,
    remove_tags: Sequence[str] | None = None,
    fill_ac_an: bool = False,
    fill_hwe: bool = False,
    fill_type: bool = False,
    script_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Annotate a VCF with the VCFtools ``vcf-annotate`` script.

    ``filters`` applies the built-in site filters to the FILTER column with the script's own syntax
    (``vcf-annotate -f``), for example ``"+/Qual=10/MinDP=5"`` (``+`` switches on every filter with its default
    value, ``-X`` excludes filter X; user-defined filter files are not supported). ``hard_filter`` removes every
    site whose FILTER is not PASS or "." (-H). ``normalize_alleles`` makes REF and ALT compact (-n, for example
    TA,TAA to T,TA). ``remove_tags`` removes columns or tags (-r; ``ID``, ``QUAL``, ``FILTER``, ``INFO/DP``,
    ``FORMAT/DP``, ``FILTER/MinDP`` or a bare tag). ``fill_ac_an`` recalculates AC and AN (--fill-AC-AN, the same
    as the ``fill-an-ac`` script), ``fill_hwe`` also adds the HWE p-value (--fill-HWE) and ``fill_type`` writes
    INFO/TYPE as snp, del, ins, mnp or complex (--fill-type, the variant-type classification of ``Vcf::event_type``).
    At least one operation must be selected. The VCF may be compressed. The annotated VCF is written to the new
    ``output_path`` and ``records`` counts its sites.
    """
    source = _cli.existing_file(vcf_path, "vcf_path")
    out = _output_path(output_path)
    arguments: list[str] = []
    if filters is not None:
        spec = _cli.text(filters, "filters")
        if _FILTER_SPEC.fullmatch(spec) is None or (out.parent / spec).exists():
            raise ConfigurationError(f"filters must be a built-in filter specification such as '+/Qual=10': {spec!r}")
        arguments += ["-f", spec]
    if _switch(hard_filter, "hard_filter"):
        arguments.append("-H")
    if _switch(normalize_alleles, "normalize_alleles"):
        arguments.append("-n")
    if remove_tags is not None:
        tags = _sequence_of_text(remove_tags, "remove_tags", pattern=_REMOVE_TAG)
        arguments += ["-r", ",".join(tags)]
    fill_counts = _switch(fill_ac_an, "fill_ac_an")
    if _switch(fill_hwe, "fill_hwe"):
        arguments.append("--fill-HWE")
    elif fill_counts:
        arguments.append("--fill-AC-AN")
    if _switch(fill_type, "fill_type"):
        arguments.append("--fill-type")
    if not arguments:
        raise ConfigurationError(
            "Select at least one operation: filters, hard_filter, normalize_alleles, remove_tags, fill_ac_an, "
            "fill_hwe or fill_type."
        )
    return _script_result(
        "vcf-annotate",
        [*arguments, str(source)],
        out,
        script_path=script_path,
        timeout_seconds=timeout_seconds,
        records=True,
    )


_SUBSET_TYPES = ("ref", "SNPs", "indels", "MNPs", "other")


def vcf_subset(
    vcf_path: str,
    output_path: str,
    *,
    samples: Sequence[str] | None = None,
    samples_file: str | None = None,
    exclude_ref_only: bool = False,
    private_only: bool = False,
    trim_alt_alleles: bool = False,
    variant_types: Sequence[str] | None = None,
    replace_with_ref: bool = False,
    keep_uncalled: bool = False,
    force: bool = False,
    script_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Select samples and variant classes of a VCF with the VCFtools ``vcf-subset`` script.

    ``samples`` (or ``samples_file`` with one sample per line) keeps those columns (-c). ``exclude_ref_only`` drops
    sites without an ALT allele in the selection (-e). ``private_only`` keeps only sites where just the selected
    samples carry an alternate allele (-p). ``trim_alt_alleles`` removes ALT alleles that the selection does not use
    (-a). ``variant_types`` keeps only the listed classes among ``ref``, ``SNPs``, ``indels``, ``MNPs`` and ``other``
    (-t); with ``replace_with_ref`` the excluded types become the reference allele instead of "." (-r).
    ``keep_uncalled`` keeps sites without calls (-u) and ``force`` proceeds when a sample is absent (-f).
    """
    source = _cli.existing_file(vcf_path, "vcf_path")
    out = _output_path(output_path)
    arguments: list[str] = []
    if samples is not None and samples_file is not None:
        raise ConfigurationError("Give samples or samples_file, not both.")
    if samples is not None:
        arguments += ["-c", ",".join(_sample_names(samples, "samples"))]
    elif samples_file is not None:
        arguments += ["-c", str(_cli.existing_file(samples_file, "samples_file"))]
    for flag, value, name in (
        ("-e", exclude_ref_only, "exclude_ref_only"),
        ("-p", private_only, "private_only"),
        ("-a", trim_alt_alleles, "trim_alt_alleles"),
        ("-r", replace_with_ref, "replace_with_ref"),
        ("-u", keep_uncalled, "keep_uncalled"),
        ("-f", force, "force"),
    ):
        if _switch(value, name):
            arguments.append(flag)
    if variant_types is not None:
        kinds = _sequence_of_text(variant_types, "variant_types")
        for kind in kinds:
            _cli.choice(kind, _SUBSET_TYPES, "variant_types")
        arguments += ["-t", ",".join(kinds)]
    return _script_result(
        "vcf-subset",
        [*arguments, str(source)],
        out,
        script_path=script_path,
        timeout_seconds=timeout_seconds,
        records=True,
    )


def vcf_convert(
    vcf_path: str,
    output_path: str,
    *,
    target_version: str = "4.2",
    reference_fasta: str | None = None,
    script_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Convert a VCF between format versions 4.0, 4.1 and 4.2 with ``vcf-convert``.

    ``reference_fasta`` (a samtools-faidx-indexed FASTA) is needed only when the file has indels (-r). The converted
    VCF is written to the new ``output_path``.
    """
    source = _cli.existing_file(vcf_path, "vcf_path")
    out = _output_path(output_path)
    arguments = ["-v", _cli.choice(target_version, ("4.0", "4.1", "4.2"), "target_version")]
    if reference_fasta is not None:
        arguments += ["-r", str(_cli.existing_file(reference_fasta, "reference_fasta"))]
    return _script_result(
        "vcf-convert",
        arguments,
        out,
        stdin=source,
        script_path=script_path,
        timeout_seconds=timeout_seconds,
        records=True,
    )


def vcf_to_tab(
    vcf_path: str,
    output_path: str,
    *,
    iupac: bool = False,
    script_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Convert a VCF to a tab-delimited table of genotypes as alleles (``vcf-to-tab``).

    Columns are CHROM, POS, REF (an ID column is not written) and one column of genotype bases per sample;
    ``iupac`` writes heterozygous genotypes as one-letter IUPAC codes (-i). The table is written to the new
    ``output_path`` and ``records`` counts its data lines.
    """
    source = _cli.existing_file(vcf_path, "vcf_path")
    out = _output_path(output_path)
    result = _script_result(
        "vcf-to-tab",
        ["-i"] if _switch(iupac, "iupac") else [],
        out,
        stdin=source,
        script_path=script_path,
        timeout_seconds=timeout_seconds,
    )
    if result["status"] == "completed":
        result["records"] = _count_records(out)
    return result


def vcf_fix_ploidy(
    vcf_path: str,
    output_path: str,
    samples_file: str,
    *,
    assumed_sex: str | None = None,
    fix_likelihoods: bool = False,
    ploidy_file: str | None = None,
    script_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Set the genotype ploidy of X, Y and MT by sample sex with ``vcf-fix-ploidy``.

    ``samples_file`` lists ``sample_name M|F`` per line (-s). ``assumed_sex`` (``"M"`` or ``"F"``) covers samples
    missing from that list (-a). ``fix_likelihoods`` adds or removes heterozygous likelihoods (-l). ``ploidy_file``
    replaces the default human ploidy definition (-p). The fixed VCF is written to the new ``output_path``.
    """
    source = _cli.existing_file(vcf_path, "vcf_path")
    out = _output_path(output_path)
    arguments = ["-s", str(_cli.existing_file(samples_file, "samples_file"))]
    if assumed_sex is not None:
        arguments += ["-a", _cli.choice(assumed_sex, ("M", "F"), "assumed_sex")]
    if _switch(fix_likelihoods, "fix_likelihoods"):
        arguments.append("-l")
    if ploidy_file is not None:
        arguments += ["-p", str(_cli.existing_file(ploidy_file, "ploidy_file"))]
    return _script_result(
        "vcf-fix-ploidy",
        arguments,
        out,
        stdin=source,
        script_path=script_path,
        timeout_seconds=timeout_seconds,
        records=True,
    )


def vcf_shuffle_columns(
    vcf_path: str,
    template_vcf_path: str,
    output_path: str,
    *,
    script_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Reorder the sample columns of a VCF to the order of a template VCF (``vcf-shuffle-cols -t``).

    ``vcf_path`` and ``template_vcf_path`` must contain the same samples. The reordered VCF is written to the new
    ``output_path``.
    """
    source = _cli.existing_file(vcf_path, "vcf_path")
    template = _cli.existing_file(template_vcf_path, "template_vcf_path")
    out = _output_path(output_path)
    return _script_result(
        "vcf-shuffle-cols",
        ["-t", str(template), str(source)],
        out,
        script_path=script_path,
        timeout_seconds=timeout_seconds,
        records=True,
    )


def vcf_contrast(
    vcf_path: str,
    output_path: str,
    variant_samples: Sequence[str],
    background_samples: Sequence[str],
    *,
    min_depth: int | None = None,
    apply_filters: bool = False,
    novel_sites_only: bool = False,
    script_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Find genotypes and alleles that are unique to some samples (``vcf-contrast``).

    ``variant_samples`` are the samples where a unique variant is expected (+list) and ``background_samples`` the
    others (-list). NOVELGT, NOVELAL and NOVELTY are added to INFO. ``min_depth`` is the minimum depth across the
    background samples (-d), ``apply_filters`` skips sites whose FILTER is not PASS or "." (-f) and
    ``novel_sites_only`` prints only records with novel genotypes (-n). The annotated VCF is written to the new
    ``output_path``.
    """
    source = _cli.existing_file(vcf_path, "vcf_path")
    out = _output_path(output_path)
    arguments = [
        "+" + ",".join(_sample_names(variant_samples, "variant_samples")),
        "-" + ",".join(_sample_names(background_samples, "background_samples")),
    ]
    if min_depth is not None:
        arguments += ["-d", str(_cli.integer(min_depth, "min_depth", minimum=0))]
    if _switch(apply_filters, "apply_filters"):
        arguments.append("-f")
    if _switch(novel_sites_only, "novel_sites_only"):
        arguments.append("-n")
    return _script_result(
        "vcf-contrast",
        [*arguments, str(source)],
        out,
        script_path=script_path,
        timeout_seconds=timeout_seconds,
        records=True,
    )


# ---------------------------------------------------------------------------
# vcf-compare, vcf-indel-stats, vcf-phased-join
# ---------------------------------------------------------------------------


def _parse_compare(text: str) -> dict[str, Any]:
    venn: list[dict[str, Any]] = []
    summary: dict[str, Any] = {}
    genotype_classes: list[dict[str, Any]] = []
    for line in text.splitlines():
        fields = line.split("\t")
        if fields[0] == "VN" and len(fields) >= 3:
            venn.append({"sites": int(fields[1]), "files": fields[2:]})
        elif fields[0] == "SN" and len(fields) >= 3:
            summary[fields[1].rstrip(":")] = _cell(fields[2])
        elif fields[0] == "GS" and len(fields) >= 5:
            genotype_classes.append(
                {
                    "genotype_class": fields[1],
                    "mismatches": _cell(fields[2]),
                    "matches": _cell(fields[3]),
                    "discordance": fields[4],
                }
            )
    return {"venn": venn, "summary": summary, "genotype_classes": genotype_classes}


def vcf_compare(
    vcf_paths: Sequence[str],
    output_path: str,
    *,
    compare_genotypes: bool = False,
    apply_filters: bool = False,
    ignore_indels: bool = False,
    regions: Sequence[str] | None = None,
    samples: Sequence[str] | None = None,
    window: int | None = None,
    reference_fasta: str | None = None,
    script_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Compare the sites and genotypes of two or more VCF files with ``vcf-compare``.

    The files must be bgzip-compressed and tabix-indexed (``tabix -p vcf``). ``compare_genotypes`` also compares
    genotypes (-g), ``ignore_indels`` excludes indel sites from that comparison, ``apply_filters`` ignores sites whose
    FILTER is not PASS or "." (-a), ``regions`` restricts to regions such as ``["chr1:1-100000"]`` (-r), ``samples``
    restricts to samples (-s), ``window`` treats records within that many bp as matching (-w) and
    ``reference_fasta`` compares the actual sequence (-R). The full report is written to the new ``output_path``; the
    result also holds ``venn`` (sites per file combination), ``summary`` (the SN lines) and ``genotype_classes``
    (the GS lines).
    """
    files = [str(_cli.existing_file(path, "vcf_paths")) for path in vcf_paths] if not isinstance(vcf_paths, str) else []
    if len(files) < 2:
        raise ConfigurationError("vcf_paths must list at least two VCF files.")
    out = _output_path(output_path)
    arguments: list[str] = []
    for flag, value, name in (
        ("-g", compare_genotypes, "compare_genotypes"),
        ("--ignore-indels", ignore_indels, "ignore_indels"),
        ("-a", apply_filters, "apply_filters"),
    ):
        if _switch(value, name):
            arguments.append(flag)
    if regions is not None:
        arguments += ["-r", ",".join(_sequence_of_text(regions, "regions", pattern=re.compile(r"[\w.|-]+(?::\d+-\d+)?")))]
    if samples is not None:
        arguments += ["-s", ",".join(_sample_names(samples, "samples"))]
    if window is not None:
        arguments += ["-w", str(_cli.integer(window, "window", minimum=0))]
    if reference_fasta is not None:
        arguments += ["-R", str(_cli.existing_file(reference_fasta, "reference_fasta"))]
    result = _script_result(
        "vcf-compare",
        [*arguments, *files],
        out,
        script_path=script_path,
        timeout_seconds=timeout_seconds,
    )
    if result["status"] == "completed":
        result.update(_parse_compare(out.read_text(encoding="utf-8", errors="replace")))
    return result


def vcf_indel_stats(
    vcf_path: str,
    output_path: str,
    exons_file: str,
    *,
    script_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """In-frame indel ratio of a VCF given exon coordinates (``vcf-indel-stats``).

    ``exons_file`` is tab-separated ``chromosome, from, to`` (1-based, inclusive). The report is written to the new
    ``output_path``; ``statistics`` holds its ``key: value`` lines, for example ``in-frame`` and ``out-frame``.
    """
    source = _cli.existing_file(vcf_path, "vcf_path")
    out = _output_path(output_path)
    result = _script_result(
        "vcf-indel-stats",
        ["-e", str(_cli.existing_file(exons_file, "exons_file"))],
        out,
        stdin=source,
        script_path=script_path,
        timeout_seconds=timeout_seconds,
    )
    if result["status"] == "completed":
        statistics: dict[str, Any] = {}
        for line in out.read_text(encoding="utf-8", errors="replace").splitlines():
            key, separator, value = line.partition("\t")
            if separator and key.strip():
                statistics[key.strip().rstrip(":")] = _cell(value.strip())
        result["statistics"] = statistics
    return result


def vcf_phased_join(
    vcf_paths: Sequence[str],
    output_path: str,
    *,
    min_join_quality: float | None = None,
    min_phase_quality: float | None = None,
    script_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Join overlapping pre-phased VCF chunks into one VCF with ``vcf-phased-join``.

    Heterozygous calls in the overlaps decide the phase of each chunk. ``vcf_paths`` lists the chunks in genomic
    order. ``min_join_quality`` is the quality threshold for gluing blocks together (-j, default 10) and
    ``min_phase_quality`` breaks pre-phased segments whose PQ value is lower (-q, default 0.6). The joined VCF is
    written to the new ``output_path``.
    """
    files = [str(_cli.existing_file(path, "vcf_paths")) for path in vcf_paths] if not isinstance(vcf_paths, str) else []
    if len(files) < 2:
        raise ConfigurationError("vcf_paths must list at least two chunks.")
    out = _output_path(output_path)
    arguments = ["-o", str(out)]
    if min_join_quality is not None:
        arguments += ["-j", format(_cli.number(min_join_quality, "min_join_quality", minimum=0.0), ".12g")]
    if min_phase_quality is not None:
        arguments += ["-q", format(_cli.number(min_phase_quality, "min_phase_quality", minimum=0.0), ".12g")]
    result = _run_script(
        "vcf-phased-join",
        [*arguments, *files],
        output_path=None,
        extra_outputs=[out],
        script_path=script_path,
        timeout_seconds=timeout_seconds,
    )
    result["output_path"] = str(out)
    if result["status"] == "completed":
        result["records"] = _count_records(out)
    return result


# ---------------------------------------------------------------------------
# fill-aa, fill-fs, fill-ref-md5
# ---------------------------------------------------------------------------

_FILL_TYPES = ("all", "indel", "ref", "snp")


def vcf_fill_ancestral_allele(
    vcf_path: str,
    output_path: str,
    ancestral_prefix: str,
    *,
    variant_types: Sequence[str] = ("all",),
    script_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Add the ancestral allele (INFO/AA) from per-chromosome ancestral FASTA files (``fill-aa``).

    ``ancestral_prefix`` is the path prefix of the gzip-compressed, ``samtools faidx``-indexed files named
    ``<prefix><chromosome>.fa.gz`` (-a). ``variant_types`` selects among ``all``, ``indel``, ``ref`` and ``snp`` (-t).
    Needs ``samtools`` on PATH. The annotated VCF is written to the new ``output_path``.
    """
    source = _cli.existing_file(vcf_path, "vcf_path")
    out = _output_path(output_path)
    prefix = _cli.text(str(ancestral_prefix), "ancestral_prefix")
    kinds = _sequence_of_text(variant_types, "variant_types")
    for kind in kinds:
        _cli.choice(kind, _FILL_TYPES, "variant_types")
    return _script_result(
        "fill-aa",
        ["-a", str(Path(prefix).expanduser().resolve()), "-t", ",".join(kinds)],
        out,
        stdin=source,
        script_path=script_path,
        timeout_seconds=timeout_seconds,
        records=True,
    )


def vcf_fill_flanking_sequence(
    vcf_path: str,
    output_path: str,
    reference_fasta: str,
    *,
    length: int = 100,
    mask_character: str = "N",
    bed_mask_files: Sequence[str] | None = None,
    vcf_mask_files: Sequence[str] | None = None,
    cluster_distance: int | None = None,
    script_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Add the flanking reference sequence (INFO/FS) around every variant (``fill-fs``).

    ``reference_fasta`` is a samtools-faidx-indexed reference (-r). ``length`` is the flank length on each side (-l,
    default 100). ``bed_mask_files`` and ``vcf_mask_files`` are tabix-indexed regions and known variants to mask in
    the flank (-b, -v) with ``mask_character`` (-m; one character or ``"lc"`` for lower case; default ``N``; it is
    passed before every mask, as the script requires, and has no effect without masks); ``cluster_distance``
    self-masks clustered variants within that distance (-c). Needs ``samtools`` on PATH. The annotated VCF is written
    to the new ``output_path``.
    """
    source = _cli.existing_file(vcf_path, "vcf_path")
    out = _output_path(output_path)
    mask = _cli.text(mask_character, "mask_character")
    if mask != "lc" and len(mask) != 1:
        raise ConfigurationError("mask_character must be one character or 'lc'.")
    arguments = ["-r", str(_cli.existing_file(reference_fasta, "reference_fasta"))]
    arguments += ["-l", str(_cli.integer(length, "length", minimum=1))]
    masks: list[str] = []
    for flag, paths, name in (("-b", bed_mask_files, "bed_mask_files"), ("-v", vcf_mask_files, "vcf_mask_files")):
        if paths is not None:
            for path in _files(paths, name):
                masks += ["-m", mask, flag, path]
    if cluster_distance is not None:
        masks += ["-m", mask, "-c", str(_cli.integer(cluster_distance, "cluster_distance", minimum=1))]
    return _script_result(
        "fill-fs",
        [*arguments, *masks, str(source)],
        out,
        script_path=script_path,
        timeout_seconds=timeout_seconds,
        records=True,
    )


def vcf_fill_reference_md5(
    vcf_path: str,
    output_path: str,
    reference_fasta: str,
    *,
    dictionary_path: str | None = None,
    assembly: str | None = None,
    species: str | None = None,
    taxonomy: str | None = None,
    script_path: str | None = None,
    timeout_seconds: float = _TIMEOUT,
) -> dict[str, Any]:
    """Insert the MD5 sums of the reference contigs into the header of a VCF (``fill-ref-md5``).

    Adds ``reference`` and ``contig`` header lines as recommended by VCF 4.1. ``vcf_path`` must be bgzip-compressed and
    tabix-indexed; ``output_path`` is a new ``.vcf.gz``. ``reference_fasta`` is a samtools-faidx-indexed FASTA (-r).
    ``dictionary_path`` reads and appends computed MD5 sums (-d; opened in append mode, so it may exist) and
    ``assembly``, ``species`` and ``taxonomy`` are the optional AS, SP and TX header values (-i). Needs ``tabix``
    and ``samtools`` on PATH.
    """
    source = _cli.existing_file(vcf_path, "vcf_path")
    if not source.name.lower().endswith((".gz", ".bgz")):
        raise ConfigurationError("vcf_path must be a bgzip-compressed, tabix-indexed VCF.")
    out = _output_path(output_path)
    arguments = ["-r", str(_cli.existing_file(reference_fasta, "reference_fasta"))]
    if dictionary_path is not None:
        arguments += ["-d", str(Path(_cli.text(str(dictionary_path), "dictionary_path")).expanduser().resolve())]
    info = []
    for key, value, name in (("AS", assembly, "assembly"), ("SP", species, "species"), ("TX", taxonomy, "taxonomy")):
        if value is not None:
            text = _cli.text(value, name)
            if "," in text or "\n" in text:
                raise ConfigurationError(f"{name} must not contain ',' or a newline.")
            info.append(f"{key}:{text}")
    if info:
        arguments += ["-i", ",".join(info)]
    result = _run_script(
        "fill-ref-md5",
        [*arguments, str(source), str(out)],
        output_path=None,
        extra_outputs=[out],
        script_path=script_path,
        timeout_seconds=timeout_seconds,
    )
    result["output_path"] = str(out)
    if result["status"] == "completed" and out.stat().st_size == 0:
        # fill-ref-md5 pipes `tabix -r` through the shell and ignores its exit code, so a crashed tabix leaves an empty file.
        result["status"] = "failed"
        result["log_tail"] = (result["log_tail"] + "\nfill-ref-md5 wrote an empty file (tabix -r failed).").strip()
    if result["status"] != "completed":
        out.unlink(missing_ok=True)
        Path(str(out) + ".header").unlink(missing_ok=True)
        result["outputs"] = []
    return result


__all__ = [
    "vcf_annotate",
    "vcf_compare",
    "vcf_contrast",
    "vcf_convert",
    "vcf_fill_ancestral_allele",
    "vcf_fill_flanking_sequence",
    "vcf_fill_reference_md5",
    "vcf_fix_ploidy",
    "vcf_indel_stats",
    "vcf_phased_join",
    "vcf_shuffle_columns",
    "vcf_subset",
    "vcf_to_tab",
    "vcftools_compare_vcfs",
    "vcftools_depth",
    "vcftools_export_genotypes",
    "vcftools_extract_fields",
    "vcftools_filter_vcf",
    "vcftools_haplotype_count",
    "vcftools_heterozygosity",
    "vcftools_indel_length_histogram",
    "vcftools_individual_burden",
    "vcftools_linkage_disequilibrium",
    "vcftools_list_sites",
    "vcftools_long_runs_of_homozygosity",
    "vcftools_missingness",
    "vcftools_pca",
    "vcftools_relatedness",
    "vcftools_singletons",
    "vcftools_site_nucleotide_diversity",
    "vcftools_site_quality",
    "vcftools_snp_density",
    "vcftools_transition_transversion",
]
_register_public_verbose(globals())
