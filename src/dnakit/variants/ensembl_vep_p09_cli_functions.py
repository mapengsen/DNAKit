"""Ensembl VEP command-line facades (pack P09): fixed argument lists, no shell, bounded execution.

The functions run the programs of the ensembl-vep distribution (116.3 is the tested release): ``vep``, ``filter_vep``,
``haplo``, ``convert_cache.pl`` (``vep_convert_cache`` in the conda package) and the Python script
``vep_input_format_validator.py``. ``variant_recoder`` is not wrapped: it needs a live Ensembl MySQL connection. The user
installs VEP (``conda install -c bioconda ensembl-vep``, or ``git clone`` of https://github.com/Ensembl/ensembl-vep with
``perl INSTALL.pl``) and either puts the programs on PATH or passes the explicit path of each. Caches, FASTA files and
plugins stay the caller's responsibility: nothing is downloaded. ``vep_annotate_variants`` and
``vep_predict_haplotypes`` never open a database connection unless ``offline=False`` (``--database``) is chosen.

Every function writes into a fresh ``output_directory`` (never overwritten) under documented file names and returns a
JSON-serialisable dictionary: status, version, command, output files, small parsed results and a bounded log tail.
Optional parameters default to ``None`` or ``False`` and are then not passed, so VEP's own defaults apply.
"""

from __future__ import annotations

import contextlib
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, Literal

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.exceptions import (
    BackendExecutionError,
    BackendTimeoutError,
    BackendUnavailableError,
    ConfigurationError,
)

_TESTED_VERSION = "116.3"
_CITATION = "https://doi.org/10.1186/s13059-016-0974-4"
_HINT = (
    "Install Ensembl VEP (for example `conda install -c bioconda ensembl-vep`, or clone "
    "https://github.com/Ensembl/ensembl-vep and run `perl INSTALL.pl`), then put its programs on PATH or pass the "
    "explicit path of the program."
)
_VALIDATOR_HINT = (
    "vep_input_format_validator.py is in the validator/ folder of https://github.com/Ensembl/ensembl-vep. "
    "Download it and pass validator_path."
)
_ROW_LIMIT = 50
_MAX_ARTIFACT_BYTES = 2_000_000_000
_MAX_LOG_BYTES = 5_000_000
_SPECIES = re.compile(r"^[A-Za-z][A-Za-z0-9_]{0,63}$")
_NAME = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9_.-]{0,63}$")
_INPUT_FORMATS = ("ensembl", "vcf", "hgvs", "id", "region", "vep", "guess")
_CUSTOM_FORMATS = ("bed", "vcf", "gff", "gtf", "bigwig")
_PHYLO_SETS = ("7", "17", "20", "30", "46", "100")
_FILTER_CHARACTERS = re.compile(r"^[A-Za-z0-9_ .,:<>=!()\"'+*/%-]+$")


# ---------------------------------------------------------------------------
# Argument helpers
# ---------------------------------------------------------------------------


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _switch(arguments: list[str], flag: str, value: object, name: str) -> None:
    if _flag(value, name):
        arguments.append(flag)


def _opt_int(
    arguments: list[str],
    flag: str,
    value: int | None,
    name: str,
    *,
    minimum: int | None = None,
    maximum: int | None = None,
) -> None:
    if value is not None:
        arguments.extend((flag, str(_cli.integer(value, name, minimum=minimum, maximum=maximum))))


def _opt_text(arguments: list[str], flag: str, value: str | None, name: str) -> None:
    if value is not None:
        arguments.extend((flag, _cli.text(value, name)))


def _binary(arguments: list[str], flag: str, value: int | bool | None, name: str) -> None:
    """A VEP option that takes 0 or 1 (``--shift_hgvs 0``)."""
    if value is None:
        return
    if isinstance(value, bool):
        value = int(value)
    arguments.extend((flag, str(_cli.integer(value, name, minimum=0, maximum=1))))


def _species(value: str) -> str:
    if not isinstance(value, str) or _SPECIES.fullmatch(value) is None:
        raise ConfigurationError("species must be a name such as 'homo_sapiens'.")
    return value


def _file(value: str | Path, name: str) -> str:
    return str(_cli.existing_file(value, name))


def _directory(value: str | Path, name: str) -> str:
    path = Path(value).expanduser().resolve()
    if not path.is_dir():
        raise ConfigurationError(f"{name} must be an existing directory: {path}")
    return str(path)


def _list_item(value: object, name: str) -> str:
    if not isinstance(value, str) or _NAME.fullmatch(value) is None:
        raise ConfigurationError(
            f"{name} entries may contain only letters, digits, '_', '.' and '-'."
        )
    return value


def _sample_list(value: str | Sequence[str], name: str) -> str:
    """``'all'`` or a list of sample names, as VEP's comma-separated value."""
    if value == "all":
        return "all"
    if isinstance(value, str) or not value:
        raise ConfigurationError(f"{name} must be 'all' or a non-empty list of sample names.")
    return ",".join(_list_item(item, name) for item in value)


def _data_source(
    arguments: list[str],
    offline: bool | None,
    cache_directory: str | Path | None,
    cache_version: int | None,
    annotated_by_file: bool,
) -> None:
    """Where VEP gets its transcripts: a cache (offline), a bgzipped GFF/GTF, or the Ensembl database.

    ``offline=None`` chooses: a cache when ``cache_directory`` is given or no GFF/GTF is given (``--offline``), otherwise
    no flag (GFF/GTF mode needs no cache). ``True`` always means ``--offline``; ``False`` means ``--database``.
    """
    if offline is not None:
        _flag(offline, "offline")
    if cache_directory is not None:
        arguments.extend(("--dir_cache", _directory(cache_directory, "cache_directory")))
    if offline is False:
        arguments.append("--database")
    elif offline is True or cache_directory is not None or not annotated_by_file:
        arguments.append("--offline")
    _opt_int(arguments, "--cache_version", cache_version, "cache_version", minimum=1)


def _custom(entry: Mapping[str, object]) -> str:
    """One ``--custom`` value ``file,short_name,format,type,coverage,fields`` from a mapping."""
    unknown = set(entry) - {"file", "short_name", "format", "type", "coverage", "fields"}
    if unknown or "file" not in entry or "short_name" not in entry or "format" not in entry:
        raise ConfigurationError(
            "each custom annotation needs file, short_name and format (optional: type, coverage, fields)."
        )
    path = _file(str(entry["file"]), "custom file")
    if "," in path:
        raise ConfigurationError("custom file paths must not contain commas.")
    short_name = _list_item(entry["short_name"], "short_name")
    file_format = _cli.choice(entry["format"], _CUSTOM_FORMATS, "custom format")
    overlap = _cli.choice(entry.get("type", "exact"), ("exact", "overlap"), "custom type")
    coverage = "1" if entry.get("coverage") else "0"
    fields = entry.get("fields")
    parts = [path, short_name, file_format, overlap, coverage]
    if fields is not None:
        if isinstance(fields, str) or not fields:
            raise ConfigurationError("custom fields must be a non-empty list of field names.")
        parts.append("%".join(_list_item(item, "custom fields") for item in fields))
    return ",".join(parts)


# ---------------------------------------------------------------------------
# Running programs
# ---------------------------------------------------------------------------


def _tail(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")[-_cli.LOG_TAIL_CHARS :]


def _stream(
    executable: str,
    arguments: Sequence[str],
    directory: Path,
    stdout_path: Path,
    stderr_path: Path,
    timeout_seconds: float,
) -> int:
    """Run a program without a shell, streaming standard output and standard error into files."""
    timeout = _cli.number(timeout_seconds, "timeout_seconds", minimum=1.0, maximum=86_400.0)
    with stdout_path.open("xb") as stdout, stderr_path.open("xb") as stderr:
        try:
            process = subprocess.Popen(
                [executable, *arguments],
                stdin=subprocess.DEVNULL,
                stdout=stdout,
                stderr=stderr,
                cwd=directory,
                shell=False,
                env={**os.environ, "LANG": "C", "LC_ALL": "C"},
                start_new_session=os.name == "posix",
            )
        except OSError as exc:
            raise BackendExecutionError(
                "External command could not be started.",
                code="VEP_START_FAILED",
                context={"error": str(exc)},
            ) from exc
        deadline = time.monotonic() + timeout
        try:
            while process.poll() is None:
                if time.monotonic() > deadline:
                    raise BackendTimeoutError(f"{Path(executable).name} exceeded its time limit.")
                if (
                    stdout_path.stat().st_size > _MAX_ARTIFACT_BYTES
                    or stderr_path.stat().st_size > _MAX_LOG_BYTES
                ):
                    raise BackendExecutionError(
                        f"{Path(executable).name} exceeded its output size limit."
                    )
                time.sleep(0.02)
        finally:
            if process.poll() is None:
                if os.name == "posix":
                    with contextlib.suppress(ProcessLookupError):
                        os.killpg(process.pid, signal.SIGKILL)
                else:
                    process.kill()
            process.wait()
    return int(process.returncode)


def _banner_version(executable: str, directory: Path) -> str | None:
    """The ensembl-vep release printed in the ``--help`` banner (``ensembl-vep : 116.3``)."""
    from dnakit.backends.external import execute_bounded_command
    from dnakit.exceptions import DNAKitError

    try:
        probe = execute_bounded_command(
            executable,
            ["--help"],
            backend_id=Path(executable).name,
            cwd=directory,
            timeout_seconds=120.0,
            max_output_bytes=2_000_000,
            verbose=False,
        )
    except DNAKitError:
        return None
    found = re.search(r"ensembl-vep\s*:\s*(\d+(?:\.\d+)*)", probe.output or "")
    return found.group(1) if found else None


def _execute(
    operation: str,
    program: str,
    explicit: str | None,
    arguments: Sequence[str],
    directory: Path,
    outputs: dict[str, Path],
    timeout_seconds: float,
    **extra: object,
) -> dict[str, Any]:
    """Run a VEP program that writes its own output files; failure raises ``BackendExecutionError``."""
    timeout = _cli.number(timeout_seconds, "timeout_seconds", minimum=1.0, maximum=86_400.0)
    executable = _cli.program(program, explicit=explicit, hint=_HINT)
    run = _cli.run(
        program,
        arguments,
        cwd=directory,
        expected=list(outputs.values()),
        explicit=executable,
        hint=_HINT,
        timeout_seconds=timeout,
    )
    if run["status"] != "completed":
        raise BackendExecutionError(
            f"{program} {operation} failed (exit code {run['return_code']}).",
            code="VEP_EXECUTION_FAILED",
            context={"operation": operation, "log_tail": run["log_tail"]},
        )
    warnings = directory / "warnings.txt"
    if warnings.is_file():
        outputs = {**outputs, "warnings": warnings}
    return _result(
        operation,
        program,
        _banner_version(executable, directory),
        run["command"],
        directory,
        outputs,
        run["log_tail"],
        **extra,
    )


def _result(
    operation: str,
    backend: str,
    version: str | None,
    command: Sequence[str],
    directory: Path,
    outputs: dict[str, Path],
    log_tail: str,
    **extra: object,
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "status": "completed",
        "backend": backend,
        "operation": operation,
        "version": version,
        "command": list(command),
        "output_directory": str(directory),
        "outputs": {name: str(path) for name, path in outputs.items()},
        "log_tail": log_tail[-_cli.LOG_TAIL_CHARS :],
        "provenance": {
            "tool": "Ensembl VEP",
            "tested_release": _TESTED_VERSION,
            "citation": _CITATION,
        },
    }
    result.update(extra)
    return result


def _parse_annotations(path: Path, output_format: str, info_field: str = "CSQ") -> dict[str, Any]:
    """Columns, first records and the number of annotated variants of a VEP result file (vcf, tab or json)."""
    records: list[dict[str, Any]] = []
    columns: list[str] = []
    count = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if output_format == "json":
                count += 1
                if len(records) < _ROW_LIMIT:
                    records.append(json.loads(line))
            elif output_format == "tab":
                if line.startswith("#Uploaded_variation"):
                    columns = line[1:].split("\t")
                elif not line.startswith("#"):
                    count += 1
                    if len(records) < _ROW_LIMIT:
                        records.append(dict(zip(columns, line.split("\t"), strict=False)))
            else:
                if line.startswith("##INFO=<ID=" + info_field + ","):
                    match = re.search(r"Format: ([^\">]+)", line)
                    columns = match.group(1).split("|") if match else []
                elif not line.startswith("#"):
                    count += 1
                    if len(records) < _ROW_LIMIT:
                        fields = line.split("\t")
                        info = next(
                            (
                                item
                                for item in fields[7].split(";")
                                if item.startswith(info_field + "=")
                            ),
                            None,
                        )
                        entries = [] if info is None else info[len(info_field) + 1 :].split(",")
                        records.append(
                            {
                                "chrom": fields[0],
                                "pos": int(fields[1]),
                                "ref": fields[3],
                                "alt": fields[4],
                                "annotations": [
                                    dict(zip(columns, entry.split("|"), strict=False))
                                    for entry in entries
                                ],
                            }
                        )
    return {"columns": columns, "records": records, "record_count": count}


# ---------------------------------------------------------------------------
# vep
# ---------------------------------------------------------------------------


def vep_annotate_variants(
    input_path: str | Path,
    output_directory: str | Path,
    *,
    output_format: Literal["vcf", "tab", "json"] = "vcf",
    input_format: Literal["ensembl", "vcf", "hgvs", "id", "region", "vep", "guess"] | None = None,
    species: str = "homo_sapiens",
    assembly: str | None = None,
    offline: bool | None = None,
    cache_directory: str | Path | None = None,
    cache_version: int | None = None,
    fasta: str | Path | None = None,
    gff: str | Path | None = None,
    gtf: str | Path | None = None,
    bigwig: str | Path | None = None,
    bam: str | Path | None = None,
    custom_annotations: Sequence[Mapping[str, object]] | None = None,
    summary: bool = False,
    variant_class: bool = False,
    hgvs: bool = False,
    hgvsg: bool = False,
    spdi: bool = False,
    shift_hgvs: int | None = None,
    shift_3prime: int | None = None,
    shift_genomic: int | None = None,
    check_ref: bool = False,
    lookup_ref: bool = False,
    individual: str | Sequence[str] | None = None,
    individual_zyg: str | Sequence[str] | None = None,
    check_existing: bool = False,
    af: bool = False,
    af_1kg: bool = False,
    af_gnomade: bool = False,
    af_gnomadg: bool = False,
    max_af: bool = False,
    check_frequency: bool = False,
    freq_pop: str | None = None,
    freq_freq: float | None = None,
    freq_gt_lt: Literal["gt", "lt"] | None = None,
    freq_filter: Literal["exclude", "include"] | None = None,
    filter_common: bool = False,
    clin_sig_allele: int | None = None,
    clinvar_somatic_classification: bool = False,
    gene_phenotype: bool = False,
    mirna: bool = False,
    regulatory: bool = False,
    cell_type: str | None = None,
    phylop: Sequence[str] | None = None,
    phastcons: Sequence[str] | None = None,
    ucsc_assembly: str | None = None,
    ucsc_data_root: str | None = None,
    domains: bool = False,
    numbers: bool = False,
    nearest: Literal["transcript", "gene", "symbol"] | None = None,
    stats: Literal["html", "text"] | None = None,
    fork: int | None = None,
    vep_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict variant consequences with Ensembl VEP (``vep``), with the annotation options as typed parameters.

    One function for the options of a ``vep`` run: ``summary`` (``--summary``, one line per variant with all
    consequences), ``variant_class``, ``hgvs``/``hgvsg``/``spdi`` notations, 3-prime shifting (``shift_hgvs``,
    ``shift_3prime``, ``shift_genomic``, each 0 or 1), ``check_ref``/``lookup_ref``, genotypes per individual
    (``individual``, ``individual_zyg``: ``'all'`` or sample names), co-located known variants (``check_existing``) with
    their allele frequencies (``af``, ``af_1kg``, ``af_gnomade``, ``af_gnomadg``, ``max_af``) and frequency filtering
    (``check_frequency`` with ``freq_pop``, ``freq_freq``, ``freq_gt_lt``, ``freq_filter``; ``filter_common``), ClinVar
    (``clin_sig_allele``, ``clinvar_somatic_classification``), ``gene_phenotype``, ``mirna``, ``regulatory`` (with
    ``cell_type``), conservation scores from UCSC bigWig files (``phylop``/``phastcons`` = sets such as ``'100'`` with
    ``ucsc_assembly``; ``ucsc_data_root`` points to a local copy), user tracks (``custom_annotations`` = mappings with
    ``file``, ``short_name``, ``format`` of bed/vcf/gff/gtf/bigwig, optional ``type``, ``coverage``, ``fields``; ``gff``,
    ``gtf`` and ``bigwig`` are VEP's shortcuts), transcript sequence correction from a ``bam`` file, ``domains``,
    ``numbers`` (exon and intron numbers), ``nearest`` and a ``stats`` report (``html`` or ``text``). ``output_format``
    is ``vcf``, ``tab`` or ``json``. Consequences need transcripts: a VEP cache (``cache_directory``, ``cache_version``;
    ``offline=True`` or the default) or a bgzipped, tabix-indexed ``gff``/``gtf`` with the reference ``fasta`` (no cache is
    needed then); ``offline=False`` reads the Ensembl database instead. VEP indexes ``fasta`` next to the file. Writes ``annotations.<vcf|tsv|json>`` and,
    with ``stats``, ``stats.<html|txt>``; the result has ``columns``, ``records`` (first rows; VCF rows carry the parsed
    CSQ entries) and ``record_count``. Options that need a database or a network connection (frequencies, ClinVar,
    UCSC scores) are passed through unchanged; they run only when the cache or the connection provides the data.
    """
    _cli.choice(output_format, ("vcf", "tab", "json"), "output_format")
    source = _file(input_path, "input_path")
    species = _species(species)
    arguments = ["--input_file", source, "--species", species]
    if input_format is not None:
        arguments.extend(("--format", _cli.choice(input_format, _INPUT_FORMATS, "input_format")))
    _opt_text(arguments, "--assembly", assembly, "assembly")
    _data_source(
        arguments, offline, cache_directory, cache_version, gff is not None or gtf is not None
    )
    for flag, value, name in (
        ("--fasta", fasta, "fasta"),
        ("--gff", gff, "gff"),
        ("--gtf", gtf, "gtf"),
        ("--bigwig", bigwig, "bigwig"),
        ("--bam", bam, "bam"),
    ):
        if value is not None:
            arguments.extend((flag, _file(value, name)))
    for entry in custom_annotations or ():
        arguments.extend(("--custom", _custom(entry)))
    for flag, value in (
        ("--summary", summary),
        ("--variant_class", variant_class),
        ("--hgvs", hgvs),
        ("--hgvsg", hgvsg),
        ("--spdi", spdi),
        ("--check_ref", check_ref),
        ("--lookup_ref", lookup_ref),
        ("--check_existing", check_existing),
        ("--af", af),
        ("--af_1kg", af_1kg),
        ("--af_gnomade", af_gnomade),
        ("--af_gnomadg", af_gnomadg),
        ("--max_af", max_af),
        ("--check_frequency", check_frequency),
        ("--filter_common", filter_common),
        ("--clinvar_somatic_classification", clinvar_somatic_classification),
        ("--gene_phenotype", gene_phenotype),
        ("--mirna", mirna),
        ("--regulatory", regulatory),
        ("--domains", domains),
        ("--numbers", numbers),
    ):
        _switch(arguments, flag, value, flag[2:])
    _binary(arguments, "--shift_hgvs", shift_hgvs, "shift_hgvs")
    _binary(arguments, "--shift_3prime", shift_3prime, "shift_3prime")
    _binary(arguments, "--shift_genomic", shift_genomic, "shift_genomic")
    _binary(arguments, "--clin_sig_allele", clin_sig_allele, "clin_sig_allele")
    if individual is not None:
        arguments.extend(("--individual", _sample_list(individual, "individual")))
    if individual_zyg is not None:
        arguments.extend(("--individual_zyg", _sample_list(individual_zyg, "individual_zyg")))
    _opt_text(arguments, "--freq_pop", freq_pop, "freq_pop")
    if freq_freq is not None:
        arguments.extend(
            (
                "--freq_freq",
                format(_cli.number(freq_freq, "freq_freq", minimum=0.0, maximum=1.0), ".12g"),
            )
        )
    if freq_gt_lt is not None:
        arguments.extend(("--freq_gt_lt", _cli.choice(freq_gt_lt, ("gt", "lt"), "freq_gt_lt")))
    if freq_filter is not None:
        arguments.extend(
            ("--freq_filter", _cli.choice(freq_filter, ("exclude", "include"), "freq_filter"))
        )
    _opt_text(arguments, "--cell_type", cell_type, "cell_type")
    for flag, values in (("--phyloP", phylop), ("--phastCons", phastcons)):
        for item in values or ():
            arguments.extend((flag, _cli.choice(str(item), _PHYLO_SETS, flag[2:])))
    if (phylop or phastcons) and ucsc_assembly is None:
        raise ConfigurationError("phylop and phastcons need ucsc_assembly (for example 'hg38').")
    _opt_text(arguments, "--ucsc_assembly", ucsc_assembly, "ucsc_assembly")
    _opt_text(arguments, "--ucsc_data_root", ucsc_data_root, "ucsc_data_root")
    if nearest is not None:
        arguments.extend(
            ("--nearest", _cli.choice(nearest, ("transcript", "gene", "symbol"), "nearest"))
        )
    _opt_int(arguments, "--fork", fork, "fork", minimum=1, maximum=64)
    _cli.program("vep", explicit=vep_path, hint=_HINT)
    directory = _cli.new_directory(output_directory, "output_directory")
    extension = {"vcf": "vcf", "tab": "tsv", "json": "json"}[output_format]
    outputs = {"output": directory / f"annotations.{extension}"}
    arguments.extend(
        (
            "--output_file",
            str(outputs["output"]),
            f"--{output_format}",
            "--warning_file",
            str(directory / "warnings.txt"),
        )
    )
    if stats is None:
        arguments.append("--no_stats")
    else:
        _cli.choice(stats, ("html", "text"), "stats")
        outputs["stats"] = directory / ("stats.html" if stats == "html" else "stats.txt")
        arguments.extend(("--stats_file", str(outputs["stats"]), f"--stats_{stats}"))
    arguments.extend(("--quiet", "--no_progress"))
    result = _execute(
        "annotate_variants",
        "vep",
        vep_path,
        arguments,
        directory,
        outputs,
        timeout_seconds,
        output_format=output_format,
    )
    result.update(_parse_annotations(outputs["output"], output_format))
    return result


def vep_filter_results(
    input_path: str | Path,
    output_directory: str | Path,
    *,
    filters: Sequence[str] | None = None,
    list_fields: bool = False,
    only_matched: bool = False,
    soft_filter: bool = False,
    count_only: bool = False,
    input_format: Literal["vcf", "tab"] | None = None,
    vcf_info_field: str | None = None,
    start: int | None = None,
    limit: int | None = None,
    gzipped: bool = False,
    filter_vep_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter VEP-annotated variant records by expression (``filter_vep``).

    ``filters`` are ``filter_vep`` expressions such as ``"Consequence is missense_variant and AF < 0.01"`` (several
    expressions are alternatives, as with repeated ``-f``); ``list_fields`` runs ``--list`` and returns the field names
    that can be filtered on instead. ``only_matched`` rewrites the VCF CSQ field with the matching entries only,
    ``soft_filter`` flags instead of dropping, ``count_only`` returns just the number of matching records (``count``),
    ``input_format`` is ``--format``, ``vcf_info_field`` is ``--vcf_info_field``, ``start`` and ``limit`` are ``--start``
    and ``--limit``, and ``gzipped`` is ``--gz``. The ontology option (which needs a database) is not exposed. Writes
    ``filtered.txt`` (VCF and tab results keep their format); with ``list_fields`` or ``count_only`` the result file is
    ``stdout.txt`` when the program printed instead of writing the file.
    """
    source = _file(input_path, "input_path")
    expressions = list(filters or ())
    if not expressions and not list_fields:
        raise ConfigurationError("Give at least one filter expression, or set list_fields=True.")
    for expression in expressions:
        if (
            not isinstance(expression, str)
            or not expression.strip()
            or _FILTER_CHARACTERS.fullmatch(expression) is None
        ):
            raise ConfigurationError(
                f"filter expression contains unsupported characters: {expression!r}"
            )
        if expression.lstrip().startswith("-"):
            raise ConfigurationError("a filter expression must not start with '-'.")
    arguments = ["--input_file", source]
    for expression in expressions:
        arguments.extend(("--filter", expression))
    _switch(arguments, "--list", list_fields, "list_fields")
    _switch(arguments, "--only_matched", only_matched, "only_matched")
    _switch(arguments, "--soft_filter", soft_filter, "soft_filter")
    _switch(arguments, "--count", count_only, "count_only")
    _switch(arguments, "--gz", gzipped, "gzipped")
    if input_format is not None:
        arguments.extend(("--format", _cli.choice(input_format, ("vcf", "tab"), "input_format")))
    _opt_text(arguments, "--vcf_info_field", vcf_info_field, "vcf_info_field")
    _opt_int(arguments, "--start", start, "start", minimum=1)
    _opt_int(arguments, "--limit", limit, "limit", minimum=1)
    executable = _cli.program("filter_vep", explicit=filter_vep_path, hint=_HINT)
    directory = _cli.new_directory(output_directory, "output_directory")
    output, log = directory / "filtered.txt", directory / "tool.log"
    arguments.extend(("--output_file", str(output)))
    return_code = _stream(
        executable, arguments, directory, directory / "stdout.txt", log, timeout_seconds
    )
    if return_code != 0:
        raise BackendExecutionError(
            f"filter_vep failed (exit code {return_code}).",
            code="VEP_EXECUTION_FAILED",
            context={"log_tail": _tail(log)},
        )
    stdout_path = directory / "stdout.txt"
    # --list prints to standard output, --count to the --output_file; the other modes write the --output_file
    source_file = output
    if (list_fields or count_only) and not output.read_text(
        encoding="utf-8", errors="replace"
    ).strip():
        source_file = stdout_path
    text = source_file.read_text(encoding="utf-8", errors="replace")
    extra: dict[str, Any] = {}
    if list_fields:
        extra["fields"] = [
            line.strip()
            for line in text.splitlines()
            if line.strip() and not line.startswith("#") and not line.startswith("Available fields")
        ]
    elif count_only:
        found = re.search(r"(\d+)", text)
        extra["count"] = int(found.group(1)) if found else None
    else:
        data_lines = [line for line in text.splitlines() if line and not line.startswith("#")]
        extra["record_count"] = len(data_lines)
        extra["records"] = [line.split("\t") for line in data_lines[:_ROW_LIMIT]]
    outputs = {"output": source_file, "log": log}
    return _result(
        "filter_results",
        "filter_vep",
        _banner_version(executable, directory),
        ["filter_vep", *arguments],
        directory,
        outputs,
        _tail(log),
        **extra,
    )


def vep_predict_haplotypes(
    input_path: str | Path,
    output_directory: str | Path,
    *,
    input_format: Literal["vcf", "ensembl", "region", "vep", "guess"] | None = None,
    species: str = "homo_sapiens",
    assembly: str | None = None,
    offline: bool | None = None,
    cache_directory: str | Path | None = None,
    cache_version: int | None = None,
    fasta: str | Path | None = None,
    gff: str | Path | None = None,
    gtf: str | Path | None = None,
    haplotype_frequencies: str | Path | None = None,
    phased: bool = False,
    chromosomes: str | None = None,
    haplo_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict the consequences of whole transcript haplotypes in phased genotypes (``haplo``).

    ``haplo`` reconstructs the transcript sequence of every sample haplotype and reports its amino-acid changes, instead
    of one variant at a time. The input is a VCF with genotypes (``phased=True`` forces phased interpretation);
    annotation comes from a VEP cache (``cache_directory``) or a bgzipped, tabix-indexed ``gff``/``gtf`` with ``fasta``.
    ``haplotype_frequencies`` is ``--haplotype_frequencies`` and ``chromosomes`` is ``--chr`` (for example ``1-5,10``).
    Writes ``haplotypes.json`` (JSON output, one object per transcript haplotype set) and parses it into ``records``.
    """
    source = _file(input_path, "input_path")
    species = _species(species)
    arguments = ["--input_file", source, "--species", species, "--json"]
    if input_format is not None:
        arguments.extend(
            (
                "--format",
                _cli.choice(
                    input_format, ("vcf", "ensembl", "region", "vep", "guess"), "input_format"
                ),
            )
        )
    _opt_text(arguments, "--assembly", assembly, "assembly")
    _data_source(
        arguments, offline, cache_directory, cache_version, gff is not None or gtf is not None
    )
    for flag, value, name in (
        ("--fasta", fasta, "fasta"),
        ("--gff", gff, "gff"),
        ("--gtf", gtf, "gtf"),
        ("--haplotype_frequencies", haplotype_frequencies, "haplotype_frequencies"),
    ):
        if value is not None:
            arguments.extend((flag, _file(value, name)))
    _switch(arguments, "--phased", phased, "phased")
    _opt_text(arguments, "--chr", chromosomes, "chromosomes")
    _cli.program("haplo", explicit=haplo_path, hint=_HINT)
    directory = _cli.new_directory(output_directory, "output_directory")
    outputs = {"output": directory / "haplotypes.json"}
    arguments.extend(
        (
            "--output_file",
            str(outputs["output"]),
            "--warning_file",
            str(directory / "warnings.txt"),
            "--quiet",
        )
    )
    result = _execute(
        "predict_haplotypes", "haplo", haplo_path, arguments, directory, outputs, timeout_seconds
    )
    rows: list[Any] = []
    count = 0
    with outputs["output"].open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.strip():
                count += 1
                if len(rows) < _ROW_LIMIT:
                    rows.append(json.loads(line))
    result["records"], result["record_count"] = rows, count
    return result


def vep_validate_input_format(
    input_path: str | Path,
    *,
    validator_path: str | Path | None = None,
    timeout_seconds: float = 600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check a file in VEP's default input format with Ensembl's ``vep_input_format_validator.py``.

    The script checks the column count, the allele characters and the ``REF/ALT`` form, the strand, that the coordinates
    match the allele length, that insertions have start greater than end, and that variants are ordered by position and
    chromosome block. It is run with the current Python interpreter. ``validator_path`` is the script (from the
    ``validator/`` folder of the ensembl-vep repository), or it is found as ``vep_input_format_validator.py`` on PATH.
    Nothing is written; the result has ``valid``, ``error_count`` and ``errors`` (line number, message, line content).
    """
    source = _file(input_path, "input_path")
    if validator_path is None:
        script = _cli.program("vep_input_format_validator.py", explicit=None, hint=_VALIDATOR_HINT)
    else:
        candidate = Path(validator_path).expanduser()
        if not candidate.is_file():
            raise BackendUnavailableError(
                f"validator script does not exist: {candidate}", hint=_VALIDATOR_HINT
            )
        script = str(candidate.resolve())
    timeout = _cli.number(timeout_seconds, "timeout_seconds", minimum=1.0, maximum=86_400.0)
    from dnakit.backends.external import execute_bounded_command

    work = Path(source).parent
    run = execute_bounded_command(
        sys.executable,
        [script, source],
        backend_id="vep_input_format_validator",
        cwd=work,
        timeout_seconds=timeout,
        max_output_bytes=_cli.MAX_OUTPUT_BYTES,
        verbose=False,
    )
    text = run.output or ""
    if run.return_code != 0:
        raise BackendExecutionError(
            f"vep_input_format_validator.py failed (exit code {run.return_code}).",
            code="VEP_EXECUTION_FAILED",
            context={"log_tail": text[-_cli.LOG_TAIL_CHARS :]},
        )
    errors = [
        {"line": int(match.group(1)), "message": match.group(2), "content": match.group(3)}
        for match in re.finditer(r'^Line (\d+): (.*?) \| Line content: "(.*)"$', text, re.MULTILINE)
    ]
    summary = re.search(r"You have (\d+) lines not valid", text)
    valid = "Your VEP input file is valid" in text
    return {
        "status": "completed",
        "backend": "vep_input_format_validator",
        "operation": "validate_input_format",
        "version": None,
        "command": [Path(sys.executable).name, Path(script).name, source],
        "valid": valid,
        "error_count": int(summary.group(1)) if summary else 0,
        "errors": errors,
        "log_tail": text[-_cli.LOG_TAIL_CHARS :],
        "provenance": {
            "tool": "Ensembl VEP",
            "tested_release": _TESTED_VERSION,
            "citation": _CITATION,
        },
    }


def vep_convert_cache(
    cache_directory: str | Path,
    *,
    species: str,
    version: str,
    compress_command: str | None = None,
    use_sereal: bool = False,
    remove_original: bool = False,
    fork: int | None = None,
    convert_cache_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a downloaded VEP cache so that it can be read through tabix (``convert_cache.pl``).

    This modifies the cache in place: it writes ``.gz`` variation files with tabix indexes next to the originals and, with
    ``remove_original=True`` (``--remove``), deletes the files it converted. ``cache_directory`` is the folder that holds
    ``<species>/<version>``; ``species`` and ``version`` (for example ``115_GRCh38``) select what is converted (a
    comma-separated list or ``all``). ``use_sereal`` is ``--sereal`` (also converts the transcript cache to Sereal),
    ``compress_command`` is ``--compress`` (default ``gzip -dc``) and ``fork`` is ``--fork``. ``bgzip`` and ``tabix`` must be
    on PATH. Nothing is returned except the tool log and the status.
    """
    cache = _directory(cache_directory, "cache_directory")
    for value, name in ((species, "species"), (version, "version")):
        if value != "all" and not all(_NAME.fullmatch(part) for part in str(value).split(",")):
            raise ConfigurationError(f"{name} must be 'all' or a comma-separated list of names.")
    arguments = ["--dir", cache, "--species", species, "--version", version, "--force_overwrite"]
    if compress_command is not None:
        if not re.fullmatch(r"[A-Za-z0-9_.-]+( -[A-Za-z]+)?", compress_command):
            raise ConfigurationError(
                "compress_command must be a program name with at most one short option, such as 'gzip -dc'."
            )
        arguments.extend(("--compress", compress_command))
    _switch(arguments, "--sereal", use_sereal, "use_sereal")
    _switch(arguments, "--remove", remove_original, "remove_original")
    _opt_int(arguments, "--fork", fork, "fork", minimum=1, maximum=64)
    timeout = _cli.number(timeout_seconds, "timeout_seconds", minimum=1.0, maximum=86_400.0)
    name = "convert_cache.pl"
    if (
        convert_cache_path is None
        and shutil.which(name) is None
        and shutil.which("vep_convert_cache") is not None
    ):
        name = "vep_convert_cache"  # the name used by the conda package
    run = _cli.run(
        name,
        arguments,
        cwd=Path(cache),
        expected=[],
        explicit=convert_cache_path,
        hint=_HINT,
        timeout_seconds=timeout,
    )
    if run["status"] != "completed":
        raise BackendExecutionError(
            f"convert_cache.pl failed (exit code {run['return_code']}).",
            code="VEP_EXECUTION_FAILED",
            context={"log_tail": run["log_tail"]},
        )
    return {
        "status": "completed",
        "backend": name,
        "operation": "convert_cache",
        "version": run["version"],
        "command": run["command"],
        "cache_directory": cache,
        "log_tail": run["log_tail"],
        "provenance": {
            "tool": "Ensembl VEP",
            "tested_release": _TESTED_VERSION,
            "citation": _CITATION,
        },
    }


__all__ = [
    "vep_annotate_variants",
    "vep_convert_cache",
    "vep_filter_results",
    "vep_predict_haplotypes",
    "vep_validate_input_format",
]

_register_public_verbose(globals(), tuple(__all__))
