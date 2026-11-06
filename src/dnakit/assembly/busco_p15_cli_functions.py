"""BUSCO command-line facade: genome, transcriptome and protein completeness with a chosen gene predictor.

BUSCO (5.x, tested against 5.8.2) scores the completeness of an assembly against the universal single-copy orthologs of a
lineage dataset. The user installs BUSCO and its programs (``conda install -c bioconda -c conda-forge busco``) and either has
BUSCO download the lineage dataset (``download_path``) or supplies a downloaded lineage folder and works ``offline``. The
function runs ``busco`` with a fixed argument whitelist and no shell in a new output folder and returns the parsed short
summary, the assembly statistics and the per-gene table.

``dnakit.assembly.external.run_busco`` runs BUSCO with its defaults; this module exposes the gene predictor, the assembly
statistics and the scaffold composition.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p15_cli as _cli
from dnakit.exceptions import ConfigurationError

_TOOL = _cli.Tool(
    name="busco",
    executable="busco",
    hint=(
        "Install BUSCO 5.x with its dependencies (for example `conda install -c bioconda -c conda-forge busco`) and put busco on PATH "
        "or pass busco_path. The lineage dataset is downloaded by BUSCO into download_path, or give a downloaded lineage folder."
    ),
    tested_version="5.8.2",
    citation="https://doi.org/10.1093/molbev/msab199",
    version_arguments=("--version",),
    version_pattern=r"BUSCO\s+(\d+\.\d+(?:\.\d+)?)",
)
_MODES = ("genome", "transcriptome", "proteins")
_PREDICTORS = ("prodigal", "augustus", "metaeuk", "miniprot")
_SCORE = re.compile(
    r"C:(?P<complete>[\d.]+)%\[S:(?P<single>[\d.]+)%,D:(?P<duplicated>[\d.]+)%\],F:(?P<fragmented>[\d.]+)%,M:(?P<missing>[\d.]+)%,n:(?P<n>\d+)"
)
_STAT_LINE = re.compile(r"^\s*(?P<value>[\d.]+%?)\s+(?P<key>[A-Za-z][^\n]*?)\s*$")


def _parse_summary(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8", errors="replace")
    result: dict[str, Any] = {}
    score = _SCORE.search(text)
    if score:
        groups = score.groupdict()
        result["scores_percent"] = {
            key: float(value) for key, value in groups.items() if key != "n"
        }
        result["marker_count"] = int(groups["n"])
    counts: dict[str, int] = {}
    statistics: dict[str, Any] = {}
    in_stats = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("Assembly Statistics"):
            in_stats = True
            continue
        match = _STAT_LINE.match(line)
        if not match:
            continue
        key = match.group("key").strip()
        value = match.group("value")
        number: Any = (
            float(value[:-1])
            if value.endswith("%")
            else (float(value) if "." in value else int(value))
        )
        if in_stats:
            statistics[re.sub(r"[^a-z0-9]+", "_", key.lower()).strip("_")] = number
        elif re.match(r"^Complete|^Fragmented|^Missing|^Total BUSCO", key):
            counts[re.sub(r"[^a-z0-9]+", "_", key.lower()).strip("_")] = number
    result["counts"] = counts
    result["assembly_statistics"] = statistics
    return result


def busco_assess_completeness(
    input_path: str,
    lineage: str,
    output_directory: str,
    *,
    mode: str = "genome",
    gene_predictor: str | None = None,
    output_name: str = "busco",
    download_path: str | None = None,
    offline: bool = False,
    auto_lineage: str | None = None,
    evalue: float = 0.001,
    limit: int = 3,
    long_mode: bool = False,
    contig_break: int = 10,
    scaffold_composition: bool = False,
    skip_assembly_statistics: bool = False,
    augustus_species: str | None = None,
    cpus: int = 1,
    busco_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Assess genome or transcriptome completeness with BUSCO markers and report assembly statistics (``busco``).

    ``mode`` is ``'genome'``, ``'transcriptome'`` or ``'proteins'``. In genome mode ``gene_predictor`` selects how genes are
    found: ``'prodigal'`` (prokaryotic lineages; BUSCO's own choice there, no flag is needed), ``'augustus'`` (eukaryotes,
    ``--augustus``, optionally with ``augustus_species``), ``'metaeuk'`` (eukaryotes, ``--metaeuk``) or ``'miniprot'``
    (``--miniprot``, protein-to-genome alignment); ``None`` leaves BUSCO's default. Transcriptome mode uses tblastn plus HMMER
    for prokaryotic lineages and MetaEuk for eukaryotic ones. ``lineage`` is a dataset name such as ``bacteria_odb10`` or a
    local folder; ``auto_lineage`` (``'all'``, ``'prok'`` or ``'euk'``) instead lets BUSCO choose the lineage. ``offline`` uses
    only datasets already in ``download_path``. ``scaffold_composition`` writes the per-scaffold A, C, G, T, N composition
    (``--scaffold_composition``); ``skip_assembly_statistics`` skips the BBTools assembly statistics (``--skip_bbtools``).
    The run goes into the new folder ``output_directory``/``output_name``. Result keys: ``summary_file``, ``scores_percent``
    (complete, single, duplicated, fragmented, missing), ``marker_count``, ``counts``, ``assembly_statistics`` (scaffold and
    contig counts, total length, percent gaps, N50 and so on), ``full_table`` rows and ``scaffold_composition`` rows.
    """
    source = _cli.existing_file(input_path, "input_path")
    chosen_mode = _cli.choice(mode, _MODES, "mode")
    out = _cli.new_directory(output_directory, "output_directory")
    name = _cli.text(output_name, "output_name")
    if "/" in name:
        raise ConfigurationError("output_name must be a simple name.")
    arguments = [
        "-i",
        str(source),
        "-m",
        chosen_mode,
        "-o",
        name,
        "--out_path",
        str(out),
        "-c",
        str(_cli.threads_value(cpus)),
        "--quiet",
    ]
    if auto_lineage is not None:
        arguments.append(
            {"all": "--auto-lineage", "prok": "--auto-lineage-prok", "euk": "--auto-lineage-euk"}[
                _cli.choice(auto_lineage, ("all", "prok", "euk"), "auto_lineage")
            ]
        )
    else:
        lineage_value = _cli.text(lineage, "lineage")
        arguments += ["-l", lineage_value]
    if gene_predictor is not None:
        predictor = _cli.choice(gene_predictor, _PREDICTORS, "gene_predictor")
        if chosen_mode != "genome" and predictor != "metaeuk":
            raise ConfigurationError(
                "gene_predictor applies to mode='genome' (MetaEuk is also used for eukaryotic transcriptomes)."
            )
        if predictor != "prodigal":
            arguments.append(f"--{predictor}")
    if augustus_species is not None:
        if gene_predictor != "augustus":
            raise ConfigurationError("augustus_species needs gene_predictor='augustus'.")
        arguments += ["--augustus_species", _cli.text(augustus_species, "augustus_species")]
    arguments += [
        "--evalue",
        format(_cli.number(evalue, "evalue", minimum=0.0), ".12g"),
        "--limit",
        str(_cli.integer(limit, "limit", minimum=1, maximum=20)),
    ]
    arguments += ["--contig_break", str(_cli.integer(contig_break, "contig_break", minimum=1))]
    if download_path is not None:
        arguments += [
            "--download_path",
            str(_cli.existing_directory(download_path, "download_path")),
        ]
    for enabled, option in (
        (offline, "--offline"),
        (long_mode, "--long"),
        (scaffold_composition, "--scaffold_composition"),
        (skip_assembly_statistics, "--skip_bbtools"),
    ):
        if _cli.flag(enabled, option):
            arguments.append(option)
    record = _cli.invoke(
        _TOOL,
        arguments,
        cwd=out,
        expected=[out / name],
        explicit=busco_path,
        timeout_seconds=timeout_seconds,
    )
    run = out / name
    summaries = sorted(run.glob("short_summary*.txt"))
    tables = sorted(run.rglob("full_table.tsv"))
    composition = sorted(run.rglob("scaffold_composition.txt"))
    summary = _parse_summary(summaries[0]) if summaries else {}
    table_rows: list[dict[str, Any]] = []
    if tables:
        lines = [
            line
            for line in tables[0].read_text(encoding="utf-8", errors="replace").splitlines()
            if line and not line.startswith("# Busco id")
        ]
        header = next(
            (
                line
                for line in tables[0].read_text(encoding="utf-8", errors="replace").splitlines()
                if line.startswith("# Busco id")
            ),
            None,
        )
        names = header.lstrip("# ").split("\t") if header else None
        table_rows = (
            _cli.read_table(tables[0], header=False, names=names, limit=20000) if names else []
        )
        table_rows = [
            row for row in table_rows if not str(next(iter(row.values()), "")).startswith("#")
        ]
        del lines
    return _cli.finish(
        record,
        "busco_completeness",
        run_directory=str(run),
        summary_file=str(summaries[0]) if summaries else None,
        full_table=table_rows,
        scaffold_composition=_cli.read_table(composition[0], limit=20000) if composition else [],
        **summary,
    )


__all__ = ["busco_assess_completeness"]
_register_public_verbose(globals())
