"""Barrnap messenger-RNA feature prediction as a typed DNAKit facade.

Barrnap 1.x (https://github.com/tseemann/barrnap, tested against 1.10.5) can annotate more than rRNA. With ``--mrna`` it
predicts protein-coding genes with Pyrodigal and reports, in GFF3, the ``mRNA`` (extended upstream over the ribosome
binding site), the ``RBS`` (ribosome binding site with its motif), the ``CDS`` (named from a DIAMOND search against the
Barrnap protein database for the chosen kingdom), the ``sig_pep`` (signal peptide, ZiggyPep) and the Rho-independent
``terminator`` (TransTermHP) features.

The user installs Barrnap with its dependencies (``conda install -c bioconda -c conda-forge barrnap``). The function runs
``barrnap`` with an explicit typed signature mapped to fixed flags; no shell is used, the program is resolved from an
explicit path, then PATH, and a missing program raises ``BackendUnavailableError`` with an install hint. Barrnap writes
GFF3 to standard output; the facade streams it to a new file, keeps standard error separate, and returns a
JSON-serialisable dictionary with the status, backend, version, command, the GFF3 path, feature counts per type, the
parsed features and a bounded log tail.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _common
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_HINT = (
    "Install Barrnap 1.10 or later with its dependencies (conda install -c bioconda -c conda-forge barrnap) and put "
    "barrnap on PATH or pass barrnap_path."
)
_KINGDOMS = ("bac", "arc", "fun")
_GFF_COLUMNS = ("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")


def _attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def _parse_gff(path: Path) -> list[dict[str, Any]]:
    """The feature lines of a GFF3 file (comments and any trailing FASTA section skipped)."""
    features: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            row: dict[str, Any] = dict(zip(_GFF_COLUMNS, fields, strict=True))
            row["start"] = int(row["start"])
            row["end"] = int(row["end"])
            row["attributes"] = _attributes(row["attributes"])
            features.append(row)
    return features


def _genetic_code(value: int | str) -> str:
    if value == "auto":
        return "auto"
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= 33:
        raise ConfigurationError("genetic_code must be 'auto' or a translation table number from 1 to 33.")
    return str(value)


def predict_mrna_features_barrnap(
    sequence_files: Sequence[str],
    output_gff: str,
    *,
    kingdom: str = "bac",
    genetic_code: int | str = "auto",
    evalue: float = 0.001,
    threads: int = 1,
    add_ids: bool = False,
    include_rrna: bool = False,
    database_dir: str | None = None,
    max_features: int = 1000,
    barrnap_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict protein-coding genes with RBS, signal peptide and terminator features in DNA with Barrnap (``--mrna``).

    ``sequence_files`` are one or more DNA files (FASTA, GenBank or FASTQ, optionally gzip-compressed); Barrnap merges
    them into one set of sequences. Barrnap skips the mRNA search with a warning when the combined input is shorter
    than about 19 kb, so use genomes or contigs, not single genes. ``kingdom`` selects the protein database used to
    name the CDS (``bac`` bacteria, ``arc`` archaea, ``fun`` fungi) and, with ``genetic_code='auto'``, the translation
    table (11, 25 and 1); give a table number to override it. ``evalue`` is the DIAMOND cut-off for naming CDS.
    ``add_ids`` adds ``ID=`` tags (``--addids``). Barrnap's default rRNA search is switched off (``--no-rrna``) unless
    ``include_rrna`` is True. ``database_dir`` overrides Barrnap's database folder (``--dbdir``).

    The GFF3 file is written to the new path ``output_gff``. Result keys: ``gff_file``, ``feature_counts`` (feature
    type to count, for example ``mRNA``, ``RBS``, ``CDS``, ``sig_pep``, ``terminator``), ``feature_count``,
    ``features`` (the first ``max_features`` GFF3 rows as dictionaries with ``seqid``, ``source``, ``type``, ``start``,
    ``end``, ``score``, ``strand``, ``phase`` and an ``attributes`` dictionary) and ``features_truncated``.
    """
    if isinstance(sequence_files, str) or not sequence_files:
        raise ConfigurationError("sequence_files must be a non-empty list of existing files.")
    inputs = [str(_cli.existing_file(path, "sequence_files")) for path in sequence_files]
    output = _cli.new_file(output_gff, "output_gff")
    chosen_kingdom = _cli.choice(kingdom, _KINGDOMS, "kingdom")
    code = _genetic_code(genetic_code)
    cutoff = _cli.number(evalue, "evalue", minimum=0.0)
    if not cutoff > 0:
        raise ConfigurationError("evalue must be greater than 0.")
    cores = _cli.integer(threads, "threads", minimum=1)
    limit = _cli.integer(max_features, "max_features", minimum=0)
    arguments: list[str] = ["--quiet", "--mrna", "--kingdom", chosen_kingdom, "--gcode", code]
    arguments += ["--evalue", format(cutoff, ".12g"), "--threads", str(cores)]
    arguments.append("--rrna" if include_rrna else "--no-rrna")
    if add_ids:
        arguments.append("--addids")
    if database_dir is not None:
        folder = Path(database_dir).expanduser().resolve()
        if not folder.is_dir():
            raise ConfigurationError(f"database_dir must be an existing folder: {folder}")
        arguments += ["--dbdir", str(folder)]
    executable = _cli.program("barrnap", explicit=barrnap_path, hint=_HINT)
    version = _cli._version_of(executable)
    arguments += inputs
    stderr_file = output.with_name(output.name + ".stderr")
    if stderr_file.exists():
        raise ConfigurationError(f"A diagnostic file already exists next to output_gff: {stderr_file}")
    config = _common.ExecutionConfig(threads=cores, timeout_seconds=_cli.DEFAULT_TIMEOUT_SECONDS)
    return_code = 0
    try:
        log = _common.run_to_file(executable, arguments, output, config)
    except BackendExecutionError as error:
        log = str(error.context.get("log", "")) if error.context else ""
        return_code = 1
    finally:
        if stderr_file.exists():
            stderr_file.unlink()
    produced = output.is_file() and return_code == 0
    result: dict[str, Any] = {
        "status": "completed" if produced else "failed",
        "backend": "barrnap",
        "version": version,
        "command": ["barrnap", *arguments],
        "return_code": return_code,
        "outputs": [str(output)] if output.is_file() else [],
        "log_tail": log[-_cli.LOG_TAIL_CHARS :],
        "gff_file": str(output) if output.is_file() else None,
    }
    features = _parse_gff(output) if produced else []
    result["feature_counts"] = dict(Counter(feature["type"] for feature in features))
    result["feature_count"] = len(features)
    result["features"] = features[:limit]
    result["features_truncated"] = len(features) > limit
    return result


__all__ = ["predict_mrna_features_barrnap"]
_register_public_verbose(globals())
