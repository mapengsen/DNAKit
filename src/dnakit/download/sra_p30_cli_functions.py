"""SRA Toolkit download and FASTA/FASTQ extraction (fixed commands, no shell).

These operations may download large sequencing runs. Ordinary tests never use the
network. Callers install SRA Toolkit; files are written to a fresh output directory.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as cli
from dnakit.exceptions import ConfigurationError, DNAKitError

_HINT = "Install SRA Toolkit: conda install -c bioconda sra-tools."
_ACCESSION = re.compile(r"(?:SRR|ERR|DRR)[0-9]+\Z")


def _accession(value: str) -> str:
    if not isinstance(value, str) or not _ACCESSION.fullmatch(value):
        raise ConfigurationError("accession must be an SRR, ERR or DRR run accession.")
    return value


def prefetch_sra_run(
    accession: str,
    output_dir: str | Path,
    *,
    executable_path: str | None = None,
    timeout_seconds: float = 3600,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download an SRA run archive with NCBI prefetch; network access is required."""
    run = _accession(accession)
    executable = cli.program("prefetch", explicit=executable_path, hint=_HINT)
    output = cli.new_directory(output_dir, "output_dir")
    archive = output / run / f"{run}.sra"
    result = cli.run(
        "prefetch",
        ["-O", str(output), run],
        cwd=output,
        expected=[archive],
        explicit=executable,
        hint=_HINT,
        timeout_seconds=timeout_seconds,
    )
    result.update(accession=run, archive=str(archive), network_operation=True)
    return result


def extract_sra_reads(
    accession: str,
    output_dir: str | Path,
    *,
    output_format: str = "fastq",
    local_sra: str | Path | None = None,
    threads: int | None = None,
    executable_path: str | None = None,
    timeout_seconds: float = 3600,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract FASTQ or FASTA with fasterq-dump; optionally use an existing local SRA archive.

    ``output_format`` is ``fastq`` (default) or ``fasta``. Without ``local_sra`` the
    tool retrieves the run from NCBI. SRA Toolkit's default read splitting is used;
    all single and paired output files are returned. No existing output is replaced.
    """
    run = _accession(accession)
    kind = cli.choice(output_format, ("fastq", "fasta"), "output_format")
    source = run if local_sra is None else str(cli.existing_file(local_sra, "local_sra"))
    count = None if threads is None else cli.integer(threads, "threads", minimum=1, maximum=256)
    executable = cli.program("fasterq-dump", explicit=executable_path, hint=_HINT)
    output = cli.new_directory(output_dir, "output_dir")
    args = ["-o", str(output / f"{run}.{kind}")]
    if kind == "fasta":
        args.append("--fasta")
    if count is not None:
        args += ["--threads", str(count)]
    args.append(source)
    result = cli.run(
        "fasterq-dump",
        args,
        cwd=output,
        expected=[],
        explicit=executable,
        hint=_HINT,
        timeout_seconds=timeout_seconds,
    )
    files = sorted(f for f in output.glob(f"{run}*.{kind}") if f.is_file())
    if not files or any(f.stat().st_size == 0 for f in files):
        raise DNAKitError("fasterq-dump finished without non-empty sequence output files.")
    result.update(
        accession=run,
        output_format=kind,
        outputs=[str(f) for f in files],
        network_operation=local_sra is None,
    )
    return result


__all__ = ["extract_sra_reads", "prefetch_sra_run"]
_register_public_verbose(globals(), tuple(__all__))
