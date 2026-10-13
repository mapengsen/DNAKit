"""Bounded facade for caller-installed bisulfite methylation extraction."""

from __future__ import annotations

from pathlib import Path

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import ExecutionConfig, artifacts, input_path, run, workspace
from dnakit.exceptions import ConfigurationError
from dnakit.tools.serialization import serialize_for_agent


def extract_methylation_methyldackel(
    reference: str | Path,
    alignment: str | Path,
    output_directory: str | Path,
    *,
    include_chg: bool = True,
    include_chh: bool = True,
    program: str = "MethylDackel",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Extract CpG/CHG/CHH methylation bedGraphs with MethylDackel."""

    reference_path = input_path(reference)
    alignment_path = input_path(alignment)
    if not isinstance(include_chg, bool) or not isinstance(include_chh, bool):
        raise ConfigurationError("include_chg and include_chh must be booleans.")
    if not isinstance(program, str) or not program.strip() or any(char in program for char in "\x00\r\n"):
        raise ConfigurationError("program must be a non-empty executable name without newlines.")
    directory = workspace(output_directory)
    prefix = directory / "methylation"
    cpg = Path(str(prefix) + "_CpG.bedGraph")
    chg = Path(str(prefix) + "_CHG.bedGraph")
    chh = Path(str(prefix) + "_CHH.bedGraph")
    outputs = [cpg]
    arguments = ["extract"]
    if include_chg:
        arguments.append("--CHG")
        outputs.append(chg)
    if include_chh:
        arguments.append("--CHH")
        outputs.append(chh)
    arguments.extend((str(reference_path), str(alignment_path), str(prefix)))
    config = ExecutionConfig(
        threads=1,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    run(program.strip(), arguments, directory, config, outputs=tuple(outputs))
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "extract_methylation_methyldackel",
            "backend": "methyldackel",
            "result": {
                "directory": directory,
                "artifacts": artifacts(directory),
                "reference": reference_path,
                "alignment": alignment_path,
                "cpg": cpg,
                "chg": chg if include_chg else None,
                "chh": chh if include_chh else None,
                "include_chg": include_chg,
                "include_chh": include_chh,
            },
            "provenance": {
                "implementation": "dnakit.epigenome.external",
                "execution": "caller-installed executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )


__all__ = ["extract_methylation_methyldackel"]

_register_public_verbose(globals(), tuple(__all__))
