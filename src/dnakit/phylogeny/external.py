"""JSON-safe facades for caller-installed phylogenetic inference programs."""

from __future__ import annotations

from pathlib import Path

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import (
    ExecutionConfig,
    artifacts,
    input_path,
    run,
    run_to_file,
    workspace,
)
from dnakit.exceptions import ConfigurationError
from dnakit.tools.serialization import serialize_for_agent


def infer_ml_tree(
    alignment: str | Path,
    output_directory: str | Path,
    *,
    model: str = "MFP",
    bootstrap: int = 0,
    threads: int = 1,
    program: str = "iqtree2",
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Infer a maximum-likelihood tree with caller-installed IQ-TREE 2.

    ``model`` is passed as one model expression to IQ-TREE's ``-m`` option;
    ``bootstrap`` is either zero or the requested ultrafast-bootstrap replicate
    count.  The wrapper keeps the command bounded and returns the generated
    ``.treefile`` plus all files written below the fresh output directory.
    """

    alignment_path = input_path(alignment)
    if not isinstance(model, str) or not model.strip() or len(model) > 200:
        raise ConfigurationError("model must be a non-empty string of at most 200 characters.")
    if isinstance(bootstrap, bool) or not isinstance(bootstrap, int) or bootstrap < 0:
        raise ConfigurationError("bootstrap must be a non-negative integer.")
    if bootstrap > 100_000:
        raise ConfigurationError("bootstrap must be at most 100000.")
    config = ExecutionConfig(
        threads=threads,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    prefix = directory / "tree"
    treefile = Path(str(prefix) + ".treefile")
    arguments = [
        "-s",
        str(alignment_path),
        "-m",
        model.strip(),
        "-nt",
        str(config.threads),
        "-pre",
        str(prefix),
    ]
    if bootstrap:
        arguments.extend(("-B", str(bootstrap)))
    log = run(program, arguments, directory, config, outputs=(treefile,))
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "infer_ml_tree",
            "backend": "iqtree2",
            "result": {
                "directory": directory,
                "treefile": treefile,
                "artifacts": artifacts(directory),
                "model": model.strip(),
                "bootstrap": bootstrap,
                "log": log,
            },
            "provenance": {
                "implementation": "dnakit.phylogeny.external",
                "execution": "caller-installed executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )


def infer_fast_tree(
    alignment: str | Path,
    output_directory: str | Path,
    *,
    nucleotide: bool = True,
    gtr: bool = False,
    gamma: bool = False,
    program: str = "FastTree",
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Infer an approximate maximum-likelihood tree with caller-installed FastTree."""

    alignment_path = input_path(alignment)
    if not isinstance(nucleotide, bool) or not isinstance(gtr, bool) or not isinstance(gamma, bool):
        raise ConfigurationError("nucleotide, gtr and gamma must be booleans.")
    config = ExecutionConfig(
        threads=1,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    treefile = directory / "tree.newick"
    arguments: list[str] = []
    if nucleotide:
        arguments.append("-nt")
    if gtr:
        arguments.append("-gtr")
    if gamma:
        arguments.append("-gamma")
    arguments.append(str(alignment_path))
    log = run_to_file(program, arguments, treefile, config)
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "infer_fast_tree",
            "backend": "fasttree",
            "result": {
                "directory": directory,
                "treefile": treefile,
                "artifacts": artifacts(directory),
                "nucleotide": nucleotide,
                "gtr": gtr,
                "gamma": gamma,
                "log": log,
            },
            "provenance": {
                "implementation": "dnakit.phylogeny.external",
                "execution": "caller-installed executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )


__all__ = ["infer_fast_tree", "infer_ml_tree"]

_register_public_verbose(globals(), tuple(__all__))
