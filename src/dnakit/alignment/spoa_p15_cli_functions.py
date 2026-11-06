"""SPOA command-line facade: partial order alignment of sequences as an MSA, a graph in GFA or DOT, with strand picking.

SPOA (4.x, tested against 4.1.6) aligns sequences into a partial order graph. The user installs it (``conda install -c bioconda
spoa`` or builds https://github.com/rvaser/spoa) and puts ``spoa`` on PATH or passes ``spoa_path``. The function runs ``spoa``
with a fixed argument whitelist and no shell and returns a JSON-serialisable dictionary (status, backend, version, command,
output files, the parsed alignment, a bounded log tail). Output files must not exist yet.

``dnakit.alignment.external.consensus_spoa`` writes only the consensus; this function exposes the alignment and the graph.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p15_cli as _cli
from dnakit.exceptions import ConfigurationError

_TOOL = _cli.Tool(
    name="spoa",
    executable="spoa",
    hint=(
        "Install SPOA 4.x (for example `conda install -c bioconda spoa`, or build https://github.com/rvaser/spoa) and put spoa on "
        "PATH or pass spoa_path."
    ),
    tested_version="4.1.6",
    citation="https://doi.org/10.1101/gr.214270.116",
    version_arguments=("--version",),
    version_pattern=r"^(\d+\.\d+\.\d+)",
)
_RESULTS = {
    "consensus": "0",
    "msa": "1",
    "consensus_and_msa": "2",
    "graph": "3",
    "consensus_and_graph": "4",
}
_ALGORITHMS = {"local": 0, "global": 1, "semi_global": 2}


def _fasta_records(path: Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            records.append({"name": line[1:].strip(), "sequence": ""})
        elif records:
            records[-1]["sequence"] += line.strip()
    return records


def spoa_partial_order_alignment(
    sequences: str,
    output_path: str,
    *,
    result: str = "msa",
    algorithm: str = "local",
    match: int = 5,
    mismatch: int = -4,
    gap_open: int = -8,
    gap_extend: int = -6,
    gap_open_second: int = -10,
    gap_extend_second: int = -4,
    strand_ambiguous: bool = False,
    min_coverage: int | None = None,
    dot_path: str | None = None,
    spoa_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align sequences into a partial order graph and write the multiple sequence alignment or the graph (``spoa``).

    ``sequences`` is a FASTA or FASTQ file (optionally gzip-compressed). ``result`` selects what goes to the new file
    ``output_path``: ``'msa'`` (multiple sequence alignment, FASTA, ``-r 1``), ``'consensus_and_msa'`` (``-r 2``, the consensus
    is the last record), ``'consensus'`` (``-r 0``), ``'graph'`` (the partial order graph in GFA, ``-r 3``) or
    ``'consensus_and_graph'`` (GFA with the consensus as a path, ``-r 4``). ``dot_path`` (optional, new file) additionally
    writes the graph in Graphviz DOT format (``-d``). ``algorithm`` is ``'local'`` (Smith-Waterman, SPOA's default),
    ``'global'`` (Needleman-Wunsch) or ``'semi_global'``; ``match``, ``mismatch``, ``gap_open``, ``gap_extend`` and the second
    affine pair ``gap_open_second`` / ``gap_extend_second`` are the scores (SPOA's defaults 5, -4, -8, -6, -10, -4; the gap
    scores must not be positive; the gap model is linear, affine or convex as in SPOA). ``strand_ambiguous`` aligns every
    sequence in the orientation that scores better and writes it that way (``-s``). ``min_coverage`` (consensus only) is the
    minimal coverage of a consensus base. With a MSA result the alignment is returned as ``alignment`` records (``name``,
    ``sequence``) and, when present, ``consensus``.
    """
    source = _cli.existing_file(sequences, "sequences")
    target = _cli.new_file(output_path, "output_path")
    mode = _cli.choice(result, tuple(_RESULTS), "result")
    arguments = [
        "-l",
        str(_ALGORITHMS[_cli.choice(algorithm, tuple(_ALGORITHMS), "algorithm")]),
        "-r",
        _RESULTS[mode],
        "-m",
        str(_cli.integer(match, "match")),
        "-n",
        str(_cli.integer(mismatch, "mismatch")),
        "-g",
        str(_cli.integer(gap_open, "gap_open", maximum=0)),
        "-e",
        str(_cli.integer(gap_extend, "gap_extend", maximum=0)),
        "-q",
        str(_cli.integer(gap_open_second, "gap_open_second", maximum=0)),
        "-c",
        str(_cli.integer(gap_extend_second, "gap_extend_second", maximum=0)),
    ]
    if _cli.flag(strand_ambiguous, "strand_ambiguous"):
        arguments.append("-s")
    if min_coverage is not None:
        if mode != "consensus":
            raise ConfigurationError("min_coverage applies to result='consensus'.")
        arguments += ["--min-coverage", str(_cli.integer(min_coverage, "min_coverage", minimum=0))]
    expected: list[Path] = []
    graph_dot: Path | None = None
    if dot_path is not None:
        graph_dot = _cli.new_file(dot_path, "dot_path")
        arguments += ["-d", str(graph_dot)]
        expected.append(graph_dot)
    arguments.append(str(source))
    record = _cli.invoke(
        _TOOL,
        arguments,
        cwd=target.parent,
        stdout_to=target,
        expected=expected,
        explicit=spoa_path,
        timeout_seconds=timeout_seconds,
    )
    extra: dict[str, Any] = {
        "output": str(target),
        "result": mode,
        "dot_file": str(graph_dot) if graph_dot else None,
    }
    if mode in {"msa", "consensus_and_msa", "consensus"}:
        rows = _fasta_records(target)
        if mode != "msa" and rows:
            extra["consensus"] = (
                rows[-1]["sequence"] if mode == "consensus_and_msa" else rows[0]["sequence"]
            )
            rows = rows[:-1] if mode == "consensus_and_msa" else []
        extra["alignment"] = rows
        extra["alignment_length"] = len(rows[0]["sequence"]) if rows else None
    else:
        text = target.read_text(encoding="utf-8", errors="replace").splitlines()
        extra["node_count"] = sum(1 for line in text if line.startswith("S\t"))
        extra["edge_count"] = sum(1 for line in text if line.startswith("L\t"))
        extra["path_count"] = sum(1 for line in text if line.startswith("P\t"))
    return _cli.finish(record, "partial_order_alignment", **extra)


__all__ = ["spoa_partial_order_alignment"]
_register_public_verbose(globals())
