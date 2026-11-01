"""VSEARCH command-line facades (pack P09): fixed argument lists, no shell, bounded execution.

Each public function runs one VSEARCH command (2.32.0 is the tested version) through the shared facade layer
``dnakit.bioinformatics._gap_cli``. The user installs VSEARCH (``conda install -c bioconda vsearch``, ``apt install
vsearch`` or a release binary) and either puts it on PATH or passes ``vsearch_path``. Every function writes into a fresh
``output_directory`` (never overwritten) under documented file names and returns a JSON-serialisable dictionary with the
status, version, command, the output files, small parsed results and a bounded log tail. Optional parameters default to
``None`` and then are not passed, so VSEARCH's own defaults apply. VSEARCH validates the options of every command, so a
parameter that the selected command does not accept makes the call fail with ``BackendExecutionError``.

Commands that are not DNA operations are not wrapped; the VSEARCH library API (``vsearch_api``) is a C++ interface of the
same commands and is covered by the corresponding function here.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Literal

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_PROGRAM = "vsearch"
_TESTED_VERSION = "2.32.0"
_CITATION = "https://doi.org/10.7717/peerj.2584"
_HINT = (
    "Install VSEARCH 2.x (for example `conda install -c bioconda vsearch` or `apt install vsearch`, or a release "
    "binary from https://github.com/torognes/vsearch/releases), then put it on PATH or pass vsearch_path."
)
_MAX_REPORT_CHARS = 200_000
_BLAST6 = (
    "query",
    "target",
    "identity",
    "alignment_length",
    "mismatches",
    "gap_opens",
    "query_start",
    "query_end",
    "target_start",
    "target_end",
    "evalue",
    "bit_score",
)
_UCHIME = (
    "score",
    "query",
    "parent_a",
    "parent_b",
    "closest_parent",
    "identity_model_query",
    "identity_query_a",
    "identity_query_b",
    "identity_a_b",
    "identity_query_top",
    "left_yes",
    "left_no",
    "left_abstain",
    "right_yes",
    "right_no",
    "right_abstain",
    "divergence",
    "chimeric",
)
_CUT_PATTERN = re.compile(r"^[ACGTUMRWSYKVHDBNacgtumrwsykvhdbn^_]+$")
_ROW_LIMIT = 1000


# ---------------------------------------------------------------------------
# Validation and argument helpers
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


def _opt_num(
    arguments: list[str],
    flag: str,
    value: float | None,
    name: str,
    *,
    minimum: float | None = None,
    maximum: float | None = None,
) -> None:
    if value is not None:
        arguments.extend(
            (flag, format(_cli.number(value, name, minimum=minimum, maximum=maximum), ".12g"))
        )


def _opt_text(arguments: list[str], flag: str, value: str | None, name: str) -> None:
    if value is not None:
        arguments.extend((flag, _cli.text(value, name)))


def _opt_choice(
    arguments: list[str], flag: str, value: str | None, options: Sequence[str], name: str
) -> None:
    if value is not None:
        arguments.extend((flag, _cli.choice(value, options, name)))


def _quality_range(
    arguments: list[str], fastq_ascii: int | None, fastq_qmin: int | None, fastq_qmax: int | None
) -> None:
    _opt_int(arguments, "--fastq_ascii", fastq_ascii, "fastq_ascii", minimum=33, maximum=64)
    _opt_int(arguments, "--fastq_qmin", fastq_qmin, "fastq_qmin", minimum=0, maximum=93)
    _opt_int(arguments, "--fastq_qmax", fastq_qmax, "fastq_qmax", minimum=0, maximum=93)


def _relabel(
    arguments: list[str],
    relabel: str | None,
    relabel_mode: str | None = None,
    *,
    name: str = "relabel",
) -> None:
    """``--relabel PREFIX`` or one of the digest modes (``md5``, ``sha1``, ``self``, ``keep`` adds to a prefix)."""
    if relabel is not None:
        arguments.extend(("--relabel", _cli.text(relabel, name)))
    if relabel_mode is not None:
        mode = _cli.choice(relabel_mode, ("keep", "md5", "sha1", "self"), "relabel_mode")
        if mode == "keep" and relabel is None:
            raise ConfigurationError("relabel_mode='keep' needs a relabel prefix.")
        arguments.append(f"--relabel_{mode}")


def _sequence_input(value: str | Path, name: str) -> str:
    return str(_cli.existing_file(value, name))


def _threads(arguments: list[str], threads: int | None) -> None:
    _opt_int(arguments, "--threads", threads, "threads", minimum=0, maximum=1024)


def _begin(vsearch_path: str | None, output_directory: str | Path) -> Path:
    """Resolve VSEARCH first (so a missing program creates nothing), then claim the fresh output directory."""
    _cli.program(_PROGRAM, explicit=vsearch_path, hint=_HINT)
    return _cli.new_directory(output_directory, "output_directory")


def _execute(
    operation: str,
    arguments: Sequence[str],
    directory: Path,
    outputs: dict[str, Path],
    *,
    vsearch_path: str | None,
    timeout_seconds: float,
    **extra: object,
) -> dict[str, Any]:
    timeout = _cli.number(timeout_seconds, "timeout_seconds", minimum=1.0, maximum=86_400.0)
    run = _cli.run(
        _PROGRAM,
        ["--quiet", *arguments],
        cwd=directory,
        expected=list(outputs.values()),
        explicit=vsearch_path,
        hint=_HINT,
        timeout_seconds=timeout,
    )
    if run["status"] != "completed":
        raise BackendExecutionError(
            f"vsearch {operation} failed (exit code {run['return_code']}).",
            code="VSEARCH_EXECUTION_FAILED",
            context={"operation": operation, "log_tail": run["log_tail"]},
        )
    result: dict[str, Any] = {
        "status": "completed",
        "backend": _PROGRAM,
        "operation": operation,
        "version": run["version"],
        "command": run["command"],
        "return_code": run["return_code"],
        "output_directory": str(directory),
        "outputs": {name: str(path) for name, path in outputs.items()},
        "log_tail": run["log_tail"],
        "provenance": {"tool": "VSEARCH", "tested_version": _TESTED_VERSION, "citation": _CITATION},
    }
    result.update(extra)
    return result


def _count(path: Path, kind: str) -> int:
    """Number of FASTA (``>``) or FASTQ (four lines per record) records of a file."""
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        if kind == "fasta":
            return sum(1 for line in handle if line.startswith(">"))
        return sum(1 for _ in handle) // 4


def _rows(
    path: Path, names: Sequence[str], limit: int = _ROW_LIMIT
) -> tuple[list[dict[str, Any]], int]:
    """Tab-separated rows as dictionaries (numbers converted); also the total row count."""
    rows: list[dict[str, Any]] = []
    total = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields == [""]:
                continue
            total += 1
            if len(rows) < limit:
                rows.append(
                    {name: _convert(value) for name, value in zip(names, fields, strict=False)}
                )
    return rows, total


def _convert(value: str) -> Any:
    for kind in (int, float):
        try:
            return kind(value)
        except ValueError:
            continue
    return value


def _fastx_kind(path: str) -> str:
    """``fastq`` when the first non-empty line starts with ``@``, otherwise ``fasta``."""
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.strip():
                return "fastq" if line.startswith("@") else "fasta"
    return "fasta"


def _output_format(value: str) -> str:
    return _cli.choice(value, ("fasta", "fastq"), "output_format")


# ---------------------------------------------------------------------------
# Alignment, search and clustering
# ---------------------------------------------------------------------------


def vsearch_allpairs_global(
    input_fasta: str | Path,
    output_directory: str | Path,
    *,
    accept_all: bool = True,
    min_identity: float | None = None,
    identity_definition: int | None = None,
    match_score: int | None = None,
    mismatch_score: int | None = None,
    gap_open: str | None = None,
    gap_extend: str | None = None,
    threads: int | None = None,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Global alignment of every pair of sequences of a FASTA file (``vsearch --allpairs_global``).

    Writes ``alignments.aln`` (human-readable alignments) and ``hits.tsv`` (BLAST-like 12-column table, parsed into
    ``records``). ``accept_all=True`` reports all pairs (``--acceptall``); with False only pairs that pass ``min_identity``
    (``--id``, 0-1) are reported. ``identity_definition`` is ``--iddef`` (0-4); ``gap_open``/``gap_extend`` are VSEARCH
    penalty strings such as ``"20I/2E"``.
    """
    source = _sequence_input(input_fasta, "input_fasta")
    directory = _begin(vsearch_path, output_directory)
    outputs = {"alignments": directory / "alignments.aln", "hits": directory / "hits.tsv"}
    arguments = [
        "--allpairs_global",
        source,
        "--alnout",
        str(outputs["alignments"]),
        "--blast6out",
        str(outputs["hits"]),
    ]
    _switch(arguments, "--acceptall", accept_all, "accept_all")
    _opt_num(arguments, "--id", min_identity, "min_identity", minimum=0.0, maximum=1.0)
    _opt_int(arguments, "--iddef", identity_definition, "identity_definition", minimum=0, maximum=4)
    _opt_int(arguments, "--match", match_score, "match_score")
    _opt_int(arguments, "--mismatch", mismatch_score, "mismatch_score")
    _opt_text(arguments, "--gapopen", gap_open, "gap_open")
    _opt_text(arguments, "--gapext", gap_extend, "gap_extend")
    _threads(arguments, threads)
    result = _execute(
        "allpairs_global",
        arguments,
        directory,
        outputs,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )
    result["records"], result["record_count"] = _rows(outputs["hits"], _BLAST6)
    return result


def vsearch_search_global(
    query_fasta: str | Path,
    database: str | Path,
    output_directory: str | Path,
    *,
    min_identity: float,
    strand: Literal["plus", "both"] | None = None,
    identity_definition: int | None = None,
    max_accepts: int | None = None,
    max_rejects: int | None = None,
    min_word_matches: int | None = None,
    word_length: int | None = None,
    query_coverage: float | None = None,
    target_coverage: float | None = None,
    top_hits_only: bool = False,
    output_no_hits: bool = False,
    query_mask: Literal["none", "dust", "soft"] | None = None,
    database_mask: Literal["none", "dust", "soft"] | None = None,
    abundance_in: bool = False,
    write_alignments: bool = False,
    write_uc: bool = False,
    write_otu_table: bool = False,
    write_sam: bool = False,
    write_matched_queries: bool = False,
    threads: int | None = None,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Global-alignment search of query sequences against a FASTA or UDB database (``vsearch --usearch_global``).

    ``min_identity`` (``--id``, 0-1) is required by VSEARCH. Always writes ``hits.tsv`` (BLAST-like table, parsed into
    ``records``); the flags ``write_alignments``, ``write_uc``, ``write_otu_table``, ``write_sam`` and
    ``write_matched_queries`` (matched and unmatched queries as FASTA) add ``alignments.aln``, ``hits.uc``,
    ``otu_table.tsv``, ``hits.sam``, ``matched.fasta`` and ``unmatched.fasta``. The VSEARCH library calls for single and
    batch search are this operation.
    """
    query = _sequence_input(query_fasta, "query_fasta")
    db = _sequence_input(database, "database")
    directory = _begin(vsearch_path, output_directory)
    outputs = {"hits": directory / "hits.tsv"}
    arguments = ["--usearch_global", query, "--db", db, "--blast6out", str(outputs["hits"])]
    _opt_num(arguments, "--id", min_identity, "min_identity", minimum=0.0, maximum=1.0)
    _opt_choice(arguments, "--strand", strand, ("plus", "both"), "strand")
    _opt_int(arguments, "--iddef", identity_definition, "identity_definition", minimum=0, maximum=4)
    _opt_int(arguments, "--maxaccepts", max_accepts, "max_accepts", minimum=0)
    _opt_int(arguments, "--maxrejects", max_rejects, "max_rejects", minimum=0)
    _opt_int(arguments, "--minwordmatches", min_word_matches, "min_word_matches", minimum=0)
    _opt_int(arguments, "--wordlength", word_length, "word_length", minimum=3, maximum=15)
    _opt_num(arguments, "--query_cov", query_coverage, "query_coverage", minimum=0.0, maximum=1.0)
    _opt_num(
        arguments, "--target_cov", target_coverage, "target_coverage", minimum=0.0, maximum=1.0
    )
    _switch(arguments, "--top_hits_only", top_hits_only, "top_hits_only")
    _switch(arguments, "--output_no_hits", output_no_hits, "output_no_hits")
    _opt_choice(arguments, "--qmask", query_mask, ("none", "dust", "soft"), "query_mask")
    _opt_choice(arguments, "--dbmask", database_mask, ("none", "dust", "soft"), "database_mask")
    _switch(arguments, "--sizein", abundance_in, "abundance_in")
    for enabled, name, flag, filename in (
        (write_alignments, "alignments", "--alnout", "alignments.aln"),
        (write_uc, "uc", "--uc", "hits.uc"),
        (write_otu_table, "otu_table", "--otutabout", "otu_table.tsv"),
        (write_sam, "sam", "--samout", "hits.sam"),
    ):
        if _flag(enabled, f"write_{name}"):
            outputs[name] = directory / filename
            arguments.extend((flag, str(outputs[name])))
    if _flag(write_matched_queries, "write_matched_queries"):
        outputs["matched"] = directory / "matched.fasta"
        outputs["unmatched"] = directory / "unmatched.fasta"
        arguments.extend(
            ("--matched", str(outputs["matched"]), "--notmatched", str(outputs["unmatched"]))
        )
    _threads(arguments, threads)
    result = _execute(
        "search_global",
        arguments,
        directory,
        outputs,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )
    result["records"], result["record_count"] = _rows(outputs["hits"], _BLAST6)
    return result


def vsearch_cluster_unoise(
    input_fasta: str | Path,
    output_directory: str | Path,
    *,
    min_abundance: int | None = None,
    unoise_alpha: float | None = None,
    identity: float | None = None,
    abundance_in: bool = False,
    abundance_out: bool = False,
    relabel: str | None = None,
    write_otu_table: bool = False,
    threads: int | None = None,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """UNOISE3 denoising of amplicon reads (``vsearch --cluster_unoise``).

    ``min_abundance`` is ``--minsize`` (VSEARCH default 8), ``unoise_alpha`` is ``--unoise_alpha`` (default 2.0). Writes
    ``centroids.fasta`` (the denoised sequences) and ``clusters.uc``; ``write_otu_table`` adds ``otu_table.tsv``.
    ``abundance_out`` writes ``;size=`` annotations; ``relabel`` renames centroids with a prefix.
    """
    source = _sequence_input(input_fasta, "input_fasta")
    directory = _begin(vsearch_path, output_directory)
    outputs = {"centroids": directory / "centroids.fasta", "clusters": directory / "clusters.uc"}
    arguments = [
        "--cluster_unoise",
        source,
        "--centroids",
        str(outputs["centroids"]),
        "--uc",
        str(outputs["clusters"]),
    ]
    _opt_int(arguments, "--minsize", min_abundance, "min_abundance", minimum=1)
    _opt_num(arguments, "--unoise_alpha", unoise_alpha, "unoise_alpha", minimum=0.0)
    _opt_num(arguments, "--id", identity, "identity", minimum=0.0, maximum=1.0)
    _switch(arguments, "--sizein", abundance_in, "abundance_in")
    _switch(arguments, "--sizeout", abundance_out, "abundance_out")
    _opt_text(arguments, "--relabel", relabel, "relabel")
    if _flag(write_otu_table, "write_otu_table"):
        outputs["otu_table"] = directory / "otu_table.tsv"
        arguments.extend(("--otutabout", str(outputs["otu_table"])))
    _threads(arguments, threads)
    result = _execute(
        "cluster_unoise",
        arguments,
        directory,
        outputs,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )
    result["centroid_count"] = _count(outputs["centroids"], "fasta")
    return result


def vsearch_detect_chimeras(
    input_fasta: str | Path,
    output_directory: str | Path,
    *,
    method: Literal[
        "uchime_denovo", "uchime2_denovo", "uchime3_denovo", "uchime_ref", "chimeras_denovo"
    ] = "uchime_denovo",
    database: str | Path | None = None,
    abundance_skew: float | None = None,
    abundance_in: bool = False,
    abundance_out: bool = False,
    min_score: float | None = None,
    min_divergence: float | None = None,
    min_diffs: int | None = None,
    no_vote_pseudocount: float | None = None,
    no_vote_weight: float | None = None,
    exclude_same_label: bool = False,
    exclude_same_sequence: bool = False,
    chimeras_diff_pct: float | None = None,
    chimeras_length_min: int | None = None,
    chimeras_parents_max: int | None = None,
    chimeras_parts: int | None = None,
    relabel: str | None = None,
    threads: int | None = None,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Chimera detection (``vsearch --uchime_denovo`` and its relatives).

    ``method`` selects the algorithm: ``uchime_denovo`` (UCHIME, abundance-ordered), ``uchime2_denovo`` and
    ``uchime3_denovo`` (denoised amplicons), ``uchime_ref`` (needs ``database``) and ``chimeras_denovo`` (the newer
    algorithm for long exact sequences, with ``chimeras_diff_pct``, ``chimeras_length_min``, ``chimeras_parents_max``
    and ``chimeras_parts``). ``abundance_skew`` is ``--abskew`` (VSEARCH defaults: 2.0, 16.0 for uchime3, 1.0 for
    chimeras_denovo). ``min_score``, ``min_divergence``, ``min_diffs``, ``no_vote_pseudocount`` and ``no_vote_weight`` are
    ``--minh``, ``--mindiv``, ``--mindiffs``, ``--dn`` and ``--xn`` of the UCHIME methods. Writes ``chimeras.fasta``,
    ``nonchimeras.fasta``, ``report.tsv`` (18-column UCHIME table, parsed into ``records``) and ``alignments.txt``; the
    UCHIME methods also write ``borderline.fasta``.
    """
    methods = ("uchime_denovo", "uchime2_denovo", "uchime3_denovo", "uchime_ref", "chimeras_denovo")
    _cli.choice(method, methods, "method")
    source = _sequence_input(input_fasta, "input_fasta")
    if method == "uchime_ref":
        if database is None:
            raise ConfigurationError("method='uchime_ref' needs a reference database.")
    elif database is not None or exclude_same_label or exclude_same_sequence:
        raise ConfigurationError(
            "database, exclude_same_label and exclude_same_sequence belong to method='uchime_ref'."
        )
    denovo_only = (chimeras_diff_pct, chimeras_length_min, chimeras_parents_max, chimeras_parts)
    uchime_only = (min_score, min_divergence, min_diffs, no_vote_pseudocount, no_vote_weight)
    if method != "chimeras_denovo" and any(value is not None for value in denovo_only):
        raise ConfigurationError("chimeras_* parameters belong to method='chimeras_denovo'.")
    if method == "chimeras_denovo" and any(value is not None for value in uchime_only):
        raise ConfigurationError(
            "min_score, min_divergence, min_diffs and the vote parameters belong to the UCHIME methods."
        )
    db = None if database is None else _sequence_input(database, "database")
    directory = _begin(vsearch_path, output_directory)
    outputs = {
        "chimeras": directory / "chimeras.fasta",
        "nonchimeras": directory / "nonchimeras.fasta",
        "report": directory / "report.tsv",
        "alignments": directory / "alignments.txt",
    }
    arguments = [
        f"--{method}",
        source,
        "--chimeras",
        str(outputs["chimeras"]),
        "--nonchimeras",
        str(outputs["nonchimeras"]),
    ]
    if method == "chimeras_denovo":
        arguments.extend(
            ("--tabbedout", str(outputs["report"]), "--alnout", str(outputs["alignments"]))
        )
    else:
        outputs["borderline"] = directory / "borderline.fasta"
        arguments.extend(
            (
                "--uchimeout",
                str(outputs["report"]),
                "--uchimealns",
                str(outputs["alignments"]),
                "--borderline",
                str(outputs["borderline"]),
            )
        )
    if db is not None:
        arguments.extend(("--db", db))
    _opt_num(arguments, "--abskew", abundance_skew, "abundance_skew", minimum=0.0)
    _switch(arguments, "--sizein", abundance_in, "abundance_in")
    _switch(arguments, "--sizeout", abundance_out, "abundance_out")
    _opt_num(arguments, "--minh", min_score, "min_score", minimum=0.0)
    _opt_num(arguments, "--mindiv", min_divergence, "min_divergence", minimum=0.0)
    _opt_int(arguments, "--mindiffs", min_diffs, "min_diffs", minimum=1)
    _opt_num(arguments, "--dn", no_vote_pseudocount, "no_vote_pseudocount", minimum=0.0)
    _opt_num(arguments, "--xn", no_vote_weight, "no_vote_weight", minimum=0.0)
    _switch(arguments, "--self", exclude_same_label, "exclude_same_label")
    _switch(arguments, "--selfid", exclude_same_sequence, "exclude_same_sequence")
    _opt_num(
        arguments,
        "--chimeras_diff_pct",
        chimeras_diff_pct,
        "chimeras_diff_pct",
        minimum=0.0,
        maximum=100.0,
    )
    _opt_int(
        arguments, "--chimeras_length_min", chimeras_length_min, "chimeras_length_min", minimum=1
    )
    _opt_int(
        arguments, "--chimeras_parents_max", chimeras_parents_max, "chimeras_parents_max", minimum=2
    )
    _opt_int(arguments, "--chimeras_parts", chimeras_parts, "chimeras_parts", minimum=2)
    _opt_text(arguments, "--relabel", relabel, "relabel")
    _threads(arguments, threads if method == "uchime_ref" else None)
    result = _execute(
        "detect_chimeras",
        arguments,
        directory,
        outputs,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
        method=method,
    )
    result["records"], result["record_count"] = _rows(outputs["report"], _UCHIME)
    result["chimera_count"] = _count(outputs["chimeras"], "fasta")
    result["nonchimera_count"] = _count(outputs["nonchimeras"], "fasta")
    return result


# ---------------------------------------------------------------------------
# Dereplication, rereplication and sorting
# ---------------------------------------------------------------------------


def _dereplicate(
    operation: str,
    command: str,
    input_fasta: str | Path,
    output_directory: str | Path,
    *,
    fastq_output: bool = False,
    strand: str | None,
    abundance_in: bool,
    abundance_out: bool,
    min_unique_size: int | None,
    max_unique_size: int | None,
    top_n: int | None,
    relabel: str | None,
    relabel_mode: str | None,
    write_uc: bool = True,
    threads: int | None = None,
    vsearch_path: str | None,
    timeout_seconds: float,
) -> dict[str, Any]:
    source = _sequence_input(input_fasta, "input_fasta")
    directory = _begin(vsearch_path, output_directory)
    if fastq_output:
        outputs = {"unique": directory / "unique.fastq"}
        output_flag = "--fastqout"
    else:
        outputs = {"unique": directory / "unique.fasta"}
        output_flag = "--fastaout" if command in {"derep_smallmem", "fastx_uniques"} else "--output"
    arguments = [f"--{command}", source, output_flag, str(outputs["unique"])]
    if _flag(write_uc, "write_uc") and command != "derep_smallmem":
        outputs["clusters"] = directory / "clusters.uc"
        arguments.extend(("--uc", str(outputs["clusters"])))
    _opt_choice(arguments, "--strand", strand, ("plus", "both"), "strand")
    _switch(arguments, "--sizein", abundance_in, "abundance_in")
    _switch(arguments, "--sizeout", abundance_out, "abundance_out")
    _opt_int(arguments, "--minuniquesize", min_unique_size, "min_unique_size", minimum=1)
    _opt_int(arguments, "--maxuniquesize", max_unique_size, "max_unique_size", minimum=1)
    _opt_int(arguments, "--topn", top_n, "top_n", minimum=1)
    _relabel(arguments, relabel, relabel_mode)
    _threads(arguments, threads)
    result = _execute(
        operation,
        arguments,
        directory,
        outputs,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )
    result["unique_count"] = _count(outputs["unique"], "fastq" if fastq_output else "fasta")
    return result


def vsearch_dereplicate_fulllength(
    input_fasta: str | Path,
    output_directory: str | Path,
    *,
    method: Literal["fulllength", "smallmem", "fastx_uniques"] = "fulllength",
    fastq_output: bool = False,
    strand: Literal["plus", "both"] | None = None,
    abundance_in: bool = False,
    abundance_out: bool = False,
    min_unique_size: int | None = None,
    max_unique_size: int | None = None,
    top_n: int | None = None,
    relabel: str | None = None,
    relabel_mode: Literal["keep", "md5", "sha1", "self"] | None = None,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Full-length dereplication: collapse identical sequences and count them (``vsearch --derep_fulllength``).

    ``method`` selects the command: ``fulllength`` (``--derep_fulllength``), ``smallmem`` (``--derep_smallmem``, less memory,
    input order kept, no UC file) or ``fastx_uniques`` (``--fastx_uniques``, FASTA or FASTQ input; ``fastq_output=True``
    keeps the qualities). Writes ``unique.fasta`` (or ``unique.fastq``) and, except for ``smallmem``, ``clusters.uc``.
    ``abundance_out`` writes ``;size=`` annotations; ``strand='both'`` also merges reverse complements; ``relabel_mode``
    is one of keep, md5, sha1, self.
    """
    _cli.choice(method, ("fulllength", "smallmem", "fastx_uniques"), "method")
    if fastq_output and method != "fastx_uniques":
        raise ConfigurationError("fastq_output needs method='fastx_uniques'.")
    command = {
        "fulllength": "derep_fulllength",
        "smallmem": "derep_smallmem",
        "fastx_uniques": "fastx_uniques",
    }[method]
    return _dereplicate(
        "dereplicate_fulllength",
        command,
        input_fasta,
        output_directory,
        fastq_output=fastq_output,
        strand=strand,
        abundance_in=abundance_in,
        abundance_out=abundance_out,
        min_unique_size=min_unique_size,
        max_unique_size=max_unique_size,
        top_n=top_n,
        relabel=relabel,
        relabel_mode=relabel_mode,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )


def vsearch_dereplicate_by_id(
    input_fasta: str | Path,
    output_directory: str | Path,
    *,
    strand: Literal["plus", "both"] | None = None,
    abundance_in: bool = False,
    abundance_out: bool = False,
    min_unique_size: int | None = None,
    max_unique_size: int | None = None,
    top_n: int | None = None,
    relabel: str | None = None,
    relabel_mode: Literal["keep", "md5", "sha1", "self"] | None = None,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Dereplicate on both the identifier and the sequence (``vsearch --derep_id``).

    Only records with the same label AND the same sequence are collapsed. Writes ``unique.fasta`` and ``clusters.uc``.
    """
    return _dereplicate(
        "dereplicate_by_id",
        "derep_id",
        input_fasta,
        output_directory,
        strand=strand,
        abundance_in=abundance_in,
        abundance_out=abundance_out,
        min_unique_size=min_unique_size,
        max_unique_size=max_unique_size,
        top_n=top_n,
        relabel=relabel,
        relabel_mode=relabel_mode,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )


def vsearch_dereplicate_prefix(
    input_fasta: str | Path,
    output_directory: str | Path,
    *,
    strand: Literal["plus", "both"] | None = None,
    abundance_in: bool = False,
    abundance_out: bool = False,
    min_unique_size: int | None = None,
    max_unique_size: int | None = None,
    top_n: int | None = None,
    relabel: str | None = None,
    relabel_mode: Literal["keep", "md5", "sha1", "self"] | None = None,
    threads: int | None = None,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Dereplicate sequences that are prefixes of longer sequences (``vsearch --derep_prefix``).

    A sequence identical to the start of a longer one is merged into the longer one. Writes ``unique.fasta`` and
    ``clusters.uc``.
    """
    return _dereplicate(
        "dereplicate_prefix",
        "derep_prefix",
        input_fasta,
        output_directory,
        strand=strand,
        abundance_in=abundance_in,
        abundance_out=abundance_out,
        min_unique_size=min_unique_size,
        max_unique_size=max_unique_size,
        top_n=top_n,
        relabel=relabel,
        relabel_mode=relabel_mode,
        threads=threads,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )


def vsearch_rereplicate(
    input_fasta: str | Path,
    output_directory: str | Path,
    *,
    strip_abundance: bool = False,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Expand each record by its ``;size=`` abundance annotation (``vsearch --rereplicate``).

    Writes ``rereplicated.fasta``; ``strip_abundance`` removes the size annotation from the labels (``--xsize``).
    """
    source = _sequence_input(input_fasta, "input_fasta")
    directory = _begin(vsearch_path, output_directory)
    outputs = {"rereplicated": directory / "rereplicated.fasta"}
    arguments = ["--rereplicate", source, "--output", str(outputs["rereplicated"])]
    _switch(arguments, "--xsize", strip_abundance, "strip_abundance")
    result = _execute(
        "rereplicate",
        arguments,
        directory,
        outputs,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )
    result["record_count"] = _count(outputs["rereplicated"], "fasta")
    return result


def vsearch_sort_by_size(
    input_fasta: str | Path,
    output_directory: str | Path,
    *,
    min_size: int | None = None,
    max_size: int | None = None,
    top_n: int | None = None,
    abundance_out: bool = False,
    relabel: str | None = None,
    relabel_mode: Literal["keep", "md5", "sha1", "self"] | None = None,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sort FASTA records by decreasing ``;size=`` abundance (``vsearch --sortbysize``).

    ``min_size``/``max_size`` drop records outside the abundance range, ``top_n`` keeps the first n. Writes
    ``sorted.fasta``.
    """
    source = _sequence_input(input_fasta, "input_fasta")
    directory = _begin(vsearch_path, output_directory)
    outputs = {"sorted": directory / "sorted.fasta"}
    arguments = ["--sortbysize", source, "--output", str(outputs["sorted"])]
    _opt_int(arguments, "--minsize", min_size, "min_size", minimum=1)
    _opt_int(arguments, "--maxsize", max_size, "max_size", minimum=1)
    _opt_int(arguments, "--topn", top_n, "top_n", minimum=1)
    _switch(arguments, "--sizeout", abundance_out, "abundance_out")
    _relabel(arguments, relabel, relabel_mode)
    result = _execute(
        "sort_by_size",
        arguments,
        directory,
        outputs,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )
    result["record_count"] = _count(outputs["sorted"], "fasta")
    return result


# ---------------------------------------------------------------------------
# FASTQ statistics, filtering, pairs
# ---------------------------------------------------------------------------


def _read_report(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")[:_MAX_REPORT_CHARS]


def vsearch_fastq_eestats(
    input_fastq: str | Path,
    output_directory: str | Path,
    *,
    fastq_ascii: int | None = None,
    fastq_qmin: int | None = None,
    fastq_qmax: int | None = None,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Per-position quality-score and expected-error statistics of a FASTQ file (``vsearch --fastq_eestats``).

    Writes ``eestats.tsv``; ``records`` has one row per read position with the number of reads and the min, low, median,
    mean, high and max of Q, Pe (error probability) and EE (cumulative expected errors).
    """
    source = _sequence_input(input_fastq, "input_fastq")
    directory = _begin(vsearch_path, output_directory)
    outputs = {"statistics": directory / "eestats.tsv"}
    arguments = ["--fastq_eestats", source, "--output", str(outputs["statistics"])]
    _quality_range(arguments, fastq_ascii, fastq_qmin, fastq_qmax)
    result = _execute(
        "fastq_eestats",
        arguments,
        directory,
        outputs,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )
    lines = outputs["statistics"].read_text(encoding="utf-8", errors="replace").splitlines()
    names = [name.lower() for name in lines[0].split("\t")] if lines else []
    result["records"] = [
        {name: _convert(value) for name, value in zip(names, line.split("\t"), strict=False)}
        for line in lines[1 : _ROW_LIMIT + 1]
    ]
    result["position_count"] = max(len(lines) - 1, 0)
    return result


def vsearch_fastq_eestats2(
    input_fastq: str | Path,
    output_directory: str | Path,
    *,
    ee_cutoffs: Sequence[float] | None = None,
    length_cutoffs: str | None = None,
    fastq_ascii: int | None = None,
    fastq_qmin: int | None = None,
    fastq_qmax: int | None = None,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read retention by length cutoff and maximum expected error (``vsearch --fastq_eestats2``).

    ``ee_cutoffs`` are the MaxEE values (VSEARCH default 0.5, 1, 2); ``length_cutoffs`` is ``"min,max,increment"`` (default
    ``"50,*,50"``; ``*`` means the longest read). Writes ``eestats2.txt``; ``records`` has one row per length cutoff with
    ``reads`` and ``percent`` for every MaxEE value.
    """
    source = _sequence_input(input_fastq, "input_fastq")
    directory = _begin(vsearch_path, output_directory)
    outputs = {"statistics": directory / "eestats2.txt"}
    arguments = ["--fastq_eestats2", source, "--output", str(outputs["statistics"])]
    if ee_cutoffs is not None:
        if isinstance(ee_cutoffs, str) or not ee_cutoffs:
            raise ConfigurationError("ee_cutoffs must be a non-empty list of numbers.")
        values = [
            format(_cli.number(item, "ee_cutoffs", minimum=0.0), ".12g") for item in ee_cutoffs
        ]
        arguments.extend(("--ee_cutoffs", ",".join(values)))
    if length_cutoffs is not None:
        if not re.fullmatch(r"\d+,(\d+|\*),\d+", length_cutoffs):
            raise ConfigurationError(
                "length_cutoffs must look like 'min,max,increment' (max may be '*')."
            )
        arguments.extend(("--length_cutoffs", length_cutoffs))
    _quality_range(arguments, fastq_ascii, fastq_qmin, fastq_qmax)
    result = _execute(
        "fastq_eestats2",
        arguments,
        directory,
        outputs,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )
    lines = outputs["statistics"].read_text(encoding="utf-8", errors="replace").splitlines()
    cutoffs: list[float] = []
    records: list[dict[str, Any]] = []
    for line in lines:
        if line.startswith("Length"):
            cutoffs = [float(item) for item in re.findall(r"MaxEE\s+([\d.]+)", line)]
        elif cutoffs and re.match(r"^\s*\d+\s+\d+\(", line):
            length = int(line.split()[0])
            cells = re.findall(r"(\d+)\(\s*([\d.]+)%\)", line)
            records.append(
                {
                    "length": length,
                    "retained": [
                        {"max_ee": ee, "reads": int(reads), "percent": float(percent)}
                        for ee, (reads, percent) in zip(cutoffs, cells, strict=False)
                    ],
                }
            )
    result["records"] = records
    summary = re.search(r"(\d+) reads, max len (\d+), avg ([\d.]+)", lines[0] if lines else "")
    if summary:
        result["reads"], result["max_length"], result["mean_length"] = (
            int(summary.group(1)),
            int(summary.group(2)),
            float(summary.group(3)),
        )
    return result


def vsearch_fastq_stats(
    input_fastq: str | Path,
    output_directory: str | Path,
    *,
    fastq_ascii: int | None = None,
    fastq_qmin: int | None = None,
    fastq_qmax: int | None = None,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Detailed read-length, quality and truncation statistics of a FASTQ file (``vsearch --fastq_stats``).

    Writes the full VSEARCH report to ``stats.txt``. ``read_length_distribution`` (length, reads, percent, cumulative
    percent) and ``quality_distribution`` (ASCII character, Q, error probability, bases, percent, cumulative percent) are
    parsed from it, with the summary ``read_count``, ``mean_length`` and ``bases``.
    """
    source = _sequence_input(input_fastq, "input_fastq")
    directory = _begin(vsearch_path, output_directory)
    outputs = {"report": directory / "stats.txt"}
    arguments = ["--fastq_stats", source, "--log", str(outputs["report"])]
    _quality_range(arguments, fastq_ascii, fastq_qmin, fastq_qmax)
    result = _execute(
        "fastq_stats",
        arguments,
        directory,
        outputs,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )
    report = outputs["report"].read_text(encoding="utf-8", errors="replace")
    lengths: list[dict[str, Any]] = []
    qualities: list[dict[str, Any]] = []
    section = ""
    for line in report.splitlines():
        if line.startswith("Read length distribution"):
            section = "length"
        elif line.startswith("Q score distribution"):
            section = "quality"
        elif line.startswith("Truncate at first Q") or line.startswith("Truncate at first EE"):
            section = ""
        elif section == "length":
            match = re.match(r"^\s*(?:>=\s*)?(\d+)\s+(\d+)\s+([\d.]+)%\s+([\d.]+)%\s*$", line)
            if match:
                lengths.append(
                    {
                        "length": int(match.group(1)),
                        "reads": int(match.group(2)),
                        "percent": float(match.group(3)),
                        "cumulative_percent": float(match.group(4)),
                    }
                )
        elif section == "quality":
            match = re.match(
                r"^\s*(\S)\s+(\d+)\s+([\d.]+)\s+(\d+)\s+([\d.]+)%\s+([\d.]+)%\s*$", line
            )
            if match:
                qualities.append(
                    {
                        "ascii": match.group(1),
                        "quality": int(match.group(2)),
                        "error_probability": float(match.group(3)),
                        "bases": int(match.group(4)),
                        "percent": float(match.group(5)),
                        "cumulative_percent": float(match.group(6)),
                    }
                )
    result["read_length_distribution"] = lengths
    result["quality_distribution"] = qualities
    summary = re.search(r"^\s*(\d+)\s+Recs", report, re.MULTILINE)
    if summary:
        result["read_count"] = int(summary.group(1))
    mean = re.search(r"^\s*([\d.]+)\s+Avg length", report, re.MULTILINE)
    if mean:
        result["mean_length"] = float(mean.group(1))
    bases = re.search(r"^\s*([\d.]+[KMG]?)\s+Bases", report, re.MULTILINE)
    if bases:
        result["bases"] = bases.group(1)
    return result


def vsearch_filter_reads(
    input_path: str | Path,
    output_directory: str | Path,
    *,
    output_format: Literal["fasta", "fastq"] | None = None,
    reverse_path: str | Path | None = None,
    min_quality: int | None = None,
    truncate_quality: int | None = None,
    truncate_ee: float | None = None,
    truncate_ee_rate: float | None = None,
    max_ee: float | None = None,
    max_ee_rate: float | None = None,
    min_length: int | None = None,
    max_length: int | None = None,
    max_ns: int | None = None,
    strip_left: int | None = None,
    strip_right: int | None = None,
    truncate_length: int | None = None,
    truncate_length_keep: int | None = None,
    min_size: int | None = None,
    max_size: int | None = None,
    write_discarded: bool = False,
    include_expected_errors: bool = False,
    relabel: str | None = None,
    relabel_mode: Literal["keep", "md5", "sha1", "self"] | None = None,
    abundance_out: bool = False,
    fastq_ascii: int | None = None,
    fastq_qmin: int | None = None,
    fastq_qmax: int | None = None,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Trim and filter FASTA or FASTQ records (``vsearch --fastx_filter``).

    One function for the filtering criteria of ``--fastx_filter``: ``min_quality`` (``--fastq_minqual``, discard a read with
    any base below it), ``truncate_quality`` (``--fastq_truncqual``), ``truncate_ee`` / ``truncate_ee_rate``
    (``--fastq_truncee`` / ``--fastq_truncee_rate``: truncate at an expected-error threshold), ``max_ee`` /
    ``max_ee_rate``, ``min_length`` / ``max_length``, ``max_ns``, ``strip_left`` / ``strip_right``, ``truncate_length``
    (discard shorter reads) / ``truncate_length_keep``, and the abundance annotation filters ``min_size`` / ``max_size``
    (``--minsize`` / ``--maxsize``, read from ``;size=``). With ``reverse_path`` paired reads are filtered together.
    ``output_format`` defaults to the format of the input. Writes ``filtered.<fasta|fastq>`` (and ``filtered_reverse``);
    ``write_discarded`` adds ``discarded`` files. Expected-error and quality criteria need FASTQ input.
    """
    source = _sequence_input(input_path, "input_path")
    kind = _output_format(output_format) if output_format is not None else _fastx_kind(source)
    reverse = None if reverse_path is None else _sequence_input(reverse_path, "reverse_path")
    directory = _begin(vsearch_path, output_directory)
    out_flag = "--fastqout" if kind == "fastq" else "--fastaout"
    extension = "fastq" if kind == "fastq" else "fasta"
    outputs = {"filtered": directory / f"filtered.{extension}"}
    arguments = ["--fastx_filter", source, out_flag, str(outputs["filtered"])]
    if reverse is not None:
        outputs["filtered_reverse"] = directory / f"filtered_reverse.{extension}"
        arguments.extend(
            ("--reverse", reverse, f"{out_flag}_rev", str(outputs["filtered_reverse"]))
        )
    if _flag(write_discarded, "write_discarded"):
        outputs["discarded"] = directory / f"discarded.{extension}"
        arguments.extend((f"{out_flag}_discarded", str(outputs["discarded"])))
        if reverse is not None:
            outputs["discarded_reverse"] = directory / f"discarded_reverse.{extension}"
            arguments.extend((f"{out_flag}_discarded_rev", str(outputs["discarded_reverse"])))
    _opt_int(arguments, "--fastq_minqual", min_quality, "min_quality", minimum=0, maximum=93)
    _opt_int(
        arguments, "--fastq_truncqual", truncate_quality, "truncate_quality", minimum=0, maximum=93
    )
    _opt_num(arguments, "--fastq_truncee", truncate_ee, "truncate_ee", minimum=0.0)
    _opt_num(arguments, "--fastq_truncee_rate", truncate_ee_rate, "truncate_ee_rate", minimum=0.0)
    _opt_num(arguments, "--fastq_maxee", max_ee, "max_ee", minimum=0.0)
    _opt_num(arguments, "--fastq_maxee_rate", max_ee_rate, "max_ee_rate", minimum=0.0)
    _opt_int(arguments, "--fastq_minlen", min_length, "min_length", minimum=0)
    _opt_int(arguments, "--fastq_maxlen", max_length, "max_length", minimum=0)
    _opt_int(arguments, "--fastq_maxns", max_ns, "max_ns", minimum=0)
    _opt_int(arguments, "--fastq_stripleft", strip_left, "strip_left", minimum=0)
    _opt_int(arguments, "--fastq_stripright", strip_right, "strip_right", minimum=0)
    _opt_int(arguments, "--fastq_trunclen", truncate_length, "truncate_length", minimum=1)
    _opt_int(
        arguments, "--fastq_trunclen_keep", truncate_length_keep, "truncate_length_keep", minimum=1
    )
    _opt_int(arguments, "--minsize", min_size, "min_size", minimum=1)
    _opt_int(arguments, "--maxsize", max_size, "max_size", minimum=1)
    _switch(arguments, "--eeout", include_expected_errors, "include_expected_errors")
    _switch(arguments, "--sizeout", abundance_out, "abundance_out")
    _relabel(arguments, relabel, relabel_mode)
    _quality_range(arguments, fastq_ascii, fastq_qmin, fastq_qmax)
    result = _execute(
        "filter_reads",
        arguments,
        directory,
        outputs,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )
    result["passed_count"] = _count(outputs["filtered"], kind)
    if "discarded" in outputs:
        result["discarded_count"] = _count(outputs["discarded"], kind)
    return result


def vsearch_sync_pairs(
    forward_path: str | Path,
    reverse_path: str | Path,
    output_directory: str | Path,
    *,
    output_format: Literal["fasta", "fastq"] | None = None,
    read_separators: str | None = None,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Reorder paired reads so that both files match and set aside orphans (``vsearch --fastx_syncpairs``).

    Writes ``synced_forward``, ``synced_reverse``, ``orphans_forward`` and ``orphans_reverse`` (extension from
    ``output_format``, which defaults to the input format). ``read_separators`` are the characters that precede the mate
    number in the labels (VSEARCH default ``/``).
    """
    forward = _sequence_input(forward_path, "forward_path")
    reverse = _sequence_input(reverse_path, "reverse_path")
    kind = _output_format(output_format) if output_format is not None else _fastx_kind(forward)
    directory = _begin(vsearch_path, output_directory)
    extension = "fastq" if kind == "fastq" else "fasta"
    prefix = "--fastqout" if kind == "fastq" else "--fastaout"
    outputs = {
        "synced_forward": directory / f"synced_forward.{extension}",
        "synced_reverse": directory / f"synced_reverse.{extension}",
        "orphans_forward": directory / f"orphans_forward.{extension}",
        "orphans_reverse": directory / f"orphans_reverse.{extension}",
    }
    arguments = [
        "--fastx_syncpairs",
        forward,
        "--reverse",
        reverse,
        prefix,
        str(outputs["synced_forward"]),
        f"{prefix}_rev",
        str(outputs["synced_reverse"]),
        f"{prefix}_orphans",
        str(outputs["orphans_forward"]),
        f"{prefix}_orphans_rev",
        str(outputs["orphans_reverse"]),
    ]
    _opt_text(arguments, "--read_separators", read_separators, "read_separators")
    result = _execute(
        "sync_pairs",
        arguments,
        directory,
        outputs,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )
    result["synced_pairs"] = _count(outputs["synced_forward"], kind)
    result["forward_orphans"] = _count(outputs["orphans_forward"], kind)
    result["reverse_orphans"] = _count(outputs["orphans_reverse"], kind)
    return result


def vsearch_merge_pairs(
    forward_fastq: str | Path,
    reverse_fastq: str | Path,
    output_directory: str | Path,
    *,
    allow_stagger: bool = False,
    max_diffs: int | None = None,
    max_diff_pct: float | None = None,
    max_ee: float | None = None,
    max_ns: int | None = None,
    min_length: int | None = None,
    min_overlap: int | None = None,
    min_merge_length: int | None = None,
    max_merge_length: int | None = None,
    truncate_quality: int | None = None,
    label_suffix: str | None = None,
    include_expected_errors: bool = False,
    write_error_table: bool = False,
    fastq_ascii: int | None = None,
    fastq_qmin: int | None = None,
    fastq_qmax: int | None = None,
    fastq_qmin_out: int | None = None,
    fastq_qmax_out: int | None = None,
    threads: int | None = None,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge overlapping paired-end reads into one sequence (``vsearch --fastq_mergepairs``).

    Writes ``merged.fastq``, ``notmerged_forward.fastq`` and ``notmerged_reverse.fastq``. Parameters map to
    ``--fastq_allowmergestagger``, ``--fastq_maxdiffs`` (default 10), ``--fastq_maxdiffpct`` (100), ``--fastq_maxee``,
    ``--fastq_maxns``, ``--fastq_minlen``, ``--fastq_minovlen`` (10), ``--fastq_minmergelen``, ``--fastq_maxmergelen`` and
    ``--fastq_truncqual``. ``write_error_table`` adds ``errors.tsv`` (``--eetabbedout``).
    """
    forward = _sequence_input(forward_fastq, "forward_fastq")
    reverse = _sequence_input(reverse_fastq, "reverse_fastq")
    directory = _begin(vsearch_path, output_directory)
    outputs = {
        "merged": directory / "merged.fastq",
        "notmerged_forward": directory / "notmerged_forward.fastq",
        "notmerged_reverse": directory / "notmerged_reverse.fastq",
    }
    arguments = [
        "--fastq_mergepairs",
        forward,
        "--reverse",
        reverse,
        "--fastqout",
        str(outputs["merged"]),
        "--fastqout_notmerged_fwd",
        str(outputs["notmerged_forward"]),
        "--fastqout_notmerged_rev",
        str(outputs["notmerged_reverse"]),
    ]
    if _flag(write_error_table, "write_error_table"):
        outputs["errors"] = directory / "errors.tsv"
        arguments.extend(("--eetabbedout", str(outputs["errors"])))
    _switch(arguments, "--fastq_allowmergestagger", allow_stagger, "allow_stagger")
    _opt_int(arguments, "--fastq_maxdiffs", max_diffs, "max_diffs", minimum=0)
    _opt_num(
        arguments, "--fastq_maxdiffpct", max_diff_pct, "max_diff_pct", minimum=0.0, maximum=100.0
    )
    _opt_num(arguments, "--fastq_maxee", max_ee, "max_ee", minimum=0.0)
    _opt_int(arguments, "--fastq_maxns", max_ns, "max_ns", minimum=0)
    _opt_int(arguments, "--fastq_minlen", min_length, "min_length", minimum=0)
    _opt_int(arguments, "--fastq_minovlen", min_overlap, "min_overlap", minimum=0)
    _opt_int(arguments, "--fastq_minmergelen", min_merge_length, "min_merge_length", minimum=0)
    _opt_int(arguments, "--fastq_maxmergelen", max_merge_length, "max_merge_length", minimum=0)
    _opt_int(
        arguments, "--fastq_truncqual", truncate_quality, "truncate_quality", minimum=0, maximum=93
    )
    _opt_text(arguments, "--label_suffix", label_suffix, "label_suffix")
    _switch(arguments, "--fastq_eeout", include_expected_errors, "include_expected_errors")
    _quality_range(arguments, fastq_ascii, fastq_qmin, fastq_qmax)
    _opt_int(arguments, "--fastq_qminout", fastq_qmin_out, "fastq_qmin_out", minimum=0, maximum=93)
    _opt_int(arguments, "--fastq_qmaxout", fastq_qmax_out, "fastq_qmax_out", minimum=0, maximum=93)
    _threads(arguments, threads)
    result = _execute(
        "merge_pairs",
        arguments,
        directory,
        outputs,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )
    result["merged_count"] = _count(outputs["merged"], "fastq")
    result["not_merged_count"] = _count(outputs["notmerged_forward"], "fastq")
    return result


def vsearch_join_pairs(
    forward_fastq: str | Path,
    reverse_fastq: str | Path,
    output_directory: str | Path,
    *,
    pad_sequence: str | None = None,
    pad_quality: str | None = None,
    write_fasta: bool = False,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Join paired-end reads into one sequence with a gap between them (``vsearch --fastq_join``).

    The reverse read is reverse-complemented and appended to the forward read after ``pad_sequence`` (``--join_padgap``,
    default ``NNNNNNNN``) with ``pad_quality`` (``--join_padgapq``, default ``IIIIIIII``). Writes ``joined.fastq``;
    ``write_fasta`` adds ``joined.fasta``.
    """
    forward = _sequence_input(forward_fastq, "forward_fastq")
    reverse = _sequence_input(reverse_fastq, "reverse_fastq")
    if pad_sequence is not None and not re.fullmatch(r"[ACGTNacgtn]+", pad_sequence):
        raise ConfigurationError("pad_sequence must consist of the letters ACGTN.")
    if pad_quality is not None and not re.fullmatch(r"[!-~]+", pad_quality):
        raise ConfigurationError("pad_quality must be printable FASTQ quality characters.")
    if pad_quality is not None and pad_quality.startswith("-"):
        raise ConfigurationError("pad_quality must not start with '-'.")
    directory = _begin(vsearch_path, output_directory)
    outputs = {"joined": directory / "joined.fastq"}
    arguments = [
        "--fastq_join",
        forward,
        "--reverse",
        reverse,
        "--fastqout",
        str(outputs["joined"]),
    ]
    if _flag(write_fasta, "write_fasta"):
        outputs["joined_fasta"] = directory / "joined.fasta"
        arguments.extend(("--fastaout", str(outputs["joined_fasta"])))
    if pad_sequence is not None:
        arguments.extend(("--join_padgap", pad_sequence))
    if pad_quality is not None:
        arguments.extend(("--join_padgapq", pad_quality))
    result = _execute(
        "join_pairs",
        arguments,
        directory,
        outputs,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )
    result["joined_count"] = _count(outputs["joined"], "fastq")
    return result


# ---------------------------------------------------------------------------
# Retrieval, orientation, cutting, scrambling, conversion
# ---------------------------------------------------------------------------


def vsearch_extract_records(
    input_path: str | Path,
    output_directory: str | Path,
    *,
    label: str | None = None,
    labels_file: str | Path | None = None,
    label_word: str | None = None,
    label_words_file: str | Path | None = None,
    label_field: str | None = None,
    substring_match: bool = False,
    output_format: Literal["fasta", "fastq"] | None = None,
    write_unmatched: bool = False,
    abundance_in: bool = False,
    abundance_out: bool = False,
    relabel: str | None = None,
    relabel_mode: Literal["keep", "md5", "sha1", "self"] | None = None,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract records by header label from FASTA or FASTQ (``vsearch --fastx_getseqs`` / ``--fastx_getseq``).

    Give at least one selector: ``label`` (one label; this is the single-record ``--fastx_getseq`` case), ``labels_file``
    (one label per line), ``label_word`` or ``label_words_file`` (match a word of the header; ``label_field`` names the
    ``key=value`` field). By default the whole header (up to the first space) must equal the label;
    ``substring_match`` allows a match anywhere. Writes ``extracted.<fasta|fastq>``; ``write_unmatched`` adds
    ``unmatched``.
    """
    selectors = (label, labels_file, label_word, label_words_file)
    if all(item is None for item in selectors):
        raise ConfigurationError(
            "Give at least one of label, labels_file, label_word or label_words_file."
        )
    source = _sequence_input(input_path, "input_path")
    kind = _output_format(output_format) if output_format is not None else _fastx_kind(source)
    list_file = None if labels_file is None else _sequence_input(labels_file, "labels_file")
    words_file = (
        None if label_words_file is None else _sequence_input(label_words_file, "label_words_file")
    )
    directory = _begin(vsearch_path, output_directory)
    extension = "fastq" if kind == "fastq" else "fasta"
    outputs = {"extracted": directory / f"extracted.{extension}"}
    arguments = [
        "--fastx_getseqs",
        source,
        "--fastqout" if kind == "fastq" else "--fastaout",
        str(outputs["extracted"]),
    ]
    if _flag(write_unmatched, "write_unmatched"):
        outputs["unmatched"] = directory / f"unmatched.{extension}"
        arguments.extend(
            ("--notmatchedfq" if kind == "fastq" else "--notmatched", str(outputs["unmatched"]))
        )
    _opt_text(arguments, "--label", label, "label")
    if list_file is not None:
        arguments.extend(("--labels", list_file))
    _opt_text(arguments, "--label_word", label_word, "label_word")
    if words_file is not None:
        arguments.extend(("--label_words", words_file))
    _opt_text(arguments, "--label_field", label_field, "label_field")
    _switch(arguments, "--label_substr_match", substring_match, "substring_match")
    _switch(arguments, "--sizein", abundance_in, "abundance_in")
    _switch(arguments, "--sizeout", abundance_out, "abundance_out")
    _relabel(arguments, relabel, relabel_mode)
    result = _execute(
        "extract_records",
        arguments,
        directory,
        outputs,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )
    result["extracted_count"] = _count(outputs["extracted"], kind)
    return result


def vsearch_extract_subsequence(
    input_path: str | Path,
    output_directory: str | Path,
    *,
    label: str,
    start: int,
    end: int,
    substring_match: bool = False,
    output_format: Literal["fasta", "fastq"] | None = None,
    label_suffix: str | None = None,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract a subsequence of one record by label and 1-based inclusive positions (``vsearch --fastx_getsubseq``).

    Writes ``subsequence.<fasta|fastq>``. ``label_suffix`` is appended to the label of the extracted record.
    """
    source = _sequence_input(input_path, "input_path")
    kind = _output_format(output_format) if output_format is not None else _fastx_kind(source)
    first = _cli.integer(start, "start", minimum=1)
    last = _cli.integer(end, "end", minimum=first)
    directory = _begin(vsearch_path, output_directory)
    outputs = {"subsequence": directory / f"subsequence.{'fastq' if kind == 'fastq' else 'fasta'}"}
    arguments = [
        "--fastx_getsubseq",
        source,
        "--label",
        _cli.text(label, "label"),
        "--subseq_start",
        str(first),
        "--subseq_end",
        str(last),
        "--fastqout" if kind == "fastq" else "--fastaout",
        str(outputs["subsequence"]),
    ]
    _switch(arguments, "--label_substr_match", substring_match, "substring_match")
    _opt_text(arguments, "--label_suffix", label_suffix, "label_suffix")
    result = _execute(
        "extract_subsequence",
        arguments,
        directory,
        outputs,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )
    result["extracted_count"] = _count(outputs["subsequence"], kind)
    return result


def vsearch_orient_sequences(
    input_path: str | Path,
    reference_fasta: str | Path,
    output_directory: str | Path,
    *,
    output_format: Literal["fasta", "fastq"] | None = None,
    word_length: int | None = None,
    query_mask: Literal["none", "dust", "soft"] | None = None,
    database_mask: Literal["none", "dust", "soft"] | None = None,
    threads: int | None = None,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Orient sequences to the strand of a reference database by k-mer matching (``vsearch --orient``).

    Reverse-complements the queries that match the reference on the opposite strand. Writes ``oriented.<fasta|fastq>``,
    ``undetermined`` (sequences that could not be oriented) and ``orientation.tsv`` (query, strand ``+``/``-``/``?``,
    forward and reverse matches; parsed into ``records``).
    """
    source = _sequence_input(input_path, "input_path")
    db = _sequence_input(reference_fasta, "reference_fasta")
    kind = _output_format(output_format) if output_format is not None else _fastx_kind(source)
    directory = _begin(vsearch_path, output_directory)
    extension = "fastq" if kind == "fastq" else "fasta"
    outputs = {
        "oriented": directory / f"oriented.{extension}",
        "undetermined": directory / f"undetermined.{extension}",
        "orientation": directory / "orientation.tsv",
    }
    arguments = [
        "--orient",
        source,
        "--db",
        db,
        "--fastqout" if kind == "fastq" else "--fastaout",
        str(outputs["oriented"]),
        "--notmatched",
        str(outputs["undetermined"]),
        "--tabbedout",
        str(outputs["orientation"]),
    ]
    _opt_int(arguments, "--wordlength", word_length, "word_length", minimum=3, maximum=15)
    _opt_choice(arguments, "--qmask", query_mask, ("none", "dust", "soft"), "query_mask")
    _opt_choice(arguments, "--dbmask", database_mask, ("none", "dust", "soft"), "database_mask")
    _threads(arguments, threads)
    result = _execute(
        "orient_sequences",
        arguments,
        directory,
        outputs,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )
    result["records"], result["record_count"] = _rows(
        outputs["orientation"], ("query", "strand", "forward_matches", "reverse_matches")
    )
    return result


def vsearch_cut_restriction_sites(
    input_fasta: str | Path,
    output_directory: str | Path,
    *,
    cut_pattern: str,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cut sequences at restriction sites (``vsearch --cut``).

    ``cut_pattern`` is a site written with ``^`` and ``_`` marking the cut on the forward and reverse strand, for example
    ``"G^AATT_C"`` (EcoRI). Writes ``fragments_forward.fasta``, ``fragments_reverse.fasta`` (fragments on the reverse
    strand), ``discarded_forward.fasta`` and ``discarded_reverse.fasta`` (sequences without the site).
    """
    pattern = _cli.text(cut_pattern, "cut_pattern")
    if _CUT_PATTERN.fullmatch(pattern) is None or "^" not in pattern or "_" not in pattern:
        raise ConfigurationError(
            "cut_pattern must use nucleotide letters and contain both '^' and '_'."
        )
    source = _sequence_input(input_fasta, "input_fasta")
    directory = _begin(vsearch_path, output_directory)
    outputs = {
        "fragments_forward": directory / "fragments_forward.fasta",
        "fragments_reverse": directory / "fragments_reverse.fasta",
        "discarded_forward": directory / "discarded_forward.fasta",
        "discarded_reverse": directory / "discarded_reverse.fasta",
    }
    arguments = [
        "--cut",
        source,
        "--cut_pattern",
        pattern,
        "--fastaout",
        str(outputs["fragments_forward"]),
        "--fastaout_rev",
        str(outputs["fragments_reverse"]),
        "--fastaout_discarded",
        str(outputs["discarded_forward"]),
        "--fastaout_discarded_rev",
        str(outputs["discarded_reverse"]),
    ]
    result = _execute(
        "cut_restriction_sites",
        arguments,
        directory,
        outputs,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )
    result["forward_fragment_count"] = _count(outputs["fragments_forward"], "fasta")
    result["reverse_fragment_count"] = _count(outputs["fragments_reverse"], "fasta")
    return result


def vsearch_scramble_sequences(
    input_path: str | Path,
    output_directory: str | Path,
    *,
    output_format: Literal["fasta", "fastq"] | None = None,
    kmer_size: int | None = None,
    random_seed: int | None = None,
    label_suffix: str | None = None,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Randomise the nucleotide order within each sequence (``vsearch --scramble``).

    ``kmer_size`` (``--scramble_kmer``, 1-9, default 1) preserves the counts of words up to that size. ``random_seed``
    (``--randseed``) makes the result reproducible (0 or omitted uses a random source). Writes
    ``scrambled.<fasta|fastq>``.
    """
    source = _sequence_input(input_path, "input_path")
    kind = _output_format(output_format) if output_format is not None else _fastx_kind(source)
    directory = _begin(vsearch_path, output_directory)
    outputs = {"scrambled": directory / f"scrambled.{'fastq' if kind == 'fastq' else 'fasta'}"}
    arguments = [
        "--scramble",
        source,
        "--fastqout" if kind == "fastq" else "--fastaout",
        str(outputs["scrambled"]),
    ]
    _opt_int(arguments, "--scramble_kmer", kmer_size, "kmer_size", minimum=1, maximum=9)
    _opt_int(arguments, "--randseed", random_seed, "random_seed", minimum=0)
    _opt_text(arguments, "--label_suffix", label_suffix, "label_suffix")
    result = _execute(
        "scramble_sequences",
        arguments,
        directory,
        outputs,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )
    result["record_count"] = _count(outputs["scrambled"], kind)
    return result


def vsearch_convert_sff(
    input_sff: str | Path,
    output_directory: str | Path,
    *,
    clip_ends: bool = False,
    fastq_ascii_out: int | None = None,
    fastq_qmin_out: int | None = None,
    fastq_qmax_out: int | None = None,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert an SFF file (454/Ion Torrent flowgram format) to FASTQ (``vsearch --sff_convert``).

    ``clip_ends`` applies the clip points stored in the file (``--sff_clip``). Writes ``converted.fastq``.
    """
    source = _sequence_input(input_sff, "input_sff")
    directory = _begin(vsearch_path, output_directory)
    outputs = {"converted": directory / "converted.fastq"}
    arguments = ["--sff_convert", source, "--fastqout", str(outputs["converted"])]
    _switch(arguments, "--sff_clip", clip_ends, "clip_ends")
    _opt_int(
        arguments, "--fastq_asciiout", fastq_ascii_out, "fastq_ascii_out", minimum=33, maximum=64
    )
    _opt_int(arguments, "--fastq_qminout", fastq_qmin_out, "fastq_qmin_out", minimum=0, maximum=93)
    _opt_int(arguments, "--fastq_qmaxout", fastq_qmax_out, "fastq_qmax_out", minimum=0, maximum=93)
    result = _execute(
        "convert_sff",
        arguments,
        directory,
        outputs,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )
    result["record_count"] = _count(outputs["converted"], "fastq")
    return result


# ---------------------------------------------------------------------------
# Taxonomy and UDB indexes
# ---------------------------------------------------------------------------


def vsearch_sintax(
    input_path: str | Path,
    database: str | Path,
    output_directory: str | Path,
    *,
    cutoff: float | None = None,
    random_seed: int | None = None,
    random_ties: bool = False,
    threads: int | None = None,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """SINTAX taxonomic classification of sequences (``vsearch --sintax``).

    ``database`` is a FASTA or UDB file whose labels carry ``;tax=d:Domain,p:Phylum,...`` annotations. ``cutoff``
    (``--sintax_cutoff``, default 0) is the minimum bootstrap confidence for the reported taxonomy; ``random_seed`` fixes
    the bootstrap; ``random_ties`` picks a random sequence among equal best hits. Writes ``classification.tsv`` and
    parses it into ``records`` (query, per-rank hits with confidence, strand, taxonomy above the cutoff).
    """
    source = _sequence_input(input_path, "input_path")
    db = _sequence_input(database, "database")
    directory = _begin(vsearch_path, output_directory)
    outputs = {"classification": directory / "classification.tsv"}
    arguments = ["--sintax", source, "--db", db, "--tabbedout", str(outputs["classification"])]
    _opt_num(arguments, "--sintax_cutoff", cutoff, "cutoff", minimum=0.0, maximum=1.0)
    _opt_int(arguments, "--randseed", random_seed, "random_seed", minimum=0)
    _switch(arguments, "--sintax_random", random_ties, "random_ties")
    _threads(arguments, threads)
    result = _execute(
        "sintax",
        arguments,
        directory,
        outputs,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )
    result["records"], result["record_count"] = _rows(
        outputs["classification"], ("query", "hits", "strand", "taxonomy")
    )
    return result


def vsearch_make_udb(
    input_fasta: str | Path,
    output_directory: str | Path,
    *,
    word_length: int | None = None,
    database_mask: Literal["none", "dust", "soft"] | None = None,
    hardmask: bool = False,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a UDB word-index database from a FASTA file for later searches (``vsearch --makeudb_usearch``).

    ``word_length`` (``--wordlength``, 3-15, default 8) and ``database_mask`` (``--dbmask``, default dust) are stored in the
    index; ``hardmask`` masks with N instead of lower case. Writes ``database.udb``.
    """
    source = _sequence_input(input_fasta, "input_fasta")
    directory = _begin(vsearch_path, output_directory)
    outputs = {"database": directory / "database.udb"}
    arguments = ["--makeudb_usearch", source, "--output", str(outputs["database"])]
    _opt_int(arguments, "--wordlength", word_length, "word_length", minimum=3, maximum=15)
    _opt_choice(arguments, "--dbmask", database_mask, ("none", "dust", "soft"), "database_mask")
    _switch(arguments, "--hardmask", hardmask, "hardmask")
    result = _execute(
        "make_udb",
        arguments,
        directory,
        outputs,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )
    result["database_bytes"] = outputs["database"].stat().st_size
    return result


def vsearch_udb_to_fasta(
    input_udb: str | Path,
    output_directory: str | Path,
    *,
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Export the sequences stored in a UDB file as FASTA (``vsearch --udb2fasta``). Writes ``sequences.fasta``."""
    source = _sequence_input(input_udb, "input_udb")
    directory = _begin(vsearch_path, output_directory)
    outputs = {"sequences": directory / "sequences.fasta"}
    arguments = ["--udb2fasta", source, "--output", str(outputs["sequences"])]
    result = _execute(
        "udb_to_fasta",
        arguments,
        directory,
        outputs,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )
    result["record_count"] = _count(outputs["sequences"], "fasta")
    return result


def vsearch_udb_report(
    input_udb: str | Path,
    output_directory: str | Path,
    *,
    report: Literal["info", "stats"] = "info",
    vsearch_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Describe a UDB file: ``report='info'`` is ``vsearch --udbinfo``, ``report='stats'`` is ``vsearch --udbstats``.

    ``info`` gives the number of sequences, alphabet, word width and index settings (parsed into ``fields``);
    ``stats`` gives the k-mer word statistics (word-size distribution and the most frequent words, kept in the report
    text). Writes ``report.txt``.
    """
    _cli.choice(report, ("info", "stats"), "report")
    source = _sequence_input(input_udb, "input_udb")
    directory = _begin(vsearch_path, output_directory)
    outputs = {"report": directory / "report.txt"}
    arguments = [f"--udb{report}", source, "--log", str(outputs["report"])]
    result = _execute(
        f"udb_{report}",
        arguments,
        directory,
        outputs,
        vsearch_path=vsearch_path,
        timeout_seconds=timeout_seconds,
    )
    text = _read_report(outputs["report"])
    body = text.split("Started", 1)[-1].split("\n", 1)[-1].split("Finished", 1)[0]
    result["report_text"] = body.strip("\n")
    if report == "info":
        fields: dict[str, Any] = {}
        for line in body.splitlines():
            match = re.match(r"^\s+(\S[\w ]*?)\s{2,}(.+)$", line)
            if match:
                fields[match.group(1)] = (
                    _convert(match.group(2).split(" (", 1)[0])
                    if match.group(2)[0].isdigit()
                    else match.group(2)
                )
        result["fields"] = fields
    return result


__all__ = [
    "vsearch_allpairs_global",
    "vsearch_cluster_unoise",
    "vsearch_convert_sff",
    "vsearch_cut_restriction_sites",
    "vsearch_dereplicate_by_id",
    "vsearch_dereplicate_fulllength",
    "vsearch_dereplicate_prefix",
    "vsearch_detect_chimeras",
    "vsearch_extract_records",
    "vsearch_extract_subsequence",
    "vsearch_fastq_eestats",
    "vsearch_fastq_eestats2",
    "vsearch_fastq_stats",
    "vsearch_filter_reads",
    "vsearch_join_pairs",
    "vsearch_make_udb",
    "vsearch_merge_pairs",
    "vsearch_orient_sequences",
    "vsearch_rereplicate",
    "vsearch_scramble_sequences",
    "vsearch_search_global",
    "vsearch_sintax",
    "vsearch_sort_by_size",
    "vsearch_sync_pairs",
    "vsearch_udb_report",
    "vsearch_udb_to_fasta",
]

_register_public_verbose(globals(), tuple(__all__))
