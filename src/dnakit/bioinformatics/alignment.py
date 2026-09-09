"""Native alignment objects with format codecs and explicit external algorithms."""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from dnakit.bioinformatics._common import (
    ExecutionConfig,
    ToolResult,
    artifacts,
    input_path,
    require,
    run,
    run_to_file,
    sequence_text,
    workspace,
)
from dnakit.exceptions import BackendExecutionError, ConfigurationError, InputFormatError


@dataclass(frozen=True)
class MultipleAlignment:
    """Aligned DNA strings with explicit '-' columns and stable unique IDs."""

    ids: tuple[str, ...]
    sequences: tuple[str, ...]

    def __post_init__(self) -> None:
        if (
            not self.ids
            or len(self.ids) != len(self.sequences)
            or len(set(self.ids)) != len(self.ids)
        ):
            raise ConfigurationError("An alignment needs matching sequences and unique IDs.")
        if any(
            not isinstance(key, str) or not key or any(c.isspace() for c in key) for key in self.ids
        ):
            raise ConfigurationError("Alignment IDs cannot be empty or contain whitespace.")
        if len({len(value) for value in self.sequences}) != 1:
            raise ConfigurationError("Aligned sequences must have equal widths.")
        object.__setattr__(
            self, "sequences", tuple(sequence_text(value, gaps=True) for value in self.sequences)
        )

    @classmethod
    def from_sequences(cls, sequences: Mapping[str, Any]) -> MultipleAlignment:
        return cls(
            tuple(sequences), tuple(sequence_text(value, gaps=True) for value in sequences.values())
        )

    def distances(self, *, gap_policy: str = "pairwise_delete") -> Any:
        """p-distance, with explicit pairwise deletion or gaps-as-a-state convention."""
        np = require("numpy", "bio-analysis")
        if gap_policy not in {"pairwise_delete", "state"}:
            raise ConfigurationError("Unknown alignment gap policy.")
        result = np.zeros((len(self.ids), len(self.ids)))
        for left in range(len(self.ids)):
            for right in range(left):
                pairs = [
                    (a, b)
                    for a, b in zip(self.sequences[left], self.sequences[right], strict=True)
                    if (a in "ACGT" and b in "ACGT")
                    or (gap_policy == "state" and a in "ACGT-" and b in "ACGT-")
                ]
                if not pairs:
                    raise ConfigurationError("No comparable alignment columns for a sequence pair.")
                result[left, right] = result[right, left] = sum(a != b for a, b in pairs) / len(
                    pairs
                )
        return result


def _biopython_alignment(alignment: MultipleAlignment) -> Any:
    seq = require("Bio.Seq", "bio-formats")
    record = require("Bio.SeqRecord", "bio-formats")
    align = require("Bio.Align", "bio-formats")
    return align.MultipleSeqAlignment(
        [
            record.SeqRecord(seq.Seq(value), id=key, description="")
            for key, value in zip(alignment.ids, alignment.sequences, strict=True)
        ]
    )


def read_alignment(path: str | Path, *, format: str = "fasta") -> MultipleAlignment:
    """Read one alignment using Biopython codecs (FASTA, Clustal, Stockholm, PHYLIP…)."""
    result = require("Bio.AlignIO", "bio-formats").read(str(input_path(path)), format)
    return MultipleAlignment(
        tuple(record.id for record in result), tuple(str(record.seq).upper() for record in result)
    )


def write_alignment(
    alignment: MultipleAlignment, path: str | Path, *, format: str = "fasta"
) -> Path:
    target = Path(path)
    with target.open("x", encoding="utf-8") as handle:
        require("Bio.AlignIO", "bio-formats").write(_biopython_alignment(alignment), handle, format)
    return target


def align_sequences(
    sequences: Mapping[str, Any],
    output_directory: str | Path,
    *,
    method: str = "mafft",
    program: str | None = None,
    config: ExecutionConfig | None = None,
) -> MultipleAlignment:
    """Run MAFFT --auto, MUSCLE 5, or Clustal Omega; retain input/output/log files."""
    if method not in {"mafft", "muscle", "clustalo"} or len(sequences) < 2:
        raise ConfigurationError("Provide two sequences and method mafft/muscle/clustalo.")
    config = config or ExecutionConfig()
    directory = workspace(output_directory)
    source, target = directory / "input.fasta", directory / "alignment.fasta"
    ids = tuple(sequences)
    if len(set(ids)) != len(ids) or any(not key or any(c.isspace() for c in key) for key in ids):
        raise ConfigurationError("Provide unique FASTA-safe IDs.")
    # Internal IDs prevent aligners from truncating or rewriting caller IDs.
    source.write_text(
        "".join(
            f">seq_{index}\n{sequence_text(value)}\n"
            for index, value in enumerate(sequences.values())
        ),
        encoding="utf-8",
    )
    executable = program or method
    if method == "mafft":
        log = run_to_file(
            executable, ("--auto", "--thread", str(config.threads), str(source)), target, config
        )
    else:
        args = (
            ("-align", str(source), "-output", str(target), "-threads", str(config.threads))
            if method == "muscle"
            else (
                "-i",
                str(source),
                "-o",
                str(target),
                "--outfmt=fasta",
                "--threads",
                str(config.threads),
            )
        )
        log = run(executable, args, directory, config, outputs=(target,))
    (directory / "alignment.log").write_text(log, encoding="utf-8")
    result = read_alignment(target)
    rows = dict(zip(result.ids, result.sequences, strict=True))
    expected = {f"seq_{index}" for index in range(len(ids))}
    if set(rows) != expected:
        raise BackendExecutionError("Aligner changed the sequence ID set.")
    for index, original in enumerate(sequences.values()):
        if rows[f"seq_{index}"].replace("-", "") != sequence_text(original):
            raise BackendExecutionError("Aligner changed sequence content.")
    return MultipleAlignment(ids, tuple(rows[f"seq_{index}"] for index in range(len(ids))))


def parse_search_results(
    path: str | Path, *, format: str = "blast-xml", **parameters: Any
) -> Iterator[Any]:
    """Stream Bio.SearchIO QueryResult objects, including BLAST and HMMER formats."""
    yield from require("Bio.SearchIO", "bio-formats").parse(
        str(input_path(path)), format, **parameters
    )


def build_blast_database(
    fasta: str | Path,
    output_directory: str | Path,
    *,
    program: str = "makeblastdb",
    config: ExecutionConfig | None = None,
) -> ToolResult:
    source = input_path(fasta)
    directory = workspace(output_directory)
    log = run(
        program,
        (
            "-in",
            str(source),
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-out",
            str(directory / "database"),
        ),
        directory,
        config or ExecutionConfig(),
    )
    return ToolResult("makeblastdb", directory, artifacts(directory), log)


def local_blast(
    query: str | Path,
    database: str | Path,
    output_directory: str | Path,
    *,
    evalue: float = 1e-5,
    program: str = "blastn",
    config: ExecutionConfig | None = None,
) -> ToolResult:
    import math

    if not math.isfinite(evalue) or evalue <= 0:
        raise ConfigurationError("evalue must be positive and finite.")
    source = input_path(query)
    database = Path(database).expanduser().resolve()
    if not database.parent.is_dir():
        raise ConfigurationError("BLAST database directory does not exist.")
    directory = workspace(output_directory)
    config = config or ExecutionConfig()
    target = directory / "results.xml"
    log = run(
        program,
        (
            "-query",
            str(source),
            "-db",
            str(database),
            "-outfmt",
            "5",
            "-out",
            str(target),
            "-evalue",
            str(evalue),
            "-num_threads",
            str(config.threads),
        ),
        directory,
        config,
        outputs=(target,),
    )
    return ToolResult("blastn", directory, artifacts(directory), log)


def build_profile(
    alignment: str | Path,
    output_directory: str | Path,
    *,
    alphabet: str = "dna",
    program: str = "hmmbuild",
    config: ExecutionConfig | None = None,
) -> ToolResult:
    if alphabet not in {"dna", "rna", "amino"}:
        raise ConfigurationError("Unknown HMM alphabet.")
    source = input_path(alignment)
    directory = workspace(output_directory)
    target = directory / "profile.hmm"
    config = config or ExecutionConfig()
    log = run(
        program,
        ("--" + alphabet, "--cpu", str(config.threads), str(target), str(source)),
        directory,
        config,
        outputs=(target,),
    )
    return ToolResult("hmmbuild", directory, artifacts(directory), log)


def search_profile(
    profile: str | Path,
    sequences: str | Path,
    output_directory: str | Path,
    *,
    method: str = "nhmmer",
    program: str | None = None,
    config: ExecutionConfig | None = None,
) -> ToolResult:
    """nhmmer for nucleotide profiles; hmmsearch for amino-acid profiles."""
    if method not in {"nhmmer", "hmmsearch"}:
        raise ConfigurationError("method must be nhmmer or hmmsearch.")
    hmm, fasta = input_path(profile), input_path(sequences)
    directory = workspace(output_directory)
    target, report = directory / "hits.tbl", directory / "report.txt"
    config = config or ExecutionConfig()
    log = run(
        program or method,
        (
            "--cpu",
            str(config.threads),
            "--tblout",
            str(target),
            "-o",
            str(report),
            str(hmm),
            str(fasta),
        ),
        directory,
        config,
        outputs=(target, report),
    )
    return ToolResult(method, directory, artifacts(directory), log)


def parse_nhmmer_table(path: str | Path) -> Iterator[dict[str, Any]]:
    """Native nhmmer --tblout parser; convert aligned/envelope spans to 0-based half-open.

    Reverse-strand spans are ordered low-to-high and carry an explicit strand.
    Protein HMMER tables instead use parse_search_results(format='hmmer3-tab').
    """
    with input_path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split(maxsplit=15)
            try:
                start, end = sorted((int(fields[6]), int(fields[7])))
                envelope_start, envelope_end = sorted((int(fields[8]), int(fields[9])))
                if fields[11] not in {"+", "-"} or start < 1 or envelope_start < 1:
                    raise ValueError("Invalid nhmmer coordinates/strand.")
                yield {
                    "target": fields[0],
                    "query": fields[2],
                    "start": start - 1,
                    "end": end,
                    "envelope_start": envelope_start - 1,
                    "envelope_end": envelope_end,
                    "strand": fields[11],
                    "profile_start": int(fields[4]) - 1,
                    "profile_end": int(fields[5]),
                    "target_length": int(fields[10]),
                    "evalue": float(fields[12]),
                    "score": float(fields[13]),
                    "bias": float(fields[14]),
                    "description": fields[15] if len(fields) > 15 else "",
                }
            except (IndexError, ValueError) as exc:
                raise InputFormatError("Malformed nhmmer table row.") from exc
