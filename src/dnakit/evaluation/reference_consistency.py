"""Reference-genome alignment summaries from FASTA or minimap2 PAF."""

from __future__ import annotations

import hashlib
import math
import re
import shutil
import tempfile
from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import cast

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.backends import execute_bounded_command
from dnakit.core import DNARecord, Provenance, ReferenceInfo, Topology
from dnakit.core._json import FrozenDict, to_json_compatible
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError
from dnakit.references import download_genome, supported_genome_aliases

from ._shared import EvaluationInput, materialize_input, mean, require_pairwise_compatible
from .config import EvaluationLimits
from .results import EvaluationEntry, EvaluationReport

_SHA256 = re.compile(r"[0-9a-f]{64}\Z")


@dataclass(frozen=True, slots=True)
class ReferenceConsistencyConfig:
    """Resource limits and minimap2 settings for reference consistency."""

    max_paf_bytes: int = 200_000_000
    max_alignments: int = 1_000_000
    minimap2_executable: str = "minimap2"
    preset: str = "sr"
    threads: int = 1
    timeout_seconds: float = 3600.0
    limits: EvaluationLimits = field(default_factory=EvaluationLimits)

    def __post_init__(self) -> None:
        for name, maximum in (
            ("max_paf_bytes", 10_000_000_000),
            ("max_alignments", 100_000_000),
            ("threads", 256),
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= maximum:
                raise ConfigurationError(f"{name} must be an integer in [1, {maximum}].")
        if (
            not isinstance(self.minimap2_executable, str)
            or not self.minimap2_executable.strip()
            or "\x00" in self.minimap2_executable
        ):
            raise ConfigurationError("minimap2_executable must be non-empty text.")
        if self.preset not in {"sr", "asm5", "asm10", "asm20", "map-hifi", "map-ont"}:
            raise ConfigurationError("Unsupported minimap2 preset.")
        if (
            isinstance(self.timeout_seconds, bool)
            or not isinstance(self.timeout_seconds, (int, float))
            or not math.isfinite(self.timeout_seconds)
            or not 0 < self.timeout_seconds <= 86_400
        ):
            raise ConfigurationError("timeout_seconds must be in (0, 86400].")
        if not isinstance(self.limits, EvaluationLimits):
            raise ConfigurationError("limits must be EvaluationLimits.")


@dataclass(frozen=True, slots=True)
class _PafHit:
    target: str
    strand: str
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    matches: int
    block_length: int
    mapq: int
    query_length: int

    @property
    def identity(self) -> float:
        return self.matches / self.block_length

    @property
    def query_coverage(self) -> float:
        return (self.query_end - self.query_start) / self.query_length

    @property
    def matched_fraction(self) -> float:
        return self.matches / self.query_length


def _parse_paf_line(
    line: str, line_number: int, lengths: dict[str, int]
) -> tuple[str, _PafHit | None]:
    columns = line.rstrip("\r\n").split("\t")
    try:
        if len(columns) < 14:
            raise ValueError("Expected PAF columns with tp:A and cg:Z tags")
        query = columns[0]
        if query not in lengths:
            raise ValueError("Query ID is absent from the input sequences")
        qlen, qs, qe = map(int, columns[1:4])
        strand, target = columns[4:6]
        tlen, ts, te, matches, block_length, mapq = map(int, columns[6:12])
        if (
            qlen != lengths[query]
            or qlen == 0
            or not 0 <= qs < qe <= qlen
            or not target
            or tlen <= 0
            or not 0 <= ts < te <= tlen
            or strand not in {"+", "-"}
            or block_length < max(qe - qs, te - ts)
            or not 0 <= matches <= min(qe - qs, te - ts, block_length)
            or not 0 <= mapq <= 255
        ):
            raise ValueError("Invalid PAF lengths, coordinates, strand, or match count")
        types = [tag[5:] for tag in columns[12:] if tag.startswith("tp:A:")]
        if len(types) != 1 or types[0] not in {"P", "I", "S", "i"}:
            raise ValueError("Expected exactly one valid minimap2 tp:A tag")
        cigars = [tag[5:] for tag in columns[12:] if tag.startswith("cg:Z:")]
        if len(cigars) != 1 or not cigars[0]:
            raise ValueError("Expected a cg:Z tag from base-level alignment (-c)")
    except ValueError as exc:
        raise ConfigurationError(
            f"Invalid PAF at line {line_number}: {exc}.",
            code="INVALID_REFERENCE_CONSISTENCY_PAF",
        ) from exc
    if types[0] in {"S", "i"}:
        return query, None
    return query, _PafHit(target, strand, qs, qe, ts, te, matches, block_length, mapq, qlen)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1_048_576):
            digest.update(block)
    return digest.hexdigest()


def _cached_hg38(directory: Path) -> tuple[Path, str] | None:
    """Reuse only an NCBI download whose FASTA matches its sidecar MD5."""

    accession = supported_genome_aliases()["hg38"]
    candidates = tuple(directory.glob(f"{accession}*_genomic.fna"))
    if not candidates:
        return None
    if len(candidates) != 1:
        raise ConfigurationError("Expected one cached hg38 genomic FASTA.")
    fasta = candidates[0]
    checksum = directory / f"{accession}_md5sum.txt"
    if not fasta.is_file() or not checksum.is_file():
        raise ConfigurationError("Cached hg38 FASTA is missing its NCBI checksum file.")
    expected: str | None = None
    for line in checksum.read_text(encoding="utf-8").splitlines():
        fields = line.split(maxsplit=1)
        if len(fields) != 2:
            continue
        candidate = Path(fields[1].lstrip("*")).name
        if candidate == fasta.name and re.fullmatch(r"[0-9a-fA-F]{32}", fields[0]):
            expected = fields[0].lower()
            break
    if expected is None:
        raise ConfigurationError("Cached hg38 checksum file has no matching FASTA entry.")
    md5 = hashlib.md5(usedforsecurity=False)
    sha256 = hashlib.sha256()
    with fasta.open("rb") as stream:
        while block := stream.read(1_048_576):
            md5.update(block)
            sha256.update(block)
    if md5.hexdigest() != expected:
        raise ConfigurationError("Cached hg38 FASTA failed its NCBI MD5 check.")
    return fasta, sha256.hexdigest()


def _evaluate_from_fasta(
    queries: EvaluationInput,
    *,
    reference_fasta: str | Path | None,
    reference_dir: str | Path,
    reference_name: str | None,
    reference_version: str | None,
    reference_sha256: str | None,
    config: ReferenceConsistencyConfig,
) -> EvaluationReport:
    items = materialize_input(queries, limits=config.limits)
    if not items:
        raise ConfigurationError("Direct reference consistency requires at least one query.")
    ids = [item.subject_id for item in items]
    if len(set(ids)) != len(ids):
        raise ConfigurationError(
            "Query IDs must be unique for PAF matching.", code="DUPLICATE_REFERENCE_CONSISTENCY_ID"
        )
    for item in items:
        if (
            not item.sequence.symbol_length
            or any(character.isspace() for character in item.subject_id)
            or ">" in item.subject_id
        ):
            raise ConfigurationError("Direct alignment requires nonempty DNA and FASTA-safe IDs.")
        require_pairwise_compatible(item.sequence, role="query")
        if item.sequence.topology is Topology.CIRCULAR:
            raise ConfigurationError("Circular queries must be linearized before alignment.")

    executable = shutil.which(config.minimap2_executable)
    if executable is None:
        raise BackendUnavailableError(
            "minimap2 is required for direct reference consistency evaluation.",
            code="MINIMAP2_UNAVAILABLE",
            hint="Install minimap2 or set config.minimap2_executable.",
        )

    if reference_fasta is None:
        directory = Path(reference_dir).expanduser().resolve()
        cached = _cached_hg38(directory) if directory.is_dir() else None
        if cached is None:
            downloaded = download_genome("hg38", directory, show_progress=True)
            fasta = Path(downloaded.fasta_path).resolve(strict=True)
            digest = downloaded.fasta_sha256
            version = reference_version or downloaded.accession
        else:
            fasta, digest = cached
            version = reference_version or supported_genome_aliases()["hg38"]
        name = reference_name or "human"
    else:
        fasta = Path(reference_fasta).expanduser().resolve(strict=True)
        if not fasta.is_file():
            raise ConfigurationError("reference_fasta must be a file.")
        digest = _sha256_file(fasta)
        name = reference_name or "local_reference"
        version = reference_version or f"local:{digest[:16]}"
    if reference_sha256 is not None and reference_sha256 != digest:
        raise ConfigurationError("reference_sha256 does not match the selected FASTA.")
    with tempfile.TemporaryDirectory(prefix="dnakit-reference-consistency-") as raw_work:
        work = Path(raw_work)
        query_path = work / "queries.fasta"
        with query_path.open("w", encoding="ascii") as stream:
            for item in items:
                stream.write(f">{item.subject_id}\n{item.sequence.symbols}\n")
        paf_path = work / "alignments.paf"
        command = execute_bounded_command(
            executable,
            (
                "-x",
                config.preset,
                "-c",
                "-t",
                str(config.threads),
                "-o",
                str(paf_path),
                str(fasta),
                str(query_path),
            ),
            backend_id="minimap2",
            cwd=work,
            timeout_seconds=config.timeout_seconds,
            max_output_bytes=1_000_000,
            monitored_output_paths=(paf_path,),
            max_monitored_output_bytes=config.max_paf_bytes,
        )
        if command.return_code or not paf_path.is_file():
            raise BackendExecutionError(
                "minimap2 did not produce a PAF alignment file.",
                code="REFERENCE_CONSISTENCY_MAPPING_FAILED",
                context={"return_code": command.return_code, "output": command.output[-2000:]},
            )
        result = evaluate_reference_consistency(
            tuple(item.record or DNARecord(item.sequence, item.subject_id) for item in items),
            paf_path,
            reference_name=name,
            reference_version=version,
            reference_sha256=digest,
            config=config,
        )
    parameters = dict(result.parameters)
    parameters.update(
        {
            "alignment_source": "minimap2-generated",
            "paf_path": None,
            "reference_fasta_path": str(fasta),
            "preset": config.preset,
            "threads": config.threads,
        }
    )
    return replace(result, parameters=FrozenDict(parameters))


def evaluate_reference_consistency(
    queries: EvaluationInput,
    paf_path: str | Path | None = None,
    *,
    reference_fasta: str | Path | None = None,
    reference_dir: str | Path = "references/hg38",
    reference_name: str | None = None,
    reference_version: str | None = None,
    reference_sha256: str | None = None,
    config: ReferenceConsistencyConfig | None = None,
    verbose: bool = True,
) -> EvaluationReport:
    """Average each query's matched-base fraction over every input query.

    Without PAF, align to a local FASTA or download and cache hg38. A supplied
    PAF is parsed directly and must have explicit reference metadata. Queries
    without a primary alignment contribute zero to the mean.
    """

    resolved = ReferenceConsistencyConfig() if config is None else config
    if not isinstance(resolved, ReferenceConsistencyConfig):
        raise TypeError("config must be ReferenceConsistencyConfig or None.")
    if paf_path is None:
        return _evaluate_from_fasta(
            queries,
            reference_fasta=reference_fasta,
            reference_dir=reference_dir,
            reference_name=reference_name,
            reference_version=reference_version,
            reference_sha256=reference_sha256,
            config=resolved,
        )
    if reference_fasta is not None:
        raise ConfigurationError("Supply either paf_path or reference_fasta, not both.")
    if not isinstance(reference_name, str) or not reference_name.strip():
        raise ConfigurationError("reference_name must be non-empty text.")
    if not isinstance(reference_version, str) or not reference_version.strip():
        raise ConfigurationError("reference_version must be non-empty text.")
    if reference_sha256 is not None and (
        not isinstance(reference_sha256, str) or _SHA256.fullmatch(reference_sha256) is None
    ):
        raise ConfigurationError("reference_sha256 must be a lowercase SHA-256 digest.")
    items = materialize_input(queries, limits=resolved.limits)
    lengths = {item.subject_id: item.sequence.symbol_length for item in items}
    if len(lengths) != len(items):
        raise ConfigurationError(
            "Query IDs must be unique for PAF matching.", code="DUPLICATE_REFERENCE_CONSISTENCY_ID"
        )
    paf = Path(paf_path).expanduser().resolve(strict=True)
    if not paf.is_file():
        raise ConfigurationError("paf_path must be a file.")
    if paf.stat().st_size > resolved.max_paf_bytes:
        raise ConfigurationError(
            "PAF exceeds max_paf_bytes.", code="REFERENCE_CONSISTENCY_PAF_LIMIT"
        )

    best: dict[str, _PafHit] = {}
    digest = hashlib.sha256()
    byte_count = 0
    alignment_count = 0
    with paf.open("rb") as stream:
        for line_number, raw_line in enumerate(stream, 1):
            byte_count += len(raw_line)
            if byte_count > resolved.max_paf_bytes:
                raise ConfigurationError(
                    "PAF exceeds max_paf_bytes.", code="REFERENCE_CONSISTENCY_PAF_LIMIT"
                )
            digest.update(raw_line)
            if not raw_line.strip():
                continue
            alignment_count += 1
            if alignment_count > resolved.max_alignments:
                raise ConfigurationError(
                    "PAF exceeds max_alignments.", code="REFERENCE_CONSISTENCY_ALIGNMENT_LIMIT"
                )
            try:
                line = raw_line.decode("utf-8")
            except UnicodeDecodeError as exc:
                raise ConfigurationError(
                    f"Invalid UTF-8 PAF at line {line_number}.",
                    code="INVALID_REFERENCE_CONSISTENCY_PAF",
                ) from exc
            query, hit = _parse_paf_line(line, line_number, lengths)
            if hit is None:
                continue
            current = best.get(query)
            rank = (hit.matched_fraction, hit.identity, hit.query_coverage, hit.mapq)
            if current is None or rank > (
                current.matched_fraction,
                current.identity,
                current.query_coverage,
                current.mapq,
            ):
                best[query] = hit

    scores = tuple(
        hit.matched_fraction if (hit := best.get(item.subject_id)) is not None else 0.0
        for item in items
    )
    result_entries = tuple(
        EvaluationEntry(item.subject_id, index, FrozenDict({"consistency_rate": scores[index]}))
        for index, item in enumerate(items)
    )
    return EvaluationReport(
        "reference_consistency",
        "primary-minimap2-paf-best-hit",
        "eval-reference-consistency-v2",
        FrozenDict(
            cast(
                Mapping[str, object],
                to_json_compatible(
                    {
                        "reference_name": reference_name,
                        "reference_version": reference_version,
                        "reference_sha256": reference_sha256,
                        "paf_path": str(paf),
                        "paf_sha256": digest.hexdigest(),
                        "paf_alignment_count": alignment_count,
                        "best_hit_rule": (
                            "max matched_fraction, identity, coverage, MAPQ; first on ties"
                        ),
                        "consistency_rate_definition": (
                            "mean(PAF nmatch / query_length for each query; no primary hit = 0)"
                        ),
                        "limits": resolved.limits,
                    }
                ),
            )
        ),
        FrozenDict({"consistency_rate": mean(scores)}),
        result_entries,
        Provenance(
            reference=ReferenceInfo(
                reference_name, version=reference_version, checksum=reference_sha256
            )
        ),
    )


__all__ = ["ReferenceConsistencyConfig", "evaluate_reference_consistency"]

_register_public_verbose(globals())
