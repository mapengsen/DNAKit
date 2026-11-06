"""Reference-genome alignment summaries from FASTA or minimap2/wfmash PAF.

Besides the matched-base consistency rate, each report carries the alignment
scores of Huang et al. (bioRxiv 10.1101/2024.09.18.612131): BLAST identity
(matches / alignment columns), gap-compressed identity (consecutive gaps count
as one difference) and the fraction of aligned queries.
"""

from __future__ import annotations

import hashlib
import math
import re
import shutil
import tempfile
from collections.abc import Iterable, Mapping
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
_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
_ALIGNERS = frozenset({"minimap2", "wfmash"})
_WFMASH_RESERVED = frozenset(
    {"-a", "--sam", "-m", "--approx-mapping", "-t", "--threads", "-i", "--align-paf", "-h"}
)


@dataclass(frozen=True, slots=True)
class ReferenceConsistencyConfig:
    """Aligner choice, aligner settings and resource limits for reference consistency.

    ``aligner`` selects minimap2 (default; suitable for short sequences) or
    wfmash (the aligner used by Huang et al.; intended for kilobase-scale and
    longer sequences). ``wfmash_args`` are extra wfmash options passed
    verbatim, e.g. ``("-w", "1000", "-S", "1000")``; wfmash options differ
    between releases, so check ``wfmash -h`` for the installed version.
    ``min_query_coverage`` is the aligned query fraction a query's best target
    needs for the query to count as aligned in ``aligned_fraction`` and the
    weighted identities; ``0.0`` counts every reported primary alignment.
    """

    max_paf_bytes: int = 200_000_000
    max_alignments: int = 1_000_000
    aligner: str = "minimap2"
    minimap2_executable: str = "minimap2"
    wfmash_executable: str = "wfmash"
    preset: str = "sr"
    wfmash_args: tuple[str, ...] = ()
    min_query_coverage: float = 0.0
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
        if self.aligner not in _ALIGNERS:
            raise ConfigurationError("aligner must be 'minimap2' or 'wfmash'.")
        for name in ("minimap2_executable", "wfmash_executable"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip() or "\x00" in value:
                raise ConfigurationError(f"{name} must be non-empty text.")
        if not isinstance(self.wfmash_args, tuple) or any(
            not isinstance(argument, str) or not argument or "\x00" in argument
            for argument in self.wfmash_args
        ):
            raise ConfigurationError("wfmash_args must be a tuple of non-empty strings.")
        if _WFMASH_RESERVED.intersection(self.wfmash_args):
            raise ConfigurationError(
                "wfmash_args must not change the PAF output, threads, or alignment mode."
            )
        if (
            isinstance(self.min_query_coverage, bool)
            or not isinstance(self.min_query_coverage, (int, float))
            or not 0.0 <= self.min_query_coverage <= 1.0
        ):
            raise ConfigurationError("min_query_coverage must be in [0, 1].")
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
    query_start: int
    query_end: int
    matches: int
    columns: int
    gap_compressed_columns: int
    mapq: int
    match_intervals: tuple[tuple[int, int], ...] | None


def _cigar_scores(
    cigar: str, strand: str, query_start: int, query_end: int, paf_matches: int
) -> tuple[int, int, int, tuple[tuple[int, int], ...] | None]:
    """Return matches, alignment columns, gap-compressed columns and '=' query intervals.

    Columns are M/=/X/I/D bases (BLAST identity denominator); gap-compressed
    columns count each I/D run once. With =/X operations the matches and
    matched query intervals come from the CIGAR (wfmash); with M operations the
    PAF match count is used and matched positions are unknown (minimap2).
    """

    operations = [(int(count), op) for count, op in _CIGAR.findall(cigar)]
    if "".join(f"{count}{op}" for count, op in operations) != cigar:
        raise ValueError("Malformed cg:Z CIGAR")
    pairs = sum(count for count, op in operations if op in "M=X")
    gaps = sum(count for count, op in operations if op in "ID")
    gap_runs = sum(1 for count, op in operations if op in "ID" and count)
    query_span = sum(count for count, op in operations if op in "M=XI")
    if pairs == 0 or query_span != query_end - query_start:
        raise ValueError("CIGAR does not span the PAF query interval")
    if not any(op in "=X" for _, op in operations):
        if paf_matches > pairs:
            raise ValueError("PAF match count exceeds CIGAR aligned pairs")
        return paf_matches, pairs + gaps, pairs + gap_runs, None
    intervals: list[tuple[int, int]] = []
    position = query_start if strand == "+" else query_end
    for count, op in operations:
        if op not in "M=XI":
            continue
        if strand == "+":
            interval = (position, position + count)
            position += count
        else:
            interval = (position - count, position)
            position -= count
        if op == "=":
            intervals.append(interval)
        elif op == "M":
            raise ValueError("CIGAR mixes M with =/X operations")
    matches = sum(stop - begin for begin, stop in intervals)
    return matches, pairs + gaps, pairs + gap_runs, tuple(intervals)


def _parse_paf_line(
    line: str, line_number: int, lengths: dict[str, int]
) -> tuple[str, _PafHit | None]:
    columns = line.rstrip("\r\n").split("\t")
    try:
        if len(columns) < 13:
            raise ValueError("Expected 12 PAF columns plus a cg:Z tag")
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
        if len(types) > 1 or (types and types[0] not in {"P", "I", "S", "i"}):
            raise ValueError("Expected at most one valid tp:A tag")
        cigars = [tag[5:] for tag in columns[12:] if tag.startswith("cg:Z:")]
        if len(cigars) != 1 or not cigars[0]:
            raise ValueError("Expected a cg:Z tag from base-level alignment")
        scores = _cigar_scores(cigars[0], strand, qs, qe, matches)
    except ValueError as exc:
        raise ConfigurationError(
            f"Invalid PAF at line {line_number}: {exc}.",
            code="INVALID_REFERENCE_CONSISTENCY_PAF",
        ) from exc
    if types and types[0] in {"S", "i"}:
        return query, None
    return query, _PafHit(target, qs, qe, *scores[:3], mapq, scores[3])


def _union_length(intervals: Iterable[tuple[int, int]]) -> int:
    total = 0
    reach = -1
    for begin, stop in sorted(intervals):
        if stop <= reach:
            continue
        total += stop - max(begin, reach)
        reach = stop
    return total


@dataclass(frozen=True, slots=True)
class _QueryAlignment:
    """All primary records of one query against its best-matching target."""

    target: str
    consistency_rate: float
    query_coverage: float
    matches: int
    columns: int
    gap_compressed_columns: int
    mapq: int

    @property
    def blast_identity(self) -> float:
        return self.matches / self.columns

    @property
    def gap_compressed_identity(self) -> float:
        return self.matches / self.gap_compressed_columns

    @classmethod
    def pool(cls, target: str, hits: list[_PafHit], query_length: int) -> _QueryAlignment:
        if all(hit.match_intervals is not None for hit in hits):
            matched = _union_length(
                interval for hit in hits for interval in hit.match_intervals or ()
            )
        else:
            matched = min(query_length, sum(hit.matches for hit in hits))
        return cls(
            target,
            matched / query_length,
            _union_length((hit.query_start, hit.query_end) for hit in hits) / query_length,
            sum(hit.matches for hit in hits),
            sum(hit.columns for hit in hits),
            sum(hit.gap_compressed_columns for hit in hits),
            max(hit.mapq for hit in hits),
        )

    @property
    def rank(self) -> tuple[float, float, float, int]:
        return (self.consistency_rate, self.blast_identity, self.query_coverage, self.mapq)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1_048_576):
            digest.update(block)
    return digest.hexdigest()


def _write_fai(fasta: Path, destination: Path) -> None:
    """Write a samtools-compatible .fai index for an uncompressed FASTA."""

    records: list[str] = []
    name: str | None = None
    length = offset = line_bases = line_width = 0
    short_line_seen = False
    position = 0

    def finish() -> None:
        if name is not None:
            records.append(f"{name}\t{length}\t{offset}\t{line_bases}\t{line_width}\n")

    with fasta.open("rb") as stream:
        for raw in stream:
            if raw.startswith(b">"):
                finish()
                fields = raw[1:].split()
                if not fields:
                    raise ConfigurationError(f"FASTA record without a name in {fasta}.")
                name = fields[0].decode("utf-8", errors="strict")
                length = line_bases = line_width = 0
                short_line_seen = False
                offset = position + len(raw)
            elif name is not None:
                bases = len(raw.rstrip(b"\r\n"))
                if bases:
                    if short_line_seen or (line_bases and bases > line_bases):
                        raise ConfigurationError(
                            f"FASTA {fasta} has uneven line lengths; index it with samtools faidx.",
                            code="REFERENCE_FASTA_NOT_INDEXABLE",
                        )
                    if not line_bases:
                        line_bases, line_width = bases, len(raw)
                    elif bases < line_bases:
                        short_line_seen = True
                    length += bases
            position += len(raw)
    finish()
    if not records:
        raise ConfigurationError(f"FASTA {fasta} contains no records.")
    temporary = destination.with_name(destination.name + ".tmp")
    temporary.write_text("".join(records), encoding="utf-8")
    temporary.replace(destination)


def _indexed_reference(fasta: Path, work: Path) -> Path:
    """Return a reference path whose .fai index wfmash can open.

    An existing ``<fasta>.fai`` is used as is; otherwise one is written next to
    the FASTA (cached for later runs) or, when that directory is read-only,
    beside a symlink inside the temporary work directory.
    """

    index = fasta.with_name(fasta.name + ".fai")
    if index.is_file() and index.stat().st_mtime >= fasta.stat().st_mtime:
        return fasta
    if fasta.suffix in {".gz", ".bgz"}:
        raise ConfigurationError(
            "wfmash needs an index for compressed FASTA; run samtools faidx on the bgzip file.",
            code="REFERENCE_FASTA_NOT_INDEXABLE",
        )
    try:
        _write_fai(fasta, index)
        return fasta
    except OSError:
        link = work / "reference.fa"
        link.symlink_to(fasta)
        _write_fai(fasta, work / "reference.fa.fai")
        return link


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

    wfmash = config.aligner == "wfmash"
    executable = shutil.which(config.wfmash_executable if wfmash else config.minimap2_executable)
    if executable is None:
        raise BackendUnavailableError(
            f"{config.aligner} is required for direct reference consistency evaluation.",
            code="WFMASH_UNAVAILABLE" if wfmash else "MINIMAP2_UNAVAILABLE",
            hint=f"Install {config.aligner} or set config.{config.aligner}_executable.",
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
        if wfmash:
            _write_fai(query_path, query_path.with_name(query_path.name + ".fai"))
            target = _indexed_reference(fasta, work)
            command = execute_bounded_command(
                executable,
                (
                    str(target),
                    str(query_path),
                    "-t",
                    str(config.threads),
                    "--quiet",
                    *config.wfmash_args,
                ),
                backend_id="wfmash",
                cwd=work,
                timeout_seconds=config.timeout_seconds,
                max_output_bytes=1_000_000,
                stdout_path=paf_path,
                max_monitored_output_bytes=config.max_paf_bytes,
            )
        else:
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
                f"{config.aligner} did not produce a PAF alignment file.",
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
            "alignment_source": f"{config.aligner}-generated",
            "aligner": config.aligner,
            "paf_path": None,
            "reference_fasta_path": str(fasta),
            "threads": config.threads,
            **({"wfmash_args": config.wfmash_args} if wfmash else {"preset": config.preset}),
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
    """Score each query's alignment to a reference and summarize the set.

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

    hits: dict[str, dict[str, list[_PafHit]]] = {}
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
            hits.setdefault(query, {}).setdefault(hit.target, []).append(hit)

    best: dict[str, _QueryAlignment] = {}
    for query, by_target in hits.items():
        for target, target_hits in by_target.items():
            candidate = _QueryAlignment.pool(target, target_hits, lengths[query])
            if query not in best or candidate.rank > best[query].rank:
                best[query] = candidate
    scores = tuple(
        alignment.consistency_rate if (alignment := best.get(item.subject_id)) else 0.0
        for item in items
    )
    aligned = {
        query: alignment
        for query, alignment in best.items()
        if alignment.query_coverage >= resolved.min_query_coverage
    }

    def _entry_metrics(subject_id: str, score: float) -> FrozenDict:
        alignment = best.get(subject_id)
        return FrozenDict(
            {
                "consistency_rate": score,
                "aligned": subject_id in aligned,
                "blast_identity": None if alignment is None else alignment.blast_identity,
                "gap_compressed_identity": (
                    None if alignment is None else alignment.gap_compressed_identity
                ),
                "query_coverage": 0.0 if alignment is None else alignment.query_coverage,
                "best_target": None if alignment is None else alignment.target,
            }
        )

    result_entries = tuple(
        EvaluationEntry(item.subject_id, index, _entry_metrics(item.subject_id, scores[index]))
        for index, item in enumerate(items)
    )
    pooled = tuple(aligned.values())
    matched_total = sum(alignment.matches for alignment in pooled)
    metrics = {
        "consistency_rate": mean(scores),
        "query_count": len(items),
        "aligned_count": len(aligned),
        "aligned_fraction": len(aligned) / len(items) if items else None,
        "weighted_blast_identity": (
            matched_total / sum(alignment.columns for alignment in pooled) if pooled else None
        ),
        "weighted_gap_compressed_identity": (
            matched_total / sum(alignment.gap_compressed_columns for alignment in pooled)
            if pooled
            else None
        ),
    }
    return EvaluationReport(
        "reference_consistency",
        "primary-paf-best-hit",
        "eval-reference-consistency-v3",
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
                            "pool each query's primary records per target; keep the target "
                            "with max consistency_rate, BLAST identity, coverage, MAPQ; "
                            "first on ties"
                        ),
                        "consistency_rate_definition": (
                            "mean(matched query bases / query_length for each query; "
                            "no primary hit = 0)"
                        ),
                        "blast_identity_definition": (
                            "matches / alignment columns (CIGAR M/=/X/I/D bases)"
                        ),
                        "gap_compressed_identity_definition": (
                            "matches / (CIGAR M/=/X bases + number of I/D runs)"
                        ),
                        "aligned_definition": ("best-target query_coverage >= min_query_coverage"),
                        "weighting": (
                            "weighted identities pool matches and columns over aligned queries"
                        ),
                        "min_query_coverage": resolved.min_query_coverage,
                        "metric_reference": "Huang et al., bioRxiv 10.1101/2024.09.18.612131",
                        "limits": resolved.limits,
                    }
                ),
            )
        ),
        FrozenDict(metrics),
        result_entries,
        Provenance(
            reference=ReferenceInfo(
                reference_name, version=reference_version, checksum=reference_sha256
            )
        ),
    )


__all__ = ["ReferenceConsistencyConfig", "evaluate_reference_consistency"]

_register_public_verbose(globals())
