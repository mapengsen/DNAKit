"""Native FASTX processing, paired reads, interval extraction and alignment summaries."""

from __future__ import annotations

import gzip
import hashlib
import itertools
import random
import re
from collections import Counter, deque
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any, TextIO, cast

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, positive, require, sequence_text
from dnakit.exceptions import ConfigurationError, InputFormatError


@dataclass(frozen=True)
class Read:
    id: str
    sequence: str
    qualities: tuple[int, ...] | None = None
    description: str = ""

    def __post_init__(self) -> None:
        if not self.id or any(c.isspace() for c in self.id):
            raise ConfigurationError("Read IDs must be non-empty without whitespace.")
        if any(c in self.description for c in "\r\n"):
            raise ConfigurationError("Read descriptions cannot contain newlines.")
        object.__setattr__(self, "sequence", sequence_text(self.sequence))
        if self.qualities is not None and (
            len(self.qualities) != len(self.sequence)
            or any(type(q) is not int or not 0 <= q <= 93 for q in self.qualities)
        ):
            raise ConfigurationError(
                "Phred qualities must match sequence length and lie in [0,93]."
            )

    @property
    def pair_id(self) -> str:
        return self.id[:-2] if self.id.endswith(("/1", "/2")) else self.id


def _open(path: Path, mode: str) -> TextIO:
    if path.suffix.lower() == ".gz":
        return cast(TextIO, gzip.open(path, mode + "t", encoding="ascii", newline=""))
    return cast(TextIO, path.open(mode, encoding="ascii", newline=""))


def read_fastx(
    path: str | Path,
    *,
    format: str = "fastq",
    phred_offset: int = 33,
    max_read_length: int = 10_000_000,
    verbose: bool = True,
) -> Iterator[Read]:
    """Stream FASTA or multiline FASTQ; quality '@' and '+' remain quality symbols."""
    if format not in {"fasta", "fastq"} or phred_offset not in {33, 64}:
        raise ConfigurationError("Use fasta/fastq and Phred offset 33/64.")
    positive(max_read_length, "max_read_length")
    with _open(input_path(path), "r") as handle:
        if format == "fasta":
            title: str | None = None
            parts: list[str] = []
            length = 0
            for line in handle:
                line = line.rstrip("\r\n")
                if line.startswith(">"):
                    if title is not None:
                        yield _read_from_title(title, "".join(parts))
                    title, parts, length = line[1:], [], 0
                elif line:
                    if title is None:
                        raise InputFormatError("FASTA sequence precedes its header.")
                    parts.append(line)
                    length += len(line)
                    if length > max_read_length:
                        raise InputFormatError("FASTA record exceeds max_read_length.")
            if title is not None:
                yield _read_from_title(title, "".join(parts))
            return
        while title_line := handle.readline():
            if not title_line.startswith("@"):
                raise InputFormatError("Expected a FASTQ @ header.")
            title = title_line[1:].rstrip("\r\n")
            parts, length = [], 0
            for line in handle:
                line = line.rstrip("\r\n")
                if line.startswith("+"):
                    if line[1:] and line[1:] != title:
                        raise InputFormatError("FASTQ repeated + title does not match @ title.")
                    break
                parts.append(line)
                length += len(line)
                if length > max_read_length:
                    raise InputFormatError("FASTQ sequence exceeds max_read_length.")
            else:
                raise InputFormatError("Truncated FASTQ sequence.")
            quality = ""
            while len(quality) < length:
                line = handle.readline()
                if not line:
                    raise InputFormatError("Truncated FASTQ quality.")
                quality += line.rstrip("\r\n")
            if len(quality) != length:
                raise InputFormatError("FASTQ quality length differs from sequence length.")
            yield _read_from_title(
                title, "".join(parts), tuple(ord(c) - phred_offset for c in quality)
            )


def _read_from_title(title: str, sequence: str, quality: tuple[int, ...] | None = None) -> Read:
    pieces = title.split(maxsplit=1)
    if not pieces:
        raise InputFormatError("Read header is empty.")
    return Read(pieces[0], sequence.upper(), quality, pieces[1] if len(pieces) > 1 else "")


def write_fastx(
    reads: Iterable[Read],
    path: str | Path,
    *,
    format: str = "fastq",
    phred_offset: int = 33,
    default_quality: int | None = None,
    verbose: bool = True,
) -> int:
    """Write fresh FASTA/FASTQ with explicit Phred conversion or synthetic quality."""
    if format not in {"fasta", "fastq"} or phred_offset not in {33, 64}:
        raise ConfigurationError("Use fasta/fastq and Phred offset 33/64.")
    if default_quality is not None and (
        type(default_quality) is not int or not 0 <= default_quality <= 93
    ):
        raise ConfigurationError("default_quality must be in [0,93].")
    count = 0
    with _open(Path(path), "x") as handle:
        for record in reads:
            title = record.id + (" " + record.description if record.description else "")
            if format == "fasta":
                handle.write(f">{title}\n{record.sequence}\n")
            else:
                qualities = record.qualities
                if qualities is None:
                    if default_quality is None:
                        raise ConfigurationError(
                            "FASTA to FASTQ conversion requires explicit default_quality."
                        )
                    qualities = (default_quality,) * len(record.sequence)
                if any(q + phred_offset > 126 for q in qualities):
                    raise ConfigurationError(
                        "Quality cannot be represented in the requested encoding."
                    )
                handle.write(
                    f"@{title}\n{record.sequence}\n+\n"
                    + "".join(chr(q + phred_offset) for q in qualities)
                    + "\n"
                )
            count += 1
    return count


@dataclass(frozen=True)
class ReadPair:
    id: str
    first: Read | None
    second: Read | None


def pair_reads(
    first: Iterable[Read],
    second: Iterable[Read],
    *,
    max_pending: int = 100_000,
    include_orphans: bool = True,
    verbose: bool = True,
) -> Iterator[ReadPair]:
    """Pair interleaved input streams by /1,/2-normalized ID with bounded pending storage.

    Output follows match completion order. Duplicate IDs are rejected per mate
    stream. The seen-ID sets grow with record count; pending read payload is bounded.
    """
    positive(max_pending, "max_pending")
    pending: tuple[dict[str, Read], dict[str, Read]] = ({}, {})
    seen: tuple[set[str], set[str]] = (set(), set())
    for items in itertools.zip_longest(first, second):
        for side, read in enumerate(items):
            if read is None:
                continue
            key = read.pair_id
            if key in seen[side]:
                raise InputFormatError(f"Duplicate read ID in mate {side + 1}: {key}")
            seen[side].add(key)
            other = pending[1 - side].pop(key, None)
            if other is not None:
                yield ReadPair(key, read if side == 0 else other, other if side == 0 else read)
            else:
                pending[side][key] = read
                if sum(len(values) for values in pending) > max_pending:
                    raise ConfigurationError(
                        "Unmatched reads exceeded max_pending; synchronize or raise the bound."
                    )
    if include_orphans:
        for side, values in enumerate(pending):
            for key, read in values.items():
                yield ReadPair(key, read if side == 0 else None, read if side == 1 else None)


@dataclass(frozen=True)
class RepairResult:
    output: Path
    recovered_reads: int
    discarded_lines: int


def repair_fastq(
    source: str | Path, destination: str | Path, *, phred_offset: int = 33, verbose: bool = True
) -> RepairResult:
    """Recover valid four-line records, dropping corrupt lines without inventing bases/qualities.

    This explicit recovery mode is for four-line FASTQ only. Multiline FASTQ
    must use the strict reader. Ambiguous synchronization cannot be guaranteed.
    """
    if phred_offset not in {33, 64}:
        raise ConfigurationError("Phred offset must be 33 or 64.")
    recovered = 0
    dropped = 0

    def records() -> Iterator[Read]:
        nonlocal recovered, dropped
        with _open(input_path(source), "r") as handle:
            buffer: deque[str] = deque()
            for line in handle:
                buffer.append(line.rstrip("\r\n"))
                while len(buffer) >= 4:
                    title, seq, plus, quality = list(itertools.islice(buffer, 4))
                    valid = (
                        title.startswith("@") and plus.startswith("+") and len(seq) == len(quality)
                    )
                    valid = valid and (plus == "+" or plus[1:] == title[1:])
                    try:
                        record = (
                            _read_from_title(
                                title[1:], seq, tuple(ord(c) - phred_offset for c in quality)
                            )
                            if valid
                            else None
                        )
                    except (ConfigurationError, InputFormatError):
                        record = None
                    if record is None:
                        buffer.popleft()
                        dropped += 1
                    else:
                        for _ in range(4):
                            buffer.popleft()
                        recovered += 1
                        yield record
            dropped += len(buffer)

    write_fastx(records(), destination, phred_offset=phred_offset)
    return RepairResult(Path(destination), recovered, dropped)


def arrange_reads(
    reads: Iterable[Read],
    *,
    operation: str = "sort",
    key: str = "id",
    seed: int = 0,
    max_records: int = 1_000_000,
    verbose: bool = True,
) -> tuple[Read, ...]:
    positive(max_records, "max_records")
    values = list(itertools.islice(reads, max_records + 1))
    if len(values) > max_records:
        raise ConfigurationError("Record collection exceeds max_records.")
    keys = {
        "id": lambda read: read.id,
        "sequence": lambda read: read.sequence,
        "length": lambda read: len(read.sequence),
        "name": lambda read: read.id + " " + read.description,
    }
    if key not in keys or operation not in {"sort", "shuffle"}:
        raise ConfigurationError("Unknown record order operation/key.")
    if operation == "sort":
        values.sort(key=keys[key])
    else:
        random.Random(seed).shuffle(values)
    return tuple(values)


def common_reads(
    collections: Iterable[Iterable[Read]], *, key: str = "id", verbose: bool = True
) -> tuple[Read, ...]:
    """Set intersection across collections, preserving first-collection records and order."""
    values = [tuple(collection) for collection in collections]
    if not values or key not in {"id", "sequence"}:
        raise ConfigurationError("Provide collections and key id/sequence.")
    shared = set(getattr(record, key) for record in values[0])
    for collection in values[1:]:
        shared &= {getattr(record, key) for record in collection}
    return tuple(record for record in values[0] if getattr(record, key) in shared)


def rename_reads(
    reads: Iterable[Read], pattern: str, replacement: str, *, verbose: bool = True
) -> Iterator[Read]:
    expression = re.compile(pattern)
    for read in reads:
        yield Read(
            expression.sub(replacement, read.id), read.sequence, read.qualities, read.description
        )


def sequence_digest(reads: Iterable[Read], *, verbose: bool = True) -> str:
    """Order-sensitive SHA256 of length-delimited uppercase sequences, independent of IDs."""
    digest = hashlib.sha256()
    for read in reads:
        encoded = read.sequence.upper().encode("ascii")
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)
    return digest.hexdigest()


@dataclass(frozen=True)
class Region:
    reference: str
    start: int
    end: int
    strand: str = "+"
    name: str = ""

    def __post_init__(self) -> None:
        if (
            not self.reference
            or type(self.start) is not int
            or type(self.end) is not int
            or not 0 <= self.start < self.end
            or self.strand not in {"+", "-"}
        ):
            raise ConfigurationError(
                "Regions use nonempty 0-based half-open intervals and +/- strand."
            )


def read_regions(
    path: str | Path, *, format: str = "bed", feature: str = "exon", verbose: bool = True
) -> Iterator[Region]:
    """BED intervals or GTF feature rows; GTF's 1-based coordinates are converted once."""
    if format not in {"bed", "gtf"}:
        raise ConfigurationError("format must be bed or gtf.")
    with _open(input_path(path), "r") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track ", "browser ")):
                continue
            fields = line.rstrip("\r\n").split("\t")
            try:
                if format == "bed":
                    yield Region(
                        fields[0],
                        int(fields[1]),
                        int(fields[2]),
                        fields[5] if len(fields) >= 6 and fields[5] != "." else "+",
                        fields[3] if len(fields) >= 4 else "",
                    )
                elif fields[2] == feature:
                    attributes: dict[str, str] = dict(re.findall(r'(\w+)\s+"([^"]*)"', fields[8]))
                    yield Region(
                        fields[0],
                        int(fields[3]) - 1,
                        int(fields[4]),
                        fields[6],
                        attributes.get("transcript_id", attributes.get("gene_id", "")),
                    )
            except (IndexError, ValueError) as exc:
                raise InputFormatError("Malformed genomic interval line.") from exc


class IndexedReference:
    """Context-managed htslib FASTA/BGZF random interval access; requires an existing .fai."""

    def __init__(self, path: str | Path) -> None:
        source = input_path(path)
        if not Path(str(source) + ".fai").is_file():
            raise ConfigurationError("Reference index missing; call index_reference first.")
        self._handle = require("pysam", "bio-reads").FastaFile(str(source))

    def __enter__(self) -> IndexedReference:
        return self

    def __exit__(self, *args: Any) -> None:
        self.close()

    def close(self) -> None:
        self._handle.close()

    def fetch(self, region: Region) -> Read:
        if (
            region.reference not in self._handle.references
            or region.end > self._handle.get_reference_length(region.reference)
        ):
            raise ConfigurationError("Interval exceeds the indexed reference.")
        text = self._handle.fetch(region.reference, region.start, region.end).upper()
        if region.strand == "-":
            text = text.translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1]
        return Read(
            region.name or f"{region.reference}:{region.start}-{region.end}:{region.strand}", text
        )

    def extract(self, regions: Iterable[Region]) -> Iterator[Read]:
        for region in regions:
            yield self.fetch(region)


def index_reference(path: str | Path, *, verbose: bool = True) -> Path:
    """Explicitly create the standard .fai (and BGZF .gzi when needed) beside the input."""
    source = input_path(path)
    if Path(str(source) + ".fai").exists():
        raise ConfigurationError(
            "Reference index already exists; validate/remove it explicitly before rebuilding."
        )
    require("pysam", "bio-reads").faidx(str(source))
    return Path(str(source) + ".fai")


def alignment_summary(path: str | Path, *, verbose: bool = True) -> dict[str, int]:
    """Stream SAM/BAM flags; secondary/supplementary records are counted explicitly."""
    pysam = require("pysam", "bio-reads")
    counts: Counter[str] = Counter(
        {
            key: 0
            for key in (
                "records",
                "mapped",
                "unmapped",
                "secondary",
                "supplementary",
                "duplicates",
                "proper_pair",
            )
        }
    )
    source = input_path(path)
    with pysam.AlignmentFile(str(source), "r" if source.suffix == ".sam" else "rb") as handle:
        for read in handle.fetch(until_eof=True):
            counts["records"] += 1
            counts["unmapped" if read.is_unmapped else "mapped"] += 1
            for key, condition in (
                ("secondary", read.is_secondary),
                ("supplementary", read.is_supplementary),
                ("duplicates", read.is_duplicate),
                ("proper_pair", read.is_proper_pair),
            ):
                counts[key] += bool(condition)
    return dict(counts)


def coverage_consensus(
    path: str | Path,
    region: Region,
    *,
    min_base_quality: int = 20,
    min_mapping_quality: int = 20,
    min_depth: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Indexed BAM regional A/C/G/T depth and majority consensus; ties/uncovered bases become N.

    Excludes unmapped, secondary, supplementary, duplicate and QC-failed reads.
    Counts aligned bases only; insertions/deletions are not inserted into consensus.
    """
    np = require("numpy", "bio-analysis")
    pysam = require("pysam", "bio-reads")
    positive(min_depth, "min_depth")
    if region.end - region.start > 10_000_000:
        raise ConfigurationError("Coverage interval exceeds 10000000 bases.")
    if min_base_quality < 0 or min_mapping_quality < 0:
        raise ConfigurationError("Quality cutoffs must be nonnegative.")
    counts = np.zeros((region.end - region.start, 4), dtype=int)
    with pysam.AlignmentFile(str(input_path(path)), "rb") as handle:
        for read in handle.fetch(region.reference, region.start, region.end):
            if (
                read.flag & (4 | 256 | 512 | 1024 | 2048)
                or read.mapping_quality < min_mapping_quality
            ):
                continue
            if read.query_sequence is None or read.query_qualities is None:
                continue
            for query_pos, reference_pos in read.get_aligned_pairs(matches_only=True):
                if (
                    region.start <= reference_pos < region.end
                    and read.query_qualities[query_pos] >= min_base_quality
                ):
                    base = read.query_sequence[query_pos].upper()
                    if base in "ACGT":
                        counts[reference_pos - region.start, "ACGT".index(base)] += 1
    depth = counts.sum(axis=1)
    symbols = [
        "ACGT"[int(row.argmax())]
        if depth[index] >= min_depth and sum(row == row.max()) == 1
        else "N"
        for index, row in enumerate(counts)
    ]
    return {"region": region, "counts": counts, "depth": depth, "consensus": "".join(symbols)}


_register_public_verbose(
    globals(),
    (
        "read_fastx",
        "write_fastx",
        "pair_reads",
        "repair_fastq",
        "arrange_reads",
        "common_reads",
        "rename_reads",
        "sequence_digest",
        "read_regions",
        "index_reference",
        "alignment_summary",
        "coverage_consensus",
    ),
)
