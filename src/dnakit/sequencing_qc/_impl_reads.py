"""Read parsing and small numeric helpers shared by the native sequencing-QC functions.

Reads are modelled by :class:`Rd` (name, upper-case sequence, optional Phred scores, header
comment). All parsers are bounded (``max_reads``, input size) and raise DNAKit exceptions with a
``code`` and a ``hint`` on malformed input. Nothing here touches the network or writes files.
"""

from __future__ import annotations

import gzip
import math
import re
from collections import Counter
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from dnakit._candidate_runtime import candidate_ok
from dnakit.exceptions import ConfigurationError, InputFormatError

MAX_PHRED = 93
DEFAULT_MAX_READS = 1_000_000
HARD_MAX_READS = 5_000_000
MAX_TEXT_BYTES = 400_000_000
_SEQ_OK = re.compile(r"^[A-Za-z.\-*]*$")
_COMPLEMENT = str.maketrans("ACGTNRYKMBDHVSWacgtnrykmbdhvsw", "TGCANYRMKVHDBSWtgcanyrmkvhdbsw")


@dataclass(frozen=True)
class Rd:
    """One read: identifier, upper-case sequence, Phred scores (or ``None``) and header comment."""

    name: str
    seq: str
    qual: tuple[int, ...] | None = None
    comment: str = ""


def revcomp(seq: str) -> str:
    """Reverse complement of IUPAC DNA text (unknown symbols are kept)."""

    return seq.translate(_COMPLEMENT)[::-1]


# --------------------------------------------------------------------------- option helpers


def opt_int(options: Mapping[str, object], key: str, default: int, *, minimum: int, maximum: int) -> int:
    value = options.get(key, default)
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(
            f"option {key!r} must be an integer in [{minimum}, {maximum}]; got {value!r}.",
            code="INVALID_OPTION",
            hint=f"Pass options={{'{key}': <int>}}.",
        )
    return value


def opt_num(
    options: Mapping[str, object], key: str, default: float | None, *, minimum: float, maximum: float
) -> float | None:
    value = options.get(key, default)
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(
            f"option {key!r} must be a finite number; got {value!r}.", code="INVALID_OPTION"
        )
    if not minimum <= value <= maximum:
        raise ConfigurationError(
            f"option {key!r} must be in [{minimum}, {maximum}]; got {value!r}.", code="INVALID_OPTION"
        )
    return float(value)


def opt_bool(options: Mapping[str, object], key: str, default: bool) -> bool:
    value = options.get(key, default)
    if not isinstance(value, bool):
        raise ConfigurationError(f"option {key!r} must be a boolean; got {value!r}.", code="INVALID_OPTION")
    return value


def opt_choice(options: Mapping[str, object], key: str, default: str, allowed: Sequence[str]) -> str:
    value = options.get(key, default)
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(
            f"option {key!r} must be one of {', '.join(allowed)}; got {value!r}.",
            code="INVALID_OPTION",
        )
    return value


def phred_offset(options: Mapping[str, object]) -> int:
    value = options.get("phred_offset", 33)
    if value not in (33, 64) or isinstance(value, bool):
        raise ConfigurationError(
            f"phred_offset must be 33 or 64; got {value!r}.",
            code="INVALID_OPTION",
            hint="Illumina 1.3-1.7 files use 64; all current platforms use 33.",
        )
    return int(value)


def max_reads(options: Mapping[str, object]) -> int:
    return opt_int(options, "max_reads", DEFAULT_MAX_READS, minimum=1, maximum=HARD_MAX_READS)


# --------------------------------------------------------------------------- parsing


def _looks_like_path(text: str) -> bool:
    return "\n" not in text and "\r" not in text and 0 < len(text) < 4096 and "\0" not in text


def _read_text(source: Path) -> str:
    if not source.is_file():
        raise InputFormatError(
            f"input file not found: {source}", code="INPUT_FILE_NOT_FOUND", hint="Pass an existing file or the text itself."
        )
    if source.stat().st_size > MAX_TEXT_BYTES:
        raise InputFormatError(f"input file exceeds {MAX_TEXT_BYTES} bytes.", code="INPUT_TOO_LARGE")
    if source.suffix == ".gz":
        with gzip.open(source, "rt", encoding="utf-8", errors="replace") as handle:
            return handle.read(MAX_TEXT_BYTES + 1)
    return source.read_text(encoding="utf-8", errors="replace")


def as_text(value: object) -> str | None:
    """Return text for ``str``/``Path``/``bytes`` inputs (reading files), else ``None``."""

    if isinstance(value, Path):
        return _read_text(value)
    if isinstance(value, (bytes, bytearray)):
        return bytes(value).decode("utf-8", errors="replace")
    if isinstance(value, str):
        if _looks_like_path(value) and not value.lstrip().startswith(("@", ">")):
            try:
                candidate = Path(value)
                if candidate.is_file():
                    return _read_text(candidate)
            except OSError:
                pass
        return value
    return None


def _clean_seq(seq: str, where: str) -> str:
    seq = seq.strip()
    if not _SEQ_OK.match(seq):
        raise InputFormatError(
            f"{where}: sequence contains characters other than letters, '.', '-' or '*'.",
            code="INVALID_SEQUENCE",
            hint="Reads must be plain nucleotide text.",
        )
    return seq.upper().replace(".", "N")


def decode_quality(qual: object, length: int, offset: int, where: str) -> tuple[int, ...]:
    """Decode a Phred string (or list of integers) to integer scores and validate its length."""

    if isinstance(qual, str):
        scores = tuple(ord(char) - offset for char in qual)
    elif isinstance(qual, Sequence) and not isinstance(qual, (bytes, bytearray)):
        try:
            scores = tuple(int(item) for item in qual)
        except (TypeError, ValueError) as exc:
            raise InputFormatError(f"{where}: quality list must contain integers.", code="INVALID_QUALITY") from exc
    else:
        raise InputFormatError(f"{where}: quality must be a string or list of integers.", code="INVALID_QUALITY")
    if len(scores) != length:
        raise InputFormatError(
            f"{where}: quality length {len(scores)} differs from sequence length {length}.",
            code="QUALITY_LENGTH_MISMATCH",
        )
    if any(not 0 <= score <= MAX_PHRED for score in scores):
        raise InputFormatError(
            f"{where}: quality values outside Phred 0-{MAX_PHRED} (offset {offset}).",
            code="INVALID_QUALITY",
            hint="Check the phred_offset option (33 or 64).",
        )
    return scores


def _split_header(header: str) -> tuple[str, str]:
    parts = header.strip().split(None, 1)
    if not parts:
        return "", ""
    return parts[0], parts[1] if len(parts) > 1 else ""


def parse_fastq_text(text: str, offset: int, limit: int) -> list[Rd]:
    lines = [line.rstrip("\r") for line in text.split("\n")]
    while lines and not lines[-1].strip():
        lines.pop()
    if len(lines) % 4:
        raise InputFormatError(
            "FASTQ text must consist of complete four-line records.", code="INVALID_FASTQ",
            hint="Each record is '@name', sequence, '+', quality.",
        )
    reads: list[Rd] = []
    for index in range(0, len(lines), 4):
        header, seq, plus, qual = lines[index : index + 4]
        where = f"FASTQ record {index // 4 + 1}"
        if not header.startswith("@") or not plus.startswith("+"):
            raise InputFormatError(f"{where}: missing '@' header or '+' separator.", code="INVALID_FASTQ")
        if len(reads) >= limit:
            raise InputFormatError(f"more than max_reads={limit} reads.", code="TOO_MANY_READS", hint="Raise options['max_reads'].")
        name, comment = _split_header(header[1:])
        seq = _clean_seq(seq, where)
        reads.append(Rd(name, seq, decode_quality(qual, len(seq), offset, where), comment))
    return reads


def parse_fasta_text(text: str, limit: int) -> list[Rd]:
    reads: list[Rd] = []
    name = comment = None
    parts: list[str] = []

    def flush() -> None:
        if name is not None:
            reads.append(Rd(name, _clean_seq("".join(parts), f"FASTA record {name}"), None, comment or ""))

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            if len(reads) >= limit:
                raise InputFormatError(f"more than max_reads={limit} reads.", code="TOO_MANY_READS")
            name, comment = _split_header(line[1:])
            parts = []
        elif name is None:
            raise InputFormatError("FASTA sequence appears before the first header.", code="INVALID_FASTA")
        else:
            parts.append(line)
    flush()
    return reads


def _record_to_read(item: object, index: int, offset: int) -> Rd:
    where = f"read {index + 1}"
    if isinstance(item, Rd):
        return item
    if isinstance(item, str):
        return Rd(f"read{index + 1}", _clean_seq(item, where))
    if isinstance(item, Mapping):
        seq_value = item.get("sequence", item.get("seq"))
        if not isinstance(seq_value, str):
            raise InputFormatError(f"{where}: mapping needs a 'sequence' string.", code="INVALID_READ")
        seq = _clean_seq(seq_value, where)
        name = str(item.get("id", item.get("name", item.get("identifier", f"read{index + 1}"))))
        qual = item.get("quality", item.get("qual"))
        scores = None if qual is None else decode_quality(qual, len(seq), offset, where)
        return Rd(name, seq, scores, str(item.get("comment", "")))
    if hasattr(item, "sequence") and hasattr(item, "identifier"):  # FastqRecord-like
        seq = _clean_seq(str(item.sequence), where)
        qual = getattr(item, "quality", None)
        scores = None if qual is None else decode_quality(qual, len(seq), offset, where)
        return Rd(str(item.identifier), seq, scores, "")
    if isinstance(item, Sequence) and not isinstance(item, (bytes, bytearray)) and 2 <= len(item) <= 3:
        seq = _clean_seq(str(item[1]), where)
        scores = decode_quality(item[2], len(seq), offset, where) if len(item) == 3 and item[2] is not None else None
        return Rd(str(item[0]), seq, scores, "")
    raise InputFormatError(
        f"{where}: unsupported record type {type(item).__name__}.",
        code="INVALID_READ",
        hint="Use FASTQ text, a file path, strings, (id, seq[, qual]) tuples or {'sequence','quality'} mappings.",
    )


_READS_KEYS = ("reads", "records", "sequences", "items")


def parse_reads(value: object, options: Mapping[str, object], *, name: str = "input") -> list[Rd]:
    """Parse reads from FASTQ/FASTA text or path, strings, tuples, mappings or record objects."""

    offset = phred_offset(options)
    limit = max_reads(options)
    if isinstance(value, Mapping):
        for key in _READS_KEYS:
            if key in value:
                return parse_reads(value[key], options, name=name)
        if "sequence" in value or "seq" in value:
            return [_record_to_read(value, 0, offset)]
        raise InputFormatError(
            f"{name}: mapping has none of the keys {', '.join(_READS_KEYS)}.", code="INVALID_INPUT",
            hint="Wrap the reads as {'reads': [...]}.",
        )
    text = as_text(value)
    if text is not None:
        stripped = text.lstrip()
        if not stripped:
            raise InputFormatError(f"{name} is empty.", code="EMPTY_INPUT", hint="Provide at least one read.")
        if stripped.startswith("@"):
            return parse_fastq_text(text, offset, limit)
        if stripped.startswith(">"):
            return parse_fasta_text(text, limit)
        lines = [line.strip() for line in text.splitlines() if line.strip()]
        if len(lines) > limit:
            raise InputFormatError(f"more than max_reads={limit} reads.", code="TOO_MANY_READS")
        return [Rd(f"read{i + 1}", _clean_seq(line, f"line {i + 1}")) for i, line in enumerate(lines)]
    if isinstance(value, Sequence) and not isinstance(value, (bytes, bytearray)):
        if not value:
            raise InputFormatError(f"{name} is empty.", code="EMPTY_INPUT", hint="Provide at least one read.")
        if len(value) > limit:
            raise InputFormatError(f"more than max_reads={limit} reads.", code="TOO_MANY_READS")
        return [_record_to_read(item, index, offset) for index, item in enumerate(value)]
    raise InputFormatError(
        f"{name}: unsupported input type {type(value).__name__}.", code="INVALID_INPUT",
        hint="Use FASTQ/FASTA text, a path, or a list of reads.",
    )


_MATE_KEYS = (("read1", "read2"), ("r1", "r2"), ("forward", "reverse"), ("R1", "R2"))


def pair_keys(value: object) -> tuple[str, str] | None:
    if isinstance(value, Mapping):
        for first, second in _MATE_KEYS:
            if first in value and second in value:
                return first, second
    return None


def parse_pair(value: object, options: Mapping[str, object]) -> tuple[list[Rd], list[Rd] | None]:
    """Parse reads; a mapping with read1/read2 (or options['read2']) yields a pair."""

    keys = pair_keys(value)
    if keys is not None:
        assert isinstance(value, Mapping)
        first = parse_reads(value[keys[0]], options, name=keys[0])
        second = parse_reads(value[keys[1]], options, name=keys[1])
        check_pair(first, second)
        return first, second
    first = parse_reads(value, options)
    if options.get("read2") is not None:
        second = parse_reads(options["read2"], options, name="read2")
        check_pair(first, second)
        return first, second
    if options.get("interleaved"):
        if len(first) % 2:
            raise InputFormatError("interleaved input has an odd number of reads.", code="UNPAIRED_READS")
        return first[::2], first[1::2]
    return first, None


def check_pair(first: Sequence[Rd], second: Sequence[Rd]) -> None:
    if len(first) != len(second):
        raise InputFormatError(
            f"read1 has {len(first)} reads but read2 has {len(second)}.",
            code="UNPAIRED_READS",
            hint="Paired files must contain the same reads in the same order.",
        )


def need_quality(reads: Sequence[Rd], what: str) -> None:
    if not any(read.qual for read in reads):
        raise InputFormatError(
            f"{what} needs Phred qualities but the reads have none.",
            code="QUALITY_REQUIRED",
            hint="Provide FASTQ input (FASTA and bare sequences carry no qualities).",
        )


# --------------------------------------------------------------------------- numeric helpers


def mean(values: Sequence[float]) -> float | None:
    return math.fsum(values) / len(values) if values else None


def stdev(values: Sequence[float]) -> float | None:
    """Sample standard deviation (n - 1 denominator); ``None`` for fewer than two values."""

    if len(values) < 2:
        return None
    centre = math.fsum(values) / len(values)
    return math.sqrt(math.fsum((v - centre) ** 2 for v in values) / (len(values) - 1))


def percentile(sorted_values: Sequence[float], fraction: float) -> float | None:
    """Percentile by linear interpolation between order statistics (numpy default)."""

    if not sorted_values:
        return None
    position = (len(sorted_values) - 1) * fraction
    low = math.floor(position)
    high = min(low + 1, len(sorted_values) - 1)
    return sorted_values[low] + (sorted_values[high] - sorted_values[low]) * (position - low)


def kth(counts: Sequence[int], k: int) -> int:
    """0-based k-th smallest value of a histogram ``counts[value]``."""

    running = 0
    for value, count in enumerate(counts):
        running += count
        if k < running:
            return value
    return len(counts) - 1


def hist_percentile(counts: Sequence[int], fraction: float) -> float | None:
    total = sum(counts)
    if total == 0:
        return None
    position = (total - 1) * fraction
    low = math.floor(position)
    low_value = kth(counts, low)
    high_value = kth(counts, min(low + 1, total - 1))
    return low_value + (high_value - low_value) * (position - low)


def nx(lengths: Sequence[int], fraction: float) -> int | None:
    """N50-style statistic: the length at which the sorted-descending running sum reaches ``fraction``."""

    if not lengths:
        return None
    target = sum(lengths) * fraction
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= target:
            return length
    return min(lengths)


def lx(lengths: Sequence[int], fraction: float) -> int | None:
    if not lengths:
        return None
    target = sum(lengths) * fraction
    running = 0
    for index, length in enumerate(sorted(lengths, reverse=True), 1):
        running += length
        if running >= target:
            return index
    return len(lengths)


def q_to_p(score: float) -> float:
    return 10.0 ** (-score / 10.0)


def p_to_q(probability: float) -> float | None:
    if probability <= 0:
        return None
    return -10.0 * math.log10(probability)


def gc_fraction(seq: str) -> float | None:
    """GC fraction among unambiguous A/C/G/T bases (N and other symbols excluded)."""

    acgt = sum(seq.count(base) for base in "ACGT")
    if acgt == 0:
        return None
    return (seq.count("G") + seq.count("C")) / acgt


def hamming(a: str, b: str) -> int:
    """Hamming distance; ``N`` in either string counts as a mismatch."""

    return sum(1 for x, y in zip(a, b, strict=False) if x != y or x == "N")


def length_counts(reads: Sequence[Rd]) -> Counter[int]:
    return Counter(len(read.seq) for read in reads)


def round6(value: float | None) -> float | None:
    return None if value is None else round(value, 6)


def run_op(module: str, operation: str, fn: Callable[[object, Mapping[str, object], list[str]], object], value: object, options: Mapping[str, object] | None) -> dict[str, object]:
    """Validate ``options``, run ``fn(value, options, warnings)`` and wrap the result in the candidate envelope."""

    if options is not None and not isinstance(options, Mapping):
        raise ConfigurationError("options must be a mapping or None.", code="INVALID_OPTION", hint="Pass options={'key': value}.")
    opts = dict(options or {})
    warnings: list[str] = []
    result = fn(value, opts, warnings)
    return candidate_ok(module, operation, value, result, options=opts, warnings=warnings)
