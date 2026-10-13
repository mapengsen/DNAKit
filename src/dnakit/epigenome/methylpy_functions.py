"""Optional methylpy wrappers for bisulfite read conversion, cytosine context and allc tables.

The functions here call methylpy's pure-Python helpers. Its shell pipelines (samtools, bowtie2,
cutadapt, tabix, bgzip, wigToBigWig) are not wrapped, so these functions need no external
executables. methylpy calls ``sys.exit`` on some bad inputs; those exits are converted into
``ConfigurationError``. Packages are imported only when a function runs.
"""

from __future__ import annotations

import contextlib
import io
import re
from collections.abc import Mapping, Sequence
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_METHYLPY_EXTRA = "bio-methylpy"
_READS_EXTRA = "bio-reads"
_MAX_TEXT = 1_000_000
_MAX_SEQUENCE = 100_000
_MAX_CONTEXT_LENGTH = 6
_MAX_FILES = 32
_MAX_CHROM_LENGTH = 300_000_000
_MAX_BUFFER_LINES = 10_000_000
_MAX_CLASSES = 64
_CONTEXT_PATTERN = re.compile(rf"[ACGTNHDBRYKMSW]{{1,{_MAX_CONTEXT_LENGTH}}}")
_SEQUENCE_PATTERN = re.compile(r"[ACGTNacgtn]+")
_ENCODING_PATTERN = re.compile(r"[!-~]*")
_BASECALL_REF = frozenset("ACGTN")


def _integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _number(value: object, name: str, minimum: float, maximum: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a number in [{minimum}, {maximum}].")
    number = float(value)
    if not minimum <= number <= maximum:
        raise ConfigurationError(f"{name} must be a number in [{minimum}, {maximum}].")
    return number


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _text(value: object, name: str, limit: int = _MAX_TEXT) -> str:
    if (
        not isinstance(value, str)
        or not value.strip()
        or len(value) > limit
        or any(char in value for char in "\x00\r\n\t")
    ):
        raise ConfigurationError(f"{name} must be a non-empty string without tabs or newlines.")
    return value


def _path(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a file path.")
    return input_path(value)


def _new_path(value: object, name: str) -> Path:
    """Output path that must not exist yet; its parent directory must exist."""

    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a file path.")
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} already exists and will not be overwritten: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The directory for {name} does not exist: {path.parent}")
    return path


def _plain(value: Any) -> Any:
    if hasattr(value, "item"):
        value = value.item()
    return value


def _provenance(statistic: str, **notes: str) -> dict[str, Any]:
    try:
        installed = version("methylpy")
    except PackageNotFoundError:
        installed = "unknown"
    return {"package": f"methylpy {installed}", "method": statistic, **notes}


def _upstream_call(function: Any, *args: Any, **kwargs: Any) -> Any:
    """Run a methylpy function that prints progress and may call sys.exit on bad input."""

    captured = io.StringIO()
    try:
        with contextlib.redirect_stdout(captured), contextlib.redirect_stderr(captured):
            return function(*args, **kwargs)
    except SystemExit as exc:
        lines = captured.getvalue().strip().splitlines()
        detail = lines[-1] if lines else str(exc.code)
        raise ConfigurationError(f"methylpy rejected the input: {detail}") from exc


def _require_fasta_index(fasta: Path) -> Path:
    index = Path(str(fasta) + ".fai")
    if not index.is_file():
        raise ConfigurationError(
            f"{fasta.name} needs a FASTA index at {index.name}; create it with "
            "samtools faidx or pysam.faidx first."
        )
    return index


def _fasta_index_length(fasta: Path, chrom: str) -> int:
    """Length of one chromosome from the samtools-style .fai index next to the FASTA file."""

    index = _require_fasta_index(fasta)
    with index.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields and fields[0] == chrom:
                return int(fields[1])
    raise ConfigurationError(f"{chrom} is not listed in {index.name}.")


def _line_count(path: Path) -> int:
    with path.open("r", encoding="utf-8") as handle:
        return sum(1 for _ in handle)


def expand_methylation_context(mc_type: str, *, verbose: bool = True) -> dict[str, Any]:
    """Expand an IUPAC cytosine-context pattern such as ``"CGN"`` into matched context strings."""

    pattern = _text(mc_type, "mc_type", limit=_MAX_CONTEXT_LENGTH)
    if not _CONTEXT_PATTERN.fullmatch(pattern):
        raise ConfigurationError(
            f"mc_type must be 1 to {_MAX_CONTEXT_LENGTH} IUPAC nucleotide codes, for example CGN."
        )
    utilities = require("methylpy.utilities", _METHYLPY_EXTRA)
    contexts = sorted(_upstream_call(utilities.expand_nucleotide_code, pattern))
    return {
        "status": "completed",
        "operation": "expand_methylation_context",
        "mc_type": pattern,
        "contexts": contexts,
        "count": len(contexts),
        "provenance": _provenance(
            "IUPAC expansion of the pattern, plus the CGN/CHG/CHH/CNN classes when it contains C",
            note="methylpy also returns single-letter codes from the pattern; they are kept",
        ),
    }


def encode_read_c_positions(
    seq: str, *, is_read2: bool = False, verbose: bool = True
) -> dict[str, Any]:
    """Encode the cytosine positions of a read as printable offset characters."""

    sequence = _text(seq, "seq", limit=_MAX_SEQUENCE)
    if not _SEQUENCE_PATTERN.fullmatch(sequence):
        raise ConfigurationError("seq must contain only A, C, G, T or N letters.")
    flag = _flag(is_read2, "is_read2")
    call_mc = require("methylpy.call_mc_se", _METHYLPY_EXTRA)
    encoding = _upstream_call(call_mc.encode_c_positions, sequence, flag)
    return {
        "status": "completed",
        "operation": "encode_read_c_positions",
        "sequence": sequence,
        "is_read2": flag,
        "encoding": encoding,
        "provenance": _provenance(
            "offset encoding of C positions; each character holds the gap to the next C",
        ),
    }


def decode_read_c_positions(
    seq: str, indexes: str, strand: str, *, is_read2: bool = False, verbose: bool = True
) -> dict[str, Any]:
    """Restore the cytosines of a converted read from its offset encoding and strand."""

    sequence = _text(seq, "seq", limit=_MAX_SEQUENCE)
    if not _SEQUENCE_PATTERN.fullmatch(sequence):
        raise ConfigurationError("seq must contain only A, C, G, T or N letters.")
    if not isinstance(indexes, str) or not _ENCODING_PATTERN.fullmatch(indexes):
        raise ConfigurationError("indexes must be printable ASCII characters, possibly empty.")
    if strand not in ("+", "-"):
        raise ConfigurationError("strand must be '+' or '-'.")
    flag = _flag(is_read2, "is_read2")
    call_mc = require("methylpy.call_mc_se", _METHYLPY_EXTRA)
    decoded = _upstream_call(call_mc.decode_c_positions, sequence, indexes, strand, flag)
    return {
        "status": "completed",
        "operation": "decode_read_c_positions",
        "sequence": decoded,
        "strand": strand,
        "is_read2": flag,
        "provenance": _provenance("restores C (or G on the reverse strand) at encoded positions"),
    }


def count_bisulfite_basecalls(
    ref_base: str, read_bases: str, *, verbose: bool = True
) -> dict[str, Any]:
    """Count consistent and inconsistent calls of one reference base from pileup base characters."""

    if not isinstance(ref_base, str) or ref_base not in _BASECALL_REF:
        raise ConfigurationError("ref_base must be one of A, C, G, T or N.")
    bases = _text(read_bases, "read_bases", limit=_MAX_TEXT)
    call_mc = require("methylpy.call_mc_se", _METHYLPY_EXTRA)
    values = list(_upstream_call(call_mc.analyze_read_basecalls, ref_base, bases))
    unconverted = int(values[2]) if len(values) > 2 else None
    converted = int(values[3]) if len(values) > 3 else None
    return {
        "status": "completed",
        "operation": "count_bisulfite_basecalls",
        "ref_base": ref_base,
        "consistent_calls": int(values[0]),
        "inconsistent_calls": int(values[1]),
        "unconverted_c": unconverted,
        "converted_c": converted,
        "provenance": _provenance(
            "pileup characters: '.' and ',' match the reference, 'T' or 'a' mark conversion",
            note="indels are skipped; unconverted_c and converted_c apply to C and G only",
        ),
    }


def read_fasta_chromosome(
    fasta_path: str | Path,
    chrom: str,
    *,
    max_length: int = 50_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read one chromosome sequence from an indexed FASTA file."""

    source = _path(fasta_path, "fasta_path")
    name = _text(chrom, "chrom")
    limit = _integer(max_length, "max_length", 1, _MAX_CHROM_LENGTH)
    length = _fasta_index_length(source, name)
    if length > limit:
        raise ConfigurationError(f"{name} has {length} bases, above max_length={limit}.")
    call_mc = require("methylpy.call_mc_se", _METHYLPY_EXTRA)
    sequence = _upstream_call(call_mc.get_chromosome_sequence, str(source), name)
    if len(sequence) != length:
        raise ConfigurationError(f"{name} read {len(sequence)} bases; the index lists {length}.")
    return {
        "status": "completed",
        "operation": "read_fasta_chromosome",
        "fasta": str(source),
        "chrom": name,
        "length": length,
        "sequence": sequence,
        "provenance": _provenance("whole-chromosome read through the FASTA index"),
    }


def convert_fastq_to_bisulfite(
    fastq_path: str | Path,
    output_path: str | Path,
    *,
    buffer_line_number: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert every C to T in FASTQ reads for alignment to a bisulfite-converted reference."""

    source = _path(fastq_path, "fastq_path")
    target = _new_path(output_path, "output_path")
    buffer_lines = _integer(buffer_line_number, "buffer_line_number", 4, _MAX_BUFFER_LINES)
    call_mc = require("methylpy.call_mc_se", _METHYLPY_EXTRA)
    _upstream_call(
        call_mc.convert_reads,
        str(source),
        str(target),
        buffer_line_number=buffer_lines,
    )
    reads = _line_count(target) // 4
    return {
        "status": "completed",
        "operation": "convert_fastq_to_bisulfite",
        "input": str(source),
        "output": str(target),
        "reads": reads,
        "provenance": _provenance(
            "C to T conversion of every read; original C positions are stored in the read name",
            note="reverse-strand reads are converted when mapped, not here",
        ),
    }


def filter_bam_by_mch_level(
    bam_path: str | Path,
    output_bam: str | Path,
    reference_fasta: str | Path,
    *,
    min_ch: int = 3,
    max_mch_level: float = 0.7,
    buffer_line_number: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Drop aligned reads whose non-CpG methylation (mCH) is above a level."""

    source = _path(bam_path, "bam_path")
    target = _new_path(output_bam, "output_bam")
    if target.suffix != ".bam":
        raise ConfigurationError("output_bam must end with .bam.")
    reference = _path(reference_fasta, "reference_fasta")
    _require_fasta_index(reference)
    informative = _integer(min_ch, "min_ch", 1, 1_000)
    level = _number(max_mch_level, "max_mch_level", 0.0, 1.0)
    buffer_lines = _integer(buffer_line_number, "buffer_line_number", 1, _MAX_BUFFER_LINES)
    pysam = require("pysam", _READS_EXTRA)
    with pysam.AlignmentFile(str(source), "rb") as handle:
        input_reads = handle.count(until_eof=True)
    call_mc = require("methylpy.call_mc_se", _METHYLPY_EXTRA)
    # methylpy appends ".bam" to the prefix it receives, so pass the output path without it.
    _upstream_call(
        call_mc.bam_quality_mch_filter,
        str(source),
        str(target)[: -len(".bam")],
        str(reference),
        min_ch=informative,
        max_mch_level=level,
        buffer_line_number=buffer_lines,
    )
    with pysam.AlignmentFile(str(target), "rb") as handle:
        kept_reads = handle.count(until_eof=True)
    return {
        "status": "completed",
        "operation": "filter_bam_by_mch_level",
        "input": str(source),
        "output": str(target),
        "input_reads": input_reads,
        "kept_reads": kept_reads,
        "removed_reads": input_reads - kept_reads,
        "min_ch": informative,
        "max_mch_level": level,
        "provenance": _provenance(
            "a read is removed when it has at least min_ch informative CH sites and "
            "mCH/(mCH+uCH) is at least max_mch_level",
            note="methylpy accepts min_mapq but does not apply it; it is not exposed here",
        ),
    }


def benjamini_hochberg_methylation_cutoffs(
    binomial_paths: Sequence[str | Path],
    mc_class_counts: Mapping[str, int],
    *,
    sig_cutoff: float = 0.01,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find the P-value cutoff per cytosine context at a Benjamini-Hochberg FDR."""

    if isinstance(binomial_paths, (str, Path)) or not isinstance(binomial_paths, Sequence):
        raise ConfigurationError("binomial_paths must be a list of binomial-test files.")
    if not 1 <= len(binomial_paths) <= _MAX_FILES:
        raise ConfigurationError(f"binomial_paths must list 1 to {_MAX_FILES} files.")
    paths = [_path(item, "binomial_paths") for item in binomial_paths]
    if len(set(paths)) != len(paths):
        raise ConfigurationError("binomial_paths must not list the same file twice.")
    if not isinstance(mc_class_counts, Mapping) or not 1 <= len(mc_class_counts) <= _MAX_CLASSES:
        raise ConfigurationError("mc_class_counts must map 1 to 64 contexts to test counts.")
    counts = {
        _text(key, "mc_class_counts key", limit=32): _integer(value, f"count for {key}", 1, 10**12)
        for key, value in mc_class_counts.items()
    }
    level = _number(sig_cutoff, "sig_cutoff", 0.0, 1.0)
    if level <= 0.0:
        raise ConfigurationError("sig_cutoff must be greater than zero.")
    sites = 0
    for path in paths:
        previous = -1.0
        with path.open("r", encoding="utf-8") as handle:
            for number, line in enumerate(handle, start=1):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 7:
                    raise ConfigurationError(f"{path.name} line {number} has fewer than 7 fields.")
                try:
                    p_value = float(fields[6])
                except ValueError as exc:
                    raise ConfigurationError(
                        f"{path.name} line {number} has a non-numeric P-value."
                    ) from exc
                if p_value < previous:
                    raise ConfigurationError(
                        f"{path.name} must be sorted by P-value; line {number} is out of order."
                    )
                if fields[3] not in counts:
                    raise ConfigurationError(f"No test count was given for context {fields[3]}.")
                previous = p_value
                sites += 1
    utilities = require("methylpy.call_mc_se", _METHYLPY_EXTRA)
    cutoffs = _upstream_call(
        utilities.benjamini_hochberg_correction_call_methylated_sites,
        [str(path) for path in paths],
        counts,
        level,
    )
    return {
        "status": "completed",
        "operation": "benjamini_hochberg_methylation_cutoffs",
        "files": [str(path) for path in paths],
        "tested_sites": sites,
        "sig_cutoff": level,
        "cutoffs": {key: _plain(cutoffs[key]) for key in sorted(cutoffs)},
        "provenance": _provenance(
            "Benjamini-Hochberg adjustment per cytosine context over sorted binomial P-values",
            note="a context with no significant site reports 0",
        ),
    }


def index_allc_table(
    allc_path: str | Path, *, reindex: bool = False, verbose: bool = True
) -> dict[str, Any]:
    """Build the byte-offset index next to an allc table whose chromosomes are contiguous."""

    source = _path(allc_path, "allc_path")
    flag = _flag(reindex, "reindex")
    utilities = require("methylpy.utilities", _METHYLPY_EXTRA)
    _upstream_call(utilities.index_allc_file, str(source), reindex=flag)
    index = Path(utilities.get_index_file_name(str(source)))
    if not index.is_file():
        raise ConfigurationError(f"methylpy did not create the index {index.name}.")
    offsets = _upstream_call(utilities.read_allc_index, str(source))
    return {
        "status": "completed",
        "operation": "index_allc_table",
        "allc": str(source),
        "index": str(index),
        "chromosomes": [
            {"chrom": chrom, "offset": int(offsets[chrom])} for chrom in sorted(offsets)
        ],
        "provenance": _provenance(
            "byte offset of each chromosome's first record in the allc table",
            note="the index is written next to the input file; reindex=True rebuilds it",
        ),
    }


def read_allc_index(allc_path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Read the chromosome offsets stored in the index of an allc table."""

    source = _path(allc_path, "allc_path")
    index = Path(str(source) + ".idx")
    if not index.is_file():
        raise ConfigurationError(f"{index.name} is missing; run index_allc_table first.")
    utilities = require("methylpy.utilities", _METHYLPY_EXTRA)
    offsets = _upstream_call(utilities.read_allc_index, str(source))
    return {
        "status": "completed",
        "operation": "read_allc_index",
        "allc": str(source),
        "index": str(index),
        "chromosomes": [
            {"chrom": chrom, "offset": int(offsets[chrom])} for chrom in sorted(offsets)
        ],
        "provenance": _provenance("reads the index file written by methylpy"),
    }


__all__ = [
    "benjamini_hochberg_methylation_cutoffs",
    "convert_fastq_to_bisulfite",
    "count_bisulfite_basecalls",
    "decode_read_c_positions",
    "encode_read_c_positions",
    "expand_methylation_context",
    "filter_bam_by_mch_level",
    "index_allc_table",
    "read_allc_index",
    "read_fasta_chromosome",
]

_register_public_verbose(globals(), tuple(__all__))
