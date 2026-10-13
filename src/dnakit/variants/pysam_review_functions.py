"""Offline pysam wrappers for alignment files (samtools) and variant files (bcftools).

Every function runs the samtools or bcftools build bundled with pysam, with stdout captured.
Inputs must already exist locally. Outputs are written only to new paths unless ``force=True``.
Commands that pysam redirects to its own temporary stdout file (view, mpileup, depad, and most
bcftools subcommands) are called without ``-o``; their captured output is written to the target
path here. Output files are never passed to the pysam dispatcher with ``catch_stdout=False``.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``bio-reads`` install hint.
"""

from __future__ import annotations

import importlib
import importlib.metadata
import re
import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require, workspace
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_EXTRA = "bio-reads"
_ALIGNMENT_INDEXES = (".bai", ".csi", ".crai")
_MAX_SOURCES = 1_000
_MAX_RECORDS = 1_000
_MAX_DEPTH_ROWS = 1_000_000
_REGION = re.compile(r"^(?P<chrom>[^\s:,]+)(?::(?P<start>\d[\d,]*)-(?P<end>\d[\d,]*))?$")
_DECODE = ("utf-8", "replace")
_PHASE_SUFFIXES = (".0.bam", ".1.bam", ".chimera.bam")


def _pysam() -> Any:
    return require("pysam", _EXTRA)


def _provenance(convention: str) -> dict[str, str]:
    return {
        "package": "pysam",
        "version": importlib.metadata.version("pysam"),
        "convention": convention,
    }


def _text(value: object, name: str) -> str:
    if isinstance(value, Path):
        value = str(value)
    if not isinstance(value, str) or not value.strip() or any(c in value for c in "\x00\r\n"):
        raise ConfigurationError(f"{name} must be a non-empty single-line string.")
    return value.strip()


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _choice(value: object, name: str, allowed: Sequence[str]) -> str:
    if value not in allowed:
        raise ConfigurationError(f"{name} must be one of: {', '.join(allowed)}.")
    return str(value)


def _sources(values: object, name: str) -> list[Path]:
    if isinstance(values, (str, bytes, Path)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list of existing local file paths.")
    if not values or len(values) > _MAX_SOURCES:
        raise ConfigurationError(f"{name} must hold between 1 and {_MAX_SOURCES} paths.")
    return [input_path(_text(item, name)) for item in values]


def _input(value: object, name: str) -> Path:
    return input_path(_text(value, name))


def _output(value: object, name: str, *, force: bool) -> Path:
    path = Path(_text(value, name)).expanduser().resolve()
    if not path.parent.is_dir():
        raise ConfigurationError(f"{name} directory does not exist: {path.parent}")
    _flag(force, "force")
    if path.is_dir():
        raise ConfigurationError(f"{name} is a directory: {path}")
    if path.exists() and not force:
        raise ConfigurationError(f"{name} already exists: {path}; pass force=True to replace it.")
    return path


def _region(value: object) -> str:
    text = _text(value, "region")
    match = _REGION.match(text)
    if match is None:
        raise ConfigurationError(
            "region must look like 'chr1' or 'chr1:1000-2000' (1-based, inclusive)."
        )
    if match["start"] is not None and int(match["start"].replace(",", "")) < 1:
        raise ConfigurationError("region positions are 1-based and must be at least 1.")
    return text


def _has_index(path: Path) -> bool:
    return any(Path(str(path) + suffix).is_file() for suffix in _ALIGNMENT_INDEXES)


def _decoded(raw: object) -> str:
    if isinstance(raw, bytes):
        return raw.decode(*_DECODE)
    return "" if raw is None else str(raw)


def _command(program: str, subcommand: str, arguments: Sequence[str]) -> object:
    """Run one samtools or bcftools subcommand through pysam with stdout captured."""

    pysam = _pysam()
    if program == "samtools":
        dispatcher = getattr(pysam.samtools, subcommand)
    else:
        dispatcher = getattr(importlib.import_module("pysam.bcftools"), subcommand)
    try:
        return dispatcher(*arguments, catch_stdout=True)
    except (pysam.SamtoolsError, OSError, RuntimeError, ValueError) as exc:
        raise BackendExecutionError(f"{program} {subcommand} failed: {exc}") from exc


def _write_raw(path: Path, raw: object) -> None:
    if isinstance(raw, bytes):
        path.write_bytes(raw)
    else:
        path.write_text(_decoded(raw), encoding="utf-8")


def _count_records(path: Path, *flags: str) -> int:
    raw = _command("samtools", "view", ("-c", *flags, str(path)))
    return int(_decoded(raw).strip())


def _int_or_float(text: str) -> int | float | str:
    for convert in (int, float):
        try:
            return convert(text)
        except ValueError:
            continue
    return text


def _lines(raw: object) -> list[str]:
    return [line for line in _decoded(raw).splitlines() if line.strip()]


def _mapping_of_rows(
    rows: list[list[str]], columns: Sequence[str]
) -> list[dict[str, int | float | str]]:
    return [
        {column: _int_or_float(value) for column, value in zip(columns, row, strict=False)}
        for row in rows
    ]


# ---------------------------------------------------------------------------
# Alignment files (samtools)
# ---------------------------------------------------------------------------


def add_bam_read_group(
    source: str | Path,
    output_path: str | Path,
    read_group: str,
    *,
    mode: str = "overwrite_all",
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Set a read-group header line and tag every read with its ID (``samtools addreplacerg``)."""

    path = _input(source, "source")
    target = _output(output_path, "output_path", force=force)
    line = _text(read_group, "read_group")
    if not line.startswith("@RG") or "ID:" not in line:
        raise ConfigurationError("read_group must be an @RG header line that contains ID:.")
    selected = _choice(mode, "mode", ("overwrite_all", "orphan_only"))
    _command(
        "samtools",
        "addreplacerg",
        ("-r", line, "-m", selected, "-o", str(target), str(path)),
    )
    read_group_id = line.split("ID:", 1)[1].split("\t", 1)[0]
    return {
        "status": "completed",
        "operation": "add_bam_read_group",
        "output": str(target),
        "read_group_id": read_group_id,
        "mode": selected,
        "records": _count_records(target),
        "provenance": _provenance("RG tag set on every read; mode overwrite_all or orphan_only"),
    }


def clip_amplicon_primers(
    source: str | Path,
    primers_bed: str | Path,
    output_path: str | Path,
    *,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Soft-clip amplicon primer regions from aligned reads (``samtools ampliconclip``)."""

    path = _input(source, "source")
    bed = _input(primers_bed, "primers_bed")
    target = _output(output_path, "output_path", force=force)
    _command(
        "samtools",
        "ampliconclip",
        ("-b", str(bed), "-o", str(target), str(path)),
    )
    return {
        "status": "completed",
        "operation": "clip_amplicon_primers",
        "output": str(target),
        "records": _count_records(target),
        "provenance": _provenance(
            "primer intervals from BED are soft-clipped; coordinates unchanged"
        ),
    }


def extract_fastq_from_bam(
    source: str | Path,
    output_directory: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write BAM reads to FASTQ files split by read role (``samtools fastq``).

    Reads without READ1/READ2 flags go to ``other.fastq``; READ1, READ2 and singleton reads
    go to ``read1.fastq``, ``read2.fastq`` and ``singletons.fastq`` in the new directory.
    """

    path = _input(source, "source")
    directory = workspace(output_directory)
    files = {
        "other": directory / "other.fastq",
        "read1": directory / "read1.fastq",
        "read2": directory / "read2.fastq",
        "singletons": directory / "singletons.fastq",
    }
    _command(
        "samtools",
        "fastq",
        (
            "-0",
            str(files["other"]),
            "-1",
            str(files["read1"]),
            "-2",
            str(files["read2"]),
            "-s",
            str(files["singletons"]),
            str(path),
        ),
    )
    counts = {
        role: sum(
            1 for line in file.read_text(encoding="utf-8").splitlines() if line.startswith("@")
        )
        if file.stat().st_size
        else 0
        for role, file in files.items()
    }
    return {
        "status": "completed",
        "operation": "extract_fastq_from_bam",
        "directory": str(directory),
        "files": {role: str(file) for role, file in files.items()},
        "records": dict(sorted(counts.items())),
        "provenance": _provenance("read order and names are kept; qualities are Phred+33"),
    }


def extract_fasta_from_bam(
    source: str | Path,
    output_directory: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write BAM sequences to FASTA files split by read role (``samtools fasta``)."""

    path = _input(source, "source")
    directory = workspace(output_directory)
    files = {
        "other": directory / "other.fasta",
        "read1": directory / "read1.fasta",
        "read2": directory / "read2.fasta",
        "singletons": directory / "singletons.fasta",
    }
    _command(
        "samtools",
        "fasta",
        (
            "-0",
            str(files["other"]),
            "-1",
            str(files["read1"]),
            "-2",
            str(files["read2"]),
            "-s",
            str(files["singletons"]),
            str(path),
        ),
    )
    counts = {
        role: sum(
            1 for line in file.read_text(encoding="utf-8").splitlines() if line.startswith(">")
        )
        if file.stat().st_size
        else 0
        for role, file in files.items()
    }
    return {
        "status": "completed",
        "operation": "extract_fasta_from_bam",
        "directory": str(directory),
        "files": {role: str(file) for role, file in files.items()},
        "records": dict(sorted(counts.items())),
        "provenance": _provenance("reverse-complemented reads are restored to reference strand"),
    }


def collate_bam_reads(
    source: str | Path,
    output_path: str | Path,
    *,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Group the reads of a BAM by name, keeping mates adjacent (``samtools collate``)."""

    path = _input(source, "source")
    target = _output(output_path, "output_path", force=force)
    _command("samtools", "collate", ("-o", str(target), str(path)))
    return {
        "status": "completed",
        "operation": "collate_bam_reads",
        "output": str(target),
        "records": _count_records(target),
        "provenance": _provenance("read names group records; order within each group is preserved"),
    }


def concatenate_bam_files(
    sources: Sequence[str | Path],
    output_path: str | Path,
    *,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Concatenate BAM files that share one header, in the order given (``samtools cat``)."""

    paths = _sources(sources, "sources")
    target = _output(output_path, "output_path", force=force)
    _command("samtools", "cat", ("-o", str(target), *map(str, paths)))
    return {
        "status": "completed",
        "operation": "concatenate_bam_files",
        "output": str(target),
        "inputs": [str(path) for path in paths],
        "records": _count_records(target),
        "provenance": _provenance("headers must agree; records keep the input order"),
    }


def count_bam_reads_per_interval(
    intervals_bed: str | Path,
    source: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count the bases covered by reads within each BED interval (``samtools bedcov``)."""

    bed = _input(intervals_bed, "intervals_bed")
    path = _input(source, "source")
    if not _has_index(path):
        raise ConfigurationError(
            "bedcov needs a BAM index (.bai or .csi); run index_alignment_file first."
        )
    rows = []
    for line in _lines(_command("samtools", "bedcov", (str(bed), str(path)))):
        fields = line.split("\t")
        if len(fields) < 4:
            raise BackendExecutionError("bedcov returned an unexpected row.")
        rows.append(
            {
                "chrom": fields[0],
                "start": int(fields[1]),
                "end": int(fields[2]),
                "covered_bases": int(fields[-1]),
            }
        )
    return {
        "status": "completed",
        "operation": "count_bam_reads_per_interval",
        "intervals": rows,
        "count": len(rows),
        "provenance": _provenance("covered bases are summed over BED intervals; BED is 0-based"),
    }


def recompute_nm_md_tags(
    source: str | Path,
    reference: str | Path,
    output_path: str | Path,
    *,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Recompute the NM and MD tags of every aligned read against a FASTA (``samtools calmd``)."""

    path = _input(source, "source")
    fasta = _input(reference, "reference")
    if not Path(str(fasta) + ".fai").is_file():
        raise ConfigurationError("reference needs a .fai index; run samtools faidx on it first.")
    target = _output(output_path, "output_path", force=force)
    # calmd runs as SAM text; with -b the pysam dispatch was observed to return nothing for
    # absolute paths, so the BAM is produced by a second view call on the text.
    text = _decoded(_command("samtools", "calmd", (str(path), str(fasta))))
    if target.suffix == ".sam":
        target.write_text(text, encoding="utf-8")
    else:
        with tempfile.TemporaryDirectory() as scratch:
            sam = Path(scratch) / "calmd.sam"
            sam.write_text(text, encoding="utf-8")
            _write_raw(target, _command("samtools", "view", ("-b", str(sam))))
    return {
        "status": "completed",
        "operation": "recompute_nm_md_tags",
        "output": str(target),
        "records": _count_records(target),
        "provenance": _provenance("NM and MD tags describe the alignment against this reference"),
    }


def summarize_bam_coverage(source: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Per-contig read count, covered bases, and mean depth and quality (``samtools coverage``)."""

    path = _input(source, "source")
    columns = (
        "rname", "startpos", "endpos", "numreads", "covbases", "coverage",
        "meandepth", "meanbaseq", "meanmapq",
    )  # fmt: skip
    rows = [line.split("\t") for line in _lines(_command("samtools", "coverage", (str(path),)))]
    contigs = _mapping_of_rows([row for row in rows if not row[0].startswith("#")], columns)
    return {
        "status": "completed",
        "operation": "summarize_bam_coverage",
        "contigs": contigs,
        "count": len(contigs),
        "provenance": _provenance("coverage excludes records filtered by samtools defaults"),
    }


def calculate_bam_read_depth(
    source: str | Path,
    *,
    region: str | None = None,
    all_positions: bool = True,
    max_positions: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Per-base read depth for a BAM (``samtools depth``).

    ``all_positions=True`` reports zero-depth positions too. A region needs a BAM index.
    """

    path = _input(source, "source")
    _flag(all_positions, "all_positions")
    limit = _integer(max_positions, "max_positions", minimum=1, maximum=_MAX_DEPTH_ROWS)
    arguments = ["-a"] if all_positions else []
    if region is not None:
        arguments.extend(["-r", _region(region)])
        if not _has_index(path):
            raise ConfigurationError("A region query needs a BAM index (.bai or .csi).")
    rows = _lines(_command("samtools", "depth", (*arguments, str(path))))
    if len(rows) > limit:
        raise ConfigurationError(f"Depth has {len(rows)} positions; raise max_positions.")
    positions = []
    for line in rows:
        contig, position, depth = line.split("\t")[:3]
        positions.append({"contig": contig, "position": int(position), "depth": int(depth)})
    return {
        "status": "completed",
        "operation": "calculate_bam_read_depth",
        "positions": positions,
        "count": len(positions),
        "provenance": _provenance(
            "1-based positions; depth counts records passing samtools defaults"
        ),
    }


def write_fasta_sequence_dictionary(
    reference: str | Path,
    output_path: str | Path,
    *,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the SAM header dictionary (``.dict``) for a FASTA reference (``samtools dict``)."""

    fasta = _input(reference, "reference")
    target = _output(output_path, "output_path", force=force)
    if target.suffix != ".dict":
        raise ConfigurationError("output_path must end with .dict.")
    raw = _command("samtools", "dict", (str(fasta),))
    text = _decoded(raw)
    target.write_text(text, encoding="utf-8")
    contigs = []
    for line in text.splitlines():
        if line.startswith("@SQ"):
            fields = dict(item.split(":", 1) for item in line.split("\t")[1:] if ":" in item)
            contigs.append(
                {"name": fields["SN"], "length": int(fields["LN"]), "md5": fields.get("M5", "")}
            )
    return {
        "status": "completed",
        "operation": "write_fasta_sequence_dictionary",
        "output": str(target),
        "contigs": contigs,
        "provenance": _provenance("MD5 values are computed from the uppercase FASTA sequence"),
    }


def fix_bam_mate_fields(
    source: str | Path,
    output_path: str | Path,
    *,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fill mate position, insert size and mate-score tags of name-sorted paired reads."""

    path = _input(source, "source")
    target = _output(output_path, "output_path", force=force)
    try:
        _command("samtools", "fixmate", ("-m", str(path), str(target)))
    except BackendExecutionError as exc:
        raise ConfigurationError(
            "fixmate needs name-sorted input; run sort_alignment_file(by_name=True) first."
        ) from exc
    return {
        "status": "completed",
        "operation": "fix_bam_mate_fields",
        "output": str(target),
        "records": _count_records(target),
        "provenance": _provenance("input must be sorted by read name; ms tags are added"),
    }


def decode_sam_flags(flags: Sequence[int], *, verbose: bool = True) -> dict[str, Any]:
    """Name the bits of SAM FLAG integers (``samtools flags``), such as 99 -> PAIRED,..."""

    if isinstance(flags, (str, bytes)) or not isinstance(flags, Sequence) or not flags:
        raise ConfigurationError("flags must be a non-empty list of integers.")
    if len(flags) > _MAX_RECORDS:
        raise ConfigurationError(f"flags may hold at most {_MAX_RECORDS} values.")
    values = [_integer(item, "flags item", minimum=0, maximum=0xFFFF) for item in flags]
    rows = _lines(_command("samtools", "flags", tuple(str(value) for value in values)))
    decoded = []
    for value, line in zip(values, rows, strict=True):
        hex_text, number, names = [*line.split("\t"), "", ""][:3]
        decoded.append(
            {
                "flag": int(number),
                "hex": hex_text,
                "names": sorted(name for name in names.split(",") if name),
                "input": value,
            }
        )
    return {
        "status": "completed",
        "operation": "decode_sam_flags",
        "flags": decoded,
        "provenance": _provenance("bit names follow the SAM specification"),
    }


def index_fastq_file(
    source: str | Path, *, force: bool = False, verbose: bool = True
) -> dict[str, Any]:
    """Write a ``.fai`` index beside a FASTQ file, one row per record (``samtools fqidx``)."""

    path = _input(source, "source")
    index = Path(str(path) + ".fai")
    _flag(force, "force")
    if index.exists() and not force:
        raise ConfigurationError(f"index already exists: {index}; pass force=True to rebuild.")
    _command("samtools", "fqidx", (str(path),))
    entries = [line for line in index.read_text(encoding="utf-8").splitlines() if line.strip()]
    return {
        "status": "completed",
        "operation": "index_fastq_file",
        "index": str(index),
        "entries": len(entries),
        "provenance": _provenance(
            "index columns: name, length, offset, bases per line, bytes per line"
        ),
    }


def import_fastq_to_bam(
    sources: Sequence[str | Path],
    output_path: str | Path,
    *,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert one FASTQ (single-end) or two FASTQ files (paired) into an unaligned BAM."""

    paths = _sources(sources, "sources")
    if len(paths) > 2:
        raise ConfigurationError("sources holds one FASTQ or a pair of FASTQ files.")
    target = _output(output_path, "output_path", force=force)
    _command("samtools", "fqimport", ("-o", str(target), *map(str, paths)))
    return {
        "status": "completed",
        "operation": "import_fastq_to_bam",
        "output": str(target),
        "layout": "paired" if len(paths) == 2 else "single",
        "records": _count_records(target),
        "provenance": _provenance("reads are unaligned; sequences and qualities are kept"),
    }


def read_bam_header_and_records(
    source: str | Path,
    *,
    records: int = 5,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the header lines and the first records of an alignment file (``samtools head``)."""

    path = _input(source, "source")
    limit = _integer(records, "records", minimum=0, maximum=_MAX_RECORDS)
    lines = _lines(_command("samtools", "head", ("-n", str(limit), str(path))))
    header = [line for line in lines if line.startswith("@")]
    rows = [line.split("\t") for line in lines if not line.startswith("@")]
    return {
        "status": "completed",
        "operation": "read_bam_header_and_records",
        "header": header,
        "records": rows,
        "provenance": _provenance("header lines first, then the first records as SAM columns"),
    }


def summarize_bam_index_stats(source: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Mapped and unmapped read counts per contig from the index (``samtools idxstats``)."""

    path = _input(source, "source")
    if not _has_index(path):
        raise ConfigurationError(
            "idxstats needs a BAM or CRAM index; run index_alignment_file first."
        )
    contigs, unplaced = [], 0
    for line in _lines(_command("samtools", "idxstats", (str(path),))):
        name, length, mapped, unmapped = line.split("\t")[:4]
        if name == "*":
            unplaced = int(unmapped)
            continue
        contigs.append(
            {
                "contig": name,
                "length": int(length),
                "mapped": int(mapped),
                "unmapped": int(unmapped),
            }
        )
    return {
        "status": "completed",
        "operation": "summarize_bam_index_stats",
        "contigs": contigs,
        "unplaced_unmapped": unplaced,
        "provenance": _provenance(
            "counts come from the index; unplaced reads are reported separately"
        ),
    }


def index_alignment_file(
    source: str | Path,
    *,
    csi: bool = False,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a BAI (or CSI with ``csi=True``) index beside a coordinate-sorted BAM or CRAM."""

    path = _input(source, "source")
    _flag(csi, "csi")
    _flag(force, "force")
    suffix = ".csi" if csi else (".crai" if path.suffix == ".cram" else ".bai")
    index = Path(str(path) + suffix)
    if index.exists() and not force:
        raise ConfigurationError(f"index already exists: {index}; pass force=True to rebuild.")
    arguments = ["-c"] if csi else []
    _command("samtools", "index", (*arguments, str(path)))
    return {
        "status": "completed",
        "operation": "index_alignment_file",
        "index": str(index),
        "format": "CSI" if csi else ("CRAI" if suffix == ".crai" else "BAI"),
        "provenance": _provenance("index requires coordinate-sorted input"),
    }


def mark_bam_duplicates(
    source: str | Path,
    output_path: str | Path,
    *,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Flag duplicate reads by 5' position and orientation (``samtools markdup``)."""

    path = _input(source, "source")
    target = _output(output_path, "output_path", force=force)
    _command("samtools", "markdup", (str(path), str(target)))
    return {
        "status": "completed",
        "operation": "mark_bam_duplicates",
        "output": str(target),
        "records": _count_records(target),
        "duplicates": _count_records(target, "-f", "1024"),
        "provenance": _provenance("duplicate flag 0x400 is set; records are not removed"),
    }


def merge_bam_files(
    sources: Sequence[str | Path],
    output_path: str | Path,
    *,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge coordinate-sorted BAM files with compatible headers (``samtools merge``)."""

    paths = _sources(sources, "sources")
    target = _output(output_path, "output_path", force=force)
    arguments = ["-f"] if force else []
    _command("samtools", "merge", (*arguments, str(target), *map(str, paths)))
    return {
        "status": "completed",
        "operation": "merge_bam_files",
        "output": str(target),
        "inputs": [str(path) for path in paths],
        "records": _count_records(target),
        "provenance": _provenance("coordinate-sorted inputs give a coordinate-sorted output"),
    }


def pileup_alignment_positions(
    source: str | Path,
    output_path: str | Path,
    *,
    reference: str | Path | None = None,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Per-position base pileup as text (``samtools mpileup``), one line per covered position."""

    path = _input(source, "source")
    target = _output(output_path, "output_path", force=force)
    arguments: list[str] = []
    if reference is not None:
        fasta = _input(reference, "reference")
        if not Path(str(fasta) + ".fai").is_file():
            raise ConfigurationError(
                "reference needs a .fai index; run samtools faidx on it first."
            )
        arguments.extend(["-f", str(fasta)])
    raw = _command("samtools", "mpileup", (*arguments, str(path)))
    _write_raw(target, raw)
    return {
        "status": "completed",
        "operation": "pileup_alignment_positions",
        "output": str(target),
        "positions": len(_lines(raw)),
        "provenance": _provenance("pileup columns follow the samtools mpileup text format"),
    }


def check_alignment_file_integrity(
    sources: Sequence[str | Path], *, verbose: bool = True
) -> dict[str, Any]:
    """Check that each BAM or CRAM opens and has a valid header (``samtools quickcheck``)."""

    paths = _sources(sources, "sources")
    files = []
    for path in paths:
        try:
            _command("samtools", "quickcheck", (str(path),))
            files.append({"path": str(path), "valid": True, "message": ""})
        except BackendExecutionError as exc:
            files.append({"path": str(path), "valid": False, "message": str(exc)})
    return {
        "status": "completed",
        "operation": "check_alignment_file_integrity",
        "files": files,
        "valid": sum(1 for item in files if item["valid"]),
        "invalid": sum(1 for item in files if not item["valid"]),
        "provenance": _provenance("quickcheck checks headers and the first records of each file"),
    }


def replace_bam_header(
    header_sam: str | Path,
    source: str | Path,
    output_path: str | Path,
    *,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a BAM with a replaced header and unchanged records (``samtools reheader``)."""

    header = _input(header_sam, "header_sam")
    path = _input(source, "source")
    target = _output(output_path, "output_path", force=force)
    if not header.read_text(encoding="utf-8", errors="replace").startswith("@"):
        raise ConfigurationError("header_sam must be a SAM header that starts with @ lines.")
    raw = _command("samtools", "reheader", (str(header), str(path)))
    _write_raw(target, raw)
    return {
        "status": "completed",
        "operation": "replace_bam_header",
        "output": str(target),
        "records": _count_records(target),
        "provenance": _provenance("only the header changes; records are copied unchanged"),
    }


def list_read_group_samples(
    sources: Sequence[str | Path], *, verbose: bool = True
) -> dict[str, Any]:
    """Sample names declared by the read groups of each alignment file (``samtools samples``)."""

    paths = _sources(sources, "sources")
    samples = []
    for line in _lines(_command("samtools", "samples", tuple(str(path) for path in paths))):
        sample, _, name = line.partition("\t")
        samples.append({"sample": sample, "source": name})
    return {
        "status": "completed",
        "operation": "list_read_group_samples",
        "samples": samples,
        "provenance": _provenance("sample names come from @RG SM fields"),
    }


def sort_alignment_file(
    source: str | Path,
    output_path: str | Path,
    *,
    by_name: bool = False,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sort reads by coordinate, or by name with ``by_name=True`` (``samtools sort``)."""

    path = _input(source, "source")
    target = _output(output_path, "output_path", force=force)
    _flag(by_name, "by_name")
    arguments = ["-n"] if by_name else []
    _command("samtools", "sort", (*arguments, "-o", str(target), str(path)))
    return {
        "status": "completed",
        "operation": "sort_alignment_file",
        "output": str(target),
        "order": "queryname" if by_name else "coordinate",
        "records": _count_records(target),
        "provenance": _provenance("sort keys follow the SAM specification"),
    }


def split_bam_by_read_group(
    source: str | Path,
    output_directory: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write one BAM per read group into a new directory (``samtools split``)."""

    path = _input(source, "source")
    directory = workspace(output_directory)
    _command("samtools", "split", ("-f", str(directory / "%!.bam"), str(path)))
    outputs = sorted(file for file in directory.glob("*.bam"))
    return {
        "status": "completed",
        "operation": "split_bam_by_read_group",
        "directory": str(directory),
        "files": [{"file": file.name, "records": _count_records(file)} for file in outputs],
        "provenance": _provenance("every read must carry an RG tag; file names are RG IDs"),
    }


def summarize_alignment_stats(source: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Summary numbers from the SN lines of ``samtools stats``: reads, mapping and length."""

    path = _input(source, "source")
    summary: dict[str, int | float | str] = {}
    for line in _lines(_command("samtools", "stats", (str(path),))):
        if not line.startswith("SN\t"):
            continue
        _, key, value = line.split("\t")[:3]
        summary[key.rstrip(":").strip()] = _int_or_float(value)
    return {
        "status": "completed",
        "operation": "summarize_alignment_stats",
        "summary": dict(sorted(summary.items())),
        "provenance": _provenance("values are the SN summary lines of samtools stats"),
    }


def select_alignment_records(
    source: str | Path,
    output_path: str | Path,
    *,
    region: str | None = None,
    required_flags: int = 0,
    excluded_flags: int = 0,
    min_mapping_quality: int = 0,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep records by region, flag masks and mapping quality (``samtools view``).

    An output ending in ``.sam`` is written as text with its header; any other name is BAM.
    A region needs an index beside the input.
    """

    path = _input(source, "source")
    target = _output(output_path, "output_path", force=force)
    req = _integer(required_flags, "required_flags", minimum=0, maximum=0xFFFF)
    exc = _integer(excluded_flags, "excluded_flags", minimum=0, maximum=0xFFFF)
    mapq = _integer(min_mapping_quality, "min_mapping_quality", minimum=0, maximum=255)
    arguments = ["-h"] if target.suffix == ".sam" else ["-b"]
    arguments.extend(["-f", str(req), "-F", str(exc), "-q", str(mapq), str(path)])
    if region is not None:
        if not _has_index(path):
            raise ConfigurationError("A region query needs an index (.bai, .csi or .crai).")
        arguments.append(_region(region))
    raw = _command("samtools", "view", tuple(arguments))
    _write_raw(target, raw)
    return {
        "status": "completed",
        "operation": "select_alignment_records",
        "output": str(target),
        "records": _count_records(target),
        "provenance": _provenance("flags use SAM numbering; mapping quality is a minimum"),
    }


def summarize_amplicon_stats(
    source: str | Path,
    primers_bed: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Amplicon layout and per-file summary numbers from ``samtools ampliconstats``."""

    path = _input(source, "source")
    bed = _input(primers_bed, "primers_bed")
    amplicons: list[dict[str, Any]] = []
    summary: list[dict[str, Any]] = []
    for line in _lines(_command("samtools", "ampliconstats", (str(bed), str(path)))):
        fields = line.split("\t")
        if fields[0] == "AMPLICON" and len(fields) >= 5:
            amplicons.append(
                {
                    "reference": fields[1],
                    "number": int(fields[2]),
                    "left": fields[3],
                    "right": fields[4],
                }
            )
        elif fields[0] == "FSS" and len(fields) >= 5:
            summary.append(
                {
                    "file": fields[1],
                    "reference": fields[2],
                    "metric": fields[3].rstrip(":"),
                    "value": _int_or_float(fields[4]),
                }
            )
    return {
        "status": "completed",
        "operation": "summarize_amplicon_stats",
        "amplicons": amplicons,
        "summary": summary,
        "count": len(amplicons),
        "provenance": _provenance("AMPLICON and FSS rows of samtools ampliconstats"),
    }


def reconstruct_reference_from_alignments(
    source: str | Path,
    output_path: str | Path,
    *,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Rebuild reference bases from aligned reads with MD tags (``samtools reference``).

    Positions no read covers are written as N. Run ``recompute_nm_md_tags`` first when the
    alignments carry no MD tags.
    """

    path = _input(source, "source")
    target = _output(output_path, "output_path", force=force)
    text = _decoded(_command("samtools", "reference", (str(path),)))
    if not text.startswith(">"):
        raise BackendExecutionError("samtools reference returned no FASTA records.")
    target.write_text(text, encoding="utf-8")
    contigs: list[dict[str, Any]] = []
    name, bases = "", ""
    for line in [*text.splitlines(), ">"]:
        if line.startswith(">"):
            if name:
                contigs.append(
                    {
                        "name": name,
                        "length": len(bases),
                        "determined_bases": sum(1 for base in bases if base in "ACGT"),
                    }
                )
            name, bases = line[1:].strip(), ""
        else:
            bases += line.strip().upper()
    return {
        "status": "completed",
        "operation": "reconstruct_reference_from_alignments",
        "output": str(target),
        "contigs": contigs,
        "provenance": _provenance("bases come from read sequences; uncovered positions are N"),
    }


def phase_bam_haplotypes(
    source: str | Path,
    output_prefix: str | Path,
    *,
    block_length: int = 13,
    verbose: bool = True,
) -> dict[str, Any]:
    """Tag reads with HP haplotypes and split them by phase (``samtools phase``).

    Writes ``<prefix>.0.bam``, ``<prefix>.1.bam`` and ``<prefix>.chimera.bam``; existing files
    under the same prefix are refused.
    """

    path = _input(source, "source")
    block = _integer(block_length, "block_length", minimum=1, maximum=1_000_000)
    prefix = Path(_text(output_prefix, "output_prefix")).expanduser().resolve()
    if not prefix.parent.is_dir():
        raise ConfigurationError(f"output_prefix directory does not exist: {prefix.parent}")
    expected = [prefix.with_name(prefix.name + suffix) for suffix in _PHASE_SUFFIXES]
    if any(item.exists() for item in expected):
        raise ConfigurationError("phase outputs already exist under this prefix; choose a new one.")
    _command("samtools", "phase", ("-k", str(block), "-b", str(prefix), str(path)))
    files = [item for item in expected if item.exists()]
    return {
        "status": "completed",
        "operation": "phase_bam_haplotypes",
        "block_length": block,
        "files": [{"file": str(item), "records": _count_records(item)} for item in files],
        "provenance": _provenance("HP tags come from heterozygous sites with read support"),
    }


def reset_alignments_to_unmapped(
    source: str | Path,
    output_path: str | Path,
    *,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Clear alignment fields so every read is unmapped (``samtools reset``).

    Sequences, qualities and read-group tags are kept.
    """

    path = _input(source, "source")
    target = _output(output_path, "output_path", force=force)
    _command("samtools", "reset", ("-o", str(target), str(path)))
    return {
        "status": "completed",
        "operation": "reset_alignments_to_unmapped",
        "output": str(target),
        "records": _count_records(target),
        "provenance": _provenance("position, CIGAR and mapping flags are cleared"),
    }


def _vcf_records(text: str) -> int:
    return sum(1 for line in text.splitlines() if line and not line.startswith("#"))


def _vcf_indexed(path: Path) -> bool:
    return any(Path(str(path) + suffix).is_file() for suffix in (".tbi", ".csi"))


def _bgzipped_vcf(path: Path, name: str) -> None:
    if not (path.name.endswith(".vcf.gz") or path.name.endswith(".bcf")):
        raise ConfigurationError(f"{name} must be a bgzipped VCF (.vcf.gz) or BCF (.bcf).")


def filter_vcf_by_expression(
    source: str | Path,
    expression: str,
    output_path: str | Path,
    *,
    mode: str = "include",
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep (``mode='include'``) or drop (``mode='exclude'``) sites matching an expression."""

    path = _input(source, "source")
    text = _text(expression, "expression")
    selected = _choice(mode, "mode", ("include", "exclude"))
    target = _output(output_path, "output_path", force=force)
    flag = "-i" if selected == "include" else "-e"
    raw = _command("bcftools", "filter", (flag, text, str(path)))
    _write_raw(target, raw)
    return {
        "status": "completed",
        "operation": "filter_vcf_by_expression",
        "output": str(target),
        "mode": selected,
        "expression": text,
        "records": _vcf_records(_decoded(raw)),
        "provenance": _provenance("expression is one bcftools filter string; no shell is involved"),
    }


def annotate_vcf_consequence(
    source: str | Path,
    gene_annotation: str | Path,
    reference: str | Path,
    output_path: str | Path,
    *,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add BCSQ consequence tags from a GFF3 gene model and a FASTA reference (``bcftools csq``).

    The input must be phased: bcftools stops on unphased heterozygous genotypes.
    """

    path = _input(source, "source")
    gff = _input(gene_annotation, "gene_annotation")
    fasta = _input(reference, "reference")
    target = _output(output_path, "output_path", force=force)
    raw = _command("bcftools", "csq", ("-g", str(gff), "-f", str(fasta), str(path)))
    text = _decoded(raw)
    _write_raw(target, raw)
    return {
        "status": "completed",
        "operation": "annotate_vcf_consequence",
        "output": str(target),
        "records": _vcf_records(text),
        "annotated_records": sum(
            1 for line in text.splitlines() if "BCSQ=" in line and not line.startswith("#")
        ),
        "provenance": _provenance("consequences use the GFF3 CDS phase and the FASTA sequence"),
    }


def check_vcf_sample_concordance(
    source: str | Path,
    *,
    genotype_source: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Genotype discordance between each query sample and each sample of a second VCF (``gtcheck``).

    Without ``genotype_source`` the VCF is compared with itself.
    """

    path = _input(source, "source")
    other = path if genotype_source is None else _input(genotype_source, "genotype_source")
    pairs = []
    for line in _lines(_command("bcftools", "gtcheck", ("-g", str(other), str(path)))):
        fields = line.split("\t")
        if fields[0] != "DCv2" or len(fields) < 7:
            continue
        pairs.append(
            {
                "query": fields[1],
                "genotyped": fields[2],
                "discordance": float(fields[3]),
                "average_log_hwe_p": float(fields[4]),
                "sites": int(fields[5]),
                "matching": int(fields[6]),
            }
        )
    return {
        "status": "completed",
        "operation": "check_vcf_sample_concordance",
        "pairs": pairs,
        "count": len(pairs),
        "provenance": _provenance(
            "discordance is the -log10 likelihood of mismatch per sample pair"
        ),
    }


def read_vcf_header(source: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Header lines and sample names of a VCF or BCF file (``bcftools head``)."""

    path = _input(source, "source")
    header = [
        line for line in _lines(_command("bcftools", "head", (str(path),))) if line.startswith("#")
    ]
    columns = next((line.split("\t") for line in header if line.startswith("#CHROM")), [])
    return {
        "status": "completed",
        "operation": "read_vcf_header",
        "header_lines": header,
        "samples": columns[9:],
        "provenance": _provenance("header lines are copied as written in the file"),
    }


def index_vcf_csi(
    source: str | Path, *, force: bool = False, verbose: bool = True
) -> dict[str, Any]:
    """Write a CSI index beside a bgzipped VCF or BCF file (``bcftools index -c``)."""

    path = _input(source, "source")
    _bgzipped_vcf(path, "source")
    _flag(force, "force")
    index = Path(str(path) + ".csi")
    if index.exists() and not force:
        raise ConfigurationError(f"index already exists: {index}; pass force=True to rebuild.")
    arguments = ["-f"] if force else []
    _command("bcftools", "index", (*arguments, "-c", str(path)))
    return {
        "status": "completed",
        "operation": "index_vcf_csi",
        "index": str(index),
        "provenance": _provenance("CSI indexes allow regions on contigs longer than 512 Mb"),
    }


def merge_indexed_vcfs(
    sources: Sequence[str | Path],
    output_path: str | Path,
    *,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Join the samples of several indexed, bgzipped VCFs into one VCF (``bcftools merge``)."""

    paths = _sources(sources, "sources")
    for path in paths:
        _bgzipped_vcf(path, "sources item")
        if not _vcf_indexed(path):
            raise ConfigurationError(
                f"{path.name} needs a .tbi or .csi index; run index_vcf_csi first."
            )
    target = _output(output_path, "output_path", force=force)
    raw = _command("bcftools", "merge", tuple(str(path) for path in paths))
    text = _decoded(raw)
    _write_raw(target, raw)
    columns = next(
        (line.split("\t") for line in text.splitlines() if line.startswith("#CHROM")), []
    )
    return {
        "status": "completed",
        "operation": "merge_indexed_vcfs",
        "output": str(target),
        "samples": columns[9:],
        "records": _vcf_records(text),
        "provenance": _provenance("only records from different files are merged"),
    }


def pileup_bam_to_vcf(
    source: str | Path,
    reference: str | Path,
    output_path: str | Path,
    *,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Genotype likelihoods of an aligned BAM as VCF, before calling (``bcftools mpileup``)."""

    path = _input(source, "source")
    fasta = _input(reference, "reference")
    target = _output(output_path, "output_path", force=force)
    raw = _command("bcftools", "mpileup", ("-f", str(fasta), "-Ov", str(path)))
    text = _decoded(raw)
    _write_raw(target, raw)
    return {
        "status": "completed",
        "operation": "pileup_bam_to_vcf",
        "output": str(target),
        "records": _vcf_records(text),
        "provenance": _provenance("PL likelihoods are written; genotype calls need a later step"),
    }


def rename_vcf_samples(
    source: str | Path,
    sample_names: str | Path,
    output_path: str | Path,
    *,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Replace sample names with those listed in a file, one per line (``bcftools reheader``)."""

    path = _input(source, "source")
    names_file = _input(sample_names, "sample_names")
    names = [
        line.strip() for line in names_file.read_text(encoding="utf-8").splitlines() if line.strip()
    ]
    if not names or len(set(names)) != len(names):
        raise ConfigurationError("sample_names must list unique, non-empty sample names.")
    target = _output(output_path, "output_path", force=force)
    raw = _command("bcftools", "reheader", ("-s", str(names_file), str(path)))
    _write_raw(target, raw)
    return {
        "status": "completed",
        "operation": "rename_vcf_samples",
        "output": str(target),
        "samples": names,
        "provenance": _provenance("only the sample names and header sample line change"),
    }


def detect_bcftools_roh(
    source: str | Path,
    *,
    sample: str,
    genotype_quality: float = 30.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Runs of homozygosity for one sample with bcftools' HMM (``bcftools roh``).

    ``genotype_quality`` is the PL value given to the two least likely genotypes, so the
    input may hold GT fields only.
    """

    path = _input(source, "source")
    name = _text(sample, "sample")
    if isinstance(genotype_quality, bool) or not isinstance(genotype_quality, (int, float)):
        raise ConfigurationError("genotype_quality must be a positive number.")
    if not 0 < float(genotype_quality) <= 1000:
        raise ConfigurationError("genotype_quality must be in (0, 1000].")
    # bcftools roh writes its table to the process stdout unless -o names a file, so the
    # result is read back from a scratch file instead of the captured stream.
    with tempfile.TemporaryDirectory() as scratch:
        table = Path(scratch) / "roh.txt"
        _command(
            "bcftools",
            "roh",
            ("-G", str(float(genotype_quality)), "-s", name, "-o", str(table), str(path)),
        )
        text = table.read_text(encoding="utf-8")
    runs = []
    for line in _lines(text):
        fields = line.split("\t")
        if fields[0] == "RG" and len(fields) >= 8:
            runs.append(
                {
                    "sample": fields[1],
                    "chromosome": fields[2],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "length": int(fields[5]),
                    "markers": int(fields[6]),
                    "quality": float(fields[7]),
                }
            )
    return {
        "status": "completed",
        "operation": "detect_bcftools_roh",
        "sample": name,
        "runs": runs,
        "count": len(runs),
        "provenance": _provenance("HMM with default transition rates; positions are 1-based"),
    }


__all__ = [
    "add_bam_read_group",
    "annotate_vcf_consequence",
    "calculate_bam_read_depth",
    "check_alignment_file_integrity",
    "check_vcf_sample_concordance",
    "clip_amplicon_primers",
    "collate_bam_reads",
    "concatenate_bam_files",
    "count_bam_reads_per_interval",
    "decode_sam_flags",
    "detect_bcftools_roh",
    "extract_fasta_from_bam",
    "extract_fastq_from_bam",
    "filter_vcf_by_expression",
    "fix_bam_mate_fields",
    "import_fastq_to_bam",
    "index_alignment_file",
    "index_fastq_file",
    "index_vcf_csi",
    "list_read_group_samples",
    "mark_bam_duplicates",
    "merge_bam_files",
    "merge_indexed_vcfs",
    "phase_bam_haplotypes",
    "pileup_alignment_positions",
    "pileup_bam_to_vcf",
    "read_bam_header_and_records",
    "read_vcf_header",
    "recompute_nm_md_tags",
    "reconstruct_reference_from_alignments",
    "rename_vcf_samples",
    "replace_bam_header",
    "reset_alignments_to_unmapped",
    "select_alignment_records",
    "sort_alignment_file",
    "split_bam_by_read_group",
    "summarize_alignment_stats",
    "summarize_amplicon_stats",
    "summarize_bam_coverage",
    "summarize_bam_index_stats",
    "write_fasta_sequence_dictionary",
]

_register_public_verbose(globals(), tuple(__all__))
