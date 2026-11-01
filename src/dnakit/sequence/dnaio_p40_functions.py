"""dnaio read/write operations beyond single and paired reading, as DNAKit functions (pack P40).

The functions call the dnaio Python package in process on explicit local files (plain or gzip/bzip2/xz
compressed, as dnaio's ``open`` supports). Outputs are written to new paths; nothing is overwritten.
"""

from __future__ import annotations

import io
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p40_helpers as _h
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-dnaio"
_LAYOUTS = ("separate", "interleaved")
_FORMATS = ("fasta", "fastq")


def _dnaio() -> Any:
    return require("dnaio", _EXTRA)


def _entry(record: Any) -> dict[str, Any]:
    out = {"name": str(record.name), "sequence": str(record.sequence)}
    if record.qualities is not None:
        out["qualities"] = str(record.qualities)
    return out


def _make(module: Any, item: Mapping[str, Any], index: int) -> Any:
    if not isinstance(item, Mapping) or "name" not in item or "sequence" not in item:
        raise ConfigurationError(f"record {index} needs 'name' and 'sequence'.")
    qualities = item.get("qualities")
    if qualities is not None and len(qualities) != len(item["sequence"]):
        raise ConfigurationError(f"record {index}: qualities and sequence differ in length.")
    return module.SequenceRecord(str(item["name"]), str(item["sequence"]), None if qualities is None else str(qualities))


def dnaio_read_parallel_files(
    sources: Sequence[str | Path],
    *,
    max_records: int = 1_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read three or more parallel FASTA/FASTQ files record by record (dnaio.open with several files).

    ``sources`` lists the files (for example R1, R2 and an index/barcode read file); the files must have the
    same number of records in the same order, otherwise dnaio raises an error that is reported as a
    ConfigurationError. Returns one tuple of records (name, sequence, qualities) per position, bounded by
    ``max_records``.
    """

    module = _dnaio()
    limit = _h.integer(max_records, "max_records", minimum=1, maximum=_h.MAX_ROWS)
    if isinstance(sources, (str, Path)) or len(sources) < 2:
        raise ConfigurationError("sources must list at least two files.")
    paths = [_h.existing(item, f"sources[{i}]") for i, item in enumerate(sources)]
    rows: list[list[dict[str, Any]]] = []
    count = 0
    try:
        with module.open(*paths, mode="r") as reader:
            for group in reader:
                count += 1
                if len(rows) < limit:
                    rows.append([_entry(record) for record in group])
    except (module.FileFormatError, module.FastaFormatError, module.FastqFormatError, ValueError) as exc:
        raise ConfigurationError(f"dnaio could not read the files in parallel: {exc}") from exc
    return {"status": "completed", "backend": "dnaio", "version": module.__version__, "file_count": len(paths),
            "record_group_count": count, "groups": rows, "truncated": count > limit}


def dnaio_write_reads(
    outputs: Sequence[str | Path],
    reads: Sequence[Sequence[Mapping[str, Any]]],
    *,
    layout: str = "separate",
    file_format: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write parallel read lists with dnaio: paired-end files, parallel FASTA/FASTQ files or interleaved reads.

    ``reads`` holds one list of records per mate (each record ``{"name", "sequence", "qualities"?}``).
    ``layout="separate"`` writes mate ``i`` to ``outputs[i]`` (two files is the paired-end writer, three or
    more the parallel FASTA/FASTQ writers); ``layout="interleaved"`` writes the first record of each mate,
    then the second, ... into the single file ``outputs[0]``. ``file_format`` (``fasta`` or ``fastq``)
    defaults to the output file extension. All mates must have the same number of records.
    """

    module = _dnaio()
    mode = _h.choice(layout, _LAYOUTS, "layout")
    fmt = None if file_format is None else _h.choice(file_format, _FORMATS, "file_format")
    if not reads or any(isinstance(r, (str, bytes)) for r in reads):
        raise ConfigurationError("reads must be a list with one list of records per mate.")
    if len({len(r) for r in reads}) != 1:
        raise ConfigurationError("All mates must have the same number of records.")
    paths = [_h.new_file(item, f"outputs[{i}]") for i, item in enumerate(outputs)]
    if mode == "separate" and len(paths) != len(reads):
        raise ConfigurationError("layout='separate' needs one output path per mate.")
    if mode == "interleaved" and (len(paths) != 1 or len(reads) != 2):
        raise ConfigurationError("layout='interleaved' needs one output path and exactly two mates.")
    if len(paths) != len({str(p) for p in paths}):
        raise ConfigurationError("outputs must be distinct.")
    groups = [[_make(module, item, i) for i, item in enumerate(mate)] for mate in reads]
    kwargs: dict[str, Any] = {"mode": "w"}
    if fmt is not None:
        kwargs["fileformat"] = fmt
    if mode == "interleaved":
        kwargs["interleaved"] = True
    try:
        with module.open(*paths, **kwargs) as writer:
            for tuple_ in zip(*groups):
                writer.write(*tuple_)
    except Exception as exc:  # noqa: BLE001 - dnaio raises several format errors
        for path in paths:
            path.unlink(missing_ok=True)
        raise ConfigurationError(f"dnaio could not write the reads: {exc}") from exc
    return {"status": "completed", "backend": "dnaio", "version": module.__version__, "layout": mode,
            "mate_count": len(groups), "record_count_per_mate": len(groups[0]),
            "outputs": [str(p) for p in paths]}


def dnaio_split_into_chunks(
    source: str | Path,
    *,
    mate_source: str | Path | None = None,
    buffer_size: int = 4_194_304,
    output_dir: str | Path | None = None,
    max_chunks: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split a FASTA/FASTQ file into chunks of complete records, without parsing them (dnaio.read_chunks).

    With ``mate_source`` the two files of a pair are split into synchronised chunks with equal record counts
    (dnaio.read_paired_chunks). ``buffer_size`` is the largest chunk in bytes. The result gives byte size and
    record count of each chunk; with ``output_dir`` (new) every chunk is also written to
    ``chunk_<n>_1<ext>`` (and ``_2`` for the mate).
    """

    module = _dnaio()
    size = _h.integer(buffer_size, "buffer_size", minimum=1024)
    chunk_limit = _h.integer(max_chunks, "max_chunks", minimum=1, maximum=_h.MAX_ROWS)
    left = _h.existing(source, "source")
    right = None if mate_source is None else _h.existing(mate_source, "mate_source")
    folder = None if output_dir is None else _h.new_dir(output_dir, "output_dir")
    chunks: list[dict[str, Any]] = []
    files: list[str] = []
    total = 0

    def count(data: memoryview) -> int:
        marker = bytes(data[:1])
        if marker == b"@":
            return bytes(data).count(b"\n") // 4
        return sum(1 for line in bytes(data).split(b"\n") if line.startswith(b">"))

    def suffix(path: Path) -> str:
        return "".join(path.suffixes[-2:]) if path.suffix in (".gz", ".bz2", ".xz") else path.suffix

    try:
        from xopen import xopen  # dnaio depends on xopen

        with xopen(left, "rb") as one:
            if right is None:
                iterator = ((chunk, None) for chunk in module.read_chunks(one, size))
                two_cm = None
            else:
                two_cm = xopen(right, "rb")
                iterator = module.read_paired_chunks(one, two_cm, size)
            try:
                for number, (first, second) in enumerate(iterator):
                    total += 1
                    if len(chunks) < chunk_limit:
                        chunks.append({"chunk": number, "bytes": len(first),
                                       "records": count(first),
                                       **({} if second is None else {"mate_bytes": len(second),
                                                                     "mate_records": count(second)})})
                    if folder is not None:
                        names = [(first, left)] + ([] if second is None else [(second, right)])
                        for k, (data, origin) in enumerate(names, 1):
                            target = folder / f"chunk_{number}_{k}{suffix(origin)}"
                            if origin.suffix in (".gz", ".bz2", ".xz"):
                                with xopen(target, "wb") as out:
                                    out.write(bytes(data))
                            else:
                                target.write_bytes(bytes(data))
                            files.append(str(target))
            finally:
                if two_cm is not None:
                    two_cm.close()
    except (module.UnknownFileFormat, module.FileFormatError, ValueError) as exc:
        raise ConfigurationError(f"dnaio could not split the file: {exc}") from exc
    return {"status": "completed", "backend": "dnaio", "version": module.__version__,
            "paired": right is not None, "chunk_count": total, "chunks": chunks,
            "truncated": total > len(chunks), "outputs": files}


def dnaio_reverse_complement_reads(
    source: str | Path,
    output: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Reverse-complement every read of a FASTA/FASTQ file (dnaio SequenceRecord.reverse_complement).

    Names stay unchanged, the sequence is reverse complemented and the quality string is reversed. The format
    of ``output`` follows the input record type (FASTQ when the input has qualities).
    """

    module = _dnaio()
    path = _h.existing(source, "source")
    target = _h.new_file(output, "output")
    count = 0
    try:
        with module.open(path) as reader:
            first = next(iter(reader), None)
        if first is None:
            raise ConfigurationError("source contains no records.")
        fmt = "fastq" if first.qualities is not None else "fasta"
        with module.open(path) as reader, module.open(target, mode="w", fileformat=fmt) as writer:
            for record in reader:
                writer.write(record.reverse_complement())
                count += 1
    except ConfigurationError:
        target.unlink(missing_ok=True)
        raise
    except (module.FileFormatError, ValueError) as exc:
        target.unlink(missing_ok=True)
        raise ConfigurationError(f"dnaio could not process the file: {exc}") from exc
    return {"status": "completed", "backend": "dnaio", "version": module.__version__, "file_format": fmt,
            "record_count": count, "outputs": [str(target)]}


__all__ = [
    "dnaio_read_parallel_files",
    "dnaio_reverse_complement_reads",
    "dnaio_split_into_chunks",
    "dnaio_write_reads",
]

_register_public_verbose(globals(), tuple(__all__))
