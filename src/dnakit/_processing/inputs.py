"""One interpretation of paths, sequence text, tables and record collections."""

from __future__ import annotations

import csv
import io
import re
from collections.abc import Callable, Generator, Iterable, Iterator, Mapping
from dataclasses import replace
from os import PathLike
from pathlib import Path
from typing import Any

from dnakit.core import DNA, DNARecord, DNASequence, DNASet
from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.io import ReadConfig, WriteConfig, read
from dnakit.io._formats import _header, _record_from_mapping
from dnakit.io.api import _infer_format, _open_read_text


def is_path(value: object) -> bool:
    return isinstance(value, PathLike) or (
        isinstance(value, str) and any(marker in value for marker in ("/", "\\", "."))
    )


def paths_in(value: object) -> tuple[Path, ...]:
    if is_path(value):
        return (Path(value).expanduser().resolve(),)  # type: ignore[arg-type]
    if isinstance(value, (list, tuple)):
        return tuple(path for item in value for path in paths_in(item))
    if isinstance(value, Mapping):
        return tuple(path for item in value.values() for path in paths_in(item))
    return ()


def input_config(
    path: Path, config: ReadConfig | None, sequence_column: str | None = None
) -> ReadConfig:
    if config is not None and not isinstance(config, ReadConfig):
        raise TypeError("read_config must be ReadConfig or None.")
    if config is None:
        size = max(1, path.stat().st_size)
        if path.suffix.lower() == ".gz":
            size = max(size, 20_000_000_000)
        config = ReadConfig(
            encoding="utf-8-sig",
            max_input_bytes=size + 3,
            max_sequence_symbols=size,
            max_record_lines=size,
            max_field_size=size,
        )
    return replace(
        config,
        preserve_extra_columns=True,
        sequence_column=config.sequence_column if sequence_column is None else sequence_column,
    )


def output_path(source: Path, operation: str, suffix: str | None = None) -> Path:
    base = source.with_suffix("") if source.suffix.lower() == ".gz" else source
    extension = base.suffix + (source.suffix if base != source else "")
    return source.with_name(base.stem + "_" + operation + (extension if suffix is None else suffix))


def protect_inputs(sources: Iterable[Path], target: Path) -> None:
    for source in sources:
        if target == source or (target.exists() and target.samefile(source)):
            raise ConfigurationError("Processing output must differ from every input file.")


def table_config(
    sources: tuple[Path, ...],
    target: Path,
    config: WriteConfig | None,
    read_config: ReadConfig | None,
    sequence_column: str | None,
) -> WriteConfig:
    resolved = WriteConfig() if config is None else config
    if not isinstance(resolved, WriteConfig):
        raise TypeError("write_config must be WriteConfig or None.")
    if _infer_format(str(target), None) not in {"csv", "tsv"} or resolved.table_columns is not None:
        return resolved
    columns: list[str] = []
    selected = sequence_column or (read_config.sequence_column if read_config else "sequence")
    id_column = read_config.id_column if read_config else "id"
    for source in sources:
        if _infer_format(str(source), None) not in {"csv", "tsv"}:
            continue
        reader_config = input_config(source, read_config, sequence_column)
        handle, close, _ = _open_read_text(source, reader_config)
        try:
            delimiter = reader_config.delimiter or (
                "," if _infer_format(str(source), None) == "csv" else "\t"
            )
            header = next(csv.reader(handle, delimiter=delimiter), [])
            if not header or len(header) != len(set(header)) or any(not name for name in header):
                raise InputFormatError("CSV/TSV headers must be non-empty and unique.")
            for column in header:
                if column not in columns:
                    columns.append(column)
        finally:
            close()
    if not columns:
        return resolved
    original = read_config or ReadConfig()
    return replace(
        resolved,
        table_columns=tuple(columns),
        sequence_column=selected,
        id_column=id_column,
        description_column=original.description_column,
        metadata_column=original.metadata_column,
        features_column=original.features_column,
        parts_column=original.parts_column,
        letter_annotations_column=original.letter_annotations_column,
    )


def _raw_rows(
    path: Path, config: ReadConfig, *, text: str | None = None
) -> Generator[Any, None, None]:
    """Normalize raw symbols before strict DNA construction, retaining record boundaries."""
    format_name = _infer_format(str(path), None)
    if text is None:
        handle, close, _ = _open_read_text(path, config)
    else:
        handle = io.StringIO(text, newline="")
        close = handle.close

    try:
        if format_name in {"csv", "tsv"}:
            previous = csv.field_size_limit(config.max_field_size)
            try:
                rows = csv.DictReader(
                    handle, delimiter=config.delimiter or ("," if format_name == "csv" else "\t")
                )
                header = rows.fieldnames
                if not header or len(header) != len(set(header)) or any(not key for key in header):
                    raise InputFormatError("CSV/TSV requires unique non-empty headers.")
                for index, row in enumerate(rows):
                    if (
                        index >= config.max_records
                        or None in row
                        or any(v is None for v in row.values())
                    ):
                        raise InputFormatError("Malformed table row or record limit exceeded.")
                    text = row.get(config.sequence_column)
                    if not isinstance(text, str):
                        raise InputFormatError(f"Missing sequence column: {config.sequence_column}")
                    yield text, None, None, row, index
            finally:
                csv.field_size_limit(previous)
        elif format_name in {"fasta", "text"}:
            identifier, description = "sequence_1", ""
            parts: list[str] = []
            count = symbols = lines = 0
            saw_header = False
            for number, raw in enumerate(handle, 1):
                line = raw.rstrip("\r\n")
                if len(line) > config.max_field_size:
                    raise InputFormatError("Sequence line exceeds max_field_size.")
                if format_name == "fasta" and line.startswith(">"):
                    if saw_header:
                        count += 1
                        if count > config.max_records:
                            raise InputFormatError("Input exceeds max_records.")
                        yield "".join(parts), identifier, description, None, count - 1
                    identifier, description = _header(line, line_number=number, marker="FASTA")
                    parts, symbols, lines, saw_header = [], 0, 1, True
                    continue
                if format_name == "fasta" and not saw_header and line.strip():
                    raise InputFormatError("FASTA sequence appeared before its header.")
                lines += 1
                if lines > config.max_record_lines:
                    raise InputFormatError("Input exceeds max_record_lines.")
                if not line.strip():
                    continue
                if format_name == "text" and config.text_mode == "lines":
                    if len(line) > config.max_sequence_symbols:
                        raise InputFormatError("Input exceeds max_sequence_symbols.")
                    count += 1
                    if count > config.max_records:
                        raise InputFormatError("Input exceeds max_records.")
                    yield line, f"sequence_{count}", "", None, count - 1
                    lines = 0
                else:
                    symbols += len(line)
                    if symbols > config.max_sequence_symbols or lines > config.max_record_lines:
                        raise InputFormatError("Input exceeds sequence record limits.")
                    parts.append(line)
            if parts or saw_header:
                if count >= config.max_records:
                    raise InputFormatError("Input exceeds max_records.")
                yield "".join(parts), identifier, description, None, count
        else:
            raise ConfigurationError(
                "Raw normalization supports FASTA/FNA, TXT, CSV and TSV inputs."
            )
    except (OSError, EOFError, UnicodeError, csv.Error) as exc:
        raise InputFormatError(
            "Input stream could not be read, decoded or decompressed.",
            code="INPUT_READ_FAILED",
            context={"source": str(path), "reason": str(exc)},
        ) from exc
    finally:
        close()


def _normalized_record(row: Any, config: ReadConfig, result: Any) -> DNARecord:
    _text, identifier, description, mapping, index = row
    sequence = result if isinstance(result, DNASequence) else result.sequence
    if sequence is None:
        raise InputFormatError("Normalization did not produce a representable DNA sequence.")
    if mapping is not None:
        mapping[config.sequence_column] = sequence.to_string()
        return _record_from_mapping(mapping, config, record_index=index)
    return DNARecord(sequence, identifier, description=description)


def _normalize_job(payload: Any) -> DNARecord:
    from dnakit._parallel import invoke

    row, config, call = payload
    module, name, _args, kwargs, _materialize = call
    if (
        module == "dnakit.standardize.normalize"
        and name == "normalize"
        and kwargs.get("config") is None
        and kwargs.get("keep_ambiguous", True) is True
        and kwargs.get("keep_u", False) is False
        and kwargs.get("keep_other", False) is False
    ):
        # File normalization exports sequences, not per-character edit audits.
        # Match the default ASCII case/deletion policies without constructing
        # millions of audit objects that this adapter immediately discards.
        text = row[0].translate(str.maketrans("acgtryswkmbdhvn", "ACGTRYSWKMBDHVN"))
        text = re.sub("[^ACGTRYSWKMBDHVN]", "", text)
        alphabet = "iupac" if re.search("[^ACGT]", text) else "strict"
        return _normalized_record(row, config, DNASequence(text, alphabet=alphabet))
    return _normalized_record(row, config, invoke(call))


def _raw_records(
    path: Path, config: ReadConfig, normalize: Callable[[str], Any], *, text: str | None = None
) -> Iterator[DNARecord]:
    from contextlib import closing

    from dnakit._parallel import ordered_map

    job_builder = getattr(normalize, "parallel_call", None)
    with closing(_raw_rows(path, config, text=text)) as rows:
        if job_builder is None:
            for row in rows:
                yield _normalized_record(row, config, normalize(row[0]))
        else:
            yield from ordered_map(
                _normalize_job, ((row, config, job_builder(row[0])) for row in rows)
            )


def iter_inputs(
    value: object,
    *,
    read_config: ReadConfig | None = None,
    sequence_column: str | None = None,
    normalize: Callable[[str], Any] | None = None,
) -> Generator[DNARecord, None, None]:
    """Stream records from files, text, DNA objects, or mixed record/file iterables."""
    if is_path(value):
        path = Path(value).expanduser().resolve()  # type: ignore[arg-type]
        config = input_config(path, read_config, sequence_column)
        if normalize is not None:
            yield from _raw_records(path, config, normalize)
        else:
            with read(path, config=config) as records:
                yield from records
    elif isinstance(value, DNARecord):
        yield value
    elif isinstance(value, (DNA, DNASet)):
        yield from value.records if isinstance(value, DNA) else value
    elif isinstance(value, DNASequence):
        yield DNARecord(value, "sequence_1")
    elif isinstance(value, str):
        yield DNARecord(DNASequence("".join(value.split()).upper(), alphabet="iupac"), "sequence_1")
    elif isinstance(value, Iterable) and not isinstance(value, (bytes, Mapping)):
        for index, item in enumerate(value, 1):
            if isinstance(item, str) and not is_path(item):
                yield DNARecord(
                    DNASequence("".join(item.split()).upper(), alphabet="iupac"),
                    f"sequence_{index}",
                )
            else:
                yield from iter_inputs(
                    item,
                    read_config=read_config,
                    sequence_column=sequence_column,
                    normalize=normalize,
                )
    else:
        raise TypeError(
            "Input must be DNA text, a file path, DNA records, or an iterable of these."
        )
