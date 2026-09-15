"""Stream sequence operations through the shared file input layer."""

from __future__ import annotations

import csv
import gzip
import io
import json
import os
import tempfile
from collections.abc import Callable, Generator, Iterator, Mapping
from contextlib import closing
from dataclasses import fields, is_dataclass
from enum import Enum
from inspect import BoundArguments
from pathlib import Path
from typing import Any, cast

from rich.console import Console
from rich.progress import BarColumn, SpinnerColumn, TextColumn, TimeElapsedColumn

from dnakit._parallel import invoke as invoke_parallel
from dnakit._parallel import ordered_map, worker_scope
from dnakit._progress import Progress
from dnakit._verbosity import progress_enabled, status_enabled
from dnakit.core import DNA, DNARecord, DNASequence
from dnakit.core._json import to_json_compatible
from dnakit.exceptions import ConfigurationError
from dnakit.io import WriteConfig, WriteResult, write
from dnakit.io.api import _artifact, _infer_format

from .inputs import input_config, iter_inputs, output_path, paths_in, protect_inputs, table_config


def _record_product(payload: Any) -> Any:
    record, call = payload
    return record, invoke_parallel(call)


class _ValidationEncoder(json.JSONEncoder):
    """Serialize generated positional reports without a metadata-sized node budget."""

    def default(self, value: Any) -> Any:
        if isinstance(value, Enum):
            return value.value
        if is_dataclass(value) and not isinstance(value, type):
            return {field.name: getattr(value, field.name) for field in fields(value)}
        if isinstance(value, Mapping):
            return {str(key): item for key, item in value.items()}
        if isinstance(value, frozenset):
            return to_json_compatible(value)
        return super().default(value)


def as_record(value: Any, source: DNARecord) -> DNARecord:
    if isinstance(value, DNARecord):
        return value
    if isinstance(value, DNA):
        return value.record
    if hasattr(value, "to_record"):
        chunk = value.to_record()
        mapped = None
        if source.features or source.letter_annotations:
            from dnakit.ops.records import trim_record

            mapped = trim_record(
                source,
                left=value.source_start,
                right=len(source.sequence) - value.source_end,
                feature_policy="truncate",
            ).record
        return DNARecord(
            chunk.sequence,
            chunk.id,
            description=chunk.description,
            metadata={**source.metadata, **chunk.metadata},
            features=mapped.features if mapped is not None else (),
            letter_annotations=mapped.letter_annotations if mapped is not None else None,
        )
    if isinstance(getattr(value, "record", None), DNARecord):
        return cast(DNARecord, value.record)
    sequence = value if isinstance(value, DNASequence) else getattr(value, "sequence", None)
    if not isinstance(sequence, DNASequence):
        raise ConfigurationError("This operation did not produce a writable DNA sequence.")
    if source.features or source.letter_annotations:
        raise ConfigurationError(
            "This sequence-only operation cannot safely remap annotations; "
            "use a record-aware operation."
        )
    return DNARecord(sequence, source.id, description=source.description, metadata=source.metadata)


def write_text_results(
    rows: Iterator[tuple[DNARecord, Any]],
    target: Path,
    *,
    config: WriteConfig,
    report: bool = False,
    _fragment: bool = False,
) -> WriteResult:
    """Atomically write RNA/protein text or report JSON without a DNA alphabet coercion."""
    fmt = _infer_format(str(target), None)
    if fmt not in ({"jsonl"} if report else {"fasta", "text", "csv", "tsv", "jsonl"}):
        raise ConfigurationError(
            "RNA/protein output supports FASTA, TXT, CSV/TSV or JSONL; reports use JSONL."
        )
    if target.exists() and not config.overwrite:
        raise FileExistsError(f"Refusing to overwrite existing output: {target}")
    if config.create_parents:
        target.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(prefix=f".{target.name}.", suffix=".tmp", dir=target.parent)
    os.close(descriptor)
    temporary = Path(name)
    count = byte_count = 0
    compressed = config.compression == "gzip" or (
        config.compression == "auto" and target.suffix.lower() == ".gz"
    )
    try:
        opener = gzip.open if compressed else open
        open_options: dict[str, Any] = {"encoding": config.encoding, "newline": ""}
        if compressed:
            open_options["compresslevel"] = config.compression_level
        with opener(temporary, "wt", **open_options) as handle:

            class BoundedWriter:
                def write(self, text: str) -> int:
                    nonlocal byte_count
                    byte_count += len(text.encode(config.encoding))
                    if config.max_output_bytes is not None and byte_count > config.max_output_bytes:
                        raise ConfigurationError("Output exceeds max_output_bytes.")
                    return cast(int, handle.write(text))

            stream = BoundedWriter()
            columns = config.table_columns or (
                config.id_column,
                config.sequence_column,
                config.metadata_column,
            )
            writer = csv.DictWriter(
                stream,
                fieldnames=columns,
                delimiter=config.delimiter or ("," if fmt == "csv" else "\t"),
                lineterminator="\n",
            )
            if fmt in {"csv", "tsv"}:
                writer.writeheader()
            for record, value in rows:
                if report:
                    from dnakit.standardize.results import ValidationReport

                    if isinstance(value, ValidationReport):
                        stream.write(
                            json.dumps(
                                {"id": record.id, "result": value},
                                ensure_ascii=False,
                                cls=_ValidationEncoder,
                            )
                            + "\n"
                        )
                        count += 1
                        continue
                    payload = (
                        value.to_dict() if hasattr(value, "to_dict") else to_json_compatible(value)
                    )
                    stream.write(
                        json.dumps({"id": record.id, "result": payload}, ensure_ascii=False) + "\n"
                    )
                elif fmt == "fasta":
                    text = str(value)
                    stream.write(
                        f">{record.id}{' ' + record.description if record.description else ''}\n"
                    )
                    block = max(config.line_width, 1024**2 // config.line_width * config.line_width)
                    for start in range(0, len(text), block):
                        end = min(len(text), start + block)
                        stream.write(
                            "\n".join(
                                text[offset : min(end, offset + config.line_width)]
                                for offset in range(start, end, config.line_width)
                            )
                            + "\n"
                        )
                elif fmt == "text":
                    stream.write(str(value) + "\n")
                else:
                    metadata = cast(dict[str, Any], to_json_compatible(record.metadata))
                    row = {
                        "id": record.id,
                        "description": record.description,
                        "sequence": str(value),
                        "metadata": metadata,
                    }
                    if fmt == "jsonl":
                        stream.write(json.dumps(row, ensure_ascii=False) + "\n")
                    else:
                        aliases = {
                            config.id_column: record.id,
                            config.sequence_column: str(value),
                            config.description_column: record.description,
                            config.metadata_column: metadata,
                        }
                        row = {**metadata, **aliases}
                        # DNA coordinates and per-base annotations do not describe RNA/protein.
                        for column in (
                            config.parts_column,
                            config.features_column,
                            config.letter_annotations_column,
                            "alphabet",
                            "topology",
                            "strandedness",
                            "schema_version",
                        ):
                            if column in columns:
                                raise ConfigurationError(
                                    f"Remove DNA-specific output column {column!r} "
                                    "for RNA/protein output."
                                )
                        writer.writerow(
                            {
                                key: json.dumps(row[key], ensure_ascii=False)
                                if isinstance(row.get(key), (dict, list))
                                else row.get(key, "")
                                for key in columns
                            }
                        )
                count += 1
        if config.overwrite:
            os.replace(temporary, target)
        else:
            os.link(temporary, target)
            temporary.unlink()
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    if _fragment:
        return WriteResult(fmt, count)
    artifact = _artifact(target, fmt)
    return WriteResult(fmt, count, byte_count=artifact.byte_size, target_artifact=artifact)


def run_sequence(
    function: Callable[..., Any],
    bound: BoundArguments,
    input_name: str,
    mode: str,
    options: dict[str, Any],
    *,
    _region_text: str | None = None,
) -> WriteResult:
    source_value = bound.arguments[input_name]
    sources = paths_in(source_value)
    if not sources:
        raise ConfigurationError("A file operation requires at least one source file.")
    target = options.get("output_path")
    if target is None:
        suffix = ".jsonl" if mode == "report" else None
        if mode == "text" and _infer_format(str(sources[0]), None) not in {"text", "csv", "tsv"}:
            suffix = ".fa"
        target = output_path(sources[0], function.__name__, suffix)
    target = Path(target).expanduser().resolve()
    protect_inputs(sources, target)
    if _region_text is None:
        if (
            function.__name__ in {"reverse", "complement"}
            and len(sources) == 1
            and isinstance(source_value, (str, os.PathLike))
            and options.get("topology") is None
            and isinstance(options.get("show_progress", True), bool)
        ):
            from dnakit._parallel import current_workers
            from dnakit.ops._fasta_direction import try_reverse_fasta

            direction = try_reverse_fasta(
                sources[0],
                target,
                read_config=options.get("read_config"),
                write_config=options.get("write_config"),
                sequence_column=options.get("sequence_column"),
                workers=current_workers(),
                show_progress=options.get("show_progress", True),
                reverse=function.__name__ == "reverse",
                complement=function.__name__ == "complement",
            )
            if direction is not None:
                return direction
        from .regions import try_regions

        fast = try_regions(function, bound, input_name, mode, options, sources, target)
        if fast is not None:
            return fast
    reader_config = options.get("read_config")
    if options.get("topology") is not None:
        from dataclasses import replace

        reader_config = replace(
            input_config(sources[0], reader_config), topology=options["topology"]
        )
    column = options.get("sequence_column")
    config = table_config(sources, target, options.get("write_config"), reader_config, column)
    if mode == "text" and config.table_columns is not None:
        from dataclasses import replace

        explicit = options.get("write_config")
        if explicit is None or explicit.table_columns is None:
            dna_columns = {
                config.parts_column,
                config.features_column,
                config.letter_annotations_column,
                "alphabet",
                "topology",
                "strandedness",
                "schema_version",
            }
            config = replace(
                config, table_columns=tuple(c for c in config.table_columns if c not in dna_columns)
            )
    if (
        mode == "chunks"
        and config.table_columns is not None
        and config.metadata_column not in config.table_columns
    ):
        from dataclasses import replace

        config = replace(config, table_columns=(*config.table_columns, config.metadata_column))
    show = options.get("show_progress", True)
    if not isinstance(show, bool):
        raise ConfigurationError("show_progress must be boolean.")
    show = progress_enabled(show)
    console = Console(stderr=True)
    if status_enabled(show):
        console.print(
            f"Starting {function.__name__.replace('_', ' ')}: {sources[0]}",
            markup=False,
            soft_wrap=True,
        )

    def invoke(value: Any) -> Any:
        bound.arguments[input_name] = value
        return function(*bound.args, **bound.kwargs)

    def parallel_call(value: Any) -> Any:
        bound.arguments[input_name] = value
        return function.__module__, function.__name__, bound.args, bound.kwargs, False

    invoke.parallel_call = parallel_call  # type: ignore[attr-defined]
    normalizer = invoke if mode == "normalize" else None
    records: Generator[DNARecord, None, None]
    if _region_text is None:
        records = iter_inputs(
            source_value, read_config=reader_config, sequence_column=column, normalize=normalizer
        )
    else:
        from dnakit.io._formats import iter_fasta

        from .inputs import _raw_records

        selected = input_config(sources[0], reader_config)
        selected_records = (
            _raw_records(sources[0], selected, normalizer, text=_region_text)
            if normalizer is not None
            else iter_fasta(io.StringIO(_region_text, newline=""), selected)
        )
        records = (record for record in selected_records)
    with (
        worker_scope(options.get("workers", -1)),
        closing(records),
        Progress(
            SpinnerColumn(),
            TextColumn("{task.description}: {task.completed:.0f} records"),
            BarColumn(),
            TimeElapsedColumn(),
            console=console,
            disable=not show,
        ) as display,
    ):
        task = display.add_task(function.__name__, total=None)

        def products() -> Iterator[tuple[DNARecord, Any]]:
            if mode in {"normalize", "chunks"}:
                # Normalization parallelizes raw rows before strict DNA construction.
                # Chunk iterators stay lazy, including chromosome-sized records.
                for record in records:
                    if mode == "normalize":
                        yield record, record
                    else:
                        for chunk in invoke(record):
                            yield record, chunk
                    display.advance(task)
                return

            def jobs() -> Iterator[Any]:
                for record in records:
                    if mode in {"sequence", "text"} and (
                        record.features or record.letter_annotations
                    ):
                        raise ConfigurationError(
                            "This sequence-only operation cannot safely remap annotated input."
                        )
                    bound.arguments[input_name] = (
                        record if mode in {"record", "report"} else record.sequence
                    )
                    yield (
                        record,
                        (function.__module__, function.__name__, bound.args, bound.kwargs, False),
                    )

            cheap = function.__name__ in {
                "reverse",
                "complement",
                "reverse_complement",
                "trim",
                "trim_record",
                "subsequence",
                "insert",
                "delete",
                "substitute",
                "mask",
                "transcribe",
            }
            for record, result in ordered_map(
                _record_product,
                jobs(),
                workers=1
                if cheap or bound.arguments.get("rng") is not None
                else options.get("workers"),
            ):
                yield record, result
                display.advance(task)

        if mode in {"text", "report"}:
            result = write_text_results(
                products(),
                target,
                config=config,
                report=mode == "report",
                _fragment=_region_text is not None,
            )
        elif _region_text is not None:
            from .fragments import write_fragment

            result = write_fragment(
                (as_record(value, record) for record, value in products()), target, config
            )
        else:
            result = write(
                (as_record(value, record) for record, value in products()), target, config=config
            )
    if status_enabled(show):
        console.print(
            f"{function.__name__.replace('_', ' ').capitalize()} complete: "
            f"{result.record_count} records -> {target}",
            markup=False,
            soft_wrap=True,
        )
    return result
