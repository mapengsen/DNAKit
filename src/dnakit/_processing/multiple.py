"""Explicit file-order concatenation and row-paired crossover."""

from __future__ import annotations

from collections.abc import Iterable
from contextlib import ExitStack, closing
from pathlib import Path
from typing import Any

from dnakit._parallel import invoke, ordered_map
from dnakit._verbosity import progress_enabled, status_enabled
from dnakit.core import DNARecord
from dnakit.exceptions import ConfigurationError
from dnakit.io import write

from .inputs import iter_inputs, output_path, paths_in, protect_inputs, table_config
from .sequence import as_record


def _crossover_job(payload: Any) -> DNARecord:
    left, call = payload
    return as_record(invoke(call), left)


def run_multiple(
    function: Any, bound: Any, input_name: str, mode: str, options: dict[str, Any]
) -> Any:
    from rich.console import Console

    names = tuple(bound.signature.parameters)
    sources = paths_in(bound.arguments[input_name])
    if mode == "pair":
        sources += paths_in(bound.arguments[names[1]])
    target = (
        Path(options.get("output_path") or output_path(sources[0], function.__name__))
        .expanduser()
        .resolve()
    )
    protect_inputs(sources, target)
    config = table_config(
        sources,
        target,
        options.get("write_config"),
        options.get("read_config"),
        options.get("sequence_column"),
    )
    console = Console(stderr=True)
    show = options.get("show_progress", True)
    if not isinstance(show, bool):
        raise ConfigurationError("show_progress must be boolean.")
    show = progress_enabled(show)
    if status_enabled(show):
        console.print(f"Starting {function.__name__}: {sources[0]}", markup=False, soft_wrap=True)
    with ExitStack() as stack:
        first = stack.enter_context(
            closing(
                iter_inputs(
                    bound.arguments[input_name],
                    read_config=options.get("read_config"),
                    sequence_column=options.get("sequence_column"),
                )
            )
        )
        products: Iterable[DNARecord]
        if mode == "concat":
            source_ids = []

            def sequences() -> Any:
                for record in first:
                    if record.features or record.letter_annotations:
                        raise ConfigurationError(
                            "Concatenation of annotated file records "
                            "requires an explicit annotation mapping."
                        )
                    source_ids.append(record.id)
                    yield record.sequence

            bound.arguments[input_name] = sequences()
            result = function(*bound.args, **bound.kwargs)
            sequence = getattr(result, "sequence", result)
            products = [DNARecord(sequence, function.__name__, metadata={"source_ids": source_ids})]
        else:
            second = stack.enter_context(
                closing(
                    iter_inputs(
                        bound.arguments[names[1]],
                        read_config=options.get("read_config"),
                        sequence_column=options.get("sequence_column"),
                    )
                )
            )

            def jobs() -> Any:
                for left, right in zip(first, second, strict=True):
                    if right.features or right.letter_annotations:
                        raise ConfigurationError(
                            "Crossover cannot automatically map parental annotations."
                        )
                    bound.arguments[input_name], bound.arguments[names[1]] = (
                        left.sequence,
                        right.sequence,
                    )
                    yield (
                        left,
                        (function.__module__, function.__name__, bound.args, bound.kwargs, False),
                    )

            products = ordered_map(
                _crossover_job,
                jobs(),
                workers=1 if bound.arguments.get("rng") is not None else None,
            )
        result = write(products, target, config=config)
    if status_enabled(show):
        console.print(
            f"{function.__name__} complete: {result.record_count} records -> {target}",
            markup=False,
            soft_wrap=True,
        )
    return result
