"""Lightweight public-entry adapters; heavy imports stay outside core initialization."""

from __future__ import annotations

from collections.abc import Callable
from functools import wraps
from inspect import Parameter, signature
from pathlib import Path
from typing import Any

_OPTIONS = {
    "output_path": None,
    "read_config": None,
    "write_config": None,
    "sequence_column": None,
    "show_progress": True,
    "topology": None,
    "workers": -1,
}


def accept_inputs(function: Callable[..., Any], *, mode: str = "sequence") -> Any:
    """Add shared file keywords without changing the native object's call semantics."""
    original = signature(function)
    input_name = next(iter(original.parameters))

    def dispatch(*args: Any, **kwargs: Any) -> Any:
        from .inputs import paths_in

        options = {
            key: kwargs.pop(key)
            for key in _OPTIONS
            if key in kwargs and key not in original.parameters
        }
        bound = original.bind(*args, **kwargs)
        value = bound.arguments[input_name]
        file_input = bool(paths_in(value))
        if mode == "normalize" and isinstance(value, str):
            # Punctuation is legitimate raw normalization input. Only file-shaped strings dispatch.
            suffix = Path(value).suffix.lower() if len(value) < 4096 and "\n" not in value else ""
            file_input = suffix in {".fna", ".fa", ".fasta", ".txt", ".csv", ".tsv", ".gz"}
        collection = mode in {"collection", "set", "mapping"}
        if collection and (
            file_input
            or isinstance(value, str)
            or (isinstance(value, (list, tuple)) and any(isinstance(v, str) for v in value))
        ):
            from .datasets import run_collection

            return run_collection(function, bound, input_name, mode, options)
        if file_input or (mode == "pair" and any(paths_in(v) for v in bound.arguments.values())):
            if mode in {"concat", "pair"}:
                from .multiple import run_multiple

                return run_multiple(function, bound, input_name, mode, options)
            from .sequence import run_sequence

            return run_sequence(
                function, bound, input_name, "text" if mode == "text_native" else mode, options
            )
        if options.get("output_path") is not None:
            from dnakit.exceptions import ConfigurationError

            raise ConfigurationError("output_path requires file input for this operation.")
        if isinstance(value, str) and mode not in {"normalize", "text_native"}:
            from dnakit.core import DNARecord, DNASequence

            sequence = DNASequence(
                "".join(value.split()).upper(),
                alphabet="iupac",
                topology=options.get("topology") or "linear",
            )
            bound.arguments[input_name] = (
                DNARecord(sequence, "sequence_1")
                if mode == "record" and function.__name__.endswith("_record")
                else sequence
            )
        if mode == "pair":
            from dnakit.core import DNASequence

            second = tuple(original.parameters)[1]
            if isinstance(bound.arguments.get(second), str):
                bound.arguments[second] = DNASequence(
                    "".join(bound.arguments[second].split()).upper(), alphabet="iupac"
                )
        return function(*bound.args, **bound.kwargs)

    @wraps(function)
    def wrapped(*args: Any, **kwargs: Any) -> Any:
        from dnakit._parallel import _workers, worker_scope

        with worker_scope(kwargs.get("workers", _workers.get())):
            return dispatch(*args, **kwargs)

    parameters = list(original.parameters.values())
    for name, default in _OPTIONS.items():
        if name not in original.parameters:
            parameters.append(Parameter(name, Parameter.KEYWORD_ONLY, default=default))
    wrapped.__signature__ = original.replace(parameters=parameters)  # type: ignore[attr-defined]
    # Python 3.14 wraps copies __annotate__, which may be unset after a prior
    # decorator assigned __annotations__. Preserve those explicit annotations too.
    wrapped.__annotations__ = dict(getattr(function, "__annotations__", {}))
    if mode in {"collection", "set", "mapping"}:
        # The native annotation still describes records, while this wrapper
        # also accepts files. Expose that contract to JSON tool adapters.
        wrapped.__dnakit_file_input__ = input_name  # type: ignore[attr-defined]
    return wrapped


def register(namespace: dict[str, Any], names: tuple[str, ...], *, mode: str) -> None:
    for name in names:
        namespace[name] = accept_inputs(namespace[name], mode=mode)
