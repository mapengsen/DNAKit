"""Binds the public ``<operation>_<format>`` functions to the format engines."""

# ruff: noqa: E501

from __future__ import annotations

import os
import re
from collections.abc import Mapping
from pathlib import Path

from dnakit._candidate_runtime import candidate_not_implemented, candidate_ok
from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.io_formats import _impl_download
from dnakit.io_formats._impl_align import FastaAlignmentFormat, NexusFormat, PhylipFormat
from dnakit.io_formats._impl_common import (
    COMMON_OPTIONS,
    DEFAULT_MAX_RECORDS,
    Format,
    Opts,
    Payload,
    bad_input,
    bad_option,
    emit,
    emit_parts,
    looks_like_path,
    page,
    read_source,
    sources_list,
)
from dnakit.io_formats._impl_tabular import MapFormat, PafFormat, PedFormat
from dnakit.io_formats._impl_trees import NewickFormat
from dnakit.io_formats._impl_twobit import TwoBitFormat
from dnakit.io_formats._impl_zarr import OptionalBackendMissing, ZarrFormat

MODULE = "io_formats"
FORMATS: dict[str, Format] = {
    f.name: f
    for f in (
        MapFormat(),
        PedFormat(),
        PafFormat(),
        NewickFormat(),
        NexusFormat(),
        PhylipFormat(),
        TwoBitFormat(),
        ZarrFormat(),
        FastaAlignmentFormat(),
    )
}
PUBLIC_FORMATS = ("map", "newick", "nexus", "paf", "ped", "phylip", "twobit", "zarr")
GENERIC_OPERATIONS = ("convert", "filter", "index", "merge", "sort", "stream", "summarize")
OPERATIONS = (
    "convert",
    "filter",
    "index",
    "merge",
    "read",
    "sort",
    "split",
    "stream",
    "summarize",
    "validate",
    "write",
)

_SUFFIX_FORMAT = {
    ".map": "map",
    ".ped": "ped",
    ".paf": "paf",
    ".nwk": "newick",
    ".newick": "newick",
    ".tre": "newick",
    ".tree": "newick",
    ".nex": "nexus",
    ".nexus": "nexus",
    ".nxs": "nexus",
    ".phy": "phylip",
    ".phylip": "phylip",
    ".2bit": "twobit",
    ".zarr": "zarr",
    ".fa": "fasta",
    ".fasta": "fasta",
    ".fna": "fasta",
}


def allowed_options(fmt: Format, op: str, generic: bool = False) -> dict[str, str]:
    """Option keys (with descriptions) accepted by ``op`` on ``fmt`` - the single source for validation and docstrings."""

    merged = dict(COMMON_OPTIONS.get(op, {}))
    if op in fmt.native_ops:  # these engines read lazily and apply their own limits
        merged.pop("max_bytes", None)
        merged.pop("max_records", None)
    merged.update(fmt.options.get(op, {}))
    if op == "convert":
        merged["to"] = "target format name (required): " + ", ".join(fmt.targets)
    if generic:
        merged["format"] = (
            "name of the input format ("
            + ", ".join(PUBLIC_FORMATS)
            + "); detected from the suffix/content when omitted"
        )
    return merged


def _head(source: object) -> bytes:
    """A short prefix of the content, to decide the format of text/bytes inputs."""

    if isinstance(source, (bytes, bytearray)):
        return bytes(source[:65536])
    if isinstance(source, os.PathLike) or (
        isinstance(source, str) and looks_like_path(source) and os.path.isfile(source)
    ):
        with open(source, "rb") as handle:  # type: ignore[arg-type]
            return handle.read(65536)
    if isinstance(source, str):
        return source[:65536].encode("utf-8")
    return b""


def detect_format(source: object) -> str:
    if isinstance(source, (os.PathLike, str)) and (
        not isinstance(source, str) or looks_like_path(source)
    ):
        path = Path(source).expanduser()
        if path.is_dir():
            return "zarr"
        if path.suffix.lower() == ".gz":
            path = path.with_suffix("")
        if path.suffix.lower() in _SUFFIX_FORMAT and (
            isinstance(source, os.PathLike) or os.path.exists(source)
        ):
            return _SUFFIX_FORMAT[path.suffix.lower()]
    head = _head(source)
    if head[:2] == b"\x1f\x8b":
        import zlib

        try:
            head = zlib.decompressobj(31).decompress(head, 65536)
        except zlib.error:
            head = b""
    if head[:4] in (b"\x43\x27\x41\x1a", b"\x1a\x41\x27\x43"):
        return "twobit"
    text = head.decode("utf-8", errors="ignore").lstrip("﻿")
    stripped = text.lstrip()
    if stripped[:6].lower() == "#nexus":
        return "nexus"
    if stripped.startswith("("):
        return "newick"
    if stripped.startswith(">"):
        return "fasta"
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if lines:
        first = lines[0]
        tabs = first.split("\t")
        if len(tabs) >= 12 and all(t.lstrip("+-").isdigit() for t in (tabs[1], tabs[2], tabs[3])):
            return "paf"
        if len(lines) > 1 and re.match(r"\s*\d+\s+\d+(\s+[A-Za-z]+)?\s*$", first):
            return "phylip"
        fields = first.split()
        if len(fields) == 4 and fields[3].lstrip("+-").isdigit():
            return "map"
        if len(fields) >= 6 and len(fields) % 2 == 0 and fields[4].lstrip("+-").isdigit():
            return "ped"
    raise bad_input(
        "Cannot detect the file format.",
        code="IO_FORMATS_UNKNOWN_FORMAT",
        hint=f"Pass options={{'format': <one of {', '.join(PUBLIC_FORMATS)}>}}.",
    )


def parse_name(function_name: str) -> tuple[str, str | None]:
    if function_name == "resume_download":
        return function_name, None
    for fmt in PUBLIC_FORMATS:
        if function_name.endswith("_" + fmt):
            return function_name[: -len(fmt) - 1], fmt
    return function_name, None


def _execute(
    function_name: str, input: object, options: Mapping[str, object] | None
) -> dict[str, object]:
    op, fmt_name = parse_name(function_name)
    where = f"io_formats.{function_name}"
    if op == "resume_download":
        opts = Opts(options, list(_impl_download.OPTIONS), where)
        opts.finish()
        return _impl_download.resume_download(input, opts)
    generic = fmt_name is None
    if generic:
        explicit = (options or {}).get("format") if isinstance(options, Mapping) else None
        if explicit is not None:
            if not isinstance(explicit, str) or explicit.lower() not in FORMATS:
                raise bad_option(f"{where}: format must be one of {list(PUBLIC_FORMATS)}.")
            fmt_name = explicit.lower()
        else:
            probe = sources_list(input)[0] if op == "merge" else input
            fmt_name = detect_format(probe)
        if fmt_name == "fasta" and op != "convert":
            raise bad_option(
                f"{where}: FASTA is only accepted as the source of convert.",
                hint="Use convert(format='fasta', to=...).",
            )
    fmt = FORMATS[fmt_name]  # type: ignore[index]
    opts = Opts(
        options,
        [*allowed_options(fmt, op, generic), *fmt.rejected_options.get(op, ())],
        where,
    )
    opts.finish()
    opts.text("format")
    native = op in fmt.native_ops
    result = _run(op, fmt, input, opts, native)
    return result


def _run(op: str, fmt: Format, input: object, opts: Opts, native: bool) -> dict[str, object]:
    if op == "read":
        if native:
            return {"format": fmt.name, **fmt.read_op(input, opts)}  # type: ignore[attr-defined]
        doc = fmt.coerce(input, opts)
        records = fmt.records(doc)
        limit = opts.integer("max_records", DEFAULT_MAX_RECORDS, 1) or DEFAULT_MAX_RECORDS
        if len(records) > limit:
            raise bad_option(
                f"The input has {len(records)} records, above max_records={limit}.",
                hint=f"Raise max_records or page with stream_{fmt.name}.",
            )
        return {
            "format": fmt.name,
            "record_count": len(records),
            **fmt.describe(doc),
            "records": records,
        }
    if op == "write":
        if native:
            return emit(fmt.write_op(input, opts), opts)  # type: ignore[attr-defined]
        doc = fmt.coerce_for_write(input, opts)
        return emit(fmt.dump(doc, opts), opts)
    if op == "validate":
        if native:
            return {"format": fmt.name, **fmt.validate_op(input, opts)}  # type: ignore[attr-defined]
        errors: list[str] = []
        warnings: list[str] = []
        stats: dict[str, object] = {}
        try:
            doc = fmt.coerce(input, opts, errors)
        except InputFormatError as exc:
            if isinstance(input, (Mapping, list, tuple)):
                errors.append(exc.message)
                doc = None
            else:
                raise
        if doc is not None and not errors:
            stats = fmt.validate_doc(doc, opts, errors, warnings)
        return {
            "format": fmt.name,
            "valid": not errors,
            "errors": errors,
            "warnings": warnings,
            **stats,
        }
    if op == "stream":
        if native:
            return {"format": fmt.name, **fmt.stream_op(input, opts)}  # type: ignore[attr-defined]
        records = fmt.records(fmt.coerce(input, opts))
        return {"format": fmt.name, **page(records, opts)}
    if op == "summarize":
        if native:
            return {"format": fmt.name, **fmt.summarize_op(input, opts)}  # type: ignore[attr-defined]
        return {"format": fmt.name, **fmt.summarize(fmt.coerce(input, opts), opts)}
    if op == "index":
        if native:
            return {"format": fmt.name, **fmt.index_op(input, opts)}  # type: ignore[attr-defined]
        return {"format": fmt.name, **fmt.index(fmt.coerce(input, opts), opts)}
    if op == "merge":
        sources = sources_list(input)
        if native:
            return emit(fmt.merge_op(sources, opts), opts)  # type: ignore[attr-defined]
        docs = [fmt.coerce(s, opts) for s in sources]
        payload = fmt.dump(fmt.merge(docs, opts), opts)
        payload.info["input_count"] = len(docs)
        return emit(payload, opts)
    if op == "split":
        if native:
            parts = fmt.split_op(input, opts)  # type: ignore[attr-defined]
        else:
            doc = fmt.coerce(input, opts)
            parts = [(label, fmt.dump(part, opts)) for label, part in fmt.split(doc, opts)]
        return emit_parts(fmt.name, parts, fmt.suffix, opts)
    if op in ("filter", "sort"):
        if native:
            return emit(getattr(fmt, f"{op}_op")(input, opts), opts)
        doc = fmt.coerce(input, opts)
        new = getattr(fmt, op)(doc, opts)
        payload = fmt.dump(new, opts)
        before, after = fmt.count(doc), fmt.count(new)
        payload.info.update({"input_count": before, "output_count": after})
        if op == "filter":
            payload.info["removed_count"] = before - after
        if isinstance(new, dict) and "kept_variant_indices" in new:
            payload.info["kept_variant_indices"] = new["kept_variant_indices"]
        return emit(payload, opts)
    if op == "convert":
        target = str(opts.text("to") or "").lower()
        if not target:
            raise bad_option(f"convert needs options['to']: one of {list(fmt.targets)}.")
        if target not in fmt.targets:
            raise bad_option(
                f"Cannot convert {fmt.name} to {target!r}; supported targets: {list(fmt.targets)}."
            )
        if native:
            return emit(fmt.convert_op(input, opts), opts)  # type: ignore[attr-defined]
        doc = fmt.coerce(input, opts)
        payload = fmt.convert(doc, target, opts)
        payload.info.setdefault("source_format", fmt.name)
        return emit(payload, opts)
    raise ConfigurationError(f"Unknown io_formats operation {op!r}.")


def run(
    function_name: str, input: object, options: Mapping[str, object] | None
) -> dict[str, object]:
    """Execute one public ``io_formats`` function and wrap the result in the candidate envelope."""

    try:
        result = _execute(function_name, input, options)
    except OptionalBackendMissing as exc:
        return candidate_not_implemented(
            MODULE,
            function_name,
            input,
            options=options,
            requires=exc.requires,
            message=exc.message,
        )
    # operations that are not validations report their caveats both in the result and in the envelope
    notes = result.get("warnings") if not function_name.startswith("validate") else None
    return candidate_ok(  # type: ignore[return-value]
        MODULE,
        function_name,
        input,
        result,
        options=options,
        warnings=notes if isinstance(notes, list) else (),
    )


__all__ = [
    "FORMATS",
    "OPERATIONS",
    "PUBLIC_FORMATS",
    "Payload",
    "allowed_options",
    "detect_format",
    "read_source",
    "run",
]
