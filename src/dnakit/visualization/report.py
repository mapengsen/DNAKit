"""Self-contained bounded HTML reports for local, read-only exploration."""

from __future__ import annotations

import html
import json
from collections.abc import Iterable, Mapping
from itertools import islice

from dnakit.core import DNARecord, DNASet
from dnakit.core._json import to_json_compatible
from dnakit.exceptions import ConfigurationError

from .results import HTMLReportArtifact

_STYLE = """
body{font-family:system-ui,sans-serif;max-width:1100px;margin:2rem auto;padding:0 1rem;
color:#172033}
table{border-collapse:collapse;width:100%}th,td{border:1px solid #d8dee9;padding:.45rem;
text-align:left}
code,pre{font-family:ui-monospace,monospace;background:#f5f7fa;padding:.15rem .3rem;overflow:auto}
details{margin:1rem 0;padding:.6rem;border:1px solid #d8dee9;border-radius:.4rem}
input{padding:.5rem;width:min(30rem,90%)}.muted{color:#5b6472}.hidden{display:none}
""".strip()
_SCRIPT = """
const q=document.getElementById('filter');
q.addEventListener('input',()=>{const v=q.value.toLowerCase();
document.querySelectorAll('tbody tr').forEach(
r=>r.classList.toggle('hidden',!r.dataset.search.includes(v)));});
""".strip()

_DEFAULT_MAX_TOTAL_SEQUENCE_SYMBOLS = 10_000_000
_DEFAULT_MAX_TOTAL_RECORD_TEXT_CHARACTERS = 5_000_000
_DEFAULT_MAX_OUTPUT_BYTES = 50_000_000
_DEFAULT_MAX_RESULT_DEPTH = 64
_DEFAULT_MAX_RESULT_NODES = 1_000_000


def _positive_integer(name: str, value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ConfigurationError(f"{name} must be a positive integer.")
    return value


class _BoundedHTMLBuilder:
    """Accumulate HTML fragments without ever accepting an oversized payload."""

    def __init__(self, max_output_bytes: int) -> None:
        self._limit = max_output_bytes
        self._byte_count = 0
        self._parts: list[str] = []

    @property
    def remaining_bytes(self) -> int:
        return self._limit - self._byte_count

    def require_capacity(self, byte_count: int) -> None:
        if byte_count > self.remaining_bytes:
            raise ConfigurationError(
                "HTML report exceeds max_output_bytes.",
                code="HTML_REPORT_OUTPUT_LIMIT",
                context={
                    "max_output_bytes": self._limit,
                    "bytes_before_fragment": self._byte_count,
                    "minimum_fragment_bytes": byte_count,
                },
            )

    def append(self, fragment: str) -> None:
        encoded_size = len(fragment.encode("utf-8"))
        self.require_capacity(encoded_size)
        self._parts.append(fragment)
        self._byte_count += encoded_size

    def finish(self) -> str:
        return "".join(self._parts)


def _escaped_utf8_size(value: str, *, quote: bool = True) -> int:
    """Return the exact UTF-8 size produced by :func:`html.escape`."""

    size = len(value.encode("utf-8"))
    replacements = {"&": "&amp;", "<": "&lt;", ">": "&gt;"}
    if quote:
        replacements.update({'"': "&quot;", "'": "&#x27;"})
    for source, replacement in replacements.items():
        size += value.count(source) * (len(replacement) - 1)
    return size


def _display_sequence(record: DNARecord, *, max_bytes: int) -> str:
    pieces: list[str] = []
    byte_count = 0
    for part in record.sequence.parts:
        if isinstance(part, str):
            piece = part
        else:
            length = "?" if part.length is None else str(part.length)
            piece = f"<gap length={length} kind={part.kind.value}>"
        byte_count += len(piece)
        if byte_count > max_bytes:
            raise ConfigurationError(
                "HTML report exceeds max_output_bytes.",
                code="HTML_REPORT_OUTPUT_LIMIT",
                context={
                    "available_output_bytes": max_bytes,
                    "minimum_sequence_display_bytes": byte_count,
                },
            )
        pieces.append(piece)
    return "".join(pieces)


def _bounded_records(
    records: DNASet | Iterable[DNARecord],
    max_records: int,
    max_total_sequence_symbols: int,
    max_total_record_text_characters: int,
) -> tuple[DNARecord, ...]:
    if isinstance(records, DNASet):
        iterator = iter(records)
    else:
        try:
            iterator = iter(records)
        except TypeError as exc:
            raise ConfigurationError("records must be DNASet or an iterable of DNARecord.") from exc
    collected: list[DNARecord] = []
    sequence_symbol_count = 0
    record_text_characters = 0
    for index, record in enumerate(islice(iterator, max_records + 1)):
        if index >= max_records:
            raise ConfigurationError(
                "HTML report exceeds max_records.",
                code="HTML_REPORT_RECORD_LIMIT",
                context={"max_records": max_records},
            )
        if not isinstance(record, DNARecord):
            raise ConfigurationError(
                "records must contain only DNARecord objects.",
                context={"index": index, "type": type(record).__name__},
            )
        collected.append(record)
        sequence_symbol_count += record.sequence.symbol_length
        if sequence_symbol_count > max_total_sequence_symbols:
            raise ConfigurationError(
                "HTML report records exceed max_total_sequence_symbols.",
                code="HTML_REPORT_SEQUENCE_LIMIT",
                context={
                    "max_total_sequence_symbols": max_total_sequence_symbols,
                    "total_sequence_symbols": sequence_symbol_count,
                    "record_index": index,
                },
            )
        record_text_characters += len(record.id) + len(record.description)
        if record_text_characters > max_total_record_text_characters:
            raise ConfigurationError(
                "HTML report record text exceeds max_total_record_text_characters.",
                code="HTML_REPORT_RECORD_TEXT_LIMIT",
                context={
                    "max_total_record_text_characters": max_total_record_text_characters,
                    "total_record_text_characters": record_text_characters,
                    "record_index": index,
                },
            )
    return tuple(collected)


def _validate_result_shape(value: object, max_depth: int, max_nodes: int) -> None:
    nodes = [0]
    active: set[int] = set()

    def visit(item: object, depth: int) -> None:
        nodes[0] += 1
        if depth > max_depth or nodes[0] > max_nodes:
            raise ConfigurationError(
                "HTML report results exceed structural limits.",
                code="HTML_REPORT_RESULT_STRUCTURE_LIMIT",
                context={"max_result_depth": max_depth, "max_result_nodes": max_nodes},
            )
        if isinstance(item, (Mapping, list, tuple)):
            identity = id(item)
            if identity in active:
                raise ConfigurationError(
                    "HTML report results contain a recursive object.",
                    code="HTML_REPORT_RESULT_STRUCTURE_LIMIT",
                )
            active.add(identity)
            try:
                if isinstance(item, Mapping):
                    for key in item:
                        visit(key, depth + 1)
                        visit(item[key], depth + 1)
                else:
                    for child in item:
                        visit(child, depth + 1)
            finally:
                active.remove(identity)

    visit(value, 0)


def _json_string_lower_bound(value: str) -> int:
    size = 2
    for character in value:
        codepoint = ord(character)
        if character in {'"', "\\", "\b", "\f", "\n", "\r", "\t"}:
            size += 2
        elif codepoint < 0x20:
            size += 6
        else:
            size += len(character.encode("utf-8"))
    return size


def _require_json_lower_bound(value: object, max_bytes: int) -> None:
    """Reject oversized JSON before conversion or encoder chunk allocation."""

    total = [0]

    def add(amount: int) -> None:
        total[0] += amount
        if total[0] > max_bytes:
            raise ConfigurationError(
                "HTML report results exceed max_result_bytes.",
                code="HTML_REPORT_RESULT_LIMIT",
                context={"max_result_bytes": max_bytes},
            )

    def visit(item: object) -> None:
        if isinstance(item, str):
            add(_json_string_lower_bound(item))
        elif item is None or item is True:
            add(4)
        elif item is False:
            add(5)
        elif isinstance(item, (int, float)):
            add(len(str(item).encode("utf-8")))
        elif isinstance(item, Mapping):
            add(2)
            first = True
            for key in item:
                if not first:
                    add(1)
                first = False
                visit(str(key))
                add(1)
                visit(item[key])
        elif isinstance(item, (list, tuple)):
            add(2)
            first = True
            for child in item:
                if not first:
                    add(1)
                first = False
                visit(child)

    visit(value)


def _bounded_json(value: object, max_bytes: int, *, result_name: str) -> str:
    """Encode JSON incrementally while retaining at most ``max_bytes`` bytes."""

    _require_json_lower_bound(value, max_bytes)
    encoder = json.JSONEncoder(
        ensure_ascii=False,
        sort_keys=True,
        indent=2,
        allow_nan=False,
    )
    pieces: list[str] = []
    total = 0
    try:
        for piece in encoder.iterencode(to_json_compatible(value)):
            total += len(piece.encode("utf-8"))
            if total > max_bytes:
                raise ConfigurationError(
                    "HTML report results exceed max_result_bytes.",
                    code="HTML_REPORT_RESULT_LIMIT",
                    context={"max_result_bytes": max_bytes, "result_name": result_name},
                )
            pieces.append(piece)
    except ConfigurationError:
        raise
    except (TypeError, ValueError, RecursionError) as exc:
        raise ConfigurationError(
            "A report result is not JSON-compatible.", context={"result_name": result_name}
        ) from exc
    return "".join(pieces)


def _result_payload(
    results: Mapping[str, object],
    max_result_bytes: int,
    max_result_depth: int,
    max_result_nodes: int,
) -> tuple[tuple[str, str], ...]:
    if not isinstance(results, Mapping):
        raise ConfigurationError("results must be a mapping.")
    if (
        isinstance(max_result_bytes, bool)
        or not isinstance(max_result_bytes, int)
        or max_result_bytes <= 0
    ):
        raise ConfigurationError("max_result_bytes must be a positive integer.")
    payloads: list[tuple[str, str]] = []
    total = 0
    if any(not isinstance(name, str) or not name.strip() for name in results):
        raise ConfigurationError("result names must be non-empty strings.")
    for name in sorted(results):
        value = results[name]
        converter = getattr(value, "to_dict", None)
        if callable(converter):
            value = converter()
        _validate_result_shape(value, max_result_depth, max_result_nodes)
        remaining = max_result_bytes - total
        encoded = _bounded_json(value, remaining, result_name=name)
        total += len(encoded.encode("utf-8"))
        payloads.append((name, encoded))
    return tuple(payloads)


def build_html_report(
    records: DNASet | Iterable[DNARecord],
    *,
    results: Mapping[str, object] | None = None,
    title: str = "DNAKit report",
    max_records: int = 1_000,
    max_result_bytes: int = 5_000_000,
    max_total_sequence_symbols: int = _DEFAULT_MAX_TOTAL_SEQUENCE_SYMBOLS,
    max_total_record_text_characters: int = _DEFAULT_MAX_TOTAL_RECORD_TEXT_CHARACTERS,
    max_output_bytes: int = _DEFAULT_MAX_OUTPUT_BYTES,
    max_result_depth: int = _DEFAULT_MAX_RESULT_DEPTH,
    max_result_nodes: int = _DEFAULT_MAX_RESULT_NODES,
) -> HTMLReportArtifact:
    """Build a self-contained searchable and expandable local HTML report.

    ``max_total_sequence_symbols`` limits the cumulative nucleotide symbols
    across all records (explicit gap spans are excluded).  The cumulative
    number of Python characters in record IDs and descriptions is limited by
    ``max_total_record_text_characters``.  ``max_output_bytes`` is an exact
    upper bound on the final UTF-8 encoded HTML document.  Every limit must be
    a positive integer.
    """

    if not isinstance(title, str) or not title.strip():
        raise ConfigurationError("title must be a non-empty string.")
    resolved_max_records = _positive_integer("max_records", max_records)
    resolved_max_result_bytes = _positive_integer("max_result_bytes", max_result_bytes)
    resolved_max_total_sequence_symbols = _positive_integer(
        "max_total_sequence_symbols", max_total_sequence_symbols
    )
    resolved_max_total_record_text_characters = _positive_integer(
        "max_total_record_text_characters", max_total_record_text_characters
    )
    resolved_max_output_bytes = _positive_integer("max_output_bytes", max_output_bytes)
    resolved_max_result_depth = _positive_integer("max_result_depth", max_result_depth)
    resolved_max_result_nodes = _positive_integer("max_result_nodes", max_result_nodes)
    materialized = _bounded_records(
        records,
        resolved_max_records,
        resolved_max_total_sequence_symbols,
        resolved_max_total_record_text_characters,
    )
    resolved_results: Mapping[str, object] = {} if results is None else results
    builder = _BoundedHTMLBuilder(resolved_max_output_bytes)
    builder.append('<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">')
    builder.append('<meta name="viewport" content="width=device-width,initial-scale=1">')
    builder.append(
        '<meta http-equiv="Content-Security-Policy" '
        "content=\"default-src 'none'; style-src 'unsafe-inline'; script-src 'unsafe-inline'\">"
    )
    escaped_title_size = _escaped_utf8_size(title)
    builder.require_capacity(escaped_title_size * 2)
    escaped_title = html.escape(title)
    builder.append(f"<title>{escaped_title}</title><style>{_STYLE}</style></head><body>")
    builder.append(f'<h1>{escaped_title}</h1><p class="muted">本报告完全自包含, 不连接网络。</p>')
    builder.append(
        '<label for="filter">筛选记录:</label> <input id="filter" type="search" '
        'placeholder="ID、描述或序列">'
    )
    builder.append(
        "<table><thead><tr><th>ID</th><th>描述</th><th>符号长度</th><th>序列</th>"
        "</tr></thead><tbody>"
    )
    for record in materialized:
        sequence = _display_sequence(record, max_bytes=builder.remaining_bytes)
        minimum_search_bytes = len(record.id) + len(record.description) + len(sequence) + 2
        builder.require_capacity(minimum_search_bytes + len(sequence))
        search_source = f"{record.id} {record.description} {sequence}".lower()
        row_size = (
            len(
                '<tr data-search=""><td><code></code></td><td></td><td></td>'
                "<td><code></code></td></tr>"
            )
            + _escaped_utf8_size(search_source)
            + _escaped_utf8_size(record.id)
            + _escaped_utf8_size(record.description)
            + len(str(record.sequence.symbol_length))
            + _escaped_utf8_size(sequence)
        )
        builder.require_capacity(row_size)
        builder.append(
            f'<tr data-search="{html.escape(search_source)}"><td><code>'
            f"{html.escape(record.id)}</code></td><td>{html.escape(record.description)}</td>"
            f"<td>{record.sequence.symbol_length}</td><td><code>{html.escape(sequence)}</code>"
            "</td></tr>"
        )
    builder.append("</tbody></table><h2>分析结果</h2>")
    result_payloads = _result_payload(
        resolved_results,
        resolved_max_result_bytes,
        resolved_max_result_depth,
        resolved_max_result_nodes,
    )
    for name, encoded in result_payloads:
        block_size = (
            len("<details><summary></summary><pre></pre></details>")
            + _escaped_utf8_size(name)
            + _escaped_utf8_size(encoded)
        )
        builder.require_capacity(block_size)
        builder.append(
            f"<details><summary>{html.escape(name)}</summary>"
            f"<pre>{html.escape(encoded)}</pre></details>"
        )
    builder.append(f"<script>{_SCRIPT}</script></body></html>")
    document = builder.finish()
    return HTMLReportArtifact(document, len(materialized), tuple(sorted(resolved_results)))


__all__ = ["build_html_report"]
