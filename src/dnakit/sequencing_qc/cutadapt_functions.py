"""Offline cutadapt helpers: neighbour sets, k-mer search sets, adapter specs and trim indices.

Each function calls the cutadapt Python API directly (``cutadapt.align``,
``cutadapt.kmer_heuristic``, ``cutadapt.qualtrim``, ``cutadapt.parser``). Results are plain
dictionaries that keep the upstream values and add a ``provenance`` entry with the package name
and version. Inputs are explicit sequences, Phred strings and local FASTA files. No reads are
trimmed and no cutadapt command line is run.

``cutadapt`` is optional. It is imported only inside each function; when it is missing,
``BackendUnavailableError`` names the ``bio-cutadapt`` extra.
"""

from __future__ import annotations

import importlib.metadata
import math
import re
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require, sequence_text
from dnakit.exceptions import ConfigurationError

_CUTADAPT_EXTRA = "bio-cutadapt"
_DNA = frozenset("ACGT")
_ADAPTER_TYPES = ("front", "back", "anywhere")
_MAX_PHRED = 93
_MAX_SEQUENCE_LENGTH = 100_000
_MAX_SPEC_LENGTH = 4096
_MAX_NEIGHBOUR_DISTANCE = 3
_DEFAULT_MAX_NEIGHBOURS = 10_000
_MAX_NEIGHBOURS_LIMIT = 1_000_000
_DEFAULT_MAX_RECORDS = 10_000
_MAX_RECORDS_LIMIT = 1_000_000
_DEFAULT_MAX_ADAPTERS = 1_000
_MAX_ADAPTERS_LIMIT = 100_000
_DEFAULT_MAX_EXPANDED = 1_000_000
_MAX_EXPANDED_LIMIT = 10_000_000
_BRACE_COUNT = re.compile(r"\{(\d+)\}")
_HIT_CHANCE = re.compile(r"([0-9]+(?:\.[0-9]+)?)%")


def _version() -> str:
    try:
        return importlib.metadata.version("cutadapt")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _bounded_int(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _rate(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a number in [0, 1).")
    number = float(value)
    if not math.isfinite(number) or not 0 <= number < 1:
        raise ConfigurationError(f"{name} must be a finite number in [0, 1).")
    return number


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _choose(value: object, allowed: tuple[str, ...], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {', '.join(allowed)}; got {value!r}.")
    return value


def _acgt(value: object, name: str) -> str:
    """Upper-case DNA restricted to A, C, G and T, as used by cutadapt's neighbour search."""

    text = sequence_text(value)
    if len(text) > _MAX_SEQUENCE_LENGTH:
        raise ConfigurationError(f"{name} exceeds {_MAX_SEQUENCE_LENGTH} bases.")
    if set(text) - _DNA:
        raise ConfigurationError(f"{name} must contain only A, C, G and T.")
    return text


def _iupac(value: object, name: str) -> str:
    text = sequence_text(value)
    if len(text) > _MAX_SEQUENCE_LENGTH:
        raise ConfigurationError(f"{name} exceeds {_MAX_SEQUENCE_LENGTH} bases.")
    return text


def _spec_text(value: object, name: str) -> str:
    if (
        not isinstance(value, str)
        or not value
        or len(value) > _MAX_SPEC_LENGTH
        or "\0" in value
        or "\n" in value
    ):
        raise ConfigurationError(f"{name} must be a non-empty single-line string.")
    return value


def _path_text(value: object, name: str) -> str:
    text = str(value) if isinstance(value, (str, Path)) else ""
    if not text or "\0" in text:
        raise ConfigurationError(f"{name} must be a non-empty path without NUL characters.")
    return text


def _qualities(value: object, sequence: str, phred_offset: int) -> str:
    if not isinstance(value, str) or not value:
        raise ConfigurationError("qualities must be a non-empty Phred quality string.")
    if len(value) != len(sequence):
        raise ConfigurationError("qualities must contain one character per base.")
    if any(not phred_offset <= ord(char) <= phred_offset + _MAX_PHRED for char in value):
        raise ConfigurationError(
            f"qualities contain characters outside the Phred+{phred_offset} range."
        )
    return value


def _provenance(method: str, **extra: str) -> dict[str, str]:
    return {"package": "cutadapt", "version": _version(), "method": method, **extra}


def hamming_neighbours(
    sequence: str,
    max_distance: int,
    *,
    exact: bool = False,
    max_neighbours: int = _DEFAULT_MAX_NEIGHBOURS,
    verbose: bool = True,
) -> dict[str, Any]:
    """List DNA strings within a Hamming distance of a sequence, using cutadapt's enumeration.

    With ``exact=True`` the strings at Hamming distance exactly ``max_distance`` are returned
    (``cutadapt.align.hamming_sphere``). With ``exact=False`` every string at distance 0 to
    ``max_distance`` is returned (``cutadapt.align.hamming_environment``). Results are sorted by
    distance and then by sequence. The count is checked against ``max_neighbours`` while the
    upstream generator runs, so an oversized request stops early.
    """

    align = require("cutadapt.align", _CUTADAPT_EXTRA)
    text = _acgt(sequence, "sequence")
    distance = _bounded_int(
        max_distance, "max_distance", minimum=0, maximum=min(_MAX_NEIGHBOUR_DISTANCE, len(text))
    )
    selected = _flag(exact, "exact")
    limit = _bounded_int(max_neighbours, "max_neighbours", minimum=1, maximum=_MAX_NEIGHBOURS_LIMIT)
    rows: list[tuple[str, int, int]] = []
    if selected:
        for item in align.hamming_sphere(text, distance):
            rows.append((str(item), distance, len(text) - distance))
            if len(rows) > limit:
                raise ConfigurationError(f"neighbours exceed max_neighbours={limit}.")
        method = "cutadapt.align.hamming_sphere (exact Hamming distance)"
    else:
        for item, found, matches in align.hamming_environment(text, distance):
            rows.append((str(item), int(found), int(matches)))
            if len(rows) > limit:
                raise ConfigurationError(f"neighbours exceed max_neighbours={limit}.")
        method = "cutadapt.align.hamming_environment (Hamming distance up to max_distance)"
    rows.sort(key=lambda row: (row[1], row[0]))
    return {
        "status": "completed",
        "operation": "hamming_neighbours",
        "sequence": text,
        "max_distance": distance,
        "exact": selected,
        "count": len(rows),
        "neighbours": [
            {"sequence": item, "distance": found, "matches": matches}
            for item, found, matches in rows
        ],
        "provenance": _provenance(method, alphabet="ACGT"),
    }


def edit_neighbours(
    sequence: str,
    max_distance: int,
    *,
    max_neighbours: int = _DEFAULT_MAX_NEIGHBOURS,
    verbose: bool = True,
) -> dict[str, Any]:
    """List DNA strings within an edit (Levenshtein) distance of a sequence.

    Uses ``cutadapt.align.edit_environment``, which yields each string with its edit distance
    and the number of matches in the optimal alignment. Each string is reported once, at its
    smallest distance. Results are sorted by distance and then by sequence.
    """

    align = require("cutadapt.align", _CUTADAPT_EXTRA)
    text = _acgt(sequence, "sequence")
    distance = _bounded_int(
        max_distance, "max_distance", minimum=0, maximum=min(_MAX_NEIGHBOUR_DISTANCE, len(text))
    )
    limit = _bounded_int(max_neighbours, "max_neighbours", minimum=1, maximum=_MAX_NEIGHBOURS_LIMIT)
    best: dict[str, tuple[int, int]] = {}
    for item, found, matches in align.edit_environment(text, distance):
        key = str(item)
        current = best.get(key)
        if current is None or int(found) < current[0]:
            best[key] = (int(found), int(matches))
        if len(best) > limit:
            raise ConfigurationError(f"neighbours exceed max_neighbours={limit}.")
    rows = sorted((found, item, matches) for item, (found, matches) in best.items())
    return {
        "status": "completed",
        "operation": "edit_neighbours",
        "sequence": text,
        "max_distance": distance,
        "count": len(rows),
        "neighbours": [
            {"sequence": item, "distance": found, "matches": matches}
            for found, item, matches in rows
        ],
        "provenance": _provenance(
            "cutadapt.align.edit_environment (edit distance up to max_distance)", alphabet="ACGT"
        ),
    }


def _kmer_inputs(
    adapter: object,
    min_overlap: object,
    max_error_rate: object,
    back_adapter: object,
    front_adapter: object,
    internal: object,
) -> tuple[str, int, float, bool, bool, bool]:
    text = _acgt(adapter, "adapter")
    overlap = _bounded_int(min_overlap, "min_overlap", minimum=1, maximum=len(text))
    rate = _rate(max_error_rate, "max_error_rate")
    return (
        text,
        overlap,
        rate,
        _flag(back_adapter, "back_adapter"),
        _flag(front_adapter, "front_adapter"),
        _flag(internal, "internal"),
    )


def adapter_kmer_positions(
    adapter: str,
    *,
    min_overlap: int = 3,
    max_error_rate: float = 0.1,
    back_adapter: bool = True,
    front_adapter: bool = False,
    internal: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build cutadapt's k-mer search sets for one adapter, one entry per search offset.

    Calls ``cutadapt.kmer_heuristic.create_positions_and_kmers``. Each entry gives a ``start``
    offset, a ``stop`` offset (``None`` means the read end) and the k-mers that must occur in the
    read for the adapter to match at that offset. Offsets follow cutadapt's convention, so a
    negative ``start`` counts back from the read end. K-mers are sorted inside each entry, and
    entries are sorted by ``(start, stop)``.
    """

    heuristic = require("cutadapt.kmer_heuristic", _CUTADAPT_EXTRA)
    text, overlap, rate, back, front, inner = _kmer_inputs(
        adapter, min_overlap, max_error_rate, back_adapter, front_adapter, internal
    )
    raw = heuristic.create_positions_and_kmers(text, overlap, rate, back, front, inner)
    entries = []
    for start, stop, kmers in raw:
        entries.append(
            {
                "start": int(start),
                "stop": None if stop is None else int(stop),
                "kmers": sorted(str(item) for item in kmers),
            }
        )
    entries.sort(
        key=lambda item: (item["start"], math.inf if item["stop"] is None else item["stop"])
    )
    return {
        "status": "completed",
        "operation": "adapter_kmer_positions",
        "adapter": text,
        "min_overlap": overlap,
        "max_error_rate": rate,
        "back_adapter": back,
        "front_adapter": front,
        "internal": inner,
        "entry_count": len(entries),
        "kmer_count": sum(len(item["kmers"]) for item in entries),
        "entries": entries,
        "provenance": _provenance("cutadapt.kmer_heuristic.create_positions_and_kmers"),
    }


def adapter_kmer_table(
    adapter: str,
    *,
    min_overlap: int = 3,
    max_error_rate: float = 0.1,
    back_adapter: bool = True,
    front_adapter: bool = False,
    internal: bool = True,
    default_length: int = 150,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate how often the adapter's k-mer search sets hit random sequence.

    Runs ``create_positions_and_kmers`` and then
    ``cutadapt.kmer_heuristic.kmer_probability_analysis`` with ``default_length`` as the read
    length. The upstream text table is parsed into rows with the k-mer, its offsets, the number
    of sites considered and the percentage of random sequences that would hit the k-mer.
    ``profile_hit_chance_percent`` is cutadapt's overall estimate for the whole search profile.
    The raw upstream text is returned as ``table_text``.
    """

    heuristic = require("cutadapt.kmer_heuristic", _CUTADAPT_EXTRA)
    text, overlap, rate, back, front, inner = _kmer_inputs(
        adapter, min_overlap, max_error_rate, back_adapter, front_adapter, internal
    )
    length = _bounded_int(default_length, "default_length", minimum=1, maximum=1_000_000)
    raw = heuristic.create_positions_and_kmers(text, overlap, rate, back, front, inner)
    table = str(heuristic.kmer_probability_analysis(raw, length))
    rows: list[dict[str, Any]] = []
    overall: float | None = None
    header_seen = False
    for line in table.splitlines():
        if not line.strip():
            continue
        if line.startswith("Chance for profile hit"):
            match = _HIT_CHANCE.search(line)
            if match is None:
                raise ConfigurationError("cutadapt returned an unexpected profile summary line.")
            overall = float(match.group(1))
            continue
        cells = [cell.strip() for cell in line.split("\t")]
        if not header_seen:
            header_seen = True
            if cells[:2] != ["kmer", "start"]:
                raise ConfigurationError("cutadapt returned an unexpected k-mer table header.")
            continue
        if len(cells) != 5:
            raise ConfigurationError("cutadapt returned a k-mer table row with unexpected columns.")
        rows.append(
            {
                "kmer": cells[0],
                "start": int(cells[1]),
                "stop": None if cells[2] == "None" else int(cells[2]),
                "considered_sites": int(cells[3]),
                "hit_chance_percent": float(cells[4]),
            }
        )
    return {
        "status": "completed",
        "operation": "adapter_kmer_table",
        "adapter": text,
        "default_length": length,
        "row_count": len(rows),
        "rows": rows,
        "profile_hit_chance_percent": overall,
        "table_text": table,
        "provenance": _provenance(
            "cutadapt.kmer_heuristic.create_positions_and_kmers + kmer_probability_analysis"
        ),
    }


def sequence_kmer_chunks(sequence: str, chunks: int, *, verbose: bool = True) -> dict[str, Any]:
    """Split a sequence into almost equal chunks and return the distinct chunk strings.

    Calls ``cutadapt.kmer_heuristic.kmer_chunks``. The first ``len % chunks`` chunks are one base
    longer than the rest. Chunks are returned sorted. For ``ACGTACGTAC`` with 3 chunks the chunks
    are ``ACGT``, ``ACG`` and ``TAC``.
    """

    heuristic = require("cutadapt.kmer_heuristic", _CUTADAPT_EXTRA)
    text = _iupac(sequence, "sequence")
    count = _bounded_int(chunks, "chunks", minimum=1, maximum=len(text))
    items = sorted(str(item) for item in heuristic.kmer_chunks(text, count))
    return {
        "status": "completed",
        "operation": "sequence_kmer_chunks",
        "sequence": text,
        "chunks": count,
        "count": len(items),
        "kmers": items,
        "provenance": _provenance("cutadapt.kmer_heuristic.kmer_chunks"),
    }


def poly_a_tail_start(
    sequence: str, *, reverse_complement: bool = False, verbose: bool = True
) -> dict[str, Any]:
    """Find the poly-A tail at the 3' end of a sequence, or the poly-T head if reverse-complemented.

    Calls ``cutadapt.qualtrim.poly_a_trim_index``. ``index`` is the 0-based start of the tail
    (or the end of the head when ``reverse_complement`` is true), and ``length`` is the number of
    bases in the tail or head. Tails shorter than three bases are ignored by cutadapt, which then
    returns the sequence length (no tail) or zero (no head).
    """

    qualtrim = require("cutadapt.qualtrim", _CUTADAPT_EXTRA)
    text = _iupac(sequence, "sequence")
    flag = _flag(reverse_complement, "reverse_complement")
    index = int(qualtrim.poly_a_trim_index(text, revcomp=flag))
    length = index if flag else len(text) - index
    return {
        "status": "completed",
        "operation": "poly_a_tail_start",
        "sequence": text,
        "sequence_length": len(text),
        "reverse_complement": flag,
        "index": index,
        "length": length,
        "provenance": _provenance(
            "cutadapt.qualtrim.poly_a_trim_index", minimum_tail_length="3 (upstream rule)"
        ),
    }


def nextseq_trim_bounds(
    sequence: str,
    qualities: str,
    *,
    cutoff: int = 20,
    phred_offset: int = 33,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the 3' trim position for NextSeq two-colour quality trimming of one read.

    Calls ``cutadapt.qualtrim.nextseq_trim_index`` on a ``dnaio`` record. The algorithm is the BWA
    running-sum rule used by ``find_quality_trim_bounds``, except that G bases are scored as if
    their quality were ``cutoff - 1``. This removes runs of G at the 3' end with no colour
    signal. ``index`` is the retained length, because only the 3' end is trimmed.
    """

    qualtrim = require("cutadapt.qualtrim", _CUTADAPT_EXTRA)
    dnaio = require("dnaio", _CUTADAPT_EXTRA)
    text = _iupac(sequence, "sequence")
    offset = _choose_offset(phred_offset)
    quals = _qualities(qualities, text, offset)
    cut = _bounded_int(cutoff, "cutoff", minimum=0, maximum=_MAX_PHRED)
    record = dnaio.SequenceRecord("read", text, quals)
    index = int(qualtrim.nextseq_trim_index(record, cut, offset))
    return {
        "status": "completed",
        "operation": "nextseq_trim_bounds",
        "length": len(text),
        "cutoff": cut,
        "phred_offset": offset,
        "index": index,
        "retained_length": index,
        "provenance": _provenance(
            "cutadapt.qualtrim.nextseq_trim_index (BWA running sum, NextSeq G rule)"
        ),
    }


def _choose_offset(value: object) -> int:
    if isinstance(value, bool) or value not in (33, 64):
        raise ConfigurationError("phred_offset must be 33 or 64.")
    return int(value)


def expand_adapter_braces(
    sequence: str,
    *,
    max_length: int = _DEFAULT_MAX_EXPANDED,
    verbose: bool = True,
) -> dict[str, Any]:
    """Expand cutadapt's ``x{n}`` repeat notation in an adapter sequence.

    Calls ``cutadapt.parser.expand_braces``, which replaces each ``x{n}`` by ``n`` copies of
    ``x``. For example ``TGA{5}CT`` becomes ``TGAAAAACT``. Before expansion the wrapper adds up the
    repeat counts, so an expansion larger than ``max_length`` is refused without building it.
    Malformed notation raises ``ConfigurationError``.
    """

    parser = require("cutadapt.parser", _CUTADAPT_EXTRA)
    text = _spec_text(sequence, "sequence")
    limit = _bounded_int(max_length, "max_length", minimum=1, maximum=_MAX_EXPANDED_LIMIT)
    estimate = len(text) + sum(int(count) for count in _BRACE_COUNT.findall(text))
    if estimate > limit:
        raise ConfigurationError(f"expanded sequence would exceed max_length={limit}.")
    try:
        expanded = str(parser.expand_braces(text))
    except ValueError as exc:
        raise ConfigurationError(f"invalid cutadapt brace notation: {exc}") from exc
    return {
        "status": "completed",
        "operation": "expand_adapter_braces",
        "input": text,
        "input_length": len(text),
        "expanded": expanded,
        "expanded_length": len(expanded),
        "provenance": _provenance("cutadapt.parser.expand_braces"),
    }


def read_adapter_fasta(
    source: str | Path,
    *,
    max_records: int = _DEFAULT_MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read adapter names and sequences from a local FASTA file, in file order.

    Calls ``cutadapt.parser.read_adapters_fasta``. The name is the first word of each header
    (``None`` for an empty header). Sequences are returned as written in the file. A file that
    is not FASTA, or that has more than ``max_records`` records, raises ``ConfigurationError``.
    """

    parser = require("cutadapt.parser", _CUTADAPT_EXTRA)
    dnaio = require("dnaio", _CUTADAPT_EXTRA)
    path = input_path(_path_text(source, "source"))
    limit = _bounded_int(max_records, "max_records", minimum=1, maximum=_MAX_RECORDS_LIMIT)
    records: list[dict[str, Any]] = []
    try:
        for name, sequence in parser.read_adapters_fasta(str(path)):
            if len(records) >= limit:
                raise ConfigurationError(f"{path.name} has more than max_records={limit} records.")
            records.append({"name": name, "sequence": str(sequence)})
    except ConfigurationError:
        raise
    except (OSError, ValueError, dnaio.exceptions.FileFormatError) as exc:
        raise ConfigurationError(
            f"cutadapt could not read adapters from {path.name}: {exc}"
        ) from exc
    return {
        "status": "completed",
        "operation": "read_adapter_fasta",
        "source": str(path),
        "count": len(records),
        "records": records,
        "provenance": _provenance("cutadapt.parser.read_adapters_fasta"),
    }


def _adapter_row(adapter: Any) -> dict[str, Any]:
    row: dict[str, Any] = {
        "name": None if adapter.name is None else str(adapter.name),
        "adapter_class": type(adapter).__name__,
        "description": str(getattr(adapter, "description", "")),
    }
    if hasattr(adapter, "front_adapter") and hasattr(adapter, "back_adapter"):
        front = adapter.front_adapter
        back = adapter.back_adapter
        row.update(
            {
                "linked": True,
                "sequence": None,
                "front_sequence": None if front is None else str(front.sequence),
                "back_sequence": None if back is None else str(back.sequence),
            }
        )
        return row
    row.update(
        {
            "linked": False,
            "sequence": str(adapter.sequence),
            "max_error_rate": float(adapter.max_error_rate),
            "min_overlap": int(adapter.min_overlap),
            "indels": bool(adapter.indels),
            "read_wildcards": bool(adapter.read_wildcards),
            "adapter_wildcards": bool(adapter.adapter_wildcards),
        }
    )
    return row


def parse_adapter_spec(
    spec: str,
    *,
    adapter_type: str = "back",
    max_errors: float = 0.1,
    min_overlap: int = 3,
    indels: bool = True,
    max_adapters: int = _DEFAULT_MAX_ADAPTERS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse one cutadapt adapter specification into adapter objects and report their settings.

    Calls ``cutadapt.parser.make_adapters_from_one_specification``. ``adapter_type`` is
    ``"front"`` (-g), ``"back"`` (-a) or ``"anywhere"`` (-b). The specification may be a sequence
    (``ACGT``), a named sequence (``name=ACGT``), a linked pair (``ACGT...TTTT``, not allowed for
    ``anywhere``) or ``file:<path>`` to read every adapter from a local FASTA file. Linked
    adapters report their front and back sequences. Returns ``adapters`` in upstream order.
    """

    parser = require("cutadapt.parser", _CUTADAPT_EXTRA)
    text = _spec_text(spec, "spec")
    selected = _choose(adapter_type, _ADAPTER_TYPES, "adapter_type")
    errors = _rate(max_errors, "max_errors")
    overlap = _bounded_int(min_overlap, "min_overlap", minimum=1, maximum=1000)
    use_indels = _flag(indels, "indels")
    limit = _bounded_int(max_adapters, "max_adapters", minimum=1, maximum=_MAX_ADAPTERS_LIMIT)
    parameters = {"max_errors": errors, "min_overlap": overlap, "indels": use_indels}
    adapters: list[dict[str, Any]] = []
    try:
        for adapter in parser.make_adapters_from_one_specification(text, selected, parameters):
            if len(adapters) >= limit:
                raise ConfigurationError(f"spec yields more than max_adapters={limit} adapters.")
            adapters.append(_adapter_row(adapter))
    except ConfigurationError:
        raise
    except (OSError, ValueError, TypeError, RuntimeError) as exc:
        raise ConfigurationError(f"cutadapt rejected the adapter specification: {exc}") from exc
    return {
        "status": "completed",
        "operation": "parse_adapter_spec",
        "spec": text,
        "adapter_type": selected,
        "count": len(adapters),
        "adapters": adapters,
        "provenance": _provenance(
            "cutadapt.parser.make_adapters_from_one_specification",
            file_notation="file:<path> reads a local FASTA file",
        ),
    }


__all__ = [
    "adapter_kmer_positions",
    "adapter_kmer_table",
    "edit_neighbours",
    "expand_adapter_braces",
    "hamming_neighbours",
    "nextseq_trim_bounds",
    "parse_adapter_spec",
    "poly_a_tail_start",
    "read_adapter_fasta",
    "sequence_kmer_chunks",
]

_register_public_verbose(globals(), tuple(__all__))
