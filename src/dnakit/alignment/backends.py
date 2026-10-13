"""Optional-backend adapters for mappy, edlib, parasail and pyhmmer.

The functions in this module are explicit adapters over optional packages:

- ``map_query_sequences_paf()``, ``build_minimap2_index()`` and
  ``fetch_minimap2_reference_region()`` use ``mappy`` (minimap2). Mapping records
  use PAF field names with zero-based half-open spans.
- ``align_dna_edlib()``, ``locate_approximate_pattern()`` and
  ``calculate_edit_distance_matrix()`` use ``edlib``. The ``iupac`` option makes
  ambiguity codes equal to the bases they cover.
- ``align_dna_parasail()``, ``score_dna_parasail()`` and
  ``score_dna_parasail_matrix()`` use ``parasail`` for global, local and
  semi-global alignment with affine gap penalties.
- ``build_dna_hmm()``, ``search_dna_hmm()`` and ``read_hmm_profile()`` use
  ``pyhmmer`` to build DNA profile HMMs from an aligned file, search DNA sequence
  files with nhmmer, and read HMMER3 profile metadata.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[...]"`` hint.
"""

from __future__ import annotations

import importlib.metadata as metadata
import re
from collections.abc import Collection, Mapping
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import (
    input_path,
    named_sequence_items,
    require,
    sequence_text,
)
from dnakit.core import DNA
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_MINIMAP2_EXTRA = "bio-mappy"
_EDLIB_EXTRA = "bio-edlib"
_PARASAIL_EXTRA = "bio-parasail"
_PYHMMER_EXTRA = "bio-pyhmmer"
_MINIMAP2_PRESETS = ("asm10", "asm20", "asm5", "map-hifi", "map-ont", "map-pb", "sr")
_PARASAIL_MODES = {"global": "nw", "local": "sw", "semiglobal": "sg"}
_PARASAIL_MATRICES = ("custom", "dnafull")
_PARASAIL_ALPHABET = "ACGTN"
_DNA_LETTERS = frozenset("ACGTN")
_IUPAC_LETTERS = frozenset("ACGTRYSWKMBDHVN")
_IUPAC_BASES = {
    "R": "AG",
    "Y": "CT",
    "S": "CG",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}
_IUPAC_EQUALITIES = [(code, base) for code, bases in _IUPAC_BASES.items() for base in bases]
_CIGAR_OPERATION = re.compile(r"(\d+)([=XID])")
_NAME_PATTERN = re.compile(r"[^\s\x00]+")
_HMM_NAME_PATTERN = re.compile(r"[A-Za-z0-9_.-]{1,128}")
_DEFAULT_MAX_REFERENCE_BYTES = 2_000_000_000
_MAX_REFERENCE_BYTES = 100_000_000_000
_SEQUENCE_LIMIT = 1_000_000
_DEFAULT_MAX_CELLS = 25_000_000
_MAX_CELLS = 1_000_000_000
_MAX_COORDINATE = 10**12


def _package(name: str) -> str:
    try:
        return f"{name} {metadata.version(name)}"
    except metadata.PackageNotFoundError:
        return name


def _choose(value: object, allowed: Collection[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        options = ", ".join(sorted(allowed))
        raise ConfigurationError(f"{name} must be one of {options}; got {value!r}.")
    return value


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _real(value: object, name: str, *, minimum: float, maximum: float) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not minimum <= float(value) <= maximum
    ):
        raise ConfigurationError(f"{name} must be a number in [{minimum:g}, {maximum:g}].")
    return float(value)


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _finite(value: object) -> float | None:
    number = float(value)
    return number if number == number and abs(number) != float("inf") else None


def _name(value: object, label: str) -> str:
    if not isinstance(value, str) or not _NAME_PATTERN.fullmatch(value):
        raise ConfigurationError(
            f"{label} must be a non-empty string without whitespace or NUL characters."
        )
    return value


def _path_argument(value: object, label: str) -> str:
    if not isinstance(value, (str, Path)) or not str(value) or "\x00" in str(value):
        raise ConfigurationError(f"{label} must be a non-empty path without NUL characters.")
    return str(value)


def _letters(value: object, label: str, *, allowed: frozenset[str], maximum: int) -> str:
    text = sequence_text(value)
    unsupported = sorted(set(text) - allowed)
    if unsupported:
        raise ConfigurationError(f"{label} contains unsupported symbols: {''.join(unsupported)}.")
    if len(text) > maximum:
        raise ConfigurationError(f"{label} is longer than max_length={maximum}.")
    return text


def _named_sequences(
    value: object,
    *,
    allowed: frozenset[str],
    max_count: int,
    limit_name: str,
    max_length: int,
) -> list[tuple[str, str]]:
    items = list(named_sequence_items(value))
    if not items:
        raise ConfigurationError("at least one named sequence is required.")
    if len(items) > max_count:
        raise ConfigurationError(f"sequence count exceeds {limit_name}={max_count}.")
    named: list[tuple[str, str]] = []
    for raw_name, raw_sequence in items:
        name = _name(raw_name, "sequence name")
        text = _letters(raw_sequence, f"sequence {name!r}", allowed=allowed, maximum=max_length)
        named.append((name, text))
    return named


def _reference_file(value: object, *, limit: int) -> Path:
    path = input_path(_path_argument(value, "reference"))
    size = path.stat().st_size
    if size > limit:
        raise ConfigurationError(
            f"reference {path.name} is {size} bytes; raise max_reference_bytes to read it."
        )
    return path


def _output_file(value: object, label: str, *, overwrite: bool) -> Path:
    path = Path(_path_argument(value, label)).expanduser().resolve()
    if not path.parent.is_dir():
        raise ConfigurationError(f"{label} directory does not exist: {path.parent}")
    if path.is_dir():
        raise ConfigurationError(f"{label} is a directory: {path}")
    if path.exists() and not overwrite:
        raise ConfigurationError(f"{label} already exists: {path}; pass overwrite=True.")
    return path


def _reference_limit(value: object) -> int:
    return _integer(value, "max_reference_bytes", minimum=1, maximum=_MAX_REFERENCE_BYTES)


# mappy (minimap2)


def _minimap2_aligner(
    mappy: Any, reference: Path, *, preset: str | None, output: Path | None = None
):
    aligner_class = getattr(mappy, "Aligner", None)
    if not callable(aligner_class):
        raise ConfigurationError("The installed mappy package lacks Aligner().")
    options: dict[str, Any] = {"n_threads": 1}
    if preset is not None:
        options["preset"] = preset
    if output is not None:
        options["fn_idx_out"] = str(output)
    try:
        aligner = aligner_class(str(reference), **options)
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"mappy could not load {reference.name}: {exc}") from exc
    if not aligner or int(aligner.n_seq) < 1:
        raise ConfigurationError(f"mappy found no reference sequences in {reference.name}.")
    return aligner


def _paf_record(name: str, query_length: int, hit: Any) -> dict[str, Any]:
    return {
        "qname": name,
        "qlen": query_length,
        "qstart": int(hit.q_st),
        "qend": int(hit.q_en),
        "strand": "+" if int(hit.strand) > 0 else "-",
        "tname": str(hit.ctg),
        "tlen": int(hit.ctg_len),
        "tstart": int(hit.r_st),
        "tend": int(hit.r_en),
        "nmatch": int(hit.mlen),
        "alen": int(hit.blen),
        "mapq": int(hit.mapq),
        "primary": bool(hit.is_primary),
        "edit_distance": int(hit.NM),
        "cigar": str(hit.cigar_str),
    }


def map_query_sequences_paf(
    reference: str | Path,
    queries: DNA | Mapping[str, str],
    *,
    preset: str = "map-ont",
    max_queries: int = 10_000,
    max_hits_per_query: int = 100,
    max_reference_bytes: int = _DEFAULT_MAX_REFERENCE_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map named DNA queries to a FASTA or .mmi reference with minimap2 as PAF-like records."""

    preset_name = _choose(preset, _MINIMAP2_PRESETS, "preset")
    query_limit = _integer(max_queries, "max_queries", minimum=1, maximum=100_000)
    hit_limit = _integer(max_hits_per_query, "max_hits_per_query", minimum=1, maximum=1_000)
    reference_path = _reference_file(reference, limit=_reference_limit(max_reference_bytes))
    named = _named_sequences(
        queries,
        allowed=_DNA_LETTERS,
        max_count=query_limit,
        limit_name="max_queries",
        max_length=_SEQUENCE_LIMIT,
    )
    mappy = require("mappy", _MINIMAP2_EXTRA)
    aligner = _minimap2_aligner(mappy, reference_path, preset=preset_name)
    records: list[dict[str, Any]] = []
    for name, sequence in named:
        hits = list(aligner.map(sequence))
        if len(hits) > hit_limit:
            raise ConfigurationError(
                f"mappy returned more than max_hits_per_query={hit_limit} for {name!r}."
            )
        records.extend(_paf_record(name, len(sequence), hit) for hit in hits)
    return {
        "status": "completed",
        "operation": "map_query_sequences_paf",
        "reference": str(reference_path),
        "preset": preset_name,
        "query_count": len(named),
        "hit_count": len(records),
        "records": records,
        "provenance": {
            "package": _package("mappy"),
            "coordinate_system": "PAF: zero-based half-open query and target spans",
            "strand": "'-' means the query aligned as its reverse complement",
            "tags": "edit_distance is the NM tag; cigar is the minimap2 CIGAR string",
            "backend_scope": "in-memory queries against one reference; no paired-end or splice",
        },
    }


def build_minimap2_index(
    reference: str | Path,
    output_path: str | Path,
    *,
    preset: str = "map-ont",
    overwrite: bool = False,
    max_reference_bytes: int = _DEFAULT_MAX_REFERENCE_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a minimap2 .mmi index for a FASTA reference and list its contigs."""

    preset_name = _choose(preset, _MINIMAP2_PRESETS, "preset")
    replace = _flag(overwrite, "overwrite")
    reference_path = _reference_file(reference, limit=_reference_limit(max_reference_bytes))
    target = _output_file(output_path, "output_path", overwrite=replace)
    if target == reference_path:
        raise ConfigurationError("output_path must differ from the reference path.")
    mappy = require("mappy", _MINIMAP2_EXTRA)
    aligner = _minimap2_aligner(mappy, reference_path, preset=preset_name, output=target)
    if not target.is_file():
        raise BackendExecutionError(f"mappy did not write the index {target}.")
    return {
        "status": "completed",
        "operation": "build_minimap2_index",
        "reference": str(reference_path),
        "output": str(target),
        "preset": preset_name,
        "contig_count": int(aligner.n_seq),
        "contigs": [str(name) for name in aligner.seq_names],
        "index_bytes": target.stat().st_size,
        "provenance": {
            "package": _package("mappy"),
            "format": "minimap2 .mmi index written by Aligner(fn_idx_out=...)",
            "coordinate_system": "contig names as listed in the FASTA headers",
        },
    }


def fetch_minimap2_reference_region(
    reference: str | Path,
    contig: str,
    start: int,
    end: int,
    *,
    max_reference_bytes: int = _DEFAULT_MAX_REFERENCE_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return one zero-based half-open reference interval through minimap2's sequence store."""

    name = _name(contig, "contig")
    start_value = _integer(start, "start", minimum=0, maximum=_MAX_COORDINATE)
    end_value = _integer(end, "end", minimum=1, maximum=_MAX_COORDINATE)
    if end_value <= start_value:
        raise ConfigurationError("end must be greater than start.")
    reference_path = _reference_file(reference, limit=_reference_limit(max_reference_bytes))
    mappy = require("mappy", _MINIMAP2_EXTRA)
    aligner = _minimap2_aligner(mappy, reference_path, preset=None)
    if name not in aligner.seq_names:
        raise ConfigurationError(f"contig {name!r} is not in the reference.")
    sequence = aligner.seq(name, start_value, end_value)
    if not isinstance(sequence, str) or len(sequence) != end_value - start_value:
        raise ConfigurationError(
            f"region {name}:{start_value}-{end_value} extends past the contig."
        )
    return {
        "status": "completed",
        "operation": "fetch_minimap2_reference_region",
        "reference": str(reference_path),
        "contig": name,
        "start": start_value,
        "end": end_value,
        "length": len(sequence),
        "sequence": sequence,
        "provenance": {
            "package": _package("mappy"),
            "coordinate_system": "zero-based half-open contig coordinates",
            "backend_scope": "loads the whole reference before slicing; a .mmi index loads faster",
        },
    }


# edlib


def _edlib_error(exc: Exception, label: str) -> ConfigurationError:
    return ConfigurationError(f"Edlib could not {label}: {exc}")


def _cigar_counts(cigar: str) -> dict[str, int]:
    counts = {"=": 0, "X": 0, "I": 0, "D": 0}
    position = 0
    for match in _CIGAR_OPERATION.finditer(cigar):
        if match.start() != position:
            raise BackendExecutionError(f"Edlib returned an unsupported CIGAR: {cigar!r}")
        counts[match.group(2)] += int(match.group(1))
        position = match.end()
    if not cigar or position != len(cigar):
        raise BackendExecutionError(f"Edlib returned an unsupported CIGAR: {cigar!r}")
    return counts


def align_dna_edlib(
    query: str | DNA,
    target: str | DNA,
    *,
    mode: str = "NW",
    k: int = -1,
    iupac: bool = False,
    max_length: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align two DNA sequences with Edlib and return an extended CIGAR and edit counts."""

    mode_name = _choose(mode, ("HW", "NW"), "mode")
    max_distance = _integer(k, "k", minimum=-1, maximum=2_147_483_647)
    ambiguous = _flag(iupac, "iupac")
    length_limit = _integer(max_length, "max_length", minimum=1, maximum=_SEQUENCE_LIMIT)
    query_text = _letters(query, "query", allowed=_IUPAC_LETTERS, maximum=length_limit)
    target_text = _letters(target, "target", allowed=_IUPAC_LETTERS, maximum=length_limit)
    edlib = require("edlib", _EDLIB_EXTRA)
    try:
        result = edlib.align(
            query_text,
            target_text,
            mode=mode_name,
            task="path",
            k=max_distance,
            additionalEqualities=_IUPAC_EQUALITIES if ambiguous else None,
        )
    except (AttributeError, TypeError, ValueError) as exc:
        raise _edlib_error(exc, "align the sequences") from exc
    distance = int(result["editDistance"])
    base = {
        "status": "completed",
        "operation": "align_dna_edlib",
        "mode": mode_name,
        "k": max_distance,
        "iupac": ambiguous,
        "query_length": len(query_text),
        "target_length": len(target_text),
        "provenance": {
            "package": _package("edlib"),
            "cigar": "extended: '=' match, 'X' mismatch, 'I' query base absent from target, "
            "'D' target base absent from query",
            "coordinate_system": "zero-based half-open spans",
            "ambiguity": "IUPAC codes equal their bases when iupac=True; literal otherwise",
        },
    }
    if distance < 0:
        return {
            **base,
            "found": False,
            "edit_distance": None,
            "query_start": None,
            "query_end": None,
            "target_start": None,
            "target_end": None,
            "cigar": None,
            "matches": None,
            "mismatches": None,
            "insertions": None,
            "deletions": None,
            "identity": None,
        }
    cigar = str(result["cigar"])
    counts = _cigar_counts(cigar)
    query_consumed = counts["="] + counts["X"] + counts["I"]
    target_consumed = counts["="] + counts["X"] + counts["D"]
    locations = result.get("locations") or []
    if query_consumed != len(query_text) or not locations:
        raise BackendExecutionError("Edlib returned an alignment that does not span the query.")
    start, end = (int(value) for value in locations[0])
    if end + 1 - start != target_consumed:
        raise BackendExecutionError("Edlib returned an alignment that does not match its span.")
    columns = sum(counts.values())
    return {
        **base,
        "found": True,
        "edit_distance": distance,
        "query_start": 0,
        "query_end": len(query_text),
        "target_start": start,
        "target_end": end + 1,
        "cigar": cigar,
        "matches": counts["="],
        "mismatches": counts["X"],
        "insertions": counts["I"],
        "deletions": counts["D"],
        "identity": counts["="] / columns if columns else None,
    }


def locate_approximate_pattern(
    pattern: str | DNA,
    text: str | DNA,
    *,
    max_edits: int = 0,
    iupac: bool = False,
    max_length: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Locate the lowest-distance approximate occurrences of a DNA pattern in a text with Edlib."""

    edits = _integer(max_edits, "max_edits", minimum=0, maximum=1_000)
    ambiguous = _flag(iupac, "iupac")
    length_limit = _integer(max_length, "max_length", minimum=1, maximum=_SEQUENCE_LIMIT)
    pattern_text = _letters(pattern, "pattern", allowed=_IUPAC_LETTERS, maximum=length_limit)
    text_text = _letters(text, "text", allowed=_IUPAC_LETTERS, maximum=length_limit)
    edlib = require("edlib", _EDLIB_EXTRA)
    try:
        result = edlib.align(
            pattern_text,
            text_text,
            mode="HW",
            task="locations",
            k=edits,
            additionalEqualities=_IUPAC_EQUALITIES if ambiguous else None,
        )
    except (AttributeError, TypeError, ValueError) as exc:
        raise _edlib_error(exc, "locate the pattern") from exc
    distance = int(result["editDistance"])
    locations = (
        [{"start": int(start), "end": int(end) + 1} for start, end in result["locations"]]
        if distance >= 0
        else []
    )
    return {
        "status": "completed",
        "operation": "locate_approximate_pattern",
        "found": distance >= 0,
        "edit_distance": distance if distance >= 0 else None,
        "max_edits": edits,
        "iupac": ambiguous,
        "pattern_length": len(pattern_text),
        "text_length": len(text_text),
        "location_count": len(locations),
        "locations": locations,
        "provenance": {
            "package": _package("edlib"),
            "reporting": "spans where the lowest edit distance is reached; other occurrences "
            "within max_edits are not listed",
            "coordinate_system": "zero-based half-open text spans",
        },
    }


def calculate_edit_distance_matrix(
    sequences: DNA | Mapping[str, str],
    *,
    iupac: bool = False,
    max_sequences: int = 100,
    max_length: int = 5_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute pairwise global edit distances between named DNA sequences with Edlib."""

    count_limit = _integer(max_sequences, "max_sequences", minimum=1, maximum=1_000)
    length_limit = _integer(max_length, "max_length", minimum=1, maximum=_SEQUENCE_LIMIT)
    ambiguous = _flag(iupac, "iupac")
    named = _named_sequences(
        sequences,
        allowed=_IUPAC_LETTERS,
        max_count=count_limit,
        limit_name="max_sequences",
        max_length=length_limit,
    )
    edlib = require("edlib", _EDLIB_EXTRA)
    size = len(named)
    matrix = [[0] * size for _ in range(size)]
    for i in range(size):
        for j in range(i + 1, size):
            try:
                result = edlib.align(
                    named[i][1],
                    named[j][1],
                    mode="NW",
                    task="distance",
                    additionalEqualities=_IUPAC_EQUALITIES if ambiguous else None,
                )
            except (AttributeError, TypeError, ValueError) as exc:
                raise _edlib_error(exc, "compare two sequences") from exc
            matrix[i][j] = matrix[j][i] = int(result["editDistance"])
    return {
        "status": "completed",
        "operation": "calculate_edit_distance_matrix",
        "names": [name for name, _ in named],
        "matrix": matrix,
        "pair_count": size * (size - 1) // 2,
        "iupac": ambiguous,
        "provenance": {
            "package": _package("edlib"),
            "distance": "global (NW) edit distance: mismatches, insertions and deletions cost 1",
        },
    }


# parasail


def _parasail_settings(
    *,
    mode: object,
    matrix: object,
    match: object,
    mismatch: object,
    gap_open: object,
    gap_extend: object,
    max_cells: object,
) -> dict[str, Any]:
    return {
        "mode": _choose(mode, _PARASAIL_MODES, "mode"),
        "matrix": _choose(matrix, _PARASAIL_MATRICES, "matrix"),
        "match": _integer(match, "match", minimum=-1_000, maximum=1_000),
        "mismatch": _integer(mismatch, "mismatch", minimum=-1_000, maximum=1_000),
        "gap_open": _integer(gap_open, "gap_open", minimum=0, maximum=1_000),
        "gap_extend": _integer(gap_extend, "gap_extend", minimum=0, maximum=1_000),
        "max_cells": _integer(max_cells, "max_cells", minimum=1, maximum=_MAX_CELLS),
    }


def _parasail_pair(query: object, target: object, *, max_cells: int) -> tuple[str, str]:
    query_text = _letters(query, "query", allowed=_DNA_LETTERS, maximum=_SEQUENCE_LIMIT)
    target_text = _letters(target, "target", allowed=_DNA_LETTERS, maximum=_SEQUENCE_LIMIT)
    if len(query_text) * len(target_text) > max_cells:
        raise ConfigurationError(f"alignment needs more than max_cells={max_cells} DP cells.")
    return query_text, target_text


def _parasail_run(
    parasail: Any, settings: Mapping[str, Any], query: str, target: str, *, traceback: bool
) -> Any:
    kind = _PARASAIL_MODES[settings["mode"]]
    name = f"{kind}_trace_scan_sat" if traceback else f"{kind}_scan_sat"
    function = getattr(parasail, name, None)
    if not callable(function):
        raise ConfigurationError(f"The installed parasail lacks {name}().")
    if settings["matrix"] == "dnafull":
        scoring = parasail.dnafull
    else:
        scoring = parasail.matrix_create(
            _PARASAIL_ALPHABET, settings["match"], settings["mismatch"]
        )
    try:
        result = function(query, target, settings["gap_open"], settings["gap_extend"], scoring)
    except (TypeError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"parasail could not align the sequences: {exc}") from exc
    if bool(result.saturated):
        raise BackendExecutionError("parasail reported a saturated score for these inputs.")
    return result


def _aligned_region(query_aligned: str, target_aligned: str) -> tuple[int, int, int, int, int, int]:
    """Half-open residue spans and column bounds where both sequences have residues."""

    both = [
        index
        for index, (q_char, t_char) in enumerate(zip(query_aligned, target_aligned, strict=True))
        if q_char != "-" and t_char != "-"
    ]
    if not both:
        return 0, 0, 0, 0, 0, -1
    first, last = both[0], both[-1]
    q_before = sum(1 for char in query_aligned[:first] if char != "-")
    t_before = sum(1 for char in target_aligned[:first] if char != "-")
    q_length = sum(1 for char in query_aligned[first : last + 1] if char != "-")
    t_length = sum(1 for char in target_aligned[first : last + 1] if char != "-")
    return q_before, q_before + q_length, t_before, t_before + t_length, first, last


def _column_counts(
    query_aligned: str, target_aligned: str, first: int, last: int
) -> dict[str, int]:
    counts = {"matches": 0, "mismatches": 0, "insertions": 0, "deletions": 0}
    for q_char, t_char in zip(
        query_aligned[first : last + 1], target_aligned[first : last + 1], strict=True
    ):
        if q_char == "-":
            counts["deletions"] += 1
        elif t_char == "-":
            counts["insertions"] += 1
        elif q_char == t_char:
            counts["matches"] += 1
        else:
            counts["mismatches"] += 1
    return counts


def align_dna_parasail(
    query: str | DNA,
    target: str | DNA,
    *,
    mode: str = "global",
    matrix: str = "custom",
    match: int = 5,
    mismatch: int = -4,
    gap_open: int = 10,
    gap_extend: int = 1,
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align two DNA sequences with parasail: global, local or semi-global, affine gaps."""

    settings = _parasail_settings(
        mode=mode,
        matrix=matrix,
        match=match,
        mismatch=mismatch,
        gap_open=gap_open,
        gap_extend=gap_extend,
        max_cells=max_cells,
    )
    query_text, target_text = _parasail_pair(query, target, max_cells=settings["max_cells"])
    parasail = require("parasail", _PARASAIL_EXTRA)
    result = _parasail_run(parasail, settings, query_text, target_text, traceback=True)
    trace = result.get_traceback()
    query_aligned = str(trace.query)
    target_aligned = str(trace.ref)
    cigar_object = result.get_cigar()
    if settings["mode"] == "global":
        first, last = 0, len(query_aligned) - 1
        q_start, q_end, t_start, t_end = 0, len(query_text), 0, len(target_text)
    else:
        # The traceback starts at the CIGAR begin offsets, so add them to the region spans.
        offset_query, offset_target = int(cigar_object.beg_query), int(cigar_object.beg_ref)
        q_begin, q_finish, t_begin, t_finish, first, last = _aligned_region(
            query_aligned, target_aligned
        )
        q_start, q_end = q_begin + offset_query, q_finish + offset_query
        t_start, t_end = t_begin + offset_target, t_finish + offset_target
    counts = (
        _column_counts(query_aligned, target_aligned, first, last)
        if last >= first
        else {"matches": 0, "mismatches": 0, "insertions": 0, "deletions": 0}
    )
    columns = last - first + 1 if last >= first else 0
    found = columns > 0
    return {
        "status": "completed",
        "operation": "align_dna_parasail",
        "mode": settings["mode"],
        "matrix": settings["matrix"],
        "score": int(result.score),
        "found": found,
        "query_length": len(query_text),
        "target_length": len(target_text),
        "query_start": q_start if found else None,
        "query_end": q_end if found else None,
        "target_start": t_start if found else None,
        "target_end": t_end if found else None,
        "query_aligned": query_aligned,
        "target_aligned": target_aligned,
        "cigar": cigar_object.decode.decode("ascii"),
        "matches": counts["matches"],
        "mismatches": counts["mismatches"],
        "insertions": counts["insertions"],
        "deletions": counts["deletions"],
        "alignment_columns": columns,
        "identity": counts["matches"] / columns if columns else None,
        "provenance": {
            "package": _package("parasail"),
            "gap_penalty": "a gap of length k costs gap_open + (k - 1) * gap_extend",
            "cigar": "'=' match, 'X' mismatch, 'I' query base absent from target, "
            "'D' target base absent from query; the CIGAR covers the whole traceback",
            "coordinate_system": "zero-based half-open spans; local and semi-global spans "
            "exclude free end gaps",
            "alphabet": "A, C, G, T and N only",
        },
    }


def _parasail_score(
    parasail: Any, settings: Mapping[str, Any], query: str, target: str
) -> tuple[int, int | None, int | None]:
    result = _parasail_run(parasail, settings, query, target, traceback=False)
    score = int(result.score)
    if settings["mode"] == "local" and score <= 0:
        return score, None, None
    return score, int(result.end_query) + 1, int(result.end_ref) + 1


def score_dna_parasail(
    query: str | DNA,
    target: str | DNA,
    *,
    mode: str = "global",
    matrix: str = "custom",
    match: int = 5,
    mismatch: int = -4,
    gap_open: int = 10,
    gap_extend: int = 1,
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Score two DNA sequences with parasail without a traceback and report the optimal end cell."""

    settings = _parasail_settings(
        mode=mode,
        matrix=matrix,
        match=match,
        mismatch=mismatch,
        gap_open=gap_open,
        gap_extend=gap_extend,
        max_cells=max_cells,
    )
    query_text, target_text = _parasail_pair(query, target, max_cells=settings["max_cells"])
    parasail = require("parasail", _PARASAIL_EXTRA)
    score, query_end, target_end = _parasail_score(parasail, settings, query_text, target_text)
    return {
        "status": "completed",
        "operation": "score_dna_parasail",
        "mode": settings["mode"],
        "matrix": settings["matrix"],
        "score": score,
        "query_end": query_end,
        "target_end": target_end,
        "query_length": len(query_text),
        "target_length": len(target_text),
        "provenance": {
            "package": _package("parasail"),
            "end_cell": "zero-based half-open end of the optimal cell; None for a local score of 0",
            "gap_penalty": "a gap of length k costs gap_open + (k - 1) * gap_extend",
        },
    }


def score_dna_parasail_matrix(
    sequences: DNA | Mapping[str, str],
    *,
    mode: str = "global",
    matrix: str = "custom",
    match: int = 5,
    mismatch: int = -4,
    gap_open: int = 10,
    gap_extend: int = 1,
    max_sequences: int = 50,
    max_length: int = 2_000,
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute a symmetric pairwise parasail score matrix for named DNA sequences."""

    settings = _parasail_settings(
        mode=mode,
        matrix=matrix,
        match=match,
        mismatch=mismatch,
        gap_open=gap_open,
        gap_extend=gap_extend,
        max_cells=max_cells,
    )
    count_limit = _integer(max_sequences, "max_sequences", minimum=1, maximum=1_000)
    length_limit = _integer(max_length, "max_length", minimum=1, maximum=_SEQUENCE_LIMIT)
    named = _named_sequences(
        sequences,
        allowed=_DNA_LETTERS,
        max_count=count_limit,
        limit_name="max_sequences",
        max_length=length_limit,
    )
    parasail = require("parasail", _PARASAIL_EXTRA)
    size = len(named)
    scores = [[0] * size for _ in range(size)]
    for i in range(size):
        for j in range(i, size):
            query_text, target_text = _parasail_pair(
                named[i][1], named[j][1], max_cells=settings["max_cells"]
            )
            score, _query_end, _target_end = _parasail_score(
                parasail, settings, query_text, target_text
            )
            scores[i][j] = scores[j][i] = score
    return {
        "status": "completed",
        "operation": "score_dna_parasail_matrix",
        "mode": settings["mode"],
        "matrix": settings["matrix"],
        "names": [name for name, _ in named],
        "scores": scores,
        "pair_count": size * (size + 1) // 2,
        "provenance": {
            "package": _package("parasail"),
            "symmetry": "only the upper triangle is computed; the matrix is mirrored",
        },
    }


# pyhmmer


def _finite_number(value: object) -> float | None:
    return _finite(value) if value is not None else None


def _hmm_name(value: object) -> str:
    if not isinstance(value, str) or not _HMM_NAME_PATTERN.fullmatch(value):
        raise ConfigurationError(
            "HMM name must be 1 to 128 ASCII letters, digits, '.', '_' or '-'; pass name=."
        )
    return value


def _read_alignment(
    pyhmmer: Any, source: Path, alphabet: Any, sequence_limit: int, column_limit: int
) -> Any:
    try:
        with pyhmmer.easel.MSAFile(str(source), digital=True, alphabet=alphabet) as handle:
            msa = handle.read()
    except (OSError, ValueError, pyhmmer.errors.EaselError) as exc:
        raise ConfigurationError(
            f"pyhmmer could not read {source.name} as an aligned DNA file: {exc}"
        ) from exc
    if msa is None or not msa.names:
        raise ConfigurationError(f"no aligned sequences found in {source.name}.")
    if len(msa.names) > sequence_limit:
        raise ConfigurationError(f"alignment has more than max_sequences={sequence_limit} rows.")
    columns = len(msa.alignment[0])
    if columns > column_limit:
        raise ConfigurationError(f"alignment has more than max_columns={column_limit} columns.")
    return msa


def build_dna_hmm(
    alignment: str | Path,
    *,
    name: str | None = None,
    output_path: str | Path | None = None,
    overwrite: bool = False,
    max_sequences: int = 10_000,
    max_columns: int = 100_000,
    max_input_bytes: int = 100_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a DNA profile HMM from an aligned DNA file with pyhmmer, optionally writing HMMER3."""

    source_text = _path_argument(alignment, "alignment")
    sequence_limit = _integer(max_sequences, "max_sequences", minimum=1, maximum=1_000_000)
    column_limit = _integer(max_columns, "max_columns", minimum=1, maximum=10_000_000)
    byte_limit = _integer(max_input_bytes, "max_input_bytes", minimum=1, maximum=10_000_000_000)
    replace = _flag(overwrite, "overwrite")
    source = input_path(source_text)
    if source.stat().st_size > byte_limit:
        raise ConfigurationError(f"{source.name} exceeds max_input_bytes={byte_limit}.")
    hmm_name = _hmm_name(name if name is not None else source.stem)
    target = (
        _output_file(output_path, "output_path", overwrite=replace)
        if output_path is not None
        else None
    )
    pyhmmer = require("pyhmmer", _PYHMMER_EXTRA)
    alphabet = pyhmmer.easel.Alphabet.dna()
    msa = _read_alignment(pyhmmer, source, alphabet, sequence_limit, column_limit)
    msa.name = hmm_name.encode("ascii")
    try:
        hmm, _profile, _optimized = pyhmmer.plan7.Builder(alphabet).build_msa(
            msa, pyhmmer.plan7.Background(alphabet)
        )
    except (ValueError, pyhmmer.errors.EaselError) as exc:
        raise ConfigurationError(
            f"pyhmmer could not build an HMM from {source.name}: {exc}"
        ) from exc
    if target is not None:
        with target.open("wb" if replace else "xb") as handle:
            hmm.write(handle, binary=False)
    return {
        "status": "completed",
        "operation": "build_dna_hmm",
        "source": str(source),
        "name": str(hmm.name),
        "alphabet": "DNA",
        "length": int(hmm.M),
        "nseq": int(hmm.nseq),
        "effective_nseq": _finite_number(hmm.nseq_effective),
        "consensus": str(hmm.consensus),
        "output": str(target) if target is not None else None,
        "provenance": {
            "package": _package("pyhmmer"),
            "builder": "pyhmmer Builder and Background with default settings",
            "coordinate_system": "model length is the number of alignment columns (1-based)",
            "consensus": "the profile consensus line, as returned by pyhmmer HMM.consensus",
        },
    }


def _single_profile(pyhmmer: Any, path: Path) -> Any:
    profiles: list[Any] = []
    try:
        with pyhmmer.plan7.HMMFile(str(path)) as handle:
            for hmm in handle:
                profiles.append(hmm)
                if len(profiles) > 1:
                    break
    except (OSError, ValueError, pyhmmer.errors.EaselError) as exc:
        raise ConfigurationError(f"pyhmmer could not read {path.name} as HMMER3: {exc}") from exc
    if len(profiles) != 1:
        raise ConfigurationError(f"{path.name} must contain exactly one HMMER profile.")
    return profiles[0]


def search_dna_hmm(
    profile: str | Path,
    sequences: str | Path,
    *,
    evalue: float = 10.0,
    strand: str = "both",
    max_hits: int = 1_000,
    max_sequences: int = 100_000,
    max_residues: int = 100_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Search a DNA sequence file with a DNA profile HMM using pyhmmer's nhmmer pipeline."""

    e_value = _real(evalue, "evalue", minimum=1e-300, maximum=1_000_000.0)
    strand_name = _choose(strand, ("both", "forward", "reverse"), "strand")
    hit_limit = _integer(max_hits, "max_hits", minimum=1, maximum=100_000)
    sequence_limit = _integer(max_sequences, "max_sequences", minimum=1, maximum=1_000_000)
    residue_limit = _integer(max_residues, "max_residues", minimum=1, maximum=1_000_000_000)
    profile_path = input_path(_path_argument(profile, "profile"))
    sequence_path = input_path(_path_argument(sequences, "sequences"))
    pyhmmer = require("pyhmmer", _PYHMMER_EXTRA)
    hmm = _single_profile(pyhmmer, profile_path)
    if not hmm.alphabet.is_dna():
        raise ConfigurationError(
            f"search_dna_hmm needs a DNA profile; {profile_path.name} uses {hmm.alphabet.type}."
        )
    targets: list[Any] = []
    residues = 0
    try:
        with pyhmmer.easel.SequenceFile(
            str(sequence_path), digital=True, alphabet=hmm.alphabet
        ) as handle:
            for sequence in handle:
                if len(targets) >= sequence_limit:
                    raise ConfigurationError(f"sequences exceed max_sequences={sequence_limit}.")
                residues += len(sequence)
                if residues > residue_limit:
                    raise ConfigurationError(f"sequences exceed max_residues={residue_limit}.")
                targets.append(sequence)
    except ConfigurationError:
        raise
    except (OSError, ValueError, pyhmmer.errors.EaselError) as exc:
        raise ConfigurationError(f"pyhmmer could not read {sequence_path.name}: {exc}") from exc
    if not targets:
        raise ConfigurationError(f"no sequences found in {sequence_path.name}.")
    options: dict[str, Any] = {"E": e_value}
    if strand_name != "both":
        options["strand"] = {"forward": "watson", "reverse": "crick"}[strand_name]
    try:
        top_hits = next(iter(pyhmmer.nhmmer(hmm, targets, cpus=1, **options)))
    except (ValueError, pyhmmer.errors.EaselError) as exc:
        raise ConfigurationError(f"pyhmmer nhmmer search failed: {exc}") from exc
    rows: list[dict[str, Any]] = []
    for hit in top_hits:
        if not hit.reported:
            continue
        for domain in hit.domains:
            if not domain.reported:
                continue
            alignment = domain.alignment
            start, end = sorted((int(alignment.target_from), int(alignment.target_to)))
            rows.append(
                {
                    "target": str(hit.name),
                    "start": start,
                    "end": end,
                    "strand": str(domain.strand),
                    "hmm_start": int(alignment.hmm_from),
                    "hmm_end": int(alignment.hmm_to),
                    "score": _finite_number(domain.score),
                    "evalue": _finite_number(domain.i_evalue),
                    "sequence_evalue": _finite_number(hit.evalue),
                }
            )
            if len(rows) > hit_limit:
                raise ConfigurationError(f"search reported more than max_hits={hit_limit} domains.")
    return {
        "status": "completed",
        "operation": "search_dna_hmm",
        "profile": str(hmm.name),
        "profile_length": int(hmm.M),
        "sequences": str(sequence_path),
        "sequences_searched": len(targets),
        "evalue_threshold": e_value,
        "strand": strand_name,
        "hit_count": len(rows),
        "hits": rows,
        "provenance": {
            "package": _package("pyhmmer"),
            "algorithm": "nhmmer long-target pipeline with an explicit E-value reporting threshold",
            "coordinate_system": "1-based inclusive target coordinates; start <= end",
            "strand": "'+' watson (top) strand, '-' crick (bottom) strand",
        },
    }


def read_hmm_profile(
    path: str | Path, *, max_profiles: int = 1_000, verbose: bool = True
) -> dict[str, Any]:
    """Read HMMER3 profile metadata (name, length, consensus, sequence counts) without searching."""

    limit = _integer(max_profiles, "max_profiles", minimum=1, maximum=100_000)
    source = input_path(_path_argument(path, "path"))
    pyhmmer = require("pyhmmer", _PYHMMER_EXTRA)
    profiles: list[dict[str, Any]] = []
    try:
        with pyhmmer.plan7.HMMFile(str(source)) as handle:
            for hmm in handle:
                if len(profiles) >= limit:
                    raise ConfigurationError(f"{source.name} has more than max_profiles={limit}.")
                profiles.append(
                    {
                        "name": str(hmm.name),
                        "accession": hmm.accession,
                        "description": hmm.description,
                        "alphabet": str(hmm.alphabet.type),
                        "length": int(hmm.M),
                        "consensus": str(hmm.consensus),
                        "nseq": int(hmm.nseq),
                        "effective_nseq": _finite_number(hmm.nseq_effective),
                    }
                )
    except ConfigurationError:
        raise
    except (OSError, ValueError, pyhmmer.errors.EaselError) as exc:
        raise ConfigurationError(f"pyhmmer could not read {source.name} as HMMER3: {exc}") from exc
    if not profiles:
        raise ConfigurationError(f"no HMMER profiles found in {source.name}.")
    return {
        "status": "completed",
        "operation": "read_hmm_profile",
        "source": str(source),
        "profile_count": len(profiles),
        "profiles": profiles,
        "provenance": {
            "package": _package("pyhmmer"),
            "format": "HMMER3 profile file",
            "coordinate_system": "model positions are 1-based match columns",
        },
    }


__all__ = [
    "align_dna_edlib",
    "align_dna_parasail",
    "build_dna_hmm",
    "build_minimap2_index",
    "calculate_edit_distance_matrix",
    "fetch_minimap2_reference_region",
    "locate_approximate_pattern",
    "map_query_sequences_paf",
    "read_hmm_profile",
    "score_dna_parasail",
    "score_dna_parasail_matrix",
    "search_dna_hmm",
]

_register_public_verbose(globals(), tuple(__all__))
