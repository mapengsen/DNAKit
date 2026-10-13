"""Offline umi-tools helpers: UMI distances and neighbours, barcode extraction and whitelists.

Each function calls umi-tools' own functions (``umi_tools.umi_methods``, ``umi_tools.network``,
``umi_tools.extract_methods``, ``umi_tools.sam_methods``, ``umi_tools.whitelist_methods`` and
``umi_tools.dedup``). Results are plain dictionaries that keep the upstream values and add a
``provenance`` entry with the package name and version. Inputs are explicit UMIs, barcodes,
reads, Phred strings, adjacency maps, count tables and local TSV or SAM/BAM files. The umi-tools
command line is not run, and no UMI-aware deduplication is performed.

``umi-tools`` is optional. It is imported only inside each function. Its modules import pysam,
matplotlib, scipy and pybktree when they are loaded, so those are checked first.
"""

from __future__ import annotations

import importlib.metadata
import re
from collections import Counter
from collections.abc import Iterable, Iterator, Mapping
from pathlib import Path
from types import SimpleNamespace
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require, sequence_text
from dnakit.exceptions import ConfigurationError

_UMI_TOOLS_EXTRA = "bio-umi-tools"
_PYSAM_EXTRA = "bio-reads"
_UMI_ALPHABET = frozenset("ACGTN")
_DNA = frozenset("ACGT")
_MAX_PHRED = 93
_QUALITY_ENCODINGS = ("phred33", "phred64", "solexa")
_KNEE_METHODS = ("distance", "density")
_ERROR_CORRECTION = ("none", "file", "derived")
_CIGAR_CODES = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5, "P": 6, "=": 7, "X": 8}
_CIGAR_PATTERN = re.compile(r"(?:\d+[MIDNSHP=X])+")
_CIGAR_OPERATION = re.compile(r"(\d+)([MIDNSHP=X])")
_MAX_UMIS = 5_000
_DEFAULT_MAX_UMIS = 2_000
_MAX_CIGAR_OPERATIONS = 10_000
_MAX_BARCODES = 1_000_000
_DEFAULT_MAX_BARCODES = 100_000
_MAX_ROWS = 10_000_000
_DEFAULT_MAX_ROWS = 1_000_000
_MAX_NODES = 1_000_000
_DEFAULT_MAX_NODES = 100_000


def _version() -> str:
    try:
        return importlib.metadata.version("umi-tools")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _umi_module(name: str) -> Any:
    """Import a umi-tools module after the packages its import chain needs."""

    require("pysam", _PYSAM_EXTRA)
    for dependency in ("matplotlib", "scipy", "pybktree"):
        require(dependency, _UMI_TOOLS_EXTRA)
    return require(name, _UMI_TOOLS_EXTRA)


def _bounded_int(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _choose(value: object, allowed: tuple[str, ...], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {', '.join(allowed)}; got {value!r}.")
    return value


def _umi(value: object, name: str, alphabet: frozenset[str] = _UMI_ALPHABET) -> str:
    if not isinstance(value, str) or not value:
        raise ConfigurationError(f"{name} must be a non-empty string of bases.")
    text = value.upper()
    if set(text) - alphabet:
        raise ConfigurationError(f"{name} must contain only {''.join(sorted(alphabet))}.")
    return text


def _umi_list(values: object, name: str, limit: int) -> list[str]:
    if isinstance(values, (str, bytes, Mapping)) or not isinstance(values, Iterable):
        raise ConfigurationError(f"{name} must be a list of UMI strings.")
    items = [_umi(item, f"{name} item") for item in values]
    if not items:
        raise ConfigurationError(f"{name} must not be empty.")
    if len(items) > limit:
        raise ConfigurationError(f"{name} exceeds max_umis={limit}; increase the explicit limit.")
    return items


def _same_length(items: Iterable[str], name: str) -> int:
    lengths = {len(item) for item in items}
    if len(lengths) != 1:
        raise ConfigurationError(f"all {name} must have the same length.")
    return lengths.pop()


def _separator(value: object) -> str:
    if (
        not isinstance(value, str)
        or not 1 <= len(value) <= 4
        or any(char.isspace() or char == "\0" for char in value)
    ):
        raise ConfigurationError("separator must be 1 to 4 non-space characters.")
    return value


def _read_name(value: object, name: str = "read_name") -> str:
    if not isinstance(value, str) or not value or any(char in value for char in ("\0", "\n", "\r")):
        raise ConfigurationError(f"{name} must be a non-empty single-line string.")
    return value


def _path_text(value: object, name: str) -> str:
    text = str(value) if isinstance(value, (str, Path)) else ""
    if not text or "\0" in text:
        raise ConfigurationError(f"{name} must be a non-empty path without NUL characters.")
    return text


def _phred(value: object, length: int, encoding: str) -> str:
    if not isinstance(value, str) or not value:
        raise ConfigurationError("qualities must be a non-empty Phred quality string.")
    if len(value) != length:
        raise ConfigurationError("qualities must contain one character per base.")
    low, high = _encoding_range(encoding)
    if any(not low <= ord(char) < high for char in value):
        raise ConfigurationError(f"qualities contain characters outside the {encoding} range.")
    return value


def _encoding_range(encoding: str) -> tuple[int, int]:
    """Printable ASCII from the encoding's lowest character up to DEL (exclusive)."""

    lowest = {"phred33": 33, "phred64": 64, "solexa": 59}
    return lowest[encoding], 127


def _provenance(method: str, **extra: str) -> dict[str, str]:
    return {"package": "umi-tools", "version": _version(), "method": method, **extra}


def average_umi_distance(
    umis: Iterable[str],
    *,
    max_umis: int = _DEFAULT_MAX_UMIS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Average pairwise Hamming distance between UMIs of one length.

    Calls ``umi_tools.umi_methods.get_average_umi_distance`` on the UMIs as bytes, which the
    upstream Hamming routine requires. Every pair is counted, so duplicate UMIs add zero
    distances. With a single UMI, upstream returns the sentinel ``-1``; this wrapper reports
    ``average_distance`` as ``None`` and ``pair_count`` as 0 instead.
    """

    umi_methods = _umi_module("umi_tools.umi_methods")
    limit = _bounded_int(max_umis, "max_umis", minimum=1, maximum=_MAX_UMIS)
    items = _umi_list(umis, "umis", limit)
    length = _same_length(items, "UMIs")
    pair_count = len(items) * (len(items) - 1) // 2
    if len(items) == 1:
        average: float | None = None
    else:
        average = float(
            umi_methods.get_average_umi_distance([item.encode("ascii") for item in items])
        )
    return {
        "status": "completed",
        "operation": "average_umi_distance",
        "umi_count": len(items),
        "umi_length": length,
        "pair_count": pair_count if len(items) > 1 else 0,
        "average_distance": average,
        "provenance": _provenance("umi_tools.umi_methods.get_average_umi_distance (Hamming)"),
    }


def find_umi_neighbours(
    umis: Iterable[str],
    *,
    max_distance: int = 1,
    max_umis: int = _DEFAULT_MAX_UMIS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find pairs of UMIs within a Hamming distance of each other.

    Duplicates are removed first. ``umi_tools.network.build_substr_idx`` splits each UMI into
    ``max_distance + 1`` substrings, so every pair within the distance shares one substring
    (pigeonhole). ``iter_nearest_neighbours`` lists the candidate pairs, and the upstream Hamming
    routine ``umi_tools.umi_methods.edit_distance`` confirms each pair. Pairs are returned with
    ``umi_a < umi_b`` and sorted.
    """

    network = _umi_module("umi_tools.network")
    umi_methods = _umi_module("umi_tools.umi_methods")
    limit = _bounded_int(max_umis, "max_umis", minimum=1, maximum=_MAX_UMIS)
    items = sorted(set(_umi_list(umis, "umis", limit)))
    length = _same_length(items, "UMIs")
    distance = _bounded_int(max_distance, "max_distance", minimum=0, maximum=min(3, length - 1))
    candidates: set[tuple[str, str]] = set()
    if len(items) > 1:
        index = network.build_substr_idx(items, length, distance)
        for first, second in network.iter_nearest_neighbours(items, index):
            candidates.add((min(first, second), max(first, second)))
    pairs = []
    for first, second in sorted(candidates):
        found = int(umi_methods.edit_distance(first.encode("ascii"), second.encode("ascii")))
        if found <= distance:
            pairs.append({"umi_a": first, "umi_b": second, "distance": found})
    return {
        "status": "completed",
        "operation": "find_umi_neighbours",
        "unique_umis": len(items),
        "umi_length": length,
        "max_distance": distance,
        "pair_count": len(pairs),
        "pairs": pairs,
        "provenance": _provenance(
            "umi_tools.network.build_substr_idx + iter_nearest_neighbours, "
            "confirmed with umi_tools.umi_methods.edit_distance (Hamming)"
        ),
    }


def umi_connected_component(
    adjacency: Mapping[str, Iterable[str]],
    node: str,
    *,
    max_nodes: int = _DEFAULT_MAX_NODES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the connected component that contains ``node`` in an undirected UMI graph.

    ``adjacency`` maps each UMI to its neighbours. Neighbours missing as keys are added with no
    edges, so the graph is complete. Calls ``umi_tools.network.breadth_first_search``, which
    collects every node reachable from ``node``. The component is returned sorted.
    """

    network = _umi_module("umi_tools.network")
    limit = _bounded_int(max_nodes, "max_nodes", minimum=1, maximum=_MAX_NODES)
    if not isinstance(adjacency, Mapping) or not adjacency:
        raise ConfigurationError("adjacency must be a non-empty mapping from UMI to neighbours.")
    if not isinstance(node, str) or not node:
        raise ConfigurationError("node must be a non-empty UMI string.")
    graph: dict[str, set[str]] = {}
    for key, neighbours in adjacency.items():
        if isinstance(neighbours, (str, bytes)) or not isinstance(neighbours, Iterable):
            raise ConfigurationError("each adjacency value must be a list of UMI strings.")
        graph.setdefault(str(key), set()).update(str(item) for item in neighbours)
    for neighbours in list(graph.values()):
        for item in neighbours:
            graph.setdefault(item, set())
    if len(graph) > limit:
        raise ConfigurationError(f"adjacency has more than max_nodes={limit} nodes.")
    if node not in graph:
        graph[node] = set()
    component = sorted(network.breadth_first_search(node, graph))
    return {
        "status": "completed",
        "operation": "umi_connected_component",
        "node": node,
        "size": len(component),
        "component": component,
        "provenance": _provenance("umi_tools.network.breadth_first_search"),
    }


def mask_low_quality_umi(
    umi: str,
    qualities: str,
    *,
    quality_encoding: str = "phred33",
    threshold: int = 10,
    verbose: bool = True,
) -> dict[str, Any]:
    """Replace UMI bases whose quality is below a threshold with ``N``.

    Calls ``umi_tools.extract_methods.mask_umi`` and ``get_below_threshold``. ``quality_encoding``
    is ``"phred33"``, ``"phred64"`` or ``"solexa"`` as in umi-tools. ``masked_positions`` lists the
    0-based positions that were below the threshold, in order.
    """

    extract = _umi_module("umi_tools.extract_methods")
    text = _umi(umi, "umi")
    encoding = _choose(quality_encoding, _QUALITY_ENCODINGS, "quality_encoding")
    quals = _phred(qualities, len(text), encoding)
    cut = _bounded_int(threshold, "threshold", minimum=0, maximum=_MAX_PHRED)
    below = extract.get_below_threshold(quals, encoding, cut)
    masked = str(extract.mask_umi(text, quals, encoding, cut))
    positions = [index for index, flag in enumerate(below) if flag]
    return {
        "status": "completed",
        "operation": "mask_low_quality_umi",
        "umi": text,
        "qualities": quals,
        "quality_encoding": encoding,
        "threshold": cut,
        "masked_umi": masked,
        "masked_positions": positions,
        "masked_count": len(positions),
        "provenance": _provenance("umi_tools.extract_methods.mask_umi + get_below_threshold"),
    }


def _positions(value: object, name: str, length: int) -> set[int]:
    if isinstance(value, (str, bytes, Mapping)) or not isinstance(value, Iterable):
        raise ConfigurationError(f"{name} must be a list of 0-based positions.")
    result: set[int] = set()
    for item in value:
        index = _bounded_int(item, f"{name} entry", minimum=0, maximum=length - 1)
        if index in result:
            raise ConfigurationError(f"{name} contains position {index} twice.")
        result.add(index)
    return result


def extract_barcode_positions(
    sequence: str,
    qualities: str,
    *,
    umi_positions: Iterable[int] = (),
    cell_positions: Iterable[int] = (),
    discard_positions: Iterable[int] = (),
    retain_umi: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove UMI, cell-barcode and discarded bases from a read by 0-based position.

    Calls ``umi_tools.extract_methods.extractSeqAndQuals``. The upstream routine returns the
    retained sequence and qualities, and the qualities of the UMI and cell bases. It does not
    return the barcode bases, so ``umi_sequence`` and ``cell_sequence`` are read from the input
    at the given positions, in read order. Position sets must not overlap.

    With ``retain_umi=True`` upstream also copies the qualities of every retained base into
    ``umi_qualities``, not only the UMI bases. The wrapper returns that value unchanged, so read
    ``umi_qualities`` with ``retain_umi=True`` with this in mind.
    """

    extract = _umi_module("umi_tools.extract_methods")
    text = sequence_text(sequence)
    quals = _phred(qualities, len(text), "phred33")
    umi = _positions(umi_positions, "umi_positions", len(text))
    cell = _positions(cell_positions, "cell_positions", len(text))
    discard = _positions(discard_positions, "discard_positions", len(text))
    if (umi & cell) or (umi & discard) or (cell & discard):
        raise ConfigurationError(
            "umi_positions, cell_positions and discard_positions must not overlap."
        )
    keep_umi = _flag(retain_umi, "retain_umi")
    new_seq, new_quals, umi_quals, cell_quals = extract.extractSeqAndQuals(
        text, quals, umi, cell, discard, keep_umi
    )
    return {
        "status": "completed",
        "operation": "extract_barcode_positions",
        "read_sequence": str(new_seq),
        "read_qualities": str(new_quals),
        "umi_sequence": "".join(text[index] for index in sorted(umi)),
        "umi_qualities": str(umi_quals),
        "cell_sequence": "".join(text[index] for index in sorted(cell)),
        "cell_qualities": str(cell_quals),
        "retain_umi": keep_umi,
        "positions": {
            "umi": sorted(umi),
            "cell": sorted(cell),
            "discard": sorted(discard),
        },
        "provenance": _provenance(
            "umi_tools.extract_methods.extractSeqAndQuals",
            barcode_sequence="read from the input at the given positions",
        ),
    }


def tag_read_name_with_umi(
    read_name: str,
    umi: str,
    *,
    cell_barcode: str | None = None,
    separator: str = "_",
    verbose: bool = True,
) -> dict[str, Any]:
    """Append a UMI, and optionally a cell barcode, to the first word of a FASTQ header.

    Calls ``umi_tools.extract_methods.addBarcodesToIdentifier``. The upstream routine splits the
    identifier on spaces and appends ``separator + [cell + separator] + umi`` to the first word.
    The rest of the header is kept unchanged.
    """

    extract = _umi_module("umi_tools.extract_methods")
    name = _read_name(read_name)
    umi_text = _umi(umi, "umi")
    cell_text = "" if cell_barcode is None else _umi(cell_barcode, "cell_barcode")
    sep = _separator(separator)
    tagged = str(
        extract.addBarcodesToIdentifier(SimpleNamespace(identifier=name), umi_text, cell_text, sep)
    )
    return {
        "status": "completed",
        "operation": "tag_read_name_with_umi",
        "read_name": name,
        "umi": umi_text,
        "cell_barcode": cell_text or None,
        "separator": sep,
        "tagged_name": tagged,
        "provenance": _provenance("umi_tools.extract_methods.addBarcodesToIdentifier"),
    }


def cigar_splice_offset(
    cigar: str,
    *,
    reverse: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the read offset of the first splice junction (N) or trailing soft clip (S) in a CIGAR.

    Calls ``umi_tools.sam_methods.find_splice`` with the operations as ``(code, length)`` pairs.
    The offset counts reference-consuming bases (M, D, =, X) and skips insertions, hard clips and
    padding. A leading soft clip is added to the offset and is not treated as a splice. Set
    ``reverse=True`` to scan from the 3' end, as umi-tools does for reverse-strand reads.
    ``splice_offset`` is ``None`` when the CIGAR has no splice.
    """

    sam = _sam_module()
    if not isinstance(cigar, str) or not _CIGAR_PATTERN.fullmatch(cigar):
        raise ConfigurationError("cigar must be a CIGAR string such as '10M100N20M'.")
    if len(cigar) > 4 * _MAX_CIGAR_OPERATIONS:
        raise ConfigurationError("cigar has too many operations.")
    operations = []
    for count, code in _CIGAR_OPERATION.findall(cigar):
        if int(count) < 1:
            raise ConfigurationError("CIGAR operation lengths must be positive.")
        operations.append((_CIGAR_CODES[code], int(count)))
    if len(operations) > _MAX_CIGAR_OPERATIONS:
        raise ConfigurationError("cigar has too many operations.")
    flag = _flag(reverse, "reverse")
    if flag:
        operations.reverse()
    offset = sam.find_splice(operations)
    splice = None if offset is False else int(offset)
    return {
        "status": "completed",
        "operation": "cigar_splice_offset",
        "cigar": cigar,
        "reverse": flag,
        "splice_offset": splice,
        "has_splice": splice is not None,
        "provenance": _provenance("umi_tools.sam_methods.find_splice"),
    }


def _sam_module() -> Any:
    require("pysam", _PYSAM_EXTRA)
    return require("umi_tools.sam_methods", _UMI_TOOLS_EXTRA)


def parse_umi_read_name(
    read_name: str,
    *,
    separator: str = "_",
    cell_barcode: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read the UMI, and optionally the cell barcode, from the end of a read name.

    Calls ``umi_tools.sam_methods.get_umi_read_string`` (UMI only) or
    ``get_cell_umi_read_string`` (cell barcode, then UMI). The upstream functions split the name
    on ``separator`` and use the last field (and the one before it). A name with too few fields
    raises ``ConfigurationError``.
    """

    sam = _sam_module()
    name = _read_name(read_name)
    sep = _separator(separator)
    flag = _flag(cell_barcode, "cell_barcode")
    try:
        if flag:
            cell_bytes, umi_bytes = sam.get_cell_umi_read_string(name, sep)
        else:
            cell_bytes, umi_bytes = sam.get_umi_read_string(name, sep)
    except ValueError as exc:
        raise ConfigurationError(f"read name does not carry the expected fields: {exc}") from exc
    return {
        "status": "completed",
        "operation": "parse_umi_read_name",
        "read_name": name,
        "separator": sep,
        "mode": "cell_and_umi" if flag else "umi_only",
        "cell_barcode": None if cell_bytes is None else bytes(cell_bytes).decode("utf-8"),
        "umi": bytes(umi_bytes).decode("utf-8"),
        "provenance": _provenance(
            "umi_tools.sam_methods.get_cell_umi_read_string"
            if flag
            else "umi_tools.sam_methods.get_umi_read_string"
        ),
    }


def _dna_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value:
        raise ConfigurationError(f"{name} must be a non-empty DNA string.")
    text = value.upper()
    if set(text) - _DNA:
        raise ConfigurationError(f"{name} must contain only A, C, G and T.")
    return text


def barcode_indel_matches(
    barcode: str,
    whitelist: Iterable[str],
    *,
    max_errors: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find whitelist barcodes within a few substitutions, insertions or deletions of a barcode.

    Calls ``umi_tools.whitelist_methods.checkError``, which uses fuzzy regex matching. The upstream
    routine skips the barcode itself and stops as soon as it has two matches. ``matches`` therefore
    holds at most two barcodes, and ``ambiguous`` is true when two were found. An exact copy of
    the barcode in the whitelist is reported by ``exact_in_whitelist`` and not as a match.
    """

    whitelist_methods = _umi_module("umi_tools.whitelist_methods")
    text = _dna_text(barcode, "barcode")
    items = sorted({_dna_text(item, "whitelist barcode") for item in whitelist})
    if not items:
        raise ConfigurationError("whitelist must not be empty.")
    errors = _bounded_int(max_errors, "max_errors", minimum=0, maximum=3)
    found = sorted(str(item) for item in whitelist_methods.checkError(text, items, errors))
    return {
        "status": "completed",
        "operation": "barcode_indel_matches",
        "barcode": text,
        "max_errors": errors,
        "whitelist_size": len(items),
        "exact_in_whitelist": text in items,
        "matches": found,
        "match_count": len(found),
        "ambiguous": len(found) > 1,
        "provenance": _provenance(
            "umi_tools.whitelist_methods.checkError (fuzzy regex; stops after two matches)"
        ),
    }


def read_barcode_whitelist(
    source: str | Path,
    *,
    error_correction: str = "none",
    threshold: int = 1,
    max_barcodes: int = _DEFAULT_MAX_BARCODES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a user-supplied cell-barcode whitelist from a local TSV file.

    Calls ``umi_tools.whitelist_methods.getUserDefinedBarcodes``. The first column of each line is
    a whitelist barcode, and ``#`` lines are comments. ``error_correction`` selects:

    - ``"none"``: the whitelist only;
    - ``"file"``: the second column lists comma-separated barcodes that correct to the first;
    - ``"derived"``: barcodes one edit away (up to ``threshold``) are mapped to the whitelist.

    The returned ``error_map`` maps each corrected barcode to its whitelist barcode.
    """

    whitelist_methods = _umi_module("umi_tools.whitelist_methods")
    path = input_path(_path_text(source, "source"))
    mode = _choose(error_correction, _ERROR_CORRECTION, "error_correction")
    distance = _bounded_int(threshold, "threshold", minimum=1, maximum=2)
    limit = _bounded_int(max_barcodes, "max_barcodes", minimum=1, maximum=_MAX_BARCODES)
    count = 0
    with path.open("r", encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            if line.startswith("#"):
                continue
            first = line.strip().split("\t")[0]
            if not first:
                raise ConfigurationError(f"line {number} of {path.name} has no barcode.")
            _dna_text(first, f"line {number} barcode")
            count += 1
            if count > limit:
                raise ConfigurationError(
                    f"{path.name} has more than max_barcodes={limit} barcodes."
                )
    whitelist, mapping = whitelist_methods.getUserDefinedBarcodes(
        str(path),
        None,
        getErrorCorrection=mode == "file",
        deriveErrorCorrection=mode == "derived",
        threshold=distance,
    )
    barcodes = sorted(str(item) for item in whitelist)
    error_map = None if mapping is None else {str(k): str(v) for k, v in sorted(mapping.items())}
    return {
        "status": "completed",
        "operation": "read_barcode_whitelist",
        "source": str(path),
        "error_correction": mode,
        "threshold": distance,
        "barcode_count": len(barcodes),
        "barcodes": barcodes,
        "error_map": error_map,
        "provenance": _provenance("umi_tools.whitelist_methods.getUserDefinedBarcodes"),
    }


def _counts(value: object, limit: int) -> Counter[str]:
    if not isinstance(value, Mapping) or not value:
        raise ConfigurationError(
            "cell_barcode_counts must be a non-empty mapping from barcode to count."
        )
    if len(value) > limit:
        raise ConfigurationError(f"cell_barcode_counts exceeds max_barcodes={limit}.")
    table: Counter[str] = Counter()
    for barcode, count in value.items():
        text = _dna_text(barcode, "cell barcode")
        if isinstance(count, bool) or not isinstance(count, int) or count < 1:
            raise ConfigurationError(f"count for barcode {text!r} must be a positive integer.")
        if text in table:
            raise ConfigurationError(f"duplicate barcode after normalization: {text!r}.")
        table[text] = count
    return table


def call_cell_whitelist(
    cell_barcode_counts: Mapping[str, int],
    *,
    knee_method: str = "distance",
    expect_cells: int | None = None,
    cell_number: int | None = None,
    error_correct_threshold: int = 0,
    max_barcodes: int = _DEFAULT_MAX_BARCODES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call true cell barcodes from barcode read counts with a knee-point method.

    Calls ``umi_tools.whitelist_methods.getCellWhitelist``. ``knee_method`` selects
    ``"distance"`` (maximum distance from the first-to-last line of the cumulative curve,
    ``getKneeEstimateDistance``) or ``"density"`` (Gaussian density, ``getKneeEstimateDensity``).
    ``expect_cells`` applies only to the density method. ``cell_number`` fixes the number of
    barcodes to accept. With ``error_correct_threshold`` above zero, the non-whitelisted barcodes
    within that Hamming distance of a called barcode are returned as ``error_map``. The upstream
    functions need a ``collections.Counter``, which this wrapper builds from the mapping.
    """

    whitelist_methods = _umi_module("umi_tools.whitelist_methods")
    limit = _bounded_int(max_barcodes, "max_barcodes", minimum=1, maximum=_MAX_BARCODES)
    counts = _counts(cell_barcode_counts, limit)
    method = _choose(knee_method, _KNEE_METHODS, "knee_method")
    if expect_cells is not None:
        if method != "density":
            raise ConfigurationError("expect_cells applies only to knee_method='density'.")
        expected = _bounded_int(expect_cells, "expect_cells", minimum=1, maximum=limit)
    else:
        expected = False
    if cell_number is None:
        fixed: int | bool = False
    else:
        fixed = _bounded_int(cell_number, "cell_number", minimum=1, maximum=len(counts))
    threshold = _bounded_int(
        error_correct_threshold, "error_correct_threshold", minimum=0, maximum=2
    )
    whitelist, mapping = whitelist_methods.getCellWhitelist(
        counts,
        knee_method=method,
        expect_cells=expected,
        cell_number=fixed,
        error_correct_threshold=threshold,
    )
    called = sorted(str(item) for item in whitelist)
    error_map = (
        None
        if mapping is None
        else {
            str(key): sorted(str(item) for item in values)
            for key, values in sorted(mapping.items())
        }
    )
    return {
        "status": "completed",
        "operation": "call_cell_whitelist",
        "knee_method": method,
        "input_barcodes": len(counts),
        "cell_count": len(called),
        "whitelist": called,
        "error_correct_threshold": threshold,
        "error_map": error_map,
        "provenance": _provenance(
            f"umi_tools.whitelist_methods.getCellWhitelist (knee_method={method!r})"
        ),
    }


def detect_bam_tags(
    source: str | Path,
    *,
    n_entries: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check whether the first mapped reads of a SAM or BAM file carry NH, X0 and XT tags.

    Calls ``umi_tools.dedup.detect_bam_features``. A tag is reported as present only if every mapped
    read examined has it. Upstream examines ``n_entries + 1`` records. The file is opened by
    upstream through pysam, which can leave a handle open until garbage collection.
    """

    require("pandas", _UMI_TOOLS_EXTRA)  # umi_tools.dedup imports pandas at module level
    dedup = _umi_module("umi_tools.dedup")
    pysam = require("pysam", _PYSAM_EXTRA)
    path = input_path(_path_text(source, "source"))
    if path.suffix.lower() not in (".bam", ".sam"):
        raise ConfigurationError("source must be a .bam or .sam file.")
    limit = _bounded_int(n_entries, "n_entries", minimum=1, maximum=1_000_000)
    try:
        found = dedup.detect_bam_features(str(path), n_entries=limit)
    except (OSError, ValueError, RuntimeError, KeyError) as exc:
        raise ConfigurationError(f"pysam could not read {path.name}: {exc}") from exc
    tags = {tag: bool(found[tag]) for tag in ("NH", "X0", "XT")}
    return {
        "status": "completed",
        "operation": "detect_bam_tags",
        "source": str(path),
        "n_entries": limit,
        "tags": tags,
        "provenance": _provenance(
            "umi_tools.dedup.detect_bam_features", pysam_version=str(pysam.__version__)
        ),
    }


def count_umis_per_gene(
    source: str | Path,
    *,
    cell_barcode: bool = False,
    separator: str = "_",
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count reads and distinct UMIs per gene from a read-name table.

    The local TSV has two columns, read name and assigned gene, with no header. Calls
    ``umi_tools.sam_methods.get_gene_count_tab`` with a getter that reads the UMI (and the cell
    barcode when ``cell_barcode=True``) from the read name. Upstream expects rows grouped by gene,
    so the wrapper merges the groups for each gene itself and does not depend on that order. For
    each gene it reports the number of reads and the number of distinct UMIs, or distinct
    cell-UMI pairs. Genes are sorted.
    """

    sam = _sam_module()
    path = input_path(_path_text(source, "source"))
    sep = _separator(separator)
    flag = _flag(cell_barcode, "cell_barcode")
    limit = _bounded_int(max_rows, "max_rows", minimum=1, maximum=_MAX_ROWS)
    rows = {"count": 0}

    def _lines() -> Iterator[str]:
        with path.open("r", encoding="utf-8") as handle:
            for number, line in enumerate(handle, start=1):
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 2 or not fields[0] or not fields[1]:
                    raise ConfigurationError(
                        f"line {number} must have two non-empty tab-separated columns: "
                        "read name and gene."
                    )
                rows["count"] += 1
                if rows["count"] > limit:
                    raise ConfigurationError(f"{path.name} has more than max_rows={limit} rows.")
                yield line

    def _getter(read_id: str) -> tuple[bytes | None, bytes]:
        try:
            if flag:
                cell, umi = sam.get_cell_umi_read_string(read_id, sep)
                return cell, umi
            _cell, umi = sam.get_umi_read_string(read_id, sep)
            return None, umi
        except ValueError as exc:
            raise ConfigurationError(
                f"read name {read_id!r} lacks the expected fields: {exc}"
            ) from exc

    totals: dict[str, dict[str, Any]] = {}
    for gene, per_cell in sam.get_gene_count_tab(_lines(), bc_getter=_getter):
        if gene is None:
            continue
        entry = totals.setdefault(str(gene), {"reads": 0, "pairs": set()})
        for cell, counter in per_cell.items():
            cell_text = None if cell is None else bytes(cell).decode("utf-8")
            for umi, reads in counter.items():
                entry["reads"] += int(reads)
                entry["pairs"].add((cell_text, bytes(umi).decode("utf-8")))
    genes = [
        {"gene": gene, "reads": entry["reads"], "unique_umis": len(entry["pairs"])}
        for gene, entry in sorted(totals.items())
    ]
    return {
        "status": "completed",
        "operation": "count_umis_per_gene",
        "source": str(path),
        "cell_barcode": flag,
        "separator": sep,
        "row_count": rows["count"],
        "gene_count": len(genes),
        "genes": genes,
        "provenance": _provenance("umi_tools.sam_methods.get_gene_count_tab"),
    }


__all__ = [
    "average_umi_distance",
    "barcode_indel_matches",
    "call_cell_whitelist",
    "cigar_splice_offset",
    "count_umis_per_gene",
    "detect_bam_tags",
    "extract_barcode_positions",
    "find_umi_neighbours",
    "mask_low_quality_umi",
    "parse_umi_read_name",
    "read_barcode_whitelist",
    "tag_read_name_with_umi",
    "umi_connected_component",
]

_register_public_verbose(globals(), tuple(__all__))
