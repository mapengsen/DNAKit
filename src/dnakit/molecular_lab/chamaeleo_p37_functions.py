"""DNA data-storage coding schemes, sequence screens and transcoding pipelines through Chamaeleo (pack P37).

Chamaeleo converts digital files to DNA and back. The bit data is a list of equal-length *bit segments*
(each a list of 0/1 integers, or a string of ``0``/``1``); DNA is a list of nucleotide strings. Every
coding scheme has an ``encode`` and a ``decode`` function with the same name stem:

- ``direct``: 2 bits per nucleotide (bit pairs to A, C, G, T).
- ``church``: Church et al. 2012, 1 bit per nucleotide with homopolymer avoidance (random choice, set ``seed``).
- ``goldman``: Goldman et al. 2013, rotating ternary Huffman code (needs 8-bit multiples).
- ``grass``: Grass et al. 2015, base-47 codewords (16-bit multiples).
- ``blawat``: Blawat et al. 2016, 8 bits to 5 nucleotides under GC and homopolymer limits.
- ``dna_fountain``: Erlich and Zielinski 2017, droplets screened by homopolymer and GC limits.
- ``yin_yang``: YinYang code, two bit segments per nucleotide sequence.

Chamaeleo is imported only when a function runs; a missing package raises ``BackendUnavailableError`` with
the ``dnakit[dna-storage]`` hint. Chamaeleo's console output is suppressed. Random schemes accept a ``seed``
and restore the global random state afterwards. Chamaeleo calls ``random.randint`` with a float bound when YinYang
adds filler segments, which Python 3.12 rejects; the call is made with an integer bound here.
"""

from __future__ import annotations

import contextlib
import copy
import io
import random
import re
import sys
from collections.abc import Mapping, Sequence
from importlib import metadata
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_EXTRA = "dna-storage"
_NUCLEOTIDES = frozenset("ACGT")
_MAX_SEGMENTS = 200_000
_MAX_BITS = 50_000_000
_SCHEMES = ("direct", "church", "goldman", "grass", "blawat", "dna_fountain", "yin_yang")
_ERROR_CORRECTIONS = (None, "hamming", "reed_solomon")
_SCHEME_OPTIONS = {
    "direct": (),
    "church": (),
    "goldman": ("support_nucleotide", "huffman_tree"),
    "grass": ("base_values",),
    "blawat": (),
    "dna_fountain": ("homopolymer", "gc_bias", "redundancy", "header_size", "c_dist", "delta", "need_pre_check", "decode_packets"),
    "yin_yang": ("yang_rule", "yin_rule", "virtual_nucleotide", "max_iterations", "max_ratio", "faster", "max_homopolymer", "max_content", "total_count", "index_length"),
}  # fmt: skip


def _provenance() -> dict[str, Any]:
    try:
        installed: str | None = metadata.version("Chamaeleo")
    except metadata.PackageNotFoundError:
        installed = None
    return {"package": "Chamaeleo", "package_version": installed}


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _integer(value: object, name: str, minimum: int, maximum: int = 10**9) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _segments(value: object, name: str = "bit_segments") -> list[list[int]]:
    """Validate bit segments: a non-empty list of equally long 0/1 lists (or '0101' strings)."""

    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(f"{name} must be a non-empty list of bit segments.")
    if len(value) > _MAX_SEGMENTS:
        raise ConfigurationError(f"{name} may hold at most {_MAX_SEGMENTS} segments.")
    rows: list[list[int]] = []
    for item in value:
        if isinstance(item, str):
            if not item or set(item) - {"0", "1"}:
                raise ConfigurationError(f"every bit segment string in {name} must use only 0 and 1.")
            rows.append([int(ch) for ch in item])
        elif isinstance(item, Sequence) and not isinstance(item, bytes) and item:
            if not all(isinstance(bit, int) and not isinstance(bit, bool) and bit in (0, 1) for bit in item):
                raise ConfigurationError(f"bit segments in {name} may contain only 0 and 1.")
            rows.append(list(item))
        else:
            raise ConfigurationError(f"every bit segment in {name} must be a non-empty list of 0/1 or a 0/1 string.")
    if len({len(row) for row in rows}) != 1:
        raise ConfigurationError(f"all bit segments in {name} must have the same length.")
    if len(rows) * len(rows[0]) > _MAX_BITS:
        raise ConfigurationError(f"{name} may hold at most {_MAX_BITS} bits.")
    return rows


def _dna(value: object, name: str = "dna_sequences") -> list[list[str]]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(f"{name} must be a non-empty list of DNA strings.")
    if len(value) > _MAX_SEGMENTS:
        raise ConfigurationError(f"{name} may hold at most {_MAX_SEGMENTS} sequences.")
    rows = []
    for item in value:
        if not isinstance(item, str) or not item or set(item) - _NUCLEOTIDES:
            raise ConfigurationError(f"every sequence in {name} must be a non-empty string of A, C, G and T (upper case).")
        rows.append(list(item))
    return rows


@contextlib.contextmanager
def _quiet(seed: int | None):
    """Silence Chamaeleo's printing and, with a seed, make random choices reproducible without leaking state."""

    state = random.getstate()
    limit = sys.getrecursionlimit()
    original_randint = random.randint

    def randint(low: float, high: float) -> int:
        # YinYang filler segments call randint with a float bound, which Python 3.12 rejects.
        return original_randint(int(low), int(high))

    if seed is not None:
        random.seed(_integer(seed, "seed", 0, 2**63))
    random.randint = randint
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            yield
    finally:
        random.randint = original_randint
        random.setstate(state)
        sys.setrecursionlimit(limit)


def _call(function: Any, *args: Any, **kwargs: Any) -> Any:
    try:
        return function(*args, **kwargs)
    except ConfigurationError:
        raise
    except (ValueError, TypeError, KeyError, IndexError, AttributeError, ZeroDivisionError) as exc:
        raise ConfigurationError(f"Chamaeleo could not process the data: {exc}") from exc


def _options(scheme: str, options: Mapping[str, Any] | None) -> dict[str, Any]:
    if options is None:
        return {}
    if not isinstance(options, Mapping):
        raise ConfigurationError("scheme_options must be a mapping.")
    unknown = set(options) - set(_SCHEME_OPTIONS[scheme])
    if unknown:
        raise ConfigurationError(f"{scheme} accepts the options {list(_SCHEME_OPTIONS[scheme])}; got {sorted(unknown)}.")
    return dict(options)


def _build(scheme: str, options: Mapping[str, Any] | None) -> Any:
    """Create a Chamaeleo coding-scheme object from a scheme name and its option mapping."""

    if scheme not in _SCHEMES:
        raise ConfigurationError(f"scheme must be one of {list(_SCHEMES)}; got {scheme!r}.")
    values = _options(scheme, options)
    default = require("Chamaeleo.methods.default", _EXTRA)
    fixed = require("Chamaeleo.methods.fixed", _EXTRA)
    flowed = require("Chamaeleo.methods.flowed", _EXTRA)
    if scheme == "direct":
        return default.BaseCodingAlgorithm()
    if scheme == "church":
        return fixed.Church()
    if scheme == "goldman":
        tree = values.pop("huffman_tree", None)
        model = _call(fixed.Goldman, fixed_huffman=True, **values)
        if tree is not None:
            if isinstance(tree, (str, bytes)) or len(tree) != 256:
                raise ConfigurationError("huffman_tree must be a list of 256 ternary code strings.")
            model.huffman_tree = list(tree)
        return model
    if scheme == "grass":
        return _call(fixed.Grass, **values)
    if scheme == "blawat":
        return fixed.Blawat()
    if scheme == "dna_fountain":
        return _call(flowed.DNAFountain, **values)
    total, index_length = values.pop("total_count", None), values.pop("index_length", None)
    model = _call(flowed.YinYangCode, **values)
    if total is not None:
        model.total_count = _integer(total, "total_count", 1)
        model.index_length = int(index_length) if index_length is not None else len(bin(total)) - 2
    return model


def _encode(scheme: str, options: Mapping[str, Any] | None, bit_segments: object, seed: int | None) -> dict[str, Any]:
    rows = _segments(bit_segments)
    model = _build(scheme, options)
    work = copy.deepcopy(rows)
    with _quiet(seed):
        result = _call(model.silicon_to_carbon, work, len(rows) * len(rows[0]))
    sequences = ["".join(item) for item in result["dna"]]
    out = {
        "status": "completed",
        "scheme": scheme,
        "sequence_count": len(sequences),
        "dna": sequences,
        "information_density": float(result["i"]),
        "bit_size": len(rows) * len(rows[0]),
        "segment_length": len(rows[0]),
        "segment_count": len(rows),
        "provenance": _provenance(),
    }
    return out


def _decode(scheme: str, options: Mapping[str, Any] | None, dna_sequences: object, segment_length: int | None, bit_size: int | None) -> dict[str, Any]:
    rows = _dna(dna_sequences)
    model = _build(scheme, options)
    model.bit_size = 0 if bit_size is None else _integer(bit_size, "bit_size", 0, _MAX_BITS)
    model.segment_length = 0 if segment_length is None else _integer(segment_length, "segment_length", 1, 10**7)
    with _quiet(None):
        result = _call(model.decode, rows)
    segments = [[int(bit) for bit in segment] for segment in result]
    if segment_length is not None:
        segments = [segment[:segment_length] for segment in segments]
    return {
        "status": "completed",
        "scheme": scheme,
        "segment_count": len(segments),
        "bit_segments": segments,
        "bit_strings": ["".join(map(str, segment)) for segment in segments],
        "provenance": _provenance(),
    }


def chamaeleo_encode_direct(bit_segments: Sequence[Sequence[int] | str], *, verbose: bool = True) -> dict[str, Any]:
    """Encode bit segments into DNA with the 2-bit direct mapping (00 to A, 01 to C, 10 to G, 11 to T).

    Segment lengths must be even. The result holds the DNA strings and the information density in bits per nucleotide.
    """

    return {"operation": "chamaeleo_encode_direct"} | _encode("direct", None, bit_segments, None)


def chamaeleo_decode_direct(dna_sequences: Sequence[str], *, verbose: bool = True) -> dict[str, Any]:
    """Decode DNA sequences back into bit segments with the 2-bit direct mapping."""

    return {"operation": "chamaeleo_decode_direct"} | _decode("direct", None, dna_sequences, None, None)


def chamaeleo_encode_church(bit_segments: Sequence[Sequence[int] | str], *, seed: int | None = None, verbose: bool = True) -> dict[str, Any]:
    """Encode bits with Church et al. 2012: bit 0 is A or C, bit 1 is G or T, and a fourth identical base in a row is avoided.

    The choice between the two bases is random; pass ``seed`` for a reproducible result.
    """

    return {"operation": "chamaeleo_encode_church"} | _encode("church", None, bit_segments, seed)


def chamaeleo_decode_church(dna_sequences: Sequence[str], *, verbose: bool = True) -> dict[str, Any]:
    """Decode Church et al. 2012 DNA (A or C is 0, G or T is 1) back into bit segments."""

    return {"operation": "chamaeleo_decode_church"} | _decode("church", None, dna_sequences, None, None)


def chamaeleo_encode_goldman(
    bit_segments: Sequence[Sequence[int] | str],
    *,
    fixed_huffman: bool = True,
    support_nucleotide: str = "A",
    verbose: bool = True,
) -> dict[str, Any]:
    """Encode bytes with Goldman et al. 2013: a ternary Huffman code whose trits pick the next base by rotation.

    Segment lengths must be multiples of 8. ``fixed_huffman=True`` uses the published code. ``False`` builds an adaptive
    tree from the data (``chamaeleo_goldman_huffman_tree``), uses it for encoding and returns it as ``huffman_tree``;
    ``chamaeleo_decode_goldman`` needs it. (Chamaeleo's own adaptive path passes the argument ``3`` as the data size
    and fails for most inputs, so the tree is built here with the intended size and multiple.)
    """

    adaptive = not _flag(fixed_huffman, "fixed_huffman")
    options: dict[str, Any] = {"support_nucleotide": support_nucleotide}
    tree = None
    if adaptive:
        tree = chamaeleo_goldman_huffman_tree(bit_segments, verbose=False)["huffman_tree"]
        options["huffman_tree"] = tree
    result = _encode("goldman", options, bit_segments, None)
    if tree is not None:
        result["huffman_tree"] = tree
    return {"operation": "chamaeleo_encode_goldman"} | result


def chamaeleo_decode_goldman(
    dna_sequences: Sequence[str],
    *,
    support_nucleotide: str = "A",
    huffman_tree: Sequence[str | None] | None = None,
    segment_length: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Decode Goldman et al. 2013 DNA into bit segments; pass the ``huffman_tree`` returned by an adaptive encoding (default: the published code)."""

    options: dict[str, Any] = {"support_nucleotide": support_nucleotide}
    if huffman_tree is not None:
        options["huffman_tree"] = list(huffman_tree)
    return {"operation": "chamaeleo_decode_goldman"} | _decode("goldman", options, dna_sequences, segment_length, None)


def chamaeleo_goldman_huffman_tree(bit_segments: Sequence[Sequence[int] | str], *, multiple: int = 3, verbose: bool = True) -> dict[str, Any]:
    """Build the adaptive ternary Huffman code of Goldman coding from the byte frequencies of the data (256 codes, ``None`` for unused bytes)."""

    rows = _segments(bit_segments)
    multiple = _integer(multiple, "multiple", 2, 16)
    model = _build("goldman", None)
    with _quiet(None):
        tree = _call(model.adaptive_huffman_tree, copy.deepcopy(rows), None, multiple)
    return {
        "status": "completed",
        "operation": "chamaeleo_goldman_huffman_tree",
        "multiple": multiple,
        "huffman_tree": list(tree),
        "used_byte_count": sum(code is not None for code in tree),
        "provenance": _provenance(),
    }


def chamaeleo_encode_grass(bit_segments: Sequence[Sequence[int] | str], *, base_values: Sequence[int] | None = None, verbose: bool = True) -> dict[str, Any]:
    """Encode 16-bit blocks into three base-47 codewords (Grass et al. 2015, 9 nucleotides per block).

    Segment lengths must be multiples of 16. ``base_values`` is an optional custom list of 48 base values.
    """

    options = {} if base_values is None else {"base_values": list(base_values)}
    return {"operation": "chamaeleo_encode_grass"} | _encode("grass", options, bit_segments, None)


def chamaeleo_decode_grass(
    dna_sequences: Sequence[str],
    *,
    base_values: Sequence[int] | None = None,
    segment_length: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Decode Grass et al. 2015 DNA (9 nucleotides per 16 bits) into bit segments."""

    options = {} if base_values is None else {"base_values": list(base_values)}
    return {"operation": "chamaeleo_decode_grass"} | _decode("grass", options, dna_sequences, segment_length, None)


def chamaeleo_encode_blawat(bit_segments: Sequence[Sequence[int] | str], *, seed: int | None = None, verbose: bool = True) -> dict[str, Any]:
    """Encode 8-bit blocks into 5 nucleotides (Blawat et al. 2016) under GC-content and homopolymer constraints.

    Segment lengths must be multiples of 8. Where several codewords satisfy the constraints one is chosen at random; set ``seed``.
    """

    return {"operation": "chamaeleo_encode_blawat"} | _encode("blawat", None, bit_segments, seed)


def chamaeleo_decode_blawat(dna_sequences: Sequence[str], *, segment_length: int | None = None, verbose: bool = True) -> dict[str, Any]:
    """Decode Blawat et al. 2016 DNA (5 nucleotides per byte) into bit segments."""

    return {"operation": "chamaeleo_decode_blawat"} | _decode("blawat", None, dna_sequences, segment_length, None)


def chamaeleo_encode_dna_fountain(
    bit_segments: Sequence[Sequence[int] | str],
    *,
    homopolymer: int = 4,
    gc_bias: float = 0.2,
    redundancy: float = 0.07,
    header_size: int = 4,
    c_dist: float = 0.1,
    delta: float = 0.05,
    need_pre_check: bool = False,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Encode bit segments as DNA Fountain droplets (Erlich and Zielinski 2017) screened by homopolymer and GC limits.

    Droplets are generated until ``ceil(segments * (1 + redundancy))`` valid sequences exist. Each starts with a
    ``header_size``-byte seed. Segment lengths must be even (odd ones get a leading 0). ``need_pre_check=True`` test-decodes first.
    Returns ``decode_packets``, the segment count that ``chamaeleo_decode_dna_fountain`` needs.
    """

    options = {
        "homopolymer": _integer(homopolymer, "homopolymer", 1, 100),
        "gc_bias": gc_bias,
        "redundancy": redundancy,
        "header_size": _integer(header_size, "header_size", 0, 16),
        "c_dist": c_dist,
        "delta": delta,
        "need_pre_check": _flag(need_pre_check, "need_pre_check"),
    }
    result = _encode("dna_fountain", options, bit_segments, seed)
    result["decode_packets"] = result["segment_count"]
    return {"operation": "chamaeleo_encode_dna_fountain"} | result


def chamaeleo_decode_dna_fountain(
    dna_sequences: Sequence[str],
    decode_packets: int,
    *,
    header_size: int = 4,
    c_dist: float = 0.1,
    delta: float = 0.05,
    segment_length: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Decode DNA Fountain droplets by peeling; ``decode_packets`` is the original segment count and must match the encoding parameters."""

    options = {
        "decode_packets": _integer(decode_packets, "decode_packets", 1, _MAX_SEGMENTS),
        "header_size": _integer(header_size, "header_size", 0, 16),
        "c_dist": c_dist,
        "delta": delta,
    }
    return {"operation": "chamaeleo_decode_dna_fountain"} | _decode("dna_fountain", options, dna_sequences, segment_length, None)


def chamaeleo_encode_yin_yang(
    bit_segments: Sequence[Sequence[int] | str],
    *,
    faster: bool = False,
    max_iterations: int = 100,
    max_ratio: float = 0.8,
    max_homopolymer: int = 4,
    max_content: float = 0.6,
    virtual_nucleotide: str = "A",
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Encode bit segments with the YinYang code, which merges two bit segments into one DNA sequence.

    ``faster=False`` is the normal setting (good and bad segments are paired deliberately); ``faster=True`` pairs segments at
    random and may add more filler sequences. Sequences obey ``max_homopolymer`` and ``max_content`` (GC limit).
    Returns ``total_count``, the segment count that ``chamaeleo_decode_yin_yang`` needs to drop filler sequences.
    """

    options = {
        "faster": _flag(faster, "faster"),
        "max_iterations": _integer(max_iterations, "max_iterations", 1, 10**6),
        "max_ratio": max_ratio,
        "max_homopolymer": _integer(max_homopolymer, "max_homopolymer", 1, 100),
        "max_content": max_content,
        "virtual_nucleotide": virtual_nucleotide,
    }
    result = _encode("yin_yang", options, bit_segments, seed)
    result["total_count"] = result["segment_count"]
    return {"operation": "chamaeleo_encode_yin_yang"} | result


def chamaeleo_decode_yin_yang(
    dna_sequences: Sequence[str],
    total_count: int,
    *,
    virtual_nucleotide: str = "A",
    segment_length: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Decode YinYang DNA into bit segments; filler segments with an index above ``total_count`` are dropped.

    The decoded segments come in pairs (upper then lower bit segment) and are not sorted; the index prefix of each segment
    (the first ``bit_length(total_count)`` bits) gives its position when the data was encoded with an index.
    """

    options = {"total_count": _integer(total_count, "total_count", 1, _MAX_SEGMENTS), "virtual_nucleotide": virtual_nucleotide}
    return {"operation": "chamaeleo_decode_yin_yang"} | _decode("yin_yang", options, dna_sequences, segment_length, None)


# Screens and files ------------------------------------------------------------------------------------


def chamaeleo_screen_dna(
    sequence: str,
    *,
    max_homopolymer: int | None = None,
    max_gc_deviation: float | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Screen one DNA sequence for homopolymer runs and GC content (``screen.homopolymer``, ``screen.gc_content``, ``screen.check``).

    The homopolymer test fails when a base repeats more than ``max_homopolymer`` times in a row. The GC test uses
    Chamaeleo's ``max_content`` ``m`` (here ``0.5 + max_gc_deviation``) and passes when ``1 - m <= GC fraction <= m``.
    Give either limit or both (the combined ``check``); a limit left as None is not tested.
    """

    if not isinstance(sequence, str) or not sequence or set(sequence) - _NUCLEOTIDES:
        raise ConfigurationError("sequence must be a non-empty upper-case DNA string of A, C, G and T.")
    screen = require("Chamaeleo.utils.screen", _EXTRA)
    homo_limit = None if max_homopolymer is None else _integer(max_homopolymer, "max_homopolymer", 1, 10**6)
    if max_gc_deviation is not None and (isinstance(max_gc_deviation, bool) or not isinstance(max_gc_deviation, (int, float)) or not 0 < max_gc_deviation <= 0.5):
        raise ConfigurationError("max_gc_deviation must be a number in (0, 0.5].")
    content_limit = None if max_gc_deviation is None else 0.5 + float(max_gc_deviation)
    homopolymer_ok = None if homo_limit is None else bool(screen.homopolymer(sequence, homo_limit))
    gc_ok = None if content_limit is None else bool(screen.gc_content(sequence, content_limit))
    return {
        "status": "completed",
        "operation": "chamaeleo_screen_dna",
        "length": len(sequence),
        "gc_fraction": (sequence.count("G") + sequence.count("C")) / len(sequence),
        "longest_homopolymer": max(len(m.group(0)) for m in re.finditer(r"A+|C+|G+|T+", sequence)),
        "homopolymer_ok": homopolymer_ok,
        "gc_ok": gc_ok,
        "passes": bool(screen.check(sequence, homo_limit, content_limit)),
        "provenance": _provenance(),
    }


def chamaeleo_read_dna_file(path: str | Path, *, max_sequences: int = 1_000_000, verbose: bool = True) -> dict[str, Any]:
    """Read DNA sequences from a text file with one sequence per line (``read_dna_file``)."""

    source = input_path(path)
    limit = _integer(max_sequences, "max_sequences", 1, 10_000_000)
    data_handle = require("Chamaeleo.utils.data_handle", _EXTRA)
    with _quiet(None):
        rows = data_handle.read_dna_file(str(source), need_logs=False)
    if len(rows) > limit:
        raise ConfigurationError(f"{source.name} has {len(rows)} sequences; max_sequences={limit}.")
    sequences = ["".join(row) for row in rows]
    return {
        "status": "completed",
        "operation": "chamaeleo_read_dna_file",
        "sequence_count": len(sequences),
        "sequences": sequences,
        "provenance": _provenance(),
    }


def chamaeleo_write_dna_file(sequences: Sequence[str], output_path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Write DNA sequences to a text file, one per line (``write_dna_file``). The file must not exist yet."""

    rows = _dna(sequences, "sequences")
    target = Path(output_path).expanduser()
    if target.exists():
        raise ConfigurationError(f"output_path already exists: {target.name}; choose a new file.")
    if not target.parent.is_dir():
        raise ConfigurationError(f"The folder of output_path does not exist: {target.parent}.")
    data_handle = require("Chamaeleo.utils.data_handle", _EXTRA)
    with _quiet(None):
        data_handle.write_dna_file(str(target), rows, need_logs=False)
    return {"status": "completed", "operation": "chamaeleo_write_dna_file", "output_path": str(target), "sequence_count": len(rows)}


# Pipelines --------------------------------------------------------------------------------------------


def _error_correction(name: str | None, check_size: int) -> Any:
    if name not in _ERROR_CORRECTIONS:
        raise ConfigurationError(f"error_correction must be one of {list(_ERROR_CORRECTIONS)}.")
    ecc = require("Chamaeleo.methods.ecc", _EXTRA)
    if name is None:
        return None
    if name == "hamming":
        return ecc.Hamming()
    return ecc.ReedSolomon(check_size=_integer(check_size, "check_size", 1, 250))


def _new(path: str | Path, name: str) -> Path:
    target = Path(path).expanduser()
    if target.exists():
        raise ConfigurationError(f"{name} already exists: {target.name}; choose a new file.")
    if not target.parent.is_dir():
        raise ConfigurationError(f"The folder of {name} does not exist: {target.parent}.")
    return target


def chamaeleo_transcode_file_to_dna(
    input_path_: str | Path,
    output_path: str | Path,
    *,
    scheme: str = "direct",
    scheme_options: Mapping[str, Any] | None = None,
    error_correction: str | None = None,
    check_size: int = 3,
    segment_length: int = 120,
    index: bool = True,
    index_length: int | None = None,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Transcode a digital file to DNA (``TranscodePipeline.transcode(direction='t_c')``).

    The file bytes are cut into ``segment_length``-bit segments (the last one zero-padded), numbered with an ``index``
    prefix, protected with ``error_correction`` (``hamming`` or ``reed_solomon`` with ``check_size`` bytes), encoded with
    ``scheme`` (one of direct, church, goldman, grass, blawat, dna_fountain, yin_yang; ``scheme_options`` are its constructor
    options) and written one sequence per line. The result holds ``decode_parameters`` for ``chamaeleo_transcode_dna_to_file``.
    """

    source = input_path(input_path_)
    target = _new(output_path, "output_path")
    seg = _integer(segment_length, "segment_length", 1, 10**6)
    index = _flag(index, "index")
    model = _build(scheme, scheme_options)
    correction = _error_correction(error_correction, check_size)
    pipelines = require("Chamaeleo.utils.pipelines", _EXTRA)
    pipeline = pipelines.TranscodePipeline(coding_scheme=model, error_correction=correction, need_logs=False)
    kwargs: dict[str, Any] = {"direction": "t_c", "input_path": str(source), "output_path": str(target), "segment_length": seg, "index": index}
    if index_length is not None:
        kwargs["index_length"] = _integer(index_length, "index_length", 1, 64)
    with _quiet(seed):
        data = _call(pipeline.transcode, **kwargs)
    payload_count = len(data["bit"])
    result = {
        "status": "completed",
        "operation": "chamaeleo_transcode_file_to_dna",
        "scheme": scheme,
        "output_path": str(target),
        "sequence_count": len(data["dna"]),
        "records": {key: value for key, value in pipeline.records.items()},
        "decode_parameters": {
            "bit_size": int(model.bit_size),
            "segment_length": int(model.segment_length),
            "payload_length": seg,
            "payload_count": payload_count,
            "index": index,
            "index_length": pipeline.records.get("index length", 0),
            "error_correction": error_correction,
            "check_size": check_size,
        },
        "provenance": _provenance(),
    }
    return result


def chamaeleo_transcode_dna_to_file(
    input_path_: str | Path,
    output_path: str | Path,
    *,
    bit_size: int,
    segment_length: int,
    scheme: str = "direct",
    scheme_options: Mapping[str, Any] | None = None,
    error_correction: str | None = None,
    check_size: int = 3,
    index: bool = True,
    index_length: int | None = None,
    payload_count: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Transcode DNA (one sequence per line) back to a digital file (``TranscodePipeline.transcode(direction='t_s')``).

    ``bit_size`` (the original file size in bits) and the other parameters are the ``decode_parameters`` returned by
    ``chamaeleo_transcode_file_to_dna``; ``segment_length`` is the length of an encoded bit segment including the index
    and error-correction bits (not the payload length). ``dna_fountain`` also needs ``decode_packets`` and ``yin_yang`` needs ``total_count`` in
    ``scheme_options`` (``payload_count`` fills these in when they are missing).
    """

    source = input_path(input_path_)
    target = _new(output_path, "output_path")
    options = dict(scheme_options or {})
    if payload_count is not None:
        if scheme == "dna_fountain":
            options.setdefault("decode_packets", payload_count)
        elif scheme == "yin_yang":
            options.setdefault("total_count", payload_count)
    model = _build(scheme, options)
    model.bit_size = _integer(bit_size, "bit_size", 1, _MAX_BITS)
    model.segment_length = _integer(segment_length, "segment_length", 1, 10**6)
    index = _flag(index, "index")
    correction = _error_correction(error_correction, check_size)
    pipelines = require("Chamaeleo.utils.pipelines", _EXTRA)
    pipeline = pipelines.TranscodePipeline(coding_scheme=model, error_correction=correction, need_logs=False)
    kwargs: dict[str, Any] = {"direction": "t_s", "input_path": str(source), "output_path": str(target), "index": index}
    if index_length is not None:
        kwargs["index_length"] = _integer(index_length, "index_length", 1, 64)
    with _quiet(None):
        data = _call(pipeline.transcode, **kwargs)
    if data["bit"] is None or not target.exists():
        raise ConfigurationError("the DNA could not be decoded back into a file (all sequences failed error correction or decoding).")
    return {
        "status": "completed",
        "operation": "chamaeleo_transcode_dna_to_file",
        "scheme": scheme,
        "output_path": str(target),
        "bytes_written": target.stat().st_size,
        "records": {key: value for key, value in pipeline.records.items()},
        "provenance": _provenance(),
    }


def _scheme_dict(schemes: Sequence[str], options: Mapping[str, Mapping[str, Any]] | None) -> dict[str, Any]:
    if isinstance(schemes, (str, bytes)) or not schemes:
        raise ConfigurationError("coding_schemes must be a non-empty list of scheme names.")
    if len(set(schemes)) != len(schemes):
        raise ConfigurationError("coding_schemes must not repeat a scheme.")
    return {name: _build(name, (options or {}).get(name)) for name in schemes}


def _file_dict(files: Mapping[str, str | Path]) -> dict[str, str]:
    if not isinstance(files, Mapping) or not files:
        raise ConfigurationError("file_paths must be a non-empty mapping from a label to a file path.")
    return {str(label): str(input_path(path)) for label, path in files.items()}


def chamaeleo_evaluate_robustness(
    file_paths: Mapping[str, str | Path],
    coding_schemes: Sequence[str],
    *,
    error_corrections: Sequence[str | None] = (None,),
    check_size: int = 3,
    needed_indices: Sequence[bool] | None = None,
    scheme_options: Mapping[str, Mapping[str, Any]] | None = None,
    nucleotide_insertion: float = 0.0,
    nucleotide_mutation: float = 0.0,
    nucleotide_deletion: float = 0.0,
    sequence_loss: float = 0.0,
    iterations: int = 1,
    segment_length: int = 120,
    index_length: int | None = None,
    seed: int | None = None,
    log_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Simulate DNA-storage errors and measure how much of each file is recovered (``RobustnessPipeline.evaluate``).

    Every file is encoded with every scheme and error-correction option, a fraction of the sequences is lost, and the
    given rates of nucleotide insertions, mutations and deletions (0 to 1, per sequence count) are applied before decoding.
    ``needed_indices`` (one flag per scheme, default all True) says whether a scheme gets an index. Each task reports the
    information density, runtimes and the recovery success rate per iteration. ``log_path`` (a new file) receives the CSV log
    that ``chamaeleo_choose_best_scheme`` reads.
    """

    files = _file_dict(file_paths)
    models = _scheme_dict(coding_schemes, scheme_options)
    if isinstance(error_corrections, (str, bytes)) or not error_corrections:
        raise ConfigurationError("error_corrections must be a non-empty list (use None for no correction).")
    corrections = {("none" if name is None else name): _error_correction(name, check_size) for name in error_corrections}
    flags = [True] * len(models) if needed_indices is None else [_flag(item, "needed_indices") for item in needed_indices]
    rates = {}
    for key, value in (("nucleotide_insertion", nucleotide_insertion), ("nucleotide_mutation", nucleotide_mutation), ("nucleotide_deletion", nucleotide_deletion), ("sequence_loss", sequence_loss)):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0 <= value <= 1:
            raise ConfigurationError(f"{key} must be a number in [0, 1].")
        rates[key] = float(value)
    target = None if log_path is None else _new(log_path, "log_path")
    pipelines = require("Chamaeleo.utils.pipelines", _EXTRA)
    pipeline = pipelines.RobustnessPipeline(
        coding_schemes=models, error_corrections=corrections, needed_indices=flags, file_paths=files,
        segment_length=_integer(segment_length, "segment_length", 1, 10**6), index_length=index_length,
        iterations=_integer(iterations, "iterations", 1, 1000), need_logs=False, **rates,
    )  # fmt: skip
    with _quiet(seed):
        _call(pipeline.evaluate)
        if target is not None:
            pipeline.output_records(type="path", path=str(target))
    tasks = []
    for task_id, task in pipeline.records.get("results", {}).items():
        tasks.append(
            {
                "task": task_id,
                "coding_scheme": task["coding scheme"],
                "error_correction": task["error-correction"],
                "file": task["file"],
                "iterations": [
                    {
                        "information_density": item.get("information density"),
                        "transcoding_state": item.get("transcoding state"),
                        "success_count": item.get("success count"),
                        "total_count": item.get("total count"),
                        "success_rate": item.get("success rate"),
                        "error_rate": item.get("error rate"),
                    }
                    for item in task["result"]
                ],
            }
        )
    return {
        "status": "completed",
        "operation": "chamaeleo_evaluate_robustness",
        "perturbation": {**rates, "iterations": iterations},
        "task_count": len(tasks),
        "tasks": tasks,
        "log_path": None if target is None else str(target),
        "provenance": _provenance(),
    }


def chamaeleo_profile_coding_features(
    file_paths: Mapping[str, str | Path],
    coding_schemes: Sequence[str],
    *,
    needed_indices: Sequence[bool] | None = None,
    scheme_options: Mapping[str, Mapping[str, Any]] | None = None,
    segment_length: int = 120,
    index_length: int | None = None,
    seed: int | None = None,
    log_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count the GC-content and homopolymer-length distributions of the DNA each coding scheme makes from each file (``BasicFeaturePipeline.calculate``).

    ``gc_distribution`` has 101 counts (sequences per rounded GC percent) and ``homopolymer_distribution`` counts
    sequences by their longest homopolymer run. ``log_path`` (a new file) receives the log for ``chamaeleo_choose_best_scheme``.
    """

    files = _file_dict(file_paths)
    models = _scheme_dict(coding_schemes, scheme_options)
    flags = [True] * len(models) if needed_indices is None else [_flag(item, "needed_indices") for item in needed_indices]
    target = None if log_path is None else _new(log_path, "log_path")
    pipelines = require("Chamaeleo.utils.pipelines", _EXTRA)
    pipeline = pipelines.BasicFeaturePipeline(
        coding_schemes=models, needed_indices=flags, file_paths=files,
        segment_length=_integer(segment_length, "segment_length", 1, 10**6), index_length=index_length, need_logs=False,
    )  # fmt: skip
    with _quiet(seed):
        _call(pipeline.calculate)
        if target is not None:
            pipeline.output_records(type="path", path=str(target))
    tasks = []
    for task_id, task in pipeline.records.get("results", {}).items():
        tasks.append(
            {
                "task": task_id,
                "coding_scheme": task["coding scheme"],
                "file": task["file"],
                "gc_distribution": [int(x) for x in task["gc"].strip("[]").split("-")],
                "homopolymer_distribution": [int(x) for x in task["homo"].strip("[]").split("-")],
            }
        )
    return {
        "status": "completed",
        "operation": "chamaeleo_profile_coding_features",
        "task_count": len(tasks),
        "tasks": tasks,
        "log_path": None if target is None else str(target),
        "provenance": _provenance(),
    }


def chamaeleo_choose_best_scheme(feature_log_path: str | Path, robustness_log_path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Compare coding schemes from the logs of ``chamaeleo_profile_coding_features`` and ``chamaeleo_evaluate_robustness`` (``OptimalChoicePipeline.calculate_best``).

    Returns one row per scheme with the average information density, the share of sequences with 40 to 60 percent GC, the
    maximum GC bias, the maximum homopolymer length and the average recovery rate, as Chamaeleo's table prints them.
    Chamaeleo's parser reads logs written with the default 101-bin GC distribution.
    """

    features, robustness = input_path(feature_log_path), input_path(robustness_log_path)
    pipelines = require("Chamaeleo.utils.pipelines", _EXTRA)
    default = require("Chamaeleo.methods.default", _EXTRA)
    pipeline = pipelines.OptimalChoicePipeline(
        coding_schemes={"direct": default.BaseCodingAlgorithm()}, feature_log_path=str(features),
        robustness_log_path=str(robustness), need_logs=False,
    )  # fmt: skip
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        _call(pipeline.calculate_best)
    lines = [line for line in buffer.getvalue().splitlines() if line.startswith("|")]
    if len(lines) < 2:
        raise ConfigurationError("Chamaeleo produced no comparison table from these logs.")
    header = [cell.strip() for cell in lines[0].strip("|").split("|")]
    rows = [dict(zip(header, (cell.strip() for cell in line.strip("|").split("|")), strict=False)) for line in lines[1:]]
    return {
        "status": "completed",
        "operation": "chamaeleo_choose_best_scheme",
        "scheme_count": len(rows),
        "schemes": rows,
        "provenance": _provenance(),
    }


__all__ = [
    "chamaeleo_choose_best_scheme",
    "chamaeleo_decode_blawat",
    "chamaeleo_decode_church",
    "chamaeleo_decode_direct",
    "chamaeleo_decode_dna_fountain",
    "chamaeleo_decode_goldman",
    "chamaeleo_decode_grass",
    "chamaeleo_decode_yin_yang",
    "chamaeleo_encode_blawat",
    "chamaeleo_encode_church",
    "chamaeleo_encode_direct",
    "chamaeleo_encode_dna_fountain",
    "chamaeleo_encode_goldman",
    "chamaeleo_encode_grass",
    "chamaeleo_encode_yin_yang",
    "chamaeleo_evaluate_robustness",
    "chamaeleo_goldman_huffman_tree",
    "chamaeleo_profile_coding_features",
    "chamaeleo_read_dna_file",
    "chamaeleo_screen_dna",
    "chamaeleo_transcode_dna_to_file",
    "chamaeleo_transcode_file_to_dna",
    "chamaeleo_write_dna_file",
]

_register_public_verbose(globals(), tuple(__all__))
