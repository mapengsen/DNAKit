"""DNA storage codecs through the separately installed Chamaeleo package.

The adapter uses the fixed Goldman codec, not the full overlapping/redundant
physical archive protocol. Cite Goldman et al. (2013), doi:10.1038/nature11875,
and Chamaeleo's authors: https://github.com/ntpz870817/Chamaeleo.
"""

from __future__ import annotations

import hashlib
import re
from importlib.metadata import version
from typing import Any

from dnakit._progress import track
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError, InputFormatError

_SCHEME = "chamaeleo-goldman-fixed-v1"
_MAX_BYTES = 1_000_000


def encode_dna_storage(
    data_hex: str, *, segment_bytes: int = 32, verbose: bool = True
) -> dict[str, Any]:
    """Encode hex-encoded bytes as indexed DNA payloads with an integrity digest.

    Segment indices and metadata are stored in the returned JSON packet, not
    synthesized into the DNA. Preserve that packet for decoding. No sequencing
    error correction, wet-lab storage, or physical indexing is performed.
    Input is limited to one MB; install ``dnakit[dna-storage]`` explicitly.
    """
    if (
        not isinstance(data_hex, str)
        or len(data_hex) > 2 * _MAX_BYTES
        or len(data_hex) % 2
        or re.fullmatch(r"[0-9a-fA-F]*", data_hex) is None
    ):
        raise ConfigurationError("data_hex must encode at most one MB using pairs of hex digits.")
    if (
        isinstance(segment_bytes, bool)
        or not isinstance(segment_bytes, int)
        or not 1 <= segment_bytes <= 4096
    ):
        raise ConfigurationError("segment_bytes must be an integer in [1, 4096].")
    codec = require("Chamaeleo.methods.fixed", "dna-storage").Goldman(
        fixed_huffman=True, need_logs=False
    )
    data = bytes.fromhex(data_hex)
    sequences = []
    for index, offset in enumerate(
        track(
            range(0, len(data), segment_bytes),
            description="DNA storage encoding",
            disable=not verbose,
        )
    ):
        chunk = data[offset : offset + segment_bytes]
        bits = [int(bit) for byte in chunk for bit in f"{byte:08b}"]
        sequence = "".join(codec.encode([bits])[0])
        sequences.append({"index": index, "byte_count": len(chunk), "sequence": sequence})
    return {
        "scheme": _SCHEME,
        "backend_version": version("Chamaeleo"),
        "byte_count": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "segments": sequences,
        "index_storage": "JSON metadata, not DNA-encoded",
        "error_correction": "none",
        "citations": [
            "https://doi.org/10.1038/nature11875",
            "https://github.com/ntpz870817/Chamaeleo",
        ],
    }


def decode_dna_storage(packet: dict[str, Any], *, verbose: bool = True) -> dict[str, Any]:
    """Decode an intact Goldman packet, rejecting missing/corrupt segments.

    A SHA-256 match checks accidental corruption; it is not authentication or
    an error-correcting code. DNA-only sequences without packet metadata are
    insufficient for this interface.
    """
    if not isinstance(packet, dict) or packet.get("scheme") != _SCHEME:
        raise InputFormatError("Unsupported DNA storage packet scheme.")
    size = packet.get("byte_count")
    if isinstance(size, bool) or not isinstance(size, int) or not 0 <= size <= _MAX_BYTES:
        raise InputFormatError("Invalid packet byte_count.")
    digest = packet.get("sha256")
    if not isinstance(digest, str) or re.fullmatch(r"[0-9a-f]{64}", digest) is None:
        raise InputFormatError("Invalid packet SHA-256.")
    segments = packet.get("segments")
    if not isinstance(segments, list) or len(segments) > size:
        raise InputFormatError("Invalid packet segments.")
    by_index = {}
    total = 0
    for segment in segments:
        if not isinstance(segment, dict):
            raise InputFormatError("Each segment must be a mapping.")
        index, count, text = (segment.get(key) for key in ("index", "byte_count", "sequence"))
        if isinstance(index, bool) or not isinstance(index, int) or index < 0 or index in by_index:
            raise InputFormatError("Invalid or duplicate segment index.")
        if isinstance(count, bool) or not isinstance(count, int) or not 1 <= count <= 4096:
            raise InputFormatError("Invalid segment byte_count.")
        if (
            not isinstance(text, str)
            or not 1 <= len(text) <= count * 8
            or re.fullmatch("[ACGT]+", text) is None
        ):
            raise InputFormatError("Invalid segment DNA sequence.")
        by_index[index] = (count, text)
        total += count
    if total != size or set(by_index) != set(range(len(segments))):
        raise InputFormatError("Missing segments or inconsistent packet byte counts.")
    codec = require("Chamaeleo.methods.fixed", "dna-storage").Goldman(
        fixed_huffman=True, need_logs=False
    )
    data = bytearray()
    for index in track(
        range(len(segments)), description="DNA storage decoding", disable=not verbose
    ):
        count, text = by_index[index]
        decoded = codec.decode([list(text)])
        # Chamaeleo may silently skip invalid segments; reject them here.
        if len(decoded) != 1 or len(decoded[0]) != count * 8:
            raise InputFormatError(f"Segment {index} cannot be decoded at the declared length.")
        bits = decoded[0]
        if "".join(codec.encode([bits])[0]) != text:
            raise InputFormatError(f"Segment {index} contains truncated or invalid codewords.")
        data.extend(
            int("".join(map(str, bits[start : start + 8])), 2) for start in range(0, len(bits), 8)
        )
    if hashlib.sha256(data).hexdigest() != digest:
        raise InputFormatError("Decoded DNA payload failed its SHA-256 integrity check.")
    return {
        "data_hex": data.hex(),
        "byte_count": len(data),
        "sha256": digest,
        "integrity": "passed",
    }
