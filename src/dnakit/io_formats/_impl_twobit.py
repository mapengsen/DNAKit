"""UCSC 2bit engine, implemented natively from the format description.

Layout (all integers 32-bit, in the byte order given by the signature 0x1A412743)::

    header   signature, version (0), sequenceCount, reserved (0)
    index    sequenceCount x { nameSize:u8, name:bytes, offset:u32 }   (offset from file start)
    record   dnaSize, nBlockCount, nBlockStarts[], nBlockSizes[], maskBlockCount, maskBlockStarts[],
             maskBlockSizes[], reserved (0), packedDna[ceil(dnaSize/4)]

Four bases are packed per byte, first base in the high bits, with T=0, C=1, A=2, G=3.  Runs of N
are stored as "N blocks" and lower-case (soft-masked) runs as "mask blocks"; any non-ACGT letter is
written as N (ambiguity codes are not preserved, like UCSC ``faToTwoBit``).  Only version 0
(32-bit offsets) is supported.
"""

# ruff: noqa: E501, RUF012, B905

from __future__ import annotations

import os
import re
import struct
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, BinaryIO

from dnakit.exceptions import InputFormatError
from dnakit.io_formats._impl_common import (
    DEFAULT_PAGE,
    Format,
    Opts,
    Payload,
    bad_input,
    bad_option,
    chrom_key,
    decode_text,
    json_payload,
    looks_like_path,
)

SIGNATURE = 0x1A412743
MAX_SEQUENCES = 5_000_000
_CODE = {"T": 0, "C": 1, "A": 2, "G": 3}
_BASE = "TCAG"
_BYTE_TO_BASES = [
    "".join(_BASE[(b >> shift) & 3] for shift in (6, 4, 2, 0)).encode("ascii") for b in range(256)
]
_KMER_TO_BYTE = {"".join(_BASE[(b >> shift) & 3] for shift in (6, 4, 2, 0)): b for b in range(256)}
# per-byte base counts, columns T C A G
_BYTE_COUNTS = [
    tuple(sum(1 for shift in (6, 4, 2, 0) if (b >> shift) & 3 == code) for code in range(4))
    for b in range(256)
]


@dataclass
class TwoBitSeq:
    name: str
    dna_size: int
    n_blocks: list[tuple[int, int]] = field(default_factory=list)
    mask_blocks: list[tuple[int, int]] = field(default_factory=list)
    packed: bytes = b""


def _runs(sequence: str, pattern: str) -> list[tuple[int, int]]:
    return [(m.start(), m.end() - m.start()) for m in re.finditer(pattern, sequence)]


def encode_sequence(name: str, sequence: str) -> TwoBitSeq:
    """Build the record of one sequence (non-ACGT -> N block, lower case -> mask block)."""

    n_blocks = _runs(sequence, r"[^ACGTacgt]+")
    mask_blocks = _runs(sequence, r"[a-z]+")
    upper = sequence.upper()
    kmers = re.sub(r"[^ACGT]", "T", upper)
    remainder = len(kmers) % 4
    if remainder:
        kmers += "T" * (4 - remainder)
    packed = bytes(_KMER_TO_BYTE[kmers[i : i + 4]] for i in range(0, len(kmers), 4))
    return TwoBitSeq(name, len(sequence), n_blocks, mask_blocks, packed)


def _pack_record(seq: TwoBitSeq, endian: str) -> bytes:
    def u32(*values: int) -> bytes:
        return struct.pack(f"{endian}{len(values)}I", *values)

    parts = [
        u32(seq.dna_size, len(seq.n_blocks)),
        u32(*[s for s, _ in seq.n_blocks]),
        u32(*[n for _, n in seq.n_blocks]),
        u32(len(seq.mask_blocks)),
        u32(*[s for s, _ in seq.mask_blocks]),
        u32(*[n for _, n in seq.mask_blocks]),
        u32(0),
        seq.packed,
    ]
    return b"".join(parts)


def build_file(seqs: Sequence[TwoBitSeq], byte_order: str = "little") -> bytes:
    endian = "<" if byte_order == "little" else ">"
    names = [s.name.encode("utf-8") for s in seqs]
    for raw in names:
        if not raw or len(raw) > 255:
            raise bad_input("2bit sequence names must be 1-255 bytes long.")
    if len(set(names)) != len(names):
        raise bad_input("2bit sequence names must be unique.", code="IO_FORMATS_DUPLICATE_ID")
    index_size = sum(1 + len(raw) + 4 for raw in names)
    offset = 16 + index_size
    records = [_pack_record(s, endian) for s in seqs]
    if offset + sum(len(r) for r in records) >= 2**32:
        raise bad_input("The 2bit file would exceed 4 GiB (version 0 limit).")
    out = [struct.pack(f"{endian}4I", SIGNATURE, 0, len(seqs), 0)]
    for raw, record in zip(names, records):
        out.append(bytes([len(raw)]) + raw + struct.pack(f"{endian}I", offset))
        offset += len(record)
    return b"".join(out + records)


def decode_slice(
    raw: bytes,
    first_byte: int,
    n_blocks: Sequence[tuple[int, int]],
    mask_blocks: Sequence[tuple[int, int]],
    start: int,
    end: int,
    soft_mask: bool = True,
) -> str:
    """Decode bases ``[start, end)`` from ``raw`` = packed bytes ``first_byte ..`` (blocks use absolute coordinates)."""

    buffer = bytearray(b"".join(_BYTE_TO_BASES[b] for b in raw))
    base = first_byte * 4
    lo, hi = start - base, end - base
    for s, n in n_blocks:
        a, b = max(s - base, lo), min(s + n - base, hi)
        if a < b:
            buffer[a:b] = b"N" * (b - a)
    if soft_mask:
        for s, n in mask_blocks:
            a, b = max(s - base, lo), min(s + n - base, hi)
            if a < b:
                buffer[a:b] = bytes(buffer[a:b]).lower()
    return buffer[lo:hi].decode("ascii")


class TwoBitFile:
    """Random-access reader over bytes or a file; nothing larger than a record header is read eagerly."""

    def __init__(self, source: bytes | Path) -> None:
        self._bytes = source if isinstance(source, bytes) else None
        self._path = None if isinstance(source, bytes) else source
        self._handle: BinaryIO | None = None
        self.size = len(source) if isinstance(source, bytes) else source.stat().st_size
        try:
            self._read_index()
        except BaseException:
            self.close()
            raise

    def _read_index(self) -> None:
        header = self.read_at(0, 16)
        if len(header) < 16:
            raise bad_input("File is too short to be a 2bit file.", code="IO_FORMATS_BAD_TWOBIT")
        sig_le = struct.unpack("<I", header[:4])[0]
        if sig_le == SIGNATURE:
            self.endian = "<"
        elif struct.unpack(">I", header[:4])[0] == SIGNATURE:
            self.endian = ">"
        else:
            raise bad_input(
                "Not a 2bit file: signature 0x1A412743 not found.",
                code="IO_FORMATS_BAD_TWOBIT",
                hint="2bit files start with the bytes 43 27 41 1A (or 1A 41 27 43 on big-endian files).",
            )
        _, self.version, self.count, self.reserved = struct.unpack(f"{self.endian}4I", header)
        if self.version != 0:
            raise bad_input(
                f"Unsupported 2bit version {self.version}; only version 0 is supported.",
                code="IO_FORMATS_BAD_TWOBIT",
                hint="Convert 64-bit-offset (version 1) files with UCSC twoBitToFa.",
            )
        if self.count > MAX_SEQUENCES:
            raise bad_input(
                f"2bit file declares {self.count} sequences (limit {MAX_SEQUENCES}).",
                code="IO_FORMATS_BAD_TWOBIT",
            )
        self.entries: list[tuple[str, int]] = []
        self.index_end = 16
        position = 16
        # the index is read in one go: at most count * (1 + 255 + 4) bytes
        blob = self.read_at(16, min(self.size - 16, self.count * 260, 64 * 1024 * 1024))
        cursor = 0
        for _ in range(self.count):
            if cursor >= len(blob):
                raise bad_input("2bit index is truncated.", code="IO_FORMATS_BAD_TWOBIT")
            length = blob[cursor]
            end = cursor + 1 + length + 4
            if end > len(blob):
                raise bad_input("2bit index is truncated.", code="IO_FORMATS_BAD_TWOBIT")
            name = blob[cursor + 1 : cursor + 1 + length].decode("utf-8", errors="replace")
            offset = struct.unpack(f"{self.endian}I", blob[end - 4 : end])[0]
            self.entries.append((name, offset))
            cursor = end
        self.index_end = position + cursor

    def read_at(self, offset: int, size: int) -> bytes:
        if self._bytes is not None:
            return self._bytes[offset : offset + size]
        if self._handle is None:
            self._handle = open(self._path, "rb")  # noqa: SIM115 - closed in close()
        self._handle.seek(offset)
        return self._handle.read(size)

    def close(self) -> None:
        if self._handle is not None:
            self._handle.close()
            self._handle = None

    def _u32s(self, offset: int, count: int) -> list[int]:
        raw = self.read_at(offset, 4 * count)
        if len(raw) != 4 * count:
            raise bad_input("2bit record is truncated.", code="IO_FORMATS_BAD_TWOBIT")
        return list(struct.unpack(f"{self.endian}{count}I", raw)) if count else []

    def meta(self, offset: int) -> dict[str, Any]:
        """Parse a record header (everything except the packed bases)."""

        dna_size, n_count = self._u32s(offset, 2)
        if n_count * 8 > self.size:
            raise bad_input(
                "2bit record declares an impossible number of N blocks.",
                code="IO_FORMATS_BAD_TWOBIT",
            )
        pos = offset + 8
        n_starts = self._u32s(pos, n_count)
        n_sizes = self._u32s(pos + 4 * n_count, n_count)
        pos += 8 * n_count
        (m_count,) = self._u32s(pos, 1)
        if m_count * 8 > self.size:
            raise bad_input(
                "2bit record declares an impossible number of mask blocks.",
                code="IO_FORMATS_BAD_TWOBIT",
            )
        pos += 4
        m_starts = self._u32s(pos, m_count)
        m_sizes = self._u32s(pos + 4 * m_count, m_count)
        pos += 8 * m_count
        (reserved,) = self._u32s(pos, 1)
        return {
            "dna_size": dna_size,
            "n_blocks": list(zip(n_starts, n_sizes)),
            "mask_blocks": list(zip(m_starts, m_sizes)),
            "reserved": reserved,
            "packed_offset": pos + 4,
            "packed_length": (dna_size + 3) // 4,
            "offset": offset,
        }

    def sequence_meta(self) -> list[tuple[str, dict[str, Any]]]:
        return [(name, self.meta(offset)) for name, offset in self.entries]

    def packed(self, meta: Mapping[str, Any]) -> bytes:
        raw = self.read_at(meta["packed_offset"], meta["packed_length"])
        if len(raw) != meta["packed_length"]:
            raise bad_input(
                "2bit record is truncated (packed bases missing).", code="IO_FORMATS_BAD_TWOBIT"
            )
        return raw

    def decode(
        self,
        meta: Mapping[str, Any],
        start: int = 0,
        end: int | None = None,
        soft_mask: bool = True,
    ) -> str:
        size = meta["dna_size"]
        end = size if end is None else end
        if not 0 <= start <= end <= size:
            raise bad_option(f"Region {start}-{end} is outside the sequence (length {size}).")
        if start == end:
            return ""
        first, last = start // 4, (end + 3) // 4
        raw = self.read_at(meta["packed_offset"] + first, last - first)
        if len(raw) != last - first:
            raise bad_input(
                "2bit record is truncated (packed bases missing).", code="IO_FORMATS_BAD_TWOBIT"
            )
        return decode_slice(
            raw, first, meta["n_blocks"], meta["mask_blocks"], start, end, soft_mask
        )

    def load(self) -> list[TwoBitSeq]:
        out = []
        for name, meta in self.sequence_meta():
            out.append(
                TwoBitSeq(
                    name, meta["dna_size"], meta["n_blocks"], meta["mask_blocks"], self.packed(meta)
                )
            )
        return out


def _count_range(packed: bytes, start: int, end: int) -> list[int]:
    """Counts of the raw 2-bit codes (T, C, A, G) of bases ``[start, end)``, ignoring N/mask blocks."""

    if end <= start:
        return [0, 0, 0, 0]
    text = b"".join(_BYTE_TO_BASES[b] for b in packed[start // 4 : (end + 3) // 4])
    text = text[start % 4 : start % 4 + (end - start)]
    return [text.count(b"T"), text.count(b"C"), text.count(b"A"), text.count(b"G")]


def composition(seq: TwoBitSeq) -> dict[str, int]:
    """Counts of A, C, G, T (excluding N) plus N and soft-masked bases of one record."""

    counts = [0, 0, 0, 0]  # T C A G
    for value, times in Counter(seq.packed).items():
        for code in range(4):
            counts[code] += _BYTE_COUNTS[value][code] * times
    for code, number in enumerate(_count_range(seq.packed, seq.dna_size, 4 * len(seq.packed))):
        counts[code] -= number  # padding bits of the last byte
    for start, length in seq.n_blocks:
        for code, number in enumerate(
            _count_range(seq.packed, start, min(start + length, seq.dna_size))
        ):
            counts[code] -= number  # positions covered by N blocks
    counts = [max(c, 0) for c in counts]
    return {
        "A": counts[2],
        "C": counts[1],
        "G": counts[3],
        "T": counts[0],
        "N": sum(n for _, n in seq.n_blocks),
        "soft_masked": sum(n for _, n in seq.mask_blocks),
    }


def _check_record(name: str, meta: Mapping[str, Any], file_size: int, errors: list[str]) -> None:
    size = meta["dna_size"]
    for label, blocks in (("N", meta["n_blocks"]), ("mask", meta["mask_blocks"])):
        previous_end = 0
        for start, length in blocks:
            if length == 0:
                errors.append(f"{name}: empty {label} block at {start}")
            if start < previous_end:
                errors.append(f"{name}: {label} blocks overlap or are not sorted at {start}")
            if start + length > size:
                errors.append(
                    f"{name}: {label} block {start}+{length} extends past the sequence (length {size})"
                )
            previous_end = max(previous_end, start + length)
    if meta["reserved"] != 0:
        errors.append(f"{name}: record reserved field is {meta['reserved']}, expected 0")
    if meta["packed_offset"] + meta["packed_length"] > file_size:
        errors.append(f"{name}: packed bases extend past the end of the file")


def _bounded_total(sizes: Sequence[int], opts: Opts) -> None:
    limit = opts.integer("max_bases", 50_000_000, 1) or 50_000_000
    if sum(sizes) > limit:
        raise bad_option(
            f"The selection has {sum(sizes)} bases, above max_bases={limit}.",
            hint="Raise max_bases or select a region.",
        )


class TwoBitFormat(Format):
    name = "twobit"
    title = "UCSC 2bit genome file"
    noun = "sequences"
    suffix = ".2bit"
    kind = "binary"
    targets = ("twobit", "fasta", "json")
    native_ops = frozenset({"read", "validate", "summarize", "index", "stream"})
    options = {
        "read": {
            "names": "list of sequence names to read (default all)",
            "start": "int: 0-based start (single name only)",
            "end": "int: 0-based exclusive end (single name only)",
            "soft_mask": "bool, default True: lower-case masked bases",
            "max_bases": "int, default 50000000: refuse to decode more bases",
        },
        "write": {"byte_order": "'little' (default) or 'big'"},
        "validate": {},
        "filter": {
            "names": "list of sequence names to keep",
            "exclude": "list of names to drop",
            "name_pattern": "regular expression (re.search, max 256 characters) the name must match",
            "min_length": "int: keep sequences with at least this many bases",
            "max_length": "int: at most this many bases",
            "max_n_fraction": "float in [0,1]: drop sequences whose N fraction exceeds it",
        },
        "sort": {
            "key": "'name' (default), 'natural' (chr1..22, X, Y, M order) or 'length'",
            "reverse": "bool, default False",
        },
        "merge": {
            "duplicates": "'error' (default) or 'first': what to do when a name occurs twice"
        },
        "split": {"sequences_per_part": "int, default 1: sequences in each output part"},
        "stream": {
            "window": "int, default 10000: bases per window (the unit that is paged)",
            "soft_mask": "bool, default True",
            "max_bases": "int, default 50000000: bases per page",
        },
        "summarize": {},
        "index": {},
        "convert": {
            "line_width": "int, default 60: FASTA line width (0 = one line)",
            "soft_mask": "bool, default True",
            "max_bases": "int, default 50000000",
            "names": "list of sequence names to convert (default all)",
        },
    }
    for _op in ("filter", "sort", "merge", "split", "convert", "write"):
        options[_op] = {
            **options[_op],
            "byte_order": "'little' (default) or 'big': byte order of the written 2bit file",
        }
    del _op
    semantics = {
        "read": "Decodes sequences (or the half-open 0-based region [start, end) of one sequence) to text: N blocks become 'N', "
        "mask blocks lower case. Result records are {name, length, sequence}.",
        "write": "Packs sequences into a version-0 2bit file. Input: FASTA text/path, a {name: sequence} mapping or a list of "
        "{'name','sequence'} records. Non-ACGT letters become N blocks; lower case becomes mask blocks.",
        "validate": "Checks the signature, version 0, reserved fields, index offsets, unique non-empty names, that N/mask blocks "
        "are sorted, non-empty, non-overlapping and inside the sequence, and that packed data fits in the file.",
        "filter": "Keeps whole sequences (records are copied bit-exactly, no re-encoding) that pass all the given criteria.",
        "sort": "Reorders sequences in the output index; the packed records are copied bit-exactly.",
        "merge": "Concatenates the sequences of several 2bit files into one version-0 file (names must be unique).",
        "split": "Writes groups of sequences to separate 2bit files; each result part is a complete 2bit file.",
        "stream": "Pages through fixed-size windows over all sequences in file order; each record is "
        "{name, start, end, sequence}. The cursor 'start' counts windows; only the windows of a page are decoded.",
        "summarize": "Per sequence: length, A/C/G/T/N counts, soft-masked bases, GC fraction = (G+C)/(A+C+G+T), and file totals.",
        "index": "Per sequence: file offset of its record, length, number of N and mask blocks and the packed-data offset.",
        "convert": "Targets: 'fasta' (soft-masking kept as lower case), 'json' (metadata and sequences) and 'twobit' "
        "(re-encode, e.g. to change byte order). FASTA -> 2bit is convert(format='fasta', to='twobit').",
    }

    # -- loading --------------------------------------------------------------------------
    def open_source(self, source: object, opts: Opts) -> TwoBitFile:
        if isinstance(source, (bytes, bytearray, memoryview)):
            return TwoBitFile(bytes(source))
        if isinstance(source, os.PathLike) or (
            isinstance(source, str) and looks_like_path(source) and os.path.isfile(source)
        ):
            return TwoBitFile(Path(source).expanduser())
        raise bad_input(
            "A 2bit input must be a path to a .2bit file or its raw bytes.",
            hint="Text/str content cannot hold a binary 2bit file; pass bytes or a path.",
        )

    def coerce(
        self, source: object, opts: Opts, errors: list[str] | None = None
    ) -> list[TwoBitSeq]:
        if isinstance(source, (Mapping, list, tuple)):
            return self.from_records(source, opts)
        handle = self.open_source(source, opts)
        try:
            self._guard_size(handle, opts)
            return handle.load()
        finally:
            handle.close()

    def _guard_size(self, handle: TwoBitFile, opts: Opts) -> None:
        limit = opts.integer("max_bytes", 256 * 1024 * 1024, 1) or 256 * 1024 * 1024
        if handle.size > limit:
            raise bad_option(
                f"The 2bit file has {handle.size} bytes, above max_bytes={limit}.",
                hint="Raise max_bytes, or use index/read(names, start, end)/stream on the large file.",
            )

    def from_records(self, records: object, opts: Opts) -> list[TwoBitSeq]:
        if isinstance(records, Mapping):
            for key in ("records", "sequences", "items"):
                if key in records:
                    return self.from_records(records[key], opts)
            pairs = [(str(k), v) for k, v in records.items()]
        elif isinstance(records, Sequence) and not isinstance(records, (str, bytes)):
            pairs = []
            for number, item in enumerate(records, 1):
                if isinstance(item, Mapping) and "name" in item and "sequence" in item:
                    pairs.append((str(item["name"]), item["sequence"]))
                elif isinstance(item, Sequence) and not isinstance(item, str) and len(item) == 2:
                    pairs.append((str(item[0]), item[1]))
                else:
                    raise bad_input(
                        f"sequence record {number} must be a {{'name','sequence'}} dict or (name, sequence) pair."
                    )
        else:
            raise bad_input(
                "Expected FASTA text, a name->sequence mapping or a list of sequence records."
            )
        out = []
        for name, seq in pairs:
            if not isinstance(seq, str):
                raise bad_input(f"Sequence {name!r} must be a string.")
            clean = "".join(seq.split())
            if not re.fullmatch(r"[A-Za-z*.\-]*", clean):
                raise bad_input(f"Sequence {name!r} contains characters that are not letters.")
            out.append(encode_sequence(name, clean))
        return out

    def coerce_for_write(self, source: object, opts: Opts) -> list[TwoBitSeq]:
        if isinstance(source, (Mapping, list, tuple)):
            return self.from_records(source, opts)
        if isinstance(source, (bytes, bytearray)) and bytes(source[:4]) in (
            b"\x43\x27\x41\x1a",
            b"\x1a\x41\x27\x43",
        ):
            return self.coerce(source, opts)
        from dnakit.io_formats._impl_align import parse_fasta
        from dnakit.io_formats._impl_common import read_source

        data, _ = read_source(
            source, opts.integer("max_bytes", 256 * 1024 * 1024, 1) or 256 * 1024 * 1024
        )
        if data[:4] in (b"\x43\x27\x41\x1a", b"\x1a\x41\x27\x43"):
            return self.coerce(data, opts)
        aln = parse_fasta(decode_text(data), None, self.bad)
        return self.from_records(list(zip(aln.names, aln.seqs)), opts)

    def parse(self, data: bytes, errors: list[str] | None, opts: Opts) -> list[TwoBitSeq]:
        return self.coerce(data, opts)

    # -- records / dump -------------------------------------------------------------------
    def decode_seq(self, seq: TwoBitSeq, soft_mask: bool = True) -> str:
        if seq.dna_size == 0:
            return ""
        return decode_slice(
            seq.packed, 0, seq.n_blocks, seq.mask_blocks, 0, seq.dna_size, soft_mask
        )

    def records(self, doc: list[TwoBitSeq]) -> list[dict[str, object]]:
        return [{"name": s.name, "length": s.dna_size, "sequence": self.decode_seq(s)} for s in doc]

    def count(self, doc: list[TwoBitSeq]) -> int:
        return len(doc)

    def dump(self, doc: list[TwoBitSeq], opts: Opts) -> Payload:
        data = build_file(doc, str(opts.text("byte_order", "little", ["little", "big"])))
        return Payload("twobit", data=data, count=len(doc), info={"sequence_count": len(doc)})

    # -- operations that stream from the file (read / validate / summarize / index / stream) --
    def read_op(self, source: object, opts: Opts) -> dict[str, object]:
        names = opts.names("names")
        start, end = opts.integer("start", None, 0), opts.integer("end", None, 0)
        soft = opts.boolean("soft_mask", True)
        handle = self.open_source(source, opts)
        try:
            metas = dict(handle.sequence_meta())
            chosen = names if names is not None else list(metas)
            for name in chosen:
                if name not in metas:
                    raise bad_option(
                        f"Sequence {name!r} is not in the 2bit file.",
                        hint=f"Available: {list(metas)[:10]}",
                    )
            if (start is not None or end is not None) and len(chosen) != 1:
                raise bad_option("start/end need exactly one name in 'names'.")
            spans = []
            for name in chosen:
                size = metas[name]["dna_size"]
                a, b = (start or 0), (size if end is None else end)
                spans.append((name, a, b))
            _bounded_total([max(0, b - a) for _, a, b in spans], opts)
            records = [
                {
                    "name": n,
                    "start": a,
                    "end": b,
                    "length": metas[n]["dna_size"],
                    "sequence": handle.decode(metas[n], a, b, soft),
                }
                for n, a, b in spans
            ]
        finally:
            handle.close()
        return {"records": records, "record_count": len(records)}

    def validate_op(self, source: object, opts: Opts) -> dict[str, object]:
        errors: list[str] = []
        warnings: list[str] = []
        try:
            handle = self.open_source(source, opts)
        except InputFormatError as exc:
            if exc.code != "IO_FORMATS_BAD_TWOBIT":
                raise
            return {"valid": False, "errors": [exc.message], "warnings": [], "sequence_count": 0}
        try:
            return self._validate(handle, errors, warnings)
        except InputFormatError as exc:
            return {
                "valid": False,
                "errors": [*errors, exc.message],
                "warnings": warnings,
                "sequence_count": handle.count,
            }
        finally:
            handle.close()

    def _validate(
        self, handle: TwoBitFile, errors: list[str], warnings: list[str]
    ) -> dict[str, object]:
        if handle.reserved != 0:
            errors.append(f"header reserved field is {handle.reserved}, expected 0")
        names = [n for n, _ in handle.entries]
        if len(set(names)) != len(names):
            errors.append("duplicate sequence names in the index")
        if any(not n for n in names):
            errors.append("empty sequence name in the index")
        offsets = [o for _, o in handle.entries]
        if len(set(offsets)) != len(offsets):
            errors.append("two index entries share one record offset")
        for name, offset in handle.entries:
            if offset < handle.index_end or offset + 16 > handle.size:
                errors.append(f"{name}: record offset {offset} is outside the data area")
                continue
            try:
                _check_record(name, handle.meta(offset), handle.size, errors)
            except InputFormatError as exc:
                errors.append(f"{name}: {exc.message}")
        return {
            "valid": not errors,
            "errors": errors,
            "warnings": warnings,
            "sequence_count": handle.count,
            "byte_order": "little" if handle.endian == "<" else "big",
            "file_size": handle.size,
        }

    def summarize_op(self, source: object, opts: Opts) -> dict[str, object]:
        handle = self.open_source(source, opts)
        try:
            rows, totals = [], Counter()
            for name, meta in handle.sequence_meta():
                seq = TwoBitSeq(
                    name,
                    meta["dna_size"],
                    meta["n_blocks"],
                    meta["mask_blocks"],
                    handle.packed(meta),
                )
                comp = composition(seq)
                acgt = comp["A"] + comp["C"] + comp["G"] + comp["T"]
                rows.append(
                    {
                        "name": name,
                        "length": seq.dna_size,
                        **comp,
                        "gc_fraction": (comp["G"] + comp["C"]) / acgt if acgt else None,
                    }
                )
                totals.update({k: v for k, v in comp.items()})
                totals["length"] += seq.dna_size
            acgt = sum(totals[k] for k in "ACGT")
            return {
                "sequence_count": len(rows),
                "total_length": totals["length"],
                "n_bases": totals["N"],
                "soft_masked_bases": totals["soft_masked"],
                "gc_fraction": (totals["G"] + totals["C"]) / acgt if acgt else None,
                "version": handle.version,
                "byte_order": "little" if handle.endian == "<" else "big",
                "sequences": rows,
            }
        finally:
            handle.close()

    def index_op(self, source: object, opts: Opts) -> dict[str, object]:
        handle = self.open_source(source, opts)
        try:
            entries = {
                name: {
                    "offset": meta["offset"],
                    "length": meta["dna_size"],
                    "n_block_count": len(meta["n_blocks"]),
                    "mask_block_count": len(meta["mask_blocks"]),
                    "packed_offset": meta["packed_offset"],
                }
                for name, meta in handle.sequence_meta()
            }
            return {"key": "sequence name", "record_count": len(entries), "entries": entries}
        finally:
            handle.close()

    def stream_op(self, source: object, opts: Opts) -> dict[str, object]:
        window = opts.integer("window", 10_000, 1) or 10_000
        start = opts.integer("start", 0, 0) or 0
        limit = opts.integer("limit", DEFAULT_PAGE, 1, 100_000) or DEFAULT_PAGE
        soft = opts.boolean("soft_mask", True)
        handle = self.open_source(source, opts)
        try:
            metas = handle.sequence_meta()
            total = sum((m["dna_size"] + window - 1) // window for _, m in metas)
            records: list[dict[str, object]] = []
            cursor = 0
            wanted_end = start + limit
            budget = 0
            max_bases = opts.integer("max_bases", 50_000_000, 1) or 50_000_000
            for name, meta in metas:
                count = (meta["dna_size"] + window - 1) // window
                lo, hi = max(start - cursor, 0), min(wanted_end - cursor, count)
                for w in range(lo, hi):
                    a, b = w * window, min((w + 1) * window, meta["dna_size"])
                    budget += b - a
                    if budget > max_bases:
                        raise bad_option(
                            f"The page would decode more than max_bases={max_bases} bases.",
                            hint="Lower limit or window.",
                        )
                    records.append(
                        {
                            "name": name,
                            "start": a,
                            "end": b,
                            "sequence": handle.decode(meta, a, b, soft),
                        }
                    )
                cursor += count
            end = start + len(records)
            return {
                "start": start,
                "limit": limit,
                "returned": len(records),
                "total": total,
                "window": window,
                "next_start": end if end < total else None,
                "done": end >= total,
                "records": records,
            }
        finally:
            handle.close()

    # -- operations on the loaded records ---------------------------------------------------
    def _n_fraction(self, seq: TwoBitSeq) -> float:
        return sum(n for _, n in seq.n_blocks) / seq.dna_size if seq.dna_size else 0.0

    def filter(self, doc: list[TwoBitSeq], opts: Opts) -> list[TwoBitSeq]:
        names, drop = opts.name_set("names"), opts.name_set("exclude")
        pattern = opts.text("name_pattern")
        min_len, max_len = opts.integer("min_length", None, 0), opts.integer("max_length", None, 0)
        max_n = opts.number("max_n_fraction", None, 0.0, 1.0)
        rx = None
        if pattern:
            if len(pattern) > 256:
                raise bad_option("name_pattern is limited to 256 characters.")
            try:
                rx = re.compile(pattern)
            except re.error as exc:
                raise bad_option(f"name_pattern is not a valid regular expression: {exc}") from exc
        out = []
        for s in doc:
            if (names is not None and s.name not in names) or (drop is not None and s.name in drop):
                continue
            if rx is not None and not rx.search(s.name):
                continue
            if (min_len is not None and s.dna_size < min_len) or (
                max_len is not None and s.dna_size > max_len
            ):
                continue
            if max_n is not None and self._n_fraction(s) > max_n:
                continue
            out.append(s)
        return out

    def sort(self, doc: list[TwoBitSeq], opts: Opts) -> list[TwoBitSeq]:
        key = opts.text("key", "name", ["name", "natural", "length"])
        reverse = opts.boolean("reverse", False)
        keys = {
            "name": lambda s: s.name,
            "natural": lambda s: (chrom_key(s.name), s.name),
            "length": lambda s: (s.dna_size, s.name),
        }
        return sorted(doc, key=keys[str(key)], reverse=reverse)

    def merge(self, docs: list[list[TwoBitSeq]], opts: Opts) -> list[TwoBitSeq]:
        policy = opts.text("duplicates", "error", ["error", "first"])
        seen: set[str] = set()
        out = []
        for doc in docs:
            for s in doc:
                if s.name in seen:
                    if policy == "first":
                        continue
                    raise bad_input(
                        f"Duplicate sequence name {s.name!r} while merging 2bit files.",
                        code="IO_FORMATS_DUPLICATE_ID",
                        hint="Pass duplicates='first' to keep the first one.",
                    )
                seen.add(s.name)
                out.append(s)
        return out

    def split(self, doc: list[TwoBitSeq], opts: Opts) -> list[tuple[str, list[TwoBitSeq]]]:
        size = opts.integer("sequences_per_part", 1, 1) or 1
        return [
            (doc[i].name if size == 1 else f"{i // size + 1:04d}", doc[i : i + size])
            for i in range(0, len(doc), size)
        ]

    def convert(self, doc: list[TwoBitSeq], target: str, opts: Opts) -> Payload:
        soft = opts.boolean("soft_mask", True)
        width = opts.integer("line_width", 60, 0) or 0
        names = opts.names("names")
        chosen = [s for s in doc if names is None or s.name in names]
        if names is not None and len(chosen) != len(set(names)):
            raise bad_option(
                f"Some of the requested names are not in the 2bit file: {sorted(set(names) - {s.name for s in chosen})[:5]}"
            )
        _bounded_total([s.dna_size for s in chosen], opts)
        if target == "twobit":
            return self.dump(chosen, opts)
        if target == "fasta":
            lines = []
            for s in chosen:
                text = self.decode_seq(s, soft)
                lines.append(f">{s.name}")
                lines.extend(
                    [text[i : i + width] for i in range(0, len(text), width)] if width else [text]
                )
            return Payload("fasta", text="".join(line + "\n" for line in lines), count=len(chosen))
        return json_payload(
            "twobit",
            [
                {
                    "name": s.name,
                    "length": s.dna_size,
                    "n_blocks": s.n_blocks,
                    "mask_blocks": s.mask_blocks,
                    "sequence": self.decode_seq(s, soft),
                }
                for s in chosen
            ],
            len(chosen),
        )
