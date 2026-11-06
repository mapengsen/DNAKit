"""Local file and metadata candidate implementations."""

from __future__ import annotations

import gzip
import hashlib
import json
import mimetypes
from collections.abc import Mapping, Sequence
from pathlib import Path


def _bytes(value: object) -> bytes:
    if isinstance(value, bytes):
        return value
    if isinstance(value, Path) or (
        isinstance(value, str) and "\n" not in value and Path(value).exists()
    ):
        return Path(value).read_bytes()
    if isinstance(value, str):
        return value.encode("utf-8")
    return json.dumps(value, ensure_ascii=False, sort_keys=True).encode("utf-8")


def _text(value: object) -> str:
    return _bytes(value).decode("utf-8")


def build_file_manifest(input: object) -> dict[str, object]:
    """Return size, checksum and media metadata for a file or virtual payload."""

    payload = _bytes(input)
    path = str(input) if isinstance(input, (str, Path)) and Path(input).exists() else None
    return {
        "path": path,
        "size_bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "media_type": mimetypes.guess_type(path or "")[0],
    }


def chunk_large_file(input: object, *, chunk_size: int = 1024) -> list[bytes]:
    """Split a payload into bounded byte chunks."""

    if chunk_size < 1:
        raise ValueError("chunk_size must be positive")
    payload = _bytes(input)
    return [payload[index : index + chunk_size] for index in range(0, len(payload), chunk_size)]


def clean_temporary_files(
    input: str | Path, *, suffixes: Sequence[str] = (".tmp", ".part")
) -> list[str]:
    """Delete matching temporary files in a directory and return deleted paths."""

    directory = Path(input)
    if not directory.exists() or not directory.is_dir():
        return []
    deleted = []
    for path in directory.iterdir():
        if path.is_file() and path.suffix in suffixes:
            path.unlink()
            deleted.append(str(path))
    return deleted


def compress(input: object, *, level: int = 6) -> bytes:
    """Compress UTF-8 or binary input with gzip."""

    return gzip.compress(_bytes(input), compresslevel=level)


def decompress(input: bytes) -> bytes:
    """Decompress gzip bytes."""

    return gzip.decompress(input)


def detect_compression(input: object) -> str:
    """Detect gzip, zip, bzip2, xz or uncompressed bytes."""

    payload = _bytes(input)
    signatures = (
        (b"\x1f\x8b", "gzip"),
        (b"PK\x03\x04", "zip"),
        (b"BZh", "bzip2"),
        (b"\xfd7zXZ", "xz"),
    )
    return next((name for signature, name in signatures if payload.startswith(signature)), "none")


def detect_file_format(input: object) -> str:
    """Detect common sequence and tabular formats from content or suffix."""

    if isinstance(input, (str, Path)) and Path(input).exists():
        suffix = Path(input).suffix.lower().lstrip(".")
        if suffix:
            return suffix
    text = _text(input).lstrip()
    if text.startswith(">"):
        return "fasta"
    if text.startswith("@") and text.count("\n") >= 3:
        return "fastq"
    if text.startswith("#NEXUS"):
        return "nexus"
    if "\t" in text:
        return "tabular"
    return "text"


def estimate_storage_size(input: object, *, compression: str = "none") -> int:
    """Estimate encoded size using the selected compression mode."""

    payload = _bytes(input)
    if compression == "gzip":
        return len(compress(payload))
    if compression != "none":
        raise ValueError("compression must be none or gzip")
    return len(payload)


def get_file_checksum(input: object, *, algorithm: str = "sha256") -> str:
    """Return a hexadecimal digest for a payload."""

    try:
        digest = hashlib.new(algorithm)
    except ValueError as exc:
        raise ValueError(f"unsupported checksum algorithm: {algorithm}") from exc
    digest.update(_bytes(input))
    return digest.hexdigest()


def list_file_records(input: object) -> list[str]:
    """Return non-empty logical lines from a text payload."""

    return [line for line in _text(input).splitlines() if line.strip()]


def read(input: object, *, encoding: str = "utf-8") -> str:
    """Read a text payload or file."""

    if isinstance(input, (str, Path)) and Path(input).exists():
        return Path(input).read_text(encoding=encoding)
    if isinstance(input, bytes):
        return input.decode(encoding)
    return str(input)


def stream_records(input: object) -> list[str]:
    """Materialize a bounded line stream for MCP-safe serialization."""

    return list_file_records(input)


def validate(input: object) -> dict[str, object]:
    """Validate that a payload is decodable UTF-8 and has no NUL bytes."""

    payload = _bytes(input)
    try:
        text = payload.decode("utf-8")
        valid = b"\x00" not in payload
        error = None if valid else "NUL byte is not allowed"
    except UnicodeDecodeError as exc:
        text, valid, error = "", False, str(exc)
    return {
        "valid": valid,
        "size_bytes": len(payload),
        "line_count": len(text.splitlines()),
        "error": error,
    }


def validate_metadata_schema(
    input: Mapping[str, object], *, required: Sequence[str] = ()
) -> dict[str, object]:
    """Validate required metadata keys and return a structured report."""

    missing = [key for key in required if key not in input]
    return {"valid": not missing, "missing": missing, "keys": sorted(input)}


def verify_checksum(input: object, *, expected: str, algorithm: str = "sha256") -> bool:
    """Compare a payload digest with an expected hexadecimal checksum."""

    return get_file_checksum(input, algorithm=algorithm).casefold() == expected.casefold()


def repack_compressed_file(input: bytes, *, level: int = 6) -> bytes:
    """Decompress and recompress gzip content using a selected level."""

    return compress(decompress(input), level=level)


def encode_metadata_categorical(input: Mapping[str, Sequence[object]]) -> dict[str, object]:
    """Encode each categorical metadata column as integer codes."""

    encoded = {}
    categories = {}
    for key, values in input.items():
        labels = sorted({str(value) for value in values})
        lookup = {label: index for index, label in enumerate(labels)}
        categories[key] = labels
        encoded[key] = [lookup[str(value)] for value in values]
    return {"data": encoded, "categories": categories}


def decode_metadata_categorical(input: Mapping[str, object]) -> dict[str, list[str]]:
    """Decode a result returned by ``encode_metadata_categorical``."""

    categories = input["categories"]
    data = input["data"]
    return {
        key: [str(categories[key][int(index)]) for index in values] for key, values in data.items()
    }  # type: ignore[index]


def normalize_sample_metadata(input: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    """Normalize sample IDs and trim string metadata values."""

    result = []
    for index, row in enumerate(input, start=1):
        normalized = {
            str(key).strip(): value.strip() if isinstance(value, str) else value
            for key, value in row.items()
        }
        normalized.setdefault("sample_id", f"sample_{index}")
        result.append(normalized)
    return result


def merge_metadata_tables(
    input: Sequence[Sequence[Mapping[str, object]]],
) -> list[dict[str, object]]:
    """Concatenate metadata tables while preserving row order."""

    return [dict(row) for table in input for row in table]


def write(input: object, *, target: str | Path | None = None, encoding: str = "utf-8") -> str:
    """Serialize text and optionally write it to a local target."""

    text = _text(input)
    if target is not None:
        path = Path(target)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding=encoding)
    return text


__all__ = [
    "build_file_manifest",
    "chunk_large_file",
    "clean_temporary_files",
    "compress",
    "decode_metadata_categorical",
    "decompress",
    "detect_compression",
    "detect_file_format",
    "encode_metadata_categorical",
    "estimate_storage_size",
    "get_file_checksum",
    "list_file_records",
    "merge_metadata_tables",
    "normalize_sample_metadata",
    "read",
    "repack_compressed_file",
    "stream_records",
    "validate",
    "validate_metadata_schema",
    "verify_checksum",
    "write",
]


from dnakit.io_formats._impl_compress import run as _run_compress
from dnakit.io_formats._impl_dispatch import run as _run_io_format


def compress_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Gzip-compress a PLINK .map variant file and return it as base64 text; the content is not parsed.

    The gzip header is fixed (modification time 0), so the same input and level always give the same output;
    decompress_map restores the exact bytes.

    input: path (str or Path) of a file (its bytes are compressed), text (str, encoded as UTF-8) or bytes (binary data
    such as a .2bit file is kept bit-exact). A directory is refused.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs larger than this
        level: int 1-9, default 6: gzip compression level
        output: path of a .gz file to write (optional; nothing is written otherwise)
        overwrite: bool, default False: allow replacing an existing output file
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'encoding': 'gzip+base64', 'data', 'original_bytes',
    'compressed_bytes', 'sha256_original'} plus 'output_path' when options['output'] is given.
    Raises InputFormatError for unusable inputs and ConfigurationError for invalid options, existing outputs or inputs above
    the size limit.
    """

    return _run_compress("compress_map", input, options)


def compress_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Gzip-compress a Newick tree file and return it as base64 text; the content is not parsed.

    The gzip header is fixed (modification time 0), so the same input and level always give the same output;
    decompress_newick restores the exact bytes.

    input: path (str or Path) of a file (its bytes are compressed), text (str, encoded as UTF-8) or bytes (binary data
    such as a .2bit file is kept bit-exact). A directory is refused.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs larger than this
        level: int 1-9, default 6: gzip compression level
        output: path of a .gz file to write (optional; nothing is written otherwise)
        overwrite: bool, default False: allow replacing an existing output file
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'encoding': 'gzip+base64', 'data', 'original_bytes',
    'compressed_bytes', 'sha256_original'} plus 'output_path' when options['output'] is given.
    Raises InputFormatError for unusable inputs and ConfigurationError for invalid options, existing outputs or inputs above
    the size limit.
    """

    return _run_compress("compress_newick", input, options)


def compress_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Gzip-compress a NEXUS file and return it as base64 text; the content is not parsed.

    The gzip header is fixed (modification time 0), so the same input and level always give the same output;
    decompress_nexus restores the exact bytes.

    input: path (str or Path) of a file (its bytes are compressed), text (str, encoded as UTF-8) or bytes (binary data
    such as a .2bit file is kept bit-exact). A directory is refused.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs larger than this
        level: int 1-9, default 6: gzip compression level
        output: path of a .gz file to write (optional; nothing is written otherwise)
        overwrite: bool, default False: allow replacing an existing output file
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'encoding': 'gzip+base64', 'data', 'original_bytes',
    'compressed_bytes', 'sha256_original'} plus 'output_path' when options['output'] is given.
    Raises InputFormatError for unusable inputs and ConfigurationError for invalid options, existing outputs or inputs above
    the size limit.
    """

    return _run_compress("compress_nexus", input, options)


def compress_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Gzip-compress a PAF alignment file and return it as base64 text; the content is not parsed.

    The gzip header is fixed (modification time 0), so the same input and level always give the same output;
    decompress_paf restores the exact bytes.

    input: path (str or Path) of a file (its bytes are compressed), text (str, encoded as UTF-8) or bytes (binary data
    such as a .2bit file is kept bit-exact). A directory is refused.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs larger than this
        level: int 1-9, default 6: gzip compression level
        output: path of a .gz file to write (optional; nothing is written otherwise)
        overwrite: bool, default False: allow replacing an existing output file
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'encoding': 'gzip+base64', 'data', 'original_bytes',
    'compressed_bytes', 'sha256_original'} plus 'output_path' when options['output'] is given.
    Raises InputFormatError for unusable inputs and ConfigurationError for invalid options, existing outputs or inputs above
    the size limit.
    """

    return _run_compress("compress_paf", input, options)


def compress_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Gzip-compress a PLINK .ped pedigree/genotype file and return it as base64 text; the content is not parsed.

    The gzip header is fixed (modification time 0), so the same input and level always give the same output;
    decompress_ped restores the exact bytes.

    input: path (str or Path) of a file (its bytes are compressed), text (str, encoded as UTF-8) or bytes (binary data
    such as a .2bit file is kept bit-exact). A directory is refused.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs larger than this
        level: int 1-9, default 6: gzip compression level
        output: path of a .gz file to write (optional; nothing is written otherwise)
        overwrite: bool, default False: allow replacing an existing output file
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'encoding': 'gzip+base64', 'data', 'original_bytes',
    'compressed_bytes', 'sha256_original'} plus 'output_path' when options['output'] is given.
    Raises InputFormatError for unusable inputs and ConfigurationError for invalid options, existing outputs or inputs above
    the size limit.
    """

    return _run_compress("compress_ped", input, options)


def compress_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Gzip-compress a PHYLIP alignment file and return it as base64 text; the content is not parsed.

    The gzip header is fixed (modification time 0), so the same input and level always give the same output;
    decompress_phylip restores the exact bytes.

    input: path (str or Path) of a file (its bytes are compressed), text (str, encoded as UTF-8) or bytes (binary data
    such as a .2bit file is kept bit-exact). A directory is refused.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs larger than this
        level: int 1-9, default 6: gzip compression level
        output: path of a .gz file to write (optional; nothing is written otherwise)
        overwrite: bool, default False: allow replacing an existing output file
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'encoding': 'gzip+base64', 'data', 'original_bytes',
    'compressed_bytes', 'sha256_original'} plus 'output_path' when options['output'] is given.
    Raises InputFormatError for unusable inputs and ConfigurationError for invalid options, existing outputs or inputs above
    the size limit.
    """

    return _run_compress("compress_phylip", input, options)


def compress_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Gzip-compress a UCSC 2bit genome file and return it as base64 text; the content is not parsed.

    The gzip header is fixed (modification time 0), so the same input and level always give the same output;
    decompress_twobit restores the exact bytes.

    input: path (str or Path) of a file (its bytes are compressed), text (str, encoded as UTF-8) or bytes (binary data
    such as a .2bit file is kept bit-exact). A directory is refused.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs larger than this
        level: int 1-9, default 6: gzip compression level
        output: path of a .gz file to write (optional; nothing is written otherwise)
        overwrite: bool, default False: allow replacing an existing output file
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'encoding': 'gzip+base64', 'data', 'original_bytes',
    'compressed_bytes', 'sha256_original'} plus 'output_path' when options['output'] is given.
    Raises InputFormatError for unusable inputs and ConfigurationError for invalid options, existing outputs or inputs above
    the size limit.
    """

    return _run_compress("compress_twobit", input, options)


def compress_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Pack a Zarr store directory as a reproducible tar archive and gzip-compress it, returned as base64 text.

    The gzip header is fixed (modification time 0), so the same input and level always give the same output;
    decompress_zarr restores the exact bytes.

    input: the path (str or Path) of the Zarr directory (archived as tar with sorted entries, no timestamps or owners, then
    gzip-compressed), the path of a file, or text (str, UTF-8) / bytes. Symbolic links inside the directory are refused.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs larger than this
        level: int 1-9, default 6: gzip compression level
        output: path of a .gz file to write (optional; nothing is written otherwise)
        overwrite: bool, default False: allow replacing an existing output file
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'encoding': 'gzip+base64', 'data', 'original_bytes',
    'compressed_bytes', 'sha256_original'} plus 'output_path' when options['output'] is given and 'archive': 'tar' for a directory.
    Raises InputFormatError for unusable inputs and ConfigurationError for invalid options, existing outputs or inputs above
    the size limit.
    """

    return _run_compress("compress_zarr", input, options)


def convert(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Convert a file of any supported format (map, ped, paf, newick, nexus, phylip, twobit, zarr; FASTA as a source) to another format.

    The format is taken from options['format'] (one of map, newick, nexus, paf, ped, phylip, twobit, zarr) or detected
    from the file suffix and content, and the call is executed by the same engine as the format-specific function (e.g.
    filter_paf); its options and result are those of that function, plus the options below.

    input: as for the format-specific function (a path, text/bytes, or - for merge - a list of them).
    options:
        max_bytes: int, default 268435456
        to: target format name (required); see the targets of the format-specific convert function
        output: path to write (optional)
        overwrite: bool, default False
        format: name of the input format; detected from the suffix/content when omitted
    Returns the DNAKit envelope (status 'ok') whose result is that of the format-specific function. Raises
    InputFormatError for malformed content and ConfigurationError for invalid options.
    """

    return _run_io_format("convert", input, options)


def convert_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Convert a PLINK .map variant file to another format (targets: map, bed, tsv, json).

    Targets: 'map' (normalised, tab-separated), 'bed' (chrom, bp-1, bp, ID; BED is 0-based half-open), 'tsv' (with header
    chrom, id, cm, pos) and 'json'.

    input: path (str or Path) of a PLINK .map variant file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        to: target format name (required): map, bed, tsv, json
        output: path to write (optional)
        overwrite: bool, default False
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("convert_map", input, options)


def convert_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Convert a Newick tree file to another format (targets: newick, nexus, tsv, json).

    Targets: 'newick' (normalised), 'nexus' (TREES block, tree names tree1..), 'tsv' (node table: tree, node_id,
    parent_id, name, length, is_tip) and 'json' (records with flat node tables).

    input: path (str or Path) of a Newick tree file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        to: target format name (required): newick, nexus, tsv, json
        output: path to write (optional)
        overwrite: bool, default False
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("convert_newick", input, options)


def convert_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Convert a NEXUS file to another format (targets: nexus, phylip, fasta, twobit, newick, json).

    Targets: 'nexus' (normalised: DATA block + TREES block), 'phylip', 'fasta', 'twobit' (nucleotides), 'newick' (one tree
    per line; needs trees) and 'json'. The phylip, fasta and twobit targets hold only the matrix: a TREES block is
    dropped and reported in 'warnings'.

    input: path (str or Path) of a NEXUS file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        to: target format name (required): nexus, phylip, fasta, twobit, newick, json
        output: path to write (optional)
        overwrite: bool, default False
        layout: 'sequential' (default) or 'interleaved' for PHYLIP/NEXUS output
        block_size: int, default 60: sites per block of an interleaved output
        strict_names: bool, default False: write strict 10-character PHYLIP names (error when longer)
        line_width: int, default 60: FASTA line width (0 = one line per sequence)
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("convert_nexus", input, options)


def convert_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Convert a PAF (Pairwise mApping Format) alignment file to another format (targets: paf, bed, tsv, json).

    Targets: 'paf' (normalised), 'bed' (target interval: target, tstart, tend, query, mapq, strand), 'tsv' (header plus
    the 12 columns and one column per tag) and 'json'.

    input: path (str or Path) of a PAF (Pairwise mApping Format) alignment file (optionally gzip-compressed), or its
    content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        to: target format name (required): paf, bed, tsv, json
        output: path to write (optional)
        overwrite: bool, default False
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("convert_paf", input, options)


def convert_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Convert a PLINK .ped pedigree/genotype file to another format (targets: ped, tsv, json).

    Targets: 'ped' (normalised), 'tsv' (header FID IID PAT MAT SEX PHENO V1..Vk, genotypes as 'a1/a2') and 'json'.

    input: path (str or Path) of a PLINK .ped pedigree/genotype file (optionally gzip-compressed), or its content as
    str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        to: target format name (required): ped, tsv, json
        output: path to write (optional)
        overwrite: bool, default False
        max_cells: int, default 20000000
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("convert_ped", input, options)


def convert_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Convert a PHYLIP alignment file to another format (targets: phylip, nexus, fasta, twobit, json).

    Targets: 'phylip' (re-layout), 'nexus', 'fasta', 'twobit' (nucleotides) and 'json'.

    input: path (str or Path) of a PHYLIP alignment file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        to: target format name (required): phylip, nexus, fasta, twobit, json
        output: path to write (optional)
        overwrite: bool, default False
        datatype: 'dna' (default), 'rna', 'protein' or 'any': alphabet used for validation and for counting missing
            characters (gap, ?, N for nucleotides, X for protein)
        input_layout: 'auto' (default), 'sequential' or 'interleaved': layout of the PHYLIP input
        input_strict_names: bool, default False: input names occupy the first 10 characters of each row (otherwise the
            first whitespace-delimited word)
        layout: 'sequential' (default) or 'interleaved' for PHYLIP/NEXUS output
        block_size: int, default 60: sites per block of an interleaved output
        strict_names: bool, default False: write strict 10-character PHYLIP names (error when longer)
        line_width: int, default 60: FASTA line width (0 = one line per sequence)
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("convert_phylip", input, options)


def convert_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Convert a UCSC 2bit genome file to another format (targets: twobit, fasta, json).

    Targets: 'fasta' (soft-masking kept as lower case), 'json' (metadata and sequences) and 'twobit' (re-encode, e.g. to
    change byte order). FASTA -> 2bit is convert(format='fasta', to='twobit').

    input: path (str or Path) of a .2bit file, or its raw bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        to: target format name (required): twobit, fasta, json
        output: path to write (optional)
        overwrite: bool, default False
        line_width: int, default 60: FASTA line width (0 = one line)
        soft_mask: bool, default True
        max_bases: int, default 50000000
        names: list of sequence names to convert (default all)
        byte_order: 'little' (default) or 'big': byte order of the written 2bit file
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("convert_twobit", input, options)


def convert_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Convert a Zarr directory store to another format (targets: zarr, json, tsv, npz).

    Targets: 'json' (metadata and values), 'tsv' (one 1-D/2-D array, tab separated), 'npz' (numpy archive of all arrays)
    and 'zarr' (a copy of the store at output, which is required).

    input: path (str or Path) of a zarr directory store (requires the optional 'zarr>=3' and 'numpy' packages).
    options (all optional unless stated):
        to: target format name (required): zarr, json, tsv, npz
        output: path to write (required for the 'zarr' target, optional otherwise)
        overwrite: bool, default False: replace an existing zarr store at output
        array: array path for the 'tsv' target (1-D or 2-D)
        max_elements: int, default 5000000: refuse to materialise more array elements than this
    Returns the DNAKit envelope (status 'ok'); result: for target 'zarr' {'format': 'zarr', 'record_count', 'output_path'}
    (the new store is written to options['output']); for 'json'/'tsv'/'npz' {'format', 'record_count', 'bytes', 'sha256'}
    plus either the content ('text', or 'data_base64' for binary formats) or, when options['output'] is given,
    'output_path' (the file is written; an existing file needs overwrite=True).
    Returns status 'not_implemented' (requires="pip install 'zarr>=3' numpy") when zarr or numpy is not installed. Raises
    InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size limits.
    """

    return _run_io_format("convert_zarr", input, options)


def decompress_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Decompress gzip data into the original bytes of a PLINK .map variant file.

    input: base64 text of a gzip stream (the 'data' of compress_map), gzip bytes, or the path (str or Path) of a .gz file.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs and decompressed data larger than this
        output: path of a file to write the decompressed bytes to (optional; they are returned inline otherwise)
        overwrite: bool, default False: allow replacing an existing output file
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'original_bytes', 'compressed_bytes', 'sha256'} plus the
    content as 'text' when the bytes are valid UTF-8, otherwise as 'data_base64', or 'output_path' when
    options['output'] is given.
    Raises InputFormatError when the input is not valid base64 or gzip and ConfigurationError for invalid options,
    existing outputs or data above the size limit.
    """

    return _run_compress("decompress_map", input, options)


def decompress_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Decompress gzip data into the original bytes of a Newick tree file.

    input: base64 text of a gzip stream (the 'data' of compress_newick), gzip bytes, or the path (str or Path) of a .gz file.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs and decompressed data larger than this
        output: path of a file to write the decompressed bytes to (optional; they are returned inline otherwise)
        overwrite: bool, default False: allow replacing an existing output file
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'original_bytes', 'compressed_bytes', 'sha256'} plus the
    content as 'text' when the bytes are valid UTF-8, otherwise as 'data_base64', or 'output_path' when
    options['output'] is given.
    Raises InputFormatError when the input is not valid base64 or gzip and ConfigurationError for invalid options,
    existing outputs or data above the size limit.
    """

    return _run_compress("decompress_newick", input, options)


def decompress_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Decompress gzip data into the original bytes of a NEXUS file.

    input: base64 text of a gzip stream (the 'data' of compress_nexus), gzip bytes, or the path (str or Path) of a .gz file.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs and decompressed data larger than this
        output: path of a file to write the decompressed bytes to (optional; they are returned inline otherwise)
        overwrite: bool, default False: allow replacing an existing output file
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'original_bytes', 'compressed_bytes', 'sha256'} plus the
    content as 'text' when the bytes are valid UTF-8, otherwise as 'data_base64', or 'output_path' when
    options['output'] is given.
    Raises InputFormatError when the input is not valid base64 or gzip and ConfigurationError for invalid options,
    existing outputs or data above the size limit.
    """

    return _run_compress("decompress_nexus", input, options)


def decompress_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Decompress gzip data into the original bytes of a PAF alignment file.

    input: base64 text of a gzip stream (the 'data' of compress_paf), gzip bytes, or the path (str or Path) of a .gz file.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs and decompressed data larger than this
        output: path of a file to write the decompressed bytes to (optional; they are returned inline otherwise)
        overwrite: bool, default False: allow replacing an existing output file
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'original_bytes', 'compressed_bytes', 'sha256'} plus the
    content as 'text' when the bytes are valid UTF-8, otherwise as 'data_base64', or 'output_path' when
    options['output'] is given.
    Raises InputFormatError when the input is not valid base64 or gzip and ConfigurationError for invalid options,
    existing outputs or data above the size limit.
    """

    return _run_compress("decompress_paf", input, options)


def decompress_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Decompress gzip data into the original bytes of a PLINK .ped pedigree/genotype file.

    input: base64 text of a gzip stream (the 'data' of compress_ped), gzip bytes, or the path (str or Path) of a .gz file.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs and decompressed data larger than this
        output: path of a file to write the decompressed bytes to (optional; they are returned inline otherwise)
        overwrite: bool, default False: allow replacing an existing output file
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'original_bytes', 'compressed_bytes', 'sha256'} plus the
    content as 'text' when the bytes are valid UTF-8, otherwise as 'data_base64', or 'output_path' when
    options['output'] is given.
    Raises InputFormatError when the input is not valid base64 or gzip and ConfigurationError for invalid options,
    existing outputs or data above the size limit.
    """

    return _run_compress("decompress_ped", input, options)


def decompress_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Decompress gzip data into the original bytes of a PHYLIP alignment file.

    input: base64 text of a gzip stream (the 'data' of compress_phylip), gzip bytes, or the path (str or Path) of a .gz file.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs and decompressed data larger than this
        output: path of a file to write the decompressed bytes to (optional; they are returned inline otherwise)
        overwrite: bool, default False: allow replacing an existing output file
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'original_bytes', 'compressed_bytes', 'sha256'} plus the
    content as 'text' when the bytes are valid UTF-8, otherwise as 'data_base64', or 'output_path' when
    options['output'] is given.
    Raises InputFormatError when the input is not valid base64 or gzip and ConfigurationError for invalid options,
    existing outputs or data above the size limit.
    """

    return _run_compress("decompress_phylip", input, options)


def decompress_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Decompress gzip data into the original bytes of a UCSC 2bit genome file.

    input: base64 text of a gzip stream (the 'data' of compress_twobit), gzip bytes, or the path (str or Path) of a .gz file.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs and decompressed data larger than this
        output: path of a file to write the decompressed bytes to (optional; they are returned inline otherwise)
        overwrite: bool, default False: allow replacing an existing output file
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'original_bytes', 'compressed_bytes', 'sha256'} plus the
    content as 'text' when the bytes are valid UTF-8, otherwise as 'data_base64', or 'output_path' when
    options['output'] is given.
    Raises InputFormatError when the input is not valid base64 or gzip and ConfigurationError for invalid options,
    existing outputs or data above the size limit.
    """

    return _run_compress("decompress_twobit", input, options)


def decompress_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Decompress gzip data made by compress_zarr: restore a Zarr store directory, or the original bytes.

    input: base64 text of a gzip stream (the 'data' of compress_zarr), gzip bytes, or the path (str or Path) of a .gz file.
    When the decompressed data is the tar archive of a directory (what compress_zarr makes from a directory path) it is
    extracted to options['output'], or returned as base64 without it; any other data is handled as by the other
    decompress_* functions (text or base64, or written to the file options['output']).
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs and decompressed data larger than this
        output: a directory to extract an archive into (it must be empty or new unless overwrite=True), or the path of a file
            for any other data (optional; without it nothing is written)
        overwrite: bool, default False: allow a non-empty directory (existing files are replaced) or an existing output file
    Extraction is safe: absolute paths, '..' components, links and special files in the archive are refused before
    anything is written.
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'original_bytes', 'compressed_bytes', 'sha256'} plus, for
    an archive, 'archive': 'tar' and either 'output_path' with 'file_count' or 'data_base64' (the tar archive); for other
    data either 'text' (valid UTF-8), 'data_base64', or 'output_path'.
    Raises InputFormatError when the input is not valid base64 or gzip or the archive has unsafe members and
    ConfigurationError for invalid options, non-empty outputs or data above the size limit.
    """

    return _run_compress("decompress_zarr", input, options)


def filter(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Filter a file of any supported format with the criteria of its format-specific filter function.

    The format is taken from options['format'] (one of map, newick, nexus, paf, ped, phylip, twobit, zarr) or detected
    from the file suffix and content, and the call is executed by the same engine as the format-specific function (e.g.
    filter_paf); its options and result are those of that function, plus the options below.

    input: as for the format-specific function (a path, text/bytes, or - for merge - a list of them).
    options:
        max_bytes: int, default 268435456
        output: path to write (optional)
        overwrite: bool, default False
        format: name of the input format; detected from the suffix/content when omitted
    Returns the DNAKit envelope (status 'ok') whose result is that of the format-specific function. Raises
    InputFormatError for malformed content and ConfigurationError for invalid options.
    """

    return _run_io_format("filter", input, options)


def filter_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Filter a PLINK .map variant file by chromosome, variant ID, position and genetic-distance criteria and return the filtered .map content.

    input: path (str or Path) of a PLINK .map variant file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output: path to write (optional)
        overwrite: bool, default False
        chromosomes: list of chromosome names to keep
        exclude_chromosomes: list of chromosome names to drop
        ids: list of variant IDs to keep
        exclude_ids: list of variant IDs to drop
        region: 'chrom:start-end' keep variants with start <= bp <= end (1-based, inclusive)
        min_position: int: keep bp >= value
        max_position: int: keep bp <= value
        min_cm: float: keep genetic distance >= value
        max_cm: float: keep genetic distance <= value
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("filter_map", input, options)


def filter_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Filter a Newick tree file by tip names, short branches and tip counts and return the filtered .nwk content.

    input: path (str or Path) of a Newick tree file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output: path to write (optional)
        overwrite: bool, default False
        keep_tips: list of tip names to keep (other tips are pruned; unary nodes are merged, summing branch lengths)
        drop_tips: list of tip names to prune
        collapse_below: float: contract internal branches shorter than this into polytomies
        min_tips: int: drop trees with fewer tips (evaluated after pruning)
        max_tips: int: drop trees with more tips (evaluated after pruning)
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("filter_newick", input, options)


def filter_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Filter a NEXUS file by taxon names and per-taxon / per-site missing data (trees are pruned to the kept taxa) and return the filtered .nex content.

    Taxon filters apply to the matrix, the TAXA labels and (by pruning tips) every tree.

    input: path (str or Path) of a NEXUS file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output: path to write (optional)
        overwrite: bool, default False
        keep: list of taxon names to keep
        exclude: list of taxon names to drop
        name_pattern: regular expression (re.search, max 256 characters) a taxon name must match
        max_missing_fraction: float in [0,1]: drop taxa whose gap/missing/ambiguous fraction exceeds it
        max_site_missing_fraction: float in [0,1]: drop alignment columns whose missing fraction (over the kept taxa)
            exceeds it
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("filter_nexus", input, options)


def filter_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Filter a PAF (Pairwise mApping Format) alignment file by mapping quality, identity, length, strand, name, region and tag criteria and return the filtered .paf content.

    input: path (str or Path) of a PAF (Pairwise mApping Format) alignment file (optionally gzip-compressed), or its
    content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output: path to write (optional)
        overwrite: bool, default False
        min_mapq: int: keep mapping quality >= value
        min_identity: float: keep matches/block_length >= value
        min_block_length: int: keep alignment block length >= value
        min_matches: int: keep matching bases >= value
        min_query_coverage: float: keep (query_end-query_start)/query_length >= value
        strand: '+' or '-'
        query_names: list of query names to keep
        exclude_query_names: list to drop
        target_names: list of target names to keep
        exclude_target_names: list to drop
        target_region: 'name:start-end' keep alignments overlapping the 0-based half-open target interval
        query_region: 'name:start-end' keep alignments overlapping the 0-based half-open query interval
        tags: mapping {tag: value}: keep alignments carrying every tag with that value (e.g. {'tp': 'P'})
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("filter_paf", input, options)


def filter_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Filter a PLINK .ped pedigree/genotype file by family, individual, sex, phenotype and missing-genotype criteria and return the filtered .ped content.

    Individual filters are combined with AND; variant filters drop genotype columns from every kept row.

    input: path (str or Path) of a PLINK .ped pedigree/genotype file (optionally gzip-compressed), or its content as
    str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output: path to write (optional)
        overwrite: bool, default False
        families: list of family IDs to keep
        exclude_families: list of family IDs to drop
        individuals: list of individual IDs to keep
        exclude_individuals: list of individual IDs to drop
        sex: keep one sex: 1/'male', 2/'female' or 0/'unknown'
        phenotype: keep one class: 'case' (2), 'control' (1), 'missing' (0 or -9), 'observed' (not missing) or
            'quantitative' (any other value)
        founders_only: bool, default False: keep individuals whose father and mother IDs are both 0
        max_missing_rate: float in [0,1]: drop individuals whose fraction of missing genotypes exceeds it
        variant_indices: list of 0-based variant columns to keep (apply the same to the .map file)
        max_variant_missing_rate: float in [0,1]: drop variant columns whose missing-call fraction exceeds it (computed
            over the individuals that pass the other filters)
        max_cells: int, default 20000000
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("filter_ped", input, options)


def filter_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Filter a PHYLIP alignment file by taxon names and per-taxon / per-site missing data and return the filtered .phy content.

    input: path (str or Path) of a PHYLIP alignment file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output: path to write (optional)
        overwrite: bool, default False
        datatype: 'dna' (default), 'rna', 'protein' or 'any': alphabet used for validation and for counting missing
            characters (gap, ?, N for nucleotides, X for protein)
        input_layout: 'auto' (default), 'sequential' or 'interleaved': layout of the PHYLIP input
        input_strict_names: bool, default False: input names occupy the first 10 characters of each row (otherwise the
            first whitespace-delimited word)
        keep: list of taxon names to keep
        exclude: list of taxon names to drop
        name_pattern: regular expression (re.search, max 256 characters) a taxon name must match
        max_missing_fraction: float in [0,1]: drop taxa whose gap/missing/ambiguous fraction exceeds it
        max_site_missing_fraction: float in [0,1]: drop alignment columns whose missing fraction (over the kept taxa)
            exceeds it
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("filter_phylip", input, options)


def filter_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Filter a UCSC 2bit genome file by sequence name, length and N fraction and return the filtered .2bit content.

    Keeps whole sequences (records are copied bit-exactly, no re-encoding) that pass all the given criteria.

    input: path (str or Path) of a .2bit file, or its raw bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output: path to write (optional)
        overwrite: bool, default False
        names: list of sequence names to keep
        exclude: list of names to drop
        name_pattern: regular expression (re.search, max 256 characters) the name must match
        min_length: int: keep sequences with at least this many bases
        max_length: int: at most this many bases
        max_n_fraction: float in [0,1]: drop sequences whose N fraction exceeds it
        byte_order: 'little' (default) or 'big': byte order of the written 2bit file
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("filter_twobit", input, options)


def filter_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Filter a Zarr directory store by array paths and rows along the first axis and return the filtered .zarr content.

    Copies the selected arrays (by path / group prefix) and optionally only the chosen rows into a new store.

    input: path (str or Path) of a zarr directory store (requires the optional 'zarr>=3' and 'numpy' packages).
    options (all optional unless stated):
        output: directory path of the new zarr store (required)
        overwrite: bool, default False: replace an existing zarr store at output
        arrays: array paths (or group prefixes) to keep
        exclude: array paths (or group prefixes) to drop
        rows: rows to keep along axis 0: 'start:stop[:step]', {start, stop, step} or a list of indices (all selected
            arrays with at least one dimension must have the same length on axis 0)
        max_elements: int, default 5000000: refuse to materialise more array elements than this
    Returns the DNAKit envelope (status 'ok'); result: {'format': 'zarr', 'record_count', 'output_path'} (the new store is
    written to options['output']).
    Returns status 'not_implemented' (requires="pip install 'zarr>=3' numpy") when zarr or numpy is not installed. Raises
    InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size limits.
    """

    return _run_io_format("filter_zarr", input, options)


def index(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build a lookup index of a file of any supported format.

    The format is taken from options['format'] (one of map, newick, nexus, paf, ped, phylip, twobit, zarr) or detected
    from the file suffix and content, and the call is executed by the same engine as the format-specific function (e.g.
    filter_paf); its options and result are those of that function, plus the options below.

    input: as for the format-specific function (a path, text/bytes, or - for merge - a list of them).
    options:
        max_bytes: int, default 268435456
        format: name of the input format; detected from the suffix/content when omitted
    Returns the DNAKit envelope (status 'ok') whose result is that of the format-specific function. Raises
    InputFormatError for malformed content and ConfigurationError for invalid options.
    """

    return _run_io_format("index", input, options)


def index_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build a lookup index of a PLINK .map variant file: variant IDs with line/byte offsets.

    Maps each variant ID to its line number, byte offset and length in the text, and lists per-chromosome counts and bp
    ranges, so a variant can be located without re-parsing.

    input: path (str or Path) of a PLINK .map variant file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
    A variant ID that occurs on several lines keeps only its last line in 'entries' (record_count still counts every
    line); such IDs are listed in 'duplicate_ids' and reported in 'warnings'.

    Returns the DNAKit envelope (status 'ok'); result: {'format', 'key', 'record_count', 'entries': {...},
    'chromosomes': {...}, 'duplicate_ids': [...]} (a lookup table), plus 'warnings' when an ID is duplicated.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("index_map", input, options)


def index_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build a lookup index of a Newick tree file: trees with byte offsets and the trees containing each tip.

    Gives each tree its byte offset and length in the text and lists, for every tip name, the trees containing it.

    input: path (str or Path) of a Newick tree file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'key', 'record_count', 'entries': {...}} (a lookup
    table).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("index_newick", input, options)


def index_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build a lookup index of a NEXUS file: taxa and trees.

    Maps taxon names to row ordinal, length and non-missing residues, and tree names to their position.

    input: path (str or Path) of a NEXUS file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'key', 'record_count', 'entries': {...}} (a lookup
    table).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("index_nexus", input, options)


def index_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build a lookup index of a PAF (Pairwise mApping Format) alignment file: query or target names with line/byte offsets and intervals.

    Maps each query (or target) name to its alignments' line number, byte offset and length, plus the interval on the
    indexed sequence, sorted by start.

    input: path (str or Path) of a PAF (Pairwise mApping Format) alignment file (optionally gzip-compressed), or its
    content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        by: 'query' (default) or 'target': the name to index alignments by
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'key', 'record_count', 'entries': {...}} (a lookup
    table).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("index_paf", input, options)


def index_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build a lookup index of a PLINK .ped pedigree/genotype file: individuals with line/byte offsets.

    Maps 'family:individual' to its line number, byte offset and length and lists the individuals of each family.

    input: path (str or Path) of a PLINK .ped pedigree/genotype file (optionally gzip-compressed), or its content as
    str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        max_cells: int, default 20000000
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'key', 'record_count', 'entries': {...}} (a lookup
    table).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("index_ped", input, options)


def index_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build a lookup index of a PHYLIP alignment file: taxa.

    Maps each taxon name to its row ordinal, sequence length and number of non-missing residues (PHYLIP rows are not
    stored one per line when interleaved, so no byte offsets are given).

    input: path (str or Path) of a PHYLIP alignment file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        datatype: 'dna' (default), 'rna', 'protein' or 'any': alphabet used for validation and for counting missing
            characters (gap, ?, N for nucleotides, X for protein)
        input_layout: 'auto' (default), 'sequential' or 'interleaved': layout of the PHYLIP input
        input_strict_names: bool, default False: input names occupy the first 10 characters of each row (otherwise the
            first whitespace-delimited word)
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'key', 'record_count', 'entries': {...}} (a lookup
    table).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("index_phylip", input, options)


def index_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build a lookup index of a UCSC 2bit genome file: sequence names with file offsets.

    Per sequence: file offset of its record, length, number of N and mask blocks and the packed-data offset.

    input: path (str or Path) of a .2bit file, or its raw bytes.
    options (all optional unless stated):
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'key', 'record_count', 'entries': {...}} (a lookup
    table).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("index_twobit", input, options)


def index_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build a lookup index of a Zarr directory store: array paths with shape and chunk grid.

    Lists every array path with kind, shape, dtype, chunk shape, chunk grid and chunk counts (no chunk data is read).

    input: path (str or Path) of a zarr directory store (requires the optional 'zarr>=3' and 'numpy' packages).
    options (all optional unless stated):
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'key', 'record_count', 'entries': {...}} (a lookup
    table).
    Returns status 'not_implemented' (requires="pip install 'zarr>=3' numpy") when zarr or numpy is not installed. Raises
    InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size limits.
    """

    return _run_io_format("index_zarr", input, options)


def merge(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Merge several files of the same supported format into one.

    The format is taken from options['format'] (one of map, newick, nexus, paf, ped, phylip, twobit, zarr) or detected
    from the file suffix and content, and the call is executed by the same engine as the format-specific function (e.g.
    filter_paf); its options and result are those of that function, plus the options below.

    input: as for the format-specific function (a path, text/bytes, or - for merge - a list of them).
    options:
        max_bytes: int, default 268435456
        output: path to write (optional)
        overwrite: bool, default False
        format: name of the input format; detected from the suffix/content when omitted
    Returns the DNAKit envelope (status 'ok') whose result is that of the format-specific function. Raises
    InputFormatError for malformed content and ConfigurationError for invalid options.
    """

    return _run_io_format("merge", input, options)


def merge_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Merge several variant lists into one PLINK .map variant file.

    input: a non-empty list of inputs (or {'inputs': [...]}), each being the path (str or Path) of a PLINK .map variant
    file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output: path to write (optional)
        overwrite: bool, default False
        duplicates: 'error' (default) or 'first': what to do when a variant ID occurs twice
        sort: bool, default False: sort the merged variants by chromosome and bp
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("merge_map", input, options)


def merge_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Merge several tree files into one multi-tree file into one Newick tree file.

    input: a non-empty list of inputs (or {'inputs': [...]}), each being the path (str or Path) of a Newick tree file
    (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output: path to write (optional)
        overwrite: bool, default False
        deduplicate: bool, default False: drop trees whose normalised Newick equals an earlier tree
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("merge_newick", input, options)


def merge_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Merge several NEXUS files (concatenate sites per taxon, or append taxa) into one NEXUS file.

    input: a non-empty list of inputs (or {'inputs': [...]}), each being the path (str or Path) of a NEXUS file
    (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output: path to write (optional)
        overwrite: bool, default False
        mode: 'concatenate' (default: join the sites of the same taxa into a supermatrix, filling taxa absent from an
            input with the missing character) or 'append' (add taxa; site counts must be equal)
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("merge_nexus", input, options)


def merge_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Merge several PAF files into one PAF (Pairwise mApping Format) alignment file.

    input: a non-empty list of inputs (or {'inputs': [...]}), each being the path (str or Path) of a PAF (Pairwise mApping
    Format) alignment file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output: path to write (optional)
        overwrite: bool, default False
        deduplicate: bool, default False: drop alignment lines identical to an earlier one
        sort_by: optional sort key as for sort_paf, applied to the merged records
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("merge_paf", input, options)


def merge_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Merge several .ped files with the same variants (rows are concatenated) into one PLINK .ped pedigree/genotype file.

    input: a non-empty list of inputs (or {'inputs': [...]}), each being the path (str or Path) of a PLINK .ped
    pedigree/genotype file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output: path to write (optional)
        overwrite: bool, default False
        duplicates: 'error' (default) or 'first': what to do when a (family, individual) pair occurs twice
        max_cells: int, default 20000000
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("merge_ped", input, options)


def merge_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Merge several PHYLIP alignments (concatenate sites per taxon, or append taxa) into one PHYLIP alignment file.

    input: a non-empty list of inputs (or {'inputs': [...]}), each being the path (str or Path) of a PHYLIP alignment file
    (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output: path to write (optional)
        overwrite: bool, default False
        datatype: 'dna' (default), 'rna', 'protein' or 'any': alphabet used for validation and for counting missing
            characters (gap, ?, N for nucleotides, X for protein)
        input_layout: 'auto' (default), 'sequential' or 'interleaved': layout of the PHYLIP input
        input_strict_names: bool, default False: input names occupy the first 10 characters of each row (otherwise the
            first whitespace-delimited word)
        mode: 'concatenate' (default: join the sites of the same taxa into a supermatrix, filling taxa absent from an
            input with the missing character) or 'append' (add taxa; site counts must be equal)
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("merge_phylip", input, options)


def merge_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Merge several 2bit files into one 2bit file into one UCSC 2bit genome file.

    Concatenates the sequences of several 2bit files into one version-0 file (names must be unique).

    input: a non-empty list of inputs (or {'inputs': [...]}), each being the path (str or Path) of a .2bit file, or its
    raw bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output: path to write (optional)
        overwrite: bool, default False
        duplicates: 'error' (default) or 'first': what to do when a name occurs twice
        byte_order: 'little' (default) or 'big': byte order of the written 2bit file
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("merge_twobit", input, options)


def merge_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Merge several zarr stores by concatenating arrays of the same path into one Zarr directory store.

    Concatenates the arrays of several stores along an axis into a new store. Every store must contain the same array
    paths; shapes must agree except along the axis; dtypes are promoted with numpy.result_type.

    input: a non-empty list of inputs (or {'inputs': [...]}), each being the path (str or Path) of a zarr directory store
    (requires the optional 'zarr>=3' and 'numpy' packages).
    options (all optional unless stated):
        output: directory path of the new zarr store (required)
        overwrite: bool, default False: replace an existing zarr store at output
        axis: int, default 0: axis to concatenate along
        max_elements: int, default 5000000: refuse to materialise more array elements than this
    Returns the DNAKit envelope (status 'ok'); result: {'format': 'zarr', 'record_count', 'output_path'} (the new store is
    written to options['output']).
    Returns status 'not_implemented' (requires="pip install 'zarr>=3' numpy") when zarr or numpy is not installed. Raises
    InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size limits.
    """

    return _run_io_format("merge_zarr", input, options)


def read_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Parse a PLINK .map variant file into structured variants.

    input: path (str or Path) of a PLINK .map variant file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs larger than this
        max_records: int, default 200000: raise instead of returning more records (use stream_* to page)
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'records': [...]} plus format-specific
    fields.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("read_map", input, options)


def read_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Parse a Newick tree file into structured trees.

    Records carry the normalised Newick string, tip names and a flat node table (id, parent, name, length in preorder).

    input: path (str or Path) of a Newick tree file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs larger than this
        max_records: int, default 200000: raise instead of returning more records (use stream_* to page)
        max_nodes: int, default 2000000: maximum nodes per tree
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'records': [...]} plus format-specific
    fields.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("read_newick", input, options)


def read_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Parse a NEXUS file into structured sequences.

    Parses the TAXA, DATA/CHARACTERS and TREES blocks (other blocks are ignored): taxon labels, the matrix (sequential or
    interleaved, matchchar expanded), the declared datatype/missing/gap symbols and the trees (TRANSLATE tables applied,
    [&R]/[&U] kept as the 'rooted' flag). Comments in [] are dropped.

    input: path (str or Path) of a NEXUS file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs larger than this
        max_records: int, default 200000: raise instead of returning more records (use stream_* to page)
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'records': [...]} plus format-specific
    fields.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("read_nexus", input, options)


def read_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Parse a PAF (Pairwise mApping Format) alignment file into structured alignments.

    input: path (str or Path) of a PAF (Pairwise mApping Format) alignment file (optionally gzip-compressed), or its
    content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs larger than this
        max_records: int, default 200000: raise instead of returning more records (use stream_* to page)
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'records': [...]} plus format-specific
    fields.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("read_paf", input, options)


def read_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Parse a PLINK .ped pedigree/genotype file into structured individuals.

    input: path (str or Path) of a PLINK .ped pedigree/genotype file (optionally gzip-compressed), or its content as
    str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs larger than this
        max_records: int, default 200000: raise instead of returning more records (use stream_* to page)
        max_cells: int, default 20000000: refuse files with more genotype cells (individuals x variants)
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'records': [...]} plus format-specific
    fields.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("read_ped", input, options)


def read_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Parse a PHYLIP alignment file into structured sequences.

    Reads relaxed PHYLIP (name = first whitespace-delimited word; blocks of residues may contain spaces) or strict PHYLIP;
    sequential and interleaved layouts are detected automatically. The header 'ntax nchar' must match.

    input: path (str or Path) of a PHYLIP alignment file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs larger than this
        max_records: int, default 200000: raise instead of returning more records (use stream_* to page)
        datatype: 'dna' (default), 'rna', 'protein' or 'any': alphabet used for validation and for counting missing
            characters (gap, ?, N for nucleotides, X for protein)
        input_layout: 'auto' (default), 'sequential' or 'interleaved': layout of the PHYLIP input
        input_strict_names: bool, default False: input names occupy the first 10 characters of each row (otherwise the
            first whitespace-delimited word)
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'records': [...]} plus format-specific
    fields.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("read_phylip", input, options)


def read_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Parse a UCSC 2bit genome file into structured sequences.

    Decodes sequences (or the half-open 0-based region [start, end) of one sequence) to text: N blocks become 'N', mask
    blocks lower case. Result records are {name, length, sequence}.

    input: path (str or Path) of a .2bit file, or its raw bytes.
    options (all optional unless stated):
        names: list of sequence names to read (default all)
        start: int: 0-based start (single name only)
        end: int: 0-based exclusive end (single name only)
        soft_mask: bool, default True: lower-case masked bases
        max_bases: int, default 50000000: refuse to decode more bases
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'records': [...]} plus format-specific
    fields.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("read_twobit", input, options)


def read_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Parse a Zarr directory store into structured arrays.

    Opens a zarr directory store read-only and returns array metadata (shape, dtype, chunks, codecs, attributes) and,
    unless include_data is false, the values as nested lists (bounded by max_elements).

    input: path (str or Path) of a zarr directory store (requires the optional 'zarr>=3' and 'numpy' packages).
    options (all optional unless stated):
        array: array path to read (default: all arrays)
        include_data: bool, default True: return element values
        start: int: first row along axis 0 (single array)
        stop: int: end row (exclusive)
        max_elements: int, default 5000000: refuse to materialise more array elements than this
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'records': [...]} plus format-specific
    fields.
    Returns status 'not_implemented' (requires="pip install 'zarr>=3' numpy") when zarr or numpy is not installed. Raises
    InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size limits.
    """

    return _run_io_format("read_zarr", input, options)


def resume_download(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Download a file over HTTP(S), resuming an interrupted transfer with a Range request.

    Download a file over HTTP(S) into 'output', resuming an interrupted transfer: bytes already stored in '<output>.part'
    are kept and requested again with 'Range: bytes=<n>-' (guarded by If-Range with the ETag/Last-Modified saved in
    '<output>.part.meta'). A server that ignores Range (HTTP 200) makes the download restart from byte 0; HTTP 416 with a
    matching size means the part file is already complete. The file is renamed to 'output' only after its size (and
    expected_sha256, if given) is verified; a transfer cut short keeps the '.part' file and raises DownloadError so that
    calling the function again resumes it.

    input: the URL (http:// or https://) as a string, or a mapping {'url': ..., 'output': ...}.
    options:
        url: http(s) URL (alternative to passing it as input)
        output: destination file path (required unless input is a mapping with 'output')
        overwrite: bool, default False: replace an existing destination file
        chunk_size: int, default 65536: bytes read per step
        timeout: number of seconds, default 30: socket timeout
        max_bytes: int, default 1073741824: refuse to download more than this
        expected_sha256: hex digest the finished file must have (a mismatch deletes the partial file)
    Returns the DNAKit envelope (status 'ok'); result: {'output_path', 'bytes', 'sha256', 'resumed_from', 'restarted',
    'complete'}. Raises DownloadError for HTTP/network failures and ConfigurationError for invalid options. Only http(s)
    URLs are accepted.
    """

    return _run_io_format("resume_download", input, options)


def sort(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Sort a file of any supported format with the ordering of its format-specific sort function.

    The format is taken from options['format'] (one of map, newick, nexus, paf, ped, phylip, twobit, zarr) or detected
    from the file suffix and content, and the call is executed by the same engine as the format-specific function (e.g.
    filter_paf); its options and result are those of that function, plus the options below.

    input: as for the format-specific function (a path, text/bytes, or - for merge - a list of them).
    options:
        max_bytes: int, default 268435456
        output: path to write (optional)
        overwrite: bool, default False
        format: name of the input format; detected from the suffix/content when omitted
    Returns the DNAKit envelope (status 'ok') whose result is that of the format-specific function. Raises
    InputFormatError for malformed content and ConfigurationError for invalid options.
    """

    return _run_io_format("sort", input, options)


def sort_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Sort variants by chromosome and position (or by ID / genetic distance) of a PLINK .map variant file and return the sorted content.

    input: path (str or Path) of a PLINK .map variant file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output: path to write (optional)
        overwrite: bool, default False
        key: 'position' (default: chromosome in natural order, then bp, then ID), 'id' or 'cm'
        reverse: bool, default False
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("sort_map", input, options)


def sort_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Sort the clades of each tree (ladderize or alphabetical order) of a Newick tree file and return the sorted content.

    input: path (str or Path) of a Newick tree file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output: path to write (optional)
        overwrite: bool, default False
        order: 'ascending' (default: clades with fewer tips first), 'descending' or 'alphabetical' (by the smallest tip
            name of each clade); ties are broken by that name
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("sort_newick", input, options)


def sort_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Sort the taxa (matrix rows) by name or missing-data fraction of a NEXUS file and return the sorted content.

    input: path (str or Path) of a NEXUS file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output: path to write (optional)
        overwrite: bool, default False
        key: 'name' (default) or 'missing'
        reverse: bool, default False
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("sort_nexus", input, options)


def sort_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Sort alignment records by target, query, mapping quality, identity, block length or matches of a PAF (Pairwise mApping Format) alignment file and return the sorted content.

    input: path (str or Path) of a PAF (Pairwise mApping Format) alignment file (optionally gzip-compressed), or its
    content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output: path to write (optional)
        overwrite: bool, default False
        key: 'target' (default: target name, start, end), 'query' (query name, start, end), 'mapq', 'identity',
            'block_length' or 'matches'
        reverse: bool, default False
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("sort_paf", input, options)


def sort_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Sort individuals by family, individual ID, sex or phenotype of a PLINK .ped pedigree/genotype file and return the sorted content.

    input: path (str or Path) of a PLINK .ped pedigree/genotype file (optionally gzip-compressed), or its content as
    str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output: path to write (optional)
        overwrite: bool, default False
        key: 'family' (default: family ID then individual ID), 'individual', 'sex' or 'phenotype'
        reverse: bool, default False
        max_cells: int, default 20000000
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("sort_ped", input, options)


def sort_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Sort the taxa (rows) by name or missing-data fraction of a PHYLIP alignment file and return the sorted content.

    input: path (str or Path) of a PHYLIP alignment file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output: path to write (optional)
        overwrite: bool, default False
        datatype: 'dna' (default), 'rna', 'protein' or 'any': alphabet used for validation and for counting missing
            characters (gap, ?, N for nucleotides, X for protein)
        input_layout: 'auto' (default), 'sequential' or 'interleaved': layout of the PHYLIP input
        input_strict_names: bool, default False: input names occupy the first 10 characters of each row (otherwise the
            first whitespace-delimited word)
        key: 'name' (default) or 'missing' (missing fraction, then name)
        reverse: bool, default False
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("sort_phylip", input, options)


def sort_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Sort the sequences of a 2bit file by name, natural chromosome order or length of a UCSC 2bit genome file and return the sorted content.

    Reorders sequences in the output index; the packed records are copied bit-exactly.

    input: path (str or Path) of a .2bit file, or its raw bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output: path to write (optional)
        overwrite: bool, default False
        key: 'name' (default), 'natural' (chr1..22, X, Y, M order) or 'length'
        reverse: bool, default False
        byte_order: 'little' (default) or 'big': byte order of the written 2bit file
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("sort_twobit", input, options)


def sort_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Sort the rows (axis 0) of all arrays by the values of a key array of a Zarr directory store and return the sorted content.

    Reorders axis 0 of every array whose first dimension equals the length of the key array (by) into a new store; other
    arrays are copied unchanged.

    input: path (str or Path) of a zarr directory store (requires the optional 'zarr>=3' and 'numpy' packages).
    options (all optional unless stated):
        output: directory path of the new zarr store (required)
        overwrite: bool, default False: replace an existing zarr store at output
        by: path of the 1-D key array (required); rows are ordered by its values (stable)
        reverse: bool, default False: reverse the ascending order
        max_elements: int, default 5000000: refuse to materialise more array elements than this
    Returns the DNAKit envelope (status 'ok'); result: {'format': 'zarr', 'record_count', 'output_path'} (the new store is
    written to options['output']).
    Returns status 'not_implemented' (requires="pip install 'zarr>=3' numpy") when zarr or numpy is not installed. Raises
    InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size limits.
    """

    return _run_io_format("sort_zarr", input, options)


def split_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Split a PLINK .map variant file per chromosome (or into blocks of N variants).

    input: path (str or Path) of a PLINK .map variant file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output_dir: directory to write the parts to as part_0001_<label><ext> (optional)
        overwrite: bool, default False
        by: 'chromosome' (default; one part per chromosome) or 'size'
        size: int: variants per part when by='size'
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'part_count', 'parts': [{'name', 'record_count', ...}]};
    every part carries its content ('text' / 'data_base64') or, with options['output_dir'], the 'output_path' of the
    written file.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("split_map", input, options)


def split_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Split a Newick tree file into files of N trees.

    input: path (str or Path) of a Newick tree file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output_dir: directory to write the parts to as part_0001_<label><ext> (optional)
        overwrite: bool, default False
        trees_per_part: int, default 1: trees in each output part
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'part_count', 'parts': [{'name', 'record_count', ...}]};
    every part carries its content ('text' / 'data_base64') or, with options['output_dir'], the 'output_path' of the
    written file.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("split_newick", input, options)


def split_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Split a NEXUS file by site blocks, taxon groups, named site partitions or tree groups.

    input: path (str or Path) of a NEXUS file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output_dir: directory to write the parts to as part_0001_<label><ext> (optional)
        overwrite: bool, default False
        site_size: int: consecutive blocks of this many sites
        taxa_size: int: consecutive groups of this many taxa
        partitions: mapping label -> 'start-end[,start-end...]' (1-based, inclusive site ranges)
        trees_per_part: int: trees per part (documents without a matrix); give exactly one of the four split options
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'part_count', 'parts': [{'name', 'record_count', ...}]};
    every part carries its content ('text' / 'data_base64') or, with options['output_dir'], the 'output_path' of the
    written file.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("split_nexus", input, options)


def split_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Split a PAF (Pairwise mApping Format) alignment file per query, target or strand (or into blocks of N alignments).

    input: path (str or Path) of a PAF (Pairwise mApping Format) alignment file (optionally gzip-compressed), or its
    content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output_dir: directory to write the parts to as part_0001_<label><ext> (optional)
        overwrite: bool, default False
        by: 'query' (default), 'target', 'strand' or 'size'
        size: int: alignments per part when by='size'
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'part_count', 'parts': [{'name', 'record_count', ...}]};
    every part carries its content ('text' / 'data_base64') or, with options['output_dir'], the 'output_path' of the
    written file.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("split_paf", input, options)


def split_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Split a PLINK .ped pedigree/genotype file per family (or into blocks of N individuals).

    input: path (str or Path) of a PLINK .ped pedigree/genotype file (optionally gzip-compressed), or its content as
    str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output_dir: directory to write the parts to as part_0001_<label><ext> (optional)
        overwrite: bool, default False
        by: 'family' (default; one part per family) or 'size'
        size: int: individuals per part when by='size'
        max_cells: int, default 20000000
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'part_count', 'parts': [{'name', 'record_count', ...}]};
    every part carries its content ('text' / 'data_base64') or, with options['output_dir'], the 'output_path' of the
    written file.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("split_ped", input, options)


def split_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Split a PHYLIP alignment file by site blocks, taxon groups or named site partitions.

    input: path (str or Path) of a PHYLIP alignment file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output_dir: directory to write the parts to as part_0001_<label><ext> (optional)
        overwrite: bool, default False
        datatype: 'dna' (default), 'rna', 'protein' or 'any': alphabet used for validation and for counting missing
            characters (gap, ?, N for nucleotides, X for protein)
        input_layout: 'auto' (default), 'sequential' or 'interleaved': layout of the PHYLIP input
        input_strict_names: bool, default False: input names occupy the first 10 characters of each row (otherwise the
            first whitespace-delimited word)
        site_size: int: consecutive blocks of this many sites
        taxa_size: int: consecutive groups of this many taxa
        partitions: mapping label -> 'start-end[,start-end...]' (1-based, inclusive site ranges)
        Give exactly one of site_size, taxa_size and partitions.
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'part_count', 'parts': [{'name', 'record_count', ...}]};
    every part carries its content ('text' / 'data_base64') or, with options['output_dir'], the 'output_path' of the
    written file.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("split_phylip", input, options)


def split_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Split a UCSC 2bit genome file into 2bit files of N sequences.

    Writes groups of sequences to separate 2bit files; each result part is a complete 2bit file.

    input: path (str or Path) of a .2bit file, or its raw bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        output_dir: directory to write the parts to as part_0001_<label><ext> (optional)
        overwrite: bool, default False
        sequences_per_part: int, default 1: sequences in each output part
        byte_order: 'little' (default) or 'big': byte order of the written 2bit file
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'part_count', 'parts': [{'name', 'record_count', ...}]};
    every part carries its content ('text' / 'data_base64') or, with options['output_dir'], the 'output_path' of the
    written file.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("split_twobit", input, options)


def split_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Split a Zarr directory store into stores of N rows (axis 0).

    Cuts axis 0 of all arrays into consecutive blocks of rows_per_part rows and writes one store per block.

    input: path (str or Path) of a zarr directory store (requires the optional 'zarr>=3' and 'numpy' packages).
    options (all optional unless stated):
        output_dir: directory that receives part_0001_<label>.zarr stores (required)
        overwrite: bool, default False
        rows_per_part: int (required): rows along axis 0 per part
        max_elements: int, default 5000000: refuse to materialise more array elements than this
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'part_count', 'parts': [{'name', 'record_count', ...}]};
    every part carries its content ('text' / 'data_base64') or, with options['output_dir'], the 'output_path' of the
    written file.
    Returns status 'not_implemented' (requires="pip install 'zarr>=3' numpy") when zarr or numpy is not installed. Raises
    InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size limits.
    """

    return _run_io_format("split_zarr", input, options)


def stream(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Return one page of records from a file of any supported format, with a cursor for the next page.

    The format is taken from options['format'] (one of map, newick, nexus, paf, ped, phylip, twobit, zarr) or detected
    from the file suffix and content, and the call is executed by the same engine as the format-specific function (e.g.
    filter_paf); its options and result are those of that function, plus the options below.

    input: as for the format-specific function (a path, text/bytes, or - for merge - a list of them).
    options:
        max_bytes: int, default 268435456
        start: int, default 0: index of the first record of the page
        limit: int, default 1000: records per page (max 100000)
        format: name of the input format; detected from the suffix/content when omitted
    Returns the DNAKit envelope (status 'ok') whose result is that of the format-specific function. Raises
    InputFormatError for malformed content and ConfigurationError for invalid options.
    """

    return _run_io_format("stream", input, options)


def stream_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Return one page of variants from a PLINK .map variant file, with a cursor for the next page.

    input: path (str or Path) of a PLINK .map variant file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        start: int, default 0: index of the first record of the page
        limit: int, default 1000: records per page (max 100000)
    Returns the DNAKit envelope (status 'ok'); result: {'start', 'limit', 'returned', 'total', 'next_start', 'done',
    'records'}; pass next_start as options['start'] for the next page.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("stream_map", input, options)


def stream_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Return one page of trees from a Newick tree file, with a cursor for the next page.

    input: path (str or Path) of a Newick tree file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        start: int, default 0: index of the first record of the page
        limit: int, default 1000: records per page (max 100000)
    Returns the DNAKit envelope (status 'ok'); result: {'start', 'limit', 'returned', 'total', 'next_start', 'done',
    'records'}; pass next_start as options['start'] for the next page.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("stream_newick", input, options)


def stream_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Return one page of sequences from a NEXUS file, with a cursor for the next page.

    input: path (str or Path) of a NEXUS file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        start: int, default 0: index of the first record of the page
        limit: int, default 1000: records per page (max 100000)
    Returns the DNAKit envelope (status 'ok'); result: {'start', 'limit', 'returned', 'total', 'next_start', 'done',
    'records'}; pass next_start as options['start'] for the next page.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("stream_nexus", input, options)


def stream_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Return one page of alignments from a PAF (Pairwise mApping Format) alignment file, with a cursor for the next page.

    input: path (str or Path) of a PAF (Pairwise mApping Format) alignment file (optionally gzip-compressed), or its
    content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        start: int, default 0: index of the first record of the page
        limit: int, default 1000: records per page (max 100000)
    Returns the DNAKit envelope (status 'ok'); result: {'start', 'limit', 'returned', 'total', 'next_start', 'done',
    'records'}; pass next_start as options['start'] for the next page.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("stream_paf", input, options)


def stream_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Return one page of individuals from a PLINK .ped pedigree/genotype file, with a cursor for the next page.

    input: path (str or Path) of a PLINK .ped pedigree/genotype file (optionally gzip-compressed), or its content as
    str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        start: int, default 0: index of the first record of the page
        limit: int, default 1000: records per page (max 100000)
        max_cells: int, default 20000000
    Returns the DNAKit envelope (status 'ok'); result: {'start', 'limit', 'returned', 'total', 'next_start', 'done',
    'records'}; pass next_start as options['start'] for the next page.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("stream_ped", input, options)


def stream_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Return one page of sequences from a PHYLIP alignment file, with a cursor for the next page.

    input: path (str or Path) of a PHYLIP alignment file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        start: int, default 0: index of the first record of the page
        limit: int, default 1000: records per page (max 100000)
        datatype: 'dna' (default), 'rna', 'protein' or 'any': alphabet used for validation and for counting missing
            characters (gap, ?, N for nucleotides, X for protein)
        input_layout: 'auto' (default), 'sequential' or 'interleaved': layout of the PHYLIP input
        input_strict_names: bool, default False: input names occupy the first 10 characters of each row (otherwise the
            first whitespace-delimited word)
    Returns the DNAKit envelope (status 'ok'); result: {'start', 'limit', 'returned', 'total', 'next_start', 'done',
    'records'}; pass next_start as options['start'] for the next page.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("stream_phylip", input, options)


def stream_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Return one page of sequences from a UCSC 2bit genome file, with a cursor for the next page.

    Pages through fixed-size windows over all sequences in file order; each record is {name, start, end, sequence}. The
    cursor 'start' counts windows; only the windows of a page are decoded.

    input: path (str or Path) of a .2bit file, or its raw bytes.
    options (all optional unless stated):
        start: int, default 0: index of the first record of the page
        limit: int, default 1000: records per page (max 100000)
        window: int, default 10000: bases per window (the unit that is paged)
        soft_mask: bool, default True
        max_bases: int, default 50000000: bases per page
    Returns the DNAKit envelope (status 'ok'); result: {'start', 'limit', 'returned', 'total', 'next_start', 'done',
    'records'}; pass next_start as options['start'] for the next page.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("stream_twobit", input, options)


def stream_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Return one page of arrays from a Zarr directory store, with a cursor for the next page.

    Pages through the rows along axis 0 of one array: records are the rows as nested lists; the cursor 'start' counts
    rows.

    input: path (str or Path) of a zarr directory store (requires the optional 'zarr>=3' and 'numpy' packages).
    options (all optional unless stated):
        start: int, default 0: index of the first record of the page
        limit: int, default 1000: records per page (max 100000)
        array: array path to page through (required when the store has several arrays)
        max_elements: int, default 5000000: maximum elements per page
    Returns the DNAKit envelope (status 'ok'); result: {'start', 'limit', 'returned', 'total', 'next_start', 'done',
    'records'}; pass next_start as options['start'] for the next page.
    Returns status 'not_implemented' (requires="pip install 'zarr>=3' numpy") when zarr or numpy is not installed. Raises
    InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size limits.
    """

    return _run_io_format("stream_zarr", input, options)


def summarize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Summarise a file of any supported format.

    The format is taken from options['format'] (one of map, newick, nexus, paf, ped, phylip, twobit, zarr) or detected
    from the file suffix and content, and the call is executed by the same engine as the format-specific function (e.g.
    filter_paf); its options and result are those of that function, plus the options below.

    input: as for the format-specific function (a path, text/bytes, or - for merge - a list of them).
    options:
        max_bytes: int, default 268435456
        format: name of the input format; detected from the suffix/content when omitted
    Returns the DNAKit envelope (status 'ok') whose result is that of the format-specific function. Raises
    InputFormatError for malformed content and ConfigurationError for invalid options.
    """

    return _run_io_format("summarize", input, options)


def summarize_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Summarise a PLINK .map variant file: variant and chromosome counts, bp ranges and duplicate IDs.

    input: path (str or Path) of a PLINK .map variant file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
    Returns the DNAKit envelope (status 'ok'); result: a dict of statistics.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("summarize_map", input, options)


def summarize_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Summarise a Newick tree file: per-tree node counts, depth, height, total branch length and resolution.

    Per tree: tip/internal/node counts, root degree, rooted (= root has exactly 2 children), fully-resolved flag, depth in
    edges, height (longest root-to-tip path length; missing lengths count as 0), total branch length.

    input: path (str or Path) of a Newick tree file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
    Returns the DNAKit envelope (status 'ok'); result: a dict of statistics.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("summarize_newick", input, options)


def summarize_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Summarise a NEXUS file: alignment statistics (constant / variable / parsimony-informative sites, composition) and tree statistics.

    Matrix statistics as for PHYLIP (constant/variable/parsimony-informative sites, composition, missing fractions) plus
    per-tree statistics when TREES blocks are present.

    input: path (str or Path) of a NEXUS file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
    Returns the DNAKit envelope (status 'ok'); result: a dict of statistics.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("summarize_nexus", input, options)


def summarize_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Summarise a PAF (Pairwise mApping Format) alignment file: alignment counts, identity, aligned bases, strands and mapping-quality histogram.

    Counts alignments, distinct queries/targets, strands, aligned bases (sum of block lengths), identity = matches /
    block_length (mean and base-weighted overall), and a mapq histogram.

    input: path (str or Path) of a PAF (Pairwise mApping Format) alignment file (optionally gzip-compressed), or its
    content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
    Returns the DNAKit envelope (status 'ok'); result: a dict of statistics.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("summarize_paf", input, options)


def summarize_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Summarise a PLINK .ped pedigree/genotype file: individual, family and variant counts, sex, phenotype classes and missing-genotype rates.

    input: path (str or Path) of a PLINK .ped pedigree/genotype file (optionally gzip-compressed), or its content as
    str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        max_cells: int, default 20000000
    Returns the DNAKit envelope (status 'ok'); result: a dict of statistics.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("summarize_ped", input, options)


def summarize_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Summarise a PHYLIP alignment file: alignment statistics (constant / variable / parsimony-informative sites, composition, missing data).

    Taxa, sites, constant/variable/parsimony-informative sites (a site is informative when at least two different
    unambiguous states each occur in at least two taxa; gaps and ambiguity codes are ignored), composition, GC content
    (nucleotides), per-taxon missing fractions and duplicate sequences.

    input: path (str or Path) of a PHYLIP alignment file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456
        datatype: 'dna' (default), 'rna', 'protein' or 'any': alphabet used for validation and for counting missing
            characters (gap, ?, N for nucleotides, X for protein)
        input_layout: 'auto' (default), 'sequential' or 'interleaved': layout of the PHYLIP input
        input_strict_names: bool, default False: input names occupy the first 10 characters of each row (otherwise the
            first whitespace-delimited word)
    Returns the DNAKit envelope (status 'ok'); result: a dict of statistics.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("summarize_phylip", input, options)


def summarize_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Summarise a UCSC 2bit genome file: per-sequence length, base composition, N and soft-masked bases, GC fraction.

    Per sequence: length, A/C/G/T/N counts, soft-masked bases, GC fraction = (G+C)/(A+C+G+T), and file totals.

    input: path (str or Path) of a .2bit file, or its raw bytes.
    options (all optional unless stated):
    Returns the DNAKit envelope (status 'ok'); result: a dict of statistics.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("summarize_twobit", input, options)


def summarize_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Summarise a Zarr directory store: array shapes, dtypes, chunking and totals.

    Counts groups/arrays and totals elements, bytes and initialised chunks; one entry per array.

    input: path (str or Path) of a zarr directory store (requires the optional 'zarr>=3' and 'numpy' packages).
    options (all optional unless stated):
    Returns the DNAKit envelope (status 'ok'); result: a dict of statistics.
    Returns status 'not_implemented' (requires="pip install 'zarr>=3' numpy") when zarr or numpy is not installed. Raises
    InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size limits.
    """

    return _run_io_format("summarize_zarr", input, options)


def validate_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Validate a PLINK .map variant file and report its errors and warnings.

    Each non-blank line needs 4 whitespace-separated columns: chromosome, variant ID, genetic distance in cM (finite
    number) and base-pair position (integer >= 0). Duplicate IDs are errors; duplicate (chromosome, bp) pairs are
    warnings. Every problem of a line is reported (line number, field name and column, offending value).

    input: path (str or Path) of a PLINK .map variant file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs larger than this
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'valid': bool, 'errors': [...], 'warnings': [...]} plus
    counts; content problems are reported, not raised.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("validate_map", input, options)


def validate_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Validate a Newick tree file and report its errors and warnings.

    Each ';'-terminated tree must have balanced parentheses, well-formed labels and finite numeric branch lengths.
    Duplicate or missing tip names are errors; negative branch lengths and single-child nodes are warnings.

    input: path (str or Path) of a Newick tree file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs larger than this
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'valid': bool, 'errors': [...], 'warnings': [...]} plus
    counts; content problems are reported, not raised.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("validate_newick", input, options)


def validate_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Validate a NEXUS file and report its errors and warnings.

    Checks the #NEXUS header, block structure, declared NTAX/NCHAR against the matrix, equal row lengths, unique taxa, the
    alphabet of nucleotide/protein data and (as warnings) tree tips absent from the matrix.

    input: path (str or Path) of a NEXUS file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs larger than this
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'valid': bool, 'errors': [...], 'warnings': [...]} plus
    counts; content problems are reported, not raised.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("validate_nexus", input, options)


def validate_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Validate a PAF (Pairwise mApping Format) alignment file and report its errors and warnings.

    Each non-blank line needs the 12 mandatory tab-separated PAF columns (query name/length/start/end, strand, target
    name/length/start/end, matching bases, block length, mapping quality) with start <= end <= length, matches <= block
    length, strand +/- and mapq in [0,255]; optional columns must look like TAG:TYPE:VALUE. Every problem of a line is
    reported (line number, field name and column, offending value).

    input: path (str or Path) of a PAF (Pairwise mApping Format) alignment file (optionally gzip-compressed), or its
    content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs larger than this
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'valid': bool, 'errors': [...], 'warnings': [...]} plus
    counts; content problems are reported, not raised.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("validate_paf", input, options)


def validate_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Validate a PLINK .ped pedigree/genotype file and report its errors and warnings.

    Rows need 6 + 2k whitespace-separated columns (family, individual, father, mother, sex, phenotype, then two alleles
    per variant) with the same k on every row, numeric phenotype and integer sex. Duplicate (family, individual) pairs are
    errors; sex outside {0,1,2} and parents absent from their family are warnings. Allele '0' marks a missing call. Every problem of a row is reported
    (line number, field name and column, offending value).

    input: path (str or Path) of a PLINK .ped pedigree/genotype file (optionally gzip-compressed), or its content as
    str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs larger than this
        variant_count: int: expected number of variants (e.g. the line count of the matching .map)
        max_cells: int, default 20000000
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'valid': bool, 'errors': [...], 'warnings': [...]} plus
    counts; content problems are reported, not raised.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("validate_ped", input, options)


def validate_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Validate a PHYLIP alignment file and report its errors and warnings.

    Checks the header counts, equal sequence lengths, unique taxon names and (by datatype) the alphabet.

    input: path (str or Path) of a PHYLIP alignment file (optionally gzip-compressed), or its content as str/bytes.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse inputs larger than this
        datatype: 'dna' (default), 'rna', 'protein' or 'any': alphabet used for validation and for counting missing
            characters (gap, ?, N for nucleotides, X for protein)
        input_layout: 'auto' (default), 'sequential' or 'interleaved': layout of the PHYLIP input
        input_strict_names: bool, default False: input names occupy the first 10 characters of each row (otherwise the
            first whitespace-delimited word)
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'valid': bool, 'errors': [...], 'warnings': [...]} plus
    counts; content problems are reported, not raised.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("validate_phylip", input, options)


def validate_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Validate a UCSC 2bit genome file and report its errors and warnings.

    Checks the signature, version 0, reserved fields, index offsets, unique non-empty names, that N/mask blocks are
    sorted, non-empty, non-overlapping and inside the sequence, and that packed data fits in the file.

    input: path (str or Path) of a .2bit file, or its raw bytes.
    options (all optional unless stated):
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'valid': bool, 'errors': [...], 'warnings': [...]} plus
    counts; content problems are reported, not raised.
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("validate_twobit", input, options)


def validate_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Validate a Zarr directory store and report its errors and warnings.

    Opens the store, checks every array's shape, chunk shape and chunk counts, and decodes the first chunk of each
    initialised array so that corrupt codecs or truncated chunk files are detected.

    input: path (str or Path) of a zarr directory store (requires the optional 'zarr>=3' and 'numpy' packages).
    options (all optional unless stated):
        max_chunk_bytes: int, default 67108864: first chunks larger than this are not decoded
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'valid': bool, 'errors': [...], 'warnings': [...]} plus
    counts; content problems are reported, not raised.
    Returns status 'not_implemented' (requires="pip install 'zarr>=3' numpy") when zarr or numpy is not installed. Raises
    InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size limits.
    """

    return _run_io_format("validate_zarr", input, options)


def write_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Serialise variants (or normalise existing PLINK .map variant file content) into a PLINK .map variant file.

    input: a list of {'chrom','id','cm','pos'} dicts or [chrom, id, cm, bp] / [chrom, id, bp] lists, or .map text/path to
    normalise.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse text/file inputs larger than this
        output: path to write (optional; the content is returned inline otherwise)
        overwrite: bool, default False: allow replacing an existing output file
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("write_map", input, options)


def write_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Serialise trees (or normalise existing Newick tree file content) into a Newick tree file.

    input: Newick text/path, a {'newick': ...} dict, a flat {'nodes': [{'id','parent','name','length'}, ...]} dict or a
    nested {'name','length','children'} dict, or a list of these.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse text/file inputs larger than this
        output: path to write (optional; the content is returned inline otherwise)
        overwrite: bool, default False: allow replacing an existing output file
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("write_newick", input, options)


def write_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Serialise sequences (or normalise existing NEXUS file content) into a NEXUS file.

    Input may be FASTA/PHYLIP/NEXUS text or path, a {name: sequence} mapping, or a list of {'name','sequence'} records;
    written as a DATA block (names quoted when needed).

    input: FASTA/PHYLIP/NEXUS text or path, a {name: sequence} mapping, or a list of {'name','sequence'} records.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse text/file inputs larger than this
        output: path to write (optional; the content is returned inline otherwise)
        overwrite: bool, default False: allow replacing an existing output file
        layout: 'sequential' (default) or 'interleaved' for PHYLIP/NEXUS output
        block_size: int, default 60: sites per block of an interleaved output
        strict_names: bool, default False: write strict 10-character PHYLIP names (error when longer)
        line_width: int, default 60: FASTA line width (0 = one line per sequence)
        datatype: 'dna' (default), 'rna', 'protein' or 'any' for records without a declared type
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("write_nexus", input, options)


def write_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Serialise alignments (or normalise existing PAF (Pairwise mApping Format) alignment file content) into a PAF (Pairwise mApping Format) alignment file.

    input: a list of dicts with the 12 PAF fields (query_name, query_length, query_start, query_end, strand, target_name,
    target_length, target_start, target_end, matches, block_length, mapq) and optional 'tags', or PAF text/path.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse text/file inputs larger than this
        output: path to write (optional; the content is returned inline otherwise)
        overwrite: bool, default False: allow replacing an existing output file
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("write_paf", input, options)


def write_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Serialise individuals (or normalise existing PLINK .ped pedigree/genotype file content) into a PLINK .ped pedigree/genotype file.

    input: a list of {'fid','iid','pat','mat','sex','pheno','genotypes': [[a1, a2], ...]} dicts, or .ped text/path to
    normalise.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse text/file inputs larger than this
        output: path to write (optional; the content is returned inline otherwise)
        overwrite: bool, default False: allow replacing an existing output file
        max_cells: int, default 20000000
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("write_ped", input, options)


def write_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Serialise sequences (or normalise existing PHYLIP alignment file content) into a PHYLIP alignment file.

    Input may be FASTA/PHYLIP/NEXUS text or path, a {name: sequence} mapping, or a list of {'name','sequence'} records;
    the alignment is written as relaxed PHYLIP (names padded to at least 10 characters) unless strict_names is true.

    input: FASTA/PHYLIP/NEXUS text or path, a {name: sequence} mapping, or a list of {'name','sequence'} records.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse text/file inputs larger than this
        output: path to write (optional; the content is returned inline otherwise)
        overwrite: bool, default False: allow replacing an existing output file
        datatype: 'dna' (default), 'rna', 'protein' or 'any': alphabet the sequences are checked against (IUPAC codes,
            gap, ?, N for nucleotides; X for protein; 'any' skips the check)
        layout: 'sequential' (default) or 'interleaved' for PHYLIP/NEXUS output
        block_size: int, default 60: sites per block of an interleaved output
        strict_names: bool, default False: write strict 10-character PHYLIP names (error when longer)
        line_width: int, default 60: FASTA line width (0 = one line per sequence)
        input_layout: 'auto' (default), 'sequential' or 'interleaved' when the input is PHYLIP text
        input_strict_names: bool, default False when the input is PHYLIP text
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content or characters outside the alphabet of datatype (the offending
    characters are listed) and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("write_phylip", input, options)


def write_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Serialise sequences (or normalise existing UCSC 2bit genome file content) into a UCSC 2bit genome file.

    Packs sequences into a version-0 2bit file. Input: FASTA text/path, a {name: sequence} mapping or a list of
    {'name','sequence'} records. Non-ACGT letters become N blocks; lower case becomes mask blocks.

    input: FASTA text/path, a {name: sequence} mapping, a list of {'name','sequence'} records, or an existing 2bit file.
    options (all optional unless stated):
        max_bytes: int, default 268435456: refuse text/file inputs larger than this
        output: path to write (optional; the content is returned inline otherwise)
        overwrite: bool, default False: allow replacing an existing output file
        byte_order: 'little' (default) or 'big': byte order of the written 2bit file
    Returns the DNAKit envelope (status 'ok'); result: {'format', 'record_count', 'bytes', 'sha256'} plus either the
    content ('text', or 'data_base64' for binary formats) or, when options['output'] is given, 'output_path' (the file is
    written; an existing file needs overwrite=True).
    Raises InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size
    limits.
    """

    return _run_io_format("write_twobit", input, options)


def write_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Create a Zarr directory store from named arrays.

    Creates a zarr v3 (or v2) directory store from a {name: array-like} mapping, a list of {'name','data'} records or a
    single array-like; requires 'output'. Supported dtypes: bool, int, uint, float.

    input: a {array name: array-like} mapping, a list of {'name','data'} records, or one array-like (stored as
    options['array']).
    options (all optional unless stated):
        output: directory path of the new zarr store (required)
        overwrite: bool, default False: replace an existing zarr store at output
        array: name of the array when the input is a single array-like (default 'data')
        chunks: int or list of ints: chunk shape (default: chosen by zarr)
        attrs: mapping stored as attributes of the root group
        zarr_format: 2 or 3 (default 3)
    Returns the DNAKit envelope (status 'ok'); result: {'format': 'zarr', 'record_count', 'output_path'} (the new store is
    written to options['output']).
    Returns status 'not_implemented' (requires="pip install 'zarr>=3' numpy") when zarr or numpy is not installed. Raises
    InputFormatError for malformed content and ConfigurationError for invalid options or inputs above the size limits.
    """

    return _run_io_format("write_zarr", input, options)


__all__ += [
    "compress_map",
    "compress_newick",
    "compress_nexus",
    "compress_paf",
    "compress_ped",
    "compress_phylip",
    "compress_twobit",
    "compress_zarr",
    "convert",
    "convert_map",
    "convert_newick",
    "convert_nexus",
    "convert_paf",
    "convert_ped",
    "convert_phylip",
    "convert_twobit",
    "convert_zarr",
    "decompress_map",
    "decompress_newick",
    "decompress_nexus",
    "decompress_paf",
    "decompress_ped",
    "decompress_phylip",
    "decompress_twobit",
    "decompress_zarr",
    "filter",
    "filter_map",
    "filter_newick",
    "filter_nexus",
    "filter_paf",
    "filter_ped",
    "filter_phylip",
    "filter_twobit",
    "filter_zarr",
    "index",
    "index_map",
    "index_newick",
    "index_nexus",
    "index_paf",
    "index_ped",
    "index_phylip",
    "index_twobit",
    "index_zarr",
    "merge",
    "merge_map",
    "merge_newick",
    "merge_nexus",
    "merge_paf",
    "merge_ped",
    "merge_phylip",
    "merge_twobit",
    "merge_zarr",
    "read_map",
    "read_newick",
    "read_nexus",
    "read_paf",
    "read_ped",
    "read_phylip",
    "read_twobit",
    "read_zarr",
    "resume_download",
    "sort",
    "sort_map",
    "sort_newick",
    "sort_nexus",
    "sort_paf",
    "sort_ped",
    "sort_phylip",
    "sort_twobit",
    "sort_zarr",
    "split_map",
    "split_newick",
    "split_nexus",
    "split_paf",
    "split_ped",
    "split_phylip",
    "split_twobit",
    "split_zarr",
    "stream",
    "stream_map",
    "stream_newick",
    "stream_nexus",
    "stream_paf",
    "stream_ped",
    "stream_phylip",
    "stream_twobit",
    "stream_zarr",
    "summarize",
    "summarize_map",
    "summarize_newick",
    "summarize_nexus",
    "summarize_paf",
    "summarize_ped",
    "summarize_phylip",
    "summarize_twobit",
    "summarize_zarr",
    "validate_map",
    "validate_newick",
    "validate_nexus",
    "validate_paf",
    "validate_ped",
    "validate_phylip",
    "validate_twobit",
    "validate_zarr",
    "write_map",
    "write_newick",
    "write_nexus",
    "write_paf",
    "write_ped",
    "write_phylip",
    "write_twobit",
    "write_zarr",
]

__candidate_adapter_operations__ = set(
    [
        "compress_map",
        "compress_newick",
        "compress_nexus",
        "compress_paf",
        "compress_ped",
        "compress_phylip",
        "compress_twobit",
        "compress_zarr",
        "convert",
        "convert_map",
        "convert_newick",
        "convert_nexus",
        "convert_paf",
        "convert_ped",
        "convert_phylip",
        "convert_twobit",
        "convert_zarr",
        "decompress_map",
        "decompress_newick",
        "decompress_nexus",
        "decompress_paf",
        "decompress_ped",
        "decompress_phylip",
        "decompress_twobit",
        "decompress_zarr",
        "filter",
        "filter_map",
        "filter_newick",
        "filter_nexus",
        "filter_paf",
        "filter_ped",
        "filter_phylip",
        "filter_twobit",
        "filter_zarr",
        "index",
        "index_map",
        "index_newick",
        "index_nexus",
        "index_paf",
        "index_ped",
        "index_phylip",
        "index_twobit",
        "index_zarr",
        "merge",
        "merge_map",
        "merge_newick",
        "merge_nexus",
        "merge_paf",
        "merge_ped",
        "merge_phylip",
        "merge_twobit",
        "merge_zarr",
        "read_map",
        "read_newick",
        "read_nexus",
        "read_paf",
        "read_ped",
        "read_phylip",
        "read_twobit",
        "read_zarr",
        "resume_download",
        "sort",
        "sort_map",
        "sort_newick",
        "sort_nexus",
        "sort_paf",
        "sort_ped",
        "sort_phylip",
        "sort_twobit",
        "sort_zarr",
        "split_map",
        "split_newick",
        "split_nexus",
        "split_paf",
        "split_ped",
        "split_phylip",
        "split_twobit",
        "split_zarr",
        "stream",
        "stream_map",
        "stream_newick",
        "stream_nexus",
        "stream_paf",
        "stream_ped",
        "stream_phylip",
        "stream_twobit",
        "stream_zarr",
        "summarize",
        "summarize_map",
        "summarize_newick",
        "summarize_nexus",
        "summarize_paf",
        "summarize_ped",
        "summarize_phylip",
        "summarize_twobit",
        "summarize_zarr",
        "validate_map",
        "validate_newick",
        "validate_nexus",
        "validate_paf",
        "validate_ped",
        "validate_phylip",
        "validate_twobit",
        "validate_zarr",
        "write_map",
        "write_newick",
        "write_nexus",
        "write_paf",
        "write_ped",
        "write_phylip",
        "write_twobit",
        "write_zarr",
    ]
)


from dnakit._verbosity import register_public as _register_public_verbose

_register_public_verbose(globals())
