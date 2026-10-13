"""Local file and metadata candidate implementations."""

from __future__ import annotations

import gzip
import hashlib
import json
import mimetypes
from pathlib import Path
from typing import Mapping, Sequence


def _bytes(value: object) -> bytes:
    if isinstance(value, bytes):
        return value
    if isinstance(value, Path) or (isinstance(value, str) and "\n" not in value and Path(value).exists()):
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
    return {"path": path, "size_bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest(), "media_type": mimetypes.guess_type(path or "")[0]}


def chunk_large_file(input: object, *, chunk_size: int = 1024) -> list[bytes]:
    """Split a payload into bounded byte chunks."""

    if chunk_size < 1:
        raise ValueError("chunk_size must be positive")
    payload = _bytes(input)
    return [payload[index : index + chunk_size] for index in range(0, len(payload), chunk_size)]


def clean_temporary_files(input: str | Path, *, suffixes: Sequence[str] = (".tmp", ".part")) -> list[str]:
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
    signatures = ((b"\x1f\x8b", "gzip"), (b"PK\x03\x04", "zip"), (b"BZh", "bzip2"), (b"\xfd7zXZ", "xz"))
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
    return {"valid": valid, "size_bytes": len(payload), "line_count": len(text.splitlines()), "error": error}


def validate_metadata_schema(input: Mapping[str, object], *, required: Sequence[str] = ()) -> dict[str, object]:
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
    return {key: [str(categories[key][int(index)]) for index in values] for key, values in data.items()}  # type: ignore[index]


def normalize_sample_metadata(input: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    """Normalize sample IDs and trim string metadata values."""

    result = []
    for index, row in enumerate(input, start=1):
        normalized = {str(key).strip(): value.strip() if isinstance(value, str) else value for key, value in row.items()}
        normalized.setdefault("sample_id", f"sample_{index}")
        result.append(normalized)
    return result


def merge_metadata_tables(input: Sequence[Sequence[Mapping[str, object]]]) -> list[dict[str, object]]:
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
    "build_file_manifest", "chunk_large_file", "clean_temporary_files", "compress", "decompress",
    "detect_compression", "detect_file_format", "encode_metadata_categorical", "decode_metadata_categorical",
    "estimate_storage_size", "get_file_checksum", "list_file_records", "merge_metadata_tables",
    "normalize_sample_metadata", "read", "repack_compressed_file", "stream_records", "validate",
    "validate_metadata_schema", "verify_checksum", "write",
]


from dnakit._candidate_runtime import run_candidate_operation

def compress_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'compress_map', input, options=options)

def compress_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'compress_newick', input, options=options)

def compress_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'compress_nexus', input, options=options)

def compress_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'compress_paf', input, options=options)

def compress_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'compress_ped', input, options=options)

def compress_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'compress_phylip', input, options=options)

def compress_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'compress_twobit', input, options=options)

def compress_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'compress_zarr', input, options=options)

def convert(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'convert', input, options=options)

def convert_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'convert_map', input, options=options)

def convert_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'convert_newick', input, options=options)

def convert_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'convert_nexus', input, options=options)

def convert_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'convert_paf', input, options=options)

def convert_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'convert_ped', input, options=options)

def convert_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'convert_phylip', input, options=options)

def convert_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'convert_twobit', input, options=options)

def convert_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'convert_zarr', input, options=options)

def decompress_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'decompress_map', input, options=options)

def decompress_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'decompress_newick', input, options=options)

def decompress_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'decompress_nexus', input, options=options)

def decompress_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'decompress_paf', input, options=options)

def decompress_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'decompress_ped', input, options=options)

def decompress_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'decompress_phylip', input, options=options)

def decompress_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'decompress_twobit', input, options=options)

def decompress_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'decompress_zarr', input, options=options)

def filter(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'filter', input, options=options)

def filter_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'filter_map', input, options=options)

def filter_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'filter_newick', input, options=options)

def filter_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'filter_nexus', input, options=options)

def filter_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'filter_paf', input, options=options)

def filter_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'filter_ped', input, options=options)

def filter_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'filter_phylip', input, options=options)

def filter_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'filter_twobit', input, options=options)

def filter_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'filter_zarr', input, options=options)

def index(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'index', input, options=options)

def index_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'index_map', input, options=options)

def index_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'index_newick', input, options=options)

def index_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'index_nexus', input, options=options)

def index_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'index_paf', input, options=options)

def index_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'index_ped', input, options=options)

def index_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'index_phylip', input, options=options)

def index_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'index_twobit', input, options=options)

def index_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'index_zarr', input, options=options)

def merge(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'merge', input, options=options)

def merge_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'merge_map', input, options=options)

def merge_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'merge_newick', input, options=options)

def merge_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'merge_nexus', input, options=options)

def merge_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'merge_paf', input, options=options)

def merge_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'merge_ped', input, options=options)

def merge_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'merge_phylip', input, options=options)

def merge_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'merge_twobit', input, options=options)

def merge_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'merge_zarr', input, options=options)

def read_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'read_map', input, options=options)

def read_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'read_newick', input, options=options)

def read_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'read_nexus', input, options=options)

def read_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'read_paf', input, options=options)

def read_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'read_ped', input, options=options)

def read_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'read_phylip', input, options=options)

def read_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'read_twobit', input, options=options)

def read_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'read_zarr', input, options=options)

def resume_download(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'resume_download', input, options=options)

def sort(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'sort', input, options=options)

def sort_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'sort_map', input, options=options)

def sort_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'sort_newick', input, options=options)

def sort_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'sort_nexus', input, options=options)

def sort_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'sort_paf', input, options=options)

def sort_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'sort_ped', input, options=options)

def sort_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'sort_phylip', input, options=options)

def sort_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'sort_twobit', input, options=options)

def sort_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'sort_zarr', input, options=options)

def split_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'split_map', input, options=options)

def split_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'split_newick', input, options=options)

def split_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'split_nexus', input, options=options)

def split_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'split_paf', input, options=options)

def split_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'split_ped', input, options=options)

def split_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'split_phylip', input, options=options)

def split_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'split_twobit', input, options=options)

def split_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'split_zarr', input, options=options)

def stream(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'stream', input, options=options)

def stream_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'stream_map', input, options=options)

def stream_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'stream_newick', input, options=options)

def stream_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'stream_nexus', input, options=options)

def stream_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'stream_paf', input, options=options)

def stream_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'stream_ped', input, options=options)

def stream_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'stream_phylip', input, options=options)

def stream_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'stream_twobit', input, options=options)

def stream_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'stream_zarr', input, options=options)

def summarize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'summarize', input, options=options)

def summarize_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'summarize_map', input, options=options)

def summarize_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'summarize_newick', input, options=options)

def summarize_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'summarize_nexus', input, options=options)

def summarize_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'summarize_paf', input, options=options)

def summarize_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'summarize_ped', input, options=options)

def summarize_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'summarize_phylip', input, options=options)

def summarize_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'summarize_twobit', input, options=options)

def summarize_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'summarize_zarr', input, options=options)

def validate_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'validate_map', input, options=options)

def validate_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'validate_newick', input, options=options)

def validate_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'validate_nexus', input, options=options)

def validate_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'validate_paf', input, options=options)

def validate_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'validate_ped', input, options=options)

def validate_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'validate_phylip', input, options=options)

def validate_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'validate_twobit', input, options=options)

def validate_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'validate_zarr', input, options=options)

def write_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'write_map', input, options=options)

def write_newick(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'write_newick', input, options=options)

def write_nexus(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'write_nexus', input, options=options)

def write_paf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'write_paf', input, options=options)

def write_ped(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'write_ped', input, options=options)

def write_phylip(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'write_phylip', input, options=options)

def write_twobit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'write_twobit', input, options=options)

def write_zarr(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('io_formats', 'write_zarr', input, options=options)

__all__ += ['compress_map', 'compress_newick', 'compress_nexus', 'compress_paf', 'compress_ped', 'compress_phylip', 'compress_twobit', 'compress_zarr', 'convert', 'convert_map', 'convert_newick', 'convert_nexus', 'convert_paf', 'convert_ped', 'convert_phylip', 'convert_twobit', 'convert_zarr', 'decompress_map', 'decompress_newick', 'decompress_nexus', 'decompress_paf', 'decompress_ped', 'decompress_phylip', 'decompress_twobit', 'decompress_zarr', 'filter', 'filter_map', 'filter_newick', 'filter_nexus', 'filter_paf', 'filter_ped', 'filter_phylip', 'filter_twobit', 'filter_zarr', 'index', 'index_map', 'index_newick', 'index_nexus', 'index_paf', 'index_ped', 'index_phylip', 'index_twobit', 'index_zarr', 'merge', 'merge_map', 'merge_newick', 'merge_nexus', 'merge_paf', 'merge_ped', 'merge_phylip', 'merge_twobit', 'merge_zarr', 'read_map', 'read_newick', 'read_nexus', 'read_paf', 'read_ped', 'read_phylip', 'read_twobit', 'read_zarr', 'resume_download', 'sort', 'sort_map', 'sort_newick', 'sort_nexus', 'sort_paf', 'sort_ped', 'sort_phylip', 'sort_twobit', 'sort_zarr', 'split_map', 'split_newick', 'split_nexus', 'split_paf', 'split_ped', 'split_phylip', 'split_twobit', 'split_zarr', 'stream', 'stream_map', 'stream_newick', 'stream_nexus', 'stream_paf', 'stream_ped', 'stream_phylip', 'stream_twobit', 'stream_zarr', 'summarize', 'summarize_map', 'summarize_newick', 'summarize_nexus', 'summarize_paf', 'summarize_ped', 'summarize_phylip', 'summarize_twobit', 'summarize_zarr', 'validate_map', 'validate_newick', 'validate_nexus', 'validate_paf', 'validate_ped', 'validate_phylip', 'validate_twobit', 'validate_zarr', 'write_map', 'write_newick', 'write_nexus', 'write_paf', 'write_ped', 'write_phylip', 'write_twobit', 'write_zarr']

__candidate_adapter_operations__ = set(['compress_map', 'compress_newick', 'compress_nexus', 'compress_paf', 'compress_ped', 'compress_phylip', 'compress_twobit', 'compress_zarr', 'convert', 'convert_map', 'convert_newick', 'convert_nexus', 'convert_paf', 'convert_ped', 'convert_phylip', 'convert_twobit', 'convert_zarr', 'decompress_map', 'decompress_newick', 'decompress_nexus', 'decompress_paf', 'decompress_ped', 'decompress_phylip', 'decompress_twobit', 'decompress_zarr', 'filter', 'filter_map', 'filter_newick', 'filter_nexus', 'filter_paf', 'filter_ped', 'filter_phylip', 'filter_twobit', 'filter_zarr', 'index', 'index_map', 'index_newick', 'index_nexus', 'index_paf', 'index_ped', 'index_phylip', 'index_twobit', 'index_zarr', 'merge', 'merge_map', 'merge_newick', 'merge_nexus', 'merge_paf', 'merge_ped', 'merge_phylip', 'merge_twobit', 'merge_zarr', 'read_map', 'read_newick', 'read_nexus', 'read_paf', 'read_ped', 'read_phylip', 'read_twobit', 'read_zarr', 'resume_download', 'sort', 'sort_map', 'sort_newick', 'sort_nexus', 'sort_paf', 'sort_ped', 'sort_phylip', 'sort_twobit', 'sort_zarr', 'split_map', 'split_newick', 'split_nexus', 'split_paf', 'split_ped', 'split_phylip', 'split_twobit', 'split_zarr', 'stream', 'stream_map', 'stream_newick', 'stream_nexus', 'stream_paf', 'stream_ped', 'stream_phylip', 'stream_twobit', 'stream_zarr', 'summarize', 'summarize_map', 'summarize_newick', 'summarize_nexus', 'summarize_paf', 'summarize_ped', 'summarize_phylip', 'summarize_twobit', 'summarize_zarr', 'validate_map', 'validate_newick', 'validate_nexus', 'validate_paf', 'validate_ped', 'validate_phylip', 'validate_twobit', 'validate_zarr', 'write_map', 'write_newick', 'write_nexus', 'write_paf', 'write_ped', 'write_phylip', 'write_twobit', 'write_zarr'])


from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
