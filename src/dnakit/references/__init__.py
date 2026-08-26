"""Reference-genome resolution and download adapters."""

from dnakit.references.models import DownloadProgress, GenomeAssembly, GenomeDownloadResult
from dnakit.references.ncbi import (
    DEFAULT_CHUNK_SIZE,
    DEFAULT_MAX_DOWNLOAD_BYTES,
    DEFAULT_NCBI_API_BASE_URL,
    DEFAULT_TIMEOUT,
    ProgressCallback,
    download_genome,
    resolve_genome_assembly,
    supported_genome_aliases,
)

__all__ = [
    "DEFAULT_CHUNK_SIZE",
    "DEFAULT_MAX_DOWNLOAD_BYTES",
    "DEFAULT_NCBI_API_BASE_URL",
    "DEFAULT_TIMEOUT",
    "DownloadProgress",
    "GenomeAssembly",
    "GenomeDownloadResult",
    "ProgressCallback",
    "download_genome",
    "resolve_genome_assembly",
    "supported_genome_aliases",
]
