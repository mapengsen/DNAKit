"""Unified, checksummed downloads from public biological databases."""

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.download.catalogs import expression, variants
from dnakit.download.ena import reads
from dnakit.download.ena import sequence as ena_sequence
from dnakit.download.encode import DEFAULT_ENCODE_DOWNLOAD_URL
from dnakit.download.encode import files as encode_files
from dnakit.download.files import ProgressCallback, dataset, download_file, download_url, tracks
from dnakit.download.indexes import IndexProgress, build_index
from dnakit.download.metadata import MetadataFormat, metadata
from dnakit.download.models import (
    DatasetDownloadResult,
    DownloadConfig,
    DownloadedFile,
    DownloadProgress,
    IndexArtifact,
    IndexBuildResult,
    RemoteFile,
)
from dnakit.download.ncbi import annotation, gene, genome, genome_package, taxonomy, virus_package
from dnakit.download.ncbi import sequence as ncbi_sequence
from dnakit.download.sequence import sequence
from dnakit.download.ucsc import DEFAULT_UCSC_DOWNLOAD_URL
from dnakit.download.ucsc import files as ucsc_files

__all__ = [
    "DEFAULT_ENCODE_DOWNLOAD_URL",
    "DEFAULT_UCSC_DOWNLOAD_URL",
    "DatasetDownloadResult",
    "DownloadConfig",
    "DownloadProgress",
    "DownloadedFile",
    "IndexArtifact",
    "IndexBuildResult",
    "IndexProgress",
    "MetadataFormat",
    "ProgressCallback",
    "RemoteFile",
    "annotation",
    "build_index",
    "dataset",
    "download_file",
    "download_url",
    "ena_sequence",
    "encode_files",
    "expression",
    "gene",
    "genome",
    "genome_package",
    "metadata",
    "ncbi_sequence",
    "reads",
    "sequence",
    "taxonomy",
    "tracks",
    "ucsc_files",
    "variants",
    "virus_package",
]

_register_public_verbose(globals())
