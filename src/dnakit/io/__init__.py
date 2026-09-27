"""Streaming DNA record input and output."""

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.io.annotations import (
    AGPComponent,
    AGPDocument,
    AGPEntry,
    AGPGap,
    AnnotationDocument,
    AnnotationEntry,
    AnnotationSource,
    AnnotationTarget,
    read_agp,
    read_bed,
    read_gff3,
    write_agp,
    write_bed,
    write_gff3,
)
from dnakit.io.api import ReadMode, read, read_one, read_set, write
from dnakit.io.artifacts import file_digest, load_arrays, load_json, save_arrays, save_json
from dnakit.io.config import FeatureWritePolicy, ReadConfig, WriteConfig
from dnakit.io.indexed import (
    FastaIndex,
    FastaIndexEntry,
    FastqIndex,
    FastqIndexEntry,
    build_fasta_index,
    build_fastq_index,
    iter_chunks,
    load_fasta_index,
    load_fastq_index,
)
from dnakit.io.inputs import read_inputs
from dnakit.io.metadata import (
    MetadataValidationIssue,
    MetadataValidationReport,
    filter_by_metadata,
    merge_metadata,
    select_metadata,
    validate_metadata,
    with_metadata,
)
from dnakit.io.parquet import (
    ParquetBackendStatus,
    parquet_backend_status,
    require_parquet_backend,
)
from dnakit.io.results import GeneratedID, WriteResult
from dnakit.io.rows import RowWriter, iter_rows, write_rows
from dnakit.io.source import RecordSource
from dnakit.io.tables import (
    ParquetCompression,
    TableColumnType,
    TableExportResult,
    TableReadResult,
    TableSchema,
    export_result,
    export_table,
    read_table,
)

__all__ = [
    "AGPComponent",
    "AGPDocument",
    "AGPEntry",
    "AGPGap",
    "AnnotationDocument",
    "AnnotationEntry",
    "AnnotationSource",
    "AnnotationTarget",
    "FastaIndex",
    "FastaIndexEntry",
    "FastqIndex",
    "FastqIndexEntry",
    "FeatureWritePolicy",
    "GeneratedID",
    "MetadataValidationIssue",
    "MetadataValidationReport",
    "ParquetBackendStatus",
    "ParquetCompression",
    "ReadConfig",
    "ReadMode",
    "RecordSource",
    "RowWriter",
    "TableColumnType",
    "TableExportResult",
    "TableReadResult",
    "TableSchema",
    "WriteConfig",
    "WriteResult",
    "build_fasta_index",
    "build_fastq_index",
    "export_result",
    "export_table",
    "file_digest",
    "filter_by_metadata",
    "iter_chunks",
    "iter_rows",
    "load_arrays",
    "load_fasta_index",
    "load_fastq_index",
    "load_json",
    "merge_metadata",
    "parquet_backend_status",
    "read",
    "read_agp",
    "read_bed",
    "read_gff3",
    "read_inputs",
    "read_one",
    "read_set",
    "read_table",
    "require_parquet_backend",
    "save_arrays",
    "save_json",
    "select_metadata",
    "validate_metadata",
    "with_metadata",
    "write",
    "write_agp",
    "write_bed",
    "write_gff3",
    "write_rows",
]

_register_public_verbose(globals())
