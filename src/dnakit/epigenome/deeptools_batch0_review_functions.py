"""deepTools BAM mapping statistics for the batch-0 coverage-gap review.

``bam_mapping_statistics()`` calls deepTools' ``plotEnrichment.openBam(returnStats=True)`` on an
indexed BAM or CRAM file and returns the mapped and unmapped read totals together with the
per-reference counts. The statistics come from the file index for BAM; for CRAM deepTools reads
the file, using ``threads`` workers. The alignment handle that ``openBam`` also returns is closed
inside the call, so the function leaves no file open.

deepTools calls ``sys.exit`` for a missing file, a missing index or a file that is not an
alignment file; those exits are converted into ``ConfigurationError``. deepTools is imported only
when the function runs.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.epigenome.deeptools_functions import _integer, _path, _provenance
from dnakit.epigenome.deeptools_review_functions import _call, _import_deeptools

_MAX_THREADS = 256
_CONVENTION = (
    "deeptools.plotEnrichment.openBam(returnStats=True): mapped and unmapped totals from the "
    "index (bam.mapped, bam.unmapped) and per-reference [mapped, unmapped] from "
    "get_index_statistics(); the file must be indexed"
)


def bam_mapping_statistics(
    bam_path: str | Path,
    *,
    threads: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count mapped and unmapped reads of an indexed BAM or CRAM file, as deepTools' openBam does.

    ``threads`` is used only for CRAM, where deepTools counts the reads itself. The result lists the
    total mapped and unmapped reads, the per-reference counts in header order and the file format.
    """

    path = _path(bam_path, "bam_path")
    workers = _integer(threads, "threads", 1, _MAX_THREADS)
    module = _import_deeptools("deeptools.plotEnrichment")
    handle, mapped, unmapped, stats = _call(
        module.openBam, str(path), returnStats=True, nThreads=workers
    )
    try:
        is_cram = bool(handle.is_cram)
    finally:
        handle.close()
    references = [
        {"name": str(name), "mapped": int(counts[0]), "unmapped": int(counts[1])}
        for name, counts in stats.items()
    ]
    return {
        "status": "completed",
        "operation": "bam_mapping_statistics",
        "bam_path": str(path),
        "format": "CRAM" if is_cram else "BAM",
        "mapped_reads": int(mapped),
        "unmapped_reads": int(unmapped),
        "references": references,
        "provenance": _provenance("deeptools", _CONVENTION),
    }


__all__ = ["bam_mapping_statistics"]

_register_public_verbose(globals())
