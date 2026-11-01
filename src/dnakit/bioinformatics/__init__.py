"""Optional bioinformatics APIs and reviewed public integration facades.

Third-party scientific packages are loaded only when an operation needs them.
Reviewed functions below are also discoverable through the Agent/MCP catalog.
See the bioinformatics API guide for dependencies.
"""

from dnakit.bioinformatics._common import ExecutionConfig, ToolResult

__all__ = ["ExecutionConfig", "ToolResult"]

# Public integration of the reviewed P00–P45 operation packs.
from .bioframe_p40_functions import (
    bioframe_fetch_ucsc_table,
    bioframe_is_complete_ucsc_string,
    bioframe_load_fasta,
)
from .cloning import (
    read_sanger,
)
from .pysam_gap_functions import (
    bcftools_call_consensus,
    samtools_ampliconstats,
    samtools_consensus,
    samtools_embedded_reference,
    samtools_flagstat,
    samtools_phase_report,
    samtools_stats,
)

__all__ = [
    *__all__,
    "bcftools_call_consensus",
    "bioframe_fetch_ucsc_table",
    "bioframe_is_complete_ucsc_string",
    "bioframe_load_fasta",
    "read_sanger",
    "samtools_ampliconstats",
    "samtools_consensus",
    "samtools_embedded_reference",
    "samtools_flagstat",
    "samtools_phase_report",
    "samtools_stats",
]
