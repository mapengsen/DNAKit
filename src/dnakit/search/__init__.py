"""Unified, bounded adapters for public biological database queries."""

from dnakit.search.blast import (
    DEFAULT_NCBI_BLAST_URL,
    BlastJob,
    BlastProgress,
    BlastStatus,
    blast_results,
    blast_status,
    identify,
    novelty,
    submit_blast,
    wait_for_blast,
)
from dnakit.search.ena import (
    DEFAULT_ENA_PORTAL_URL,
    DEFAULT_READ_FIELDS,
    ena_search,
    reads,
)
from dnakit.search.ena import (
    project as ena_project,
)
from dnakit.search.ena import (
    sample as ena_sample,
)
from dnakit.search.encode import DEFAULT_ENCODE_URL, encode_search
from dnakit.search.ensembl import (
    DEFAULT_ENSEMBL_REST_URL,
    annotation,
    chromosome,
    comparative_alignment,
    genome_info,
    homology,
    id_convert,
    map_coordinates,
    nearby,
    regulation,
    sequence,
    sequence_by_id,
    transcripts,
    variant,
)
from dnakit.search.ensembl import (
    lookup as ensembl_lookup,
)
from dnakit.search.models import QueryProgress, QueryResult, SearchConfig
from dnakit.search.ncbi import (
    DEFAULT_NCBI_DATASETS_URL,
    DEFAULT_NCBI_EUTILS_URL,
    accession,
    assembly,
    clinical_variant,
    database_version,
    entrez,
    expression,
    gene,
    literature,
    ncbi_orthologs,
    project,
    sample,
    taxonomy,
    virus,
)
from dnakit.search.ucsc import DEFAULT_UCSC_API_URL
from dnakit.search.ucsc import chromosomes as ucsc_chromosomes
from dnakit.search.ucsc import files as ucsc_files
from dnakit.search.ucsc import genomes as ucsc_genomes
from dnakit.search.ucsc import sequence as ucsc_sequence
from dnakit.search.ucsc import track_data as ucsc_track_data
from dnakit.search.ucsc import tracks as ucsc_tracks

__all__ = [
    "DEFAULT_ENA_PORTAL_URL",
    "DEFAULT_ENCODE_URL",
    "DEFAULT_ENSEMBL_REST_URL",
    "DEFAULT_NCBI_BLAST_URL",
    "DEFAULT_NCBI_DATASETS_URL",
    "DEFAULT_NCBI_EUTILS_URL",
    "DEFAULT_READ_FIELDS",
    "DEFAULT_UCSC_API_URL",
    "BlastJob",
    "BlastProgress",
    "BlastStatus",
    "QueryProgress",
    "QueryResult",
    "SearchConfig",
    "accession",
    "annotation",
    "assembly",
    "blast_results",
    "blast_status",
    "chromosome",
    "clinical_variant",
    "comparative_alignment",
    "database_version",
    "ena_project",
    "ena_sample",
    "ena_search",
    "encode_search",
    "ensembl_lookup",
    "entrez",
    "expression",
    "gene",
    "genome_info",
    "homology",
    "id_convert",
    "identify",
    "literature",
    "map_coordinates",
    "ncbi_orthologs",
    "nearby",
    "novelty",
    "project",
    "reads",
    "regulation",
    "sample",
    "sequence",
    "sequence_by_id",
    "submit_blast",
    "taxonomy",
    "transcripts",
    "ucsc_chromosomes",
    "ucsc_files",
    "ucsc_genomes",
    "ucsc_sequence",
    "ucsc_track_data",
    "ucsc_tracks",
    "variant",
    "virus",
    "wait_for_blast",
]
