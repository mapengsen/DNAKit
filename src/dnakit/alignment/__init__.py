"""Deterministic bounded native pairwise DNA alignment."""

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.alignment._candidates import *  # noqa: F403
from dnakit.alignment._candidates import __all__ as _candidate_all
from dnakit.alignment.external import (
    align_sequences,
    build_blast_database,
    cluster_mmseqs_protein,
    consensus_spoa,
    intersect_bed,
    local_blast,
    map_long_reads,
    search_diamond_blastp,
    search_diamond_blastx,
    search_mmseqs_protein,
    sort_bed_bedops,
)
from dnakit.alignment.pairwise import AlignmentConfig, align_pairwise
from dnakit.alignment.results import AlignmentColumn, AlignmentResult

__all__ = [
    "AlignmentColumn",
    "AlignmentConfig",
    "AlignmentResult",
    "align_pairwise",
    "align_sequences",
    "build_blast_database",
    "cluster_mmseqs_protein",
    "consensus_spoa",
    "intersect_bed",
    "local_blast",
    "map_long_reads",
    "search_diamond_blastp",
    "search_diamond_blastx",
    "search_mmseqs_protein",
    "sort_bed_bedops",
    *_candidate_all,
]

# Optional-backend exports (see docs/api/features) 
from .backends import (
    align_dna_edlib,
    align_dna_parasail,
    build_dna_hmm,
    build_minimap2_index,
    calculate_edit_distance_matrix,
    fetch_minimap2_reference_region,
    locate_approximate_pattern,
    map_query_sequences_paf,
    read_hmm_profile,
    score_dna_parasail,
    score_dna_parasail_matrix,
    search_dna_hmm,
)

__all__ = [*__all__, "align_dna_edlib", "align_dna_parasail", "build_dna_hmm", "build_minimap2_index", "calculate_edit_distance_matrix", "fetch_minimap2_reference_region", "locate_approximate_pattern", "map_query_sequences_paf", "read_hmm_profile", "score_dna_parasail", "score_dna_parasail_matrix", "search_dna_hmm"]

# Full-pass exports (see docs/api/features)
from .edlib_functions import (
    edlib_nice_alignment,
)

from .parasail_functions import (
    parasail_nw_banded_score,
    parasail_sg_db_align,
    parasail_sg_de_align,
    parasail_sg_dx_align,
    parasail_sg_qb_align,
    parasail_sg_qb_db_align,
    parasail_sg_qb_de_align,
    parasail_sg_qe_align,
    parasail_sg_qe_db_align,
    parasail_sg_qe_de_align,
    parasail_sg_qx_align,
    parasail_ssw_align,
)

from .parasail_family_functions import (
    parasail_dp_variant,
    parasail_ssw_profile_align,
)

from .pyhmmer_functions import (
    pyhmmer_align_to_hmm,
)

__all__ = [*__all__, "edlib_nice_alignment", "parasail_nw_banded_score", "parasail_sg_db_align", "parasail_sg_de_align", "parasail_sg_dx_align", "parasail_sg_qb_align", "parasail_sg_qb_db_align", "parasail_sg_qb_de_align", "parasail_sg_qe_align", "parasail_sg_qe_db_align", "parasail_sg_qe_de_align", "parasail_sg_qx_align", "parasail_ssw_align", "pyhmmer_align_to_hmm"]
__all__ = [*__all__, "parasail_dp_variant", "parasail_ssw_profile_align"]

_register_public_verbose(globals())
