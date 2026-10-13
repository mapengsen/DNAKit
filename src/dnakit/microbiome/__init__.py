"""Local candidate adapters for the microbiome capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation
from dnakit.microbiome.external import (
    classify_centrifuge,
    classify_gtdbtk,
    classify_kaiju,
    classify_kraken2,
    estimate_bracken_abundance,
    profile_humann,
)

def annotate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('microbiome', 'annotate', input, options=options)

def assemble(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('microbiome', 'assemble', input, options=options)

def bin(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('microbiome', 'bin', input, options=options)

def build_cooccurrence_network(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('microbiome', 'build_cooccurrence_network', input, options=options)

def build_microbial_gene_catalog(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('microbiome', 'build_microbial_gene_catalog', input, options=options)

def build_pathogen_transmission_tree(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('microbiome', 'build_pathogen_transmission_tree', input, options=options)

def calculate_ani(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('microbiome', 'calculate_ani', input, options=options)

def calculate_differential_abundance(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('microbiome', 'calculate_differential_abundance', input, options=options)

def calculate_digital_ddh(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('microbiome', 'calculate_digital_ddh', input, options=options)

def classify(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('microbiome', 'classify', input, options=options)

def compare(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('microbiome', 'compare', input, options=options)

def dereplicate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('microbiome', 'dereplicate', input, options=options)

def detect(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('microbiome', 'detect', input, options=options)

def export(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('microbiome', 'export', input, options=options)

def infer_viral_quasispecies(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('microbiome', 'infer_viral_quasispecies', input, options=options)

def plot(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('microbiome', 'plot', input, options=options)

def predict(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('microbiome', 'predict', input, options=options)

def qc(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('microbiome', 'qc', input, options=options)

def quantify(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('microbiome', 'quantify', input, options=options)

def read(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('microbiome', 'read', input, options=options)

def track(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('microbiome', 'track', input, options=options)

__all__ = ['annotate', 'assemble', 'bin', 'build_cooccurrence_network', 'build_microbial_gene_catalog', 'build_pathogen_transmission_tree', 'calculate_ani', 'calculate_differential_abundance', 'calculate_digital_ddh', 'classify', 'classify_centrifuge', 'classify_gtdbtk', 'classify_kaiju', 'classify_kraken2', 'compare', 'dereplicate', 'detect', 'estimate_bracken_abundance', 'export', 'infer_viral_quasispecies', 'plot', 'predict', 'profile_humann', 'qc', 'quantify', 'read', 'track']

__candidate_adapter_operations__ = set(['annotate', 'assemble', 'bin', 'build_cooccurrence_network', 'build_microbial_gene_catalog', 'build_pathogen_transmission_tree', 'calculate_ani', 'calculate_differential_abundance', 'calculate_digital_ddh', 'classify', 'compare', 'dereplicate', 'detect', 'export', 'infer_viral_quasispecies', 'plot', 'predict', 'qc', 'quantify', 'read', 'track'])


from dnakit._verbosity import register_public as _register_public_verbose
# Optional-backend exports (see docs/api/features) 
from .backends import (
    compare_minhash_sketches,
    estimate_unique_kmers,
    kmer_abundance_histogram,
    normalize_reads_by_median,
    query_kmer_counts,
    read_sourmash_signatures,
    sketch_fasta_minhash,
    write_sourmash_signatures,
)

__all__ = [*__all__, "compare_minhash_sketches", "estimate_unique_kmers", "kmer_abundance_histogram", "normalize_reads_by_median", "query_kmer_counts", "read_sourmash_signatures", "sketch_fasta_minhash", "write_sourmash_signatures"]

# Full-pass exports (see docs/api/features)
from .khmer_functions import (
    decode_kmer_2bit,
    encode_kmer_2bit,
    estimate_kmer_table_size,
    expected_kmer_table_collisions,
    hash_kmer_murmur3,
    inspect_khmer_graph,
    query_saved_countgraph,
    trim_low_abundance_reads,
)

from .sourmash_functions import (
    compare_minhash_matrix,
    downsample_minhash_sketch,
    estimate_mutation_rate,
    hash_canonical_kmer,
    intersect_minhash_sketches,
    lowest_common_lineage,
    search_minhash_database,
    set_size_error_probability,
)

__all__ = [*__all__, "decode_kmer_2bit", "encode_kmer_2bit", "estimate_kmer_table_size", "expected_kmer_table_collisions", "hash_kmer_murmur3", "inspect_khmer_graph", "query_saved_countgraph", "trim_low_abundance_reads", "compare_minhash_matrix", "downsample_minhash_sketch", "estimate_mutation_rate", "hash_canonical_kmer", "intersect_minhash_sketches", "lowest_common_lineage", "search_minhash_database", "set_size_error_probability"]

# Review-pass exports (see docs/api/features)
from .khmer_review_functions import (
    clean_khmer_read_sequences,
    read_paired_fastx,
)

from .sourmash_review_functions import (
    bp_threshold_to_hashes,
    count_kmer_lca_assignments,
    expected_mutated_kmers,
    kmer_mutation_probability,
    probability_no_shared_kmers,
    variance_mutated_kmers,
)

__all__ = [*__all__, "clean_khmer_read_sequences", "read_paired_fastx", "bp_threshold_to_hashes", "count_kmer_lca_assignments", "expected_mutated_kmers", "kmer_mutation_probability", "probability_no_shared_kmers", "variance_mutated_kmers"]

# Family-pass exports (see docs/api/features)
from .khmer_family_functions import (
    query_saved_khmer_table,
)

from .sourmash_family_functions import (
    compare_minhash_abundance_matrix,
    downsample_minhash_num_sketch,
    parse_sourmash_signature_json,
    search_minhash_database_ranked,
    search_sourmash_sbt_index,
    select_sourmash_signatures,
    write_sourmash_signature_archive,
)

__all__ = [*__all__, "query_saved_khmer_table", "compare_minhash_abundance_matrix", "downsample_minhash_num_sketch", "parse_sourmash_signature_json", "search_minhash_database_ranked", "search_sourmash_sbt_index", "select_sourmash_signatures", "write_sourmash_signature_archive"]

_register_public_verbose(globals())
