"""Local candidate adapters for the phylogeny capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation
from dnakit._verbosity import register_public as _register_public_verbose

from .external import infer_fast_tree, infer_ml_tree


def annotate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('phylogeny', 'annotate', input, options=options)

def build(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('phylogeny', 'build', input, options=options)

def calculate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('phylogeny', 'calculate', input, options=options)

def compare(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('phylogeny', 'compare', input, options=options)

def export(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('phylogeny', 'export', input, options=options)

def infer(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('phylogeny', 'infer', input, options=options)

def plot(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('phylogeny', 'plot', input, options=options)

def test(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('phylogeny', 'test', input, options=options)

__all__ = [
    'annotate', 'build', 'calculate', 'compare', 'export', 'infer',
    'infer_fast_tree', 'infer_ml_tree', 'plot', 'test'
]

__candidate_adapter_operations__ = set([
    'annotate', 'build', 'calculate', 'compare', 'export', 'infer', 'plot', 'test'
])


# Optional-backend exports (see docs/api/features) 
from .backends import (
    calculate_branch_length_distance,
    calculate_node_coordinates,
    calculate_patristic_distances,
    calculate_quartet_distance,
    read_tree_summary,
    root_tree_on_outgroup,
    simulate_nucleotide_alignment,
    simulate_nucleotide_sequences,
    simulate_substitution_counts,
)

__all__ = [*__all__, "calculate_branch_length_distance", "calculate_node_coordinates", "calculate_patristic_distances", "calculate_quartet_distance", "read_tree_summary", "root_tree_on_outgroup", "simulate_nucleotide_alignment", "simulate_nucleotide_sequences", "simulate_substitution_counts"]

# Full-pass exports (see docs/api/features)
from .dendropy_functions import (
    build_derived_state_matrix,
    calculate_b1_index,
    calculate_colless_tree_imbalance,
    calculate_mean_pairwise_differences,
    calculate_monophyletic_partition_discordance,
    calculate_nucleotide_diversity,
    calculate_parsimony_score,
    calculate_pybus_harvey_gamma,
    calculate_reconciliation_discordance,
    calculate_sackin_index,
    calculate_tajimas_d,
    calculate_treeness,
    calculate_unfolded_site_frequency_spectrum,
    calculate_wattersons_theta,
    calculate_weighted_robinson_foulds_distance,
    count_false_positive_and_negative_splits,
    count_segregating_sites,
)

from .toytree_functions import (
    calculate_equal_splits,
    calculate_information_rf_distance,
    calculate_matching_split_distance,
    calculate_matching_split_information_distance,
    calculate_mutual_clustering_information_distance,
    calculate_relative_evolutionary_divergence,
    calculate_shared_phylogenetic_information_distance,
    calculate_tree_correlation_matrix,
    calculate_tree_vcv_matrix,
    collapse_weak_nodes,
    extract_subtree,
    ladderize_tree,
    list_tree_bipartitions,
    remove_unary_nodes,
    root_tree_on_balanced_midpoint,
    root_tree_on_midpoint,
    root_tree_on_minimal_ancestor_deviation,
    scale_tree_to_root_height,
    unroot_tree,
)

__all__ = [*__all__, "build_derived_state_matrix", "calculate_b1_index", "calculate_colless_tree_imbalance", "calculate_mean_pairwise_differences", "calculate_monophyletic_partition_discordance", "calculate_nucleotide_diversity", "calculate_parsimony_score", "calculate_pybus_harvey_gamma", "calculate_reconciliation_discordance", "calculate_sackin_index", "calculate_tajimas_d", "calculate_treeness", "calculate_unfolded_site_frequency_spectrum", "calculate_wattersons_theta", "calculate_weighted_robinson_foulds_distance", "count_false_positive_and_negative_splits", "count_segregating_sites", "calculate_equal_splits", "calculate_information_rf_distance", "calculate_matching_split_distance", "calculate_matching_split_information_distance", "calculate_mutual_clustering_information_distance", "calculate_relative_evolutionary_divergence", "calculate_shared_phylogenetic_information_distance", "calculate_tree_correlation_matrix", "calculate_tree_vcv_matrix", "collapse_weak_nodes", "extract_subtree", "ladderize_tree", "list_tree_bipartitions", "remove_unary_nodes", "root_tree_on_balanced_midpoint", "root_tree_on_midpoint", "root_tree_on_minimal_ancestor_deviation", "scale_tree_to_root_height", "unroot_tree"]

# Review-pass exports (see docs/api/features)
from .dendropy_review_functions import (
    calculate_coalescence_ages,
    calculate_divergence_times,
    calculate_n_bar_index,
    calculate_node_ages,
    calculate_node_depths,
    list_missing_bipartitions,
)

from .toytree_review_functions import (
    add_tree_child_node,
    add_tree_internal_node,
    add_tree_internal_node_and_child,
    add_tree_internal_node_and_subtree,
    add_tree_sister_node,
    bisect_tree,
    calculate_all_node_distance_matrix,
    calculate_descendant_distances,
    calculate_farthest_node_distance,
    calculate_internal_node_distance_matrix,
    calculate_tip_distance_matrix_from_vcv,
    calculate_tip_diversification_rates,
    calculate_tree_node_distance,
    drop_tree_tips,
    extend_tree_tips_to_align,
    find_tree_farthest_node,
    list_tree_node_path,
    list_tree_quadripartitions,
    list_tree_quartets,
    map_consensus_tree_features,
    merge_tree_nodes,
    parse_tree_network,
    remove_tree_nodes,
    resolve_tree_polytomies,
    rotate_tree_node,
    set_tree_node_heights,
)

__all__ = [*__all__, "calculate_coalescence_ages", "calculate_divergence_times", "calculate_n_bar_index", "calculate_node_ages", "calculate_node_depths", "list_missing_bipartitions", "add_tree_child_node", "add_tree_internal_node", "add_tree_internal_node_and_child", "add_tree_internal_node_and_subtree", "add_tree_sister_node", "bisect_tree", "calculate_all_node_distance_matrix", "calculate_descendant_distances", "calculate_farthest_node_distance", "calculate_internal_node_distance_matrix", "calculate_tip_distance_matrix_from_vcv", "calculate_tip_diversification_rates", "calculate_tree_node_distance", "drop_tree_tips", "extend_tree_tips_to_align", "find_tree_farthest_node", "list_tree_node_path", "list_tree_quadripartitions", "list_tree_quartets", "map_consensus_tree_features", "merge_tree_nodes", "parse_tree_network", "remove_tree_nodes", "resolve_tree_polytomies", "rotate_tree_node", "set_tree_node_heights"]

_register_public_verbose(globals())
