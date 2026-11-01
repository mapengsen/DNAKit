"""Phylogeny operations: native Newick parsing, distance trees (NJ/UPGMA), Robinson-Foulds and bootstrap, plus backend wrappers."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import candidate_not_implemented, candidate_ok
from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import ConfigurationError, InputFormatError

from . import _impl_phylo as _impl

from .external import infer_fast_tree, infer_ml_tree


_MODULE = "phylogeny"


def _ok(operation: str, input: object, result: object, options: Mapping[str, object] | None, warnings: tuple[str, ...] = ()) -> object:
    return candidate_ok(_MODULE, operation, input, result, options=options, warnings=warnings)


def _is_tree_text(value: object) -> bool:
    text = _impl.read_text(value)
    return text is not None and text.lstrip().startswith("(")


def _distance_kind(settings: Mapping[str, object]) -> str:
    kind = _impl.get_choice(settings, "distance", "jc69", ("jc69", "jc", "p"))
    return "p" if kind == "p" else "jc69"


def _saturation(settings: Mapping[str, object]) -> float | None:
    value = settings.get("saturation_distance")
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or value <= 0:
        raise ConfigurationError("options.saturation_distance must be a positive number.", code="PHYLOGENY_OPTION")
    return float(value)


def _matrix_result(names: list[str], matrix: list[list[float]], kind: str) -> dict[str, object]:
    return {"names": names, "matrix": matrix, "distance": kind}


def _tree_result(root: _impl.Node, method: str, names: list[str], matrix: list[list[float]]) -> dict[str, object]:
    stats = _impl.tree_statistics(root)
    return {
        "newick": _impl.to_newick(root),
        "method": method,
        "rooted": method == "upgma",
        "leaf_count": len(names),
        "total_branch_length": stats["total_branch_length"],
        "negative_branch_count": stats["negative_branch_count"],
        "note": "UPGMA gives a rooted ultrametric tree; neighbour-joining gives an unrooted tree drawn with a trifurcating root and may contain negative branch lengths",
    }


def annotate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Annotate every node of a Newick tree with structural fields (depth, leaf count, root distance, clade) and optional leaf metadata.

    input: Newick text or file path. Unnamed internal nodes are labelled N<preorder index> in the returned
    annotated_newick. options: metadata (optional {leaf name: value}; each internal node then lists the distinct
    values below it and whether they are homogeneous); max_clade (int, default 50; clades with more leaves are
    not listed). root_distance is the sum of branch lengths from the root, None if a length is missing.
    """

    settings = _impl.options_dict(options)
    metadata = settings.get("metadata")
    if metadata is not None and not isinstance(metadata, Mapping):
        raise ConfigurationError("options.metadata must map leaf names to values.", code="PHYLOGENY_OPTION")
    max_clade = _impl.get_int(settings, "max_clade", 50, low=1, high=1_000_000)
    root = _impl.load_tree(input)
    return _ok("annotate", input, _impl.annotate_tree(root, metadata, max_clade), options)


def build(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build a tree from a distance matrix by neighbour-joining (Saitou-Nei) or UPGMA (average linkage).

    input: {'names': [...], 'matrix': [[...]]} (symmetric, zero diagonal, non-negative; at most 400 taxa) or a
    nested {a: {b: distance}} mapping. options: method ('nj' default or 'upgma'). Ties are broken by input order.
    UPGMA returns a rooted ultrametric tree; NJ returns an unrooted tree written with a trifurcating root and can
    produce negative branch lengths (they are not clipped). Returns the Newick string.
    """

    settings = _impl.options_dict(options)
    method = _impl.get_choice(settings, "method", "nj", ("nj", "upgma"))
    names, matrix = _impl.parse_distance_matrix(input)
    root = _impl.build_tree(names, matrix, method)
    return _ok("build", input, _tree_result(root, method, names, matrix), options)


def calculate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Calculate tree statistics from a Newick tree, or a pairwise distance matrix from an aligned set of sequences.

    options.what: 'auto' (default: tree_statistics for a Newick tree, distance_matrix for an alignment),
    'tree_statistics', 'patristic_matrix' (tree) or 'distance_matrix' (alignment).
    Tree input: Newick text/path. tree_statistics reports leaf/internal counts, polytomies, depth in edges, Sackin
    index, total branch length, root-to-tip distances, ultrametricity, treeness and the Colless index (binary trees).
    Alignment input: FASTA text/path, {name: sequence} or [{'name','sequence'}]. options.distance 'jc69' (default,
    d = -3/4 ln(1 - 4p/3)) or 'p' (mismatches / compared sites; gaps and non-ACGT skipped pairwise); saturated pairs
    (p >= 0.75) raise unless options.saturation_distance is set.
    """

    settings = _impl.options_dict(options)
    what = _impl.get_choice(settings, "what", "auto", ("auto", "tree_statistics", "patristic_matrix", "distance_matrix"))
    if what == "auto":
        what = "tree_statistics" if _is_tree_text(input) else "distance_matrix"
    if what == "distance_matrix":
        kind = _distance_kind(settings)
        names, seqs = _impl.parse_alignment(input)
        matrix, compared = _impl.distance_matrix_from_alignment(names, seqs, kind, _saturation(settings))
        result = {"what": what, **_matrix_result(names, matrix, kind), "compared_sites": compared, "alignment_length": len(seqs[0])}
    else:
        root = _impl.load_tree(input)
        if what == "tree_statistics":
            result = {"what": what, **_impl.tree_statistics(root)}
        else:
            names, matrix = _impl.patristic_matrix(root)
            result = {"what": what, "names": names, "matrix": matrix, "distance": "patristic"}
    return _ok("calculate", input, result, options)


def _tree_pair(input: object, settings: Mapping[str, object]) -> tuple[object, object]:
    if settings.get("other") is not None:
        return input, settings["other"]
    if isinstance(input, Mapping):
        for left, right in (("a", "b"), ("left", "right"), ("tree1", "tree2")):
            if left in input and right in input:
                return input[left], input[right]
    if isinstance(input, (list, tuple)) and len(input) == 2:
        return input[0], input[1]
    raise InputFormatError("compare needs two trees.", code="TREE_COMPARE_INPUT", hint="Pass [newick_a, newick_b], {'a': .., 'b': ..} or options={'other': newick_b}.")


def compare(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compute the Robinson-Foulds distance between two trees on the same leaf set.

    input: [newick_a, newick_b], {'a': .., 'b': ..} or one Newick plus options.other (text or paths). RF is the
    number of non-trivial splits found in only one tree; unrooted splits by default. options: rooted (bool,
    default False; True compares rooted clades). Also returns normalized RF = RF / (splits_a + splits_b), shared
    splits and the maximum RF of binary trees. Different leaf sets raise an error. Branch lengths are ignored.
    """

    settings = _impl.options_dict(options)
    first, second = _tree_pair(input, settings)
    result = _impl.robinson_foulds(_impl.load_tree(first, "first tree"), _impl.load_tree(second, "second tree"), bool(settings.get("rooted", False)))
    return _ok("compare", input, result, options)


def export(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Serialise a Newick tree to Newick (normalised), NEXUS, nested JSON or an edge list.

    input: Newick text or path. options: format ('newick' default, 'nexus', 'json', 'edgelist'); output
    (optional file path; the text or JSON is written there; JSON is dumped as text). Newick output re-quotes
    labels when needed and writes branch lengths with 12 significant digits.
    """

    import json

    settings = _impl.options_dict(options)
    fmt = _impl.get_choice(settings, "format", "newick", ("newick", "nexus", "json", "edgelist"))
    root = _impl.load_tree(input)
    if fmt == "newick":
        payload: object = _impl.to_newick(root)
        text = str(payload)
    elif fmt == "nexus":
        payload = _impl.to_nexus(root)
        text = str(payload)
    elif fmt == "json":
        payload = _impl.tree_to_json(root)
        text = json.dumps(payload, indent=2)
    else:
        payload = _impl.edge_list(root)
        text = json.dumps(payload, indent=2)
    return _ok("export", input, {"format": fmt, "content": payload, "output_path": _impl.write_output(settings, text)}, options)


def infer(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Infer a distance-based tree (neighbour-joining or UPGMA on Jukes-Cantor or p distances) from aligned DNA sequences.

    input: FASTA text/path, {name: sequence} or [{'name','sequence'}]; equal lengths, at most 400 sequences. options:
    method ('nj' default or 'upgma'; maximum-likelihood or Bayesian methods are NOT available here - use
    infer_fast_tree / infer_ml_tree, which call external programs); distance ('jc69' default or 'p');
    saturation_distance (value used for saturated pairs instead of raising). Returns the Newick tree and the matrix.
    """

    settings = _impl.options_dict(options)
    requested = str(settings.get("method", "nj")).lower()
    if requested not in ("nj", "upgma"):
        return candidate_not_implemented(
            _MODULE, "infer", input, options=options, requires="an external phylogenetic program (FastTree, IQ-TREE, RAxML, MrBayes)",
            message=f"Method {requested!r} is not available natively; use infer_fast_tree or infer_ml_tree. Native methods: nj, upgma.",
        )
    kind = _distance_kind(settings)
    names, seqs = _impl.parse_alignment(input)
    matrix, _ = _impl.distance_matrix_from_alignment(names, seqs, kind, _saturation(settings))
    root = _impl.build_tree(names, matrix, requested)
    result = {**_tree_result(root, requested, names, matrix), "distance": kind, "distance_matrix": matrix, "names": names}
    return _ok("infer", input, result, options)


def plot(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Draw a Newick tree as a left-to-right rectangular SVG phylogram (or cladogram).

    input: Newick text or path (at most 2,000 leaves). options: branch_lengths (bool, default True: x position =
    root distance; False: x = depth in edges, requires no lengths); width (int pixels, 200-4000, default 640);
    output (optional path where the SVG is written). Leaves are drawn in tree order; no tree is reordered.
    """

    settings = _impl.options_dict(options)
    width = _impl.get_int(settings, "width", 640, low=200, high=4000)
    root = _impl.load_tree(input)
    use_lengths = bool(settings.get("branch_lengths", True))
    svg, leaf_count = _impl.tree_svg(root, use_lengths, width)
    return _ok("plot", input, {"svg": svg, "leaf_count": leaf_count, "branch_lengths": use_lengths, "output_path": _impl.write_output(settings, svg)}, options)


def test(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Estimate nonparametric bootstrap support for the clades of a distance-based tree built from an alignment.

    input: aligned DNA sequences (FASTA text/path, {name: sequence}, [{'name','sequence'}]; at most 400 sequences).
    The reference tree is built with the chosen method on the original alignment; each replicate resamples
    alignment columns with replacement, rebuilds the tree and counts which reference splits recur (Felsenstein 1985).
    options: replicates (int 1-10000, default 100); seed (int, default 0; same seed gives the same result);
    method ('nj' default or 'upgma'); distance ('jc69' default or 'p'); saturation_distance (distance used for
    saturated/undefined pairs in the replicates, default 10.0 for jc69 and 1.0 for p). Support is a fraction in
    [0, 1]; the returned Newick carries percentages as internal node labels. Work is capped (pairs x patterns x replicates).
    """

    settings = _impl.options_dict(options)
    replicates = _impl.get_int(settings, "replicates", 100, low=1, high=10_000)
    seed = _impl.get_int(settings, "seed", 0, low=0, high=2**32 - 1)
    method = _impl.get_choice(settings, "method", "nj", ("nj", "upgma"))
    kind = _distance_kind(settings)
    saturation = _saturation(settings) or (10.0 if kind == "jc69" else 1.0)
    names, seqs = _impl.parse_alignment(input)
    outcome = _impl.bootstrap_support(names, seqs, method, kind, replicates, seed, saturation)
    support = {split: outcome["split_hits"][split] / replicates for split in outcome["ref_splits"]}
    reference = outcome["reference"]
    _impl.label_support(reference, support, rooted=(method == "upgma"))
    result = {
        "newick_with_support": _impl.to_newick(reference),
        "clades": outcome["clades"],
        "replicates": replicates,
        "seed": seed,
        "method": method,
        "distance": kind,
        "mean_support": sum(support.values()) / len(support) if support else None,
    }
    return _ok("test", input, result, options)


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


# Coverage-gap exports (see paper/write_paper/figure/fig2/task_coverage/coverage_gap/waves/CANONICAL_MODULES.md)
from .dendropy_gap_wrappers import (
    align_with_muscle,
    convert_format,
    fetch_entrez_records,
    format_to_fasta,
    format_to_newick,
    format_to_nexml,
    format_to_nexus,
    format_to_phylip,
    paup_estimate_model,
    paup_estimate_tree,
    paup_estimate_ultrametric_tree,
    paup_prune_taxa,
    paup_pscore_trees,
    run_dendropy_format,
    summarize_tree_samples,
    transfer_tree_labels,
)
from .toytree_gap_functions import (
    draw_tree_cloud_figure,
    draw_tree_figure,
    draw_tree_grid_figure,
    draw_tree_node_feature_figure,
    export_tree_node_data,
    map_node_feature_values,
    prune_tree_tips,
    root_tree_on_outgroup as root_toytree_on_outgroup,
)

__all__ = [*__all__, "align_with_muscle", "convert_format", "fetch_entrez_records", "format_to_fasta", "format_to_newick", "format_to_nexml", "format_to_nexus", "format_to_phylip", "paup_estimate_model", "paup_estimate_tree", "paup_estimate_ultrametric_tree", "paup_prune_taxa", "paup_pscore_trees", "run_dendropy_format", "summarize_tree_samples", "transfer_tree_labels", "draw_tree_cloud_figure", "draw_tree_figure", "draw_tree_grid_figure", "draw_tree_node_feature_figure", "export_tree_node_data", "map_node_feature_values", "prune_tree_tips", "root_toytree_on_outgroup"]

# Public integration of the reviewed P00–P45 operation packs.
from .biopython_p18_tree_functions import (
    format_phylip_distance_matrix,
    write_tree_file,
)
from .biopython_p19_phylo_functions import (
    biopython_alignment_distance_matrix,
    codon_dn_ds_estimates,
    codon_dn_ds_trees,
    collapse_tree_genera,
    convert_tree_support,
    find_tree_clades,
    nexus_bootstrap_alignment,
    nexus_combine_matrices,
    nexus_terminal_gaps_to_missing,
    nexus_weighted_stepmatrix,
    nexus_write_partition_files,
    parsimony_tree_search,
    random_tree,
    read_paml_results,
    run_paml,
    split_tree_clade,
    strict_consensus_tree,
    tree_is_bifurcating,
    tree_is_parent_of,
    trees_are_compatible,
    trees_are_identical,
)
from .dendropy_cli_functions import (
    dendropy_convert_format,
)
from .dendropy_p38_functions import (
    calculate_folded_site_frequency_spectrum,
    calculate_hky85_rate_matrix,
    calculate_population_pair_statistics,
    concatenate_alignments,
    extract_alignment_columns,
    genbank_xml_to_alignment,
    parse_genbank_xml_records,
    read_multi_phylip_alignments,
    sample_stationary_nucleotide_sequence,
    simulate_fragmented_population_sequences,
)
from .iqtree_p11_cli_functions import (
    alignment_site_statistics,
    compare_tree_topologies,
    compute_branch_support,
    compute_site_concordance_factors,
    date_phylogeny,
    evaluate_tree_likelihood,
    infer_pathogen_tree,
    infer_polymorphism_aware_tree,
    likelihood_mapping,
    reconstruct_ancestral_sequences,
    run_symmetry_tests,
    select_substitution_model,
    simulate_alignment,
)
from .p16_fasttree_cli_functions import (
    compute_fasttree_distance_matrix,
    infer_fasttree_tree,
)
from .pyvolve_p41_functions import (
    pyvolve_state_frequencies,
    simulate_pyvolve_sequences,
)
from .skbio_p20_tree_functions import (
    build_taxonomy_tree,
    build_tree_from_taxdump,
    export_tree_lineages,
    extract_tree_supports,
    propagate_tip_counts,
    reroot_tree,
    shuffle_tree_tips,
    tree_diameter,
    tree_from_linkage_matrix,
    tree_tip_subsets,
)
from .toytree_p39_functions import (
    toytree_consensus,
    toytree_k80_matrix,
    toytree_local_genealogy,
    toytree_monophyly,
    toytree_mrca,
    toytree_neighbor_joining,
    toytree_rf_distance,
    toytree_tn93_matrix,
    toytree_unique_topologies,
    toytree_upgma,
)

__all__ = [
    *__all__,
    "alignment_site_statistics",
    "biopython_alignment_distance_matrix",
    "build_taxonomy_tree",
    "build_tree_from_taxdump",
    "calculate_folded_site_frequency_spectrum",
    "calculate_hky85_rate_matrix",
    "calculate_population_pair_statistics",
    "codon_dn_ds_estimates",
    "codon_dn_ds_trees",
    "collapse_tree_genera",
    "compare_tree_topologies",
    "compute_branch_support",
    "compute_fasttree_distance_matrix",
    "compute_site_concordance_factors",
    "concatenate_alignments",
    "convert_tree_support",
    "date_phylogeny",
    "dendropy_convert_format",
    "evaluate_tree_likelihood",
    "export_tree_lineages",
    "extract_alignment_columns",
    "extract_tree_supports",
    "find_tree_clades",
    "format_phylip_distance_matrix",
    "genbank_xml_to_alignment",
    "infer_fasttree_tree",
    "infer_pathogen_tree",
    "infer_polymorphism_aware_tree",
    "likelihood_mapping",
    "nexus_bootstrap_alignment",
    "nexus_combine_matrices",
    "nexus_terminal_gaps_to_missing",
    "nexus_weighted_stepmatrix",
    "nexus_write_partition_files",
    "parse_genbank_xml_records",
    "parsimony_tree_search",
    "propagate_tip_counts",
    "pyvolve_state_frequencies",
    "random_tree",
    "read_multi_phylip_alignments",
    "read_paml_results",
    "reconstruct_ancestral_sequences",
    "reroot_tree",
    "run_paml",
    "run_symmetry_tests",
    "sample_stationary_nucleotide_sequence",
    "select_substitution_model",
    "shuffle_tree_tips",
    "simulate_alignment",
    "simulate_fragmented_population_sequences",
    "simulate_pyvolve_sequences",
    "split_tree_clade",
    "strict_consensus_tree",
    "toytree_consensus",
    "toytree_k80_matrix",
    "toytree_local_genealogy",
    "toytree_monophyly",
    "toytree_mrca",
    "toytree_neighbor_joining",
    "toytree_rf_distance",
    "toytree_tn93_matrix",
    "toytree_unique_topologies",
    "toytree_upgma",
    "tree_diameter",
    "tree_from_linkage_matrix",
    "tree_is_bifurcating",
    "tree_is_parent_of",
    "tree_tip_subsets",
    "trees_are_compatible",
    "trees_are_identical",
    "write_tree_file",
]

_register_public_verbose(globals())
