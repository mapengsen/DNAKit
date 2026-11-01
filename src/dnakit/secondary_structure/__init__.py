"""Secondary-structure parsing and explicit NUPACK integration."""

from dnakit._verbosity import register_public as _register_public_verbose

from .dotbracket import (
    analyze_dot_bracket,
    ensemble_defect_from_probabilities,
    pair_probability_metrics,
    target_structure_probability,
)
from .nupack import NupackAdapter, probe_nupack
from .results import (
    AccessibilityWindow,
    ComplexConcentration,
    NupackComplexResult,
    NupackTubeResult,
    PairProbabilityResult,
    PredictedSecondaryStructure,
    SecondaryBasePair,
    SecondaryStem,
    SecondaryStructureSummary,
)

__all__ = [
    "AccessibilityWindow",
    "ComplexConcentration",
    "NupackAdapter",
    "NupackComplexResult",
    "NupackTubeResult",
    "PairProbabilityResult",
    "PredictedSecondaryStructure",
    "SecondaryBasePair",
    "SecondaryStem",
    "SecondaryStructureSummary",
    "analyze_dot_bracket",
    "ensemble_defect_from_probabilities",
    "pair_probability_metrics",
    "probe_nupack",
    "target_structure_probability",
]

# Public integration of the reviewed P00–P45 operation packs.
from .biotite_p30_structure_functions import (
    base_pairs_to_dot_bracket,
    nucleic_acid_base_pairs,
    nucleic_acid_base_stacking,
    nucleic_acid_dot_bracket,
    plot_nucleotide_secondary_structure,
    pseudoknot_order,
)
from .nupack_p39_functions import (
    nupack_complex_concentrations,
    nupack_design_complex,
    nupack_design_tubes,
    nupack_ensemble_size,
    nupack_evaluate_design,
    nupack_loop_energy,
    nupack_loop_stacking_energies,
    nupack_mfe_structures,
    nupack_pair_probabilities,
    nupack_partition_function,
    nupack_sample_structures,
    nupack_sequence_distance,
    nupack_structure_distance,
    nupack_structure_energy,
    nupack_suboptimal_structures,
    nupack_tube_analysis,
)
from .viennarna_p30_functions import (
    cofold_rna,
    duplex_rna,
    duplex_rna_alignments,
    find_rna_interactions,
    fold_rna,
    fold_rna_alignment,
    fold_rna_alignment_local,
    fold_rna_complex,
    fold_rna_local,
    fold_rna_pseudoknots,
    rna_structure_coordinates,
)

__all__ = [
    *__all__,
    "base_pairs_to_dot_bracket",
    "cofold_rna",
    "duplex_rna",
    "duplex_rna_alignments",
    "find_rna_interactions",
    "fold_rna",
    "fold_rna_alignment",
    "fold_rna_alignment_local",
    "fold_rna_complex",
    "fold_rna_local",
    "fold_rna_pseudoknots",
    "nucleic_acid_base_pairs",
    "nucleic_acid_base_stacking",
    "nucleic_acid_dot_bracket",
    "nupack_complex_concentrations",
    "nupack_design_complex",
    "nupack_design_tubes",
    "nupack_ensemble_size",
    "nupack_evaluate_design",
    "nupack_loop_energy",
    "nupack_loop_stacking_energies",
    "nupack_mfe_structures",
    "nupack_pair_probabilities",
    "nupack_partition_function",
    "nupack_sample_structures",
    "nupack_sequence_distance",
    "nupack_structure_distance",
    "nupack_structure_energy",
    "nupack_suboptimal_structures",
    "nupack_tube_analysis",
    "plot_nucleotide_secondary_structure",
    "pseudoknot_order",
    "rna_structure_coordinates",
]

_register_public_verbose(globals())
