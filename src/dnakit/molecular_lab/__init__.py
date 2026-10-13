"""Local candidate adapters for the molecular lab capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation

def analyze(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('molecular_lab', 'analyze', input, options=options)

def assemble_sanger_reads(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('molecular_lab', 'assemble_sanger_reads', input, options=options)

def calculate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('molecular_lab', 'calculate', input, options=options)

def call_bases_from_sanger(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('molecular_lab', 'call_bases_from_sanger', input, options=options)

def export(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('molecular_lab', 'export', input, options=options)

def optimize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('molecular_lab', 'optimize', input, options=options)

def parse(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('molecular_lab', 'parse', input, options=options)

def plan(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('molecular_lab', 'plan', input, options=options)

def plot(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('molecular_lab', 'plot', input, options=options)

def validate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('molecular_lab', 'validate', input, options=options)

__all__ = ['analyze', 'assemble_sanger_reads', 'calculate', 'call_bases_from_sanger', 'export', 'optimize', 'parse', 'plan', 'plot', 'validate']

__candidate_adapter_operations__ = set(['analyze', 'assemble_sanger_reads', 'calculate', 'call_bases_from_sanger', 'export', 'optimize', 'parse', 'plan', 'plot', 'validate'])


from dnakit._verbosity import register_public as _register_public_verbose
# Optional-backend exports (see docs/api/features) 
from .backends import (
    build_feature_map,
    build_genbank_feature_map,
    build_sbol3_document,
    calculate_primer3_dimer,
    calculate_primer3_hairpin,
    calculate_primer3_tm,
    convert_snapgene_to_genbank,
    design_assembly_primers,
    design_pcr_primers,
    evaluate_dnachisel_specs,
    find_recombination_sites,
    optimize_dnachisel_record,
    read_dnacauldron_parts,
    read_genbank_genome_tracks,
    read_sbol3_sequences,
    read_snapgene_features,
    read_snapgene_sequence,
    render_feature_map_svg,
    render_genome_tracks_svg,
    reverse_translate_protein,
    validate_sbol3_document,
    write_dnacauldron_constructs,
)

__all__ = [*__all__, "build_feature_map", "build_genbank_feature_map", "build_sbol3_document", "calculate_primer3_dimer", "calculate_primer3_hairpin", "calculate_primer3_tm", "convert_snapgene_to_genbank", "design_assembly_primers", "design_pcr_primers", "evaluate_dnachisel_specs", "find_recombination_sites", "optimize_dnachisel_record", "read_dnacauldron_parts", "read_genbank_genome_tracks", "read_sbol3_sequences", "read_snapgene_features", "read_snapgene_sequence", "render_feature_map_svg", "render_genome_tracks_svg", "reverse_translate_protein", "validate_sbol3_document", "write_dnacauldron_constructs"]

# Full-pass exports (see docs/api/features)
from .dna_features_viewer_functions import (
    annotate_sequence_feature,
    graphical_translation,
)

from .dnacauldron_functions import (
    autoselect_golden_gate_enzyme,
    crop_dna_record_with_features,
    list_type2s_digestion_overhangs,
)

from .dnachisel_functions import (
    find_sequence_difference_segments,
    generate_constrained_dna_sequence,
    group_nearby_positions,
    list_common_restriction_enzymes,
)

from .primer3_functions import (
    calculate_primer3_end_stability,
    design_primer3_pairs,
    sanitize_primer3_sequence,
)

from .pydna_functions import (
    assemble_pydna_fragments,
    fetch_genbank_record,
    find_common_substrings,
    find_terminal_overlaps,
    parse_primer_sequences,
    primer_melting_temperature,
    primer_pair_hybridization,
    simulate_pcr_amplicon,
)

from .pygenomeviz_functions import (
    fetch_genbank_by_accession,
    summarize_genbank_records,
)

from .sbol3_functions import (
    check_sbol3_display_id,
    copy_sbol3_design,
    extract_sbol3_display_id,
    to_sbol3_display_id,
)

from .snapgene_reader_functions import (
    snapgene_seqrecord_summary,
)

__all__ = [*__all__, "annotate_sequence_feature", "graphical_translation", "autoselect_golden_gate_enzyme", "crop_dna_record_with_features", "list_type2s_digestion_overhangs", "find_sequence_difference_segments", "generate_constrained_dna_sequence", "group_nearby_positions", "list_common_restriction_enzymes", "calculate_primer3_end_stability", "design_primer3_pairs", "sanitize_primer3_sequence", "assemble_pydna_fragments", "fetch_genbank_record", "find_common_substrings", "find_terminal_overlaps", "parse_primer_sequences", "primer_melting_temperature", "primer_pair_hybridization", "simulate_pcr_amplicon", "fetch_genbank_by_accession", "summarize_genbank_records", "check_sbol3_display_id", "copy_sbol3_design", "extract_sbol3_display_id", "to_sbol3_display_id", "snapgene_seqrecord_summary"]

# Review-pass exports (see docs/api/features)
from .dnacauldron_review_functions import (
    assign_parts_to_backbones,
    clone_into_donor_vector,
    record_has_backbone,
    simulate_type2s_assembly,
    substitute_type2s_overhangs,
)

from .dnachisel_review_functions import (
    codon_backtranslation_table,
    dna_window_overlap,
    group_nearby_dna_segments,
    subdivide_dna_window,
)

from .primer3_review_functions import (
    validate_acgt_sequence,
)

from .pydna_review_functions import (
    annotate_primer_binding_sites,
    assemble_recombinase_fragments,
    design_user_cloning_primers,
    dna_molecules_equal,
    extract_protospacers,
    gateway_recombination_products,
    generate_random_orf,
    ligate_blunt_fragments,
    parse_sequence_text,
    rare_codon_positions,
    simulate_genome_edit,
    smallest_circular_rotation,
    suggest_pcr_program,
)

__all__ = [*__all__, "assign_parts_to_backbones", "clone_into_donor_vector", "record_has_backbone", "simulate_type2s_assembly", "substitute_type2s_overhangs", "codon_backtranslation_table", "dna_window_overlap", "group_nearby_dna_segments", "subdivide_dna_window", "validate_acgt_sequence", "annotate_primer_binding_sites", "assemble_recombinase_fragments", "design_user_cloning_primers", "dna_molecules_equal", "extract_protospacers", "gateway_recombination_products", "generate_random_orf", "ligate_blunt_fragments", "parse_sequence_text", "rare_codon_positions", "simulate_genome_edit", "smallest_circular_rotation", "suggest_pcr_program"]

# Family-pass exports (see docs/api/features)
from .dnacauldron_family_functions import (
    parse_dnacauldron_sequence_text,
    read_dnacauldron_records,
)

from .dnachisel_family_functions import (
    annotate_sequence_differences,
    find_sequence_difference_positions,
)

from .pydna_family_functions import (
    amplicon_melting_temperature,
    find_common_substrings_with_topology,
    pcr_annealing_temperature,
    primer_melting_temperature_breslauer,
    primer_melting_temperature_dbd,
    simulate_pcr_assembly,
)

__all__ = [*__all__, "parse_dnacauldron_sequence_text", "read_dnacauldron_records", "annotate_sequence_differences", "find_sequence_difference_positions", "amplicon_melting_temperature", "find_common_substrings_with_topology", "pcr_annealing_temperature", "primer_melting_temperature_breslauer", "primer_melting_temperature_dbd", "simulate_pcr_assembly"]

_register_public_verbose(globals())
