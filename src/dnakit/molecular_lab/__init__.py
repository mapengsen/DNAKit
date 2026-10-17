"""Local candidate adapters for the molecular lab capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation

from . import _impl_lab as _impl

def analyze(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Analyze oligonucleotide sequences: composition, GC fraction, reverse complement, molecular weight, A260 extinction coefficient and melting temperature.

    input: a DNA string (A, C, G, T), FASTA text, a list of sequences / {'id', 'sequence'} mappings, or a mapping name -> sequence (up to 1000 sequences of 100 kb).
    options: tm_method ('nearest_neighbor' default = SantaLucia 1998 via dnakit.thermodynamics, 'wallace' for 2-13 nt, 'basic' = 64.9 + 41 (G+C-16.4)/N for >= 14 nt),
    sodium_mm (default 50), primer_nm (default 250).
    Result per sequence: length, base counts, gc_fraction, longest homopolymer, 3' GC clamp, molecular_weight_da (anhydrous, 5'-OH), extinction_coefficient_260nm,
    nmol_per_od260, ug_per_od260, tm_celsius.
    """
    return _impl.run_op('analyze', _impl.op_analyze, input, options)

def assemble_sanger_reads(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Assemble overlapping Sanger reads into contigs by greedy suffix/prefix overlap merging (both strands).

    input: FASTA text, a list of sequences / {'id', 'sequence'} mappings or a mapping name -> sequence (up to 100 reads of 20 kb; IUPAC letters allowed).
    options: min_overlap (default 20 bases), min_identity (default 0.95, ungapped), both_strands (default True: reads may be reverse complemented).
    The pair of contigs with the longest overlap is merged repeatedly; the consensus keeps the earlier contig's base in overlaps (no quality weighting, no gap
    handling), so reads should be quality-trimmed first. Result: contigs (sequence, reads with strand and offset) and singlets.
    """
    return _impl.run_op('assemble_sanger_reads', _impl.op_assemble_sanger, input, options)

def calculate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Perform a bench calculation selected by options['calculation']: dilution, molarity, copy_number, ligation, a260, tm, pcr_yield or gel_size.

    Parameters are read from the input mapping and from options (options win).
    dilution: three of c1, v1, c2, v2 (C1*V1 = C2*V2); returns the fourth plus the diluent volume.
    molarity: length (bp, or sequence for the exact weight), strand ('ds' 660 g/mol per bp, 'ss' 330 per nt), and one of mass_concentration_ng_per_ul or molar_concentration_nm.
    copy_number: mass_ng, length or sequence, strand; copies = mass / MW * 6.02214076e23.
    ligation: vector_bp, insert_bp, vector_ng, insert_to_vector_ratio (default 3), optional vector_conc_ng_per_ul / insert_conc_ng_per_ul for volumes.
    a260: a260, nucleic_acid ('dsDNA' 50, 'ssDNA' 33, 'RNA' 40 ng/uL per A260) or sequence (oligo extinction coefficient), dilution_factor, path_length_cm.
    tm: sequence, method ('nearest_neighbor' default via dnakit.thermodynamics, 'wallace', 'basic'), sodium_mm (50), primer_nm (250).
    pcr_yield: template_copies, cycles, efficiency (default 1.0): N0 * (1+E)^cycles.
    gel_size: ladder [{'bp', 'distance'}], distances: log-linear least-squares fit of size against migration distance.
    """
    return _impl.run_op('calculate', _impl.op_calculate, input, options)

def call_bases_from_sanger(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Call bases from Sanger chromatogram traces by picking the strongest channel at each peak, with IUPAC codes for mixed peaks.

    input: an ABIF .ab1 file path (channels DATA9-12 ordered by FWO_1, peaks from PLOC2) or {'traces': {'A': [...], 'C': [...], 'G': [...], 'T': [...]}, 'peak_locations': [...] (optional)}.
    options: mixed_peak_ratio (default 0.5: second/first channel height at or above this gives the two-base IUPAC code), min_peak_height (default 0.1 of the global maximum, used when peaks are detected),
    min_peak_spacing (default 3 samples).
    Result: sequence, per-peak secondary/primary ratio, mixed positions and, for .ab1 files, agreement with the instrument's own calls. No Phred-style quality values are estimated.
    """
    return _impl.run_op('call_bases_from_sanger', _impl.op_call_bases, input, options)

def export(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Serialise lab records (samples, primers, sequences) to CSV, TSV, JSON or FASTA text, optionally writing a file.

    input: a list of record mappings (or sequence strings), one record, or a mapping id -> sequence.
    options: format ('csv' default, 'tsv', 'json', 'fasta'), columns (field order; default = fields in order of appearance), line_width (FASTA, default 80),
    output_path (optional file; refuses to overwrite unless overwrite=True). FASTA needs a 'sequence' (or 'seq') field and uses 'id' or 'name' as the header.
    """
    return _impl.run_op('export', _impl.op_export, input, options)

def optimize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Optimize primer length to reach a target melting temperature.

    input: template DNA string; the primer is the first n bases (direction 'forward') or the reverse complement of the last n bases ('reverse').
    options: target_tm (default 60), min_length (18), max_length (30), direction ('forward' default), require_gc_clamp (default False: the 3' base must be G or C),
    tm_method ('nearest_neighbor' default via dnakit.thermodynamics with sodium_mm=50, primer_nm=250, or 'basic').
    The length whose Tm is closest to target_tm is chosen (ties: shorter). Secondary structure and specificity are not evaluated.
    """
    return _impl.run_op('optimize', _impl.op_optimize, input, options)

def parse(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Parse laboratory quantity strings ('50 ng/uL', '10 mM', '200 uL') into SI values, or expand plate well ranges ('A1:B3').

    input: a quantity string or list of them (kind 'quantity') or well text (kind 'wells').
    options: kind ('quantity' default or 'wells'), plate (96 default or 384) and order ('row' default or 'column') for wells.
    Units: g, L, mol, M with prefixes k, m, u (or micro sign), n, p, f; ratios mass/volume and amount/volume; 'X' (fold). Result: value, unit, dimension, si_value and si_unit (g, L, mol, g/L, mol/L).
    """
    return _impl.run_op('parse', _impl.op_parse, input, options)

def plan(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Plan a bench protocol with exact volumes: PCR master mix, serial dilution or ligation setup.

    Parameters are read from the input mapping and options (options win); options['protocol'] selects 'pcr' (default), 'serial_dilution' or 'ligation'.
    pcr: reaction_volume_ul (25), n_reactions (1), overage_fraction (0.1), template_volume_ul (1), components [{'name', 'stock', 'final', 'unit'}] (default 1x buffer from 5x, 0.2 mM dNTP,
    0.5 uM primers, 0.025 U/uL polymerase); volume = final / stock * reaction volume, water to volume, master mix scaled by n * (1 + overage).
    serial_dilution: start_concentration, dilution_factor (10), steps, final_volume_ul (100).
    ligation: vector_bp, insert_bp, vector_ng, insert_to_vector_ratio (3), optional concentrations, reaction_volume_ul (20).
    """
    return _impl.run_op('plan', _impl.op_plan, input, options)

def plot(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Draw a schematic virtual agarose gel as an SVG string.

    input: mapping lane name -> list of fragment sizes in bp (or a list of such lists), e.g. {'ladder': [1000, 500, 100], 'PCR': [650]}.
    options: min_bp (default 100), max_bp (default 10000), output_path (optional .svg; no overwrite unless overwrite=True).
    Migration is linear in log10(size) between min_bp and max_bp; brightness grows with fragment mass. It is a drawing, not a simulation of a particular gel.
    """
    return _impl.run_op('plot', _impl.op_plot, input, options)

def validate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Validate PCR primers against rule-of-thumb design criteria.

    input: a primer sequence, list/mapping of primers, or a pair {'forward': ..., 'reverse': ...}.
    options: min_length (18), max_length (30), min_gc (0.40), max_gc (0.60), min_tm (52), max_tm (65), max_homopolymer (4), max_tm_difference (pair, 5), tm_method ('nearest_neighbor' default), sodium_mm, primer_nm.
    Checks length, GC content, 3' G/C clamp (3' base G/C, at most 3 G/C in the last 5), homopolymer runs, dinucleotide repeats and Tm. Passing does not predict amplification;
    secondary structure and off-target binding are not assessed.
    """
    return _impl.run_op('validate', _impl.op_validate, input, options)

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
