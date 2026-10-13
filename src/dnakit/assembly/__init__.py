"""Local candidate adapters for the assembly capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation
from dnakit._verbosity import register_public as _register_public_verbose

from .external import assemble_flye, assemble_hifiasm, assemble_megahit, assemble_shasta, assemble_unicycler, evaluate_merqury, run_busco
from .local import (
    calculate_l50 as _calculate_l50_native,
    calculate_n50 as _calculate_n50_native,
    calculate_ng50 as _calculate_ng50_native,
    calculate_nx as _calculate_nx_native,
    calculate_stats as _calculate_stats_native,
)

def add_genome_to_pangenome(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'add_genome_to_pangenome', input, options=options)

def build_chromosome_level_assembly(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'build_chromosome_level_assembly', input, options=options)

def build_command(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'build_command', input, options=options)

def build_graph_genome(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'build_graph_genome', input, options=options)

def build_haplotype_assembly(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'build_haplotype_assembly', input, options=options)

def build_pangenome(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'build_pangenome', input, options=options)

def calculate_assembly_completeness(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'calculate_assembly_completeness', input, options=options)

def calculate_busco_completeness(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'calculate_busco_completeness', input, options=options)

def calculate_pangenome_openness(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'calculate_pangenome_openness', input, options=options)

def check_input(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'check_input', input, options=options)

def compare_assemblies(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'compare_assemblies', input, options=options)

def compare_assembly_graphs(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'compare_assembly_graphs', input, options=options)

def detect_assembly_contamination(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'detect_assembly_contamination', input, options=options)

def detect_chimeric_contigs(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'detect_chimeric_contigs', input, options=options)

def detect_misassemblies(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'detect_misassemblies', input, options=options)

def estimate_heterozygosity(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'estimate_heterozygosity', input, options=options)

def estimate_repeats(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'estimate_repeats', input, options=options)

def estimate_size(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'estimate_size', input, options=options)

def extract_accessory_pangenome(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'extract_accessory_pangenome', input, options=options)

def extract_core_pangenome(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'extract_core_pangenome', input, options=options)

def identify_centromeres(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'identify_centromeres', input, options=options)

def identify_repeats_in_assembly(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'identify_repeats_in_assembly', input, options=options)

def identify_telomeres(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'identify_telomeres', input, options=options)

def parse_log(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'parse_log', input, options=options)

def parse_output(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'parse_output', input, options=options)

def phase_assembly(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'phase_assembly', input, options=options)

def plot_graph(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'plot_graph', input, options=options)

def polish(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'polish', input, options=options)

def purge_haplotigs(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'purge_haplotigs', input, options=options)

def remove_contamination(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'remove_contamination', input, options=options)

def run(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'run', input, options=options)

def scaffold(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('assembly', 'scaffold', input, options=options)


# Promote deterministic contiguity statistics to the formal local namespace;
# all other assembly operations remain dependency-checked candidate adapters.
calculate_l50 = _calculate_l50_native
calculate_n50 = _calculate_n50_native
calculate_ng50 = _calculate_ng50_native
calculate_nx = _calculate_nx_native
calculate_stats = _calculate_stats_native

__all__ = ['add_genome_to_pangenome', 'assemble_flye', 'assemble_hifiasm', 'assemble_megahit', 'assemble_shasta', 'assemble_unicycler', 'build_chromosome_level_assembly', 'build_command', 'build_graph_genome', 'build_haplotype_assembly', 'build_pangenome', 'calculate_assembly_completeness', 'calculate_busco_completeness', 'calculate_l50', 'calculate_n50', 'calculate_ng50', 'calculate_nx', 'calculate_pangenome_openness', 'calculate_stats', 'check_input', 'compare_assemblies', 'compare_assembly_graphs', 'detect_assembly_contamination', 'detect_chimeric_contigs', 'detect_misassemblies', 'estimate_heterozygosity', 'estimate_repeats', 'estimate_size', 'extract_accessory_pangenome', 'extract_core_pangenome', 'evaluate_merqury', 'identify_centromeres', 'identify_repeats_in_assembly', 'identify_telomeres', 'parse_log', 'parse_output', 'phase_assembly', 'plot_graph', 'polish', 'purge_haplotigs', 'remove_contamination', 'run', 'run_busco', 'scaffold']

__candidate_adapter_operations__ = set(['add_genome_to_pangenome', 'build_chromosome_level_assembly', 'build_command', 'build_graph_genome', 'build_haplotype_assembly', 'build_pangenome', 'calculate_assembly_completeness', 'calculate_busco_completeness', 'calculate_l50', 'calculate_n50', 'calculate_ng50', 'calculate_nx', 'calculate_pangenome_openness', 'calculate_stats', 'check_input', 'compare_assemblies', 'compare_assembly_graphs', 'detect_assembly_contamination', 'detect_chimeric_contigs', 'detect_misassemblies', 'estimate_heterozygosity', 'estimate_repeats', 'estimate_size', 'extract_accessory_pangenome', 'extract_core_pangenome', 'identify_centromeres', 'identify_repeats_in_assembly', 'identify_telomeres', 'parse_log', 'parse_output', 'phase_assembly', 'plot_graph', 'polish', 'purge_haplotigs', 'remove_contamination', 'run', 'run_busco', 'scaffold']) - {'calculate_l50', 'calculate_n50', 'calculate_ng50', 'calculate_nx', 'calculate_stats'}
_register_public_verbose(globals())
