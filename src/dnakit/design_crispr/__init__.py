"""Local candidate adapters for the design crispr capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation

def assemble_dna_fragments(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('design_crispr', 'assemble_dna_fragments', input, options=options)

def calculate_primer_tm(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('design_crispr', 'calculate_primer_tm', input, options=options)

def check_primer_dimer(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('design_crispr', 'check_primer_dimer', input, options=options)

def check_primer_secondary_structure(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('design_crispr', 'check_primer_secondary_structure', input, options=options)

def check_primer_specificity(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('design_crispr', 'check_primer_specificity', input, options=options)

def design(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('design_crispr', 'design', input, options=options)

def encode_data_into_dna(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('design_crispr', 'encode_data_into_dna', input, options=options)

def evaluate_primer_pair(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('design_crispr', 'evaluate_primer_pair', input, options=options)

def export(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('design_crispr', 'export', input, options=options)

def filter(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('design_crispr', 'filter', input, options=options)

def plot(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('design_crispr', 'plot', input, options=options)

def predict_editing_outcomes(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('design_crispr', 'predict_editing_outcomes', input, options=options)

def rank(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('design_crispr', 'rank', input, options=options)

def remove_restriction_sites(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('design_crispr', 'remove_restriction_sites', input, options=options)

def run_crispr_screen_design(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('design_crispr', 'run_crispr_screen_design', input, options=options)

def scan_grna_offtargets(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('design_crispr', 'scan_grna_offtargets', input, options=options)

def score(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('design_crispr', 'score', input, options=options)

def simulate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('design_crispr', 'simulate', input, options=options)

def validate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('design_crispr', 'validate', input, options=options)

__all__ = ['assemble_dna_fragments', 'calculate_primer_tm', 'check_primer_dimer', 'check_primer_secondary_structure', 'check_primer_specificity', 'design', 'encode_data_into_dna', 'evaluate_primer_pair', 'export', 'filter', 'plot', 'predict_editing_outcomes', 'rank', 'remove_restriction_sites', 'run_crispr_screen_design', 'scan_grna_offtargets', 'score', 'simulate', 'validate']

__candidate_adapter_operations__ = set(['assemble_dna_fragments', 'calculate_primer_tm', 'check_primer_dimer', 'check_primer_secondary_structure', 'check_primer_specificity', 'design', 'encode_data_into_dna', 'evaluate_primer_pair', 'export', 'filter', 'plot', 'predict_editing_outcomes', 'rank', 'remove_restriction_sites', 'run_crispr_screen_design', 'scan_grna_offtargets', 'score', 'simulate', 'validate'])


from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
