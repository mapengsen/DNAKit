"""Local candidate adapters for the clinical capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation

def analyze_trio_variants(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('clinical', 'analyze_trio_variants', input, options=options)

def annotate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('clinical', 'annotate', input, options=options)

def build_cancer_evolution_tree(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('clinical', 'build_cancer_evolution_tree', input, options=options)

def calculate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('clinical', 'calculate', input, options=options)

def call(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('clinical', 'call', input, options=options)

def classify(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('clinical', 'classify', input, options=options)

def compare(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('clinical', 'compare', input, options=options)

def detect_circulating_tumor_dna(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('clinical', 'detect_circulating_tumor_dna', input, options=options)

def detect_de_novo_variants(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('clinical', 'detect_de_novo_variants', input, options=options)

def detect_minimal_residual_disease(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('clinical', 'detect_minimal_residual_disease', input, options=options)

def generate_clinical_variant_report(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('clinical', 'generate_clinical_variant_report', input, options=options)

def infer_tumor_clones(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('clinical', 'infer_tumor_clones', input, options=options)

def predict(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('clinical', 'predict', input, options=options)

def prioritize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('clinical', 'prioritize', input, options=options)

def qc(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('clinical', 'qc', input, options=options)

def read(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('clinical', 'read', input, options=options)

def report(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('clinical', 'report', input, options=options)

def type_hla(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('clinical', 'type_hla', input, options=options)

__all__ = ['analyze_trio_variants', 'annotate', 'build_cancer_evolution_tree', 'calculate', 'call', 'classify', 'compare', 'detect_circulating_tumor_dna', 'detect_de_novo_variants', 'detect_minimal_residual_disease', 'generate_clinical_variant_report', 'infer_tumor_clones', 'predict', 'prioritize', 'qc', 'read', 'report', 'type_hla']

__candidate_adapter_operations__ = set(['analyze_trio_variants', 'annotate', 'build_cancer_evolution_tree', 'calculate', 'call', 'classify', 'compare', 'detect_circulating_tumor_dna', 'detect_de_novo_variants', 'detect_minimal_residual_disease', 'generate_clinical_variant_report', 'infer_tumor_clones', 'predict', 'prioritize', 'qc', 'read', 'report', 'type_hla'])


from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
