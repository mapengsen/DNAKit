"""Local candidate adapters for the haplotypes capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation

def analyze_trio_inheritance(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('haplotypes', 'analyze_trio_inheritance', input, options=options)

def calculate_pedigree_likelihood(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('haplotypes', 'calculate_pedigree_likelihood', input, options=options)

def cluster(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('haplotypes', 'cluster', input, options=options)

def compare(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('haplotypes', 'compare', input, options=options)

def construct(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('haplotypes', 'construct', input, options=options)

def count(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('haplotypes', 'count', input, options=options)

def detect_compound_heterozygotes(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('haplotypes', 'detect_compound_heterozygotes', input, options=options)

def detect_denovo_variants(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('haplotypes', 'detect_denovo_variants', input, options=options)

def frequency(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('haplotypes', 'frequency', input, options=options)

def impute(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('haplotypes', 'impute', input, options=options)

def infer_parentage(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('haplotypes', 'infer_parentage', input, options=options)

def phase(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('haplotypes', 'phase', input, options=options)

def reconstruct(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('haplotypes', 'reconstruct', input, options=options)

def select(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('haplotypes', 'select', input, options=options)

def simulate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('haplotypes', 'simulate', input, options=options)

def tag(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('haplotypes', 'tag', input, options=options)

def visualize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('haplotypes', 'visualize', input, options=options)

__all__ = ['analyze_trio_inheritance', 'calculate_pedigree_likelihood', 'cluster', 'compare', 'construct', 'count', 'detect_compound_heterozygotes', 'detect_denovo_variants', 'frequency', 'impute', 'infer_parentage', 'phase', 'reconstruct', 'select', 'simulate', 'tag', 'visualize']

__candidate_adapter_operations__ = set(['analyze_trio_inheritance', 'calculate_pedigree_likelihood', 'cluster', 'compare', 'construct', 'count', 'detect_compound_heterozygotes', 'detect_denovo_variants', 'frequency', 'impute', 'infer_parentage', 'phase', 'reconstruct', 'select', 'simulate', 'tag', 'visualize'])


from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
