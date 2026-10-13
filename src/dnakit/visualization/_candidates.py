"""Local candidate adapters for the visualization capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation

def export_figure(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('visualization', 'export_figure', input, options=options)

def generate_analysis_summary(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('visualization', 'generate_analysis_summary', input, options=options)

def generate_clinical_pdf_report(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('visualization', 'generate_clinical_pdf_report', input, options=options)

def generate_html_genome_report(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('visualization', 'generate_html_genome_report', input, options=options)

def generate_qc_html_report(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('visualization', 'generate_qc_html_report', input, options=options)

def save_figure(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('visualization', 'save_figure', input, options=options)

def plot(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    """Render a generic visualization request through the selected backend."""

    return run_candidate_operation('visualization', 'plot', input, options=options)

__all__ = ['export_figure', 'generate_analysis_summary', 'generate_clinical_pdf_report', 'generate_html_genome_report', 'generate_qc_html_report', 'plot']

from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
