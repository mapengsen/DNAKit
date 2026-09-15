# Adapted from ClawBio/ClawBio @ 1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d.
# Original file: clawbio/common/__init__.py. See LICENSE and UPSTREAM.json at the package root.
# DNAKit changes: package namespace, isolated cache, lazy imports and native execution.
"""Lazy shared scientific workflow helpers; see UPSTREAM.json."""
from importlib import import_module
_EXPORTS = {'detect_format': ('dnakit.workflows.common.parsers', 'detect_format'), 'parse_genetic_file': ('dnakit.workflows.common.parsers', 'parse_genetic_file'), 'GenotypeRecord': ('dnakit.workflows.common.parsers', 'GenotypeRecord'), 'sha256_file': ('dnakit.workflows.common.checksums', 'sha256_file'), 'sha256_hex': ('dnakit.workflows.common.checksums', 'sha256_hex'), 'generate_report_header': ('dnakit.workflows.common.report', 'generate_report_header'), 'generate_report_footer': ('dnakit.workflows.common.report', 'generate_report_footer'), 'DISCLAIMER': ('dnakit.workflows.common.report', 'DISCLAIMER'), 'PatientProfile': ('dnakit.workflows.common.profile', 'PatientProfile'), 'HtmlReportBuilder': ('dnakit.workflows.common.html_report', 'HtmlReportBuilder'), 'write_html_report': ('dnakit.workflows.common.html_report', 'write_html_report'), 'compute_input_checksum': ('dnakit.workflows.common.scrna_io', 'compute_input_checksum'), 'detect_processed_input_reason': ('dnakit.workflows.common.scrna_io', 'detect_processed_input_reason'), 'load_count_adata': ('dnakit.workflows.common.scrna_io', 'load_count_adata'), 'load_10x_mtx_data': ('dnakit.workflows.common.scrna_io', 'load_10x_mtx_data'), 'resolve_input_source': ('dnakit.workflows.common.scrna_io', 'resolve_input_source'), 'write_checksums': ('dnakit.workflows.common.reproducibility', 'write_checksums'), 'write_environment_yml': ('dnakit.workflows.common.reproducibility', 'write_environment_yml'), 'write_commands_sh': ('dnakit.workflows.common.reproducibility', 'write_commands_sh'), 'write_conda_lock': ('dnakit.workflows.common.reproducibility', 'write_conda_lock'), 'SarekConfig': ('dnakit.workflows.common.sarek', 'SarekConfig'), 'SarekSample': ('dnakit.workflows.common.sarek', 'SarekSample'), 'SarekWrapper': ('dnakit.workflows.common.sarek', 'SarekWrapper'), 'build_samplesheet': ('dnakit.workflows.common.sarek', 'build_samplesheet'), 'QcConfig': ('dnakit.workflows.common.vcf_qc', 'QcConfig'), 'QcResult': ('dnakit.workflows.common.vcf_qc', 'QcResult'), 'VcfQC': ('dnakit.workflows.common.vcf_qc', 'VcfQC')}
__all__ = list(_EXPORTS)
def __getattr__(name):
    if name not in _EXPORTS:
        raise AttributeError(name)
    module, attribute = _EXPORTS[name]
    if not module.startswith("dnakit."):
        module = __name__ + "." + module
    return getattr(import_module(module), attribute)
