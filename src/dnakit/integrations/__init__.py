"""Reviewed external genomics analyses with optional, caller-installed backends."""

from dnakit._verbosity import register_public as _register_public_verbose

from .catalog import describe_analysis, list_analyses, search_projects
from .execution import plan_analysis, run_analysis
from .mcp import call_mcp_tool, call_mcp_workflow, list_mcp_servers, list_mcp_tools
from .pip_tools import (
    align_biotite,
    align_edlib,
    beta_diversity_scikit_bio,
    compare_sourmash_sequences,
    diff_geneblocks_sequences,
    evaluate_seqme,
    fetch_pyfaidx_region,
    intersect_pybedtools_bed,
    map_mappy_sequences,
    merge_pyranges_bed,
    optimize_with_dnachisel,
    overlap_polars_bio_bed,
    permanova_scikit_bio,
    query_gffutils_database,
    simulate_dnacauldron_assembly,
    summarize_dna_features_viewer,
    summarize_dnaio_fastx,
    summarize_gtars_regions,
    summarize_pybigwig_region,
    summarize_pyfastx,
    summarize_pyranges_bed,
    summarize_sbol3_document,
    summarize_scadnano_design,
    summarize_vcf_cyvcf2,
    summarize_vcf_vcfpy,
)
from .publication import (
    build_software_citation,
    describe_publication_record,
    list_publication_records,
    publication_readiness,
    publication_registry_summary,
)
from .publication_audit import (
    audit_publication_evidence,
    describe_publication_audit,
    publication_audit_summary,
)
from .third_party import (
    describe_third_party_tool,
    describe_unified_capability,
    integration_inventory,
    list_third_party_tools,
    list_unified_capabilities,
    run_unified_capability,
    third_party_coverage,
    third_party_route,
    unified_capability_summary,
)

__all__ = [
    "align_biotite",
    "align_edlib",
    "audit_publication_evidence",
    "beta_diversity_scikit_bio",
    "build_software_citation",
    "call_mcp_tool",
    "call_mcp_workflow",
    "compare_sourmash_sequences",
    "describe_analysis",
    "describe_publication_record",
    "describe_publication_audit",
    "describe_third_party_tool",
    "describe_unified_capability",
    "diff_geneblocks_sequences",
    "evaluate_seqme",
    "fetch_pyfaidx_region",
    "integration_inventory",
    "intersect_pybedtools_bed",
    "list_analyses",
    "list_publication_records",
    "list_mcp_servers",
    "list_mcp_tools",
    "list_third_party_tools",
    "list_unified_capabilities",
    "map_mappy_sequences",
    "merge_pyranges_bed",
    "optimize_with_dnachisel",
    "overlap_polars_bio_bed",
    "permanova_scikit_bio",
    "publication_readiness",
    "publication_audit_summary",
    "publication_registry_summary",
    "plan_analysis",
    "query_gffutils_database",
    "run_analysis",
    "run_unified_capability",
    "search_projects",
    "simulate_dnacauldron_assembly",
    "summarize_dna_features_viewer",
    "summarize_dnaio_fastx",
    "summarize_gtars_regions",
    "summarize_pybigwig_region",
    "summarize_pyfastx",
    "summarize_pyranges_bed",
    "summarize_sbol3_document",
    "summarize_scadnano_design",
    "summarize_vcf_cyvcf2",
    "summarize_vcf_vcfpy",
    "third_party_coverage",
    "third_party_route",
    "unified_capability_summary",
]

# Optional-backend exports (see docs/api/features) 
from .service_backends import (
    analyze_fasta_with_clawbio,
    compute_toil_sequence_stats,
    describe_smartapi_api,
    get_biomcp_gene,
    get_biomcp_variant,
    list_galaxy_tools,
    search_biomcp_genes,
    search_smartapi_registry,
    show_galaxy_tool,
)

__all__ = [*__all__, "analyze_fasta_with_clawbio", "compute_toil_sequence_stats", "describe_smartapi_api", "get_biomcp_gene", "get_biomcp_variant", "list_galaxy_tools", "search_biomcp_genes", "search_smartapi_registry", "show_galaxy_tool"]

_register_public_verbose(globals(), tuple(__all__))
