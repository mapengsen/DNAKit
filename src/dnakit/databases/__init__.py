"""DNAKit scientific database APIs with local, attributed provider implementations."""

from ._candidates import *  # noqa: F403
from ._candidates import __all__ as _candidate_all
from ._generated import *  # noqa: F403
from ._generated import __all__ as _functions
from ._routes import (
    cadd_get_scores,
    orphanet_get_disease_genes,
    pdb_get_structure,
    string_get_functional_enrichment,
)
from ._runtime import call, describe, execute, list_operations, search

__all__ = [
    "call",
    "describe",
    "search",
    "execute",
    "list_operations",
    *_functions,
    "pdb_get_structure",
    "orphanet_get_disease_genes",
    "cadd_get_scores",
    "string_get_functional_enrichment",
    *_candidate_all,
]

from dnakit._verbosity import register_public as _register_public_verbose

# Optional-backend exports (see docs/api/features) 
from .gene_backends import (
    apply_cdna_mutations,
    fetch_ensembl_gene_info,
    fetch_ensembl_sequences,
    fetch_gene_sequence,
    fetch_mygene_gene,
    fetch_ncbi_gene_by_id,
    fetch_ncbi_gene_by_symbol,
    fetch_ncbi_gene_orthologs,
    fetch_spliced_sequence,
    fetch_variant_annotations,
    format_hgvs_id,
    hgvs_ids_from_vcf,
    index_gene_models,
    install_genome_assembly,
    list_biomart_attributes,
    list_biomart_datasets,
    map_gene_identifiers,
    query_biomart_dataset,
    query_genes_at_locus,
    query_genes_by_name,
    query_mygene_genes,
    query_variant_annotations,
    read_genome_summary,
    search_ensembl_genes,
    write_bed_sequences,
)

__all__ = [*__all__, "apply_cdna_mutations", "fetch_ensembl_gene_info", "fetch_ensembl_sequences", "fetch_gene_sequence", "fetch_mygene_gene", "fetch_ncbi_gene_by_id", "fetch_ncbi_gene_by_symbol", "fetch_ncbi_gene_orthologs", "fetch_spliced_sequence", "fetch_variant_annotations", "format_hgvs_id", "hgvs_ids_from_vcf", "index_gene_models", "install_genome_assembly", "list_biomart_attributes", "list_biomart_datasets", "map_gene_identifiers", "query_biomart_dataset", "query_genes_at_locus", "query_genes_by_name", "query_mygene_genes", "query_variant_annotations", "read_genome_summary", "search_ensembl_genes", "write_bed_sequences"]

# Full-pass exports (see docs/api/features)
from .biothings_client_functions import (
    biothings_client_get_homologs,
)

from .genomepy_functions import (
    filter_annotation_contigs,
    filter_fasta_records,
    read_annotation_records,
    sample_random_regions,
    write_annotation_records,
    write_fasta_sizes,
    write_gap_bed,
)

from .gget_functions import (
    gget_archs4,
    gget_bgee,
    gget_blast,
    gget_blat,
    gget_cosmic_local,
    gget_enrichr,
    gget_gene_expression,
    gget_opentargets,
    gget_psi_block,
    gget_specificity,
)

from .mygene_functions import (
    mygene_getgenes,
)

from .myvariant_functions import (
    myvariant_querymany,
)

from .ncbi_datasets_functions import (
    fetch_ncbi_assemblies_by_bioproject,
    fetch_ncbi_assemblies_by_taxon,
    fetch_ncbi_assembly_metadata,
    select_ncbi_assembly_fields,
    select_ncbi_gene_fields,
)

from .pyensembl_functions import (
    ensembl_genome_fasta_url,
    ensembl_gtf_url,
    ensembl_reference_name,
    ensembl_sequence_fasta_url,
    find_nearest_genomic_locus,
    lookup_sequence_version_fallback,
    match_identifier_version,
    merge_genomic_intervals,
    normalize_chromosome_name,
    normalize_gene_aliases,
    normalize_genome_reference_name,
    normalize_strand_symbol,
    read_fasta_records,
    read_first_fasta_id,
    stable_ensembl_id,
)

__all__ = [*__all__, "biothings_client_get_homologs", "filter_annotation_contigs", "filter_fasta_records", "read_annotation_records", "sample_random_regions", "write_annotation_records", "write_fasta_sizes", "write_gap_bed", "gget_archs4", "gget_bgee", "gget_blast", "gget_blat", "gget_cosmic_local", "gget_enrichr", "gget_gene_expression", "gget_opentargets", "gget_psi_block", "gget_specificity", "mygene_getgenes", "myvariant_querymany", "fetch_ncbi_assemblies_by_bioproject", "fetch_ncbi_assemblies_by_taxon", "fetch_ncbi_assembly_metadata", "select_ncbi_assembly_fields", "select_ncbi_gene_fields", "ensembl_genome_fasta_url", "ensembl_gtf_url", "ensembl_reference_name", "ensembl_sequence_fasta_url", "find_nearest_genomic_locus", "lookup_sequence_version_fallback", "match_identifier_version", "merge_genomic_intervals", "normalize_chromosome_name", "normalize_gene_aliases", "normalize_genome_reference_name", "normalize_strand_symbol", "read_fasta_records", "read_first_fasta_id", "stable_ensembl_id"]

# Review-pass exports (see docs/api/features)
from .genomepy_review_functions import (
    map_genome_contig_names,
)

from .pyensembl_review_functions import (
    ensembl_dated_genome_fasta_mask,
    ensembl_dated_release_urls,
    ensembl_genome_fasta_source_options,
    ensembl_species_record,
)

__all__ = [*__all__, "map_genome_contig_names", "ensembl_dated_genome_fasta_mask", "ensembl_dated_release_urls", "ensembl_genome_fasta_source_options", "ensembl_species_record"]

# Family-pass exports (see docs/api/features)
from .pyensembl_family_functions import (
    check_ensembl_canonical_source,
    ensembl_gtf_filename,
    ensembl_sequence_fasta_filename,
)

__all__ = [*__all__, "check_ensembl_canonical_source", "ensembl_gtf_filename", "ensembl_sequence_fasta_filename"]


# Coverage-gap exports (see paper/write_paper/figure/fig2/task_coverage/coverage_gap/waves/CANONICAL_MODULES.md)
from .genomepy_gap_functions import (
    genomepy_bowtie2_index_command,
    genomepy_bwa_index_command,
    genomepy_download_annotation,
    genomepy_download_annotation_head,
    genomepy_download_assembly_report,
    genomepy_ensembl_genomes,
    genomepy_find_annotation_urls,
    genomepy_gencode_genomes,
    genomepy_gencode_releases,
    genomepy_gencode_to_ucsc_names,
    genomepy_generate_annotation,
    genomepy_genome_size,
    genomepy_gmap_index_command,
    genomepy_head_annotations,
    genomepy_hisat2_index_command,
    genomepy_list_available_genomes,
    genomepy_minimap2_index_command,
    genomepy_ncbi_genomes,
    genomepy_run_index_cmd,
    genomepy_search_genomes,
    genomepy_star_index_command,
    genomepy_ucsc_download_annotation,
    genomepy_ucsc_genomes,
)
from .genomepy_verify_functions import (
    genomepy_ucsc_assembly_accessions,
    genomepy_ucsc_sourcename_accessions,
    genomepy_ucsc_annotation_links,
    genomepy_ensembl_add_grch37,
    genomepy_gencode_add_grch37,
)
from .gget_gap_functions import (
    gget_alphainfluenza_download,
    gget_cbio_download,
    gget_cbio_plot,
    gget_diamond_align,
    gget_genbank_metadata,
    gget_muscle_align,
    gget_ref,
    gget_sars_cov2_download,
    gget_virus,
    gget_virus_metadata,
    gget_virus_metadata_chunked,
    gget_virus_sequences_by_accession,
)
from .gget_precondition_functions import (
    gget_cellxgene_query,
    gget_cosmic_download_reference,
    gget_cosmic_select_reference,
)

__all__ = [*__all__, "genomepy_bowtie2_index_command", "genomepy_bwa_index_command", "genomepy_download_annotation", "genomepy_download_annotation_head", "genomepy_download_assembly_report", "genomepy_ensembl_genomes", "genomepy_find_annotation_urls", "genomepy_gencode_genomes", "genomepy_gencode_releases", "genomepy_gencode_to_ucsc_names", "genomepy_generate_annotation", "genomepy_genome_size", "genomepy_gmap_index_command", "genomepy_head_annotations", "genomepy_hisat2_index_command", "genomepy_list_available_genomes", "genomepy_minimap2_index_command", "genomepy_ncbi_genomes", "genomepy_run_index_cmd", "genomepy_search_genomes", "genomepy_star_index_command", "genomepy_ucsc_download_annotation", "genomepy_ucsc_genomes", "genomepy_ucsc_assembly_accessions", "genomepy_ucsc_sourcename_accessions", "genomepy_ucsc_annotation_links", "genomepy_ensembl_add_grch37", "genomepy_gencode_add_grch37", "gget_alphainfluenza_download", "gget_cbio_download", "gget_cbio_plot", "gget_diamond_align", "gget_genbank_metadata", "gget_muscle_align", "gget_ref", "gget_sars_cov2_download", "gget_virus", "gget_virus_metadata", "gget_virus_metadata_chunked", "gget_virus_sequences_by_accession", "gget_cellxgene_query", "gget_cosmic_download_reference", "gget_cosmic_select_reference"]

# Public integration of the reviewed P00–P45 operation packs.
from .biopython_p19_biosql_functions import (
    biosql_get_record,
    biosql_get_subsequence,
    biosql_list_entries,
    biosql_load_records,
)
from .genomepy_p40_functions import (
    genomepy_annotation_attributes,
    genomepy_annotation_dictionary,
    genomepy_annotation_genes,
    genomepy_annotation_lengths,
    genomepy_download_blacklist,
    genomepy_import_local_annotation,
    genomepy_list_installed,
    genomepy_sanitize_annotation,
)
from .gget_p41_functions import (
    gget_download_reference,
    gget_isoform_structure_map,
)
from .ncbi_eutils_p43_functions import (
    ncbi_fetch_dna_records,
    ncbi_link_dna_records,
    ncbi_search_dna_database,
    ncbi_summarize_dna_records,
)
from .public_apis_p44_functions import (
    datasets_gene_id_download_summary,
    datasets_gene_id_expression_report,
    datasets_gene_id_links,
    datasets_gene_product_report,
    datasets_gene_taxon_annotation_chromosome_summary,
    datasets_gene_taxon_counts,
    datasets_gene_taxon_dataset_report,
    datasets_gene_taxon_product_report,
    datasets_genome_accession_annotation_report,
    datasets_genome_accession_annotation_report_download_summary,
    datasets_genome_accession_annotation_summary,
    datasets_genome_accession_check,
    datasets_genome_accession_download_summary,
    datasets_genome_accession_links,
    datasets_genome_accession_revision_history,
    datasets_genome_report,
    datasets_genome_sequence_accession_sequence_assemblies,
    datasets_genome_taxon_checkm_histogram,
    datasets_organelle_accession_download,
    datasets_organelle_accessions_dataset_report,
    datasets_organelle_taxon_dataset_report,
    datasets_protein_accession_download,
    datasets_virus_accession_annotation_report,
    datasets_virus_accession_check,
    datasets_virus_download_summary,
    datasets_virus_taxon_annotation_report,
    ensembl_ga4gh_beacon_query,
    ensembl_ga4gh_callsets_record,
    ensembl_ga4gh_callsets_search,
    ensembl_ga4gh_datasets_record,
    ensembl_ga4gh_datasets_search,
    ensembl_ga4gh_features_record,
    ensembl_ga4gh_features_search,
    ensembl_ga4gh_featuresets_record,
    ensembl_ga4gh_featuresets_search,
    ensembl_ga4gh_references_record,
    ensembl_ga4gh_references_search,
    ensembl_ga4gh_referencesets_record,
    ensembl_ga4gh_referencesets_search,
    ensembl_ga4gh_variantannotations_search,
    ensembl_ga4gh_variantannotationsets_record,
    ensembl_ga4gh_variantannotationsets_search,
    ensembl_ga4gh_variants_record,
    ensembl_ga4gh_variants_search,
    ensembl_ga4gh_variantsets_record,
    ensembl_ga4gh_variantsets_search,
    ensembl_info_analysis,
    ensembl_info_assembly,
    ensembl_info_biotypes,
    ensembl_info_biotypes_groups,
    ensembl_info_biotypes_name,
    ensembl_info_compara_methods,
    ensembl_info_compara_species_sets,
    ensembl_info_external_dbs,
    ensembl_info_genomes_accession,
    ensembl_info_genomes_assembly,
    ensembl_info_genomes_division,
    ensembl_info_genomes_taxonomy,
    ensembl_info_variation,
    ensembl_info_variation_consequence_types,
    ensembl_info_variation_populations,
    ensembl_info_variation_populations_list,
    ensembl_map_cds,
    ensembl_map_translation,
    ensembl_overlap_translation,
    ensembl_phenotype_accession,
    ensembl_phenotype_gene,
    ensembl_phenotype_term,
    ensembl_variant_recoder,
    ensembl_variation_publication,
    ensembl_vep_region as ensembl_vep_region_query,
)
from .public_apis_p45_functions import (
    encode_body_system_summary,
    encode_download_urls,
    encode_experiment_matrix,
    encode_object_metadata,
    encode_tabular_report,
    rcsb_nucleotide_motif_search,
    rcsb_nucleotide_sequence_search,
    ucsc_search_assembly,
    ucsc_track_schema,
)

__all__ = [
    *__all__,
    "biosql_get_record",
    "biosql_get_subsequence",
    "biosql_list_entries",
    "biosql_load_records",
    "datasets_gene_id_download_summary",
    "datasets_gene_id_expression_report",
    "datasets_gene_id_links",
    "datasets_gene_product_report",
    "datasets_gene_taxon_annotation_chromosome_summary",
    "datasets_gene_taxon_counts",
    "datasets_gene_taxon_dataset_report",
    "datasets_gene_taxon_product_report",
    "datasets_genome_accession_annotation_report",
    "datasets_genome_accession_annotation_report_download_summary",
    "datasets_genome_accession_annotation_summary",
    "datasets_genome_accession_check",
    "datasets_genome_accession_download_summary",
    "datasets_genome_accession_links",
    "datasets_genome_accession_revision_history",
    "datasets_genome_report",
    "datasets_genome_sequence_accession_sequence_assemblies",
    "datasets_genome_taxon_checkm_histogram",
    "datasets_organelle_accession_download",
    "datasets_organelle_accessions_dataset_report",
    "datasets_organelle_taxon_dataset_report",
    "datasets_protein_accession_download",
    "datasets_virus_accession_annotation_report",
    "datasets_virus_accession_check",
    "datasets_virus_download_summary",
    "datasets_virus_taxon_annotation_report",
    "encode_body_system_summary",
    "encode_download_urls",
    "encode_experiment_matrix",
    "encode_object_metadata",
    "encode_tabular_report",
    "ensembl_ga4gh_beacon_query",
    "ensembl_ga4gh_callsets_record",
    "ensembl_ga4gh_callsets_search",
    "ensembl_ga4gh_datasets_record",
    "ensembl_ga4gh_datasets_search",
    "ensembl_ga4gh_features_record",
    "ensembl_ga4gh_features_search",
    "ensembl_ga4gh_featuresets_record",
    "ensembl_ga4gh_featuresets_search",
    "ensembl_ga4gh_references_record",
    "ensembl_ga4gh_references_search",
    "ensembl_ga4gh_referencesets_record",
    "ensembl_ga4gh_referencesets_search",
    "ensembl_ga4gh_variantannotations_search",
    "ensembl_ga4gh_variantannotationsets_record",
    "ensembl_ga4gh_variantannotationsets_search",
    "ensembl_ga4gh_variants_record",
    "ensembl_ga4gh_variants_search",
    "ensembl_ga4gh_variantsets_record",
    "ensembl_ga4gh_variantsets_search",
    "ensembl_info_analysis",
    "ensembl_info_assembly",
    "ensembl_info_biotypes",
    "ensembl_info_biotypes_groups",
    "ensembl_info_biotypes_name",
    "ensembl_info_compara_methods",
    "ensembl_info_compara_species_sets",
    "ensembl_info_external_dbs",
    "ensembl_info_genomes_accession",
    "ensembl_info_genomes_assembly",
    "ensembl_info_genomes_division",
    "ensembl_info_genomes_taxonomy",
    "ensembl_info_variation",
    "ensembl_info_variation_consequence_types",
    "ensembl_info_variation_populations",
    "ensembl_info_variation_populations_list",
    "ensembl_map_cds",
    "ensembl_map_translation",
    "ensembl_overlap_translation",
    "ensembl_phenotype_accession",
    "ensembl_phenotype_gene",
    "ensembl_phenotype_term",
    "ensembl_variant_recoder",
    "ensembl_variation_publication",
    "ensembl_vep_region_query",
    "genomepy_annotation_attributes",
    "genomepy_annotation_dictionary",
    "genomepy_annotation_genes",
    "genomepy_annotation_lengths",
    "genomepy_download_blacklist",
    "genomepy_import_local_annotation",
    "genomepy_list_installed",
    "genomepy_sanitize_annotation",
    "gget_download_reference",
    "gget_isoform_structure_map",
    "ncbi_fetch_dna_records",
    "ncbi_link_dna_records",
    "ncbi_search_dna_database",
    "ncbi_summarize_dna_records",
    "rcsb_nucleotide_motif_search",
    "rcsb_nucleotide_sequence_search",
    "ucsc_search_assembly",
    "ucsc_track_schema",
]

_register_public_verbose(globals())
