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

_register_public_verbose(globals())
