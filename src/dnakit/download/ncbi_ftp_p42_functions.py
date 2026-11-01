"""Curated NCBI FTP biological data collections via HTTPS, with real streamed downloads.

``collection`` selects one of the named datasets in ``NCBI_FTP_COLLECTIONS``.
Paths are relative to that dataset root; callers select release-specific filenames
from the listing. Archives are downloaded intact and are not silently extracted.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import PurePosixPath
from typing import Any
from urllib.parse import quote, unquote, urljoin, urlsplit

from dnakit.download.files import dataset
from dnakit.download.models import DownloadConfig, RemoteFile
from dnakit.exceptions import ConfigurationError

NCBI_FTP_COLLECTIONS = {
    "cds_nucleotide_sequences_fasta": {
        "root": "genomes/all/",
        "operation": "download CDS nucleotide sequences (FASTA)",
        "group": "G4011",
    },
    "rna_sequences_fasta": {
        "root": "genomes/all/",
        "operation": "download RNA sequences (FASTA)",
        "group": "G4012",
    },
    "repeatmasker_repeat_annotation_rm_out": {
        "root": "genomes/all/",
        "operation": "download RepeatMasker repeat annotation (rm.out)",
        "group": "G4013",
    },
    "genome_feature_tables": {
        "root": "genomes/all/",
        "operation": "download genome feature tables",
        "group": "G4014",
    },
    "genome_feature_count_summaries": {
        "root": "genomes/all/",
        "operation": "download genome feature count summaries",
        "group": "G4015",
    },
    "assembly_structure_and_region_descriptors_agp": {
        "root": "genomes/all/",
        "operation": "download assembly structure and region descriptors (AGP)",
        "group": "G4016",
    },
    "mapview_genome_map_data_sequence_and_feature_objects": {
        "root": "genomes/MapView/",
        "operation": "download MapView genome map data (sequence and feature objects)",
        "group": "G4017",
    },
    "influenza_nucleotide_sequences_and_annotation_tables": {
        "root": "genomes/INFLUENZA/",
        "operation": "download influenza nucleotide sequences and annotation tables",
        "group": "G4018",
    },
    "rrna_targeted_locus_sequences_16s_18s_23s_28s_5s_its": {
        "root": "refseq/TargetedLoci/",
        "operation": "download rRNA targeted-locus sequences (16S/18S/23S/28S/5S/ITS)",
        "group": "G4019",
    },
    "prokaryote_ani_and_checkm_genome_quality_reports": {
        "root": "genomes/ASSEMBLY_REPORTS/",
        "operation": "download prokaryote ANI and CheckM genome-quality reports",
        "group": "G4020",
    },
    "pairwise_genome_alignments_comparative_genome_viewer": {
        "root": "comparative-genome-viewer/",
        "operation": "download pairwise genome alignments (Comparative Genome Viewer)",
        "group": "G4021",
    },
    "genome_reference_consortium_assembly_sequences": {
        "root": "pub/grc/",
        "operation": "download Genome Reference Consortium assembly sequences",
        "group": "G4022",
    },
    "nist_pathogen_reference_genome_assemblies_fasta": {
        "root": "nist-immsa/NISTPathogen/",
        "operation": "download NIST pathogen reference genome assemblies (FASTA)",
        "group": "G4023",
    },
    "nucleotide_blast_database_archives": {
        "root": "blast/db/",
        "operation": "download nucleotide BLAST database archives",
        "group": "G4024",
    },
    "windowmasker_genome_masking_files": {
        "root": "blast/windowmasker_files/",
        "operation": "download WindowMasker genome masking files",
        "group": "G4025",
    },
    "genbank_nucleotide_release_flat_files": {
        "root": "genbank/",
        "operation": "download GenBank nucleotide release flat files",
        "group": "G4026",
    },
    "genbank_est_gss_and_other_record_catalogs_and_gene_lists": {
        "root": "genbank/catalog/",
        "operation": "list GenBank EST, GSS and other record catalogs and gene lists",
        "group": "G4027",
    },
    "third_party_annotation_tpa_records_in_genbank_flat_file_format": {
        "root": "tpa/release/gbff/",
        "operation": "download third-party annotation (TPA) records in GenBank flat-file format",
        "group": "G4028",
    },
    "third_party_annotation_tpa_nucleotide_sequences_fasta": {
        "root": "tpa/release/",
        "operation": "download third-party annotation (TPA) nucleotide sequences (FASTA)",
        "group": "G4029",
    },
    "refseq_wgs_master_contig_records": {
        "root": "refseq/wgs/",
        "operation": "download RefSeq WGS master contig records",
        "group": "G4030",
    },
    "refseq_functional_element_biological_region_sequences": {
        "root": "refseq/FunctionalElements/",
        "operation": "download RefSeq functional-element (biological region) sequences",
        "group": "G4031",
    },
    "mane_select_transcript_set": {
        "root": "refseq/MANE/MANE_human/",
        "operation": "download MANE Select transcript set",
        "group": "G4032",
    },
    "ncbi_gene_records_table_gene_info": {
        "root": "gene/DATA/",
        "operation": "download NCBI Gene records table (gene_info)",
        "group": "G4033",
    },
    "gene_to_nucleotide_accession_mapping_gene2accession_gene2refseq": {
        "root": "gene/DATA/",
        "operation": "download gene-to-nucleotide accession mapping (gene2accession, gene2refseq)",
        "group": "G4034",
    },
    "gene_to_ensembl_identifier_mapping_gene2ensembl": {
        "root": "gene/DATA/",
        "operation": "download gene-to-Ensembl identifier mapping (gene2ensembl)",
        "group": "G4035",
    },
    "gene_to_go_annotation_table_gene2go": {
        "root": "gene/DATA/",
        "operation": "download gene-to-GO annotation table (gene2go)",
        "group": "G4036",
    },
    "gene_orthology_table_gene_orthologs": {
        "root": "gene/DATA/",
        "operation": "download gene orthology table (gene_orthologs)",
        "group": "G4037",
    },
    "gene_genomic_neighbourhood_table_gene_neighbors": {
        "root": "gene/DATA/",
        "operation": "download gene genomic neighbourhood table (gene_neighbors)",
        "group": "G4038",
    },
    "gene_group_table_gene_group": {
        "root": "gene/DATA/",
        "operation": "download gene group table (gene_group)",
        "group": "G4039",
    },
    "discontinued_gene_identifier_history_gene_history": {
        "root": "gene/DATA/",
        "operation": "download discontinued gene identifier history (gene_history)",
        "group": "G4040",
    },
    "gene_summary_text_gene_summary": {
        "root": "gene/DATA/",
        "operation": "download gene summary text (gene_summary)",
        "group": "G4041",
    },
    "gene_refseq_uniprotkb_mapping_gene_refseq_uniprotkb_collab": {
        "root": "gene/DATA/",
        "operation": "download gene-RefSeq-UniProtKB mapping (gene_refseq_uniprotkb_collab)",
        "group": "G4042",
    },
    "ncbi_gene_expression_data_mammalia": {
        "root": "gene/DATA/expression/",
        "operation": "download NCBI gene expression data (Mammalia)",
        "group": "G4043",
    },
    "ccds_consensus_coding_sequence_records": {
        "root": "pub/CCDS/",
        "operation": "download CCDS consensus coding sequence records",
        "group": "G4044",
    },
    "homologene_gene_homology_group_tables": {
        "root": "pub/HomoloGene/",
        "operation": "download HomoloGene gene homology group tables",
        "group": "G4045",
    },
    "dbsnp_variant_calls_vcf": {
        "root": "snp/",
        "operation": "download dbSNP variant calls (VCF)",
        "group": "G4046",
    },
    "dbsnp_refsnp_records_json": {
        "root": "snp/latest_release/JSON/",
        "operation": "download dbSNP refSNP records (JSON)",
        "group": "G4047",
    },
    "dbsnp_coordinate_bed_tracks": {
        "root": "snp/organisms/human_9606/BED/",
        "operation": "download dbSNP coordinate BED tracks",
        "group": "G4048",
    },
    "dbsnp_rs_flanking_sequences_fasta": {
        "root": "snp/organisms/human_9606/rs_fasta/",
        "operation": "download dbSNP rs flanking sequences (FASTA)",
        "group": "G4049",
    },
    "dbsnp_genotype_tables": {
        "root": "snp/organisms/human_9606/",
        "operation": "download dbSNP genotype tables",
        "group": "G4050",
    },
    "dbsnp_population_allele_frequencies": {
        "root": "snp/population_frequency/",
        "operation": "download dbSNP population allele frequencies",
        "group": "G4051",
    },
    "dbsnp_submitted_variant_records": {
        "root": "snp/submission/",
        "operation": "download dbSNP submitted variant records",
        "group": "G4052",
    },
    "dbvar_structural_variant_records": {
        "root": "pub/dbVar/data/",
        "operation": "download dbVar structural variant records",
        "group": "G4054",
    },
    "exome_validated_variant_calls_vcf_tsv_bed": {
        "root": "variation/get-rm/",
        "operation": "download exome-validated variant calls (VCF, TSV, BED)",
        "group": "G4055",
    },
    "molecular_qtl_association_dumps_eqtl_meqtl": {
        "root": "eqtl/",
        "operation": "download molecular QTL association dumps (eQTL, meQTL)",
        "group": "G4057",
    },
    "1000_genomes_variant_call_sets_vcf": {
        "root": "1000genomes/ftp/",
        "operation": "download 1000 Genomes variant call sets (VCF)",
        "group": "G4058",
    },
    "1000_genomes_structural_variant_calls": {
        "root": "1000genomes/ftp/",
        "operation": "download 1000 Genomes structural variant calls",
        "group": "G4059",
    },
    "1000_genomes_alignments_and_sequence_reads_bam_cram_fastq": {
        "root": "1000genomes/ftp/",
        "operation": "download 1000 Genomes alignments and sequence reads (BAM, CRAM, FASTQ)",
        "group": "G4060",
    },
    "giab_benchmark_variant_calls_vcf_bed": {
        "root": "giab/ftp/release/",
        "operation": "download GIAB benchmark variant calls (VCF, BED)",
        "group": "G4061",
    },
    "giab_sample_sequencing_reads_and_alignments_fastq_bam": {
        "root": "giab/ftp/data/",
        "operation": "download GIAB sample sequencing reads and alignments (FASTQ, BAM)",
        "group": "G4062",
    },
    "giab_rna_seq_data": {
        "root": "giab/ftp/data_RNAseq/",
        "operation": "download GIAB RNA-seq data",
        "group": "G4063",
    },
    "giab_somatic_benchmark_data": {
        "root": "giab/ftp/data_somatic/",
        "operation": "download GIAB somatic benchmark data",
        "group": "G4064",
    },
    "seqc_somatic_mutation_benchmark_data": {
        "root": "seqc/ftp/Somatic_Mutation_WG/",
        "operation": "download SEQC somatic mutation benchmark data",
        "group": "G4065",
    },
    "sample_genome_variant_analysis_sets": {
        "root": "SampleData/Genomes/Genomic_Analysis_Packages/",
        "operation": "download sample genome variant analysis sets",
        "group": "G4066",
    },
    "analyzed_rna_seq_coverage_tracks": {
        "root": "SampleData/Genomes/Analyzed_RNAseq/",
        "operation": "download analyzed RNA-seq coverage tracks",
        "group": "G4067",
    },
    "hapmap_genotype_data": {
        "root": "hapmap/",
        "operation": "download HapMap genotype data",
        "group": "G4068",
    },
    "hapmap_inferred_genotypes": {
        "root": "hapmap/inferred_genotypes/",
        "operation": "download HapMap inferred genotypes",
        "group": "G4069",
    },
    "hapmap_raw_array_genotype_data": {
        "root": "hapmap/raw_data/",
        "operation": "download HapMap raw array genotype data",
        "group": "G4070",
    },
    "hapmap_snp_genotyping_assay_designs": {
        "root": "hapmap/assays/",
        "operation": "download HapMap SNP genotyping assay designs",
        "group": "G4071",
    },
    "hapmap_allele_frequency_tables": {
        "root": "hapmap/frequencies/",
        "operation": "download HapMap allele frequency tables",
        "group": "G4072",
    },
    "hapmap_linkage_disequilibrium_tables": {
        "root": "hapmap/ld_data/",
        "operation": "download HapMap linkage disequilibrium tables",
        "group": "G4073",
    },
    "hapmap_phased_haplotypes": {
        "root": "hapmap/phasing/",
        "operation": "download HapMap phased haplotypes",
        "group": "G4074",
    },
    "hapmap_recombination_maps": {
        "root": "hapmap/recombination/",
        "operation": "download HapMap recombination maps",
        "group": "G4075",
    },
    "hapmap_copy_number_variant_data": {
        "root": "hapmap/cnv_data/",
        "operation": "download HapMap copy-number variant data",
        "group": "G4076",
    },
    "hapmap_mtdna_and_chry_haplogroup_assignments": {
        "root": "hapmap/mtDNA_and_chrY_haplogroups/",
        "operation": "download HapMap mtDNA and chrY haplogroup assignments",
        "group": "G4077",
    },
    "hapmap_allocated_snp_lists": {
        "root": "hapmap/allocated_snps/",
        "operation": "download HapMap allocated SNP lists",
        "group": "G4078",
    },
    "perlegen_resequencing_amplicon_sequences": {
        "root": "hapmap/perlegen_amplicons/",
        "operation": "download Perlegen resequencing amplicon sequences",
        "group": "G4079",
    },
    "j_watson_genome_snp_calls_and_blocks": {
        "root": "hapmap/jimwatsonsequence/",
        "operation": "download J. Watson genome SNP calls and blocks",
        "group": "G4080",
    },
    "hapmap_tag_snp_tsnp_data": {
        "root": "hapmap/tscsnp-data/",
        "operation": "download HapMap tag-SNP (tSNP) data",
        "group": "G4081",
    },
    "epigenome_track_files_chip_seq_dna_methylation": {
        "root": "epigenomics/",
        "operation": "download epigenome track files (ChIP-seq, DNA methylation)",
        "group": "G4082",
    },
    "genome_browser_bigwig_tracks": {
        "root": "bigwig/",
        "operation": "download genome browser bigWig tracks",
        "group": "G4083",
    },
    "sra_sequencing_run_data_reads": {
        "root": "sra/data/",
        "operation": "download SRA sequencing run data (reads)",
        "group": "G4084",
    },
    "pathogen_isolate_assemblies_and_results": {
        "root": "pathogen/Results/",
        "operation": "download pathogen isolate assemblies and results",
        "group": "G4085",
    },
    "pathogen_reference_genome_sets": {
        "root": "pathogen/Reference/",
        "operation": "download pathogen reference genome sets",
        "group": "G4086",
    },
    "pathogen_reference_phylogenetic_trees": {
        "root": "pathogen/Reference_Trees/",
        "operation": "download pathogen reference phylogenetic trees",
        "group": "G4087",
    },
    "amrfinderplus_reference_database": {
        "root": "pathogen/Antimicrobial_resistance/Data/",
        "operation": "download AMRFinderPlus reference database",
        "group": "G4088",
    },
    "geo_datasets_gds": {
        "root": "geo/datasets/",
        "operation": "download GEO DataSets (GDS)",
        "group": "G4090",
    },
    "geo_platform_annotation_gpl": {
        "root": "geo/platforms/",
        "operation": "download GEO platform annotation (GPL)",
        "group": "G4091",
    },
}


def _relative(value: str, *, directory: bool = False) -> str:
    if not isinstance(value, str) or (not value and not directory):
        raise ConfigurationError("A non-empty relative file path is required.")
    decoded = unquote(value)
    if (
        decoded.startswith("/")
        or any(p in {".", ".."} for p in decoded.split("/"))
        or any(c in decoded for c in ("\\", "\x00", "?", "#", ":"))
    ):
        raise ConfigurationError("Paths must stay within the selected collection root.")
    if not directory and decoded.endswith("/"):
        raise ConfigurationError("A file path, rather than a directory, is required.")
    return quote(decoded, safe="/._-")


def _root(collection: str) -> str:
    if collection not in NCBI_FTP_COLLECTIONS:
        raise ConfigurationError("Unknown collection; inspect NCBI_FTP_COLLECTIONS.")
    return "https://ftp.ncbi.nlm.nih.gov/" + NCBI_FTP_COLLECTIONS[collection]["root"]


class _Links(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag == "a":
            self.links.extend(value for key, value in attrs if key == "href" and value)


def list_ncbi_ftp_collection(
    collection: str,
    *,
    directory: str = "",
    timeout: float = 60,
    max_bytes: int = 2_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """List release directories or files of a curated NCBI collection, including GenBank
    catalogs.
    """
    root = _root(collection)
    relative = _relative(directory, directory=True).rstrip("/")
    url = root + (relative + "/" if relative else "")
    # Directory listings are HTML by design, unlike data-record API responses.
    from dnakit.search._http import request_text
    from dnakit.search.models import SearchConfig

    text = request_text(
        url, SearchConfig(timeout=timeout, max_response_bytes=max_bytes), provider="NCBI FTP"
    )
    parser = _Links()
    parser.feed(text)
    files = []
    for href in parser.links:
        full = urljoin(url, href)
        parts = urlsplit(full)
        if not full.startswith(root) or parts.query or parts.fragment or full == url:
            continue
        path = unquote(parts.path[len(urlsplit(root).path) :])
        _relative(path, directory=full.endswith("/"))
        files.append({"path": path, "url": full, "is_directory": full.endswith("/")})
    return {
        "status": "completed",
        "collection": collection,
        "files": files,
        "provenance": {"source": url, "retrieved_at": datetime.now(timezone.utc).isoformat()},
    }


def download_ncbi_ftp_collection(
    collection: str,
    paths: Sequence[str],
    output_dir: str,
    *,
    expected_md5: Mapping[str, str] | None = None,
    max_file_bytes: int = 20_000_000_000,
    max_total_bytes: int = 50_000_000_000,
    overwrite: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download selected dataset files atomically with progress, MD5/SHA256 and a manifest.

    ``collection`` parameterises release, file format and dataset selections from
    the same upstream FTP provider. ``paths`` are relative paths returned by
    ``list_ncbi_ftp_collection``. This function writes bytes, not just URLs.
    Optional MD5 expectations use the caller's unencoded relative paths as keys.
    """
    root = _root(collection)
    if (
        isinstance(paths, (str, bytes))
        or not isinstance(paths, Sequence)
        or not 1 <= len(paths) <= 100
    ):
        raise ConfigurationError("paths must contain 1-100 relative file paths.")
    hashes = dict(expected_md5 or {})
    if set(hashes) - set(paths):
        raise ConfigurationError("expected_md5 contains a path not requested.")
    resources = []
    for path in paths:
        encoded = _relative(path)
        filename = PurePosixPath(unquote(encoded)).name
        resources.append(
            RemoteFile(root + encoded, filename=filename, expected_md5=hashes.get(path))
        )
    if len({r.filename for r in resources}) != len(resources):
        raise ConfigurationError(
            "Requested paths have duplicate basenames; use separate output directories."
        )
    config = DownloadConfig(
        max_file_bytes=max_file_bytes,
        max_total_bytes=max_total_bytes,
        max_files=100,
        overwrite=overwrite,
        show_progress=verbose,
    )
    result = dataset(
        resources,
        output_dir,
        kind=collection,
        source="NCBI FTP",
        config=config,
        metadata={
            "collection_root": root,
            "operation": NCBI_FTP_COLLECTIONS[collection]["operation"],
            "retrieved_at": datetime.now(timezone.utc).isoformat(),
        },
    )
    return result.to_dict()


__all__ = ["download_ncbi_ftp_collection", "list_ncbi_ftp_collection"]
