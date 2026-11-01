# P42 公共数据库操作

通过固定的官方端点执行真实检索和下载。普通测试与下方回放示例不联网；未声称已经逐条完成线上数据验证。

下载继承字节限制、原子写入、校验和与进度条。返回 JSON 可序列化结果，保留服务响应、查询日期和分页信息；分页需按上游令牌显式继续。旧数据集文件名与版本需先通过目录列表确认。

## `download_ncbi_ftp_collection()`

- **作用:** Download selected dataset files atomically with progress, MD5/SHA256 and a manifest.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `download_ncbi_ftp_collection(collection: 'str', paths: 'Sequence[str]', output_dir: 'str', *, expected_md5: 'Mapping[str, str] | None' = None, max_file_bytes: 'int' = 20000000000, max_total_bytes: 'int' = 50000000000, overwrite: 'bool' = False, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `list_ncbi_ftp_collection()`

- **作用:** List release directories or files of a curated NCBI collection, including GenBank
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `list_ncbi_ftp_collection(collection: 'str', *, directory: 'str' = '', timeout: 'float' = 60, max_bytes: 'int' = 2000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## 数据集选择

| collection | 来源目录 |
|---|---|
| `cds_nucleotide_sequences_fasta` | `genomes/all/` |
| `rna_sequences_fasta` | `genomes/all/` |
| `repeatmasker_repeat_annotation_rm_out` | `genomes/all/` |
| `genome_feature_tables` | `genomes/all/` |
| `genome_feature_count_summaries` | `genomes/all/` |
| `assembly_structure_and_region_descriptors_agp` | `genomes/all/` |
| `mapview_genome_map_data_sequence_and_feature_objects` | `genomes/MapView/` |
| `influenza_nucleotide_sequences_and_annotation_tables` | `genomes/INFLUENZA/` |
| `rrna_targeted_locus_sequences_16s_18s_23s_28s_5s_its` | `refseq/TargetedLoci/` |
| `prokaryote_ani_and_checkm_genome_quality_reports` | `genomes/ASSEMBLY_REPORTS/` |
| `pairwise_genome_alignments_comparative_genome_viewer` | `comparative-genome-viewer/` |
| `genome_reference_consortium_assembly_sequences` | `pub/grc/` |
| `nist_pathogen_reference_genome_assemblies_fasta` | `nist-immsa/NISTPathogen/` |
| `nucleotide_blast_database_archives` | `blast/db/` |
| `windowmasker_genome_masking_files` | `blast/windowmasker_files/` |
| `genbank_nucleotide_release_flat_files` | `genbank/` |
| `genbank_est_gss_and_other_record_catalogs_and_gene_lists` | `genbank/catalog/` |
| `third_party_annotation_tpa_records_in_genbank_flat_file_format` | `tpa/release/gbff/` |
| `third_party_annotation_tpa_nucleotide_sequences_fasta` | `tpa/release/` |
| `refseq_wgs_master_contig_records` | `refseq/wgs/` |
| `refseq_functional_element_biological_region_sequences` | `refseq/FunctionalElements/` |
| `mane_select_transcript_set` | `refseq/MANE/MANE_human/` |
| `ncbi_gene_records_table_gene_info` | `gene/DATA/` |
| `gene_to_nucleotide_accession_mapping_gene2accession_gene2refseq` | `gene/DATA/` |
| `gene_to_ensembl_identifier_mapping_gene2ensembl` | `gene/DATA/` |
| `gene_to_go_annotation_table_gene2go` | `gene/DATA/` |
| `gene_orthology_table_gene_orthologs` | `gene/DATA/` |
| `gene_genomic_neighbourhood_table_gene_neighbors` | `gene/DATA/` |
| `gene_group_table_gene_group` | `gene/DATA/` |
| `discontinued_gene_identifier_history_gene_history` | `gene/DATA/` |
| `gene_summary_text_gene_summary` | `gene/DATA/` |
| `gene_refseq_uniprotkb_mapping_gene_refseq_uniprotkb_collab` | `gene/DATA/` |
| `ncbi_gene_expression_data_mammalia` | `gene/DATA/expression/` |
| `ccds_consensus_coding_sequence_records` | `pub/CCDS/` |
| `homologene_gene_homology_group_tables` | `pub/HomoloGene/` |
| `dbsnp_variant_calls_vcf` | `snp/` |
| `dbsnp_refsnp_records_json` | `snp/latest_release/JSON/` |
| `dbsnp_coordinate_bed_tracks` | `snp/organisms/human_9606/BED/` |
| `dbsnp_rs_flanking_sequences_fasta` | `snp/organisms/human_9606/rs_fasta/` |
| `dbsnp_genotype_tables` | `snp/organisms/human_9606/` |
| `dbsnp_population_allele_frequencies` | `snp/population_frequency/` |
| `dbsnp_submitted_variant_records` | `snp/submission/` |
| `dbvar_structural_variant_records` | `pub/dbVar/data/` |
| `exome_validated_variant_calls_vcf_tsv_bed` | `variation/get-rm/` |
| `molecular_qtl_association_dumps_eqtl_meqtl` | `eqtl/` |
| `1000_genomes_variant_call_sets_vcf` | `1000genomes/ftp/` |
| `1000_genomes_structural_variant_calls` | `1000genomes/ftp/` |
| `1000_genomes_alignments_and_sequence_reads_bam_cram_fastq` | `1000genomes/ftp/` |
| `giab_benchmark_variant_calls_vcf_bed` | `giab/ftp/release/` |
| `giab_sample_sequencing_reads_and_alignments_fastq_bam` | `giab/ftp/data/` |
| `giab_rna_seq_data` | `giab/ftp/data_RNAseq/` |
| `giab_somatic_benchmark_data` | `giab/ftp/data_somatic/` |
| `seqc_somatic_mutation_benchmark_data` | `seqc/ftp/Somatic_Mutation_WG/` |
| `sample_genome_variant_analysis_sets` | `SampleData/Genomes/Genomic_Analysis_Packages/` |
| `analyzed_rna_seq_coverage_tracks` | `SampleData/Genomes/Analyzed_RNAseq/` |
| `hapmap_genotype_data` | `hapmap/` |
| `hapmap_inferred_genotypes` | `hapmap/inferred_genotypes/` |
| `hapmap_raw_array_genotype_data` | `hapmap/raw_data/` |
| `hapmap_snp_genotyping_assay_designs` | `hapmap/assays/` |
| `hapmap_allele_frequency_tables` | `hapmap/frequencies/` |
| `hapmap_linkage_disequilibrium_tables` | `hapmap/ld_data/` |
| `hapmap_phased_haplotypes` | `hapmap/phasing/` |
| `hapmap_recombination_maps` | `hapmap/recombination/` |
| `hapmap_copy_number_variant_data` | `hapmap/cnv_data/` |
| `hapmap_mtdna_and_chry_haplogroup_assignments` | `hapmap/mtDNA_and_chrY_haplogroups/` |
| `hapmap_allocated_snp_lists` | `hapmap/allocated_snps/` |
| `perlegen_resequencing_amplicon_sequences` | `hapmap/perlegen_amplicons/` |
| `j_watson_genome_snp_calls_and_blocks` | `hapmap/jimwatsonsequence/` |
| `hapmap_tag_snp_tsnp_data` | `hapmap/tscsnp-data/` |
| `epigenome_track_files_chip_seq_dna_methylation` | `epigenomics/` |
| `genome_browser_bigwig_tracks` | `bigwig/` |
| `sra_sequencing_run_data_reads` | `sra/data/` |
| `pathogen_isolate_assemblies_and_results` | `pathogen/Results/` |
| `pathogen_reference_genome_sets` | `pathogen/Reference/` |
| `pathogen_reference_phylogenetic_trees` | `pathogen/Reference_Trees/` |
| `amrfinderplus_reference_database` | `pathogen/Antimicrobial_resistance/Data/` |
| `geo_datasets_gds` | `geo/datasets/` |
| `geo_platform_annotation_gpl` | `geo/platforms/` |
## 离线输入检查示例

**示例代码：**

```python
from dnakit.download.ncbi_ftp_p42_functions import NCBI_FTP_COLLECTIONS
print(len(NCBI_FTP_COLLECTIONS))
```

**示例输出：**

```text
78
```

## 上游文档

- [官方文档 1](https://ftp.ncbi.nlm.nih.gov/)
