"""Canonical tutorial locations for the migrated public APIs.

Source-provider groups describe provenance; routes describe the user's task.
Keep individual overrides explicit so a provider's unrelated operations do not
end up together merely because they share an implementation module.
"""

SECTIONS = {
    "discovery": (
        "22_database_query_discovery",
        "查找与调用数据库接口",
        "Discover and call database APIs",
    ),
    "sequences": (
        "22_database_query_sequences",
        "序列、组装与测序数据查询",
        "Sequence, assembly and sequencing queries",
    ),
    "phenotypes": (
        "22_database_query_phenotypes",
        "疾病、表型与本体查询",
        "Disease, phenotype and ontology queries",
    ),
    "regulation": (
        "22_database_query_regulation",
        "表达、调控与肿瘤数据查询",
        "Expression, regulation and cancer queries",
    ),
    "organisms": (
        "22_database_query_organisms",
        "模式生物、物种与微生物组查询",
        "Organism, taxonomy and microbiome queries",
    ),
    "pathways": (
        "22_database_query_pathways",
        "通路、互作、蛋白与结构查询",
        "Pathway, interaction, protein and structure queries",
    ),
    "targets": (
        "22_database_query_targets",
        "靶点、药物与关联证据查询",
        "Target, drug and association queries",
    ),
    "literature": (
        "22_database_query_literature",
        "文献、预印本与临床研究查询",
        "Literature, preprint and clinical study queries",
    ),
    "encode": (
        "22_database_query_encode",
        "ENCODE 实验与文件查询",
        "ENCODE experiment and file queries",
    ),
    "gene-annotation": (
        "26_gene_annotation",
        "基因标识、转录本与同源注释查询",
        "Gene, transcript and homology annotation queries",
    ),
    "variant-records": (
        "22_database_query_variant_frequency",
        "变异记录与群体频率查询",
        "Variant Records and Population Frequencies",
    ),
    "structural-variants": (
        "22_database_query_structural_variants",
        "结构变异查询",
        "Structural Variant Queries",
    ),
    "clinical-variants": (
        "22_database_query_clinical_variants",
        "ClinVar 临床意义查询",
        "ClinVar Clinical Significance",
    ),
    "genetic-evidence": (
        "22_database_query_genetic_evidence",
        "遗传病关联、剂量与可干预性",
        "Gene–Disease Evidence, Dosage and Actionability",  # noqa: RUF001
    ),
    "cancer-variants": (
        "22_database_query_cancer_variants",
        "肿瘤变异与治疗证据查询",
        "Cancer Variants and Treatment Evidence",
    ),
    "gwas": (
        "22_database_query_gwas",
        "GWAS 研究与关联查询",
        "GWAS Studies and Associations",
    ),
    "phewas": (
        "22_database_query_phewas",
        "队列表型关联与精细定位查询",
        "Cohort Phenotype Associations and Fine-mapping",
    ),
    "population-resources": (
        "22_database_query_population_resources",
        "人群、样本与连锁查询",
        "Populations, Samples and Linkage Disequilibrium",
    ),
    "variant-annotation": (
        "35_variant_annotation",
        "变异注释与标识转换",
        "Variant Annotation and Identifier Conversion",
    ),
    "variant-scores": (
        "36_variant_effect_scores",
        "变异效应评分与实验数据",
        "Variant Effect Scores and Experimental Data",
    ),
    "variant-splicing": (
        "23_prediction_variant_effects",
        "变异效应与蛋白质稳定性",
        "Variant Effects and Protein Stability",
    ),
    "selection": (
        "30_phylogeny",
        "系统发育与分子进化",
        "Phylogenetics and Molecular Evolution",
    ),
    "variant-statistics": (
        "33_variant_calling",
        "变异检测与统计",
        "Variant Calling and Statistics",
    ),
    "enrichment": (
        "37_functional_enrichment",
        "基因集与通路富集",
        "Gene Set and Pathway Enrichment",
    ),
    "model-tracks": (
        "23_prediction_regulatory_tracks",
        "表达、调控与接触图预测",
        "Expression, regulatory tracks and contact maps",
    ),
    "model-variants": (
        "23_prediction_variant_effects",
        "变异效应与蛋白质稳定性",
        "Variant effects and protein stability",
    ),
    "model-proteins": (
        "23_prediction_protein_features",
        "蛋白质特征与结构预测",
        "Protein features and structure",
    ),
    "downloads": (
        "15_download",
        "按编号与链接下载文件的参数",
        "Accession and URL download parameters",
    ),
    "encode-downloads": (
        "15_download",
        "ENCODE 文件与批量下载接口",
        "ENCODE file and batch download APIs",
    ),
    "encode-collection": (
        "10_datasets",
        "ENCODE 实验整理、引用与来源记录",
        "ENCODE collections, citations and provenance",
    ),
    "encode-credentials": ("14_engineering", "ENCODE 凭据管理", "ENCODE credential management"),
    "execution": (
        "14_engineering",
        "随包脚本、执行结果与报告接口",
        "Bundled scripts, execution results and report APIs",
    ),
    "blast": ("09_search", "远程蛋白序列 BLAST 检索", "Remote protein BLAST search"),
    "sequence-direction": (
        "04_ops_001_sequence_direction",
        "反向互补计算接口",
        "Reverse complement API",
    ),
    "translation": (
        "04_ops_002_transcription_translation",
        "多阅读框翻译接口",
        "Reading frame translation API",
    ),
    "sequence-descriptors": (
        "05_all_descriptors",
        "序列统计与组成计算接口",
        "Sequence statistics and composition APIs",
    ),
    "digest": ("04_ops_004_sequence_editing", "虚拟限制性酶切", "Virtual restriction digestion"),
    "cloning": (
        "04_ops_006_sequence_concatenation",
        "Golden Gate 与 Gibson 拼接设计",
        "Golden Gate and Gibson assembly design",
    ),
    "vcf": ("33_variant_calling", "VCF 规范化与统计接口", "VCF normalization and statistics APIs"),
    "inheritance": (
        "34_inheritance",
        "独立位点遗传模拟接口",
        "Independent-locus inheritance API",
    ),
}

FUNCTION_ROUTES = {}


def _assign(namespace: str, section: str, names: str) -> None:
    for name in names.split():
        qualified = f"dnakit.{namespace}.{name}"
        if qualified in FUNCTION_ROUTES:
            raise ValueError(f"Duplicate tutorial route: {qualified}")
        FUNCTION_ROUTES[qualified] = section


_assign("databases", "discovery", "call describe execute execute_tool list_operations search")
_assign("databases", "downloads", "download_file")
_assign("download", "downloads", "ncbi.sequence ena.sequence")
_assign(
    "databases",
    "gene-annotation",
    """
    ensembl_compara_get_gene_tree ensembl_compara_get_orthologues
    ensembl_compara_get_paralogues ensembl_get_homology ensembl_get_overlap_features
    ensembl_get_regulatory_features ensembl_lookup_gene hgnc_fetch_gene_by_symbol
    my_gene_batch_query my_gene_get_gene_annotation my_gene_query_genes
    ncbidatasets_get_gene ncbidatasets_get_gene_by_symbol ncbidatasets_get_orthologs
    ncbigene_get_summary ncbigene_search tark_get_mane_transcripts tark_get_transcript
    variant_validator_gene2transcripts
""",
)
_assign(
    "databases",
    "variant-records",
    """
    dbsnp_get_frequencies
    dbsnp_get_variant_by_rsid
    gnomad_get_gene_constraints
    gnomad_get_variant
    gnomad_get_variant_populations
    gnomad_search_variants
""",
)
_assign(
    "databases",
    "structural-variants",
    """
    ensembl_get_structural_variants
    ensembl_get_sv_detail
    gnomad_get_sv_by_gene
    gnomad_get_sv_by_region
    gnomad_get_sv_detail
""",
)
_assign(
    "databases",
    "clinical-variants",
    """
    clin_var_get_clinical_significance
    clin_var_get_variant_details
    clin_var_search_variants
""",
)
_assign(
    "databases",
    "genetic-evidence",
    """
    clin_gen_dosage_by_gene
    clin_gen_dosage_region_search
    clin_gen_get_actionability_adult
    clin_gen_get_actionability_pediatric
    clin_gen_get_dosage_sensitivity
    clin_gen_get_gene_validity
    clin_gen_get_variant_classifications
    clin_gen_search_actionability
    clin_gen_search_dosage_sensitivity
    clin_gen_search_gene_validity
""",
)
_assign(
    "databases",
    "cancer-variants",
    """
    civic_get_molecular_profile
    civic_get_variant
    civic_get_variants_by_gene
    civic_search_assertions
    civic_search_evidence_items
    civic_search_genes
    civic_search_therapies
    cosmic_get_mutations_by_gene
    cosmic_search_mutations
    onco_kb_annotate_variant
    progenetix_cnv_search
    progenetix_list_cohorts
    progenetix_list_filtering_terms
    progenetix_search_biosamples
""",
)
_assign(
    "databases",
    "gwas",
    """
    gwas_get_associations_for_snp
    gwas_get_associations_for_study
    gwas_get_associations_for_trait
    gwas_get_snp_by_id
    gwas_get_snps_for_gene
    gwas_get_studies_for_trait
    gwas_get_study_by_id
    gwas_get_variants_for_trait
    gwas_search_associations
    gwas_search_associations_by_gene
    gwas_search_snps
    gwas_search_studies
""",
)
_assign(
    "databases",
    "phewas",
    """
    bio_bank_japan_phewas_by_variant
    finn_gen_get_variant_finemapping
    genebass_gene_burden_phewas
    tpmi_phewas_by_variant
    ukbtopmed_phewas_by_variant
""",
)
_assign(
    "databases",
    "population-resources",
    """
    igsr_list_data_collections
    igsr_search_populations
    igsr_search_samples
    ldlink_get_proxies
""",
)
_assign(
    "databases",
    "variant-annotation",
    """
    ensembl_var_get_variant_consequences
    ensembl_vep_annotate_hgvs
    ensembl_vep_annotate_rsid
    ensembl_vep_region
    ensembl_vep_variant_recoder
    favor_annotate_variant
    my_variant_get_variant_annotation
    my_variant_query_variants
    variant_validator_validate_variant
""",
)
_assign(
    "databases",
    "variant-scores",
    """
    alpha_missense_get_protein_scores
    alpha_missense_get_residue_scores
    alpha_missense_get_variant_score
    cadd_get_position_scores
    cadd_get_range_scores
    cadd_get_scores
    cadd_get_variant_score
    eve_get_gene_info
    eve_get_variant_score
    mave_db_get_effect_matrix
    my_variant_get_pathogenicity_scores
""",
)
_assign(
    "databases",
    "variant-splicing",
    """
    splice_ai_get_max_delta
    splice_ai_predict_pangolin
    splice_ai_predict_splice
""",
)
_assign("databases", "variant-statistics", "coding_variant_fraction")
_assign(
    "databases",
    "enrichment",
    """
    enrichr_enrich enrichr_gene_enrichment_analysis panther_enrichment
    reactome_analysis_pathway_enrichment string_functional_enrichment
    string_get_functional_enrichment
""",
)
_assign("databases", "model-tracks", "alpha_genome_predict_interval")
_assign(
    "databases",
    "model-variants",
    """
    alpha_genome_score_variant dyna_mut2_predict_stability esm2_score_missense_variant
    esm_explain_variant_mechanism esm_score_variant_sae_batch esm_score_variant_sae_disruption
    evo2_score_variant
""",
)
_assign(
    "databases",
    "model-proteins",
    """
    esm_describe_sae_feature esm_get_region_sae_features esm_get_sae_features
    esm_score_sequence nvidia_nim_alphafold2 nvidia_nim_esmfold
""",
)
_assign("databases", "blast", "blast_protein_search")
_assign("databases", "sequence-direction", "sequence_reverse_complement")
_assign("databases", "translation", "dna_translate_reading_frames")
_assign(
    "databases",
    "sequence-descriptors",
    "sequence_count_residues sequence_gc_content sequence_stats",
)
_assign("databases", "digest", "dna_virtual_digest")
_assign("databases", "cloning", "dna_gibson_design dna_golden_gate_assemble dna_golden_gate_design")
_assign("databases", "selection", "sequence_dn_ds")
_assign("databases", "vcf", "vcf_count_variants vcf_normalize vcf_summary_stats")
_assign("workflows", "inheritance", "simulate_inheritance")
_assign("encode", "encode-downloads", "batch_download download_files")
_assign(
    "encode",
    "encode-collection",
    """
    compare_experiments export_data get_citations get_provenance get_references
    link_reference list_tracked log_derived_file summarize_collection track_experiment
""",
)
_assign("encode", "encode-credentials", "manage_credentials")


def section_for(record: dict) -> str:
    """Return one canonical tutorial section for a public function."""
    section = FUNCTION_ROUTES.get(record["name"], record["group"])
    if section not in SECTIONS:
        raise ValueError(f"Missing tutorial route: {record['name']}")
    return section


def page_for(record: dict) -> str:
    return SECTIONS[section_for(record)][0]
