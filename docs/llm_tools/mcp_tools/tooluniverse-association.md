# ToolUniverse: association (`tooluniverse`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **845** 个 DNA/基因相关工具（服务器共列出 2927 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

本页列出类别 `association` 的 138 个工具，其他类别：[`sequence`](tooluniverse-sequence.md) (44), [`genome_annotation`](tooluniverse-genome_annotation.md) (253), [`variant`](tooluniverse-variant.md) (155), [`expression`](tooluniverse-expression.md) (56), [`regulation_epigenome`](tooluniverse-regulation_epigenome.md) (54), [`sequencing_data`](tooluniverse-sequencing_data.md) (91), [`prediction_model`](tooluniverse-prediction_model.md) (17), [`molecular_tools`](tooluniverse-molecular_tools.md) (23), [`breeding_population`](tooluniverse-breeding_population.md) (14)。

| 项目 | 内容 |
| --- | --- |
| 源码 | [mims-harvard/ToolUniverse @ `9959ebe`](https://github.com/mims-harvard/ToolUniverse/tree/9959ebee7638773314362f1746980a15f5b17e5b) |
| 许可证 | Apache-2.0 |
| 启动命令 | `tooluniverse` |
| 传输方式 | stdio |
| 安装 | Install the pinned ToolUniverse distribution and required provider dependencies separately. Use the installed tooluniverse console script; unlike the README uvx --refresh example, this bridge does not refresh or install anything. |
| 副作用 | Starts the compact ToolUniverse MCP catalog. Selected tools may contact services, load models, execute programs or write files according to their provider definitions. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`Alliance_get_disease_genes`](#Alliance_get_disease_genes) | `association` | Lists genes linked to a Disease Ontology term across human and model organisms in the Alliance. |
| [`Alliance_get_gene_phenotypes`](#Alliance_get_gene_phenotypes) | `association` | Retrieves phenotype annotations for a gene from model organism databases via the Alliance. |
| [`BioBankJapan_phewas_by_variant`](#BioBankJapan_phewas_by_variant) | `association` | Returns phenotypes associated with a variant in BioBank Japan PheWeb, sorted by p-value. |
| [`cancer_biomarkers_disease_target_score`](#cancer_biomarkers_disease_target_score) | `association` | Gets Open Targets target-disease scores from the cancer biomarkers source for a disease EFO ID. |
| [`cancer_gene_census_disease_target_score`](#cancer_gene_census_disease_target_score) | `association` | Gets Open Targets target-disease scores from Cancer Gene Census for a disease EFO ID. |
| [`chembl_disease_target_score`](#chembl_disease_target_score) | `association` | Gets Open Targets target-disease scores from the ChEMBL source for a disease EFO ID. |
| [`ClinGen_get_actionability_adult`](#ClinGen_get_actionability_adult) | `association` | Returns ClinGen adult clinical actionability curations for genes. |
| [`ClinGen_get_actionability_pediatric`](#ClinGen_get_actionability_pediatric) | `association` | Returns ClinGen pediatric clinical actionability curations for genes. |
| [`ClinGen_get_gene_validity`](#ClinGen_get_gene_validity) | `association` | Returns ClinGen gene-disease validity classifications, optionally filtered by gene. |
| [`ClinGen_search_actionability`](#ClinGen_search_actionability) | `association` | Searches ClinGen actionability curations for a gene across adult and pediatric contexts. |
| [`ClinGen_search_gene_validity`](#ClinGen_search_gene_validity) | `association` | Searches ClinGen gene-disease validity curations by gene symbol. |
| [`Coloc_abf_test`](#Coloc_abf_test) | `association` | Runs coloc.abf Bayesian colocalization on two GWAS/QTL summary statistics over a shared SNP set, returning PP0-PP4. |
| [`CPIC_get_gene_drug_pairs`](#CPIC_get_gene_drug_pairs) | `association` | Gets all CPIC gene-drug pairs for a pharmacogene with CPIC and PharmGKB evidence levels. |
| [`CPIC_list_guidelines`](#CPIC_list_guidelines) | `association` | Lists CPIC pharmacogenomic guidelines, optionally filtered by gene or drug. |
| [`CPIC_search_gene_drug_pairs`](#CPIC_search_gene_drug_pairs) | `association` | Searches CPIC gene-drug pairs by gene symbol and/or CPIC evidence level. |
| [`disease_target_score`](#disease_target_score) | `association` | Gets Open Targets target-disease scores for a disease from a specified datasource. |
| [`DisGeNET_get_disease_genes`](#DisGeNET_get_disease_genes) | `association` | Gets DisGeNET genes associated with a disease, ranked by score. |
| [`DisGeNET_get_gda`](#DisGeNET_get_gda) | `association` | Gets DisGeNET gene-disease associations filtered by source and minimum score. |
| [`DisGeNET_get_vda`](#DisGeNET_get_vda) | `association` | Gets DisGeNET variant-disease associations by rsID or gene. |
| [`DisGeNET_search_disease`](#DisGeNET_search_disease) | `association` | Searches DisGeNET for genes associated with a disease given a UMLS CUI. |
| [`DisGeNET_search_gene`](#DisGeNET_search_gene) | `association` | Searches DisGeNET for diseases associated with a gene. |
| [`ensembl_get_variation_phenotypes`](#ensembl_get_variation_phenotypes) | `association` | Gets phenotype associations for a variant from Ensembl. |
| [`EnsemblLD_get_ld_pairwise`](#EnsemblLD_get_ld_pairwise) | `association` | Gets pairwise LD (r2, D') between two variants across 1000 Genomes populations. |
| [`EnsemblLD_get_ld_region`](#EnsemblLD_get_ld_region) | `association` | Gets all pairwise LD values among variants in a region for one population. |
| [`EnsemblLD_get_ld_variants`](#EnsemblLD_get_ld_variants) | `association` | Gets variants in LD with a variant in a 1000 Genomes population. |
| [`EnsemblPheno_get_by_gene`](#EnsemblPheno_get_by_gene) | `association` | Gets gene-phenotype associations from Ensembl (OMIM, ClinVar, GWAS catalog, Orphanet). |
| [`EnsemblPheno_get_by_region`](#EnsemblPheno_get_by_region) | `association` | Gets phenotype associations for features in a genomic region from Ensembl. |
| [`EnsemblPheno_get_by_term`](#EnsemblPheno_get_by_term) | `association` | Gets all variants and genes associated with a trait name or EFO/HP/MONDO accession. |
| [`EnsemblPheno_get_by_variant`](#EnsemblPheno_get_by_variant) | `association` | Gets phenotype associations for an rsID from GWAS catalog, ClinVar and others. |
| [`EpiGraphDB_get_gene_drug_associations`](#EpiGraphDB_get_gene_drug_associations) | `association` | Gets CPIC and PharmGKB pharmacogenomic drug associations for a gene via EpiGraphDB. |
| [`EpiGraphDB_get_genetic_correlations`](#EpiGraphDB_get_genetic_correlations) | `association` | Gets strong genetic correlations between a GWAS trait and other traits from OpenGWAS. |
| [`EpiGraphDB_get_mendelian_randomization`](#EpiGraphDB_get_mendelian_randomization) | `association` | Gets Mendelian randomization estimates between exposure and outcome traits from OpenGWAS. |
| [`EpiGraphDB_search_opengwas`](#EpiGraphDB_search_opengwas) | `association` | Searches IEU OpenGWAS for GWAS studies by trait name, returning study IDs for MR analyses. |
| [`eQTL_get_associations`](#eQTL_get_associations) | `association` | Gets eQTL Catalogue associations between variants and gene/exon/transcript expression for a dataset. |
| [`eQTL_list_datasets`](#eQTL_list_datasets) | `association` | Lists EBI eQTL Catalogue datasets by tissue and quantification method. |
| [`europepmc_disease_target_score`](#europepmc_disease_target_score) | `association` | Gets all Open Targets target-disease scores from Europe PMC text mining for a disease. |
| [`eva_disease_target_score`](#eva_disease_target_score) | `association` | Gets Open Targets target-disease scores from EVA (ClinVar) germline variants for a disease. |
| [`eva_somatic_disease_target_score`](#eva_somatic_disease_target_score) | `association` | Gets Open Targets target-disease scores from EVA somatic mutations for a disease. |
| [`expression_atlas_disease_target_score`](#expression_atlas_disease_target_score) | `association` | Gets Open Targets target-disease scores from Expression Atlas for a disease. |
| [`fda_pharmacogenomic_biomarkers`](#fda_pharmacogenomic_biomarkers) | `association` | Fetches the FDA Table of Pharmacogenomic Biomarkers in Drug Labeling, filterable by drug name or gene biomarker. |
| [`Finemap_credible_set`](#Finemap_credible_set) | `association` | Computes per-SNP posterior inclusion probabilities via Wakefield approximate Bayes factors and a credible set from effect sizes and SEs. |
| [`FinnGen_get_phenotype`](#FinnGen_get_phenotype) | `association` | Gets a FinnGen disease endpoint's case/control counts, genomic control lambda and number of genome-wide significant loci. |
| [`FinnGen_get_region_associations`](#FinnGen_get_region_associations) | `association` | Gets regional association and fine-mapping credible sets (lead variants, PIPs, MAF) for a FinnGen phenotype in a region. |
| [`FinnGen_get_variant_finemapping`](#FinnGen_get_variant_finemapping) | `association` | Returns FinnGen fine-mapping regions and phenotypes in which a GRCh38 chr:pos:ref:alt variant appears in a credible set. |
| [`FinnGen_list_phenotypes`](#FinnGen_list_phenotypes) | `association` | Searches/lists FinnGen disease endpoints by keyword or category with case/control counts and genome-wide significant loci. |
| [`FlyBase_get_gene_disease_models`](#FlyBase_get_gene_disease_models) | `association` | Returns genotype-based Drosophila disease models for a gene with associated human disease names and DOIDs. |
| [`G2P_download_panel`](#G2P_download_panel) | `association` | Downloads all curated gene-disease entries of a Gene2Phenotype panel with genotype, mechanism and evidence. |
| [`G2P_get_panel`](#G2P_get_panel) | `association` | Returns summary statistics (counts, description, update date) for a Gene2Phenotype clinical panel. |
| [`G2P_get_record`](#G2P_get_record) | `association` | Gets a G2P gene-disease record by stable ID with locus, inheritance, variant consequences and evidence. |
| [`G2P_search`](#G2P_search) | `association` | Searches Gene2Phenotype for curated gene-disease pairs with genotype and confidence by gene or disease. |
| [`gather_disease_profile`](#gather_disease_profile) | `association` | Aggregates a disease profile from Orphanet, OMIM, DisGeNET, Open Targets and OLS, including IDs, associated genes, inheritance and prevalence. |
| [`gather_gene_disease_associations`](#gather_gene_disease_associations) | `association` | Queries DisGeNET, OMIM, Open Targets, GenCC and ClinVar for gene-disease associations and reports cross-source concordance. |
| [`GenCC_get_classifications`](#GenCC_get_classifications) | `association` | Returns counts of GenCC gene-disease validity classifications, filterable by gene, disease or submitter. |
| [`GenCC_search_disease`](#GenCC_search_disease) | `association` | Finds genes with curated gene-disease validity classifications for a disease name or MONDO/OMIM ID. |
| [`GenCC_search_gene`](#GenCC_search_gene) | `association` | Gets curated gene-disease validity classifications (Definitive to Refuted) for a gene from GenCC submitters. |
| [`Genebass_gene_burden_phewas`](#Genebass_gene_burden_phewas) | `association` | Returns rare-variant gene burden associations across phenotypes from Genebass (UK Biobank exomes). |
| [`GeneNetwork_get_trait_info`](#GeneNetwork_get_trait_info) | `association` | Gets trait (probe/gene) details in a GeneNetwork dataset: gene symbol, position, mean expression, LRS, locus and p-value. |
| [`genomics_england_disease_target_score`](#genomics_england_disease_target_score) | `association` | Gets Open Targets target-disease scores from Genomics England PanelApp for a disease. |
| [`GTEx_calculate_eqtl`](#GTEx_calculate_eqtl) | `association` | Calculates eQTL statistics (p-value, NES, t-stat, genotypes) for any gene-variant pair in a GTEx tissue. |
| [`GTEx_get_eqtl_genes`](#GTEx_get_eqtl_genes) | `association` | Gets genes with significant cis-eQTLs (eGenes) per tissue with allelic fold change and q-values. |
| [`GTEx_get_finemapping_and_independent_eqtl`](#GTEx_get_finemapping_and_independent_eqtl) | `association` | Gets GTEx DAP-G fine-mapping credible sets or conditionally independent eQTLs for a gene and tissue. |
| [`GTEx_get_multi_tissue_eqtls`](#GTEx_get_multi_tissue_eqtls) | `association` | Gets Metasoft multi-tissue eQTL results (m-values, NES, p-values) for a gene-variant pair. |
| [`GTEx_get_single_tissue_eqtls`](#GTEx_get_single_tissue_eqtls) | `association` | Gets significant single-tissue gene-variant eQTL associations with p-values and NES. |
| [`GTEx_get_single_tissue_sqtls`](#GTEx_get_single_tissue_sqtls) | `association` | Gets GTEx single-tissue splicing QTLs or sGenes for a gene/variant and tissue. |
| [`GTEx_query_eqtl`](#GTEx_query_eqtl) | `association` | Queries GTEx single-tissue eQTL associations (variant, p-value, slope) for a gene. |
| [`GtoPdb_get_disease_associations`](#GtoPdb_get_disease_associations) | `association` | Gets GtoPdb curated disease-target and disease-ligand associations for a disease ID. |
| [`gwas_get_association_by_id`](#gwas_get_association_by_id) | `association` | Gets a specific GWAS Catalog association by its identifier. |
| [`gwas_get_associations_for_snp`](#gwas_get_associations_for_snp) | `association` | Gets all GWAS Catalog associations for an rsID with optional sorting. |
| [`gwas_get_associations_for_study`](#gwas_get_associations_for_study) | `association` | Gets all associations for a GWAS Catalog study accession sorted by p-value. |
| [`gwas_get_associations_for_trait`](#gwas_get_associations_for_trait) | `association` | Gets GWAS associations for a trait (resolved to EFO/MONDO) sorted by p-value. |
| [`gwas_get_snp_by_id`](#gwas_get_snp_by_id) | `association` | Gets a GWAS Catalog SNP record by rsID with location and mapped genes. |
| [`gwas_get_snps_for_gene`](#gwas_get_snps_for_gene) | `association` | Gets all GWAS Catalog SNPs mapped to a gene with rsIDs, positions and functional annotations. |
| [`gwas_get_studies_for_trait`](#gwas_get_studies_for_trait) | `association` | Gets GWAS studies for a trait with cohort, GxE and summary-statistics filters. |
| [`gwas_get_study_by_id`](#gwas_get_study_by_id) | `association` | Gets a GWAS Catalog study by its GCST accession. |
| [`gwas_get_variants_for_trait`](#gwas_get_variants_for_trait) | `association` | Searches the GWAS Catalog for variants associated with a disease/trait with rsIDs and p-values. |
| [`gwas_search_associations`](#gwas_search_associations) | `association` | Searches GWAS Catalog associations by trait, rsID, EFO ID, study accession or p-value threshold. |
| [`GWAS_search_associations_by_gene`](#GWAS_search_associations_by_gene) | `association` | Searches GWAS Catalog associations by gene name returning strongest risk allele and p-values. |
| [`gwas_search_snps`](#gwas_search_snps) | `association` | Searches GWAS Catalog SNPs by rsID or mapped gene returning location, genes and functional class. |
| [`gwas_search_studies`](#gwas_search_studies) | `association` | Searches GWAS studies by trait, cohort, GxE and summary-statistics availability. |
| [`GWASSumStats_get_region_associations`](#GWASSumStats_get_region_associations) | `association` | Queries full GWAS summary statistics for variants in a chromosomal region with p-values, betas and allele frequencies. |
| [`GWASSumStats_get_trait_studies`](#GWASSumStats_get_trait_studies) | `association` | Gets GWAS studies with deposited full summary statistics for an EFO trait ID. |
| [`GWASSumStats_list_studies`](#GWASSumStats_list_studies) | `association` | Lists GWAS Catalog studies with complete deposited variant-level summary statistics. |
| [`HPO_get_disease_annotations`](#HPO_get_disease_annotations) | `association` | Gets a disease's HPO annotation profile: MAxO medical actions, phenotypes by body system and causative genes. |
| [`HPO_get_genes_by_phenotype`](#HPO_get_genes_by_phenotype) | `association` | Gets genes in which variants cause an HPO phenotype term, with NCBI gene IDs and symbols. |
| [`IMPC_get_gene_phenotype_hits`](#IMPC_get_gene_phenotype_hits) | `association` | Gets IMPC phenotyping statistical results (p-values, effect sizes) for a mouse gene knockout. |
| [`IMPC_get_gene_summary`](#IMPC_get_gene_summary) | `association` | Gets IMPC mouse gene summary: phenotyping status, viability, orthologs and MP/HP annotations. |
| [`IMPC_get_phenotypes_by_gene`](#IMPC_get_phenotypes_by_gene) | `association` | Gets MP phenotype annotations from IMPC knockout screens for a gene, grouped by system. |
| [`IntOGen_get_drivers`](#IntOGen_get_drivers) | `association` | Returns cancer driver genes for a cancer type identified by IntOGen methods on somatic mutation data. |
| [`IntOGen_get_gene_info`](#IntOGen_get_gene_info) | `association` | Returns cancer types where a gene is an IntOGen cancer driver and the detection methods supporting it. |
| [`KEGG_get_disease`](#KEGG_get_disease) | `association` | Returns a KEGG disease entry with category, associated genes and mutations, pathways, drugs and cross-references. |
| [`KEGG_get_disease_genes`](#KEGG_get_disease_genes) | `association` | Returns human KEGG gene IDs linked to a KEGG disease through mutations, fusions or other mechanisms. |
| [`LDlink_get_proxies`](#LDlink_get_proxies) | `association` | Returns variants in linkage disequilibrium with a query SNP for a population via LDlink LDproxy. |
| [`MARRVEL_get_omim_phenotypes`](#MARRVEL_get_omim_phenotypes) | `association` | Returns OMIM phenotype associations for a human gene via MARRVEL with MIM numbers and inheritance. |
| [`MedGen_get_condition`](#MedGen_get_condition) | `association` | Returns a MedGen condition by UID or CUI with definition, associated genes, OMIM IDs, inheritance and HPO features. |
| [`MedGen_search_conditions`](#MedGen_search_conditions) | `association` | Searches MedGen genetic conditions by keyword and returns names, definitions, OMIM IDs and associated genes. |
| [`MedlinePlus_get_genetics_condition_by_name`](#MedlinePlus_get_genetics_condition_by_name) | `association` | Returns MedlinePlus Genetics information for a condition including related genes and inheritance pattern. |
| [`MGI_get_phenotypes`](#MGI_get_phenotypes) | `association` | Returns Mammalian Phenotype annotations for a mouse gene from mutant allele studies via Alliance/MGI. |
| [`Monarch_get_associations`](#Monarch_get_associations) | `association` | Queries Monarch gene-phenotype, disease-phenotype, gene-disease or gene-pathway associations by CURIE. |
| [`Monarch_get_gene_diseases`](#Monarch_get_gene_diseases) | `association` | Returns gene-disease associations with evidence for a gene CURIE from Monarch. |
| [`Monarch_get_gene_phenotypes`](#Monarch_get_gene_phenotypes) | `association` | Returns human and model organism phenotypes associated with a gene CURIE from Monarch. |
| [`MonarchV3_get_associations`](#MonarchV3_get_associations) | `association` | Queries cross-species Monarch associations by subject/object CURIE and Biolink category. |
| [`MR_estimate`](#MR_estimate) | `association` | Computes two-sample Mendelian randomization estimates (IVW, MR-Egger, weighted median) from SNP exposure/outcome betas. |
| [`MyDisease_get_disease`](#MyDisease_get_disease) | `association` | Returns MyDisease.info disease annotations from MONDO, DO, CTD, HPO and DisGeNET. |
| [`OMIM_get_entry`](#OMIM_get_entry) | `association` | Returns an OMIM entry by MIM number with description, clinical features, inheritance and molecular genetics. |
| [`OMIM_get_gene_map`](#OMIM_get_gene_map) | `association` | Returns OMIM gene map entries with chromosomal location, phenotypes and inheritance by MIM number or chromosome. |
| [`OMIM_search`](#OMIM_search) | `association` | Searches OMIM for Mendelian diseases and genes by keyword. |
| [`OncoKB_get_cancer_genes`](#OncoKB_get_cancer_genes) | `association` | Returns all OncoKB cancer genes classified as oncogenes or tumor suppressors. |
| [`OncoKB_get_gene_info`](#OncoKB_get_gene_info) | `association` | Returns OncoKB gene-level information on oncogene or tumor suppressor status with background. |
| [`OpenGWAS_get_mr_instruments`](#OpenGWAS_get_mr_instruments) | `association` | Assembles harmonized LD-clumped SNP instruments for an exposure/outcome GWAS pair from IEU OpenGWAS. |
| [`OpenTargets_get_associated_targets_by_disease_efoId`](#OpenTargets_get_associated_targets_by_disease_efoId) | `association` | Lists targets (genes) associated with a disease or phenotype by EFO ID, paginated with association scores. |
| [`OpenTargets_get_credible_set_colocalisation`](#OpenTargets_get_credible_set_colocalisation) | `association` | Retrieves colocalisation results (COLOC, eCAVIAR, H4) for a GWAS/QTL credible set by studyLocusId. |
| [`OpenTargets_get_credible_set_detail`](#OpenTargets_get_credible_set_detail) | `association` | Returns detail for a GWAS credible set: lead variant, study, fine-mapping method, significance and L2G predictions. |
| [`OpenTargets_get_diseases_phenotypes_by_target_ensembl`](#OpenTargets_get_diseases_phenotypes_by_target_ensembl) | `association` | Lists diseases and phenotypes associated with a target Ensembl gene ID, paginated with association scores. |
| [`OpenTargets_get_evidence_by_datasource`](#OpenTargets_get_evidence_by_datasource) | `association` | Retrieves Open Targets target-disease evidence filtered by chosen datasources such as GWAS credible sets, EVA or gene burden. |
| [`OpenTargets_get_gwas_study`](#OpenTargets_get_gwas_study) | `association` | Retrieves GWAS study metadata (traits, sample sizes, cohorts, summary-stat availability) by GCST accession. |
| [`OpenTargets_get_study_credible_sets`](#OpenTargets_get_study_credible_sets) | `association` | Lists credible sets for a GWAS study with lead variants, p-values, effect sizes and L2G scores. |
| [`OpenTargets_get_variant_credible_sets`](#OpenTargets_get_variant_credible_sets) | `association` | Lists GWAS credible sets in which a variant appears, with L2G predictions and effect sizes. |
| [`OpenTargets_search_gwas_studies_by_disease`](#OpenTargets_search_gwas_studies_by_disease) | `association` | Searches Open Targets GWAS studies for a disease or trait by EFO/MONDO IDs. |
| [`OpenTargets_target_disease_evidence`](#OpenTargets_target_disease_evidence) | `association` | Retrieves IntOGen somatic cancer driver evidence for a target-disease association. |
| [`Orphanet_get_gene_diseases`](#Orphanet_get_gene_diseases) | `association` | Lists Orphanet rare diseases associated with a gene symbol or gene name, with association details. |
| [`Orphanet_get_genes`](#Orphanet_get_genes) | `association` | Lists genes associated with an Orphanet rare disease with association types (causative, modifier, susceptibility). |
| [`PanelApp_get_panel`](#PanelApp_get_panel) | `association` | Returns all genes in a Genomics England PanelApp panel with confidence level, inheritance mode and phenotypes. |
| [`PanelApp_search_genes`](#PanelApp_search_genes) | `association` | Lists PanelApp panels containing a gene, with confidence level, inheritance and phenotypes. |
| [`PanelApp_search_panels`](#PanelApp_search_panels) | `association` | Searches Genomics England PanelApp gene panels by name or disease group. |
| [`PharmacoDB_get_biomarker_assoc`](#PharmacoDB_get_biomarker_assoc) | `association` | Returns gene-compound biomarker associations linking gene expression or mutation to drug sensitivity. |
| [`Pharos_get_disease_targets`](#Pharos_get_disease_targets) | `association` | Returns Pharos targets associated with a disease with Target Development Level filtering. |
| [`PomBase_get_gene_phenotypes`](#PomBase_get_gene_phenotypes) | `association` | Returns FYPO phenotype annotations and deletion viability for an S. pombe gene. |
| [`reactome_disease_target_score`](#reactome_disease_target_score) | `association` | Gets Open Targets target-disease scores from Reactome for a disease. |
| [`RGD_get_qtls_in_region`](#RGD_get_qtls_in_region) | `association` | Returns RGD QTLs overlapping a genomic region with LOD score and p-value. |
| [`SGD_get_disease_annotations`](#SGD_get_disease_annotations) | `association` | Gets curated Disease Ontology associations for a yeast gene from SGD, inferred through orthology to human disease genes. |
| [`SGD_get_phenotypes`](#SGD_get_phenotypes) | `association` | Gets curated phenotype annotations for a yeast gene from SGD, including mutant types, methods, strain backgrounds and conditions. |
| [`SPrediXcan_associate`](#SPrediXcan_associate) | `association` | Runs a summary-based S-PrediXcan TWAS test for a gene by combining eQTL prediction weights with GWAS per-SNP z-scores, returning z-score and p-value. |
| [`TIMER2_survival_association`](#TIMER2_survival_association) | `association` | Tests overall-survival association of high vs low expression of a gene in a TCGA cancer type with a log-rank test. |
| [`TPMI_phewas_by_variant`](#TPMI_phewas_by_variant) | `association` | Returns phenotypes associated with a variant in the Taiwan Precision Medicine Initiative PheWeb. |
| [`UKBTOPMed_phewas_by_variant`](#UKBTOPMed_phewas_by_variant) | `association` | Returns phenotypes associated with a variant in the UKB-TOPMed PheWeb. |
| [`WormBase_get_human_diseases`](#WormBase_get_human_diseases) | `association` | Gets human disease models (DOID) for a C. elegans gene inferred by orthology from WormBase. |
| [`WormBase_get_phenotypes`](#WormBase_get_phenotypes) | `association` | Gets observed and not-observed phenotype annotations for a C. elegans gene from RNAi, mutation and transgene experiments. |
| [`ZFIN_get_gene_phenotypes`](#ZFIN_get_gene_phenotypes) | `association` | Gets phenotype annotations observed when a zebrafish gene is mutated or knocked down. |

## `Alliance_get_disease_genes` {#Alliance_get_disease_genes}

Lists genes linked to a Disease Ontology term across human and model organisms in the Alliance.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`disease_id`
- **关键词：**disease genes, model organism, gene-disease association, doid
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `disease_id` | string | 是 | Disease Ontology (DO) ID. Examples: 'DOID:162' (cancer), 'DOID:10652' (Alzheimer's disease), 'DOID:9351' (diabetes mellitus), 'DOID:1612' (breast cancer). |
| `limit` | integer | 否 | Maximum number of gene associations to return (1-100). Default: 20. |
| `page` | integer | 否 | Page number for pagination. Default: 1. |

??? note "上游说明（英文原文）"

    Get genes associated with a disease across all model organisms from the Alliance of Genome Resources. Uses Disease Ontology (DO) IDs to find genes linked to a disease from human, mouse, rat, zebrafish, fly, worm, and other organisms. Includes both direct annotations and inferred associations from model organism phenotypes. Example: 'DOID:162' (cancer) returns TP53, BRCA1, and thousands more.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Alliance_get_disease_genes',
     'arguments': {'disease_id': '...'}},
)
```

## `Alliance_get_gene_phenotypes` {#Alliance_get_gene_phenotypes}

Retrieves phenotype annotations for a gene from model organism databases via the Alliance.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`gene_id`
- **关键词：**phenotype, knockout, mutant, gene-phenotype
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_id` | string | 是 | Alliance gene ID with source prefix. Examples: 'HGNC:6081' (human INS), 'MGI:98834' (mouse Trp53), 'FB:FBgn0003996' (fly w). |
| `limit` | integer | 否 | Maximum number of phenotypes to return (1-100). Default: 20. |
| `page` | integer | 否 | Page number for pagination. Default: 1. |

??? note "上游说明（英文原文）"

    Get phenotype annotations for a gene from the Alliance of Genome Resources. Returns phenotypes observed when the gene is mutated or knocked out, aggregated from model organism databases. Supports all model organisms (yeast, fly, worm, zebrafish, mouse, rat, frog). Example: 'HGNC:6081' (INS) returns phenotypes like 'Abnormal circulating insulin concentration'.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Alliance_get_gene_phenotypes',
     'arguments': {'gene_id': '...'}},
)
```

## `BioBankJapan_phewas_by_variant` {#BioBankJapan_phewas_by_variant}

Returns phenotypes associated with a variant in BioBank Japan PheWeb, sorted by p-value.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**phewas, biobank japan, rsid, p-value
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `rsid` | string | 否 | dbSNP rsID, e.g. 'rs7903146'. Auto-resolved to GRCh37 coordinates via Ensembl. Provide this OR 'variant'. |
| `variant` | string | 否 | Variant in chr:pos:ref:alt or chr-pos-ref-alt format, in GRCh37, e.g. '10:114758349:C:T'. Provide this OR 'rsid'. A bare rsID is also accepted here. |
| `limit` | integer | 否 | Maximum number of phenotype associations to return (sorted by ascending p-value). 默认 `20` |
| `max_pval` | number \| null | 否 | Optional: only return associations with p-value &lt;= this threshold (e.g. 5e-8 for genome-wide significance). |

??? note "上游说明（英文原文）"

    Phenome-wide association study (PheWAS) lookup for a single variant in BioBank Japan (BBJ), the largest Japanese-ancestry biobank. Returns every phenotype associated with the variant (the inverse of a GWAS), sorted by p-value. Accepts an rsID (auto-resolved to GRCh37 via Ensembl) or an explicit chr:pos:ref:alt variant in GRCh37. Use for cross-ancestry replication of a genetic association in an East Asian (Japanese) population.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'BioBankJapan_phewas_by_variant', 'arguments': {}},
)
```

## `cancer_biomarkers_disease_target_score` {#cancer_biomarkers_disease_target_score}

Gets Open Targets target-disease scores from the cancer biomarkers source for a disease EFO ID.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`efoId`
- **关键词：**open targets, target-disease, cancer biomarkers, efo
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `efoId` | string | 是 | The EFO (Experimental Factor Ontology) ID of the disease, e.g., 'MONDO_0011996' for chronic myelogenous leukemia |
| `pageSize` | integer | 否 | Number of results per page (default: 100, max: 100) 默认 `100` |

??? note "上游说明（英文原文）"

    Extract disease-target association scores from cancer biomarkers data. This includes known cancer biomarkers.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'cancer_biomarkers_disease_target_score',
     'arguments': {'efoId': '...'}},
)
```

## `cancer_gene_census_disease_target_score` {#cancer_gene_census_disease_target_score}

Gets Open Targets target-disease scores from Cancer Gene Census for a disease EFO ID.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`efoId`
- **关键词：**open targets, cancer gene census, target-disease, efo
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `efoId` | string | 是 | The EFO (Experimental Factor Ontology) ID of the disease, e.g., 'MONDO_0011996' for chronic myelogenous leukemia |
| `pageSize` | integer | 否 | Number of results per page (default: 100, max: 100) 默认 `100` |

??? note "上游说明（英文原文）"

    Extract disease-target association scores from Cancer Gene Census. This provides curated cancer gene data.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'cancer_gene_census_disease_target_score',
     'arguments': {'efoId': '...'}},
)
```

## `chembl_disease_target_score` {#chembl_disease_target_score}

Gets Open Targets target-disease scores from the ChEMBL source for a disease EFO ID.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`efoId`
- **关键词：**open targets, chembl, target-disease, efo
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `efoId` | string | 是 | The EFO (Experimental Factor Ontology) ID of the disease, e.g., 'MONDO_0011996' for chronic myelogenous leukemia |
| `pageSize` | integer | 否 | Number of results per page (default: 100, max: 100) 默认 `100` |

??? note "上游说明（英文原文）"

    Extract disease-target association scores specifically from ChEMBL database. ChEMBL provides bioactivity data for drug-target interactions.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'chembl_disease_target_score',
     'arguments': {'efoId': '...'}},
)
```

## `ClinGen_get_actionability_adult` {#ClinGen_get_actionability_adult}

Returns ClinGen adult clinical actionability curations for genes.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**actionability, incidental findings, brca1, clingen
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 否 | Optional: Filter by gene symbol |

??? note "上游说明（英文原文）"

    Get ClinGen clinical actionability curations for adult context. Returns genes with actionable findings, intervention scores, and outcome scores. Actionability scores help prioritize which incidental findings to return. Example: BRCA1/2 have high actionability for cancer screening.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ClinGen_get_actionability_adult', 'arguments': {}},
)
```

## `ClinGen_get_actionability_pediatric` {#ClinGen_get_actionability_pediatric}

Returns ClinGen pediatric clinical actionability curations for genes.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**actionability, pediatric, clingen, genes
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 否 | Optional: Filter by gene symbol |

??? note "上游说明（英文原文）"

    Get ClinGen clinical actionability curations for pediatric context. Returns genes with actionable findings in children, intervention recommendations, and evidence scores. Pediatric curations may differ from adult context due to age-specific interventions.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ClinGen_get_actionability_pediatric', 'arguments': {}},
)
```

## `ClinGen_get_gene_validity` {#ClinGen_get_gene_validity}

Returns ClinGen gene-disease validity classifications, optionally filtered by gene.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**gene-disease validity, definitive, clingen
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 否 | Optional: Filter by gene symbol |

??? note "上游说明（英文原文）"

    Get all ClinGen gene-disease validity curations. Returns comprehensive list of gene-disease relationships with classifications (Definitive, Strong, Moderate, Limited). Useful for batch analysis or exploring curated genes. Optional gene filter available.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ClinGen_get_gene_validity', 'arguments': {}},
)
```

## `ClinGen_search_actionability` {#ClinGen_search_actionability}

Searches ClinGen actionability curations for a gene across adult and pediatric contexts.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`gene`
- **关键词：**actionability, gene, clingen
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 是 | Gene symbol to search (e.g., BRCA1, MLH1, LDLR) |

??? note "上游说明（英文原文）"

    Search ClinGen clinical actionability across both adult and pediatric contexts by gene. Returns actionability data from both contexts for comparison. Useful for comprehensive actionability assessment of a gene across age groups.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ClinGen_search_actionability',
     'arguments': {'gene': '...'}},
)
```

## `ClinGen_search_gene_validity` {#ClinGen_search_gene_validity}

Searches ClinGen gene-disease validity curations by gene symbol.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`gene`
- **关键词：**gene-disease validity, inheritance, clingen
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 是 | Gene symbol to search (e.g., BRCA1, TP53, CFTR) |

??? note "上游说明（英文原文）"

    Search ClinGen gene-disease validity curations by gene symbol. Returns classification (Definitive, Strong, Moderate, Limited, Disputed, Refuted), disease associations, and inheritance patterns. Essential for ACMG variant classification. Example: BRCA1 returns validity scores for hereditary breast/ovarian cancer.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ClinGen_search_gene_validity',
     'arguments': {'gene': '...'}},
)
```

## `Coloc_abf_test` {#Coloc_abf_test}

Runs coloc.abf Bayesian colocalization on two GWAS/QTL summary statistics over a shared SNP set, returning PP0-PP4.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`beta1`, `se1`, `beta2`, `se2`
- **关键词：**coloc, colocalization, gwas, eqtl, causal variant, posterior probability
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `beta1` | array of number | 是 | Trait 1 per-SNP effect sizes over the region (same SNP order as trait 2). |
| `se1` | array of number | 是 | Standard errors of beta1 (positive). |
| `beta2` | array of number | 是 | Trait 2 per-SNP effect sizes (same SNP order/alleles). |
| `se2` | array of number | 是 | Standard errors of beta2 (positive). |
| `snp` | array \| null | 否 | Optional SNP identifiers (same order) to label the best shared variant. |
| `sd_prior` | number \| null | 否 | Prior SD of the effect size (default 0.15; 0.2 common for case/control log-OR). |
| `p1` | number \| null | 否 | Prior prob a SNP is causal for trait 1 (default 1e-4). |
| `p2` | number \| null | 否 | Prior prob a SNP is causal for trait 2 (default 1e-4). |
| `p12` | number \| null | 否 | Prior prob a SNP is causal for both traits (default 1e-5). |

??? note "上游说明（英文原文）"

    Bayesian colocalization (coloc.abf, Giambartolomei 2014) of two GWAS/QTL signals over the same SNP set in a region. Tests whether the two traits share a single causal variant. Returns the five posterior probabilities PP0-PP4; PP4 (&gt;0.8 = strong) is the probability of a shared causal variant - the standard way to link a GWAS locus to a gene via an eQTL. Pure-compute (no R, no LD matrix); pairs with eQTL_get_associations + GWAS summary tools.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Coloc_abf_test',
     'arguments': {'beta1': [], 'se1': [], 'beta2': [], 'se2': []}},
)
```

## `CPIC_get_gene_drug_pairs` {#CPIC_get_gene_drug_pairs}

Gets all CPIC gene-drug pairs for a pharmacogene with CPIC and PharmGKB evidence levels.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**pharmacogenomics, gene-drug pair, cpic level, pharmgkb
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `genesymbol` | string \| null | 否 | Gene symbol (e.g., 'CYP2D6', 'CYP2C19', 'SLCO1B1', 'TPMT', 'DPYD', 'VKORC1') |
| `gene_symbol` | string \| null | 否 | Alias for genesymbol (e.g., 'CYP2D6') |
| `gene` | string \| null | 否 | Alias for genesymbol (e.g., 'CYP2D6') |

??? note "上游说明（英文原文）"

    Get all gene-drug pairs for a pharmacogene from CPIC. Shows which drugs are affected by variants in a specific gene, with CPIC evidence level and PharmGKB clinical annotation level. Evidence levels: A (strongest) to D (weakest). Useful for understanding the clinical significance of a pharmacogene. Provide gene as 'genesymbol', 'gene_symbol', or 'gene'.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CPIC_get_gene_drug_pairs', 'arguments': {}},
)
```

## `CPIC_list_guidelines` {#CPIC_list_guidelines}

Lists CPIC pharmacogenomic guidelines, optionally filtered by gene or drug.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**cpic, pharmacogenomics, guidelines, gene-drug
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 否 | Filter by gene symbol (e.g., CYP2D6, TPMT). Returns only guidelines involving this gene. |
| `gene_symbol` | string | 否 | Alias for gene. Filter by gene symbol (e.g., CYP2D6, TPMT). |
| `drug` | string | 否 | Filter by drug name (e.g., 'codeine', 'warfarin', 'clopidogrel'). Case-insensitive substring match against drug names in the guideline. |
| `drug_name` | string | 否 | Alias for drug. |

??? note "上游说明（英文原文）"

    List all CPIC pharmacogenomic guidelines. Returns ~29 evidence-based guidelines (count may vary as guidelines are added/updated) for using pharmacogenomic test results to optimize drug therapy. Each guideline covers one or more gene-drug pairs with specific dosing recommendations based on genotype/phenotype. Includes guideline ID, name, URL, associated genes, and PharmGKB IDs. Use the returned guideline ID with CPIC_get_recommendations to retrieve specific dosing recommendations.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CPIC_list_guidelines', 'arguments': {}},
)
```

## `CPIC_search_gene_drug_pairs` {#CPIC_search_gene_drug_pairs}

Searches CPIC gene-drug pairs by gene symbol and/or CPIC evidence level.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**pharmacogenomics, gene-drug pair, cpic level, hgnc symbol
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `genesymbol` | string \| null | 否 | Gene symbol to filter by (e.g., 'CYP2D6', 'DPYD', 'TPMT'). Omit to search all genes. |
| `cpiclevel` | string \| null | 否 | CPIC evidence level to filter by (e.g., 'A', 'B', 'B/C', 'C', 'D'). Omit to include all levels. |
| `limit` | integer \| null | 否 | Maximum number of pairs to return. Omit to get every matching pair (CPIC publishes 635 in total, 96 at level A); there is no default cap, so 'count' equals the full match count unless you set this. |
| `gene_symbol` | string \| null | 否 | Gene symbol alias (e.g., 'CYP2D6', 'VKORC1') — alternative to genesymbol |
| `gene` | string \| null | 否 | Gene symbol alias (e.g., 'CYP2D6') — alternative to genesymbol |

??? note "上游说明（英文原文）"

    Search CPIC gene-drug pairs by gene symbol and/or CPIC evidence level. Gene symbol should be the standard HGNC symbol (e.g., "CYP2D6", "DPYD", "TPMT"). CPIC level indicates evidence strength: A=highest evidence, B=moderate, C=limited, D=for information only. Returns gene-drug pairs with guideline IDs for further lookup.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CPIC_search_gene_drug_pairs', 'arguments': {}},
)
```

## `disease_target_score` {#disease_target_score}

Gets Open Targets target-disease scores for a disease from a specified datasource.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`efoId`, `datasourceId`
- **关键词：**open targets, target-disease, datasource, efo
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `efoId` | string | 是 | The EFO (Experimental Factor Ontology) ID of the disease, e.g., 'MONDO_0011996' for chronic myelogenous leukemia |
| `datasourceId` | string | 是 | The datasource ID to extract scores from. Available options: 'gwas_credible_sets' (genetic/GWAS association), 'gene_burden', 'clinical_precedence', 'eva', 'eva_somatic', 'cancer_gene_census', 'cancer_biomarkers', 'europepmc', 'expression_atlas', 'genomics_england', 'impc', 'reactome', 'uniprot_literature', 'uniprot_variants'. Two IDs were retired upstream and now match nothing rather than errorin… |
| `pageSize` | integer | 否 | Number of results per page (default: 100, max: 100) 默认 `100` |

??? note "上游说明（英文原文）"

    Extract disease-target association scores from a specific data source using GraphQL API. This tool retrieves all targets associated with a disease and their scores from a specified datasource (e.g., clinical_precedence, eva, cancer_gene_census, etc.).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'disease_target_score',
     'arguments': {'efoId': '...', 'datasourceId': '...'}},
)
```

## `DisGeNET_get_disease_genes` {#DisGeNET_get_disease_genes}

Gets DisGeNET genes associated with a disease, ranked by score.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`disease`
- **关键词：**disgenet, disease genes, gene prioritization, gda score
- **需要 API 密钥：**`DISGENET_API_KEY`
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type (fixed: get_disease_genes) 默认 `get_disease_genes` |
| `disease` | string | 是 | Disease ID (UMLS CUI) or disease name |
| `min_score` | number | 否 | Minimum association score (0-1). Default: no filter |
| `limit` | integer | 否 | Maximum genes to return (default: 50) 默认 `50` |

??? note "上游说明（英文原文）"

    Get all genes associated with a disease from DisGeNET. Returns ranked gene list with association scores and evidence counts. Filter by minimum score to get high-confidence associations. Useful for disease gene prioritization.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DisGeNET_get_disease_genes',
     'arguments': {'disease': '...'}},
)
```

## `DisGeNET_get_gda` {#DisGeNET_get_gda}

Gets DisGeNET gene-disease associations filtered by source and minimum score.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**disgenet, gda, curated, evidence
- **需要 API 密钥：**`DISGENET_API_KEY`
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type (fixed: get_gda) 默认 `get_gda` |
| `gene` | string | 否 | Gene symbol (optional if disease provided) |
| `disease` | string | 否 | Disease ID/UMLS CUI (optional if gene provided) |
| `source` | string | 否 | Data source filter: CURATED, ANIMAL_MODELS, LITERATURE, etc. 可选值：`CURATED`, `ANIMAL_MODELS`, `LITERATURE`, `ALL` |
| `min_score` | number | 否 | Minimum GDA score (0-1). Higher = stronger evidence. |
| `limit` | integer | 否 | Maximum results (default: 25) 默认 `25` |

??? note "上游说明（英文原文）"

    Get gene-disease associations with filtering. Filter by data source (CURATED, LITERATURE, ANIMAL_MODELS) and minimum score. Returns comprehensive GDA data with evidence levels. Requires DISGENET_API_KEY.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DisGeNET_get_gda', 'arguments': {}},
)
```

## `DisGeNET_get_vda` {#DisGeNET_get_vda}

Gets DisGeNET variant-disease associations by rsID or gene.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**disgenet, vda, rsid, variant-disease
- **需要 API 密钥：**`DISGENET_API_KEY`
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type (fixed: get_vda) 默认 `get_vda` |
| `variant` | string | 否 | Variant rsID (e.g., rs1234) - optional if gene provided |
| `gene` | string | 否 | Gene symbol to get all variants - optional if variant provided |
| `limit` | integer | 否 | Maximum results (default: 25) 默认 `25` |

??? note "上游说明（英文原文）"

    Get variant-disease associations from DisGeNET. Query by rsID or gene symbol to find disease-associated variants. Returns variant IDs, diseases, and evidence scores. Useful for pharmacogenomics and variant interpretation.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DisGeNET_get_vda', 'arguments': {}},
)
```

## `DisGeNET_search_disease` {#DisGeNET_search_disease}

Searches DisGeNET for genes associated with a disease given a UMLS CUI.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`disease`
- **关键词：**disgenet, umls cui, disease genes
- **需要 API 密钥：**`DISGENET_API_KEY`
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type (fixed: search_disease) 默认 `search_disease` |
| `disease` | string | 是 | UMLS CUI for the disease, e.g. 'C0006142' (breast cancer). Free-text disease names are NOT accepted -- resolve one via umls_search_concepts first. |
| `limit` | integer | 否 | Maximum associations to return (default: 10) 默认 `10` |

??? note "上游说明（英文原文）"

    Search DisGeNET for genes associated with a disease. Requires a UMLS CUI, NOT a free-text disease name (e.g. 'C0006142' for breast cancer, not 'breast cancer') -- resolve a disease name to its CUI first via umls_search_concepts. Returns associated genes with scores and evidence. Requires DISGENET_API_KEY.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DisGeNET_search_disease',
     'arguments': {'disease': '...'}},
)
```

## `DisGeNET_search_gene` {#DisGeNET_search_gene}

Searches DisGeNET for diseases associated with a gene.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`gene`
- **关键词：**disgenet, gene-disease, association score
- **需要 API 密钥：**`DISGENET_API_KEY`
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type (fixed: search_gene) 默认 `search_gene` |
| `gene` | string | 是 | Gene symbol (e.g., BRCA1, TP53, EGFR) |
| `limit` | integer | 否 | Maximum associations to return (default: 10) 默认 `10` |

??? note "上游说明（英文原文）"

    Search DisGeNET for diseases associated with a gene. DisGeNET is the largest public database of gene-disease associations. Returns diseases, association scores, and evidence counts. Use the 'gene' parameter (e.g., 'BRCA1', 'TP53') — note: parameter name is 'gene', not 'gene_symbol'. Requires DISGENET_API_KEY env var (register free at disgenet.org).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DisGeNET_search_gene', 'arguments': {'gene': '...'}},
)
```

## `ensembl_get_variation_phenotypes` {#ensembl_get_variation_phenotypes}

Gets phenotype associations for a variant from Ensembl.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`id`
- **关键词：**variant phenotype, trait, gwas, effect size
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `id` | string | 是 | Variation ID (e.g., 'rs699'). Use ensembl_get_variants to find variant IDs in a genomic region. |
| `species` | string | 否 | Species name (default 'human') 默认 `human` |
| `phenotypes` | integer | 否 | Include phenotype associations (must be 1 for this tool) 默认 `1` |

??? note "上游说明（英文原文）"

    Get phenotype associations for a specific genetic variation. Returns detailed phenotype information including traits, studies, and effect sizes. Use this to understand the phenotypic consequences of a variant.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ensembl_get_variation_phenotypes',
     'arguments': {'id': '...'}},
)
```

## `EnsemblLD_get_ld_pairwise` {#EnsemblLD_get_ld_pairwise}

Gets pairwise LD (r2, D') between two variants across 1000 Genomes populations.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`variant1`, `variant2`
- **关键词：**linkage disequilibrium, r2, d prime, 1000 genomes
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant1` | string | 是 | rs ID of the first variant. Examples: 'rs6792369', 'rs429358'. |
| `variant2` | string | 是 | rs ID of the second variant. Examples: 'rs1042779', 'rs7412'. |

??? note "上游说明（英文原文）"

    Get pairwise linkage disequilibrium (LD) statistics between two specific variants from the Ensembl REST API. Returns r-squared and D' across ALL 1000 Genomes populations, allowing comparison of LD patterns across ancestries. Useful for assessing whether two GWAS signals are correlated or independent. Example: rs6792369 vs rs1042779 shows r2=1.0 in East Asian populations but r2=0.84 in African populations.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblLD_get_ld_pairwise',
     'arguments': {'variant1': '...', 'variant2': '...'}},
)
```

## `EnsemblLD_get_ld_region` {#EnsemblLD_get_ld_region}

Gets all pairwise LD values among variants in a region for one population.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`region`, `population`
- **关键词：**ld matrix, fine-mapping, r2, population
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `region` | string | 是 | Chromosomal window as 'chr:start..end' (1-based, GRCh38), max 1 Mb. Example: '6:25837556..25840000'. |
| `population` | string | 是 | 1000 Genomes population. Format: '1000GENOMES:phase_3:&lt;POP&gt;'. Common populations: 'CEU' (CEPH Utah), 'YRI' (Yoruba), 'CHB' (Han Chinese Beijing), 'GBR' (British), 'JPT' (Japanese Tokyo), 'KHV' (Kinh Vietnamese). Example: '1000GENOMES:phase_3:CEU'. |
| `r2_threshold` | number \| null | 否 | Minimum r-squared threshold to report. Set higher (e.g., 0.8) for strong LD only. |
| `d_prime_threshold` | number \| null | 否 | Minimum D' threshold to report. Set e.g., 0.8 for high D' only. |

??? note "上游说明（英文原文）"

    Get region-wide pairwise linkage disequilibrium (LD) from the Ensembl REST API: every r-squared and D' pair among ALL variants inside a chromosomal window for one 1000 Genomes Phase 3 population. This is the LD matrix used as input to statistical fine-mapping (e.g. SuSiE, FINEMAP) and LD-aware clumping at a GWAS locus. The window must be &lt;= 1 Mb (Ensembl limit). Region format: 'chr:start..end'. Example: region '6:25837556..25840000' in CEU returns 4 LD pairs, e.g. rs62394274/rs62394275 with r2=1.0 and d_prime=1.0.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblLD_get_ld_region',
     'arguments': {'region': '...', 'population': '...'}},
)
```

## `EnsemblLD_get_ld_variants` {#EnsemblLD_get_ld_variants}

Gets variants in LD with a variant in a 1000 Genomes population.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`variant_id`, `population`
- **关键词：**linkage disequilibrium, proxy snps, r2
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_id` | string | 是 | rs ID of the variant. Examples: 'rs1042779', 'rs429358', 'rs7903146', 'rs6792369'. |
| `population` | string | 是 | 1000 Genomes population. Format: '1000GENOMES:phase_3:&lt;POP&gt;'. Common populations: 'KHV' (Kinh Vietnamese), 'CEU' (CEPH Utah), 'YRI' (Yoruba), 'CHB' (Han Chinese Beijing), 'GBR' (British), 'JPT' (Japanese Tokyo), 'ASW' (African SW US). Example: '1000GENOMES:phase_3:CEU'. |
| `r2_threshold` | number \| null | 否 | Minimum r-squared threshold to report. Default: 0.05. Set higher (e.g., 0.8) for strong LD only. |
| `d_prime_threshold` | number \| null | 否 | Minimum D' threshold. Default: none. Set e.g., 0.8 for high D' only. |
| `limit` | integer \| null | 否 | Maximum number of LD partners to return, strongest r2 first. Default: 200. Increase to retrieve every partner when `truncated` is true. |

??? note "上游说明（英文原文）"

    Get linkage disequilibrium (LD) data for a variant from the Ensembl REST API. Returns variants in LD with the query variant in a specified 1000 Genomes population, with r-squared and D' statistics. LD measures the non-random association of alleles at different loci, critical for GWAS interpretation and fine-mapping. Results are sorted by r2 (strongest first) and capped at `limit` (default 200); `total_ld_count` always reports how many partners passed the thresholds and `truncated` flags whether any were withheld -- raise `limit` to retrieve them all. Example: rs1042779 in KHV population returns 708 LD partners with r2 and D' values.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblLD_get_ld_variants',
     'arguments': {'variant_id': '...', 'population': '...'}},
)
```

## `EnsemblPheno_get_by_gene` {#EnsemblPheno_get_by_gene}

Gets gene-phenotype associations from Ensembl (OMIM, ClinVar, GWAS catalog, Orphanet).

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`gene`
- **关键词：**gene phenotype, omim, gwas catalog, orphanet
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `species` | string | 否 | Species name. Use 'homo_sapiens' for human. Default: 'homo_sapiens'. 默认 `homo_sapiens` |
| `gene` | string | 是 | Gene symbol. Examples: 'BRCA1', 'TP53', 'EGFR', 'KRAS', 'BRAF'. |

??? note "上游说明（英文原文）"

    Get phenotype associations for a gene from the Ensembl REST API. Returns diseases, traits, and clinical conditions linked to the gene, aggregated from multiple sources including Cancer Gene Census, OMIM, ClinVar, NHGRI-EBI GWAS catalog, and Orphanet. Each entry includes the phenotype description, data source, genomic location, and ontology accessions (EFO, HP, Orphanet). Example: BRCA1 returns 162 phenotype associations across multiple cancer types and other conditions.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblPheno_get_by_gene', 'arguments': {'gene': '...'}},
)
```

## `EnsemblPheno_get_by_region` {#EnsemblPheno_get_by_region}

Gets phenotype associations for features in a genomic region from Ensembl.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`region`
- **关键词：**region phenotype, gwas, clinvar
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `species` | string | 否 | Species name. Use 'homo_sapiens' for human. Default: 'homo_sapiens'. 默认 `homo_sapiens` |
| `region` | string | 是 | Genomic region in format 'chromosome:start-end'. Example: '17:7661779-7687538' for TP53, '17:43044295-43125482' for BRCA1. Keep regions under 5Mb for reasonable response times. |
| `feature_type` | string \| null | 否 | Filter by feature type: 'Variation' (variants only), 'Gene' (genes only), 'QTL' (quantitative trait loci). Default: all types. |

??? note "上游说明（英文原文）"

    Get phenotype associations for all variants and genes in a genomic region from the Ensembl REST API. Returns disease/trait associations from GWAS catalog, ClinVar, Cancer Gene Census, etc. for features overlapping the specified region. Useful for understanding the phenotypic landscape of a genomic interval, especially for interpreting GWAS loci or candidate regions. Example: TP53 region (chr17:7661779-7687538) returns phenotype associations from the NHGRI-EBI GWAS catalog.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblPheno_get_by_region',
     'arguments': {'region': '...'}},
)
```

## `EnsemblPheno_get_by_term` {#EnsemblPheno_get_by_term}

Gets all variants and genes associated with a trait name or EFO/HP/MONDO accession.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**trait, gwas, risk allele, efo
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `species` | string | 否 | Species name. Use 'homo_sapiens' for human. Default: 'homo_sapiens'. 默认 `homo_sapiens` |
| `term` | string \| null | 否 | Trait or disease name. Examples: 'Alzheimer disease', 'Coronary artery disease', 'Type 2 diabetes mellitus'. If you pass an ontology accession (e.g. 'EFO:0000249') here it is auto-detected and routed correctly. |
| `accession` | string \| null | 否 | Ontology accession (EFO, HP, MONDO, Orphanet). Examples: 'EFO:0000249' (Alzheimer disease), 'HP:0002511' (Alzheimer disease phenotype), 'MONDO:0004975'. Mutually exclusive with 'term'. |

??? note "上游说明（英文原文）"

    Reverse phenotype lookup from the Ensembl REST API: given a trait/disease NAME (e.g. 'Alzheimer disease') OR an ontology accession (EFO/HP/MONDO, e.g. 'EFO:0000249', 'HP:0002511'), return ALL associated variants and genes with their risk allele, p-value, beta, clinical significance, and source (NHGRI-EBI GWAS catalog, ClinVar, dbGaP, etc.). This is the inverse of the gene/region/variant -&gt; phenotype tools: it starts from a phenotype and finds the underlying genetics. Provide the query under 'term' (name) or 'accession' (ontology ID); an ontology-looking value passed as 'term' is auto-routed to the accession endpoint. Example: 'Alzheimer disease' returns 4450 records (e.g. rs201377759 in APP); 'EFO:0000249' returns 3111 records (e.g. rs7810606 in EPHA1, GWAS p=4e-11). Up to 500 associations are returned.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblPheno_get_by_term', 'arguments': {}},
)
```

## `EnsemblPheno_get_by_variant` {#EnsemblPheno_get_by_variant}

Gets phenotype associations for an rsID from GWAS catalog, ClinVar and others.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`variant_id`
- **关键词：**rsid, phenotype, odds ratio, gwas
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `species` | string | 否 | Species name. Use 'homo_sapiens' for human. Default: 'homo_sapiens'. 默认 `homo_sapiens` |
| `variant_id` | string | 是 | Variant rsID. Examples: 'rs429358' (APOE e4), 'rs7903146' (TCF7L2 diabetes), 'rs1801133' (MTHFR). |

??? note "上游说明（英文原文）"

    Get phenotype associations for a specific variant (by rsID) from the Ensembl REST API. Returns traits, diseases, and phenotypes associated with the variant from NHGRI-EBI GWAS catalog, ClinVar, and other sources. Includes risk alleles, beta coefficients/odds ratios, p-values, and study references. Example: rs429358 (APOE epsilon 4 variant) returns 1004 phenotype associations including Alzheimer disease, LDL cholesterol levels, and metabolite measurements.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblPheno_get_by_variant',
     'arguments': {'variant_id': '...'}},
)
```

## `EpiGraphDB_get_gene_drug_associations` {#EpiGraphDB_get_gene_drug_associations}

Gets CPIC and PharmGKB pharmacogenomic drug associations for a gene via EpiGraphDB.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**pharmacogenomics, cpic, pharmgkb, gene-drug
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_name` | string | 否 | Gene symbol to find drug associations for (e.g., 'TP53', 'EGFR', 'DPYD', 'CYP2D6'). |
| `gene_symbol` | string | 否 | Alias for gene_name. Gene symbol (e.g., 'CYP2D6', 'BRCA1'). |
| `pval_threshold` | number | 否 | P-value threshold (default 0.05). 默认 `0.05` |

??? note "上游说明（英文原文）"

    Get pharmacogenomics drug-gene associations from CPIC and PharmGKB databases via EpiGraphDB. Returns drugs associated with a gene with clinical evidence levels. Example: 'DPYD' returns fluorouracil with CPIC Level A evidence.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EpiGraphDB_get_gene_drug_associations', 'arguments': {}},
)
```

## `EpiGraphDB_get_genetic_correlations` {#EpiGraphDB_get_genetic_correlations}

Gets strong genetic correlations between a GWAS trait and other traits from OpenGWAS.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`trait`
- **关键词：**genetic correlation, rg, gwas, opengwas
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `trait` | string | 是 | Exact, case-sensitive GWAS trait label (e.g., 'Waist circumference', 'Hip circumference'). Many labels have no node in the sparse genetic-correlation graph and return empty. |
| `pval_threshold` | number | 否 | Accepted for compatibility but NOT applied — the endpoint always returns the fixed \|rg\| &gt; 0.8 set. 默认 `0.05` |

??? note "上游说明（英文原文）"

    Get strong genetic correlations (\|rg\| &gt; 0.8) between a GWAS trait and other traits in the IEU OpenGWAS database, indicating shared genetic architecture. NOTE: the genetic-correlation graph is sparse and matches exact, case-sensitive trait labels, so many common traits (e.g. 'Body mass index') return no rows even though they exist in OpenGWAS — an empty result means 'no strong-rg edge in the graph', not 'no shared genetics'. Genetic correlation is NOT causal evidence; use EpiGraphDB_get_mendelian_randomization for causality. Example traits with data: 'Waist circumference', 'Hip circumference', 'Body fat percentage'.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EpiGraphDB_get_genetic_correlations',
     'arguments': {'trait': '...'}},
)
```

## `EpiGraphDB_get_mendelian_randomization` {#EpiGraphDB_get_mendelian_randomization}

Gets Mendelian randomization estimates between exposure and outcome traits from OpenGWAS.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`exposure_trait`, `outcome_trait`
- **关键词：**mendelian randomization, causal inference, gwas, opengwas
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `exposure_trait` | string | 是 | GWAS trait name for the exposure (e.g., 'Body mass index', 'LDL cholesterol', 'Type 2 diabetes'). Use EpiGraphDB_search_opengwas to find exact trait names. |
| `outcome_trait` | string | 是 | GWAS trait name for the outcome (e.g., 'Coronary heart disease', 'Stroke', 'Alzheimer disease'). Use EpiGraphDB_search_opengwas to find exact trait names. |
| `pval_threshold` | number | 否 | P-value threshold for instrument selection (default 1e-5). 默认 `1e-05` |

??? note "上游说明（英文原文）"

    Get Mendelian Randomization (MR) evidence for a causal relationship between an exposure trait and an outcome trait using IEU OpenGWAS data. Returns MR estimates (beta, SE, p-value) from multiple instruments/methods. Example: exposure='Body mass index', outcome='Coronary heart disease' returns causal effect estimates.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EpiGraphDB_get_mendelian_randomization',
     'arguments': {'exposure_trait': '...', 'outcome_trait': '...'}},
)
```

## `EpiGraphDB_search_opengwas` {#EpiGraphDB_search_opengwas}

Searches IEU OpenGWAS for GWAS studies by trait name, returning study IDs for MR analyses.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`query`
- **关键词：**opengwas, gwas study, trait, ieu
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Trait or phenotype name to search for in OpenGWAS (e.g., 'body mass index', 'coronary heart disease', 'type 2 diabetes', 'blood pressure'). |

??? note "上游说明（英文原文）"

    Search the IEU OpenGWAS database for GWAS studies by trait name. Returns matching GWAS study IDs and trait names with similarity scores. OpenGWAS IDs (e.g., 'ieu-a-2') are needed as inputs for Mendelian Randomization analyses. Example: 'body mass index' returns 'ieu-a-2', 'ieu-b-40', etc.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EpiGraphDB_search_opengwas',
     'arguments': {'query': '...'}},
)
```

## `eQTL_get_associations` {#eQTL_get_associations}

Gets eQTL Catalogue associations between variants and gene/exon/transcript expression for a dataset.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`dataset_id`
- **关键词：**eqtl, eqtl catalogue, molecular trait, splicing qtl
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `dataset_id` | string | 是 | Dataset identifier from eQTL_list_datasets (e.g., 'QTD000001', 'QTD000584') |
| `size` | integer | 否 | Number of associations to return (default 20) 默认 `20` |
| `gene_id` | string \| null | 否 | Filter by Ensembl gene ID (e.g., 'ENSG00000141510' for TP53). Returns all eQTLs affecting this gene. |
| `variant` | string \| null | 否 | Filter by variant in chr_pos_ref_alt format (e.g., 'chr1_791100_G_GGGA'). Returns all genes affected by this variant. |
| `molecular_trait_id` | string \| null | 否 | Filter by molecular trait ID (usually same as gene_id for 'ge' datasets) |

??? note "上游说明（英文原文）"

    Get eQTL associations (variant-gene expression correlations) from a specific dataset in the EBI eQTL Catalogue. Returns statistical associations between genetic variants and molecular traits (gene/exon/transcript expression or splicing). Each result includes beta effect size, p-value, allele frequency, and variant coordinates. Use eQTL_list_datasets first to find available dataset IDs. Can filter by gene_id (Ensembl gene ID) or variant (chr_pos_ref_alt format).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'eQTL_get_associations',
     'arguments': {'dataset_id': '...'}},
)
```

## `eQTL_list_datasets` {#eQTL_list_datasets}

Lists EBI eQTL Catalogue datasets by tissue and quantification method.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**eqtl catalogue, datasets, tissue, quantification
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `size` | integer | 否 | Number of datasets to return (default 20) 默认 `20` |
| `quant_method` | string \| null | 否 | Filter by quantification method: 'ge' (gene expression), 'exon', 'tx' (transcript), 'txrev' (transcript ratio), 'leafcutter' (splicing) |
| `tissue_id` | string \| null | 否 | Filter by tissue ontology ID (e.g., 'CL_0000235' for macrophage, 'UBERON_0002107' for liver) |

??? note "上游说明（英文原文）"

    List available eQTL (expression Quantitative Trait Loci) datasets from the EBI eQTL Catalogue. Each dataset represents a study-tissue-quantification combination with measured genetic associations to gene expression. Returns dataset IDs, study labels, tissue types, sample sizes, and quantification methods (ge=gene expression, exon, tx=transcript, txrev=transcript ratio, leafcutter=splicing). Use dataset IDs with eQTL_get_associations for specific variant-gene associations.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'eQTL_list_datasets', 'arguments': {}},
)
```

## `europepmc_disease_target_score` {#europepmc_disease_target_score}

Gets all Open Targets target-disease scores from Europe PMC text mining for a disease.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`efoId`
- **关键词：**open targets, europe pmc, text mining, target-disease
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `efoId` | string | 是 | The EFO (Experimental Factor Ontology) ID of the disease, e.g., 'MONDO_0011996' for chronic myelogenous leukemia |
| `pageSize` | integer | 否 | Upstream fetch batch size (max 100), NOT a limit on results returned. The tool pages through all of the disease's associated targets regardless of this value -- 'pageSize': 10 on MONDO_0005011 (Crohn disease) still returns 390 scored targets. Lower it only to make each upstream request smaller. 默认 `100` |

??? note "上游说明（英文原文）"

    Extract disease-target association scores from Europe PMC literature. This includes literature-based evidence. Returns EVERY target with a Europe PMC score for the disease, sorted strongest-first -- it is not a paged or top-N view. Scanning is capped at 25 seconds; when the cap is hit the response sets 'truncated' and a 'note' saying how many of the disease's associated targets were scanned.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'europepmc_disease_target_score',
     'arguments': {'efoId': '...'}},
)
```

## `eva_disease_target_score` {#eva_disease_target_score}

Gets Open Targets target-disease scores from EVA (ClinVar) germline variants for a disease.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`efoId`
- **关键词：**open targets, eva, clinvar, germline, target-disease
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `efoId` | string | 是 | The EFO (Experimental Factor Ontology) ID of the disease, e.g., 'MONDO_0011996' for chronic myelogenous leukemia |
| `pageSize` | integer | 否 | Number of results per page (default: 100, max: 100) 默认 `100` |

??? note "上游说明（英文原文）"

    Extract disease-target association scores from EVA (European Variation Archive). EVA provides genetic variant data.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'eva_disease_target_score', 'arguments': {'efoId': '...'}},
)
```

## `eva_somatic_disease_target_score` {#eva_somatic_disease_target_score}

Gets Open Targets target-disease scores from EVA somatic mutations for a disease.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`efoId`
- **关键词：**open targets, eva, somatic, target-disease
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `efoId` | string | 是 | The EFO (Experimental Factor Ontology) ID of the disease, e.g., 'MONDO_0011996' for chronic myelogenous leukemia |
| `pageSize` | integer | 否 | Number of results per page (default: 100, max: 100) 默认 `100` |

??? note "上游说明（英文原文）"

    Extract disease-target association scores from EVA somatic mutations. This includes somatic variant data.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'eva_somatic_disease_target_score',
     'arguments': {'efoId': '...'}},
)
```

## `expression_atlas_disease_target_score` {#expression_atlas_disease_target_score}

Gets Open Targets target-disease scores from Expression Atlas for a disease.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`efoId`
- **关键词：**open targets, expression atlas, target-disease
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `efoId` | string | 是 | The EFO (Experimental Factor Ontology) ID of the disease, e.g., 'MONDO_0011996' for chronic myelogenous leukemia |
| `pageSize` | integer | 否 | Number of results per page (default: 100, max: 100) 默认 `100` |

??? note "上游说明（英文原文）"

    Extract disease-target association scores from Expression Atlas. This provides gene expression data.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'expression_atlas_disease_target_score',
     'arguments': {'efoId': '...'}},
)
```

## `fda_pharmacogenomic_biomarkers` {#fda_pharmacogenomic_biomarkers}

Fetches the FDA Table of Pharmacogenomic Biomarkers in Drug Labeling, filterable by drug name or gene biomarker.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**pharmacogenomics, pgx, biomarker gene, cyp2d6, hla, drug labeling
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `drug_name` | string | 否 | Filter by the name of the drug (e.g., 'Sivextro', 'Abacavir'). Case-insensitive partial match. |
| `biomarker` | string | 否 | Filter by the specific biomarker (e.g., 'CYP2D6', 'HLA-B*5701'). Case-insensitive partial match. |
| `limit` | integer | 否 | Maximum number of results to return. 默认 `10` |

??? note "上游说明（英文原文）"

    Retrieve pharmacogenomic biomarkers from FDA drug labels. This tool fetches the Table of Pharmacogenomic Biomarkers in Drug Labeling from the FDA website and allows filtering by drug name or biomarker.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'fda_pharmacogenomic_biomarkers', 'arguments': {}},
)
```

## `Finemap_credible_set` {#Finemap_credible_set}

Computes per-SNP posterior inclusion probabilities via Wakefield approximate Bayes factors and a credible set from effect sizes and SEs.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`beta`, `se`
- **关键词：**fine-mapping, pip, credible set, abf, summary statistics, causal variant
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `beta` | array of number | 是 | Per-SNP effect sizes over the region. |
| `se` | array of number | 是 | Standard errors of beta (positive). |
| `snp` | array \| null | 否 | Optional SNP identifiers (same order) to label PIPs and the credible set. |
| `coverage` | number \| null | 否 | Credible-set coverage in (0,1\] (default 0.95). |
| `sd_prior` | number \| null | 否 | Prior SD of the effect size (default 0.15). |

??? note "上游说明（英文原文）"

    Single-causal-variant fine-mapping from GWAS/QTL summary statistics (per-SNP effect size + SE over a region). Computes each SNP's posterior inclusion probability (PIP) via Wakefield approximate Bayes factors and the 95% credible set (smallest SNP set capturing the causal variant). LD-free (assumes one causal signal); pure-compute, no R. For multi-signal loci use SuSiE/FINEMAP with LD.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Finemap_credible_set',
     'arguments': {'beta': [], 'se': []}},
)
```

## `FinnGen_get_phenotype` {#FinnGen_get_phenotype}

Gets a FinnGen disease endpoint's case/control counts, genomic control lambda and number of genome-wide significant loci.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`phenocode`
- **关键词：**finngen, gwas endpoint, phenocode, lambda gc, finnish biobank
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `phenocode` | string | 是 | FinnGen phenotype code (e.g., 'T2D', 'I9_CHD', 'C3_BREAST', 'J10_ASTHMA'). Get valid codes from FinnGen_list_phenotypes. |

??? note "上游说明（英文原文）"

    Get detailed information about a specific FinnGen disease endpoint (phenotype). Returns case/control counts, genomic control lambda values, and genome-wide significant loci count for a given phenocode. Use phenocodes from FinnGen_list_phenotypes (e.g., 'T2D' for type 2 diabetes, 'I9_CHD' for coronary heart disease, 'C3_BREAST' for breast cancer).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'FinnGen_get_phenotype',
     'arguments': {'phenocode': '...'}},
)
```

## `FinnGen_get_region_associations` {#FinnGen_get_region_associations}

Gets regional association and fine-mapping credible sets (lead variants, PIPs, MAF) for a FinnGen phenotype in a region.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`phenocode`, `region`
- **关键词：**finngen, gwas, credible set, locus, fine-mapping, maf
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `phenocode` | string | 是 | FinnGen phenotype code (e.g., 'T2D', 'I9_CHD'). Get codes from FinnGen_list_phenotypes. |
| `region` | string | 是 | Genomic region in chr:start-end format (GRCh38). Example: '9:22000000-22200000' for the CDKN2A/B locus on chromosome 9. |

??? note "上游说明（英文原文）"

    Get regional association and fine-mapping credible set data for a FinnGen phenotype in a specific genomic region. Returns credible sets with lead variants, posterior probabilities, and minor allele frequencies. Use to examine GWAS signals at a specific locus for a Finnish population disease endpoint.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'FinnGen_get_region_associations',
     'arguments': {'phenocode': '...', 'region': '...'}},
)
```

## `FinnGen_get_variant_finemapping` {#FinnGen_get_variant_finemapping}

Returns FinnGen fine-mapping regions and phenotypes in which a GRCh38 chr:pos:ref:alt variant appears in a credible set.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`variant`
- **关键词：**finngen, fine-mapping, credible set, variant phenotype, grch38
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant` | string | 是 | Genomic variant in chr:pos:ref:alt format (GRCh38). Examples: '19:44908684:T:C' (APOE rs429358), '10:112998590:C:T' (TCF7L2 rs7903146). Also accepts chr-pos-ref-alt format. |

??? note "上游说明（英文原文）"

    Get fine-mapping regions associated with a genomic variant in FinnGen. Returns phenotypes where this variant appears in a credible set from statistical fine-mapping. Input variant in chr:pos:ref:alt format (GRCh38). Useful for understanding which disease associations a variant contributes to in the Finnish population.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'FinnGen_get_variant_finemapping',
     'arguments': {'variant': '...'}},
)
```

## `FinnGen_list_phenotypes` {#FinnGen_list_phenotypes}

Searches/lists FinnGen disease endpoints by keyword or category with case/control counts and genome-wide significant loci.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**finngen, gwas endpoints, phenocode, biobank, disease phenotypes
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string \| null | 否 | Search query to filter phenotypes by code, name, or category (e.g., 'diabetes', 'breast cancer', 'asthma'). Case-insensitive. |
| `category` | string \| null | 否 | Filter by disease category (e.g., 'Diabetes', 'circulatory', 'Neoplasms'). Case-insensitive substring match. |
| `min_cases` | integer \| null | 否 | Minimum number of cases required. Use to filter for well-powered phenotypes (e.g., 1000). |
| `limit` | integer \| null | 否 | Maximum number of results to return (default 50). |

??? note "上游说明（英文原文）"

    Search and list FinnGen disease endpoints (phenotypes). FinnGen is a Finnish population genomics study with 486,367 participants and 2,470 phenotypes. Search by keyword (e.g., 'diabetes', 'asthma', 'cancer') or browse by category. Returns phenotype codes, case/control counts, and number of genome-wide significant loci. Use phenocode values from results with FinnGen_get_phenotype and FinnGen_get_region_associations.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'FinnGen_list_phenotypes', 'arguments': {}},
)
```

## `FlyBase_get_gene_disease_models` {#FlyBase_get_gene_disease_models}

Returns genotype-based Drosophila disease models for a gene with associated human disease names and DOIDs.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`gene_id`
- **关键词：**drosophila, disease model, doid, model organism, human disease
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_id` | string | 是 | FlyBase gene ID with 'FB:' prefix. Examples: 'FB:FBgn0000490' (dpp), 'FB:FBgn0004647' (Notch/N), 'FB:FBgn0003996' (white/w). |
| `limit` | integer | 否 | Maximum number of disease models to return (1-100). Default: 20. |
| `page` | integer | 否 | Page number for pagination. Default: 1. |

??? note "上游说明（英文原文）"

    Get disease models involving a Drosophila gene from the Alliance of Genome Resources. Drosophila is a powerful model organism for human disease research. Returns genotype-based disease models with associated human disease names and Disease Ontology IDs. Example: 'FB:FBgn0000490' (dpp) has 946 disease model entries including models for ALS, Parkinson's disease, and other neurodegenerative conditions.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'FlyBase_get_gene_disease_models',
     'arguments': {'gene_id': '...'}},
)
```

## `G2P_download_panel` {#G2P_download_panel}

Downloads all curated gene-disease entries of a Gene2Phenotype panel with genotype, mechanism and evidence.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`panel`
- **关键词：**gene2phenotype, g2p panel, developmental disorders, allelic requirement, curation
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `panel` | string | 是 | Panel name (case-insensitive). Examples: 'Cancer', 'DD' (developmental disorders), 'Skeletal', 'Eye', 'Cardiac', 'Ear', 'Skin'. |

??? note "上游说明（英文原文）"

    Download the FULL content of a Gene2Phenotype (G2P) clinical panel as structured gene-disease records. Unlike G2P_get_panel (which returns only summary counts), this returns EVERY curated gene-disease entry in the panel with: G2P stable id, gene symbol, HGNC id, gene MIM, disease name, disease MIM, disease MONDO id, allelic requirement (inheritance), confidence (definitive/strong/moderate/limited), variant consequence, variant types, molecular mechanism, phenotypes (HPA/HP terms), supporting publications (PubMed IDs), panels, and review dates. Use this to enumerate all genes in a panel or to build clinical gene lists. Panels: DD (developmental disorders, ~2856 records), Cancer (~125 records), Skeletal, Eye, Cardiac, Ear, Skin, Hearing loss, Mitochondrial.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'G2P_download_panel', 'arguments': {'panel': '...'}},
)
```

## `G2P_get_panel` {#G2P_get_panel}

Returns summary statistics (counts, description, update date) for a Gene2Phenotype clinical panel.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`panel`
- **关键词：**gene2phenotype, clinical panel, dd panel, gene-disease
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `panel` | string | 是 | Panel name: 'cancer', 'dd' (developmental disorders), 'skeletal', 'eye', 'cardiac', 'ear', 'skin', 'prenatal' |

??? note "上游说明（英文原文）"

    Get summary statistics for a Gene2Phenotype clinical panel. G2P organizes gene-disease associations into expert-curated panels: DD (developmental disorders), Cancer, Skeletal, Eye, Cardiac, Ear, Skin, and Prenatal. Returns panel name, description, last update date, total records count, total genes, and breakdown by confidence level (definitive, strong, moderate, limited).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'G2P_get_panel', 'arguments': {'panel': '...'}},
)
```

## `G2P_get_record` {#G2P_get_record}

Gets a G2P gene-disease record by stable ID with locus, inheritance, variant consequences and evidence.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`stable_id`
- **关键词：**gene2phenotype, inheritance, variant consequence, mechanism, clinical genetics
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `stable_id` | string | 是 | G2P stable identifier (e.g., 'G2P01796', 'G2P00119', 'G2P00979'). Get IDs from G2P_search results. |

??? note "上游说明（英文原文）"

    Get detailed Gene2Phenotype (G2P) record for a specific gene-disease association by its stable ID. Returns comprehensive clinical genetics information including locus details (gene symbol, genomic coordinates, cross-references), genotype/inheritance pattern, variant consequences (with Sequence Ontology terms), molecular mechanism with evidence, disease details with ontology terms, confidence level, and supporting publications. Use G2P_search first to find stable IDs.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'G2P_get_record', 'arguments': {'stable_id': '...'}},
)
```

## `G2P_search` {#G2P_search}

Searches Gene2Phenotype for curated gene-disease pairs with genotype and confidence by gene or disease.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`query`
- **关键词：**gene2phenotype, gene-disease, monoallelic, biallelic, confidence
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Search query - gene symbol (e.g., 'BRCA1', 'TP53') or disease name (e.g., 'epilepsy', 'cardiomyopathy', 'cancer') |
| `page` | integer \| null | 否 | Page number for pagination (default 1) |

??? note "上游说明（英文原文）"

    Search the Gene2Phenotype (G2P) database for gene-disease associations curated by clinical geneticists. G2P is maintained by EBI and provides expert-curated evidence for gene-disease relationships used in clinical diagnosis. Returns gene-disease pairs with genotype (monoallelic/biallelic), molecular mechanism (loss/gain of function), clinical panel assignment (DD, Cancer, Skeletal, etc.), and confidence level (definitive, strong, moderate, limited). Useful for clinical variant interpretation and genetic diagnosis.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'G2P_search', 'arguments': {'query': '...'}},
)
```

## `gather_disease_profile` {#gather_disease_profile}

Aggregates a disease profile from Orphanet, OMIM, DisGeNET, Open Targets and OLS, including IDs, associated genes, inheritance and prevalence.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`disease`
- **关键词：**disease profile, orphanet, omim, disease genes, inheritance, mondo
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `disease` | string | 是 | Disease name (e.g., 'breast cancer', 'cystic fibrosis', 'distal renal tubular acidosis') |

??? note "上游说明（英文原文）"

    Gather a comprehensive disease profile from Orphanet, OMIM, DisGeNET, OpenTargets, and OLS in a single call. Returns disease identifiers (ORPHA, MIM, EFO, MONDO), associated genes with scores, inheritance pattern, prevalence, and ontology cross-references. Use for disease research requiring multi-database evidence.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gather_disease_profile', 'arguments': {'disease': '...'}},
)
```

## `gather_gene_disease_associations` {#gather_gene_disease_associations}

Queries DisGeNET, OMIM, Open Targets, GenCC and ClinVar for gene-disease associations and reports cross-source concordance.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**gene-disease, disgenet, omim, gencc, clinvar, concordance
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string \| null | 否 | Gene symbol (e.g., 'TP53', 'BRCA1'). At least one of gene or disease required. |
| `disease` | string \| null | 否 | Disease name (e.g., 'breast cancer', 'distal renal tubular acidosis'). At least one of gene or disease required. |

??? note "上游说明（英文原文）"

    Query multiple gene-disease databases (DisGeNET, OMIM, OpenTargets, GenCC, ClinVar) in a single call. Returns a cross-referenced comparison table showing which sources agree on each association, with concordance scores. Use when you need to compare gene-disease evidence across databases or find associations unique to one source. Each association carries 'concordance' (how many sources named it) with two denominators: 'total_sources_queried' (all five, including any that failed) and 'total_sources_with_data' (those that actually contributed) -- use the latter, since a source with no API key configured otherwise looks like a source that disagreed. 'associations' is capped at 50 and 'per_source_results' at 10 rows per source; 'num_associations', 'truncated' and 'per_source_result_counts' give the untruncated figures and 'notes' names any truncation that applied.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gather_gene_disease_associations', 'arguments': {}},
)
```

## `GenCC_get_classifications` {#GenCC_get_classifications}

Returns counts of GenCC gene-disease validity classifications, filterable by gene, disease or submitter.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**gencc, gene-disease validity, clingen, curation, definitive
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type (fixed: get_classifications) 默认 `get_classifications` |
| `gene_symbol` | string | 否 | Optional HGNC gene symbol to filter classification stats (e.g., 'GBA', 'BRCA2') |
| `disease` | string | 否 | Optional disease name substring to filter classification stats (e.g., 'Gaucher', 'breast cancer') |
| `submitter` | string | 否 | Optional filter by submitting organization name (e.g., 'ClinGen', 'Ambry', 'Genomics England') 默认 `` |

??? note "上游说明（英文原文）"

    Get summary statistics of gene-disease validity classifications from GenCC. Returns counts per classification level (Definitive, Strong, Moderate, Limited, etc.). Optionally filter by gene_symbol (e.g., 'GBA'), disease name substring, or submitter organization. Without filters, returns global stats across all 40,000+ curated submissions.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GenCC_get_classifications', 'arguments': {}},
)
```

## `GenCC_search_disease` {#GenCC_search_disease}

Finds genes with curated gene-disease validity classifications for a disease name or MONDO/OMIM ID.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**gencc, disease genes, mondo, omim, clinical validity
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type (fixed: search_disease) 默认 `search_disease` |
| `disease` | string | 否 | Disease name or ID to search (e.g., 'Marfan syndrome', 'MONDO:0007947', 'breast cancer') |
| `disease_name` | string | 否 | Alias for disease. Disease name or ID to search. |
| `classification` | string | 否 | Optional filter by classification level (e.g., Definitive, Strong) 默认 `` |

??? note "上游说明（英文原文）"

    Find genes with validity evidence for a disease from GenCC (Gene Curation Coalition). Search by disease name or MONDO/OMIM identifier. Returns gene-disease validity classifications from multiple expert curators. Useful for identifying established gene-disease relationships for clinical interpretation.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GenCC_search_disease', 'arguments': {}},
)
```

## `GenCC_search_gene` {#GenCC_search_gene}

Gets curated gene-disease validity classifications (Definitive to Refuted) for a gene from GenCC submitters.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`gene_symbol`
- **关键词：**gencc, gene-disease validity, clingen, panelapp, orphanet
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type (fixed: search_gene) 默认 `search_gene` |
| `gene_symbol` | string | 是 | HGNC gene symbol (e.g., BRCA2, TP53, FBN1, CFTR) |
| `classification` | string | 否 | Optional filter by classification level (e.g., Definitive, Strong, Moderate, Limited) 默认 `` |

??? note "上游说明（英文原文）"

    Get gene-disease validity classifications for a gene from GenCC (Gene Curation Coalition). Returns curated evidence levels (Definitive, Strong, Moderate, Limited, Disputed, Refuted) from multiple expert panels including ClinGen, Ambry Genetics, Genomics England PanelApp, and Orphanet. Essential for assessing clinical validity of gene-disease associations in genetic testing.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GenCC_search_gene', 'arguments': {'gene_symbol': '...'}},
)
```

## `Genebass_gene_burden_phewas` {#Genebass_gene_burden_phewas}

Returns rare-variant gene burden associations across phenotypes from Genebass (UK Biobank exomes).

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`gene`
- **关键词：**gene burden, exome, saige-gene, phewas
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 是 | Ensembl gene ID (e.g. 'ENSG00000148737') or gene symbol (e.g. 'TCF7L2'). Symbols are auto-resolved via Ensembl. |
| `burden_set` | string | 否 | Variant category aggregated in the burden test. One of: 'pLoF' (default), 'missense\|LC', 'synonymous'. Aliases accepted: lof/plof -&gt; pLoF, missense -&gt; missense\|LC. 默认 `pLoF` |
| `limit` | integer | 否 | Maximum number of phenotype associations to return (sorted by ascending burden p-value). 默认 `25` |
| `max_pval` | number \| null | 否 | Optional: only return associations with burden p-value &lt;= this threshold (e.g. 2.5e-6 for exome-wide significance). |

??? note "上游说明（英文原文）"

    Exome-wide gene-level burden PheWAS from Genebass (UK Biobank ~394K exomes, SAIGE-GENE+). Given one gene, returns its rare-variant burden association across ~4,500 phenotypes (GRCh38), sorted by p-value. Aggregates rare coding variants into a single gene-level test, complementing single-variant PheWAS. Accepts an Ensembl gene ID (ENSG...) or a gene symbol (auto-resolved via Ensembl). Choose a burden_set: pLoF (predicted loss-of-function, default), missense\|LC, or synonymous.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Genebass_gene_burden_phewas',
     'arguments': {'gene': '...'}},
)
```

## `GeneNetwork_get_trait_info` {#GeneNetwork_get_trait_info}

Gets trait (probe/gene) details in a GeneNetwork dataset: gene symbol, position, mean expression, LRS, locus and p-value.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`dataset_name`, `trait_name`
- **关键词：**genenetwork, qtl, lrs, eqtl, probe, linkage
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `dataset_name` | string | 是 | Dataset name/abbreviation. Examples: 'HC_M2_0606_P' (Hippocampus Consortium M430v2), 'BR_M2_1106_R' (Brain mRNA). |
| `trait_name` | string | 是 | Trait/probe identifier. Examples: '1436869_at' (Shh probe), '1457211_at', '1418701_at'. These are microarray probe IDs for expression datasets. |

??? note "上游说明（英文原文）"

    Get detailed information about a specific trait (probe/gene) in a GeneNetwork dataset. Returns trait metadata including gene symbol, description, chromosome, megabase position, mean expression value, linkage statistics (LRS score, locus), additive effect, and p-value. Requires a dataset name and trait name (probe ID). Use GeneNetwork_list_datasets to find available datasets and their names. Common datasets include BXD hippocampus, liver, and brain expression arrays.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GeneNetwork_get_trait_info',
     'arguments': {'dataset_name': '...', 'trait_name': '...'}},
)
```

## `genomics_england_disease_target_score` {#genomics_england_disease_target_score}

Gets Open Targets target-disease scores from Genomics England PanelApp for a disease.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`efoId`
- **关键词：**open targets, genomics england, panelapp, target-disease
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `efoId` | string | 是 | The EFO (Experimental Factor Ontology) ID of the disease, e.g., 'MONDO_0011996' for chronic myelogenous leukemia |
| `pageSize` | integer | 否 | Number of results per page (default: 100, max: 100) 默认 `100` |

??? note "上游说明（英文原文）"

    Extract disease-target association scores from Genomics England data. This includes clinical genomics evidence.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'genomics_england_disease_target_score',
     'arguments': {'efoId': '...'}},
)
```

## `GTEx_calculate_eqtl` {#GTEx_calculate_eqtl}

Calculates eQTL statistics (p-value, NES, t-stat, genotypes) for any gene-variant pair in a GTEx tissue.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`gencode_id`, `variant_id`, `tissue_site_detail_id`
- **关键词：**eqtl, nes, genotype, gtex, association test
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`calculate_eqtl` 默认 `calculate_eqtl` |
| `gencode_id` | string | 是 | Required: Versioned GENCODE ID (e.g. ENSG00000141510.18) |
| `variant_id` | string | 是 | Required: GTEx variant ID |
| `tissue_site_detail_id` | string | 是 | Required: Tissue ID (e.g. 'Whole_Blood') |
| `dataset_id` | string | 否 | 可选值：`gtex_v8`, `gtex_v10`, `gtex_snrnaseq_pilot` 默认 `gtex_v8` |

??? note "上游说明（英文原文）"

    Calculate custom eQTL for any gene-variant pair in any tissue. Dynamically calculates gene-variant association statistics including p-value, normalized effect size (NES), t-statistic, and genotype data. Works for both significant and non-significant associations. Returns detailed statistics even if association not in precomputed results. Use for: testing specific hypotheses, custom variant-gene pairs, dynamic eQTL analysis, validating findings. Requires all three: gencode_id, variant_id, and tissue_site_detail_id.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GTEx_calculate_eqtl',
     'arguments': {'gencode_id': '...',
                   'variant_id': '...',
                   'tissue_site_detail_id': '...'}},
)
```

## `GTEx_get_eqtl_genes` {#GTEx_get_eqtl_genes}

Gets genes with significant cis-eQTLs (eGenes) per tissue with allelic fold change and q-values.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**egenes, cis-eqtl, gtex, afc, fdr
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`get_eqtl_genes` 默认 `get_eqtl_genes` |
| `tissue_site_detail_id` | array of string | 否 | Optional: Filter by tissue IDs. Omit for all tissues |
| `dataset_id` | string | 否 | GTEx dataset version 可选值：`gtex_v8`, `gtex_v10`, `gtex_snrnaseq_pilot` 默认 `gtex_v8` |
| `page` | integer | 否 | 默认 `0` |
| `items_per_page` | integer | 否 | 默认 `250` |

??? note "上游说明（英文原文）"

    Get eQTL genes (eGenes) with significant cis-eQTLs. Returns genes with at least one significant expression quantitative trait locus (eQTL) in specified tissues. Includes allelic fold change, p-values, q-values (FDR &lt; 0.05), and effect sizes. Filter by tissue or get all eGenes. Use for: identifying genes with genetic regulation, finding tissue-specific eQTLs, studying regulatory variants, genetic association studies. Essential for understanding genetic control of gene expression.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GTEx_get_eqtl_genes', 'arguments': {}},
)
```

## `GTEx_get_finemapping_and_independent_eqtl` {#GTEx_get_finemapping_and_independent_eqtl}

Gets GTEx DAP-G fine-mapping credible sets or conditionally independent eQTLs for a gene and tissue.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**fine-mapping, dap-g, independent eqtl, credible set, gtex
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`get_finemapping_and_independent_eqtl` 默认 `get_finemapping_and_independent_eqtl` |
| `result_type` | string | 否 | 'finemapping' for DAP-G credible sets with PIP (default), 'independent' for rank-ordered conditionally-independent eQTL signals. 可选值：`finemapping`, `independent` 默认 `finemapping` |
| `gencode_id` | string \| array | 否 | Gene identifier: gene symbol (e.g. 'ERAP2'), Ensembl ID, or versioned GENCODE ID (e.g. 'ENSG00000164308.16'). Auto-resolved. |
| `gene_symbol` | string | 否 | Gene symbol alias for gencode_id. |
| `tissue_site_detail_id` | array of string | 否 | Optional tissue ID(s) to filter (e.g. \['Adipose_Subcutaneous'\]). Omit for all tissues. |
| `page` | integer | 否 | Page number (0-based) 默认 `0` |
| `items_per_page` | integer | 否 | Results per page 默认 `250` |
| `dataset_id` | string | 否 | GTEx dataset version (default gtex_v8) 可选值：`gtex_v8`, `gtex_v10`, `gtex_snrnaseq_pilot` 默认 `gtex_v8` |

??? note "上游说明（英文原文）"

    Get GTEx statistical fine-mapping credible sets and conditionally-independent eQTLs - identifying likely causal variants and multiple independent regulatory signals per gene beyond the marginal single-tissue eQTL list. result_type='finemapping' (default) returns DAP-G credible sets with per-variant PIP (posterior inclusion probability) and setId. result_type='independent' returns rank-ordered conditionally-independent (secondary) eQTL signals. Based on gtex_v8. Example: ERAP2 (ENSG00000164308.16) Adipose_Subcutaneous -&gt; method DAP-G, variant chr5_96916728_G_A_b38, pip 0.9468, setId 1; independent -&gt; rs2927608 rank 1 pValue 2.88e-123 nes 1.076.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GTEx_get_finemapping_and_independent_eqtl',
     'arguments': {}},
)
```

## `GTEx_get_multi_tissue_eqtls` {#GTEx_get_multi_tissue_eqtls}

Gets Metasoft multi-tissue eQTL results (m-values, NES, p-values) for a gene-variant pair.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`gencode_id`
- **关键词：**multi-tissue eqtl, metasoft, m-value, gtex
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`get_multi_tissue_eqtls` 默认 `get_multi_tissue_eqtls` |
| `gencode_id` | string | 是 | Required: Versioned GENCODE ID |
| `variant_id` | string | 否 | Optional: GTEx variant ID to filter specific variant |
| `dataset_id` | string | 否 | 可选值：`gtex_v8`, `gtex_v10`, `gtex_snrnaseq_pilot` 默认 `gtex_v8` |
| `page` | integer | 否 | 默认 `0` |
| `items_per_page` | integer | 否 | 默认 `250` |

??? note "上游说明（英文原文）"

    Get multi-tissue eQTL meta-analysis results (Metasoft). Returns m-values (posterior probability of eQTL effect in each tissue), normalized effect sizes, p-values, and standard errors across all tissues for a gene-variant pair. Shows tissue-sharing patterns of eQTL effects. Use for: understanding tissue-specificity of eQTLs, finding shared vs tissue-specific effects, cross-tissue meta-analysis, regulatory mechanism studies. Requires gencode_id, optionally filter by variant_id.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GTEx_get_multi_tissue_eqtls',
     'arguments': {'gencode_id': '...'}},
)
```

## `GTEx_get_single_tissue_eqtls` {#GTEx_get_single_tissue_eqtls}

Gets significant single-tissue gene-variant eQTL associations with p-values and NES.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**eqtl, gtex, variant, tissue, nes
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`get_single_tissue_eqtls` 默认 `get_single_tissue_eqtls` |
| `gencode_id` | array of string | 否 | Optional: Versioned GENCODE ID(s) to query |
| `variant_id` | array of string | 否 | Optional: GTEx variant ID(s) to query |
| `tissue_site_detail_id` | array of string | 否 | Optional: Tissue ID(s) to filter. At least one of gencode_id, variant_id, or tissue_site_detail_id required |
| `dataset_id` | string | 否 | 可选值：`gtex_v8`, `gtex_v10`, `gtex_snrnaseq_pilot` 默认 `gtex_v8` |
| `page` | integer | 否 | 默认 `0` |
| `items_per_page` | integer | 否 | 默认 `250` |

??? note "上游说明（英文原文）"

    Get significant single-tissue eQTL associations. Returns precomputed gene-variant associations with p-values, effect sizes (NES), and genomic positions. Query by gene ID, variant ID, tissue, or combination. Results show which genetic variants affect gene expression in specific tissues. Use for: finding eQTLs for candidate genes, looking up variant effects, tissue-specific regulatory analysis, GWAS follow-up. Must provide at least one filter: gencode_id, variant_id, or tissue_site_detail_id.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GTEx_get_single_tissue_eqtls', 'arguments': {}},
)
```

## `GTEx_get_single_tissue_sqtls` {#GTEx_get_single_tissue_sqtls}

Gets GTEx single-tissue splicing QTLs or sGenes for a gene/variant and tissue.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**sqtl, splicing qtl, sgenes, alternative splicing, gtex
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`get_single_tissue_sqtls` 默认 `get_single_tissue_sqtls` |
| `result_type` | string | 否 | 'sqtl' for single-tissue sQTL associations (default), 'sgene' for genes with significant sQTLs in a tissue. 可选值：`sqtl`, `sgene` 默认 `sqtl` |
| `gencode_id` | string \| array | 否 | Gene identifier(s) for sqtl mode: gene symbol (e.g. 'ERAP2'), Ensembl ID, or versioned GENCODE ID (e.g. 'ENSG00000164308.16'). Auto-resolved. |
| `gene_symbol` | string | 否 | Gene symbol alias for gencode_id. |
| `variant_id` | string \| array | 否 | Optional GTEx variant ID(s) to filter sqtl results (e.g. 'chr5_96659855_A_G_b38'). |
| `tissue_site_detail_id` | array of string | 否 | Tissue ID(s) (e.g. \['Whole_Blood'\]). Recommended for both modes. See GTEx_get_tissue_sites. |
| `page` | integer | 否 | Page number (0-based) 默认 `0` |
| `items_per_page` | integer | 否 | Results per page 默认 `250` |
| `dataset_id` | string | 否 | GTEx dataset version (default gtex_v8) 可选值：`gtex_v8`, `gtex_v10`, `gtex_snrnaseq_pilot` 默认 `gtex_v8` |

??? note "上游说明（英文原文）"

    Get GTEx single-tissue splicing QTLs (sQTLs) and sGenes - a molecular layer distinct from eQTLs that captures genetic regulation of alternative splicing. result_type='sqtl' (default) returns significant single-tissue sQTL associations for a gene/tissue (variant, pValue, nes, and phenotypeId which encodes the LeafCutter intron-excision cluster, e.g. chr5:96900189:96901506:clu_31758). result_type='sgene' returns genes with at least one significant sQTL in a tissue, with qValue. Based on gtex_v8. Example: ERAP2 (ENSG00000164308.16) in Whole_Blood -&gt; rs72772072 / chr5_96659855_A_G_b38, pValue 1.2e-06, nes 0.439.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GTEx_get_single_tissue_sqtls', 'arguments': {}},
)
```

## `GTEx_query_eqtl` {#GTEx_query_eqtl}

Queries GTEx single-tissue eQTL associations (variant, p-value, slope) for a gene.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**eqtl, gtex, regulatory variant, tissue
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 否 | Gene symbol (e.g., 'TP53', 'BRCA1'). Auto-resolved to versioned GENCODE ID. |
| `ensembl_gene_id` | string | 否 | Ensembl gene identifier (e.g., 'ENSG00000141510'). Use gene_symbol instead if you have a gene name. |
| `gene_id` | string | 否 | Alias for ensembl_gene_id. Ensembl gene ID (e.g., 'ENSG00000141510'). |
| `gene_input` | string | 否 | Alias for gene_symbol. Gene symbol or Ensembl ID (e.g., 'TP53', 'ENSG00000141510'). |
| `tissue_id` | string | 否 | GTEx tissue ID to filter eQTL results (e.g., 'Brain_Cortex', 'Whole_Blood', 'Liver'). Case-sensitive. |
| `tissue` | string | 否 | Alias for tissue_id. GTEx tissue name (e.g., 'Brain_Cortex', 'Whole_Blood'). |
| `page` | integer | 否 | Page number (1-based). 默认 `1` |
| `size` | integer | 否 | Number of records per page (1–100). 默认 `10` |

??? note "上游说明（英文原文）"

    Query GTEx single-tissue eQTL associations for a gene. Accepts gene symbols (TP53, BRCA1) or Ensembl IDs. Use to identify regulatory variants (variantId, pValue, slope) relevant to expression regulation.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GTEx_query_eqtl', 'arguments': {}},
)
```

## `GtoPdb_get_disease_associations` {#GtoPdb_get_disease_associations}

Gets GtoPdb curated disease-target and disease-ligand associations for a disease ID.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**gtopdb, drug target, iuphar, disease target
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `disease_id` | integer \| null | 否 | GtoPdb disease ID. Get from GtoPdb_search_diseases. Examples: 1161 (non-allergic asthma), 34 (acute myeloid leukemia). Alias: diseaseId also accepted. |
| `diseaseId` | integer \| null | 否 | Alias for disease_id. GtoPdb disease ID. |

??? note "上游说明（英文原文）"

    Get curated disease-to-target and disease-to-ligand associations from the Guide to Pharmacology (GtoPdb / IUPHAR-BPS) for a given disease ID. Merges two keyless endpoints into one record: diseaseTargets (drug targets linked to the disease, with roles, references, and ligand-target interactions) and diseaseLigands (approved/experimental ligands linked to the disease, with approval status and clinical-use notes). An empty list is a valid result (some diseases have targets but no curated ligands). Use GtoPdb_search_diseases first to find a disease_id (e.g., non-allergic asthma = 1161). No authentication required.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GtoPdb_get_disease_associations', 'arguments': {}},
)
```

## `gwas_get_association_by_id` {#gwas_get_association_by_id}

Gets a specific GWAS Catalog association by its identifier.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`association_id`
- **关键词：**gwas catalog, association id, risk allele
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `association_id` | string | 是 | GWAS association identifier |

??? note "上游说明（英文原文）"

    Get a specific GWAS association by its unique identifier.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gwas_get_association_by_id',
     'arguments': {'association_id': '...'}},
)
```

## `gwas_get_associations_for_snp` {#gwas_get_associations_for_snp}

Gets all GWAS Catalog associations for an rsID with optional sorting.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`rs_id`
- **关键词：**rsid, gwas catalog, snp associations, traits
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `rs_id` | string | 是 | dbSNP rs identifier |
| `sort` | string | 否 | Sort field (e.g., 'p_value', 'or_value') |
| `direction` | string | 否 | Sort direction ('asc' or 'desc') |
| `size` | integer | 否 | Number of results to return per page |
| `page` | integer | 否 | Page number for pagination |

??? note "上游说明（英文原文）"

    Get all associations for a specific SNP with optional sorting.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gwas_get_associations_for_snp',
     'arguments': {'rs_id': '...'}},
)
```

## `gwas_get_associations_for_study` {#gwas_get_associations_for_study}

Gets all associations for a GWAS Catalog study accession sorted by p-value.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`accession_id`
- **关键词：**gcst, gwas study, associations, p-value
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `accession_id` | string | 是 | Study accession identifier |
| `size` | integer | 否 | Number of results to return per page |
| `page` | integer | 否 | Page number for pagination |

??? note "上游说明（英文原文）"

    Get all associations for a specific study, sorted by p-value.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gwas_get_associations_for_study',
     'arguments': {'accession_id': '...'}},
)
```

## `gwas_get_associations_for_trait` {#gwas_get_associations_for_trait}

Gets GWAS associations for a trait (resolved to EFO/MONDO) sorted by p-value.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**efo, trait, gwas catalog, disease associations
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `disease_trait` | string | 否 | Disease or trait name for text-based search (e.g., 'diabetes', 'hypertension') |
| `efo_uri` | string | 否 | Full EFO ontology URI (e.g., 'http://www.ebi.ac.uk/efo/EFO_0001645') |
| `efo_id` | string | 否 | EFO/OBA term ID (e.g., 'EFO_0001645', 'OBA_2050062'). Recommended for reliable trait filtering. |
| `efo_trait` | string | 否 | Exact EFO trait label. Use when you know the canonical trait string. |
| `size` | integer | 否 | Number of results to return per page |
| `page` | integer | 否 | Page number for pagination |

??? note "上游说明（英文原文）"

    Get all associations for a specific trait, sorted by p-value (most significant first). A free-text disease_trait is resolved to an EFO/MONDO term first, and that resolution is disclosed in the response: 'query_trait' (what you asked for), 'resolved_efo_id', 'resolved_efo_label' (the term's label in GWAS Catalog) and 'trait_resolution_source', plus a 'trait_resolution_note' whenever the resolved label is not an exact match for your query. Check those before trusting the results: e.g. 'cardiorespiratory fitness' resolves to EFO_0009184 'heart rate response to exercise', which is a different phenotype from VO2max (EFO_0004887 'maximal oxygen uptake measurement'). Pass efo_id explicitly to bypass the resolver.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gwas_get_associations_for_trait', 'arguments': {}},
)
```

## `gwas_get_snp_by_id` {#gwas_get_snp_by_id}

Gets a GWAS Catalog SNP record by rsID with location and mapped genes.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`rs_id`
- **关键词：**rsid, snp, gwas catalog, mapped gene
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `rs_id` | string | 是 | dbSNP rs identifier |

??? note "上游说明（英文原文）"

    Get a specific GWAS SNP by its rs ID.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gwas_get_snp_by_id', 'arguments': {'rs_id': '...'}},
)
```

## `gwas_get_snps_for_gene` {#gwas_get_snps_for_gene}

Gets all GWAS Catalog SNPs mapped to a gene with rsIDs, positions and functional annotations.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**gwas catalog, snps, mapped gene, rsid
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `size` | integer | 否 | Number of results to return per page |
| `page` | integer | 否 | Page number for pagination |
| `gene_symbol` | string | 否 | Gene name or symbol (e.g., 'BRCA1', 'TP53', 'GBA') |
| `gene` | string | 否 | Alias for gene_symbol |
| `mapped_gene` | string | 否 | Alias for gene_symbol |

??? note "上游说明（英文原文）"

    Get all GWAS Catalog SNPs mapped to a specific gene by gene name (e.g., 'BRCA1', 'TP53'). Returns variant details including rs IDs, chromosomal positions, and functional annotations. Use this when you want to find all catalogued variants within or near a specific gene.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gwas_get_snps_for_gene', 'arguments': {}},
)
```

## `gwas_get_studies_for_trait` {#gwas_get_studies_for_trait}

Gets GWAS studies for a trait with cohort, GxE and summary-statistics filters.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**gwas studies, efo, cohort, gxe, summary statistics
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `disease_trait` | string | 否 | Disease or trait name for text-based search (e.g., 'diabetes', 'alzheimer disease') |
| `efo_uri` | string | 否 | Full EFO ontology URI (e.g., 'http://www.ebi.ac.uk/efo/EFO_0001645') |
| `efo_id` | string | 否 | EFO/OBA term ID (e.g., 'EFO_0001645', 'OBA_2050062'). Recommended for reliable trait filtering. |
| `efo_trait` | string | 否 | Exact EFO trait label. Use when you know the canonical trait string. |
| `cohort` | string | 否 | Cohort name (e.g., 'UKB' for UK Biobank) |
| `gxe` | boolean | 否 | Filter for Gene-by-Environment interaction studies |
| `full_pvalue_set` | boolean | 否 | Filter for studies with full summary statistics |
| `size` | integer | 否 | Number of results to return per page |
| `page` | integer | 否 | Page number for pagination |

??? note "上游说明（英文原文）"

    Get studies for a specific trait with optional filters for cohort, GxE interactions, and summary statistics.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gwas_get_studies_for_trait', 'arguments': {}},
)
```

## `gwas_get_study_by_id` {#gwas_get_study_by_id}

Gets a GWAS Catalog study by its GCST accession.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**gcst, gwas study, study metadata
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `study_id` | string | 否 | GWAS study identifier (GCST accession, e.g. 'GCST000392') |
| `accession_id` | string | 否 | Alias for study_id -- the same GCST accession identifier, accepted so results from other gwas_* tools (which use the field name 'accession_id') can be passed straight in. |

??? note "上游说明（英文原文）"

    Get a specific GWAS study by its unique identifier. Accepts either 'study_id' or 'accession_id' (the same GCST identifier field name used by GWASAssociationsForStudy and returned on every association record) so IDs can be chained directly between gwas_* tools without renaming.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gwas_get_study_by_id', 'arguments': {}},
)
```

## `gwas_get_variants_for_trait` {#gwas_get_variants_for_trait}

Searches the GWAS Catalog for variants associated with a disease/trait with rsIDs and p-values.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**gwas catalog, trait variants, rsid, efo
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `disease_trait` | string | 否 | Disease or trait name for text-based search (e.g., 'diabetes', 'breast cancer') |
| `trait` | string | 否 | Alias for disease_trait. Disease or trait name (e.g., 'type 2 diabetes', 'breast cancer') |
| `efo_uri` | string | 否 | Full EFO ontology URI (e.g., 'http://www.ebi.ac.uk/efo/EFO_0001645') |
| `efo_id` | string | 否 | EFO/OBA term ID (e.g., 'EFO_0001645', 'OBA_2050062'). Recommended for reliable trait filtering. |
| `efo_trait` | string | 否 | Exact EFO trait label. Use when you know the canonical trait string. |
| `size` | integer | 否 | Number of results to return per page |
| `limit` | integer | 否 | Alias for size. Number of results to return per page |
| `page` | integer | 否 | Page number for pagination |
| `sort` | string | 否 | Sort field (e.g., 'p_value', 'or_value'). Defaults to 'p_value' so the most significant associations appear first. |
| `direction` | string | 否 | Sort direction ('asc' or 'desc'). Defaults to 'asc'. |

??? note "上游说明（英文原文）"

    Search the GWAS Catalog for all genetic variants (SNPs) linked to a specific disease or trait. Accepts free-text disease/trait names (e.g., 'diabetes', 'breast cancer'), EFO ontology IDs, or EFO trait labels. Returns variant details including rsIDs, p-values, mapped genes, and genomic locations. Results are sorted by p-value ascending (most significant first) by default, so the first page holds the strongest loci; override with the 'sort' and 'direction' parameters. This is the best tool for finding all GWAS variants for a disease when you do NOT already have a specific variant ID. A free-text disease_trait is resolved to an EFO/MONDO term first, and that resolution is disclosed: 'query_trait', 'resolved_efo_id', 'resolved_efo_label' and 'trait_resolution_source' are returned, plus a 'trait_resolution_note' whenever the resolved label is not an exact match for your query - check them before trusting the results, and pass efo_id explicitly to bypass the resolver.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gwas_get_variants_for_trait', 'arguments': {}},
)
```

## `gwas_search_associations` {#gwas_search_associations}

Searches GWAS Catalog associations by trait, rsID, EFO ID, study accession or p-value threshold.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**gwas catalog, association search, efo, rsid, p-value
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `disease_trait` | string | 否 | Disease or trait name for text-based search (e.g., 'diabetes', 'coronary artery disease', 'breast cancer'). NOTE: GWAS associations are indexed by trait/disease, not by gene or drug name — passing a gene symbol or drug name here will fail. Use efo_id for precise filtering. |
| `query` | string | 否 | Alias for disease_trait. Disease or trait name (e.g., 'type 2 diabetes'). NOTE: gene symbols and drug names are not valid GWAS trait queries. |
| `efo_uri` | string | 否 | Full EFO ontology URI (e.g., 'http://www.ebi.ac.uk/efo/EFO_0001645') |
| `efo_id` | string | 否 | EFO/OBA term ID (e.g., 'EFO_0001645', 'OBA_2050062'). Recommended for reliable trait filtering. |
| `efo_trait` | string | 否 | Exact EFO trait label (e.g., 'serum alanine aminotransferase amount'). Use when you know the canonical trait string. |
| `rs_id` | string | 否 | dbSNP rs identifier |
| `accession_id` | string | 否 | Study accession identifier |
| `p_value` | number | 否 | Maximum p-value threshold for filtering (e.g., 5e-8 for genome-wide significance). IMPORTANT: the GWAS Catalog REST API v2 has NO server-side p-value filter (its OpenAPI document lists only pubmed_id, rs_id, full_pvalue_set, accession_id, efo_trait, efo_id, show_child_trait, mapped_gene, extended_geneset, sort, direction, page and size on GET /v2/associations), so this threshold is applied CLIENT… |
| `p_value_threshold` | number | 否 | Alias for p_value: maximum p-value threshold (also applied client-side; see p_value). |
| `sort` | string | 否 | Sort field (e.g., 'p_value', 'or_value'). Defaults to 'p_value'. When a p_value threshold is supplied and sort is 'p_value', the upstream page is always fetched ascending so the threshold sees the significant end; metadata.sort_override records this and the returned rows are ordered per your 'direction'. |
| `direction` | string | 否 | Sort direction ('asc' or 'desc'). Defaults to 'asc'. Combined with a p_value threshold this orders the RETURNED rows; see metadata.sort_override. |
| `size` | integer | 否 | Number of results to return |
| `page` | integer | 否 | Page number for pagination |

??? note "上游说明（英文原文）"

    Search the GWAS Catalog for genome-wide association study results by keyword. Accepts free-text disease/trait names (e.g., 'melanoma', 'diabetes', 'breast cancer'), SNP rs IDs, EFO ontology IDs, or study accession IDs. This is the primary tool for discovering GWAS associations when you do NOT already have a specific association ID. Returns SNPs, p-values, odds ratios, mapped genes, and study metadata with pagination support. Response envelope keys beyond 'data': 'metadata.pagination' (upstream page block - its totalElements/totalPages always describe the UNFILTERED result set); when a p_value threshold is given, 'metadata.p_value_filter', 'metadata.p_value_filter_scope' and 'metadata.p_value_filter_note' state that the threshold was applied CLIENT-SIDE (GWAS Catalog REST v2 has no server-side p-value filter), 'metadata.pagination_totals_are_prefilter' is true, 'metadata.filtered_total' / 'filtered_total_is_exact' / 'filtered_total_at_least' / 'filtered_total_note' give the post-filter denominator, and 'metadata.sort_override' records that the page was fetched p_value-ascending so the threshold could match anything. When a disease_trait is resolved to an ontology term, 'query_trait', 'resolved_efo_id', 'resolved_efo_label' and 'trait_resolution_source' are returned, plus a 'trait_resolution_note' whenever the resolved label is not an exact match for the trait you asked for. An empty 'data' comes with a 'note' that names the real cause: a client-side filter that removed every row is reported as such, and only a genuinely empty unfiltered result set suggests trying a trait synonym or a gene-based tool.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gwas_search_associations', 'arguments': {}},
)
```

## `GWAS_search_associations_by_gene` {#GWAS_search_associations_by_gene}

Searches GWAS Catalog associations by gene name returning strongest risk allele and p-values.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`gene_name`
- **关键词：**gwas catalog, risk allele, p-value, gene associations
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_name` | string | 是 | Gene symbol (e.g., BRCA1). |
| `size` | integer | 否 | Max associations to return (default 100; results are sorted by p-value, most significant first). 默认 `100` |

??? note "上游说明（英文原文）"

    Search GWAS Catalog associations by gene name (returns strongest risk allele and p-value fields).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GWAS_search_associations_by_gene',
     'arguments': {'gene_name': '...'}},
)
```

## `gwas_search_snps` {#gwas_search_snps}

Searches GWAS Catalog SNPs by rsID or mapped gene returning location, genes and functional class.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**rsid, snp search, mapped gene, functional class
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `rs_id` | string | 否 | dbSNP rs identifier |
| `mapped_gene` | string | 否 | Gene name or symbol |
| `size` | integer | 否 | Number of results to return |
| `page` | integer | 否 | Page number for pagination |

??? note "上游说明（英文原文）"

    Search the GWAS Catalog for single nucleotide polymorphisms (SNPs) by rs ID (e.g., 'rs7903146') or mapped gene name (e.g., 'TCF7L2'). Returns SNP details including chromosomal location, mapped genes, and functional class. Use this to look up specific variants or find all GWAS-catalogued SNPs for a gene.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gwas_search_snps', 'arguments': {}},
)
```

## `gwas_search_studies` {#gwas_search_studies}

Searches GWAS studies by trait, cohort, GxE and summary-statistics availability.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**gwas studies, cohort, gxe, summary statistics
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `disease_trait` | string | 否 | Disease or trait name for text-based search (e.g., 'diabetes', 'breast cancer') |
| `efo_uri` | string | 否 | Full EFO ontology URI (e.g., 'http://www.ebi.ac.uk/efo/EFO_0001645') |
| `efo_id` | string | 否 | EFO/OBA term ID (e.g., 'EFO_0001645', 'OBA_2050062'). Recommended for reliable trait filtering. |
| `efo_trait` | string | 否 | Exact EFO trait label. Use when you know the canonical trait string. |
| `cohort` | string | 否 | Cohort name (e.g., 'UKB' for UK Biobank) |
| `gxe` | boolean | 否 | Filter for Gene-by-Environment interaction studies |
| `full_pvalue_set` | boolean | 否 | Filter for studies with full summary statistics |
| `size` | integer | 否 | Number of results to return |
| `page` | integer | 否 | Page number for pagination |

??? note "上游说明（英文原文）"

    Search for GWAS studies by various criteria including disease trait, cohort, GxE interactions, and summary statistics availability.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gwas_search_studies', 'arguments': {}},
)
```

## `GWASSumStats_get_region_associations` {#GWASSumStats_get_region_associations}

Queries full GWAS summary statistics for variants in a chromosomal region with p-values, betas and allele frequencies.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`chromosome`, `bp_lower`, `bp_upper`
- **关键词：**summary statistics, region, beta, odds ratio, allele frequency
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chromosome` | integer | 是 | Chromosome number (1-22, or 23 for X). |
| `bp_lower` | integer | 是 | Start position of the region (GRCh38 coordinates). |
| `bp_upper` | integer | 是 | End position of the region (GRCh38 coordinates). |
| `p_upper` | number \| null | 否 | Maximum p-value threshold (default 5e-8 for genome-wide significance). Use larger values like 0.001 to include suggestive hits. |
| `study_accession` | string \| null | 否 | Filter by specific GWAS study accession (e.g., 'GCST002245'). Get study IDs from GWASSumStats_list_studies or GWASSumStats_get_trait_studies. |
| `size` | integer \| null | 否 | Maximum number of variants to return (default 50, max 1000). |

??? note "上游说明（英文原文）"

    Query full GWAS summary statistics for variants in a chromosomal region. Returns variant-level association data including p-values, effect sizes (beta/OR), allele frequencies, and study/trait information. Unlike gwas_search_associations (which returns curated top hits), this endpoint provides complete summary statistics from deposited studies. Filter by chromosome region, p-value threshold, and optionally by specific study. Results are sorted most-significant-first. Some deposited studies encode an unreported p-value as the GWAS Catalog missing-value sentinel -99; such rows cannot satisfy p_upper and are excluded when p_upper is set (the count is reported in metadata.sentinel_rows_excluded). Omit p_upper to retain them, flagged with p_value_reported = false.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GWASSumStats_get_region_associations',
     'arguments': {'chromosome': 0, 'bp_lower': 0, 'bp_upper': 0}},
)
```

## `GWASSumStats_get_trait_studies` {#GWASSumStats_get_trait_studies}

Gets GWAS studies with deposited full summary statistics for an EFO trait ID.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`trait_id`
- **关键词：**summary statistics, efo, gwas studies
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `trait_id` | string | 是 | EFO trait ontology ID (e.g., 'EFO_0000249' for Alzheimer disease, 'EFO_0001645' for coronary heart disease, 'EFO_0001359' for type 2 diabetes). |

??? note "上游说明（英文原文）"

    Get GWAS studies with deposited summary statistics for a specific EFO trait. Returns study accession IDs for studies that have full variant-level summary statistics. Use EFO trait IDs like 'EFO_0000249' (Alzheimer disease), 'EFO_0001645' (coronary heart disease), 'EFO_0000253' (asthma), 'EFO_0001359' (type 2 diabetes). Get EFO IDs from ols_search_efo_terms or gwas_search_studies.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GWASSumStats_get_trait_studies',
     'arguments': {'trait_id': '...'}},
)
```

## `GWASSumStats_list_studies` {#GWASSumStats_list_studies}

Lists GWAS Catalog studies with complete deposited variant-level summary statistics.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**summary statistics, gwas catalog, studies
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `size` | integer \| null | 否 | Number of studies to return (default 20, max 100). |

??? note "上游说明（英文原文）"

    List GWAS studies that have deposited full summary statistics with the EBI GWAS Catalog. Unlike the main GWAS Catalog (which stores curated top hits), studies listed here have complete variant-level summary statistics available for download and query. Returns study accession IDs that can be used as filters in GWASSumStats_get_region_associations.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GWASSumStats_list_studies', 'arguments': {}},
)
```

## `HPO_get_disease_annotations` {#HPO_get_disease_annotations}

Gets a disease's HPO annotation profile: MAxO medical actions, phenotypes by body system and causative genes.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`disease_id`
- **关键词：**hpo, disease annotation, causative genes, maxo, omim
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `disease_id` | string | 是 | Disease identifier as an OMIM, ORPHA, or DECIPHER CURIE (e.g., 'OMIM:154700' for Marfan syndrome, 'OMIM:277900' for Wilson disease, 'ORPHA:558'). |

??? note "上游说明（英文原文）"

    Get the full curated HPO annotation profile for a disease, keyed by a disease ID (OMIM / ORPHA / DECIPHER), from the JAX HPO network annotations. Returns three things: (1) medical_actions — MAxO curated treatments / management actions, each with TREATS / PREVENTS relations and the specific phenotypes they target; (2) categories — the disease's phenotypes grouped by body system (e.g., Cardiovascular, Nervous System, Eye), each with frequency / onset; (3) genes — the causative genes (NCBI gene IDs and symbols). Use to retrieve evidence-based management options and the organ-system phenotype breakdown for a known disease. Example: OMIM:154700 (Marfan syndrome) returns 3 medical actions including MAXO:0000653 'angiotensin receptor blocker therapy' (PREVENTS/TREATS aortic root aneurysm), 13 body systems, and gene NCBIGene:2200 (FBN1). Other examples: OMIM:277900 (Wilson disease) → copper chelator / zinc therapy; OMIM:261600 (phenylketonuria) → dietary phenylalanine avoidance.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'HPO_get_disease_annotations',
     'arguments': {'disease_id': '...'}},
)
```

## `HPO_get_genes_by_phenotype` {#HPO_get_genes_by_phenotype}

Gets genes in which variants cause an HPO phenotype term, with NCBI gene IDs and symbols.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`term_id`
- **关键词：**hpo, phenotype genes, rare disease, ncbi gene
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `term_id` | string | 是 | HPO term identifier (e.g., 'HP:0001250'). Must start with 'HP:' followed by 7 digits. |
| `limit` | integer \| null | 否 | Maximum number of results to return (1-500, default 50). Phenotype terms can be linked to thousands of genes. |
| `offset` | integer \| null | 否 | Number of results to skip before returning 'limit' results (default 0). Use with the response's metadata.has_more/metadata.total to page through phenotypes with more genes than fit in one call. |

??? note "上游说明（英文原文）"

    Get the genes associated with an HPO phenotype term (genes in which variants cause this phenotype), from the JAX HPO network annotations. Returns NCBI gene IDs and symbols. Example: HP:0001250 (Seizure) is linked to genes such as RELN, SCN1A, and ~2000 others. Use to go from an observed patient phenotype to candidate disease genes for variant prioritization. A phenotype term can be linked to thousands of genes but 'limit' caps a single response at 500 -- check metadata.has_more in the response and pass 'offset' to page through the rest; a well-known gene for a common phenotype may not appear in the first page.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'HPO_get_genes_by_phenotype',
     'arguments': {'term_id': '...'}},
)
```

## `IMPC_get_gene_phenotype_hits` {#IMPC_get_gene_phenotype_hits}

Gets IMPC phenotyping statistical results (p-values, effect sizes) for a mouse gene knockout.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**impc, knockout mouse, phenotyping, effect size
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 否 | Mouse gene symbol (e.g., 'Trp53', 'Brca1') |
| `mgi_id` | string | 否 | MGI accession ID (e.g., 'MGI:98834') |
| `significant_only` | boolean | 否 | Return only significant results (default: true) 默认 `True` |
| `limit` | integer | 否 | Maximum results (default: 100) 默认 `100` |

??? note "上游说明（英文原文）"

    Get detailed statistical results from IMPC phenotyping for a gene. Returns p-values, effect sizes, and classification tags for each phenotyping parameter. Filter by significant results only. Useful for quantifying phenotype severity.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'IMPC_get_gene_phenotype_hits', 'arguments': {}},
)
```

## `IMPC_get_gene_summary` {#IMPC_get_gene_summary}

Gets IMPC mouse gene summary: phenotyping status, viability, orthologs and MP/HP annotations.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**impc, mgi, knockout, viability, mammalian phenotype
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 否 | Mouse gene symbol (e.g., 'Trp53', 'Brca1', 'Wdr7'). Case-sensitive. Use mouse symbol, not human. |
| `mgi_id` | string | 否 | MGI accession ID (e.g., 'MGI:98834'). Use if gene symbol is ambiguous. |

??? note "上游说明（英文原文）"

    Get mouse gene summary from IMPC including phenotyping status, viability, ortholog info, and MP/HP term annotations. Query by human gene symbol or MGI ID. No API key needed. Essential for safety/liability triangulation via mouse KO phenotypes.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'IMPC_get_gene_summary', 'arguments': {}},
)
```

## `IMPC_get_phenotypes_by_gene` {#IMPC_get_phenotypes_by_gene}

Gets MP phenotype annotations from IMPC knockout screens for a gene, grouped by system.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**impc, knockout, mp terms, mouse phenotype
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 否 | Mouse gene symbol (e.g., 'Trp53', 'Brca1') |
| `mgi_id` | string | 否 | MGI accession ID (e.g., 'MGI:98834') |
| `limit` | integer | 否 | Maximum phenotype records to return (default: 100) 默认 `100` |

??? note "上游说明（英文原文）"

    Get all phenotype annotations for a mouse gene from IMPC knockout screens. Returns MP ontology terms grouped by physiological system (e.g., nervous, cardiovascular, immune). Essential for understanding gene function and predicting drug side effects.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'IMPC_get_phenotypes_by_gene', 'arguments': {}},
)
```

## `IntOGen_get_drivers` {#IntOGen_get_drivers}

Returns cancer driver genes for a cancer type identified by IntOGen methods on somatic mutation data.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`cancer_type`
- **关键词：**cancer driver genes, somatic mutation, dndscv, tumor type
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `cancer_type` | string | 是 | IntOGen cancer type code. Examples: 'BRCA' (breast), 'LUAD' (lung adenocarcinoma), 'COADREAD' (colorectal), 'AML' (acute myeloid leukemia), 'SKCM' (melanoma), 'GBM' (glioblastoma). Use IntOGen_list_cancer_types to see all available codes. |

??? note "上游说明（英文原文）"

    Get cancer driver genes for a specific cancer type from IntOGen. Returns genes identified as drivers by computational methods (dNdScv, OncodriveFML, OncodriveCLUSTL, HotMAPS, smRegions, CBaSE, MutPanning) applied to somatic mutation data from thousands of tumors. Each result includes the gene symbol, number of mutations and mutated samples across cohorts. Use IntOGen_list_cancer_types to see available cancer type codes. Example: cancer_type='BRCA' returns 118 driver genes for breast cancer including TP53, PIK3CA, GATA3.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'IntOGen_get_drivers',
     'arguments': {'cancer_type': '...'}},
)
```

## `IntOGen_get_gene_info` {#IntOGen_get_gene_info}

Returns cancer types where a gene is an IntOGen cancer driver and the detection methods supporting it.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`gene`
- **关键词：**driver gene, cancer type, somatic mutation, oncodrive
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 是 | Gene symbol (HUGO nomenclature). Examples: 'TP53', 'KRAS', 'PIK3CA', 'BRAF', 'EGFR', 'PTEN', 'APC', 'BRCA1', 'IDH1', 'NRAS'. |

??? note "上游说明（英文原文）"

    Get driver gene classification for a specific gene across all cancer types from IntOGen. Returns the list of cancer types where the gene is identified as a cancer driver, along with the detection methods used (dNdScv, OncodriveFML, OncodriveCLUSTL, HotMAPS, smRegions, CBaSE, MutPanning) and mutation statistics. Example: gene='TP53' is a driver in 74 cancer types; gene='KRAS' is a driver in 22 cancer types.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'IntOGen_get_gene_info', 'arguments': {'gene': '...'}},
)
```

## `KEGG_get_disease` {#KEGG_get_disease}

Returns a KEGG disease entry with category, associated genes and mutations, pathways, drugs and cross-references.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`disease_id`
- **关键词：**kegg disease, disease genes, mutations, omim
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `disease_id` | string | 是 | KEGG disease ID in H\#\#\#\#\# format (e.g., 'H00001' for B-cell ALL, 'H00408' for type 1 diabetes, 'H00409' for type 2 diabetes, 'H00031' for breast cancer). |

??? note "上游说明（英文原文）"

    Get detailed KEGG disease information including disease name, category, description, associated genes (with mutations/translocations), pathways, therapeutic drugs, database cross-references (ICD-10, MeSH, OMIM), and references. Example: H00001 returns B-cell acute lymphoblastic leukemia with BCR-ABL, MLL-AF4, E2A-PBX1 gene fusions and 10+ associated pathways. Use KEGG_search_disease to find disease IDs.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'KEGG_get_disease', 'arguments': {'disease_id': '...'}},
)
```

## `KEGG_get_disease_genes` {#KEGG_get_disease_genes}

Returns human KEGG gene IDs linked to a KEGG disease through mutations, fusions or other mechanisms.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`disease_id`
- **关键词：**disease genes, gene fusion, mutation, kegg disease
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `disease_id` | string | 是 | KEGG disease ID (e.g., 'H00001' for B-cell ALL, 'H00031' for breast cancer). |
| `organism` | string | 否 | KEGG organism code (default: 'hsa' for human). Other examples: 'mmu' (mouse). 默认 `hsa` |

??? note "上游说明（英文原文）"

    Get all human genes linked to a KEGG disease. Returns KEGG gene IDs (hsa:XXXXX format) for genes associated with the disease through mutations, gene fusions, or other mechanisms. Example: H00001 (B-cell ALL) returns hsa:25 (ABL1), hsa:4297 (KMT2A), hsa:4299 (AFF1), etc. Useful for building disease-gene association networks. Optionally specify a different organism.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'KEGG_get_disease_genes',
     'arguments': {'disease_id': '...'}},
)
```

## `LDlink_get_proxies` {#LDlink_get_proxies}

Returns variants in linkage disequilibrium with a query SNP for a population via LDlink LDproxy.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`variant`
- **关键词：**linkage disequilibrium, ld proxy, r2, rsid, population
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant` | string | 是 | Query SNP rsID, e.g. 'rs7903146'. |
| `population` | string \| null | 否 | 1000 Genomes population code (default 'CEU'). Examples: 'CEU' (European), 'YRI' (African), 'CHB' (East Asian), or super-populations 'EUR'/'AFR'/'EAS'/'SAS'/'AMR'. |
| `r2_threshold` | number \| null | 否 | Minimum R2 to report (default 0.8; high LD). |
| `genome_build` | string \| null | 否 | 'grch38' (default) or 'grch37'. |
| `limit` | integer \| null | 否 | Max proxies to return (1-500, default 50). |

??? note "上游说明（英文原文）"

    Find variants in linkage disequilibrium (LD proxies) with a query SNP, population-specific, via NIH LDlink's LDproxy module. Essential for GWAS interpretation: a lead SNP is often just a tag for the true causal variant, so high-R2 proxies (especially in coding/regulatory regions) are better mechanistic candidates. Returns proxy rsIDs with R2, coordinates, alleles, MAF, and distance. Requires a free LDlink token (LDLINK_TOKEN env var; register at ldlink.nih.gov).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'LDlink_get_proxies', 'arguments': {'variant': '...'}},
)
```

## `MARRVEL_get_omim_phenotypes` {#MARRVEL_get_omim_phenotypes}

Returns OMIM phenotype associations for a human gene via MARRVEL with MIM numbers and inheritance.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**omim, mendelian disease, inheritance, phenotype
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `symbol` | string \| null | 否 | HGNC gene symbol, e.g. 'CFTR', 'BRCA1'. |
| `gene_symbol` | string \| null | 否 | Alias for symbol. |
| `gene` | string \| null | 否 | Alias for symbol. |

??? note "上游说明（英文原文）"

    Get OMIM phenotype/disease associations for a human gene by symbol via MARRVEL: each association's phenotype name, gene and phenotype MIM numbers, mode of inheritance, and phenotypic series. Use to enumerate the Mendelian diseases linked to a gene and their inheritance patterns. MARRVEL indexes only current approved HGNC symbols, so a retired symbol or alias (e.g. 'RP20'/'LCA2' for RPE65) is resolved once via HGNC and re-queried; when that happens the substitution is disclosed in metadata.resolved_symbol and metadata.resolved_note. No API key required.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MARRVEL_get_omim_phenotypes', 'arguments': {}},
)
```

## `MedGen_get_condition` {#MedGen_get_condition}

Returns a MedGen condition by UID or CUI with definition, associated genes, OMIM IDs, inheritance and HPO features.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**medgen, genetic condition, omim, inheritance
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `uid` | string | 否 | MedGen UID (e.g., '41393' for cystic fibrosis). Get from MedGen_search_conditions. |
| `cui` | string | 否 | UMLS Concept Unique Identifier (e.g., 'C0010674' for cystic fibrosis). Alternative to uid. |
| `concept_id` | string \| null | 否 | Alias for cui. UMLS CUI identifier (e.g., "C0017205" for Gaucher disease). |

??? note "上游说明（英文原文）"

    Get detailed information about a specific genetic condition from NCBI MedGen by UID or UMLS CUI. Returns full definition, associated genes, OMIM IDs, synonyms, modes of inheritance, and clinical features (HPO phenotypes). Example: uid='41393' returns cystic fibrosis with CFTR/TGFB1/FCGR2A genes, autosomal recessive inheritance, 30+ clinical features.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MedGen_get_condition', 'arguments': {}},
)
```

## `MedGen_search_conditions` {#MedGen_search_conditions}

Searches MedGen genetic conditions by keyword and returns names, definitions, OMIM IDs and associated genes.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`query`
- **关键词：**medgen search, genetic condition, omim, disease genes
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Search term (e.g., 'cystic fibrosis', 'BRCA1', 'Marfan syndrome', 'autosomal recessive'). |
| `max_results` | integer | 否 | Maximum number of results (default 10, max 50). 默认 `10` |

??? note "上游说明（英文原文）"

    Search NCBI MedGen for genetic conditions, diseases, or syndromes by keyword. Returns condition names, definitions, OMIM IDs, and associated genes. MedGen aggregates data from OMIM, Orphanet, ClinVar, HPO, GTR, and GeneReviews. Example: query='cystic fibrosis' returns CF entry with CFTR gene, OMIM 219700.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MedGen_search_conditions', 'arguments': {'query': '...'}},
)
```

## `MedlinePlus_get_genetics_condition_by_name` {#MedlinePlus_get_genetics_condition_by_name}

Returns MedlinePlus Genetics information for a condition including related genes and inheritance pattern.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`condition`
- **关键词：**genetic condition, inheritance, related genes
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `condition` | string | 是 | URL slug of genetic condition, e.g., "alzheimer-disease", must match MedlinePlus page path. |
| `format` | string | 否 | Format parameter for the endpoint (note: API always returns XML regardless of this parameter). 可选值：`json`, `xml` 默认 `json` |

??? note "上游说明（英文原文）"

    Get detailed information from MedlinePlus Genetics corresponding to genetic condition name. Returns condition description, related genes, inheritance pattern, synonyms, and other details.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MedlinePlus_get_genetics_condition_by_name',
     'arguments': {'condition': '...'}},
)
```

## `MGI_get_phenotypes` {#MGI_get_phenotypes}

Returns Mammalian Phenotype annotations for a mouse gene from mutant allele studies via Alliance/MGI.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`gene_id`
- **关键词：**mammalian phenotype, mouse mutant, mp ontology, gene-phenotype
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_id` | string | 是 | MGI gene identifier with 'MGI:' prefix. Examples: 'MGI:98834' (Trp53), 'MGI:104537' (Brca1), 'MGI:96677' (Kit). |
| `limit` | integer | 否 | Number of phenotype annotations to return per page (default: 20, max: 100). 默认 `20` |
| `page` | integer | 否 | Page number for pagination (default: 1). 默认 `1` |

??? note "上游说明（英文原文）"

    Get phenotype annotations for a mouse gene from MGI (Mouse Genome Informatics). Returns Mammalian Phenotype (MP) ontology annotations including phenotype statements derived from mutant allele studies. Uses Alliance of Genome Resources API which integrates MGI phenotype curation data. Useful for understanding mouse gene function through loss-of-function and gain-of-function phenotypes. Example: 'MGI:98834' returns Trp53 phenotypes including tumor susceptibility, apoptosis defects, etc.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MGI_get_phenotypes', 'arguments': {'gene_id': '...'}},
)
```

## `Monarch_get_associations` {#Monarch_get_associations}

Queries Monarch gene-phenotype, disease-phenotype, gene-disease or gene-pathway associations by CURIE.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**biolink, gene-disease, gene-phenotype, knowledge graph
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `subject` | string | 否 | Subject entity CURIE. Examples: 'HGNC:11998' (TP53), 'MONDO:0005148' (type 2 diabetes). |
| `object` | string | 否 | Object entity CURIE. Use instead of subject to find entities associated TO this object. |
| `category` | string | 否 | Optional Biolink association category filter. Options: 'biolink:GeneToPhenotypicFeatureAssociation', 'biolink:DiseaseToPhenotypicFeatureAssociation', 'biolink:CorrelatedGeneToDiseaseAssociation', 'biolink:CausalGeneToDiseaseAssociation', 'biolink:GeneToPathwayAssociation'. |
| `limit` | integer | 否 | Maximum associations to return (default: 20, max: 200). |

??? note "上游说明（英文原文）"

    Query disease-gene-phenotype associations in the Monarch Initiative knowledge graph. Find gene-phenotype, disease-phenotype, gene-disease, or gene-pathway associations. Provide a subject CURIE (entity of interest) or object CURIE, with optional category filter. Uses the /association/all endpoint. Example: subject=HGNC:11998 returns all TP53 associations including phenotypes and diseases.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Monarch_get_associations', 'arguments': {}},
)
```

## `Monarch_get_gene_diseases` {#Monarch_get_gene_diseases}

Returns gene-disease associations with evidence for a gene CURIE from Monarch.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`subject`
- **关键词：**gene-disease, monarch, hgnc curie
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `subject` | string | 是 | Gene CURIE identifier (e.g., 'HGNC:11998' for TP53). Use Monarch_search_gene to find this. |
| `limit` | integer | 否 | Maximum disease associations to return (default: 100) 默认 `100` |
| `offset` | integer | 否 | Number of initial entries to skip. 默认 `0` |

??? note "上游说明（英文原文）"

    Get diseases associated with a gene from Monarch Initiative. Input is a gene CURIE (e.g., 'HGNC:11998' for TP53). Returns disease associations with evidence. Bridges gene→disease across species for hypothesis generation.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Monarch_get_gene_diseases',
     'arguments': {'subject': '...'}},
)
```

## `Monarch_get_gene_phenotypes` {#Monarch_get_gene_phenotypes}

Returns human and model organism phenotypes associated with a gene CURIE from Monarch.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`subject`
- **关键词：**gene-phenotype, hpo, mp, cross-species
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `subject` | string | 是 | Gene CURIE identifier (e.g., 'HGNC:11998' for TP53, 'HGNC:1100' for BRCA1). Use Monarch_search_gene to find this. |
| `limit` | integer | 否 | Maximum phenotype associations to return (default: 100) 默认 `100` |
| `offset` | integer | 否 | Number of initial entries to skip. 默认 `0` |

??? note "上游说明（英文原文）"

    Get phenotypes associated with a gene from Monarch Initiative (cross-species). Input is a gene CURIE (e.g., 'HGNC:11998' for TP53). Returns human and model organism phenotypes (HP/MP terms). Use Monarch_search_gene first to find the CURIE.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Monarch_get_gene_phenotypes',
     'arguments': {'subject': '...'}},
)
```

## `MonarchV3_get_associations` {#MonarchV3_get_associations}

Queries cross-species Monarch associations by subject/object CURIE and Biolink category.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`category`
- **关键词：**biolink, gene-disease, cross-species, associations
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `subject` | string | 否 | Subject-side entity CURIE, i.e. the entity on the LEFT of the association (the arrow's source). For gene-&gt;disease categories (CausalGeneToDiseaseAssociation, CorrelatedGeneToDiseaseAssociation) the subject is the GENE; for gene-&gt;phenotype and gene-&gt;pathway the subject is the gene; for DiseaseToPhenotypicFeatureAssociation the subject is the disease. Examples: 'HGNC:11998' (TP53), 'MONDO:0005148'… |
| `object` | string | 否 | Object-side entity CURIE, i.e. the entity on the RIGHT of the association (the arrow's target). Use this to query the association backwards. For gene-&gt;disease categories the object is the DISEASE, so 'which genes cause this disease?' is object='MONDO:0006559' with category='biolink:CausalGeneToDiseaseAssociation'; for gene/disease-&gt;phenotype categories the object is the HP phenotype term. Optiona… |
| `category` | string | 是 | Biolink association category, written as subject-&gt;object. Options: 'biolink:GeneToPhenotypicFeatureAssociation' (gene-&gt;phenotypes), 'biolink:DiseaseToPhenotypicFeatureAssociation' (disease-&gt;phenotypes), 'biolink:CorrelatedGeneToDiseaseAssociation' (gene-&gt;diseases, correlated), 'biolink:CausalGeneToDiseaseAssociation' (gene-&gt;diseases, causal), 'biolink:GeneToPathwayAssociation' (gene-&gt;pathways), '… |
| `limit` | integer | 否 | Maximum number of associations to return (default: 20, max: 200). |

??? note "上游说明（英文原文）"

    Query cross-species associations between biomedical entities in the Monarch Initiative knowledge graph. Find gene-phenotype, disease-phenotype, gene-disease, or gene-pathway associations. Uses Biolink model categories to specify association type. Integrates data from OMIM, ClinVar, HPO, MGI, ZFIN, FlyBase, WormBase. Associations are DIRECTED: every category has a fixed subject side and object side, so anchor your entity on the correct side. Supply 'subject', 'object', or both (at least one is required). Example: HGNC:11998 (TP53) as subject with category GeneToPhenotypicFeatureAssociation returns phenotypes like hepatocellular carcinoma; MONDO:0006559 as OBJECT with category CausalGeneToDiseaseAssociation returns the genes that cause that disease (NCSTN, PSENEN, PSEN1).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MonarchV3_get_associations',
     'arguments': {'category': '...'}},
)
```

## `MR_estimate` {#MR_estimate}

Computes two-sample Mendelian randomization estimates (IVW, MR-Egger, weighted median) from SNP exposure/outcome betas.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`beta_exposure`, `se_exposure`, `beta_outcome`, `se_outcome`
- **关键词：**mendelian randomization, ivw, mr-egger, instruments, causal inference
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `beta_exposure` | array of number | 是 | Per-instrument SNP-exposure effect sizes (harmonized to the same effect allele; non-zero). |
| `se_exposure` | array of number | 是 | Standard errors of beta_exposure (positive). |
| `beta_outcome` | array of number | 是 | Per-instrument SNP-outcome effect sizes (same SNP order / effect allele). |
| `se_outcome` | array of number | 是 | Standard errors of beta_outcome (positive). |

??? note "上游说明（英文原文）"

    Two-sample Mendelian randomization causal estimate from harmonized instrument summary statistics (SNP-exposure + SNP-outcome betas and SEs, e.g. from OpenGWAS_get_mr_instruments). Returns IVW (primary), MR-Egger (with directional-pleiotropy intercept), weighted-median estimates, and Cochran's Q / I^2 heterogeneity. Pure-compute (no R, no API key); reimplements TwoSampleMR's core estimators.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MR_estimate',
     'arguments': {'beta_exposure': [],
                   'se_exposure': [],
                   'beta_outcome': [],
                   'se_outcome': []}},
)
```

## `MyDisease_get_disease` {#MyDisease_get_disease}

Returns MyDisease.info disease annotations from MONDO, DO, CTD, HPO and DisGeNET.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`disease_id`
- **关键词：**mydisease, disgenet, disease annotation
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `disease_id` | string | 是 | Disease identifier. Supports MONDO IDs (e.g., 'MONDO:0005148' for T2DM), DOID IDs (e.g., 'DOID:9352'), OMIM IDs (e.g., 'OMIM:125853'), or MESH IDs (e.g., 'MESH:D003924'). MONDO IDs provide the most comprehensive results. |
| `fields` | string | 否 | Comma-separated list of data fields to return. Options: 'mondo', 'disease_ontology', 'ctd', 'hpo', 'disgenet', 'all'. Default: 'mondo,disease_ontology,ctd,hpo'. Use 'all' for complete data. |

??? note "上游说明（英文原文）"

    Get comprehensive disease annotations from MyDisease.info, a BioThings disease aggregator that integrates data from MONDO, Disease Ontology, CTD (chemicals and pathways related to disease), HPO (phenotypes), and DisGeNET (gene-disease associations). Returns standardized disease identifiers, definitions, cross-references, associated chemicals, pathways, and phenotypes. Example: MONDO:0005148 (type 2 diabetes) returns 160 CTD chemicals, 984 related pathways, and 3 HPO phenotypes including late onset and autosomal dominant inheritance.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MyDisease_get_disease',
     'arguments': {'disease_id': '...'}},
)
```

## `OMIM_get_entry` {#OMIM_get_entry}

Returns an OMIM entry by MIM number with description, clinical features, inheritance and molecular genetics.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`mim_number`
- **关键词：**omim entry, mim number, mendelian, molecular genetics
- **需要 API 密钥：**`OMIM_API_KEY`
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type (fixed: get_entry) 默认 `get_entry` |
| `mim_number` | string | 是 | OMIM MIM number (e.g., 164730, 219700). Can include OMIM: prefix. |
| `include` | string | 否 | Data sections to include (default: text,clinicalSynopsis,geneMap) 默认 `text,clinicalSynopsis,geneMap` |

??? note "上游说明（英文原文）"

    Get detailed OMIM entry by MIM number. Returns comprehensive information including text description, clinical features, inheritance, molecular genetics, and references. MIM numbers are 6-digit identifiers (e.g., 164730 for BRAF gene). Requires OMIM_API_KEY.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OMIM_get_entry', 'arguments': {'mim_number': '...'}},
)
```

## `OMIM_get_gene_map` {#OMIM_get_gene_map}

Returns OMIM gene map entries with chromosomal location, phenotypes and inheritance by MIM number or chromosome.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**omim gene map, chromosome, disease genes
- **需要 API 密钥：**`OMIM_API_KEY`
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type (fixed: get_gene_map) 默认 `get_gene_map` |
| `mim_number` | string | 否 | OMIM MIM number for specific gene (optional if chromosome provided) |
| `chromosome` | string | 否 | Chromosome number (1-22, X, Y) to get all genes (optional if mim_number provided) |
| `limit` | integer | 否 | Maximum results for chromosome query (default: 50, max: 100) 默认 `50` |

??? note "上游说明（英文原文）"

    Get OMIM gene-disease mapping information. Returns genes mapped to chromosomal locations with associated phenotypes and inheritance patterns. Query by MIM number or chromosome. Useful for finding all disease genes on a chromosome.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OMIM_get_gene_map', 'arguments': {}},
)
```

## `OMIM_search` {#OMIM_search}

Searches OMIM for Mendelian diseases and genes by keyword.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`query`
- **关键词：**omim search, mendelian disease, mim
- **需要 API 密钥：**`OMIM_API_KEY`
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type (fixed: search) 默认 `search` |
| `query` | string | 是 | Search query - gene name (BRCA1), disease (cystic fibrosis), or phenotype (intellectual disability) |
| `limit` | integer | 否 | Maximum results (default: 10, max: 20) 默认 `10` |
| `start` | integer | 否 | Offset for pagination (default: 0) 默认 `0` |

??? note "上游说明（英文原文）"

    Search OMIM database for Mendelian diseases and genes. OMIM is the authoritative source for human genetic disorders. Search by gene name, disease name, or phenotype keywords. Returns MIM numbers, titles, and entry types. Requires OMIM_API_KEY env var (register at omim.org/api).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OMIM_search', 'arguments': {'query': '...'}},
)
```

## `OncoKB_get_cancer_genes` {#OncoKB_get_cancer_genes}

Returns all OncoKB cancer genes classified as oncogenes or tumor suppressors.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**cancer genes, oncogene, tumor suppressor
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type (fixed: get_cancer_genes) 默认 `get_cancer_genes` |

??? note "上游说明（英文原文）"

    Get list of all cancer genes curated in OncoKB. Returns genes classified as oncogenes and/or tumor suppressors. Useful for checking if a gene is a known cancer driver. Requires ONCOKB_API_TOKEN env var.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OncoKB_get_cancer_genes', 'arguments': {}},
)
```

## `OncoKB_get_gene_info` {#OncoKB_get_gene_info}

Returns OncoKB gene-level information on oncogene or tumor suppressor status with background.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**oncogene, tumor suppressor, gene background
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type (fixed: get_gene_info) 默认 `get_gene_info` |
| `gene` | string | 否 | Gene symbol (e.g., BRAF, TP53, ROS1). Without ONCOKB_API_TOKEN only BRAF, TP53, and ROS1 are available (demo mode). |
| `gene_symbol` | string | 否 | Gene symbol alias — alternative to gene parameter |

??? note "上游说明（英文原文）"

    Get gene-level oncogenic information from OncoKB. Returns whether gene is an oncogene, tumor suppressor, or both. Includes curated background information. Requires ONCOKB_API_TOKEN env var (demo mode available for BRAF, TP53, ROS1).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OncoKB_get_gene_info', 'arguments': {}},
)
```

## `OpenGWAS_get_mr_instruments` {#OpenGWAS_get_mr_instruments}

Assembles harmonized LD-clumped SNP instruments for an exposure/outcome GWAS pair from IEU OpenGWAS.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`exposure_id`
- **关键词：**opengwas, mendelian randomization, snp instruments, ld clumping
- **需要 API 密钥：**`OPENGWAS_JWT`
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `exposure_id` | string | 是 | IEU OpenGWAS study ID for the exposure (e.g. 'ieu-a-2'). Find IDs via EpiGraphDB_search_opengwas. |
| `outcome_id` | string | 否 | IEU OpenGWAS study ID for the outcome (e.g. 'ieu-a-7'). Optional — omit to return only the exposure instruments. |
| `pval` | number | 否 | Instrument p-value threshold for the exposure (default 5e-8, genome-wide significance). 默认 `5e-08` |
| `clump` | integer | 否 | 1 to LD-clump instruments to independent SNPs (default), 0 to skip clumping. 默认 `1` |
| `r2` | number | 否 | LD r2 threshold for clumping (default 0.001). 默认 `0.001` |
| `kb` | integer | 否 | Clumping window in kb (default 10000). 默认 `10000` |
| `pop` | string | 否 | Reference population for LD clumping (default 'EUR'). 默认 `EUR` |

??? note "上游说明（英文原文）"

    Assemble harmonized two-sample Mendelian randomization (MR) instruments from the IEU OpenGWAS database for an arbitrary exposure/outcome GWAS pair. Fetches the genome-wide significant, LD-clumped SNP instruments for the exposure (/tophits) and those SNPs' effects in the outcome GWAS (/associations), then aligns both onto a common effect allele so the result is ready to feed into an IVW or MR-Egger estimator. Use this when EpiGraphDB's pre-computed MR-EvE does not cover the trait pair, or when you need custom instruments (own p-value/clumping thresholds). Requires a free OpenGWAS JWT token. Pass OpenGWAS study IDs (e.g. 'ieu-a-2' for BMI, 'ieu-a-7' for coronary heart disease); use EpiGraphDB_search_opengwas to look IDs up by trait name.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenGWAS_get_mr_instruments',
     'arguments': {'exposure_id': '...'}},
)
```

## `OpenTargets_get_associated_targets_by_disease_efoId` {#OpenTargets_get_associated_targets_by_disease_efoId}

Lists targets (genes) associated with a disease or phenotype by EFO ID, paginated with association scores.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`efoId`
- **关键词：**target-disease, association score, efo, gene prioritization
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `efoId` | string | 是 | The efoId of a disease or phenotype. |
| `size` | integer | 否 | Number of top-scored associated targets to return (default 50). A disease can have thousands of associated targets; raise this to retrieve more than the first page. 默认 `50` |
| `index` | integer | 否 | Zero-based page index for paginating beyond the first `size` targets (default 0). |

??? note "上游说明（英文原文）"

    Find targets associated with a specific disease or phenotype based on efoId. Results are PAGINATED (50 per call by default): `associatedTargets.count` is the total number of associated targets, `associatedTargets.returned` is how many are in `rows`, and `associatedTargets.truncated` is true when targets remain unfetched, in which case `metadata.truncation_note` names the `size`/`index` values to pass. Raise `size` (Open Targets rejects a page above 3000) and step `index` to page through; a common disease has thousands of associated targets, so the first page is NOT the full list.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_associated_targets_by_disease_efoId',
     'arguments': {'efoId': '...'}},
)
```

## `OpenTargets_get_credible_set_colocalisation` {#OpenTargets_get_credible_set_colocalisation}

Retrieves colocalisation results (COLOC, eCAVIAR, H4) for a GWAS/QTL credible set by studyLocusId.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`studyLocusId`
- **关键词：**coloc, h4, credible set, qtl
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `studyLocusId` | string | 是 | The credible-set / study-locus ID (e.g., '0004aaf4841475f0caf5bcde04432bcf'). |
| `size` | integer | 否 | Maximum number of colocalisation rows to return (1-500). 默认 `20` |

??? note "上游说明（英文原文）"

    Retrieve GWAS/QTL colocalisation evidence for a credible set (by studyLocusId) from OpenTargets. Returns the colocalisation count and rows, each with the method (COLOC, COLOC_PIP_ECAVIAR, etc.), the posterior probability h4 of a shared causal variant, the colocalisation posterior probability clpp, numberColocalisingVariants, and the other colocalising study locus (studyLocusId, study type/id). Use to link a GWAS signal to a molecular QTL.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_credible_set_colocalisation',
     'arguments': {'studyLocusId': '...'}},
)
```

## `OpenTargets_get_credible_set_detail` {#OpenTargets_get_credible_set_detail}

Returns detail for a GWAS credible set: lead variant, study, fine-mapping method, significance and L2G predictions.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`studyLocusId`
- **关键词：**fine-mapping, l2g, lead variant, credible set
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `studyLocusId` | string | 是 | Credible set identifier (studyLocusId), a 32-character hash (e.g., 'b758d8fb10924f5338cbad8d27c7dee8'). Obtain from variant credible sets or study credible sets queries. |

??? note "上游说明（英文原文）"

    Get detailed information about a specific GWAS credible set (fine-mapped locus) from Open Targets. Returns the lead variant, associated study, fine-mapping method, statistical significance, locus-to-gene (L2G) predictions mapping the signal to likely causal genes, and colocalization evidence with QTLs.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_credible_set_detail',
     'arguments': {'studyLocusId': '...'}},
)
```

## `OpenTargets_get_diseases_phenotypes_by_target_ensembl` {#OpenTargets_get_diseases_phenotypes_by_target_ensembl}

Lists diseases and phenotypes associated with a target Ensembl gene ID, paginated with association scores.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`ensemblId`
- **关键词：**gene-disease, association score, ensembl, phenotype
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ensemblId` | string | 是 | The ensemblId of a target. |
| `page` | object | 否 | Pagination parameters. `count` in the response is the total across all pages; omit this to get the first page only. |

??? note "上游说明（英文原文）"

    Find diseases or phenotypes associated with a specific target using ensemblId. `associatedDiseases.count` is the total number of associations; `rows` returns only one page of them (25 by default) -- pass `page` to retrieve more (e.g. page 2 of a target with thousands of associations like SIRT1's 3432). `associatedDiseases.returned` reports how many rows this response actually contains and `associatedDiseases.truncated` is true when associations remain unfetched, in which case `metadata.truncation_note` names the exact `page` argument to pass.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_diseases_phenotypes_by_target_ensembl',
     'arguments': {'ensemblId': '...'}},
)
```

## `OpenTargets_get_evidence_by_datasource` {#OpenTargets_get_evidence_by_datasource}

Retrieves Open Targets target-disease evidence filtered by chosen datasources such as GWAS credible sets, EVA or gene burden.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**target-disease evidence, gwas, eva, gene burden, somatic
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `efoId` | string | 否 | Disease EFO/MONDO ID (e.g., 'MONDO_0005011' for Crohn disease). Alternative to disease_name. |
| `ensemblId` | string | 否 | Target Ensembl gene ID (e.g., 'ENSG00000141510' for TP53). Alternative to gene_symbol. |
| `gene_symbol` | string | 否 | HGNC gene symbol (e.g., 'TP53', 'BRCA1'). Auto-resolved to ensemblId. |
| `disease_name` | string | 否 | Disease or phenotype name (e.g., 'Crohn disease'). Auto-resolved to efoId. |
| `datasourceIds` | array of string | 否 | List of datasource IDs to filter evidence. Examples: \['gwas_credible_sets', 'europepmc'\], \['eva', 'clinical_precedence'\], \['intogen', 'cancer_gene_census'\]. For genetic (GWAS) association evidence use 'gwas_credible_sets' — the retired 'ot_genetics_portal' name is remapped automatically. Omit or pass empty array for all sources. |
| `size` | integer | 否 | Maximum evidence rows to return (default: 50) 默认 `50` |

??? note "上游说明（英文原文）"

    Get target-disease evidence from Open Targets filtered by configurable data sources. Unlike OpenTargets_target_disease_evidence (intogen-only), this tool accepts any datasourceIds. Valid sources include: gwas_credible_sets, eva, eva_somatic, gene_burden, genomics_england, intogen, cancer_gene_census, cancer_biomarkers, clinical_precedence, crispr, europepmc, expression_atlas, gene2phenotype, impc, orphanet, uniprot_literature, uniprot_variants, clingen, phewas_catalog, reactome, progeny, slapenrich, sysbio. Note: two source IDs were renamed upstream and the retired names now match nothing rather than erroring — 'ot_genetics_portal' became 'gwas_credible_sets' (confirmed live: IL23R/Crohn disease returns 0 rows for the old name and 43 for the new one) and 'chembl' became 'clinical_precedence'. This tool remaps both retired names automatically and reports the substitution in metadata.datasource_rename_note. Returns evidence rows with scores, literature references, and source attribution.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_evidence_by_datasource', 'arguments': {}},
)
```

## `OpenTargets_get_gwas_study` {#OpenTargets_get_gwas_study}

Retrieves GWAS study metadata (traits, sample sizes, cohorts, summary-stat availability) by GCST accession.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`studyId`
- **关键词：**gcst, gwas catalog, summary statistics
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `studyId` | string | 是 | GWAS study accession ID (e.g., 'GCST000392' for T1D Barrett et al., 'GCST90029024' for T2D Loh et al.) |

??? note "上游说明（英文原文）"

    Get detailed metadata for a GWAS study from Open Targets Platform, including trait information, sample sizes, publication details, disease mappings, cohorts, LD population structure, and availability of summary statistics. Use GCST accession IDs.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_gwas_study',
     'arguments': {'studyId': '...'}},
)
```

## `OpenTargets_get_study_credible_sets` {#OpenTargets_get_study_credible_sets}

Lists credible sets for a GWAS study with lead variants, p-values, effect sizes and L2G scores.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`studyIds`
- **关键词：**fine-mapping, l2g, gcst, loci
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `studyIds` | array of string | 是 | GWAS study accession IDs (e.g., \['GCST000392'\]) |
| `size` | integer | 否 | Number of credible sets to return (default 20) 默认 `20` |
| `index` | integer | 否 | Page index for pagination (0-based, default 0) 默认 `0` |

??? note "上游说明（英文原文）"

    Get all credible sets (fine-mapped GWAS loci) for a specific GWAS study from Open Targets. Returns the lead variants, fine-mapping details, p-values, effect sizes, and locus-to-gene (L2G) predictions for each locus in the study. Use this to explore all significant loci discovered in a GWAS.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_study_credible_sets',
     'arguments': {'studyIds': []}},
)
```

## `OpenTargets_get_variant_credible_sets` {#OpenTargets_get_variant_credible_sets}

Lists GWAS credible sets in which a variant appears, with L2G predictions and effect sizes.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`variantId`
- **关键词：**fine-mapping, l2g, variant id, gwas
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variantId` | string | 是 | Open Targets variant ID in format chr_position_ref_alt (e.g., '10_112998590_C_T' for rs7903146) |
| `size` | integer | 否 | Number of credible sets to return (default 10) 默认 `10` |

??? note "上游说明（英文原文）"

    Get GWAS credible sets (fine-mapped loci) associated with a genetic variant from Open Targets. Returns studies where this variant appears in a fine-mapped credible set, with locus-to-gene (L2G) predictions, p-values, effect sizes, and the fine-mapping method used. Essential for understanding which GWAS signals a variant is part of and which genes it likely affects.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_variant_credible_sets',
     'arguments': {'variantId': '...'}},
)
```

## `OpenTargets_search_gwas_studies_by_disease` {#OpenTargets_search_gwas_studies_by_disease}

Searches Open Targets GWAS studies for a disease or trait by EFO/MONDO IDs.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**gwas, efo, mondo, trait
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `diseaseIds` | array of string | 否 | Disease ontology IDs (e.g., \['MONDO_0005148'\] for type 2 diabetes, \['MONDO_0004975'\] for Alzheimer disease) |
| `enableIndirect` | boolean | 否 | Include studies for child disease terms (default true) 默认 `True` |
| `disease_name` | string | 否 | Disease or trait name for auto-resolution (e.g., 'type 2 diabetes'). Resolved to diseaseIds automatically. |
| `disease` | string | 否 | Alias for disease_name. |
| `trait` | string | 否 | Alias for disease_name. |

??? note "上游说明（英文原文）"

    Search for GWAS studies associated with a disease or trait in Open Targets. Takes disease ontology IDs (EFO, MONDO) and returns matching GWAS studies with their metadata including sample sizes, traits, publications, and associated diseases. Use enableIndirect=true to include child disease terms.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_search_gwas_studies_by_disease',
     'arguments': {}},
)
```

## `OpenTargets_target_disease_evidence` {#OpenTargets_target_disease_evidence}

Retrieves IntOGen somatic cancer driver evidence for a target-disease association.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**intogen, somatic mutation, cancer driver
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `efoId` | string | 否 | EFO/MONDO disease ID (e.g., MONDO_0005011). Alternative to disease_name. |
| `ensemblId` | string | 否 | Ensembl gene ID (e.g., ENSG00000141510). Alternative to gene_symbol. |
| `gene_symbol` | string | 否 | HGNC gene symbol (e.g., 'TP53', 'BRCA1'). Auto-resolved to ensemblId. |
| `disease_name` | string | 否 | Disease or phenotype name (e.g., 'Crohn disease', 'breast carcinoma'). Auto-resolved to efoId. |

??? note "上游说明（英文原文）"

    Explore IntOGen somatic driver evidence for a target-disease association. IMPORTANT: IntOGen only covers somatic tumor driver mutations — this tool returns empty data (\{\}) for non-cancer diseases (e.g., type 2 diabetes, hypertension). Empty results are expected and normal when querying non-oncology phenotypes. For non-cancer diseases or to query other data sources, use OpenTargets_get_evidence_by_datasource instead. Accepts gene_symbol (e.g. 'TP53') and disease_name (e.g. 'Crohn disease') as alternatives to ensemblId/efoId — IDs are auto-resolved via OpenTargets search. Results are PAGINATED and Open Targets returns only 25 evidence rows: `evidences.count` is the total, `evidences.returned` is how many are in `rows`, and `evidences.truncated` is true when evidence remains unfetched, in which case `metadata.truncation_note` says so. This tool exposes no pagination parameter, so use OpenTargets_get_evidence_by_datasource (which takes `size`) when the full evidence set is needed.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_target_disease_evidence', 'arguments': {}},
)
```

## `Orphanet_get_gene_diseases` {#Orphanet_get_gene_diseases}

Lists Orphanet rare diseases associated with a gene symbol or gene name, with association details.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**gene-disease, rare disease, orpha code
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type (fixed: get_gene_diseases) 默认 `get_gene_diseases` |
| `gene_name` | string | 否 | Gene symbol or name keyword (e.g., 'FBN1', 'BRCA1', 'huntingtin'). Use gene_symbol for standard HGNC symbols. Note: Orphanet uses current HGNC symbols (e.g., 'GBA1' not 'GBA'). |
| `gene_symbol` | string \| null | 否 | Alternative to gene_name: standard HGNC gene symbol (e.g., 'GBA1', 'FBN1'). Note: Orphanet uses current HGNC symbols (e.g., 'GBA1' not 'GBA'). |

??? note "上游说明（英文原文）"

    Get rare diseases associated with a gene from Orphanet. Accepts gene symbols (e.g., 'FBN1', 'BRCA1', 'HTT') or gene name keywords (e.g., 'fibrillin', 'huntingtin', 'dystrophin'). Gene symbols are auto-resolved to full names. Returns disease ORPHA codes, gene symbols, association types, and genomic locus information.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Orphanet_get_gene_diseases', 'arguments': {}},
)
```

## `Orphanet_get_genes` {#Orphanet_get_genes}

Lists genes associated with an Orphanet rare disease with association types (causative, modifier, susceptibility).

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**disease genes, causative gene, rare disease, genetic diagnosis
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type (fixed: get_genes) 默认 `get_genes` |
| `orpha_code` | string \| integer | 否 | Orphanet ORPHA code for the disease |
| `orphacode` | integer \| string | 否 | Alias for orpha_code (the Orphanet disease code). |
| `orpha_id` | integer \| string | 否 | Alias for orpha_code (the Orphanet disease code). |

??? note "上游说明（英文原文）"

    Get genes associated with a rare disease from Orphanet. Returns gene symbols, names, and association types (causative, modifier, susceptibility). Essential for rare disease genetic diagnosis.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Orphanet_get_genes', 'arguments': {}},
)
```

## `PanelApp_get_panel` {#PanelApp_get_panel}

Returns all genes in a Genomics England PanelApp panel with confidence level, inheritance mode and phenotypes.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`panel_id`
- **关键词：**gene panel, genomics england, diagnostic genes, inheritance
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `panel_id` | integer | 是 | PanelApp panel ID (e.g., 504 for 'Inherited polyposis and early onset colorectal cancer'). Find IDs using PanelApp_search_panels. |

??? note "上游说明（英文原文）"

    Get the full gene list and details for a specific Genomics England PanelApp gene panel by panel ID. Returns all genes in the panel with their confidence levels (3=green/diagnostic, 2=amber/borderline, 1=red/limited), mode of inheritance, mode of pathogenicity, penetrance, phenotypes, publications, and evidence sources. This is the definitive resource for which genes are included in clinical genetic testing panels used by the NHS. Use PanelApp_search_panels to find panel IDs.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PanelApp_get_panel', 'arguments': {'panel_id': 0}},
)
```

## `PanelApp_search_genes` {#PanelApp_search_genes}

Lists PanelApp panels containing a gene, with confidence level, inheritance and phenotypes.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`gene_symbol`
- **关键词：**gene panel, clinical genetic testing, genomics england
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 是 | Gene symbol to search across all panels (e.g., 'TP53', 'BRCA1', 'EGFR', 'MLH1') |

??? note "上游说明（英文原文）"

    Search for a gene across all Genomics England PanelApp panels to find which clinical testing panels include it. Returns all panels containing the gene, along with the confidence level, mode of inheritance, and associated phenotypes for each panel. This is invaluable for understanding the clinical significance of a gene - if it appears in multiple high-confidence panels, it has strong clinical evidence. Use gene symbols (e.g., 'TP53', 'BRCA1').

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PanelApp_search_genes',
     'arguments': {'gene_symbol': '...'}},
)
```

## `PanelApp_search_panels` {#PanelApp_search_panels}

Searches Genomics England PanelApp gene panels by name or disease group.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`search`
- **关键词：**gene panel, genetic testing, nhs
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `search` | string | 是 | Search query for panel name or disease (e.g., 'breast cancer', 'epilepsy', 'hearing loss', 'intellectual disability') |

??? note "上游说明（英文原文）"

    Search Genomics England PanelApp for gene panels used in clinical genetic testing. PanelApp is the authoritative source for curated gene panels used by the NHS and global genomic medicine services. Returns panel IDs, names, disease groups, status, version, and gene count statistics. Use this to find relevant gene panels for a disease area, then use PanelApp_get_panel to get the full gene list. Panels contain genes with confidence levels (green=3 for diagnostic grade, amber=2 for borderline, red=1 for limited evidence). Matches only against panel-level metadata (name/disease group), not gene-level phenotype text, so some disease names (e.g. 'haemophilia') won't match any panel directly even though a relevant panel exists -- if a disease-name search returns nothing, try PanelApp_search_genes with the causal gene symbol instead (e.g. F8, F9 for haemophilia).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PanelApp_search_panels', 'arguments': {'search': '...'}},
)
```

## `PharmacoDB_get_biomarker_assoc` {#PharmacoDB_get_biomarker_assoc}

Returns gene-compound biomarker associations linking gene expression or mutation to drug sensitivity.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**biomarker, drug sensitivity, pharmacogenomics, aac
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`get_biomarker_associations` 默认 `get_biomarker_associations` |
| `compound_name` | string \| null | 否 | Compound name (e.g., 'Paclitaxel', 'Erlotinib'). Either compound_name or compound_id required. |
| `compound_id` | integer \| null | 否 | PharmacoDB compound ID. Either compound_name or compound_id required. |
| `gene_name` | string \| null | 否 | Gene Ensembl ID to filter (e.g., 'ENSG00000141510' for TP53) |
| `tissue_name` | string \| null | 否 | Tissue type to filter (e.g., 'Breast', 'Lung', 'Skin') |
| `mdata_type` | string \| null | 否 | Molecular data type: 'rna' (gene expression), 'cnv' (copy number), 'mutation' |
| `per_page` | integer | 否 | Number of results per page (default 20, max 100) 默认 `20` |

??? note "上游说明（英文原文）"

    Get gene-compound biomarker associations from PharmacoDB showing how gene expression or other molecular features correlate with drug sensitivity. Returns statistical associations including effect estimate, p-values, FDR, sensitivity statistic (AAC/IC50), molecular data type (rna/cnv/mutation), sample size, and tissue/dataset context. Filter by compound, gene (Ensembl ID), tissue, and molecular data type. Useful for pharmacogenomics biomarker discovery and validating drug sensitivity predictors. Example: find genes associated with Paclitaxel sensitivity in breast cancer.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PharmacoDB_get_biomarker_assoc', 'arguments': {}},
)
```

## `Pharos_get_disease_targets` {#Pharos_get_disease_targets}

Returns Pharos targets associated with a disease with Target Development Level filtering.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`disease`
- **关键词：**tdl, tdark, disease targets
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `disease` | string | 是 | Disease name (e.g., 'breast cancer', 'Alzheimer disease', 'diabetes mellitus') |
| `tdl` | string | 否 | Optional filter by Target Development Level 可选值：`Tdark`, `Tbio`, `Tchem`, `Tclin` |
| `top` | integer | 否 | Maximum number of results (1-100) 默认 `20` |

??? note "上游说明（英文原文）"

    Get drug targets associated with a disease from Pharos. Returns targets with TDL classification to identify druggable vs understudied proteins for the disease. Supports TDL filtering to find Tdark opportunities.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Pharos_get_disease_targets',
     'arguments': {'disease': '...'}},
)
```

## `PomBase_get_gene_phenotypes` {#PomBase_get_gene_phenotypes}

Returns FYPO phenotype annotations and deletion viability for an S. pombe gene.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`gene_id`
- **关键词：**fypo, deletion viability, s. pombe
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_id` | string | 是 | PomBase systematic gene identifier. Examples: 'SPBC11B10.09' (cdc2), 'SPAC2F7.03c' (pom1). |

??? note "上游说明（英文原文）"

    Get phenotype annotations for a fission yeast (S. pombe) gene from PomBase. Returns deletion viability, FYPO (Fission Yeast Phenotype Ontology) term annotations with evidence codes. S. pombe phenotypes are a key resource for understanding gene function in cell biology. Example: SPBC11B10.09 (cdc2) is inviable when deleted, reflecting its essential role in cell division.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PomBase_get_gene_phenotypes',
     'arguments': {'gene_id': '...'}},
)
```

## `reactome_disease_target_score` {#reactome_disease_target_score}

Gets Open Targets target-disease scores from Reactome for a disease.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`efoId`
- **关键词：**open targets, reactome, pathway, target-disease
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `efoId` | string | 是 | The EFO (Experimental Factor Ontology) ID of the disease, e.g., 'MONDO_0011996' for chronic myelogenous leukemia |
| `pageSize` | integer | 否 | Number of results per page (default: 100, max: 100) 默认 `100` |

??? note "上游说明（英文原文）"

    Extract disease-target association scores from Reactome pathway data. This includes pathway-based evidence.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'reactome_disease_target_score',
     'arguments': {'efoId': '...'}},
)
```

## `RGD_get_qtls_in_region` {#RGD_get_qtls_in_region}

Returns RGD QTLs overlapping a genomic region with LOD score and p-value.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`chromosome`, `start`, `stop`, `map_key`
- **关键词：**qtl, lod score, rat, genomic region
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chromosome` | string | 是 | Chromosome name, e.g. '10' or 'chr10'. |
| `start` | integer | 是 | Region start coordinate (1-based). |
| `stop` | integer | 是 | Region stop coordinate. |
| `map_key` | integer | 是 | Assembly map key. 360=rat rn7/GRCr8, 372=rat mRatBN7.2, 38=human GRCh38, 17=human GRCh37, 35=mouse GRCm39, 18=mouse GRCm38. |

??? note "上游说明（英文原文）"

    Get QTLs (Quantitative Trait Loci) overlapping a genomic region from the Rat Genome Database (RGD's flagship data type). Returns QTL symbol, name, RGD ID, LOD score, p-value, and inheritance type. Requires a chromosome, start/stop coordinates, and a map_key identifying the assembly (e.g., 360 for rat rn7/GRCr8, 38 for human GRCh38, 35 for mouse GRCm39). Example: chromosome='10', start=1, stop=50000000, map_key=360 returns QTL Uae41 'Urinary albumin excretion QTL 41' lod=5.26 rgdId=7387235.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'RGD_get_qtls_in_region',
     'arguments': {'chromosome': '...',
                   'start': 0,
                   'stop': 0,
                   'map_key': 0}},
)
```

## `SGD_get_disease_annotations` {#SGD_get_disease_annotations}

Gets curated Disease Ontology associations for a yeast gene from SGD, inferred through orthology to human disease genes.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`sgd_id`
- **关键词：**yeast, saccharomyces cerevisiae, disease ontology, doid, orthology, gene-disease
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sgd_id` | string | 是 | SGD locus identifier. Examples: 'S000000364' (CDC28), 'S000000259' (GAL4). |

??? note "上游说明（英文原文）"

    Get curated gene-disease (Disease Ontology / DOID) associations for a yeast gene from SGD, inferred by orthology to human genes. Returns disease name, DOID, annotation type, evidence, and references. Example: S000000364 (CDC28) is associated with 'cancer' (DOID:162), manually curated.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SGD_get_disease_annotations',
     'arguments': {'sgd_id': '...'}},
)
```

## `SGD_get_phenotypes` {#SGD_get_phenotypes}

Gets curated phenotype annotations for a yeast gene from SGD, including mutant types, methods, strain backgrounds and conditions.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`sgd_id`
- **关键词：**yeast, mutant phenotype, gene-phenotype, strain background, fitness
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sgd_id` | string | 是 | SGD locus identifier. Examples: 'S000003219' (RMR1), 'S000000259' (GAL4), 'S000005564' (ACT1). |

??? note "上游说明（英文原文）"

    Get curated phenotype annotations for a yeast gene from SGD. Returns observable phenotypes (growth, morphology, fitness), mutant types, experimental methods, strain backgrounds, and chemical conditions. Yeast phenotypes are a key resource for functional genomics. Example: S000003219 (RMR1) has 29 phenotype annotations including competitive fitness, respiratory growth defects.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SGD_get_phenotypes', 'arguments': {'sgd_id': '...'}},
)
```

## `SPrediXcan_associate` {#SPrediXcan_associate}

Runs a summary-based S-PrediXcan TWAS test for a gene by combining eQTL prediction weights with GWAS per-SNP z-scores, returning z-score and p-value.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`weight`, `gwas_z`
- **关键词：**twas, eqtl weights, gwas summary statistics, snp z-score, gene-trait association
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `weight` | array of number | 是 | Per-SNP eQTL prediction weights for the gene (same SNP order as gwas_z). |
| `gwas_z` | array of number | 是 | Per-SNP GWAS z-scores (beta/se), same SNP order. |
| `snp_sd` | array \| null | 否 | Per-SNP standard deviation (~sqrt(2*MAF*(1-MAF))); default 1.0 (standardized). |
| `covariance` | array \| null | 否 | Optional SNP covariance matrix (n x n, same SNP order). If omitted, SNPs are assumed independent. |
| `snp` | array \| null | 否 | Optional SNP identifiers (same order), for reference. |

??? note "上游说明（英文原文）"

    Summary-based transcriptome-wide association (S-PrediXcan, Barbeira 2018): test a gene for trait association by combining its eQTL prediction weights with GWAS per-SNP z-scores. Returns the gene TWAS z-score and p-value. Pure-compute (no individual-level data, no LD panel needed when SNPs are independent/LD-pruned; optionally supply a SNP covariance matrix). Pairs with eQTL_get_associations + GWAS summary tools.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SPrediXcan_associate',
     'arguments': {'weight': [], 'gwas_z': []}},
)
```

## `TIMER2_survival_association` {#TIMER2_survival_association}

Tests overall-survival association of high vs low expression of a gene in a TCGA cancer type with a log-rank test.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`cancer`, `gene`
- **关键词：**survival, log-rank, gene expression, tcga, prognostic gene
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`survival_association` 默认 `survival_association` |
| `cancer` | string | 是 | TCGA cancer type (e.g., 'BRCA', 'LUAD', 'COAD') |
| `gene` | string | 是 | Gene symbol for survival analysis (e.g., 'CD8A', 'TP53', 'BRCA1') |

??? note "上游说明（英文原文）"

    Analyze overall survival association of gene expression in TCGA cancer patients using cBioPortal data. Performs log-rank test comparing high vs. low gene expression groups (median split). Returns log-rank p-value, median survival times, and event counts. Note: TIMER2.0/3.0 REST API is no longer public; uses cBioPortal TCGA data.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'TIMER2_survival_association',
     'arguments': {'cancer': '...', 'gene': '...'}},
)
```

## `TPMI_phewas_by_variant` {#TPMI_phewas_by_variant}

Returns phenotypes associated with a variant in the Taiwan Precision Medicine Initiative PheWeb.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**phewas, tpmi, rsid
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `rsid` | string | 否 | dbSNP rsID, e.g. 'rs7903146'. Auto-resolved to GRCh38 coordinates via Ensembl. Provide this OR 'variant'. |
| `variant` | string | 否 | Variant in chr:pos:ref:alt or chr-pos-ref-alt format, in GRCh38, e.g. '10:112998590:C:T'. Provide this OR 'rsid'. A bare rsID is also accepted here. |
| `limit` | integer | 否 | Maximum number of phenotype associations to return (sorted by ascending p-value). 默认 `20` |
| `max_pval` | number \| null | 否 | Optional: only return associations with p-value &lt;= this threshold (e.g. 5e-8 for genome-wide significance). |

??? note "上游说明（英文原文）"

    Phenome-wide association study (PheWAS) lookup for a single variant in the Taiwan Precision Medicine Initiative (TPMI) PheWeb, a large Taiwanese (Han Chinese) biobank. Returns every phenotype associated with the variant, sorted by p-value. Accepts an rsID (auto-resolved to GRCh38 via Ensembl) or an explicit chr:pos:ref:alt variant in GRCh38. Use for cross-ancestry replication in an East Asian (Han Chinese) population.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'TPMI_phewas_by_variant', 'arguments': {}},
)
```

## `UKBTOPMed_phewas_by_variant` {#UKBTOPMed_phewas_by_variant}

Returns phenotypes associated with a variant in the UKB-TOPMed PheWeb.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**无
- **关键词：**phewas, uk biobank, topmed, rsid
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `rsid` | string | 否 | dbSNP rsID, e.g. 'rs7903146'. Auto-resolved to GRCh38 coordinates via Ensembl. Provide this OR 'variant'. |
| `variant` | string | 否 | Variant in chr:pos:ref:alt or chr-pos-ref-alt format, in GRCh38, e.g. '10:112998590:C:T'. Provide this OR 'rsid'. A bare rsID is also accepted here. |
| `limit` | integer | 否 | Maximum number of phenotype associations to return (sorted by ascending p-value). 默认 `20` |
| `max_pval` | number \| null | 否 | Optional: only return associations with p-value &lt;= this threshold (e.g. 5e-8 for genome-wide significance). |

??? note "上游说明（英文原文）"

    Phenome-wide association study (PheWAS) lookup for a single variant in the UKB-TOPMed PheWeb (UK Biobank, ~400K European-ancestry samples, 1,400+ phecode phenotypes). Returns every phenotype associated with the variant, sorted by p-value. Accepts an rsID (auto-resolved to GRCh38 via Ensembl) or an explicit chr:pos:ref:alt variant in GRCh38. Use for variant pleiotropy assessment and replication in a European population.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'UKBTOPMed_phewas_by_variant', 'arguments': {}},
)
```

## `WormBase_get_human_diseases` {#WormBase_get_human_diseases}

Gets human disease models (DOID) for a C. elegans gene inferred by orthology from WormBase.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`gene_id`
- **关键词：**c. elegans, disease ontology, orthology, gene-disease, disease model
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_id` | string | 是 | WBGene identifier or gene name. Examples: 'WBGene00006818', 'unc-86'. |

??? note "上游说明（英文原文）"

    Get human disease models (Disease Ontology / DOID) for a C. elegans gene from WormBase, inferred by orthology to human disease genes. Accepts a WBGene ID or gene name. Example: WBGene00006818 (unc-86) is a potential model for DOID:0110546 'autosomal dominant nonsyndromic deafness 15', inferred by orthology to human HGNC:9220.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'WormBase_get_human_diseases',
     'arguments': {'gene_id': '...'}},
)
```

## `WormBase_get_phenotypes` {#WormBase_get_phenotypes}

Gets observed and not-observed phenotype annotations for a C. elegans gene from RNAi, mutation and transgene experiments.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`gene_id`
- **关键词：**c. elegans, rnai, mutant phenotype, gene-phenotype
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_id` | string | 是 | WormBase gene identifier. Examples: 'WBGene00006763' (unc-26), 'WBGene00000898' (daf-2). |

??? note "上游说明（英文原文）"

    Get phenotype annotations for a C. elegans gene from WormBase. Returns observed phenotypes (from RNAi, mutation, transgene experiments), phenotypes not observed, and total counts. Phenotype data uses the WormBase Phenotype Ontology (WPO). Example: WBGene00006763 (unc-26) has 125+ observed phenotypes including locomotion defects, synaptic vesicle recycling abnormalities, and pharyngeal pumping defects.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'WormBase_get_phenotypes',
     'arguments': {'gene_id': '...'}},
)
```

## `ZFIN_get_gene_phenotypes` {#ZFIN_get_gene_phenotypes}

Gets phenotype annotations observed when a zebrafish gene is mutated or knocked down.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`gene_id`
- **关键词：**zebrafish, mutant phenotype, morpholino, gene-phenotype
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_id` | string | 是 | ZFIN gene ID with 'ZFIN:' prefix. Examples: 'ZFIN:ZDB-GENE-990415-8' (pax2a), 'ZFIN:ZDB-GENE-990415-72' (tp53), 'ZFIN:ZDB-GENE-980526-166' (shha). |
| `limit` | integer | 否 | Maximum number of phenotypes to return (1-100). Default: 20. |
| `page` | integer | 否 | Page number for pagination. Default: 1. |

??? note "上游说明（英文原文）"

    Get phenotype annotations for a zebrafish gene from ZFIN via the Alliance of Genome Resources. Returns phenotypes observed when the gene is mutated or knocked out in zebrafish, including affected anatomical structures and phenotype statements. Zebrafish is widely used for forward genetic screens and disease modeling. Example: 'ZFIN:ZDB-GENE-990415-8' (pax2a) returns phenotypes related to eye, kidney, and nervous system development.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ZFIN_get_gene_phenotypes',
     'arguments': {'gene_id': '...'}},
)
```
