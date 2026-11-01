# ToolUniverse: variant (`tooluniverse`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **845** 个 DNA/基因相关工具（服务器共列出 2927 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

本页列出类别 `variant` 的 155 个工具，其他类别：[`sequence`](tooluniverse-sequence.md) (44), [`genome_annotation`](tooluniverse-genome_annotation.md) (253), [`association`](tooluniverse-association.md) (138), [`expression`](tooluniverse-expression.md) (56), [`regulation_epigenome`](tooluniverse-regulation_epigenome.md) (54), [`sequencing_data`](tooluniverse-sequencing_data.md) (91), [`prediction_model`](tooluniverse-prediction_model.md) (17), [`molecular_tools`](tooluniverse-molecular_tools.md) (23), [`breeding_population`](tooluniverse-breeding_population.md) (14)。

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
| [`Alliance_get_alleles_and_models`](#Alliance_get_alleles_and_models) | `variant` | Retrieves alleles, variants and affected genetic models for a gene from the Alliance of Genome Resources. |
| [`AlphaMissense_get_protein_scores`](#AlphaMissense_get_protein_scores) | `variant` | Retrieves AlphaMissense pathogenicity scores for every residue of a protein by UniProt ID. |
| [`AlphaMissense_get_residue_scores`](#AlphaMissense_get_residue_scores) | `variant` | Returns AlphaMissense benign/ambiguous/pathogenic class breakdown for substitutions at one residue. |
| [`AlphaMissense_get_variant_score`](#AlphaMissense_get_variant_score) | `variant` | Classifies a missense variant (e.g. p.R123H) with AlphaMissense. |
| [`annotate_variant_multi_source`](#annotate_variant_multi_source) | `variant` | Annotates a variant from ClinVar, gnomAD, CIViC and UniProt, returning pathogenicity, population data and clinical evidence. |
| [`CADD_get_position_scores`](#CADD_get_position_scores) | `variant` | Returns CADD PHRED scores for all possible substitutions at a genomic position. |
| [`CADD_get_range_scores`](#CADD_get_range_scores) | `variant` | Returns CADD scores for all variants in a genomic range up to 100 bp. |
| [`CADD_get_variant_score`](#CADD_get_variant_score) | `variant` | Returns the CADD PHRED deleteriousness score for a variant given chrom, pos, ref and alt. |
| [`CancerVar_classify_variant`](#CancerVar_classify_variant) | `variant` | Classifies a somatic variant into AMP/ASCO/CAP Tier I-IV using CancerVar. |
| [`cBioPortal_get_copy_number_alterations`](#cBioPortal_get_copy_number_alterations) | `variant` | Returns GISTIC copy-number alteration calls per sample for genes in a cancer study. |
| [`cBioPortal_get_mutations`](#cBioPortal_get_mutations) | `variant` | Retrieves mutation data for given genes in a cBioPortal cancer study. |
| [`civic_get_assertion`](#civic_get_assertion) | `variant` | Retrieves a CIViC assertion integrating evidence on a variant's clinical actionability. |
| [`civic_get_evidence_item`](#civic_get_evidence_item) | `variant` | Retrieves a CIViC evidence item linking a variant or molecular profile to clinical outcomes. |
| [`civic_get_molecular_profile`](#civic_get_molecular_profile) | `variant` | Retrieves a CIViC molecular profile (variant combination) with its evidence. |
| [`civic_get_variant`](#civic_get_variant) | `variant` | Retrieves a CIViC variant with coordinates, gene and clinical significance. |
| [`civic_get_variants_by_gene`](#civic_get_variants_by_gene) | `variant` | Lists all CIViC variants for a given gene. |
| [`civic_search_assertions`](#civic_search_assertions) | `variant` | Searches CIViC assertions by therapy, disease or variant name. |
| [`civic_search_evidence_items`](#civic_search_evidence_items) | `variant` | Searches CIViC evidence items by disease, therapy, molecular profile and significance. |
| [`civic_search_molecular_profiles`](#civic_search_molecular_profiles) | `variant` | Searches CIViC molecular profiles by name. |
| [`civic_search_variants`](#civic_search_variants) | `variant` | Searches CIViC variants by variant name or gene. |
| [`ClinGen_dosage_by_gene`](#ClinGen_dosage_by_gene) | `variant` | Returns ClinGen haploinsufficiency and triplosensitivity scores with coordinates for a gene. |
| [`ClinGen_dosage_region_search`](#ClinGen_dosage_region_search) | `variant` | Searches ClinGen dosage sensitivity curations overlapping chromosomal coordinates. |
| [`ClinGen_get_dosage_sensitivity`](#ClinGen_get_dosage_sensitivity) | `variant` | Returns ClinGen haploinsufficiency and triplosensitivity curations, optionally with CNV regions. |
| [`ClinGen_get_variant_classifications`](#ClinGen_get_variant_classifications) | `variant` | Returns ClinGen Evidence Repository variant classifications filtered by gene or variant. |
| [`ClinGen_search_dosage_sensitivity`](#ClinGen_search_dosage_sensitivity) | `variant` | Searches ClinGen dosage sensitivity curations by gene symbol. |
| [`ClinGenAllele_get_allele`](#ClinGenAllele_get_allele) | `variant` | Retrieves a ClinGen Allele Registry canonical allele with HGVS and database cross-references. |
| [`ClinGenAllele_lookup_hgvs`](#ClinGenAllele_lookup_hgvs) | `variant` | Looks up an HGVS expression in the ClinGen Allele Registry to get its CA ID and cross-references. |
| [`ClinGenAR_get_external_records`](#ClinGenAR_get_external_records) | `variant` | Returns external database records (ClinVar, dbSNP, gnomAD, etc.) for a ClinGen allele. |
| [`ClinGenAR_lookup_allele`](#ClinGenAR_lookup_allele) | `variant` | Looks up a variant by HGVS in the ClinGen Allele Registry and returns its CA ID. |
| [`ClinGenAR_lookup_by_external_id`](#ClinGenAR_lookup_by_external_id) | `variant` | Resolves a dbSNP rsID or ClinVar VariationID to a ClinGen canonical allele with cross-references. |
| [`ClinVar_get_clinical_significance`](#ClinVar_get_clinical_significance) | `variant` | Retrieves pathogenicity classification and clinical interpretations for a ClinVar variant. |
| [`ClinVar_get_submitted_records`](#ClinVar_get_submitted_records) | `variant` | Retrieves individual submitter SCV assertions for a ClinVar variant. |
| [`ClinVar_get_variant_details`](#ClinVar_get_variant_details) | `variant` | Retrieves ClinVar variant summary including genes, location and review status. |
| [`ClinVar_search_by_region`](#ClinVar_search_by_region) | `variant` | Finds ClinVar variants overlapping a genomic region with optional significance filter. |
| [`ClinVar_search_variants`](#ClinVar_search_variants) | `variant` | Searches ClinVar variants by gene, condition, ID or variant name. |
| [`coding_variant_fraction`](#coding_variant_fraction) | `variant` | Computes the fraction of coding variants in an Excel/CSV table matching a Sequence Ontology consequence and below a VAF threshold. |
| [`COSMIC_get_mutations_by_gene`](#COSMIC_get_mutations_by_gene) | `variant` | Gets COSMIC somatic mutations for a gene, with CDS and amino-acid changes, cancer types and genomic positions. |
| [`COSMIC_search_mutations`](#COSMIC_search_mutations) | `variant` | Searches COSMIC somatic mutations by gene, mutation ID or description, returning IDs and CDS/AA changes. |
| [`CPIC_get_alleles`](#CPIC_get_alleles) | `variant` | Gets CPIC alleles for a pharmacogene with activity values and functional status. |
| [`CPIC_get_recommendations`](#CPIC_get_recommendations) | `variant` | Gets CPIC guideline recommendations linking pharmacogene genotype or phenotype to prescribing actions. |
| [`dbsnp_get_frequencies`](#dbsnp_get_frequencies) | `variant` | Gets population-specific allele frequencies for an rsID from dbSNP. |
| [`dbsnp_get_variant_by_rsid`](#dbsnp_get_variant_by_rsid) | `variant` | Gets dbSNP variant information by rsID with GRCh38/GRCh37 coordinates and alleles. |
| [`dbsnp_search_by_gene`](#dbsnp_search_by_gene) | `variant` | Searches dbSNP for variants in a gene by symbol. |
| [`EBIProteins_get_variation`](#EBIProteins_get_variation) | `variant` | Gets known sequence variants for a UniProt protein from COSMIC, ClinVar, gnomAD/ExAC and UniProt. |
| [`EBIProteins_get_variation_by_dbsnp`](#EBIProteins_get_variation_by_dbsnp) | `variant` | Returns UniProt proteins and residue changes carrying a given dbSNP rsID. |
| [`EBIProteins_get_variation_by_hgvs`](#EBIProteins_get_variation_by_hgvs) | `variant` | Returns overlapping UniProt proteins and protein-level consequences for an HGVS genomic expression. |
| [`ensembl_get_structural_variants`](#ensembl_get_structural_variants) | `variant` | Gets known structural variants overlapping a region from Ensembl (DGVa, ClinGen). |
| [`ensembl_get_sv_detail`](#ensembl_get_sv_detail) | `variant` | Gets details of an Ensembl structural variant by nsv/esv accession. |
| [`Ensembl_get_transcript_haplotypes`](#Ensembl_get_transcript_haplotypes) | `variant` | Gets 1000 Genomes protein and CDS haplotypes for an Ensembl transcript with frequencies and variants. |
| [`ensembl_get_variants`](#ensembl_get_variants) | `variant` | Gets SNPs and indels overlapping a genomic region from Ensembl. |
| [`ensembl_get_variation`](#ensembl_get_variation) | `variant` | Gets Ensembl details for a variant ID such as rs699, including alleles, frequencies and consequences. |
| [`ensembl_vep_region`](#ensembl_vep_region) | `variant` | Runs Ensembl VEP on a region and allele, returning consequences and SIFT/PolyPhen/CADD scores. |
| [`EnsemblVar_get_population_frequencies`](#EnsemblVar_get_population_frequencies) | `variant` | Gets gnomAD and 1000 Genomes allele frequencies for a variant via Ensembl. |
| [`EnsemblVar_get_variant_consequences`](#EnsemblVar_get_variant_consequences) | `variant` | Gets an Ensembl variant record with consequence, clinical significance and synonyms. |
| [`EnsemblVEP_annotate_hgvs`](#EnsemblVEP_annotate_hgvs) | `variant` | Runs Ensembl VEP on an HGVS notation, returning consequences and pathogenicity scores. |
| [`EnsemblVEP_annotate_rsid`](#EnsemblVEP_annotate_rsid) | `variant` | Runs Ensembl VEP on a dbSNP rsID, returning consequences and pathogenicity scores. |
| [`EnsemblVEP_variant_recoder`](#EnsemblVEP_variant_recoder) | `variant` | Converts a variant ID to HGVS genomic/coding/protein and SPDI notations. |
| [`EVA_get_clustered_variant_by_rs`](#EVA_get_clustered_variant_by_rs) | `variant` | Resolves an EVA clustered-variant RS accession to assembly, contig and position. |
| [`EVA_get_variants_by_gene`](#EVA_get_variants_by_gene) | `variant` | Gets EVA variants for a gene with allele frequencies and consequences. |
| [`EVA_get_variants_by_region`](#EVA_get_variants_by_region) | `variant` | Gets EVA variants within a genomic region. |
| [`EVA_list_studies`](#EVA_list_studies) | `variant` | Lists EVA variant studies available for a species. |
| [`EVE_get_gene_info`](#EVE_get_gene_info) | `variant` | Checks whether EVE pathogenicity scores exist for a gene and returns its Ensembl ID. |
| [`EVE_get_variant_score`](#EVE_get_variant_score) | `variant` | Gets the EVE pathogenicity score for a variant via Ensembl VEP. |
| [`FAVOR_annotate_variant`](#FAVOR_annotate_variant) | `variant` | Annotates a GRCh38 variant via FAVOR with frequencies, consequences and functional scores. |
| [`FlyBase_get_gene_alleles`](#FlyBase_get_gene_alleles) | `variant` | Lists alleles/variants for a Drosophila gene with categories, phenotypes, disease annotations and genomic coordinates. |
| [`GDC_get_cnv_data`](#GDC_get_cnv_data) | `variant` | Queries GDC/TCGA copy number variation segment and gene-level files by project and gene symbol. |
| [`GDC_get_mutation_frequency`](#GDC_get_mutation_frequency) | `variant` | Returns total pan-cancer somatic mutation (SSM) occurrence count for a gene across GDC/TCGA projects. |
| [`GDC_get_mutation_frequency_by_project`](#GDC_get_mutation_frequency_by_project) | `variant` | Returns per-project mutated case counts and total case counts for a gene, giving per-cancer mutation rates. |
| [`GDC_get_ssm_by_gene`](#GDC_get_ssm_by_gene) | `variant` | Gets somatic simple mutations for a gene across GDC/TCGA projects with type, genomic coordinates and amino acid change. |
| [`GeneBe_classify_variant`](#GeneBe_classify_variant) | `variant` | Classifies a germline variant (chr, pos, ref, alt; hg38/hg19) under ACMG/AMP returning verdict, score and criteria. |
| [`GeneBe_classify_variants_batch`](#GeneBe_classify_variants_batch) | `variant` | Batch-classifies up to 1000 germline variants under ACMG/AMP guidelines via GeneBe. |
| [`GenomeNexus_annotate_dbsnp`](#GenomeNexus_annotate_dbsnp) | `variant` | Annotates a dbSNP rsID via Genome Nexus with VEP, SIFT, PolyPhen-2, AlphaMissense and cancer hotspot data. |
| [`GenomeNexus_annotate_mutation`](#GenomeNexus_annotate_mutation) | `variant` | Annotates a mutation given chromosome, start, end, ref and alt alleles (GRCh37) via Genome Nexus. |
| [`GenomeNexus_annotate_variant`](#GenomeNexus_annotate_variant) | `variant` | Annotates a GRCh37 HGVSg variant with aggregated VEP, SIFT, PolyPhen-2, AlphaMissense and hotspot annotations. |
| [`GenomeNexus_get_cancer_hotspots`](#GenomeNexus_get_cancer_hotspots) | `variant` | Returns cancer mutation hotspot annotations (residue, tumor count, type) for an HGVSg variant. |
| [`GIAB_list_directory`](#GIAB_list_directory) | `variant` | Browses NIST Genome in a Bottle release tree of high-confidence benchmark VCF and BED files for HG001-HG007. |
| [`gnomad_get_constraint`](#gnomad_get_constraint) | `variant` | Gets a curated gnomAD constraint summary (pLI, oe_lof, LOEUF, missense Z) for a gene symbol or Ensembl ID. |
| [`gnomad_get_gene_constraints`](#gnomad_get_gene_constraints) | `variant` | Gets gnomAD/ExAC gene constraint metrics (pLI, oe_lof, oe_mis, expected/observed counts) for a gene. |
| [`gnomad_get_region`](#gnomad_get_region) | `variant` | Returns genes overlapping a genomic interval and gnomAD variant IDs in it for a chosen dataset. |
| [`gnomad_get_sv_by_gene`](#gnomad_get_sv_by_gene) | `variant` | Gets gnomAD structural variants for a gene from the v4 (GRCh38) or v2.1 (GRCh37) SV callset. |
| [`gnomad_get_sv_by_region`](#gnomad_get_sv_by_region) | `variant` | Gets gnomAD structural variants overlapping a genomic region in the matching assembly callset. |
| [`gnomad_get_sv_detail`](#gnomad_get_sv_detail) | `variant` | Gets detailed info for a gnomAD structural variant by ID, auto-selecting v4 or v2.1 callset. |
| [`gnomad_get_variant`](#gnomad_get_variant) | `variant` | Gets gnomAD variant metadata and aggregate allele frequency by variant ID (chrom-pos-ref-alt). |
| [`gnomad_get_variant_populations`](#gnomad_get_variant_populations) | `variant` | Gets population-stratified allele frequencies (afr, nfe, eas, sas, etc.) for a gnomAD variant. |
| [`gnomad_search_variants`](#gnomad_search_variants) | `variant` | Searches gnomAD variants by free text such as an rsID and returns matching variant IDs. |
| [`gProfiler_annotate_snps`](#gProfiler_annotate_snps) | `variant` | Maps SNP rsIDs to positions, Ensembl genes and predicted variant consequences using g:SNPense. |
| [`InterVar_classify_variant`](#InterVar_classify_variant) | `variant` | Classifies a germline variant given chrom/pos/ref/alt with InterVar ACMG/AMP 2015 rules and reports evidence criteria. |
| [`IPD_get_hla_allele`](#IPD_get_hla_allele) | `variant` | Fetches a full IPD-IMGT/HLA allele record by accession including name, locus, class, sequence and release history. |
| [`IPD_search_cells`](#IPD_search_cells) | `variant` | Searches IPD HLA reference cell lines by name or other fields and returns cell IDs and names. |
| [`IPD_search_hla_alleles`](#IPD_search_hla_alleles) | `variant` | Searches IPD-IMGT/HLA for human HLA alleles by allele name prefix or substring and returns accessions and names. |
| [`KEGG_get_variant`](#KEGG_get_variant) | `variant` | Returns a KEGG variant entry with gene, mutation type and ClinVar/dbSNP/COSMIC cross-references, networks and drugs. |
| [`KEGG_search_variant`](#KEGG_search_variant) | `variant` | Searches the KEGG VARIANT database for human gene variants associated with disease. |
| [`LitVar_get_variant_details`](#LitVar_get_variant_details) | `variant` | Returns the LitVar2 record for an rsID with ClinGen allele IDs, associated genes, HGVS names and clinical annotations. |
| [`LitVar_search_variants`](#LitVar_search_variants) | `variant` | Searches LitVar2 variants by rsID, gene name or HGVS and returns identifiers, genes, HGVS and publication counts. |
| [`LOVD_get_variants`](#LOVD_get_variants) | `variant` | Returns all curated LOVD variants for a gene with HGVS DNA/RNA/protein notation and hg19 positions. |
| [`LOVD_search_variants`](#LOVD_search_variants) | `variant` | Searches LOVD variants by variant DBID or HGVS DNA notation. |
| [`MaveDB_get_clinical_controls`](#MaveDB_get_clinical_controls) | `variant` | Returns ClinVar pathogenic/benign/VUS variants overlapping a MaveDB score set for score calibration. |
| [`MaveDB_get_effect_matrix`](#MaveDB_get_effect_matrix) | `variant` | Loads a MaveDB score set and reshapes single missense variant scores into an amino acid by position effect matrix. |
| [`MaveDB_get_gnomad_variants`](#MaveDB_get_gnomad_variants) | `variant` | Returns gnomAD variants with allele counts and frequencies observed within a MaveDB score set. |
| [`MaveDB_get_mapped_variants`](#MaveDB_get_mapped_variants) | `variant` | Returns GA4GH VRS genomic allele mappings and ClinGen allele IDs for variants in a MaveDB score set. |
| [`MaveDB_get_score_set`](#MaveDB_get_score_set) | `variant` | Returns MaveDB score set details including target gene, methods, variant count and publications. |
| [`MaveDB_get_variant_scores`](#MaveDB_get_variant_scores) | `variant` | Returns HGVS-annotated variant functional scores from a MaveDB score set. |
| [`MaveDB_search_experiments`](#MaveDB_search_experiments) | `variant` | Searches MaveDB experiments grouping related variant effect score sets. |
| [`MaveDB_search_score_sets`](#MaveDB_search_score_sets) | `variant` | Searches MaveDB score sets from MAVE and deep mutational scanning experiments. |
| [`Mutalyzer_back_translate`](#Mutalyzer_back_translate) | `variant` | Converts a protein-level HGVS variant to possible coding DNA variant descriptions with Mutalyzer. |
| [`Mutalyzer_normalize_variant`](#Mutalyzer_normalize_variant) | `variant` | Normalizes an HGVS variant with Mutalyzer and returns corrected notation, protein effect and gene. |
| [`Mutalyzer_parse_hgvs`](#Mutalyzer_parse_hgvs) | `variant` | Parses an HGVS variant into reference, coordinate system, variant type and positions. |
| [`MyVariant_get_pathogenicity_scores`](#MyVariant_get_pathogenicity_scores) | `variant` | Retrieves dbNSFP pathogenicity scores (REVEL, CADD, SIFT, PolyPhen-2, etc.) for a variant. |
| [`MyVariant_get_variant_annotation`](#MyVariant_get_variant_annotation) | `variant` | Retrieves MyVariant.info annotation for a variant by HGVS ID or rsID. |
| [`MyVariant_query_variants`](#MyVariant_query_variants) | `variant` | Searches MyVariant.info for variants by rsID, gene, position or clinical significance. |
| [`NCBIVariation_alfa_frequencies_by_rsid`](#NCBIVariation_alfa_frequencies_by_rsid) | `variant` | Returns ALFA per-population allele counts and frequencies for a dbSNP rsID. |
| [`NCBIVariation_hgvs_to_spdi`](#NCBIVariation_hgvs_to_spdi) | `variant` | Converts and validates an HGVS variant description to SPDI using NCBI Variation Services. |
| [`NCBIVariation_rsid_lookup`](#NCBIVariation_rsid_lookup) | `variant` | Returns dbSNP rsID details with GRCh38 placements, genes, ClinVar significance and MANE transcripts. |
| [`NCBIVariation_spdi_canonical`](#NCBIVariation_spdi_canonical) | `variant` | Returns the canonical normalized SPDI form of a variant. |
| [`NCBIVariation_spdi_equivalents`](#NCBIVariation_spdi_equivalents) | `variant` | Returns all equivalent SPDI representations of a variant across GRCh37/38, RefSeqGene and transcripts. |
| [`NCBIVariation_spdi_to_hgvs`](#NCBIVariation_spdi_to_hgvs) | `variant` | Converts an SPDI variant to HGVS notation using NCBI Variation Services. |
| [`NCBIVariation_spdi_to_rsids`](#NCBIVariation_spdi_to_rsids) | `variant` | Returns dbSNP rsIDs co-located with a normalized SPDI variant. |
| [`NCBIVariation_vcf_to_spdi`](#NCBIVariation_vcf_to_spdi) | `variant` | Converts a VCF chrom/pos/ref/alt variant to a normalized contextual SPDI on RefSeq. |
| [`OncoKB_annotate_copy_number`](#OncoKB_annotate_copy_number) | `variant` | Annotates gene amplification or deletion for oncogenicity and therapies with OncoKB. |
| [`OncoKB_annotate_mutations`](#OncoKB_annotate_mutations) | `variant` | Annotates a cancer mutation such as BRAF V600E for oncogenicity, effect and therapies with OncoKB. |
| [`OncoKB_annotate_variant`](#OncoKB_annotate_variant) | `variant` | Annotates a cancer variant with evidence level, oncogenicity, mutation effect and FDA-approved therapies. |
| [`OpenCRAVAT_annotate_variant`](#OpenCRAVAT_annotate_variant) | `variant` | Annotates a single genomic variant with OpenCRAVAT sources such as ClinVar, gnomAD, REVEL, AlphaMissense and SpliceAI. |
| [`OpenTargets_drug_pharmacogenomics_data`](#OpenTargets_drug_pharmacogenomics_data) | `variant` | Retrieves Open Targets pharmacogenomics records for a drug by ChEMBL ID, including genotype annotations and evidence levels. |
| [`OpenTargets_get_target_constraint_info_by_ensemblID`](#OpenTargets_get_target_constraint_info_by_ensemblID) | `variant` | Retrieves genetic constraint metrics (observed/expected, LOEUF-type scores) for a gene by Ensembl ID. |
| [`OpenTargets_get_variant_effect_predictions`](#OpenTargets_get_variant_effect_predictions) | `variant` | Retrieves variant effect predictions (VEP, SIFT, PolyPhen, AlphaMissense, CADD, GERP) for a chrom_pos_ref_alt variant. |
| [`OpenTargets_get_variant_info`](#OpenTargets_get_variant_info) | `variant` | Returns variant details including coordinates, rsIDs, gnomAD allele frequencies, consequences and HGVS. |
| [`OpenTargets_get_variant_pharmacogenomics`](#OpenTargets_get_variant_pharmacogenomics) | `variant` | Retrieves PharmGKB-derived pharmacogenomic annotations (drug, effect, evidence level, genotype) for a variant. |
| [`OpenTargets_get_variant_transcript_consequences`](#OpenTargets_get_variant_transcript_consequences) | `variant` | Returns per-transcript consequences of a variant: gene, amino acid change, impact and consequence score. |
| [`Pathoplexus_get_mutations`](#Pathoplexus_get_mutations) | `variant` | Returns characteristic nucleotide or amino-acid mutations for a Pathoplexus organism above a proportion threshold. |
| [`PGSCatalog_get_performance_metrics`](#PGSCatalog_get_performance_metrics) | `variant` | Returns published evaluation metrics (OR, HR, AUC, sample sets) for a PGS Catalog score. |
| [`PGSCatalog_get_score`](#PGSCatalog_get_score) | `variant` | Returns metadata of a PGS Catalog polygenic score: traits, variant count, method, build and ancestry. |
| [`PGSCatalog_get_scores_by_trait`](#PGSCatalog_get_scores_by_trait) | `variant` | Lists published polygenic scores for an EFO/MONDO trait from the PGS Catalog. |
| [`PGSCatalog_search_traits`](#PGSCatalog_search_traits) | `variant` | Searches PGS Catalog traits by term and returns trait IDs with polygenic score counts. |
| [`PharmGKB_get_clinical_annotations`](#PharmGKB_get_clinical_annotations) | `variant` | Returns PharmGKB clinical annotations with variant-level gene-drug-phenotype significance by annotation ID. |
| [`PharmGKB_get_dosing_guidelines`](#PharmGKB_get_dosing_guidelines) | `variant` | Returns CPIC/DPWG pharmacogenetic dosing guidelines from PharmGKB by guideline ID. |
| [`PharmGKB_get_variant_annotations`](#PharmGKB_get_variant_annotations) | `variant` | Returns PharmGKB variant annotations (allele/genotype, phenotype, significance) filtered by gene or chemical. |
| [`PharmGKB_search_variants`](#PharmGKB_search_variants) | `variant` | Searches PharmGKB genetic variants by rsID or name with gene info. |
| [`PharmVar_get_alleles`](#PharmVar_get_alleles) | `variant` | Returns PharmVar star allele definitions, function, activity scores and variant details for a pharmacogene. |
| [`PharmVar_list_genes`](#PharmVar_list_genes) | `variant` | Lists pharmacogenes curated by PharmVar. |
| [`Progenetix_cnv_search`](#Progenetix_cnv_search) | `variant` | Searches cancer biosamples with CNVs (gains/losses) in a genomic region and disease filter. |
| [`Progenetix_get_cnv_frequencies`](#Progenetix_get_cnv_frequencies) | `variant` | Returns genome-wide gain and loss frequencies per cytoband for a cancer type from Progenetix. |
| [`proteins_api_get_variants`](#proteins_api_get_variants) | `variant` | Returns natural, disease-associated and population variants for UniProt accessions, paginated. |
| [`ProtVar_get_function`](#ProtVar_get_function) | `variant` | Returns UniProt functional features and structural context at an amino-acid position, optionally for a variant residue. |
| [`ProtVar_get_population`](#ProtVar_get_population) | `variant` | Returns co-located variants with gnomAD/1000G frequencies, ClinVar significance and SIFT/PolyPhen at a protein position. |
| [`ProtVar_map_variant`](#ProtVar_map_variant) | `variant` | Maps a protein change, rsID or VCF-style variant to genomic coordinates with AlphaMissense, EVE and ESM scores. |
| [`StarAlleles_search`](#StarAlleles_search) | `variant` | Autocompletes pharmacogenomic star alleles (e.g. CYP2D6*4) with nucleotide and protein changes. |
| [`UniProt_get_disease_variants_by_accession`](#UniProt_get_disease_variants_by_accession) | `variant` | Extracts all VARIANT features and their disease annotations from a UniProtKB entry. |
| [`VariantValidator_format_genomic_to_transcripts`](#VariantValidator_format_genomic_to_transcripts) | `variant` | Projects a genomic HGVS variant onto every overlapping RefSeq transcript with c./p. descriptions via VariantFormatter. |
| [`VariantValidator_validate_variant`](#VariantValidator_validate_variant) | `variant` | Validates an HGVS variant description and returns genomic, transcript and protein-level representations. |
| [`VCF_count_variants`](#VCF_count_variants) | `variant` | Counts records in a local VCF/BCF with bcftools after filters, with SNP/indel/MNP breakdown. |
| [`VCF_normalize`](#VCF_normalize) | `variant` | Normalizes a local VCF/BCF with bcftools norm and reports counts before and after. |
| [`VCF_summary_stats`](#VCF_summary_stats) | `variant` | Computes bcftools stats for a local VCF/BCF: SNPs, indels, ts/tv ratio and per-sample genotype counts. |
| [`VEP_predict_pathogenicity`](#VEP_predict_pathogenicity) | `variant` | Predicts missense variant pathogenicity with AlphaMissense, SIFT and PolyPhen-2 via Ensembl VEP for an rsID, HGVS or genomic variant. |
| [`ZFIN_get_allele`](#ZFIN_get_allele) | `variant` | Gets a zebrafish allele record from ZFIN via the Alliance API with associated gene, alteration type and synonyms. |
| [`ZFIN_get_gene_alleles`](#ZFIN_get_gene_alleles) | `variant` | Gets mutant alleles and variants of a zebrafish gene with variant types, molecular consequences and phenotypes. |

## `Alliance_get_alleles_and_models` {#Alliance_get_alleles_and_models}

Retrieves alleles, variants and affected genetic models for a gene from the Alliance of Genome Resources.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`gene_id`
- **关键词：**allele, mutant strain, hgvs, model organism, phenotype
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_id` | string | 是 | Alliance gene ID with member prefix, e.g. 'MGI:97490' (Pax6), 'HGNC:8620'. |
| `limit` | integer \| null | 否 | Max results per category (default 20, max 100). |

??? note "上游说明（英文原文）"

    Get alleles/variants and affected genomic models (mutant strains) for a gene via the Alliance of Genome Resources (mouse/MGI, frog/Xenbase, and other members). Returns allele/variant records with HGVS genomic locations and disease/phenotype flags, plus affected genomic models (genotypes/strains). Example: MGI:97490 (Pax6) returns variant_summary records such as rs231853465 (NC_000068.8:g.105516528T&gt;C) and models such as AffectedGenomicModel MGI:2680573.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Alliance_get_alleles_and_models',
     'arguments': {'gene_id': '...'}},
)
```

## `AlphaMissense_get_protein_scores` {#AlphaMissense_get_protein_scores}

Retrieves AlphaMissense pathogenicity scores for every residue of a protein by UniProt ID.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`uniprot_id`
- **关键词：**missense, pathogenicity, alphamissense, residue scores
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `uniprot_id` | string | 是 | UniProt accession ID (e.g., 'P00533' for EGFR, 'P01116' for KRAS, 'P04637' for TP53). |
| `max_residues` | integer \| null | 否 | Optional cap on the number of residues fetched. 0 / null / omitted = no cap (fetch every residue of the protein, recommended for whole-protein analysis). Positive integer = fetch only the first N residues (useful for fast diagnostic calls on very long proteins). 默认 `0` |

??? note "上游说明（英文原文）"

    Get AlphaMissense pathogenicity scores for EVERY residue of a protein in one call. The upstream hegelab API only exposes per-residue queries; this tool hides that by looking up the protein length from UniProt and concurrently fetching all residues via a thread pool. Returns a list of per-position records (one per residue) plus protein length, count of successful fetches, and the pathogenic/ambiguous/benign thresholds. Use this when you want whole-protein vulnerability assessment, hotspot detection, or to build a (20 × n_positions) AlphaMissense score matrix for benchmarking. Set max_residues &gt; 0 to cap the loop for fast diagnostic calls.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'AlphaMissense_get_protein_scores',
     'arguments': {'uniprot_id': '...'}},
)
```

## `AlphaMissense_get_residue_scores` {#AlphaMissense_get_residue_scores}

Returns AlphaMissense benign/ambiguous/pathogenic class breakdown for substitutions at one residue.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`uniprot_id`, `position`
- **关键词：**missense, substitution, pathogenicity class, snv-accessible
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `uniprot_id` | string | 是 | UniProt accession ID (e.g., 'P00533' for EGFR) |
| `position` | integer | 是 | Amino acid position in the protein (1-indexed) |

??? note "上游说明（英文原文）"

    Get the AlphaMissense class breakdown for substitutions at a specific protein position: which residues fall in each class (benign/ambiguous/pathogenic), reported separately for the SNV-accessible substitutions and for all 19, plus the residue mean scores. NOTE: the upstream hotspot API does not expose a numeric score per individual substitution, so this returns class membership rather than 20 scores.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'AlphaMissense_get_residue_scores',
     'arguments': {'uniprot_id': '...', 'position': 0}},
)
```

## `AlphaMissense_get_variant_score` {#AlphaMissense_get_variant_score}

Classifies a missense variant (e.g. p.R123H) with AlphaMissense.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`uniprot_id`, `variant`
- **关键词：**missense variant, pathogenicity, protein change, hgvs p
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `uniprot_id` | string | 是 | UniProt accession ID (e.g., 'P00533' for EGFR) |
| `variant` | string | 是 | Variant in protein notation: p.X123Y or X123Y where X is reference amino acid, 123 is position, Y is variant (e.g., 'p.R123H', 'V600E') |

??? note "上游说明（英文原文）"

    Classify a specific missense variant with AlphaMissense. Input: UniProt ID and variant (e.g., 'p.R123H'). Returns the variant's AlphaMissense class (pathogenic/ambiguous/benign) and whether it came from the SNV-accessible or the all-substitutions set. NOTE: the upstream hotspot API publishes per-residue class membership, not a numeric score per substitution, so `pathogenicity_score` is null and `score_available` is false; the `residue_mean_score_*` fields are residue-level averages, NOT this variant's score. For a numeric per-variant score use GenomeNexus_annotate_variant.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'AlphaMissense_get_variant_score',
     'arguments': {'uniprot_id': '...', 'variant': '...'}},
)
```

## `annotate_variant_multi_source` {#annotate_variant_multi_source}

Annotates a variant from ClinVar, gnomAD, CIViC and UniProt, returning pathogenicity, population data and clinical evidence.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**variant annotation, clinvar, gnomad, civic, pathogenicity, rsid
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant` | string \| null | 否 | Variant name (e.g., 'BRAF V600E', 'TP53 R175H'). Gene symbol + amino acid change. |
| `gene` | string \| null | 否 | Gene symbol (e.g., 'BRCA1'). Used for gnomAD/CIViC lookup. |
| `rsid` | string \| null | 否 | dbSNP rsID (e.g., 'rs113488022'). |

??? note "上游说明（英文原文）"

    Annotate a genetic variant from ClinVar, gnomAD, CIViC, and UniProt in a single call. Returns pathogenicity classification, population constraints, clinical evidence, and protein function. Use when you need comprehensive variant interpretation from multiple databases at once.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'annotate_variant_multi_source', 'arguments': {}},
)
```

## `CADD_get_position_scores` {#CADD_get_position_scores}

Returns CADD PHRED scores for all possible substitutions at a genomic position.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`chrom`, `pos`
- **关键词：**cadd, phred, deleteriousness, saturation
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chrom` | string | 是 | Chromosome (1-22, X, Y, MT) |
| `pos` | integer | 是 | Genomic position (1-based) |
| `version` | string | 否 | CADD version and genome build 默认 `GRCh38-v1.7` |

??? note "上游说明（英文原文）"

    Get CADD scores for all possible substitutions at a genomic position. Returns PHRED scores for A, C, G, T variants at the position. Useful for saturation mutagenesis analysis.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CADD_get_position_scores',
     'arguments': {'chrom': '...', 'pos': 0}},
)
```

## `CADD_get_range_scores` {#CADD_get_range_scores}

Returns CADD scores for all variants in a genomic range up to 100 bp.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`chrom`, `start`, `end`
- **关键词：**cadd, region scan, deleteriousness
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chrom` | string | 是 | Chromosome (1-22, X, Y, MT) |
| `start` | integer | 是 | Start position (1-based, inclusive) |
| `end` | integer | 是 | End position (1-based, inclusive). Max range: 100bp. |
| `version` | string | 否 | CADD version and genome build 默认 `GRCh38-v1.7` |

??? note "上游说明（英文原文）"

    Get CADD scores for all variants in a genomic range (max 100bp). Returns all pre-computed deleteriousness scores in the region. Use for hotspot analysis or exon scanning.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CADD_get_range_scores',
     'arguments': {'chrom': '...', 'start': 0, 'end': 0}},
)
```

## `CADD_get_variant_score` {#CADD_get_variant_score}

Returns the CADD PHRED deleteriousness score for a variant given chrom, pos, ref and alt.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`chrom`, `pos`, `ref`, `alt`
- **关键词：**cadd, phred score, pathogenicity, grch38
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chrom` | string | 是 | Chromosome (1-22, X, Y, MT). Can include 'chr' prefix. |
| `pos` | integer | 是 | Genomic position (1-based) |
| `ref` | string | 是 | Reference allele (e.g., 'A', 'G') |
| `alt` | string | 是 | Alternate allele (e.g., 'T', 'C') |
| `version` | string | 否 | CADD version and genome build 可选值：`GRCh38-v1.7`, `GRCh37-v1.7`, `GRCh38-v1.6`, `GRCh37-v1.6` 默认 `GRCh38-v1.7` |
| `include_annotations` | boolean | 否 | Include full annotation details in response 默认 `False` |

??? note "上游说明（英文原文）"

    Get CADD deleteriousness score for a specific variant. PHRED scores: &gt;=20 top 1% deleterious (pathogenic threshold ~15-20), &gt;=30 top 0.1%. Input: chromosome, position, ref, alt alleles. Supports GRCh37 and GRCh38.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CADD_get_variant_score',
     'arguments': {'chrom': '...', 'pos': 0, 'ref': '...', 'alt': '...'}},
)
```

## `CancerVar_classify_variant` {#CancerVar_classify_variant}

Classifies a somatic variant into AMP/ASCO/CAP Tier I-IV using CancerVar.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`chrom`, `pos`, `ref`, `alt`
- **关键词：**somatic variant, amp asco cap, tier, cancer
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chrom` | string | 是 | Chromosome number or name. Can include 'chr' prefix (e.g., '7', 'chr7', 'X'). |
| `pos` | integer | 是 | Genomic position (1-based). |
| `ref` | string | 是 | Reference allele (e.g., 'A', 'G', 'ATCG'). |
| `alt` | string | 是 | Alternate allele (e.g., 'T', 'C'). |
| `build` | string | 否 | Genome build. Default: hg19. 可选值：`hg19`, `hg38` 默认 `hg19` |

??? note "上游说明（英文原文）"

    Classify a somatic cancer variant using AMP/ASCO/CAP 2017 guidelines (CancerVar). Returns a Tier I–IV assignment: Tier I = Strong Clinical Significance (FDA-approved therapy association), Tier II = Potential Clinical Significance (investigational/preclinical evidence), Tier III = Unknown Clinical Significance, Tier IV = Benign/Likely Benign. Also returns 12 CBP evidence scores and an Oncogenicity Pathogenicity Index (OPAI, 0–1). Input requires genomic coordinates in hg19 or hg38. No API key needed. Example: BRAF V600E (chr7:140453136 A&gt;T hg19) is Tier I (strong) with OPAI=0.99.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CancerVar_classify_variant',
     'arguments': {'chrom': '...', 'pos': 0, 'ref': '...', 'alt': '...'}},
)
```

## `cBioPortal_get_copy_number_alterations` {#cBioPortal_get_copy_number_alterations}

Returns GISTIC copy-number alteration calls per sample for genes in a cancer study.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`study_id`, `gene_list`
- **关键词：**cna, gistic, amplification, deletion
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `study_id` | string | 是 | Cancer study ID (e.g., 'brca_tcga_pan_can_atlas_2018'). The discrete GISTIC profile is resolved automatically. |
| `gene_list` | string | 是 | Comma-separated gene symbols (e.g., 'ERBB2' or 'TP53,ERBB2') |
| `gene` | string | 否 | Alias for gene_list: comma-separated gene symbols |
| `alteration_type` | string | 否 | CNA category filter: 'AMP' (amplification), 'GAIN', 'DIPLOID', 'HETLOSS' (shallow loss), 'HOMDEL' (deep deletion), or 'ALL' for every altered sample (default 'ALL'). 可选值：`AMP`, `GAIN`, `DIPLOID`, `HETLOSS`, `HOMDEL`, `ALL` 默认 `ALL` |
| `sample_list_id` | string | 否 | Optional sample list ID. Defaults to '\{study_id\}_all'. |
| `molecular_profile_id` | string | 否 | Optional explicit GISTIC molecular profile ID (overrides auto-resolution). |

??? note "上游说明（英文原文）"

    Get discrete copy-number alteration (CNA) calls per sample from GISTIC profiles for a gene in a cancer study. Returns per-sample alteration values: -2 (deep deletion), -1 (shallow loss), 0 (neutral), 1 (gain), 2 (amplification), plus a count breakdown by category. Use alteration_type to filter to a single category (AMP, GAIN, DIPLOID, HETLOSS, HOMDEL) or ALL for every altered sample. Example: ERBB2 (entrez 2064) AMP in brca_tcga_pan_can_atlas_2018 returns 123 amplified samples.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'cBioPortal_get_copy_number_alterations',
     'arguments': {'study_id': '...', 'gene_list': '...'}},
)
```

## `cBioPortal_get_mutations` {#cBioPortal_get_mutations}

Retrieves mutation data for given genes in a cBioPortal cancer study.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`study_id`, `gene_list`
- **关键词：**somatic mutations, tumor, cbioportal, maf
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `study_id` | string | 是 | Cancer study ID (e.g., 'brca_tcga') |
| `gene_list` | string | 是 | Comma-separated gene symbols (e.g., 'BRCA1,BRCA2') |
| `sample_list_id` | string | 否 | Optional sample list ID. If not provided, uses all samples in the study. |

??? note "上游说明（英文原文）"

    Get mutation data for specific genes in a cancer study. This uses the updated cBioPortal API that requires molecular profiles.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'cBioPortal_get_mutations',
     'arguments': {'study_id': '...', 'gene_list': '...'}},
)
```

## `civic_get_assertion` {#civic_get_assertion}

Retrieves a CIViC assertion integrating evidence on a variant's clinical actionability.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`assertion_id`
- **关键词：**civic, assertion, clinical actionability, cancer variant
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `assertion_id` | integer | 是 | CIViC assertion ID (e.g., 101) |

??? note "上游说明（英文原文）"

    Get detailed information about a specific assertion in CIViC database by assertion ID. Assertions are higher-level clinical interpretations that integrate multiple evidence items into formal statements about clinical actionability. Returns assertion description, evidence summary, and associated molecular profile or variant information.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'civic_get_assertion', 'arguments': {'assertion_id': 0}},
)
```

## `civic_get_evidence_item` {#civic_get_evidence_item}

Retrieves a CIViC evidence item linking a variant or molecular profile to clinical outcomes.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`evidence_id`
- **关键词：**civic, evidence item, biomarker, therapy response
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `evidence_id` | integer | 是 | CIViC evidence item ID (e.g., 116) |

??? note "上游说明（英文原文）"

    Get detailed information about a specific evidence item in CIViC database by evidence ID. Evidence items link molecular profiles or variants to clinical outcomes, therapies, or disease states with curated literature support. Returns evidence description, level (A-E), type (PREDICTIVE, DIAGNOSTIC, PROGNOSTIC), significance, direction, status, disease, therapies, molecular profile, and citation.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'civic_get_evidence_item',
     'arguments': {'evidence_id': 0}},
)
```

## `civic_get_molecular_profile` {#civic_get_molecular_profile}

Retrieves a CIViC molecular profile (variant combination) with its evidence.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`molecular_profile_id`
- **关键词：**molecular profile, biomarker, braf v600e, civic
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `molecular_profile_id` | integer | 是 | CIViC molecular profile ID (e.g., 12 for BRAF V600E) |
| `evidence_limit` | integer | 否 | Maximum number of evidence items to return (default: 10) 默认 `10` |

??? note "上游说明（英文原文）"

    Get detailed information about a specific molecular profile in CIViC database by molecular profile ID. Molecular profiles represent combinations of variants or features (e.g., 'BRAF V600E') that serve as biomarkers for clinical interpretation. Returns profile name, description, clinical significance score, and associated evidence items with diseases, therapies, and evidence levels.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'civic_get_molecular_profile',
     'arguments': {'molecular_profile_id': 0}},
)
```

## `civic_get_variant` {#civic_get_variant}

Retrieves a CIViC variant with coordinates, gene and clinical significance.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant_id`
- **关键词：**civic variant, snv, indel, coordinates
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_id` | integer | 是 | CIViC variant ID (e.g., 4170) |

??? note "上游说明（英文原文）"

    Get detailed information about a specific variant in CIViC database by variant ID. Variants represent specific genetic alterations (SNVs, indels, structural variants, etc.) with clinical significance in cancer. Returns variant name, coordinates, gene association, and related information.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'civic_get_variant', 'arguments': {'variant_id': 0}},
)
```

## `civic_get_variants_by_gene` {#civic_get_variants_by_gene}

Lists all CIViC variants for a given gene.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**gene variants, civic, cancer mutations
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_id` | integer | 否 | CIViC gene ID (e.g., 19 for EGFR, 12 for BRAF). Find gene IDs using civic_search_genes. |
| `gene_name` | string | 否 | Gene symbol (e.g., 'EGFR', 'BRAF', 'TP53'). Will be looked up automatically. Aliases: gene, gene_symbol, query. |
| `gene_symbol` | string | 否 | Alias for gene_name. Standard gene symbol (e.g., 'KRAS', 'BRCA1', 'EGFR'). |
| `gene` | string | 否 | Alias for gene_name. Gene symbol (e.g., 'EGFR', 'KRAS'). |
| `limit` | integer | 否 | Maximum number of variants to return (default: 500, uses cursor pagination to bypass CIViC's 100/page server cap) 默认 `500` |

??? note "上游说明（英文原文）"

    Get all variants associated with a specific gene in CIViC database. Returns variant information including names, coordinates, and associated evidence. Use gene_name or gene_symbol (e.g., 'EGFR', 'BRAF') for convenience, or gene_id if already known. Uses cursor-based pagination to return up to 500 variants per call.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'civic_get_variants_by_gene', 'arguments': {}},
)
```

## `civic_search_assertions` {#civic_search_assertions}

Searches CIViC assertions by therapy, disease or variant name.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**assertions, clinical interpretation, civic
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `limit` | integer | 否 | Maximum number of assertions to return (default: 20, recommended max: 100) 默认 `20` |
| `therapy` | string | 否 | Filter by therapy/drug name (e.g., 'imatinib', 'ponatinib'). Alias: therapy_name. |
| `therapy_name` | string | 否 | Alias for therapy. Filter by therapy/drug name. |
| `disease` | string | 否 | Filter by disease name (e.g., 'leukemia', 'melanoma'). Alias: disease_name. |
| `disease_name` | string | 否 | Alias for disease. Filter by disease name. |
| `variant_name` | string | 否 | Filter by variant name (e.g., 'V600E', 'T315I'). |

??? note "上游说明（英文原文）"

    Search for assertions in CIViC database. Assertions are higher-level clinical interpretations that integrate multiple evidence items into formal statements about clinical actionability. Returns a list of assertions with descriptions and associated molecular profiles or variants. Filter by therapy (e.g., 'ponatinib'), disease (e.g., 'leukemia'), or variant name.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'civic_search_assertions', 'arguments': {}},
)
```

## `civic_search_evidence_items` {#civic_search_evidence_items}

Searches CIViC evidence items by disease, therapy, molecular profile and significance.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**evidence search, predictive, prognostic, civic
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `limit` | integer | 否 | Maximum number of evidence items to return (default: 20, recommended max: 100) 默认 `20` |
| `status` | string \| null | 否 | Filter by curation status. Default: ACCEPTED (peer-reviewed). Options: ACCEPTED, SUBMITTED, REJECTED, ALL (returns all statuses combined). |
| `therapy` | string | 否 | Filter by therapy/drug name (e.g., 'imatinib', 'pembrolizumab'). Alias: therapy_name. Note: matches any evidence item where the therapy appears, including combination regimens — results may include multi-drug combinations. |
| `therapy_name` | string | 否 | Alias for therapy. Filter by therapy/drug name. |
| `disease` | string | 否 | Filter by disease name (e.g., 'leukemia', 'melanoma', 'lung cancer'). Alias: disease_name. Note: CIViC uses specific disease names (e.g., 'Lung Non-small Cell Carcinoma', not 'NSCLC'); try partial names or multiple searches if results are empty. |
| `disease_name` | string | 否 | Alias for disease. Filter by disease name. |
| `molecular_profile` | string | 否 | Filter by molecular profile name (e.g., 'BRAF V600E', 'EGFR T790M', 'KRAS G12C'). Uses substring matching — 'FLT3 ITD' will also match 'FLT3 ITD AND FLT3 D835Y'. For gene fusions, CIViC uses double-colon notation: 'GENE::PARTNER Fusion' (e.g., 'FGFR2::BICC1 Fusion', 'ALK::EML4 Fusion'). Use civic_search_molecular_profiles to discover exact profile names. Gene-prefixed clinical spellings ('EGFRvII… |
| `evidence_type` | string \| null | 否 | Filter by evidence type. Values: PREDICTIVE (drug response), DIAGNOSTIC (disease diagnosis), PROGNOSTIC (patient outcomes), PREDISPOSING (disease risk), ONCOGENIC (variant pathogenicity), FUNCTIONAL (molecular function). |
| `significance` | string \| null | 否 | Filter by clinical significance. Values depend on evidence_type, e.g. SENSITIVITYRESPONSE, RESISTANCE (predictive); BETTER_OUTCOME, POOR_OUTCOME (prognostic); POSITIVE, NEGATIVE (diagnostic); ONCOGENIC, PROTECTIVE. Use to separate, e.g., resistance from sensitivity evidence for a variant. |

??? note "上游说明（英文原文）"

    Search for evidence items in CIViC database. Evidence items are curated statements linking variants or molecular profiles to clinical outcomes, therapies, or disease states. Returns evidence with disease, therapies, significance, direction, and citation. Filter by therapy name (e.g., 'imatinib'), disease name (e.g., 'leukemia'), or molecular profile name (e.g., 'BRAF V600E').

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'civic_search_evidence_items', 'arguments': {}},
)
```

## `civic_search_molecular_profiles` {#civic_search_molecular_profiles}

Searches CIViC molecular profiles by name.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**molecular profiles, biomarker search, civic
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `limit` | integer | 否 | Maximum number of molecular profiles to return (default: 20, recommended max: 100) 默认 `20` |
| `query` | string | 否 | Filter by molecular profile name (e.g., 'BRAF V600E', 'EGFR T790M', 'FLT3 ITD'). Returns profiles whose name contains this string. Alias: 'name'. |
| `name` | string | 否 | Alias for query. Filter by molecular profile name (e.g., 'PD-L1', 'KRAS G12C'). |

??? note "上游说明（英文原文）"

    Search for molecular profiles in CIViC database. Molecular profiles represent combinations of variants or features (e.g., 'BRAF V600E', 'EGFR T790M') that serve as biomarkers for clinical interpretation. Returns a list of molecular profiles with IDs and names. Use the 'query' parameter to filter by profile name (alias: 'name').

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'civic_search_molecular_profiles', 'arguments': {}},
)
```

## `civic_search_variants` {#civic_search_variants}

Searches CIViC variants by variant name or gene.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**variant search, v600e, civic, gene
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `limit` | integer | 否 | Maximum number of variants to return (default: 20, recommended max: 100) 默认 `20` |
| `query` | string | 否 | Variant name to search for (e.g., "T790M", "V600E", "exon 19 deletion"). Returns name-matching variants. |
| `gene` | string | 否 | Gene symbol to filter variants by (e.g., 'EGFR', 'BRAF', 'KRAS'). Returns all variants for that gene. |
| `gene_name` | string | 否 | Alias for gene. Gene symbol (e.g., 'TP53', 'BRCA1'). |
| `variant_name` | string | 否 | Specific variant name to filter within a gene (e.g., 'L858R', 'V600E', 'T790M'). Use together with gene_name — CIViC stores variants without the gene prefix (use 'L858R' not 'EGFR L858R'). Gene-prefixed clinical spellings ('EGFRvIII', 'BRAFV600E') are still handled: the gene-stripped designation is queried as well and both result sets are merged, reported in normalization_note. |

??? note "上游说明（英文原文）"

    Search for variants in CIViC database. Returns a list of variants with their IDs and names. Use 'query' to search by variant name (e.g., 'V600E'), or 'gene'/'gene_name' to get all variants for a specific gene (e.g., 'EGFR'). To find a specific variant in a gene, combine gene_name + variant_name (e.g., gene_name='EGFR', variant_name='L858R'). Note: CIViC stores variants without gene prefix — do not combine gene and variant in a single string.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'civic_search_variants', 'arguments': {}},
)
```

## `ClinGen_dosage_by_gene` {#ClinGen_dosage_by_gene}

Returns ClinGen haploinsufficiency and triplosensitivity scores with coordinates for a gene.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`gene`
- **关键词：**dosage sensitivity, grch38, haploinsufficiency
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 是 | Gene symbol to search for (e.g., 'BRCA1', 'MECP2', 'TP53'). Case-insensitive. |

??? note "上游说明（英文原文）"

    Get ClinGen dosage sensitivity curation for a specific gene using the JSON API. Returns haploinsufficiency (HI) and triplosensitivity (TS) assertion scores with genomic coordinates in both GRCh37 and GRCh38. Scores: 0=No evidence, 1=Little evidence, 2=Emerging evidence, 3=Sufficient evidence, 30=Gene associated with autosomal recessive phenotype, 40=Dosage sensitivity unlikely. Also includes pLI score, HI prediction score, and associated diseases with MONDO IDs.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ClinGen_dosage_by_gene', 'arguments': {'gene': '...'}},
)
```

## `ClinGen_dosage_region_search` {#ClinGen_dosage_region_search}

Searches ClinGen dosage sensitivity curations overlapping chromosomal coordinates.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`chromosome`, `start`, `end`
- **关键词：**genomic region, cnv, dosage, overlap
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chromosome` | string | 是 | Chromosome (e.g., '17', 'X'). |
| `start` | integer | 是 | Start position. |
| `end` | integer | 是 | End position. |
| `assembly` | string | 否 | Genome assembly version. 可选值：`GRCh37`, `GRCh38` 默认 `GRCh38` |

??? note "上游说明（英文原文）"

    Search ClinGen dosage sensitivity curations by genomic region. Returns all genes and recurrent CNV regions with dosage sensitivity curations overlapping the specified chromosomal coordinates. Supports both GRCh37 and GRCh38. Results include haploinsufficiency and triplosensitivity scores. Useful for interpreting the clinical significance of deletions and duplications detected by chromosomal microarray.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ClinGen_dosage_region_search',
     'arguments': {'chromosome': '...', 'start': 0, 'end': 0}},
)
```

## `ClinGen_get_dosage_sensitivity` {#ClinGen_get_dosage_sensitivity}

Returns ClinGen haploinsufficiency and triplosensitivity curations, optionally with CNV regions.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**haploinsufficiency, triplosensitivity, cnv, dosage
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 否 | Optional: Filter by gene symbol |
| `include_regions` | boolean | 否 | Include recurrent CNV regions (default: false) 默认 `False` |

??? note "上游说明（英文原文）"

    Get all ClinGen dosage sensitivity curations. Returns haploinsufficiency and triplosensitivity scores for curated genes. Set include_regions=true to also get recurrent CNV regions. Scores: 0=No evidence, 1=Little, 2=Emerging, 3=Sufficient evidence.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ClinGen_get_dosage_sensitivity', 'arguments': {}},
)
```

## `ClinGen_get_variant_classifications` {#ClinGen_get_variant_classifications}

Returns ClinGen Evidence Repository variant classifications filtered by gene or variant.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**acmg, pathogenicity, vcep, variant classification
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 否 | Optional: Filter by gene symbol |
| `variant` | string | 否 | Optional: filter to one variant. Accepts a ClinGen canonical allele ID ('CA16020993', with or without a 'CAR:' prefix), a ClinVar VariationID ('370074'), full HGVS ('NM_000277.3:c.1315+1G&gt;T'), or a protein change ('p.Arg408Trp'). The form supplied selects the matching server-side filter. HGVS and protein-change strings match on substring, so a short fragment can return many rows -- combine with '… |

??? note "上游说明（英文原文）"

    Get variant pathogenicity classifications from ClinGen Evidence Repository. Returns expert-curated variant classifications with supporting evidence. Optional filters by gene or variant. Classifications follow ACMG/AMP guidelines with expert review.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ClinGen_get_variant_classifications', 'arguments': {}},
)
```

## `ClinGen_search_dosage_sensitivity` {#ClinGen_search_dosage_sensitivity}

Searches ClinGen dosage sensitivity curations by gene symbol.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`gene`
- **关键词：**haploinsufficiency, triplosensitivity, cnv interpretation
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 是 | Gene symbol to search (e.g., MECP2, PMP22, RAI1) |

??? note "上游说明（英文原文）"

    Search ClinGen dosage sensitivity curations by gene symbol. Returns haploinsufficiency (HI) and triplosensitivity (TS) scores (0-3 scale). Score 3 = Sufficient evidence for dosage pathogenicity. Critical for CNV interpretation. Example: MECP2 has HI score of 3.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ClinGen_search_dosage_sensitivity',
     'arguments': {'gene': '...'}},
)
```

## `ClinGenAllele_get_allele` {#ClinGenAllele_get_allele}

Retrieves a ClinGen Allele Registry canonical allele with HGVS and database cross-references.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**ca id, allele registry, hgvs, dbsnp
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ca_id` | string \| null | 否 | ClinGen canonical allele identifier (e.g., 'CA000387'). Obtain from ClinGenAllele_lookup_hgvs or ClinVar. |
| `allele_id` | string \| null | 否 | Alias for ca_id. ClinGen canonical allele identifier (e.g., 'CA000387'). |

??? note "上游说明（英文原文）"

    Get detailed information for a canonical allele by its ClinGen Allele Registry CA ID. Returns cross-references to ClinVar, dbSNP, COSMIC, gnomAD, genomic coordinates, and transcript-level HGVS expressions. COORDINATES: every genomic_alleles row reports start/end in 0-based interbase (half-open) form and carries a 'coordinate_system' field saying so; the 1-based genomic position used by HGVS, VCF and gnomAD is start+1 (e.g. start 32332711 = g.32332712). Each row also carries 'reference_sequence', the registry RefSeq URI (e.g. http://reg.genome.network/refseq/RS000061) identifying the sequence the coordinates are counted against - for RefSeqGene/LRG rows that have no reference_genome and no chromosome, this is the only field identifying the coordinate frame.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ClinGenAllele_get_allele', 'arguments': {}},
)
```

## `ClinGenAllele_lookup_hgvs` {#ClinGenAllele_lookup_hgvs}

Looks up an HGVS expression in the ClinGen Allele Registry to get its CA ID and cross-references.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`hgvs`
- **关键词：**hgvs, canonical allele, normalization, clinvar
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `hgvs` | string | 是 | HGVS expression for the variant. Examples: 'NC_000017.11:g.7674220C&gt;T' (genomic), 'NM_000546.6:c.743G&gt;A' (coding). Must use a valid reference sequence accession. |

??? note "上游说明（英文原文）"

    Look up a genetic variant by HGVS expression in the ClinGen Allele Registry to get its canonical allele identifier (CA ID) and cross-references to ClinVar, dbSNP, COSMIC, gnomAD, and other databases. Accepts genomic (NC_*), coding (NM_*), or protein (NP_*) HGVS expressions. COORDINATES: every genomic_alleles row reports start/end in 0-based interbase (half-open) form and carries a 'coordinate_system' field saying so; the 1-based genomic position used by HGVS, VCF and gnomAD is start+1 (e.g. start 32332711 = g.32332712). Each row also carries 'reference_sequence', the registry RefSeq URI (e.g. http://reg.genome.network/refseq/RS000061) identifying the sequence the coordinates are counted against - for RefSeqGene/LRG rows that have no reference_genome and no chromosome, this is the only field identifying the coordinate frame.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ClinGenAllele_lookup_hgvs', 'arguments': {'hgvs': '...'}},
)
```

## `ClinGenAR_get_external_records` {#ClinGenAR_get_external_records}

Returns external database records (ClinVar, dbSNP, gnomAD, etc.) for a ClinGen allele.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`allele_id`
- **关键词：**cross-references, ca id, gnomad, cosmic
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `allele_id` | string | 是 | ClinGen canonical allele ID. Examples: 'CA000387' (TP53 R248Q), 'CA006116' (BRCA1 variant). Use ClinGenAR_lookup_allele first to get the CA ID. |

??? note "上游说明（英文原文）"

    Get detailed cross-references and external database records for a ClinGen canonical allele. Takes a ClinGen allele ID (CA ID) and returns links to ClinVar, dbSNP, COSMIC, gnomAD, ExAC, and other databases. Useful for finding all known database entries for a given variant. Example: CA000387 (TP53 R248Q) returns ClinVar variation IDs, COSMIC IDs, dbSNP rsIDs, and gnomAD allele frequencies.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ClinGenAR_get_external_records',
     'arguments': {'allele_id': '...'}},
)
```

## `ClinGenAR_lookup_allele` {#ClinGenAR_lookup_allele}

Looks up a variant by HGVS in the ClinGen Allele Registry and returns its CA ID.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`hgvs`
- **关键词：**hgvs, allele registry, canonical allele
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `hgvs` | string | 是 | HGVS expression for the variant. Supports coding (c.), genomic (g.), and protein (p.) notations. Examples: 'NM_000546.6:c.743G&gt;A' (TP53 R248Q), 'NC_000017.11:g.7674220C&gt;T', 'NM_000059.4:c.5946delT' (BRCA2). |

??? note "上游说明（英文原文）"

    Look up a genetic variant in the ClinGen Allele Registry by HGVS notation. Returns the canonical allele identifier (CA ID), community standard title, and cross-references to ClinVar, dbSNP, COSMIC, gnomAD, and other databases. The Allele Registry normalizes variant representations across different coordinate systems and nomenclatures. Example: 'NM_000546.6:c.743G&gt;A' returns CA000387 (TP53 p.Arg248Gln) with links to ClinVar, COSMIC, and gnomAD.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ClinGenAR_lookup_allele', 'arguments': {'hgvs': '...'}},
)
```

## `ClinGenAR_lookup_by_external_id` {#ClinGenAR_lookup_by_external_id}

Resolves a dbSNP rsID or ClinVar VariationID to a ClinGen canonical allele with cross-references.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**rsid, clinvar variation id, reverse lookup
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `dbsnp_rs` | string \| integer \| null | 否 | dbSNP rsID to resolve. The 'rs' prefix is optional (e.g. 113488022 or 'rs113488022'). Mutually exclusive with clinvar_variation_id. |
| `clinvar_variation_id` | string \| integer \| null | 否 | ClinVar VariationID to resolve (e.g. 13961). Mutually exclusive with dbsnp_rs. |

??? note "上游说明（英文原文）"

    Reverse-lookup the ClinGen Allele Registry by an external identifier. Given a dbSNP rsID (dbsnp_rs) or a ClinVar VariationID (clinvar_variation_id), resolves it to the canonical ClinGen allele (CA id) and returns ALL cross-references (ClinVar alleles/variations, dbSNP, COSMIC, gnomAD, MyVariantInfo) plus the community standard HGVS title and genomic allele HGVS. This complements ClinGenAR_lookup_allele (which takes an HGVS expression): use this when you start from an rsID or ClinVar variation id instead of HGVS. Example: dbsnp_rs=113488022 -&gt; CA123643, 'NM_004333.6(BRAF):c.1799T&gt;A (p.Val600Glu)', incl COSM476; clinvar_variation_id=13961 resolves to the same CA123643 record. Public API, no key required.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ClinGenAR_lookup_by_external_id', 'arguments': {}},
)
```

## `ClinVar_get_clinical_significance` {#ClinVar_get_clinical_significance}

Retrieves pathogenicity classification and clinical interpretations for a ClinVar variant.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant_id`
- **关键词：**clinvar, pathogenicity, clinical significance
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_id` | string | 是 | ClinVar variant ID (e.g., '12345', '123456') |

??? note "上游说明（英文原文）"

    Get clinical significance information for a variant from ClinVar. Returns pathogenicity classification and clinical interpretations.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ClinVar_get_clinical_significance',
     'arguments': {'variant_id': '...'}},
)
```

## `ClinVar_get_submitted_records` {#ClinVar_get_submitted_records}

Retrieves individual submitter SCV assertions for a ClinVar variant.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant_id`
- **关键词：**scv, submitter, clinical assertion, clinvar
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_id` | string | 是 | ClinVar variant identifier. Accepts a VCV accession (e.g. 'VCV000013961') OR a bare ClinVar variation id (e.g. '13961'); a numeric id is normalized to VCV%09d. Example: 'VCV000013961' (BRAF V600E). |

??? note "上游说明（英文原文）"

    Retrieve the individual per-submitter assertions (SCV / ClinicalAssertion records) for a ClinVar variant. Unlike ClinVar_search_variants / ClinVar_get_variant_details / ClinVar_get_clinical_significance (which use esummary and return only the single AGGREGATE germline classification), this tool calls NCBI eutils efetch (rettype=vcv) and parses the full VariationArchive XML to expose every submitter's own classification, review status, condition and last-evaluated date. Use this when you need to see who classified a variant and how (e.g. to inspect conflicting interpretations behind a 'Conflicting classifications of pathogenicity' aggregate). Keyless NCBI endpoint.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ClinVar_get_submitted_records',
     'arguments': {'variant_id': '...'}},
)
```

## `ClinVar_get_variant_details` {#ClinVar_get_variant_details}

Retrieves ClinVar variant summary including genes, location and review status.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant_id`
- **关键词：**clinvar, variant summary, review status
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_id` | string | 是 | ClinVar variant ID (e.g., '12345', '123456') |

??? note "上游说明（英文原文）"

    Get variant summary information from ClinVar by variant ID. Returns accession, title, genes, clinical significance, review status, chromosome location, and variation name. For the same variant, this provides a different field view than ClinVar_get_clinical_significance.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ClinVar_get_variant_details',
     'arguments': {'variant_id': '...'}},
)
```

## `ClinVar_search_by_region` {#ClinVar_search_by_region}

Finds ClinVar variants overlapping a genomic region with optional significance filter.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**genomic region, overlap, cnv, clinvar
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `region` | string | 否 | Region as written, e.g. 'chr7:155593770-155593780'. Alternative to chrom/start/end. |
| `chrom` | string | 否 | Chromosome, e.g. '7' or 'chr7'. |
| `start` | integer | 否 | Region start (1-based). |
| `end` | integer | 否 | Region end (1-based, inclusive). |
| `assembly` | string | 否 | Assembly the coordinates are in. Default GRCh37. 可选值：`GRCh37`, `GRCh38` 默认 `GRCh37` |
| `clinical_significance` | string | 否 | Optional filter, e.g. 'pathogenic'. |
| `margin` | integer | 否 | How far upstream to look for variants that start before the region and span into it. Default 2000000. Raising it costs few extra candidates because single-base records upstream are filtered out, so it is a safe way to reach larger spanning CNVs -- unless 'truncated' is already true, in which case a wider window only makes the truncation worse. 默认 `2000000` |
| `max_results` | integer | 否 | Maximum overlapping variants to return; the full number found is reported as 'n_overlapping'. Default 50. 默认 `50` |

??? note "上游说明（英文原文）"

    Find ClinVar variants that OVERLAP a genomic region, including large copy-number variants that begin far outside it. Use this for 'which variant is located at chr:start-end' questions: Entrez's chrpos search matches a variant's START position, so a narrow window misses a pathogenic CNV that spans the region but starts megabases upstream (an 11 bp search at chr7:155593770-155593780 returns nothing, though a 1.5 Mb pathogenic loss covers it). The search therefore widens the window upstream by 'margin' and keeps only records whose span genuinely overlaps; records Entrez indexes as single-base are excluded from the upstream part of that window, since a one-base variant starting outside the region cannot reach into it. Results are sorted smallest span first, so the most specific overlapping variant comes first. Completeness is reported, not assumed: 'total_available' is Entrez's own count for the search, 'n_candidates' is how many of those records were actually fetched and checked for overlap, and 'truncated' is true (with a 'truncation_note') when the second number is smaller than the first -- in that case an empty or short 'variants' list is NOT evidence that nothing overlaps the region.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ClinVar_search_by_region', 'arguments': {}},
)
```

## `ClinVar_search_variants` {#ClinVar_search_variants}

Searches ClinVar variants by gene, condition, ID or variant name.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**clinvar search, hgvs, gene, condition
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 否 | Gene name or symbol (e.g., 'BRCA1', 'BRCA2') At least one of gene, condition, or variant_id must be provided. |
| `condition` | string | 否 | Disease or condition name (e.g., 'breast cancer', 'diabetes') At least one of gene, condition, or variant_id must be provided. Matched against ClinVar's disease index, which carries MedGen's own condition names exactly -- a synonym, a locus name or a misspelling ('Lebers congenital amaurosis' for the indexed 'Leber congenital amaurosis') matches nothing there, so the search falls back to an unres… |
| `variant_id` | string | 否 | ClinVar variant ID (e.g., '12345') At least one of gene, condition, or variant_id must be provided. |
| `variant_name` | string \| array \| null | 否 | Protein change or HGVS notation to search for (e.g., 'Glu6Val', 'V600E', 'c.20A&gt;T'), matched against ClinVar's own \[Variant name\] index -- combine with 'gene' for an exact-match lookup of a specific known variant instead of browsing all of a gene's ClinVar entries. A gene-level browse can silently miss well-known but old/low-ID records once a gene has more entries than max_results (e.g. HBB's can… |
| `max_results` | integer | 否 | Maximum number of results to return (default: 20). Alias: limit. 默认 `20` |
| `limit` | integer | 否 | Alias for max_results: maximum number of results to return. |
| `clinical_significance` | string | 否 | Filter by clinical significance. Must be one of the classes ClinVar's index actually recognizes: 'Pathogenic', 'Likely pathogenic', 'Uncertain significance', 'VUS', 'Likely benign', 'Benign', 'Established risk allele', 'Likely risk allele', 'Uncertain risk allele', 'drug response', 'risk factor', 'not provided', 'other'. Case-insensitive; underscores/hyphens are treated as spaces. Combine two cla… |
| `gene_symbol` | string \| null | 否 | Alias for gene. HGNC gene symbol (e.g., "DPYD", "CYP2C19"). |
| `significance` | string \| null | 否 | Alias for clinical_significance (e.g., "pathogenic", "benign", "uncertain_significance"). Same accepted value set and same validation as clinical_significance. |
| `query` | string \| null | 否 | Alias for condition. Free-text search mapped to condition/disease field. |

??? note "上游说明（英文原文）"

    Search for variants in ClinVar database by gene name, condition, variant ID, or variant name (protein change / HGVS notation). Returns variant identifiers and basic information. At least one of gene, condition, variant_id, or variant_name must be provided. Combine gene + variant_name for an exact-match lookup of one specific known variant instead of browsing all of a gene's ClinVar entries.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ClinVar_search_variants', 'arguments': {}},
)
```

## `coding_variant_fraction` {#coding_variant_fraction}

Computes the fraction of coding variants in an Excel/CSV table matching a Sequence Ontology consequence and below a VAF threshold.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`file`
- **关键词：**vaf, variant allele frequency, sequence ontology, missense, coding variants
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `file` | string | 是 | Path to variant Excel (.xlsx) or CSV file. |
| `vaf_threshold` | number | 否 | Keep variants with VAF strictly less than this value. Default 0.3. |
| `annotation` | string | 否 | Sequence Ontology term to count in numerator. Default 'synonymous_variant'. |
| `header_rows` | integer | 否 | 1 for flat headers, 2 for MultiIndex (common in Excel variant reports). Default 1. 可选值：`1`, `2` |
| `vaf_col` | string \| null | 否 | Explicit VAF column name. If empty, auto-detect from 'variant allele freq' or 'vaf'. |
| `so_col` | string \| null | 否 | Explicit Sequence Ontology column name. If empty, auto-detect from 'sequence ontology'. |

??? note "上游说明（英文原文）"

    Compute the fraction of CODING variants in an Excel/CSV variant table that match a given Sequence Ontology annotation and fall below a VAF threshold. Uses the CODING-only denominator convention (synonymous_variant, missense_variant, splice_region_variant, stop_gained, stop_lost, start_lost, frameshift_variant, inframe_insertion, inframe_deletion) - intronic, UTR, intergenic, and up/downstream records are excluded. Auto-detects 'variant allele frequency'/'VAF' and 'Sequence Ontology' columns if not specified. Supports flat headers (header_rows=1) and MultiIndex headers (header_rows=2) common in cancer variant reports. Returns the fraction with coding denominator plus the 'naive' fraction (over all VAF-filtered records) for comparison.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'coding_variant_fraction', 'arguments': {'file': '...'}},
)
```

## `COSMIC_get_mutations_by_gene` {#COSMIC_get_mutations_by_gene}

Gets COSMIC somatic mutations for a gene, with CDS and amino-acid changes, cancer types and genomic positions.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**somatic mutation, cosmic, cancer, cds change, mutational landscape
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type (fixed: get_by_gene) 默认 `get_by_gene` |
| `gene` | string | 否 | Gene symbol (e.g., BRAF, TP53, EGFR, KRAS, PIK3CA). Alias: gene_name also accepted. |
| `gene_name` | string | 否 | Alias for gene parameter. Gene symbol (e.g., FLT3, BRAF, TP53). |
| `max_results` | integer | 否 | Maximum number of mutations to return (default: 100, max: 500) 默认 `100` |
| `genome_build` | integer | 否 | Genome build version: 37 (GRCh37/hg19) or 38 (GRCh38/hg38). Default: 37 可选值：`37`, `38` 默认 `37` |

??? note "上游说明（英文原文）"

    Get all somatic mutations for a specific gene from COSMIC database. Returns comprehensive mutation data for a gene symbol including CDS changes, amino acid changes, cancer types, and genomic positions. Useful for exploring the mutational landscape of cancer genes like BRAF, TP53, EGFR, KRAS, etc.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'COSMIC_get_mutations_by_gene', 'arguments': {}},
)
```

## `COSMIC_search_mutations` {#COSMIC_search_mutations}

Searches COSMIC somatic mutations by gene, mutation ID or description, returning IDs and CDS/AA changes.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**somatic mutation, cosmic, cancer, mutation id
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type (fixed: search) 默认 `search` |
| `terms` | string | 否 | Search query - gene name (e.g., BRAF), mutation (e.g., V600E), or mutation ID (e.g., COSM476). Aliases: query or gene also accepted. |
| `query` | string | 否 | Alias for terms. Search query - gene name, mutation, or COSMIC ID. |
| `max_results` | integer | 否 | Maximum number of results to return (default: 20, max: 500) 默认 `20` |
| `genome_build` | integer | 否 | Genome build version: 37 (GRCh37/hg19) or 38 (GRCh38/hg38). Default: 37 可选值：`37`, `38` 默认 `37` |

??? note "上游说明（英文原文）"

    Search COSMIC database for somatic mutations in cancer. COSMIC is the world's largest expert-curated database of somatic mutations in human cancer. Search by gene name, mutation ID, or mutation description. Returns mutation IDs, gene names, CDS changes, and amino acid changes. No authentication required.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'COSMIC_search_mutations', 'arguments': {}},
)
```

## `CPIC_get_alleles` {#CPIC_get_alleles}

Gets CPIC alleles for a pharmacogene with activity values and functional status.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`genesymbol`
- **关键词：**pharmacogenomics, star allele, cpic, cyp2d6, allele function
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `genesymbol` | string | 是 | Gene symbol (e.g., 'CYP2D6', 'CYP2C19', 'TPMT', 'DPYD') |
| `limit` | integer \| null | 否 | Maximum number of alleles to return per call (default 50, max 1000). Compare with 'total_count' in the response to see how many exist in total. 默认 `50` |
| `offset` | integer \| null | 否 | Zero-based index of the first allele to return, for paging through genes with more alleles than one page (e.g. offset=50 returns alleles 51-100). Default: 0. 默认 `0` |

??? note "上游说明（英文原文）"

    Get pharmacogenomics alleles and their functional status for a gene from CPIC. Returns allele names, activity values, and functional status, plus 'total_count' -- the number of alleles CPIC curates for the gene, which is often far larger than one page (CYP2D6 has 208). Returns 50 per call by default; use 'limit' (max 1000) and 'offset' to retrieve the rest. Note: use 'clinicalfunctionalstatus' field for the active clinical classification (e.g., 'No function', 'Normal function'); the 'functionalstatus' field may be null for some genes. Critical for translating genotype to phenotype (e.g., CYP2D6 *4 = no function, CYP2D6 *1 = normal function). Each allele also carries 'frequency', CPIC's per-ancestry star-allele frequencies: a map keyed by CPIC's own biogeographic-group labels ('European', 'East Asian', 'Sub-Saharan African', 'Central/South Asian', 'African American/Afro-Caribbean', 'Latino', 'American', 'Near Eastern', 'Oceanian'), passed through verbatim and deliberately NOT remapped to gnomAD population codes. CPIC does not publish a frequency for every allele -- 'frequency' is null on 13 of CYP2C19's 49 alleles and 25 of CYP2D6's 208 -- and null means 'not published', never zero: a published 0.0 for a group is real data and appears as 0.0, and single groups inside the map can be null while others carry values. Alongside it, 'strength' is CPIC's evidence level for the function assignment, 'inferredfrequency' is true when CPIC inferred the frequency rather than measuring it, and 'citations' lists supporting PMIDs; these are null on the alleles that lack a frequency. CPIC's narrative 'findings' text is not returned (it restates 'strength' at ~78 KB per gene).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CPIC_get_alleles', 'arguments': {'genesymbol': '...'}},
)
```

## `CPIC_get_recommendations` {#CPIC_get_recommendations}

Gets CPIC guideline recommendations linking pharmacogene genotype or phenotype to prescribing actions.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**pharmacogenomics, dosing guideline, genotype, metabolizer phenotype, cpic
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `guideline_id` | integer | 否 | CPIC guideline numeric ID. Alternative to drug/drug_name. Use CPIC_list_guidelines to browse IDs. |
| `drug` | string | 否 | Drug name to auto-resolve guideline_id (e.g., 'codeine', 'abacavir', 'tamoxifen'). |
| `drug_name` | string | 否 | Alias for drug. |
| `gene` | string \| null | 否 | Filter to recommendation rows for this pharmacogene (e.g. 'CYP2C19', 'CYP2D6'), matched case-insensitively against each row's own phenotypes dict. Combine with phenotype for an exact single-row lookup. |
| `phenotype` | string \| null | 否 | Filter to recommendation rows whose phenotype exactly matches this CPIC term (case-insensitive), e.g. 'Poor Metabolizer', 'Likely Poor Metabolizer', 'Normal Metabolizer', 'Ultrarapid Metabolizer'. Note 'Poor Metabolizer' and 'Likely Poor Metabolizer' are distinct CPIC categories, not synonyms -- use the exact term. If gene is also given, matches only that gene's phenotype value; otherwise matches… |
| `limit` | integer \| null | 否 | Maximum number of recommendations to return (default 50). Applied after gene/phenotype filtering, if any. |
| `offset` | integer \| null | 否 | Number of recommendations to skip for pagination (default 0). Applied after gene/phenotype filtering, if any. |

??? note "上游说明（英文原文）"

    Get drug dosing recommendations from a CPIC pharmacogenomic guideline. Returns clinically actionable recommendations linking genotype/phenotype to prescribing actions (classification strength, recommendation text, allele status). Accepts drug name (e.g. 'codeine', 'abacavir') for auto-resolution to guideline_id, or provide guideline_id directly. Without gene/phenotype, returns ALL recommendation rows for every genotype/phenotype combination the guideline covers -- pass gene (e.g. 'CYP2C19') and/or phenotype (e.g. 'Poor Metabolizer', exact CPIC term) to filter to one specific recommendation. Key guideline IDs: 100421 (HLA-B/abacavir), 100416 (CYP2D6,OPRM1,COMT/opioids-codeine), 100414 (CYP2D6,CYP2C19/tricyclic-antidepressants), 100415 (CYP2D6/tamoxifen), 100412 (CYP2C9,HLA-B/phenytoin). Note: warfarin (guideline 100425) uses a dosing calculator and returns 0 rows from this endpoint.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CPIC_get_recommendations', 'arguments': {}},
)
```

## `dbsnp_get_frequencies` {#dbsnp_get_frequencies}

Gets population-specific allele frequencies for an rsID from dbSNP.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`rsid`
- **关键词：**dbsnp, allele frequency, rsid, population
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `rsid` | string | 是 | dbSNP rsID (e.g., 'rs12345', '12345') |

??? note "上游说明（英文原文）"

    Get allele frequencies for a variant from dbSNP. Returns population-specific allele frequency data. Each study reports an 'allele_count' -- the number of ALLELES (chromosomes) observed carrying the allele, NOT a number of individuals. Do not read it as study size: for a diploid cohort the number of people is at most half of it, and the denominator is the cohort's total allele count, which dbSNP does not report here.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'dbsnp_get_frequencies', 'arguments': {'rsid': '...'}},
)
```

## `dbsnp_get_variant_by_rsid` {#dbsnp_get_variant_by_rsid}

Gets dbSNP variant information by rsID with GRCh38/GRCh37 coordinates and alleles.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`rsid`
- **关键词：**dbsnp, rsid, coordinates, alleles, grch38
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `rsid` | string | 是 | dbSNP rsID (e.g., 'rs12345', '12345') |

??? note "上游说明（英文原文）"

    Get variant information from dbSNP by rsID. Returns genomic coordinates, alleles, and basic variant information. Coordinates are reported on BOTH assemblies: 'position'/'position_grch38' are on GRCh38 (dbSNP's current assembly, named by the 'assembly' field) and 'position_grch37' is the same variant on GRCh37. Always pick the field matching your pipeline's build -- the two can differ by hundreds of kilobases.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'dbsnp_get_variant_by_rsid', 'arguments': {'rsid': '...'}},
)
```

## `dbsnp_search_by_gene` {#dbsnp_search_by_gene}

Searches dbSNP for variants in a gene by symbol.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`gene_symbol`
- **关键词：**dbsnp, gene variants, snp, rsid
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 是 | Gene symbol (e.g., 'BRCA1', 'TP53', 'APOE') |
| `limit` | integer | 否 | Maximum number of results to return 默认 `20` |

??? note "上游说明（英文原文）"

    Search for variants in a specific gene. Returns variants associated with the gene symbol.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'dbsnp_search_by_gene',
     'arguments': {'gene_symbol': '...'}},
)
```

## `EBIProteins_get_variation` {#EBIProteins_get_variation}

Gets known sequence variants for a UniProt protein from COSMIC, ClinVar, gnomAD/ExAC and UniProt.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`accession`
- **关键词：**clinvar, cosmic, gnomad, missense, variants
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `accession` | string | 是 | UniProt accession. Examples: 'P04637' (TP53), 'P00533' (EGFR), 'P38398' (BRCA1). |
| `source_type` | string \| null | 否 | Filter by variant source. Options: 'large_scale_study' (COSMIC, gnomAD), 'mixed' (multiple sources), 'uniprot' (curated). Null returns all. |
| `disease_only` | boolean | 否 | If true, return only variants with disease associations. Default: false. |

??? note "上游说明（英文原文）"

    Get all known protein sequence variants for a UniProt protein from the EBI Proteins API. Returns variants from multiple sources including COSMIC (somatic mutations in cancer), ClinVar (clinical significance), gnomAD/ExAC (population frequencies), and UniProt (curated). Includes disease associations, clinical significance classifications, and cross-references. Essential for variant impact analysis and clinical interpretation. Example: P04637 (TP53) returns 3469 variants with cancer associations and clinical classifications.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EBIProteins_get_variation',
     'arguments': {'accession': '...'}},
)
```

## `EBIProteins_get_variation_by_dbsnp` {#EBIProteins_get_variation_by_dbsnp}

Returns UniProt proteins and residue changes carrying a given dbSNP rsID.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`dbsnp_id`
- **关键词：**rsid, dbsnp, missense, protein consequence
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `dbsnp_id` | string | 是 | dbSNP reference SNP ID (rsID). Examples: 'rs1042522' (TP53 Pro72Arg), 'rs429358' (APOE). |

??? note "上游说明（英文原文）"

    Reverse variation lookup by dbSNP rsID from the EBI Proteins API. Given a dbSNP rsID, returns every UniProt protein (including isoforms) that carries that variant, each with the mapped natural-variant features (wild-type and alternative residue, protein position, genomic location, clinical significances, cross-references to ClinVar/COSMIC/dbSNP). Inverse of looking up variants by protein accession. Example: 'rs1042522' (TP53 codon 72 Pro/Arg) returns 11 UniProt entries including TP53 isoforms and A0AAQ5BHX5, each with the position-72 variant feature.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EBIProteins_get_variation_by_dbsnp',
     'arguments': {'dbsnp_id': '...'}},
)
```

## `EBIProteins_get_variation_by_hgvs` {#EBIProteins_get_variation_by_hgvs}

Returns overlapping UniProt proteins and protein-level consequences for an HGVS genomic expression.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`hgvs`
- **关键词：**hgvs, variant consequence, genomic variant
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `hgvs` | string | 是 | HGVS genomic expression using a RefSeq NC_ accession. Example: 'NC_000017.11:g.7676154G&gt;A' (TP53 region). |

??? note "上游说明（英文原文）"

    Reverse variation lookup by HGVS genomic expression from the EBI Proteins API. Given an HGVS genomic-level expression (e.g. 'NC_000017.11:g.7676154G&gt;A'), returns every overlapping UniProt protein with the resulting variant consequence projected onto the protein sequence (wild-type and alternative residue, protein position, genomic location). Use to translate a genomic DNA change into its protein consequence across all affected UniProt entries. Example: 'NC_000017.11:g.7676154G&gt;A' returns 4 UniProt entries; the TP53 (P04637) entry reports VARIANT P-&gt;L at residue 72 (the Pro72Leu protein consequence).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EBIProteins_get_variation_by_hgvs',
     'arguments': {'hgvs': '...'}},
)
```

## `ensembl_get_structural_variants` {#ensembl_get_structural_variants}

Gets known structural variants overlapping a region from Ensembl (DGVa, ClinGen).

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`species`, `region`
- **关键词：**structural variant, cnv, dgva, clingen
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `species` | string | 是 | Species name (e.g., 'human', 'homo_sapiens'). |
| `region` | string | 是 | Genomic region in format 'chr:start-end' (e.g., '17:43044295-43125370'). Maximum 5Mb. |
| `feature` | string | 否 | Feature type to retrieve. Fixed to 'structural_variation' for this tool. 默认 `structural_variation` |

??? note "上游说明（英文原文）"

    Get structural variants overlapping a genomic region from Ensembl. Returns known SVs from DGVa (Database of Genomic Variants archive), ClinGen, and other sources. Each SV includes its ID (nsv/esv accession), coordinates, and feature type. Use ensembl_get_sv_detail to get clinical significance, variant class (CNV, deletion, duplication, inversion), and supporting evidence for specific variants. Region must be &lt; 5Mb.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ensembl_get_structural_variants',
     'arguments': {'species': '...', 'region': '...'}},
)
```

## `ensembl_get_sv_detail` {#ensembl_get_sv_detail}

Gets details of an Ensembl structural variant by nsv/esv accession.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`species`, `id`
- **关键词：**structural variant, cnv, deletion, clinical significance
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `species` | string | 是 | Species name (e.g., 'human', 'homo_sapiens'). |
| `id` | string | 是 | Structural variant accession from Ensembl/DGVa (e.g., 'nsv2769779', 'esv3647486'). Use ensembl_get_structural_variants to find IDs. |

??? note "上游说明（英文原文）"

    Get detailed information for a specific structural variant from Ensembl by its accession (nsv/esv ID). Returns variant class (CNV, deletion, duplication, inversion, gain, loss), clinical significance (pathogenic, benign, etc.), source database, supporting evidence with individual sample observations, study information, and genomic mappings. Essential for CNV clinical interpretation.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ensembl_get_sv_detail',
     'arguments': {'species': '...', 'id': '...'}},
)
```

## `Ensembl_get_transcript_haplotypes` {#Ensembl_get_transcript_haplotypes}

Gets 1000 Genomes protein and CDS haplotypes for an Ensembl transcript with frequencies and variants.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`id`
- **关键词：**haplotype, 1000 genomes, transcript, population frequency
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `id` | string | 是 | Ensembl transcript stable ID. Examples: 'ENST00000288602' (BRAF), 'ENST00000269305' (TP53), 'ENST00000357654' (BRCA1). |
| `species` | string | 否 | Species name (default 'human'). Use 'human' or 'homo_sapiens'. 默认 `human` |

??? note "上游说明（英文原文）"

    Get population protein and CDS haplotypes for a transcript from the Ensembl REST API (transcript_haplotypes endpoint). Returns the real observed haplotypes across 1000 Genomes phase 3 populations, each with per-population frequencies, contributing variants, sequence diffs, an indel flag, and observed counts/frequency, plus total_population_counts and total_haplotype_count. Use this to assess how a coding transcript varies across human populations. Example: ENST00000288602 (BRAF) returns 11 protein_haplotypes and 26 cds_haplotypes; the most common protein haplotype has 1000GENOMES:phase_3 population frequencies up to ~0.99.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Ensembl_get_transcript_haplotypes',
     'arguments': {'id': '...'}},
)
```

## `ensembl_get_variants` {#ensembl_get_variants}

Gets SNPs and indels overlapping a genomic region from Ensembl.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`region`
- **关键词：**snp, indel, genomic region, consequence
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `region` | string | 是 | Genomic region in format 'chromosome:start-end' (e.g., '13:32315086-32400268' or '1:1000000-2000000') |
| `species` | string | 否 | Species name (default 'human'). Examples: 'human', 'mouse', 'rat' 默认 `human` |
| `feature` | string | 否 | Feature type to retrieve (default 'variation') 默认 `variation` |

??? note "上游说明（英文原文）"

    Get genetic variants (SNPs, indels) in a genomic region. Returns variant information including positions, alleles, and consequences. Use this to find variants overlapping a specific genomic region.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ensembl_get_variants', 'arguments': {'region': '...'}},
)
```

## `ensembl_get_variation` {#ensembl_get_variation}

Gets Ensembl details for a variant ID such as rs699, including alleles, frequencies and consequences.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**rsid, alleles, consequence, maf
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `id` | string | 否 | Variation ID (e.g., 'rs699', 'rs1421085'). Use ensembl_get_variants to find variant IDs in a genomic region. |
| `variant_id` | string | 否 | Alias for id. dbSNP rsID or Ensembl variation ID (e.g., 'rs7903146', 'rs699'). |
| `species` | string | 否 | Species name (default 'human') 默认 `human` |
| `phenotypes` | integer | 否 | Include phenotype associations (0=no, 1=yes) 默认 `0` |

??? note "上游说明（英文原文）"

    Get detailed information about a specific genetic variation by its ID (e.g., rs699). Returns variation details including alleles, frequencies, consequences, and evidence. Use ensembl_get_variants to find variant IDs in a region.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ensembl_get_variation', 'arguments': {}},
)
```

## `ensembl_vep_region` {#ensembl_vep_region}

Runs Ensembl VEP on a region and allele, returning consequences and SIFT/PolyPhen/CADD scores.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`species`, `region`, `allele`
- **关键词：**vep, sift, polyphen, cadd, consequence
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `species` | string | 是 | Species (e.g., 'human', 'homo_sapiens') |
| `region` | string | 是 | Genomic region with variant (e.g., '21:25891796-25891796') |
| `allele` | string | 是 | Alternative allele (e.g., 'T', 'C') |
| `AncestralAllele` | integer | 否 | Include ancestral allele information 默认 `0` |
| `Blosum62` | integer | 否 | Include BLOSUM62 protein conservation score 默认 `0` |
| `CADD` | integer | 否 | Include CADD scores 默认 `0` |
| `Conservation` | integer | 否 | Include conservation scores 默认 `0` |
| `DisGeNET` | integer | 否 | Include DisGeNET disease associations 默认 `0` |
| `GeneSplicer` | integer | 否 | Include GeneSplicer predictions 默认 `0` |
| `GO` | integer | 否 | Include Gene Ontology terms 默认 `0` |
| `LoF` | integer | 否 | Include Loss-of-Function predictions 默认 `0` |
| `MaxEntScan` | integer | 否 | Include MaxEntScan splice predictions 默认 `0` |
| `Phenotypes` | integer | 否 | Include phenotype annotations 默认 `0` |
| `SIFT` | string | 否 | SIFT prediction: p=prediction, s=score, b=both 可选值：`p`, `s`, `b` |
| `PolyPhen` | string | 否 | PolyPhen prediction: p=prediction, s=score, b=both 可选值：`p`, `s`, `b` |

??? note "上游说明（英文原文）"

    Variant Effect Predictor (VEP) for genomic variants. Predicts functional consequences of variants including impact on genes, transcripts, and proteins. Returns SIFT/PolyPhen scores, consequence types, and affected transcripts. Essential for variant interpretation.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ensembl_vep_region',
     'arguments': {'species': '...', 'region': '...', 'allele': '...'}},
)
```

## `EnsemblVar_get_population_frequencies` {#EnsemblVar_get_population_frequencies}

Gets gnomAD and 1000 Genomes allele frequencies for a variant via Ensembl.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant_id`
- **关键词：**allele frequency, gnomad, 1000 genomes, population
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_id` | string | 是 | dbSNP rsID of the variant. Examples: 'rs429358' (APOE), 'rs7903146' (TCF7L2), 'rs1042779'. |
| `species` | string | 否 | Species name. Default 'human'. Examples: 'human', 'homo_sapiens'. |

??? note "上游说明（英文原文）"

    Get allele frequency data for a variant (SNP/indel) across global populations from gnomAD and 1000 Genomes via Ensembl. Returns allele frequencies broken down by population (e.g., gnomADe:ALL, gnomADe:nfe, 1000GENOMES:phase_3:CEU, etc.). Essential for assessing variant rarity in different ethnic groups, interpreting clinical significance, and population genetics studies. Example: rs429358 (APOE epsilon-4 allele) shows C allele frequency of 14.8% globally (gnomAD), 22.7% in African populations, 9.8% in East Asian populations. Each population entry carries population, allele, frequency and (where the source reports it) allele_count; Ensembl publishes no sample size or allele-number denominator, so none is returned.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblVar_get_population_frequencies',
     'arguments': {'variant_id': '...'}},
)
```

## `EnsemblVar_get_variant_consequences` {#EnsemblVar_get_variant_consequences}

Gets an Ensembl variant record with consequence, clinical significance and synonyms.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant_id`
- **关键词：**variant consequence, clinical significance, synonyms
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_id` | string | 是 | dbSNP rsID of the variant. Examples: 'rs429358' (APOE), 'rs7903146' (TCF7L2), 'rs1042779'. |
| `species` | string | 否 | Species name. Default 'human'. Examples: 'human', 'homo_sapiens'. |

??? note "上游说明（英文原文）"

    Get detailed variant information including consequence type, clinical significance, synonyms, evidence attributes, and allele mappings from Ensembl. More comprehensive than VEP - returns the full variant record with all known synonyms, mapped positions, ancestral allele, and minor allele data. Example: rs429358 returns most_severe_consequence='missense_variant', 27 synonyms, mapped to chr19:44908684 (GRCh38), and is classified as a missense_variant affecting APOE.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblVar_get_variant_consequences',
     'arguments': {'variant_id': '...'}},
)
```

## `EnsemblVEP_annotate_hgvs` {#EnsemblVEP_annotate_hgvs}

Runs Ensembl VEP on an HGVS notation, returning consequences and pathogenicity scores.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`hgvs_notation`
- **关键词：**vep, hgvs, sift, polyphen
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `hgvs_notation` | string | 是 | Variant in HGVS notation. Supports protein (p.), coding DNA (c.), and genomic (g.) notation. Examples: 'BRAF:p.Val600Glu', 'ENST00000366667:c.803C&gt;T', 'NC_000007.14:g.140753336A&gt;T'. |
| `species` | string | 否 | Species name. Default: 'human'. Other options: 'mouse', 'rat', 'zebrafish', etc. 默认 `human` |

??? note "上游说明（英文原文）"

    Predict functional consequences of a genetic variant using HGVS notation via the Ensembl Variant Effect Predictor (VEP). Returns impact on genes, transcripts, and proteins including SIFT/PolyPhen pathogenicity scores, consequence types, amino acid changes, and colocated known variants. HGVS notation examples: 'BRAF:p.Val600Glu' (protein), 'ENST00000366667:c.803C&gt;T' (coding DNA). This complements the existing region-based VEP tool by accepting standard clinical variant nomenclature.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblVEP_annotate_hgvs',
     'arguments': {'hgvs_notation': '...'}},
)
```

## `EnsemblVEP_annotate_rsid` {#EnsemblVEP_annotate_rsid}

Runs Ensembl VEP on a dbSNP rsID, returning consequences and pathogenicity scores.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**vep, rsid, sift, polyphen
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_id` | string | 否 | dbSNP rs identifier. Examples: 'rs7903146' (TCF7L2, type 2 diabetes), 'rs429358' (APOE, Alzheimer's), 'rs1801133' (MTHFR C677T). Aliases: rsid, rs_id. |
| `rsid` | string | 否 | Alias for variant_id. dbSNP rs identifier (e.g., 'rs7903146'). |
| `rs_id` | string | 否 | Alias for variant_id. dbSNP rs identifier (e.g., 'rs7903146'). |
| `species` | string | 否 | Species name. Default: 'human'. 默认 `human` |

??? note "上游说明（英文原文）"

    Predict functional consequences of a genetic variant using its dbSNP rs identifier via the Ensembl Variant Effect Predictor (VEP). Returns impact on genes, transcripts, and proteins including SIFT/PolyPhen scores, consequence types, amino acid changes, and colocated variants. Useful when you have a variant's rs number from GWAS, clinical reports, or literature. Example: 'rs7903146' (TCF7L2 T2D variant), 'rs429358' (APOE Alzheimer variant).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblVEP_annotate_rsid', 'arguments': {}},
)
```

## `EnsemblVEP_variant_recoder` {#EnsemblVEP_variant_recoder}

Converts a variant ID to HGVS genomic/coding/protein and SPDI notations.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant_id`
- **关键词：**variant recoder, hgvs, spdi, rsid
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_id` | string | 是 | Variant identifier to recode. Accepts: rsID (e.g., 'rs429358'), HGVS notation (e.g., 'AGT:c.803T&gt;C'), SPDI (e.g., 'NC_000019.10:44908683:T:C'). |
| `species` | string | 否 | Species name. Default: 'human'. 默认 `human` |

??? note "上游说明（英文原文）"

    Convert genetic variant identifiers between different nomenclature formats using the Ensembl Variant Recoder. Given a variant ID (e.g., rs number), returns equivalent representations in HGVS genomic (hgvsg), HGVS coding (hgvsc), HGVS protein (hgvsp), and SPDI format. Essential for cross-referencing variants across databases that use different notation systems. Example: 'rs429358' returns HGVS coding notation for all APOE transcripts.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblVEP_variant_recoder',
     'arguments': {'variant_id': '...'}},
)
```

## `EVA_get_clustered_variant_by_rs` {#EVA_get_clustered_variant_by_rs}

Resolves an EVA clustered-variant RS accession to assembly, contig and position.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`accession`
- **关键词：**eva, rs accession, variant position, accessioning
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `accession` | string | 是 | EVA clustered-variant RS accession (numeric, no 'rs' prefix). Examples: '429358' (human APOE region, GRCh38 contig CM000681.2 start 44908684), '3091781465' (horse clustered variant) |

??? note "上游说明（英文原文）"

    Resolve an EVA clustered-variant RS accession to its canonical assembly/contig/position record using the EVA accessioning (identifiers) service. This is distinct from EVA's variant-DATA API (webservices/rest) used by the other EVA tools: the accessioning API (webservices/identifiers) maps an RS accession to its definitive coordinate record (assemblyAccession, taxonomyAccession, contig, start, variant type). The RS accession is the numeric clustered-variant id (e.g. 429358), not the 'rs' prefixed dbSNP string. Returns a list of accession records (usually one).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EVA_get_clustered_variant_by_rs',
     'arguments': {'accession': '...'}},
)
```

## `EVA_get_variants_by_gene` {#EVA_get_variants_by_gene}

Gets EVA variants for a gene with allele frequencies and consequences.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`gene`
- **关键词：**eva, gene variants, allele frequency
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 是 | Gene symbol (HGNC format, e.g., 'BRCA1', 'TP53', 'EGFR', 'BRCA2', 'APOE') |
| `species` | string | 否 | Species and assembly (e.g., 'hsapiens_grch38', 'hsapiens_grch37', 'mmusculus_grcm38') 默认 `hsapiens_grch38` |
| `limit` | integer | 否 | Maximum number of variants to return (default 20, max ~1000) 默认 `20` |

??? note "上游说明（英文原文）"

    Get genetic variants from the European Variation Archive (EVA) for a specific gene. EVA is hosted at EBI and aggregates variant data from multiple studies including ClinVar, dbSNP, and population sequencing projects. Returns variants with allele frequencies, study sources, and variant type. The gene name uses HGNC symbols (e.g., 'BRCA1', 'TP53', 'EGFR'). Species defaults to human GRCh38.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EVA_get_variants_by_gene', 'arguments': {'gene': '...'}},
)
```

## `EVA_get_variants_by_region` {#EVA_get_variants_by_region}

Gets EVA variants within a genomic region.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`region`
- **关键词：**eva, region variants, snv, indel
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `region` | string | 是 | Genomic region in 'chr:start-end' format (e.g., '17:41200000-41210000', '7:140753300-140753400', '19:44905700-44910000') |
| `species` | string | 否 | Species and assembly (e.g., 'hsapiens_grch38', 'hsapiens_grch37') 默认 `hsapiens_grch38` |
| `limit` | integer | 否 | Maximum number of variants to return 默认 `20` |

??? note "上游说明（英文原文）"

    Get genetic variants from the European Variation Archive (EVA) within a specific genomic region. Returns all known variants in the specified chromosomal coordinate range, including SNVs, indels, and structural variants with allele frequencies from multiple studies. Useful for scanning variant hotspots or checking variants in regulatory regions. The endpoint URL uses 'chromosome:start-end' format.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EVA_get_variants_by_region',
     'arguments': {'region': '...'}},
)
```

## `EVA_list_studies` {#EVA_list_studies}

Lists EVA variant studies available for a species.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**eva, variant studies, prjeb, species
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `species` | string | 否 | Species and assembly (e.g., 'hsapiens_grch38', 'hsapiens_grch37', 'mmusculus_grcm38') 默认 `hsapiens_grch38` |

??? note "上游说明（英文原文）"

    List all variant studies available in the European Variation Archive (EVA) for a given species. Returns study accessions (e.g., PRJEB84021) that have submitted variant data. Useful for discovering available datasets before querying variants. Each study typically corresponds to a population sequencing project or variant calling pipeline.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EVA_list_studies', 'arguments': {}},
)
```

## `EVE_get_gene_info` {#EVE_get_gene_info}

Checks whether EVE pathogenicity scores exist for a gene and returns its Ensembl ID.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`gene_symbol`
- **关键词：**eve, gene coverage, pathogenicity
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 是 | Gene symbol (e.g., 'TP53', 'BRCA1', 'EGFR') |

??? note "上游说明（英文原文）"

    Check if EVE scores are available for a gene and get gene information. EVE covers ~3,000 disease-related human genes. Returns Ensembl ID and link to EVE website.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EVE_get_gene_info', 'arguments': {'gene_symbol': '...'}},
)
```

## `EVE_get_variant_score` {#EVE_get_variant_score}

Gets the EVE pathogenicity score for a variant via Ensembl VEP.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**eve, pathogenicity, deep learning, missense
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant` | string | 否 | Variant in HGVS format (e.g., 'ENST00000269305.4:c.100G&gt;A', 'NM_000546.5:c.215C&gt;G') |
| `chrom` | string | 否 | Chromosome (1-22, X, Y). Use with pos, ref, alt instead of variant. |
| `pos` | integer | 否 | Genomic position (GRCh38) |
| `ref` | string | 否 | Reference allele |
| `alt` | string | 否 | Alternate allele |
| `species` | string | 否 | Species (default: human) 默认 `human` |

??? note "上游说明（英文原文）"

    Get EVE (Evolutionary Variant Effect) pathogenicity score for a variant. EVE uses unsupervised deep learning on evolutionary data. Scores: 0-1 where &gt;0.5 is likely pathogenic. Accessed via Ensembl VEP API.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EVE_get_variant_score', 'arguments': {}},
)
```

## `FAVOR_annotate_variant` {#FAVOR_annotate_variant}

Annotates a GRCh38 variant via FAVOR with frequencies, consequences and functional scores.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant`
- **关键词：**favor, variant annotation, topmed, cadd, gencode
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant` | string | 是 | GRCh38 SNV/indel as chr-pos-ref-alt. Accepts '19-44908822-C-T', 'chr19:44908822:C:T', etc. |

??? note "上游说明（英文原文）"

    Comprehensive functional annotation of a single GRCh38 variant via FAVOR (Functional Annotation of Variants Online Resource, Harvard). One call returns: allele frequencies (BRAVO/TOPMed, gnomAD by ancestry, 1000 Genomes), gene and consequence (GENCODE category, exonic effect, protein variant, HGVS), deleteriousness scores (CADD, SIFT, PolyPhen-2, AlphaMissense, MetaSVM, MutationTaster/Assessor, FATHMM-XF, Grantham, LINSIGHT, FunSeq), conservation (GERP, phyloP mammalian/primate/vertebrate, phastCons), ClinVar clinical significance/disease, and regulatory/epigenomic context (CAGE promoter/enhancer, GeneHancer, ENCODE DNase, ReMap TF overlap). Use for variant interpretation across both coding and non-coding regions in one step. Input is a GRCh38 (hg38) coordinate string; resolve rsIDs to GRCh38 coordinates first if needed. No API key required.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'FAVOR_annotate_variant', 'arguments': {'variant': '...'}},
)
```

## `FlyBase_get_gene_alleles` {#FlyBase_get_gene_alleles}

Lists alleles/variants for a Drosophila gene with categories, phenotypes, disease annotations and genomic coordinates.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`gene_id`
- **关键词：**drosophila, alleles, mutations, insertion, fbgn
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_id` | string | 是 | FlyBase gene ID with 'FB:' prefix. Examples: 'FB:FBgn0003996' (white/w), 'FB:FBgn0000490' (dpp), 'FB:FBgn0004647' (Notch/N). |
| `limit` | integer | 否 | Maximum number of alleles to return (1-100). Default: 20. |
| `page` | integer | 否 | Page number for pagination. Default: 1. |

??? note "上游说明（英文原文）"

    Get alleles and variants for a Drosophila gene from FlyBase via the Alliance of Genome Resources. Returns allele symbols, categories (insertion, point mutation, etc.), associated phenotypes, disease annotations, and variant details including genomic coordinates. Example: 'FB:FBgn0000490' (dpp) returns alleles like dpp\[MI03752\] with variant location and molecular consequence.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'FlyBase_get_gene_alleles',
     'arguments': {'gene_id': '...'}},
)
```

## `GDC_get_cnv_data` {#GDC_get_cnv_data}

Queries GDC/TCGA copy number variation segment and gene-level files by project and gene symbol.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`project_id`
- **关键词：**cnv, copy number, tcga, amplification, deletion
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `project_id` | string | 是 | GDC project (e.g., 'TCGA-BRCA') |
| `gene_symbol` | string | 否 | Optional: Gene symbol to focus analysis |
| `size` | integer | 否 | Number of results 默认 `10` |

??? note "上游说明（英文原文）"

    Query copy number variation (CNV) files from GDC/TCGA projects. Returns CNV segment and gene-level copy number files. Use for amplification/deletion analysis.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GDC_get_cnv_data', 'arguments': {'project_id': '...'}},
)
```

## `GDC_get_mutation_frequency` {#GDC_get_mutation_frequency}

Returns total pan-cancer somatic mutation (SSM) occurrence count for a gene across GDC/TCGA projects.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`gene_symbol`
- **关键词：**somatic mutation, ssm, pan-cancer, tcga, mutation frequency
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 是 | Gene symbol (e.g., 'TP53', 'KRAS', 'EGFR') |
| `gene` | string | 否 | Gene symbol alias — alternative to gene_symbol |

??? note "上游说明（英文原文）"

    Get pan-cancer mutation frequency statistics for a gene across all TCGA/GDC cancer projects. Returns total SSM occurrence count (no per-project breakdown — use GDC_get_ssm_by_gene with project_id for cancer-specific data). Answers 'How often is this gene mutated across all cancers?'

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GDC_get_mutation_frequency',
     'arguments': {'gene_symbol': '...'}},
)
```

## `GDC_get_mutation_frequency_by_project` {#GDC_get_mutation_frequency_by_project}

Returns per-project mutated case counts and total case counts for a gene, giving per-cancer mutation rates.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**somatic mutation, mutation rate, tcga, cancer type
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 否 | Gene symbol (e.g., 'KRAS', 'TP53', 'EGFR') |
| `gene` | string | 否 | Gene symbol alias — alternative to gene_symbol |
| `size` | integer | 否 | Maximum number of projects to return (default 100) 默认 `100` |

??? note "上游说明（英文原文）"

    Get per-project (per-cancer-type) somatic mutation frequency for a gene from NCI GDC. For each TCGA/GDC project, returns the numerator (mutated case count) AND denominator (total case count), enabling true per-cancer mutation rates. Answers 'In which cancer types is this gene most frequently mutated?' Example: gene KRAS returns CPTAC-3 mutated=510/total=3521, ALCHEMIST-ALCH mutated=391/total=3038.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GDC_get_mutation_frequency_by_project', 'arguments': {}},
)
```

## `GDC_get_ssm_by_gene` {#GDC_get_ssm_by_gene}

Gets somatic simple mutations for a gene across GDC/TCGA projects with type, genomic coordinates and amino acid change.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`gene_symbol`
- **关键词：**somatic mutation, ssm, tcga, amino acid change, pan-cancer
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 是 | Gene symbol (e.g., 'TP53', 'EGFR', 'BRAF', 'KRAS') |
| `project_id` | string | 否 | Optional: Filter by project (e.g., 'TCGA-BRCA', 'TCGA-LUAD') |
| `size` | integer | 否 | Number of results (1–100) 默认 `20` |

??? note "上游说明（英文原文）"

    Get somatic mutations (SSMs) for a gene across TCGA/GDC projects. Returns mutation type, genomic coordinates, and amino acid changes. Essential for pan-cancer mutation analysis.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GDC_get_ssm_by_gene',
     'arguments': {'gene_symbol': '...'}},
)
```

## `GeneBe_classify_variant` {#GeneBe_classify_variant}

Classifies a germline variant (chr, pos, ref, alt; hg38/hg19) under ACMG/AMP returning verdict, score and criteria.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`chr`, `pos`, `ref`, `alt`
- **关键词：**acmg, pathogenicity, germline, variant classification, hg38
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chr` | string | 是 | Chromosome, e.g. '7', 'X' (a 'chr' prefix is also accepted). |
| `pos` | integer | 是 | 1-based genomic position of the variant, e.g. 140753336. |
| `ref` | string | 是 | Reference allele, e.g. 'A'. |
| `alt` | string | 是 | Alternate allele, e.g. 'T'. |
| `genome` | string | 否 | Genome build: 'hg38' (default) / 'GRCh38' or 'hg19' / 'GRCh37'. 默认 `hg38` |

??? note "上游说明（英文原文）"

    Classify a germline variant under the ACMG/AMP guidelines using GeneBe (genebe.net). Given chromosome, position, ref and alt alleles on genome build hg38 (default) or hg19, returns the ACMG verdict (e.g. Pathogenic), the numeric ACMG score, and the exact triggered criteria (e.g. PS3,PM1,PM2,PM5,PP2,PP3_Moderate,PP5), plus gene symbol, transcript, HGVS c./p., predicted effect, dbSNP id, gnomAD allele frequency, ClinVar classification, and an AlphaMissense score. This is an independent second opinion complementary to InterVar_classify_variant (a different group's ACMG implementation) — cross-check both for borderline variants. Public API, no key required (a free account raises the rate limit).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GeneBe_classify_variant',
     'arguments': {'chr': '...', 'pos': 0, 'ref': '...', 'alt': '...'}},
)
```

## `GeneBe_classify_variants_batch` {#GeneBe_classify_variants_batch}

Batch-classifies up to 1000 germline variants under ACMG/AMP guidelines via GeneBe.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variants`
- **关键词：**acmg, batch, pathogenicity, germline, variant classification
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variants` | array of object | 是 | List of variants to classify (max 1000). Each item is an object with chr, pos, ref, alt. |
| `genome` | string | 否 | Genome build for ALL variants in the batch: 'hg38' (default) / 'GRCh38' or 'hg19' / 'GRCh37'. 默认 `hg38` |

??? note "上游说明（英文原文）"

    Batch-classify up to 1000 germline variants under the ACMG/AMP guidelines in a single GeneBe (genebe.net) request. Provide 'variants' as a JSON array of objects, each with chr, pos, ref, alt on genome build hg38 (default) or hg19. For every variant GeneBe returns the full ACMG output: the verdict (e.g. Pathogenic), the numeric ACMG score, the exact triggered criteria (e.g. PS3,PM1,PM2,PM5,PP2,PP3_Moderate,PP5), per-gene HGVS c./p. (acmg_by_gene), AlphaMissense score/prediction, gnomAD allele frequency, and ClinVar classification. Use this instead of calling GeneBe_classify_variant once per variant when you have many variants. Public API, no key required (a free account raises the rate limit). Example: \[\{"chr":"7","pos":140753336,"ref":"A","alt":"T"\},\{"chr":"22","pos":28695868,"ref":"AG","alt":"A"\}\] returns BRAF p.Val600Glu Pathogenic score 14 plus the second variant.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GeneBe_classify_variants_batch',
     'arguments': {'variants': []}},
)
```

## `GenomeNexus_annotate_dbsnp` {#GenomeNexus_annotate_dbsnp}

Annotates a dbSNP rsID via Genome Nexus with VEP, SIFT, PolyPhen-2, AlphaMissense and cancer hotspot data.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`rsid`
- **关键词：**rsid, dbsnp, vep, alphamissense, polyphen, sift
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `rsid` | string | 是 | dbSNP reference SNP ID. Examples: 'rs121913529' (KRAS G12), 'rs113488022' (BRAF V600E), 'rs28934578' (TP53 R273H). A bare numeric id (e.g. '121913529') is accepted and prefixed with 'rs'. |

??? note "上游说明（英文原文）"

    Annotate a variant directly by its dbSNP rsID using Genome Nexus (Memorial Sloan Kettering). Genome Nexus resolves the rsID to genomic coordinates and returns aggregated annotation from VEP, SIFT, PolyPhen-2, AlphaMissense, and cancer hotspots, plus an annotation summary with the genomic location (chromosome/start/end) and HGVS protein notation. Use this when you have a dbSNP rsID rather than HGVS coordinates. Example: 'rs121913529' (KRAS G12 codon) returns 12 transcript consequences, gene KRAS, SIFT 'deleterious', PolyPhen 'probably_damaging', AlphaMissense score 0.9949 (pathogenic), genomic location chr12:25398284.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GenomeNexus_annotate_dbsnp',
     'arguments': {'rsid': '...'}},
)
```

## `GenomeNexus_annotate_mutation` {#GenomeNexus_annotate_mutation}

Annotates a mutation given chromosome, start, end, ref and alt alleles (GRCh37) via Genome Nexus.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`chromosome`, `start`, `end`, `reference_allele`, `variant_allele`
- **关键词：**mutation annotation, grch37, vep, cancer, coordinates
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chromosome` | string | 是 | Chromosome (without 'chr' prefix). Examples: '7', '17', '12', 'X'. |
| `start` | integer | 是 | Start position (GRCh37/hg19, 1-based). |
| `end` | integer | 是 | End position (GRCh37/hg19, 1-based). Same as start for SNVs. |
| `reference_allele` | string | 是 | Reference allele. Example: 'A'. |
| `variant_allele` | string | 是 | Variant (alternate) allele. Example: 'T'. |

??? note "上游说明（英文原文）"

    Annotate a cancer mutation using Genome Nexus (MSK) by specifying chromosome, start, end, reference allele, and variant allele separately. This is an alternative to the HGVS-based annotation for when you have genomic coordinates directly (GRCh37/hg19). Returns the same comprehensive annotation as GenomeNexus_annotate_variant including VEP consequences, SIFT, PolyPhen, AlphaMissense scores, cancer hotspots, and mutation assessor. Example: chr=7, start=140453136, end=140453136, ref=A, alt=T for BRAF V600E.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GenomeNexus_annotate_mutation',
     'arguments': {'chromosome': '...',
                   'start': 0,
                   'end': 0,
                   'reference_allele': '...',
                   'variant_allele': '...'}},
)
```

## `GenomeNexus_annotate_variant` {#GenomeNexus_annotate_variant}

Annotates a GRCh37 HGVSg variant with aggregated VEP, SIFT, PolyPhen-2, AlphaMissense and hotspot annotations.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`hgvsg`
- **关键词：**hgvs, hgvsg, vep, alphamissense, variant effect
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `hgvsg` | string | 是 | Genomic variant in HGVS notation using GRCh37/hg19 coordinates. Format: 'chr:g.posRef&gt;Alt'. Examples: '7:g.140453136A&gt;T' (BRAF V600E), '17:g.7577120C&gt;T' (TP53 R273H), '12:g.25398284C&gt;T' (KRAS G12D). |

??? note "上游说明（英文原文）"

    Annotate a genomic variant using Genome Nexus (Memorial Sloan Kettering), which aggregates variant annotations from multiple sources including VEP, SIFT, PolyPhen-2, AlphaMissense, cancer hotspots, and mutation assessor. Input is a GRCh37/hg19 HGVS genomic notation (e.g., '7:g.140453136A&gt;T' for BRAF V600E, '17:g.7577120C&gt;T' for TP53 R273H). Returns the most severe consequence, transcript consequences with pathogenicity scores (SIFT, PolyPhen, AlphaMissense), cancer hotspot information, and annotation summary with HGVS protein notation. This is the primary variant annotation tool for cancer genomics. Note: coordinates must be in GRCh37/hg19 assembly.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GenomeNexus_annotate_variant',
     'arguments': {'hgvsg': '...'}},
)
```

## `GenomeNexus_get_cancer_hotspots` {#GenomeNexus_get_cancer_hotspots}

Returns cancer mutation hotspot annotations (residue, tumor count, type) for an HGVSg variant.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`hgvsg`
- **关键词：**cancer hotspots, recurrent mutation, hgvsg, 3d cluster
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `hgvsg` | string | 是 | Genomic variant in HGVS notation using GRCh37/hg19 coordinates. Format: 'chr:g.posRef&gt;Alt'. Examples: '7:g.140453136A&gt;T' (BRAF V600E), '17:g.7577120C&gt;T' (TP53 R273H). |

??? note "上游说明（英文原文）"

    Get cancer mutation hotspot annotations from Genome Nexus (MSK) for a specific genomic variant. Cancer hotspots are recurrently mutated residues in cancer, curated by the Chang et al. study. Returns hotspot residue, tumor count, hotspot type (single residue or 3D cluster), and the gene symbol. Use this to check if a variant falls in a known cancer hotspot. Input is a GRCh37/hg19 HGVS genomic variant. Example: '7:g.140453136A&gt;T' (BRAF V600E) returns V600 hotspot with 897 tumors.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GenomeNexus_get_cancer_hotspots',
     'arguments': {'hgvsg': '...'}},
)
```

## `GIAB_list_directory` {#GIAB_list_directory}

Browses NIST Genome in a Bottle release tree of high-confidence benchmark VCF and BED files for HG001-HG007.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**giab, benchmark, vcf, high-confidence calls, variant calling, na12878
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `path` | string \| null | 否 | Path relative to the GIAB release root, e.g. '' (top level), 'AshkenazimTrio', or 'AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38'. Default '' lists top-level trios/samples. |

??? note "上游说明（英文原文）"

    Browse NIST Genome in a Bottle's (GIAB) benchmark file release tree -- the reference-standard high-confidence variant (VCF) and callable-region (BED) sets used to validate variant-calling pipelines for samples HG001-HG007. There is no JSON API for this data, only a directory listing, so this tool navigates it: an empty path lists top-level trios/samples; e.g. path='AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38' lists that release's benchmark VCF/BED files with direct download URLs.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GIAB_list_directory', 'arguments': {}},
)
```

## `gnomad_get_constraint` {#gnomad_get_constraint}

Gets a curated gnomAD constraint summary (pLI, oe_lof, LOEUF, missense Z) for a gene symbol or Ensembl ID.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**constraint, pli, loeuf, missense z, loss-of-function intolerance
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 否 | HGNC gene symbol, e.g. 'BRCA1', 'TP53', 'GBA1'. Provide this OR gene_id. |
| `gene_id` | string | 否 | Ensembl gene ID, e.g. 'ENSG00000012048' (BRCA1). Provide this OR gene_symbol. |
| `dataset` | string | 否 | gnomAD dataset release selecting the constraint table. 'gnomad_r4' (default) and 'gnomad_r3' -&gt; GRCh38 / gnomAD v4 constraint; 'gnomad_r2_1' and 'exac' -&gt; GRCh37 / gnomAD v2.1.1 era constraint. 可选值：`gnomad_r4`, `gnomad_r3`, `gnomad_r2_1`, `exac` 默认 `gnomad_r4` |

??? note "上游说明（英文原文）"

    Get gene-level constraint metrics from gnomAD (Genome Aggregation Database) for a single gene, given a gene symbol (e.g. BRCA1) OR an Ensembl gene ID (e.g. ENSG00000012048). Returns a flat, curated summary: pLI (probability of loss-of-function intolerance), oe_lof (observed/expected LoF ratio) with its lower/upper confidence bounds, LOEUF (= oe_lof_upper, the recommended constraint metric; lower = more constrained), missense Z-score (mis_z), synonymous Z-score (syn_z), and expected/observed LoF counts (exp_lof, obs_lof). The `dataset` parameter selects the constraint release: gnomad_r4 (default) and gnomad_r3 use GRCh38 (gnomAD v4 constraint); gnomad_r2_1 and exac use GRCh37. Note: gnomAD uses current HGNC symbols (e.g. GBA1, not GBA). No API key required.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gnomad_get_constraint', 'arguments': {}},
)
```

## `gnomad_get_gene_constraints` {#gnomad_get_gene_constraints}

Gets gnomAD/ExAC gene constraint metrics (pLI, oe_lof, oe_mis, expected/observed counts) for a gene.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`gene_symbol`
- **关键词：**constraint, pli, loeuf, loss-of-function, oe ratio
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 是 | Gene symbol (e.g., 'BRCA1', 'TP53') |
| `reference_genome` | string | 否 | Reference genome. 可选值：`GRCh37`, `GRCh38` 默认 `GRCh38` |

??? note "上游说明（英文原文）"

    Get gene constraint metrics from gnomAD for a gene symbol (GRCh38). Returns constraint fields under `exac_constraint` and `gnomad_constraint` (e.g., pLI, oe_lof/oe_mis/oe_syn, exp/obs counts). Mitochondrially-encoded genes (MT-*) are scored under `mitochondrial_constraint` instead -- their `gnomad_constraint` and `exac_constraint` are always null, which does not mean unconstrained. Note: does NOT include LOEUF, missense Z-score, or synonymous Z-score -- for the full constraint panel (pLI + LOEUF + mis_z + syn_z), use gnomad_get_constraint instead.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gnomad_get_gene_constraints',
     'arguments': {'gene_symbol': '...'}},
)
```

## `gnomad_get_region` {#gnomad_get_region}

Returns genes overlapping a genomic interval and gnomAD variant IDs in it for a chosen dataset.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`chrom`, `start`, `stop`
- **关键词：**gnomad, genomic region, interval, variant ids
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chrom` | string | 是 | Chromosome (e.g., '19'). |
| `start` | integer | 是 | 1-based start position. |
| `stop` | integer | 是 | 1-based stop position. |
| `reference_genome` | string | 否 | Reference genome. 可选值：`GRCh37`, `GRCh38` 默认 `GRCh38` |
| `dataset` | string | 否 | gnomAD dataset ID used for `variants(dataset: ...)`. Allowed values: gnomad_r4, gnomad_r4_non_ukb, gnomad_r3, gnomad_r3_controls_and_biobanks, gnomad_r3_non_cancer, gnomad_r3_non_neuro, gnomad_r3_non_topmed, gnomad_r3_non_v2, gnomad_r2_1, gnomad_r2_1_controls, gnomad_r2_1_non_neuro, gnomad_r2_1_non_cancer, gnomad_r2_1_non_topmed, exac. Coverage: gnomad_r4* = GRCh38 with genomes + exomes; gnomad_r… 默认 `gnomad_r3` |

??? note "上游说明（英文原文）"

    Get basic regional information from gnomAD by genomic interval. Returns genes overlapping the region and a list of `variant_id` values in the region for the selected dataset. Use those IDs with `gnomad_get_variant`. Keep regions small to avoid large payloads. `dataset` defaults to `gnomad_r3` (GRCh38, genomes only), so the variant list is the one that callset contains. Defaults deliberately differ across this gnomAD family: gnomad_get_variant, gnomad_search_variants and gnomad_get_region default to `gnomad_r3` while gnomad_get_variant_populations and gnomad_get_constraint default to `gnomad_r4`, so gnomad_get_variant and gnomad_get_variant_populations report different numbers for the same variant unless you pass the same `dataset` to both. Every response echoes `dataset` and `reference_genome` (gnomad_r2_1* and exac = GRCh37; gnomad_r3* and gnomad_r4* = GRCh38) so the callset that answered is never a guess.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gnomad_get_region',
     'arguments': {'chrom': '...', 'start': 0, 'stop': 0}},
)
```

## `gnomad_get_sv_by_gene` {#gnomad_get_sv_by_gene}

Gets gnomAD structural variants for a gene from the v4 (GRCh38) or v2.1 (GRCh37) SV callset.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`gene_symbol`
- **关键词：**structural variant, sv, deletion, duplication, gnomad sv
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 是 | Gene symbol (e.g., 'BRCA1', 'TP53', 'MECP2') |
| `reference_genome` | string | 否 | Reference genome assembly. 可选值：`GRCh37`, `GRCh38` 默认 `GRCh38` |

??? note "上游说明（英文原文）"

    Get structural variants (SVs) from gnomAD for a gene. The SV callset follows the requested assembly: GRCh38 (default) uses the gnomAD v4 SV callset (gnomad_sv_r4), GRCh37 uses the gnomAD v2.1 SV callset (gnomad_sv_r2_1), so SV coordinates always match the gene coordinates returned. Returns deletions, duplications, inversions, translocations, and complex SVs with population allele frequencies, SV type, consequence, and homozygote counts (null for multi-allelic CNVs). Essential for interpreting CNVs in clinical genomics.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gnomad_get_sv_by_gene',
     'arguments': {'gene_symbol': '...'}},
)
```

## `gnomad_get_sv_by_region` {#gnomad_get_sv_by_region}

Gets gnomAD structural variants overlapping a genomic region in the matching assembly callset.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`chrom`, `start`, `stop`
- **关键词：**structural variant, sv, genomic region, gnomad sv
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chrom` | string | 是 | Chromosome (e.g., '17', 'X'). |
| `start` | integer | 是 | Start position (1-based, GRCh38 by default). |
| `stop` | integer | 是 | Stop position. |
| `reference_genome` | string | 否 | Reference genome assembly. 可选值：`GRCh37`, `GRCh38` 默认 `GRCh38` |

??? note "上游说明（英文原文）"

    Get structural variants (SVs) from gnomAD overlapping a genomic region. The SV callset follows the requested assembly: GRCh38 (default) uses the gnomAD v4 SV callset (gnomad_sv_r4), GRCh37 uses the gnomAD v2.1 SV callset (gnomad_sv_r2_1), so SV coordinates are always in the same assembly as the requested region. Returns DEL/DUP/INV/BND/CPX variants with allele frequencies and consequences. Use for CNV interpretation in specific chromosomal regions. Region should be &lt; 5Mb.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gnomad_get_sv_by_region',
     'arguments': {'chrom': '...', 'start': 0, 'stop': 0}},
)
```

## `gnomad_get_sv_detail` {#gnomad_get_sv_detail}

Gets detailed info for a gnomAD structural variant by ID, auto-selecting v4 or v2.1 callset.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant_id`
- **关键词：**structural variant, sv detail, allele frequency, gnomad sv
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_id` | string | 是 | gnomAD structural variant ID. GRCh38/v4 IDs look like 'DEL_chr17_24e4872b'; GRCh37/v2.1 IDs look like 'DEL_1_1937'. |
| `dataset` | string | 否 | Optional SV callset override. Defaults to the callset auto-detected from the variant ID format (gnomad_sv_r2_1 for GRCh37-style IDs, gnomad_sv_r4 otherwise). Set explicitly only to query a v2.1 subset such as gnomad_sv_r2_1_controls. 可选值：`gnomad_sv_r2_1`, `gnomad_sv_r2_1_controls`, `gnomad_sv_r2_1_non_neuro`, `gnomad_sv_r4` |

??? note "上游说明（英文原文）"

    Get detailed information for a specific gnomAD structural variant by its ID. Resolves IDs from both SV callsets, selected automatically from the ID format: GRCh38 IDs such as 'DEL_chr17_24e4872b' resolve against the gnomAD v4 SV callset (gnomad_sv_r4), while GRCh37 IDs such as 'DEL_1_1937' resolve against the gnomAD v2.1 SV callset (gnomad_sv_r2_1). Returns SV coordinates, type, length, population allele frequency, consequence, and filter status (consequence and homozygote/hemizygote counts are null for multi-allelic CNVs and intergenic SVs). Use gnomad_get_sv_by_gene or gnomad_get_sv_by_region to discover variant IDs first.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gnomad_get_sv_detail',
     'arguments': {'variant_id': '...'}},
)
```

## `gnomad_get_variant` {#gnomad_get_variant}

Gets gnomAD variant metadata and aggregate allele frequency by variant ID (chrom-pos-ref-alt).

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant_id`
- **关键词：**allele frequency, gnomad, population variant, exome, genome
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_id` | string | 是 | Variant ID (e.g., '19-44908822-C-T'). |
| `dataset` | string | 否 | gnomAD dataset ID. Allowed values: gnomad_r4, gnomad_r4_non_ukb, gnomad_r3, gnomad_r3_controls_and_biobanks, gnomad_r3_non_cancer, gnomad_r3_non_neuro, gnomad_r3_non_topmed, gnomad_r3_non_v2, gnomad_r2_1, gnomad_r2_1_controls, gnomad_r2_1_non_neuro, gnomad_r2_1_non_cancer, gnomad_r2_1_non_topmed, exac. Coverage: gnomad_r4* = GRCh38 with genomes + exomes; gnomad_r3* = GRCh38 with genomes ONLY (no… 默认 `gnomad_r3` |

??? note "上游说明（英文原文）"

    Get basic variant metadata from gnomAD by `variant_id` (format like '19-44908822-C-T'). Use `gnomad_search_variants` to discover valid `variant_id` values from an rsID query. `dataset` is required by gnomAD and defaults to `gnomad_r3`. `gnomad_r3` is a GENOMES-ONLY callset with no exome component, so under that default `exome` is null for every variant -- that means 'not present in this callset', NOT 'not observed in exomes'; responses carry a `dataset_note` saying so. Pass dataset='gnomad_r4' for exome allele counts (on GJB2 c.35delG that is a further 10,423 allele observations out of 1,461,600 alleles). Defaults deliberately differ across this gnomAD family: gnomad_get_variant, gnomad_search_variants and gnomad_get_region default to `gnomad_r3` while gnomad_get_variant_populations and gnomad_get_constraint default to `gnomad_r4`, so gnomad_get_variant and gnomad_get_variant_populations report different numbers for the same variant unless you pass the same `dataset` to both. Every response echoes `dataset` and `reference_genome` (gnomad_r2_1* and exac = GRCh37; gnomad_r3* and gnomad_r4* = GRCh38) so the callset that answered is never a guess.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gnomad_get_variant', 'arguments': {'variant_id': '...'}},
)
```

## `gnomad_get_variant_populations` {#gnomad_get_variant_populations}

Gets population-stratified allele frequencies (afr, nfe, eas, sas, etc.) for a gnomAD variant.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant_id`
- **关键词：**allele frequency, ancestry, population, gnomad, af
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_id` | string | 是 | Variant ID (format 'chrom-pos-ref-alt', e.g., '1-55051215-G-GA'). |
| `dataset` | string | 否 | gnomAD dataset ID. Allowed values: gnomad_r4, gnomad_r4_non_ukb, gnomad_r3, gnomad_r3_controls_and_biobanks, gnomad_r3_non_cancer, gnomad_r3_non_neuro, gnomad_r3_non_topmed, gnomad_r3_non_v2, gnomad_r2_1, gnomad_r2_1_controls, gnomad_r2_1_non_neuro, gnomad_r2_1_non_cancer, gnomad_r2_1_non_topmed, exac. Coverage: gnomad_r4* = GRCh38 with genomes + exomes; gnomad_r3* = GRCh38 with genomes ONLY (no… 默认 `gnomad_r4` |

??? note "上游说明（英文原文）"

    Get per-ancestry (population-stratified) allele frequencies for a variant from gnomAD by `variant_id` (format like '1-55051215-G-GA'). Unlike `gnomad_get_variant` (aggregate AF only), this returns per-ancestry rows (`amr`, `nfe`, `fin`, `afr`, `eas`, `sas`, `asj`, `ami`, `mid`, `remaining`, plus `_XX`/`_XY` sex splits) with computed `af = ac / an` (af is null when an = 0), separated into `genome` and `exome` callsets. `dataset` defaults to `gnomad_r4`. Use `gnomad_search_variants` to find a `variant_id` from an rsID. `gnomad_r4` carries both genomes and exomes; `gnomad_r3` is genomes-only, so selecting it makes `exome` null for every variant (a property of the callset, not of the variant) and the response says so in `dataset_note`. Defaults deliberately differ across this gnomAD family: gnomad_get_variant, gnomad_search_variants and gnomad_get_region default to `gnomad_r3` while gnomad_get_variant_populations and gnomad_get_constraint default to `gnomad_r4`, so gnomad_get_variant and gnomad_get_variant_populations report different numbers for the same variant unless you pass the same `dataset` to both. Every response echoes `dataset` and `reference_genome` (gnomad_r2_1* and exac = GRCh37; gnomad_r3* and gnomad_r4* = GRCh38) so the callset that answered is never a guess.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gnomad_get_variant_populations',
     'arguments': {'variant_id': '...'}},
)
```

## `gnomad_search_variants` {#gnomad_search_variants}

Searches gnomAD variants by free text such as an rsID and returns matching variant IDs.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`query`
- **关键词：**rsid, variant search, gnomad, variant id
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Variant search query (e.g., 'rs7412'). |
| `dataset` | string | 否 | gnomAD dataset ID. Allowed values: gnomad_r4, gnomad_r4_non_ukb, gnomad_r3, gnomad_r3_controls_and_biobanks, gnomad_r3_non_cancer, gnomad_r3_non_neuro, gnomad_r3_non_topmed, gnomad_r3_non_v2, gnomad_r2_1, gnomad_r2_1_controls, gnomad_r2_1_non_neuro, gnomad_r2_1_non_cancer, gnomad_r2_1_non_topmed, exac. Coverage: gnomad_r4* = GRCh38 with genomes + exomes; gnomad_r3* = GRCh38 with genomes ONLY (no… 默认 `gnomad_r3` |

??? note "上游说明（英文原文）"

    Search for variants in gnomAD by free-text query (commonly an rsID like 'rs7412'). Returns matching `variant_id` strings. Note: gnomAD's `variant_search` response only exposes `variant_id`. Use the returned `variant_id` with `gnomad_get_variant`. `dataset` defaults to `gnomad_r3` (GRCh38, genomes only); which variants exist depends on the callset searched. Defaults deliberately differ across this gnomAD family: gnomad_get_variant, gnomad_search_variants and gnomad_get_region default to `gnomad_r3` while gnomad_get_variant_populations and gnomad_get_constraint default to `gnomad_r4`, so gnomad_get_variant and gnomad_get_variant_populations report different numbers for the same variant unless you pass the same `dataset` to both. Every response echoes `dataset` and `reference_genome` (gnomad_r2_1* and exac = GRCh37; gnomad_r3* and gnomad_r4* = GRCh38) so the callset that answered is never a guess.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gnomad_search_variants', 'arguments': {'query': '...'}},
)
```

## `gProfiler_annotate_snps` {#gProfiler_annotate_snps}

Maps SNP rsIDs to positions, Ensembl genes and predicted variant consequences using g:SNPense.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`snp_list`
- **关键词：**rsid, snp annotation, variant consequence, g:snpense
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `snp_list` | string | 是 | Comma-separated list of SNP rsIDs. Examples: 'rs11540652,rs429358,rs7903146' or 'rs1042779,rs6792369'. Maximum ~200 SNPs per query. |
| `organism` | string | 否 | Organism identifier. Default: 'hsapiens'. Examples: 'hsapiens' (human), 'mmusculus' (mouse). |

??? note "上游说明（英文原文）"

    Map SNP rsIDs to genes and annotate their functional consequences using g:Profiler (g:SNPense) from the University of Tartu. Given a list of SNP rsIDs, returns the chromosome, position, strand, associated Ensembl gene IDs, gene names, and predicted variant consequences (missense_variant, synonymous_variant, intron_variant, etc.) with counts. Useful for annotating GWAS hits, linking variants to genes, or understanding the functional impact of SNPs. Example: 'rs11540652' maps to chromosome 17 position 7674220, gene TP53 (ENSG00000141510), with 18 missense variant consequences; 'rs429358' maps to APOE on chromosome 19.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gProfiler_annotate_snps',
     'arguments': {'snp_list': '...'}},
)
```

## `InterVar_classify_variant` {#InterVar_classify_variant}

Classifies a germline variant given chrom/pos/ref/alt with InterVar ACMG/AMP 2015 rules and reports evidence criteria.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`chrom`, `pos`, `ref`, `alt`
- **关键词：**acmg, pathogenicity, germline variant, clinical significance, pvs1
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chrom` | string | 是 | Chromosome number or name. Can include 'chr' prefix (e.g., '17', 'chr17', 'X'). |
| `pos` | integer | 是 | Genomic position (1-based). |
| `ref` | string | 是 | Reference allele (e.g., 'G', 'A', 'ATCG'). |
| `alt` | string | 是 | Alternate allele (e.g., 'A', 'T', 'C'). |
| `build` | string | 否 | Genome build. Default: hg19. 可选值：`hg19`, `hg38` 默认 `hg19` |

??? note "上游说明（英文原文）"

    Classify a germline variant using ACMG/AMP 2015 standards (InterVar). Returns the clinical significance category — Pathogenic, Likely Pathogenic, Uncertain Significance, Likely Benign, or Benign — plus the activation status of all 28 evidence criteria (PVS1, PS1-PS4, PM1-PM6, PP1-PP5, BA1, BP1-BP7, BS1-BS4). Input requires genomic coordinates in hg19 or hg38. No API key needed. Example: BRCA1 c.5266dupC (chr17:41245466 G&gt;A in hg19) is Likely Pathogenic with PM2+PP3+PP5.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'InterVar_classify_variant',
     'arguments': {'chrom': '...', 'pos': 0, 'ref': '...', 'alt': '...'}},
)
```

## `IPD_get_hla_allele` {#IPD_get_hla_allele}

Fetches a full IPD-IMGT/HLA allele record by accession including name, locus, class, sequence and release history.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`accession`
- **关键词：**hla allele, imgt, locus, allele sequence, hla typing
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `accession` | string | 是 | IPD-IMGT/HLA allele accession (e.g. 'HLA00001'). Obtain via IPD_search_hla_alleles. |

??? note "上游说明（英文原文）"

    Fetch the full IPD-IMGT/HLA record for a single HLA allele by its IPD accession (e.g. 'HLA00001'). Returns the complete record including allele name, locus, class, sequence, confirmation status, citations, CWD (common/well-documented) status, release history, and cross-references. Use IPD_search_hla_alleles first to obtain an accession. Keyless EMBL-EBI API.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'IPD_get_hla_allele', 'arguments': {'accession': '...'}},
)
```

## `IPD_search_cells` {#IPD_search_cells}

Searches IPD HLA reference cell lines by name or other fields and returns cell IDs and names.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`query`
- **关键词：**hla reference cells, hla typing, cell line, genotype
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Text to search for (e.g. a cell primary name, or a lab name when field='lab_of_origin'). The filter expression is built for you. |
| `field` | string | 否 | Cell field to search (default 'primary_name'). Other useful fields include 'lab_of_origin'. 默认 `primary_name` |
| `match` | string | 否 | How to match: 'contains' (default, substring), 'startsWith' (prefix), or 'eq' (exact). 可选值：`startsWith`, `contains`, `eq` 默认 `contains` |
| `limit` | integer | 否 | Maximum number of cells to return (default 10, max 100). 默认 `10` |

??? note "上游说明（英文原文）"

    Search IPD cell lines (HLA reference cells with their HLA typing) by a field value. By default searches the cell's primary name; set 'field' to search other fields such as 'lab_of_origin'. Returns matching cells with their IPD cell id and primary name. Useful for finding HLA-typed reference cell lines. Keyless EMBL-EBI API.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'IPD_search_cells', 'arguments': {'query': '...'}},
)
```

## `IPD_search_hla_alleles` {#IPD_search_hla_alleles}

Searches IPD-IMGT/HLA for human HLA alleles by allele name prefix or substring and returns accessions and names.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`name`
- **关键词：**hla allele, hla typing, allele name, imgt
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `name` | string | 是 | HLA allele name or partial name to search (e.g. 'A*01:01', 'DRB1*15:01', 'B*07'). The filter expression is built for you. |
| `match` | string | 否 | How to match the name: 'startsWith' (default, prefix match), 'contains' (substring), or 'eq' (exact full name). 可选值：`startsWith`, `contains`, `eq` 默认 `startsWith` |
| `limit` | integer | 否 | Maximum number of alleles to return (default 10, max 100). 默认 `10` |

??? note "上游说明（英文原文）"

    Search the IPD-IMGT/HLA database for human HLA alleles by allele name prefix or substring (e.g. 'A*01:01', 'DRB1*15', 'B*07'). Returns matching alleles with their IPD accession (e.g. HLA00001) and full allele name. Use this to resolve an HLA type or partial type into specific allele accessions for follow-up detail lookups. Keyless EMBL-EBI API.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'IPD_search_hla_alleles', 'arguments': {'name': '...'}},
)
```

## `KEGG_get_variant` {#KEGG_get_variant}

Returns a KEGG variant entry with gene, mutation type and ClinVar/dbSNP/COSMIC cross-references, networks and drugs.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant_id`
- **关键词：**kegg variant, gain of function, clinvar, cosmic, dbsnp
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_id` | string | 是 | KEGG variant ID. Can be with prefix 'hsa_var:673v1' or without '673v1'. Get IDs from KEGG_search_variant. Format: \{entrez_gene_id\}v\{variant_number\}. |

??? note "上游说明（英文原文）"

    Get detailed KEGG variant information including gene name, mutation type (gain/loss of function), specific mutations with ClinVar/dbSNP/COSMIC cross-references, associated signaling networks, diseases, and targeted drugs. Example: hsa_var:673v1 (BRAF) returns V600E and V600K mutations with ClinVar:13961, dbSNP:rs113488022, COSM:476, linked to melanoma and thyroid cancer, targeted by dabrafenib and encorafenib.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'KEGG_get_variant', 'arguments': {'variant_id': '...'}},
)
```

## `KEGG_search_variant` {#KEGG_search_variant}

Searches the KEGG VARIANT database for human gene variants associated with disease.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`keyword`
- **关键词：**variant search, mutation, loss of function, cancer
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `keyword` | string | 是 | Gene name to search (e.g., 'BRAF', 'TP53', 'EGFR', 'KRAS', 'PIK3CA'). Case-insensitive. |
| `max_results` | integer | 否 | Maximum number of results (default: 25). 默认 `25` |

??? note "上游说明（英文原文）"

    Search the KEGG VARIANT database for human gene variants associated with disease. KEGG VARIANT catalogs gain-of-function mutations, loss-of-function mutations, and other variants in human genes linked to cancer and other diseases. Returns variant IDs and descriptions. Example: searching 'BRAF' returns 673v1 (BRAF mutation) and 673v2 (BRAF mutation). Use KEGG_get_variant for mutation details, ClinVar/dbSNP cross-references, and associated diseases/drugs.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'KEGG_search_variant', 'arguments': {'keyword': '...'}},
)
```

## `LitVar_get_variant_details` {#LitVar_get_variant_details}

Returns the LitVar2 record for an rsID with ClinGen allele IDs, associated genes, HGVS names and clinical annotations.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`rsid`
- **关键词：**rsid, dbsnp, hgvs, clingen allele, litvar
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `rsid` | string | 是 | dbSNP rsID of the variant (e.g., 'rs113488022' = BRAF V600E, 'rs121913529' = KRAS G12, 'rs328', 'rs7903146'). Provide the rsID only; the litvar@...\#\# wrapper is added automatically. |

??? note "上游说明（英文原文）"

    Get the full structured variant record from NCBI LitVar2 by rsID (the 'variant/get' endpoint, WITHOUT the publications list). Given a dbSNP rsID, returns a single normalized record combining: ClinGen canonical allele IDs (CA...), associated gene(s), variant name and HGVS notation, genomic chromosome:position (data_chromosome_base_position), allele, SNP class (e.g. snv), species and NCBI tax id, the literature-derived clinical-significance set (data_clinical_significance), and, when available, ALFA minor allele frequency (data_maf). This is the variant 'profile' record - use LitVar_get_variant_publications for the list of citing PMIDs, or LitVar_search_variants to find an rsID from a gene or HGVS query.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'LitVar_get_variant_details',
     'arguments': {'rsid': '...'}},
)
```

## `LitVar_search_variants` {#LitVar_search_variants}

Searches LitVar2 variants by rsID, gene name or HGVS and returns identifiers, genes, HGVS and publication counts.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`query`
- **关键词：**rsid, hgvs, variant search, gene name
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Search query: rsID (e.g., 'rs328', 'rs7903146'), gene name (e.g., 'BRCA1', 'TP53'), or HGVS notation (e.g., 'p.V600E') |

??? note "上游说明（英文原文）"

    Search for genetic variants in the NCBI LitVar2 database by rsID, gene name, or HGVS notation. LitVar2 links genetic variants to their mentions in PubMed literature. Returns variant identifiers, associated gene names, HGVS nomenclature, publication counts, and clinical significance when available. Useful for finding which variants have been studied for a gene or identifying literature-linked variants by rsID.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'LitVar_search_variants', 'arguments': {'query': '...'}},
)
```

## `LOVD_get_variants` {#LOVD_get_variants}

Returns all curated LOVD variants for a gene with HGVS DNA/RNA/protein notation and hg19 positions.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`gene_symbol`
- **关键词：**lovd, hgvs, curated variants, hg19
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 是 | HGNC gene symbol (e.g., 'TP53', 'BRCA1'). |

??? note "上游说明（英文原文）"

    Get all curated variants for a gene from LOVD. Returns variant details including HGVS DNA/RNA/protein notation, genomic positions (hg19), LOVD variant DBID, and reporting count. Note: some genes may have hundreds or thousands of variants.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'LOVD_get_variants', 'arguments': {'gene_symbol': '...'}},
)
```

## `LOVD_search_variants` {#LOVD_search_variants}

Searches LOVD variants by variant DBID or HGVS DNA notation.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`gene_symbol`
- **关键词：**lovd, hgvs, variant dbid, c. notation
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 是 | HGNC gene symbol (e.g., 'TP53', 'BRCA1'). |
| `variant_dbid` | string | 否 | LOVD variant DBID (e.g., 'TP53_010464'). Unique identifier assigned by LOVD. |
| `dna_notation` | string | 否 | HGVS DNA notation without RefSeq prefix (e.g., 'c.*2609C&gt;A', 'c.742C&gt;T'). Use coding 'c.' notation. |

??? note "上游说明（英文原文）"

    Search variants in LOVD by DBID or HGVS DNA notation. Use variant_dbid for exact LOVD variant IDs (e.g., 'TP53_010464') or dna_notation for HGVS notation without RefSeq prefix (e.g., 'c.*2609C&gt;A'). At least one search parameter is required.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'LOVD_search_variants',
     'arguments': {'gene_symbol': '...'}},
)
```

## `MaveDB_get_clinical_controls` {#MaveDB_get_clinical_controls}

Returns ClinVar pathogenic/benign/VUS variants overlapping a MaveDB score set for score calibration.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`urn`
- **关键词：**clinvar, mave, calibration, clinical significance
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `urn` | string | 是 | MaveDB score set URN (e.g., 'urn:mavedb:00001263-a-2' for BRCA2 SGE). Obtain from MaveDB_search_score_sets. |
| `clinical_significance` | string \| null | 否 | Optional exact-match filter on ClinVar clinical significance (case-insensitive), e.g. 'Pathogenic', 'Benign', 'Uncertain significance'. |
| `limit` | integer \| null | 否 | Maximum number of clinical controls to return (client-side truncation). Set to 0 or null to return all. 默认 `0` |

??? note "上游说明（英文原文）"

    Get clinical-control variants (ClinVar pathogenic/benign/VUS) that overlap a MaveDB score set, used to calibrate MAVE functional scores into clinically-interpretable thresholds. Returns ClinVar entries with dbIdentifier, gene symbol, clinical significance (e.g. 'Uncertain significance', 'Pathogenic'), review status, ClinVar version (e.g. '01_2025'), and the MaveDB variant URN(s) sharing the same genomic coordinate. Optionally filter by clinical_significance. This is the reference set for building functional-score evidence (e.g. ACMG PS3/BS3). Obtain the URN from MaveDB_search_score_sets.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MaveDB_get_clinical_controls',
     'arguments': {'urn': '...'}},
)
```

## `MaveDB_get_effect_matrix` {#MaveDB_get_effect_matrix}

Loads a MaveDB score set and reshapes single missense variant scores into an amino acid by position effect matrix.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`urn`
- **关键词：**deep mutational scanning, missense, effect matrix, hgvs
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `urn` | string | 是 | MaveDB score set URN (e.g., 'urn:mavedb:00000115-a-7'). Obtain from MaveDB_search_score_sets. |
| `uniprot_id` | string \| null | 否 | Optional UniProt accession (e.g., 'P01116' for KRAS). If supplied, the tool fetches the canonical sequence and verifies the MaveDB position numbering against it via a landmark check; sets numbering_offset and surfaces a warning in numbering_check if the positions don't match (silent off-by-N is the \#1 silent bug in DMS analysis). |
| `score_field` | string \| null | 否 | Optional explicit column name to use as the variant score (e.g., 'score', 'ddG', 'fitness'). If omitted, tries 'score' → 'ddG' → 'fitness' → first numeric non-HGVS field in that order. The actual field used is reported in data.score_field_used. |

??? note "上游说明（英文原文）"

    One-shot DMS data loader: fetch a MaveDB score set, parse HGVS to (ref_aa, position, alt_aa), filter to single missense, optionally verify position numbering against a UniProt canonical sequence, and reshape into a (20 amino acids × n_positions) effect matrix ready for downstream analysis. Hides the boilerplate that variant-effect-prediction benchmarking and residue-mechanism interpretation skills used to inline. Returns the matrix, the position list, the standard amino acid order, and explicit numbering / drop / score-field metadata so the caller can audit what was loaded.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MaveDB_get_effect_matrix', 'arguments': {'urn': '...'}},
)
```

## `MaveDB_get_gnomad_variants` {#MaveDB_get_gnomad_variants}

Returns gnomAD variants with allele counts and frequencies observed within a MaveDB score set.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`urn`
- **关键词：**gnomad, allele frequency, mave, population
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `urn` | string | 是 | MaveDB score set URN (e.g., 'urn:mavedb:00001263-a-2' for BRCA2 SGE). Obtain from MaveDB_search_score_sets. |
| `limit` | integer \| null | 否 | Maximum number of gnomAD variants to return (client-side truncation). Set to 0 or null to return all. 默认 `0` |

??? note "上游说明（英文原文）"

    Get gnomAD population-frequency variants observed within a MaveDB score set. Returns gnomAD variants with dbIdentifier (e.g. '13-32319097-A-C'), gnomAD version (e.g. 'v4.1'), allele count, allele number, and allele frequency, each linked to the MaveDB mapped variant via a genomic VRS Allele. Cross-reference allele frequency with the measured MAVE functional score: rare variants with strong functional impact are candidate disease alleles. Obtain the URN from MaveDB_search_score_sets.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MaveDB_get_gnomad_variants', 'arguments': {'urn': '...'}},
)
```

## `MaveDB_get_mapped_variants` {#MaveDB_get_mapped_variants}

Returns GA4GH VRS genomic allele mappings and ClinGen allele IDs for variants in a MaveDB score set.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`urn`
- **关键词：**vrs allele, genomic coordinates, clingen, mave
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `urn` | string | 是 | MaveDB score set URN (e.g., 'urn:mavedb:00001263-a-2' for BRCA2 SGE). Obtain from MaveDB_search_score_sets. |
| `limit` | integer \| null | 否 | Maximum number of mapped variants to return (client-side truncation). Set to 0 or null to return all (the API returns the full list in one call). 默认 `0` |

??? note "上游说明（英文原文）"

    Get genomically-mapped variant coordinates (GA4GH VRS Allele + ClinGen Allele ID) for a MaveDB score set. Each MAVE variant is mapped to a postMapped VRS Allele carrying a genomic SequenceLocation (chromosome + start/end on the reference assembly, e.g. NC_000013.11 for chr13) plus the ClinGen canonical Allele ID (e.g. CA387754009), which links the functional score to ClinVar, dbSNP, and gnomAD. Intronic or unmappable variants have post_mapped=null and an error_message. Use this to translate raw MAVE scores into genome coordinates for downstream variant annotation. Obtain the URN from MaveDB_search_score_sets.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MaveDB_get_mapped_variants', 'arguments': {'urn': '...'}},
)
```

## `MaveDB_get_score_set` {#MaveDB_get_score_set}

Returns MaveDB score set details including target gene, methods, variant count and publications.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`urn`
- **关键词：**score set, target gene, saturation mutagenesis, mavedb
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `urn` | string | 是 | MaveDB score set URN (e.g., 'urn:mavedb:00000001-a-1'). Obtain URNs from MaveDB_search_score_sets. |

??? note "上游说明（英文原文）"

    Get detailed information about a specific MaveDB score set by its URN identifier. Returns the target gene, experimental methods, variant count, publications, and license. Score sets contain functional impact scores from deep mutational scanning or saturation genome editing experiments.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MaveDB_get_score_set', 'arguments': {'urn': '...'}},
)
```

## `MaveDB_get_variant_scores` {#MaveDB_get_variant_scores}

Returns HGVS-annotated variant functional scores from a MaveDB score set.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`urn`
- **关键词：**variant effect, functional score, hgvs, dms
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `urn` | string | 是 | MaveDB score set URN (e.g., 'urn:mavedb:00001234-c-1'). Obtain from MaveDB_search_score_sets. |
| `hgvs_pro` | string | 否 | Optional filter: HGVS protein notation substring (e.g., 'Arg175', 'p.Arg175His'). Case-insensitive partial match on protein-level HGVS. |
| `limit` | integer \| null | 否 | Maximum variants to return after client-side filtering. Set to 0 or null to return all variants in the score set (default — needed for whole-protein DMS workflows). Set to a positive integer to truncate. The MaveDB API returns the full CSV in one HTTP call regardless, so 'unlimited' has no extra cost. 默认 `0` |

??? note "上游说明（英文原文）"

    Get functional variant effect scores from a MaveDB score set. Returns HGVS-annotated variants with their functional scores from deep mutational scanning or other MAVE assays. Scores indicate variant effect on protein function (interpretation depends on the specific assay - check score set metadata). Optionally filter by HGVS protein notation (e.g., 'Arg175' to find all R175 substitutions). The /scores endpoint returns every variant in one CSV; `limit` is client-side truncation only — set `limit=0` (or omit) to receive ALL variants in one call (use this for whole-protein DMS analyses, e.g. KRAS folding ΔΔG has ~2200 variants).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MaveDB_get_variant_scores', 'arguments': {'urn': '...'}},
)
```

## `MaveDB_search_experiments` {#MaveDB_search_experiments}

Searches MaveDB experiments grouping related variant effect score sets.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`query`
- **关键词：**mavedb experiment, mave, score sets
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Search query: gene symbol (e.g., 'TP53', 'BRCA1'), protein name, or keywords describing the experiment. |

??? note "上游说明（英文原文）"

    Search MaveDB for experiments (groups of related score sets from the same study). Each experiment may contain multiple score sets from different replicates, conditions, or analysis methods. Returns experiment URNs, titles, associated score sets, and publications. Search by gene name, protein, or assay keywords.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MaveDB_search_experiments',
     'arguments': {'query': '...'}},
)
```

## `MaveDB_search_score_sets` {#MaveDB_search_score_sets}

Searches MaveDB score sets from MAVE and deep mutational scanning experiments.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`query`
- **关键词：**mave, deep mutational scanning, score set search, gene
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Search text - gene name (e.g., 'BRCA1', 'TP53'), protein name, or keyword. Searches titles and descriptions of score sets. |
| `limit` | integer | 否 | Maximum number of results to return. 默认 `20` |

??? note "上游说明（英文原文）"

    Search MaveDB for variant effect score sets from Multiplexed Assays of Variant Effect (MAVE) experiments, including deep mutational scanning. Returns score set URNs, titles, variant counts, and publication dates. Use gene names or protein names as search terms.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MaveDB_search_score_sets', 'arguments': {'query': '...'}},
)
```

## `Mutalyzer_back_translate` {#Mutalyzer_back_translate}

Converts a protein-level HGVS variant to possible coding DNA variant descriptions with Mutalyzer.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant`
- **关键词：**hgvs, back translation, p. to c., mutalyzer
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant` | string | 是 | Protein-level HGVS variant description. Examples: 'NP_002993.1:p.Asp92Tyr' (SDHD), 'NP_000537.3:p.Pro72Arg' (TP53), 'NP_000050.2:p.Tyr32Cys' (BRCA2). Format: NP_accession:p.AminoAcidChange |

??? note "上游说明（英文原文）"

    Convert a protein-level HGVS variant to possible DNA-level variant descriptions using Mutalyzer. Given a protein change like NP_002993.1:p.Asp92Tyr, returns the corresponding coding DNA change(s) like NM_003002.4:c.(274G&gt;T). Useful for finding the DNA-level representation of a known amino acid change, accounting for codon degeneracy.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Mutalyzer_back_translate',
     'arguments': {'variant': '...'}},
)
```

## `Mutalyzer_normalize_variant` {#Mutalyzer_normalize_variant}

Normalizes an HGVS variant with Mutalyzer and returns corrected notation, protein effect and gene.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant`
- **关键词：**hgvs normalization, variant validation, mutalyzer
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant` | string | 是 | HGVS variant description. Examples: 'NM_003002.4:c.274G&gt;T' (SDHD missense), 'NM_000546.6:c.215C&gt;G' (TP53 p.Pro72Arg), 'NM_000059.4:c.7397C&gt;T' (BRCA2). Format: reference:coordinate_system.change |

??? note "上游说明（英文原文）"

    Validate and normalize an HGVS variant description using Mutalyzer. Returns the corrected/normalized HGVS notation, predicted protein effect, RNA description, and gene symbol. Essential for clinical variant interpretation - ensures HGVS nomenclature compliance and predicts amino acid changes. Accepts coding DNA (c.), genomic (g.), protein (p.), or RNA (r.) variant descriptions. Example: 'NM_003002.4:c.274G&gt;T' normalizes to a SDHD missense variant with protein prediction p.(Asp92Tyr).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Mutalyzer_normalize_variant',
     'arguments': {'variant': '...'}},
)
```

## `Mutalyzer_parse_hgvs` {#Mutalyzer_parse_hgvs}

Parses an HGVS variant into reference, coordinate system, variant type and positions.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant`
- **关键词：**hgvs parsing, variant type, coordinates
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant` | string | 是 | HGVS variant description to parse. Examples: 'NM_003002.4:c.274G&gt;T', 'NC_000017.11:g.7674220C&gt;G', 'NM_000546.6:c.215C&gt;G'. Supports coding (c.), genomic (g.), protein (p.), and RNA (r.) coordinates. |

??? note "上游说明（英文原文）"

    Parse an HGVS variant description into a structured model using Mutalyzer. Returns the reference sequence ID, coordinate system (c./g./p./r.), variant type (substitution, deletion, insertion, etc.), and exact positions. Useful for programmatically extracting components of HGVS notation without regex parsing.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Mutalyzer_parse_hgvs', 'arguments': {'variant': '...'}},
)
```

## `MyVariant_get_pathogenicity_scores` {#MyVariant_get_pathogenicity_scores}

Retrieves dbNSFP pathogenicity scores (REVEL, CADD, SIFT, PolyPhen-2, etc.) for a variant.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant_id`
- **关键词：**dbnsfp, revel, sift, polyphen, conservation
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_id` | string | 是 | Variant ID: rsID (e.g., rs45478192) or HGVS genomic (e.g., chr16:g.23635348A&gt;C). Both formats work; an HGVS id must be expressed in the assembly given by the 'assembly' parameter. |
| `fields` | string | 否 | Fields to return (pre-configured for pathogenicity scores) 默认 `dbnsfp.revel.score,cadd.phred,dbnsfp.alphamissense.score,db…` |
| `assembly` | string | 否 | Reference genome assembly for both the HGVS variant_id and the returned coordinates: 'hg19' (GRCh37, the default) or 'hg38' (GRCh38). 可选值：`hg19`, `hg38` 默认 `hg19` |

??? note "上游说明（英文原文）"

    Get comprehensive pathogenicity prediction scores for a variant from dbNSFP. Returns REVEL, CADD, AlphaMissense, SIFT, PolyPhen-2, MetaRNN, GERP, PhyloP, and PhastCons scores in a single call. Use this when MyVariant_query_variants does not return dbnsfp data, or when you specifically need REVEL/AlphaMissense scores for ACMG PP3/BP4 classification. Coordinates default to hg19/GRCh37 -- set assembly='hg38' for hg38/GRCh38 (e.g. rs4244285 is chr10:g.96541616G&gt;C on hg19 and chr10:g.94781859G&gt;C on hg38). An HGVS variant_id must be written in the SAME assembly you request: this endpoint looks ids up verbatim and never lifts coordinates over, so an hg19 id with assembly='hg38' is a 404, not a translation. An rsID works in either assembly. The response reports the frame in a top-level 'coordinate_assembly' key. MyVariant sometimes files several records under one rsID and resolves the rsID to the highest-scoring one, which can be an identifier-only stub; when that happens the sibling records that do carry the requested fields are reported in the top-level 'sibling_records_with_requested_fields', 'sibling_variant_ids_with_requested_fields' and 'sibling_record_note' keys, and 'data' still holds the record MyVariant resolved.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MyVariant_get_pathogenicity_scores',
     'arguments': {'variant_id': '...'}},
)
```

## `MyVariant_get_variant_annotation` {#MyVariant_get_variant_annotation}

Retrieves MyVariant.info annotation for a variant by HGVS ID or rsID.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant_id`
- **关键词：**rsid, hgvs, dbsnp, clinvar, variant annotation
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_id` | string | 是 | HGVS variant ID (e.g., 'chr1:g.2186318G&gt;A') or rsID (e.g., 'rs123'). An HGVS id must be expressed in the assembly given by the 'assembly' parameter. |
| `fields` | string | 否 | Fields to return. Useful fields: dbsnp, clinvar, cadd, gnomad_genome, gnomad_exome, dbnsfp (functional predictions). 默认 `dbsnp,clinvar,cadd,gnomad_genome,dbnsfp` |
| `assembly` | string | 否 | Reference genome assembly for both the HGVS variant_id and the returned coordinates: 'hg19' (GRCh37, the default) or 'hg38' (GRCh38). 可选值：`hg19`, `hg38` 默认 `hg19` |

??? note "上游说明（英文原文）"

    Get detailed annotations for a specific variant by its HGVS ID (e.g., 'chr1:g.2186318G&gt;A') or rsID. Returns information from dbSNP, ClinVar, CADD, and other sources. Coordinates default to hg19/GRCh37 -- set assembly='hg38' for hg38/GRCh38 (e.g. rs4244285 is chr10:g.96541616G&gt;C on hg19 and chr10:g.94781859G&gt;C on hg38). An HGVS variant_id must be written in the SAME assembly you request: this endpoint looks ids up verbatim and never lifts coordinates over, so an hg19 id with assembly='hg38' is a 404, not a translation. An rsID works in either assembly. The response reports the frame in a top-level 'coordinate_assembly' key. MyVariant sometimes files several records under one rsID and resolves the rsID to the highest-scoring one, which can be an identifier-only stub; when that happens the sibling records that do carry the requested fields are reported in the top-level 'sibling_records_with_requested_fields', 'sibling_variant_ids_with_requested_fields' and 'sibling_record_note' keys, and 'data' still holds the record MyVariant resolved.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MyVariant_get_variant_annotation',
     'arguments': {'variant_id': '...'}},
)
```

## `MyVariant_query_variants` {#MyVariant_query_variants}

Searches MyVariant.info for variants by rsID, gene, position or clinical significance.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`query`
- **关键词：**variant search, gnomad, clinical significance, rsid
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Search query. Examples: 'rs58991260' (rsID), 'chr7:g.55249071G&gt;A' (HGVS), 'clinvar.gene.symbol:BRCA1' (gene in ClinVar), 'cadd.phred:&gt;30' (high CADD score). |
| `fields` | string | 否 | Comma-separated fields to return. Common: dbsnp, clinvar, cadd, gnomad_genome, dbnsfp. 默认 `dbsnp.rsid,clinvar.rcv.clinical_significance,cadd.phred,gno…` |
| `size` | integer | 否 | Maximum number of results (1-100). 默认 `10` |
| `assembly` | string | 否 | Reference genome assembly for the coordinates in the query and in the results: 'hg19' (GRCh37, the default) or 'hg38' (GRCh38). 可选值：`hg19`, `hg38` 默认 `hg19` |

??? note "上游说明（英文原文）"

    Search for genetic variants by various criteria including rsID, gene, chromosome position, or clinical significance. Returns variant annotations from sources like dbSNP, ClinVar, gnomAD, and CADD. Coordinates default to hg19/GRCh37 -- set assembly='hg38' for hg38/GRCh38 (e.g. rs4244285 is chr10:g.96541616G&gt;C on hg19 and chr10:g.94781859G&gt;C on hg38). Any HGVS coordinate in the query must be written in the assembly you request; an rsID works in either. The response reports the frame in a top-level 'coordinate_assembly' key.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MyVariant_query_variants', 'arguments': {'query': '...'}},
)
```

## `NCBIVariation_alfa_frequencies_by_rsid` {#NCBIVariation_alfa_frequencies_by_rsid}

Returns ALFA per-population allele counts and frequencies for a dbSNP rsID.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`rsid`
- **关键词：**alfa, allele frequency, ancestry, dbsnp
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `rsid` | string | 是 | dbSNP reference SNP ID. Examples: 'rs429358' (APOE, per-ancestry C allele frequency), 'rs334' (sickle-cell HBB, per-ancestry A/T counts), 'rs7412' (APOE epsilon-2) |

??? note "上游说明（英文原文）"

    Get ALFA (Allele Frequency Aggregator) per-ancestry allele frequencies for a dbSNP rsID from NCBI Variation Services. Unlike dbSNP esummary global_mafs (one aggregate global frequency per study), this returns allele COUNTS and computed frequencies broken down by population/ancestry: European, African, African American, East Asian, Latin American, South Asian, Other Asian, etc. SAMN biosample IDs are resolved to ancestry names via the live ALFA frequency metadata, aggregating &gt;1M dbGaP subjects. Accepts rsID with or without 'rs' prefix.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NCBIVariation_alfa_frequencies_by_rsid',
     'arguments': {'rsid': '...'}},
)
```

## `NCBIVariation_hgvs_to_spdi` {#NCBIVariation_hgvs_to_spdi}

Converts and validates an HGVS variant description to SPDI using NCBI Variation Services.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`hgvs`
- **关键词：**hgvs, spdi, variant notation
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `hgvs` | string | 是 | HGVS variant description. Examples: 'NC_000001.11:g.230710048A&gt;G' (genomic), 'NM_000059.4:c.68_69del' (coding deletion in BRCA2), 'NM_000546.6:c.215C&gt;G' (TP53 Pro72Arg) |

??? note "上游说明（英文原文）"

    Convert an HGVS variant description to SPDI notation using NCBI Variation Services. Accepts any valid HGVS description (genomic g., coding c., or RNA r.) and returns the equivalent SPDI representation(s). Also validates the input HGVS expression. Useful for normalizing variant representations and converting between notation systems.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NCBIVariation_hgvs_to_spdi',
     'arguments': {'hgvs': '...'}},
)
```

## `NCBIVariation_rsid_lookup` {#NCBIVariation_rsid_lookup}

Returns dbSNP rsID details with GRCh38 placements, genes, ClinVar significance and MANE transcripts.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`rsid`
- **关键词：**dbsnp, rsid, grch38, clinvar, mane
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `rsid` | string | 是 | dbSNP reference SNP ID. Examples: 'rs429358' (APOE epsilon-4 allele), 'rs7412' (APOE epsilon-2), 'rs699' (AGT M235T), 'rs1801133' (MTHFR C677T) |

??? note "上游说明（英文原文）"

    Look up a dbSNP rsID and retrieve variant details including GRCh38 genomic coordinates, associated genes, clinical significance (ClinVar), citations, and MANE Select transcript IDs. Each entry in 'grch38_placements' reports 'position' as a 1-based GRCh38 coordinate (matching ClinVar, Ensembl VEP and MyVariant) and keeps dbSNP's 0-based interbase SPDI offset in 'spdi_position'. Accepts rsID with or without 'rs' prefix (e.g., 'rs429358' or '429358'). Returns structured summary of the variant.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NCBIVariation_rsid_lookup', 'arguments': {'rsid': '...'}},
)
```

## `NCBIVariation_spdi_canonical` {#NCBIVariation_spdi_canonical}

Returns the canonical normalized SPDI form of a variant.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`spdi`
- **关键词：**spdi, normalization, indel, canonical
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `spdi` | string | 是 | Variant in SPDI notation to normalize. Examples: 'NC_000001.11:230710047:A:G', 'NM_000546.6:214:C:G' (TP53) |

??? note "上游说明（英文原文）"

    Get the canonical (normalized) SPDI representation for a variant. Returns the single canonical form of a variant after normalization (right-shifting for indels, using the top-level chromosome sequence). Useful for variant deduplication and consistent storage.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NCBIVariation_spdi_canonical',
     'arguments': {'spdi': '...'}},
)
```

## `NCBIVariation_spdi_equivalents` {#NCBIVariation_spdi_equivalents}

Returns all equivalent SPDI representations of a variant across GRCh37/38, RefSeqGene and transcripts.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`spdi`
- **关键词：**spdi equivalents, liftover, grch37, transcript
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `spdi` | string | 是 | Variant in SPDI notation. Examples: 'NC_000001.11:230710047:A:G' (GRCh38 AGT variant), 'NM_000059.4:266:AT:' (BRCA2 transcript-level) |

??? note "上游说明（英文原文）"

    Get all equivalent SPDI representations of a variant across different reference sequences and genome assemblies (GRCh37/GRCh38, RefSeqGene, transcript). Given one SPDI, returns the same variant mapped to all available coordinate systems. Useful for liftover between assemblies and cross-referencing variants described on different reference sequences.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NCBIVariation_spdi_equivalents',
     'arguments': {'spdi': '...'}},
)
```

## `NCBIVariation_spdi_to_hgvs` {#NCBIVariation_spdi_to_hgvs}

Converts an SPDI variant to HGVS notation using NCBI Variation Services.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`spdi`
- **关键词：**spdi, hgvs, notation conversion
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `spdi` | string | 是 | Variant in SPDI notation (SeqID:Position:Deleted:Inserted). Examples: 'NC_000001.11:230710047:A:G' (AGT M259T), 'NC_000017.11:7674220:G:C' (TP53 region), 'NM_000059.4:266:AT:' (BRCA2 c.68_69del) |

??? note "上游说明（英文原文）"

    Convert a variant in SPDI notation to HGVS notation using NCBI Variation Services. SPDI (Sequence Position Deletion Insertion) is NCBI's canonical variant representation. Input format: 'SeqID:Position:DeletedSequence:InsertedSequence' (e.g., 'NC_000001.11:230710047:A:G'). Returns the equivalent HGVS genomic description (e.g., 'NC_000001.11:g.230710048A&gt;G'). Note: SPDI uses 0-based interbase coordinates, while HGVS uses 1-based.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NCBIVariation_spdi_to_hgvs',
     'arguments': {'spdi': '...'}},
)
```

## `NCBIVariation_spdi_to_rsids` {#NCBIVariation_spdi_to_rsids}

Returns dbSNP rsIDs co-located with a normalized SPDI variant.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`spdi`
- **关键词：**rsid lookup, spdi, dbsnp
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `spdi` | string | 是 | Variant in SPDI notation (SeqID:Position:Deleted:Inserted). Examples: 'NC_000019.10:44908683:T:C' (recovers rs429358 APOE), 'NC_000001.11:230710047:A:G' (AGT) |

??? note "上游说明（英文原文）"

    Reverse lookup from a normalized SPDI variant to the co-located dbSNP rsID(s), using NCBI Variation Services /spdi/\{spdi\}/rsids. Complements the forward direction (rsid_lookup) and the SPDI conversion tools, letting you recover an rsID from a genomic SPDI/coordinate. SPDI format: 'SeqID:Position:DeletedSequence:InsertedSequence' (0-based interbase). Returns a list of integer rsIDs.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NCBIVariation_spdi_to_rsids',
     'arguments': {'spdi': '...'}},
)
```

## `NCBIVariation_vcf_to_spdi` {#NCBIVariation_vcf_to_spdi}

Converts a VCF chrom/pos/ref/alt variant to a normalized contextual SPDI on RefSeq.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`chrom`, `pos`, `ref`, `alt`
- **关键词：**vcf, spdi, ref alt, normalization
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chrom` | string | 是 | Chromosome RefSeq accession (e.g., 'NC_000007.14' for chr7 GRCh38, 'NC_000019.10' for chr19) |
| `pos` | string | 是 | 1-based VCF position (e.g., '140753336' for BRAF V600E, '44908684' for APOE rs429358) |
| `ref` | string | 是 | Reference allele (e.g., 'A', 'T', 'G', 'C') |
| `alt` | string | 是 | Alternate allele (e.g., 'T', 'C'). Single allele; comma-separated alts also accepted by the API. |

??? note "上游说明（英文原文）"

    Convert a raw VCF four-field variant (chromosome RefSeq accession, position, reference allele, alternate allele) directly to a normalized contextual SPDI on a RefSeq sequence, using NCBI Variation Services /vcf/\{chrom\}/\{pos\}/\{ref\}/\{alts\}/contextuals. Lets you enter the SPDI pipeline from VCF coordinates without first building an HGVS string. Position is 1-based VCF coordinate; the returned SPDI position is 0-based interbase.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NCBIVariation_vcf_to_spdi',
     'arguments': {'chrom': '...',
                   'pos': '...',
                   'ref': '...',
                   'alt': '...'}},
)
```

## `OncoKB_annotate_copy_number` {#OncoKB_annotate_copy_number}

Annotates gene amplification or deletion for oncogenicity and therapies with OncoKB.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`gene`
- **关键词：**cna, amplification, deletion, oncogenicity
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type (fixed: annotate_copy_number) 默认 `annotate_copy_number` |
| `gene` | string | 是 | Gene symbol (e.g., ERBB2, MYC, CDKN2A) |
| `copy_number_type` | string | 否 | Type of copy number alteration (case-insensitive: AMPLIFICATION or DELETION) 可选值：`AMPLIFICATION`, `DELETION`, `Amplification`, `Deletion`, `amplification`, `deletion` |
| `tumor_type` | string | 否 | Optional OncoTree tumor type code |
| `copy_number_alteration` | string | 否 | Alias for copy_number_type. Amplification or Deletion. |

??? note "上游说明（英文原文）"

    Annotate copy number alterations (amplification/deletion) for oncogenic potential. Returns evidence for gene amplification or deletion in cancer. Example: ERBB2 amplification in breast cancer, CDKN2A deletion.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OncoKB_annotate_copy_number',
     'arguments': {'gene': '...'}},
)
```

## `OncoKB_annotate_mutations` {#OncoKB_annotate_mutations}

Annotates a cancer mutation such as BRAF V600E for oncogenicity, effect and therapies with OncoKB.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**cancer mutation, oncogenic, mutation effect, oncokb
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type (fixed: annotate_variant) |
| `query` | string | 否 | Mutation query string, e.g. 'BRAF V600E' or 'EGFR T790M'. Gene and variant separated by space. |
| `gene` | string | 否 | Gene symbol (e.g., BRAF, EGFR, TP53) |
| `variant` | string | 否 | Protein change (e.g., V600E, T790M) |
| `tumor_type` | string | 否 | Optional OncoTree tumor type code (e.g., MEL, LUAD) |

??? note "上游说明（英文原文）"

    Annotate a cancer mutation using OncoKB. Accepts query string like 'BRAF V600E' or separate gene/variant parameters. Returns oncogenicity, mutation effect, and FDA-approved therapies. Alias for OncoKB_annotate_variant with query-string parsing. Demo mode available for BRAF, TP53, ROS1.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OncoKB_annotate_mutations', 'arguments': {}},
)
```

## `OncoKB_annotate_variant` {#OncoKB_annotate_variant}

Annotates a cancer variant with evidence level, oncogenicity, mutation effect and FDA-approved therapies.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**oncokb, actionability, evidence level, variant annotation
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type (fixed: annotate_variant) 默认 `annotate_variant` |
| `gene` | string | 否 | Gene symbol (e.g., BRAF, EGFR, TP53, KRAS) |
| `gene_symbol` | string | 否 | Alias for gene. Gene symbol (e.g., BRAF, EGFR, TP53, KRAS). |
| `variant` | string | 否 | Variant notation - protein change (e.g., V600E, T790M, G12D) |
| `tumor_type` | string | 否 | Optional OncoTree tumor type code (e.g., MEL for melanoma, LUAD for lung adenocarcinoma) |
| `alteration` | string | 否 | Alias for variant. Protein change notation (e.g., "V600E", "T790M", "G12D"). |

??? note "上游说明（英文原文）"

    Annotate a cancer variant for oncogenic potential and treatment implications using OncoKB. Returns evidence level, oncogenicity, mutation effect, and FDA-approved therapies. Requires ONCOKB_API_TOKEN env var (demo mode available for BRAF, TP53, ROS1). Example: BRAF V600E in melanoma.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OncoKB_annotate_variant', 'arguments': {}},
)
```

## `OpenCRAVAT_annotate_variant` {#OpenCRAVAT_annotate_variant}

Annotates a single genomic variant with OpenCRAVAT sources such as ClinVar, gnomAD, REVEL, AlphaMissense and SpliceAI.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**opencravat, variant annotation, gnomad, alphamissense, spliceai
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chrom` | string \| null | 否 | Chromosome (e.g., 'chr17', 'chr7', '17'). 'chr' prefix added automatically if missing. |
| `chromosome` | string \| null | 否 | Alias for chrom. Chromosome (e.g., 'chr17', 'chr7', '17'). |
| `pos` | integer \| null | 否 | Genomic position (1-based, GRCh38 coordinates) |
| `position` | integer \| null | 否 | Alias for pos. Genomic position (1-based, GRCh38 coordinates). |
| `ref_base` | string \| null | 否 | Reference allele (e.g., 'C', 'A', 'G', 'T') |
| `ref` | string \| null | 否 | Alias for ref_base. Reference allele. |
| `alt_base` | string \| null | 否 | Alternate allele (e.g., 'T', 'A', 'G', 'C') |
| `alt` | string \| null | 否 | Alias for alt_base. Alternate allele. |
| `annotators` | string \| null | 否 | Comma-separated annotator names (e.g., 'clinvar,gnomad3,sift,polyphen2,revel,alphamissense'). If omitted, returns only basic variant annotation (gene, consequence). Common annotators: clinvar, gnomad3, sift, polyphen2, cadd_exome, revel, alphamissense, spliceai, dann, fathmm, gerp, phastcons, dbsnp, mutation_assessor, provean, vest, primateai, lrt. |

??? note "上游说明（英文原文）"

    Annotate a single genomic variant using OpenCRAVAT with 182+ annotation sources including ClinVar, gnomAD, SIFT, PolyPhen-2, REVEL, AlphaMissense, SpliceAI, GERP, PhastCons, DANN, FATHMM, and more. Returns gene name, amino acid change, consequence type, and all requested annotator results. Input: chromosome, position (1-based, GRCh38), reference and alternate alleles. Optionally specify which annotators to use (comma-separated).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenCRAVAT_annotate_variant', 'arguments': {}},
)
```

## `OpenTargets_drug_pharmacogenomics_data` {#OpenTargets_drug_pharmacogenomics_data}

Retrieves Open Targets pharmacogenomics records for a drug by ChEMBL ID, including genotype annotations and evidence levels.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`chemblId`
- **关键词：**pgx, genotype, drug response, chembl, evidence level
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chemblId` | string | 是 | The ChEMBL ID of the drug (e.g., 'CHEMBL25' for aspirin, 'CHEMBL1201583' for pembrolizumab). |

??? note "上游说明（英文原文）"

    Retrieve pharmacogenomics data for a specific drug, including evidence levels and genotype annotations.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_drug_pharmacogenomics_data',
     'arguments': {'chemblId': '...'}},
)
```

## `OpenTargets_get_target_constraint_info_by_ensemblID` {#OpenTargets_get_target_constraint_info_by_ensemblID}

Retrieves genetic constraint metrics (observed/expected, LOEUF-type scores) for a gene by Ensembl ID.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`ensemblId`
- **关键词：**gnomad, loeuf, loss of function, constraint
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ensemblId` | string | 是 | The Ensembl ID of the target. |

??? note "上游说明（英文原文）"

    Retrieve genetic constraint information for a specific target ensemblID, including expected and observed values, and scores.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_target_constraint_info_by_ensemblID',
     'arguments': {'ensemblId': '...'}},
)
```

## `OpenTargets_get_variant_effect_predictions` {#OpenTargets_get_variant_effect_predictions}

Retrieves variant effect predictions (VEP, SIFT, PolyPhen, AlphaMissense, CADD, GERP) for a chrom_pos_ref_alt variant.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variantId`
- **关键词：**sift, polyphen, alphamissense, cadd, vep
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variantId` | string | 是 | Variant ID in chrom_pos_ref_alt format (e.g., '19_44908822_C_T' for the APOE rs7412 variant). |

??? note "上游说明（英文原文）"

    Retrieve in-silico functional effect predictions for a variant (chrom_pos_ref_alt) from OpenTargets. Returns one record per prediction method (VEP, SIFT, PolyPhen, AlphaMissense, FoldX, GERP, CADD, etc.), each with an assessment (e.g. 'missense_variant', 'deleterious', 'tolerated') and a numeric score. Use to judge whether a variant is likely damaging.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_variant_effect_predictions',
     'arguments': {'variantId': '...'}},
)
```

## `OpenTargets_get_variant_info` {#OpenTargets_get_variant_info}

Returns variant details including coordinates, rsIDs, gnomAD allele frequencies, consequences and HGVS.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variantId`
- **关键词：**rsid, gnomad, allele frequency, hgvs
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variantId` | string | 是 | Open Targets variant ID in format chr_position_ref_alt (e.g., '10_112998590_C_T' for rs7903146, '19_44908684_T_C' for rs429358). Use the search tool to find variant IDs from rsIDs. |

??? note "上游说明（英文原文）"

    Get detailed information about a genetic variant from Open Targets Platform, including genomic coordinates, rsIDs, allele frequencies across populations (gnomAD), functional consequence, transcript consequences, and HGVS notation. Useful for GWAS variant annotation and population genetics.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_variant_info',
     'arguments': {'variantId': '...'}},
)
```

## `OpenTargets_get_variant_pharmacogenomics` {#OpenTargets_get_variant_pharmacogenomics}

Retrieves PharmGKB-derived pharmacogenomic annotations (drug, effect, evidence level, genotype) for a variant.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variantId`
- **关键词：**pgx, pharmgkb, genotype, drug response
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variantId` | string | 是 | Variant ID in chrom_pos_ref_alt format (e.g., '19_44908822_C_T' for the APOE rs7412 variant). |

??? note "上游说明（英文原文）"

    Retrieve variant-level pharmacogenomics (PGx) annotations for a variant (chrom_pos_ref_alt) from OpenTargets (PharmGKB-derived). Returns one record per PGx association with the drug(s), phenotypeText (clinical effect), evidenceLevel (PharmGKB 1A-4), genotypeAnnotationText, phenotype/PGx category, and source. Use to find drug-response associations for a variant.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_variant_pharmacogenomics',
     'arguments': {'variantId': '...'}},
)
```

## `OpenTargets_get_variant_transcript_consequences` {#OpenTargets_get_variant_transcript_consequences}

Returns per-transcript consequences of a variant: gene, amino acid change, impact and consequence score.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variantId`
- **关键词：**vep, transcript, impact, amino acid change
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variantId` | string | 是 | Variant ID in chrom_pos_ref_alt format (e.g., '19_44908822_C_T' for the APOE rs7412 variant). |

??? note "上游说明（英文原文）"

    Retrieve per-transcript consequences for a variant (chrom_pos_ref_alt) from OpenTargets. Returns one record per affected transcript with transcriptId, the target gene, aminoAcidChange (when coding), consequenceScore, impact (HIGH/MODERATE/LOW/MODIFIER), and whether it is the canonical/distance transcript. Use to see how a variant affects each transcript/isoform.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_variant_transcript_consequences',
     'arguments': {'variantId': '...'}},
)
```

## `Pathoplexus_get_mutations` {#Pathoplexus_get_mutations}

Returns characteristic nucleotide or amino-acid mutations for a Pathoplexus organism above a proportion threshold.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`organism`
- **关键词：**mutation prevalence, pathogen, lineage, genomic surveillance
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `organism` | string | 是 | Pathoplexus organism slug, e.g. 'west-nile', 'mpox'. |
| `country` | string \| null | 否 | Filter by country (geoLocCountry). |
| `lineage` | string \| null | 否 | Filter by lineage. |
| `min_proportion` | number \| null | 否 | Minimum proportion of sequences carrying the mutation (default 0.8). |
| `mutation_type` | string \| null | 否 | 'aminoAcid' (default) or 'nucleotide'. |
| `limit` | integer \| null | 否 | Max mutations (default 50, max 500). |

??? note "上游说明（英文原文）"

    Get characteristic amino-acid (default) or nucleotide mutations for a Pathoplexus organism above a proportion threshold, optionally filtered by country/lineage — for surveillance of mutation prevalence in emerging pathogens. Returns each mutation with gene, position, from/to residue, proportion (fraction of sequences carrying it), count, and coverage. Known organisms: 'west-nile', 'ebola-zaire', 'ebola-sudan', 'cchf', 'mpox'. No API key required.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Pathoplexus_get_mutations',
     'arguments': {'organism': '...'}},
)
```

## `PGSCatalog_get_performance_metrics` {#PGSCatalog_get_performance_metrics}

Returns published evaluation metrics (OR, HR, AUC, sample sets) for a PGS Catalog score.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`pgs_id`
- **关键词：**polygenic score, prs, auc, validation
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `pgs_id` | string | 是 | PGS Catalog score id, e.g. 'PGS000001' or 'PGS000004'. |

??? note "上游说明（英文原文）"

    Retrieve the published validation/evaluation performance metrics for a polygenic score from the PGS Catalog by its PGS id (e.g. 'PGS000001'). Returns one record per evaluation sample set, each with: effect sizes (Odds Ratio OR, Hazard Ratio HR, beta), classification accuracy (AUROC, Concordance/C-index), and other metrics (R2, top-vs-bottom quantile OR) — each as an estimate with 95% CI (ci_lower/ci_upper) — plus the evaluation cohort ancestry (ancestry_broad/free/country) and sample size (sample_number, cases, controls), covariates, and source publication. Use this to judge how well a score predicts and in which ancestries it was validated. Find PGS ids with PGSCatalog_get_scores_by_trait. No API key.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PGSCatalog_get_performance_metrics',
     'arguments': {'pgs_id': '...'}},
)
```

## `PGSCatalog_get_score` {#PGSCatalog_get_score}

Returns metadata of a PGS Catalog polygenic score: traits, variant count, method, build and ancestry.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`pgs_id`
- **关键词：**polygenic score, prs, genome build, ancestry
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `pgs_id` | string | 是 | PGS Catalog score id, e.g. 'PGS000001'. |

??? note "上游说明（英文原文）"

    Retrieve the full metadata for a specific polygenic score from the PGS Catalog by its PGS id (e.g. 'PGS000001'): reported and EFO/MONDO traits, number of variants, development method, genome build, weight type, publication, training samples, ancestry distribution, and the downloadable scoring-file URL. Find PGS ids with PGSCatalog_get_scores_by_trait. No API key.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PGSCatalog_get_score', 'arguments': {'pgs_id': '...'}},
)
```

## `PGSCatalog_get_scores_by_trait` {#PGSCatalog_get_scores_by_trait}

Lists published polygenic scores for an EFO/MONDO trait from the PGS Catalog.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`trait_id`
- **关键词：**polygenic score, prs, efo trait
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `trait_id` | string | 是 | EFO or MONDO trait id, e.g. 'MONDO_0004989' (breast carcinoma) or 'EFO_0001645' (coronary artery disease). |

??? note "上游说明（英文原文）"

    List the published polygenic scores (PGS) for a trait in the PGS Catalog, given its EFO/MONDO trait id (find one with PGSCatalog_search_traits). Each score summary includes the PGS id, name, reported trait, number of variants, development method, genome build, and publication (author/journal/year/PMID/DOI) plus the scoring-file URL. Use PGSCatalog_get_score for the full record of one score. No API key.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PGSCatalog_get_scores_by_trait',
     'arguments': {'trait_id': '...'}},
)
```

## `PGSCatalog_search_traits` {#PGSCatalog_search_traits}

Searches PGS Catalog traits by term and returns trait IDs with polygenic score counts.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`query`
- **关键词：**polygenic score, trait search, efo
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Trait or phenotype term, e.g. 'coronary artery disease', 'type 2 diabetes', 'breast cancer'. PGS Catalog's own search can return zero results for some multi-word phrases that don't exactly match an indexed label even when a near-identical trait exists (e.g. 'familial hypercholesterolemia' finds nothing while 'hypercholesterolemia' alone does) -- if a natural-language query returns no results, ret… |

??? note "上游说明（英文原文）"

    Search the PGS Catalog (pgscatalog.org, EMBL-EBI — the open database of published polygenic scores) for traits/phenotypes matching a term. Returns each matching trait's id (EFO/MONDO, e.g. 'MONDO_0004989'), label, description, and the number of polygenic scores published for it. Use this first to find the trait_id, then PGSCatalog_get_scores_by_trait to list its scores. No API key.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PGSCatalog_search_traits', 'arguments': {'query': '...'}},
)
```

## `PharmGKB_get_clinical_annotations` {#PharmGKB_get_clinical_annotations}

Returns PharmGKB clinical annotations with variant-level gene-drug-phenotype significance by annotation ID.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**pharmgkb, pgx, clinical annotation, genotype
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `annotation_id` | string | 否 | PharmGKB clinical annotation ID (e.g., '1447954390'). Required for reliable results. |
| `gene_id` | string | 否 | PharmGKB Gene Accession ID (e.g., "PA128"). NOTE: Gene-based lookup is not supported — will return an error with instructions to use annotation_id instead. |
| `gene` | string | 否 | NOT SUPPORTED: PharmGKB API requires a specific annotation_id (e.g. "1447954390"), not a gene symbol. Browse https://www.pharmgkb.org/clinicalAnnotation to find annotation IDs. For gene-drug dosing use CPIC_list_guidelines instead. |
| `gene_symbol` | string | 否 | Alias for gene — NOT SUPPORTED. See gene parameter. |

??? note "上游说明（英文原文）"

    Get clinical annotations showing gene-drug-phenotype relationships. Returns variant-level clinical significance data. Always provide annotation_id — gene_id lookup is not supported by the API and returns an error with guidance.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PharmGKB_get_clinical_annotations', 'arguments': {}},
)
```

## `PharmGKB_get_dosing_guidelines` {#PharmGKB_get_dosing_guidelines}

Returns CPIC/DPWG pharmacogenetic dosing guidelines from PharmGKB by guideline ID.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`guideline_id`
- **关键词：**cpic, dpwg, dosing, pharmacogene
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `guideline_id` | string | 是 | PharmGKB ClinPGx guideline ID from CPIC_list_guidelines 'clinpgxid' field (e.g., 'PA166251465' for warfarin, 'PA166251454' for opioids/codeine, 'PA166251458' for tamoxifen). Use clinpgxid, NOT pharmgkbid. |
| `gene` | string | 否 | Gene symbol (e.g., 'CYP2D6'). NOTE: Filtering by gene symbol is unreliable and may return a generic prompt instead of actual guidelines. Use guideline_id instead. |

??? note "上游说明（英文原文）"

    Get pharmacogenetic dosing guidelines (CPIC/DPWG) from PharmGKB by guideline_id. IMPORTANT: use the 'clinpgxid' field from CPIC_list_guidelines (e.g., 'PA166251454'), NOT the 'pharmgkbid' field — old-format pharmgkbid values (e.g., 'PA166104996') return HTTP 404. The test example PA166251465 is for the warfarin guideline.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PharmGKB_get_dosing_guidelines',
     'arguments': {'guideline_id': '...'}},
)
```

## `PharmGKB_get_variant_annotations` {#PharmGKB_get_variant_annotations}

Returns PharmGKB variant annotations (allele/genotype, phenotype, significance) filtered by gene or chemical.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**pgx, allele, genotype, pharmgkb
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_id` | string | 否 | PharmGKB Gene Accession ID (e.g., 'PA126' for CYP2C9). Queried via location.genes.accessionId. Provide either gene_id or chemical_id. |
| `chemical_id` | string | 否 | PharmGKB Chemical/drug Accession ID (e.g., 'PA451906' for warfarin). Queried via relatedChemicals.accessionId. Provide either gene_id or chemical_id. |
| `drug_id` | string | 否 | Alias for chemical_id. |
| `limit` | integer | 否 | Maximum number of variant annotations to return (default 50, max 500). CYP2C9 alone has over 1000 annotations. |

??? note "上游说明（英文原文）"

    Get variant-level literature annotations from PharmGKB/ClinPGx, filtered by gene or chemical PharmGKB Accession ID. Provide gene_id (e.g., 'PA126' for CYP2C9) OR chemical_id (e.g., 'PA451906' for warfarin). Each annotation gives allele/genotype, phenotype category, drug context, significance, and source literature. Results are capped by 'limit'. Keyless. Resolve gene/drug names to PA IDs first via PharmGKB_search_genes / PharmGKB_search_drugs.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PharmGKB_get_variant_annotations', 'arguments': {}},
)
```

## `PharmGKB_search_variants` {#PharmGKB_search_variants}

Searches PharmGKB genetic variants by rsID or name with gene info.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`query`
- **关键词：**rsid, pharmgkb, variant search
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Variant name or rsID (e.g., 'rs1799853'). |

??? note "上游说明（英文原文）"

    Search for genetic variants in PharmGKB by rsID or name. Returns variant IDs and associated gene info.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PharmGKB_search_variants', 'arguments': {'query': '...'}},
)
```

## `PharmVar_get_alleles` {#PharmVar_get_alleles}

Returns PharmVar star allele definitions, function, activity scores and variant details for a pharmacogene.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`gene`
- **关键词：**star allele, cyp2d6, hgvs, activity score
- **需要 API 密钥：**`PHARMVAR_API_KEY`
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 是 | Gene symbol (e.g., 'CYP2D6', 'CYP2C9', 'CYP2C19', 'CYP2B6', 'DPYD', 'SLCO1B1', 'UGT1A1', 'VKORC1') |

??? note "上游说明（英文原文）"

    Get all allele definitions for a pharmacogene from PharmVar (Pharmacogene Variation Consortium). Returns star allele names (e.g., CYP2D6*1, *4, *17), functional status, activity scores, variant details (HGVS notation, rsIDs, genomic positions), and evidence levels. PharmVar is the authoritative source for pharmacogenomic star allele nomenclature used in CPIC guidelines.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PharmVar_get_alleles', 'arguments': {'gene': '...'}},
)
```

## `PharmVar_list_genes` {#PharmVar_list_genes}

Lists pharmacogenes curated by PharmVar.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**pharmacogene, star allele, cyp
- **需要 API 密钥：**`PHARMVAR_API_KEY`
- **只读：**是（服务器声明）

??? note "上游说明（英文原文）"

    List all pharmacogenes curated by PharmVar (Pharmacogene Variation Consortium). PharmVar is the authoritative repository for pharmacogene allele nomenclature and variant data, defining star (*) alleles for genes like CYP2D6, CYP2C9, CYP2C19, CYP2B6, and DPYD. Returns 15+ genes with their allele counts and descriptions.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PharmVar_list_genes', 'arguments': {}},
)
```

## `Progenetix_cnv_search` {#Progenetix_cnv_search}

Searches cancer biosamples with CNVs (gains/losses) in a genomic region and disease filter.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`reference_name`, `start`, `end`
- **关键词：**cnv, copy number, genomic region, beacon
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `reference_name` | string | 是 | RefSeq chromosome accession. Examples: 'refseq:NC_000007.14' (chr7/GRCh38), 'refseq:NC_000017.11' (chr17/GRCh38), 'refseq:NC_000001.11' (chr1/GRCh38). |
| `start` | integer | 是 | Start position (1-based, GRCh38). Example: 55019017 for EGFR start. |
| `end` | integer | 是 | End position (1-based, GRCh38). Example: 55211628 for EGFR end. |
| `variant_type` | string | 否 | CNV type: 'DUP' for amplification/duplication, 'DEL' for deletion. Leave empty for any CNV type. 默认 `` |
| `filters` | string | 否 | Optional NCIt ontology code to filter by cancer type. Example: 'NCIT:C4017' for breast cancer. 默认 `` |
| `limit` | integer | 否 | Maximum number of biosamples to return (default: 10). 默认 `10` |

??? note "上游说明（英文原文）"

    Search for cancer biosamples with copy number variations (CNVs) in a specific genomic region. Combines disease filter with genomic coordinates to find samples with amplifications or deletions in genes of interest. Uses RefSeq chromosome accessions for coordinates. Example: searching for EGFR amplifications (refseq:NC_000007.14:55019017-55211628, DUP) in breast cancer (NCIT:C4017) returns 9 matching biosamples.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Progenetix_cnv_search',
     'arguments': {'reference_name': '...', 'start': 0, 'end': 0}},
)
```

## `Progenetix_get_cnv_frequencies` {#Progenetix_get_cnv_frequencies}

Returns genome-wide gain and loss frequencies per cytoband for a cancer type from Progenetix.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`filters`
- **关键词：**cnv frequency, cytoband, cancer
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `filters` | string | 是 | NCIt ontology code identifying the cancer-type collation. Examples: 'NCIT:C4017' (Breast Ductal Carcinoma), 'NCIT:C3058' (Glioblastoma), 'NCIT:C9145' (AML), 'NCIT:C3224' (Melanoma), 'NCIT:C4872' (Lung Adenocarcinoma). |
| `dataset_ids` | string | 否 | Progenetix dataset to query (default 'progenetix'). 默认 `progenetix` |
| `max_intervals` | integer | 否 | Optional cap on the number of interval bins returned (default: all ~3000 bins). |

??? note "上游说明（英文原文）"

    Get the genome-wide aggregate CNV frequency profile for a cancer type from Progenetix — the signature Progenetix output. For each genomic interval bin (cytoband), returns the gain frequency and loss frequency aggregated across all samples of a Progenetix collation, selected by NCIt code. Example: NCIT:C4017 (Breast Ductal Carcinoma) returns sampleCount 13869 across 3106 interval bins; NCIT:C3058 (Glioblastoma) bin 1p36.33 has gain_frequency 4.258 and loss_frequency 5.278.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Progenetix_get_cnv_frequencies',
     'arguments': {'filters': '...'}},
)
```

## `proteins_api_get_variants` {#proteins_api_get_variants}

Returns natural, disease-associated and population variants for UniProt accessions, paginated.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`accession`
- **关键词：**missense, disease variant, population frequency, uniprot
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `accession` | string \| array | 是 | UniProt protein accession(s). Can be a single accession (e.g., 'P05067'), comma-separated string (e.g., 'P05067,P04637'), or list (e.g., \['P05067', 'P04637'\]). Use UniProt search or ebi_search_domain with domain='uniprot' to find valid accessions. |
| `size` | integer | 否 | Maximum number of variants to return per protein (default: 100) 默认 `100` |
| `offset` | integer | 否 | Offset for pagination (default: 0) 默认 `0` |
| `format` | string | 否 | Response format. JSON is recommended for most use cases. 可选值：`json`, `xml` 默认 `json` |

??? note "上游说明（英文原文）"

    Get variation data for a protein including disease-associated variants, natural variants, and population data. Supports batch operations: pass a list of accessions or comma-separated string to retrieve multiple proteins in parallel. Returns variant information with positions, amino acid changes, disease associations, and evidence. Note: Not all proteins have annotated variants.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'proteins_api_get_variants',
     'arguments': {'accession': '...'}},
)
```

## `ProtVar_get_function` {#ProtVar_get_function}

Returns UniProt functional features and structural context at an amino-acid position, optionally for a variant residue.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`accession`, `position`
- **关键词：**missense, functional annotation, active site
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `accession` | string | 是 | UniProt accession (e.g. 'P04637' for TP53). |
| `position` | integer | 是 | Amino acid position in the protein (1-based). |
| `variant_aa` | string | 否 | Single-letter code for the variant amino acid (e.g. 'H' for histidine). Optional but recommended. |

??? note "上游说明（英文原文）"

    Get functional annotations at a specific amino acid position from ProtVar. Returns annotated UniProt features overlapping that position (domains, active sites, binding sites, PTMs), the protein's function description, and structural context. Use to understand the biological importance of a variant site before interpreting pathogenicity scores. Example: P04637 position 175 returns DNA-binding domain overlap and TP53 tumor suppressor function.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ProtVar_get_function',
     'arguments': {'accession': '...', 'position': 0}},
)
```

## `ProtVar_get_population` {#ProtVar_get_population}

Returns co-located variants with gnomAD/1000G frequencies, ClinVar significance and SIFT/PolyPhen at a protein position.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`accession`, `position`, `genomic_location`
- **关键词：**gnomad, clinvar, allele frequency, sift
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `accession` | string | 是 | UniProt accession (e.g. 'P22304' for IDS). |
| `position` | integer | 是 | Amino acid position in the protein (1-based). |
| `genomic_location` | integer | 是 | Genomic coordinate (GRCh38) for the variant. Obtain from ProtVar_map_variant output. |

??? note "上游说明（英文原文）"

    Get population observation data for a protein variant position from ProtVar. Returns co-located variants with population allele frequencies (gnomAD, 1000Genomes), clinical significance (ClinVar), and computational predictions (SIFT, PolyPhen).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ProtVar_get_population',
     'arguments': {'accession': '...',
                   'position': 0,
                   'genomic_location': 0}},
)
```

## `ProtVar_map_variant` {#ProtVar_map_variant}

Maps a protein change, rsID or VCF-style variant to genomic coordinates with AlphaMissense, EVE and ESM scores.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant`
- **关键词：**alphamissense, rsid, vcf, pathogenicity
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant` | string | 是 | Variant identifier in one of three formats: (1) 'UniProtAccession SingleLetterChange' e.g. 'P04637 R175H' (TP53 R175H); (2) dbSNP rsID e.g. 'rs1799966'; (3) VCF genomic 'chr17 43057065 . REF ALT'. Note: colon-separated 'chr:pos:ref:alt' is auto-converted to VCF format. |

??? note "上游说明（英文原文）"

    Map a human protein variant to genomic coordinates and get pathogenicity predictions (AlphaMissense, EVE, ESM, conservation). Accepts three input formats: (1) protein 'ACCESSION CHANGE' e.g. 'P04637 R175H'; (2) dbSNP rsID e.g. 'rs1799966'; (3) VCF-style genomic 'chr17 43057065 . T G'. Returns isoform mappings, consequence type, and variant effect scores.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ProtVar_map_variant', 'arguments': {'variant': '...'}},
)
```

## `StarAlleles_search` {#StarAlleles_search}

Autocompletes pharmacogenomic star alleles (e.g. CYP2D6*4) with nucleotide and protein changes.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`terms`
- **关键词：**star allele, pharmacogenomics, cyp2d6, haplotype
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `terms` | string | 是 | Gene or allele name/prefix to search, e.g. 'CYP2D6', 'CYP2D6*4'. |
| `max_results` | integer | 否 | Maximum number of matches to return (default 20, max 500). 默认 `20` |

??? note "上游说明（英文原文）"

    Pharmacogenomic star allele autocomplete via the NLM Clinical Table Search Service, e.g. CYP2D6*4. Given a gene or allele prefix, returns matching star allele designations with their defining nucleotide change, alternate/historical name, and resulting protein change. Complements the existing CPIC tools (which give dosing recommendations per diplotype but do not look up individual allele definitions). No API key required.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'StarAlleles_search', 'arguments': {'terms': '...'}},
)
```

## `UniProt_get_disease_variants_by_accession` {#UniProt_get_disease_variants_by_accession}

Extracts all VARIANT features and their disease annotations from a UniProtKB entry.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`accession`
- **关键词：**missense variants, disease variants, natural variants, uniprotkb
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `accession` | string | 是 | UniProtKB accession, e.g., P05067. |

??? note "上游说明（英文原文）"

    Extract all variants (feature type = VARIANT) and their related annotations from UniProtKB entry.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'UniProt_get_disease_variants_by_accession',
     'arguments': {'accession': '...'}},
)
```

## `VariantValidator_format_genomic_to_transcripts` {#VariantValidator_format_genomic_to_transcripts}

Projects a genomic HGVS variant onto every overlapping RefSeq transcript with c./p. descriptions via VariantFormatter.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant_description`
- **关键词：**hgvs, refseq transcripts, c. notation, variant projection
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `genome_build` | string | 否 | Reference genome assembly: 'GRCh37' (hg19) or 'GRCh38' (hg38). 默认 `GRCh38` |
| `variant_description` | string | 是 | Genomic-level variant description. Accepts genomic HGVS (e.g. 'NC_000017.11:g.50198002C&gt;A') or a pseudo-VCF string (e.g. '17-50198002-C-A'). |

??? note "上游说明（英文原文）"

    Project a genomic variant onto EVERY overlapping RefSeq transcript in one call using the VariantValidator VariantFormatter endpoint. Given a genomic-level variant (g.HGVS like 'NC_000017.11:g.50198002C&gt;A' or a pseudo-VCF string like '17-50198002-C-A') on GRCh37 or GRCh38, returns the per-transcript coding (c.) HGVS, predicted protein (p.) consequence in three- and single-letter codes, gene info (HGNC id + symbol), the pseudo-VCF representation, and MANE Select / RefSeq Select flags for each transcript. Use this when you need all transcript projections for a genomic variant; the standard VariantValidator_validate_variant endpoint no longer accepts select_transcripts='all' for genomic input (it returns HTTP 404 directing callers to VariantFormatter), so this is the only way to reach this capability in ToolUniverse. Public API, no key required.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'VariantValidator_format_genomic_to_transcripts',
     'arguments': {'variant_description': '...'}},
)
```

## `VariantValidator_validate_variant` {#VariantValidator_validate_variant}

Validates an HGVS variant description and returns genomic, transcript and protein-level representations.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`genome_build`, `variant_description`, `select_transcripts`
- **关键词：**hgvs, variant validation, nomenclature, liftover
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `genome_build` | string | 是 | Reference genome assembly: 'GRCh37' (hg19) or 'GRCh38' (hg38) |
| `variant_description` | string | 是 | HGVS variant description (e.g., 'NM_007294.4:c.5266dup' for BRCA1 c.5266dupC, 'NM_000179.3:c.943C&gt;T', 'chr17:g.43092919del') |
| `select_transcripts` | string | 是 | Transcript to validate against (e.g., 'NM_007294.4'). Use 'all' to get all transcripts for a genomic variant. |

??? note "上游说明（英文原文）"

    Validate and convert genetic variant descriptions using VariantValidator. Takes HGVS-formatted variant notation and returns validated variant descriptions with genomic context, transcript context, protein consequence, and cross-references. Validates variant nomenclature against current RefSeq sequences and LOVD syntax checker. Supports substitutions, deletions, insertions, duplications, inversions. Reference assemblies: GRCh37 (hg19) or GRCh38 (hg38).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'VariantValidator_validate_variant',
     'arguments': {'genome_build': '...',
                   'variant_description': '...',
                   'select_transcripts': '...'}},
)
```

## `VCF_count_variants` {#VCF_count_variants}

Counts records in a local VCF/BCF with bcftools after filters, with SNP/indel/MNP breakdown.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`vcf_path`
- **关键词：**vcf, bcftools, snp, indel, filtering
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | any | 否 | Operation (fixed) |
| `vcf_path` | string | 是 | Path to a local .vcf, .vcf.gz, or .bcf file |
| `pass_only` | boolean \| null | 否 | If true, count only FILTER=PASS or '.' records |
| `regions` | string \| null | 否 | Restrict to a region, e.g. 'chr1' or 'chr1:1-1000000' |
| `min_qual` | number \| null | 否 | Keep only records with QUAL &gt;= this value |
| `include_expr` | string \| null | 否 | Extra bcftools -i include expression, e.g. 'INFO/AF&gt;0.01' |

??? note "上游说明（英文原文）"

    Count variants in a local VCF/BCF after applying filters, via bcftools. Returns total records plus SNP/indel/MNP breakdown. Supports PASS-only, region, min QUAL, and a custom bcftools include expression so SNP/indel counts under a filter are deterministic. Requires bcftools on PATH.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'VCF_count_variants', 'arguments': {'vcf_path': '...'}},
)
```

## `VCF_normalize` {#VCF_normalize}

Normalizes a local VCF/BCF with bcftools norm and reports counts before and after.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`vcf_path`
- **关键词：**vcf, bcftools norm, left-align, multiallelic, indels
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | any | 否 | Operation (fixed) |
| `vcf_path` | string | 是 | Path to a local .vcf, .vcf.gz, or .bcf file |
| `multiallelics` | string \| null | 否 | 'split' (default) breaks multiallelic records into biallelic; 'join' merges; 'none' leaves as-is 可选值：`split`, `join`, `none`, `None` |
| `reference_fasta` | string \| null | 否 | Optional reference FASTA path; enables indel left-alignment (-f) |

??? note "上游说明（英文原文）"

    Normalize a local VCF/BCF with bcftools norm (split or join multiallelics; left-align indels against a reference FASTA when given) and report SNP/indel/record counts before vs after. Resolves the common cause of disagreeing indel counts: un-split multiallelics and non-left-aligned indels. Requires bcftools on PATH.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'VCF_normalize', 'arguments': {'vcf_path': '...'}},
)
```

## `VCF_summary_stats` {#VCF_summary_stats}

Computes bcftools stats for a local VCF/BCF: SNPs, indels, ts/tv ratio and per-sample genotype counts.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`vcf_path`
- **关键词：**vcf, bcftools stats, ts/tv, genotypes, heterozygosity
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | any | 否 | Operation (fixed) |
| `vcf_path` | string | 是 | Path to a local .vcf, .vcf.gz, or .bcf file |
| `pass_only` | boolean \| null | 否 | If true, count only FILTER=PASS or '.' records |
| `regions` | string \| null | 否 | Restrict to a region, e.g. 'chr1' or 'chr1:1-1000000' (streamed, no index needed) |
| `min_qual` | number \| null | 否 | Keep only records with QUAL &gt;= this value |
| `include_expr` | string \| null | 否 | Extra bcftools -i include expression, e.g. 'INFO/DP&gt;10' |

??? note "上游说明（英文原文）"

    Canonical variant statistics for a local VCF/BCF file via bcftools stats. Returns records, SNPs, indels, MNPs, multiallelic sites, ts/tv ratio, and per-sample het/hom/missing counts. Optional PASS-only, region, QUAL, and bcftools-expression filters. Use this instead of writing an ad-hoc parser so counts are reproducible. Requires bcftools on PATH.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'VCF_summary_stats', 'arguments': {'vcf_path': '...'}},
)
```

## `VEP_predict_pathogenicity` {#VEP_predict_pathogenicity}

Predicts missense variant pathogenicity with AlphaMissense, SIFT and PolyPhen-2 via Ensembl VEP for an rsID, HGVS or genomic variant.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**alphamissense, sift, polyphen, missense, rsid, hgvs
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `rsid` | string \| null | 否 | dbSNP rsID, e.g. 'rs699' (one input mode) |
| `hgvs_notation` | string \| null | 否 | HGVS notation, e.g. 'ENST00000269305.9:c.524G&gt;A' (one input mode) |
| `chrom` | string \| null | 否 | Chromosome for region input, e.g. '17' (with pos+alt) |
| `pos` | integer \| null | 否 | 1-based position for region input |
| `alt` | string \| null | 否 | Alternate allele for region input, e.g. 'T' |
| `genome_build` | string \| null | 否 | GRCh38 (default) or GRCh37 |

??? note "上游说明（英文原文）"

    Predict missense variant pathogenicity via the Ensembl VEP REST API with AlphaMissense (DeepMind 2023) + SIFT + PolyPhen-2 per transcript, in one call. Provide ONE of: a dbSNP rsID ('rsid'), HGVS notation ('hgvs_notation'), or a genomic locus ('chrom'+'pos'+'alt'). Returns per-transcript predictions (gene, amino_acids, alphamissense_class/score, sift, polyphen, most_severe_consequence). No API key.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'VEP_predict_pathogenicity', 'arguments': {}},
)
```

## `ZFIN_get_allele` {#ZFIN_get_allele}

Gets a zebrafish allele record from ZFIN via the Alliance API with associated gene, alteration type and synonyms.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`allele_id`
- **关键词：**zebrafish, danio rerio, mutant allele, alteration type
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `allele_id` | string | 是 | ZFIN allele ID with 'ZFIN:' prefix. Examples: 'ZFIN:ZDB-ALT-210706-1' (szu2 allele of pax2a), 'ZFIN:ZDB-ALT-980203-1590' (tm243a allele of pax2a). |

??? note "上游说明（英文原文）"

    Get detailed information about a specific zebrafish allele from ZFIN via the Alliance of Genome Resources. Returns allele symbol, associated gene, alteration type, species, and synonyms. Useful for looking up specific zebrafish mutant lines. Example: 'ZFIN:ZDB-ALT-210706-1' returns the szu2 allele of pax2a with its deletion variant details.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ZFIN_get_allele', 'arguments': {'allele_id': '...'}},
)
```

## `ZFIN_get_gene_alleles` {#ZFIN_get_gene_alleles}

Gets mutant alleles and variants of a zebrafish gene with variant types, molecular consequences and phenotypes.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`gene_id`
- **关键词：**zebrafish, mutant alleles, molecular consequence, insertion deletion
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_id` | string | 是 | ZFIN gene ID with 'ZFIN:' prefix. Examples: 'ZFIN:ZDB-GENE-990415-8' (pax2a), 'ZFIN:ZDB-GENE-990415-72' (tp53), 'ZFIN:ZDB-GENE-980526-166' (shha). |
| `limit` | integer | 否 | Maximum number of alleles to return (1-100). Default: 20. |
| `page` | integer | 否 | Page number for pagination. Default: 1. |

??? note "上游说明（英文原文）"

    Get alleles and variants for a zebrafish gene from ZFIN via the Alliance of Genome Resources. Returns mutant alleles with variant types (deletion, insertion, point mutation), molecular consequences, and phenotype/disease associations. Zebrafish genetics benefits from large-scale mutagenesis screens producing thousands of alleles. Example: 'ZFIN:ZDB-GENE-990415-8' (pax2a) returns alleles like szu2 (deletion) with start_lost consequence.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ZFIN_get_gene_alleles', 'arguments': {'gene_id': '...'}},
)
```
