# ToolUniverse: genome_annotation (`tooluniverse`)

DNAKit's MCP tool index keeps **845** DNA/gene tools of this server (2927 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

This page lists the 253 tools in category `genome_annotation`; other categories: [`sequence`](tooluniverse-sequence.md) (44), [`variant`](tooluniverse-variant.md) (155), [`association`](tooluniverse-association.md) (138), [`expression`](tooluniverse-expression.md) (56), [`regulation_epigenome`](tooluniverse-regulation_epigenome.md) (54), [`sequencing_data`](tooluniverse-sequencing_data.md) (91), [`prediction_model`](tooluniverse-prediction_model.md) (17), [`molecular_tools`](tooluniverse-molecular_tools.md) (23), [`breeding_population`](tooluniverse-breeding_population.md) (14).

| Item | Value |
| --- | --- |
| Source | [mims-harvard/ToolUniverse @ `9959ebe`](https://github.com/mims-harvard/ToolUniverse/tree/9959ebee7638773314362f1746980a15f5b17e5b) |
| License | Apache-2.0 |
| Command | `tooluniverse` |
| Transports | stdio |
| Installation | Install the pinned ToolUniverse distribution and required provider dependencies separately. Use the installed tooluniverse console script; unlike the README uvx --refresh example, this bridge does not refresh or install anything. |
| Effects | Starts the compact ToolUniverse MCP catalog. Selected tools may contact services, load models, execute programs or write files according to their provider definitions. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`AllenBrain_search_genes`](#AllenBrain_search_genes) | `genome_annotation` | Searches genes in the Allen Brain Atlas by symbol or name, returning gene metadata. |
| [`Alliance_get_gene`](#Alliance_get_gene) | `genome_annotation` | Retrieves an integrated gene record (symbol, name, synopsis, location) from the Alliance of Genome Resources. |
| [`Alliance_get_gene_orthologs`](#Alliance_get_gene_orthologs) | `genome_annotation` | Retrieves orthologs and paralogs of a gene across model organisms with prediction methods. |
| [`Alliance_get_molecular_interactions`](#Alliance_get_molecular_interactions) | `genome_annotation` | Retrieves physical molecular interactions for a gene with interactor genes and evidence from the Alliance. |
| [`Alliance_search_genes`](#Alliance_search_genes) | `genome_annotation` | Searches genes by symbol or synonym across human and model organism databases in the Alliance. |
| [`APPRIS_get_functional_annotations`](#APPRIS_get_functional_annotations) | `genome_annotation` | Retrieves APPRIS functional annotations (domains, residues, conservation) for a gene's isoforms. |
| [`APPRIS_get_isoforms`](#APPRIS_get_isoforms) | `genome_annotation` | Lists all transcript isoforms of a gene with their APPRIS classification. |
| [`APPRIS_get_principal_isoform`](#APPRIS_get_principal_isoform) | `genome_annotation` | Returns the APPRIS principal transcript isoform for an Ensembl gene. |
| [`BAR_get_gene_info`](#BAR_get_gene_info) | `genome_annotation` | Retrieves plant gene annotation (position, strand, aliases, function) from BAR. |
| [`Bgee_search_genes`](#Bgee_search_genes) | `genome_annotation` | Searches genes in Bgee by symbol, name or description across species. |
| [`BioGRID_get_interactions`](#BioGRID_get_interactions) | `genome_annotation` | Queries curated protein and genetic interactions for genes from BioGRID with evidence. |
| [`BVBRC_get_genome`](#BVBRC_get_genome) | `genome_annotation` | Retrieves BV-BRC genome metadata including taxonomy, size, GC content and assembly status. |
| [`BVBRC_search_genome_features`](#BVBRC_search_genome_features) | `genome_annotation` | Searches BV-BRC genome features such as genes and CDS across pathogen genomes. |
| [`BVBRC_search_genomes`](#BVBRC_search_genomes) | `genome_annotation` | Searches BV-BRC pathogen genome assemblies by organism or keyword. |
| [`BVBRC_search_pathways`](#BVBRC_search_pathways) | `genome_annotation` | Searches KEGG pathway assignments with EC numbers in BV-BRC pathogen genomes. |
| [`BVBRC_search_specialty_genes`](#BVBRC_search_specialty_genes) | `genome_annotation` | Searches BV-BRC specialty genes such as virulence factors and antibiotic resistance genes. |
| [`BVBRC_search_subsystems`](#BVBRC_search_subsystems) | `genome_annotation` | Searches BV-BRC functional subsystems and roles in pathogen genomes. |
| [`cBioPortal_get_gene_info`](#cBioPortal_get_gene_info) | `genome_annotation` | Returns the HUGO symbol and gene type for an Entrez Gene ID. |
| [`cBioPortal_get_gene_panel_genes`](#cBioPortal_get_gene_panel_genes) | `genome_annotation` | Lists all genes covered by a cBioPortal gene panel. |
| [`cBioPortal_get_genes`](#cBioPortal_get_genes) | `genome_annotation` | Searches cBioPortal genes by symbol or alias. |
| [`civic_search_genes`](#civic_search_genes) | `genome_annotation` | Searches CIViC genes by name. |
| [`CPIC_get_gene_info`](#CPIC_get_gene_info) | `genome_annotation` | Gets CPIC gene information including coordinates, sequence IDs and PharmGKB/HGNC/NCBI/Ensembl xrefs. |
| [`CPIC_list_pgx_genes`](#CPIC_list_pgx_genes) | `genome_annotation` | Lists all CPIC-curated pharmacogenomic genes with coordinates and cross-references. |
| [`DepMap_get_gene_dependencies`](#DepMap_get_gene_dependencies) | `genome_annotation` | Gets DepMap CRISPR knockout gene dependency (effect) scores for a gene, optionally per cell line. |
| [`DepMap_search_genes`](#DepMap_search_genes) | `genome_annotation` | Searches DepMap/Sanger genes by symbol or name, returning Ensembl IDs. |
| [`Dfam_get_annotations`](#Dfam_get_annotations) | `genome_annotation` | Gets Dfam transposable element hits for a chromosomal region of an assembly. |
| [`Dfam_get_family`](#Dfam_get_family) | `genome_annotation` | Gets a Dfam TE family by accession with consensus sequence, classification and clade distribution. |
| [`Dfam_search_families`](#Dfam_search_families) | `genome_annotation` | Searches Dfam TE families by name prefix, clade or repeat type. |
| [`DGIdb_get_gene_druggability`](#DGIdb_get_gene_druggability) | `genome_annotation` | Gets DGIdb druggability categories for genes. |
| [`DGIdb_get_gene_info`](#DGIdb_get_gene_info) | `genome_annotation` | Gets DGIdb gene information including aliases, categories and interaction counts. |
| [`EBIProteins_get_coordinate_mapping`](#EBIProteins_get_coordinate_mapping) | `genome_annotation` | Maps UniProt protein residue positions to chromosomal coordinates per exon for all isoforms. |
| [`EBIProteins_get_coordinates`](#EBIProteins_get_coordinates) | `genome_annotation` | Gets genomic coordinates for a UniProt protein including Ensembl gene/transcript IDs and exon structure. |
| [`EBIProteins_get_proteins_by_genomic_loc`](#EBIProteins_get_proteins_by_genomic_loc) | `genome_annotation` | Returns UniProt proteins whose coding sequence overlaps a genomic position, with Ensembl gene/transcript mapping. |
| [`EBIProteins_get_rna_editing`](#EBIProteins_get_rna_editing) | `genome_annotation` | Gets A-to-I RNA editing (recoding) sites for a UniProt protein with positions. |
| [`Enrichr_enrich`](#Enrichr_enrich) | `genome_annotation` | Runs Enrichr enrichment for a gene list against a chosen library, returning ranked terms and overlapping genes. |
| [`enrichr_gene_enrichment_analysis`](#enrichr_gene_enrichment_analysis) | `genome_annotation` | Runs Enrichr enrichment analysis on a gene list against selected pathway and ontology libraries. |
| [`Enrichr_gene_to_genesets`](#Enrichr_gene_to_genesets) | `genome_annotation` | Returns all Enrichr gene sets across libraries that contain a given gene. |
| [`Enrichr_get_top_enriched`](#Enrichr_get_top_enriched) | `genome_annotation` | Gets top enriched terms for a gene list across several Enrichr libraries at once. |
| [`ensembl_get_archive`](#ensembl_get_archive) | `genome_annotation` | Gets archived Ensembl gene information by stable ID across releases. |
| [`Ensembl_get_assembly_info`](#Ensembl_get_assembly_info) | `genome_annotation` | Gets Ensembl genome assembly metadata including karyotype and chromosome sizes. |
| [`Ensembl_get_cross_references`](#Ensembl_get_cross_references) | `genome_annotation` | Gets external cross-references (HGNC, Entrez, UniProt, OMIM) for an Ensembl stable ID. |
| [`Ensembl_get_gene_overlapping_features`](#Ensembl_get_gene_overlapping_features) | `genome_annotation` | Gets genomic features overlapping an Ensembl gene. |
| [`ensembl_get_genetree`](#ensembl_get_genetree) | `genome_annotation` | Gets an Ensembl gene tree for a gene family. |
| [`ensembl_get_homology`](#ensembl_get_homology) | `genome_annotation` | Gets Ensembl homologues of a gene symbol across species. |
| [`ensembl_get_overlap_features`](#ensembl_get_overlap_features) | `genome_annotation` | Gets genes, transcripts and other features overlapping a genomic region from Ensembl. |
| [`Ensembl_get_region_features`](#Ensembl_get_region_features) | `genome_annotation` | Gets genes, transcripts, regulatory and constrained elements in a region. |
| [`ensembl_get_species`](#ensembl_get_species) | `genome_annotation` | Gets Ensembl species information with assembly and database version. |
| [`Ensembl_get_species_info`](#Ensembl_get_species_info) | `genome_annotation` | Searches Ensembl species and their current genome assemblies. |
| [`ensembl_get_xrefs`](#ensembl_get_xrefs) | `genome_annotation` | Gets external database IDs for an Ensembl object. |
| [`ensembl_get_xrefs_by_name`](#ensembl_get_xrefs_by_name) | `genome_annotation` | Gets external database cross-references for a gene name or symbol. |
| [`ensembl_lookup_gene`](#ensembl_lookup_gene) | `genome_annotation` | Looks up an Ensembl gene by stable ID or symbol with location and description. |
| [`Ensembl_lookup_gene_by_symbol`](#Ensembl_lookup_gene_by_symbol) | `genome_annotation` | Looks up Ensembl gene and transcript IDs for a gene symbol. |
| [`EnsemblArchive_batch_lookup`](#EnsemblArchive_batch_lookup) | `genome_annotation` | Looks up version histories and current status for up to 50 Ensembl stable IDs. |
| [`EnsemblArchive_get_id_history`](#EnsemblArchive_get_id_history) | `genome_annotation` | Gets the version history and current status of an Ensembl stable ID. |
| [`EnsemblCompara_get_cafe_tree`](#EnsemblCompara_get_cafe_tree) | `genome_annotation` | Gets the Ensembl Compara CAFE tree with ancestral gene copy numbers for a gene family. |
| [`EnsemblCompara_get_gene_tree`](#EnsemblCompara_get_gene_tree) | `genome_annotation` | Gets the Ensembl Compara gene tree for a gene. |
| [`EnsemblCompara_get_orthologues`](#EnsemblCompara_get_orthologues) | `genome_annotation` | Finds orthologues of a gene using Ensembl Compara. |
| [`EnsemblCompara_get_paralogues`](#EnsemblCompara_get_paralogues) | `genome_annotation` | Finds paralogues of a gene using Ensembl Compara. |
| [`EnsemblMap_convert_coordinates`](#EnsemblMap_convert_coordinates) | `genome_annotation` | Converts genomic coordinates between assemblies such as GRCh37 and GRCh38. |
| [`EnsemblMap_translate_coordinates`](#EnsemblMap_translate_coordinates) | `genome_annotation` | Maps cDNA or protein positions of an Ensembl ID to chromosomal coordinates. |
| [`EnsemblReg_get_constrained_elements`](#EnsemblReg_get_constrained_elements) | `genome_annotation` | Gets constrained elements in a genomic region from Ensembl Compara. |
| [`EpiGraphDB_search_gene`](#EpiGraphDB_search_gene) | `genome_annotation` | Searches EpiGraphDB genes by symbol or Ensembl ID. |
| [`FlyBase_get_gene`](#FlyBase_get_gene) | `genome_annotation` | Gets Drosophila gene symbol, name, synopsis, genomic location, synonyms and cross-references from FlyBase via the Alliance. |
| [`FlyBase_get_gene_interactions`](#FlyBase_get_gene_interactions) | `genome_annotation` | Returns molecular (protein-protein) or genetic (suppression/enhancement) interactions for a Drosophila gene. |
| [`FlyBase_get_gene_orthologs`](#FlyBase_get_gene_orthologs) | `genome_annotation` | Lists orthologs of a Drosophila gene in human, mouse, zebrafish, worm, yeast and others with prediction methods. |
| [`FlyMine_search`](#FlyMine_search) | `genome_annotation` | Searches FlyMine across Drosophila genes, proteins, publications, protein domains and GO terms, with category facets. |
| [`FlyMine_search_genes`](#FlyMine_search_genes) | `genome_annotation` | Searches FlyMine for Drosophila genes and returns symbols, FBgn IDs, CG numbers and organism. |
| [`G2P_get_gene`](#G2P_get_gene) | `genome_annotation` | Gets Gene2Phenotype gene info: location, HGNC/Ensembl/OMIM cross-references, synonyms and coordinates. |
| [`GenomeNexus_get_canonical_transcript`](#GenomeNexus_get_canonical_transcript) | `genome_annotation` | Gets a gene's canonical Ensembl transcript, protein ID, RefSeq, CCDS and Pfam domains from Genome Nexus. |
| [`gnomad_get_gene`](#gnomad_get_gene) | `genome_annotation` | Gets gnomAD gene metadata (Ensembl ID, coordinates, canonical transcript) by symbol or Ensembl ID. |
| [`gnomad_get_transcript`](#gnomad_get_transcript) | `genome_annotation` | Gets gnomAD transcript metadata and coordinates with linked gene by Ensembl transcript ID. |
| [`gnomad_search_genes`](#gnomad_search_genes) | `genome_annotation` | Searches gnomAD genes by free text (usually symbol) and returns Ensembl IDs and symbols. |
| [`GO_get_annotations_for_gene`](#GO_get_annotations_for_gene) | `genome_annotation` | Returns a page of GO annotations for a gene/protein with evidence codes via GOlr. |
| [`GO_get_genes_for_term`](#GO_get_genes_for_term) | `genome_annotation` | Finds distinct genes annotated to a GO term and its descendants, optionally restricted to a taxon. |
| [`GOAPI_get_gene_functions`](#GOAPI_get_gene_functions) | `genome_annotation` | Gets GO biological process, molecular function and cellular component annotations for a gene with evidence types. |
| [`GOAPI_get_genes_by_function`](#GOAPI_get_genes_by_function) | `genome_annotation` | Gets genes/proteins annotated with a GO term, optionally filtered by taxon. |
| [`GoaT_get_species`](#GoaT_get_species) | `genome_annotation` | Looks up genome sequencing status, assembly level, genome size, chromosome number and ploidy via GoaT. |
| [`gProfiler_convert_ids`](#gProfiler_convert_ids) | `genome_annotation` | Converts gene identifiers (symbols, Ensembl, Entrez, UniProt, RefSeq) between namespaces with g:Convert. |
| [`gProfiler_enrichment`](#gProfiler_enrichment) | `genome_annotation` | Runs g:GOSt functional enrichment of a gene list against GO, KEGG, Reactome and other sources. |
| [`gProfiler_find_orthologs`](#gProfiler_find_orthologs) | `genome_annotation` | Maps genes from a source organism to orthologs in a target organism with g:Orth (Ensembl orthology). |
| [`GSEA_prerank`](#GSEA_prerank) | `genome_annotation` | Runs pre-ranked GSEA on a ranked gene list against gene sets returning ES, NES, p-values and leading-edge genes. |
| [`GTDB_get_genome`](#GTDB_get_genome) | `genome_annotation` | Gets GTDB genome metadata by accession: size, GC%, N50, CheckM completeness, gene counts and taxonomy. |
| [`GTDB_get_species`](#GTDB_get_species) | `genome_annotation` | Gets a GTDB species cluster with its assigned genome accessions, NCBI names and representative status. |
| [`GTDB_search_genomes`](#GTDB_search_genomes) | `genome_annotation` | Searches GTDB prokaryotic genomes by organism name and returns genome accessions and taxonomy. |
| [`Harmonizome_get_dataset`](#Harmonizome_get_dataset) | `genome_annotation` | Gets a Harmonizome dataset's description, association type, measurement and catalog of gene sets. |
| [`Harmonizome_get_gene`](#Harmonizome_get_gene) | `genome_annotation` | Gets Harmonizome gene info (symbol, name, Entrez ID, description) integrated from 100+ datasets. |
| [`Harmonizome_get_gene_set_members`](#Harmonizome_get_gene_set_members) | `genome_annotation` | Returns member genes of a Harmonizome attribute gene set with values, or a gene's associated attributes. |
| [`Harmonizome_list_datasets`](#Harmonizome_list_datasets) | `genome_annotation` | Lists 100+ integrated genomics datasets in Harmonizome (expression, pathways, regulation, disease). |
| [`Harmonizome_search`](#Harmonizome_search) | `genome_annotation` | Searches Harmonizome genes, datasets or attributes by keyword. |
| [`Harmonizome_search_genes`](#Harmonizome_search_genes) | `genome_annotation` | Searches Harmonizome genes by keyword returning matching gene symbols. |
| [`HGNC_fetch_gene_by_id`](#HGNC_fetch_gene_by_id) | `genome_annotation` | Fetches the HGNC gene record (symbol, name, location, cross-references) by HGNC ID. |
| [`HGNC_fetch_gene_by_symbol`](#HGNC_fetch_gene_by_symbol) | `genome_annotation` | Fetches HGNC gene data (ID, name, location, Ensembl/Entrez/UniProt/OMIM/RefSeq links) by symbol. |
| [`HGNC_fetch_gene_family_members`](#HGNC_fetch_gene_family_members) | `genome_annotation` | Lists all member genes of an HGNC gene group by gene_group_id with full records. |
| [`HGNC_search_by_location`](#HGNC_search_by_location) | `genome_annotation` | Searches HGNC for genes mapped to a cytogenetic location such as 17p13.1. |
| [`HGNC_search_genes`](#HGNC_search_genes) | `genome_annotation` | Searches HGNC human genes by symbol, name or alias with wildcard support. |
| [`HPA_generic_search`](#HPA_generic_search) | `genome_annotation` | Runs a custom Human Protein Atlas search query and returns selected columns. |
| [`HPA_get_biological_processes_by_gene`](#HPA_get_biological_processes_by_gene) | `genome_annotation` | Gets biological process annotations for a gene from HPA, highlighting key processes. |
| [`HPA_get_comprehensive_gene_details_by_ensembl_id`](#HPA_get_comprehensive_gene_details_by_ensembl_id) | `genome_annotation` | Gets full HPA gene page data by Ensembl ID: expression, antibodies, IHC and immunofluorescence images. |
| [`HPA_get_gene_basic_info_by_ensembl_id`](#HPA_get_gene_basic_info_by_ensembl_id) | `genome_annotation` | Gets basic gene info and expression data from the Human Protein Atlas by Ensembl ID. |
| [`HPA_get_gene_tsv_data_by_ensembl_id`](#HPA_get_gene_tsv_data_by_ensembl_id) | `genome_annotation` | Gets detailed HPA gene data in TSV format by Ensembl ID. |
| [`HPA_search_genes_by_query`](#HPA_search_genes_by_query) | `genome_annotation` | Searches HPA by gene name, keyword or cell line and returns Ensembl gene IDs. |
| [`humanbase_ppi_analysis`](#humanbase_ppi_analysis) | `genome_annotation` | Retrieves tissue-specific interaction networks and associated GO biological processes for a gene list from HumanBase. |
| [`HumanMine_search`](#HumanMine_search) | `genome_annotation` | Searches HumanMine (human, mouse, rat) across genes, proteins, pathways and diseases with category facets. |
| [`HumanMine_search_genes`](#HumanMine_search_genes) | `genome_annotation` | Searches HumanMine genes by symbol, name or ID for human, mouse and rat. |
| [`IMGT_search_genes`](#IMGT_search_genes) | `genome_annotation` | Searches IMGT immunoglobulin/TCR germline genes by type (IGHV, TRAV, etc.) and species. |
| [`IMPC_search_genes`](#IMPC_search_genes) | `genome_annotation` | Searches IMPC mouse genes by symbol, name, synonym or MGI ID with phenotyping status. |
| [`InterMine_run_pathquery`](#InterMine_run_pathquery) | `genome_annotation` | Runs an InterMine XML PathQuery against HumanMine and returns rows (e.g. gene pathways, domains, region features). |
| [`KEGG_convert_ids`](#KEGG_convert_ids) | `genome_annotation` | Converts KEGG gene or compound identifiers to UniProt, NCBI Gene, NCBI Protein, ChEBI or PubChem IDs. |
| [`kegg_find_genes`](#kegg_find_genes) | `genome_annotation` | Finds KEGG genes by keyword across all organisms or within one organism. |
| [`kegg_get_gene_info`](#kegg_get_gene_info) | `genome_annotation` | Returns KEGG gene information by gene ID including sequence, function and pathway associations. |
| [`KEGG_get_gene_pathways`](#KEGG_get_gene_pathways) | `genome_annotation` | Returns all KEGG pathways a gene participates in for a KEGG gene ID or symbol. |
| [`KEGG_get_network`](#KEGG_get_network) | `genome_annotation` | Returns a KEGG network entry with signaling definition, disease classes, drug targets and gene elements. |
| [`KEGG_get_pathway_genes`](#KEGG_get_pathway_genes) | `genome_annotation` | Returns gene IDs for all member genes of a KEGG pathway. |
| [`kegg_get_pathway_info`](#kegg_get_pathway_info) | `genome_annotation` | Returns KEGG pathway information by pathway ID including genes, compounds and reactions. |
| [`KEGG_link_entries`](#KEGG_link_entries) | `genome_annotation` | Finds linked entries between KEGG databases, e.g. diseases for a gene or genes in a pathway. |
| [`KEGG_search_network`](#KEGG_search_network) | `genome_annotation` | Searches KEGG NETWORK for disease/drug-associated signaling networks by gene or keyword. |
| [`LipidMaps_get_gene`](#LipidMaps_get_gene) | `genome_annotation` | Resolves a lipid-metabolism gene in LMPD to its record with names, chromosome, cytogenetic location and summary. |
| [`LipidMaps_get_protein`](#LipidMaps_get_protein) | `genome_annotation` | Resolves a lipid-metabolism protein in LMPD to gene symbol, species and UniProt/RefSeq cross-references. |
| [`LNCipedia_get_lncrna`](#LNCipedia_get_lncrna) | `genome_annotation` | Returns lncRNA details from RNAcentral by URS ID including RNA sequence, species, gene associations and sources. |
| [`LNCipedia_get_lncrna_xrefs`](#LNCipedia_get_lncrna_xrefs) | `genome_annotation` | Returns cross-references for a lncRNA from LNCipedia, NONCODE, HGNC, RefSeq, Ensembl and other databases. |
| [`LNCipedia_search_lncrna`](#LNCipedia_search_lncrna) | `genome_annotation` | Searches lncRNAs by name, accession or keyword over RNAcentral with optional species filter. |
| [`LNCipedia_search_ncrna_by_type`](#LNCipedia_search_ncrna_by_type) | `genome_annotation` | Searches RNAcentral non-coding RNAs by type such as miRNA, lncRNA, rRNA, tRNA or snoRNA, with optional species. |
| [`LOVD_get_gene`](#LOVD_get_gene) | `genome_annotation` | Returns LOVD gene information including HGNC/Entrez IDs, chromosomal location, RefSeq transcripts and build. |
| [`MARRVEL_get_gene`](#MARRVEL_get_gene) | `genome_annotation` | Returns MARRVEL aggregated gene identity for a human gene symbol with HGNC/Ensembl/Entrez IDs, location and summary. |
| [`MedlinePlus_get_genetics_gene_by_name`](#MedlinePlus_get_genetics_gene_by_name) | `genome_annotation` | Returns MedlinePlus Genetics information for a gene including function, related conditions and synonyms. |
| [`MetabolomicsWorkbench_get_gene_protein`](#MetabolomicsWorkbench_get_gene_protein) | `genome_annotation` | Looks up a Metabolomics Workbench gene or protein (MGP) annotation record by gene or protein identifier. |
| [`MGI_get_gene`](#MGI_get_gene) | `genome_annotation` | Returns MGI mouse gene details by MGI ID including synopsis, synonyms, GRCm39 location and cross-references. |
| [`MGI_search_genes`](#MGI_search_genes) | `genome_annotation` | Searches MGI mouse genes by symbol, name or keyword and returns MGI IDs. |
| [`MGnify_get_genome`](#MGnify_get_genome) | `genome_annotation` | Returns MGnify genome catalog details with completeness, contamination, N50, taxonomy and pangenome statistics. |
| [`MGnify_get_go_terms`](#MGnify_get_go_terms) | `genome_annotation` | Returns GO term annotations with gene counts from a MGnify metagenomics analysis. |
| [`MGnify_search_genomes`](#MGnify_search_genomes) | `genome_annotation` | Searches MGnify genome catalog for MAGs and isolate genomes by taxonomy lineage. |
| [`MIBiG_list_clusters`](#MIBiG_list_clusters) | `genome_annotation` | Lists MIBiG biosynthetic gene clusters matching a search term. |
| [`miRBase_get_mirna`](#miRBase_get_mirna) | `genome_annotation` | Returns miRNA details from RNAcentral including sequence, species, RNA type and gene associations. |
| [`miRBase_get_mirna_xrefs`](#miRBase_get_mirna_xrefs) | `genome_annotation` | Returns miRNA cross-references from miRBase, TarBase, ENA, RefSeq and other sources. |
| [`miRBase_search_mirna`](#miRBase_search_mirna) | `genome_annotation` | Searches miRNAs by name, accession or keyword over RNAcentral with species filter. |
| [`Monarch_search_gene`](#Monarch_search_gene) | `genome_annotation` | Searches Monarch for a gene by symbol or name and returns HGNC/NCBIGene CURIEs and orthologs. |
| [`MouseMine_search`](#MouseMine_search) | `genome_annotation` | Searches MouseMine (MGI) records across genes, alleles, proteins, GO terms and publications. |
| [`MouseMine_search_genes`](#MouseMine_search_genes) | `genome_annotation` | Searches MouseMine for mouse protein-coding genes and returns symbols and MGI IDs. |
| [`MSigDB_check_gene_in_set`](#MSigDB_check_gene_in_set) | `genome_annotation` | Checks whether a gene belongs to an MSigDB gene set such as TF target or hallmark sets. |
| [`MSigDB_get_gene_set_members`](#MSigDB_get_gene_set_members) | `genome_annotation` | Returns the parsed gene list of an MSigDB gene set by exact name. |
| [`MSigDB_get_geneset`](#MSigDB_get_geneset) | `genome_annotation` | Returns an MSigDB gene set by exact name from the API. |
| [`MSigDB_get_hallmark_geneset`](#MSigDB_get_hallmark_geneset) | `genome_annotation` | Returns an MSigDB Hallmark gene set by pathway name. |
| [`MyGene_batch_query`](#MyGene_batch_query) | `genome_annotation` | Queries MyGene.info for annotations of multiple gene IDs, supporting ID mapping. |
| [`MyGene_get_gene_annotation`](#MyGene_get_gene_annotation) | `genome_annotation` | Retrieves MyGene.info annotation for a gene by Entrez or Ensembl ID. |
| [`MyGene_query_genes`](#MyGene_query_genes) | `genome_annotation` | Searches MyGene.info for genes by keyword, symbol, name or identifier. |
| [`NCBIDatasets_get_gene`](#NCBIDatasets_get_gene) | `genome_annotation` | Returns NCBI gene information by Gene ID including symbol, chromosomal location, genomic coordinates and cross-references. |
| [`NCBIDatasets_get_gene_by_symbol`](#NCBIDatasets_get_gene_by_symbol) | `genome_annotation` | Looks up NCBI Gene records by official gene symbol and species (name or taxonomy ID). |
| [`NCBIDatasets_get_genome_assembly`](#NCBIDatasets_get_genome_assembly) | `genome_annotation` | Returns NCBI genome assembly metadata for a GCF_ or GCA_ accession including level, status and statistics. |
| [`NCBIDatasets_get_orthologs`](#NCBIDatasets_get_orthologs) | `genome_annotation` | Returns NCBI orthologs of a gene across species with symbols and gene types. |
| [`NCBIDatasets_get_sequence_reports`](#NCBIDatasets_get_sequence_reports) | `genome_annotation` | Returns per-sequence reports (chromosome, plasmid, scaffold) for a genome assembly with accessions, length and GC%. |
| [`NCBIDatasets_list_genomes_by_taxon`](#NCBIDatasets_list_genomes_by_taxon) | `genome_annotation` | Lists NCBI genome assemblies for a taxon with accession, level, N50, GC% and RefSeq category. |
| [`NCBIGene_get_summary`](#NCBIGene_get_summary) | `genome_annotation` | Returns an NCBI Gene summary by Entrez ID with function, chromosomal location, coordinates, exon count and MIM IDs. |
| [`NCBIGene_search`](#NCBIGene_search) | `genome_annotation` | Searches NCBI Gene by symbol, name or keyword with organism filter and returns Entrez Gene IDs. |
| [`OMA_get_genome_pair_orthologs`](#OMA_get_genome_pair_orthologs) | `genome_annotation` | Returns all pairwise orthologs between two whole genomes from OMA. |
| [`OMA_get_group`](#OMA_get_group) | `genome_annotation` | Returns an OMA Group of strict 1:1 orthologs with member proteins. |
| [`OMA_get_hog`](#OMA_get_hog) | `genome_annotation` | Returns OMA Hierarchical Orthologous Group information covering gene family evolution and duplications. |
| [`OMA_get_orthologs`](#OMA_get_orthologs) | `genome_annotation` | Returns OMA pairwise orthologs for a protein with relationship type and evolutionary distance. |
| [`OMA_get_protein`](#OMA_get_protein) | `genome_annotation` | Returns OMA protein entry details with species, genomic location and group membership. |
| [`OMA_get_protein_go`](#OMA_get_protein_go) | `genome_annotation` | Returns GO annotations for an OMA entry with aspect, evidence code and information content. |
| [`OMA_resolve_xref`](#OMA_resolve_xref) | `genome_annotation` | Resolves a gene name, symbol or cross-reference ID to OMA protein entries. |
| [`OmniPath_get_annotation_resource_geneset`](#OmniPath_get_annotation_resource_geneset) | `genome_annotation` | Returns the full gene/protein membership of an OmniPath annotation resource. |
| [`OpenGenes_get_gene`](#OpenGenes_get_gene) | `genome_annotation` | Returns the Open Genes aging/longevity profile of a human gene with mechanisms, clusters and disease categories. |
| [`OpenGenes_search_genes`](#OpenGenes_search_genes) | `genome_annotation` | Returns a page of Open Genes aging-gene summaries with mechanisms and confidence. |
| [`OpenTargets_get_target_cancer_hallmarks`](#OpenTargets_get_target_cancer_hallmarks) | `genome_annotation` | Retrieves COSMIC Cancer Gene Census hallmark annotations (promotes/suppresses) for a gene by Ensembl ID. |
| [`OpenTargets_get_target_classes_by_ensemblID`](#OpenTargets_get_target_classes_by_ensemblID) | `genome_annotation` | Returns target class annotations (e.g. enzyme, receptor families) for an Ensembl gene ID. |
| [`OpenTargets_get_target_depmap_essentiality`](#OpenTargets_get_target_depmap_essentiality) | `genome_annotation` | Retrieves DepMap CRISPR essentiality and Chronos gene effect scores per tissue/cell line for an Ensembl gene ID. |
| [`OpenTargets_get_target_gene_ontology_by_ensemblID`](#OpenTargets_get_target_gene_ontology_by_ensemblID) | `genome_annotation` | Retrieves Gene Ontology annotations for a target by Ensembl gene ID. |
| [`OpenTargets_get_target_genomic_location_by_ensemblID`](#OpenTargets_get_target_genomic_location_by_ensemblID) | `genome_annotation` | Returns chromosome, start, end and strand for a target Ensembl gene ID. |
| [`OpenTargets_get_target_homologues_by_ensemblID`](#OpenTargets_get_target_homologues_by_ensemblID) | `genome_annotation` | Retrieves homologues (orthologs/paralogs) of a target by Ensembl gene ID. |
| [`OpenTargets_get_target_id_description_by_name`](#OpenTargets_get_target_id_description_by_name) | `genome_annotation` | Looks up the Ensembl gene ID and description for a target name or symbol. |
| [`OpenTargets_get_target_info_by_ensemblID`](#OpenTargets_get_target_info_by_ensemblID) | `genome_annotation` | Returns core target information for an Ensembl gene ID: symbol, name, biotype, UniProt IDs, location and function. |
| [`OpenTargets_get_target_interactions_by_ensemblID`](#OpenTargets_get_target_interactions_by_ensemblID) | `genome_annotation` | Returns interaction partners and evidence for a target Ensembl gene ID, paginated. |
| [`OpenTargets_get_target_pathways_by_ensemblID`](#OpenTargets_get_target_pathways_by_ensemblID) | `genome_annotation` | Lists Reactome pathways a gene participates in by Ensembl ID, with top-level terms. |
| [`OpenTargets_get_target_synonyms_by_ensemblID`](#OpenTargets_get_target_synonyms_by_ensemblID) | `genome_annotation` | Retrieves synonyms and alternative symbols for a target Ensembl gene ID. |
| [`OrthoDB_get_group_details`](#OrthoDB_get_group_details) | `genome_annotation` | Returns KEGG, GO and InterPro annotations for an OrthoDB orthologous group. |
| [`OrthoDB_get_orthologs`](#OrthoDB_get_orthologs) | `genome_annotation` | Lists member genes (IDs, organisms, descriptions) of an OrthoDB orthologous group, optionally by species. |
| [`OrthoDB_search_groups`](#OrthoDB_search_groups) | `genome_annotation` | Searches OrthoDB v12 orthologous groups by gene/protein name and taxonomic level. |
| [`PANTHER_enrichment`](#PANTHER_enrichment) | `genome_annotation` | Runs PANTHER overrepresentation analysis on a gene list for GO terms, pathways or protein classes. |
| [`PANTHER_gene_info`](#PANTHER_gene_info) | `genome_annotation` | Returns PANTHER family/subfamily, GO annotations and protein class for a gene. |
| [`PANTHER_ortholog`](#PANTHER_ortholog) | `genome_annotation` | Finds orthologs of a gene across species using PANTHER, with ortholog type. |
| [`PathwayCommons_get_neighborhood`](#PathwayCommons_get_neighborhood) | `genome_annotation` | Returns known interactions (PPI, regulatory, catalytic) for a gene from Pathway Commons. |
| [`PathwayCommons_get_pathway`](#PathwayCommons_get_pathway) | `genome_annotation` | Returns a Pathway Commons pathway by URI with its name, participants and interaction processes. |
| [`PathwayCommons_paths_between`](#PathwayCommons_paths_between) | `genome_annotation` | Finds mechanistic interaction paths between two or more gene symbols across Pathway Commons sources. |
| [`PDC_get_gene_protein`](#PDC_get_gene_protein) | `genome_annotation` | Returns gene description, NCBI gene ID, protein accessions and PDC study coverage for a gene symbol. |
| [`PharmGKB_get_gene_details`](#PharmGKB_get_gene_details) | `genome_annotation` | Returns chromosomal location and cross-references for a PharmGKB gene accession. |
| [`PharmGKB_search_genes`](#PharmGKB_search_genes) | `genome_annotation` | Searches PharmGKB genes by name, symbol or ID. |
| [`Pharos_get_target`](#Pharos_get_target) | `genome_annotation` | Returns Pharos target details (TDL, family, novelty, diseases) by gene symbol or UniProt ID. |
| [`Planteome_search_annotations`](#Planteome_search_annotations) | `genome_annotation` | Returns GO, PO and TO annotations for a plant gene from Planteome. |
| [`PlantReactome_get_pathway_participants`](#PlantReactome_get_pathway_participants) | `genome_annotation` | Lists gene and protein participants with UniProt accessions and gene names for a Plant Reactome pathway. |
| [`PomBase_get_gene`](#PomBase_get_gene) | `genome_annotation` | Returns PomBase gene details: name, product, domains, deletion viability and UniProt cross-reference. |
| [`PomBase_get_go_annotations`](#PomBase_get_go_annotations) | `genome_annotation` | Returns GO annotations by aspect for an S. pombe gene from PomBase. |
| [`PomBase_get_interactions`](#PomBase_get_interactions) | `genome_annotation` | Returns physical and genetic interactions with evidence for an S. pombe gene. |
| [`PomBase_get_orthologs`](#PomBase_get_orthologs) | `genome_annotation` | Returns human and budding yeast orthologs of an S. pombe gene. |
| [`PomBase_search_genes`](#PomBase_search_genes) | `genome_annotation` | Searches PomBase S. pombe genes by name, systematic ID or product. |
| [`proteins_api_get_genome_mappings`](#proteins_api_get_genome_mappings) | `genome_annotation` | Returns genome mappings (Ensembl, RefSeq, EMBL cross-references and coordinates) for UniProt accessions. |
| [`proteins_api_get_xrefs`](#proteins_api_get_xrefs) | `genome_annotation` | Returns database cross-references (Ensembl, RefSeq, EMBL, PDB, Reactome) for UniProt accessions. |
| [`QuickGO_annotations_by_gene`](#QuickGO_annotations_by_gene) | `genome_annotation` | Returns GO annotations with evidence codes for a gene product from QuickGO. |
| [`QuickGO_annotations_by_goterm`](#QuickGO_annotations_by_goterm) | `genome_annotation` | Lists gene products annotated to a GO term with evidence, optionally by taxon. |
| [`Reactome_get_participant_reference_entities`](#Reactome_get_participant_reference_entities) | `genome_annotation` | Returns reference entities (UniProt proteins, ChEBI chemicals) for participants of a Reactome event. |
| [`Reactome_get_participants`](#Reactome_get_participants) | `genome_annotation` | Returns physical entities (proteins, complexes, small molecules) participating in a Reactome event. |
| [`Reactome_get_pathways_low_entity`](#Reactome_get_pathways_low_entity) | `genome_annotation` | Returns the most specific Reactome pathways containing a given entity. |
| [`Reactome_map_uniprot_to_pathways`](#Reactome_map_uniprot_to_pathways) | `genome_annotation` | Maps a UniProt accession to all Reactome pathways containing it. |
| [`Reactome_map_uniprot_to_reactions`](#Reactome_map_uniprot_to_reactions) | `genome_annotation` | Maps a UniProt accession to Reactome reactions involving it. |
| [`ReactomeAnalysis_pathway_enrichment`](#ReactomeAnalysis_pathway_enrichment) | `genome_annotation` | Runs Reactome overrepresentation analysis on gene symbols, UniProt or Ensembl IDs. |
| [`ReactomeAnalysis_pathway_found_entities`](#ReactomeAnalysis_pathway_found_entities) | `genome_annotation` | Lists submitted identifiers matching a hit pathway in a Reactome analysis, with cross-references. |
| [`ReactomeAnalysis_species_comparison`](#ReactomeAnalysis_species_comparison) | `genome_annotation` | Projects gene/protein identifiers from another species onto human Reactome pathways for enrichment. |
| [`ReactomeAnalysis_token_result`](#ReactomeAnalysis_token_result) | `genome_annotation` | Retrieves previously computed Reactome analysis results by token with pagination. |
| [`ReactomeInteractors_get_entity_pathways`](#ReactomeInteractors_get_entity_pathways) | `genome_annotation` | Returns Reactome pathways containing a given entity stable ID, flagging disease pathways. |
| [`resolve_identifier_for_gene_or_protein`](#resolve_identifier_for_gene_or_protein) | `genome_annotation` | Detects the namespace of a gene or protein identifier (symbol, Ensembl, UniProt, Entrez, RefSeq, HGNC) and maps it to others. |
| [`Rfam_accession_to_id`](#Rfam_accession_to_id) | `genome_annotation` | Converts an Rfam accession to its family ID/name. |
| [`Rfam_get_family`](#Rfam_get_family) | `genome_annotation` | Returns Rfam family metadata: RNA type, sequence counts, species coverage and search parameters. |
| [`Rfam_get_sequence_regions`](#Rfam_get_sequence_regions) | `genome_annotation` | Lists genomic regions (sequence IDs, coordinates, scores) containing Rfam family members. |
| [`Rfam_id_to_accession`](#Rfam_id_to_accession) | `genome_annotation` | Converts an Rfam family ID/name to its RF accession. |
| [`RGD_get_annotations`](#RGD_get_annotations) | `genome_annotation` | Returns RGD disease, phenotype, pathway and GO annotations for a rat gene. |
| [`RGD_get_gene`](#RGD_get_gene) | `genome_annotation` | Returns RGD rat gene symbol, name, description, type and Ensembl cross-reference. |
| [`RGD_get_orthologs`](#RGD_get_orthologs) | `genome_annotation` | Returns orthologs of a rat gene across species from RGD. |
| [`RGD_resolve_symbol_or_region`](#RGD_resolve_symbol_or_region) | `genome_annotation` | Resolves a rat gene symbol to its RGD record or lists genes overlapping a genomic region. |
| [`RGD_search_genes`](#RGD_search_genes) | `genome_annotation` | Searches RGD genes by symbol or keyword for rat, human or mouse. |
| [`RNAcentral_get_by_accession`](#RNAcentral_get_by_accession) | `genome_annotation` | Returns an RNAcentral ncRNA entry with annotations and source cross-references. |
| [`RNAcentral_get_genome_locations`](#RNAcentral_get_genome_locations) | `genome_annotation` | Returns chromosome, strand and coordinates of an RNAcentral ncRNA in a species. |
| [`RNAcentral_get_region_ncRNAs`](#RNAcentral_get_region_ncRNAs) | `genome_annotation` | Lists RNAcentral non-coding RNAs overlapping a genomic region for a species. |
| [`RNAcentral_get_xrefs_and_pubs`](#RNAcentral_get_xrefs_and_pubs) | `genome_annotation` | Returns per-source cross-references (NCBI gene, Ensembl/GENCODE transcripts) and publications for an RNAcentral entry. |
| [`RNAcentral_search`](#RNAcentral_search) | `genome_annotation` | Searches RNAcentral ncRNA records (miRNA, lncRNA, rRNA) by query. |
| [`SGD_get_gene`](#SGD_get_gene) | `genome_annotation` | Gets a yeast gene/locus overview from SGD by SGD ID, including names, description, aliases, UniProt ID and functional annotations. |
| [`SGD_get_go_annotations`](#SGD_get_go_annotations) | `genome_annotation` | Gets Gene Ontology annotations (process, function, component) with evidence codes for a yeast gene from SGD. |
| [`SGD_get_interactions`](#SGD_get_interactions) | `genome_annotation` | Gets genetic and physical interactions for a yeast gene from SGD, with experiment types, interacting genes and references. |
| [`SGD_search`](#SGD_search) | `genome_annotation` | Searches the Saccharomyces Genome Database across genes, GO terms, phenotypes, alleles, diseases and more, with category filtering. |
| [`STRING_functional_enrichment`](#STRING_functional_enrichment) | `genome_annotation` | Runs STRING enrichment of a gene/protein set against GO, KEGG and Reactome. |
| [`STRING_get_enrichment`](#STRING_get_enrichment) | `genome_annotation` | Runs STRING functional enrichment on a list of gene or protein identifiers and returns enriched GO, pathway and other terms with FDR. |
| [`STRING_get_functional_annotations`](#STRING_get_functional_annotations) | `genome_annotation` | Returns all STRING functional annotations (GO, KEGG, Reactome, domains, phenotypes, diseases, tissues) for a single gene or protein. |
| [`TargetMine_search_genes`](#TargetMine_search_genes) | `genome_annotation` | Searches TargetMine for genes and returns gene symbols, NCBI Gene IDs, names and organisms. |
| [`Tark_get_mane_transcripts`](#Tark_get_mane_transcripts) | `genome_annotation` | Looks up MANE Select and MANE Plus Clinical transcripts for a gene, pairing Ensembl ENST with RefSeq NM IDs. |
| [`Tark_get_transcript`](#Tark_get_transcript) | `genome_annotation` | Gets the archived Ensembl Tark record for a transcript ID including assembly, biotype, genomic location and sequence checksum. |
| [`TogoID_convert`](#TogoID_convert) | `genome_annotation` | Converts biological identifiers between databases (e.g. Ensembl gene to UniProt, NCBI Gene to HGNC) via TogoID. |
| [`UCSC_get_track`](#UCSC_get_track) | `genome_annotation` | Gets UCSC Genome Browser annotation track records (genes, regulatory elements, conservation, variants, repeats) for a genomic region. |
| [`UCSC_list_tracks`](#UCSC_list_tracks) | `genome_annotation` | Lists available UCSC Genome Browser annotation tracks for a genome assembly or shows a track's column schema. |
| [`UCSC_search`](#UCSC_search) | `genome_annotation` | Searches the UCSC Genome Browser for genes, transcripts or features by name and returns genomic positions. |
| [`UniProt_id_mapping`](#UniProt_id_mapping) | `genome_annotation` | Maps identifiers between databases (e.g. Ensembl or gene name to UniProt) with batch and async support. |
| [`UniProtIDMap_convert_ids`](#UniProtIDMap_convert_ids) | `genome_annotation` | Converts protein or gene identifiers between 100+ databases via the UniProt ID Mapping service. |
| [`UniProtIDMap_gene_to_uniprot`](#UniProtIDMap_gene_to_uniprot) | `genome_annotation` | Converts gene names or symbols to UniProt accessions for a given organism via UniProt ID Mapping. |
| [`VariantValidator_gene2transcripts`](#VariantValidator_gene2transcripts) | `genome_annotation` | Gets all transcripts for a gene with MANE annotations, Ensembl/RefSeq IDs and genomic coordinates via VariantValidator. |
| [`VEuPathDB_get_gene_record`](#VEuPathDB_get_gene_record) | `genome_annotation` | Retrieves a gene record from a VEuPathDB pathogen-genome site (e.g. PlasmoDB) by gene ID. |
| [`VEuPathDB_search_genes_by_organism`](#VEuPathDB_search_genes_by_organism) | `genome_annotation` | Searches a VEuPathDB site for all genes of a given organism or strain. |
| [`WikiPathways_find_pathways_by_gene`](#WikiPathways_find_pathways_by_gene) | `genome_annotation` | Finds all WikiPathways pathways containing a gene given as HGNC symbol, Entrez or Ensembl ID. |
| [`WikiPathways_get_pathway_genes`](#WikiPathways_get_pathway_genes) | `genome_annotation` | Gets the genes (gene products) participating in a WikiPathways pathway. |
| [`WormBase_get_gene`](#WormBase_get_gene) | `genome_annotation` | Gets a C. elegans gene overview from WormBase by gene ID, including names, description and type. |
| [`WormBase_get_interactions`](#WormBase_get_interactions) | `genome_annotation` | Gets physical and genetic interactions for a C. elegans gene from WormBase. |
| [`WormBase_get_orthologs`](#WormBase_get_orthologs) | `genome_annotation` | Gets orthologs, homologs and paralogs of a C. elegans gene across species from WormBase. |
| [`Xenbase_get_gene`](#Xenbase_get_gene) | `genome_annotation` | Gets Xenopus gene details from Xenbase via the Alliance of Genome Resources, including synopsis and genomic location. |
| [`Xenbase_search_genes`](#Xenbase_search_genes) | `genome_annotation` | Searches Xenopus tropicalis and laevis genes in Xenbase via the Alliance API. |
| [`ZFIN_get_gene`](#ZFIN_get_gene) | `genome_annotation` | Gets zebrafish gene details from ZFIN via the Alliance API, including synopsis, genomic location and cross-references. |
| [`ZFIN_get_gene_orthologs`](#ZFIN_get_gene_orthologs) | `genome_annotation` | Gets orthologs of a zebrafish gene in human, mouse, fly, worm, yeast and other species. |
| [`ZFIN_search`](#ZFIN_search) | `genome_annotation` | Searches ZFIN zebrafish genes by symbol or name and returns gene IDs. |

## `AllenBrain_search_genes` {#AllenBrain_search_genes}

Searches genes in the Allen Brain Atlas by symbol or name, returning gene metadata.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** gene symbol, entrez id, chromosome, homologene
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_acronym` | string \| null | no | Gene symbol/acronym for exact match. Examples: 'Gad1', 'Pvalb', 'Sst', 'Bdnf', 'Th'. |
| `gene_name` | string \| null | no | Gene name for partial match search (alternative to gene_acronym). Example: 'glutamate decarboxylase'. |
| `num_rows` | integer | no | Max results to return. Default: 50. Default `50` |

??? note "Upstream description"

    Search for genes in the Allen Brain Atlas by gene symbol (acronym) or name. Returns gene metadata including Entrez ID, chromosome, homologene ID. Useful for finding gene IDs needed to query expression data. Example: 'Gad1' returns glutamate decarboxylase 1 records for mouse and rat.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'AllenBrain_search_genes', 'arguments': {}},
)
```

## `Alliance_get_gene` {#Alliance_get_gene}

Retrieves an integrated gene record (symbol, name, synopsis, location) from the Alliance of Genome Resources.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** gene record, model organism, mgi, flybase, genomic location
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | Alliance gene ID with source prefix. Examples: 'HGNC:6081' (human INS), 'MGI:98834' (mouse Trp53), 'FB:FBgn0003996' (fly w), 'SGD:S000000259' (yeast GAL4), 'ZFIN:ZDB-GENE-990415-8' (zebrafish pax2a), 'WB:WBGene00006763' (worm unc-13). |

??? note "Upstream description"

    Get detailed gene information from the Alliance of Genome Resources, which integrates data from 7 model organism databases (SGD, FlyBase, WormBase, ZFIN, RGD, MGI, Xenbase) plus human. Returns gene symbol, name, species, synopsis, synonyms, genomic location, and cross-references. Example: 'HGNC:6081' returns insulin gene info.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Alliance_get_gene', 'arguments': {'gene_id': '...'}},
)
```

## `Alliance_get_gene_orthologs` {#Alliance_get_gene_orthologs}

Retrieves orthologs and paralogs of a gene across model organisms with prediction methods.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** orthologs, paralogs, homology, model organism
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | Alliance gene ID with member prefix, e.g. 'MGI:97490' (Pax6), 'HGNC:8620', 'RGD:3889', 'Xenbase:XB-GENE-865965'. |
| `stringency` | string \| null | no | Ortholog stringency filter: 'stringent' (default), 'moderate', or 'all'. |
| `limit` | integer \| null | no | Max results per category (default 20, max 100). |

??? note "Upstream description"

    Get orthologs and paralogs of a gene across model organisms via the Alliance of Genome Resources (used by MGI, Xenbase, and other members). Returns the stringent ortholog set with per-ortholog disease/expression annotation flags and prediction methods, plus the paralog set. Gene IDs use member prefixes (MGI:, HGNC:, RGD:, ZFIN:, FB:, WB:, SGD:, Xenbase:); bare WBGene/FBgn/MGI numeric IDs are auto-prefixed. Example: MGI:97490 (Pax6) returns gene_to_gene_orthology orthologs (e.g. Drosophila ey) and paralogs.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Alliance_get_gene_orthologs',
     'arguments': {'gene_id': '...'}},
)
```

## `Alliance_get_molecular_interactions` {#Alliance_get_molecular_interactions}

Retrieves physical molecular interactions for a gene with interactor genes and evidence from the Alliance.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** molecular interactions, interactor, gene network, evidence
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | Alliance gene ID with member prefix, e.g. 'MGI:97490' (Pax6), 'HGNC:8620'. |
| `limit` | integer \| null | no | Max results (default 20, max 100). |
| `page` | integer \| null | no | Page number (default 1). |

??? note "Upstream description"

    Get molecular (physical) interactions for a gene via the Alliance of Genome Resources (mouse/MGI, frog/Xenbase, and other members). Returns interactor gene IDs/symbols/species, interaction type, and supporting evidence references (PubMed / member curie). Gene IDs use member prefixes. Example: MGI:97490 (Pax6) returns gene_molecular_interaction records with evidence references such as MGI:6822438 (PMID:22863532).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Alliance_get_molecular_interactions',
     'arguments': {'gene_id': '...'}},
)
```

## `Alliance_search_genes` {#Alliance_search_genes}

Searches genes by symbol or synonym across human and model organism databases in the Alliance.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** gene search, symbol, synonym, model organism
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Gene symbol or synonym to search for (not a descriptive name). Examples: 'INS', 'TP53', 'shh', 'notch1'. |
| `limit` | integer | no | Maximum number of results to return (1-50). Default: 10. |

??? note "Upstream description"

    Search for genes across all model organisms in the Alliance of Genome Resources. Matches by gene symbol or synonym (type-ahead) and returns genes from human, mouse, rat, zebrafish, fly, worm, yeast, and frog databases. Use the official gene symbol rather than a descriptive name: 'INS' returns INS (human), Ins1/Ins2 (rat), ins (zebrafish), whereas 'insulin' matches GO terms/diseases, not genes.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Alliance_search_genes', 'arguments': {'query': '...'}},
)
```

## `APPRIS_get_functional_annotations` {#APPRIS_get_functional_annotations}

Retrieves APPRIS functional annotations (domains, residues, conservation) for a gene's isoforms.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** isoforms, transcript, functional residues, conservation
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | Ensembl gene ID (e.g., ENSG00000141510 for TP53). |
| `species` | string | no | Species name (e.g., homo_sapiens). Default `homo_sapiens` |
| `methods` | string | no | Comma-separated analysis methods to retrieve (e.g., 'firestar,spade,corsair'). Leave empty for all. |
| `transcript_id` | string | no | Optional Ensembl transcript ID to filter results for a specific isoform (e.g., ENST00000269305). |

??? note "Upstream description"

    Get detailed functional annotations for gene isoforms from APPRIS, including structural domains (spade), functional residues (firestar), 3D structure coverage (matador3d), cross-species conservation (corsair), transmembrane helices (thump), signal peptides (crash), and proteomics evidence (proteo).

    Annotations are grouped by type. Optionally filter by transcript_id to focus on a specific isoform.

    Available methods: appris, firestar, matador3d, spade, corsair, thump, crash, proteo, tsl

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'APPRIS_get_functional_annotations',
     'arguments': {'gene_id': '...'}},
)
```

## `APPRIS_get_isoforms` {#APPRIS_get_isoforms}

Lists all transcript isoforms of a gene with their APPRIS classification.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** transcript isoforms, ensembl gene, principal isoform
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | Ensembl gene ID (e.g., ENSG00000141510 for TP53, ENSG00000012048 for BRCA1). |
| `species` | string | no | Species name in lowercase with underscore (e.g., homo_sapiens, mus_musculus, danio_rerio). Default `homo_sapiens` |

??? note "Upstream description"

    Get all transcript isoforms and their APPRIS annotations for a gene. APPRIS classifies isoforms as PRINCIPAL (most functionally important) or ALTERNATIVE based on structural, functional, and evolutionary evidence.

    Requires an Ensembl gene ID (ENSG...). Use HGNC_fetch_gene_by_symbol or Ensembl tools to convert gene symbols to Ensembl IDs first.

    Reliability levels: PRINCIPAL:1 (highest) through PRINCIPAL:5, MINOR (alternative), or 'No Principal Isoform'.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'APPRIS_get_isoforms', 'arguments': {'gene_id': '...'}},
)
```

## `APPRIS_get_principal_isoform` {#APPRIS_get_principal_isoform}

Returns the APPRIS principal transcript isoform for an Ensembl gene.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** principal isoform, canonical transcript, ensembl
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | Ensembl gene ID (e.g., ENSG00000141510 for TP53). |
| `species` | string | no | Species name (e.g., homo_sapiens, mus_musculus). Default `homo_sapiens` |

??? note "Upstream description"

    Get the principal (most functionally important) transcript isoform for a gene as determined by the APPRIS database.

    APPRIS selects the principal isoform based on protein structure (matador3d, spade), functional residues (firestar), cross-species conservation (corsair), transmembrane topology (thump, crash), and proteomics evidence (proteo).

    Requires an Ensembl gene ID (ENSG...). Use HGNC_fetch_gene_by_symbol to convert gene symbols first.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'APPRIS_get_principal_isoform',
     'arguments': {'gene_id': '...'}},
)
```

## `BAR_get_gene_info` {#BAR_get_gene_info}

Retrieves plant gene annotation (position, strand, aliases, function) from BAR.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** plant gene, arabidopsis, agi locus, annotation
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | Gene locus identifier, e.g. 'AT1G01010' (Arabidopsis). |
| `species` | string \| null | no | Species, e.g. 'arabidopsis' (default). |

??? note "Upstream description"

    Get plant gene annotation from the BAR (Bio-Analytic Resource for Plant Biology, a Global Core Biodata Resource): genomic position, strand, aliases, and functional annotation. Example: gene_id='AT1G01010' (species='arabidopsis', the default) returns NAC domain protein 1's Chr1 coordinates and aliases ANAC001/NAC001/NTL10. Distinct from ToolUniverse's plant_reactome tools, which cover pathways, not per-gene annotation.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'BAR_get_gene_info', 'arguments': {'gene_id': '...'}},
)
```

## `Bgee_search_genes` {#Bgee_search_genes}

Searches genes in Bgee by symbol, name or description across species.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** gene search, ensembl id, species
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Gene name, symbol, or description to search for. Examples: 'TP53', 'insulin', 'BRCA1', 'sonic hedgehog'. |

??? note "Upstream description"

    Search for genes in the Bgee comparative gene expression database across 29+ animal species. Bgee integrates RNA-Seq, Affymetrix, EST, and in situ hybridization data to provide curated expression calls. Search by gene symbol, name, or description. Returns matching genes with species information, biotype, and expression summaries. Example: searching 'TP53' returns tumor protein p53 across human, mouse, zebrafish, and other model organisms.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Bgee_search_genes', 'arguments': {'query': '...'}},
)
```

## `BioGRID_get_interactions` {#BioGRID_get_interactions}

Queries curated protein and genetic interactions for genes from BioGRID with evidence.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** genetic interaction, protein interaction, biogrid, evidence
- **Needs API keys:** `BIOGRID_ACCESS_KEY`
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_names` | array of string | no | List of gene names or protein identifiers (e.g., \['TP53', 'BRCA1', 'MYC'\]). Accepts official gene symbols. |
| `gene_name` | string | no | Alias for gene_names: single gene symbol (e.g., 'TP53'). Converted to gene_names list automatically. |
| `organism` | string | no | Organism name (e.g., 'Homo sapiens', 'Mus musculus') or NCBI taxonomy ID (e.g., '9606' for human, '10090' for mouse). Default: 9606 (human) Default `9606` |
| `interaction_type` | string | no | Type of interaction: 'physical' (protein-protein), 'genetic' (epistasis, synthetic lethality), 'both' (all interactions) One of `physical`, `genetic`, `both` Default `both` |
| `evidence_types` | array of string | no | Filter by evidence types (e.g., \['Affinity Capture-MS', 'Two-hybrid'\] for physical, \['Synthetic Lethality', 'Dosage Rescue'\] for genetic). Leave empty for all evidence types. |
| `limit` | integer | no | Maximum number of interactions to return (default: 100, max: 10000) Default `100` |
| `throughput` | string | no | Filter by throughput: 'low' (low-throughput studies), 'high' (high-throughput screens), or leave empty for all One of `low`, `high`, `` |

??? note "Upstream description"

    Query experimentally validated protein and genetic interactions from BioGRID (Biological General Repository for Interaction Datasets). Returns curated interactions from published studies with evidence methods, PubMed citations, and throughput information. **Prerequisites**: Requires BIOGRID_API_KEY environment variable (free academic API key at https://webservice.thebiogrid.org/). BioGRID contains 2.3M+ interactions from 80+ organisms, all experimentally validated (no predictions). Use for: finding experimentally proven interactions, getting literature evidence for interactions, validating predicted interactions, finding interaction methods used, accessing high-confidence curated data.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'BioGRID_get_interactions', 'arguments': {}},
)
```

## `BVBRC_get_genome` {#BVBRC_get_genome}

Retrieves BV-BRC genome metadata including taxonomy, size, GC content and assembly status.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `genome_id`
- **Keywords:** pathogen genome, assembly, gc content, patric
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `genome_id` | string | yes | BV-BRC genome identifier (taxon_id.version format). Examples: '83332.12' (M. tuberculosis H37Rv), '511145.12' (E. coli K-12 MG1655), '93061.5' (S. aureus NCTC 8325). |

??? note "Upstream description"

    Get detailed genome information from BV-BRC (Bacterial and Viral Bioinformatics Resource Center, formerly PATRIC). Returns genome metadata including organism name, taxonomy, genome size, GC content, assembly status, host, disease, and isolation information. BV-BRC is the primary NIAID-funded resource for bacterial and viral pathogen genomics. Example: '83332.12' returns Mycobacterium tuberculosis H37Rv.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'BVBRC_get_genome', 'arguments': {'genome_id': '...'}},
)
```

## `BVBRC_search_genome_features` {#BVBRC_search_genome_features}

Searches BV-BRC genome features such as genes and CDS across pathogen genomes.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** genome features, cds, amr genes, virulence factor
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string \| null | no | Gene name to search for. Examples: 'mecA' (methicillin resistance), 'blaTEM' (beta-lactamase), 'katG' (isoniazid resistance). |
| `product` | string \| null | no | Gene product keyword to search. Examples: 'beta-lactamase', 'penicillin-binding protein', 'DNA gyrase'. |
| `genome_id` | string \| null | no | Restrict search to a specific genome. Example: '83332.12'. |
| `limit` | integer \| null | no | Maximum number of results. Default: 10. Max: 100. |
| `keyword` | string \| null | no | Alias for product: gene product keyword to search (e.g., 'spike protein', 'beta-lactamase'). |

??? note "Upstream description"

    Search for genome features (genes, proteins, CDS) in BV-BRC pathogen genomes. Find specific genes like antimicrobial resistance genes (mecA, blaTEM), virulence factors, or any annotated feature across thousands of pathogen genomes. Example: search for 'mecA' gene (methicillin resistance) across all genomes.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'BVBRC_search_genome_features', 'arguments': {}},
)
```

## `BVBRC_search_genomes` {#BVBRC_search_genomes}

Searches BV-BRC pathogen genome assemblies by organism or keyword.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `keyword`
- **Keywords:** genome search, bacterial genome, viral genome, assembly
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `keyword` | string | yes | Search keyword - organism name, strain, or pathogen. Examples: 'Staphylococcus aureus', 'Escherichia coli', 'Salmonella enterica', 'SARS-CoV-2'. |
| `limit` | integer \| null | no | Maximum number of results to return. Default: 10. Max: 100. |

??? note "Upstream description"

    Search for pathogen genomes in BV-BRC by organism name or keyword. Returns a list of matching genomes with metadata. Useful for finding bacterial or viral genome assemblies. Example: 'Staphylococcus aureus' returns all S. aureus genome assemblies.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'BVBRC_search_genomes', 'arguments': {'keyword': '...'}},
)
```

## `BVBRC_search_pathways` {#BVBRC_search_pathways}

Searches KEGG pathway assignments with EC numbers in BV-BRC pathogen genomes.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** kegg pathway, ec number, pathogen metabolism
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `taxon_id` | string \| null | no | NCBI Taxonomy ID. Examples: '1773' (M. tuberculosis), '562' (E. coli), '1280' (S. aureus). |
| `pathway_name` | string \| null | no | Pathway name keyword. Examples: 'Glycolysis', 'Fatty acid', 'Amino acid', 'Sulfur metabolism'. |
| `ec_number` | string \| null | no | Enzyme Commission number. Examples: '1.11.1.6' (catalase), '2.7.1.1' (hexokinase). |
| `genome_id` | string \| null | no | Restrict to a specific genome. Example: '83332.12' (M. tuberculosis H37Rv). |
| `limit` | integer \| null | no | Maximum number of results. Default: 25. Max: 100. |

??? note "Upstream description"

    Search for metabolic pathways in pathogen genomes from BV-BRC. Returns KEGG pathway assignments with EC numbers, enzyme descriptions, and genome associations. Useful for understanding pathogen metabolism and identifying potential drug targets. Example: search for pathways in M. tuberculosis (taxon_id 1773).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'BVBRC_search_pathways', 'arguments': {}},
)
```

## `BVBRC_search_specialty_genes` {#BVBRC_search_specialty_genes}

Searches BV-BRC specialty genes such as virulence factors and antibiotic resistance genes.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** virulence genes, resistance genes, card, vfdb
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string \| null | no | Gene name to search. Examples: 'mecA' (methicillin resistance), 'katG' (isoniazid resistance), 'drrC' (efflux pump). |
| `property` | string \| null | no | Specialty gene property category. Options: 'Antibiotic Resistance', 'Virulence Factor', 'Drug Target', 'Transporter', 'Human Homolog'. |
| `source` | string \| null | no | Source database for annotation. Options: 'CARD', 'NDARO', 'PATRIC_VF', 'VFDB', 'DrugBank', 'TTD', 'Human'. |
| `taxon_id` | string \| null | no | NCBI Taxonomy ID. Examples: '1280' (S. aureus), '1773' (M. tuberculosis), '562' (E. coli). |
| `limit` | integer \| null | no | Maximum number of results. Default: 25. Max: 100. |

??? note "Upstream description"

    Search for specialty genes in BV-BRC including virulence factors, antibiotic resistance genes, drug targets, and human homologs across pathogen genomes. Data sourced from CARD, NDARO, VFDB, PATRIC_VF, and other curated databases. Essential for identifying key pathogenicity and resistance determinants. Example: search for 'mecA' (methicillin resistance gene in Staphylococcus).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'BVBRC_search_specialty_genes', 'arguments': {}},
)
```

## `BVBRC_search_subsystems` {#BVBRC_search_subsystems}

Searches BV-BRC functional subsystems and roles in pathogen genomes.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** subsystems, functional roles, genome annotation
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `taxon_id` | string \| null | no | NCBI Taxonomy ID. Examples: '1773' (M. tuberculosis), '1280' (S. aureus), '562' (E. coli). |
| `superclass` | string \| null | no | Top-level functional category. Options include: 'METABOLISM', 'STRESS RESPONSE, DEFENSE, VIRULENCE', 'ENERGY', 'CELLULAR PROCESSES', 'DNA PROCESSING', 'PROTEIN PROCESSING', 'CELL ENVELOPE', 'RNA PROCESSING'. |
| `subsystem_name` | string \| null | no | Subsystem name keyword. Examples: 'mycothiol', 'iron acquisition', 'beta-lactam resistance'. |
| `role_name` | string \| null | no | Specific functional role keyword. Examples: 'catalase', 'penicillin-binding protein'. |
| `genome_id` | string \| null | no | Restrict to a specific genome. Example: '83332.12'. |
| `limit` | integer \| null | no | Maximum number of results. Default: 25. Max: 100. |

??? note "Upstream description"

    Search for functional subsystems (curated groups of functionally related proteins) in BV-BRC pathogen genomes. Returns subsystem classifications including superclass, class, subclass, and specific functional roles. Useful for understanding pathogen functional organization and comparing metabolic capabilities. Example: search for virulence-related subsystems in M. tuberculosis.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'BVBRC_search_subsystems', 'arguments': {}},
)
```

## `cBioPortal_get_gene_info` {#cBioPortal_get_gene_info}

Returns the HUGO symbol and gene type for an Entrez Gene ID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `entrez_gene_id`
- **Keywords:** entrez id, hugo symbol, gene type
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `entrez_gene_id` | integer | yes | Entrez Gene ID (e.g., 672 for BRCA1) |

??? note "Upstream description"

    Get basic information about a specific gene by Entrez Gene ID. Returns the HUGO gene symbol and gene type (e.g., 'protein-coding'); does not return aliases -- use HGNC_fetch_gene_by_symbol or NCBIGene_get_summary for gene aliases/chromosome location.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'cBioPortal_get_gene_info',
     'arguments': {'entrez_gene_id': 0}},
)
```

## `cBioPortal_get_gene_panel_genes` {#cBioPortal_get_gene_panel_genes}

Lists all genes covered by a cBioPortal gene panel.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_panel_id`
- **Keywords:** gene panel, targeted sequencing, panel genes
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_panel_id` | string | yes | Gene panel ID (e.g., 'IMPACT468') |

??? note "Upstream description"

    Get all genes in a specific gene panel. Essential for understanding what genes are covered when analyzing mutation frequency data from targeted sequencing studies.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'cBioPortal_get_gene_panel_genes',
     'arguments': {'gene_panel_id': '...'}},
)
```

## `cBioPortal_get_genes` {#cBioPortal_get_genes}

Searches cBioPortal genes by symbol or alias.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** gene search, alias, symbol
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `keyword` | string | no | Gene symbol or alias to search for (e.g., 'BRCA1', 'TP53') |
| `query` | string | no | Alias for keyword: gene symbol or alias to search for |

??? note "Upstream description"

    Search for genes by keyword or query (gene symbol or alias). Alias: cBioPortal_search.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'cBioPortal_get_genes', 'arguments': {}},
)
```

## `civic_search_genes` {#civic_search_genes}

Searches CIViC genes by name.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** gene search, civic, cancer genes
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `limit` | integer | no | Maximum number of genes to return (default: 10, recommended max: 100) Default `10` |
| `name` | string | no | Gene symbol to search for (e.g., "EGFR", "BRAF", "BRCA1"). Alias: use 'query' or 'gene_name' instead. |
| `query` | string | no | Gene symbol to search for (e.g., "FLT3", "KRAS", "TP53"). Alias for 'name'. |
| `gene_name` | string | no | Gene symbol to search for. Alias for 'name'. |

??? note "Upstream description"

    Search for genes in CIViC (Clinical Interpretation of Variants in Cancer) database. CIViC is a community knowledgebase for expert-curated interpretations of variants in cancer. Returns genes with their IDs, names, and descriptions. Use this to find genes of interest before querying variants or evidence.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'civic_search_genes', 'arguments': {}},
)
```

## `CPIC_get_gene_info` {#CPIC_get_gene_info}

Gets CPIC gene information including coordinates, sequence IDs and PharmGKB/HGNC/NCBI/Ensembl xrefs.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** pharmacogene, cpic, gene coordinates, hgnc, pharmgkb
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `symbol` | string \| null | no | Gene symbol (e.g., 'CYP2D6', 'CYP2C19', 'SLCO1B1', 'TPMT', 'DPYD') |
| `gene` | string \| null | no | Gene symbol alias -- alternative to symbol. |
| `gene_symbol` | string \| null | no | Gene symbol alias -- alternative to symbol. |

??? note "Upstream description"

    Get pharmacogenomics information for a gene from the Clinical Pharmacogenomics Implementation Consortium (CPIC) database. Returns gene coordinates, sequence IDs, and cross-references to PharmGKB, HGNC, NCBI, and Ensembl. Use this to find which genes are pharmacogenomically relevant. Example genes: CYP2D6, CYP2C19, CYP2C9, SLCO1B1, TPMT, DPYD, UGT1A1, VKORC1, G6PD.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CPIC_get_gene_info', 'arguments': {}},
)
```

## `CPIC_list_pgx_genes` {#CPIC_list_pgx_genes}

Lists all CPIC-curated pharmacogenomic genes with coordinates and cross-references.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** pharmacogenes, cpic, cyp enzymes, transporters
- **Read-only:** yes (declared by the server)

??? note "Upstream description"

    List all pharmacogenomics genes curated by the Clinical Pharmacogenomics Implementation Consortium (CPIC). Returns 130+ genes with known drug interactions including CYP enzymes, transporters, and other pharmacogenes. Includes gene coordinates, sequence IDs, and database cross-references.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CPIC_list_pgx_genes', 'arguments': {}},
)
```

## `DepMap_get_gene_dependencies` {#DepMap_get_gene_dependencies}

Gets DepMap CRISPR knockout gene dependency (effect) scores for a gene, optionally per cell line.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_symbol`
- **Keywords:** gene essentiality, crispr screen, dependency, gene effect
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_symbol` | string | yes | Gene symbol (e.g., 'EGFR', 'KRAS', 'TP53') |
| `model_id` | string | no | Optional: Filter by specific cell line |

??? note "Upstream description"

    Get gene dependency information from CRISPR knockout screens. Returns gene IDs and annotations. Negative gene effect scores indicate essentiality (cell death upon knockout).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DepMap_get_gene_dependencies',
     'arguments': {'gene_symbol': '...'}},
)
```

## `DepMap_search_genes` {#DepMap_search_genes}

Searches DepMap/Sanger genes by symbol or name, returning Ensembl IDs.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** gene symbol, ensembl id, depmap, gene search
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Gene symbol or name to search (e.g., 'WDR7', 'EGFR', 'TP53') |

??? note "Upstream description"

    Search for genes in DepMap/Sanger by symbol or name. Returns matching genes with Ensembl IDs and exact-match flags. Use this helper to find the correct gene symbol before querying dependencies, especially when DepMap_get_gene_dependencies returns mismatched genes.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DepMap_search_genes', 'arguments': {'query': '...'}},
)
```

## `Dfam_get_annotations` {#Dfam_get_annotations}

Gets Dfam transposable element hits for a chromosomal region of an assembly.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `chrom`, `start`, `end`
- **Keywords:** transposable element, repeat annotation, dfam, genomic region
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `assembly` | string | no | Genome assembly name (default 'hg38'). Options: 'hg38', 'hg19', 'mm10', 'mm39', 'dm6', etc. |
| `chrom` | string | yes | Chromosome name (e.g., 'chr1', 'chr17', 'chrX'). |
| `start` | integer | yes | Start position (0-based). |
| `end` | integer | yes | End position. |
| `nrph` | boolean | no | Non-redundant profile hits only (default true). Set false for all hits. |

??? note "Upstream description"

    Get transposable element (TE) annotations for a genomic region from Dfam. Returns all TE/repeat element hits in the specified chromosomal region including family name, type (SINE, LINE, DNA, LTR, Satellite), strand, bit score, E-value, and alignment coordinates. Useful for analyzing repeat content of a genomic locus, identifying TEs near genes, or studying TE distribution. Example: chr1:10000-50000 in hg38 returns hits including TAR1 satellites, L1MC4a LINEs, MER5B DNA transposons, MIR1 SINEs.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Dfam_get_annotations',
     'arguments': {'chrom': '...', 'start': 0, 'end': 0}},
)
```

## `Dfam_get_family` {#Dfam_get_family}

Gets a Dfam TE family by accession with consensus sequence, classification and clade distribution.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `accession`
- **Keywords:** transposable element, consensus sequence, repeat family, dfam
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | Dfam family accession (e.g., 'DF000000003' for AluSc, 'DF000000053' for AluYa5). |

??? note "Upstream description"

    Get detailed information for a specific Dfam transposable element family by accession. Returns consensus sequence, taxonomic classification, repeat type/subtype, clade distribution, author, citations, and curation state. Useful for studying TE biology, understanding repeat element structure, or getting consensus sequences for bioinformatics analysis. Example: DF000000003 returns AluSc (309bp SINE/Alu, consensus sequence, found in primates, published Jurka 1988).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Dfam_get_family', 'arguments': {'accession': '...'}},
)
```

## `Dfam_search_families` {#Dfam_search_families}

Searches Dfam TE families by name prefix, clade or repeat type.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** transposable element, alu, line, ltr, repeat family
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `name_prefix` | string | no | TE family name prefix to search (e.g., 'AluS', 'L1M', 'SVA', 'MER'). Case-sensitive. |
| `clade` | string | no | Taxonomic clade filter as NCBI taxon ID (e.g., '9606' for human, '40674' for Mammalia, '7227' for Drosophila). |
| `repeat_type` | string | no | Repeat type filter. Options: 'SINE', 'LINE', 'LTR', 'DNA', 'Satellite', 'RC', 'Other'. |
| `limit` | integer | no | Maximum number of results (default 20, max 50). |

??? note "Upstream description"

    Search Dfam for transposable element (TE) / repeat element families by name prefix, taxonomic clade, or repeat type. Dfam is the authoritative database of repetitive DNA, including SINEs (Alu), LINEs (L1), DNA transposons, LTR retrotransposons, and satellites. Returns family accession, name, length, classification, and repeat type. Useful for identifying repeat families in a genome or studying transposon evolution. Example: searching name_prefix='AluS' returns Alu SINE subfamilies (AluSc, AluSx, etc.) found in primates; clade=9606 with repeat_type='SINE' returns human-specific SINEs.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Dfam_search_families', 'arguments': {}},
)
```

## `DGIdb_get_gene_druggability` {#DGIdb_get_gene_druggability}

Gets DGIdb druggability categories for genes.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** druggable genome, gene category, kinase, gpcr
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `genes` | array of string | no | List of gene symbols to check druggability. Aliases: gene_name, gene. |
| `gene_name` | string | no | Alias for genes. Single gene symbol (e.g., 'EGFR'). |
| `gene` | string | no | Alias for genes. Single gene symbol (e.g., 'EGFR'). |

??? note "Upstream description"

    Get druggability information for genes. Returns gene categories indicating if a gene is druggable (kinase, GPCR, ion channel, etc.). Gene symbols DGIdb does not recognize are omitted from the results and named in metadata.unmatched_genes.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DGIdb_get_gene_druggability', 'arguments': {}},
)
```

## `DGIdb_get_gene_info` {#DGIdb_get_gene_info}

Gets DGIdb gene information including aliases, categories and interaction counts.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** dgidb, gene aliases, gene category
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `genes` | array of string | no | List of gene symbols. Aliases: gene_name, gene. |
| `gene_name` | string | no | Alias for genes. Single gene symbol (e.g., 'EGFR'). |
| `gene` | string | no | Alias for genes. Single gene symbol (e.g., 'EGFR'). |

??? note "Upstream description"

    Get detailed gene information from DGIdb including aliases, categories, and interaction counts. Gene symbols DGIdb does not recognize are omitted from the results and named in metadata.unmatched_genes.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DGIdb_get_gene_info', 'arguments': {}},
)
```

## `EBIProteins_get_coordinate_mapping` {#EBIProteins_get_coordinate_mapping}

Maps UniProt protein residue positions to chromosomal coordinates per exon for all isoforms.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `accession`
- **Keywords:** protein-to-genome, exon mapping, chromosome, isoforms
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | UniProt accession. Examples: 'P04637' (TP53), 'P00533' (EGFR), 'P38398' (BRCA1), 'P01308' (Insulin), 'P69905' (HBA1 hemoglobin alpha). |

??? note "Upstream description"

    Map UniProt protein positions to genomic coordinates at exon-level resolution. Returns chromosome, strand, and per-exon mappings between protein residue positions and chromosomal coordinates for all transcript isoforms. Essential for translating protein variant positions (e.g., 'p.R273H at position 273') to genomic coordinates for variant analysis, connecting protein annotations to genome browsers, or understanding exon-intron structure at the protein level. Example: 'P04637' (TP53) returns 6 transcript mappings on chromosome 17 (reverse strand) with 10 exons spanning positions 7669612-7676594; 'P00533' (EGFR) maps to chromosome 7.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EBIProteins_get_coordinate_mapping',
     'arguments': {'accession': '...'}},
)
```

## `EBIProteins_get_coordinates` {#EBIProteins_get_coordinates}

Gets genomic coordinates for a UniProt protein including Ensembl gene/transcript IDs and exon structure.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `accession`
- **Keywords:** genomic coordinates, exon structure, ensembl gene, strand
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | UniProt accession. Examples: 'P04637' (TP53), 'P00533' (EGFR), 'P38398' (BRCA1). |

??? note "Upstream description"

    Get genomic coordinate mappings for a protein from the EBI Proteins API. Maps UniProt protein to genomic coordinates including Ensembl gene/transcript IDs, chromosome location, strand, and exon structure. Essential for cross-referencing protein features with genomic data, variant mapping between coordinate systems, and genome browser integration. Example: P04637 (TP53) maps to ENSG00000141510 on chr17.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EBIProteins_get_coordinates',
     'arguments': {'accession': '...'}},
)
```

## `EBIProteins_get_proteins_by_genomic_loc` {#EBIProteins_get_proteins_by_genomic_loc}

Returns UniProt proteins whose coding sequence overlaps a genomic position, with Ensembl gene/transcript mapping.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** genome-to-protein, genomic position, cds, ensembl transcript
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `taxonomy` | string \| integer | no | NCBI taxonomy ID. Default '9606' (human). Example: 9606 (Homo sapiens), 10090 (mouse). |
| `location` | string \| null | no | Genomic location as 'chromosome:position' (1-based). Example: '17:7676154'. Provide this OR both 'chromosome' and 'position'. |
| `chromosome` | string \| integer \| null | no | Chromosome (used with 'position' if 'location' not given). Example: '17', 'X'. |
| `position` | string \| integer \| null | no | 1-based genomic position (used with 'chromosome' if 'location' not given). Example: 7676154. |

??? note "Upstream description"

    Reverse genome-to-protein coordinate lookup from the EBI Proteins API: given a genomic position (taxonomy + chromosome:position), return all UniProt proteins whose coding sequence overlaps that position, with exon-level mapping (accession, Ensembl gene/transcript/translation IDs, protein start/end residue, amino acid, chromosome, strand, assembly). Inverse of EBIProteins_get_coordinate_mapping. Use to find which proteins a chromosomal variant affects and at which residue. Example: taxonomy 9606, location '17:7676154' returns 18 entries including TP53 (P04637) at protein residue 72 (Pro) on ENST00000923569.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EBIProteins_get_proteins_by_genomic_loc',
     'arguments': {}},
)
```

## `EBIProteins_get_rna_editing` {#EBIProteins_get_rna_editing}

Gets A-to-I RNA editing (recoding) sites for a UniProt protein with positions.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `accession`
- **Keywords:** rna editing, a-to-i, recoding, adar
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | UniProt accession. Examples: 'P42262' (GRIA2, canonical Q/R editing site), 'P28335' (HTR2C, 5 editing sites). |

??? note "Upstream description"

    Get RNA editing annotations (A-to-I / Q-&gt;R recoding sites) for a UniProt protein from the EBI Proteins API. Returns post-transcriptional RNA editing sites where a genomic adenosine is read as guanosine, recoding the protein. Each site includes the protein position, genomic location, codon, consequence type, wild-type and edited amino acid, and tissue/sample editing statistics. Distinct from DNA-level variation. Example: P42262 (GRIA2) returns 3 sites including the canonical Q/R site at residue 607 (Q-&gt;R, genomic NC_000004.12:g.157336723A&gt;G, missense). P28335 (HTR2C) returns 5 sites.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EBIProteins_get_rna_editing',
     'arguments': {'accession': '...'}},
)
```

## `Enrichr_enrich` {#Enrichr_enrich}

Runs Enrichr enrichment for a gene list against a chosen library, returning ranked terms and overlapping genes.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_list`
- **Keywords:** enrichr, gene set enrichment, pathway, go, kegg
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type (fixed: enrich) Default `enrich` |
| `gene_list` | array of string | yes | List of gene symbols (e.g., \['TP53', 'BRCA1', 'EGFR'\]). Minimum 1 gene. |
| `library` | string | no | Enrichr library name (default: GO_Biological_Process_2023). Common: GO_Biological_Process_2023, KEGG_2021_Human, Reactome_2022, WikiPathways_2024_Human, MSigDB_Hallmark_2020 Default `GO_Biological_Process_2023` |
| `top_n` | integer | no | Number of top enriched terms to return (default: 10) Default `10` |
| `species` | string \| null | no | Organism instance to query: 'human' (default), 'mouse', 'fly' (FlyEnrichr), 'worm' (WormEnrichr), 'yeast' (YeastEnrichr), or 'fish' (FishEnrichr). Gene set libraries differ per organism, and gene symbols must match that organism's nomenclature. One of `None`, `human`, `mouse`, `fly`, `worm`, `yeast`, `fish` |

??? note "Upstream description"

    Perform gene set enrichment analysis using Enrichr. Submit a list of genes and get enriched terms from a specific library. Returns ranked terms with p-values, z-scores, combined scores, and overlapping genes. No authentication required. Set species to query the fly, worm, yeast, or fish Enrichr instances; libraries and gene symbols differ per organism.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Enrichr_enrich', 'arguments': {'gene_list': []}},
)
```

## `enrichr_gene_enrichment_analysis` {#enrichr_gene_enrichment_analysis}

Runs Enrichr enrichment analysis on a gene list against selected pathway and ontology libraries.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_list`
- **Keywords:** enrichment, gene set, pathway, go terms, over-representation
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_list` | array of string | yes | List of gene names or symbols to analyze. At least 2 genes are required for path ranking analysis. |
| `libs` | array of string | no | List of Enrichr enrichment libraries to query (e.g. 'KEGG_2021_Human', 'GO_Biological_Process_2023'). Omitting this or passing an empty list uses a broad default set of pathway/ontology libraries. Default `\['WikiPathways_2024_Human', 'Reactome_Pathways_2024', 'MSig…` |

??? note "Upstream description"

    Perform gene enrichment analysis using Enrichr to find biological pathways, processes, and molecular functions associated with a list of genes. Returns connectivity paths between genes and enrichment terms.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'enrichr_gene_enrichment_analysis',
     'arguments': {'gene_list': []}},
)
```

## `Enrichr_gene_to_genesets` {#Enrichr_gene_to_genesets}

Returns all Enrichr gene sets across libraries that contain a given gene.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene`
- **Keywords:** enrichr, gene set membership, pathways, tf targets
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type (fixed: gene_to_genesets) Default `gene_to_genesets` |
| `gene` | string | yes | A single official HGNC gene symbol (e.g., 'STAT3', 'BRCA1'). Case-sensitive; use the canonical symbol. |
| `include_metadata` | boolean | no | If true, also return library category/description metadata (Enrichr setup=true). Default: false. Default `False` |
| `max_terms_per_library` | integer | no | Optional cap on the number of term names returned per library (0 = no cap, return all). Default: 0. Default `0` |
| `species` | string \| null | no | Organism instance to query: 'human' (default), 'mouse', 'fly' (FlyEnrichr), 'worm' (WormEnrichr), 'yeast' (YeastEnrichr), or 'fish' (FishEnrichr). Gene set libraries differ per organism, and gene symbols must match that organism's nomenclature. One of `None`, `human`, `mouse`, `fly`, `worm`, `yeast`, `fish` |

??? note "Upstream description"

    Reverse gene-to-geneset membership lookup. Given ONE gene symbol, return every Enrichr term/gene-set (by name) that contains it, broken down across all ~243 Enrichr libraries (KEGG/Reactome/WikiPathways pathways, GO terms, transcription-factor target sets, disease/drug signatures, GeneSigDB publication tables, etc.). This is the inverse of forward enrichment (which takes a gene list): use it to ask 'what known programs/signatures is gene X a member of?' without running an enrichment. Returns a dict of library_name -&gt; list of term/gene-set names containing the gene. Set include_metadata=true to also return per-library category/description metadata. No authentication required. Set species to query the fly, worm, yeast, or fish Enrichr instances; libraries and gene symbols differ per organism.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Enrichr_gene_to_genesets', 'arguments': {'gene': '...'}},
)
```

## `Enrichr_get_top_enriched` {#Enrichr_get_top_enriched}

Gets top enriched terms for a gene list across several Enrichr libraries at once.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_list`
- **Keywords:** enrichr, functional annotation, go, reactome, kegg
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type (fixed: get_top_enriched) Default `get_top_enriched` |
| `gene_list` | array of string | yes | List of gene symbols (e.g., \['TP53', 'BRCA1', 'EGFR'\]) |
| `libraries` | array of string | no | Libraries to query (default: GO_Biological_Process_2023, KEGG_2021_Human, Reactome_2022, WikiPathways_2024_Human) |
| `top_n` | integer | no | Number of top terms per library (default: 5) Default `5` |
| `species` | string \| null | no | Organism instance to query: 'human' (default), 'mouse', 'fly' (FlyEnrichr), 'worm' (WormEnrichr), 'yeast' (YeastEnrichr), or 'fish' (FishEnrichr). Gene set libraries differ per organism, and gene symbols must match that organism's nomenclature. One of `None`, `human`, `mouse`, `fly`, `worm`, `yeast`, `fish` |

??? note "Upstream description"

    Get top enriched terms across multiple libraries simultaneously. Submit a gene list and get the most significant pathways, processes, and functions from GO, KEGG, Reactome, and WikiPathways in one call. Ideal for quick functional annotation overview. Set species to query the fly, worm, yeast, or fish Enrichr instances; libraries and gene symbols differ per organism.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Enrichr_get_top_enriched',
     'arguments': {'gene_list': []}},
)
```

## `ensembl_get_archive` {#ensembl_get_archive}

Gets archived Ensembl gene information by stable ID across releases.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `id`
- **Keywords:** ensembl archive, stable id, gene history
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `id` | string | yes | Ensembl stable ID (gene, transcript, or protein, e.g., 'ENSG00000139618') |

??? note "Upstream description"

    Get historical (archived) gene information by stable ID. Returns gene location and status across all Ensembl releases. Use this to track gene ID changes, assembly updates, or retrieve information from older Ensembl versions.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ensembl_get_archive', 'arguments': {'id': '...'}},
)
```

## `Ensembl_get_assembly_info` {#Ensembl_get_assembly_info}

Gets Ensembl genome assembly metadata including karyotype and chromosome sizes.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `species`
- **Keywords:** genome assembly, karyotype, chromosome length, genebuild
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `species` | string | yes | Ensembl species name (lowercase with underscores). Examples: 'homo_sapiens' (human), 'mus_musculus' (mouse), 'danio_rerio' (zebrafish), 'drosophila_melanogaster' (fruit fly), 'saccharomyces_cerevisiae' (yeast). |

??? note "Upstream description"

    Get genome assembly metadata for a species from the Ensembl REST API. Returns assembly name, accession, version history, chromosome karyotype with sizes, genebuild date, golden path length, and available coordinate system versions. Essential for interpreting genomic coordinates and understanding assembly composition. Example: 'homo_sapiens' returns GRCh38.p14, assembly date 2013-12, genebuild May 2025, 194 top-level regions, karyotype chromosomes 1-22 plus X and Y.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Ensembl_get_assembly_info',
     'arguments': {'species': '...'}},
)
```

## `Ensembl_get_cross_references` {#Ensembl_get_cross_references}

Gets external cross-references (HGNC, Entrez, UniProt, OMIM) for an Ensembl stable ID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensembl_id`
- **Keywords:** xref, hgnc, entrez, omim
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensembl_id` | string | yes | Ensembl stable identifier. Can be gene (ENSG*), transcript (ENST*), or translation (ENSP*). Examples: 'ENSG00000141510' (TP53 gene), 'ENST00000269305' (TP53 canonical transcript). |
| `external_db` | string \| null | no | Optional: filter by external database name. Examples: 'HGNC', 'EntrezGene', 'Uniprot_gn', 'MIM_GENE', 'Reactome_gene', 'GeneCards'. If not specified, returns all cross-references. |

??? note "Upstream description"

    Get all external database cross-references for an Ensembl stable identifier. Returns linked records from HGNC, EntrezGene, UniProt, GeneCards, Reactome, MIM/OMIM, ArrayExpress, and other databases. Essential for mapping between Ensembl and other identifier systems, discovering related database entries, or building cross-database links for a gene. Can filter by specific external database. Example: 'ENSG00000141510' (TP53) returns 157 cross-references across 11 databases including HGNC:11998, EntrezGene:7157, GeneCards:TP53, MIM_GENE:191170, and Reactome pathways.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Ensembl_get_cross_references',
     'arguments': {'ensembl_id': '...'}},
)
```

## `Ensembl_get_gene_overlapping_features` {#Ensembl_get_gene_overlapping_features}

Gets genomic features overlapping an Ensembl gene.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** overlapping genes, antisense, regulatory elements
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | Ensembl gene ID. Examples: 'ENSG00000141510' (TP53), 'ENSG00000012048' (BRCA1), 'ENSG00000157764' (BRAF). |
| `feature_types` | string | no | Comma-separated feature types. Options: 'gene', 'transcript', 'regulatory', 'constrained', 'variation', 'repeat'. Default: 'gene'. |

??? note "Upstream description"

    Get features overlapping an Ensembl gene by gene ID. Returns all genomic features co-located with the gene, including neighboring genes, overlapping transcripts, regulatory elements, and antisense transcripts. Useful for finding overlapping genes, identifying all transcript isoforms in a locus, or discovering regulatory features at a gene. Example: 'ENSG00000141510' (TP53) returns 2 overlapping genes (TP53 and WRAP53), 42 transcripts across 4 biotypes (protein_coding, nonsense_mediated_decay, retained_intron, protein_coding_CDS_not_defined), and 9 regulatory features.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Ensembl_get_gene_overlapping_features',
     'arguments': {'gene_id': '...'}},
)
```

## `ensembl_get_genetree` {#ensembl_get_genetree}

Gets an Ensembl gene tree for a gene family.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `id`
- **Keywords:** gene tree, orthologues, paralogues, phylogeny
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `id` | string | yes | Gene tree ID (e.g., 'ENSGT00390000003602'). Gene tree IDs can be found in gene lookup responses with expand=1, or by querying Ensembl Compara database. |
| `prune_species` | string | no | Prune tree to specific species (optional, e.g., 'homo_sapiens,mus_musculus') |
| `prune_taxon` | string | no | Prune tree to specific taxon (optional) |
| `aligned` | boolean | no | Return aligned sequences Default `False` |
| `sequence` | string | no | Include sequences in tree (cdna, protein, or none) One of `cdna`, `protein`, `none` Default `none` |

??? note "Upstream description"

    Get gene tree (phylogenetic tree) for a gene. Returns evolutionary relationships showing orthologues and paralogues across species. Use this to study gene evolution and find related genes in other species. Gene tree IDs can be found in gene lookup responses or by querying the Compara database.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ensembl_get_genetree', 'arguments': {'id': '...'}},
)
```

## `ensembl_get_homology` {#ensembl_get_homology}

Gets Ensembl homologues of a gene symbol across species.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `species`, `symbol`
- **Keywords:** orthologues, paralogues, homology, comparative genomics
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `species` | string | yes | Source species (e.g., 'human', 'homo_sapiens') |
| `symbol` | string | yes | Gene symbol (e.g., 'BRCA2', 'TP53') |
| `target_species` | string | no | Target species to find homologues in (optional, e.g., 'mouse', 'mus_musculus'). Leave empty to search all species. |
| `target_taxon` | string | no | Target taxonomic group (optional, e.g., 'Mammalia', 'Primates') |
| `type` | string | no | Homology type: 'orthologues' (different species), 'paralogues' (same species), 'all' (both) One of `orthologues`, `paralogues`, `all` Default `all` |
| `sequence` | string | no | Include sequences in response One of `cdna`, `protein`, `none` Default `none` |
| `aligned` | boolean | no | Return aligned sequences Default `False` |

??? note "Upstream description"

    Get homology (orthologues and paralogues) for a gene by symbol across species. Returns evolutionary relationships with gene IDs, protein IDs, chromosomal locations, and homology types (1-to-1, 1-to-many, etc). Essential for comparative genomics studies.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ensembl_get_homology',
     'arguments': {'species': '...', 'symbol': '...'}},
)
```

## `ensembl_get_overlap_features` {#ensembl_get_overlap_features}

Gets genes, transcripts and other features overlapping a genomic region from Ensembl.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `species`, `region`
- **Keywords:** overlap, genomic region, genes, transcripts, biotype
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `species` | string | yes | Species (e.g., 'human', 'homo_sapiens') |
| `region` | string | yes | Genomic region (e.g., '7:140424943-140624564') |
| `feature` | string | no | Feature type to retrieve One of `gene`, `transcript`, `cds`, `exon`, `repeat`, `simple_feature`, `misc_feature`, `constrained_element` Default `gene` |
| `biotype` | string | no | Filter by biotype (optional, e.g., 'protein_coding', 'lincRNA', 'miRNA') |

??? note "Upstream description"

    Get genomic features (genes, transcripts, repeats) overlapping a region. Returns comprehensive annotations including coding genes, ncRNAs, pseudogenes, and structural features. Use this for genome annotation and feature discovery.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ensembl_get_overlap_features',
     'arguments': {'species': '...', 'region': '...'}},
)
```

## `Ensembl_get_region_features` {#Ensembl_get_region_features}

Gets genes, transcripts, regulatory and constrained elements in a region.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `region`
- **Keywords:** genomic region, features, constrained elements
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `species` | string | no | Species name. Default: 'human'. Examples: 'human', 'mouse', 'rat', 'zebrafish'. |
| `region` | string | yes | Genomic region in format 'chr:start-end'. Examples: '17:7661779-7687546' (TP53), '7:140424943-140524564' (BRAF), '13:32315086-32400266' (BRCA2). Max region size ~5Mb. |
| `feature_types` | string | no | Comma-separated feature types to retrieve. Options: 'gene', 'transcript', 'regulatory', 'constrained', 'variation', 'repeat', 'motif'. Default: 'gene'. Example: 'gene,regulatory,constrained'. |

??? note "Upstream description"

    Get all genomic features overlapping a specified chromosomal region. Returns genes, transcripts, regulatory elements, constrained elements, and other features within the coordinates. Essential for understanding the functional landscape of a genomic region, finding genes near a variant, or identifying regulatory elements in a locus. Supports multiple feature types in one query. Example: region '7:140424943-140524564' (BRAF locus) returns 5 genes, 20 regulatory elements, and 82 constrained elements; region '17:7661779-7687546' (TP53 locus) returns genes, enhancers, CTCF binding sites, and promoters.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Ensembl_get_region_features',
     'arguments': {'region': '...'}},
)
```

## `ensembl_get_species` {#ensembl_get_species}

Gets Ensembl species information with assembly and database version.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** ensembl species, assembly, genome build
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `division` | string | no | Filter by Ensembl division (optional) One of `EnsemblVertebrates`, `EnsemblPlants`, `EnsemblFungi`, `EnsemblBacteria`, `EnsemblProtists`, `EnsemblMetazoa` |

??? note "Upstream description"

    Get information about a species or list all available species in Ensembl. Returns scientific name, common name, assembly details, and database version. Use without ID to list all species (200+ available).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ensembl_get_species', 'arguments': {}},
)
```

## `Ensembl_get_species_info` {#Ensembl_get_species_info}

Searches Ensembl species and their current genome assemblies.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** ensembl species, assembly, taxid
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `search` | string | no | Search term to filter species. Can be a species name ('mus_musculus'), common name ('mouse'), display name ('Mouse'), or taxonomy ID ('10090'). Leave empty to list all 348+ species. |

??? note "Upstream description"

    Search and list species available in the Ensembl database with their genome assembly information. Returns species name, common name, taxonomy ID, current assembly version, and Ensembl division. Supports filtering by species name, common name, or taxonomy ID. Useful for discovering which genomes are available and their assembly versions before querying other Ensembl endpoints. Example: searching 'mouse' returns Mus musculus (taxon 10090, assembly GRCm39) plus related species.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Ensembl_get_species_info', 'arguments': {}},
)
```

## `ensembl_get_xrefs` {#ensembl_get_xrefs}

Gets external database IDs for an Ensembl object.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `id`
- **Keywords:** xref, id mapping, uniprot, refseq, hgnc
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `id` | string | yes | Ensembl ID (gene, transcript, or protein ID, e.g., 'ENSG00000139618', 'ENST00000380152', 'ENSP00000369497') |
| `external_db` | string | no | Filter by external database name (optional, e.g., 'UniProt', 'RefSeq', 'HGNC') |
| `object_type` | string | no | Object type filter (optional, e.g., 'gene', 'transcript', 'translation') |

??? note "Upstream description"

    Get cross-references (external database IDs) for an Ensembl gene, transcript, or protein. Returns mappings to databases like UniProt, RefSeq, HGNC, MGI, and many others. Use this to convert between different database identifiers.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ensembl_get_xrefs', 'arguments': {'id': '...'}},
)
```

## `ensembl_get_xrefs_by_name` {#ensembl_get_xrefs_by_name}

Gets external database cross-references for a gene name or symbol.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `name`
- **Keywords:** xref, gene symbol, id mapping
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `name` | string | yes | Gene symbol or name (e.g., 'BRCA1', 'TP53') |
| `species` | string | no | Species name (default 'homo_sapiens') Default `homo_sapiens` |
| `external_db` | string | no | Filter by external database name (optional) |

??? note "Upstream description"

    Get cross-references for an Ensembl object by its name/symbol. Returns mappings to external databases. Use this when you have a gene symbol or name instead of an Ensembl ID.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ensembl_get_xrefs_by_name', 'arguments': {'name': '...'}},
)
```

## `ensembl_lookup_gene` {#ensembl_lookup_gene}

Looks up an Ensembl gene by stable ID or symbol with location and description.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** gene lookup, ensembl id, location
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | Ensembl gene ID or symbol (e.g., 'ENSG00000139618' or 'BRCA1'). If using a stable ID, the tool will automatically route to /lookup/id endpoint. |
| `species` | string | no | Species name required for gene symbols (default 'homo_sapiens'). Examples: 'homo_sapiens', 'mus_musculus', 'rattus_norvegicus' Default `homo_sapiens` |

??? note "Upstream description"

    Lookup gene information by Ensembl gene ID or symbol. Returns gene details including location, description, and external references. Supports both stable IDs (e.g., ENSG00000139618) and gene symbols (e.g., BRCA1). For gene symbols, species parameter is required.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ensembl_lookup_gene', 'arguments': {'gene_id': '...'}},
)
```

## `Ensembl_lookup_gene_by_symbol` {#Ensembl_lookup_gene_by_symbol}

Looks up Ensembl gene and transcript IDs for a gene symbol.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `symbol`
- **Keywords:** gene symbol, ensembl id, transcript id
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `symbol` | string | yes | Gene symbol to look up. Examples: 'TP53', 'BRCA1', 'EGFR', 'BRAF', 'KRAS'. |
| `species` | string | no | Species name. Default: 'human'. Examples: 'human', 'mouse', 'rat', 'zebrafish'. |
| `external_db` | string \| null | no | Optional: filter by external database. Examples: 'HGNC', 'EntrezGene'. |

??? note "Upstream description"

    Look up Ensembl gene IDs for a gene symbol across external databases. Given a gene symbol (e.g., 'TP53'), returns all matching Ensembl gene and transcript IDs. Useful for converting gene symbols to Ensembl stable identifiers for use with other Ensembl APIs. Example: symbol 'TP53' in human returns ENSG00000141510 (gene) and related transcript/translation IDs; 'BRCA1' returns ENSG00000012048.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Ensembl_lookup_gene_by_symbol',
     'arguments': {'symbol': '...'}},
)
```

## `EnsemblArchive_batch_lookup` {#EnsemblArchive_batch_lookup}

Looks up version histories and current status for up to 50 Ensembl stable IDs.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensembl_ids`
- **Keywords:** ensembl archive, deprecated id, id history
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensembl_ids` | string | yes | Comma-separated Ensembl stable IDs (max 50). Example: 'ENSG00000141510,ENSG00000012048' for TP53 and BRCA1. |

??? note "Upstream description"

    Batch lookup version histories for multiple Ensembl stable IDs. Returns current status, release, and assembly for each ID. Supports up to 50 IDs per request. Useful for checking a list of gene IDs for currency, identifying deprecated IDs, and finding their replacements. Example: ENSG00000141510 (TP53) and ENSG00000012048 (BRCA1) both return as current at release 115 on GRCh38.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblArchive_batch_lookup',
     'arguments': {'ensembl_ids': '...'}},
)
```

## `EnsemblArchive_get_id_history` {#EnsemblArchive_get_id_history}

Gets the version history and current status of an Ensembl stable ID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensembl_id`
- **Keywords:** ensembl archive, stable id, release
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensembl_id` | string | yes | Ensembl stable ID to look up (gene: ENSG*, transcript: ENST*, protein: ENSP*, e.g., 'ENSG00000141510' for TP53). |

??? note "Upstream description"

    Get the version history and current status of an Ensembl stable ID. Returns the latest version string, current release number, assembly, and whether the ID is still active. Useful for checking if a gene/transcript/protein ID is current or deprecated, and what replaced it. Example: ENSG00000141510 (TP53) is current at release 115 on GRCh38 with latest version ENSG00000141510.20.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblArchive_get_id_history',
     'arguments': {'ensembl_id': '...'}},
)
```

## `EnsemblCompara_get_cafe_tree` {#EnsemblCompara_get_cafe_tree}

Gets the Ensembl Compara CAFE tree with ancestral gene copy numbers for a gene family.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** cafe, gene family, expansion, contraction
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_tree_id` | string \| null | no | Ensembl gene tree (family) stable ID, e.g. 'ENSGT00390000003602'. Use this OR 'gene'. |
| `gene` | string \| null | no | Gene symbol (e.g. 'BRCA2') or Ensembl gene ID (e.g. 'ENSG00000139618'). The CAFE tree for the family containing this gene is returned. Use this OR 'gene_tree_id'. |
| `species` | string \| null | no | Species for symbol lookup when 'gene' is a symbol (default: 'human'). Examples: 'human', 'mouse', 'homo_sapiens'. |

??? note "Upstream description"

    Retrieve the Ensembl Compara CAFE tree for a gene family: the gene gain/loss (expansion/contraction) analysis that infers ancestral gene-copy counts (n_members) at every node of the species tree, with the birth-death rate (lambda) and per-node p-values. Unlike EnsemblCompara_get_gene_tree (which returns the phylogenetic gene tree topology + members), this returns the CAFE family-size dynamics used to detect lineage-specific gene-family expansions/contractions. Provide either gene_tree_id (e.g. 'ENSGT00390000003602') or gene (symbol e.g. 'BRCA2', or Ensembl gene ID). Example: gene_tree_id='ENSGT00390000003602' returns type 'cafe tree', rooted, pvalue_avg, lambda ~5.08e-08, and per-node n_members / pvalue along the species tree.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblCompara_get_cafe_tree', 'arguments': {}},
)
```

## `EnsemblCompara_get_gene_tree` {#EnsemblCompara_get_gene_tree}

Gets the Ensembl Compara gene tree for a gene.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene`
- **Keywords:** gene tree, compara, duplication, speciation
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string | yes | Gene symbol or Ensembl gene ID. Examples: 'TP53', 'BRCA1', 'ENSG00000141510'. Ensembl IDs are more reliable than symbols for gene tree retrieval. |
| `species` | string | no | Species (default: 'human'). Used when providing a gene symbol. Examples: 'human', 'mouse'. |

??? note "Upstream description"

    Get the phylogenetic gene tree showing evolutionary relationships of a gene family across species from Ensembl Compara. Gene trees trace the history of speciation and duplication events for a gene family. Returns tree members (genes across species), tree structure, and optionally Newick format. Useful for evolutionary analysis and understanding gene family expansion. Example: TP53 gene tree shows conservation across vertebrates.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblCompara_get_gene_tree',
     'arguments': {'gene': '...'}},
)
```

## `EnsemblCompara_get_orthologues` {#EnsemblCompara_get_orthologues}

Finds orthologues of a gene using Ensembl Compara.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene`
- **Keywords:** orthologues, compara, mouse, zebrafish
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string | yes | Gene symbol or Ensembl gene ID. Examples: 'BRCA1', 'TP53', 'ENSG00000141510', 'ENSG00000012048'. |
| `species` | string | no | Source species (default: 'human'). Examples: 'human', 'mouse', 'zebrafish', 'drosophila_melanogaster'. |
| `target_species` | string | no | Target species to limit orthologue search. Examples: 'mouse', 'zebrafish', 'rat', 'drosophila_melanogaster'. Default: all species. |
| `target_taxon` | integer | no | NCBI taxon ID to limit orthologue search. Examples: 10090 (mouse), 7955 (zebrafish), 7227 (fruit fly). Alternative to target_species. |

??? note "Upstream description"

    Find orthologous genes (genes in different species that evolved from a common ancestor) using Ensembl Compara. Essential for comparative genomics - finding the mouse, zebrafish, or fly equivalent of a human gene. Returns orthologue gene IDs, species, homology type (one-to-one, one-to-many), and taxonomy level of divergence. Example: TP53 in human has a one-to-one orthologue (Trp53) in mouse.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblCompara_get_orthologues',
     'arguments': {'gene': '...'}},
)
```

## `EnsemblCompara_get_paralogues` {#EnsemblCompara_get_paralogues}

Finds paralogues of a gene using Ensembl Compara.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene`
- **Keywords:** paralogues, gene duplication, compara
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string | yes | Gene symbol or Ensembl gene ID. Examples: 'TP53', 'BRCA1', 'ENSG00000141510'. |
| `species` | string | no | Species (default: 'human'). Examples: 'human', 'mouse', 'zebrafish'. |

??? note "Upstream description"

    Find paralogous genes (genes within the same species that arose from gene duplication events) using Ensembl Compara. Paralogues often have related but diverged functions. Returns paralogue gene IDs, duplication type, and evolutionary context. Example: TP53 has paralogues TP63 (ENSG00000073282) and TP73 (ENSG00000078900) in human, from ancient vertebrate gene duplication.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblCompara_get_paralogues',
     'arguments': {'gene': '...'}},
)
```

## `EnsemblMap_convert_coordinates` {#EnsemblMap_convert_coordinates}

Converts genomic coordinates between assemblies such as GRCh37 and GRCh38.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `source_assembly`, `chromosome`, `start`, `end`, `target_assembly`
- **Keywords:** liftover, grch37, grch38, assembly mapping
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `species` | string | no | Species name. Default 'human'. Examples: 'human', 'mouse', 'homo_sapiens', 'mus_musculus'. |
| `source_assembly` | string | yes | Source genome assembly name. Examples: 'GRCh37', 'GRCh38', 'GRCm38', 'GRCm39'. |
| `chromosome` | string | yes | Chromosome name (e.g., '1', '17', 'X', 'MT'). |
| `start` | integer | yes | Start position on the chromosome (1-based, inclusive). |
| `end` | integer | yes | End position on the chromosome (1-based, inclusive). |
| `target_assembly` | string | yes | Target genome assembly name. Examples: 'GRCh38', 'GRCh37'. |

??? note "Upstream description"

    Convert genomic coordinates between human genome assemblies (e.g., GRCh37/hg19 to GRCh38/hg38 or vice versa). Essential for migrating variant positions, GWAS loci, or BED intervals between assembly versions. Uses the Ensembl REST API assembly mapping endpoint. Example: converting BRAF V600E position chr7:140453136 from GRCh37 maps to chr7:140753336 in GRCh38. Supports any species with multiple assemblies.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblMap_convert_coordinates',
     'arguments': {'source_assembly': '...',
                   'chromosome': '...',
                   'start': 0,
                   'end': 0,
                   'target_assembly': '...'}},
)
```

## `EnsemblMap_translate_coordinates` {#EnsemblMap_translate_coordinates}

Maps cDNA or protein positions of an Ensembl ID to chromosomal coordinates.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensembl_id`, `start`, `end`
- **Keywords:** coordinate mapping, cdna, protein position, genomic
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensembl_id` | string | yes | Ensembl transcript ID (ENST*) for cDNA mapping, or Ensembl protein ID (ENSP*) for protein mapping. Examples: 'ENST00000269305' (TP53 transcript), 'ENSP00000269305' (TP53 protein). |
| `start` | integer | yes | Start position in transcript/protein coordinates (1-based). For cDNA: nucleotide position. For protein: amino acid position. |
| `end` | integer | yes | End position in transcript/protein coordinates (1-based). |

??? note "Upstream description"

    Map protein or cDNA positions to genomic coordinates using Ensembl. Converts positions from transcript coordinate space (cDNA nucleotide positions) or protein coordinate space (amino acid positions) to chromosomal positions. Essential for mapping protein mutations (e.g., p53 amino acids 100-200) or transcript features to their genomic locations. Example: mapping protein positions 100-200 of ENSP00000269305 (p53) returns genomic coordinates chr17:7674931-7676071 across 3 exons.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblMap_translate_coordinates',
     'arguments': {'ensembl_id': '...', 'start': 0, 'end': 0}},
)
```

## `EnsemblReg_get_constrained_elements` {#EnsemblReg_get_constrained_elements}

Gets constrained elements in a genomic region from Ensembl Compara.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `region`
- **Keywords:** constrained elements, conservation, purifying selection
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `species` | string | no | Species name. Use 'homo_sapiens' for human. Default: 'homo_sapiens'. Default `homo_sapiens` |
| `region` | string | yes | Genomic region in format 'chromosome:start-end'. Example: '17:7661779-7687538' for TP53, '7:140424943-140624564' for BRAF. |

??? note "Upstream description"

    Get evolutionarily constrained elements in a genomic region from the Ensembl Compara database. Constrained elements are regions under purifying selection across species, indicating functional importance. Higher scores indicate stronger conservation. Use this to identify regulatory regions, coding exons, or other functional elements that are preserved across evolution. Essential for prioritizing non-coding variants and understanding regulatory conservation. Example: TP53 locus (chr17:7661779-7687538) returns 20 constrained elements with conservation scores up to 47.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblReg_get_constrained_elements',
     'arguments': {'region': '...'}},
)
```

## `EpiGraphDB_search_gene` {#EpiGraphDB_search_gene}

Searches EpiGraphDB genes by symbol or Ensembl ID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** gene search, druggability, chromosome, epigraphdb
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_name` | string \| null | no | Gene symbol to search for (e.g., 'BRCA1', 'TP53', 'EGFR'). Either gene_name or gene_id required. |
| `gene_id` | string \| null | no | Ensembl gene ID (e.g., 'ENSG00000012048'). Either gene_name or gene_id required. |
| `limit` | integer | no | Maximum results to return (default 10, max 50). Default `10` |

??? note "Upstream description"

    Search for gene information in EpiGraphDB by gene symbol or Ensembl ID. Returns gene metadata including druggability classification, chromosomal location, and gene type. Example: 'BRCA1' returns Tier 2 druggability, chr17 location.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EpiGraphDB_search_gene', 'arguments': {}},
)
```

## `FlyBase_get_gene` {#FlyBase_get_gene}

Gets Drosophila gene symbol, name, synopsis, genomic location, synonyms and cross-references from FlyBase via the Alliance.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** drosophila, fbgn, fruit fly, alliance genome, gene record
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | FlyBase gene ID with 'FB:' prefix. Examples: 'FB:FBgn0003996' (white/w), 'FB:FBgn0000490' (decapentaplegic/dpp), 'FB:FBgn0004647' (Notch/N), 'FB:FBgn0000546' (eyeless/ey), 'FB:FBgn0003731' (hedgehog/hh). |

??? note "Upstream description"

    Get detailed Drosophila melanogaster (fruit fly) gene information from FlyBase via the Alliance of Genome Resources. Returns gene symbol, name, synopsis, genomic location, synonyms, and cross-references. FlyBase is the primary database for Drosophila genetics. Example: 'FB:FBgn0003996' returns the white (w) gene with its function as an ABCG2-class transporter.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'FlyBase_get_gene', 'arguments': {'gene_id': '...'}},
)
```

## `FlyBase_get_gene_interactions` {#FlyBase_get_gene_interactions}

Returns molecular (protein-protein) or genetic (suppression/enhancement) interactions for a Drosophila gene.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** drosophila, genetic interaction, suppressor, enhancer, ppi
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | FlyBase gene ID with 'FB:' prefix. Examples: 'FB:FBgn0003996' (white/w), 'FB:FBgn0000490' (dpp), 'FB:FBgn0004647' (Notch/N). |
| `interaction_type` | string | no | Type of interaction: 'molecular' (physical PPI) or 'genetic' (phenotype-based). Default: 'molecular'. |
| `limit` | integer | no | Maximum number of interactions to return (1-100). Default: 20. |
| `page` | integer | no | Page number for pagination. Default: 1. |

??? note "Upstream description"

    Get molecular or genetic interactions for a Drosophila gene from the Alliance of Genome Resources. Molecular interactions are physical protein-protein interactions. Genetic interactions are phenotype-based (suppression, enhancement, etc.). Returns interacting gene pairs with species and interaction type. Example: 'FB:FBgn0000490' (dpp) has 169 genetic interactions including genetic suppression with anchor.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'FlyBase_get_gene_interactions',
     'arguments': {'gene_id': '...'}},
)
```

## `FlyBase_get_gene_orthologs` {#FlyBase_get_gene_orthologs}

Lists orthologs of a Drosophila gene in human, mouse, zebrafish, worm, yeast and others with prediction methods.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** drosophila, orthologs, homology, cross-species, model organism
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | FlyBase gene ID with 'FB:' prefix. Examples: 'FB:FBgn0004647' (Notch/N), 'FB:FBgn0003996' (white/w), 'FB:FBgn0000546' (eyeless/ey). |
| `limit` | integer | no | Maximum number of orthologs to return (1-100). Default: 20. |
| `page` | integer | no | Page number for pagination. Default: 1. |
| `stringency` | string | no | Stringency filter for orthology calls: 'stringent' (high confidence) or 'moderate'. Default: 'stringent'. |

??? note "Upstream description"

    Get orthologs of a Drosophila gene across other species from the Alliance of Genome Resources. Returns orthologous genes in human, mouse, zebrafish, worm, yeast, and other model organisms with prediction method details and best-score indicators. Essential for translational research from fly models to human disease genes. Example: 'FB:FBgn0004647' (Notch) returns orthologs including human NOTCH1-4.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'FlyBase_get_gene_orthologs',
     'arguments': {'gene_id': '...'}},
)
```

## `FlyMine_search` {#FlyMine_search}

Searches FlyMine across Drosophila genes, proteins, publications, protein domains and GO terms, with category facets.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `q`
- **Keywords:** flymine, intermine, drosophila, genomics warehouse, keyword search
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `q` | string | yes | Search query: gene symbol (e.g., 'dpp', 'Notch', 'hedgehog'), protein name, domain name, or free text. |
| `size` | integer | no | Number of results to return (default: 10, max: 100). Default `10` |
| `format` | string | no | Response format. Must be 'json'. Default `json` |
| `facet_Category` | string | no | Filter by entity category. Options: 'Gene', 'Protein', 'Publication', 'ProteinDomain', 'CDNAClone', 'GOTerm', 'InteractionExperiment', 'RNAiScreen', 'UniProtFeature'. |

??? note "Upstream description"

    Search FlyMine, the InterMine data warehouse for Drosophila melanogaster (fruit fly) genomics data from FlyBase. Searches across genes, proteins, publications, protein domains, gene ontology terms, and more. FlyMine integrates data from FlyBase, UniProt, Reactome, KEGG, TreeFam, and other model organism databases. Returns matching records with category facets and relevance scores. Use facet_Category parameter to filter by entity type (e.g., 'Gene', 'Protein', 'Publication', 'ProteinDomain', 'CDNAClone', 'GOTerm', 'InteractionExperiment').

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'FlyMine_search', 'arguments': {'q': '...'}},
)
```

## `FlyMine_search_genes` {#FlyMine_search_genes}

Searches FlyMine for Drosophila genes and returns symbols, FBgn IDs, CG numbers and organism.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `q`
- **Keywords:** flymine, drosophila, fbgn, cg number, gene search
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `q` | string | yes | Gene search query: gene symbol (e.g., 'dpp', 'Notch', 'wg', 'hh'), FlyBase ID (e.g., 'FBgn0000490'), or keyword. |
| `size` | integer | no | Number of results to return (default: 10). Default `10` |

??? note "Upstream description"

    Search FlyMine specifically for Drosophila genes. This is a convenience wrapper that filters results to the Gene category only. Returns gene symbols, FlyBase gene IDs (FBgn numbers), secondary identifiers (CG numbers), and organism information. FlyMine covers 12+ Drosophila species plus other model organisms. For broader searches across all entity types, use FlyMine_search instead.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'FlyMine_search_genes', 'arguments': {'q': '...'}},
)
```

## `G2P_get_gene` {#G2P_get_gene}

Gets Gene2Phenotype gene info: location, HGNC/Ensembl/OMIM cross-references, synonyms and coordinates.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_symbol`
- **Keywords:** gene2phenotype, hgnc, omim, gene location
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_symbol` | string | yes | Gene symbol (e.g., 'BRCA1', 'TP53', 'EGFR', 'SCN1A') |

??? note "Upstream description"

    Get Gene2Phenotype information for a specific gene by symbol. Returns gene location, cross-references to HGNC/Ensembl/OMIM, gene synonyms, and genomic coordinates. This provides the gene-level data without disease associations. For disease associations, use G2P_search with the gene name instead.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'G2P_get_gene', 'arguments': {'gene_symbol': '...'}},
)
```

## `GenomeNexus_get_canonical_transcript` {#GenomeNexus_get_canonical_transcript}

Gets a gene's canonical Ensembl transcript, protein ID, RefSeq, CCDS and Pfam domains from Genome Nexus.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_symbol`
- **Keywords:** canonical transcript, ensembl, refseq, ccds, pfam
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_symbol` | string | yes | HUGO gene symbol. Examples: 'TP53', 'BRAF', 'BRCA1', 'EGFR', 'KRAS'. |

??? note "Upstream description"

    Get the canonical transcript for a gene from Genome Nexus (MSK). Returns the Ensembl transcript ID, protein ID, protein length, RefSeq mRNA ID, CCDS ID, HUGO symbol, and Pfam protein domain annotations. Useful for determining the standard transcript to use for variant annotation, and for getting a quick overview of protein domain architecture. Input is a HUGO gene symbol (e.g., 'TP53', 'BRAF', 'BRCA1').

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GenomeNexus_get_canonical_transcript',
     'arguments': {'gene_symbol': '...'}},
)
```

## `gnomad_get_gene` {#gnomad_get_gene}

Gets gnomAD gene metadata (Ensembl ID, coordinates, canonical transcript) by symbol or Ensembl ID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** gnomad, ensembl gene, coordinates, canonical transcript
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_symbol` | string | no | Gene symbol (e.g., 'BRCA1'). |
| `gene_id` | string | no | Ensembl gene ID (e.g., 'ENSG00000012048'). |
| `reference_genome` | string | no | Reference genome. One of `GRCh37`, `GRCh38` Default `GRCh38` |

??? note "Upstream description"

    Get basic gene metadata from gnomAD by `gene_symbol` or `gene_id` (Ensembl gene ID). Use `gnomad_search_genes` to discover `gene_id`. Use the returned `canonical_transcript_id` with `gnomad_get_transcript`. The returned `chrom`/`start`/`stop` are in the assembly named by `reference_genome`, which gnomAD reports itself on both the `gene` object and the top level of the response (GRCh38 unless the `reference_genome` argument says otherwise). Coordinates are not comparable across assemblies: BRCA2 is 13:32,315,086-32,400,268 on GRCh38 and 13:32,889,611-32,973,805 on GRCh37.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gnomad_get_gene', 'arguments': {}},
)
```

## `gnomad_get_transcript` {#gnomad_get_transcript}

Gets gnomAD transcript metadata and coordinates with linked gene by Ensembl transcript ID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `transcript_id`
- **Keywords:** ensembl transcript, enst, coordinates, gnomad
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `transcript_id` | string | yes | Ensembl transcript ID (e.g., 'ENST00000357654'). |
| `reference_genome` | string | no | Reference genome. One of `GRCh37`, `GRCh38` Default `GRCh38` |

??? note "Upstream description"

    Get basic transcript metadata from gnomAD by Ensembl transcript ID (e.g., ENST...). The response includes the linked `gene` object (gene ID + symbol). Use `gnomad_get_gene` to find a gene's `canonical_transcript_id`. The returned `chrom`/`start`/`stop` are in the assembly named by `reference_genome`, which gnomAD reports itself on both the `transcript` object and the top level of the response (GRCh38 unless the `reference_genome` argument says otherwise). Coordinates are not comparable across assemblies.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gnomad_get_transcript',
     'arguments': {'transcript_id': '...'}},
)
```

## `gnomad_search_genes` {#gnomad_search_genes}

Searches gnomAD genes by free text (usually symbol) and returns Ensembl IDs and symbols.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** gene search, ensembl id, gene symbol, gnomad
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Gene search string (e.g., 'BRCA1'). |
| `reference_genome` | string | no | Reference genome for the search. One of `GRCh37`, `GRCh38` Default `GRCh38` |

??? note "Upstream description"

    Search for genes in gnomAD by free-text query (typically a gene symbol). Returns matching `ensembl_id` (Ensembl gene ID) and `symbol`. Use `ensembl_id` as `gene_id` for `gnomad_get_gene`. The response reports `reference_genome`, the assembly the search ran against (GRCh38 unless `reference_genome` says otherwise); results are assembly-specific, e.g. BRCA2 resolves to `ensembl_version` 17 on GRCh38 and 10 on GRCh37.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gnomad_search_genes', 'arguments': {'query': '...'}},
)
```

## `GO_get_annotations_for_gene` {#GO_get_annotations_for_gene}

Returns a page of GO annotations for a gene/protein with evidence codes via GOlr.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** gene ontology, go annotation, evidence code, golr
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | A gene identifier such as gene symbol (e.g., 'TP53') or database ID. |
| `rows` | integer | no | Maximum number of annotation rows to return. Default: 100. Raise it (e.g. 2000) for a heavily annotated gene such as TP53, which has 1812 annotation rows; lower it to keep responses small. Default `100` |

??? note "Upstream description"

    Finds GO annotations for a specific gene/protein using GOlr search. Returns one row per annotation, so a gene with several evidence codes for the same GO term appears several times. This returns a PAGE, not the complete set: heavily annotated genes have far more annotations than the default page holds (measured: 'TP53' matches 1812 annotation rows), and the response carries no total, so a short list is NOT evidence the gene has no further annotations -- raise 'rows' if completeness matters.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GO_get_annotations_for_gene',
     'arguments': {'gene_id': '...'}},
)
```

## `GO_get_genes_for_term` {#GO_get_genes_for_term}

Finds distinct genes annotated to a GO term and its descendants, optionally restricted to a taxon.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `id`
- **Keywords:** gene ontology, gene set, go term, taxon
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `id` | string | yes | The standard GO term ID, e.g., 'GO:0006915'. |
| `taxon` | string | no | Optional species filter using an NCBI taxon ID, e.g. 'NCBITaxon:9606' (human) or 'NCBITaxon:10090' (mouse). Omit to return genes from every species. |
| `rows` | integer | no | Maximum number of annotation rows to scan. Distinct genes returned may be fewer, since one gene can carry several annotations. Default is 100. Default `100` |

??? note "Upstream description"

    Finds genes/proteins annotated to a Gene Ontology term, including annotations to its descendant terms (is_a/part_of closure), via the GO Solr index. Rows are collapsed to distinct genes, each carrying its evidence codes and annotation count. Supply `taxon` to restrict to one species -- omit it and genes from all species are returned together. Example: id='GO:0006915' with taxon='NCBITaxon:9606' returns human apoptosis genes.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GO_get_genes_for_term', 'arguments': {'id': '...'}},
)
```

## `GOAPI_get_gene_functions` {#GOAPI_get_gene_functions}

Gets GO biological process, molecular function and cellular component annotations for a gene with evidence types.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** gene ontology, gene function, evidence code, biological process
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | Gene identifier as a CURIE. Examples: 'HGNC:11998' (TP53), 'UniProtKB:P04637' (p53 protein), 'HGNC:1100' (BRCA1), 'MGI:MGI:98834' (Trp53 mouse). For human genes, use HGNC prefix. |
| `rows` | integer | no | Maximum number of annotations to return (default: 20, max: 100). |
| `aspect` | string | no | Filter by GO aspect: 'P' (Biological Process), 'F' (Molecular Function), 'C' (Cellular Component). Default: all aspects. |

??? note "Upstream description"

    Get Gene Ontology (GO) annotations for a specific gene, showing its known biological processes, molecular functions, and cellular components. Returns GO terms with evidence types (IDA=direct assay, IMP=mutant phenotype, ISS=sequence similarity, etc.) and literature references. Input is a gene CURIE identifier. Example: HGNC:11998 (TP53) returns annotations like 'apoptotic process', 'DNA binding', 'nucleus'.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GOAPI_get_gene_functions',
     'arguments': {'gene_id': '...'}},
)
```

## `GOAPI_get_genes_by_function` {#GOAPI_get_genes_by_function}

Gets genes/proteins annotated with a GO term, optionally filtered by taxon.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `go_id`
- **Keywords:** gene ontology, gene set, go term, taxon
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `go_id` | string | yes | Gene Ontology term ID. Examples: 'GO:0006915' (apoptotic process), 'GO:0003677' (DNA binding), 'GO:0005634' (nucleus). |
| `rows` | integer | no | Maximum number of gene annotations to return (default: 20, max: 100). |
| `taxon` | string | no | Filter by taxon. Format: NCBITaxon:ID. Examples: 'NCBITaxon:9606' (human), 'NCBITaxon:10090' (mouse), 'NCBITaxon:7227' (Drosophila). |

??? note "Upstream description"

    Get genes annotated with a specific Gene Ontology (GO) term. Find all genes/proteins that have been experimentally or computationally associated with a biological process, molecular function, or cellular component. Optionally filter by taxon. Example: GO:0006915 (apoptotic process) returns genes like TP53, BAX, BCL2 that function in apoptosis.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GOAPI_get_genes_by_function',
     'arguments': {'go_id': '...'}},
)
```

## `GoaT_get_species` {#GoaT_get_species}

Looks up genome sequencing status, assembly level, genome size, chromosome number and ploidy via GoaT.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `taxon`
- **Keywords:** goat, genome size, assembly level, earth biogenome, ploidy, chromosome number
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `taxon` | string | yes | Scientific name (e.g. 'Panthera leo') or NCBI taxon id (e.g. '9689'). |
| `include_descendants` | boolean \| null | no | If true and taxon is a numeric id, include all descendant taxa (e.g. every species in a genus). Default false. |
| `limit` | integer \| null | no | Max taxa to return, 1-50. Default 10. |

??? note "Upstream description"

    Look up genome sequencing status via GoaT (Genomes on a Tree), which indexes the Earth BioGenome Project, Darwin Tree of Life, Vertebrate Genomes Project, and ERGA against the full NCBI taxonomy. Returns assembly level, genome size, chromosome number, and ploidy. Example: taxon='Panthera leo' (African lion) returns a chromosome-level assembly with 38 chromosomes. Accepts a scientific name or NCBI taxon id; set include_descendants for genus- or family-level coverage.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GoaT_get_species', 'arguments': {'taxon': '...'}},
)
```

## `gProfiler_convert_ids` {#gProfiler_convert_ids}

Converts gene identifiers (symbols, Ensembl, Entrez, UniProt, RefSeq) between namespaces with g:Convert.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_list`
- **Keywords:** id conversion, ensembl, entrez, refseq, g:convert
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_list` | string | yes | Comma-separated list of gene identifiers to convert. Examples: 'TP53,BRCA1,EGFR' (symbols), 'ENSG00000141510,ENSG00000012048' (Ensembl IDs). |
| `target_namespace` | string | no | Target namespace for conversion. Default: 'ENSG'. Options: ENSG (Ensembl Gene), ENST (Ensembl Transcript), ENSP (Ensembl Protein), UNIPROTSWISSPROT (UniProt), ENTREZGENE_ACC (Entrez), WIKIGENE_ACC (WikiGene), REFSEQ_MRNA_ACC, REFSEQ_PEPTIDE_ACC. |
| `organism` | string | no | Organism identifier. Default: 'hsapiens'. Examples: 'hsapiens', 'mmusculus'. |

??? note "Upstream description"

    Convert gene identifiers between different namespaces using g:Profiler (g:Convert) from the University of Tartu. Supports conversion between gene symbols, Ensembl IDs, Entrez IDs, UniProt accessions, RefSeq IDs, and many other identifiers. Input is a list of gene identifiers and a target namespace. Supports 500+ organisms. Example: convert 'TP53,BRCA1,EGFR' to Ensembl gene IDs (ENSG). Available target namespaces: ENSG (Ensembl Gene), ENST (Ensembl Transcript), ENSP (Ensembl Protein), UNIPROTSWISSPROT (UniProt reviewed), ENTREZGENE_ACC (Entrez Gene), WIKIGENE_ACC (WikiGene).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gProfiler_convert_ids',
     'arguments': {'gene_list': '...'}},
)
```

## `gProfiler_enrichment` {#gProfiler_enrichment}

Runs g:GOSt functional enrichment of a gene list against GO, KEGG, Reactome and other sources.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_list`
- **Keywords:** enrichment, g:gost, gene list, go, pathways
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_list` | string | yes | Comma-separated list of gene symbols. Examples: 'TP53,BRCA1,EGFR,KRAS,PTEN' or 'Trp53,Brca1,Egfr' (mouse). |
| `organism` | string | no | Organism identifier. Default: 'hsapiens'. Examples: 'hsapiens' (human), 'mmusculus' (mouse), 'dmelanogaster' (fly), 'scerevisiae' (yeast), 'drerio' (zebrafish). |
| `sources` | string | no | Comma-separated data sources to query. Default: 'GO:BP,GO:MF,GO:CC,KEGG,REAC,WP,HP'. Options: GO:BP, GO:MF, GO:CC, KEGG, REAC, WP, MIRNA, HPA, CORUM, HP, TF. |

??? note "Upstream description"

    Perform functional enrichment analysis on a gene list using g:Profiler (g:GOSt) from the University of Tartu. Identifies statistically overrepresented biological functions, pathways, and phenotypes in your gene list. Supports multiple data sources: Gene Ontology (GO:BP, GO:MF, GO:CC), KEGG pathways, Reactome, WikiPathways, Human Phenotype Ontology (HP), miRNA targets (MIRNA), CORUM protein complexes, and Transcription Factor targets (TF). Returns enriched terms with p-values (g:SCS multiple testing correction), term sizes, intersection genes, precision, and recall. Input is a comma-separated string of gene symbols (e.g., 'TP53,BRCA1,EGFR,KRAS,PTEN') -- a JSON array is not accepted. Supports 500+ organisms (use NCBI taxonomy prefix, e.g., 'hsapiens', 'mmusculus', 'dmelanogaster').

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gProfiler_enrichment', 'arguments': {'gene_list': '...'}},
)
```

## `gProfiler_find_orthologs` {#gProfiler_find_orthologs}

Maps genes from a source organism to orthologs in a target organism with g:Orth (Ensembl orthology).

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_list`
- **Keywords:** orthologs, g:orth, cross-species, gene list
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_list` | string | yes | Comma-separated list of gene identifiers from the source organism. Examples: 'TP53,BRCA1,EGFR' (human genes to find mouse orthologs). |
| `source_organism` | string | no | Source organism. Default: 'hsapiens'. Examples: 'hsapiens' (human), 'mmusculus' (mouse). |
| `target_organism` | string | no | Target organism to find orthologs in. Default: 'mmusculus'. Examples: 'mmusculus' (mouse), 'rnorvegicus' (rat), 'drerio' (zebrafish), 'dmelanogaster' (fly), 'celegans' (worm). |

??? note "Upstream description"

    Find orthologous genes across species using g:Profiler (g:Orth) from the University of Tartu. Maps genes from a source organism to their orthologs in a target organism using Ensembl orthology data. Useful for translating gene lists between model organisms. Example: find mouse orthologs of human TP53 and BRCA1. Supports 500+ organisms in the Ensembl database. Input is a list of gene identifiers from the source organism, and the target organism to map to.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'gProfiler_find_orthologs',
     'arguments': {'gene_list': '...'}},
)
```

## `GSEA_prerank` {#GSEA_prerank}

Runs pre-ranked GSEA on a ranked gene list against gene sets returning ES, NES, p-values and leading-edge genes.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** gsea, gene set enrichment, ranked list, nes, leading edge
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ranked_genes` | object \| null | no | \{gene: score\} mapping (score = DE statistic / fold change). Genes are ranked by score descending. Alternative to genes+scores. |
| `genes` | array \| null | no | Gene identifiers (use with `scores`, same order). |
| `scores` | array \| null | no | Per-gene ranking metric (same order as `genes`). |
| `gene_sets` | object \| null | no | \{set_name: \[gene, ...\]\} collection to test. Alternative to gene_set. |
| `gene_set` | array \| null | no | A single gene set (list of gene symbols). |
| `weight` | number \| null | no | Score weighting exponent for the running sum (default 1.0; classic GSEA). |
| `n_permutations` | integer \| null | no | Gene-label permutations for the p-value (default 1000). |

??? note "Upstream description"

    Pre-ranked Gene Set Enrichment Analysis (Subramanian 2005) on a ranked gene list (e.g. genes ranked by a DESeq2/edgeR statistic) against one or more gene sets. Returns per-set enrichment score (ES), normalized ES (NES), permutation p-value, and leading-edge genes. Pure-compute (no R, no bundled MSigDB - supply gene sets, e.g. from Enrichr). Complements over-representation tools by using the full ranking.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GSEA_prerank', 'arguments': {}},
)
```

## `GTDB_get_genome` {#GTDB_get_genome}

Gets GTDB genome metadata by accession: size, GC%, N50, CheckM completeness, gene counts and taxonomy.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `accession`
- **Keywords:** gtdb, genome assembly, checkm, n50, gc content, prokaryote
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type (fixed: get_genome) One of `get_genome` Default `get_genome` |
| `accession` | string | yes | NCBI genome accession (e.g., 'GCF_000005845.2' for E. coli K-12, 'GCF_000009045.1' for B. subtilis) |

??? note "Upstream description"

    Get comprehensive metadata for a specific genome by NCBI accession from GTDB. Returns nucleotide statistics (genome size, GC%, contig count, N50), gene-level metrics (completeness via CheckM/CheckM2, contamination, SSU/LSU counts, protein count, coding density), NCBI metadata (assembly level, BioProject, BioSample, isolation source, geographic coordinates), and taxonomic history across GTDB releases. Useful for evaluating genome quality, finding genome properties, and tracking taxonomic changes across GTDB versions.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GTDB_get_genome', 'arguments': {'accession': '...'}},
)
```

## `GTDB_get_species` {#GTDB_get_species}

Gets a GTDB species cluster with its assigned genome accessions, NCBI names and representative status.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `species`
- **Keywords:** gtdb, species cluster, ani, genome accessions, prokaryote
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type (fixed: get_species) One of `get_species` Default `get_species` |
| `species` | string | yes | Species name WITHOUT the s__ prefix (e.g., 'Escherichia coli', 'Hydrogenothermus marinus', 'Bacillus subtilis') |
| `max_genomes` | integer | no | Maximum number of genomes to return (default: 20). Species clusters can have thousands of genomes. |

??? note "Upstream description"

    Get detailed information about a GTDB species cluster, including all genomes assigned to that species. A GTDB species cluster groups genomes by average nucleotide identity (ANI). Returns the full list of genomes with their accessions, NCBI organism names, NCBI taxonomy, and whether they are the GTDB species representative or NCBI type material. Useful for exploring species diversity, finding representative genomes, and comparing GTDB vs NCBI classifications.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GTDB_get_species', 'arguments': {'species': '...'}},
)
```

## `GTDB_search_genomes` {#GTDB_search_genomes}

Searches GTDB prokaryotic genomes by organism name and returns genome accessions and taxonomy.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** gtdb, genome search, prokaryote, accession
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type (fixed: search_genomes) One of `search_genomes` Default `search_genomes` |
| `query` | string | yes | Search query for organism name (e.g., 'Escherichia', 'Staphylococcus aureus', 'Lactobacillus') |
| `page` | integer | no | Page number for pagination (default: 1) |
| `items_per_page` | integer | no | Number of results per page (default: 10, max: 50) |

??? note "Upstream description"

    Search for prokaryotic genomes in GTDB by organism name. A multi-word query (e.g. a binomial like 'Mycobacterium tuberculosis') is first matched precisely, requiring EVERY term of the query to appear in the genome's name or taxonomy; 'match_type' is then 'exact'. If nothing matches the whole query, the tool falls back to GTDB's loose search, which matches ANY term of the query and therefore returns a broader taxon (often just the genus) or a different species; that case is reported as 'match_type': 'broad' with a top-level 'note', and 'total_results' then counts the broader match rather than the queried organism -- always read 'match_type' and 'total_results_scope' before quoting a count or trusting a result. Results also match against the whole taxonomic clade rather than organism-name text alone, so exact/prefix name matches are ranked first within a page. Returns genome accessions with both NCBI and GTDB taxonomy classifications, allowing comparison between the two taxonomic systems. Each result shows the genome accession (GCA/GCF), NCBI organism name, NCBI taxonomy, GTDB taxonomy, and whether the genome is a GTDB species representative. Supports pagination. Useful for finding genome accessions, exploring taxonomic reclassifications, and linking genomes to the GTDB taxonomy.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GTDB_search_genomes', 'arguments': {'query': '...'}},
)
```

## `Harmonizome_get_dataset` {#Harmonizome_get_dataset}

Gets a Harmonizome dataset's description, association type, measurement and catalog of gene sets.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `dataset_name`
- **Keywords:** harmonizome, gene sets, dataset, ma'ayan lab
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `dataset_name` | string | yes | Exact dataset name from Harmonizome. Examples: 'CTD Gene-Disease Associations', 'KEGG Pathways', 'BioGRID Protein-Protein Interactions'. Use Harmonizome_list_datasets to find names. |
| `gene_set_limit` | integer | no | Maximum number of gene set names to return (default: 50). |

??? note "Upstream description"

    Get detailed information about a specific dataset integrated in Harmonizome (Ma'ayan Lab). Returns the dataset description, association type, measurement method, attribute type, PubMed references, and a catalog of gene sets within the dataset. Use Harmonizome_list_datasets first to find available dataset names. Example: 'CTD Gene-Disease Associations' returns 5218 gene sets with curated disease-gene associations.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Harmonizome_get_dataset',
     'arguments': {'dataset_name': '...'}},
)
```

## `Harmonizome_get_gene` {#Harmonizome_get_gene}

Gets Harmonizome gene info (symbol, name, Entrez ID, description) integrated from 100+ datasets.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_symbol`
- **Keywords:** harmonizome, gene summary, entrez, integrated data
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_symbol` | string | yes | Official gene symbol (case-insensitive). Examples: 'TP53', 'BRCA1', 'EGFR', 'KRAS', 'CDK2'. |

??? note "Upstream description"

    Get comprehensive gene information from Harmonizome (Ma'ayan Lab), which integrates data from 100+ genomics datasets covering gene expression, protein interactions, pathways, diseases, drug targets, and more. Returns gene symbol, name, NCBI Entrez ID, description, synonyms, and linked protein identifiers. Useful as a quick gene summary with curated descriptions from RefSeq. Example: TP53 returns description about tumor suppressor protein, NCBI ID 7157, synonyms \[TRP53, LFS1, P53, BMFS5, BCC7\], and linked UniProt protein P53_HUMAN.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Harmonizome_get_gene',
     'arguments': {'gene_symbol': '...'}},
)
```

## `Harmonizome_get_gene_set_members` {#Harmonizome_get_gene_set_members}

Returns member genes of a Harmonizome attribute gene set with values, or a gene's associated attributes.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** harmonizome, gene set members, attributes, association values
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `mode` | string | no | 'gene_set' for member genes of an attribute set (requires attribute + dataset); 'gene' for a gene's full association table (requires gene_symbol). One of `gene_set`, `gene` Default `gene_set` |
| `attribute` | string | no | Attribute/gene-set name for mode='gene_set' (e.g. 'heart', 'liver', a cell type or disease). Spaces are auto-encoded. |
| `dataset` | string | no | Dataset name for mode='gene_set' (e.g. 'GTEx Tissue Gene Expression Profiles'). Use Harmonizome_list_datasets to find names. |
| `gene_symbol` | string | no | Gene symbol for mode='gene' (e.g. 'DTD2', 'TP53'). |
| `limit` | integer | no | Maximum number of association rows to return (default 100). Sets can have thousands of members. Default `100` |

??? note "Upstream description"

    Retrieve the actual Harmonizome association payload (which Harmonizome_get_dataset and Harmonizome_get_gene deliberately omit). mode='gene_set' (default) returns the member genes of a curated attribute set with thresholdValue + standardizedValue - requires `attribute` (e.g. 'heart') and `dataset` (e.g. 'GTEx Tissue Gene Expression Profiles'). mode='gene' returns the full per-gene association table across all datasets - requires `gene_symbol` (e.g. 'DTD2'). Example: gene_set heart / GTEx Tissue Gene Expression Profiles -&gt; 6473 members incl. \{gene: DTD2, threshold_value: -1.0, standardized_value: -1.05248\}. gene DTD2 -&gt; 3787 associations across datasets.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Harmonizome_get_gene_set_members', 'arguments': {}},
)
```

## `Harmonizome_list_datasets` {#Harmonizome_list_datasets}

Lists 100+ integrated genomics datasets in Harmonizome (expression, pathways, regulation, disease).

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** harmonizome, datasets, gene sets, catalog
- **Read-only:** yes (declared by the server)

??? note "Upstream description"

    List all 100+ integrated genomics datasets available in Harmonizome (Ma'ayan Lab). Returns dataset names that can be used to explore gene-gene set associations. Covers gene expression (Allen Brain Atlas, GTEx, BioGPS, CCLE, TCGA), pathways (BioCarta, KEGG, Reactome, WikiPathways), protein interactions (BioGRID, PPI Hub), gene regulation (ChEA, ENCODE, JASPAR), drug targets (DrugBank, CTD), disease associations (OMIM, ClinVar, GAD, GWAS), phenotypes (MGI, HPO), and more.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Harmonizome_list_datasets', 'arguments': {}},
)
```

## `Harmonizome_search` {#Harmonizome_search}

Searches Harmonizome genes, datasets or attributes by keyword.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** harmonizome, gene search, dataset search, attributes
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search query. Examples: 'TP53', 'BRCA1', 'kinase', 'expression', 'pathway'. |
| `entity_type` | string | no | What to search for: 'gene' (default), 'dataset', or 'attribute'. |
| `limit` | integer | no | Maximum results to return (default: 20). |

??? note "Upstream description"

    Search Harmonizome across genes, datasets, or attributes by keyword. For gene search, returns matching gene symbols with names and NCBI IDs via exact lookup and suggest API. For dataset search, filters the 100+ dataset catalog by keyword. For attribute search, returns auto-complete suggestions. Use entity_type to control what is searched.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Harmonizome_search', 'arguments': {'query': '...'}},
)
```

## `Harmonizome_search_genes` {#Harmonizome_search_genes}

Searches Harmonizome genes by keyword returning matching gene symbols.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** harmonizome, gene search, gene symbol
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search keyword. Examples: 'kinase', 'tumor suppressor', 'receptor'. |
| `limit` | integer | no | Maximum results to return (default: 20). |

??? note "Upstream description"

    Search for genes by keyword in Harmonizome (Ma'ayan Lab). Returns matching gene symbols that can be used with Harmonizome_get_gene for detailed information. Note: Harmonizome API does not support server-side text filtering; this endpoint returns a paginated list from the gene catalog. Use Harmonizome_search with entity_type='gene' for keyword-based gene search.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Harmonizome_search_genes', 'arguments': {'query': '...'}},
)
```

## `HGNC_fetch_gene_by_id` {#HGNC_fetch_gene_by_id}

Fetches the HGNC gene record (symbol, name, location, cross-references) by HGNC ID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `hgnc_id`
- **Keywords:** hgnc, gene nomenclature, official symbol
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `hgnc_id` | string | yes | HGNC identifier, with or without 'HGNC:' prefix. Examples: 'HGNC:11998' (TP53), 'HGNC:1100' (BRCA1), '1101' (BRCA2). |

??? note "Upstream description"

    Fetch gene information from HGNC by HGNC ID. Useful when you have a specific HGNC identifier and want the full gene record. Returns official symbol, name, location, and cross-references to other databases. Example: 'HGNC:11998' returns TP53.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'HGNC_fetch_gene_by_id', 'arguments': {'hgnc_id': '...'}},
)
```

## `HGNC_fetch_gene_by_symbol` {#HGNC_fetch_gene_by_symbol}

Fetches HGNC gene data (ID, name, location, Ensembl/Entrez/UniProt/OMIM/RefSeq links) by symbol.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `symbol`
- **Keywords:** hgnc, gene symbol, nomenclature, cross-references
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `symbol` | string | yes | Official HGNC gene symbol. Examples: 'TP53', 'BRCA1', 'EGFR', 'INS', 'ACE2'. |

??? note "Upstream description"

    Fetch detailed gene information from HGNC (HUGO Gene Nomenclature Committee) by official gene symbol. Returns comprehensive gene data including HGNC ID, full name, chromosomal location, cross-references to Ensembl, NCBI Entrez, UniProt, OMIM, and RefSeq. HGNC is the authoritative source for human gene nomenclature. Example: 'TP53' returns tumor protein p53 with all its identifiers. Retired symbols are resolved automatically: if the symbol is not a current approved symbol, it is looked up as a previous symbol and then as an alias symbol, and the current gene's full record is returned with the substitution disclosed in metadata ('resolved_from', 'resolved_symbol', 'resolution_relation', 'symbol_resolution_note'). Example: 'CYP2E' returns CYP2E1 (renamed in 2002), 'DIA4' returns NQO1. If the retired symbol is ambiguous (shared by several genes) an error naming every candidate is returned instead of an arbitrary pick; if the symbol is unknown to HGNC entirely, metadata carries 'no_results_note' pointing at HGNC_search_genes.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'HGNC_fetch_gene_by_symbol',
     'arguments': {'symbol': '...'}},
)
```

## `HGNC_fetch_gene_family_members` {#HGNC_fetch_gene_family_members}

Lists all member genes of an HGNC gene group by gene_group_id with full records.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_group_id`
- **Keywords:** hgnc, gene family, gene group, members
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_group_id` | string | yes | Numeric HGNC gene group/family ID. Examples: '366' (Ubiquitin specific peptidases), '588' (BCL2 family), '475' (Solute carrier family). Accepts a string or integer. |

??? note "Upstream description"

    Enumerate ALL member genes of an HGNC gene family/group by its numeric gene_group_id, returning the full HGNC record for every member gene. Unlike HGNC_fetch_gene_by_symbol/by_id (single gene) this lists the entire family. Each member record includes hgnc_id, approved symbol, full name, chromosomal location, locus type, and cross-references (Ensembl, Entrez, UniProt, OMIM, RefSeq). Example: gene_group_id '366' returns the 56 members of the 'Ubiquitin specific peptidases' family (USP1-USP54, USPL1, CYLD, PAN2). Find a gene_group_id from the 'gene_group_id' field of any gene record returned by HGNC_fetch_gene_by_symbol.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'HGNC_fetch_gene_family_members',
     'arguments': {'gene_group_id': '...'}},
)
```

## `HGNC_search_by_location` {#HGNC_search_by_location}

Searches HGNC for genes mapped to a cytogenetic location such as 17p13.1.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `location`
- **Keywords:** cytoband, chromosomal location, hgnc, genes in region
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `location` | string | yes | Cytogenetic band location. Examples: '17p13.1' (TP53 region), '17q21.31' (BRCA1 region), '7p11.2' (EGFR region), '13q14.2' (RB1 region). |

??? note "Upstream description"

    Search HGNC for genes at a specific chromosomal location. Returns all genes mapped to that cytogenetic band. Useful for identifying genes in a chromosomal region of interest. Example: '17p13.1' returns genes including TP53.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'HGNC_search_by_location',
     'arguments': {'location': '...'}},
)
```

## `HGNC_search_genes` {#HGNC_search_genes}

Searches HGNC human genes by symbol, name or alias with wildcard support.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** hgnc, gene search, alias, wildcard
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search query - gene symbol, partial name, or keyword. Supports wildcards (*). Examples: 'BRCA*', 'tumor protein', 'insulin receptor'. |
| `search_field` | string \| null | no | Field to search in. Options: 'symbol', 'name', 'alias_symbol', 'prev_symbol', 'entrez_id', 'ensembl_gene_id', 'uniprot_ids', 'location'. Default searches across symbol and name. |

??? note "Upstream description"

    Search for human genes in HGNC by symbol, name, or alias. Supports wildcards (e.g., 'BRCA*' finds BRCA1, BRCA2). Searches across approved symbols and names. Returns a list of matching genes with basic info. Example: 'kinase' searches gene names containing kinase.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'HGNC_search_genes', 'arguments': {'query': '...'}},
)
```

## `HPA_generic_search` {#HPA_generic_search}

Runs a custom Human Protein Atlas search query and returns selected columns.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `search_query`
- **Keywords:** human protein atlas, custom query, columns, gene search
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `search_query` | string | yes | Search term for the query |
| `columns` | string | no | Comma-separated list of columns to retrieve. Defaults to 'g,gs,gd'. Common options: 'g' (Gene), 'gs' (Synonym), 'e' (Ensembl), 'u' (UniProt), 'pe' (Protein existence), 's' (Subcellular location), 'scml' (Main location), 'scal' (Addl location), 'rnat' (RNA tissue specificity), 'rnatsm' (RNA tissue nTPM), 'rnablm' (RNA blood nTPM), 'rnascm' (RNA single cell nTPM). |
| `format` | string | no | Response format (json or tsv). Defaults to 'json'. |

??? note "Upstream description"

    Generic search tool for Human Protein Atlas. Allows custom search queries and retrieval of specific columns.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'HPA_generic_search',
     'arguments': {'search_query': '...'}},
)
```

## `HPA_get_biological_processes_by_gene` {#HPA_get_biological_processes_by_gene}

Gets biological process annotations for a gene from HPA, highlighting key processes.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_name`, `filter_processes`
- **Keywords:** human protein atlas, biological process, gene function
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_name` | string | yes | Gene name or gene symbol, e.g., 'TP53', 'CDK1', 'CASP3', etc. |
| `filter_processes` | boolean | yes | Whether to filter and focus on specific biological processes (apoptosis, cell cycle, transcription regulation, etc.), defaults to true. |

??? note "Upstream description"

    Get biological process information for a gene, with special focus on key processes like apoptosis, cell cycle, etc.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'HPA_get_biological_processes_by_gene',
     'arguments': {'gene_name': '...', 'filter_processes': False}},
)
```

## `HPA_get_comprehensive_gene_details_by_ensembl_id` {#HPA_get_comprehensive_gene_details_by_ensembl_id}

Gets full HPA gene page data by Ensembl ID: expression, antibodies, IHC and immunofluorescence images.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensembl_id`
- **Keywords:** human protein atlas, ensembl id, immunohistochemistry, antibody, images
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensembl_id` | string | yes | Ensembl Gene ID, e.g., 'ENSG00000064787' (BCAS1), 'ENSG00000141510' (TP53), etc. Usually obtained through HPA_search_genes_by_query tool. |
| `include_images` | boolean | no | Whether to include image URL information (immunofluorescence, cell line images, etc.), defaults to true. Default `True` |
| `include_antibodies` | boolean | no | Whether to include detailed antibody information (validation status, Western blot data, etc.), defaults to true. Default `True` |
| `include_expression` | boolean | no | Whether to include detailed expression data (tissue specificity, subcellular localization, etc.), defaults to true. Default `True` |

??? note "Upstream description"

    Get detailed in-depth information from gene page using Ensembl Gene ID, including image URLs, antibody data, protein expression, and comprehensive information. This is the core tool for retrieving all images (tissue immunohistochemistry, subcellular immunofluorescence).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'HPA_get_comprehensive_gene_details_by_ensembl_id',
     'arguments': {'ensembl_id': '...'}},
)
```

## `HPA_get_gene_basic_info_by_ensembl_id` {#HPA_get_gene_basic_info_by_ensembl_id}

Gets basic gene info and expression data from the Human Protein Atlas by Ensembl ID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensembl_id`
- **Keywords:** human protein atlas, ensembl id, gene info
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensembl_id` | string | yes | Ensembl Gene ID, e.g., 'ENSG00000134057', 'ENSG00000141510', etc. |

??? note "Upstream description"

    Get gene basic information and expression data from Human Protein Atlas using Ensembl Gene ID. Enhanced version now uses efficient JSON API.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'HPA_get_gene_basic_info_by_ensembl_id',
     'arguments': {'ensembl_id': '...'}},
)
```

## `HPA_get_gene_tsv_data_by_ensembl_id` {#HPA_get_gene_tsv_data_by_ensembl_id}

Gets detailed HPA gene data in TSV format by Ensembl ID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensembl_id`
- **Keywords:** human protein atlas, tsv, ensembl id
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensembl_id` | string | yes | Ensembl Gene ID, e.g., 'ENSG00000134057', 'ENSG00000141510', etc. |

??? note "Upstream description"

    Get detailed gene data in TSV format from Human Protein Atlas using Ensembl Gene ID (backward compatibility tool).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'HPA_get_gene_tsv_data_by_ensembl_id',
     'arguments': {'ensembl_id': '...'}},
)
```

## `HPA_search_genes_by_query` {#HPA_search_genes_by_query}

Searches HPA by gene name, keyword or cell line and returns Ensembl gene IDs.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `search_query`
- **Keywords:** human protein atlas, gene search, ensembl id
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `search_query` | string | yes | Gene name, alias, keyword, or cell line name to search for, e.g., 'EGFR', 'TP53', or 'MCF7'. |
| `max_results` | integer \| null | no | Maximum number of gene matches to return (default: 50). This search does a broad full-text match with no server-side limit, so a short/common query can otherwise return thousands of results; exact gene-symbol matches are ranked first. |

??? note "Upstream description"

    Search for matching genes by gene name, keywords, or cell line names and return Ensembl ID list. This is the entry tool for many HPA query workflows.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'HPA_search_genes_by_query',
     'arguments': {'search_query': '...'}},
)
```

## `humanbase_ppi_analysis` {#humanbase_ppi_analysis}

Retrieves tissue-specific interaction networks and associated GO biological processes for a gene list from HumanBase.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** tissue-specific network, gene network, go process, functional interaction
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_list` | array of string | no | List of gene names or symbols to analyze for protein-protein interactions. The gene name should be the official gene symbol, not the synonym. |
| `tissue` | string | no | Tissue type for tissue-specific interactions. Use simple names like 'brain', 'heart', 'liver', 'kidney', 'lung', 'blood'. For v3 networks with more cell-type specificity, append '-v3' (e.g., 'brain-v3', 'liver-v3'). Default `brain` |
| `max_node` | integer | no | Maximum number of nodes to retrieve in the interaction network. Warning: the more nodes, the more time it takes to retrieve the data. Default is 10 (~30 seconds). Default `10` |
| `interaction` | string | no | Specific interaction type to filter by. Available types: 'co-expression', 'interaction', 'tf-binding', 'gsea-microrna-targets', 'gsea-perturbations'. If not specified, all types will be included. One of `co-expression`, `interaction`, `tf-binding`, `gsea-microrna-targets`, `gsea-perturbations` |
| `string_mode` | boolean | no | Whether to return the result in string mode. If True, the result will be a string of the network graph and the biological processes. If False, the result will be a NetworkX graph and a list of biological processes. Default `True` |
| `genes` | array of string | no | Alias for gene_list. List of gene symbols to analyze. |
| `top_n` | integer | no | Alias for max_node. Maximum number of network nodes to retrieve. |

??? note "Upstream description"

    Retrieve tissue-specific protein-protein interactions and biological processes from HumanBase. Returns a NetworkX graph of tissue specific protein-protein interactions and a list of associated biological processes involeed by the given genes from Gene Ontology.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'humanbase_ppi_analysis', 'arguments': {}},
)
```

## `HumanMine_search` {#HumanMine_search}

Searches HumanMine (human, mouse, rat) across genes, proteins, pathways and diseases with category facets.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `q`
- **Keywords:** humanmine, intermine, genomics warehouse, keyword search
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `q` | string | yes | Search query: gene symbol (e.g., 'TP53', 'BRCA1'), gene name, disease name, pathway name, or any biological term |
| `size` | integer \| null | no | Maximum number of results to return (default 10, max 100) |
| `facet_Category` | string \| null | no | Filter by entity type: 'Gene', 'Protein', 'Pathway', 'Disease', 'Publication', 'GWASResult', 'GOTerm', 'ProteinDomain', etc. |

??? note "Upstream description"

    Search HumanMine (InterMine) for genes, proteins, pathways, diseases, and other biological entities related to human, mouse, and rat. HumanMine integrates data from 30+ sources including NCBI Gene, UniProt, Ensembl, GO, KEGG, Reactome, OMIM, and GWAS Catalog. Returns entity IDs, types, names, symbols, and organism information. Use the facet_Category parameter to filter results to a specific type (Gene, Protein, Pathway, Disease, Publication, etc.).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'HumanMine_search', 'arguments': {'q': '...'}},
)
```

## `HumanMine_search_genes` {#HumanMine_search_genes}

Searches HumanMine genes by symbol, name or ID for human, mouse and rat.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `q`
- **Keywords:** humanmine, gene search, entrez, ensembl
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `q` | string | yes | Gene search query: gene symbol (e.g., 'EGFR', 'TP53'), gene name (e.g., 'tumor protein p53'), or Entrez Gene ID |
| `size` | integer \| null | no | Maximum number of results to return (default 10, max 100) |

??? note "Upstream description"

    Search HumanMine specifically for genes by symbol, name, or ID across human, mouse, and rat. HumanMine integrates gene data from NCBI Gene, Ensembl, UniProt, and other sources. Returns gene symbols, Entrez IDs, Ensembl IDs, gene names, and organism information. More focused than HumanMine_search as it pre-filters to Gene category only.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'HumanMine_search_genes', 'arguments': {'q': '...'}},
)
```

## `IMGT_search_genes` {#IMGT_search_genes}

Searches IMGT immunoglobulin/TCR germline genes by type (IGHV, TRAV, etc.) and species.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** imgt, ighv, trav, germline genes, immunogenetics
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Default `search_genes` |
| `query` | string | no | Search query (gene name) |
| `gene_type` | string | no | Gene type: IGHV, IGHD, IGHJ, IGKV, IGLV, TRAV, TRBV, etc. |
| `species` | string | no | Species name (default: Homo sapiens) Default `Homo sapiens` |

??? note "Upstream description"

    Search IMGT for immunoglobulin/TCR germline genes. Filter by gene type (IGHV, IGKV, TRAV, etc.) and species. Returns search URLs and gene type information for navigating IMGT GENE-DB.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'IMGT_search_genes', 'arguments': {}},
)
```

## `IMPC_search_genes` {#IMPC_search_genes}

Searches IMPC mouse genes by symbol, name, synonym or MGI ID with phenotyping status.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** impc, mgi id, mouse gene, search
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search query: gene symbol, name fragment, or MGI ID (e.g., 'Wdr7', 'WD repeat', 'MGI:2442934') |
| `limit` | integer | no | Maximum results (default: 20) Default `20` |

??? note "Upstream description"

    Search IMPC for mouse genes by symbol, name, synonym, or MGI ID. Returns matching genes with their phenotyping and production status. Use to find the correct MGI ID for a gene before querying phenotypes.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'IMPC_search_genes', 'arguments': {'query': '...'}},
)
```

## `InterMine_run_pathquery` {#InterMine_run_pathquery}

Runs an InterMine XML PathQuery against HumanMine and returns rows (e.g. gene pathways, domains, region features).

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** intermine, pathquery, humanmine, genomic model, xml query
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | InterMine XML PathQuery. Must include model, view (space-separated paths), and constraints. Example: '&lt;query model="genomic" view="Gene.symbol Gene.pathways.name"&gt;&lt;constraint path="Gene.symbol" op="=" value="PAX6"/&gt;&lt;/query&gt;'. |

??? note "Upstream description"

    Run an arbitrary InterMine PathQuery against HumanMine (an InterMine data warehouse) to traverse the data model graph: gene-&gt;pathways, gene-&gt;protein domains, region-&gt;features, etc. The 'query' argument is an InterMine XML PathQuery string (with model, view columns, and constraints) — results come back as JSON rows in the order of the 'view' columns. Example: query '&lt;query model="genomic" view="Gene.symbol Gene.pathways.name Gene.pathways.dataSets.name"&gt;&lt;constraint path="Gene.symbol" op="=" value="PAX6"/&gt;&lt;/query&gt;' returns PAX6 Reactome pathways like 'Activation of HOX genes during differentiation'. A second view 'Gene.proteins.proteinDomainRegions.proteinDomain.name' returns protein domains such as 'Homeobox domain'.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'InterMine_run_pathquery', 'arguments': {'query': '...'}},
)
```

## `KEGG_convert_ids` {#KEGG_convert_ids}

Converts KEGG gene or compound identifiers to UniProt, NCBI Gene, NCBI Protein, ChEBI or PubChem IDs.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `kegg_id`, `target_db`
- **Keywords:** id conversion, ncbi gene, uniprot, cross-reference
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `kegg_id` | string | yes | KEGG entry ID. For genes use format 'organism:gene_id' (e.g., 'hsa:7157' for human TP53, 'mmu:22059' for mouse Trp53). For compounds use 'cpd:C00002' (ATP). Multiple IDs separated by '+' (e.g., 'hsa:7157+hsa:672'). |
| `target_db` | string | yes | External database to convert to. Options: 'uniprot' (UniProt accessions), 'ncbi-geneid' (NCBI Gene/Entrez IDs), 'ncbi-proteinid' (NCBI protein accessions), 'chebi' (ChEBI compound IDs), 'pubchem' (PubChem substance IDs). |

??? note "Upstream description"

    Convert KEGG gene/compound identifiers to external database IDs using the KEGG /conv API. Maps KEGG gene IDs (e.g., hsa:7157 for TP53) to UniProt, NCBI Gene, NCBI Protein, ChEBI, or PubChem IDs. Essential for cross-database integration in bioinformatics workflows. Example: converting hsa:7157 to UniProt returns P04637 (TP53_HUMAN).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'KEGG_convert_ids',
     'arguments': {'kegg_id': '...', 'target_db': '...'}},
)
```

## `kegg_find_genes` {#kegg_find_genes}

Finds KEGG genes by keyword across all organisms or within one organism.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `keyword`
- **Keywords:** kegg gene, gene search, organism code
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `keyword` | string | yes | Search keyword for gene names or descriptions |
| `organism` | string | no | Organism code (e.g., 'hsa' for human, 'mmu' for mouse). Optional - searches all organisms if not specified Default `` |

??? note "Upstream description"

    Find genes in KEGG database by keyword. Can search across all organisms or within a specific organism.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'kegg_find_genes', 'arguments': {'keyword': '...'}},
)
```

## `kegg_get_gene_info` {#kegg_get_gene_info}

Returns KEGG gene information by gene ID including sequence, function and pathway associations.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** kegg gene id, gene function, nucleotide sequence, pathways
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | KEGG gene identifier (e.g., 'hsa:348', 'hsa:3480') |

??? note "Upstream description"

    Get detailed gene information from KEGG by gene ID. Returns gene data including sequence, function, and pathway associations.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'kegg_get_gene_info', 'arguments': {'gene_id': '...'}},
)
```

## `KEGG_get_gene_pathways` {#KEGG_get_gene_pathways}

Returns all KEGG pathways a gene participates in for a KEGG gene ID or symbol.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** gene pathways, kegg, tp53, pathway membership
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | no | KEGG gene identifier in organism:id format. Examples: 'hsa:7157' (human TP53), 'hsa:672' (human BRCA1), 'mmu:22059' (mouse Trp53). Use organism prefix: hsa=human, mmu=mouse, rno=rat, dme=fly, sce=yeast. |
| `gene_symbol` | string | no | Alias for gene_id: human gene symbol (e.g. 'TP53', 'BRCA1'). Auto-resolved to KEGG gene ID using organism='hsa'. |
| `gene` | string | no | Alias for gene_id: human gene symbol (e.g. 'TP53', 'EGFR'). |
| `organism` | string | no | KEGG organism code for gene symbol resolution (default: 'hsa' for human). Used only when gene_symbol is provided instead of gene_id. Default `hsa` |

??? note "Upstream description"

    Get all KEGG pathways that a gene participates in. Returns pathway IDs and names for a given KEGG gene identifier. Example: hsa:7157 (human TP53) participates in 51 pathways including p53 signaling pathway (hsa04115), apoptosis (hsa04210), cell cycle (hsa04110), and 48 disease pathways (various cancers, viral infections). Requires KEGG organism:id format.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'KEGG_get_gene_pathways', 'arguments': {}},
)
```

## `KEGG_get_network` {#KEGG_get_network}

Returns a KEGG network entry with signaling definition, disease classes, drug targets and gene elements.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `network_id`
- **Keywords:** kegg network, signaling pathway, gene elements, disease
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `network_id` | string | yes | KEGG network ID in N\#\#\#\#\# format (e.g., 'N00001' for EGF-EGFR-RAS-ERK, 'N00010' for WNT signaling). |

??? note "Upstream description"

    Get detailed KEGG network information including signaling pathway definition (e.g., EGF -&gt; EGFR -&gt; RAS -&gt; RAF -&gt; MEK -&gt; ERK), associated disease classes, linked diseases, drug targets, and pathway elements. KEGG networks model specific molecular interaction cascades perturbed in disease. Example: N00001 returns EGF-EGFR-RAS-ERK signaling pathway linked to colorectal cancer, melanoma, breast cancer, and 7+ other cancer types. For N\#\#\#\#\#-style perturbed-network IDs, `genes` gives the gene symbols (with Entrez IDs and descriptions) that map the numeric IDs in `expanded`, and `pathways` cross-references the canonical KEGG PATHWAY map(s) (e.g. hsa04668) the network belongs to. For nt\#\#\#\#\#\#-style disease/map-class network IDs, `elements` instead lists the child N\#\#\#\#\# networks that compose it (`genes`/`pathways` are empty for that ID family).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'KEGG_get_network', 'arguments': {'network_id': '...'}},
)
```

## `KEGG_get_pathway_genes` {#KEGG_get_pathway_genes}

Returns gene IDs for all member genes of a KEGG pathway.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `pathway_id`
- **Keywords:** pathway genes, gene set, kegg pathway
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `pathway_id` | string | yes | KEGG pathway identifier (e.g., 'hsa04115' for p53 signaling pathway, 'hsa04110' for cell cycle, 'hsa04210' for apoptosis). Use organism prefix: hsa=human, mmu=mouse. |

??? note "Upstream description"

    Get all genes in a KEGG pathway. Returns gene IDs for all member genes of a specific pathway. Example: hsa04115 (p53 signaling pathway) contains genes including CDK2 (hsa:1017), CDK4 (hsa:1019), CDK6 (hsa:1021), CDKN1A (hsa:1026), CDKN2A (hsa:1029), and many more involved in p53-mediated cell cycle arrest, apoptosis, and senescence.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'KEGG_get_pathway_genes',
     'arguments': {'pathway_id': '...'}},
)
```

## `kegg_get_pathway_info` {#kegg_get_pathway_info}

Returns KEGG pathway information by pathway ID including genes, compounds and reactions.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `pathway_id`
- **Keywords:** kegg pathway, pathway genes, reactions
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `pathway_id` | string | yes | KEGG pathway identifier (e.g., 'hsa00010', 'path:hsa00010') |

??? note "Upstream description"

    Get detailed pathway information from KEGG by pathway ID. Returns pathway data including genes, compounds, and reactions.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'kegg_get_pathway_info',
     'arguments': {'pathway_id': '...'}},
)
```

## `KEGG_link_entries` {#KEGG_link_entries}

Finds linked entries between KEGG databases, e.g. diseases for a gene or genes in a pathway.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `source`, `target`
- **Keywords:** kegg link, cross-reference, gene-disease, pathway genes
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `source` | string | yes | KEGG entry ID to query. Examples: 'hsa:7157' (gene), 'hsa05200' (pathway), 'H00004' (disease), 'D00279' (drug), 'C00002' (compound). |
| `target` | string | yes | Target KEGG database to search for links. Options: 'pathway', 'disease', 'drug', 'compound', 'ko' (orthology), 'enzyme', 'reaction', 'module'. Can also use organism prefix for genes (e.g., 'hsa'). |

??? note "Upstream description"

    Find cross-references between KEGG databases using the KEGG /link API. Given a KEGG entry (gene, pathway, disease, drug, compound), find all linked entries in another KEGG database. Examples: find all diseases for gene hsa:7157 (TP53), all genes in a pathway, or all pathways for a compound. Complements KEGG_get_gene_pathways with more flexible cross-database queries.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'KEGG_link_entries',
     'arguments': {'source': '...', 'target': '...'}},
)
```

## `KEGG_search_network` {#KEGG_search_network}

Searches KEGG NETWORK for disease/drug-associated signaling networks by gene or keyword.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `keyword`
- **Keywords:** network search, signaling, perturbed pathway
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `keyword` | string | yes | Search keyword for network/pathway name or gene (e.g., 'EGFR', 'RAS', 'p53', 'WNT', 'NOTCH'). Case-insensitive. |
| `max_results` | integer | no | Maximum number of results (default: 25). Default `25` |

??? note "Upstream description"

    Search the KEGG NETWORK database for signaling pathway networks by gene or keyword. KEGG NETWORK contains 800+ disease/drug-associated molecular interaction networks representing signal transduction pathways perturbed in cancer, infection, and other diseases. Returns network IDs (N\#\#\#\#\#) and pathway names. Example: searching 'EGFR' returns N00001 EGF-EGFR-RAS-ERK signaling pathway, N00006 Amplified EGFR, etc.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'KEGG_search_network', 'arguments': {'keyword': '...'}},
)
```

## `LipidMaps_get_gene` {#LipidMaps_get_gene}

Resolves a lipid-metabolism gene in LMPD to its record with names, chromosome, cytogenetic location and summary.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `input_value`
- **Keywords:** lmpd, lipid metabolism gene, chromosome, cytogenetic band
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `input_value` | string | yes | The gene identifier value to resolve, e.g. 'FASN' (gene symbol), '2194' (NCBI Gene ID), or 'LMP002507' (LIPID MAPS protein id). |
| `input_item` | string | no | Namespace of input_value. One of: 'gene_symbol', 'gene_id', 'gene_name', 'lmp_id'. One of `gene_symbol`, `gene_id`, `gene_name`, `lmp_id` Default `gene_symbol` |
| `output_item` | string | no | Which fields to return. 'all' (default) returns the complete record; a single field name (e.g. 'gene_name', 'summary', 'chromosome') returns just that field. Default `all` |

??? note "Upstream description"

    Look up a lipid-metabolism gene in the LIPID MAPS Proteome Database (LMPD). Resolves a gene to its LMPD record (gene name, synonyms, alternative names, chromosome, cytogenetic map location, NCBI summary, and species). Set input_item to the namespace of input_value: 'gene_symbol' (e.g. 'FASN'), 'gene_id' (NCBI Gene ID, e.g. '2194'), 'gene_name' (exact gene-name substring), or 'lmp_id' (LIPID MAPS protein id, e.g. 'LMP002507'). A gene_symbol may match multiple species (e.g. 'FASN' returns mouse LMP001577 and human LMP002507); all matches are returned as a list. Use this for lipid-related enzymes such as fatty acid synthase (FASN) or stearoyl-CoA desaturase (SCD). Complements the LipidMaps_get_compound_* structure tools, which cover lipids, not the genes/enzymes that make them.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'LipidMaps_get_gene', 'arguments': {'input_value': '...'}},
)
```

## `LipidMaps_get_protein` {#LipidMaps_get_protein}

Resolves a lipid-metabolism protein in LMPD to gene symbol, species and UniProt/RefSeq cross-references.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `input_value`
- **Keywords:** lmpd, refseq, uniprot, gene symbol
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `input_value` | string | yes | The protein identifier value to resolve, e.g. 'P49327' (UniProt), 'NP_004095' (RefSeq protein), 'SCD' (gene symbol), or 'LMP002507' (LIPID MAPS protein id). |
| `input_item` | string | no | Namespace of input_value. One of: 'uniprot_id', 'gene_symbol', 'gene_id', 'gene_name', 'refseq_id', 'mrna_id', 'protein_entry', 'lmp_id'. One of `uniprot_id`, `gene_symbol`, `gene_id`, `gene_name`, `refseq_id`, `mrna_id`, `protein_entry`, `lmp_id` Default `uniprot_id` |
| `output_item` | string | no | Which fields to return. 'all' (default) returns the complete record; a single field name (e.g. 'uniprot_id', 'protein_name', 'seq') returns just that field. Default `all` |

??? note "Upstream description"

    Look up a lipid-metabolism protein/enzyme in the LIPID MAPS Proteome Database (LMPD). Resolves a protein to its LMPD record including gene symbol/name, species, and cross-references (UniProt accession, RefSeq protein id, mRNA RefSeq id, protein GI, UniProt entry name) plus the amino-acid sequence. Set input_item to the namespace of input_value: 'uniprot_id' (e.g. 'P49327'), 'gene_symbol' (e.g. 'SCD'), 'gene_id' (NCBI Gene ID), 'refseq_id' (RefSeq protein id, e.g. 'NP_004095'), 'mrna_id' (RefSeq mRNA id), 'protein_entry' (UniProt entry name, e.g. 'FAS_HUMAN'), or 'lmp_id' (LIPID MAPS protein id). A gene_symbol may match multiple species; all matches are returned as a list. Use this to map a lipid-related enzyme (e.g. fatty acid synthase, stearoyl-CoA desaturase) to its sequence and database cross-references.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'LipidMaps_get_protein',
     'arguments': {'input_value': '...'}},
)
```

## `LNCipedia_get_lncrna` {#LNCipedia_get_lncrna}

Returns lncRNA details from RNAcentral by URS ID including RNA sequence, species, gene associations and sources.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `rnacentral_id`
- **Keywords:** lncrna, rnacentral, urs id, ncrna sequence
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `rnacentral_id` | string | yes | RNAcentral URS identifier. Use species-specific format 'URS_TAXID' (e.g., 'URS0002A146E4_9606' for human MALAT1) or base URS ID with taxid parameter. |
| `taxid` | integer | no | NCBI Taxonomy ID for species-specific lookup (e.g., 9606 for human). Optional if rnacentral_id already includes taxid. |

??? note "Upstream description"

    Get detailed lncRNA information from RNAcentral including full RNA sequence, species, gene associations, publication count, and source databases. Use a species-specific RNAcentral ID (format: URS_TAXID) or base URS ID with optional taxid. Key lncRNA IDs: URS0002A146E4_9606 (MALAT1, 7472 nt), URS00026A23F2_9606 (HOTAIR, 2273 nt), URS000075DAEC_9606 (NEAT1, 3756 nt). Returns the full nucleotide sequence which can be used for downstream analysis (structure prediction, conservation, etc.).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'LNCipedia_get_lncrna',
     'arguments': {'rnacentral_id': '...'}},
)
```

## `LNCipedia_get_lncrna_xrefs` {#LNCipedia_get_lncrna_xrefs}

Returns cross-references for a lncRNA from LNCipedia, NONCODE, HGNC, RefSeq, Ensembl and other databases.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `rnacentral_id`
- **Keywords:** lncrna xref, hgnc, refseq, ensembl, noncode
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `rnacentral_id` | string | yes | RNAcentral URS identifier (e.g., 'URS0002A146E4' for MALAT1). Use base URS ID without taxid suffix. |
| `page_size` | integer | no | Number of cross-references to return (1-50). Default: 10. |

??? note "Upstream description"

    Get cross-references for a long non-coding RNA from RNAcentral. Returns database links from LNCipedia, NONCODE, LncBook, HGNC, RefSeq, GeneCards, MalaCards, Ensembl, and other sources. Includes external accession IDs, expert database URLs, and species information. Essential for mapping lncRNA IDs across databases and finding functional annotation data.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'LNCipedia_get_lncrna_xrefs',
     'arguments': {'rnacentral_id': '...'}},
)
```

## `LNCipedia_search_lncrna` {#LNCipedia_search_lncrna}

Searches lncRNAs by name, accession or keyword over RNAcentral with optional species filter.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** lncrna search, rnacentral, lncipedia, species
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search query: lncRNA name (e.g., 'MALAT1', 'HOTAIR'), gene symbol, accession, or keyword. Supports AND/OR operators. |
| `species` | string | no | Species name filter (e.g., 'Homo sapiens', 'Mus musculus'). Leave empty for all species. |
| `size` | integer | no | Number of results to return (1-100). Default: 10. |

??? note "Upstream description"

    Search for long non-coding RNAs (lncRNAs) by name, accession, or keyword using EBI Search over RNAcentral. RNAcentral aggregates lncRNA data from LNCipedia, NONCODE, LncBook, HGNC, RefSeq, GeneCards, and 50+ other databases. Supports filtering by species. Returns RNAcentral IDs, descriptions, sequence lengths, and source databases. Use to search for lncRNAs like MALAT1, HOTAIR, NEAT1, XIST, H19, or KCNQ1OT1. The human genome encodes 50,000+ lncRNA transcripts. Example: search 'MALAT1' returns URS0002A146E4_9606 (7472 nt, metastasis-associated lung adenocarcinoma transcript 1).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'LNCipedia_search_lncrna', 'arguments': {'query': '...'}},
)
```

## `LNCipedia_search_ncrna_by_type` {#LNCipedia_search_ncrna_by_type}

Searches RNAcentral non-coding RNAs by type such as miRNA, lncRNA, rRNA, tRNA or snoRNA, with optional species.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** ncrna, mirna, trna, snorna, rna type
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search query: gene name, keyword, or accession. Use '*' for all entries of the specified RNA type. |
| `rna_type` | string | no | RNA type filter: miRNA, lncRNA, rRNA, tRNA, snRNA, snoRNA, piRNA, pre_miRNA, SRP_RNA, RNase_P_RNA, antisense_RNA, ribozyme, etc. |
| `species` | string | no | Species name filter (e.g., 'Homo sapiens'). Leave empty for all species. |
| `size` | integer | no | Number of results to return (1-100). Default: 10. |

??? note "Upstream description"

    Search for any type of non-coding RNA by RNA type classification. Supports all RNA types in RNAcentral: miRNA, lncRNA, rRNA, tRNA, snRNA, snoRNA, piRNA, pre_miRNA, SRP_RNA, RNase_P_RNA, antisense_RNA, and more. Useful for browsing or discovering ncRNAs of a specific type, optionally filtered by species or keyword. Example: search for all human snoRNAs, or find piRNAs related to a specific gene.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'LNCipedia_search_ncrna_by_type',
     'arguments': {'query': '...'}},
)
```

## `LOVD_get_gene` {#LOVD_get_gene}

Returns LOVD gene information including HGNC/Entrez IDs, chromosomal location, RefSeq transcripts and build.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_symbol`
- **Keywords:** lovd, hgnc id, refseq transcript, chromosomal location
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_symbol` | string | yes | HGNC gene symbol (e.g., 'TP53', 'BRCA1', 'BRCA2', 'EGFR'). |

??? note "Upstream description"

    Get gene information from LOVD (Leiden Open Variation Database). Returns HGNC ID, Entrez ID, chromosomal location, RefSeq transcript(s), genomic build, curation details, and dates. LOVD hosts curated variant databases for over 23,000 genes used in clinical diagnostics.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'LOVD_get_gene', 'arguments': {'gene_symbol': '...'}},
)
```

## `MARRVEL_get_gene` {#MARRVEL_get_gene}

Returns MARRVEL aggregated gene identity for a human gene symbol with HGNC/Ensembl/Entrez IDs, location and summary.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** marrvel, gene identity, hgnc, cytogenetic location
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `symbol` | string \| null | no | HGNC gene symbol, e.g. 'CFTR', 'BRCA1'. |
| `gene_symbol` | string \| null | no | Alias for symbol. |
| `gene` | string \| null | no | Alias for symbol. |

??? note "Upstream description"

    Get aggregated identity and annotation for a human gene by symbol from MARRVEL: cross-references (HGNC, Ensembl, Entrez, UniProt), chromosome/cytogenetic location, gene type, aliases, previous symbols, and the Entrez gene summary. Use as a one-call gene-identity resolver in rare-disease / Mendelian variant triage workflows. Does NOT include an OMIM ID -- MARRVEL's own gene-lookup endpoint returns an unreliable placeholder value for it (confirmed live), use MARRVEL_get_omim_phenotypes for OMIM associations instead. MARRVEL indexes only current approved HGNC symbols, so a retired symbol or alias (e.g. 'RP20'/'LCA2' for RPE65) is resolved once via HGNC and re-queried; when that happens the substitution is disclosed in metadata.resolved_symbol and metadata.resolved_note, and metadata.query_symbol keeps the symbol you asked for. No API key required.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MARRVEL_get_gene', 'arguments': {}},
)
```

## `MedlinePlus_get_genetics_gene_by_name` {#MedlinePlus_get_genetics_gene_by_name}

Returns MedlinePlus Genetics information for a gene including function, related conditions and synonyms.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene`
- **Keywords:** gene function, medlineplus genetics, gene summary
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string | yes | URL slug of gene name, e.g., "BRCA1", must match MedlinePlus page path. |
| `format` | string | no | Format parameter for the endpoint (note: API always returns XML regardless of this parameter). One of `json`, `xml` Default `json` |

??? note "Upstream description"

    Get detailed information from MedlinePlus Genetics corresponding to gene name. Returns XML data parsed into a dictionary with gene function, related health conditions, synonyms, and other details.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MedlinePlus_get_genetics_gene_by_name',
     'arguments': {'gene': '...'}},
)
```

## `MetabolomicsWorkbench_get_gene_protein` {#MetabolomicsWorkbench_get_gene_protein}

Looks up a Metabolomics Workbench gene or protein (MGP) annotation record by gene or protein identifier.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `input_value`
- **Keywords:** mgp, gene record, metabolic gene, entrez
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `input_value` | string | yes | The identifier to look up (e.g. gene symbol 'acaca', UniProt id 'Q13085', mgp_id 'MGP000016'). |
| `entity` | string \| null | no | Which record type to fetch: 'gene' (default) or 'protein'. |
| `id_type` | string \| null | no | Namespace of input_value. For genes: gene_symbol (default), gene_id, gene_name, mgp_id, refseq_id. For proteins: uniprot_id (default), gene_symbol, mgp_id, refseq_id, gene_id. |

??? note "Upstream description"

    Look up a Metabolomics Workbench gene or protein (MGP) annotation record linking gene/protein identifiers to the Metabolomics Workbench. Set entity='gene' (default) or entity='protein'; id_type selects the lookup namespace. Gene records (entity='gene') return gene_symbol, gene_id, gene_name, synonyms, chromosome/map_location, MGP id, and a summary. Protein records (entity='protein') return UniProt id, gene_symbol, species, mRNA/RefSeq ids, protein entry name, sequence length and the amino-acid sequence. Examples: entity='gene', id_type='gene_symbol', input_value='acaca' -&gt; ACACA 'acetyl-CoA carboxylase alpha', gene_id 31, MGP000016; entity='protein', id_type='uniprot_id', input_value='Q13085' -&gt; ACACA_HUMAN, NM_198839, NP_942136, 2346 aa.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MetabolomicsWorkbench_get_gene_protein',
     'arguments': {'input_value': '...'}},
)
```

## `MGI_get_gene` {#MGI_get_gene}

Returns MGI mouse gene details by MGI ID including synopsis, synonyms, GRCm39 location and cross-references.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** mgi, mouse gene, grcm39, gene synopsis
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | MGI gene identifier with 'MGI:' prefix. Examples: 'MGI:98834' (Trp53/p53), 'MGI:104537' (Brca1), 'MGI:96677' (Kit), 'MGI:88039' (Apoe), 'MGI:88276' (Braf). |

??? note "Upstream description"

    Get detailed information about a mouse gene from MGI (Mouse Genome Informatics) by its MGI identifier. Returns gene symbol, full name, species info, gene synopsis (curated and automated), synonyms, genomic location (GRCm39 assembly), SO term, and cross-references to UniProt, NCBI, Ensembl, etc. Powered by the Alliance of Genome Resources API which integrates MGI curated data. Example: 'MGI:98834' returns Trp53 (transformation related protein 53) details.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MGI_get_gene', 'arguments': {'gene_id': '...'}},
)
```

## `MGI_search_genes` {#MGI_search_genes}

Searches MGI mouse genes by symbol, name or keyword and returns MGI IDs.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** mgi, mouse gene search, mus musculus
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Gene search query: mouse gene symbol (e.g., 'Trp53', 'Brca1', 'Sox9'), gene name keyword (e.g., 'insulin', 'tumor'), or MGI ID (e.g., 'MGI:98834'). |
| `limit` | integer | no | Maximum number of results to return (default: 10, max: 50). Default `10` |

??? note "Upstream description"

    Search for mouse (Mus musculus) genes in MGI (Mouse Genome Informatics) by gene symbol, name, or keyword. Returns matching mouse gene symbols, names, and MGI identifiers. Powered by the Alliance of Genome Resources API which integrates MGI data. Use this to find mouse gene IDs before calling MGI_get_gene or MGI_get_phenotypes. Examples: 'Trp53' (tumor protein p53), 'Brca1', 'insulin', 'Notch'.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MGI_search_genes', 'arguments': {'query': '...'}},
)
```

## `MGnify_get_genome` {#MGnify_get_genome}

Returns MGnify genome catalog details with completeness, contamination, N50, taxonomy and pangenome statistics.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `genome_id`
- **Keywords:** mag, genome assembly, n50, completeness, pangenome
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `genome_id` | string | yes | MGnify genome accession. Examples: 'MGYG000000001', 'MGYG000535630'. |

??? note "Upstream description"

    Get detailed genome information from the MGnify genome catalog. Returns assembly quality metrics (completeness, contamination, N50), taxonomy, geographic origin, rRNA counts, protein count, pangenome statistics, and cross-references to ENA. Example: MGYG000000001 returns a high-quality Firmicutes isolate genome with 98.6% completeness.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MGnify_get_genome', 'arguments': {'genome_id': '...'}},
)
```

## `MGnify_get_go_terms` {#MGnify_get_go_terms}

Returns GO term annotations with gene counts from a MGnify metagenomics analysis.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `analysis_id`
- **Keywords:** go terms, gene counts, functional profile, metagenome
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `analysis_id` | string | yes | MGnify analysis accession (e.g., 'MGYA00585482'). Find IDs via MGnify_list_analyses. |
| `page_size` | integer | no | Number of GO terms to return per page (default 25, max 100). Default `25` |

??? note "Upstream description"

    Get Gene Ontology (GO) functional annotations from a MGnify metagenomics analysis. Returns GO terms with their accession, description, gene count, and category (biological_process, molecular_function, or cellular_component). Use to understand the functional capacity of a microbial community. Example: MGYA00585482 returns GO:0055114 oxidation-reduction process (count=1507).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MGnify_get_go_terms',
     'arguments': {'analysis_id': '...'}},
)
```

## `MGnify_search_genomes` {#MGnify_search_genomes}

Searches MGnify genome catalog for MAGs and isolate genomes by taxonomy lineage.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** mag, isolate genome, taxonomy lineage, reference genome
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `taxonomy` | string | no | GTDB taxonomy lineage filter, matched anywhere in the lineage. Examples: 'Bacteroidota' (phylum), 'Bacteroides' (genus), 'Bacillota', 'Pseudomonadota'. NCBI-style phylum names such as 'Bacteroidetes' or 'Firmicutes' return few or no hits -- use the GTDB spelling. |
| `page` | integer | no | Page number (default 1). |
| `page_size` | integer | no | Results per page (default 25, max 100). |

??? note "Upstream description"

    Search the MGnify genome catalog for metagenome-assembled genomes (MAGs) and isolate genomes. Filter by taxonomy lineage. Returns quality metrics, taxonomy, and protein counts. Useful for finding reference genomes from specific microbial lineages. Taxonomy names follow GTDB, not NCBI: use 'Bacteroidota' (not 'Bacteroidetes'), 'Bacillota' (not 'Firmicutes'), 'Pseudomonadota' (not 'Proteobacteria') -- an NCBI-style name usually returns 0 hits rather than an error.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MGnify_search_genomes', 'arguments': {}},
)
```

## `MIBiG_list_clusters` {#MIBiG_list_clusters}

Lists MIBiG biosynthetic gene clusters matching a search term.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** biosynthetic gene cluster, bgc, natural products, mibig
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `search_term` | string | no | Optional text to filter results by product name or organism. The API returns all entries; filtering is applied client-side. Examples: 'erythromycin', 'vancomycin', 'Streptomyces', 'PKS'. Leave empty to get all 3,013 BGC entries. |

??? note "Upstream description"

    List biosynthetic gene clusters (BGCs) from the MIBiG repository. MIBiG (Minimum Information about a Biosynthetic Gene Cluster) is the standard repository for experimentally characterized BGCs producing natural products like antibiotics, antifungals, and anticancer compounds. Returns all BGC entries with accession, products, biosynthetic class (PKS, NRPS, Terpene, etc.), organism, quality, and completeness. Contains 3,013 BGC entries from bacteria, fungi, and plants. Use search_term to filter results client-side by product name or organism.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MIBiG_list_clusters', 'arguments': {}},
)
```

## `miRBase_get_mirna` {#miRBase_get_mirna}

Returns miRNA details from RNAcentral including sequence, species, RNA type and gene associations.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `rnacentral_id`
- **Keywords:** mirna, rnacentral, mature sequence, urs id
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `rnacentral_id` | string | yes | RNAcentral URS identifier. Use species-specific format 'URS_TAXID' (e.g., 'URS000039ED8D_9606') or base URS ID (e.g., 'URS000039ED8D') with taxid parameter. |
| `taxid` | integer | no | NCBI Taxonomy ID for species-specific lookup (e.g., 9606 for human, 10090 for mouse). Optional if rnacentral_id already includes taxid. |

??? note "Upstream description"

    Get detailed miRNA information from RNAcentral including full RNA sequence, species, RNA type, gene associations, publication count, and source databases. Use a species-specific RNAcentral ID (format: URS_TAXID, e.g., 'URS000039ED8D_9606' for human hsa-miR-21-5p) or a base URS ID with optional taxid parameter. Key miRNA IDs: URS000039ED8D_9606 (hsa-miR-21-5p), URS0000D54CAD_9606 (hsa-miR-155-5p), URS0000416056_9606 (hsa-let-7a-5p), URS000039ED8D_10090 (mmu-miR-21a-5p).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'miRBase_get_mirna',
     'arguments': {'rnacentral_id': '...'}},
)
```

## `miRBase_get_mirna_xrefs` {#miRBase_get_mirna_xrefs}

Returns miRNA cross-references from miRBase, TarBase, ENA, RefSeq and other sources.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `rnacentral_id`
- **Keywords:** mirbase, tarbase, xref, refseq
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `rnacentral_id` | string | yes | RNAcentral URS identifier (e.g., 'URS000039ED8D' for hsa-miR-21-5p). Use base URS ID without taxid suffix. |
| `page_size` | integer | no | Number of cross-references to return (1-50). Default: 10. |

??? note "Upstream description"

    Get cross-references (xrefs) for a microRNA from RNAcentral. Returns database links from miRBase, TarBase, ENA, RefSeq, piRBase, LncBase, and other sources. Includes external accession IDs, expert database URLs, species information, and miRBase-specific fields (mature products, precursor links). Essential for mapping miRNA IDs across databases and finding target interaction data.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'miRBase_get_mirna_xrefs',
     'arguments': {'rnacentral_id': '...'}},
)
```

## `miRBase_search_mirna` {#miRBase_search_mirna}

Searches miRNAs by name, accession or keyword over RNAcentral with species filter.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** mirna search, mirbase, species
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search query: miRNA name (e.g., 'hsa-miR-21-5p'), miRBase accession (e.g., 'MIMAT0000076'), gene name, or keyword. Supports AND/OR operators. |
| `species` | string | no | Species name filter (e.g., 'Homo sapiens', 'Mus musculus', 'Drosophila melanogaster'). Leave empty for all species. |
| `size` | integer | no | Number of results to return (1-100). Default: 10. |

??? note "Upstream description"

    Search for microRNAs (miRNAs) by name, accession, or keyword using EBI Search over RNAcentral. RNAcentral aggregates miRNA data from miRBase, TarBase, LncBase, and 50+ other databases. Supports filtering by species (e.g., 'Homo sapiens', 'Mus musculus') and RNA type. Returns RNAcentral IDs, descriptions, sequence lengths, and source databases. Use this to find miRNAs like hsa-miR-21-5p, let-7, miR-155, or search by miRBase accession (e.g., MIMAT0000076). Examples: search 'hsa-miR-21' for human miR-21 family, 'let-7 AND TAXONOMY:Homo sapiens' for human let-7 family members.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'miRBase_search_mirna', 'arguments': {'query': '...'}},
)
```

## `Monarch_search_gene` {#Monarch_search_gene}

Searches Monarch for a gene by symbol or name and returns HGNC/NCBIGene CURIEs and orthologs.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** gene curie, orthologs, monarch, gene search
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Gene symbol or name (e.g., 'TP53', 'BRCA1', 'WDR7') |
| `limit` | integer | no | Number of results to return (default: 10) Default `10` |
| `offset` | integer | no | Number of initial entries to skip. Default `0` |

??? note "Upstream description"

    Search Monarch Initiative knowledge graph for a gene by symbol or name. Returns HGNC/NCBIGene CURIEs, cross-species orthologs, and gene metadata. Use this first to get the gene CURIE needed for Monarch_get_gene_phenotypes and Monarch_get_gene_diseases.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Monarch_search_gene', 'arguments': {'query': '...'}},
)
```

## `MouseMine_search` {#MouseMine_search}

Searches MouseMine (MGI) records across genes, alleles, proteins, GO terms and publications.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `q`
- **Keywords:** mousemine, mgi, intermine, mouse genomics
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `q` | string | yes | Search query: gene symbol (e.g., 'Brca1', 'Tp53', 'Sox9'), protein name, disease term, or free text. |
| `size` | integer | no | Number of results to return (default: 10, max: 100). Default `10` |
| `format` | string | no | Response format. Must be 'json'. Default `json` |
| `facet_Category` | string | no | Filter results by entity category. Options include: 'ProteinCodingGene', 'Gene', 'Allele', 'Protein', 'Publication', 'Genotype', 'GOTerm', 'ProteinDomain', 'GXDExpression', 'NcRNAGene', 'AntisenseLncRNAGene', 'Strain', 'HTExperiment', 'InteractionExperiment', 'DOTerm', 'UniProtFeature'. |

??? note "Upstream description"

    Search MouseMine, the InterMine data warehouse for mouse (Mus musculus) genomics data from MGI (Mouse Genome Informatics). Searches across genes, alleles, proteins, publications, gene ontology terms, protein domains, and more. Returns matching records with category facets and relevance scores. MouseMine integrates data from MGI including gene annotations, phenotypes, alleles, expression data, pathways, and disease models. Use facet_Category parameter to filter by specific entity type (e.g., 'ProteinCodingGene', 'Allele', 'Publication', 'Protein').

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MouseMine_search', 'arguments': {'q': '...'}},
)
```

## `MouseMine_search_genes` {#MouseMine_search_genes}

Searches MouseMine for mouse protein-coding genes and returns symbols and MGI IDs.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `q`
- **Keywords:** mouse genes, mgi id, protein-coding
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `q` | string | yes | Gene search query: gene symbol (e.g., 'Brca1', 'Tp53'), gene name, or keyword. |
| `size` | integer | no | Number of results to return (default: 10). Default `10` |

??? note "Upstream description"

    Search MouseMine specifically for mouse protein-coding genes. This is a convenience wrapper that filters results to ProteinCodingGene category only. Returns gene symbols, MGI identifiers, and organism information from Mouse Genome Informatics (MGI). For broader searches across all entity types, use MouseMine_search instead.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MouseMine_search_genes', 'arguments': {'q': '...'}},
)
```

## `MSigDB_check_gene_in_set` {#MSigDB_check_gene_in_set}

Checks whether a gene belongs to an MSigDB gene set such as TF target or hallmark sets.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_set_name`, `gene`
- **Keywords:** gene set membership, msigdb, tf targets, mirna targets
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | any | no |  |
| `gene_set_name` | string | yes | Exact gene set name (e.g., 'ZNF549_TARGET_GENES', 'ESC_V6.5_UP_EARLY.V1_DN', 'HALLMARK_APOPTOSIS') |
| `gene` | string | yes | Gene symbol to check (e.g., 'SELENOP', 'TP53', 'CD37') |
| `species` | string | no | MSigDB collection to query. 'auto' (default) tries human then mouse. Use 'mouse' for mouse-only sets such as the MP_* mammalian-phenotype sets and other MH/M1-M8 collections. One of `auto`, `human`, `mouse` Default `auto` |

??? note "Upstream description"

    Check if a specific gene is a member of an MSigDB gene set. Covers GTRD transcription factor targets (e.g., ZNF549_TARGET_GENES), miRDB miRNA targets, oncogenic signatures (C6), hallmark sets (H), and 33,000+ other sets. Use to verify gene membership in TF target lists, miRNA target lists, or any MSigDB collection. Supports human and mouse collections (species parameter).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MSigDB_check_gene_in_set',
     'arguments': {'gene_set_name': '...', 'gene': '...'}},
)
```

## `MSigDB_get_gene_set_members` {#MSigDB_get_gene_set_members}

Returns the parsed gene list of an MSigDB gene set by exact name.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_set_name`
- **Keywords:** gene set, msigdb, gsea, members
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | any | no |  |
| `gene_set_name` | string | yes | Exact MSigDB gene set name (e.g., 'ZNF549_TARGET_GENES', 'MIR675_3P_TARGET_GENES', 'ESC_V6.5_UP_EARLY.V1_DN') |
| `species` | string | no | MSigDB collection to query. 'auto' (default) tries human then mouse. Use 'mouse' for mouse-only sets such as the MP_* mammalian-phenotype sets and other MH/M1-M8 collections. One of `auto`, `human`, `mouse` Default `auto` |

??? note "Upstream description"

    Get all gene members of an MSigDB gene set by exact name. Returns parsed gene list (unlike MSigDB_get_geneset which returns raw API response). Covers 33,000+ gene sets including GTRD TF targets, miRDB miRNA targets, oncogenic signatures, hallmark sets, GO terms, and pathway gene sets. Supports human and mouse collections (species parameter).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MSigDB_get_gene_set_members',
     'arguments': {'gene_set_name': '...'}},
)
```

## `MSigDB_get_geneset` {#MSigDB_get_geneset}

Returns an MSigDB gene set by exact name from the API.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `geneSetName`
- **Keywords:** msigdb, gsea, gene set collection
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `geneSetName` | string | yes | Exact MSigDB gene set name (case-sensitive). Examples: 'HALLMARK_APOPTOSIS', 'KEGG_P53_SIGNALING_PATHWAY', 'GOBP_APOPTOTIC_PROCESS', 'HALLMARK_HYPOXIA', 'KRAS.LUNG_UP.V1_UP'. Find names at https://www.gsea-msigdb.org/gsea/msigdb/ |
| `species` | string | no | MSigDB collection to query. 'auto' (default) tries human then mouse. Use 'mouse' for mouse-only sets such as the MP_* mammalian-phenotype sets and other MH/M1-M8 collections. One of `auto`, `human`, `mouse` Default `auto` |

??? note "Upstream description"

    Retrieve a specific gene set from MSigDB (Molecular Signatures Database) by its exact name. MSigDB is the definitive database for gene sets used in Gene Set Enrichment Analysis (GSEA), containing 33,000+ gene sets organized into 9 major collections: H (hallmark), C1 (positional), C2 (curated pathways from KEGG/Reactome/BioCarta/WikiPathways), C3 (regulatory motifs), C4 (computational cancer), C5 (Gene Ontology), C6 (oncogenic signatures), C7 (immunologic signatures), C8 (cell type signatures). Returns the gene symbols in the set, systematic name, PMID, collection, and URLs. Supports human and mouse collections (species parameter).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MSigDB_get_geneset', 'arguments': {'geneSetName': '...'}},
)
```

## `MSigDB_get_hallmark_geneset` {#MSigDB_get_hallmark_geneset}

Returns an MSigDB Hallmark gene set by pathway name.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `pathway_name`
- **Keywords:** hallmark, gene set, gsea
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `pathway_name` | string | yes | Hallmark pathway name without the 'HALLMARK_' prefix (will be auto-prepended). Examples: 'APOPTOSIS', 'P53_PATHWAY', 'HYPOXIA', 'DNA_REPAIR', 'INFLAMMATORY_RESPONSE', 'EPITHELIAL_MESENCHYMAL_TRANSITION', 'TNFA_SIGNALING_VIA_NFKB', 'KRAS_SIGNALING_UP', 'OXIDATIVE_PHOSPHORYLATION', 'GLYCOLYSIS', 'MYC_TARGETS_V1' |
| `species` | string | no | MSigDB collection to query. 'auto' (default) tries human then mouse. Use 'mouse' for mouse-only sets such as the MP_* mammalian-phenotype sets and other MH/M1-M8 collections. One of `auto`, `human`, `mouse` Default `auto` |

??? note "Upstream description"

    Retrieve a Hallmark gene set from MSigDB. The Hallmark collection (H) contains 50 gene sets representing well-defined biological states and processes with coherent expression, distilled from multiple founder sets. These are the most commonly used gene sets for pathway analysis and GSEA. Hallmark gene set names are prefixed with 'HALLMARK_' followed by the pathway name in uppercase with underscores (e.g., HALLMARK_APOPTOSIS, HALLMARK_HYPOXIA, HALLMARK_P53_PATHWAY, HALLMARK_INFLAMMATORY_RESPONSE, HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION). Supports human and mouse collections (species parameter).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MSigDB_get_hallmark_geneset',
     'arguments': {'pathway_name': '...'}},
)
```

## `MyGene_batch_query` {#MyGene_batch_query}

Queries MyGene.info for annotations of multiple gene IDs, supporting ID mapping.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_ids`
- **Keywords:** gene id mapping, batch, entrez, ensembl
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_ids` | array of string | yes | List of gene IDs to query. Can mix different ID types (Entrez, Ensembl, symbols). |
| `species` | string | no | Species filter for the query. Default `human` |
| `fields` | string | no | Comma-separated list of fields to return for each gene. Default `symbol,name,entrezgene,ensembl.gene` |

??? note "Upstream description"

    Query multiple genes at once using a list of gene identifiers. Efficient for retrieving annotations for a gene set or performing ID mapping between different gene ID types.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MyGene_batch_query', 'arguments': {'gene_ids': []}},
)
```

## `MyGene_get_gene_annotation` {#MyGene_get_gene_annotation}

Retrieves MyGene.info annotation for a gene by Entrez or Ensembl ID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** gene annotation, gene ontology, pathways, entrez
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | Gene ID to query. Can be Entrez Gene ID (e.g., '1017' for CDK2) or Ensembl ID (e.g., 'ENSG00000123374'). |
| `fields` | string | no | Comma-separated list of fields to return. Available fields include: symbol, name, summary, go (gene ontology), pathway (KEGG/Reactome), interpro (protein domains), generif (gene references), pdb, uniprot. Default `symbol,name,entrezgene,ensembl,summary,go,pathway,interpro` |

??? note "Upstream description"

    Get detailed annotation for a specific gene by its ID (Entrez or Ensembl). Returns comprehensive information including gene ontology, pathways, protein domains, and literature references.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MyGene_get_gene_annotation',
     'arguments': {'gene_id': '...'}},
)
```

## `MyGene_query_genes` {#MyGene_query_genes}

Searches MyGene.info for genes by keyword, symbol, name or identifier.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** gene search, hgnc, entrez, ensembl, symbol lookup
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search query. Can be gene symbol (e.g., 'CDK2'), name ('cyclin dependent kinase'), Entrez ID ('1017'), or Ensembl ID ('ENSG00000123374'). Supports wildcards (*) and boolean operators (AND, OR, NOT). |
| `species` | string | no | Species filter. Use common name or NCBI taxonomy ID. Examples: 'human', 'mouse', '9606' (human), 'all'. Default `human` |
| `fields` | string | no | Comma-separated list of fields to return. Common fields: symbol, name, entrezgene, ensembl.gene, summary, go, pathway, interpro. Default `symbol,name,entrezgene,ensembl.gene,summary` |
| `size` | integer | no | Maximum number of results to return (1-100). Default `10` |

??? note "Upstream description"

    Search for genes by keyword, symbol, name, or other identifiers. Returns gene annotations from 30+ sources including Entrez Gene, Ensembl, UniProt, HGNC. Useful for finding gene IDs, symbols, descriptions, and basic annotations.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MyGene_query_genes', 'arguments': {'query': '...'}},
)
```

## `NCBIDatasets_get_gene` {#NCBIDatasets_get_gene}

Returns NCBI gene information by Gene ID including symbol, chromosomal location, genomic coordinates and cross-references.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** ncbi gene id, entrez, genomic coordinates, gene record
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | NCBI Gene ID (numeric). Examples: '7157' (TP53), '672' (BRCA1), '3630' (INS), '22059' (mouse Trp53). |

??? note "Upstream description"

    Get comprehensive gene information from the National Center for Biotechnology Information (NCBI) by Gene ID. Returns gene symbol, description, chromosomal location, genomic coordinates, cross-references to UniProt/Ensembl/OMIM, synonyms, and annotation details. Covers all organisms in the NCBI Gene database. No API key is required; optionally set NCBI_API_KEY for a higher request-rate limit. Example: gene_id 7157 returns TP53 (tumor protein p53) on chromosome 17.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NCBIDatasets_get_gene', 'arguments': {'gene_id': '...'}},
)
```

## `NCBIDatasets_get_gene_by_symbol` {#NCBIDatasets_get_gene_by_symbol}

Looks up NCBI Gene records by official gene symbol and species (name or taxonomy ID).

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `symbol`
- **Keywords:** gene symbol, taxon, ncbi gene, species
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `symbol` | string | yes | Gene symbol. Examples: 'TP53', 'BRCA1', 'INS', 'EGFR'. |
| `taxon` | string | no | Organism as common name or taxonomy ID. Examples: 'human', 'mouse', 'rat', '9606', '10090'. Default: 'human' (reported as metadata.query_taxon_defaulted). Synonyms `organism` and `species` are accepted. |
| `organism` | string | no | Synonym for `taxon`. Used when `taxon` is not supplied. |
| `species` | string | no | Synonym for `taxon`. Used when neither `taxon` nor `organism` is supplied. |

??? note "Upstream description"

    Look up gene information by gene symbol and organism. Searches NCBI Gene database using official gene symbols. Specify the species with `taxon` (or its synonyms `organism`/`species`) as a common name (human, mouse, rat) or taxonomy ID. IMPORTANT -- the species defaults to HUMAN when none is given, so a gene symbol shared across species returns the human gene unless you say otherwise; `metadata.query_taxon_defaulted` is true whenever that default was applied, and `metadata.query_taxon` echoes the species actually used. Every returned row also carries its own `tax_id`/`taxname`, which is the authoritative check. Returns gene ID, description, chromosomal location, and database cross-references.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NCBIDatasets_get_gene_by_symbol',
     'arguments': {'symbol': '...'}},
)
```

## `NCBIDatasets_get_genome_assembly` {#NCBIDatasets_get_genome_assembly}

Returns NCBI genome assembly metadata for a GCF_ or GCA_ accession including level, status and statistics.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `accession`
- **Keywords:** genome assembly, gcf, gca, refseq category, bioproject
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | Assembly accession, e.g. 'GCF_000005845.2' (E. coli K-12) or 'GCF_000001405.40' (human GRCh38.p14). |

??? note "Upstream description"

    Get genome assembly metadata for an assembly accession (RefSeq GCF_ or GenBank GCA_) from the NCBI Datasets v2 API. Returns organism/strain, assembly name/level/status, RefSeq category (e.g. reference genome), release date, submitter, BioProject, and assembly statistics (total length, chromosome count, contig/scaffold N50, GC%, annotation provider). Use to characterize a specific genome assembly across any organism (bacteria, fungi, plant, animal). No API key required.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NCBIDatasets_get_genome_assembly',
     'arguments': {'accession': '...'}},
)
```

## `NCBIDatasets_get_orthologs` {#NCBIDatasets_get_orthologs}

Returns NCBI orthologs of a gene across species with symbols and gene types.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** orthologs, comparative genomics, ncbi gene id
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | NCBI Gene ID (numeric). Examples: '7157' (TP53), '672' (BRCA1), '3630' (INS). |
| `page_size` | integer | no | Maximum number of orthologs to return (1-100). Default: 20. |

??? note "Upstream description"

    Get orthologous genes across species for a given NCBI Gene ID. Returns orthologs identified by NCBI's Ortholog pipeline, including gene symbols, species, and gene types for each ortholog. Useful for comparative genomics and cross-species gene function studies. Example: gene_id 7157 (TP53) returns mouse Trp53, rat Tp53, zebrafish tp53, etc.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NCBIDatasets_get_orthologs',
     'arguments': {'gene_id': '...'}},
)
```

## `NCBIDatasets_get_sequence_reports` {#NCBIDatasets_get_sequence_reports}

Returns per-sequence reports (chromosome, plasmid, scaffold) for a genome assembly with accessions, length and GC%.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `accession`
- **Keywords:** chromosome, scaffold, gc content, assembly report
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | Assembly accession, e.g. 'GCF_000005845.2'. |
| `page_size` | integer \| null | no | Maximum sequence records in this page (default 100, maximum 1000). |
| `page_token` | string \| null | no | Opaque next_page_token returned by a previous call. Omit for the first page. |

??? note "Upstream description"

    Get a page of per-sequence (chromosome/plasmid/scaffold) reports for a genome assembly via the NCBI Datasets v2 API. Returns each molecule's chromosome name, role, location type, RefSeq and GenBank accessions, length, and GC%, plus total_available and next_page_token metadata. Use next_page_token to retrieve additional pages for assemblies with many sequences. Use to map an assembly to its individual sequence accessions (e.g. resolve E. coli K-12 to NC_000913.3, or human chromosomes to their RefSeq NC_ accessions). No API key required.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NCBIDatasets_get_sequence_reports',
     'arguments': {'accession': '...'}},
)
```

## `NCBIDatasets_list_genomes_by_taxon` {#NCBIDatasets_list_genomes_by_taxon}

Lists NCBI genome assemblies for a taxon with accession, level, N50, GC% and RefSeq category.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `taxon`
- **Keywords:** genome assemblies, reference genome, n50, taxon
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `taxon` | string | yes | NCBI tax id or scientific name, e.g. '562', 'Escherichia coli', 'Saccharomyces cerevisiae'. |
| `limit` | integer \| null | no | Max assemblies to return (default 20, max 100). |
| `reference_only` | boolean \| null | no | If true, return only reference/representative genomes. |

??? note "Upstream description"

    List genome assemblies available for a taxon (NCBI tax id or scientific name) via the NCBI Datasets v2 API. Returns a summary list of assemblies (accession, assembly name/level, RefSeq category, release date, total length, N50, GC%) plus the total count available. Use to discover what genomes exist for an organism or clade, e.g. all Escherichia coli or Mycobacterium tuberculosis assemblies. Set reference_only to restrict to reference/representative genomes. No API key required.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NCBIDatasets_list_genomes_by_taxon',
     'arguments': {'taxon': '...'}},
)
```

## `NCBIGene_get_summary` {#NCBIGene_get_summary}

Returns an NCBI Gene summary by Entrez ID with function, chromosomal location, coordinates, exon count and MIM IDs.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `id`
- **Keywords:** entrez gene, esummary, gene summary, map location
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `id` | string | yes | NCBI Gene ID (Entrez Gene ID). Can be a single ID or comma-separated list. Examples: '7157' (TP53), '672' (BRCA1), '1956' (EGFR), '3845' (KRAS). Get IDs from NCBIGene_search. |

??? note "Upstream description"

    Get a detailed gene summary from NCBI Gene by Entrez Gene ID using the E-utilities ESummary API. Returns gene name, symbol, description, function summary, chromosomal location, map position, aliases, genomic coordinates, exon count, and MIM (OMIM) IDs. The most authoritative source for human gene information. Use NCBIGene_search first to find Gene IDs. Example: Gene ID 7157 returns TP53 on chromosome 17p13.1 with a detailed functional summary.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NCBIGene_get_summary', 'arguments': {'id': '...'}},
)
```

## `NCBIGene_search` {#NCBIGene_search}

Searches NCBI Gene by symbol, name or keyword with organism filter and returns Entrez Gene IDs.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `term`
- **Keywords:** esearch, entrez gene id, gene search
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `term` | string | yes | Search query. Use gene symbol with organism filter for best results. Format: 'GENE\[Symbol\] AND Organism\[Organism\]'. Examples: 'TP53\[Symbol\] AND Homo sapiens\[Organism\]', 'BRCA1\[Symbol\] AND Homo sapiens\[Organism\]', 'insulin AND Homo sapiens\[Organism\]', 'EGFR AND Mus musculus\[Organism\]'. |
| `retmax` | integer | no | Maximum number of results to return (default 10, max 100). |

??? note "Upstream description"

    Search the NCBI Gene database using the E-utilities API. Find gene IDs (Entrez Gene IDs) by gene symbol, name, or keyword combined with an organism filter. Returns a list of matching Gene IDs that can be used with NCBIGene_get_summary. Essential for resolving gene names to standardized NCBI Gene IDs. Example: searching 'TP53' in Homo sapiens returns Gene ID 7157; 'BRCA1' returns 672.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NCBIGene_search', 'arguments': {'term': '...'}},
)
```

## `OMA_get_genome_pair_orthologs` {#OMA_get_genome_pair_orthologs}

Returns all pairwise orthologs between two whole genomes from OMA.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `genome1`, `genome2`
- **Keywords:** orthologs, genome pair, oma, comparative genomics
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `genome1` | string | yes | First genome: UniProt species code (e.g. 'HUMAN', 'PANTR') or NCBI taxon ID (e.g. '9606'). |
| `genome2` | string | yes | Second genome: UniProt species code (e.g. 'MOUSE', 'PANTR') or NCBI taxon ID (e.g. '10090'). |
| `per_page` | integer \| null | no | Number of ortholog pairs to return per page (default: 20, max: 100). The full proteome-vs-proteome set is large; paginate with 'page'. |
| `page` | integer \| null | no | Page number for paginated results (default: 1). |

??? note "Upstream description"

    Retrieve the complete set of pairwise orthologs between TWO whole genomes (proteome-vs-proteome) in OMA. Unlike OMA_get_orthologs (which returns orthologs of a single query protein), this returns every orthologous protein PAIR between two species genomes -- the data used to build species-pair ortholog tables, dN/dS panels, and comparative-genomics / synteny analyses. Each pair reports entry_1 and entry_2 (omaid, canonical_id, species, HOG/group), the relationship type (1:1, 1:n, n:1, n:m), evolutionary distance, and alignment score. Provide species codes (e.g. genome1='HUMAN', genome2='MOUSE') or NCBI taxon IDs. Results are paginated; use per_page and page. Example: genome1='HUMAN', genome2='MOUSE' returns pairs such as HUMAN01436 (NK2R_HUMAN) &lt;-&gt; MOUSE00001 (NK2R_MOUSE), rel_type 1:1.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OMA_get_genome_pair_orthologs',
     'arguments': {'genome1': '...', 'genome2': '...'}},
)
```

## `OMA_get_group` {#OMA_get_group}

Returns an OMA Group of strict 1:1 orthologs with member proteins.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `group_id`
- **Keywords:** oma group, 1:1 orthologs, gene family
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `group_id` | string | yes | OMA Group number (numeric ID). Examples: '1388790' (p53 group), '839588' (EGFR group). Use OMA_get_protein to find the group ID for a specific protein. |

??? note "Upstream description"

    Get OMA Group details including member proteins. OMA Groups contain sets of strict 1:1 orthologs across all genomes, ensuring each species contributes at most one gene. This is the most stringent orthology definition in OMA. Example: group 1388790 returns the p53 OMA Group with fingerprint 'HKKGEPC', description 'cellular tumor antigen p53', and members from fish, mammals, and other vertebrates.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OMA_get_group', 'arguments': {'group_id': '...'}},
)
```

## `OMA_get_hog` {#OMA_get_hog}

Returns OMA Hierarchical Orthologous Group information covering gene family evolution and duplications.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `hog_id`
- **Keywords:** hog, gene family, duplication, orthology
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `hog_id` | string | yes | Hierarchical Orthologous Group ID. Examples: 'HOG:F0782425' (p53 family), 'HOG:F0798498' (insulin family). Use OMA_get_protein to find the current HOG ID for a specific protein (the 'oma_hog_id' field). Note: OMA reassigns HOG IDs between releases, so IDs from older releases may return 'no longer valid'. |

??? note "Upstream description"

    Get Hierarchical Orthologous Group (HOG) information from OMA. HOGs represent groups of genes that descended from a single ancestral gene within a given taxonomic range. HOGs capture the full evolutionary history of gene families including gene duplications. Example: HOG:F0782425 returns the p53 HOG at the Euteleostomi level with description 'cellular tumor antigen p53', completeness score 0.61, and child HOGs for fish, mammals, etc.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OMA_get_hog', 'arguments': {'hog_id': '...'}},
)
```

## `OMA_get_orthologs` {#OMA_get_orthologs}

Returns OMA pairwise orthologs for a protein with relationship type and evolutionary distance.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `protein_id`
- **Keywords:** orthologs, 1:many, oma, homologs
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `protein_id` | string | yes | UniProt accession (e.g., 'P04637' for human p53) or OMA ID. |
| `rel_type` | string \| null | no | Filter by orthology relationship type. Options: '1:1' (one-to-one), '1:n' (one-to-many), 'n:1' (many-to-one), 'n:m' (many-to-many). Leave empty for all types. |
| `per_page` | integer \| null | no | Number of orthologs to return (default: 20). OMA does not paginate this endpoint -- it returns the complete ortholog set in one response -- so this limit is applied by ToolUniverse after fetching, and there is no maximum: set it to 'total_count' to receive every ortholog. |

??? note "Upstream description"

    Get pairwise orthologs for a protein from the OMA (Orthologous MAtrix) Browser. Returns orthologous proteins across species with relationship type (1:1, 1:many, many:1, many:many), evolutionary distance, and alignment score. OMA orthologs are inferred using a rigorous algorithm that provides high-confidence predictions. Example: P04637 (human p53) returns orthologs in 200+ species including Mus musculus, Danio rerio, and Drosophila melanogaster. The response reports 'count' (orthologs in this response) alongside 'total_count' (all orthologs OMA holds for the protein); these differ whenever 'per_page' is smaller than the full set, which is usual -- P04637 has 157. When they differ 'truncated' is true and 'truncation_note' says how to widen the request. Treat 'count' as a page size and never as the number of species with an ortholog.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OMA_get_orthologs', 'arguments': {'protein_id': '...'}},
)
```

## `OMA_get_protein` {#OMA_get_protein}

Returns OMA protein entry details with species, genomic location and group membership.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `protein_id`
- **Keywords:** oma entry, genomic location, uniprot
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `protein_id` | string | yes | UniProt accession (e.g., 'P04637' for human p53) or OMA ID (e.g., 'HUMAN31534'). UniProt accessions are recommended for cross-database compatibility. |

??? note "Upstream description"

    Get protein information from the OMA (Orthologous MAtrix) Browser by UniProt accession or OMA ID. OMA is a comprehensive orthology database covering 2,600+ complete genomes. Returns protein details including species, genomic location, OMA group membership, and Hierarchical Orthologous Group (HOG) assignment. Example: P04637 returns human TP53 with OMA ID HUMAN31534, HOG:F0782425.2c.7a, and genomic locus on chromosome 17.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OMA_get_protein', 'arguments': {'protein_id': '...'}},
)
```

## `OMA_get_protein_go` {#OMA_get_protein_go}

Returns GO annotations for an OMA entry with aspect, evidence code and information content.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `protein_id`
- **Keywords:** go annotation, evidence code, oma
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `protein_id` | string | yes | OMA ID (e.g. 'HUMAN17018') or UniProt accession (e.g. 'P04637'). |
| `aspect` | string \| null | no | Optional GO aspect filter: 'biological_process', 'molecular_function', or 'cellular_component'. Leave empty for all aspects. |

??? note "Upstream description"

    Get Gene Ontology (GO) functional annotations for a specific OMA protein. OMA attaches per-entry GO annotations with the aspect (biological_process / molecular_function / cellular_component), information content (ic), evidence code (e.g. IDA, NAS, IEA), and a literature reference (PMID). Other OMA tools return orthology / HOG / group membership but not these functional annotations. Accepts an OMA ID (e.g. 'HUMAN17018') or UniProt accession (e.g. 'P04637'); use OMA_resolve_xref or OMA_get_protein first if you only have a gene symbol. Example: protein_id='HUMAN17018' returns GO:0060261 'positive regulation of transcription initiation by RNA polymerase II', aspect biological_process, ic 8.27, evidence IDA, reference PMID:12218053.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OMA_get_protein_go', 'arguments': {'protein_id': '...'}},
)
```

## `OMA_resolve_xref` {#OMA_resolve_xref}

Resolves a gene name, symbol or cross-reference ID to OMA protein entries.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `search`
- **Keywords:** gene symbol, xref resolution, oma
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `search` | string | yes | Gene symbol (e.g. 'BRCA2'), UniProt entry name (e.g. 'MED4_HUMAN'), or any cross-reference identifier to resolve to OMA entries. |
| `limit` | integer \| null | no | Maximum number of cross-reference matches to return (default: 25, max: 100). |

??? note "Upstream description"

    Resolve a gene name, gene symbol, UniProt entry name, or any cross-reference identifier to OMA protein entries. This is the entry point into the OMA toolchain when you only have a common identifier (e.g. 'BRCA2' or 'MED4_HUMAN') rather than a UniProt accession or OMA ID. OMA_get_protein only accepts a UniProt accession (P04637) or OMA ID (HUMAN31534); use this tool first to obtain the entry_nr / oma_id / genome for the gene of interest, then feed the oma_id into OMA_get_protein, OMA_get_orthologs, OMA_get_protein_go, etc. Example: search='MED4_HUMAN' returns entry_nr 24423700, oma_id HUMAN17018, species Homo sapiens (taxon 9606); search='BRCA2' returns matching OMA entries across species.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OMA_resolve_xref', 'arguments': {'search': '...'}},
)
```

## `OmniPath_get_annotation_resource_geneset` {#OmniPath_get_annotation_resource_geneset}

Returns the full gene/protein membership of an OmniPath annotation resource.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `resource`
- **Keywords:** omnipath, gene set, annotation resource
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `resource` | string | yes | Name of the OmniPath annotation resource whose full gene/protein set to retrieve. Examples: 'CancerGeneCensus' (cancer driver genes), 'Surfaceome' (cell-surface proteins), 'TFcensus' (transcription factors), 'kinase.com' (kinases), 'Phosphatome' (phosphatases). Many others available (80+ resources). Verify exact names via https://omnipathdb.org/queries/annotations. |
| `entity_types` | string \| null | no | Optional filter on entity type, comma-separated. Options: 'protein', 'complex', 'mirna'. Default: all entity types in the resource. |

??? note "Upstream description"

    Retrieve the ENTIRE annotated gene/protein set (resource-wide membership) for an OmniPath annotation resource, with NO protein filter. This is the reverse of per-protein annotation lookup: instead of 'what does OmniPath say about gene X', it answers 'list every gene/protein in resource Y'. Use it to fetch curated functional gene sets such as all cell-surface proteins (Surfaceome), all kinases (kinase.com), all transcription factors (TFcensus), all phosphatases (Phosphatome), or all cancer driver genes (CancerGeneCensus). OmniPath integrates 80+ such annotation resources. Each returned member lists its label/value annotations from the resource (grouped from OmniPath's long-format rows). Example: resource='CancerGeneCensus' returns 752 unique driver genes from 9755 records, including KRAS (UniProt P01116, tier=1, hallmark=True); resource='Surfaceome' returns the full surfaceome including SLC12A8 (A0AV02, mainclass=Transporters); resource='kinase.com' returns PRKCB (P05771, group=AGC); resource='TFcensus' returns HIVEP2 (P31629, tfcensus_class=a).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OmniPath_get_annotation_resource_geneset',
     'arguments': {'resource': '...'}},
)
```

## `OpenGenes_get_gene` {#OpenGenes_get_gene}

Returns the Open Genes aging/longevity profile of a human gene with mechanisms, clusters and disease categories.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `symbol`
- **Keywords:** aging gene, longevity, open genes
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `symbol` | string | yes | Human gene symbol, e.g. 'GHR', 'FOXO3', 'TP53', 'SIRT1'. |

??? note "Upstream description"

    Get the aging/longevity profile of a human gene by symbol from Open Genes, a manually-curated aging-gene database. Returns the gene's aging mechanisms (e.g. INS/IGF-1 pathway dysregulation), functional clusters, associated disease categories, curator confidence level, expression change with age, and an evidence-count summary across study types (lifespan-increase studies, longevity associations, age-related changes, progeria associations, interventions). Use to assess whether and how a gene is implicated in aging. No API key required.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenGenes_get_gene', 'arguments': {'symbol': '...'}},
)
```

## `OpenGenes_search_genes` {#OpenGenes_search_genes}

Returns a page of Open Genes aging-gene summaries with mechanisms and confidence.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** aging genes, longevity, gene catalog
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `limit` | integer \| null | no | Genes per page (default 20, max 100). |
| `page` | integer \| null | no | Page number (default 1). |

??? note "Upstream description"

    Browse the Open Genes catalog of human aging/longevity genes (paginated). Returns a page of aging-gene summaries (symbol, name, aging mechanisms, functional clusters, confidence level, expression change) plus the total count of curated aging genes. Use to enumerate aging-associated genes or sample the catalog. No API key required.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenGenes_search_genes', 'arguments': {}},
)
```

## `OpenTargets_get_target_cancer_hallmarks` {#OpenTargets_get_target_cancer_hallmarks}

Retrieves COSMIC Cancer Gene Census hallmark annotations (promotes/suppresses) for a gene by Ensembl ID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensemblId`
- **Keywords:** cosmic, cancer gene census, oncogene, tumor suppressor
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensemblId` | string | yes | Ensembl gene ID of the target (e.g., 'ENSG00000146648' for EGFR). |

??? note "Upstream description"

    Retrieve cancer hallmarks annotation for a target by Ensembl gene ID (from the COSMIC Cancer Gene Census / Hanahan-Weinberg hallmarks). Returns cancerHallmarks (each with hallmark label, impact = 'promotes' or 'suppresses', and a mechanistic description with PMID) plus role-in-cancer attributes (name + description). Use to summarise a target's oncogenic mechanisms.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_target_cancer_hallmarks',
     'arguments': {'ensemblId': '...'}},
)
```

## `OpenTargets_get_target_classes_by_ensemblID` {#OpenTargets_get_target_classes_by_ensemblID}

Returns target class annotations (e.g. enzyme, receptor families) for an Ensembl gene ID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensemblId`
- **Keywords:** protein class, chembl target class, gene family
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensemblId` | string | yes | The Ensembl ID of the target. |

??? note "Upstream description"

    Retrieve the target classes associated with a specific target ensemblID.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_target_classes_by_ensemblID',
     'arguments': {'ensemblId': '...'}},
)
```

## `OpenTargets_get_target_depmap_essentiality` {#OpenTargets_get_target_depmap_essentiality}

Retrieves DepMap CRISPR essentiality and Chronos gene effect scores per tissue/cell line for an Ensembl gene ID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensemblId`
- **Keywords:** depmap, crispr screen, chronos, essential gene
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensemblId` | string | yes | Ensembl gene ID of the target (e.g., 'ENSG00000146648' for EGFR). |

??? note "Upstream description"

    Retrieve DepMap CRISPR gene essentiality for a target by Ensembl gene ID. Returns the overall isEssential flag plus per-tissue essentiality: for each tissue, the cell-line screens with depmapId, cell line name, geneEffect (Chronos score; more negative = more essential, roughly &lt; -0.5 = dependency), and the cell line's expression of the gene. Use to assess cancer cell-line dependency on a target.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_target_depmap_essentiality',
     'arguments': {'ensemblId': '...'}},
)
```

## `OpenTargets_get_target_gene_ontology_by_ensemblID` {#OpenTargets_get_target_gene_ontology_by_ensemblID}

Retrieves Gene Ontology annotations for a target by Ensembl gene ID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensemblId`
- **Keywords:** go annotation, biological process, molecular function
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensemblId` | string | yes | The Ensembl ID of the target for which to retrieve Gene Ontology annotations. |

??? note "Upstream description"

    Retrieve Gene Ontology annotations for a specific target by Ensembl ID.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_target_gene_ontology_by_ensemblID',
     'arguments': {'ensemblId': '...'}},
)
```

## `OpenTargets_get_target_genomic_location_by_ensemblID` {#OpenTargets_get_target_genomic_location_by_ensemblID}

Returns chromosome, start, end and strand for a target Ensembl gene ID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensemblId`
- **Keywords:** chromosome, coordinates, strand, locus
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensemblId` | string | yes | The Ensembl ID of the target for which to retrieve genomic location information. |

??? note "Upstream description"

    Retrieve genomic location data for a specific target, including chromosome, start, end, and strand.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_target_genomic_location_by_ensemblID',
     'arguments': {'ensemblId': '...'}},
)
```

## `OpenTargets_get_target_homologues_by_ensemblID` {#OpenTargets_get_target_homologues_by_ensemblID}

Retrieves homologues (orthologs/paralogs) of a target by Ensembl gene ID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensemblId`
- **Keywords:** ortholog, paralog, homology
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensemblId` | string | yes | The Ensembl ID of the target for which to retrieve homologues. |

??? note "Upstream description"

    Fetch homologues for a specific target by Ensembl ID.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_target_homologues_by_ensemblID',
     'arguments': {'ensemblId': '...'}},
)
```

## `OpenTargets_get_target_id_description_by_name` {#OpenTargets_get_target_id_description_by_name}

Looks up the Ensembl gene ID and description for a target name or symbol.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `targetName`
- **Keywords:** gene symbol, ensembl id, id mapping
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `targetName` | string | yes | The name of the target for which the ID is required. |

??? note "Upstream description"

    Get the ensemblId and description based on the target name.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_target_id_description_by_name',
     'arguments': {'targetName': '...'}},
)
```

## `OpenTargets_get_target_info_by_ensemblID` {#OpenTargets_get_target_info_by_ensemblID}

Returns core target information for an Ensembl gene ID: symbol, name, biotype, UniProt IDs, location and function.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensemblId`
- **Keywords:** gene symbol, biotype, uniprot, gene function
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensemblId` | string | yes | Ensembl gene ID of the target (e.g., 'ENSG00000141510' for TP53). Resolve a gene symbol to an Ensembl ID with OpenTargets_get_target_id_description_by_name. |

??? note "Upstream description"

    Get core information about a drug target from Open Targets by its Ensembl gene ID: approved symbol and name, biotype, UniProt/protein cross-reference IDs, genomic location, and curated function descriptions. Use this as the primary lookup for a target before fetching its disease associations, tractability, or safety. Example: ENSG00000141510 returns TP53 (tumor protein p53), protein_coding, chr17.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_target_info_by_ensemblID',
     'arguments': {'ensemblId': '...'}},
)
```

## `OpenTargets_get_target_interactions_by_ensemblID` {#OpenTargets_get_target_interactions_by_ensemblID}

Returns interaction partners and evidence for a target Ensembl gene ID, paginated.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensemblId`
- **Keywords:** ppi, intact, string, interaction network
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensemblId` | string | yes | The Ensembl ID of the target. |
| `page` | object | no | Pagination parameters. |

??? note "Upstream description"

    Retrieve interaction data for a specific target ensemblID, including interaction partners and evidence. Results are PAGINATED and Open Targets returns only 25 rows when no `page` is given: `interactions.count` is the total number of interactions (a well-studied target has thousands), `interactions.returned` is how many are in `rows`, and `interactions.truncated` is true when interactions remain unfetched, in which case `metadata.truncation_note` names the exact `page` argument to pass. Pass page: \{"index": 0, "size": N\} (Open Targets rejects a size above 3000) and step `index` to page through.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_target_interactions_by_ensemblID',
     'arguments': {'ensemblId': '...'}},
)
```

## `OpenTargets_get_target_pathways_by_ensemblID` {#OpenTargets_get_target_pathways_by_ensemblID}

Lists Reactome pathways a gene participates in by Ensembl ID, with top-level terms.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensemblId`
- **Keywords:** reactome, pathway membership, gene to pathway
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensemblId` | string | yes | Ensembl gene ID of the target (e.g., 'ENSG00000146648' for EGFR). |

??? note "Upstream description"

    Retrieve Reactome pathway membership for a target by Ensembl gene ID. Returns each pathway the gene participates in with its Reactome pathwayId (e.g., R-HSA-...), pathway name, and top-level Reactome term (e.g., 'Signal Transduction'). Use to understand the biological processes a target is involved in.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_target_pathways_by_ensemblID',
     'arguments': {'ensemblId': '...'}},
)
```

## `OpenTargets_get_target_synonyms_by_ensemblID` {#OpenTargets_get_target_synonyms_by_ensemblID}

Retrieves synonyms and alternative symbols for a target Ensembl gene ID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensemblId`
- **Keywords:** gene synonyms, alias, symbol
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensemblId` | string | yes | The Ensembl ID of the target. |

??? note "Upstream description"

    Retrieve synonyms for specified target, including alternative names and symbols, using given ensemblID.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_target_synonyms_by_ensemblID',
     'arguments': {'ensemblId': '...'}},
)
```

## `OrthoDB_get_group_details` {#OrthoDB_get_group_details}

Returns KEGG, GO and InterPro annotations for an OrthoDB orthologous group.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `group_id`
- **Keywords:** orthologous group, go, kegg, interpro
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `group_id` | string | yes | OrthoDB group ID. Format: 'numberATtaxid'. Examples: '727649at7742' (Vertebrata level), '59506at9443' (Primates level). Get IDs from OrthoDB_search_groups. |

??? note "Upstream description"

    Get detailed information about an OrthoDB orthologous group including functional annotations. Returns KEGG pathway associations, GO terms, and InterPro domain annotations for the group. Useful for understanding the evolutionary conservation and function of gene families. Example: group '727649at7742' shows olfactory receptor group at Vertebrata level with associated KEGG pathways for signal transduction.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OrthoDB_get_group_details',
     'arguments': {'group_id': '...'}},
)
```

## `OrthoDB_get_orthologs` {#OrthoDB_get_orthologs}

Lists member genes (IDs, organisms, descriptions) of an OrthoDB orthologous group, optionally by species.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `group_id`
- **Keywords:** ortholog, gene id, species, conservation
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `group_id` | string | yes | OrthoDB group ID. Format: 'numberATtaxid'. Get IDs from OrthoDB_search_groups. |
| `species` | string \| null | no | Comma-separated NCBI taxonomy IDs to filter. Examples: '9606' (human only), '9606,10090' (human+mouse). Null for all species. |

??? note "Upstream description"

    Get orthologous genes from an OrthoDB group, optionally filtered by species. Returns gene IDs, organism names, and descriptions for all members of an orthologous group. Useful for finding orthologs across species, evolutionary conservation analysis, and functional annotation transfer. Example: group '727649at7742' with species=9606 returns human genes in that orthologous group.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OrthoDB_get_orthologs', 'arguments': {'group_id': '...'}},
)
```

## `OrthoDB_search_groups` {#OrthoDB_search_groups}

Searches OrthoDB v12 orthologous groups by gene/protein name and taxonomic level.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** ortholog, orthologous group, taxonomic level
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Gene or protein name to search. Examples: 'BRCA1', 'TP53', 'insulin', 'EGFR'. |
| `species` | integer \| null | no | NCBI taxonomy ID to filter by species. Examples: 9606 (human), 10090 (mouse), 7227 (fly). Null for all species. |
| `level` | integer \| null | no | Taxonomic level for group definition. Examples: 7742 (Vertebrata), 33208 (Metazoa), 2759 (Eukaryota). Null for all levels. |
| `limit` | integer | no | Maximum number of groups to return (1-50). Default: 10. |

??? note "Upstream description"

    Search for orthologous groups by gene/protein name in OrthoDB v12. OrthoDB catalogs orthologous groups at different taxonomic levels (vertebrates, metazoa, eukaryotes, etc.), enabling evolutionary analysis of gene conservation. Returns group IDs that can be queried for details. Complements OMA and Ensembl Compara with group-level orthology at deep evolutionary scales. Example: searching 'BRCA1' with species=9606 returns orthologous groups containing BRCA1 at various taxonomic levels.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OrthoDB_search_groups', 'arguments': {'query': '...'}},
)
```

## `PANTHER_enrichment` {#PANTHER_enrichment}

Runs PANTHER overrepresentation analysis on a gene list for GO terms, pathways or protein classes.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_list`
- **Keywords:** enrichment, overrepresentation, go, fisher test
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_list` | string | yes | Comma-separated list of gene symbols, UniProt IDs, or Ensembl IDs. Example: 'TP53,BRCA1,EGFR,KRAS,PTEN,APC,RB1'. |
| `organism` | integer \| null | no | NCBI taxonomy ID. Default: 9606 (human). Others: 10090 (mouse), 7227 (fly), 559292 (yeast). |
| `annotation_dataset` | string \| null | no | Annotation dataset for enrichment. Options: 'GO:0008150' (biological_process, default), 'GO:0003674' (molecular_function), 'GO:0005575' (cellular_component), 'ANNOT_TYPE_ID_PANTHER_PATHWAY' (PANTHER pathway), 'ANNOT_TYPE_ID_REACTOME_PATHWAY' (Reactome pathway), 'ANNOT_TYPE_ID_PANTHER_PC' (protein class). |

??? note "Upstream description"

    Perform gene set enrichment (overrepresentation) analysis using PANTHER. Given a gene list, identifies statistically enriched GO terms, PANTHER pathways, Reactome pathways, or protein classes using Fisher's exact test with FDR correction. PANTHER is one of the most widely used tools for functional enrichment in bioinformatics. Example: gene list 'TP53,BRCA1,EGFR,KRAS,PTEN' enriched for GO biological_process returns 'regulation of mitotic cell cycle phase transition' (FDR=0.003).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PANTHER_enrichment', 'arguments': {'gene_list': '...'}},
)
```

## `PANTHER_gene_info` {#PANTHER_gene_info}

Returns PANTHER family/subfamily, GO annotations and protein class for a gene.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** panther family, go annotation, gene function
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | Gene identifier - UniProt accession, Ensembl ID, or gene symbol. Examples: 'P04637' (TP53), 'P38398' (BRCA1), 'ENSG00000141510'. |
| `organism` | integer \| null | no | NCBI taxonomy ID. Default: 9606 (human). Others: 10090 (mouse), 10116 (rat), 7227 (fly), 6239 (worm), 559292 (yeast). |

??? note "Upstream description"

    Get gene classification and functional annotation from PANTHER (Protein ANalysis THrough Evolutionary Relationships). Returns PANTHER family/subfamily, Gene Ontology annotations (molecular function, biological process, cellular component), PANTHER protein class, and pathway membership. PANTHER classifies genes by evolutionary relationships across 144 organisms. Example: P04637 returns TP53 classified in family PTHR11447 (tumor protein p53) with GO annotations for DNA-binding transcription factor activity.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PANTHER_gene_info', 'arguments': {'gene_id': '...'}},
)
```

## `PANTHER_ortholog` {#PANTHER_ortholog}

Finds orthologs of a gene across species using PANTHER, with ortholog type.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** ortholog, ldo, cross-species
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | Gene identifier - UniProt accession, Ensembl ID, or gene symbol. Examples: 'P04637' (TP53), 'P38398' (BRCA1). |
| `organism` | integer \| null | no | Source organism NCBI taxonomy ID. Default: 9606 (human). |
| `target_organism` | integer \| null | no | Target organism NCBI taxonomy ID. Default: 10090 (mouse). Others: 10116 (rat), 7227 (fly), 6239 (worm). |
| `ortholog_type` | string \| null | no | Ortholog type filter. Options: 'LDO' (least diverged ortholog, default), 'O' (all orthologs), 'P' (paralogs), 'X' (xenologs), 'LDX' (least diverged xenolog), 'all' (all types). |

??? note "Upstream description"

    Find orthologs of a gene across species using PANTHER evolutionary classification. Returns EVERY matching ortholog in 'mappings' (count in 'total_mappings'), each with its persistent ID, gene symbol, and ortholog type (LDO=Least Diverged Ortholog, O=ortholog, P=paralog, X=xenolog, LDX=least diverged xenolog); 'mapping' repeats the first entry for backward compatibility and is not the whole answer. A gene often has several orthologs in one target species -- human ABCB1 maps to both mouse Abcb1a (LDO) and Abcb1b (O) -- so read 'mappings', especially with ortholog_type='O'. PANTHER classifies orthologs by evolutionary divergence across 144 fully sequenced genomes. Example: P04637 (human TP53) maps to Tp53 (mouse P02340) as LDO.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PANTHER_ortholog', 'arguments': {'gene_id': '...'}},
)
```

## `PathwayCommons_get_neighborhood` {#PathwayCommons_get_neighborhood}

Returns known interactions (PPI, regulatory, catalytic) for a gene from Pathway Commons.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene`
- **Keywords:** interaction network, biopax, partner genes
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string | yes | Gene symbol (e.g., BRCA1, TP53, EGFR) |
| `limit` | integer \| null | no | Graph search depth (1 = direct neighbors, 2 = neighbors of neighbors). Default 1. Higher values return much more data. Default `1` |
| `datasource` | string \| null | no | Filter by data source: reactome, kegg, biogrid, intact, hprd, pid, corum, ctd, drugbank |
| `max_results` | integer \| null | no | Maximum number of interactions to return (default: 100) Default `100` |

??? note "Upstream description"

    Get the interaction neighborhood for a gene from Pathway Commons. Returns all known interactions (protein-protein, regulatory, catalytic, complex membership) for the query gene from 22 integrated databases. Each interaction includes the partner gene, interaction type (interacts-with, controls-state-change-of, controls-expression-of, in-complex-with, etc.), data source, and PubMed references. Useful for building gene interaction networks and finding functional partners.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PathwayCommons_get_neighborhood',
     'arguments': {'gene': '...'}},
)
```

## `PathwayCommons_get_pathway` {#PathwayCommons_get_pathway}

Returns a Pathway Commons pathway by URI with its name, participants and interaction processes.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `uri`
- **Keywords:** biopax, pathway participants, reactome, kegg
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `uri` | string | yes | Pathway Commons URI from search results, e.g., http://bioregistry.io/reactome:R-HSA-5693532 (DNA Double-Strand Break Repair) |

??? note "Upstream description"

    Get detailed pathway information from Pathway Commons by URI. Returns the pathway name, participants (proteins, small molecules, complexes), and processes (interactions, catalysis, controls). Use URIs from PathwayCommons_search results. Covers Reactome, KEGG, WikiPathways, PID, and other sources. Response is parsed from BioPAX/JSON-LD format.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PathwayCommons_get_pathway', 'arguments': {'uri': '...'}},
)
```

## `PathwayCommons_paths_between` {#PathwayCommons_paths_between}

Finds mechanistic interaction paths between two or more gene symbols across Pathway Commons sources.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `genes`
- **Keywords:** pathsbetween, gene network, biopax, interaction path
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `genes` | array of string | yes | Two or more gene symbols / identifiers to connect, e.g. \["BRCA1", "BRCA2", "TP53"\]. A comma-separated string is also accepted. |
| `limit` | integer \| null | no | Path length limit (graph search distance). Default 1 (direct mechanistic links). Higher values include longer indirect paths and much more data. Default `1` |
| `datasource` | string \| null | no | Filter by data source: reactome, kegg, biogrid, intact, hprd, pid, corum, ctd, drugbank |
| `max_results` | integer \| null | no | Maximum number of edges to return (default: 200) Default `200` |

??? note "Upstream description"

    Find mechanistic interaction paths connecting a set of genes/proteins in Pathway Commons (PC2 kind=PATHSBETWEEN across 22 integrated databases: Reactome, KEGG, WikiPathways, PID, BioGRID, IntAct, etc.). Given two or more gene symbols, returns the directed/undirected interaction edges that link them, parsed from SIF triples into \{source, relation, target\}. Relations include in-complex-with, controls-state-change-of, controls-phosphorylation-of, controls-expression-of, chemical-affects, used-to-produce. Use to reveal how a gene set is mechanistically connected (e.g. a DNA-repair module BRCA1/BRCA2/TP53). Requires at least 2 distinct genes. Note: an uncached query can be slow; if it times out, retry (results are cached server-side).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PathwayCommons_paths_between',
     'arguments': {'genes': []}},
)
```

## `PDC_get_gene_protein` {#PDC_get_gene_protein}

Returns gene description, NCBI gene ID, protein accessions and PDC study coverage for a gene symbol.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_name`
- **Keywords:** gene to protein, ncbi gene id, cptac
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `get_gene_protein` Default `get_gene_protein` |
| `gene_name` | string | yes | Gene symbol to look up (e.g., 'TP53', 'EGFR', 'BRCA1', 'MYC') |

??? note "Upstream description"

    Get protein information and proteomics study coverage for a gene from the NCI Proteomics Data Commons (PDC). Given a gene symbol (e.g., TP53, EGFR, MYC), returns the gene description, NCBI gene ID, HGNC authority, associated protein accessions (UniProt, RefSeq, Ensembl), and spectral count data across all PDC studies where the protein was detected. Spectral counts indicate protein abundance evidence in each study. Useful for understanding which cancer proteomics datasets have measured a specific protein.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PDC_get_gene_protein', 'arguments': {'gene_name': '...'}},
)
```

## `PharmGKB_get_gene_details` {#PharmGKB_get_gene_details}

Returns chromosomal location and cross-references for a PharmGKB gene accession.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** pharmgkb, pharmacogene, chromosome
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | PharmGKB Gene Accession ID (e.g., 'PA128'). |

??? note "Upstream description"

    Get detailed information for a gene using its PharmGKB Accession ID. Returns chromosomal location, cross-references, and other metadata.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PharmGKB_get_gene_details',
     'arguments': {'gene_id': '...'}},
)
```

## `PharmGKB_search_genes` {#PharmGKB_search_genes}

Searches PharmGKB genes by name, symbol or ID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** pharmgkb, gene search, pharmacogene
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Gene name, symbol, or PharmGKB Accession ID (e.g., 'CYP2D6', 'PA128'). |

??? note "Upstream description"

    Search for genes in PharmGKB by name, symbol, or ID. Returns PharmGKB Gene ID and basic gene metadata.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PharmGKB_search_genes', 'arguments': {'query': '...'}},
)
```

## `Pharos_get_target` {#Pharos_get_target}

Returns Pharos target details (TDL, family, novelty, diseases) by gene symbol or UniProt ID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** tdl, tcrd, druggability, novelty
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string | no | Gene symbol (e.g., 'EGFR', 'TP53', 'BRCA1'). Use either gene or uniprot. |
| `uniprot` | string | no | UniProt accession (e.g., 'P00533'). Use either gene or uniprot. |

??? note "Upstream description"

    Get detailed drug target information from Pharos/TCRD by gene symbol or UniProt ID. Returns Target Development Level (Tdark/Tbio/Tchem/Tclin -- Pharos' druggability classification), protein family, novelty score, publication count, associated diseases (full per-disease association counts plus the top 10 diseases with MONDO IDs), and ligand/drug counts plus the top 10 ligands. Essential for identifying understudied proteins. For ligand bioactivity values use Pharos_get_target_ligands.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Pharos_get_target', 'arguments': {}},
)
```

## `Planteome_search_annotations` {#Planteome_search_annotations}

Returns GO, PO and TO annotations for a plant gene from Planteome.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** plant gene, go annotation, trait ontology
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Gene id or accession, e.g. 'AT4G32150'. |
| `limit` | integer \| null | no | Max annotations to return, 1-100. Default 20. |

??? note "Upstream description"

    Search Planteome's gene-to-ontology-term annotation index by gene id or accession. Resolves a plant gene to every ontology term (GO/PO/TO/crop ontologies) it has been annotated with, across evidence types and source databases. Example: query='AT4G32150' (Arabidopsis VAMP711) returns its GO annotations including 'vesicle-mediated transport'.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Planteome_search_annotations',
     'arguments': {'query': '...'}},
)
```

## `PlantReactome_get_pathway_participants` {#PlantReactome_get_pathway_participants}

Lists gene and protein participants with UniProt accessions and gene names for a Plant Reactome pathway.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `pathway_id`
- **Keywords:** pathway genes, plant, uniprot
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `pathway_id` | string | yes | Plant Reactome stable identifier of a pathway or reaction. Format: R-\{species_code\}-\{number\}. Example: 'R-OSA-1119263' (Arginine biosynthesis in rice). |

??? note "Upstream description"

    List the gene/protein/molecule participants of a Plant Reactome pathway or reaction (the macromolecules that are IN the pathway), with UniProt accessions and gene names. PlantReactome_get_pathway returns only pathway metadata and 'hasEvent' sub-event IDs, not the participant gene products; use this tool to get the actual genes/proteins. Each participant carries its refEntities, where 'identifier' is the UniProt accession and 'gene_names' lists the gene symbols. Example: R-OSA-1119263 (Arginine biosynthesis, rice) returns 36 participants; the first refEntity is UniProt Q10MK5 with gene names including LOC_OS03G19280. No API key required.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PlantReactome_get_pathway_participants',
     'arguments': {'pathway_id': '...'}},
)
```

## `PomBase_get_gene` {#PomBase_get_gene}

Returns PomBase gene details: name, product, domains, deletion viability and UniProt cross-reference.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** s. pombe, fission yeast, systematic id
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | PomBase systematic gene identifier. Examples: 'SPBC11B10.09' (cdc2), 'SPAC2F7.03c' (pom1), 'SPAC1002.17c' (ase1). |

??? note "Upstream description"

    Get detailed gene information from PomBase for fission yeast (S. pombe) by systematic ID. Returns gene name, protein product, InterPro domain annotations, deletion viability, UniProt cross-reference, and transmembrane domain coordinates. S. pombe is a key model organism for cell cycle, chromosome biology, and RNA processing research. Example: SPBC11B10.09 returns cdc2 (cyclin-dependent kinase), essential for cell division.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PomBase_get_gene', 'arguments': {'gene_id': '...'}},
)
```

## `PomBase_get_go_annotations` {#PomBase_get_go_annotations}

Returns GO annotations by aspect for an S. pombe gene from PomBase.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** go annotation, s. pombe, biological process
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | PomBase systematic ID. Examples: 'SPBC11B10.09' (cdc2). |

??? note "Upstream description"

    Get Gene Ontology annotations (biological_process / molecular_function / cellular_component) for a fission-yeast (S. pombe) gene from PomBase, with per-aspect counts and resolved term names. Example: SPBC11B10.09 (cdc2) returns BP terms like GO:0000082 'G1/S transition of mitotic cell cycle', plus MF and CC terms.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PomBase_get_go_annotations',
     'arguments': {'gene_id': '...'}},
)
```

## `PomBase_get_interactions` {#PomBase_get_interactions}

Returns physical and genetic interactions with evidence for an S. pombe gene.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** genetic interaction, physical interaction, s. pombe
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | PomBase systematic ID. Examples: 'SPBC11B10.09' (cdc2). |

??? note "Upstream description"

    Get physical and genetic interactions for a fission-yeast (S. pombe) gene from PomBase, each with evidence method, reference (PMID), throughput, and source database. Example: SPBC11B10.09 (cdc2) has ~100 physical interactions (e.g. SPAC1006.03c/red1 via Affinity Capture-MS, PMID:24713849, BIOGRID) and ~118 genetic interactions.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PomBase_get_interactions',
     'arguments': {'gene_id': '...'}},
)
```

## `PomBase_get_orthologs` {#PomBase_get_orthologs}

Returns human and budding yeast orthologs of an S. pombe gene.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** ortholog, hgnc, s. cerevisiae
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | PomBase systematic ID. Examples: 'SPBC11B10.09' (cdc2), 'SPAC1006.03c' (red1). |

??? note "Upstream description"

    Get human and S. cerevisiae orthologs of a fission-yeast (S. pombe) gene from PomBase. Returns ortholog IDs with taxon (9606=human via HGNC, 4932=budding yeast via systematic name). Example: SPBC11B10.09 (cdc2) returns human CDKs HGNC:1722/1771/1772 (taxon 9606) and S. cerevisiae CDC28 / YBR160W (taxon 4932).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PomBase_get_orthologs', 'arguments': {'gene_id': '...'}},
)
```

## `PomBase_search_genes` {#PomBase_search_genes}

Searches PomBase S. pombe genes by name, systematic ID or product.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** s. pombe, gene search, systematic id
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search query - gene name, systematic ID prefix, or product keyword. Examples: 'cdc2', 'kinase', 'SPAC', 'pom1'. |
| `limit` | integer \| null | no | Maximum results to return (1-50, default 10). |

??? note "Upstream description"

    Search PomBase for fission yeast (S. pombe) genes by name, systematic ID, or protein product keyword. Searches across 12,600+ gene entries. Returns matching genes with systematic IDs, gene names, and protein products. Example: searching 'kinase' returns CDK and other kinase genes; searching 'cdc' returns cell division cycle genes.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PomBase_search_genes', 'arguments': {'query': '...'}},
)
```

## `proteins_api_get_genome_mappings` {#proteins_api_get_genome_mappings}

Returns genome mappings (Ensembl, RefSeq, EMBL cross-references and coordinates) for UniProt accessions.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `accession`
- **Keywords:** genome coordinates, ensembl, refseq, protein to gene
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string \| array | yes | UniProt protein accession(s). Can be a single accession (e.g., 'P05067'), comma-separated string (e.g., 'P05067,P04637'), or list (e.g., \['P05067', 'P04637'\]). Use UniProt search or ebi_search_domain with domain='uniprot' to find valid accessions. |
| `format` | string | no | Response format. JSON is recommended for most use cases. One of `json`, `xml` Default `json` |

??? note "Upstream description"

    Get reference genome mappings for a protein using the Proteins API. Supports batch operations: pass a list of accessions or comma-separated string to retrieve multiple proteins in parallel. This extracts genome-related cross-references (Ensembl, RefSeq, EMBL, GenBank) from the main protein endpoint, which may contain genomic location information. Note: The separate genome endpoint is not available, so this tool extracts genome-related cross-references from the main protein data. Use UniProt or EBI Search (domain='uniprot') to find accessions such as 'P05067' or 'P04637'.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'proteins_api_get_genome_mappings',
     'arguments': {'accession': '...'}},
)
```

## `proteins_api_get_xrefs` {#proteins_api_get_xrefs}

Returns database cross-references (Ensembl, RefSeq, EMBL, PDB, Reactome) for UniProt accessions.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `accession`
- **Keywords:** cross reference, id mapping, ensembl, refseq
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string \| array | yes | UniProt protein accession(s). Can be a single accession (e.g., 'P05067'), comma-separated string (e.g., 'P05067,P04637'), or list (e.g., \['P05067', 'P04637'\]). Use UniProt search or ebi_search_domain with domain='uniprot' to find valid accessions. |
| `format` | string | no | Response format. JSON is recommended for most use cases. One of `json`, `xml` Default `json` |

??? note "Upstream description"

    Get cross-references for a protein from the Proteins API, including links to other databases (Ensembl, PDB, Reactome, EMBL, RefSeq, etc.). Supports batch operations: pass a list of accessions or comma-separated string to retrieve multiple proteins in parallel. Note: This tool extracts cross-references from the main protein endpoint (the separate xrefs endpoint is not available). Returns database cross-reference records.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'proteins_api_get_xrefs',
     'arguments': {'accession': '...'}},
)
```

## `QuickGO_annotations_by_gene` {#QuickGO_annotations_by_gene}

Returns GO annotations with evidence codes for a gene product from QuickGO.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_product_id`
- **Keywords:** go annotation, evidence code, quickgo
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_product_id` | string | yes | Gene product identifier, typically a UniProtKB accession. Format: 'UniProtKB:P04637' (TP53), 'UniProtKB:P38398' (BRCA1), 'UniProtKB:P00533' (EGFR). |
| `aspect` | string \| null | no | Filter by GO aspect: 'biological_process', 'molecular_function', or 'cellular_component'. Default: all aspects. |
| `taxon_id` | integer \| null | no | Filter by NCBI taxonomy ID. Example: 9606 (human), 10090 (mouse). |
| `limit` | integer | no | Maximum number of annotations to return. Default: 25, max: 100. Default `25` |

??? note "Upstream description"

    Search Gene Ontology (GO) annotations for a specific gene product using the EBI QuickGO browser. Returns GO terms annotated to the protein/gene with evidence codes, qualifiers, and annotation sources. Accepts UniProtKB accessions (e.g., 'UniProtKB:P04637' for TP53). Provides more detailed annotation data than the standard GO API, including evidence provenance and assigned-by information. Useful for understanding gene function, localization, and biological processes.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'QuickGO_annotations_by_gene',
     'arguments': {'gene_product_id': '...'}},
)
```

## `QuickGO_annotations_by_goterm` {#QuickGO_annotations_by_goterm}

Lists gene products annotated to a GO term with evidence, optionally by taxon.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `go_id`
- **Keywords:** go term, gene products, quickgo
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `go_id` | string | yes | Gene Ontology term accession. Examples: 'GO:0006915' (apoptotic process), 'GO:0008150' (biological process), 'GO:0005634' (nucleus), 'GO:0003677' (DNA binding). |
| `taxon_id` | integer \| null | no | Filter by NCBI taxonomy ID. Example: 9606 (human), 10090 (mouse), 7227 (fly). |
| `limit` | integer | no | Maximum number of annotations to return. Default: 25, max: 100. Default `25` |

??? note "Upstream description"

    Search for all gene products annotated with a specific Gene Ontology (GO) term using the EBI QuickGO browser. Returns the proteins/genes associated with the GO term, their evidence codes, and annotation sources. Can be filtered by species (taxon_id). Useful for finding all genes involved in a biological process, molecular function, or cellular component. Example: GO:0006915 (apoptotic process) returns 6800+ human annotations.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'QuickGO_annotations_by_goterm',
     'arguments': {'go_id': '...'}},
)
```

## `Reactome_get_participant_reference_entities` {#Reactome_get_participant_reference_entities}

Returns reference entities (UniProt proteins, ChEBI chemicals) for participants of a Reactome event.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `stId`
- **Keywords:** participants, uniprot, reaction members
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `stId` | string | yes | Event Stable ID (reaction or pathway, e.g., 'R-HSA-111289'). To find reaction IDs, use Reactome_get_pathway_reactions with a pathway ID. To find pathway IDs, use Reactome_list_top_pathways. |

??? note "Upstream description"

    Get reference entities for participants of a reaction or event. Returns list of reference entities (UniProt proteins, ChEBI chemicals, etc.) associated with the participants, including database identifiers and names.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Reactome_get_participant_reference_entities',
     'arguments': {'stId': '...'}},
)
```

## `Reactome_get_participants` {#Reactome_get_participants}

Returns physical entities (proteins, complexes, small molecules) participating in a Reactome event.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `stId`
- **Keywords:** participants, physical entities, reactome
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `stId` | string | yes | Event Stable ID (reaction or pathway, e.g., 'R-HSA-111289'). To find reaction IDs, use Reactome_get_pathway_reactions with a pathway ID. To find pathway IDs, use Reactome_list_top_pathways. |

??? note "Upstream description"

    Get all participants (physical entities) of a reaction or event. Returns list of physical entities (proteins, complexes, small molecules) participating in the event, including their reference entities.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Reactome_get_participants', 'arguments': {'stId': '...'}},
)
```

## `Reactome_get_pathways_low_entity` {#Reactome_get_pathways_low_entity}

Returns the most specific Reactome pathways containing a given entity.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `id`
- **Keywords:** entity to pathway, lowest level, reactome
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `id` | string | yes | Entity Stable ID (pathway, reaction, or physical entity, e.g., 'R-HSA-73817'). To find pathway IDs, use Reactome_list_top_pathways. To find reaction IDs, use Reactome_get_pathway_reactions. To find entity IDs, use Reactome_get_participants. |

??? note "Upstream description"

    Get low-level pathways (most specific pathways) containing a specific entity. Returns pathways at the lowest level of the hierarchy.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Reactome_get_pathways_low_entity',
     'arguments': {'id': '...'}},
)
```

## `Reactome_map_uniprot_to_pathways` {#Reactome_map_uniprot_to_pathways}

Maps a UniProt accession to all Reactome pathways containing it.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `uniprot_id`
- **Keywords:** uniprot, gene to pathway, reactome
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `uniprot_id` | string | yes | UniProt protein accession (e.g., 'P04637' for TP53, 'P00533' for EGFR) |

??? note "Upstream description"

    Map a UniProt protein identifier to Reactome pathways. Returns all pathways that contain this protein.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Reactome_map_uniprot_to_pathways',
     'arguments': {'uniprot_id': '...'}},
)
```

## `Reactome_map_uniprot_to_reactions` {#Reactome_map_uniprot_to_reactions}

Maps a UniProt accession to Reactome reactions involving it.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `uniprot_id`
- **Keywords:** uniprot, reactions, mapping
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `uniprot_id` | string | yes | UniProt protein accession (e.g., 'P04637' for TP53, 'P00533' for EGFR) |

??? note "Upstream description"

    Map a UniProt protein identifier to Reactome reactions. Returns all reactions that involve this protein.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Reactome_map_uniprot_to_reactions',
     'arguments': {'uniprot_id': '...'}},
)
```

## `ReactomeAnalysis_pathway_enrichment` {#ReactomeAnalysis_pathway_enrichment}

Runs Reactome overrepresentation analysis on gene symbols, UniProt or Ensembl IDs.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `identifiers`
- **Keywords:** enrichment, overrepresentation, fdr, gene list
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `identifiers` | string | yes | Newline-separated list of gene/protein identifiers. Supports gene symbols (TP53, BRCA1), UniProt IDs (P04637, P38398), or Ensembl IDs. Example: 'TP53\\nBRCA1\\nCDH1\\nEGFR\\nKRAS'. |
| `page_size` | integer \| null | no | Number of pathways to return (default 20, max 50). |
| `include_disease` | boolean \| null | no | Include disease pathways in results (default true). |
| `projection` | boolean \| null | no | Project identifiers to human Reactome pathways for cross-species analysis (default true). |

??? note "Upstream description"

    Perform pathway overrepresentation (enrichment) analysis using the Reactome Analysis Service. Submit a list of gene symbols, UniProt IDs, or Ensembl IDs and identify which Reactome pathways are statistically enriched. Returns pathways ranked by p-value with FDR (False Discovery Rate) correction, entity counts, and reaction coverage. Supports projection to human pathways for cross-species analysis. INTERPRETING RESULTS: rank pathways by 'fdr' or 'p_value' - those are the enrichment statistics. 'entities_coverage' (= entities_found / entities_total) is the fraction of a pathway your identifiers cover. 'pathway_size_fraction_of_reactome' (deprecated alias 'entities_ratio') is Reactome's raw 'ratio': pathway size over Reactome's entity universe for the species, independent of the identifiers you submit - never rank by it. Example: submitting TP53, BRCA1, CDH1, EGFR, KRAS returns 394 enriched pathways including 'Developmental Lineage of Mammary Gland Myoepithelial Cells' (p=3.94e-08) and 'Signaling by ERBB2 KD Mutants'.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ReactomeAnalysis_pathway_enrichment',
     'arguments': {'identifiers': '...'}},
)
```

## `ReactomeAnalysis_pathway_found_entities` {#ReactomeAnalysis_pathway_found_entities}

Lists submitted identifiers matching a hit pathway in a Reactome analysis, with cross-references.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `token`, `pathway`
- **Keywords:** found entities, gene list, pathway hit
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `token` | string | yes | Analysis token from a previous ReactomeAnalysis_pathway_enrichment or ReactomeAnalysis_expression_analysis call (the 'token' field in its data). Tokens expire after a period, so run the analysis first. |
| `pathway` | string | yes | Reactome pathway stable id of a hit pathway from the analysis result (e.g. 'R-HSA-3700989'). |
| `resource` | string \| null | no | Resource to map entities against (default 'TOTAL'). Other options: 'UNIPROT', 'ENSEMBL'. |

??? note "Upstream description"

    Per-pathway found-entities drill-down: for an analysis token and a specific HIT pathway, list exactly which submitted identifiers matched that pathway and their Reactome cross-references (UniProt/ENSEMBL ids each input mapsTo). Standard token result tools return only the pathway summary table - they never reveal WHICH of your genes landed in a given pathway. Obtain the token from ReactomeAnalysis_pathway_enrichment / _expression_analysis (data.token), then pass that token plus a pathway stable id. Example: token from submitting PTEN,TP53,EGFR,BRCA1 with pathway 'R-HSA-3700989' returns identifiers including \{id:'BRCA1', mapsTo:\[\{resource:'UNIPROT', ids:\['P38398'\]\},\{resource:'ENSEMBL', ids:\['ENST00000357654','ENSG00000012048'\]\}\]\}.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ReactomeAnalysis_pathway_found_entities',
     'arguments': {'token': '...', 'pathway': '...'}},
)
```

## `ReactomeAnalysis_species_comparison` {#ReactomeAnalysis_species_comparison}

Projects gene/protein identifiers from another species onto human Reactome pathways for enrichment.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `identifiers`
- **Keywords:** orthology projection, enrichment, cross-species
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `identifiers` | string | yes | Newline-separated list of gene/protein identifiers from any species. Supports UniProt, Ensembl, or gene symbols. Example: 'Trp53\\nBrca1\\nCdh1' (mouse gene symbols). |
| `species` | integer \| null | no | NCBI taxonomy ID of the source species (default 9606 for human). Examples: 10090 (mouse), 7955 (zebrafish), 7227 (fly). |
| `page_size` | integer \| null | no | Number of pathways to return (default 20, max 50). |

??? note "Upstream description"

    Perform Reactome pathway analysis with cross-species projection. Submit protein/gene identifiers from any supported species and Reactome will project them to human pathways, enabling comparison of pathway enrichment across species. This uses the same overrepresentation analysis but ensures identifiers are mapped via orthology to human pathways. Example: submitting mouse protein IDs returns enriched human Reactome pathways with mapping statistics. INTERPRETING RESULTS: rank pathways by 'fdr' or 'p_value' - those are the enrichment statistics. 'entities_coverage' (= entities_found / entities_total) is the fraction of a pathway your identifiers cover. 'pathway_size_fraction_of_reactome' (deprecated alias 'entities_ratio') is Reactome's raw 'ratio': pathway size over Reactome's entity universe for the species, independent of the identifiers you submit - never rank by it.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ReactomeAnalysis_species_comparison',
     'arguments': {'identifiers': '...'}},
)
```

## `ReactomeAnalysis_token_result` {#ReactomeAnalysis_token_result}

Retrieves previously computed Reactome analysis results by token with pagination.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `token`
- **Keywords:** analysis token, enrichment results, pagination
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `token` | string | yes | Analysis token from a previous ReactomeAnalysis_pathway_enrichment or ReactomeAnalysis_species_comparison call. |
| `page_size` | integer \| null | no | Number of pathways per page (default 20, max 50). |

??? note "Upstream description"

    Retrieve previously computed Reactome pathway analysis results using an analysis token. Tokens are returned from pathway_enrichment or species_comparison analyses and can be used to access results without rerunning the analysis. Useful for paginating through large result sets or sharing analysis results. The token contains the full analysis result and is valid for a limited time. INTERPRETING RESULTS: rank pathways by 'fdr' or 'p_value' - those are the enrichment statistics. 'entities_coverage' (= entities_found / entities_total) is the fraction of a pathway your identifiers cover. 'pathway_size_fraction_of_reactome' (deprecated alias 'entities_ratio') is Reactome's raw 'ratio': pathway size over Reactome's entity universe for the species, independent of the identifiers you submit - never rank by it.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ReactomeAnalysis_token_result',
     'arguments': {'token': '...'}},
)
```

## `ReactomeInteractors_get_entity_pathways` {#ReactomeInteractors_get_entity_pathways}

Returns Reactome pathways containing a given entity stable ID, flagging disease pathways.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `entity_id`
- **Keywords:** entity pathways, disease pathway, reactome
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `entity_id` | string | yes | Reactome stable identifier for the entity. Use ReactomeInteractors_search_entity to find IDs. Examples: 'R-HSA-199420' (PTEN), 'R-HSA-69541' (TP53 tetramer), 'R-HSA-8869337' (TP53). |
| `species` | integer \| null | no | NCBI taxonomy ID for species filter (default: 9606 for human). Examples: 9606 (human), 10090 (mouse). |

??? note "Upstream description"

    Find Reactome biological pathways that contain a specific molecular entity. Returns pathway stable IDs, names, and whether they are disease-associated. Use Reactome stable IDs (e.g., R-HSA-199420 for PTEN) obtained from ReactomeInteractors_search_entity. Example: R-HSA-199420 (PTEN) participates in 10 pathways including PI3K/AKT signaling, PTEN regulation, and TCR signaling.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ReactomeInteractors_get_entity_pathways',
     'arguments': {'entity_id': '...'}},
)
```

## `resolve_identifier_for_gene_or_protein` {#resolve_identifier_for_gene_or_protein}

Detects the namespace of a gene or protein identifier (symbol, Ensembl, UniProt, Entrez, RefSeq, HGNC) and maps it to others.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `identifier`
- **Keywords:** id mapping, hgnc, entrez, ensembl id, refseq, uniprot
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `identifier` | string | yes | The identifier to resolve, e.g. 'TP53', 'ENSG00000141510', 'P04637', '7157', 'NM_000546', 'HGNC:11998'. |
| `namespace` | string \| null | no | Override namespace detection. Use when the identifier is ambiguous — a bare integer is read as an Entrez Gene ID but could be another namespace. One of `gene_symbol`, `entrez_gene`, `ensembl_gene`, `ensembl_transcript`, `ensembl_protein`, `uniprot`, `refseq_rna`, `refseq_protein`, `hgnc`, `None` |
| `species` | string \| null | no | Species for symbol lookups, common name or NCBI taxonomy ID. Default 'human'. Examples: 'human', 'mouse', '10090'. |
| `tax_id` | integer \| null | no | NCBI taxonomy ID restricting UniProt ID mapping, e.g. 9606 (human), 10090 (mouse). Optional. |
| `sources` | array \| null | no | Restrict which resolvers are queried. Default: all. UniProt ID Mapping submits a job and is the slowest; drop 'uniprot_idmap' for a faster answer. |

??? note "Upstream description"

    Resolve a gene or protein identifier across namespaces when you do NOT know which namespace it belongs to. Detects whether a string like 'TP53', 'ENSG00000141510', 'P04637', '7157' or 'NM_000546' is a symbol, Ensembl, UniProt, Entrez, RefSeq or HGNC id, then queries MyGene, UniProt ID Mapping, g:Profiler, BridgeDb and Ensembl and returns every equivalent id with the sources that agree on it. Use this to obtain the id another tool requires. Use UniProtIDMap_convert_ids, BridgeDb_xrefs or gProfiler_convert_ids directly instead when you already know the source namespace and want one authority. Genes and proteins only: for compounds use PubChem/ChEMBL tools, for ontology terms use ols_get_term_xrefs.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'resolve_identifier_for_gene_or_protein',
     'arguments': {'identifier': '...'}},
)
```

## `Rfam_accession_to_id` {#Rfam_accession_to_id}

Converts an Rfam accession to its family ID/name.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `accession`
- **Keywords:** rfam, ncrna family, id conversion
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `get_family_id` Default `get_family_id` |
| `accession` | string | yes | Required: Rfam accession (e.g., RF00360) |

??? note "Upstream description"

    Convert Rfam accession to family ID/name. Takes RF accession (e.g., RF00360) and returns human-readable family name (e.g., snoZ107_R87). Simple ID conversion utility. Use for: displaying readable names, documentation, cross-referencing, data integration. Quick lookup with plain text response.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Rfam_accession_to_id', 'arguments': {'accession': '...'}},
)
```

## `Rfam_get_family` {#Rfam_get_family}

Returns Rfam family metadata: RNA type, sequence counts, species coverage and search parameters.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `family_id`
- **Keywords:** ncrna, mirna, snorna, rfam family
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `get_family` Default `get_family` |
| `family_id` | string | yes | Required: Rfam accession (e.g., RF00360) or family name (e.g., snoZ107_R87) |
| `format` | string | no | Output format One of `json`, `xml` Default `json` |

??? note "Upstream description"

    Get comprehensive RNA family information from Rfam database. Returns family metadata including description, curation details, RNA type (snoRNA, miRNA, tRNA, etc.), number of sequences, species coverage, search parameters, and structure source. Rfam v15.1 contains 4,227 RNA families. Query by RF accession (e.g., RF00360) or family name (e.g., snoZ107_R87). Use for: RNA family discovery, understanding RNA function, checking curation quality, identifying RNA types. Available formats: json, xml.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Rfam_get_family', 'arguments': {'family_id': '...'}},
)
```

## `Rfam_get_sequence_regions` {#Rfam_get_sequence_regions}

Lists genomic regions (sequence IDs, coordinates, scores) containing Rfam family members.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `family_id`
- **Keywords:** genomic regions, coordinates, ncrna
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `get_sequence_regions` Default `get_sequence_regions` |
| `family_id` | string | yes | Required: Rfam accession or family name |
| `format` | string | no | Output format One of `text`, `json`, `xml` Default `text` |

??? note "Upstream description"

    Get list of all genomic regions containing RNA family members. Returns tab-delimited data with sequence IDs, coordinates, scores, and species. Shows where family is found across genomes. Warning: Some large families return 403 error due to too many regions. Use for: genome annotation, identifying RNA locations, cross-species comparisons, coverage analysis. Available formats: text (tab-delimited), json, xml.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Rfam_get_sequence_regions',
     'arguments': {'family_id': '...'}},
)
```

## `Rfam_id_to_accession` {#Rfam_id_to_accession}

Converts an Rfam family ID/name to its RF accession.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `family_id`
- **Keywords:** rfam, accession, id conversion
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `get_family_accession` Default `get_family_accession` |
| `family_id` | string | yes | Required: Family name (e.g., snoZ107_R87) |

??? note "Upstream description"

    Convert Rfam family ID/name to accession. Takes family name (e.g., snoZ107_R87) and returns RF accession (e.g., RF00360). Simple ID conversion utility. Use for: getting stable identifiers, cross-referencing, data integration, creating accession-based links. Accessions are stable across Rfam versions.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Rfam_id_to_accession', 'arguments': {'family_id': '...'}},
)
```

## `RGD_get_annotations` {#RGD_get_annotations}

Returns RGD disease, phenotype, pathway and GO annotations for a rat gene.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `rgd_id`
- **Keywords:** rgd, rat gene, go, disease annotation
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `rgd_id` | string | yes | RGD gene ID (e.g., '2218' for Brca1) |
| `aspect` | string | no | Filter by annotation type: 'D' (Disease), 'P' (Pathway), 'W' (Phenotype), 'F' (Function), 'C' (Component), 'E' (Expression). Leave empty for all. Default `` |

??? note "Upstream description"

    Get disease, phenotype, pathway, and GO annotations for a rat gene from RGD. Returns curated annotations with ontology terms, evidence codes, and qualifiers. Filter by aspect: D=Disease, P=Pathway, W=Phenotype, F=Molecular Function, C=Cellular Component, E=Expression. Essential for rat disease model research and toxicogenomics.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'RGD_get_annotations', 'arguments': {'rgd_id': '...'}},
)
```

## `RGD_get_gene` {#RGD_get_gene}

Returns RGD rat gene symbol, name, description, type and Ensembl cross-reference.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `rgd_id`
- **Keywords:** rgd, rat gene, ensembl
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `rgd_id` | string | yes | RGD gene ID (numeric, e.g., '2218' for Brca1). Can include 'RGD:' prefix. |

??? note "Upstream description"

    Get detailed rat gene information from the Rat Genome Database by RGD ID. Returns gene symbol, name, description (including AGR-aggregated functional summary with disease associations, pathways, and interactions), gene type, and Ensembl cross-reference. Use when looking up specific rat genes for toxicology, pharmacology, or disease model research.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'RGD_get_gene', 'arguments': {'rgd_id': '...'}},
)
```

## `RGD_get_orthologs` {#RGD_get_orthologs}

Returns orthologs of a rat gene across species from RGD.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `rgd_id`
- **Keywords:** ortholog, rat, cross-species
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `rgd_id` | string | yes | RGD gene ID (e.g., '2218' for Brca1) |

??? note "Upstream description"

    Get orthologs of a rat gene across species (human, mouse, chinchilla, dog, pig, etc.) from RGD. Returns ortholog gene symbols, RGD IDs, names, and species. Use for cross-species comparison in model organism genetics, finding human disease gene orthologs in rat.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'RGD_get_orthologs', 'arguments': {'rgd_id': '...'}},
)
```

## `RGD_resolve_symbol_or_region` {#RGD_resolve_symbol_or_region}

Resolves a rat gene symbol to its RGD record or lists genes overlapping a genomic region.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** rgd, gene symbol, genomic region
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `symbol` | string \| null | no | Rat gene symbol, e.g. 'Tp53'. Use this for symbol mode. |
| `species_type_key` | integer \| null | no | Species type key for symbol mode (default 3=rat, 1=human, 2=mouse). |
| `chromosome` | string \| null | no | Chromosome for region mode, e.g. '10'. |
| `start` | integer \| null | no | Region start (region mode). |
| `stop` | integer \| null | no | Region stop (region mode). |
| `map_key` | integer \| null | no | Assembly map key for region mode (e.g. 360 for rat rn7). |

??? note "Upstream description"

    Resolve a rat gene symbol to its native RGD record/ID, or list genes overlapping a genomic region. Symbol mode: provide 'symbol' (and optional species_type_key, default 3=rat) -&gt; returns RGD key, rgd_id, name, type. Region mode: provide chromosome + start + stop + map_key -&gt; returns genes in that interval. Example: symbol='Tp53' returns key=1887 rgd_id=3889 name='tumor protein p53'; region chromosome='10', start=1, stop=5000000, map_key=360 returns Mir484 rgdId=2325341.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'RGD_resolve_symbol_or_region', 'arguments': {}},
)
```

## `RGD_search_genes` {#RGD_search_genes}

Searches RGD genes by symbol or keyword for rat, human or mouse.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** rgd, gene search, rat
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Gene symbol or keyword (e.g., 'Brca1', 'Tp53', 'insulin') |
| `species` | string | no | Species to search: 'rat' (default), 'human', 'mouse' Default `rat` |

??? note "Upstream description"

    Search for rat genes by symbol or keyword in the Rat Genome Database. Returns matching genes with RGD IDs, symbols, names, types, and descriptions. Supports rat, human, and mouse species. Use for finding rat orthologs, toxicology studies, or rat disease model genes.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'RGD_search_genes', 'arguments': {'query': '...'}},
)
```

## `RNAcentral_get_by_accession` {#RNAcentral_get_by_accession}

Returns an RNAcentral ncRNA entry with annotations and source cross-references.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `accession`
- **Keywords:** ncrna, rnacentral, cross reference
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | RNAcentral accession (e.g., 'URS000075C808'). |

??? note "Upstream description"

    Retrieve a single RNAcentral entry by accession for detailed annotations and source cross-references. Use for cross-database ID mapping and metadata.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'RNAcentral_get_by_accession',
     'arguments': {'accession': '...'}},
)
```

## `RNAcentral_get_genome_locations` {#RNAcentral_get_genome_locations}

Returns chromosome, strand and coordinates of an RNAcentral ncRNA in a species.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `urs_id`, `taxid`
- **Keywords:** genomic location, ncrna, strand, urs
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Fixed to 'genome_locations' for this tool. One of `genome_locations` Default `genome_locations` |
| `urs_id` | string | yes | RNAcentral URS identifier, e.g. 'URS00003B7674'. A trailing '_taxid' suffix is accepted and stripped. |
| `taxid` | integer | yes | NCBI taxonomy id of the organism whose genome mapping is requested, e.g. 9606 (human), 10090 (mouse). Required. |

??? note "Upstream description"

    Get the genomic coordinates of a non-coding RNA in a specific organism from RNAcentral. Provide an RNAcentral URS identifier (e.g. 'URS00003B7674' for let-7a) and an NCBI taxonomy id (e.g. 9606 for human). Returns each genomic mapping: chromosome, strand, start, end, UCSC/Ensembl assembly (GRCh38), and GCA accession. RNAcentral aggregates ncRNAs (miRNA, lncRNA, rRNA, snoRNA, tRNA, etc.) from 40+ member databases.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'RNAcentral_get_genome_locations',
     'arguments': {'urs_id': '...', 'taxid': 0}},
)
```

## `RNAcentral_get_region_ncRNAs` {#RNAcentral_get_region_ncRNAs}

Lists RNAcentral non-coding RNAs overlapping a genomic region for a species.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** ncrna, genomic region, lncrna, mirna
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `region` | string \| null | no | Locus as 'chr:start-end' (1-based), e.g. '2:39745816-39826679'. Omit if giving chromosome/start/end. |
| `species` | string \| null | no | Ensembl species slug (default 'homo_sapiens'). E.g. 'mus_musculus', 'danio_rerio'. |
| `chromosome` | string \| null | no | Chromosome name without 'chr' prefix, e.g. '2'. Used with start/end when 'region' is omitted. |
| `start` | integer \| null | no | Region start (1-based). Used with chromosome/end when 'region' is omitted. |
| `end` | integer \| null | no | Region end (1-based). Used with chromosome/start when 'region' is omitted. |

??? note "Upstream description"

    List all RNAcentral non-coding RNAs overlapping a genome locus, aggregated from 50+ member databases (miRBase, snoRNA, tRNA, lncRNA, piRNA, etc.). Provide a species (default 'homo_sapiens') and either a 'region' string ('2:39745816-39826679') or discrete 'chromosome'/'start'/'end'. Returns overlapping ncRNA transcripts and their exon features with biotype, strand and coordinates, split into data.transcripts and data.exons. Use to annotate what ncRNAs sit in a region of interest. Public, no authentication.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'RNAcentral_get_region_ncRNAs', 'arguments': {}},
)
```

## `RNAcentral_get_xrefs_and_pubs` {#RNAcentral_get_xrefs_and_pubs}

Returns per-source cross-references (NCBI gene, Ensembl/GENCODE transcripts) and publications for an RNAcentral entry.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `accession`
- **Keywords:** ncrna, gencode, ncbi gene, xrefs
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | RNAcentral URS identifier (e.g., 'URS000063A371', 'URS000075C808'). |

??? note "Upstream description"

    Get the resolved per-source cross-references AND literature for an RNAcentral entry in one call. Returns the full xref list (with database name, taxid, accession id, NCBI gene id, miRBase mature/precursor products, Ensembl/GENCODE transcript ids, RefSeq, etc.) under data.xrefs, and the publications list (title, authors, PMID, DOI) under data.publications. Use after RNAcentral_search or RNAcentral_get_by_accession to map an URS id to external databases and supporting papers. Example: URS000063A371 returns 8 xrefs (database=Rfam, taxid=7955, accession id 'CR381647.10:84955..85035:rfam') and 8 publications (e.g. PMID 12844358, DOI 10.1186/gb-2003-4-7-r42).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'RNAcentral_get_xrefs_and_pubs',
     'arguments': {'accession': '...'}},
)
```

## `RNAcentral_search` {#RNAcentral_search}

Searches RNAcentral ncRNA records (miRNA, lncRNA, rRNA) by query.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** ncrna, mirna, lncrna, search
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Keyword, accession, or sequence-based query (per RNAcentral API). |
| `page_size` | integer | no | Number of records per page (1–100). Default `10` |

??? note "Upstream description"

    Search aggregated ncRNA records (miRNA, rRNA, lncRNA, etc.) across sources via RNAcentral. Use to find accessions, types, species, and descriptions for ncRNA analysis.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'RNAcentral_search', 'arguments': {'query': '...'}},
)
```

## `SGD_get_gene` {#SGD_get_gene}

Gets a yeast gene/locus overview from SGD by SGD ID, including names, description, aliases, UniProt ID and functional annotations.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `sgd_id`
- **Keywords:** yeast, saccharomyces cerevisiae, locus, systematic name, gene description
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sgd_id` | string | yes | SGD locus identifier. Examples: 'S000000259' (GAL4), 'S000003219' (RMR1), 'S000005564' (ACT1). |

??? note "Upstream description"

    Get yeast gene/locus overview from SGD (Saccharomyces Genome Database) by SGD identifier. Returns gene name, systematic name, description, UniProt ID, aliases, and functional annotations. S. cerevisiae is the primary model organism for eukaryotic cell biology. Example: S000000259 returns GAL4, the galactose metabolism transcription factor.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SGD_get_gene', 'arguments': {'sgd_id': '...'}},
)
```

## `SGD_get_go_annotations` {#SGD_get_go_annotations}

Gets Gene Ontology annotations (process, function, component) with evidence codes for a yeast gene from SGD.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `sgd_id`
- **Keywords:** yeast, gene ontology, go terms, evidence codes, gene function
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sgd_id` | string | yes | SGD locus identifier. Examples: 'S000003219' (RMR1), 'S000005564' (ACT1), 'S000000259' (GAL4). |

??? note "Upstream description"

    Get Gene Ontology annotations for a yeast gene from SGD. Returns GO terms (biological process, molecular function, cellular component), evidence codes, qualifiers, and references. SGD GO annotations are among the most thoroughly curated in any model organism. Example: S000003219 (RMR1) has annotations for meiotic recombination and cytoplasm localization.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SGD_get_go_annotations', 'arguments': {'sgd_id': '...'}},
)
```

## `SGD_get_interactions` {#SGD_get_interactions}

Gets genetic and physical interactions for a yeast gene from SGD, with experiment types, interacting genes and references.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `sgd_id`
- **Keywords:** yeast, genetic interaction, synthetic lethality, two-hybrid, interactome
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sgd_id` | string | yes | SGD locus identifier. Examples: 'S000003219' (RMR1), 'S000005564' (ACT1), 'S000000259' (GAL4). |

??? note "Upstream description"

    Get genetic and physical interactions for a yeast gene from SGD. Returns interaction type (Physical/Genetic), experiment type (Two-hybrid, Affinity Capture-MS, Synthetic Lethality, etc.), interacting gene names, and references. Yeast interaction data is the most comprehensive of any organism. Example: S000003219 (RMR1) has 178 interactions.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SGD_get_interactions', 'arguments': {'sgd_id': '...'}},
)
```

## `SGD_search` {#SGD_search}

Searches the Saccharomyces Genome Database across genes, GO terms, phenotypes, alleles, diseases and more, with category filtering.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** yeast, saccharomyces genome database, gene search, alleles, locus
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search query string. Examples: 'GAL4', 'cell cycle', 'sporulation'. |
| `category` | string | no | Filter by category. Options: 'locus', 'reference', 'biological_process', 'molecular_function', 'cellular_component', 'phenotype', 'disease', 'allele', 'complex', 'pathway', 'chemical', 'strain', 'download'. |
| `limit` | integer | no | Maximum results to return (default 10, max 50). |
| `offset` | integer | no | Pagination offset (default 0). |

??? note "Upstream description"

    Search SGD across all data types (genes, GO terms, phenotypes, references, alleles, diseases, complexes, pathways). Returns matching entries with categories, descriptions, and links. Supports category filtering. Example: searching 'GAL4' finds the transcription factor gene, related references, alleles, and GO terms.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SGD_search', 'arguments': {'query': '...'}},
)
```

## `STRING_functional_enrichment` {#STRING_functional_enrichment}

Runs STRING enrichment of a gene/protein set against GO, KEGG and Reactome.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `protein_ids`
- **Keywords:** enrichment, go, kegg, fdr
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `protein_ids` | array of string | yes | List of protein identifiers (UniProt IDs, gene names, Ensembl IDs). Minimum 3 proteins recommended for meaningful enrichment. |
| `species` | integer | no | NCBI taxonomy ID (default: 9606 for human) Default `9606` |
| `category` | string | no | Enrichment category: 'Process' (GO Biological Process), 'Component' (GO Cellular Component), 'Function' (GO Molecular Function), 'KEGG', 'Reactome', 'WikiPathways', 'COMPARTMENTS', 'TISSUES', 'DISEASES' One of `Process`, `Component`, `Function`, `KEGG`, `Reactome`, `WikiPathways`, `COMPARTMENTS`, `TISSUES`, `DISEASES` Default `Process` |

??? note "Upstream description"

    Identify enriched biological functions, pathways, and processes for a protein set using STRING (Search Tool for Retrieval of Interacting Genes/Proteins) database. Performs statistical enrichment analysis against GO (Gene Ontology) terms, KEGG pathways, Reactome pathways, and other annotation databases. Returns significantly enriched terms with FDR-corrected p-values. Minimum 3-4 proteins required for meaningful statistical analysis. No API key required. Use for: discovering pathways enriched in protein network, identifying biological processes, finding shared functions among proteins, interpreting omics data (proteomics, interactomics), hypothesis generation.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'STRING_functional_enrichment',
     'arguments': {'protein_ids': []}},
)
```

## `STRING_get_enrichment` {#STRING_get_enrichment}

Runs STRING functional enrichment on a list of gene or protein identifiers and returns enriched GO, pathway and other terms with FDR.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `identifiers`
- **Keywords:** gene set enrichment, go enrichment, kegg, reactome, gene list, fdr
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `identifiers` | string | yes | Newline-separated list of protein identifiers (gene names or STRING IDs). Example: 'BRCA1\\nBRCA2\\nTP53\\nATM'. Minimum 2 proteins, recommended 5-100 for meaningful enrichment. |
| `species` | integer | no | NCBI taxonomy ID (default 9606 for human). Examples: 9606 (human), 10090 (mouse), 7227 (fly), 6239 (worm). Default `9606` |
| `background_string_identifiers` | string \| null | no | Optional newline-separated background gene set. If not provided, the entire genome proteome is used as background. |

??? note "Upstream description"

    Perform functional enrichment analysis on a set of proteins using the STRING database. Tests whether the input protein list is statistically enriched for Gene Ontology terms, KEGG pathways, Reactome pathways, WikiPathways, InterPro domains, TISSUES, DISEASES, COMPARTMENTS, and other functional categories. Returns enriched terms with p-values, FDR-corrected q-values, and the matching proteins. Useful for interpreting gene lists from RNA-seq, GWAS, proteomics experiments.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'STRING_get_enrichment',
     'arguments': {'identifiers': '...'}},
)
```

## `STRING_get_functional_annotations` {#STRING_get_functional_annotations}

Returns all STRING functional annotations (GO, KEGG, Reactome, domains, phenotypes, diseases, tissues) for a single gene or protein.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `identifiers`
- **Keywords:** gene ontology, pathways, gene function, hpo, tissues, annotation
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `identifiers` | string | yes | Protein identifier (gene name or STRING ID). Examples: 'TP53', 'BRCA1', 'EGFR', '9606.ENSP00000269305'. |
| `species` | integer | no | NCBI taxonomy ID (default 9606 for human). Examples: 9606 (human), 10090 (mouse), 7227 (fly), 6239 (worm). Default `9606` |
| `category` | string | no | Filter by annotation category. Options: 'Process' (GO BP), 'Function' (GO MF), 'Component' (GO CC), 'KEGG', 'RCTM' (Reactome), 'WikiPathways', 'InterPro', 'Pfam', 'HPO', 'DISEASES', 'TISSUES', 'COMPARTMENTS', 'Keyword'. Leave null for all categories. |

??? note "Upstream description"

    Get comprehensive functional annotations for a protein from STRING database. Returns per-protein annotations from Gene Ontology (Process, Function, Component), KEGG pathways, Reactome, WikiPathways, InterPro domains, Pfam, HPO phenotypes, DISEASES, TISSUES expression, and COMPARTMENTS subcellular localization. Different from STRING_functional_enrichment which tests for statistical enrichment of a protein SET - this tool returns all known annotations for individual proteins. Example: 'TP53' returns 1,877 annotations across 13 categories including 655 biological processes, 462 HPO phenotypes, 212 disease associations.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'STRING_get_functional_annotations',
     'arguments': {'identifiers': '...'}},
)
```

## `TargetMine_search_genes` {#TargetMine_search_genes}

Searches TargetMine for genes and returns gene symbols, NCBI Gene IDs, names and organisms.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `q`
- **Keywords:** gene search, ncbi gene id, intermine, drug target
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `q` | string | yes | Gene search query: gene symbol (e.g., 'TP53', 'EGFR', 'BRCA1'), Entrez Gene ID, or keyword. |
| `size` | integer | no | Number of results to return (default: 10). Default `10` |

??? note "Upstream description"

    Search TargetMine specifically for genes relevant to drug target discovery. This convenience wrapper filters results to the Gene category. Returns gene symbols, NCBI Gene IDs, organism information, and gene names. TargetMine covers human, rat, mouse, and other model organism genes with extensive drug target annotations. For broader searches across all entity types, use TargetMine_search instead.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'TargetMine_search_genes', 'arguments': {'q': '...'}},
)
```

## `Tark_get_mane_transcripts` {#Tark_get_mane_transcripts}

Looks up MANE Select and MANE Plus Clinical transcripts for a gene, pairing Ensembl ENST with RefSeq NM IDs.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** mane select, enst, refseq nm, transcript, ensembl
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string \| null | no | HGNC gene symbol, e.g. 'BRCA2', 'CFTR'. |
| `ensembl_id` | string \| null | no | Ensembl transcript id, e.g. 'ENST00000380152' (version optional). |
| `refseq_id` | string \| null | no | RefSeq transcript id, e.g. 'NM_000059' (version optional). |

??? note "Upstream description"

    Look up the MANE (Matched Annotation from NCBI and EMBL-EBI) transcript(s) for a gene, mapping the Ensembl transcript (ENST) to its equivalent RefSeq transcript (NM_). Returns MANE Select and, where defined, MANE Plus Clinical pairings from the Ensembl Tark transcript archive. Query by HGNC gene symbol, an Ensembl transcript id, or a RefSeq transcript id (versions optional — they are matched on the base id). Use to find the clinically-recommended representative transcript for a gene and to convert between RefSeq and Ensembl transcript namespaces. No API key required.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Tark_get_mane_transcripts', 'arguments': {}},
)
```

## `Tark_get_transcript` {#Tark_get_transcript}

Gets the archived Ensembl Tark record for a transcript ID including assembly, biotype, genomic location and sequence checksum.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `stable_id`
- **Keywords:** ensembl transcript, enst, genomic location, biotype, release history
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `stable_id` | string | yes | Ensembl transcript stable id, e.g. 'ENST00000380152' (version optional). |

??? note "Upstream description"

    Get the archived transcript record for an Ensembl transcript id (ENST) from the Ensembl Tark transcript archive: assembly, biotype, genomic location (region/start/end/strand), the transcript sequence checksum, and the Ensembl release versions in which the transcript appears. Use to resolve transcript versions and verify transcript identity/equivalence across releases. No API key required.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Tark_get_transcript', 'arguments': {'stable_id': '...'}},
)
```

## `TogoID_convert` {#TogoID_convert}

Converts biological identifiers between databases (e.g. Ensembl gene to UniProt, NCBI Gene to HGNC) via TogoID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ids`, `source`, `target`
- **Keywords:** id conversion, ensembl gene, ncbi gene, hgnc, identifier mapping
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ids` | string | yes | Source identifier(s), comma-separated, in the source dataset's native format (e.g. 'ENSG00000139618', '675', '2244'). |
| `source` | string | yes | Source TogoID dataset name, e.g. 'ensembl_gene', 'ncbigene', 'uniprot', 'pubchem_compound', 'hgnc'. |
| `target` | string | yes | Target TogoID dataset name, e.g. 'uniprot', 'hgnc', 'chebi', 'pdb', 'ensembl_transcript'. |

??? note "Upstream description"

    Convert biological identifiers between databases using TogoID (DBCLS), which links 117 ID types (genes, proteins, transcripts, compounds, diseases, pathways, variants, cell lines, organisms, etc.). Provide the source ID(s), the source dataset name, and the target dataset name. Works for directly-related dataset pairs (e.g. ensembl_gene-&gt;uniprot, ncbigene-&gt;hgnc, pubchem_compound-&gt;chebi, uniprot-&gt;pdb); unrelated pairs return a 'no route' error — convert via an intermediate. Call TogoID_list_datasets for the exact dataset names. No API key required.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'TogoID_convert',
     'arguments': {'ids': '...', 'source': '...', 'target': '...'}},
)
```

## `UCSC_get_track` {#UCSC_get_track}

Gets UCSC Genome Browser annotation track records (genes, regulatory elements, conservation, variants, repeats) for a genomic region.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `genome`, `track`, `chrom`, `start`, `end`
- **Keywords:** annotation track, knowngene, genomic region, conservation, dbsnp
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `genome` | string | yes | UCSC genome assembly. Examples: 'hg38', 'hg19', 'mm39'. |
| `track` | string | yes | Track name. Common tracks: 'knownGene' (GENCODE genes), 'ncbiRefSeq' (RefSeq), 'snp151' (dbSNP), 'encodeCcreCombined' (ENCODE cCREs), 'phastCons100way' (conservation), 'rmsk' (repeats). |
| `chrom` | string | yes | Chromosome name. Examples: 'chr1', 'chr17', 'chrX'. |
| `start` | integer | yes | 0-based start position (inclusive). |
| `end` | integer | yes | End position (exclusive). Must be &gt; start. |
| `maxItemsOutput` | integer \| null | no | Maximum number of items to return. Default: 100. Set to limit large result sets. |

??? note "Upstream description"

    Get annotation track data from the UCSC Genome Browser for a specified genomic region. Retrieves features like genes, regulatory elements, conservation scores, variants, or repeats from any of UCSC's 370+ tracks. Use UCSC_search first to find the genomic coordinates of your gene of interest. Example: track 'knownGene' on hg38 chr17:7668421-7687490 returns GENCODE gene annotations for the TP53 region.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'UCSC_get_track',
     'arguments': {'genome': '...',
                   'track': '...',
                   'chrom': '...',
                   'start': 0,
                   'end': 0}},
)
```

## `UCSC_list_tracks` {#UCSC_list_tracks}

Lists available UCSC Genome Browser annotation tracks for a genome assembly or shows a track's column schema.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `genome`
- **Keywords:** annotation tracks, genome assembly, track schema, hg38
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `genome` | string | yes | UCSC genome assembly identifier. Examples: 'hg38', 'hg19', 'mm39', 'mm10', 'dm6', 'danRer11'. |
| `track` | string \| null | no | Optional. When set, return this track's column schema instead of the track list. Example: 'knownGene', 'dbSnp155', 'ncbiRefSeq'. |
| `name_filter` | string \| null | no | Optional substring filter (case-insensitive) applied to track name/shortLabel/longLabel when listing tracks. Example: 'dbSnp', 'gencode', 'clinvar'. |
| `max_tracks` | integer \| null | no | Maximum number of tracks to return when listing (default 500). The total matched count is always reported in track_count. |

??? note "Upstream description"

    Discover available UCSC Genome Browser annotation tracks for a genome assembly, and inspect a track's column schema. With only 'genome' set, returns the list of leaf tracks (each with track name, type, shortLabel, longLabel, parent) so you can find valid track names to pass to UCSC_get_track (there are tens of thousands, e.g. dbSnp155, dbSnp155ClinVar, dbSnp155Common, knownGene, ncbiRefSeq, encodeCcreCombined). Use 'name_filter' to narrow the list by substring (matched against track name, shortLabel, longLabel). When 'track' is also provided, returns that track's column schema instead (column name, sqlType, jsonType, description) from the list/schema endpoint. Example: genome 'hg38' returns 23000+ tracks; genome 'hg38' + track 'knownGene' returns its 27 columns (chrom, chromStart, ...).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'UCSC_list_tracks', 'arguments': {'genome': '...'}},
)
```

## `UCSC_search` {#UCSC_search}

Searches the UCSC Genome Browser for genes, transcripts or features by name and returns genomic positions.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `search_term`
- **Keywords:** gene position, genomic coordinates, refseq, gencode, mane
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `search_term` | string | yes | Gene symbol, transcript ID, SNP rsID, or genomic feature to search. Examples: 'TP53', 'BRCA1', 'rs429358', 'ENST00000269305'. |
| `genome` | string | no | UCSC genome assembly identifier. Default: 'hg38'. Others: 'hg19', 'mm10', 'mm39', 'dm6', 'ce11', 'danRer11', 'sacCer3'. |

??? note "Upstream description"

    Search the UCSC Genome Browser for genes, transcripts, or genomic features by name. Returns position matches across tracks (MANE, RefSeq, GENCODE, HGNC, etc.) for a given genome assembly. Supports searching for gene symbols, transcript IDs, SNP rsIDs, or chromosomal bands. Example: searching 'TP53' in hg38 returns chr17:7668421-7687490 from MANE Select and other tracks.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'UCSC_search', 'arguments': {'search_term': '...'}},
)
```

## `UniProt_id_mapping` {#UniProt_id_mapping}

Maps identifiers between databases (e.g. Ensembl or gene name to UniProt) with batch and async support.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ids`, `from_db`
- **Keywords:** id mapping, ensembl, gene name, identifier conversion
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ids` | any | yes | ID(s) to map. Can be single string or array of strings, e.g., 'ENSG00000141510' or \['MEIOB', 'TP53'\] |
| `from_db` | string | yes | Source database. Examples: 'Ensembl', 'Gene_Name', 'RefSeq_Protein', 'PDB', 'EMBL' |
| `to_db` | string | no | Target database (default: 'UniProtKB') |
| `max_wait_time` | integer | no | Maximum time to wait for async task completion in seconds (default: 30) |

??? note "Upstream description"

    Map IDs between different databases (e.g., Ensembl to UniProt, Gene Name to UniProt). Supports batch mapping and async task processing. Use this to convert between different database identifiers.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'UniProt_id_mapping',
     'arguments': {'ids': '...', 'from_db': '...'}},
)
```

## `UniProtIDMap_convert_ids` {#UniProtIDMap_convert_ids}

Converts protein or gene identifiers between 100+ databases via the UniProt ID Mapping service.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ids`, `from_db`
- **Keywords:** id mapping, ensembl, refseq, gene names, identifier conversion
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ids` | string | yes | Comma-separated identifiers to convert. Examples: 'TP53,BRCA1,EGFR' (gene names), 'P04637,P38398' (UniProt ACs), 'ENSG00000141510' (Ensembl). |
| `from_db` | string | yes | Source database. Common values: 'Gene_Name', 'UniProtKB_AC-ID', 'Ensembl', 'Ensembl_Genomes', 'RefSeq_Protein', 'PDB', 'ChEMBL', 'GeneCards', 'HGNC'. See UniProt ID Mapping docs for full list. |
| `to_db` | string | no | Target database. Common values: 'UniProtKB', 'UniProtKB-Swiss-Prot', 'Gene_Name', 'Ensembl', 'RefSeq_Protein', 'PDB', 'ChEMBL', 'GeneCards'. Default: 'UniProtKB'. Default `UniProtKB` |
| `tax_id` | integer \| null | no | NCBI Taxonomy ID to restrict results. Examples: 9606 (human), 10090 (mouse), 10116 (rat). Optional but recommended for gene names. |

??? note "Upstream description"

    Convert protein/gene identifiers between databases using the UniProt ID Mapping service. Submits a mapping job and polls for results. Supports 100+ databases including UniProtKB, Gene Names, Ensembl, RefSeq, PDB, ChEMBL, and more. The canonical service for cross-database identifier conversion in proteomics and genomics. Example: convert gene names TP53,BRCA1,EGFR to UniProtKB accessions (P04637, P38398, P00533). Returns mapped IDs with from/to pairs.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'UniProtIDMap_convert_ids',
     'arguments': {'ids': '...', 'from_db': '...'}},
)
```

## `UniProtIDMap_gene_to_uniprot` {#UniProtIDMap_gene_to_uniprot}

Converts gene names or symbols to UniProt accessions for a given organism via UniProt ID Mapping.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_names`
- **Keywords:** gene symbol, gene to protein, id mapping, accession
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_names` | string | yes | Comma-separated gene symbols. Examples: 'TP53', 'BRCA1,TP53,EGFR,KRAS', 'INS,GCG'. |
| `tax_id` | integer | no | NCBI Taxonomy ID restricting the species searched. Optional; when omitted it defaults to 9606 (human), so a non-human query MUST pass this explicitly or it will silently return human proteins. Examples: 9606 = human, 10090 = mouse, 10116 = rat, 7227 = fruit fly, 6239 = C. elegans. There is no species-name parameter -- give the numeric taxonomy ID. Default `9606` |
| `reviewed_only` | boolean \| null | no | If true, return only Swiss-Prot (reviewed) entries. Default: false (all UniProtKB entries). |

??? note "Upstream description"

    Convert gene names or symbols to UniProt accessions using the UniProt ID Mapping service. Shortcut for the most common ID mapping task in bioinformatics: going from gene symbols to UniProt protein IDs. Restricts to a species via taxonomy ID to avoid ambiguous results. Example: 'TP53' in human (9606) maps to P04637 (canonical) plus isoforms.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'UniProtIDMap_gene_to_uniprot',
     'arguments': {'gene_names': '...'}},
)
```

## `VariantValidator_gene2transcripts` {#VariantValidator_gene2transcripts}

Gets all transcripts for a gene with MANE annotations, Ensembl/RefSeq IDs and genomic coordinates via VariantValidator.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_symbol`
- **Keywords:** transcripts, mane select, refseq, ensembl, gene coordinates
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_symbol` | string | yes | HGNC gene symbol (e.g., 'TP53', 'BRCA1', 'EGFR'). Aliases: gene, gene_name. |
| `gene` | string | no | Alias for gene_symbol. |
| `gene_name` | string | no | Alias for gene_symbol. |
| `transcript_set` | string | no | Transcript filter: 'mane' for MANE Select/Plus Clinical only, 'refseq' for RefSeq transcripts, 'ensembl' for Ensembl, 'all' for everything Default `mane` |
| `genome_build` | string | no | Reference genome assembly: 'GRCh37' or 'GRCh38' (default: GRCh38) Default `GRCh38` |

??? note "Upstream description"

    Get all transcripts for a gene using VariantValidator. Returns MANE Select and MANE Plus Clinical transcript annotations, Ensembl and RefSeq transcript identifiers, genomic coordinates, and cross-references. Useful for identifying the correct transcript for variant annotation. Supports filter for MANE transcripts only or all transcripts.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'VariantValidator_gene2transcripts',
     'arguments': {'gene_symbol': '...'}},
)
```

## `VEuPathDB_get_gene_record` {#VEuPathDB_get_gene_record}

Retrieves a gene record from a VEuPathDB pathogen-genome site (e.g. PlasmoDB) by gene ID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** pathogen genomics, plasmodb, gene record, parasite
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | Gene identifier / primary key (source_id) as used by VEuPathDB. Examples: 'PF3D7_0417200' (PlasmoDB), 'TGME49_200010' (ToxoDB). |
| `project` | string \| null | no | VEuPathDB member site that owns the gene. One of: plasmodb (default), toxodb, fungidb, vectorbase, cryptodb, giardiadb, microsporidiadb, piroplasmadb, trichdb, tritrypdb, amoebadb. |
| `attributes` | array \| string \| null | no | Optional list of gene attribute columns to return. Defaults to primary_key, product, organism, gene_type, location_text, source_id. Other valid columns include chromosome, transcript_count, is_pseudo. |

??? note "Upstream description"

    Retrieve the record for a single gene from a VEuPathDB pathogen-genome database by its gene identifier (primary key), via the EuPathDB WDK gene records endpoint (POST query). Pick the member site with the 'project' parameter (plasmodb, toxodb, fungidb, vectorbase, cryptodb, giardiadb, microsporidiadb, piroplasmadb, trichdb, tritrypdb, amoebadb). Returns gene attributes: protein product, organism, gene type, genomic location, chromosome and transcript count. Example: gene_id='PF3D7_0417200', project='plasmodb' returns the bifunctional dihydrofolate reductase-thymidylate synthase (DHFR-TS) gene of Plasmodium falciparum 3D7.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'VEuPathDB_get_gene_record',
     'arguments': {'gene_id': '...'}},
)
```

## `VEuPathDB_search_genes_by_organism` {#VEuPathDB_search_genes_by_organism}

Searches a VEuPathDB site for all genes of a given organism or strain.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `organism`
- **Keywords:** pathogen genomics, gene list, organism, plasmodium
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `organism` | string | yes | Full organism/strain name as used by VEuPathDB. Examples: 'Plasmodium falciparum 3D7', 'Toxoplasma gondii ME49', 'Aspergillus fumigatus Af293'. Must match the site's organism vocabulary (genus + species + strain). |
| `project` | string \| null | no | VEuPathDB member site to query. One of: plasmodb (default), toxodb, fungidb, vectorbase, cryptodb, giardiadb, microsporidiadb, piroplasmadb, trichdb, tritrypdb, amoebadb. |
| `attributes` | array \| string \| null | no | Optional list of gene attribute columns to return. Defaults to primary_key, product, organism, gene_type, location_text, source_id. Other valid columns include chromosome, transcript_count, is_pseudo. |
| `limit` | integer \| null | no | Maximum number of genes to return (1-500, default 25). |

??? note "Upstream description"

    Search a VEuPathDB pathogen-genome database for all genes belonging to a given organism/strain, using the EuPathDB WDK GenesByTaxonGene search (POST query). VEuPathDB is a family of sites that share one API: PlasmoDB (malaria, Plasmodium), ToxoDB (Toxoplasma), FungiDB (fungi), VectorBase (disease vectors), CryptoDB, GiardiaDB, MicrosporidiaDB, PiroplasmaDB, TrichDB, TriTrypDB and AmoebaDB. Choose the site with the 'project' parameter and give the full organism name (genus + species + strain). Returns one row per gene with primary key (gene ID), protein product, organism, gene type and genomic location. Example: organism='Plasmodium falciparum 3D7', project='plasmodb' returns the ~5,720 P. falciparum 3D7 genes.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'VEuPathDB_search_genes_by_organism',
     'arguments': {'organism': '...'}},
)
```

## `WikiPathways_find_pathways_by_gene` {#WikiPathways_find_pathways_by_gene}

Finds all WikiPathways pathways containing a gene given as HGNC symbol, Entrez or Ensembl ID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene`
- **Keywords:** gene to pathway, pathway membership, entrez, ensembl
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string | yes | Gene identifier to search. Examples: 'TP53', 'BRCA1', 'EGFR', '7157' (Entrez), 'ENSG00000141510' (Ensembl). |
| `species` | string | no | Species to filter results. Default: 'Homo sapiens'. Other options: 'Mus musculus', 'Rattus norvegicus', 'Danio rerio'. Default `Homo sapiens` |

??? note "Upstream description"

    Find all WikiPathways pathways containing a specific gene. Takes a gene identifier (HGNC symbol, Entrez ID, or Ensembl ID) and returns all pathways in which that gene appears. Useful for understanding the biological context of a gene across curated pathway databases. Example: searching for 'TP53' returns multiple pathways including WP254 (Apoptosis), WP179 (Cell Cycle), and others. Supports multiple species.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'WikiPathways_find_pathways_by_gene',
     'arguments': {'gene': '...'}},
)
```

## `WikiPathways_get_pathway_genes` {#WikiPathways_get_pathway_genes}

Gets the genes (gene products) participating in a WikiPathways pathway.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** pathway genes, gene set, wpid
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `pathway_id` | string \| null | no | WikiPathways pathway identifier. Examples: 'WP254' (Apoptosis, 87 gene products), 'WP179' (Cell cycle), 'WP4216' (Chromosomal and microsatellite instability), 'WP1584' (Type II diabetes mellitus). Required unless `wpid` is given instead. |
| `wpid` | string \| null | no | Alias for `pathway_id`, matching the parameter name used by the sibling WikiPathways_get_pathway tool. Provide either `pathway_id` or `wpid`. |
| `code` | string \| null | no | Optional filter restricting results to gene products annotated from one identifier system: 'H' (HGNC), 'En' (Ensembl), 'S' (UniProt), 'L' (Entrez Gene), 'Ce' (ChEBI, metabolites). Omit to return every gene product in the pathway. Most pathways are annotated from Entrez Gene or Ensembl and carry no HGNC-sourced entries, so filtering can legitimately return 0. |

??? note "Upstream description"

    Get the genes involved in a WikiPathways pathway. `gene_count` counts distinct gene products; `genes` lists every symbol label WikiPathways records for them, which is larger because one gene product carries several aliases (AKT, AKT1, Akt1 are one node) and occasionally an unrelated symbol -- use `gene_products` for node-level identifiers and their exact labels. By default every gene product is returned regardless of which database annotated it; set `code` only when you need genes cross-referenced to one identifier system, and note most pathways are annotated from Entrez Gene or Ensembl rather than HGNC. Useful for extracting gene sets from curated pathways for enrichment analysis or network building. Example: WP254 (Apoptosis) has 87 gene products (139 labels) including AKT1, APAF1, BAD, BAK1, BAX, BCL2, BCL2L1, BID, CASP3, CASP8, CASP9, CFLAR, CYCS, DIABLO, FADD, FAS, NFKB1, TP53, TRAF1, XIAP.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'WikiPathways_get_pathway_genes', 'arguments': {}},
)
```

## `WormBase_get_gene` {#WormBase_get_gene}

Gets a C. elegans gene overview from WormBase by gene ID, including names, description and type.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** c. elegans, nematode, wbgene, gene description
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | WormBase gene identifier (WBGene format). Examples: 'WBGene00006763' (unc-26), 'WBGene00006789' (unc-52), 'WBGene00000898' (daf-2). |

??? note "Upstream description"

    Get detailed gene overview from WormBase for C. elegans (nematode worm) by WormBase gene ID. Returns gene name, sequence name, species, description, gene type, and status. C. elegans is one of the most important model organisms for genetics, neuroscience, and developmental biology. Example: WBGene00006763 returns unc-26 (synaptojanin, a polyphosphoinositide phosphatase required for synaptic vesicle recycling).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'WormBase_get_gene', 'arguments': {'gene_id': '...'}},
)
```

## `WormBase_get_interactions` {#WormBase_get_interactions}

Gets physical and genetic interactions for a C. elegans gene from WormBase.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** c. elegans, genetic interaction, physical interaction, interactors
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | WBGene identifier or gene name. Examples: 'WBGene00006818', 'unc-86'. |

??? note "Upstream description"

    Get physical and genetic interactions for a C. elegans gene from WormBase, split into physical and genetic lists with interactor names, IDs, interaction type, and citation. Accepts a WBGene ID or gene name. Example: WBGene00006818 (unc-86) has a physical interactor die-1 (WBGene00000995) and a genetic interactor ttx-3 (WBGene00006654).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'WormBase_get_interactions',
     'arguments': {'gene_id': '...'}},
)
```

## `WormBase_get_orthologs` {#WormBase_get_orthologs}

Gets orthologs, homologs and paralogs of a C. elegans gene across species from WormBase.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** c. elegans, orthologs, paralogs, homology, inparanoid
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | WBGene identifier or gene name. Examples: 'WBGene00006818', 'unc-86'. |

??? note "Upstream description"

    Get orthologs and homologs of a C. elegans gene across species from WormBase (Inparanoid/EnsEMBL-Compara/OrthoFinder etc.), plus nematode orthologs and within-worm paralogs. Accepts a WBGene ID or gene name (auto-resolved). Example: WBGene00006818 (unc-86) returns cross-species orthologs like zebrafish pou4f4 and nematode orthologs like CSP40.g8217 (c_tribulationis).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'WormBase_get_orthologs', 'arguments': {'gene_id': '...'}},
)
```

## `Xenbase_get_gene` {#Xenbase_get_gene}

Gets Xenopus gene details from Xenbase via the Alliance of Genome Resources, including synopsis and genomic location.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** xenopus, frog, alliance of genome resources, genomic location
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | Xenbase gene ID. Can be with or without 'Xenbase:' prefix. Examples: 'XB-GENE-484286' (tp53, X. tropicalis), 'XB-GENE-6252168' (tp53.L, X. laevis), 'XB-GENE-487192' (pax6, X. tropicalis). |

??? note "Upstream description"

    Get detailed Xenopus gene information from Xenbase via the Alliance of Genome Resources. Returns gene symbol, name, synopsis, genomic location, synonyms, cross-references, and species details. Supports both Xenopus tropicalis and Xenopus laevis genes. Gene IDs use Xenbase format (e.g., 'XB-GENE-484286' or 'Xenbase:XB-GENE-484286').

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Xenbase_get_gene', 'arguments': {'gene_id': '...'}},
)
```

## `Xenbase_search_genes` {#Xenbase_search_genes}

Searches Xenopus tropicalis and laevis genes in Xenbase via the Alliance API.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** xenopus, frog, gene search, alliance of genome resources
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Gene symbol or name to search. Examples: 'tp53', 'pax6', 'sox2', 'bmp4', 'wnt3a', 'fgf8'. |
| `species` | string | no | Optional species filter: 'tropicalis' for X. tropicalis, 'laevis' for X. laevis. Leave empty for both species. Default `` |
| `limit` | integer | no | Maximum number of results to return. Default: 10. Default `10` |

??? note "Upstream description"

    Search for Xenopus (frog) genes in Xenbase via the Alliance of Genome Resources API. Returns gene IDs, symbols, names, and species for both Xenopus tropicalis and Xenopus laevis. Xenbase is the primary database for Xenopus genetics and developmental biology. Optionally filter by species (tropicalis or laevis). Useful for finding amphibian homologs of human disease genes.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Xenbase_search_genes', 'arguments': {'query': '...'}},
)
```

## `ZFIN_get_gene` {#ZFIN_get_gene}

Gets zebrafish gene details from ZFIN via the Alliance API, including synopsis, genomic location and cross-references.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** zebrafish, danio rerio, alliance of genome resources, genomic location
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | ZFIN gene ID with 'ZFIN:' prefix. Examples: 'ZFIN:ZDB-GENE-990415-8' (pax2a), 'ZFIN:ZDB-GENE-980526-166' (shha/sonic hedgehog a), 'ZFIN:ZDB-GENE-990415-72' (tp53/tumor protein p53), 'ZFIN:ZDB-GENE-030131-7740' (bmp4/bone morphogenetic protein 4). |

??? note "Upstream description"

    Get detailed Danio rerio (zebrafish) gene information from ZFIN via the Alliance of Genome Resources. Returns gene symbol, name, synopsis, genomic location, synonyms, and cross-references. ZFIN is the primary database for zebrafish genetics. Example: 'ZFIN:ZDB-GENE-990415-8' returns pax2a (paired box 2a) with its role in nervous system development.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ZFIN_get_gene', 'arguments': {'gene_id': '...'}},
)
```

## `ZFIN_get_gene_orthologs` {#ZFIN_get_gene_orthologs}

Gets orthologs of a zebrafish gene in human, mouse, fly, worm, yeast and other species.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** zebrafish, orthologs, human ortholog, alliance
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | ZFIN gene ID with 'ZFIN:' prefix. Examples: 'ZFIN:ZDB-GENE-990415-8' (pax2a), 'ZFIN:ZDB-GENE-990415-72' (tp53), 'ZFIN:ZDB-GENE-980526-166' (shha). |
| `limit` | integer | no | Maximum number of orthologs to return (1-100). Default: 20. |
| `page` | integer | no | Page number for pagination. Default: 1. |
| `stringency` | string | no | Stringency filter for orthology calls: 'stringent' (high confidence) or 'moderate'. Default: 'stringent'. |

??? note "Upstream description"

    Get orthologs of a zebrafish gene across other species from the Alliance of Genome Resources. Returns orthologous genes in human, mouse, fly, worm, yeast, and other model organisms. Essential for translating zebrafish findings to human gene function. Zebrafish share approximately 70% of genes with humans. Example: 'ZFIN:ZDB-GENE-990415-8' (pax2a) returns human PAX2 as ortholog.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ZFIN_get_gene_orthologs',
     'arguments': {'gene_id': '...'}},
)
```

## `ZFIN_search` {#ZFIN_search}

Searches ZFIN zebrafish genes by symbol or name and returns gene IDs.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** zebrafish, gene search, zfin id
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Gene symbol or name to search for (e.g., 'tp53', 'sox9a', 'pax2a', 'fgf8a') |
| `limit` | integer | no | Maximum number of results (default: 10) Default `10` |

??? note "Upstream description"

    Search for zebrafish (Danio rerio) genes in ZFIN by keyword (gene symbol or name). Returns matching genes with their ZFIN IDs, symbols, and names. Use the returned gene_id with ZFIN_get_gene for full details. Example: query='tp53' returns zebrafish tumor protein p53.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ZFIN_search', 'arguments': {'query': '...'}},
)
```
