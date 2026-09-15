"""Generated DNAKit database APIs; implementations are loaded only when called."""
from __future__ import annotations
from builtins import object as _Object

from ._runtime import call as _call


def alpha_genome_predict_interval(*,
    chromosome: str,
    start: int,
    end: int,
    output_types: list[_Object] | None | None = None,
    ontology_terms: list[_Object] | None | None = None,
    organism: str | None | None = None,
    sequence_length: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    'Predict multimodal genomic tracks for a genomic interval with DeepMind AlphaGenome (single DNA-sequence model; up to 1 Mb at single-base resolution). Returns a compact per-modality summary (track counts and shapes) for the requested outputs. Hosted API via the alphagenome SDK; free for non-commercial use (set ALPHA_GENOME_API_KEY).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['chromosome'] = chromosome
    _dnakit_arguments['start'] = start
    _dnakit_arguments['end'] = end
    if output_types is not None:
        _dnakit_arguments['output_types'] = output_types
    if ontology_terms is not None:
        _dnakit_arguments['ontology_terms'] = ontology_terms
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if sequence_length is not None:
        _dnakit_arguments['sequence_length'] = sequence_length
    return _call('alpha_genome_predict_interval', _dnakit_arguments, verbose=verbose)

def alpha_genome_score_variant(*,
    chromosome: str,
    position: int,
    reference_bases: str,
    alternate_bases: str,
    output_type: str | None | None = None,
    organism: str | None | None = None,
    sequence_length: str | None | None = None,
    top_n: int | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Score a regulatory variant's effect with DeepMind AlphaGenome (Avsec, Nature 2026), the hosted successor to Enformer/Borzoi. Predicts ref vs alt over up to 1 Mb at single-base resolution and returns the recommended per-track effect scores (sorted by |effect|). Hosted API via the alphagenome SDK; free for non-commercial use (set ALPHA_GENOME_API_KEY).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['chromosome'] = chromosome
    _dnakit_arguments['position'] = position
    _dnakit_arguments['reference_bases'] = reference_bases
    _dnakit_arguments['alternate_bases'] = alternate_bases
    if output_type is not None:
        _dnakit_arguments['output_type'] = output_type
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if sequence_length is not None:
        _dnakit_arguments['sequence_length'] = sequence_length
    if top_n is not None:
        _dnakit_arguments['top_n'] = top_n
    return _call('alpha_genome_score_variant', _dnakit_arguments, verbose=verbose)

def alpha_missense_get_protein_scores(*,
    uniprot_id: str,
    max_residues: int | None | None = None,
    verbose: bool = True,
) -> _Object:
    'Get AlphaMissense pathogenicity scores for EVERY residue of a protein in one call. The upstream hegelab API only exposes per-residue queries; this tool hides that by looking up the protein length from UniProt and concurrently fetching all residues via a thread pool. Returns a list of per-position records (one per residue) plus protein length, count of successful fetches, and the pathogenic/ambiguous/benign thresholds. Use this when you want whole-protein vulnerability assessment, hotspot detection, or to build a (20 × n_positions) AlphaMissense score matrix for benchmarking. Set max_residues > 0 to cap the loop for fast diagnostic calls.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['uniprot_id'] = uniprot_id
    if max_residues is not None:
        _dnakit_arguments['max_residues'] = max_residues
    return _call('alpha_missense_get_protein_scores', _dnakit_arguments, verbose=verbose)

def alpha_missense_get_residue_scores(*,
    uniprot_id: str,
    position: int,
    verbose: bool = True,
) -> _Object:
    'Get the AlphaMissense class breakdown for substitutions at a specific protein position: which residues fall in each class (benign/ambiguous/pathogenic), reported separately for the SNV-accessible substitutions and for all 19, plus the residue mean scores. NOTE: the upstream hotspot API does not expose a numeric score per individual substitution, so this returns class membership rather than 20 scores.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['uniprot_id'] = uniprot_id
    _dnakit_arguments['position'] = position
    return _call('alpha_missense_get_residue_scores', _dnakit_arguments, verbose=verbose)

def alpha_missense_get_variant_score(*,
    uniprot_id: str,
    variant: str,
    verbose: bool = True,
) -> _Object:
    "Classify a specific missense variant with AlphaMissense. Input: UniProt ID and variant (e.g., 'p.R123H'). Returns the variant's AlphaMissense class (pathogenic/ambiguous/benign) and whether it came from the SNV-accessible or the all-substitutions set. NOTE: the upstream hotspot API publishes per-residue class membership, not a numeric score per substitution, so `pathogenicity_score` is null and `score_available` is false; the `residue_mean_score_*` fields are residue-level averages, NOT this variant's score. For a numeric per-variant score use GenomeNexus_annotate_variant.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['uniprot_id'] = uniprot_id
    _dnakit_arguments['variant'] = variant
    return _call('alpha_missense_get_variant_score', _dnakit_arguments, verbose=verbose)

def ar_xiv_search_papers(*,
    query: str,
    limit: int | None = None,
    sort_by: str | None = None,
    sort_order: str | None = None,
    date_from: str | None = None,
    date_to: str | None = None,
    category: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Search arXiv for papers by keyword using the public arXiv API. Returns papers with title, abstract, authors, publication date, category, and URL.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if sort_by is not None:
        _dnakit_arguments['sort_by'] = sort_by
    if sort_order is not None:
        _dnakit_arguments['sort_order'] = sort_order
    if date_from is not None:
        _dnakit_arguments['date_from'] = date_from
    if date_to is not None:
        _dnakit_arguments['date_to'] = date_to
    if category is not None:
        _dnakit_arguments['category'] = category
    return _call('ar_xiv_search_papers', _dnakit_arguments, verbose=verbose)

def blast_protein_search(*,
    sequence: str,
    database: str | None = None,
    expect: float | None = None,
    hitlist_size: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Search protein sequences using NCBI BLAST blastp against protein databases. Requires Biopython (pip install biopython). IMPORTANT: This tool uses NCBI's remote BLAST service which is inherently SLOW - expect 5-30 minutes per query depending on server load. This is normal NCBI behavior, not a bug. For faster results, use 'swissprot' or 'pdb' databases instead of 'nr'.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['sequence'] = sequence
    if database is not None:
        _dnakit_arguments['database'] = database
    if expect is not None:
        _dnakit_arguments['expect'] = expect
    if hitlist_size is not None:
        _dnakit_arguments['hitlist_size'] = hitlist_size
    return _call('blast_protein_search', _dnakit_arguments, verbose=verbose)

def bvbrc_search_epitopes(*,
    taxon_id: str | None | None = None,
    protein_name: str | None | None = None,
    epitope_type: str | None | None = None,
    organism: str | None | None = None,
    limit: int | None | None = None,
    verbose: bool = True,
) -> _Object:
    'Search for pathogen epitopes (B-cell and T-cell) in BV-BRC. Returns epitope sequences, types, protein associations, and assay results for vaccine and immunology research. Covers epitopes from IEDB for bacteria and viruses including SARS-CoV-2, influenza, tuberculosis, and other pathogens. Example: search taxon_id 2697049 for SARS-CoV-2 spike glycoprotein epitopes.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if taxon_id is not None:
        _dnakit_arguments['taxon_id'] = taxon_id
    if protein_name is not None:
        _dnakit_arguments['protein_name'] = protein_name
    if epitope_type is not None:
        _dnakit_arguments['epitope_type'] = epitope_type
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('bvbrc_search_epitopes', _dnakit_arguments, verbose=verbose)

def bio_bank_japan_phewas_by_variant(*,
    rsid: str | None = None,
    variant: str | None = None,
    limit: int | None = None,
    max_pval: float | None | None = None,
    verbose: bool = True,
) -> _Object:
    'Phenome-wide association study (PheWAS) lookup for a single variant in BioBank Japan (BBJ), the largest Japanese-ancestry biobank. Returns every phenotype associated with the variant (the inverse of a GWAS), sorted by p-value. Accepts an rsID (auto-resolved to GRCh37 via Ensembl) or an explicit chr:pos:ref:alt variant in GRCh37. Use for cross-ancestry replication of a genetic association in an East Asian (Japanese) population.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if rsid is not None:
        _dnakit_arguments['rsid'] = rsid
    if variant is not None:
        _dnakit_arguments['variant'] = variant
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if max_pval is not None:
        _dnakit_arguments['max_pval'] = max_pval
    return _call('bio_bank_japan_phewas_by_variant', _dnakit_arguments, verbose=verbose)

def bio_grid_get_interactions(*,
    gene_names: list[str] | None = None,
    gene_name: str | None = None,
    organism: str | None = None,
    interaction_type: str | None = None,
    evidence_types: list[str] | None = None,
    limit: int | None = None,
    throughput: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Query experimentally validated protein and genetic interactions from BioGRID (Biological General Repository for Interaction Datasets). Returns curated interactions from published studies with evidence methods, PubMed citations, and throughput information. **Prerequisites**: Requires BIOGRID_API_KEY environment variable (free academic API key at https://webservice.thebiogrid.org/). BioGRID contains 2.3M+ interactions from 80+ organisms, all experimentally validated (no predictions). Use for: finding experimentally proven interactions, getting literature evidence for interactions, validating predicted interactions, finding interaction methods used, accessing high-confidence curated data.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if gene_names is not None:
        _dnakit_arguments['gene_names'] = gene_names
    if gene_name is not None:
        _dnakit_arguments['gene_name'] = gene_name
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if interaction_type is not None:
        _dnakit_arguments['interaction_type'] = interaction_type
    if evidence_types is not None:
        _dnakit_arguments['evidence_types'] = evidence_types
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if throughput is not None:
        _dnakit_arguments['throughput'] = throughput
    return _call('bio_grid_get_interactions', _dnakit_arguments, verbose=verbose)

def bio_rxiv_get_preprint(*,
    doi: str,
    server: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Get full metadata for a specific bioRxiv or medRxiv preprint by its DOI. Returns comprehensive metadata including title, authors, abstract, publication date, version, category, corresponding author, funding information, and URLs (HTML, PDF, XML). Response time: <1 second. Use when you have a bioRxiv DOI (format: 10.1101/YYYY.MM.DD.XXXXXX). To SEARCH for preprints by keywords, use EuropePMC_search_articles with 'SRC:PPR' filter (fast, relevance-ranked), or web_search with 'site:biorxiv.org'. bioRxiv API does not support text search.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['doi'] = doi
    if server is not None:
        _dnakit_arguments['server'] = server
    return _call('bio_rxiv_get_preprint', _dnakit_arguments, verbose=verbose)

def bio_rxiv_list_recent_preprints(*,
    start_date: str | None = None,
    end_date: str | None = None,
    server: str | None = None,
    cursor: int | None = None,
    query: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    "List recent bioRxiv or medRxiv preprints by date range. Returns up to 100 preprints posted between the given start and end dates, with title, authors, abstract, DOI, category, and publication status. The bioRxiv API returns results in batches of up to 100; use the cursor parameter for pagination (0 for first page, 100 for second, etc.). Use this to monitor new preprints in a date range or discover recent work in a field. For searching by DOI, use BioRxiv_get_preprint instead. For searching by keywords, use EuropePMC_search_articles with 'SRC:PPR' filter.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if start_date is not None:
        _dnakit_arguments['start_date'] = start_date
    if end_date is not None:
        _dnakit_arguments['end_date'] = end_date
    if server is not None:
        _dnakit_arguments['server'] = server
    if cursor is not None:
        _dnakit_arguments['cursor'] = cursor
    if query is not None:
        _dnakit_arguments['query'] = query
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('bio_rxiv_list_recent_preprints', _dnakit_arguments, verbose=verbose)

def cadd_get_position_scores(*,
    chrom: str,
    pos: int,
    version: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get CADD scores for all possible substitutions at a genomic position. Returns PHRED scores for A, C, G, T variants at the position. Useful for saturation mutagenesis analysis.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['chrom'] = chrom
    _dnakit_arguments['pos'] = pos
    if version is not None:
        _dnakit_arguments['version'] = version
    return _call('cadd_get_position_scores', _dnakit_arguments, verbose=verbose)

def cadd_get_range_scores(*,
    chrom: str,
    start: int,
    end: int,
    version: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get CADD scores for all variants in a genomic range (max 100bp). Returns all pre-computed deleteriousness scores in the region. Use for hotspot analysis or exon scanning.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['chrom'] = chrom
    _dnakit_arguments['start'] = start
    _dnakit_arguments['end'] = end
    if version is not None:
        _dnakit_arguments['version'] = version
    return _call('cadd_get_range_scores', _dnakit_arguments, verbose=verbose)

def cadd_get_variant_score(*,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    version: str | None = None,
    include_annotations: bool | None = None,
    verbose: bool = True,
) -> _Object:
    'Get CADD deleteriousness score for a specific variant. PHRED scores: >=20 top 1% deleterious (pathogenic threshold ~15-20), >=30 top 0.1%. Input: chromosome, position, ref, alt alleles. Supports GRCh37 and GRCh38.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['chrom'] = chrom
    _dnakit_arguments['pos'] = pos
    _dnakit_arguments['ref'] = ref
    _dnakit_arguments['alt'] = alt
    if version is not None:
        _dnakit_arguments['version'] = version
    if include_annotations is not None:
        _dnakit_arguments['include_annotations'] = include_annotations
    return _call('cadd_get_variant_score', _dnakit_arguments, verbose=verbose)

def cellx_gene_download_h5ad(*,
    operation: str | None = None,
    dataset_id: str,
    output_path: str | None = None,
    census_version: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Download original H5AD (HDF5-based AnnData) files from CELLxGENE datasets or get their URIs. Access unprocessed data from individual studies with complete metadata and analysis results. Requires dataset_id from CELLxGENE Discover. Note: Files can be large (GBs). Use for: accessing raw data, reproducibility, offline analysis, custom processing pipelines. Note: requires the cellxgene-census Python package; install with `pip install cellxgene-census` for live access.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['dataset_id'] = dataset_id
    if output_path is not None:
        _dnakit_arguments['output_path'] = output_path
    if census_version is not None:
        _dnakit_arguments['census_version'] = census_version
    return _call('cellx_gene_download_h5ad', _dnakit_arguments, verbose=verbose)

def cellx_gene_get_cell_metadata(*,
    operation: str | None = None,
    organism: str | None = None,
    obs_value_filter: str | None = None,
    column_names: list[str] | None = None,
    census_version: str | None = None,
    gene: str | None = None,
    tissue: str | None = None,
    max_cells: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Query cell metadata from CELLxGENE Census (50M+ human/mouse single cells). CRITICAL: obs_value_filter is REQUIRED - unfiltered queries will timeout. Common filters: tissue_general (brain, heart, lung, liver, blood, kidney, skin, eye, intestine), cell_type (T cell, B cell, neuron, macrophage, fibroblast), disease (normal, COVID-19, Alzheimer disease, cancer). Start with: obs_value_filter=\'tissue_general == "lung" and disease == "normal"\'. Returns cell type, tissue, disease, donor info. Use for: finding cells by criteria, cohort selection, exploring cell populations. Note: requires the cellxgene-census Python package; install with `pip install cellxgene-census` for live access.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if obs_value_filter is not None:
        _dnakit_arguments['obs_value_filter'] = obs_value_filter
    if column_names is not None:
        _dnakit_arguments['column_names'] = column_names
    if census_version is not None:
        _dnakit_arguments['census_version'] = census_version
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    if tissue is not None:
        _dnakit_arguments['tissue'] = tissue
    if max_cells is not None:
        _dnakit_arguments['max_cells'] = max_cells
    return _call('cellx_gene_get_cell_metadata', _dnakit_arguments, verbose=verbose)

def cellx_gene_get_embeddings(*,
    operation: str | None = None,
    organism: str | None = None,
    embedding_name: str | None = None,
    census_version: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Access pre-calculated cell embeddings (scVI, Geneformer) from CELLxGENE Census. Returns available embedding names or specific embedding data. Embeddings enable rapid visualization and analysis without recomputation. Use for: UMAP/t-SNE visualization, cell similarity analysis, transfer learning. Note: requires the cellxgene-census Python package; install with `pip install cellxgene-census` for live access.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if embedding_name is not None:
        _dnakit_arguments['embedding_name'] = embedding_name
    if census_version is not None:
        _dnakit_arguments['census_version'] = census_version
    return _call('cellx_gene_get_embeddings', _dnakit_arguments, verbose=verbose)

def cellx_gene_get_expression_data(*,
    operation: str | None = None,
    organism: str | None = None,
    obs_value_filter: str | None = None,
    var_value_filter: str | None = None,
    obs_column_names: list[str] | None = None,
    var_column_names: list[str] | None = None,
    census_version: str | None = None,
    gene: str | None = None,
    tissue: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Query gene expression data from CELLxGENE Census (50M+ cells, 60K+ genes). CRITICAL: At least one filter REQUIRED - unfiltered queries timeout. Start with: obs_value_filter=\'tissue_general == "lung" and disease == "normal"\' AND var_value_filter=\'feature_name in ["TP53", "EGFR"]\'. Returns AnnData summary (dimensions, metadata). For full matrix data, use Python API with sufficient memory. Note: requires the cellxgene-census Python package; install with `pip install cellxgene-census` for live access.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if obs_value_filter is not None:
        _dnakit_arguments['obs_value_filter'] = obs_value_filter
    if var_value_filter is not None:
        _dnakit_arguments['var_value_filter'] = var_value_filter
    if obs_column_names is not None:
        _dnakit_arguments['obs_column_names'] = obs_column_names
    if var_column_names is not None:
        _dnakit_arguments['var_column_names'] = var_column_names
    if census_version is not None:
        _dnakit_arguments['census_version'] = census_version
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    if tissue is not None:
        _dnakit_arguments['tissue'] = tissue
    return _call('cellx_gene_get_expression_data', _dnakit_arguments, verbose=verbose)

def cosmic_get_mutations_by_gene(*,
    operation: str | None = None,
    gene: str | None = None,
    gene_name: str | None = None,
    max_results: int | None = None,
    genome_build: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get all somatic mutations for a specific gene from COSMIC database. Returns comprehensive mutation data for a gene symbol including CDS changes, amino acid changes, cancer types, and genomic positions. Useful for exploring the mutational landscape of cancer genes like BRAF, TP53, EGFR, KRAS, etc.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    if gene_name is not None:
        _dnakit_arguments['gene_name'] = gene_name
    if max_results is not None:
        _dnakit_arguments['max_results'] = max_results
    if genome_build is not None:
        _dnakit_arguments['genome_build'] = genome_build
    return _call('cosmic_get_mutations_by_gene', _dnakit_arguments, verbose=verbose)

def cosmic_search_mutations(*,
    operation: str | None = None,
    terms: str | None = None,
    query: str | None = None,
    max_results: int | None = None,
    genome_build: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Search COSMIC database for somatic mutations in cancer. COSMIC is the world's largest expert-curated database of somatic mutations in human cancer. Search by gene name, mutation ID, or mutation description. Returns mutation IDs, gene names, CDS changes, and amino acid changes. No authentication required.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if terms is not None:
        _dnakit_arguments['terms'] = terms
    if query is not None:
        _dnakit_arguments['query'] = query
    if max_results is not None:
        _dnakit_arguments['max_results'] = max_results
    if genome_build is not None:
        _dnakit_arguments['genome_build'] = genome_build
    return _call('cosmic_search_mutations', _dnakit_arguments, verbose=verbose)

def ctd_get_chemical_gene_interactions(*,
    input_terms: str,
    verbose: bool = True,
) -> _Object:
    "PERMANENTLY UNAVAILABLE: this tool always returns an error. CTD's native batchQuery.go is CAPTCHA-blocked for programmatic clients, and the RENCI Automat mirror this tool used to fall back on has been fully decommissioned (its registry no longer lists a 'ctd' backend at all). No other free, live source of CTD's curated chemical-gene interactions was found. Use DGIdb_get_drug_gene_interactions or ChEMBL/OpenTargets tools for target-compound bioactivity data instead.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['input_terms'] = input_terms
    return _call('ctd_get_chemical_gene_interactions', _dnakit_arguments, verbose=verbose)

def ch_embl_get_drug_mechanisms(*,
    limit: int | None = None,
    offset: int | None = None,
    drug_chembl_id: str | None = None,
    drug_name: str | None = None,
    molecule_chembl_id: str | None = None,
    chembl_id: str | None = None,
    drug_chembl_id__exact: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Get mechanisms of action for a drug by ChEMBL drug ID or drug name. Accepts drug_chembl_id (e.g., 'CHEMBL1201581') or drug_name (e.g., 'trastuzumab') — drug_name triggers automatic ChEMBL ID lookup.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    if drug_chembl_id is not None:
        _dnakit_arguments['drug_chembl_id'] = drug_chembl_id
    if drug_name is not None:
        _dnakit_arguments['drug_name'] = drug_name
    if molecule_chembl_id is not None:
        _dnakit_arguments['molecule_chembl_id'] = molecule_chembl_id
    if chembl_id is not None:
        _dnakit_arguments['chembl_id'] = chembl_id
    if drug_chembl_id__exact is not None:
        _dnakit_arguments['drug_chembl_id__exact'] = drug_chembl_id__exact
    return _call('ch_embl_get_drug_mechanisms', _dnakit_arguments, verbose=verbose)

def ch_embl_get_molecule(*,
    chembl_id: str | None = None,
    format: str | None = None,
    molecule_chembl_id: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get detailed information about a molecule by its ChEMBL ID. Returns molecule properties, structures, synonyms, and related data.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if chembl_id is not None:
        _dnakit_arguments['chembl_id'] = chembl_id
    if format is not None:
        _dnakit_arguments['format'] = format
    if molecule_chembl_id is not None:
        _dnakit_arguments['molecule_chembl_id'] = molecule_chembl_id
    return _call('ch_embl_get_molecule', _dnakit_arguments, verbose=verbose)

def ch_embl_get_target_activities(*,
    target_chembl_id__exact: str | None = None,
    limit: int | None = None,
    offset: int | None = None,
    target_chembl_id: str | None = None,
    chembl_id: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Get all activity data for a target by ChEMBL target ID. Returns bioactivity measurements (IC50, Ki, etc.) for molecules tested against this target. Note: activity values may be in different units (nM or uM) depending on the assay — check the 'standard_units' field in each result. To find a target ID, use ChEMBL_search_targets with a target name or gene symbol.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if target_chembl_id__exact is not None:
        _dnakit_arguments['target_chembl_id__exact'] = target_chembl_id__exact
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    if target_chembl_id is not None:
        _dnakit_arguments['target_chembl_id'] = target_chembl_id
    if chembl_id is not None:
        _dnakit_arguments['chembl_id'] = chembl_id
    return _call('ch_embl_get_target_activities', _dnakit_arguments, verbose=verbose)

def ch_embl_search_drugs(*,
    query: str | None = None,
    max_phase: int | None = None,
    limit: int | None = None,
    offset: int | None = None,
    molecule_chembl_id: str | None = None,
    q: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Search drugs by name, approval status, or other criteria. Use the `query` parameter to search by drug name (e.g., {"query": "sotorasib"}). Use `max_phase` to filter by development phase.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if query is not None:
        _dnakit_arguments['query'] = query
    if max_phase is not None:
        _dnakit_arguments['max_phase'] = max_phase
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    if molecule_chembl_id is not None:
        _dnakit_arguments['molecule_chembl_id'] = molecule_chembl_id
    if q is not None:
        _dnakit_arguments['q'] = q
    return _call('ch_embl_search_drugs', _dnakit_arguments, verbose=verbose)

def ch_embl_search_targets(*,
    target_chembl_id: str | None = None,
    pref_name__contains: str | None = None,
    organism: str | None = None,
    target_type: str | None = None,
    fields: list[str] | None = None,
    limit: int | None = None,
    offset: int | None = None,
    q: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Search targets by name, organism, target type, or other criteria. Use this tool to find target ChEMBL IDs needed for other tools.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if target_chembl_id is not None:
        _dnakit_arguments['target_chembl_id'] = target_chembl_id
    if pref_name__contains is not None:
        _dnakit_arguments['pref_name__contains'] = pref_name__contains
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if target_type is not None:
        _dnakit_arguments['target_type'] = target_type
    if fields is not None:
        _dnakit_arguments['fields'] = fields
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    if q is not None:
        _dnakit_arguments['q'] = q
    return _call('ch_embl_search_targets', _dnakit_arguments, verbose=verbose)

def ch_ipatlas_enrichment_analysis(*,
    operation: str | None = None,
    bed_data: str | None = None,
    motif: str | None = None,
    gene_list: list[_Object] | str | None = None,
    genome: str | None = None,
    antigen_class: str | None = None,
    cell_type_class: str | None = None,
    threshold: str | None = None,
    distance: str | None = None,
    gene: list[_Object] | str | None = None,
    cell_type: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Perform enrichment analysis to identify transcription factors and histone modifications enriched in your data. **Required: Provide ONE input type** - (1) BED genomic regions, (2) DNA sequence motif (IUPAC notation: A/T/G/C/W/S/M/K/R/Y/B/D/H/V/N), or (3) gene symbol list. Compares your input against 433,000+ ChIP-seq/ATAC-seq/Bisulfite-seq experiments to identify significant enrichment. Returns ranked list of proteins bound to your regions/motif or regulating your genes. Note: Returns submission URL (web form-based analysis). Use for: identifying regulators of genomic regions, finding proteins bound to sequence motifs, discovering transcription factors regulating gene sets.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if bed_data is not None:
        _dnakit_arguments['bed_data'] = bed_data
    if motif is not None:
        _dnakit_arguments['motif'] = motif
    if gene_list is not None:
        _dnakit_arguments['gene_list'] = gene_list
    if genome is not None:
        _dnakit_arguments['genome'] = genome
    if antigen_class is not None:
        _dnakit_arguments['antigen_class'] = antigen_class
    if cell_type_class is not None:
        _dnakit_arguments['cell_type_class'] = cell_type_class
    if threshold is not None:
        _dnakit_arguments['threshold'] = threshold
    if distance is not None:
        _dnakit_arguments['distance'] = distance
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    if cell_type is not None:
        _dnakit_arguments['cell_type'] = cell_type
    return _call('ch_ipatlas_enrichment_analysis', _dnakit_arguments, verbose=verbose)

def ch_ipatlas_get_experiments(*,
    operation: str | None = None,
    genome: str | None = None,
    antigen: str | None = None,
    cell_type: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search ChIP-Atlas experiment metadata including experiment IDs, antigens, cell types, and processing statistics. Filter by genome, antigen (TF/histone), or cell type. Returns experiment details with SRX/ERX/DRX IDs. Use for: finding relevant experiments, getting experiment metadata, dataset discovery, batch analysis planning.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if genome is not None:
        _dnakit_arguments['genome'] = genome
    if antigen is not None:
        _dnakit_arguments['antigen'] = antigen
    if cell_type is not None:
        _dnakit_arguments['cell_type'] = cell_type
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('ch_ipatlas_get_experiments', _dnakit_arguments, verbose=verbose)

def ch_ipatlas_get_peak_data(*,
    operation: str | None = None,
    experiment_id: str | None = None,
    genome: str | None = None,
    format: str | None = None,
    threshold: str | None = None,
    gene: str | None = None,
    experiment_type: str | None = None,
    limit: int | None = None,
    peak_file: str | None = None,
    flank_bp: int | None = None,
    max_file_bytes: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get download URLs for ChIP-Atlas peak-call data in BigWig, BED, or BigBed format. BigWig contains coverage scores in RPM (Reads Per Million). BED/BigBed contain peak regions with MACS2 (Model-based Analysis of ChIP-Seq) peak caller scores. Use for: downloading raw data, integrating with genome browsers (UCSC, IGV), custom analysis, visualization.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if experiment_id is not None:
        _dnakit_arguments['experiment_id'] = experiment_id
    if genome is not None:
        _dnakit_arguments['genome'] = genome
    if format is not None:
        _dnakit_arguments['format'] = format
    if threshold is not None:
        _dnakit_arguments['threshold'] = threshold
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    if experiment_type is not None:
        _dnakit_arguments['experiment_type'] = experiment_type
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if peak_file is not None:
        _dnakit_arguments['peak_file'] = peak_file
    if flank_bp is not None:
        _dnakit_arguments['flank_bp'] = flank_bp
    if max_file_bytes is not None:
        _dnakit_arguments['max_file_bytes'] = max_file_bytes
    return _call('ch_ipatlas_get_peak_data', _dnakit_arguments, verbose=verbose)

def ch_ipatlas_search_datasets(*,
    operation: str | None = None,
    antigen: str | None = None,
    cell_type: str | None = None,
    genome: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Search ChIP-Atlas by antigen or cell type to find available datasets. Returns number of experiments and experiment IDs. Use antigenList for finding TF/histone datasets across cell types. Use celltypeList for finding all factors in specific cell type. Use for: dataset discovery, planning comparative analysis, checking data availability.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if antigen is not None:
        _dnakit_arguments['antigen'] = antigen
    if cell_type is not None:
        _dnakit_arguments['cell_type'] = cell_type
    if genome is not None:
        _dnakit_arguments['genome'] = genome
    return _call('ch_ipatlas_search_datasets', _dnakit_arguments, verbose=verbose)

def clin_gen_dosage_by_gene(*,
    gene: str | None = None,
    gene_symbol: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get ClinGen dosage sensitivity curation for a specific gene using the JSON API. Returns haploinsufficiency (HI) and triplosensitivity (TS) assertion scores with genomic coordinates in both GRCh37 and GRCh38. Scores: 0=No evidence, 1=Little evidence, 2=Emerging evidence, 3=Sufficient evidence, 30=Gene associated with autosomal recessive phenotype, 40=Dosage sensitivity unlikely. Also includes pLI score, HI prediction score, and associated diseases with MONDO IDs.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    if gene_symbol is not None:
        _dnakit_arguments['gene_symbol'] = gene_symbol
    return _call('clin_gen_dosage_by_gene', _dnakit_arguments, verbose=verbose)

def clin_gen_dosage_region_search(*,
    chromosome: str,
    start: int,
    end: int,
    assembly: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Search ClinGen dosage sensitivity curations by genomic region. Returns all genes and recurrent CNV regions with dosage sensitivity curations overlapping the specified chromosomal coordinates. Supports both GRCh37 and GRCh38. Results include haploinsufficiency and triplosensitivity scores. Useful for interpreting the clinical significance of deletions and duplications detected by chromosomal microarray.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['chromosome'] = chromosome
    _dnakit_arguments['start'] = start
    _dnakit_arguments['end'] = end
    if assembly is not None:
        _dnakit_arguments['assembly'] = assembly
    return _call('clin_gen_dosage_region_search', _dnakit_arguments, verbose=verbose)

def clin_gen_get_actionability_adult(*,
    gene: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get ClinGen clinical actionability curations for adult context. Returns genes with actionable findings, intervention scores, and outcome scores. Actionability scores help prioritize which incidental findings to return. Example: BRCA1/2 have high actionability for cancer screening.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    return _call('clin_gen_get_actionability_adult', _dnakit_arguments, verbose=verbose)

def clin_gen_get_actionability_pediatric(*,
    gene: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get ClinGen clinical actionability curations for pediatric context. Returns genes with actionable findings in children, intervention recommendations, and evidence scores. Pediatric curations may differ from adult context due to age-specific interventions.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    return _call('clin_gen_get_actionability_pediatric', _dnakit_arguments, verbose=verbose)

def clin_gen_get_dosage_sensitivity(*,
    gene: str | None = None,
    include_regions: bool | None = None,
    verbose: bool = True,
) -> _Object:
    'Get all ClinGen dosage sensitivity curations. Returns haploinsufficiency and triplosensitivity scores for curated genes. Set include_regions=true to also get recurrent CNV regions. Scores: 0=No evidence, 1=Little, 2=Emerging, 3=Sufficient evidence.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    if include_regions is not None:
        _dnakit_arguments['include_regions'] = include_regions
    return _call('clin_gen_get_dosage_sensitivity', _dnakit_arguments, verbose=verbose)

def clin_gen_get_gene_validity(*,
    gene: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get all ClinGen gene-disease validity curations. Returns comprehensive list of gene-disease relationships with classifications (Definitive, Strong, Moderate, Limited). Useful for batch analysis or exploring curated genes. Optional gene filter available.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    return _call('clin_gen_get_gene_validity', _dnakit_arguments, verbose=verbose)

def clin_gen_get_variant_classifications(*,
    gene: str | None = None,
    variant: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get variant pathogenicity classifications from ClinGen Evidence Repository. Returns expert-curated variant classifications with supporting evidence. Optional filters by gene or variant. Classifications follow ACMG/AMP guidelines with expert review.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    if variant is not None:
        _dnakit_arguments['variant'] = variant
    return _call('clin_gen_get_variant_classifications', _dnakit_arguments, verbose=verbose)

def clin_gen_search_actionability(*,
    gene: str,
    verbose: bool = True,
) -> _Object:
    'Search ClinGen clinical actionability across both adult and pediatric contexts by gene. Returns actionability data from both contexts for comparison. Useful for comprehensive actionability assessment of a gene across age groups.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene'] = gene
    return _call('clin_gen_search_actionability', _dnakit_arguments, verbose=verbose)

def clin_gen_search_dosage_sensitivity(*,
    gene: str,
    verbose: bool = True,
) -> _Object:
    'Search ClinGen dosage sensitivity curations by gene symbol. Returns haploinsufficiency (HI) and triplosensitivity (TS) scores (0-3 scale). Score 3 = Sufficient evidence for dosage pathogenicity. Critical for CNV interpretation. Example: MECP2 has HI score of 3.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene'] = gene
    return _call('clin_gen_search_dosage_sensitivity', _dnakit_arguments, verbose=verbose)

def clin_gen_search_gene_validity(*,
    gene: str,
    verbose: bool = True,
) -> _Object:
    'Search ClinGen gene-disease validity curations by gene symbol. Returns classification (Definitive, Strong, Moderate, Limited, Disputed, Refuted), disease associations, and inheritance patterns. Essential for ACMG variant classification. Example: BRCA1 returns validity scores for hereditary breast/ovarian cancer.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene'] = gene
    return _call('clin_gen_search_gene_validity', _dnakit_arguments, verbose=verbose)

def clin_var_get_clinical_significance(*,
    variant_id: str,
    verbose: bool = True,
) -> _Object:
    'Get clinical significance information for a variant from ClinVar. Returns pathogenicity classification and clinical interpretations.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['variant_id'] = variant_id
    return _call('clin_var_get_clinical_significance', _dnakit_arguments, verbose=verbose)

def clin_var_get_variant_details(*,
    variant_id: str,
    verbose: bool = True,
) -> _Object:
    'Get variant summary information from ClinVar by variant ID. Returns accession, title, genes, clinical significance, review status, chromosome location, and variation name. For the same variant, this provides a different field view than ClinVar_get_clinical_significance.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['variant_id'] = variant_id
    return _call('clin_var_get_variant_details', _dnakit_arguments, verbose=verbose)

def clin_var_search_variants(*,
    gene: str | None = None,
    condition: str | None = None,
    variant_id: str | None = None,
    variant_name: str | list[_Object] | None | None = None,
    max_results: int | None = None,
    limit: int | None = None,
    clinical_significance: str | None = None,
    gene_symbol: str | None | None = None,
    significance: str | None | None = None,
    query: str | None | None = None,
    variant: str | list[_Object] | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Search for variants in ClinVar database by gene name, condition, variant ID, or variant name (protein change / HGVS notation). Returns variant identifiers and basic information. At least one of gene, condition, variant_id, or variant_name must be provided. Combine gene + variant_name for an exact-match lookup of one specific known variant instead of browsing all of a gene's ClinVar entries.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    if condition is not None:
        _dnakit_arguments['condition'] = condition
    if variant_id is not None:
        _dnakit_arguments['variant_id'] = variant_id
    if variant_name is not None:
        _dnakit_arguments['variant_name'] = variant_name
    if max_results is not None:
        _dnakit_arguments['max_results'] = max_results
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if clinical_significance is not None:
        _dnakit_arguments['clinical_significance'] = clinical_significance
    if gene_symbol is not None:
        _dnakit_arguments['gene_symbol'] = gene_symbol
    if significance is not None:
        _dnakit_arguments['significance'] = significance
    if query is not None:
        _dnakit_arguments['query'] = query
    if variant is not None:
        _dnakit_arguments['variant'] = variant
    return _call('clin_var_search_variants', _dnakit_arguments, verbose=verbose)

def clinical_trials_search_studies(*,
    operation: str | None = None,
    query_cond: str | None | None = None,
    query_intr: str | None | None = None,
    query_term: str | None | None = None,
    filter_status: str | None | None = None,
    filter_phase: str | None | None = None,
    filter_study_type: str | None | None = None,
    page_size: int | None = None,
    next_page_token: str | None | None = None,
    query: str | None | None = None,
    condition: str | None | None = None,
    status: str | None | None = None,
    max_results: int | None | None = None,
    intervention: str | None | None = None,
    sponsor: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Search ClinicalTrials.gov for clinical trial studies by condition, intervention, sponsor, or other criteria. Returns NCT IDs, titles, status, phase, enrollment, and key trial metadata. Example: search for 'breast cancer' trials with 'pembrolizumab' that are RECRUITING in Phase 3. Query rewriting: multi-word intervention/condition/sponsor values are sent to ClinicalTrials.gov as an EXACT PHRASE, and intervention/sponsor values are further restricted to specific registry fields, so a phrase no registrant recorded verbatim can match 0 studies even when many trials mention it. The query actually executed is always returned in `executed_query` (with per-parameter detail in `query_rewrites`); when a rewritten query matches nothing, `relaxed_match_check` reports how many studies your original wording would have matched and how to retrieve them. To opt out of the rewrite, wrap a value in parentheses (e.g. '(donor derived cell therapy)') and it is sent verbatim as a loose keyword match.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if query_cond is not None:
        _dnakit_arguments['query_cond'] = query_cond
    if query_intr is not None:
        _dnakit_arguments['query_intr'] = query_intr
    if query_term is not None:
        _dnakit_arguments['query_term'] = query_term
    if filter_status is not None:
        _dnakit_arguments['filter_status'] = filter_status
    if filter_phase is not None:
        _dnakit_arguments['filter_phase'] = filter_phase
    if filter_study_type is not None:
        _dnakit_arguments['filter_study_type'] = filter_study_type
    if page_size is not None:
        _dnakit_arguments['page_size'] = page_size
    if next_page_token is not None:
        _dnakit_arguments['next_page_token'] = next_page_token
    if query is not None:
        _dnakit_arguments['query'] = query
    if condition is not None:
        _dnakit_arguments['condition'] = condition
    if status is not None:
        _dnakit_arguments['status'] = status
    if max_results is not None:
        _dnakit_arguments['max_results'] = max_results
    if intervention is not None:
        _dnakit_arguments['intervention'] = intervention
    if sponsor is not None:
        _dnakit_arguments['sponsor'] = sponsor
    return _call('clinical_trials_search_studies', _dnakit_arguments, verbose=verbose)

def crossref_search_works(*,
    query: str,
    limit: int | None = None,
    offset: int | None = None,
    filter: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Search Crossref Works API for scholarly articles, publications, and research outputs by keyword. Returns metadata including title, abstract, journal, authors, publication year, DOI, URL, citations, and publisher. Supports filtering by publication type (journal article, conference paper, etc.) and date range. Use this to discover academic literature, verify citations, or build bibliographies. Find DOIs to use with Crossref_get_work for complete metadata. The response reports 'count' (works in this page) alongside 'total_count' (works matching the query in all of Crossref) -- these differ on almost every real search, since 'limit' defaults to 10 and caps at 100, so treat 'count' as a page size and never as the size of the literature. When they differ, 'truncated' is true and 'truncation_note' gives the 'offset' for the next page. Note: Crossref's own API weakens keyword relevance-matching whenever a non-default 'sort' query param is combined with 'query' (e.g. sort=is-referenced-by-count), returning globally top-cited works largely unrelated to the search term -- omit 'sort' (the default relevance ranking) unless you have separately verified the combination for your query.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    if filter is not None:
        _dnakit_arguments['filter'] = filter
    return _call('crossref_search_works', _dnakit_arguments, verbose=verbose)

def dgidb_get_drug_gene_interactions(*,
    genes: list[str] | None = None,
    gene_name: str | None = None,
    gene: str | None = None,
    interaction_sources: list[str] | None = None,
    interaction_types: list[str] | None = None,
    verbose: bool = True,
) -> _Object:
    'Get drug-gene interactions for a list of genes from DGIdb. Returns drugs targeting each gene, interaction types, and evidence sources. Gene symbols DGIdb does not recognize are omitted from the results and named in metadata.unmatched_genes, so a missing gene is distinguishable from a gene with no interactions.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if genes is not None:
        _dnakit_arguments['genes'] = genes
    if gene_name is not None:
        _dnakit_arguments['gene_name'] = gene_name
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    if interaction_sources is not None:
        _dnakit_arguments['interaction_sources'] = interaction_sources
    if interaction_types is not None:
        _dnakit_arguments['interaction_types'] = interaction_types
    return _call('dgidb_get_drug_gene_interactions', _dnakit_arguments, verbose=verbose)

def dgidb_get_gene_druggability(*,
    genes: list[str] | None = None,
    gene_name: str | None = None,
    gene: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get druggability information for genes. Returns gene categories indicating if a gene is druggable (kinase, GPCR, ion channel, etc.). Gene symbols DGIdb does not recognize are omitted from the results and named in metadata.unmatched_genes.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if genes is not None:
        _dnakit_arguments['genes'] = genes
    if gene_name is not None:
        _dnakit_arguments['gene_name'] = gene_name
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    return _call('dgidb_get_gene_druggability', _dnakit_arguments, verbose=verbose)

def dna_gibson_design(*,
    operation: str | None = None,
    fragments: list[str],
    overlap_length: int | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Design Gibson Assembly overlaps for seamless DNA fragment assembly. For each fragment, computes the 5' overhang (from the previous fragment's 3' end) and 3' overhang (from the next fragment's 5' end) based on the specified overlap length (default 20 bp). Returns the modified fragment sequences ready for Gibson Assembly cloning. Use for: seamless cloning, multi-fragment assembly, vector backbone insertion.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['fragments'] = fragments
    if overlap_length is not None:
        _dnakit_arguments['overlap_length'] = overlap_length
    return _call('dna_gibson_design', _dnakit_arguments, verbose=verbose)

def dna_golden_gate_assemble(*,
    operation: str | None = None,
    fragments: list[str],
    enzyme: str | None | None = None,
    circular: bool | None | None = None,
    labels: list[_Object] | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Work out the product of a Golden Gate reaction that has already been performed: given the plasmids/parts that were combined and the Type IIS enzyme used (BsaI, BbsI, Esp3I/BsmBI, SapI), digest each input, discard fragments that keep a recognition site (they are re-cut in the reaction), and ligate the rest by matching 4-bp overhangs into the circular product. Returns the assembled sequence, its length, and the fragment order with the overhang at each junction. Use to answer 'I combined these plasmids with <enzyme> — what does the resulting plasmid contain/express?'. This is the reverse of DNA_golden_gate_design, which designs parts rather than analysing a finished reaction.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['fragments'] = fragments
    if enzyme is not None:
        _dnakit_arguments['enzyme'] = enzyme
    if circular is not None:
        _dnakit_arguments['circular'] = circular
    if labels is not None:
        _dnakit_arguments['labels'] = labels
    return _call('dna_golden_gate_assemble', _dnakit_arguments, verbose=verbose)

def dna_golden_gate_design(*,
    operation: str | None = None,
    parts: list[str],
    enzyme: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    'Design Golden Gate Assembly parts with Type IIS restriction enzyme sites (BsaI or BbsI). Assigns unique non-palindromic 4-bp overhangs to each junction and builds the full synthetic part sequences including enzyme recognition sites. Use for: scarless multi-part assembly, modular cloning (MoClo), Golden Gate cloning workflows.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['parts'] = parts
    if enzyme is not None:
        _dnakit_arguments['enzyme'] = enzyme
    return _call('dna_golden_gate_design', _dnakit_arguments, verbose=verbose)

def dna_translate_reading_frames(*,
    operation: str | None = None,
    sequence: str,
    frame: str | None | None = None,
    genetic_code: int | None | None = None,
    verbose: bool = True,
) -> _Object:
    'Translate a DNA sequence to protein in all 3 forward reading frames and identify the longest ORF. Uses the standard genetic code (NCBI Code 1). Returns protein sequences for each frame, longest ORF per frame, and the best frame with longest uninterrupted protein. Use for: reading frame analysis, finding the correct ORF in an unknown sequence, comparing all three frames.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['sequence'] = sequence
    if frame is not None:
        _dnakit_arguments['frame'] = frame
    if genetic_code is not None:
        _dnakit_arguments['genetic_code'] = genetic_code
    return _call('dna_translate_reading_frames', _dnakit_arguments, verbose=verbose)

def dna_virtual_digest(*,
    operation: str | None = None,
    sequence: str,
    enzymes: list[_Object] | None | None = None,
    circular: bool | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Perform a virtual restriction enzyme digest on a DNA sequence using the NEB enzyme library (25 curated enzymes plus Biopython's ~600-enzyme library, including Type IIS enzymes such as BsaI, BbsI, Esp3I/BsmBI and SapI, which cut outside their recognition site). Returns all resulting fragments with their sequences, lengths, and positions, plus a list of cut sites. Supports both linear and circular DNA. Use for: cloning planning, gel electrophoresis prediction, verifying construct integrity.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['sequence'] = sequence
    if enzymes is not None:
        _dnakit_arguments['enzymes'] = enzymes
    if circular is not None:
        _dnakit_arguments['circular'] = circular
    return _call('dna_virtual_digest', _dnakit_arguments, verbose=verbose)

def dep_map_get_gene_dependencies(*,
    gene_symbol: str,
    model_id: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get gene dependency information from CRISPR knockout screens. Returns gene IDs and annotations. Negative gene effect scores indicate essentiality (cell death upon knockout).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_symbol'] = gene_symbol
    if model_id is not None:
        _dnakit_arguments['model_id'] = model_id
    return _call('dep_map_get_gene_dependencies', _dnakit_arguments, verbose=verbose)

def dep_map_search_cell_lines(*,
    query: str,
    verbose: bool = True,
) -> _Object:
    'Search cancer cell lines by name. Returns matching cell lines with IDs and cancer types. Use for finding cell models by common names (A549, HeLa, MCF7, etc.).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    return _call('dep_map_search_cell_lines', _dnakit_arguments, verbose=verbose)

def dis_ge_net_get_disease_genes(*,
    operation: str | None = None,
    disease: str,
    min_score: float | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get all genes associated with a disease from DisGeNET. Returns ranked gene list with association scores and evidence counts. Filter by minimum score to get high-confidence associations. Useful for disease gene prioritization.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['disease'] = disease
    if min_score is not None:
        _dnakit_arguments['min_score'] = min_score
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('dis_ge_net_get_disease_genes', _dnakit_arguments, verbose=verbose)

def dis_ge_net_get_gda(*,
    operation: str | None = None,
    gene: str | None = None,
    disease: str | None = None,
    source: str | None = None,
    min_score: float | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get gene-disease associations with filtering. Filter by data source (CURATED, LITERATURE, ANIMAL_MODELS) and minimum score. Returns comprehensive GDA data with evidence levels. Requires DISGENET_API_KEY.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    if disease is not None:
        _dnakit_arguments['disease'] = disease
    if source is not None:
        _dnakit_arguments['source'] = source
    if min_score is not None:
        _dnakit_arguments['min_score'] = min_score
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('dis_ge_net_get_gda', _dnakit_arguments, verbose=verbose)

def dis_ge_net_get_vda(*,
    operation: str | None = None,
    variant: str | None = None,
    gene: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get variant-disease associations from DisGeNET. Query by rsID or gene symbol to find disease-associated variants. Returns variant IDs, diseases, and evidence scores. Useful for pharmacogenomics and variant interpretation.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if variant is not None:
        _dnakit_arguments['variant'] = variant
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('dis_ge_net_get_vda', _dnakit_arguments, verbose=verbose)

def dis_ge_net_search_disease(*,
    operation: str | None = None,
    disease: str,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Search DisGeNET for genes associated with a disease. Requires a UMLS CUI, NOT a free-text disease name (e.g. 'C0006142' for breast cancer, not 'breast cancer') -- resolve a disease name to its CUI first via umls_search_concepts. Returns associated genes with scores and evidence. Requires DISGENET_API_KEY.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['disease'] = disease
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('dis_ge_net_search_disease', _dnakit_arguments, verbose=verbose)

def dis_ge_net_search_gene(*,
    operation: str | None = None,
    gene: str,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Search DisGeNET for diseases associated with a gene. DisGeNET is the largest public database of gene-disease associations. Returns diseases, association scores, and evidence counts. Use the 'gene' parameter (e.g., 'BRCA1', 'TP53') — note: parameter name is 'gene', not 'gene_symbol'. Requires DISGENET_API_KEY env var (register free at disgenet.org).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['gene'] = gene
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('dis_ge_net_search_gene', _dnakit_arguments, verbose=verbose)

def dyna_mut2_predict_stability(*,
    operation: str | None = None,
    pdb_id: str,
    chain: str,
    mutation: str,
    verbose: bool = True,
) -> _Object:
    'Predict the effect of a single-point amino acid mutation on protein stability and dynamics using DynaMut2. Provide a PDB ID (4-character code from RCSB PDB), chain identifier, and mutation in the format WtPosNew (e.g., V1A means Valine at position 1 mutated to Alanine). Returns the predicted change in Gibbs free energy (ddG in kcal/mol): positive = stabilizing, negative = destabilizing. The tool downloads the PDB structure from RCSB and submits it to DynaMut2 for analysis. Computation typically takes 15-30 seconds. Based on Normal Mode Analysis and graph-based signatures.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['pdb_id'] = pdb_id
    _dnakit_arguments['chain'] = chain
    _dnakit_arguments['mutation'] = mutation
    return _call('dyna_mut2_predict_stability', _dnakit_arguments, verbose=verbose)

def enaportal_search_samples(*,
    query: str,
    limit: int | None | None = None,
    fields: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    'Search the European Nucleotide Archive (ENA) for biological samples by taxonomy or text query. Returns sample accessions, descriptions, taxonomy, and publication dates. Useful for finding source material metadata across the world\'s nucleotide archives (INSDC). Example: \'tax_tree(562)\' finds all E. coli samples; \'description="liver"\' finds samples with liver in the description.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if fields is not None:
        _dnakit_arguments['fields'] = fields
    return _call('enaportal_search_samples', _dnakit_arguments, verbose=verbose)

def enaportal_search_studies(*,
    query: str,
    limit: int | None | None = None,
    fields: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    'Search the European Nucleotide Archive (ENA) for sequencing studies using text queries or taxonomy filters. Returns study accessions, titles, submitting centers, and publication dates. ENA is the world\'s largest repository of nucleotide sequences (>3 billion records). Supports powerful query syntax: use \'description="keyword"\' for text search, \'tax_tree(TAXID)\' for taxonomy-based queries, or combine with AND/OR. Example: searching \'description="SARS-CoV-2"\' returns COVID-19 sequencing studies.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if fields is not None:
        _dnakit_arguments['fields'] = fields
    return _call('enaportal_search_studies', _dnakit_arguments, verbose=verbose)

def encode_get_biosample(*,
    accession: str,
    verbose: bool = True,
) -> _Object:
    'Get detailed metadata for a specific ENCODE biosample by its accession ID. Returns comprehensive information including biosample type, organism, life stage, age, sex, tissue/cell line name, donor information, treatments, genetic modifications, culture conditions, and characterizations. Essential for understanding sample provenance, experimental conditions, and biological context of ENCODE datasets. Find biosample accessions using ENCODE_search_biosamples (format: ENCBS######).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['accession'] = accession
    return _call('encode_get_biosample', _dnakit_arguments, verbose=verbose)

def encode_get_chromatin_state(*,
    biosample_term_name: str | None | None = None,
    organism: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search ENCODE chromatin state annotations (ChromHMM segmentations) for specific biosamples or tissues. Chromatin states are genome-wide annotations that classify each genomic position into a functional state (e.g., active TSS, flanking active TSS, strong transcription, weak transcription, genic enhancers, active enhancers, bivalent enhancer, bivalent/poised TSS, flanking bivalent TSS, repressed polycomb, weak repressed polycomb, quiescent, heterochromatin) based on combinatorial patterns of histone marks. Use this to understand the chromatin landscape of specific cell types.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if biosample_term_name is not None:
        _dnakit_arguments['biosample_term_name'] = biosample_term_name
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('encode_get_chromatin_state', _dnakit_arguments, verbose=verbose)

def encode_get_experiment(*,
    accession: str,
    verbose: bool = True,
) -> _Object:
    'Get comprehensive metadata for a specific ENCODE experiment by its accession ID. Returns detailed information including assay type, target protein, biosample details, replicates, files, quality metrics, analysis pipeline, lab, award, and publication references. Use this to understand experimental design, data quality, and available files before downloading. Find experiment accessions using ENCODE_search_experiments (format: ENCSR######).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['accession'] = accession
    return _call('encode_get_experiment', _dnakit_arguments, verbose=verbose)

def encode_list_files(*,
    file_type: str | None = None,
    assay_title: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'List ENCODE data files with filters by file format, output type, or assay. Returns file metadata including accessions, formats (FASTQ, BAM, bigWig), output types (alignments, peaks, signal), and file sizes. Use this to programmatically discover downloadable artifacts for specific experiments or data types. Find file accessions to use with ENCODE_get_file for detailed file information and download URLs.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if file_type is not None:
        _dnakit_arguments['file_type'] = file_type
    if assay_title is not None:
        _dnakit_arguments['assay_title'] = assay_title
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('encode_list_files', _dnakit_arguments, verbose=verbose)

def encode_search_annotations(*,
    annotation_type: str | None = None,
    biosample_term_name: str | None | None = None,
    organism: str | None = None,
    assembly: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search ENCODE annotations including candidate cis-Regulatory Elements (cCREs), chromatin states, and imputed signals. cCREs are classified as promoter-like (PLS), proximal enhancer-like (pELS), distal enhancer-like (dELS), DNase-H3K4me3 (DNase-H3K4me3), and CTCF-only. Chromatin state annotations use ChromHMM models to segment the genome into functional states. Use this to find regulatory element annotations for specific biosamples or genome assemblies.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if annotation_type is not None:
        _dnakit_arguments['annotation_type'] = annotation_type
    if biosample_term_name is not None:
        _dnakit_arguments['biosample_term_name'] = biosample_term_name
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if assembly is not None:
        _dnakit_arguments['assembly'] = assembly
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('encode_search_annotations', _dnakit_arguments, verbose=verbose)

def encode_search_biosamples(*,
    organism: str | None = None,
    biosample_type: str | None = None,
    treatment: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search ENCODE biosamples (cell lines, tissues, primary cells) by organism, biosample type, or treatment. Returns biosample accessions, types, organism, life stage, and descriptions. Use this to discover biological samples used in ENCODE experiments, find cell lines or tissues of interest, or identify treated vs. control samples. Find biosample accessions to use with ENCODE_get_biosample for detailed characterization data.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if biosample_type is not None:
        _dnakit_arguments['biosample_type'] = biosample_type
    if treatment is not None:
        _dnakit_arguments['treatment'] = treatment
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('encode_search_biosamples', _dnakit_arguments, verbose=verbose)

def encode_search_chromatin_accessibility(*,
    assay_type: str | None = None,
    biosample_term_name: str | None | None = None,
    organism: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search ENCODE chromatin accessibility experiments (ATAC-seq and DNase-seq). ATAC-seq and DNase-seq identify open chromatin regions including promoters, enhancers, and other regulatory elements. Returns experiment accessions, biosample details, and accessibility assay metadata. Open chromatin regions are sites of active gene regulation where transcription factors can bind. Use this to find chromatin accessibility data for specific tissues or cell types.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if assay_type is not None:
        _dnakit_arguments['assay_type'] = assay_type
    if biosample_term_name is not None:
        _dnakit_arguments['biosample_term_name'] = biosample_term_name
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('encode_search_chromatin_accessibility', _dnakit_arguments, verbose=verbose)

def encode_search_experiments(*,
    assay_title: str | None = None,
    target: str | None = None,
    organism: str | None = None,
    status: str | None = None,
    limit: int | None = None,
    biosample: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Search ENCODE experiments\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if assay_title is not None:
        _dnakit_arguments['assay_title'] = assay_title
    if target is not None:
        _dnakit_arguments['target'] = target
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if status is not None:
        _dnakit_arguments['status'] = status
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if biosample is not None:
        _dnakit_arguments['biosample'] = biosample
    return _call('encode_search_experiments', _dnakit_arguments, verbose=verbose)

def encode_search_histone_experiments(*,
    histone_mark: str | None | None = None,
    target: str | None | None = None,
    biosample_term_name: str | None | None = None,
    biosample: str | None | None = None,
    organism: str | None = None,
    limit: int | None = None,
    biosample_term: str | None | None = None,
    cell_type: str | None | None = None,
    tissue: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    'Search ENCODE histone ChIP-seq experiments by histone modification mark, biosample, or organism. Returns experiment accessions, histone marks (H3K4me3, H3K27ac, H3K27me3, H3K36me3, H3K4me1, H3K9me3), biosample summaries, and metadata. Use this to find histone modification profiling data for specific marks or tissues. Common marks: H3K4me3 (active promoters), H3K27ac (active enhancers), H3K27me3 (polycomb repression), H3K4me1 (poised enhancers), H3K36me3 (gene bodies), H3K9me3 (heterochromatin).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if histone_mark is not None:
        _dnakit_arguments['histone_mark'] = histone_mark
    if target is not None:
        _dnakit_arguments['target'] = target
    if biosample_term_name is not None:
        _dnakit_arguments['biosample_term_name'] = biosample_term_name
    if biosample is not None:
        _dnakit_arguments['biosample'] = biosample
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if biosample_term is not None:
        _dnakit_arguments['biosample_term'] = biosample_term
    if cell_type is not None:
        _dnakit_arguments['cell_type'] = cell_type
    if tissue is not None:
        _dnakit_arguments['tissue'] = tissue
    return _call('encode_search_histone_experiments', _dnakit_arguments, verbose=verbose)

def encode_search_rnaseq_experiments(*,
    biosample_term_name: str | None | None = None,
    biosample: str | None | None = None,
    cell_type: str | None | None = None,
    tissue: str | None | None = None,
    organism: str | None = None,
    assay_type: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search ENCODE RNA-seq experiments by biosample (cell type or tissue), organism, and assay type (total RNA-seq, polyA plus RNA-seq, small RNA-seq, microRNA-seq). Returns experiment accessions, biosample summaries, lab names, and metadata from the ENCODE project. Use this to find ENCODE-curated RNA-seq data for specific cell lines or tissues.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if biosample_term_name is not None:
        _dnakit_arguments['biosample_term_name'] = biosample_term_name
    if biosample is not None:
        _dnakit_arguments['biosample'] = biosample
    if cell_type is not None:
        _dnakit_arguments['cell_type'] = cell_type
    if tissue is not None:
        _dnakit_arguments['tissue'] = tissue
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if assay_type is not None:
        _dnakit_arguments['assay_type'] = assay_type
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('encode_search_rnaseq_experiments', _dnakit_arguments, verbose=verbose)

def esm2_score_missense_variant(*,
    sequence: str,
    position: int,
    mutant: str,
    wild_type: str | None | None = None,
    model_id: str | None | None = None,
    wait_for_model: bool | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Score a missense protein variant with ESM-2 masked-marginal log-likelihood ratio (Meier 2021) over HuggingFace's keyless hf-inference. Input wild-type 1-letter sequence + 1-based position + mutant residue; returns LLR = logP(mut) - logP(wt). Negative => mutant disfavored (candidate deleterious). No API key. For ESMC/SAE scoring with a key use ESM_score_sequence.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['sequence'] = sequence
    _dnakit_arguments['position'] = position
    _dnakit_arguments['mutant'] = mutant
    if wild_type is not None:
        _dnakit_arguments['wild_type'] = wild_type
    if model_id is not None:
        _dnakit_arguments['model_id'] = model_id
    if wait_for_model is not None:
        _dnakit_arguments['wait_for_model'] = wait_for_model
    return _call('esm2_score_missense_variant', _dnakit_arguments, verbose=verbose)

def esm_describe_sae_feature(*,
    feature_id: int,
    sae_model: str | None = None,
    model: str | None = None,
    n_proteins: int | None = None,
    top_residues_per_protein: int | None = None,
    use_cache: bool | None = None,
    verbose: bool = True,
) -> _Object:
    "Label a single SAE feature_id with its dominant biological category by aggregating UniProt feature-type overlap at SAE-activating residues across a curated 10-protein panel. Returns one of: catalytic, ligand-binding, ptm, domain, motif, structural-stability, secondary-structure, transmembrane, signal-peptide, propeptide, or 'uncategorized' (no informative UniProt overlap found). Use when ESM_get_sae_features / ESM_score_variant_sae_disruption returned a feature_id that needs biological interpretation. Cost: 10 Forge credits FIRST time (1 per panel protein). Subsequent calls for the same feature_id are FREE (filesystem cache at ~/.cache/tooluniverse/sae_labels/{sae_model}/feature_{id}.json). Latency: ~30-60s first call (panel of 10 SAE inferences), <1s on cache hit. Prerequisites: same as ESM_get_sae_features. License: Cambrian Inference License — non-commercial.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['feature_id'] = feature_id
    if sae_model is not None:
        _dnakit_arguments['sae_model'] = sae_model
    if model is not None:
        _dnakit_arguments['model'] = model
    if n_proteins is not None:
        _dnakit_arguments['n_proteins'] = n_proteins
    if top_residues_per_protein is not None:
        _dnakit_arguments['top_residues_per_protein'] = top_residues_per_protein
    if use_cache is not None:
        _dnakit_arguments['use_cache'] = use_cache
    return _call('esm_describe_sae_feature', _dnakit_arguments, verbose=verbose)

def esm_explain_variant_mechanism(*,
    sequence: str,
    position: int,
    ref_aa: str,
    alt_aa: str,
    window: int | None = None,
    top_k_features: int | None = None,
    include_descriptions: bool | None = None,
    model: str | None = None,
    sae_model: str | None = None,
    verbose: bool = True,
) -> _Object:
    "One-call composite for variant mechanism: runs ESMC-6B SAE variant disruption + describe_sae_feature on each top affected feature + composes a category-level mechanism summary (e.g. 'Disrupted feature categories (lost): catalytic=2, ligand-binding=1'). Use when a variant-interpretation skill wants the mechanism in a single tool call rather than orchestrating disruption + N description calls. Set include_descriptions=false to skip the category labeling step (2 Forge calls only). Otherwise: 2 disruption calls + up to (2 * top_k_features) describe calls (most cached after first use). Prerequisites: pip install 'esm @ git+https://github.com/evolutionaryscale/esm@ee891c52'; set ESM_API_KEY. Outputs governed by EvolutionaryScale Cambrian Inference License — non-commercial use only unless covered by a separate commercial agreement.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['sequence'] = sequence
    _dnakit_arguments['position'] = position
    _dnakit_arguments['ref_aa'] = ref_aa
    _dnakit_arguments['alt_aa'] = alt_aa
    if window is not None:
        _dnakit_arguments['window'] = window
    if top_k_features is not None:
        _dnakit_arguments['top_k_features'] = top_k_features
    if include_descriptions is not None:
        _dnakit_arguments['include_descriptions'] = include_descriptions
    if model is not None:
        _dnakit_arguments['model'] = model
    if sae_model is not None:
        _dnakit_arguments['sae_model'] = sae_model
    return _call('esm_explain_variant_mechanism', _dnakit_arguments, verbose=verbose)

def esm_get_region_sae_features(*,
    sequence: str,
    start_position: int,
    end_position: int,
    top_k_features: int | None = None,
    model: str | None = None,
    sae_model: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Aggregate ESMC-6B SAE features over a contiguous residue range to characterize the region's biological identity. Returns features ranked by total |activation| across the region, with each feature's per-residue hit pattern. Use for: interpreting a known domain (residues X-Y), profiling an antibody epitope, characterizing a binding pocket, comparing two regions of the same protein, or seeding ESM_describe_sae_feature calls on the top-K to get category labels. Max region length 500 residues. Prerequisites: pip install 'esm @ git+https://github.com/evolutionaryscale/esm@ee891c52'; set ESM_API_KEY. Outputs governed by EvolutionaryScale Cambrian Inference License — non-commercial use only unless covered by a separate commercial agreement.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['sequence'] = sequence
    _dnakit_arguments['start_position'] = start_position
    _dnakit_arguments['end_position'] = end_position
    if top_k_features is not None:
        _dnakit_arguments['top_k_features'] = top_k_features
    if model is not None:
        _dnakit_arguments['model'] = model
    if sae_model is not None:
        _dnakit_arguments['sae_model'] = sae_model
    return _call('esm_get_region_sae_features', _dnakit_arguments, verbose=verbose)

def esm_get_sae_features(*,
    sequence: str,
    model: str | None = None,
    sae_model: str | None = None,
    position: int | None | None = None,
    window: int | None = None,
    top_k_per_residue: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Run a protein sequence through an ESMC Sparse Autoencoder (SAE) and return sparse feature activations per residue. SAEs decompose the protein language model's hidden states into a 16,384-feature sparse codebook (top-k=64 active per residue). Each feature is an interpretable latent — many activate on specific biological categories (catalytic site, ligand-binding region, PTM sequon, structural motif, etc.) that can be looked up via ESM_describe_sae_feature. Use for: variant interpretation (compare ref vs mutant activations), mechanism inference (which functional features the protein engages), protein-LM interpretability research. Prerequisites: (1) pip install 'esm @ git+https://github.com/evolutionaryscale/esm@ee891c52' — the PyPI esm release does NOT yet include SAE support, which lives on an upstream feature branch. (2) ESM_API_KEY env var with a Forge token (https://forge.evolutionaryscale.ai). License note: SAE outputs are governed by the EvolutionaryScale Cambrian Inference Clickthrough License — non-commercial / academic use only.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['sequence'] = sequence
    if model is not None:
        _dnakit_arguments['model'] = model
    if sae_model is not None:
        _dnakit_arguments['sae_model'] = sae_model
    if position is not None:
        _dnakit_arguments['position'] = position
    if window is not None:
        _dnakit_arguments['window'] = window
    if top_k_per_residue is not None:
        _dnakit_arguments['top_k_per_residue'] = top_k_per_residue
    return _call('esm_get_sae_features', _dnakit_arguments, verbose=verbose)

def esm_score_sequence(*,
    sequence: str,
    model: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Score a protein sequence using ESMC logits to compute per-residue log-probabilities and mean pseudo-log-likelihood. Higher (less negative) mean log-probability indicates the sequence is more natural/likely according to the model's learned protein distribution. Use for: evaluating mutant fitness, comparing sequence variants, identifying low-confidence residues. Mean log-probability near 0 is highest confidence; values below -4 suggest unusual or low-quality residues. Prerequisites: pip install esm; set ESM_API_KEY (from https://forge.evolutionaryscale.ai).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['sequence'] = sequence
    if model is not None:
        _dnakit_arguments['model'] = model
    return _call('esm_score_sequence', _dnakit_arguments, verbose=verbose)

def esm_score_variant_sae_batch(*,
    sequence: str,
    variants: list[dict[str, _Object]],
    window: int | None = None,
    top_k_features: int | None = None,
    model: str | None = None,
    sae_model: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Score many missense variants against one reference protein using ESMC-6B SAE. More Forge-efficient than calling ESM_score_variant_sae_disruption N times: runs the reference SAE once and reuses it across every variant, giving N+1 Forge calls instead of 2N. Use for saturation mutagenesis at a single position (all 19 alternative residues), DMS-style sweeps, clinical variant panels, or any workflow scoring multiple variants on the same protein. Max 100 variants per call. Prerequisites: pip install 'esm @ git+https://github.com/evolutionaryscale/esm@ee891c52' (SAE support is on an unmerged feature branch — PyPI esm 3.2.x does NOT include SAEConfig). Set ESM_API_KEY environment variable. Outputs governed by EvolutionaryScale Cambrian Inference License — non-commercial use only unless covered by a separate commercial agreement.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['sequence'] = sequence
    _dnakit_arguments['variants'] = variants
    if window is not None:
        _dnakit_arguments['window'] = window
    if top_k_features is not None:
        _dnakit_arguments['top_k_features'] = top_k_features
    if model is not None:
        _dnakit_arguments['model'] = model
    if sae_model is not None:
        _dnakit_arguments['sae_model'] = sae_model
    return _call('esm_score_variant_sae_batch', _dnakit_arguments, verbose=verbose)

def esm_score_variant_sae_disruption(*,
    sequence: str,
    position: int,
    ref_aa: str,
    alt_aa: str,
    window: int | None = None,
    top_k_features: int | None = None,
    model: str | None = None,
    sae_model: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Composite SAE-based variant scoring. Given a protein sequence and a missense variant (position + ref_aa + alt_aa), runs ESMC-6B SAE on both reference and mutant sequences, computes per-feature activation deltas summed over a window centered on the mutation site, and returns ranked top-K features LOST and GAINED in the mutant. This is the convenience layer over ESM_get_sae_features — use this for one-shot variant interpretation; use ESM_get_sae_features directly only when you need raw per-residue features for non-variant analyses. Validates ref_aa matches the position in the supplied sequence (returns clear error if not, so you can detect transcript / isoform mismatches). Cost: 2 Forge API credits (1 ref + 1 mut). Latency: ~3-6s for typical-length human proteins. Prerequisites: same as ESM_get_sae_features (pip install 'esm @ git+https://github.com/evolutionaryscale/esm@ee891c52', ESM_API_KEY env var). License: SAE outputs governed by Cambrian Inference License — non-commercial / academic only.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['sequence'] = sequence
    _dnakit_arguments['position'] = position
    _dnakit_arguments['ref_aa'] = ref_aa
    _dnakit_arguments['alt_aa'] = alt_aa
    if window is not None:
        _dnakit_arguments['window'] = window
    if top_k_features is not None:
        _dnakit_arguments['top_k_features'] = top_k_features
    if model is not None:
        _dnakit_arguments['model'] = model
    if sae_model is not None:
        _dnakit_arguments['sae_model'] = sae_model
    return _call('esm_score_variant_sae_disruption', _dnakit_arguments, verbose=verbose)

def eve_get_gene_info(*,
    gene_symbol: str,
    verbose: bool = True,
) -> _Object:
    'Check if EVE scores are available for a gene and get gene information. EVE covers ~3,000 disease-related human genes. Returns Ensembl ID and link to EVE website.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_symbol'] = gene_symbol
    return _call('eve_get_gene_info', _dnakit_arguments, verbose=verbose)

def eve_get_variant_score(*,
    variant: str | None = None,
    chrom: str | None = None,
    pos: int | None = None,
    ref: str | None = None,
    alt: str | None = None,
    species: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get EVE (Evolutionary Variant Effect) pathogenicity score for a variant. EVE uses unsupervised deep learning on evolutionary data. Scores: 0-1 where >0.5 is likely pathogenic. Accessed via Ensembl VEP API.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if variant is not None:
        _dnakit_arguments['variant'] = variant
    if chrom is not None:
        _dnakit_arguments['chrom'] = chrom
    if pos is not None:
        _dnakit_arguments['pos'] = pos
    if ref is not None:
        _dnakit_arguments['ref'] = ref
    if alt is not None:
        _dnakit_arguments['alt'] = alt
    if species is not None:
        _dnakit_arguments['species'] = species
    return _call('eve_get_variant_score', _dnakit_arguments, verbose=verbose)

def enrichr_enrich(*,
    operation: str | None = None,
    gene_list: list[str],
    library: str | None = None,
    top_n: int | None = None,
    species: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    'Perform gene set enrichment analysis using Enrichr. Submit a list of genes and get enriched terms from a specific library. Returns ranked terms with p-values, z-scores, combined scores, and overlapping genes. No authentication required. Set species to query the fly, worm, yeast, or fish Enrichr instances; libraries and gene symbols differ per organism.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['gene_list'] = gene_list
    if library is not None:
        _dnakit_arguments['library'] = library
    if top_n is not None:
        _dnakit_arguments['top_n'] = top_n
    if species is not None:
        _dnakit_arguments['species'] = species
    return _call('enrichr_enrich', _dnakit_arguments, verbose=verbose)

def ensembl_compara_get_gene_tree(*,
    gene: str,
    species: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get the phylogenetic gene tree showing evolutionary relationships of a gene family across species from Ensembl Compara. Gene trees trace the history of speciation and duplication events for a gene family. Returns tree members (genes across species), tree structure, and optionally Newick format. Useful for evolutionary analysis and understanding gene family expansion. Example: TP53 gene tree shows conservation across vertebrates.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene'] = gene
    if species is not None:
        _dnakit_arguments['species'] = species
    return _call('ensembl_compara_get_gene_tree', _dnakit_arguments, verbose=verbose)

def ensembl_compara_get_orthologues(*,
    gene: str,
    species: str | None = None,
    target_species: str | None = None,
    target_taxon: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Find orthologous genes (genes in different species that evolved from a common ancestor) using Ensembl Compara. Essential for comparative genomics - finding the mouse, zebrafish, or fly equivalent of a human gene. Returns orthologue gene IDs, species, homology type (one-to-one, one-to-many), and taxonomy level of divergence. Example: TP53 in human has a one-to-one orthologue (Trp53) in mouse.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene'] = gene
    if species is not None:
        _dnakit_arguments['species'] = species
    if target_species is not None:
        _dnakit_arguments['target_species'] = target_species
    if target_taxon is not None:
        _dnakit_arguments['target_taxon'] = target_taxon
    return _call('ensembl_compara_get_orthologues', _dnakit_arguments, verbose=verbose)

def ensembl_compara_get_paralogues(*,
    gene: str,
    species: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Find paralogous genes (genes within the same species that arose from gene duplication events) using Ensembl Compara. Paralogues often have related but diverged functions. Returns paralogue gene IDs, duplication type, and evolutionary context. Example: TP53 has paralogues TP63 (ENSG00000073282) and TP73 (ENSG00000078900) in human, from ancient vertebrate gene duplication.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene'] = gene
    if species is not None:
        _dnakit_arguments['species'] = species
    return _call('ensembl_compara_get_paralogues', _dnakit_arguments, verbose=verbose)

def ensembl_seq_get_id_sequence(*,
    ensembl_id: str,
    type: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Get the amino acid or nucleotide sequence for an Ensembl ID (protein, transcript, or gene). Retrieves protein sequences from ENSP IDs, cDNA/CDS from ENST IDs, or genomic DNA from ENSG IDs. Example: ENSP00000269305 returns the 393 amino acid protein sequence of TP53 canonical isoform. ENST00000269305 with type='cdna' returns the 2512 bp cDNA sequence of TP53-201.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['ensembl_id'] = ensembl_id
    if type is not None:
        _dnakit_arguments['type'] = type
    return _call('ensembl_seq_get_id_sequence', _dnakit_arguments, verbose=verbose)

def ensembl_seq_get_region_sequence(*,
    region: str,
    species: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Get the DNA nucleotide sequence for a specific genomic region from Ensembl. Retrieves the reference genome sequence for given chromosome coordinates. Useful for examining regulatory regions, intergenic sequences, or any arbitrary genomic interval. Example: region '17:7668421-7668520' in homo_sapiens returns 100 bp of the TP53 locus. Supports any Ensembl species and coordinate strand.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['region'] = region
    if species is not None:
        _dnakit_arguments['species'] = species
    return _call('ensembl_seq_get_region_sequence', _dnakit_arguments, verbose=verbose)

def ensembl_vep_annotate_hgvs(*,
    hgvs_notation: str,
    species: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Predict functional consequences of a genetic variant using HGVS notation via the Ensembl Variant Effect Predictor (VEP). Returns impact on genes, transcripts, and proteins including SIFT/PolyPhen pathogenicity scores, consequence types, amino acid changes, and colocated known variants. HGVS notation examples: 'BRAF:p.Val600Glu' (protein), 'ENST00000366667:c.803C>T' (coding DNA). This complements the existing region-based VEP tool by accepting standard clinical variant nomenclature.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['hgvs_notation'] = hgvs_notation
    if species is not None:
        _dnakit_arguments['species'] = species
    return _call('ensembl_vep_annotate_hgvs', _dnakit_arguments, verbose=verbose)

def ensembl_vep_annotate_rsid(*,
    variant_id: str | None = None,
    rsid: str | None = None,
    rs_id: str | None = None,
    species: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Predict functional consequences of a genetic variant using its dbSNP rs identifier via the Ensembl Variant Effect Predictor (VEP). Returns impact on genes, transcripts, and proteins including SIFT/PolyPhen scores, consequence types, amino acid changes, and colocated variants. Useful when you have a variant's rs number from GWAS, clinical reports, or literature. Example: 'rs7903146' (TCF7L2 T2D variant), 'rs429358' (APOE Alzheimer variant).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if variant_id is not None:
        _dnakit_arguments['variant_id'] = variant_id
    if rsid is not None:
        _dnakit_arguments['rsid'] = rsid
    if rs_id is not None:
        _dnakit_arguments['rs_id'] = rs_id
    if species is not None:
        _dnakit_arguments['species'] = species
    return _call('ensembl_vep_annotate_rsid', _dnakit_arguments, verbose=verbose)

def ensembl_vep_variant_recoder(*,
    variant_id: str,
    species: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Convert genetic variant identifiers between different nomenclature formats using the Ensembl Variant Recoder. Given a variant ID (e.g., rs number), returns equivalent representations in HGVS genomic (hgvsg), HGVS coding (hgvsc), HGVS protein (hgvsp), and SPDI format. Essential for cross-referencing variants across databases that use different notation systems. Example: 'rs429358' returns HGVS coding notation for all APOE transcripts.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['variant_id'] = variant_id
    if species is not None:
        _dnakit_arguments['species'] = species
    return _call('ensembl_vep_variant_recoder', _dnakit_arguments, verbose=verbose)

def ensembl_var_get_variant_consequences(*,
    variant_id: str,
    species: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Get detailed variant information including consequence type, clinical significance, synonyms, evidence attributes, and allele mappings from Ensembl. More comprehensive than VEP - returns the full variant record with all known synonyms, mapped positions, ancestral allele, and minor allele data. Example: rs429358 returns most_severe_consequence='missense_variant', 27 synonyms, mapped to chr19:44908684 (GRCh38), and is classified as a missense_variant affecting APOE.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['variant_id'] = variant_id
    if species is not None:
        _dnakit_arguments['species'] = species
    return _call('ensembl_var_get_variant_consequences', _dnakit_arguments, verbose=verbose)

def europe_pmc_search_articles(*,
    query: str,
    limit: int | None = None,
    require_has_ft: bool | None = None,
    fulltext_terms: list[str] | None = None,
    enrich_missing_abstract: bool | None = None,
    extract_terms_from_fulltext: list[str] | None = None,
    page_size: int | None = None,
    pageSize: int | None = None,
    source: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Search for articles on Europe PMC including abstracts and metadata. Returns the top `limit` results ranked by relevance, which is normally a SLICE of the matches: `metadata.total_results` carries Europe PMC\'s own hitCount (every article matching the query) while `metadata.count` is only how many were returned, and a top-level `truncated: true` plus `truncation_note` appear whenever matches were left behind. Do not read a small `count` as \'few articles exist\' — check `total_results`. Europe PMC supports fielded queries for indexed full text (e.g., BODY:"term") but only for records where Europe PMC has full text indexed (HAS_FT:Y). Use require_has_ft/fulltext_terms to avoid confusing \'has a full-text link somewhere\' with \'full text searchable in Europe PMC\'.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if require_has_ft is not None:
        _dnakit_arguments['require_has_ft'] = require_has_ft
    if fulltext_terms is not None:
        _dnakit_arguments['fulltext_terms'] = fulltext_terms
    if enrich_missing_abstract is not None:
        _dnakit_arguments['enrich_missing_abstract'] = enrich_missing_abstract
    if extract_terms_from_fulltext is not None:
        _dnakit_arguments['extract_terms_from_fulltext'] = extract_terms_from_fulltext
    if page_size is not None:
        _dnakit_arguments['page_size'] = page_size
    if pageSize is not None:
        _dnakit_arguments['pageSize'] = pageSize
    if source is not None:
        _dnakit_arguments['source'] = source
    return _call('europe_pmc_search_articles', _dnakit_arguments, verbose=verbose)

def evo2_score_variant(*,
    ref_sequence: str | None | None = None,
    alt_sequence: str | None | None = None,
    sequence: str | None | None = None,
    position: int | None | None = None,
    reference: str | None | None = None,
    alternate: str | None | None = None,
    model: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Zero-shot variant-effect scoring with NVIDIA-hosted Evo 2 (Arc Institute genome foundation model). Runs the model's forward pass on the reference and alternate DNA windows and returns the delta log-likelihood (alt - ref). Negative = variant disfavored by the genome model (candidate deleterious). Two hosted sizes are selectable via `model`: evo2-40b (default, most accurate) and evo2-7b (smaller, faster). Hosted NIM API; requires NVIDIA_API_KEY. Complements NvidiaNIM_evo2 (which only generates sequences).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if ref_sequence is not None:
        _dnakit_arguments['ref_sequence'] = ref_sequence
    if alt_sequence is not None:
        _dnakit_arguments['alt_sequence'] = alt_sequence
    if sequence is not None:
        _dnakit_arguments['sequence'] = sequence
    if position is not None:
        _dnakit_arguments['position'] = position
    if reference is not None:
        _dnakit_arguments['reference'] = reference
    if alternate is not None:
        _dnakit_arguments['alternate'] = alternate
    if model is not None:
        _dnakit_arguments['model'] = model
    return _call('evo2_score_variant', _dnakit_arguments, verbose=verbose)

def expression_atlas_search_differential(*,
    gene: str | None = None,
    condition: str | None = None,
    species: str | None = None,
    limit: int | None | None = None,
    offset: int | None | None = None,
    verbose: bool = True,
) -> _Object:
    'Search for differential expression EXPERIMENTS from EBI Expression Atlas by condition text and/or a gene mentioned in the experiment description. Supports disease-context expression analysis. Caveat: `condition` does real text filtering, but `gene` is a text-search coincidence (experiment descriptions rarely name individual genes) and does not return actual per-gene differential-expression values -- treat gene-tagged results as a starting point for manual review, not a validated gene-specific result set.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    if condition is not None:
        _dnakit_arguments['condition'] = condition
    if species is not None:
        _dnakit_arguments['species'] = species
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    return _call('expression_atlas_search_differential', _dnakit_arguments, verbose=verbose)

def expression_atlas_search_experiments(*,
    gene: str | None = None,
    condition: str | None = None,
    species: str | None = None,
    limit: int | None | None = None,
    offset: int | None | None = None,
    verbose: bool = True,
) -> _Object:
    'Search all EBI Expression Atlas experiments (baseline + differential) by condition text and/or a gene mentioned in the experiment description. Returns experiment accessions, types, descriptions, and assay counts. Caveat: `condition` does real text filtering, but `gene` is a text-search coincidence (experiment descriptions rarely name individual genes), so a gene-only query effectively returns the unfiltered catalog -- use to find relevant expression datasets, not as a validated per-gene result set.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    if condition is not None:
        _dnakit_arguments['condition'] = condition
    if species is not None:
        _dnakit_arguments['species'] = species
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    return _call('expression_atlas_search_experiments', _dnakit_arguments, verbose=verbose)

def favor_annotate_variant(*,
    variant: str,
    verbose: bool = True,
) -> _Object:
    'Comprehensive functional annotation of a single GRCh38 variant via FAVOR (Functional Annotation of Variants Online Resource, Harvard). One call returns: allele frequencies (BRAVO/TOPMed, gnomAD by ancestry, 1000 Genomes), gene and consequence (GENCODE category, exonic effect, protein variant, HGVS), deleteriousness scores (CADD, SIFT, PolyPhen-2, AlphaMissense, MetaSVM, MutationTaster/Assessor, FATHMM-XF, Grantham, LINSIGHT, FunSeq), conservation (GERP, phyloP mammalian/primate/vertebrate, phastCons), ClinVar clinical significance/disease, and regulatory/epigenomic context (CAGE promoter/enhancer, GeneHancer, ENCODE DNase, ReMap TF overlap). Use for variant interpretation across both coding and non-coding regions in one step. Input is a GRCh38 (hg38) coordinate string; resolve rsIDs to GRCh38 coordinates first if needed. No API key required.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['variant'] = variant
    return _call('favor_annotate_variant', _dnakit_arguments, verbose=verbose)

def fda_get_active_ingredient_info_by_drug_name(*,
    drug_name: str,
    limit: int | None = None,
    skip: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Fetch a list of active ingredients in a specific drug product.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['drug_name'] = drug_name
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if skip is not None:
        _dnakit_arguments['skip'] = skip
    return _call('fda_get_active_ingredient_info_by_drug_name', _dnakit_arguments, verbose=verbose)

def fda_get_adverse_reactions_by_drug_name(*,
    drug_name: str,
    limit: int | None = None,
    skip: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Retrieve adverse reactions information based on the drug name. IMPORTANT -- FDA labels file safety content under different section names depending on the label's FORMAT, not the drug: modern PLR-format labels (Physician Labeling Rule, 2006 onward) use 'warnings_and_cautions' and usually carry NO 'warnings' section, while legacy and OTC-monograph labels use 'warnings'/'precautions' and usually carry no 'warnings_and_cautions'. Measured against openFDA drug/label on 2026-08-11 (search=_exists_:<field>): 261,639 labels total; 208,073 have 'warnings'; 46,925 have 'warnings_and_cautions'; 43,486 have 'warnings_and_cautions' but no 'warnings' at all; 33,054 have 'boxed_warning'. A null section therefore does NOT mean the drug lacks that content. When every requested section comes back null for a matched label but that same label carries a related section, the related section's content is added to the row under ITS OWN key (never rewritten into the requested key), its name is listed in 'related_sections_present', and a top-level 'section_note' explains the split and names the sibling tool that returns it. Example: GIAPREZA (angiotensin II) has warnings=null and boxed_warning=null but documents its thrombosis risk in 'warnings_and_cautions'. Related sections consulted for this tool: warnings, boxed_warning, precautions, general_precautions (tools: FDA_get_warnings_by_drug_name, FDA_get_boxed_warning_info_by_drug_name, FDA_get_precautions_by_drug_name, FDA_get_general_precautions_by_drug_name).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['drug_name'] = drug_name
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if skip is not None:
        _dnakit_arguments['skip'] = skip
    return _call('fda_get_adverse_reactions_by_drug_name', _dnakit_arguments, verbose=verbose)

def fda_get_boxed_warning_info_by_drug_name(*,
    drug_name: str,
    limit: int | None = None,
    skip: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Retrieve boxed warning and adverse effects information for a specific drug. IMPORTANT -- FDA labels file safety content under different section names depending on the label's FORMAT, not the drug: modern PLR-format labels (Physician Labeling Rule, 2006 onward) use 'warnings_and_cautions' and usually carry NO 'warnings' section, while legacy and OTC-monograph labels use 'warnings'/'precautions' and usually carry no 'warnings_and_cautions'. Measured against openFDA drug/label on 2026-08-11 (search=_exists_:<field>): 261,639 labels total; 208,073 have 'warnings'; 46,925 have 'warnings_and_cautions'; 43,486 have 'warnings_and_cautions' but no 'warnings' at all; 33,054 have 'boxed_warning'. A null section therefore does NOT mean the drug lacks that content. When every requested section comes back null for a matched label but that same label carries a related section, the related section's content is added to the row under ITS OWN key (never rewritten into the requested key), its name is listed in 'related_sections_present', and a top-level 'section_note' explains the split and names the sibling tool that returns it. Example: GIAPREZA (angiotensin II) has warnings=null and boxed_warning=null but documents its thrombosis risk in 'warnings_and_cautions'. Related sections consulted for this tool: warnings_and_cautions, warnings (tools: FDA_get_warnings_and_cautions_by_drug_name, FDA_get_warnings_by_drug_name).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['drug_name'] = drug_name
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if skip is not None:
        _dnakit_arguments['skip'] = skip
    return _call('fda_get_boxed_warning_info_by_drug_name', _dnakit_arguments, verbose=verbose)

def fda_get_clinical_studies_info_by_drug_name(*,
    drug_name: str,
    limit: int | None = None,
    skip: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Retrieve clinical studies information based on the drug name.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['drug_name'] = drug_name
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if skip is not None:
        _dnakit_arguments['skip'] = skip
    return _call('fda_get_clinical_studies_info_by_drug_name', _dnakit_arguments, verbose=verbose)

def fda_get_indications_by_drug_name(*,
    drug_name: str,
    limit: int | None = None,
    skip: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Retrieve indications and usage information based on a specific drug name.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['drug_name'] = drug_name
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if skip is not None:
        _dnakit_arguments['skip'] = skip
    return _call('fda_get_indications_by_drug_name', _dnakit_arguments, verbose=verbose)

def fda_get_mechanism_of_action_by_drug_name(*,
    drug_name: str,
    limit: int | None = None,
    skip: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Retrieve the mechanism of action information for a specific drug.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['drug_name'] = drug_name
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if skip is not None:
        _dnakit_arguments['skip'] = skip
    return _call('fda_get_mechanism_of_action_by_drug_name', _dnakit_arguments, verbose=verbose)

def finn_gen_get_variant_finemapping(*,
    variant: str,
    verbose: bool = True,
) -> _Object:
    'Get fine-mapping regions associated with a genomic variant in FinnGen. Returns phenotypes where this variant appears in a credible set from statistical fine-mapping. Input variant in chr:pos:ref:alt format (GRCh38). Useful for understanding which disease associations a variant contributes to in the Finnish population.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['variant'] = variant
    return _call('finn_gen_get_variant_finemapping', _dnakit_arguments, verbose=verbose)

def fly_base_get_gene(*,
    gene_id: str,
    verbose: bool = True,
) -> _Object:
    "Get detailed Drosophila melanogaster (fruit fly) gene information from FlyBase via the Alliance of Genome Resources. Returns gene symbol, name, synopsis, genomic location, synonyms, and cross-references. FlyBase is the primary database for Drosophila genetics. Example: 'FB:FBgn0003996' returns the white (w) gene with its function as an ABCG2-class transporter.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_id'] = gene_id
    return _call('fly_base_get_gene', _dnakit_arguments, verbose=verbose)

def fly_base_get_gene_alleles(*,
    gene_id: str,
    limit: int | None = None,
    page: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Get alleles and variants for a Drosophila gene from FlyBase via the Alliance of Genome Resources. Returns allele symbols, categories (insertion, point mutation, etc.), associated phenotypes, disease annotations, and variant details including genomic coordinates. Example: 'FB:FBgn0000490' (dpp) returns alleles like dpp[MI03752] with variant location and molecular consequence.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_id'] = gene_id
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if page is not None:
        _dnakit_arguments['page'] = page
    return _call('fly_base_get_gene_alleles', _dnakit_arguments, verbose=verbose)

def fly_base_get_gene_disease_models(*,
    gene_id: str,
    limit: int | None = None,
    page: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Get disease models involving a Drosophila gene from the Alliance of Genome Resources. Drosophila is a powerful model organism for human disease research. Returns genotype-based disease models with associated human disease names and Disease Ontology IDs. Example: 'FB:FBgn0000490' (dpp) has 946 disease model entries including models for ALS, Parkinson's disease, and other neurodegenerative conditions.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_id'] = gene_id
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if page is not None:
        _dnakit_arguments['page'] = page
    return _call('fly_base_get_gene_disease_models', _dnakit_arguments, verbose=verbose)

def fly_base_get_gene_expression(*,
    gene_id: str,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Get expression annotations for a Drosophila gene from the Alliance of Genome Resources. Returns per-annotation records of developmental stage (whenExpressed) and anatomical location (whereExpressed) for the gene, plus the total annotation count. Example: 'FB:FBgn0000490' (dpp) has 317 expression annotations such as 'embryonic stage 14' in 'abdominal transverse muscle cell'. Useful for understanding when and where a fly gene is expressed.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_id'] = gene_id
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('fly_base_get_gene_expression', _dnakit_arguments, verbose=verbose)

def fly_base_get_gene_interactions(*,
    gene_id: str,
    interaction_type: str | None = None,
    limit: int | None = None,
    page: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Get molecular or genetic interactions for a Drosophila gene from the Alliance of Genome Resources. Molecular interactions are physical protein-protein interactions. Genetic interactions are phenotype-based (suppression, enhancement, etc.). Returns interacting gene pairs with species and interaction type. Example: 'FB:FBgn0000490' (dpp) has 169 genetic interactions including genetic suppression with anchor.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_id'] = gene_id
    if interaction_type is not None:
        _dnakit_arguments['interaction_type'] = interaction_type
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if page is not None:
        _dnakit_arguments['page'] = page
    return _call('fly_base_get_gene_interactions', _dnakit_arguments, verbose=verbose)

def fly_base_get_gene_orthologs(*,
    gene_id: str,
    limit: int | None = None,
    page: int | None = None,
    stringency: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Get orthologs of a Drosophila gene across other species from the Alliance of Genome Resources. Returns orthologous genes in human, mouse, zebrafish, worm, yeast, and other model organisms with prediction method details and best-score indicators. Essential for translational research from fly models to human disease genes. Example: 'FB:FBgn0004647' (Notch) returns orthologs including human NOTCH1-4.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_id'] = gene_id
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if page is not None:
        _dnakit_arguments['page'] = page
    if stringency is not None:
        _dnakit_arguments['stringency'] = stringency
    return _call('fly_base_get_gene_orthologs', _dnakit_arguments, verbose=verbose)

def fly_mine_search(*,
    q: str,
    size: int | None = None,
    format: str | None = None,
    facet_Category: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Search FlyMine, the InterMine data warehouse for Drosophila melanogaster (fruit fly) genomics data from FlyBase. Searches across genes, proteins, publications, protein domains, gene ontology terms, and more. FlyMine integrates data from FlyBase, UniProt, Reactome, KEGG, TreeFam, and other model organism databases. Returns matching records with category facets and relevance scores. Use facet_Category parameter to filter by entity type (e.g., 'Gene', 'Protein', 'Publication', 'ProteinDomain', 'CDNAClone', 'GOTerm', 'InteractionExperiment').\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['q'] = q
    if size is not None:
        _dnakit_arguments['size'] = size
    if format is not None:
        _dnakit_arguments['format'] = format
    if facet_Category is not None:
        _dnakit_arguments['facet_Category'] = facet_Category
    return _call('fly_mine_search', _dnakit_arguments, verbose=verbose)

def four_dn_search_data(*,
    operation: str | None = None,
    query: str | None = None,
    item_type: str | None = None,
    file_type: str | None = None,
    assay_title: str | None = None,
    biosource_name: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search 4DN Data Portal for Hi-C chromatin conformation files and experiments. Access 350+ uniformly processed Hi-C contact files from multiple technologies (in situ Hi-C, Capture Hi-C, Micro-C XL). Search by assay type, cell type, or keywords. Returns file metadata with download URLs. Note: Downloads require free 4DN account (create at data.4dnucleome.org). Use for: finding Hi-C datasets, discovering 3D genome data, accessing TAD (Topologically Associating Domain) boundaries, A/B compartments, chromatin loop calls.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if query is not None:
        _dnakit_arguments['query'] = query
    if item_type is not None:
        _dnakit_arguments['item_type'] = item_type
    if file_type is not None:
        _dnakit_arguments['file_type'] = file_type
    if assay_title is not None:
        _dnakit_arguments['assay_title'] = assay_title
    if biosource_name is not None:
        _dnakit_arguments['biosource_name'] = biosource_name
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('four_dn_search_data', _dnakit_arguments, verbose=verbose)

def gbif_get_species(*,
    speciesKey: int,
    verbose: bool = True,
) -> _Object:
    'Get detailed species information from the Global Biodiversity Information Facility (GBIF) by its taxon key. Returns the full taxonomic record including scientific name, authorship, rank, full classification (kingdom through species), taxonomic status, basionym, vernacular name, and number of descendant taxa. Use GBIF_search_species first to find taxon keys. Example: key 2435099 returns Puma concolor (puma) with full classification from Animalia to Puma concolor.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['speciesKey'] = speciesKey
    return _call('gbif_get_species', _dnakit_arguments, verbose=verbose)

def gbif_search_species(*,
    query: str,
    limit: int | None = None,
    offset: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search species via GBIF species/search\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    return _call('gbif_search_species', _dnakit_arguments, verbose=verbose)

def gdc_get_clinical_data(*,
    project_id: str | None = None,
    primary_site: str | None = None,
    disease_type: str | None = None,
    vital_status: str | None = None,
    gender: str | None = None,
    size: int | None = None,
    offset: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get detailed clinical data for cancer cases from NCI GDC/TCGA. Returns demographics (gender, race, vital_status, age_at_index), diagnoses (primary_diagnosis, tumor_stage, age_at_diagnosis, days_to_last_follow_up), and treatments (therapeutic_agents, treatment_type). Filter by project, primary_site, disease_type, or vital_status.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if project_id is not None:
        _dnakit_arguments['project_id'] = project_id
    if primary_site is not None:
        _dnakit_arguments['primary_site'] = primary_site
    if disease_type is not None:
        _dnakit_arguments['disease_type'] = disease_type
    if vital_status is not None:
        _dnakit_arguments['vital_status'] = vital_status
    if gender is not None:
        _dnakit_arguments['gender'] = gender
    if size is not None:
        _dnakit_arguments['size'] = size
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    return _call('gdc_get_clinical_data', _dnakit_arguments, verbose=verbose)

def gdc_get_mutation_frequency(*,
    gene_symbol: str | None = None,
    gene: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get pan-cancer mutation frequency statistics for a gene across all TCGA projects. Returns overall and per-project mutation rates. Note: this tool is pan-cancer only and does not support filtering by cancer type.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if gene_symbol is not None:
        _dnakit_arguments['gene_symbol'] = gene_symbol
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    return _call('gdc_get_mutation_frequency', _dnakit_arguments, verbose=verbose)

def gdc_get_ssm_by_gene(*,
    gene_symbol: str,
    project_id: str | None = None,
    size: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get somatic mutations (SSMs) for a gene across TCGA/GDC projects\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_symbol'] = gene_symbol
    if project_id is not None:
        _dnakit_arguments['project_id'] = project_id
    if size is not None:
        _dnakit_arguments['size'] = size
    return _call('gdc_get_ssm_by_gene', _dnakit_arguments, verbose=verbose)

def gdc_get_survival(*,
    project_id: str,
    gene_symbol: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get Kaplan-Meier survival data for a GDC/TCGA cancer cohort. Returns time-to-event data with censoring status and survival estimates for each patient. Filter by project and optionally by gene mutation status. Use for overall survival analysis of TCGA cancer types.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['project_id'] = project_id
    if gene_symbol is not None:
        _dnakit_arguments['gene_symbol'] = gene_symbol
    return _call('gdc_get_survival', _dnakit_arguments, verbose=verbose)

def gdc_list_files(*,
    data_type: str | None = None,
    size: int | None = None,
    offset: int | None = None,
    verbose: bool = True,
) -> _Object:
    'List NCI GDC files via /files with optional data_type filter\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if data_type is not None:
        _dnakit_arguments['data_type'] = data_type
    if size is not None:
        _dnakit_arguments['size'] = size
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    return _call('gdc_list_files', _dnakit_arguments, verbose=verbose)

def gdc_list_projects(*,
    program: str | None = None,
    size: int | None = None,
    verbose: bool = True,
) -> _Object:
    'List GDC projects (TCGA, TARGET, etc.) with summary statistics\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if program is not None:
        _dnakit_arguments['program'] = program
    if size is not None:
        _dnakit_arguments['size'] = size
    return _call('gdc_list_projects', _dnakit_arguments, verbose=verbose)

def gdc_search_cases(*,
    project_id: str | None = None,
    size: int | None = None,
    offset: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search NCI GDC cases via /cases\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if project_id is not None:
        _dnakit_arguments['project_id'] = project_id
    if size is not None:
        _dnakit_arguments['size'] = size
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    return _call('gdc_search_cases', _dnakit_arguments, verbose=verbose)

def geo_search_atacseq_datasets(*,
    query: str | None | None = None,
    organism: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search NCBI GEO for ATAC-seq chromatin accessibility datasets by query keyword, tissue, or cell type. Returns dataset accessions (GSE IDs), titles, summaries, organism, and sample counts. ATAC-seq measures open chromatin regions and is used to map regulatory elements, enhancers, and TF binding sites genome-wide.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if query is not None:
        _dnakit_arguments['query'] = query
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('geo_search_atacseq_datasets', _dnakit_arguments, verbose=verbose)

def geo_search_chipseq_datasets(*,
    query: str,
    organism: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search NCBI GEO for ChIP-seq (Chromatin Immunoprecipitation followed by sequencing) datasets. ChIP-seq profiles genome-wide binding of proteins (transcription factors, histones) to DNA. Returns dataset accessions, titles, summaries, and sample counts. Use this to find ChIP-seq studies for specific transcription factors, histone marks, or in specific cell types/diseases.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('geo_search_chipseq_datasets', _dnakit_arguments, verbose=verbose)

def geo_search_rnaseq_datasets(*,
    query: str | None | None = None,
    organism: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search NCBI GEO for RNA-seq gene expression datasets by query keyword, disease, tissue, or gene of interest. Returns dataset accessions (GSE IDs), titles, summaries, organism, and sample counts. Use this to find transcriptomics experiments relevant to a biological question. Covers bulk RNA-seq and total RNA-seq studies across all tissues and diseases.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if query is not None:
        _dnakit_arguments['query'] = query
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('geo_search_rnaseq_datasets', _dnakit_arguments, verbose=verbose)

def go_get_annotations_for_gene(*,
    gene_id: str,
    rows: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Finds GO annotations for a specific gene/protein using GOlr search. Returns one row per annotation, so a gene with several evidence codes for the same GO term appears several times. This returns a PAGE, not the complete set: heavily annotated genes have far more annotations than the default page holds (measured: 'TP53' matches 1812 annotation rows), and the response carries no total, so a short list is NOT evidence the gene has no further annotations -- raise 'rows' if completeness matters.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_id'] = gene_id
    if rows is not None:
        _dnakit_arguments['rows'] = rows
    return _call('go_get_annotations_for_gene', _dnakit_arguments, verbose=verbose)

def go_search_terms(*,
    query: str,
    verbose: bool = True,
) -> _Object:
    'Searches for Gene Ontology (GO) terms by a keyword using the GOlr search engine. Returns GO terms and related biological entities.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    return _call('go_search_terms', _dnakit_arguments, verbose=verbose)

def gtdb_get_species(*,
    operation: str | None = None,
    species: str,
    max_genomes: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get detailed information about a GTDB species cluster, including all genomes assigned to that species. A GTDB species cluster groups genomes by average nucleotide identity (ANI). Returns the full list of genomes with their accessions, NCBI organism names, NCBI taxonomy, and whether they are the GTDB species representative or NCBI type material. Useful for exploring species diversity, finding representative genomes, and comparing GTDB vs NCBI classifications.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['species'] = species
    if max_genomes is not None:
        _dnakit_arguments['max_genomes'] = max_genomes
    return _call('gtdb_get_species', _dnakit_arguments, verbose=verbose)

def gtdb_get_taxon_info(*,
    operation: str | None = None,
    taxon: str,
    verbose: bool = True,
) -> _Object:
    'Get information about a GTDB taxon including its rank, total genome count, higher-rank lineage, and release history. GTDB uses rank prefixes: d__ (domain), p__ (phylum), c__ (class), o__ (order), f__ (family), g__ (genus), s__ (species). For example, f__Lachnospiraceae returns its rank (Family), genome count, and lineage (Bacteria > Bacillota > Clostridia > Lachnospirales). Also returns the full lineage breakdown from domain to species. Useful for understanding where a taxon fits in the GTDB tree and how many genomes support it.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['taxon'] = taxon
    return _call('gtdb_get_taxon_info', _dnakit_arguments, verbose=verbose)

def gtdb_search_genomes(*,
    operation: str | None = None,
    query: str,
    page: int | None = None,
    items_per_page: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Search for prokaryotic genomes in GTDB by organism name. A multi-word query (e.g. a binomial like 'Mycobacterium tuberculosis') is first matched precisely, requiring EVERY term of the query to appear in the genome's name or taxonomy; 'match_type' is then 'exact'. If nothing matches the whole query, the tool falls back to GTDB's loose search, which matches ANY term of the query and therefore returns a broader taxon (often just the genus) or a different species; that case is reported as 'match_type': 'broad' with a top-level 'note', and 'total_results' then counts the broader match rather than the queried organism -- always read 'match_type' and 'total_results_scope' before quoting a count or trusting a result. Results also match against the whole taxonomic clade rather than organism-name text alone, so exact/prefix name matches are ranked first within a page. Returns genome accessions with both NCBI and GTDB taxonomy classifications, allowing comparison between the two taxonomic systems. Each result shows the genome accession (GCA/GCF), NCBI organism name, NCBI taxonomy, GTDB taxonomy, and whether the genome is a GTDB species representative. Supports pagination. Useful for finding genome accessions, exploring taxonomic reclassifications, and linking genomes to the GTDB taxonomy.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['query'] = query
    if page is not None:
        _dnakit_arguments['page'] = page
    if items_per_page is not None:
        _dnakit_arguments['items_per_page'] = items_per_page
    return _call('gtdb_search_genomes', _dnakit_arguments, verbose=verbose)

def gtdb_search_taxon(*,
    operation: str | None = None,
    query: str,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search the Genome Taxonomy Database (GTDB) for prokaryotic taxa by partial name. GTDB provides a standardized genome-based taxonomy for Bacteria and Archaea, resolving polyphyletic groups in NCBI taxonomy. Returns matching taxon names with GTDB rank prefixes (d__=domain, p__=phylum, c__=class, o__=order, f__=family, g__=genus, s__=species). Search is case-insensitive and matches partial names. Use this for looking up standardized prokaryotic taxonomy names, finding correct GTDB identifiers, or exploring taxonomic groups like Lachnospiraceae, Enterobacteriaceae, or Clostridia.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['query'] = query
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('gtdb_search_taxon', _dnakit_arguments, verbose=verbose)

def gtex_calculate_eqtl(*,
    operation: str | None = None,
    gencode_id: str,
    variant_id: str,
    tissue_site_detail_id: str,
    dataset_id: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Calculate custom eQTL for any gene-variant pair in any tissue. Dynamically calculates gene-variant association statistics including p-value, normalized effect size (NES), t-statistic, and genotype data. Works for both significant and non-significant associations. Returns detailed statistics even if association not in precomputed results. Use for: testing specific hypotheses, custom variant-gene pairs, dynamic eQTL analysis, validating findings. Requires all three: gencode_id, variant_id, and tissue_site_detail_id.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['gencode_id'] = gencode_id
    _dnakit_arguments['variant_id'] = variant_id
    _dnakit_arguments['tissue_site_detail_id'] = tissue_site_detail_id
    if dataset_id is not None:
        _dnakit_arguments['dataset_id'] = dataset_id
    return _call('gtex_calculate_eqtl', _dnakit_arguments, verbose=verbose)

def gtex_get_expression_summary(*,
    gene_symbol: str | None = None,
    ensembl_gene_id: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get GTEx expression summary for a gene via /expression/geneExpression\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if gene_symbol is not None:
        _dnakit_arguments['gene_symbol'] = gene_symbol
    if ensembl_gene_id is not None:
        _dnakit_arguments['ensembl_gene_id'] = ensembl_gene_id
    return _call('gtex_get_expression_summary', _dnakit_arguments, verbose=verbose)

def gtex_get_median_gene_expression(*,
    operation: str | None = None,
    gencode_id: str | list[_Object] | None = None,
    tissue_site_detail_id: list[str] | None = None,
    dataset_id: str | None = None,
    page: int | None = None,
    items_per_page: int | None = None,
    gene_symbol: str | None = None,
    gene: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Get median gene expression levels across GTEx tissues. Returns median expression in TPM (Transcripts Per Million) for one or more genes across 54 tissue sites. Filter by specific tissues or get expression across all tissues. Defaults to GTEx V8 (2020) -- confirmed live no V11 dataset is queryable via this API (dataset_id accepts gtex_v7/gtex_v8/gtex_v10). Set dataset_id='gtex_v10' for the newer release; gene IDs are resolved to each dataset's own GENCODE version automatically, so V10 returns full results. Use for: tissue-specific expression profiling, identifying highly expressed genes, comparing expression patterns across tissues. Example: Get brain vs liver expression for TP53 gene.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if gencode_id is not None:
        _dnakit_arguments['gencode_id'] = gencode_id
    if tissue_site_detail_id is not None:
        _dnakit_arguments['tissue_site_detail_id'] = tissue_site_detail_id
    if dataset_id is not None:
        _dnakit_arguments['dataset_id'] = dataset_id
    if page is not None:
        _dnakit_arguments['page'] = page
    if items_per_page is not None:
        _dnakit_arguments['items_per_page'] = items_per_page
    if gene_symbol is not None:
        _dnakit_arguments['gene_symbol'] = gene_symbol
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    return _call('gtex_get_median_gene_expression', _dnakit_arguments, verbose=verbose)

def gtex_get_multi_tissue_eqtls(*,
    operation: str | None = None,
    gencode_id: str,
    variant_id: str | None = None,
    dataset_id: str | None = None,
    page: int | None = None,
    items_per_page: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get multi-tissue eQTL meta-analysis results (Metasoft). Returns m-values (posterior probability of eQTL effect in each tissue), normalized effect sizes, p-values, and standard errors across all tissues for a gene-variant pair. Shows tissue-sharing patterns of eQTL effects. Use for: understanding tissue-specificity of eQTLs, finding shared vs tissue-specific effects, cross-tissue meta-analysis, regulatory mechanism studies. Requires gencode_id, optionally filter by variant_id.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['gencode_id'] = gencode_id
    if variant_id is not None:
        _dnakit_arguments['variant_id'] = variant_id
    if dataset_id is not None:
        _dnakit_arguments['dataset_id'] = dataset_id
    if page is not None:
        _dnakit_arguments['page'] = page
    if items_per_page is not None:
        _dnakit_arguments['items_per_page'] = items_per_page
    return _call('gtex_get_multi_tissue_eqtls', _dnakit_arguments, verbose=verbose)

def gtex_get_single_tissue_eqtls(*,
    operation: str | None = None,
    gencode_id: list[str] | None = None,
    variant_id: list[str] | None = None,
    tissue_site_detail_id: list[str] | None = None,
    dataset_id: str | None = None,
    page: int | None = None,
    items_per_page: int | None = None,
    gene_symbol: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get significant single-tissue eQTL associations. Returns precomputed gene-variant associations with p-values, effect sizes (NES), and genomic positions. Query by gene ID, variant ID, tissue, or combination. Results show which genetic variants affect gene expression in specific tissues. Use for: finding eQTLs for candidate genes, looking up variant effects, tissue-specific regulatory analysis, GWAS follow-up. Must provide at least one filter: gencode_id, variant_id, or tissue_site_detail_id.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if gencode_id is not None:
        _dnakit_arguments['gencode_id'] = gencode_id
    if variant_id is not None:
        _dnakit_arguments['variant_id'] = variant_id
    if tissue_site_detail_id is not None:
        _dnakit_arguments['tissue_site_detail_id'] = tissue_site_detail_id
    if dataset_id is not None:
        _dnakit_arguments['dataset_id'] = dataset_id
    if page is not None:
        _dnakit_arguments['page'] = page
    if items_per_page is not None:
        _dnakit_arguments['items_per_page'] = items_per_page
    if gene_symbol is not None:
        _dnakit_arguments['gene_symbol'] = gene_symbol
    return _call('gtex_get_single_tissue_eqtls', _dnakit_arguments, verbose=verbose)

def gtex_get_tissue_sites(*,
    operation: str | None = None,
    dataset_id: str | None = None,
    page: int | None = None,
    items_per_page: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get all available GTEx tissue sites with metadata. Returns tissue IDs, names, sample counts, eGene counts, and color codes for 54 non-diseased tissue sites. Includes RNA-seq sample summaries with donor age and sex statistics. Use for: discovering available tissues, finding tissue IDs for queries, checking sample sizes, understanding tissue coverage. Essential first step before querying tissue-specific data.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if dataset_id is not None:
        _dnakit_arguments['dataset_id'] = dataset_id
    if page is not None:
        _dnakit_arguments['page'] = page
    if items_per_page is not None:
        _dnakit_arguments['items_per_page'] = items_per_page
    return _call('gtex_get_tissue_sites', _dnakit_arguments, verbose=verbose)

def gtex_query_eqtl(*,
    gene_symbol: str | None = None,
    ensembl_gene_id: str | None = None,
    page: int | None = None,
    size: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Query GTEx single-tissue eQTL via /association/singleTissueEqtl\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if gene_symbol is not None:
        _dnakit_arguments['gene_symbol'] = gene_symbol
    if ensembl_gene_id is not None:
        _dnakit_arguments['ensembl_gene_id'] = ensembl_gene_id
    if page is not None:
        _dnakit_arguments['page'] = page
    if size is not None:
        _dnakit_arguments['size'] = size
    return _call('gtex_query_eqtl', _dnakit_arguments, verbose=verbose)

def gwas_search_associations_by_gene(*,
    gene_name: str,
    size: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search GWAS Catalog associations by gene name (returns strongest risk allele and p-value fields).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_name'] = gene_name
    if size is not None:
        _dnakit_arguments['size'] = size
    return _call('gwas_search_associations_by_gene', _dnakit_arguments, verbose=verbose)

def gen_cc_get_classifications(*,
    operation: str | None = None,
    gene_symbol: str | None = None,
    disease: str | None = None,
    submitter: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Get summary statistics of gene-disease validity classifications from GenCC. Returns counts per classification level (Definitive, Strong, Moderate, Limited, etc.). Optionally filter by gene_symbol (e.g., 'GBA'), disease name substring, or submitter organization. Without filters, returns global stats across all 40,000+ curated submissions.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if gene_symbol is not None:
        _dnakit_arguments['gene_symbol'] = gene_symbol
    if disease is not None:
        _dnakit_arguments['disease'] = disease
    if submitter is not None:
        _dnakit_arguments['submitter'] = submitter
    return _call('gen_cc_get_classifications', _dnakit_arguments, verbose=verbose)

def gen_cc_search_disease(*,
    operation: str | None = None,
    disease: str | None = None,
    disease_name: str | None = None,
    classification: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Find genes with validity evidence for a disease from GenCC (Gene Curation Coalition). Search by disease name or MONDO/OMIM identifier. Returns gene-disease validity classifications from multiple expert curators. Useful for identifying established gene-disease relationships for clinical interpretation.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if disease is not None:
        _dnakit_arguments['disease'] = disease
    if disease_name is not None:
        _dnakit_arguments['disease_name'] = disease_name
    if classification is not None:
        _dnakit_arguments['classification'] = classification
    return _call('gen_cc_search_disease', _dnakit_arguments, verbose=verbose)

def gen_cc_search_gene(*,
    operation: str | None = None,
    gene_symbol: str,
    classification: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get gene-disease validity classifications for a gene from GenCC (Gene Curation Coalition). Returns curated evidence levels (Definitive, Strong, Moderate, Limited, Disputed, Refuted) from multiple expert panels including ClinGen, Ambry Genetics, Genomics England PanelApp, and Orphanet. Essential for assessing clinical validity of gene-disease associations in genetic testing.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['gene_symbol'] = gene_symbol
    if classification is not None:
        _dnakit_arguments['classification'] = classification
    return _call('gen_cc_search_gene', _dnakit_arguments, verbose=verbose)

def genebass_gene_burden_phewas(*,
    gene: str,
    burden_set: str | None = None,
    limit: int | None = None,
    max_pval: float | None | None = None,
    verbose: bool = True,
) -> _Object:
    'Exome-wide gene-level burden PheWAS from Genebass (UK Biobank ~394K exomes, SAIGE-GENE+). Given one gene, returns its rare-variant burden association across ~4,500 phenotypes (GRCh38), sorted by p-value. Aggregates rare coding variants into a single gene-level test, complementing single-variant PheWAS. Accepts an Ensembl gene ID (ENSG...) or a gene symbol (auto-resolved via Ensembl). Choose a burden_set: pLoF (predicted loss-of-function, default), missense|LC, or synonymous.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene'] = gene
    if burden_set is not None:
        _dnakit_arguments['burden_set'] = burden_set
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if max_pval is not None:
        _dnakit_arguments['max_pval'] = max_pval
    return _call('genebass_gene_burden_phewas', _dnakit_arguments, verbose=verbose)

def hgnc_fetch_gene_by_symbol(*,
    symbol: str,
    verbose: bool = True,
) -> _Object:
    "Fetch detailed gene information from HGNC (HUGO Gene Nomenclature Committee) by official gene symbol. Returns comprehensive gene data including HGNC ID, full name, chromosomal location, cross-references to Ensembl, NCBI Entrez, UniProt, OMIM, and RefSeq. HGNC is the authoritative source for human gene nomenclature. Example: 'TP53' returns tumor protein p53 with all its identifiers. Retired symbols are resolved automatically: if the symbol is not a current approved symbol, it is looked up as a previous symbol and then as an alias symbol, and the current gene's full record is returned with the substitution disclosed in metadata ('resolved_from', 'resolved_symbol', 'resolution_relation', 'symbol_resolution_note'). Example: 'CYP2E' returns CYP2E1 (renamed in 2002), 'DIA4' returns NQO1. If the retired symbol is ambiguous (shared by several genes) an error naming every candidate is returned instead of an arbitrary pick; if the symbol is unknown to HGNC entirely, metadata carries 'no_results_note' pointing at HGNC_search_genes.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['symbol'] = symbol
    return _call('hgnc_fetch_gene_by_symbol', _dnakit_arguments, verbose=verbose)

def hmdb_get_diseases(*,
    operation: str | None = None,
    hmdb_id: str | None = None,
    compound_name: str | None = None,
    pubchem_cid: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Alias for Metabolite_get_diseases. Get disease associations for a metabolite via CTD database.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if hmdb_id is not None:
        _dnakit_arguments['hmdb_id'] = hmdb_id
    if compound_name is not None:
        _dnakit_arguments['compound_name'] = compound_name
    if pubchem_cid is not None:
        _dnakit_arguments['pubchem_cid'] = pubchem_cid
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('hmdb_get_diseases', _dnakit_arguments, verbose=verbose)

def hmdb_get_metabolite(*,
    operation: str | None = None,
    hmdb_id: str | None = None,
    compound_name: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Alias for Metabolite_get_info. Get compound information by HMDB ID, compound name, or PubChem CID.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if hmdb_id is not None:
        _dnakit_arguments['hmdb_id'] = hmdb_id
    if compound_name is not None:
        _dnakit_arguments['compound_name'] = compound_name
    return _call('hmdb_get_metabolite', _dnakit_arguments, verbose=verbose)

def hmdb_search(*,
    operation: str | None = None,
    query: str,
    verbose: bool = True,
) -> _Object:
    "Alias for Metabolite_search. Search for metabolites by compound name or molecular formula via PubChem. Accepts common names, synonyms (e.g., 'glucosylceramide', 'glucose', 'ATP').\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['query'] = query
    return _call('hmdb_search', _dnakit_arguments, verbose=verbose)

def hpa_get_cancer_prognostics_by_gene(*,
    ensembl_id: str | None = None,
    gene_name: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Retrieve prognostic value of a gene across various cancer types, indicating if its expression level correlates with patient survival outcomes. Reports both the significant results ('prognostic_summary') and the cancer types that were assessed but found NOT prognostic ('non_prognostic_summary'), plus 'cancers_assessed_count' as the denominator, so that a negative result can be distinguished from missing data.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if ensembl_id is not None:
        _dnakit_arguments['ensembl_id'] = ensembl_id
    if gene_name is not None:
        _dnakit_arguments['gene_name'] = gene_name
    return _call('hpa_get_cancer_prognostics_by_gene', _dnakit_arguments, verbose=verbose)

def hpa_get_protein_interactions_by_gene(*,
    gene_name: str,
    verbose: bool = True,
) -> _Object:
    "DEPRECATED / non-functional: the Human Protein Atlas search API no longer exposes protein-protein interaction (ppi) data, so this tool returns an error. For PPI partners use STRING_get_protein_interactions (by gene symbol, parameter 'protein_ids') or EBIProteins_get_interactions (by UniProt accession) instead.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_name'] = gene_name
    return _call('hpa_get_protein_interactions_by_gene', _dnakit_arguments, verbose=verbose)

def hpa_get_rna_expression_by_source(*,
    gene_name: str,
    source_type: str,
    source_name: str,
    verbose: bool = True,
) -> _Object:
    "Get RNA expression for a gene in a specific biological source, read from HPA's dedicated per-source download columns (t_RNA_*, blood_RNA_*, brain_RNA_*, sc_RNA_*). Supports source_type 'tissue', 'blood', 'brain' and 'single_cell'. Units differ by source: single_cell values are nCPM, all others nTPM -- the returned expression_unit is read from HPA's own column header, never assumed.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_name'] = gene_name
    _dnakit_arguments['source_type'] = source_type
    _dnakit_arguments['source_name'] = source_name
    return _call('hpa_get_rna_expression_by_source', _dnakit_arguments, verbose=verbose)

def hpa_get_rna_expression_in_specific_tissues(*,
    ensembl_id: str | None = None,
    tissue_names: list[str] | None = None,
    tissues: list[str] | None = None,
    gene_name: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Query RNA expression levels (nTPM) for a specific gene in one or more user-specified tissues. Reads HPA's full per-tissue nTPM panel (the 'Tissue RNA - <tissue> [nTPM]' columns), so any of HPA's 51 published consensus tissues and brain regions can be queried, e.g. 'cerebral cortex', 'basal ganglia', 'midbrain', 'liver'. Common names are mapped to HPA's actual column names, including the ones HPA suffixes with '1' ('skin' -> 'Skin 1', 'stomach' -> 'Stomach 1', 'endometrium' -> 'Endometrium 1'); when a name resolves to a differently-named HPA column the response says so in a 'note'. A tissue reported as 'Not found' means HPA publishes no RNA column for that name (it is a naming problem or a tissue HPA covers only with protein/IHC data), not that the gene lacks expression there.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if ensembl_id is not None:
        _dnakit_arguments['ensembl_id'] = ensembl_id
    if tissue_names is not None:
        _dnakit_arguments['tissue_names'] = tissue_names
    if tissues is not None:
        _dnakit_arguments['tissues'] = tissues
    if gene_name is not None:
        _dnakit_arguments['gene_name'] = gene_name
    return _call('hpa_get_rna_expression_in_specific_tissues', _dnakit_arguments, verbose=verbose)

def hpa_get_subcellular_location(*,
    gene_name: str,
    verbose: bool = True,
) -> _Object:
    'Get annotated subcellular locations for a protein using optimized columns parameter. Retrieves both main and additional subcellular locations efficiently.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_name'] = gene_name
    return _call('hpa_get_subcellular_location', _dnakit_arguments, verbose=verbose)

def hpa_search_genes_by_query(*,
    search_query: str,
    max_results: int | None | None = None,
    verbose: bool = True,
) -> _Object:
    'Search for matching genes by gene name, keywords, or cell line names and return Ensembl ID list. This is the entry tool for many HPA query workflows.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['search_query'] = search_query
    if max_results is not None:
        _dnakit_arguments['max_results'] = max_results
    return _call('hpa_search_genes_by_query', _dnakit_arguments, verbose=verbose)

def hpo_get_diseases_by_phenotype(*,
    term_id: str | None = None,
    limit: int | None | None = None,
    offset: int | None | None = None,
    hp_id: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Get the diseases that present with an HPO phenotype term, from the JAX HPO network annotations. Returns disease IDs (ORPHA/OMIM/DECIPHER), names, and the cross-referenced MONDO ID. Example: HP:0001250 (Seizure) is a feature of thousands of diseases. Use to build a differential diagnosis from a set of observed phenotypes. A phenotype term can be linked to thousands of diseases but 'limit' caps a single response at 500 -- check metadata.has_more in the response and pass 'offset' to page through the rest; a well-known disease for a common phenotype may not appear in the first page.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if term_id is not None:
        _dnakit_arguments['term_id'] = term_id
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    if hp_id is not None:
        _dnakit_arguments['hp_id'] = hp_id
    return _call('hpo_get_diseases_by_phenotype', _dnakit_arguments, verbose=verbose)

def hpo_get_genes_by_phenotype(*,
    term_id: str,
    limit: int | None | None = None,
    offset: int | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Get the genes associated with an HPO phenotype term (genes in which variants cause this phenotype), from the JAX HPO network annotations. Returns NCBI gene IDs and symbols. Example: HP:0001250 (Seizure) is linked to genes such as RELN, SCN1A, and ~2000 others. Use to go from an observed patient phenotype to candidate disease genes for variant prioritization. A phenotype term can be linked to thousands of genes but 'limit' caps a single response at 500 -- check metadata.has_more in the response and pass 'offset' to page through the rest; a well-known gene for a common phenotype may not appear in the first page.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['term_id'] = term_id
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    return _call('hpo_get_genes_by_phenotype', _dnakit_arguments, verbose=verbose)

def hpo_get_term(*,
    term_id: str,
    verbose: bool = True,
) -> _Object:
    "Get detailed information about a Human Phenotype Ontology (HPO) term by its ID. HPO provides standardized descriptions of phenotypic abnormalities in human disease with over 18,000 terms. Returns the term name, definition, synonyms, cross-references to UMLS/SNOMED, descendant count, and translations. Example: HP:0001250 (Seizure) returns definition, 346 descendants, synonyms including 'Epileptic seizure' and 'Epilepsy', and cross-references to UMLS and SNOMED-CT.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['term_id'] = term_id
    return _call('hpo_get_term', _dnakit_arguments, verbose=verbose)

def hpo_search_terms(*,
    query: str,
    max_results: int | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Search for Human Phenotype Ontology (HPO) terms by keyword. HPO contains over 18,000 standardized phenotype terms used in clinical genetics and rare disease research. Returns matching terms with definitions, synonym counts, and descendant counts. Useful for finding the correct HPO ID for a clinical phenotype before querying disease or gene associations. Example: searching 'seizure' returns terms like 'Bilateral tonic-clonic seizure' (HP:0007334), 'Generalized-onset seizure' (HP:0002197), 'Focal-onset seizure' (HP:0007359).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if max_results is not None:
        _dnakit_arguments['max_results'] = max_results
    return _call('hpo_search_terms', _dnakit_arguments, verbose=verbose)

def harmonizome_get_gene(*,
    gene_symbol: str,
    verbose: bool = True,
) -> _Object:
    "Get comprehensive gene information from Harmonizome (Ma'ayan Lab), which integrates data from 100+ genomics datasets covering gene expression, protein interactions, pathways, diseases, drug targets, and more. Returns gene symbol, name, NCBI Entrez ID, description, synonyms, and linked protein identifiers. Useful as a quick gene summary with curated descriptions from RefSeq. Example: TP53 returns description about tumor suppressor protein, NCBI ID 7157, synonyms [TRP53, LFS1, P53, BMFS5, BCC7], and linked UniProt protein P53_HUMAN.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_symbol'] = gene_symbol
    return _call('harmonizome_get_gene', _dnakit_arguments, verbose=verbose)

def igsr_list_data_collections(*,
    operation: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'List available data collections from the International Genome Sample Resource (IGSR). Data collections include the original 1000 Genomes Project phases, the 30x high-coverage resequencing, Human Genome Structural Variation Consortium (HGSVC), Human Pangenome Reference Consortium (HPRC), and other related studies. Returns collection codes, titles, sample/population counts, and available data types (sequence, alignment, variants). Use this to discover what genomics datasets are available through IGSR.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('igsr_list_data_collections', _dnakit_arguments, verbose=verbose)

def igsr_search_populations(*,
    operation: str | None = None,
    superpopulation: str | None | None = None,
    query: str | None | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search 1000 Genomes Project populations from the International Genome Sample Resource (IGSR). The 1000 Genomes Project characterized genetic variation across 26 populations grouped into 5 superpopulations: AFR (African), AMR (Admixed American), EAS (East Asian), EUR (European), SAS (South Asian). Returns population codes, names, geographic coordinates, sample counts, and superpopulation assignments. Use this to find populations for allele frequency comparisons, population stratification analysis, or ancestry-specific variant studies.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if superpopulation is not None:
        _dnakit_arguments['superpopulation'] = superpopulation
    if query is not None:
        _dnakit_arguments['query'] = query
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('igsr_search_populations', _dnakit_arguments, verbose=verbose)

def igsr_search_samples(*,
    operation: str | None = None,
    population: str | None | None = None,
    data_collection: str | None | None = None,
    sample_name: str | None | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search samples from the 1000 Genomes Project by population, data collection, or sample name. Returns sample names (e.g., HG00096, NA12878), sex, BioSample IDs, population assignments, and available data collections. Use this to find samples from specific populations for variant analysis, to identify samples with particular sequencing data types (WGS, exome, RNA-seq), or to look up specific well-characterized samples like NA12878 (reference genome sample).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if population is not None:
        _dnakit_arguments['population'] = population
    if data_collection is not None:
        _dnakit_arguments['data_collection'] = data_collection
    if sample_name is not None:
        _dnakit_arguments['sample_name'] = sample_name
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('igsr_search_samples', _dnakit_arguments, verbose=verbose)

def imgt_get_gene_info(*,
    operation: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get information about IMGT databases, gene nomenclature, and analysis tools. Returns descriptions of IMGT resources (LIGM-DB, GENE-DB, V-QUEST) and gene naming conventions. Good starting point for immunogenetics research.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    return _call('imgt_get_gene_info', _dnakit_arguments, verbose=verbose)

def imgt_search_genes(*,
    operation: str | None = None,
    query: str | None = None,
    gene_type: str | None = None,
    species: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Search IMGT for immunoglobulin/TCR germline genes. Filter by gene type (IGHV, IGKV, TRAV, etc.) and species. Returns search URLs and gene type information for navigating IMGT GENE-DB.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if query is not None:
        _dnakit_arguments['query'] = query
    if gene_type is not None:
        _dnakit_arguments['gene_type'] = gene_type
    if species is not None:
        _dnakit_arguments['species'] = species
    return _call('imgt_search_genes', _dnakit_arguments, verbose=verbose)

def imm_port_search_studies(*,
    query: str,
    condition_or_disease: str | None = None,
    assay_method: str | None = None,
    research_focus: str | None = None,
    species: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search the ImmPort immunology database for studies by keyword, disease condition, assay method, or research focus. ImmPort is a NIAID-funded repository with 900+ open studies covering vaccine trials, flow cytometry, ELISPOT, RNA-seq, and clinical immunology data. Assay methods include: Flow Cytometry, ELISPOT, ELISA, RNA sequencing, Transcription profiling by array, Hemagglutination Inhibition, Luminex xMAP, Q-PCR, CyTOF, and more. Research focuses include: Vaccine Response, Infection Response, Immune Response, Transplantation, Autoimmune, Allergy, Preterm Birth. The species filter is applied server-side; condition_or_disease, assay_method, and research_focus filters are applied client-side.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if condition_or_disease is not None:
        _dnakit_arguments['condition_or_disease'] = condition_or_disease
    if assay_method is not None:
        _dnakit_arguments['assay_method'] = assay_method
    if research_focus is not None:
        _dnakit_arguments['research_focus'] = research_focus
    if species is not None:
        _dnakit_arguments['species'] = species
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('imm_port_search_studies', _dnakit_arguments, verbose=verbose)

def inter_pro_get_entries_for_protein(*,
    accession: str,
    verbose: bool = True,
) -> _Object:
    'Get all InterPro domain and family entries annotated on a specific protein. Performs a reverse lookup from a UniProt protein accession to all InterPro entries (domains, families, homologous superfamilies, conserved sites) that have been detected on that protein. Example: P04637 (TP53) has 9 InterPro entries including IPR002117 (p53 family), IPR011615 (p53 DNA-binding domain), IPR010991 (p53 tetramerisation domain), IPR013872 (p53 transactivation domain).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['accession'] = accession
    return _call('inter_pro_get_entries_for_protein', _dnakit_arguments, verbose=verbose)

def inter_pro_get_protein_domains(*,
    protein_id: str | None = None,
    accession: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get protein domain annotations from InterPro database using UniProt protein ID. Returns domain families, signatures, and functional annotations.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if protein_id is not None:
        _dnakit_arguments['protein_id'] = protein_id
    if accession is not None:
        _dnakit_arguments['accession'] = accession
    return _call('inter_pro_get_protein_domains', _dnakit_arguments, verbose=verbose)

def jaspar_get_transcription_factors(*,
    collection: str | None = None,
    page: int | None = None,
    page_size: int | None = None,
    verbose: bool = True,
) -> _Object:
    'List all JASPAR transcription factor binding profiles. Returns all matrices; does NOT support filtering by TF name. To search by name, use jaspar_search_matrices with {"name": "FOXA1"}. Supports pagination via page/page_size.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if collection is not None:
        _dnakit_arguments['collection'] = collection
    if page is not None:
        _dnakit_arguments['page'] = page
    if page_size is not None:
        _dnakit_arguments['page_size'] = page_size
    return _call('jaspar_get_transcription_factors', _dnakit_arguments, verbose=verbose)

def kegg_get_gene_pathways(*,
    gene_id: str | None = None,
    gene_symbol: str | None = None,
    gene: str | None = None,
    organism: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get all KEGG pathways that a gene participates in. Returns pathway IDs and names for a given KEGG gene identifier. Example: hsa:7157 (human TP53) participates in 51 pathways including p53 signaling pathway (hsa04115), apoptosis (hsa04210), cell cycle (hsa04110), and 48 disease pathways (various cancers, viral infections). Requires KEGG organism:id format.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if gene_id is not None:
        _dnakit_arguments['gene_id'] = gene_id
    if gene_symbol is not None:
        _dnakit_arguments['gene_symbol'] = gene_symbol
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    return _call('kegg_get_gene_pathways', _dnakit_arguments, verbose=verbose)

def kegg_get_pathway_genes(*,
    pathway_id: str,
    verbose: bool = True,
) -> _Object:
    'Get all genes in a KEGG pathway. Returns gene IDs for all member genes of a specific pathway. Example: hsa04115 (p53 signaling pathway) contains genes including CDK2 (hsa:1017), CDK4 (hsa:1019), CDK6 (hsa:1021), CDKN1A (hsa:1026), CDKN2A (hsa:1029), and many more involved in p53-mediated cell cycle arrest, apoptosis, and senescence.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['pathway_id'] = pathway_id
    return _call('kegg_get_pathway_genes', _dnakit_arguments, verbose=verbose)

def ldlink_get_proxies(*,
    variant: str,
    population: str | None | None = None,
    r2_threshold: float | None | None = None,
    genome_build: str | None | None = None,
    limit: int | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Find variants in linkage disequilibrium (LD proxies) with a query SNP, population-specific, via NIH LDlink's LDproxy module. Essential for GWAS interpretation: a lead SNP is often just a tag for the true causal variant, so high-R2 proxies (especially in coding/regulatory regions) are better mechanistic candidates. Returns proxy rsIDs with R2, coordinates, alleles, MAF, and distance. Requires a free LDlink token (LDLINK_TOKEN env var; register at ldlink.nih.gov).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['variant'] = variant
    if population is not None:
        _dnakit_arguments['population'] = population
    if r2_threshold is not None:
        _dnakit_arguments['r2_threshold'] = r2_threshold
    if genome_build is not None:
        _dnakit_arguments['genome_build'] = genome_build
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('ldlink_get_proxies', _dnakit_arguments, verbose=verbose)

def marrvel_get_gene(*,
    symbol: str | None | None = None,
    gene_symbol: str | None | None = None,
    gene: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Get aggregated identity and annotation for a human gene by symbol from MARRVEL: cross-references (HGNC, Ensembl, Entrez, UniProt), chromosome/cytogenetic location, gene type, aliases, previous symbols, and the Entrez gene summary. Use as a one-call gene-identity resolver in rare-disease / Mendelian variant triage workflows. Does NOT include an OMIM ID -- MARRVEL's own gene-lookup endpoint returns an unreliable placeholder value for it (confirmed live), use MARRVEL_get_omim_phenotypes for OMIM associations instead. MARRVEL indexes only current approved HGNC symbols, so a retired symbol or alias (e.g. 'RP20'/'LCA2' for RPE65) is resolved once via HGNC and re-queried; when that happens the substitution is disclosed in metadata.resolved_symbol and metadata.resolved_note, and metadata.query_symbol keeps the symbol you asked for. No API key required.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if symbol is not None:
        _dnakit_arguments['symbol'] = symbol
    if gene_symbol is not None:
        _dnakit_arguments['gene_symbol'] = gene_symbol
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    return _call('marrvel_get_gene', _dnakit_arguments, verbose=verbose)

def marrvel_get_omim_phenotypes(*,
    symbol: str | None | None = None,
    gene_symbol: str | None | None = None,
    gene: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Get OMIM phenotype/disease associations for a human gene by symbol via MARRVEL: each association's phenotype name, gene and phenotype MIM numbers, mode of inheritance, and phenotypic series. Use to enumerate the Mendelian diseases linked to a gene and their inheritance patterns. MARRVEL indexes only current approved HGNC symbols, so a retired symbol or alias (e.g. 'RP20'/'LCA2' for RPE65) is resolved once via HGNC and re-queried; when that happens the substitution is disclosed in metadata.resolved_symbol and metadata.resolved_note. No API key required.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if symbol is not None:
        _dnakit_arguments['symbol'] = symbol
    if gene_symbol is not None:
        _dnakit_arguments['gene_symbol'] = gene_symbol
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    return _call('marrvel_get_omim_phenotypes', _dnakit_arguments, verbose=verbose)

def mgi_get_gene(*,
    gene_id: str,
    verbose: bool = True,
) -> _Object:
    "Get detailed information about a mouse gene from MGI (Mouse Genome Informatics) by its MGI identifier. Returns gene symbol, full name, species info, gene synopsis (curated and automated), synonyms, genomic location (GRCm39 assembly), SO term, and cross-references to UniProt, NCBI, Ensembl, etc. Powered by the Alliance of Genome Resources API which integrates MGI curated data. Example: 'MGI:98834' returns Trp53 (transformation related protein 53) details.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_id'] = gene_id
    return _call('mgi_get_gene', _dnakit_arguments, verbose=verbose)

def mgi_get_phenotypes(*,
    gene_id: str,
    limit: int | None = None,
    page: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Get phenotype annotations for a mouse gene from MGI (Mouse Genome Informatics). Returns Mammalian Phenotype (MP) ontology annotations including phenotype statements derived from mutant allele studies. Uses Alliance of Genome Resources API which integrates MGI phenotype curation data. Useful for understanding mouse gene function through loss-of-function and gain-of-function phenotypes. Example: 'MGI:98834' returns Trp53 phenotypes including tumor susceptibility, apoptosis defects, etc.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_id'] = gene_id
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if page is not None:
        _dnakit_arguments['page'] = page
    return _call('mgi_get_phenotypes', _dnakit_arguments, verbose=verbose)

def mgi_search_genes(*,
    query: str,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Search for mouse (Mus musculus) genes in MGI (Mouse Genome Informatics) by gene symbol, name, or keyword. Returns matching mouse gene symbols, names, and MGI identifiers. Powered by the Alliance of Genome Resources API which integrates MGI data. Use this to find mouse gene IDs before calling MGI_get_gene or MGI_get_phenotypes. Examples: 'Trp53' (tumor protein p53), 'Brca1', 'insulin', 'Notch'.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('mgi_search_genes', _dnakit_arguments, verbose=verbose)

def mgnify_get_genome(*,
    genome_id: str,
    verbose: bool = True,
) -> _Object:
    'Get detailed genome information from the MGnify genome catalog. Returns assembly quality metrics (completeness, contamination, N50), taxonomy, geographic origin, rRNA counts, protein count, pangenome statistics, and cross-references to ENA. Example: MGYG000000001 returns a high-quality Firmicutes isolate genome with 98.6% completeness.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['genome_id'] = genome_id
    return _call('mgnify_get_genome', _dnakit_arguments, verbose=verbose)

def mgnify_get_go_terms(*,
    analysis_id: str,
    page_size: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get Gene Ontology (GO) functional annotations from a MGnify metagenomics analysis. Returns GO terms with their accession, description, gene count, and category (biological_process, molecular_function, or cellular_component). Use to understand the functional capacity of a microbial community. Example: MGYA00585482 returns GO:0055114 oxidation-reduction process (count=1507).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['analysis_id'] = analysis_id
    if page_size is not None:
        _dnakit_arguments['page_size'] = page_size
    return _call('mgnify_get_go_terms', _dnakit_arguments, verbose=verbose)

def mgnify_get_interpro(*,
    analysis_id: str,
    page_size: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get InterPro protein domain annotations from a MGnify metagenomics analysis. Returns InterPro identifiers with description and gene counts, ranked by abundance. Useful for understanding the protein family composition of a microbial community. Example: MGYA00585482 returns IPR036388 Winged helix-like DNA-binding domain (count=360).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['analysis_id'] = analysis_id
    if page_size is not None:
        _dnakit_arguments['page_size'] = page_size
    return _call('mgnify_get_interpro', _dnakit_arguments, verbose=verbose)

def mgnify_get_study_detail(*,
    study_accession: str,
    verbose: bool = True,
) -> _Object:
    'Get detailed information about a specific MGnify metagenomics study. Returns study name, abstract, associated BioProject, centre, public/private status, sample count, and biomes. Complementary to MGnify_search_studies which lists studies; use MGnify_list_analyses / MGnify_list_analysis_downloads for exact analyses/downloads counts (this endpoint does not provide them). Example: MGYS00002008 is a Tara Oceans marine metagenomics study.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['study_accession'] = study_accession
    return _call('mgnify_get_study_detail', _dnakit_arguments, verbose=verbose)

def mgnify_get_taxonomy(*,
    analysis_id: str,
    rna_type: str | None = None,
    page_size: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get taxonomic composition from a MGnify metagenomics analysis. Returns organisms detected with their abundance counts, taxonomic lineage (domain/phylum/class/order/family/genus/species), and rank. Supports SSU (16S/18S) and LSU (23S/28S) rRNA marker genes. Use MGnify_list_analyses to find analysis IDs for a study. Example: MGYA00585482 returns Bacteria, Proteobacteria, Firmicutes lineages with counts.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['analysis_id'] = analysis_id
    if rna_type is not None:
        _dnakit_arguments['rna_type'] = rna_type
    if page_size is not None:
        _dnakit_arguments['page_size'] = page_size
    return _call('mgnify_get_taxonomy', _dnakit_arguments, verbose=verbose)

def mgnify_list_analyses(*,
    study_accession: str,
    size: int | None = None,
    verbose: bool = True,
) -> _Object:
    'List MGnify analyses via /analyses for a given study_accession\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['study_accession'] = study_accession
    if size is not None:
        _dnakit_arguments['size'] = size
    return _call('mgnify_list_analyses', _dnakit_arguments, verbose=verbose)

def mgnify_list_biomes(*,
    depth: int | None = None,
    page: int | None = None,
    page_size: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Browse the MGnify biome hierarchy. Returns biome identifiers, names, and sample counts. Biomes range from engineered environments (bioreactors) to natural habitats (ocean, soil, host-associated). Useful for discovering what metagenomic data is available for specific environments.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if depth is not None:
        _dnakit_arguments['depth'] = depth
    if page is not None:
        _dnakit_arguments['page'] = page
    if page_size is not None:
        _dnakit_arguments['page_size'] = page_size
    return _call('mgnify_list_biomes', _dnakit_arguments, verbose=verbose)

def mgnify_search_genomes(*,
    taxonomy: str | None = None,
    page: int | None = None,
    page_size: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Search the MGnify genome catalog for metagenome-assembled genomes (MAGs) and isolate genomes. Filter by taxonomy lineage. Returns quality metrics, taxonomy, and protein counts. Useful for finding reference genomes from specific microbial lineages. Taxonomy names follow GTDB, not NCBI: use 'Bacteroidota' (not 'Bacteroidetes'), 'Bacillota' (not 'Firmicutes'), 'Pseudomonadota' (not 'Proteobacteria') -- an NCBI-style name usually returns 0 hits rather than an error.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if taxonomy is not None:
        _dnakit_arguments['taxonomy'] = taxonomy
    if page is not None:
        _dnakit_arguments['page'] = page
    if page_size is not None:
        _dnakit_arguments['page_size'] = page_size
    return _call('mgnify_search_genomes', _dnakit_arguments, verbose=verbose)

def mgnify_search_studies(*,
    biome: str | None = None,
    search: str | None = None,
    size: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search MGnify studies via /studies with optional biome/search filters\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if biome is not None:
        _dnakit_arguments['biome'] = biome
    if search is not None:
        _dnakit_arguments['search'] = search
    if size is not None:
        _dnakit_arguments['size'] = size
    return _call('mgnify_search_studies', _dnakit_arguments, verbose=verbose)

def mave_db_get_effect_matrix(*,
    urn: str,
    uniprot_id: str | None | None = None,
    score_field: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    'One-shot DMS data loader: fetch a MaveDB score set, parse HGVS to (ref_aa, position, alt_aa), filter to single missense, optionally verify position numbering against a UniProt canonical sequence, and reshape into a (20 amino acids × n_positions) effect matrix ready for downstream analysis. Hides the boilerplate that variant-effect-prediction benchmarking and residue-mechanism interpretation skills used to inline. Returns the matrix, the position list, the standard amino acid order, and explicit numbering / drop / score-field metadata so the caller can audit what was loaded.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['urn'] = urn
    if uniprot_id is not None:
        _dnakit_arguments['uniprot_id'] = uniprot_id
    if score_field is not None:
        _dnakit_arguments['score_field'] = score_field
    return _call('mave_db_get_effect_matrix', _dnakit_arguments, verbose=verbose)

def med_rxiv_get_preprint(*,
    doi: str | None = None,
    server: str | None = None,
    query: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Get full metadata for a specific medRxiv preprint by its DOI. Returns comprehensive metadata including title, authors, abstract, publication date, version, category, and URLs. Response time: <1 second. Use when you have a medRxiv DOI (format: 10.1101/YYYY.MM.DD.XXXXXXX). To SEARCH for medRxiv preprints by keywords, use EuropePMC_search_articles with 'SRC:PPR' filter (fast, relevance-ranked), or web_search with 'site:medrxiv.org'. medRxiv API does not support text search.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if doi is not None:
        _dnakit_arguments['doi'] = doi
    if server is not None:
        _dnakit_arguments['server'] = server
    if query is not None:
        _dnakit_arguments['query'] = query
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('med_rxiv_get_preprint', _dnakit_arguments, verbose=verbose)

def monarch_v3_get_associations(*,
    subject: str | None = None,
    object: str | None = None,
    category: str,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Query cross-species associations between biomedical entities in the Monarch Initiative knowledge graph. Find gene-phenotype, disease-phenotype, gene-disease, or gene-pathway associations. Uses Biolink model categories to specify association type. Integrates data from OMIM, ClinVar, HPO, MGI, ZFIN, FlyBase, WormBase. Associations are DIRECTED: every category has a fixed subject side and object side, so anchor your entity on the correct side. Supply 'subject', 'object', or both (at least one is required). Example: HGNC:11998 (TP53) as subject with category GeneToPhenotypicFeatureAssociation returns phenotypes like hepatocellular carcinoma; MONDO:0006559 as OBJECT with category CausalGeneToDiseaseAssociation returns the genes that cause that disease (NCSTN, PSENEN, PSEN1).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if subject is not None:
        _dnakit_arguments['subject'] = subject
    if object is not None:
        _dnakit_arguments['object'] = object
    _dnakit_arguments['category'] = category
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('monarch_v3_get_associations', _dnakit_arguments, verbose=verbose)

def monarch_v3_get_entity(*,
    entity_id: str,
    verbose: bool = True,
) -> _Object:
    'Look up detailed information about any biomedical entity in the Monarch Initiative knowledge graph by its CURIE identifier. Supports genes (HGNC), diseases (MONDO, OMIM), phenotypes (HP), variants, and more. Returns the entity name, category, description, synonyms, cross-references, and species. Example: HGNC:11998 returns TP53 gene with cross-refs to ENSEMBL, OMIM; MONDO:0005148 returns type 2 diabetes mellitus with detailed description.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['entity_id'] = entity_id
    return _call('monarch_v3_get_entity', _dnakit_arguments, verbose=verbose)

def monarch_v3_get_histopheno(*,
    entity_id: str,
    verbose: bool = True,
) -> _Object:
    'Get the phenotype profile of a disease grouped by body system from the Monarch Initiative. Returns counts of phenotypic features (HPO terms) organized by anatomical system (e.g., nervous_system, skeletal_system, eye, cardiovascular_system). Only systems with at least one phenotype are returned, sorted by count descending. Useful for understanding which organ systems are most affected by a disease. Example: MONDO:0007947 (Marfan syndrome) returns high counts for head_neck (58), musculature (58), skeletal_system (56), eye (31).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['entity_id'] = entity_id
    return _call('monarch_v3_get_histopheno', _dnakit_arguments, verbose=verbose)

def monarch_v3_get_mappings(*,
    entity_id: str,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get cross-ontology mappings for a biomedical entity from the Monarch Initiative. Returns equivalence mappings between ontology identifiers, such as MONDO to OMIM, ICD10, SNOMED CT, NCIT, MESH, DOID, and Orphanet. Useful for translating disease identifiers across coding systems (e.g., for clinical data integration, ICD-10 billing codes, or linking to OMIM entries). Example: MONDO:0005148 (type 2 diabetes) maps to ICD10CM:E11, DOID:9352, OMIM:125853, MESH:D003924.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['entity_id'] = entity_id
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('monarch_v3_get_mappings', _dnakit_arguments, verbose=verbose)

def monarch_v3_phenotype_similarity_search(*,
    phenotypes: list[str],
    group: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Find diseases that best match a set of phenotypes (HPO terms) using semantic similarity from the Monarch Initiative. Given a patient phenotype profile as a list of HPO IDs, returns diseases ranked by phenotype overlap score. Uses information content-based similarity (AncestorInformationContent metric). Essential for rare disease differential diagnosis: input observed phenotypes, get ranked candidate diseases. Example: [HP:0001250 (seizure), HP:0001249 (intellectual disability), HP:0002353 (EEG abnormality)] returns diseases like intellectual disability AR 31, Dravet syndrome ranked by semantic similarity score.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['phenotypes'] = phenotypes
    if group is not None:
        _dnakit_arguments['group'] = group
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('monarch_v3_phenotype_similarity_search', _dnakit_arguments, verbose=verbose)

def monarch_v3_search(*,
    query: str | None = None,
    limit: int | None = None,
    category: str | None = None,
    q: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Search the Monarch Initiative knowledge graph for biomedical entities by name or keyword. Returns genes, diseases, phenotypes, and other entities matching the query, with their CURIE identifiers and categories. Useful for finding entity IDs to use with MonarchV3_get_entity and MonarchV3_get_associations. Example: searching 'BRCA1' returns HGNC:1100 (gene), 'Alzheimer' returns MONDO disease entries.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if query is not None:
        _dnakit_arguments['query'] = query
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if category is not None:
        _dnakit_arguments['category'] = category
    if q is not None:
        _dnakit_arguments['q'] = q
    return _call('monarch_v3_search', _dnakit_arguments, verbose=verbose)

def monarch_get_gene_diseases(*,
    subject: str,
    limit: int | None = None,
    offset: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Get diseases associated with a gene from Monarch Initiative. Input is a gene CURIE (e.g., 'HGNC:11998' for TP53). Returns disease associations with evidence. Bridges gene→disease across species for hypothesis generation.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['subject'] = subject
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    return _call('monarch_get_gene_diseases', _dnakit_arguments, verbose=verbose)

def monarch_get_gene_phenotypes(*,
    subject: str,
    limit: int | None = None,
    offset: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Get phenotypes associated with a gene from Monarch Initiative (cross-species). Input is a gene CURIE (e.g., 'HGNC:11998' for TP53). Returns human and model organism phenotypes (HP/MP terms). Use Monarch_search_gene first to find the CURIE.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['subject'] = subject
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    return _call('monarch_get_gene_phenotypes', _dnakit_arguments, verbose=verbose)

def monarch_search(*,
    query: str,
    category: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search the Monarch Initiative knowledge graph for biomedical entities by name or keyword. Returns genes, diseases, phenotypes, and other entities matching the query, with CURIE identifiers and categories. Filter by category (biolink:Gene, biolink:Disease, biolink:PhenotypicFeature) to narrow results. Use returned CURIEs with Monarch_get_entity and Monarch_get_associations for detailed lookups.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if category is not None:
        _dnakit_arguments['category'] = category
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('monarch_search', _dnakit_arguments, verbose=verbose)

def monarch_search_gene(*,
    query: str,
    limit: int | None = None,
    offset: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search Monarch Initiative knowledge graph for a gene by symbol or name. Returns HGNC/NCBIGene CURIEs, cross-species orthologs, and gene metadata. Use this first to get the gene CURIE needed for Monarch_get_gene_phenotypes and Monarch_get_gene_diseases.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    return _call('monarch_search_gene', _dnakit_arguments, verbose=verbose)

def my_gene_batch_query(*,
    gene_ids: list[str],
    species: str | None = None,
    fields: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Query multiple genes at once using a list of gene identifiers. Efficient for retrieving annotations for a gene set or performing ID mapping between different gene ID types.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_ids'] = gene_ids
    if species is not None:
        _dnakit_arguments['species'] = species
    if fields is not None:
        _dnakit_arguments['fields'] = fields
    return _call('my_gene_batch_query', _dnakit_arguments, verbose=verbose)

def my_gene_get_gene_annotation(*,
    gene_id: str,
    fields: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get detailed annotation for a specific gene by its ID (Entrez or Ensembl). Returns comprehensive information including gene ontology, pathways, protein domains, and literature references.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_id'] = gene_id
    if fields is not None:
        _dnakit_arguments['fields'] = fields
    return _call('my_gene_get_gene_annotation', _dnakit_arguments, verbose=verbose)

def my_gene_query_genes(*,
    query: str | None = None,
    species: str | None = None,
    fields: str | None = None,
    size: int | None = None,
    q: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Search for genes by keyword, symbol, name, or other identifiers. Returns gene annotations from 30+ sources including Entrez Gene, Ensembl, UniProt, HGNC. Useful for finding gene IDs, symbols, descriptions, and basic annotations.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if query is not None:
        _dnakit_arguments['query'] = query
    if species is not None:
        _dnakit_arguments['species'] = species
    if fields is not None:
        _dnakit_arguments['fields'] = fields
    if size is not None:
        _dnakit_arguments['size'] = size
    if q is not None:
        _dnakit_arguments['q'] = q
    return _call('my_gene_query_genes', _dnakit_arguments, verbose=verbose)

def my_variant_get_pathogenicity_scores(*,
    variant_id: str,
    fields: str | None = None,
    assembly: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Get comprehensive pathogenicity prediction scores for a variant from dbNSFP. Returns REVEL, CADD, AlphaMissense, SIFT, PolyPhen-2, MetaRNN, GERP, PhyloP, and PhastCons scores in a single call. Use this when MyVariant_query_variants does not return dbnsfp data, or when you specifically need REVEL/AlphaMissense scores for ACMG PP3/BP4 classification. Coordinates default to hg19/GRCh37 -- set assembly='hg38' for hg38/GRCh38 (e.g. rs4244285 is chr10:g.96541616G>C on hg19 and chr10:g.94781859G>C on hg38). An HGVS variant_id must be written in the SAME assembly you request: this endpoint looks ids up verbatim and never lifts coordinates over, so an hg19 id with assembly='hg38' is a 404, not a translation. An rsID works in either assembly. The response reports the frame in a top-level 'coordinate_assembly' key. MyVariant sometimes files several records under one rsID and resolves the rsID to the highest-scoring one, which can be an identifier-only stub; when that happens the sibling records that do carry the requested fields are reported in the top-level 'sibling_records_with_requested_fields', 'sibling_variant_ids_with_requested_fields' and 'sibling_record_note' keys, and 'data' still holds the record MyVariant resolved.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['variant_id'] = variant_id
    if fields is not None:
        _dnakit_arguments['fields'] = fields
    if assembly is not None:
        _dnakit_arguments['assembly'] = assembly
    return _call('my_variant_get_pathogenicity_scores', _dnakit_arguments, verbose=verbose)

def my_variant_get_variant_annotation(*,
    variant_id: str,
    fields: str | None = None,
    assembly: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Get detailed annotations for a specific variant by its HGVS ID (e.g., 'chr1:g.2186318G>A') or rsID. Returns information from dbSNP, ClinVar, CADD, and other sources. Coordinates default to hg19/GRCh37 -- set assembly='hg38' for hg38/GRCh38 (e.g. rs4244285 is chr10:g.96541616G>C on hg19 and chr10:g.94781859G>C on hg38). An HGVS variant_id must be written in the SAME assembly you request: this endpoint looks ids up verbatim and never lifts coordinates over, so an hg19 id with assembly='hg38' is a 404, not a translation. An rsID works in either assembly. The response reports the frame in a top-level 'coordinate_assembly' key. MyVariant sometimes files several records under one rsID and resolves the rsID to the highest-scoring one, which can be an identifier-only stub; when that happens the sibling records that do carry the requested fields are reported in the top-level 'sibling_records_with_requested_fields', 'sibling_variant_ids_with_requested_fields' and 'sibling_record_note' keys, and 'data' still holds the record MyVariant resolved.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['variant_id'] = variant_id
    if fields is not None:
        _dnakit_arguments['fields'] = fields
    if assembly is not None:
        _dnakit_arguments['assembly'] = assembly
    return _call('my_variant_get_variant_annotation', _dnakit_arguments, verbose=verbose)

def my_variant_query_variants(*,
    query: str | None = None,
    fields: str | None = None,
    size: int | None = None,
    assembly: str | None = None,
    variant_id: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Search for genetic variants by various criteria including rsID, gene, chromosome position, or clinical significance. Returns variant annotations from sources like dbSNP, ClinVar, gnomAD, and CADD. Coordinates default to hg19/GRCh37 -- set assembly='hg38' for hg38/GRCh38 (e.g. rs4244285 is chr10:g.96541616G>C on hg19 and chr10:g.94781859G>C on hg38). Any HGVS coordinate in the query must be written in the assembly you request; an rsID works in either. The response reports the frame in a top-level 'coordinate_assembly' key.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if query is not None:
        _dnakit_arguments['query'] = query
    if fields is not None:
        _dnakit_arguments['fields'] = fields
    if size is not None:
        _dnakit_arguments['size'] = size
    if assembly is not None:
        _dnakit_arguments['assembly'] = assembly
    if variant_id is not None:
        _dnakit_arguments['variant_id'] = variant_id
    return _call('my_variant_query_variants', _dnakit_arguments, verbose=verbose)

def ncbidatasets_get_gene(*,
    gene_id: str,
    verbose: bool = True,
) -> _Object:
    'Get comprehensive gene information from the National Center for Biotechnology Information (NCBI) by Gene ID. Returns gene symbol, description, chromosomal location, genomic coordinates, cross-references to UniProt/Ensembl/OMIM, synonyms, and annotation details. Covers all organisms in the NCBI Gene database. No API key is required; optionally set NCBI_API_KEY for a higher request-rate limit. Example: gene_id 7157 returns TP53 (tumor protein p53) on chromosome 17.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_id'] = gene_id
    return _call('ncbidatasets_get_gene', _dnakit_arguments, verbose=verbose)

def ncbidatasets_get_gene_by_symbol(*,
    symbol: str,
    taxon: str | None = None,
    organism: str | None = None,
    species: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Look up gene information by gene symbol and organism. Searches NCBI Gene database using official gene symbols. Specify the species with `taxon` (or its synonyms `organism`/`species`) as a common name (human, mouse, rat) or taxonomy ID. IMPORTANT -- the species defaults to HUMAN when none is given, so a gene symbol shared across species returns the human gene unless you say otherwise; `metadata.query_taxon_defaulted` is true whenever that default was applied, and `metadata.query_taxon` echoes the species actually used. Every returned row also carries its own `tax_id`/`taxname`, which is the authoritative check. Returns gene ID, description, chromosomal location, and database cross-references.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['symbol'] = symbol
    if taxon is not None:
        _dnakit_arguments['taxon'] = taxon
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if species is not None:
        _dnakit_arguments['species'] = species
    return _call('ncbidatasets_get_gene_by_symbol', _dnakit_arguments, verbose=verbose)

def ncbidatasets_get_genome_assembly(*,
    accession: str,
    verbose: bool = True,
) -> _Object:
    'Get genome assembly metadata for an assembly accession (RefSeq GCF_ or GenBank GCA_) from the NCBI Datasets v2 API. Returns organism/strain, assembly name/level/status, RefSeq category (e.g. reference genome), release date, submitter, BioProject, and assembly statistics (total length, chromosome count, contig/scaffold N50, GC%, annotation provider). Use to characterize a specific genome assembly across any organism (bacteria, fungi, plant, animal). No API key required.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['accession'] = accession
    return _call('ncbidatasets_get_genome_assembly', _dnakit_arguments, verbose=verbose)

def ncbidatasets_get_orthologs(*,
    gene_id: str,
    page_size: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Get orthologous genes across species for a given NCBI Gene ID. Returns orthologs identified by NCBI's Ortholog pipeline, including gene symbols, species, and gene types for each ortholog. Useful for comparative genomics and cross-species gene function studies. Example: gene_id 7157 (TP53) returns mouse Trp53, rat Tp53, zebrafish tp53, etc.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_id'] = gene_id
    if page_size is not None:
        _dnakit_arguments['page_size'] = page_size
    return _call('ncbidatasets_get_orthologs', _dnakit_arguments, verbose=verbose)

def ncbidatasets_get_sequence_reports(*,
    accession: str,
    page_size: int | None | None = None,
    page_token: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Get a page of per-sequence (chromosome/plasmid/scaffold) reports for a genome assembly via the NCBI Datasets v2 API. Returns each molecule's chromosome name, role, location type, RefSeq and GenBank accessions, length, and GC%, plus total_available and next_page_token metadata. Use next_page_token to retrieve additional pages for assemblies with many sequences. Use to map an assembly to its individual sequence accessions (e.g. resolve E. coli K-12 to NC_000913.3, or human chromosomes to their RefSeq NC_ accessions). No API key required.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['accession'] = accession
    if page_size is not None:
        _dnakit_arguments['page_size'] = page_size
    if page_token is not None:
        _dnakit_arguments['page_token'] = page_token
    return _call('ncbidatasets_get_sequence_reports', _dnakit_arguments, verbose=verbose)

def ncbidatasets_get_taxonomy(*,
    tax_id: str,
    verbose: bool = True,
) -> _Object:
    'Get detailed taxonomy information from NCBI by taxonomy ID. Returns organism name, rank, full lineage (as bare ancestor tax_ids in `lineage`, and as resolved organism_name/rank pairs in `lineage_names`), current counts of genes/assemblies/RNA types, and child taxa. Covers all organisms in NCBI Taxonomy including bacteria, archaea, eukaryotes, and viruses. Example: tax_id 9606 returns Homo sapiens and its root-to-parent lineage.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['tax_id'] = tax_id
    return _call('ncbidatasets_get_taxonomy', _dnakit_arguments, verbose=verbose)

def ncbidatasets_list_genomes_by_taxon(*,
    taxon: str,
    limit: int | None | None = None,
    reference_only: bool | None | None = None,
    verbose: bool = True,
) -> _Object:
    'List genome assemblies available for a taxon (NCBI tax id or scientific name) via the NCBI Datasets v2 API. Returns a summary list of assemblies (accession, assembly name/level, RefSeq category, release date, total length, N50, GC%) plus the total count available. Use to discover what genomes exist for an organism or clade, e.g. all Escherichia coli or Mycobacterium tuberculosis assemblies. Set reference_only to restrict to reference/representative genomes. No API key required.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['taxon'] = taxon
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if reference_only is not None:
        _dnakit_arguments['reference_only'] = reference_only
    return _call('ncbidatasets_list_genomes_by_taxon', _dnakit_arguments, verbose=verbose)

def ncbidatasets_suggest_taxonomy(*,
    query: str,
    verbose: bool = True,
) -> _Object:
    "Suggest taxonomic names matching a query string. Searches NCBI Taxonomy for organisms by partial name, returning matching scientific names, taxonomy IDs, common names, and ranks. Useful for finding taxonomy IDs from partial organism names. Example: 'drosophila' returns D. melanogaster (7227), D. simulans (7240), etc.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    return _call('ncbidatasets_suggest_taxonomy', _dnakit_arguments, verbose=verbose)

def ncbigene_get_summary(*,
    id: str,
    verbose: bool = True,
) -> _Object:
    'Get a detailed gene summary from NCBI Gene by Entrez Gene ID using the E-utilities ESummary API. Returns gene name, symbol, description, function summary, chromosomal location, map position, aliases, genomic coordinates, exon count, and MIM (OMIM) IDs. The most authoritative source for human gene information. Use NCBIGene_search first to find Gene IDs. Example: Gene ID 7157 returns TP53 on chromosome 17p13.1 with a detailed functional summary.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['id'] = id
    return _call('ncbigene_get_summary', _dnakit_arguments, verbose=verbose)

def ncbigene_search(*,
    term: str,
    retmax: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Search the NCBI Gene database using the E-utilities API. Find gene IDs (Entrez Gene IDs) by gene symbol, name, or keyword combined with an organism filter. Returns a list of matching Gene IDs that can be used with NCBIGene_get_summary. Essential for resolving gene names to standardized NCBI Gene IDs. Example: searching 'TP53' in Homo sapiens returns Gene ID 7157; 'BRCA1' returns 672.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['term'] = term
    if retmax is not None:
        _dnakit_arguments['retmax'] = retmax
    return _call('ncbigene_search', _dnakit_arguments, verbose=verbose)

def ncbiprotein_get_summary(*,
    id: str,
    verbose: bool = True,
) -> _Object:
    "Get protein record summary from the NCBI Protein database by GI number or accession using the E-utilities ESummary API. Returns protein title, organism, accession version, sequence length, molecular type, source database, and creation/update dates. NCBI Protein aggregates sequences from GenBank, RefSeq, Swiss-Prot, PDB, and other sources. Use for resolving NCBI protein GI numbers to names, accessions, and species. Example: GI 119395750 returns 'keratin, type II cytoskeletal 1 [Homo sapiens]' with accession NP_006112.3 and length 644 aa.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['id'] = id
    return _call('ncbiprotein_get_summary', _dnakit_arguments, verbose=verbose)

def ncbi_sra_search_runs(*,
    operation: _Object | None = None,
    study: str | None = None,
    organism: str | None = None,
    strategy: str | None = None,
    platform: str | None = None,
    source: str | None = None,
    query: str | None = None,
    limit: int | None = None,
    sort: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Search NCBI Sequence Read Archive (SRA) for NGS/RNA-seq sequencing runs by study accession, organism, sequencing strategy (RNA-Seq, WGS, ChIP-Seq), platform (ILLUMINA, OXFORD_NANOPORE), or library source. Returns SRA UIDs which can be used with other SRA tools to retrieve metadata and download links. Essential for discovering high-throughput sequencing data from public repositories. Use this to find relevant sequencing experiments before accessing the actual FASTQ data files.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if study is not None:
        _dnakit_arguments['study'] = study
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if strategy is not None:
        _dnakit_arguments['strategy'] = strategy
    if platform is not None:
        _dnakit_arguments['platform'] = platform
    if source is not None:
        _dnakit_arguments['source'] = source
    if query is not None:
        _dnakit_arguments['query'] = query
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if sort is not None:
        _dnakit_arguments['sort'] = sort
    return _call('ncbi_sra_search_runs', _dnakit_arguments, verbose=verbose)

def ncbi_fetch_accessions(*,
    operation: _Object | None = None,
    uids: _Object,
    verbose: bool = True,
) -> _Object:
    'Convert GenBank UIDs to accession numbers (U00096, NC_045512, etc.). Takes UIDs from NCBI_search_nucleotide and returns actual accession numbers that can be used with database tools. Returns both EMBL/GenBank accessions (U*, M*, AJ*) that work with ENA Browser tools, and RefSeq accessions (NC_*, NM_*) for NCBI tools. Use this as second step after NCBI_search_nucleotide to get usable accession numbers.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['uids'] = uids
    return _call('ncbi_fetch_accessions', _dnakit_arguments, verbose=verbose)

def ncbi_get_sequence(*,
    operation: _Object | None = None,
    accession: str,
    format: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Retrieve DNA/RNA sequence data from NCBI by accession number. Returns sequences in specified format (fasta, gb, embl). Works with both EMBL/GenBank accessions (U00096, M15390) and RefSeq accessions (NC_045512, NM_000546). Use this to get actual sequence data after finding accessions with NCBI_search_nucleotide. Returns full sequence with annotations depending on format.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['accession'] = accession
    if format is not None:
        _dnakit_arguments['format'] = format
    return _call('ncbi_get_sequence', _dnakit_arguments, verbose=verbose)

def ncbi_search_nucleotide(*,
    operation: _Object | None = None,
    organism: str | None = None,
    gene: str | None = None,
    strain: str | None = None,
    keywords: str | None = None,
    seq_type: str | None = None,
    query: str | None = None,
    limit: int | None = None,
    sort: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Search NCBI Nucleotide database (GenBank/EMBL/RefSeq) by organism, gene name, or keywords. Returns GenBank UIDs which can be used to fetch accession numbers. Supports filtering by organism, gene, strain, sequence type (complete genome, mRNA, RefSeq). Use this to discover DNA/RNA sequences for specific organisms or genes. Returns both EMBL/GenBank accessions (U00096, M15390) and RefSeq accessions (NC_*, NM_*). Follow up with NCBI_fetch_accessions to get actual accession numbers.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    if strain is not None:
        _dnakit_arguments['strain'] = strain
    if keywords is not None:
        _dnakit_arguments['keywords'] = keywords
    if seq_type is not None:
        _dnakit_arguments['seq_type'] = seq_type
    if query is not None:
        _dnakit_arguments['query'] = query
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if sort is not None:
        _dnakit_arguments['sort'] = sort
    return _call('ncbi_search_nucleotide', _dnakit_arguments, verbose=verbose)

def nvidia_nim_alphafold2(*,
    sequence: str,
    algorithm: str | None = None,
    databases: list[str] | None = None,
    e_value: float | None = None,
    iterations: int | None = None,
    relax_prediction: bool | None = None,
    skip_template_search: bool | None = None,
    verbose: bool = True,
) -> _Object:
    "Predict protein 3D structure from amino acid sequence using DeepMind's AlphaFold2 via NVIDIA NIM. Returns PDB format structure. Requires NVIDIA_API_KEY. Input: sequence (amino acids), algorithm (jackhmmer/mmseqs2), databases (uniref90, mgnify, small_bfd). Long-running async operation.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['sequence'] = sequence
    if algorithm is not None:
        _dnakit_arguments['algorithm'] = algorithm
    if databases is not None:
        _dnakit_arguments['databases'] = databases
    if e_value is not None:
        _dnakit_arguments['e_value'] = e_value
    if iterations is not None:
        _dnakit_arguments['iterations'] = iterations
    if relax_prediction is not None:
        _dnakit_arguments['relax_prediction'] = relax_prediction
    if skip_template_search is not None:
        _dnakit_arguments['skip_template_search'] = skip_template_search
    return _call('nvidia_nim_alphafold2', _dnakit_arguments, verbose=verbose)

def nvidia_nim_esmfold(*,
    sequence: str,
    verbose: bool = True,
) -> _Object:
    'Fast alignment-free protein structure prediction using ESMFold via NVIDIA NIM. No MSA required - predicts directly from sequence. Input: sequence (max 1024 amino acids). Returns PDB structure. Fast synchronous response.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['sequence'] = sequence
    return _call('nvidia_nim_esmfold', _dnakit_arguments, verbose=verbose)

def omim_get_clinical_synopsis(*,
    operation: str | None = None,
    mim_number: str,
    verbose: bool = True,
) -> _Object:
    'Get clinical synopsis (phenotype features) for an OMIM phenotype entry. Returns structured clinical features organized by body system (e.g., neurologic, cardiovascular). Essential for phenotype-based disease matching. Requires OMIM_API_KEY.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['mim_number'] = mim_number
    return _call('omim_get_clinical_synopsis', _dnakit_arguments, verbose=verbose)

def omim_get_entry(*,
    operation: str | None = None,
    mim_number: str,
    include: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get detailed OMIM entry by MIM number. Returns comprehensive information including text description, clinical features, inheritance, molecular genetics, and references. MIM numbers are 6-digit identifiers (e.g., 164730 for BRAF gene). Requires OMIM_API_KEY.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['mim_number'] = mim_number
    if include is not None:
        _dnakit_arguments['include'] = include
    return _call('omim_get_entry', _dnakit_arguments, verbose=verbose)

def omim_get_gene_map(*,
    operation: str | None = None,
    mim_number: str | None = None,
    chromosome: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get OMIM gene-disease mapping information. Returns genes mapped to chromosomal locations with associated phenotypes and inheritance patterns. Query by MIM number or chromosome. Useful for finding all disease genes on a chromosome.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if mim_number is not None:
        _dnakit_arguments['mim_number'] = mim_number
    if chromosome is not None:
        _dnakit_arguments['chromosome'] = chromosome
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('omim_get_gene_map', _dnakit_arguments, verbose=verbose)

def omim_search(*,
    operation: str | None = None,
    query: str,
    limit: int | None = None,
    start: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search OMIM database for Mendelian diseases and genes. OMIM is the authoritative source for human genetic disorders. Search by gene name, disease name, or phenotype keywords. Returns MIM numbers, titles, and entry types. Requires OMIM_API_KEY env var (register at omim.org/api).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['query'] = query
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if start is not None:
        _dnakit_arguments['start'] = start
    return _call('omim_search', _dnakit_arguments, verbose=verbose)

def osl_get_efo_id_by_disease_name(*,
    disease: str,
    verbose: bool = True,
) -> _Object:
    'Legacy helper to find an EFO term ID from a disease name using EMBL-EBI OLS search.\n\nReturns the best match as `{efo_id, name}` (or `null` if none). Note: the tool name contains a legacy typo (`OSL_` vs `OLS_`). If you need richer ontology data (IRI, synonyms, children), use the newer `ols_*` tools in `efo_tools.json`.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['disease'] = disease
    return _call('osl_get_efo_id_by_disease_name', _dnakit_arguments, verbose=verbose)

def onco_kb_annotate_variant(*,
    operation: str | None = None,
    gene: str | None = None,
    gene_symbol: str | None = None,
    variant: str | None = None,
    tumor_type: str | None = None,
    alteration: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Annotate a cancer variant for oncogenic potential and treatment implications using OncoKB. Returns evidence level, oncogenicity, mutation effect, and FDA-approved therapies. Requires ONCOKB_API_TOKEN env var (demo mode available for BRAF, TP53, ROS1). Example: BRAF V600E in melanoma.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    if gene_symbol is not None:
        _dnakit_arguments['gene_symbol'] = gene_symbol
    if variant is not None:
        _dnakit_arguments['variant'] = variant
    if tumor_type is not None:
        _dnakit_arguments['tumor_type'] = tumor_type
    if alteration is not None:
        _dnakit_arguments['alteration'] = alteration
    return _call('onco_kb_annotate_variant', _dnakit_arguments, verbose=verbose)

def open_targets_get_associated_drugs_by_disease_efo_id(*,
    efoId: str,
    size: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Retrieve known drugs associated with a specific disease by disease efoId.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['efoId'] = efoId
    if size is not None:
        _dnakit_arguments['size'] = size
    return _call('open_targets_get_associated_drugs_by_disease_efo_id', _dnakit_arguments, verbose=verbose)

def open_targets_get_associated_drugs_by_target_ensembl_id(*,
    ensemblId: str,
    size: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get known drugs and information (e.g. id, name, MoA) associated with a specific target ensemblID, including clinical trial phase and mechanism of action of the drugs.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['ensemblId'] = ensemblId
    if size is not None:
        _dnakit_arguments['size'] = size
    return _call('open_targets_get_associated_drugs_by_target_ensembl_id', _dnakit_arguments, verbose=verbose)

def open_targets_get_associated_phenotypes_by_disease_efo_id(*,
    efoId: str,
    verbose: bool = True,
) -> _Object:
    'Find HPO phenotypes asosciated with the specified disease efoId.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['efoId'] = efoId
    return _call('open_targets_get_associated_phenotypes_by_disease_efo_id', _dnakit_arguments, verbose=verbose)

def open_targets_get_associated_targets_by_disease_efo_id(*,
    efoId: str,
    size: int | None = None,
    index: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Find targets associated with a specific disease or phenotype based on efoId. Results are PAGINATED (50 per call by default): `associatedTargets.count` is the total number of associated targets, `associatedTargets.returned` is how many are in `rows`, and `associatedTargets.truncated` is true when targets remain unfetched, in which case `metadata.truncation_note` names the `size`/`index` values to pass. Raise `size` (Open Targets rejects a page above 3000) and step `index` to page through; a common disease has thousands of associated targets, so the first page is NOT the full list.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['efoId'] = efoId
    if size is not None:
        _dnakit_arguments['size'] = size
    if index is not None:
        _dnakit_arguments['index'] = index
    return _call('open_targets_get_associated_targets_by_disease_efo_id', _dnakit_arguments, verbose=verbose)

def open_targets_get_credible_set_detail(*,
    studyLocusId: str,
    verbose: bool = True,
) -> _Object:
    'Get detailed information about a specific GWAS credible set (fine-mapped locus) from Open Targets. Returns the lead variant, associated study, fine-mapping method, statistical significance, locus-to-gene (L2G) predictions mapping the signal to likely causal genes, and colocalization evidence with QTLs.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['studyLocusId'] = studyLocusId
    return _call('open_targets_get_credible_set_detail', _dnakit_arguments, verbose=verbose)

def open_targets_get_disease_ancestors_parents_by_efo_id(*,
    efoId: str,
    verbose: bool = True,
) -> _Object:
    'Retrieve the disease ancestors and parents in the ontology using the disease EFO ID.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['efoId'] = efoId
    return _call('open_targets_get_disease_ancestors_parents_by_efo_id', _dnakit_arguments, verbose=verbose)

def open_targets_get_disease_descendants_children_by_efo_id(*,
    efoId: str,
    verbose: bool = True,
) -> _Object:
    'Retrieve the disease descendants and children in the ontology using the disease EFO ID.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['efoId'] = efoId
    return _call('open_targets_get_disease_descendants_children_by_efo_id', _dnakit_arguments, verbose=verbose)

def open_targets_get_disease_description_by_efo_id(*,
    efoId: str,
    verbose: bool = True,
) -> _Object:
    "Retrieve disease description, name, database cros references, obsolete terms, and whether it's a therapeutic area, all using the specified efoId.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['efoId'] = efoId
    return _call('open_targets_get_disease_description_by_efo_id', _dnakit_arguments, verbose=verbose)

def open_targets_get_disease_id_description_by_name(*,
    diseaseName: str | None = None,
    query: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Retrieve the efoId and additional details of a disease based on its name.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if diseaseName is not None:
        _dnakit_arguments['diseaseName'] = diseaseName
    if query is not None:
        _dnakit_arguments['query'] = query
    return _call('open_targets_get_disease_id_description_by_name', _dnakit_arguments, verbose=verbose)

def open_targets_get_disease_synonyms_by_efo_id(*,
    efoId: str,
    verbose: bool = True,
) -> _Object:
    'Retrieve disease synonyms by its EFO ID.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['efoId'] = efoId
    return _call('open_targets_get_disease_synonyms_by_efo_id', _dnakit_arguments, verbose=verbose)

def open_targets_get_disease_therapeutic_areas_by_efo_id(*,
    efoId: str,
    verbose: bool = True,
) -> _Object:
    'Retrieve the therapeutic areas associated with a specific disease efoId.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['efoId'] = efoId
    return _call('open_targets_get_disease_therapeutic_areas_by_efo_id', _dnakit_arguments, verbose=verbose)

def open_targets_get_diseases_phenotypes_by_target_ensembl(*,
    ensemblId: str,
    page: dict[str, _Object] | None = None,
    verbose: bool = True,
) -> _Object:
    "Find diseases or phenotypes associated with a specific target using ensemblId. `associatedDiseases.count` is the total number of associations; `rows` returns only one page of them (25 by default) -- pass `page` to retrieve more (e.g. page 2 of a target with thousands of associations like SIRT1's 3432). `associatedDiseases.returned` reports how many rows this response actually contains and `associatedDiseases.truncated` is true when associations remain unfetched, in which case `metadata.truncation_note` names the exact `page` argument to pass.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['ensemblId'] = ensemblId
    if page is not None:
        _dnakit_arguments['page'] = page
    return _call('open_targets_get_diseases_phenotypes_by_target_ensembl', _dnakit_arguments, verbose=verbose)

def open_targets_get_drug_adverse_events_by_chembl_id(*,
    chemblId: str,
    page: dict[str, _Object] | None = None,
    verbose: bool = True,
) -> _Object:
    'Retrieve significant adverse events reported for a specific drug chemblId, ranked by log-likelihood ratio. Results are PAGINATED and Open Targets returns only 25 rows when no `page` is given: `adverseEvents.count` is the total number of signals, `adverseEvents.returned` is how many are present in `rows`, and `adverseEvents.truncated` is true when signals remain unfetched, in which case a `truncation_note` naming the exact `page` argument to pass is added to `metadata`. Pass page: {"index": 0, "size": <count>} to retrieve every signal before drawing any safety conclusion — a default page is NOT the full list (ibuprofen, CHEMBL521, has 55 signals; the default 25 omit toxic epidermal necrolysis, urticaria, systemic lupus erythematosus and hepatic enzyme increased).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['chemblId'] = chemblId
    if page is not None:
        _dnakit_arguments['page'] = page
    return _call('open_targets_get_drug_adverse_events_by_chembl_id', _dnakit_arguments, verbose=verbose)

def open_targets_get_drug_chemb_id_by_generic_name(*,
    drugName: str | None = None,
    query: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Fetch the drug chemblId and description based on the drug generic name.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if drugName is not None:
        _dnakit_arguments['drugName'] = drugName
    if query is not None:
        _dnakit_arguments['query'] = query
    return _call('open_targets_get_drug_chemb_id_by_generic_name', _dnakit_arguments, verbose=verbose)

def open_targets_get_drug_indications_by_chembl_id(*,
    chemblId: str,
    verbose: bool = True,
) -> _Object:
    'Fetch indications (treatable phenotypes/diseases) for a given drug chemblId.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['chemblId'] = chemblId
    return _call('open_targets_get_drug_indications_by_chembl_id', _dnakit_arguments, verbose=verbose)

def open_targets_get_drug_mechanisms_of_action_by_chembl_id(*,
    chemblId: str,
    verbose: bool = True,
) -> _Object:
    'Retrieve the mechanisms of action associated with a specific drug using chemblId.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['chemblId'] = chemblId
    return _call('open_targets_get_drug_mechanisms_of_action_by_chembl_id', _dnakit_arguments, verbose=verbose)

def open_targets_get_drug_warnings_by_chembl_id(*,
    chemblId: str,
    verbose: bool = True,
) -> _Object:
    'Retrieve warnings for a specific drug using ChEMBL ID.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['chemblId'] = chemblId
    return _call('open_targets_get_drug_warnings_by_chembl_id', _dnakit_arguments, verbose=verbose)

def open_targets_get_evidence_by_datasource(*,
    efoId: str | None = None,
    ensemblId: str | None = None,
    gene_symbol: str | None = None,
    disease_name: str | None = None,
    datasourceIds: list[str] | None = None,
    size: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Get target-disease evidence from Open Targets filtered by configurable data sources. Unlike OpenTargets_target_disease_evidence (intogen-only), this tool accepts any datasourceIds. Valid sources include: gwas_credible_sets, eva, eva_somatic, gene_burden, genomics_england, intogen, cancer_gene_census, cancer_biomarkers, clinical_precedence, crispr, europepmc, expression_atlas, gene2phenotype, impc, orphanet, uniprot_literature, uniprot_variants, clingen, phewas_catalog, reactome, progeny, slapenrich, sysbio. Note: two source IDs were renamed upstream and the retired names now match nothing rather than erroring — 'ot_genetics_portal' became 'gwas_credible_sets' (confirmed live: IL23R/Crohn disease returns 0 rows for the old name and 43 for the new one) and 'chembl' became 'clinical_precedence'. This tool remaps both retired names automatically and reports the substitution in metadata.datasource_rename_note. Returns evidence rows with scores, literature references, and source attribution.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if efoId is not None:
        _dnakit_arguments['efoId'] = efoId
    if ensemblId is not None:
        _dnakit_arguments['ensemblId'] = ensemblId
    if gene_symbol is not None:
        _dnakit_arguments['gene_symbol'] = gene_symbol
    if disease_name is not None:
        _dnakit_arguments['disease_name'] = disease_name
    if datasourceIds is not None:
        _dnakit_arguments['datasourceIds'] = datasourceIds
    if size is not None:
        _dnakit_arguments['size'] = size
    return _call('open_targets_get_evidence_by_datasource', _dnakit_arguments, verbose=verbose)

def open_targets_get_gwas_study(*,
    studyId: str,
    verbose: bool = True,
) -> _Object:
    'Get detailed metadata for a GWAS study from Open Targets Platform, including trait information, sample sizes, publication details, disease mappings, cohorts, LD population structure, and availability of summary statistics. Use GCST accession IDs.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['studyId'] = studyId
    return _call('open_targets_get_gwas_study', _dnakit_arguments, verbose=verbose)

def open_targets_get_publications_by_target_ensembl_id(*,
    entityId: str,
    additionalIds: list[str] | None = None,
    startYear: int | None = None,
    startMonth: int | None = None,
    endYear: int | None = None,
    endMonth: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Retrieve publications related to a target ensemblID, including PubMed IDs and publication dates. Only the first 25 publications are returned and this tool exposes no pagination parameter: `literatureOcurrences.count` is the total (often hundreds of thousands), `literatureOcurrences.returned` is how many are in `rows`, and `literatureOcurrences.truncated` is true when publications remain unfetched, in which case `metadata.truncation_note` says so. Narrow the query with startYear/endYear/additionalIds rather than treating `rows` as the literature for this target.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['entityId'] = entityId
    if additionalIds is not None:
        _dnakit_arguments['additionalIds'] = additionalIds
    if startYear is not None:
        _dnakit_arguments['startYear'] = startYear
    if startMonth is not None:
        _dnakit_arguments['startMonth'] = startMonth
    if endYear is not None:
        _dnakit_arguments['endYear'] = endYear
    if endMonth is not None:
        _dnakit_arguments['endMonth'] = endMonth
    return _call('open_targets_get_publications_by_target_ensembl_id', _dnakit_arguments, verbose=verbose)

def open_targets_get_similar_entities_by_disease_efo_id(*,
    efoId: str,
    threshold: float,
    size: int,
    verbose: bool = True,
) -> _Object:
    'Retrieve similar entities for a given disease efoId using a model trained with PubMed.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['efoId'] = efoId
    _dnakit_arguments['threshold'] = threshold
    _dnakit_arguments['size'] = size
    return _call('open_targets_get_similar_entities_by_disease_efo_id', _dnakit_arguments, verbose=verbose)

def open_targets_get_study_credible_sets(*,
    studyIds: list[str],
    size: int | None = None,
    index: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get all credible sets (fine-mapped GWAS loci) for a specific GWAS study from Open Targets. Returns the lead variants, fine-mapping details, p-values, effect sizes, and locus-to-gene (L2G) predictions for each locus in the study. Use this to explore all significant loci discovered in a GWAS.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['studyIds'] = studyIds
    if size is not None:
        _dnakit_arguments['size'] = size
    if index is not None:
        _dnakit_arguments['index'] = index
    return _call('open_targets_get_study_credible_sets', _dnakit_arguments, verbose=verbose)

def open_targets_get_target_classes_by_ensembl_id(*,
    ensemblId: str,
    verbose: bool = True,
) -> _Object:
    'Retrieve the target classes associated with a specific target ensemblID.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['ensemblId'] = ensemblId
    return _call('open_targets_get_target_classes_by_ensembl_id', _dnakit_arguments, verbose=verbose)

def open_targets_get_target_gene_ontology_by_ensembl_id(*,
    ensemblId: str,
    verbose: bool = True,
) -> _Object:
    'Retrieve Gene Ontology annotations for a specific target by Ensembl ID.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['ensemblId'] = ensemblId
    return _call('open_targets_get_target_gene_ontology_by_ensembl_id', _dnakit_arguments, verbose=verbose)

def open_targets_get_target_genomic_location_by_ensembl_id(*,
    ensemblId: str,
    verbose: bool = True,
) -> _Object:
    'Retrieve genomic location data for a specific target, including chromosome, start, end, and strand.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['ensemblId'] = ensemblId
    return _call('open_targets_get_target_genomic_location_by_ensembl_id', _dnakit_arguments, verbose=verbose)

def open_targets_get_target_homologues_by_ensembl_id(*,
    ensemblId: str,
    verbose: bool = True,
) -> _Object:
    'Fetch homologues for a specific target by Ensembl ID.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['ensemblId'] = ensemblId
    return _call('open_targets_get_target_homologues_by_ensembl_id', _dnakit_arguments, verbose=verbose)

def open_targets_get_target_id_description_by_name(*,
    targetName: str,
    verbose: bool = True,
) -> _Object:
    'Get the ensemblId and description based on the target name.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['targetName'] = targetName
    return _call('open_targets_get_target_id_description_by_name', _dnakit_arguments, verbose=verbose)

def open_targets_get_target_interactions_by_ensembl_id(*,
    ensemblId: str,
    page: dict[str, _Object] | None = None,
    verbose: bool = True,
) -> _Object:
    'Retrieve interaction data for a specific target ensemblID, including interaction partners and evidence. Results are PAGINATED and Open Targets returns only 25 rows when no `page` is given: `interactions.count` is the total number of interactions (a well-studied target has thousands), `interactions.returned` is how many are in `rows`, and `interactions.truncated` is true when interactions remain unfetched, in which case `metadata.truncation_note` names the exact `page` argument to pass. Pass page: {"index": 0, "size": N} (Open Targets rejects a size above 3000) and step `index` to page through.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['ensemblId'] = ensemblId
    if page is not None:
        _dnakit_arguments['page'] = page
    return _call('open_targets_get_target_interactions_by_ensembl_id', _dnakit_arguments, verbose=verbose)

def open_targets_get_target_safety_profile_by_ensembl_id(*,
    ensemblId: str,
    verbose: bool = True,
) -> _Object:
    'Retrieve known target safety liabilities for a specific target Ensembl ID.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['ensemblId'] = ensemblId
    return _call('open_targets_get_target_safety_profile_by_ensembl_id', _dnakit_arguments, verbose=verbose)

def open_targets_get_target_tractability_by_ensembl_id(*,
    ensemblId: str,
    verbose: bool = True,
) -> _Object:
    'Retrieve tractability assessments, including modality and values, for a specific target ensembl ID.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['ensemblId'] = ensemblId
    return _call('open_targets_get_target_tractability_by_ensembl_id', _dnakit_arguments, verbose=verbose)

def open_targets_get_variant_credible_sets(*,
    variantId: str | None = None,
    size: int | None = None,
    variant_id: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get GWAS credible sets (fine-mapped loci) associated with a genetic variant from Open Targets. Returns studies where this variant appears in a fine-mapped credible set, with locus-to-gene (L2G) predictions, p-values, effect sizes, and the fine-mapping method used. Essential for understanding which GWAS signals a variant is part of and which genes it likely affects.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if variantId is not None:
        _dnakit_arguments['variantId'] = variantId
    if size is not None:
        _dnakit_arguments['size'] = size
    if variant_id is not None:
        _dnakit_arguments['variant_id'] = variant_id
    return _call('open_targets_get_variant_credible_sets', _dnakit_arguments, verbose=verbose)

def open_targets_get_variant_info(*,
    variantId: str,
    verbose: bool = True,
) -> _Object:
    'Get detailed information about a genetic variant from Open Targets Platform, including genomic coordinates, rsIDs, allele frequencies across populations (gnomAD), functional consequence, transcript consequences, and HGVS notation. Useful for GWAS variant annotation and population genetics.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['variantId'] = variantId
    return _call('open_targets_get_variant_info', _dnakit_arguments, verbose=verbose)

def open_targets_map_any_disease_id_to_all_other_ids(*,
    inputId: str,
    verbose: bool = True,
) -> _Object:
    'Given any known disease or phenotype ID (EFO, OMIM, MONDO, UMLS, ICD10, MedDRA, etc.), return all known cross-referenced IDs including the EFO ID.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['inputId'] = inputId
    return _call('open_targets_map_any_disease_id_to_all_other_ids', _dnakit_arguments, verbose=verbose)

def open_targets_multi_entity_search_by_query_string(*,
    queryString: str,
    entityNames: list[str] | None = None,
    page: dict[str, _Object] | None = None,
    verbose: bool = True,
) -> _Object:
    'Perform a multi-entity search based on a query string, filtering by entity names and pagination settings.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['queryString'] = queryString
    if entityNames is not None:
        _dnakit_arguments['entityNames'] = entityNames
    if page is not None:
        _dnakit_arguments['page'] = page
    return _call('open_targets_multi_entity_search_by_query_string', _dnakit_arguments, verbose=verbose)

def open_targets_search_gwas_studies_by_disease(*,
    diseaseIds: list[str] | None = None,
    enableIndirect: bool | None = None,
    disease_name: str | None = None,
    disease: str | None = None,
    trait: str | None = None,
    size: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search for GWAS studies associated with a disease or trait in Open Targets. Takes disease ontology IDs (EFO, MONDO) and returns matching GWAS studies with their metadata including sample sizes, traits, publications, and associated diseases. Use enableIndirect=true to include child disease terms.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if diseaseIds is not None:
        _dnakit_arguments['diseaseIds'] = diseaseIds
    if enableIndirect is not None:
        _dnakit_arguments['enableIndirect'] = enableIndirect
    if disease_name is not None:
        _dnakit_arguments['disease_name'] = disease_name
    if disease is not None:
        _dnakit_arguments['disease'] = disease
    if trait is not None:
        _dnakit_arguments['trait'] = trait
    if size is not None:
        _dnakit_arguments['size'] = size
    return _call('open_targets_search_gwas_studies_by_disease', _dnakit_arguments, verbose=verbose)

def open_targets_target_disease_evidence(*,
    efoId: str | None = None,
    ensemblId: str | None = None,
    gene_symbol: str | None = None,
    disease_name: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Explore IntOGen somatic driver evidence for a target-disease association. IMPORTANT: IntOGen only covers somatic tumor driver mutations — this tool returns empty data ({}) for non-cancer diseases (e.g., type 2 diabetes, hypertension). Empty results are expected and normal when querying non-oncology phenotypes. For non-cancer diseases or to query other data sources, use OpenTargets_get_evidence_by_datasource instead. Accepts gene_symbol (e.g. 'TP53') and disease_name (e.g. 'Crohn disease') as alternatives to ensemblId/efoId — IDs are auto-resolved via OpenTargets search. Results are PAGINATED and Open Targets returns only 25 evidence rows: `evidences.count` is the total, `evidences.returned` is how many are in `rows`, and `evidences.truncated` is true when evidence remains unfetched, in which case `metadata.truncation_note` says so. This tool exposes no pagination parameter, so use OpenTargets_get_evidence_by_datasource (which takes `size`) when the full evidence set is needed.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if efoId is not None:
        _dnakit_arguments['efoId'] = efoId
    if ensemblId is not None:
        _dnakit_arguments['ensemblId'] = ensemblId
    if gene_symbol is not None:
        _dnakit_arguments['gene_symbol'] = gene_symbol
    if disease_name is not None:
        _dnakit_arguments['disease_name'] = disease_name
    return _call('open_targets_target_disease_evidence', _dnakit_arguments, verbose=verbose)

def orphanet_get_classification(*,
    operation: str | None = None,
    orpha_code: str | int | None = None,
    lang: str | None = None,
    orphacode: int | str | None = None,
    orpha_id: int | str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get disease classification hierarchy from Orphanet. Returns parent and child disease categories showing where a disease fits in the rare disease taxonomy. Useful for understanding disease relationships.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if orpha_code is not None:
        _dnakit_arguments['orpha_code'] = orpha_code
    if lang is not None:
        _dnakit_arguments['lang'] = lang
    if orphacode is not None:
        _dnakit_arguments['orphacode'] = orphacode
    if orpha_id is not None:
        _dnakit_arguments['orpha_id'] = orpha_id
    return _call('orphanet_get_classification', _dnakit_arguments, verbose=verbose)

def orphanet_get_disease(*,
    operation: str | None = None,
    orpha_code: str | int | None = None,
    lang: str | None = None,
    orphacode: int | str | None = None,
    orpha_id: int | str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get detailed information about a rare disease by ORPHA code. Returns disease name, definition, synonyms, and classification. ORPHA codes are unique Orphanet identifiers (e.g., 558 for Marfan syndrome). No authentication required.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if orpha_code is not None:
        _dnakit_arguments['orpha_code'] = orpha_code
    if lang is not None:
        _dnakit_arguments['lang'] = lang
    if orphacode is not None:
        _dnakit_arguments['orphacode'] = orphacode
    if orpha_id is not None:
        _dnakit_arguments['orpha_id'] = orpha_id
    return _call('orphanet_get_disease', _dnakit_arguments, verbose=verbose)

def orphanet_get_epidemiology(*,
    operation: str | None = None,
    orpha_code: str | int | None = None,
    orphacode: int | str | None = None,
    orpha_id: int | str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get epidemiology data (prevalence, incidence) for a rare disease from Orphanet. Returns prevalence estimates by geographic region and type (point prevalence, birth prevalence, annual incidence). Data includes validation status and literature sources.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if orpha_code is not None:
        _dnakit_arguments['orpha_code'] = orpha_code
    if orphacode is not None:
        _dnakit_arguments['orphacode'] = orphacode
    if orpha_id is not None:
        _dnakit_arguments['orpha_id'] = orpha_id
    return _call('orphanet_get_epidemiology', _dnakit_arguments, verbose=verbose)

def orphanet_get_gene_diseases(*,
    operation: str | None = None,
    gene_name: str | None = None,
    gene_symbol: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Get rare diseases associated with a gene from Orphanet. Accepts gene symbols (e.g., 'FBN1', 'BRCA1', 'HTT') or gene name keywords (e.g., 'fibrillin', 'huntingtin', 'dystrophin'). Gene symbols are auto-resolved to full names. Returns disease ORPHA codes, gene symbols, association types, and genomic locus information.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if gene_name is not None:
        _dnakit_arguments['gene_name'] = gene_name
    if gene_symbol is not None:
        _dnakit_arguments['gene_symbol'] = gene_symbol
    return _call('orphanet_get_gene_diseases', _dnakit_arguments, verbose=verbose)

def orphanet_get_genes(*,
    operation: str | None = None,
    orpha_code: str | int | None = None,
    orphacode: int | str | None = None,
    orpha_id: int | str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get genes associated with a rare disease from Orphanet. Returns gene symbols, names, and association types (causative, modifier, susceptibility). Essential for rare disease genetic diagnosis.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if orpha_code is not None:
        _dnakit_arguments['orpha_code'] = orpha_code
    if orphacode is not None:
        _dnakit_arguments['orphacode'] = orphacode
    if orpha_id is not None:
        _dnakit_arguments['orpha_id'] = orpha_id
    return _call('orphanet_get_genes', _dnakit_arguments, verbose=verbose)

def orphanet_get_icd_mapping(*,
    operation: str | None = None,
    orpha_code: str | int | None = None,
    coding_system: str | None = None,
    lang: str | None = None,
    orphacode: int | str | None = None,
    orpha_id: int | str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get cross-references between Orphanet ORPHA codes and other medical coding systems (ICD-10, ICD-11, OMIM, SNOMED-CT). Essential for mapping rare disease codes across clinical and research nomenclatures. Returns mapping relationships and validation status.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if orpha_code is not None:
        _dnakit_arguments['orpha_code'] = orpha_code
    if coding_system is not None:
        _dnakit_arguments['coding_system'] = coding_system
    if lang is not None:
        _dnakit_arguments['lang'] = lang
    if orphacode is not None:
        _dnakit_arguments['orphacode'] = orphacode
    if orpha_id is not None:
        _dnakit_arguments['orpha_id'] = orpha_id
    return _call('orphanet_get_icd_mapping', _dnakit_arguments, verbose=verbose)

def orphanet_get_natural_history(*,
    operation: str | None = None,
    orpha_code: str | int | None = None,
    orphacode: int | str | None = None,
    orpha_id: int | str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get natural history data for a rare disease from Orphanet. Returns average age of onset and type of inheritance. Useful for understanding disease progression patterns and genetic counseling.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if orpha_code is not None:
        _dnakit_arguments['orpha_code'] = orpha_code
    if orphacode is not None:
        _dnakit_arguments['orphacode'] = orphacode
    if orpha_id is not None:
        _dnakit_arguments['orpha_id'] = orpha_id
    return _call('orphanet_get_natural_history', _dnakit_arguments, verbose=verbose)

def orphanet_get_phenotypes(*,
    operation: str | None = None,
    orpha_code: str | int | None = None,
    orpha_id: str | int | None = None,
    orphacode: int | str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get HPO (Human Phenotype Ontology) phenotypes associated with a rare disease from Orphanet. Returns HPO terms with their frequency (Very frequent, Frequent, Occasional, Very rare, Excluded) and whether they are diagnostic criteria. Essential for phenotype-driven rare disease diagnosis.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if orpha_code is not None:
        _dnakit_arguments['orpha_code'] = orpha_code
    if orpha_id is not None:
        _dnakit_arguments['orpha_id'] = orpha_id
    if orphacode is not None:
        _dnakit_arguments['orphacode'] = orphacode
    return _call('orphanet_get_phenotypes', _dnakit_arguments, verbose=verbose)

def orphanet_search_by_name(*,
    operation: str | None = None,
    name: str | None = None,
    exact: bool | None = None,
    lang: str | None = None,
    limit: int | None = None,
    query: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Search for rare diseases by exact or partial name match in Orphanet. More precise than keyword search. Returns ORPHA codes and disease names. Supports exact matching mode for disambiguation.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if name is not None:
        _dnakit_arguments['name'] = name
    if exact is not None:
        _dnakit_arguments['exact'] = exact
    if lang is not None:
        _dnakit_arguments['lang'] = lang
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if query is not None:
        _dnakit_arguments['query'] = query
    return _call('orphanet_search_by_name', _dnakit_arguments, verbose=verbose)

def orphanet_search_diseases(*,
    operation: str | None = None,
    query: str,
    limit: int | None = None,
    lang: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Search Orphanet for rare diseases by keyword. Orphanet is the reference portal for rare diseases and orphan drugs. Returns disease names, ORPHA codes, and synonyms. No authentication required. Supports multiple languages (en, fr, de, es, it, pt, pl, nl).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['query'] = query
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if lang is not None:
        _dnakit_arguments['lang'] = lang
    return _call('orphanet_search_diseases', _dnakit_arguments, verbose=verbose)

def panther_enrichment(*,
    gene_list: str,
    organism: int | None | None = None,
    annotation_dataset: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Perform gene set enrichment (overrepresentation) analysis using PANTHER. Given a gene list, identifies statistically enriched GO terms, PANTHER pathways, Reactome pathways, or protein classes using Fisher's exact test with FDR correction. PANTHER is one of the most widely used tools for functional enrichment in bioinformatics. Example: gene list 'TP53,BRCA1,EGFR,KRAS,PTEN' enriched for GO biological_process returns 'regulation of mitotic cell cycle phase transition' (FDR=0.003).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_list'] = gene_list
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if annotation_dataset is not None:
        _dnakit_arguments['annotation_dataset'] = annotation_dataset
    return _call('panther_enrichment', _dnakit_arguments, verbose=verbose)

def panther_ortholog(*,
    gene_id: str,
    organism: int | None | None = None,
    target_organism: int | None | None = None,
    ortholog_type: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Find orthologs of a gene across species using PANTHER evolutionary classification. Returns EVERY matching ortholog in 'mappings' (count in 'total_mappings'), each with its persistent ID, gene symbol, and ortholog type (LDO=Least Diverged Ortholog, O=ortholog, P=paralog, X=xenolog, LDX=least diverged xenolog); 'mapping' repeats the first entry for backward compatibility and is not the whole answer. A gene often has several orthologs in one target species -- human ABCB1 maps to both mouse Abcb1a (LDO) and Abcb1b (O) -- so read 'mappings', especially with ortholog_type='O'. PANTHER classifies orthologs by evolutionary divergence across 144 fully sequenced genomes. Example: P04637 (human TP53) maps to Tp53 (mouse P02340) as LDO.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_id'] = gene_id
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if target_organism is not None:
        _dnakit_arguments['target_organism'] = target_organism
    if ortholog_type is not None:
        _dnakit_arguments['ortholog_type'] = ortholog_type
    return _call('panther_ortholog', _dnakit_arguments, verbose=verbose)

def pdbe_sifts_get_best_structures(*,
    uniprot_accession: str,
    limit: int | None | None = None,
    offset: int | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Get the best available PDB structures for a UniProt protein, ranked by coverage and resolution from PDBe SIFTS. Returns PDB IDs with chain mapping, resolution, experimental method (X-ray, cryo-EM, NMR), and UniProt sequence coverage. Essential for structural biologists to find the most suitable structure for a protein of interest. Example: P04637 (TP53) returns 641 structure entries, P00533 (EGFR) returns 633 entries, with X-ray and cryo-EM structures ranked by quality. Returns 50 entries per call by default; use 'limit' (max 500) and 'offset' to page through the full ranked list, whose size is reported in the response.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['uniprot_accession'] = uniprot_accession
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    return _call('pdbe_sifts_get_best_structures', _dnakit_arguments, verbose=verbose)

def pdbe_get_uniprot_mappings(*,
    pdb_id: str | None = None,
    uniprot_id: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get UniProt-to-PDB chain mappings for a PDB structure from the PDBe Graph API. For each UniProt accession represented in the structure, returns all chain mappings with residue-level alignment (UniProt start/end, PDB start/end, identity, coverage). Essential for understanding which proteins are in a crystal structure and how they map to UniProt reference sequences. Example: 4HHB returns P69905 (HBA_HUMAN) mapped to chains A,C and P68871 (HBB_HUMAN) mapped to chains B,D.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if pdb_id is not None:
        _dnakit_arguments['pdb_id'] = pdb_id
    if uniprot_id is not None:
        _dnakit_arguments['uniprot_id'] = uniprot_id
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('pdbe_get_uniprot_mappings', _dnakit_arguments, verbose=verbose)

def powo_search_plants(*,
    q: str,
    f: str | None | None = None,
    perPage: int | None | None = None,
    cursor: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    'Search the Plants of the World Online (POWO) database maintained by the Royal Botanic Gardens, Kew. POWO is the most comprehensive global checklist of vascular plants with 1.3+ million names and 350,000+ accepted species. Search by plant name, genus, family, or common name. Returns taxonomic information including accepted names, synonyms, distribution, and images. Essential for plant taxonomy, botany, and biodiversity research.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['q'] = q
    if f is not None:
        _dnakit_arguments['f'] = f
    if perPage is not None:
        _dnakit_arguments['perPage'] = perPage
    if cursor is not None:
        _dnakit_arguments['cursor'] = cursor
    return _call('powo_search_plants', _dnakit_arguments, verbose=verbose)

def pfam_get_protein_annotations(*,
    accession: str,
    verbose: bool = True,
) -> _Object:
    'Get all Pfam domain annotations for a specific protein by its UniProt accession. Returns every Pfam domain found in the protein with exact start/end residue positions, E-values (scores), short names, and integrated InterPro entries. This provides the complete Pfam domain architecture of a protein. For example, P04637 (TP53) has 4 Pfam domains: PF08563 (P53_TAD, 6-30), PF18521 (P53_TAD2, 35-59), PF00870 (P53, 100-288), PF07710 (P53_tetramer, 319-357). Use this to understand the domain composition of any protein.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['accession'] = accession
    return _call('pfam_get_protein_annotations', _dnakit_arguments, verbose=verbose)

def plant_reactome_get_pathway(*,
    pathway_id: str,
    verbose: bool = True,
) -> _Object:
    "Get detailed information about a specific Plant Reactome pathway by stable ID. Returns pathway name, species, sub-events (reactions), compartments, GO annotations, and cross-references. Example: 'R-OSA-1119312' returns the Photorespiration pathway in Oryza sativa.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['pathway_id'] = pathway_id
    return _call('plant_reactome_get_pathway', _dnakit_arguments, verbose=verbose)

def plant_reactome_list_species(*,
    verbose: bool = True,
) -> _Object:
    'List all plant species available in Plant Reactome with their taxonomy IDs. Returns the complete catalog of 140+ plant species including model organisms and major crop species. Useful for discovering which species have curated pathway data available.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    return _call('plant_reactome_list_species', _dnakit_arguments, verbose=verbose)

def plant_reactome_search_pathways(*,
    query: str,
    species: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Search for plant biological pathways in Plant Reactome (Gramene). Plant Reactome is a curated resource for plant metabolic and regulatory pathways, covering model and crop species including rice (Oryza sativa), Arabidopsis thaliana, maize, wheat, and 140+ plant species. Returns pathways matching the query with their stable IDs and species. Example: 'photosynthesis' finds photosynthesis-related pathways.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if species is not None:
        _dnakit_arguments['species'] = species
    return _call('plant_reactome_search_pathways', _dnakit_arguments, verbose=verbose)

def progenetix_cnv_search(*,
    reference_name: str,
    start: int,
    end: int,
    variant_type: str | None = None,
    filters: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search for cancer biosamples with copy number variations (CNVs) in a specific genomic region. Combines disease filter with genomic coordinates to find samples with amplifications or deletions in genes of interest. Uses RefSeq chromosome accessions for coordinates. Example: searching for EGFR amplifications (refseq:NC_000007.14:55019017-55211628, DUP) in breast cancer (NCIT:C4017) returns 9 matching biosamples.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['reference_name'] = reference_name
    _dnakit_arguments['start'] = start
    _dnakit_arguments['end'] = end
    if variant_type is not None:
        _dnakit_arguments['variant_type'] = variant_type
    if filters is not None:
        _dnakit_arguments['filters'] = filters
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('progenetix_cnv_search', _dnakit_arguments, verbose=verbose)

def progenetix_list_cohorts(*,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'List available cohorts in the Progenetix cancer CNV database. Returns cohort IDs, names, sizes, and data types. Progenetix organizes samples into cohorts like TCGA, DIPG studies, and curated cancer sets. Example: returns 6 cohorts including TCGA (517 samples) and DIPG (16 samples).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('progenetix_list_cohorts', _dnakit_arguments, verbose=verbose)

def progenetix_list_filtering_terms(*,
    prefixes: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    "List available filtering terms (ontology codes) in the Progenetix cancer CNV database. Returns NCIt cancer type codes, EDAM operation codes, or other ontology terms with the number of matching samples. Useful for discovering which cancer types are available before querying biosamples or individuals. Example: prefixes='NCIT' returns 890 NCIt terms across all cancer types.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if prefixes is not None:
        _dnakit_arguments['prefixes'] = prefixes
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('progenetix_list_filtering_terms', _dnakit_arguments, verbose=verbose)

def progenetix_search_biosamples(*,
    filters: str,
    limit: int | None = None,
    skip: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search cancer tumor biosamples in the Progenetix database by NCIt disease ontology code. Progenetix contains genome-wide copy number variation (CNV) profiles from over 100,000 tumor samples. Returns biosample IDs, histological diagnosis, pathological stage, and external references. Uses the GA4GH Beacon v2 protocol. Example: NCIT:C3058 returns 28 glioblastoma biosamples with histological classification and study references.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['filters'] = filters
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if skip is not None:
        _dnakit_arguments['skip'] = skip
    return _call('progenetix_search_biosamples', _dnakit_arguments, verbose=verbose)

def pub_chem_get_cid_by_compound_name(*,
    name: str | None = None,
    compound_name: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Search PubChem by compound name to get CID(s). This is the PRIMARY tool for looking up a drug or chemical by name (e.g., 'cytarabine', 'aspirin', 'metformin', 'venetoclax') to obtain PubChem Compound IDs. Use the returned CID with other PubChem tools for properties, patents, etc. Only works with chemical/drug names, not disease names.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if name is not None:
        _dnakit_arguments['name'] = name
    if compound_name is not None:
        _dnakit_arguments['compound_name'] = compound_name
    return _call('pub_chem_get_cid_by_compound_name', _dnakit_arguments, verbose=verbose)

def pub_chem_get_compound_properties_by_cid(*,
    cid: int,
    properties: list[_Object] | str | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Get a set of specified molecular properties through CID (Compound ID), such as molecular weight, IUPAC name, and SMILES. PubChem's PUG-REST renamed its SMILES properties: ask for 'ConnectivitySMILES' (formerly 'CanonicalSMILES') and 'SMILES' (formerly 'IsomericSMILES'). Requesting a legacy name returns the value under the CURRENT key, so a caller that looks up the key it asked for will not find it.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['cid'] = cid
    if properties is not None:
        _dnakit_arguments['properties'] = properties
    return _call('pub_chem_get_compound_properties_by_cid', _dnakit_arguments, verbose=verbose)

def pub_med_get_article(*,
    pmid: str | int,
    verbose: bool = True,
) -> _Object:
    'Get complete metadata for a specific PubMed article by its PMID (PubMed ID) using efetch. Returns comprehensive article information including title, full abstract, complete author list with affiliations, journal details (name, volume, issue, pages), publication dates, DOI, MeSH terms, article type, grant information, and references. Essential for accessing full abstracts, detailed metadata, and citation information. Find PMIDs using PubMed_search_articles (format: integer like 12345678).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['pmid'] = pmid
    return _call('pub_med_get_article', _dnakit_arguments, verbose=verbose)

def pub_med_get_cited_by(*,
    pmid: str,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get a list of PubMed articles that cite a specific PMID using elink. Returns PMIDs of articles that have cited the input article in their reference list. Use this to track citation impact, find more recent work building on a paper, or explore how research has been applied. Essential for literature reviews, impact analysis, and staying current with research developments. Find PMIDs using PubMed_search_articles.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['pmid'] = pmid
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('pub_med_get_cited_by', _dnakit_arguments, verbose=verbose)

def pub_med_get_related(*,
    pmid: str,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get related PubMed articles for a specific PMID using elink. Returns a list of computationally similar articles based on title, abstract, and MeSH terms. Use this to discover additional relevant literature, find similar studies, or explore a research topic more comprehensively. The algorithm identifies articles with similar content, methodology, or subject matter. Find PMIDs using PubMed_search_articles.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['pmid'] = pmid
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('pub_med_get_related', _dnakit_arguments, verbose=verbose)

def pub_med_search_articles(*,
    query: str,
    limit: int | None = None,
    mindate: str | None = None,
    maxdate: str | None = None,
    datetype: str | None = None,
    include_abstract: bool | None = None,
    sort: str | None = None,
    max_results: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search PubMed biomedical literature database using NCBI E-utilities (esearch + esummary). Returns articles with PMID, title, authors, journal, publication year, DOI, article type, and PubMed URL. Automatically fetches rich metadata for each article. Use this to discover medical/life sciences literature, find articles by topic, or search for publications by author or journal. Rate limiting is automatically handled (3 req/sec without API key, 10 req/sec with key set via NCBI_API_KEY environment variable).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if mindate is not None:
        _dnakit_arguments['mindate'] = mindate
    if maxdate is not None:
        _dnakit_arguments['maxdate'] = maxdate
    if datetype is not None:
        _dnakit_arguments['datetype'] = datetype
    if include_abstract is not None:
        _dnakit_arguments['include_abstract'] = include_abstract
    if sort is not None:
        _dnakit_arguments['sort'] = sort
    if max_results is not None:
        _dnakit_arguments['max_results'] = max_results
    return _call('pub_med_search_articles', _dnakit_arguments, verbose=verbose)

def quick_go_annotations_by_gene(*,
    gene_product_id: str,
    aspect: str | None | None = None,
    taxon_id: int | None | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Search Gene Ontology (GO) annotations for a specific gene product using the EBI QuickGO browser. Returns GO terms annotated to the protein/gene with evidence codes, qualifiers, and annotation sources. Accepts UniProtKB accessions (e.g., 'UniProtKB:P04637' for TP53). Provides more detailed annotation data than the standard GO API, including evidence provenance and assigned-by information. Useful for understanding gene function, localization, and biological processes.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_product_id'] = gene_product_id
    if aspect is not None:
        _dnakit_arguments['aspect'] = aspect
    if taxon_id is not None:
        _dnakit_arguments['taxon_id'] = taxon_id
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('quick_go_annotations_by_gene', _dnakit_arguments, verbose=verbose)

def rcsbdata_get_entry(*,
    pdb_id: str | None = None,
    entry_id: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Get comprehensive structure metadata for a PDB entry from the RCSB Data API. Returns experimental method, resolution, deposition/release dates, polymer/nonpolymer entity counts, molecular weight, assembly count, space group, and unit cell dimensions. This is the direct REST endpoint for PDB entry details, complementing GraphQL-based tools. Example: '4HHB' returns hemoglobin at 1.74A by X-ray diffraction, deposited 1984-03-07, 2 polymer entities.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if pdb_id is not None:
        _dnakit_arguments['pdb_id'] = pdb_id
    if entry_id is not None:
        _dnakit_arguments['entry_id'] = entry_id
    return _call('rcsbdata_get_entry', _dnakit_arguments, verbose=verbose)

def rgd_get_annotations(*,
    rgd_id: str,
    aspect: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get disease, phenotype, pathway, and GO annotations for a rat gene from RGD. Returns curated annotations with ontology terms, evidence codes, and qualifiers. Filter by aspect: D=Disease, P=Pathway, W=Phenotype, F=Molecular Function, C=Cellular Component, E=Expression. Essential for rat disease model research and toxicogenomics.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['rgd_id'] = rgd_id
    if aspect is not None:
        _dnakit_arguments['aspect'] = aspect
    return _call('rgd_get_annotations', _dnakit_arguments, verbose=verbose)

def rgd_get_gene(*,
    rgd_id: str,
    verbose: bool = True,
) -> _Object:
    'Get detailed rat gene information from the Rat Genome Database by RGD ID. Returns gene symbol, name, description (including AGR-aggregated functional summary with disease associations, pathways, and interactions), gene type, and Ensembl cross-reference. Use when looking up specific rat genes for toxicology, pharmacology, or disease model research.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['rgd_id'] = rgd_id
    return _call('rgd_get_gene', _dnakit_arguments, verbose=verbose)

def rgd_get_orthologs(*,
    rgd_id: str,
    verbose: bool = True,
) -> _Object:
    'Get orthologs of a rat gene across species (human, mouse, chinchilla, dog, pig, etc.) from RGD. Returns ortholog gene symbols, RGD IDs, names, and species. Use for cross-species comparison in model organism genetics, finding human disease gene orthologs in rat.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['rgd_id'] = rgd_id
    return _call('rgd_get_orthologs', _dnakit_arguments, verbose=verbose)

def rgd_get_qtls_in_region(*,
    chromosome: str,
    start: int,
    stop: int,
    map_key: int,
    verbose: bool = True,
) -> _Object:
    "Get QTLs (Quantitative Trait Loci) overlapping a genomic region from the Rat Genome Database (RGD's flagship data type). Returns QTL symbol, name, RGD ID, LOD score, p-value, and inheritance type. Requires a chromosome, start/stop coordinates, and a map_key identifying the assembly (e.g., 360 for rat rn7/GRCr8, 38 for human GRCh38, 35 for mouse GRCm39). Example: chromosome='10', start=1, stop=50000000, map_key=360 returns QTL Uae41 'Urinary albumin excretion QTL 41' lod=5.26 rgdId=7387235.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['chromosome'] = chromosome
    _dnakit_arguments['start'] = start
    _dnakit_arguments['stop'] = stop
    _dnakit_arguments['map_key'] = map_key
    return _call('rgd_get_qtls_in_region', _dnakit_arguments, verbose=verbose)

def rgd_get_strain(*,
    symbol: str,
    verbose: bool = True,
) -> _Object:
    "Get the full Rat Genome Database (RGD) record for a single rat strain by its exact strain symbol (e.g. 'SHR', 'BN', 'SS', 'WKY', 'GK'). Returns the strain RGD ID, name, full strain name, substrain, strain type, genetics formula, origin, background strain, and description. Use RGD_search_strains first if you only know a keyword.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['symbol'] = symbol
    return _call('rgd_get_strain', _dnakit_arguments, verbose=verbose)

def rgd_get_strain_annotations(*,
    symbol: str,
    category: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Get curated disease, phenotype, and variant-trait ontology annotations for a rat strain from the Rat Genome Database (RGD), by strain symbol (e.g. 'SHR', 'BN', 'SS'). This exposes RGD's rat disease-MODEL annotations - for example SHR is annotated to 'Left Ventricular Hypertrophy' (disease, MODEL: spontaneous) and 'arterial blood pressure trait' (phenotype). Each annotation includes the ontology term, accession (DOID/VT/CMO/MP/RS), category, qualifier, and evidence code. Optional 'category' filter: disease, phenotype, strain, pathway, biological_process, molecular_function, cellular_component.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['symbol'] = symbol
    if category is not None:
        _dnakit_arguments['category'] = category
    return _call('rgd_get_strain_annotations', _dnakit_arguments, verbose=verbose)

def rgd_search_genes(*,
    query: str,
    species: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Search for rat genes by symbol or keyword in the Rat Genome Database. Returns matching genes with RGD IDs, symbols, names, types, and descriptions. Supports rat, human, and mouse species. Use for finding rat orthologs, toxicology studies, or rat disease model genes.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if species is not None:
        _dnakit_arguments['species'] = species
    return _call('rgd_search_genes', _dnakit_arguments, verbose=verbose)

def rgd_search_strains(*,
    query: str | None | None = None,
    strain_type: str | None | None = None,
    limit: int | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Search the Rat Genome Database (RGD) rat strain catalog by keyword and/or strain type. Rat inbred/congenic/transgenic strains are RGD's flagship rat-specific resource (not in Alliance of Genome Resources or MGI). Provide 'query' (a keyword matched against symbol, name, full strain name, genetics, and description) and/or 'strain_type' (one of: inbred, congenic, mutant, transgenic, recombinant_inbred, outbred, consomic, coisogenic, segregating_inbred, hybrid, wild, conplastic, advanced_intercross_line, transchromosomal). At least one of the two must be supplied. Returns matching strains with their RGD ID, symbol, name, type, genetics, and origin.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if query is not None:
        _dnakit_arguments['query'] = query
    if strain_type is not None:
        _dnakit_arguments['strain_type'] = strain_type
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('rgd_search_strains', _dnakit_arguments, verbose=verbose)

def re_map_get_transcription_factor_binding(*,
    gene_name: str,
    cell_type: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get transcription factor binding sites from ReMap database for specific genes and cell types\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_name'] = gene_name
    if cell_type is not None:
        _dnakit_arguments['cell_type'] = cell_type
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('re_map_get_transcription_factor_binding', _dnakit_arguments, verbose=verbose)

def reactome_analysis_pathway_enrichment(*,
    identifiers: str,
    page_size: int | None | None = None,
    include_disease: bool | None | None = None,
    projection: bool | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Perform pathway overrepresentation (enrichment) analysis using the Reactome Analysis Service. Submit a list of gene symbols, UniProt IDs, or Ensembl IDs and identify which Reactome pathways are statistically enriched. Returns pathways ranked by p-value with FDR (False Discovery Rate) correction, entity counts, and reaction coverage. Supports projection to human pathways for cross-species analysis. INTERPRETING RESULTS: rank pathways by 'fdr' or 'p_value' - those are the enrichment statistics. 'entities_coverage' (= entities_found / entities_total) is the fraction of a pathway your identifiers cover. 'pathway_size_fraction_of_reactome' (deprecated alias 'entities_ratio') is Reactome's raw 'ratio': pathway size over Reactome's entity universe for the species, independent of the identifiers you submit - never rank by it. Example: submitting TP53, BRCA1, CDH1, EGFR, KRAS returns 394 enriched pathways including 'Developmental Lineage of Mammary Gland Myoepithelial Cells' (p=3.94e-08) and 'Signaling by ERBB2 KD Mutants'.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['identifiers'] = identifiers
    if page_size is not None:
        _dnakit_arguments['page_size'] = page_size
    if include_disease is not None:
        _dnakit_arguments['include_disease'] = include_disease
    if projection is not None:
        _dnakit_arguments['projection'] = projection
    return _call('reactome_analysis_pathway_enrichment', _dnakit_arguments, verbose=verbose)

def reactome_content_search(*,
    query: str,
    species: str | None = None,
    types: str | None = None,
    rows: int | None | None = None,
    start: int | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Search Reactome for pathways, reactions, and biological entities by keyword. Performs free-text search across Reactome's knowledge base of biological pathways. Supports filtering by species and entity type. Example: searching 'apoptosis' in Homo sapiens returns pathways like R-HSA-109581 (Apoptosis), R-HSA-169911 (Regulation of Apoptosis), R-HSA-9635465 (Suppression of apoptosis). Use this to discover relevant pathways before querying their details. Results are paged: 'count' is how many entries came back, 'total_results' is how many matched in Reactome overall, and 'has_more' flags that the list is a slice -- use 'rows' and 'start' to page through the rest.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if species is not None:
        _dnakit_arguments['species'] = species
    if types is not None:
        _dnakit_arguments['types'] = types
    if rows is not None:
        _dnakit_arguments['rows'] = rows
    if start is not None:
        _dnakit_arguments['start'] = start
    return _call('reactome_content_search', _dnakit_arguments, verbose=verbose)

def reactome_get_participants(*,
    stId: str,
    verbose: bool = True,
) -> _Object:
    'Get all participants (physical entities) of a reaction or event. Returns list of physical entities (proteins, complexes, small molecules) participating in the event, including their reference entities.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['stId'] = stId
    return _call('reactome_get_participants', _dnakit_arguments, verbose=verbose)

def reactome_get_pathway(*,
    stId: str,
    verbose: bool = True,
) -> _Object:
    'Get detailed information about a specific pathway by its Stable ID. Returns comprehensive pathway data including metadata, compartments, literature references, orthologous events, and related pathways. Uses the unified query endpoint that works for all entity types.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['stId'] = stId
    return _call('reactome_get_pathway', _dnakit_arguments, verbose=verbose)

def reactome_get_pathway_reactions(*,
    stId: str,
    verbose: bool = True,
) -> _Object:
    'Query all reactions and subpathways contained under a Pathway using Pathway Stable ID. Returns list of events (reactions and subpathways) contained in the pathway, including their metadata and relationships.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['stId'] = stId
    return _call('reactome_get_pathway_reactions', _dnakit_arguments, verbose=verbose)

def reactome_map_uniprot_to_pathways(*,
    uniprot_id: str | None = None,
    id: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Map a UniProt protein identifier to Reactome pathways. Returns all pathways that contain this protein.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if uniprot_id is not None:
        _dnakit_arguments['uniprot_id'] = uniprot_id
    if id is not None:
        _dnakit_arguments['id'] = id
    return _call('reactome_map_uniprot_to_pathways', _dnakit_arguments, verbose=verbose)

def regulome_db_query_variant(*,
    rsid: str | None = None,
    variant: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Query RegulomeDB for regulatory variant annotations using rsID\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if rsid is not None:
        _dnakit_arguments['rsid'] = rsid
    if variant is not None:
        _dnakit_arguments['variant'] = variant
    return _call('regulome_db_query_variant', _dnakit_arguments, verbose=verbose)

def screen_get_regulatory_elements(*,
    gene_name: str | None = None,
    query: str | None = None,
    element_type: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get candidate cis-regulatory elements (cCREs) from SCREEN database for specific genomic regions\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if gene_name is not None:
        _dnakit_arguments['gene_name'] = gene_name
    if query is not None:
        _dnakit_arguments['query'] = query
    if element_type is not None:
        _dnakit_arguments['element_type'] = element_type
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('screen_get_regulatory_elements', _dnakit_arguments, verbose=verbose)

def sgd_get_gene(*,
    sgd_id: str,
    verbose: bool = True,
) -> _Object:
    'Get yeast gene/locus overview from SGD (Saccharomyces Genome Database) by SGD identifier. Returns gene name, systematic name, description, UniProt ID, aliases, and functional annotations. S. cerevisiae is the primary model organism for eukaryotic cell biology. Example: S000000259 returns GAL4, the galactose metabolism transcription factor.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['sgd_id'] = sgd_id
    return _call('sgd_get_gene', _dnakit_arguments, verbose=verbose)

def sgd_get_go_annotations(*,
    sgd_id: str,
    verbose: bool = True,
) -> _Object:
    'Get Gene Ontology annotations for a yeast gene from SGD. Returns GO terms (biological process, molecular function, cellular component), evidence codes, qualifiers, and references. SGD GO annotations are among the most thoroughly curated in any model organism. Example: S000003219 (RMR1) has annotations for meiotic recombination and cytoplasm localization.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['sgd_id'] = sgd_id
    return _call('sgd_get_go_annotations', _dnakit_arguments, verbose=verbose)

def sgd_get_interactions(*,
    sgd_id: str,
    verbose: bool = True,
) -> _Object:
    'Get genetic and physical interactions for a yeast gene from SGD. Returns interaction type (Physical/Genetic), experiment type (Two-hybrid, Affinity Capture-MS, Synthetic Lethality, etc.), interacting gene names, and references. Yeast interaction data is the most comprehensive of any organism. Example: S000003219 (RMR1) has 178 interactions.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['sgd_id'] = sgd_id
    return _call('sgd_get_interactions', _dnakit_arguments, verbose=verbose)

def sgd_get_phenotypes(*,
    sgd_id: str,
    verbose: bool = True,
) -> _Object:
    'Get curated phenotype annotations for a yeast gene from SGD. Returns observable phenotypes (growth, morphology, fitness), mutant types, experimental methods, strain backgrounds, and chemical conditions. Yeast phenotypes are a key resource for functional genomics. Example: S000003219 (RMR1) has 29 phenotype annotations including competitive fitness, respiratory growth defects.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['sgd_id'] = sgd_id
    return _call('sgd_get_phenotypes', _dnakit_arguments, verbose=verbose)

def sgd_search(*,
    query: str,
    category: str | None = None,
    limit: int | None = None,
    offset: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Search SGD across all data types (genes, GO terms, phenotypes, references, alleles, diseases, complexes, pathways). Returns matching entries with categories, descriptions, and links. Supports category filtering. Example: searching 'GAL4' finds the transcription factor gene, related references, alleles, and GO terms.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if category is not None:
        _dnakit_arguments['category'] = category
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    return _call('sgd_search', _dnakit_arguments, verbose=verbose)

def sra_get_experiment(*,
    accession: str,
    verbose: bool = True,
) -> _Object:
    'Get detailed metadata for a specific SRA experiment by accession number. Returns experiment title, organism, sequencing platform, library information, study details, and run statistics. Accepts SRA experiment accessions (SRX/ERX/DRX), study accessions (SRP/ERP/DRP), or sample accessions (SRS/ERS/DRS).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['accession'] = accession
    return _call('sra_get_experiment', _dnakit_arguments, verbose=verbose)

def sra_search_experiments(*,
    query: str | None = None,
    organism: str | None = None,
    library_strategy: str | None = None,
    platform: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search the NCBI Sequence Read Archive (SRA) for sequencing experiments by keyword, organism, library strategy, or sequencing platform. SRA is the largest repository of raw sequencing data, containing data from Illumina, PacBio, Oxford Nanopore, and other platforms. Library strategies include: RNA-Seq, WGS (whole genome), WXS (whole exome), ChIP-Seq, ATAC-seq, Bisulfite-Seq, AMPLICON, scRNA-Seq, and more.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if query is not None:
        _dnakit_arguments['query'] = query
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if library_strategy is not None:
        _dnakit_arguments['library_strategy'] = library_strategy
    if platform is not None:
        _dnakit_arguments['platform'] = platform
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('sra_search_experiments', _dnakit_arguments, verbose=verbose)

def string_functional_enrichment(*,
    protein_ids: list[str],
    species: int | None = None,
    category: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Identify enriched biological functions, pathways, and processes for a protein set using STRING (Search Tool for Retrieval of Interacting Genes/Proteins) database. Performs statistical enrichment analysis against GO (Gene Ontology) terms, KEGG pathways, Reactome pathways, and other annotation databases. Returns significantly enriched terms with FDR-corrected p-values. Minimum 3-4 proteins required for meaningful statistical analysis. No API key required. Use for: discovering pathways enriched in protein network, identifying biological processes, finding shared functions among proteins, interpreting omics data (proteomics, interactomics), hypothesis generation.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['protein_ids'] = protein_ids
    if species is not None:
        _dnakit_arguments['species'] = species
    if category is not None:
        _dnakit_arguments['category'] = category
    return _call('string_functional_enrichment', _dnakit_arguments, verbose=verbose)

def string_get_functional_annotations(*,
    identifiers: str,
    species: int | None = None,
    category: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Get comprehensive functional annotations for a protein from STRING database. Returns per-protein annotations from Gene Ontology (Process, Function, Component), KEGG pathways, Reactome, WikiPathways, InterPro domains, Pfam, HPO phenotypes, DISEASES, TISSUES expression, and COMPARTMENTS subcellular localization. Different from STRING_functional_enrichment which tests for statistical enrichment of a protein SET - this tool returns all known annotations for individual proteins. Example: 'TP53' returns 1,877 annotations across 13 categories including 655 biological processes, 462 HPO phenotypes, 212 disease associations.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['identifiers'] = identifiers
    if species is not None:
        _dnakit_arguments['species'] = species
    if category is not None:
        _dnakit_arguments['category'] = category
    return _call('string_get_functional_annotations', _dnakit_arguments, verbose=verbose)

def string_get_interaction_partners(*,
    identifiers: str | None = None,
    species: int | None = None,
    limit: int | None = None,
    required_score: int | None | None = None,
    protein_ids: list[_Object] | str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get interaction partners for a single protein from STRING database. Returns a ranked list of proteins that interact with the query protein, along with detailed confidence scores by evidence type (experimental, database, text-mining, etc.). Different from STRING_get_network which returns all edges in a network - this returns a simple partner list focused on one protein. Useful for finding what proteins a given protein interacts with.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if identifiers is not None:
        _dnakit_arguments['identifiers'] = identifiers
    if species is not None:
        _dnakit_arguments['species'] = species
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if required_score is not None:
        _dnakit_arguments['required_score'] = required_score
    if protein_ids is not None:
        _dnakit_arguments['protein_ids'] = protein_ids
    return _call('string_get_interaction_partners', _dnakit_arguments, verbose=verbose)

def string_get_network(*,
    identifiers: str | None = None,
    species: int | None = None,
    limit: int | None = None,
    required_score: int | None | None = None,
    protein_ids: list[_Object] | str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get the protein-protein interaction network for one or more proteins from STRING database. Returns all interactions (edges) between the query proteins and their neighbors, including confidence scores broken down by evidence type (experiments, co-expression, databases, text-mining, etc.). Different from STRING_get_functional_annotations (which returns GO/pathway annotations) and STRING_get_enrichment (which tests enrichment of a set). Use this to build or visualize PPI networks. Multiple proteins can be separated by carriage return characters (\\r).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if identifiers is not None:
        _dnakit_arguments['identifiers'] = identifiers
    if species is not None:
        _dnakit_arguments['species'] = species
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if required_score is not None:
        _dnakit_arguments['required_score'] = required_score
    if protein_ids is not None:
        _dnakit_arguments['protein_ids'] = protein_ids
    return _call('string_get_network', _dnakit_arguments, verbose=verbose)

def string_ppi_enrichment(*,
    protein_ids: list[str],
    species: int | None = None,
    confidence_score: float | None = None,
    verbose: bool = True,
) -> _Object:
    'Test if your protein set has more interactions than expected by chance (PPI = Protein-Protein Interaction enrichment). Compares observed interactions in your network to random expectation - significant enrichment suggests proteins work together functionally. Returns enrichment p-value and expected vs observed interaction counts. Useful for validating that your protein list forms a real functional module rather than random proteins. No API key required. Use for: validating protein complex predictions, testing if proteins form functional module, quality control for network analysis, distinguishing real complexes from random protein lists, assessing functional coherence.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['protein_ids'] = protein_ids
    if species is not None:
        _dnakit_arguments['species'] = species
    if confidence_score is not None:
        _dnakit_arguments['confidence_score'] = confidence_score
    return _call('string_ppi_enrichment', _dnakit_arguments, verbose=verbose)

def semantic_scholar_search_papers(*,
    query: str,
    limit: int | None = None,
    year: str | None = None,
    sort: str | None = None,
    include_abstract: bool | None = None,
    verbose: bool = True,
) -> _Object:
    "Search for papers on Semantic Scholar including abstracts and AI-generated TLDR summaries. This tool queries the Semantic Scholar API using natural language keywords and returns papers with details such as title, abstract, TLDR summary, fields of study, publication year, journal (venue), and URL. Returns only one page of results: `metadata.returned` is how many papers are in `data`, `metadata.total` is Semantic Scholar's approximate count of ALL papers matching the query in the corpus (null when the API does not report it, see `metadata.total_source`), and the top-level `truncated` flag plus `truncation_note` say when more matches exist than were returned. Do not treat `metadata.returned` as the size of the literature. Rate limits: 1 req/sec without API key, 100 req/sec with SEMANTIC_SCHOLAR_API_KEY environment variable set.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if year is not None:
        _dnakit_arguments['year'] = year
    if sort is not None:
        _dnakit_arguments['sort'] = sort
    if include_abstract is not None:
        _dnakit_arguments['include_abstract'] = include_abstract
    return _call('semantic_scholar_search_papers', _dnakit_arguments, verbose=verbose)

def sequence_count_residues(*,
    operation: str | None = None,
    sequence: str | None = None,
    uniprot_id: str | None = None,
    residue: str,
    verbose: bool = True,
) -> _Object:
    'Count occurrences of a specific amino acid or nucleotide residue in a sequence. Accepts either a raw sequence string or a UniProt accession (fetches protein sequence from UniProt REST API). Reports count, fraction, percentage, and 1-based positions. Useful for cysteine counting (disulfide bond estimation), charged residue analysis, or nucleotide composition.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if sequence is not None:
        _dnakit_arguments['sequence'] = sequence
    if uniprot_id is not None:
        _dnakit_arguments['uniprot_id'] = uniprot_id
    _dnakit_arguments['residue'] = residue
    return _call('sequence_count_residues', _dnakit_arguments, verbose=verbose)

def sequence_dn_ds(*,
    seq1: str | None | None = None,
    seq2: str | None | None = None,
    fasta1_path: str | None | None = None,
    fasta2_path: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Compute dN/dS (Ka/Ks) between two coding sequences via the Nei-Gojobori (1986) estimator with Jukes-Cantor correction — the standard test for selection: dN/dS > 1 = positive/diversifying, << 1 = purifying, ~ 1 = near-neutral. Give two in-frame, codon-aligned coding sequences inline ('seq1'/'seq2') or as single-record FASTA files. Returns dN, dS, the ratio, per-site counts, and an interpretation. Pure Python, no API key. Works for any organisms/genes.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if seq1 is not None:
        _dnakit_arguments['seq1'] = seq1
    if seq2 is not None:
        _dnakit_arguments['seq2'] = seq2
    if fasta1_path is not None:
        _dnakit_arguments['fasta1_path'] = fasta1_path
    if fasta2_path is not None:
        _dnakit_arguments['fasta2_path'] = fasta2_path
    return _call('sequence_dn_ds', _dnakit_arguments, verbose=verbose)

def sequence_gc_content(*,
    operation: str | None = None,
    sequence: str,
    verbose: bool = True,
) -> _Object:
    'Calculate GC content of a DNA or RNA sequence. Reports GC count, GC fraction, GC percentage, and full nucleotide composition. Ambiguous bases (N) are excluded from the denominator. Useful for primer design, genome characterization, and codon usage analysis.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['sequence'] = sequence
    return _call('sequence_gc_content', _dnakit_arguments, verbose=verbose)

def sequence_reverse_complement(*,
    operation: str | None = None,
    sequence: str,
    verbose: bool = True,
) -> _Object:
    'Return the reverse complement of a DNA sequence. Supports A, T, C, G, and N bases. Useful for primer design, probe design, and strand orientation analysis.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['sequence'] = sequence
    return _call('sequence_reverse_complement', _dnakit_arguments, verbose=verbose)

def sequence_stats(*,
    operation: str | None = None,
    sequence: str | None = None,
    uniprot_id: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Compute basic statistics for a DNA, RNA, or protein sequence. Auto-detects sequence type, reports length, full composition, and type-specific metrics (GC% for nucleotides, estimated monoisotopic molecular weight for proteins -- for the conventional average mass reported by ExPASy ProtParam/UniProt, use ProtParam_calculate instead). Accepts sequence string or UniProt accession.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if sequence is not None:
        _dnakit_arguments['sequence'] = sequence
    if uniprot_id is not None:
        _dnakit_arguments['uniprot_id'] = uniprot_id
    return _call('sequence_stats', _dnakit_arguments, verbose=verbose)

def splice_ai_get_max_delta(*,
    variant: str,
    genome: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get simplified maximum SpliceAI delta score and pathogenicity interpretation for a variant. Returns single max score across all splice effects (AG/AL/DG/DL) with interpretation: benign (<0.2), low (0.2-0.5), moderate (0.5-0.8), high (≥0.8). Best for quick variant triage.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['variant'] = variant
    if genome is not None:
        _dnakit_arguments['genome'] = genome
    return _call('splice_ai_get_max_delta', _dnakit_arguments, verbose=verbose)

def splice_ai_predict_pangolin(*,
    variant: str,
    genome: str | None = None,
    distance: int | None = None,
    mask: bool | None = None,
    verbose: bool = True,
) -> _Object:
    'Predict splice-altering effects using Pangolin model (alternative to SpliceAI). Pangolin uses a different architecture and may provide complementary predictions. Returns splice effect scores per transcript. Useful for validation or when SpliceAI gives borderline results.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['variant'] = variant
    if genome is not None:
        _dnakit_arguments['genome'] = genome
    if distance is not None:
        _dnakit_arguments['distance'] = distance
    if mask is not None:
        _dnakit_arguments['mask'] = mask
    return _call('splice_ai_predict_pangolin', _dnakit_arguments, verbose=verbose)

def splice_ai_predict_splice(*,
    variant: str,
    genome: str | None = None,
    distance: int | None = None,
    mask: bool | None = None,
    verbose: bool = True,
) -> _Object:
    'Predict splice-altering effects using SpliceAI deep learning model. Returns delta scores for acceptor gain (DS_AG), acceptor loss (DS_AL), donor gain (DS_DG), donor loss (DS_DL). Scores ≥0.2 suggest splice impact, ≥0.5 moderate, ≥0.8 high pathogenicity. Critical for intronic/exonic variant interpretation.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['variant'] = variant
    if genome is not None:
        _dnakit_arguments['genome'] = genome
    if distance is not None:
        _dnakit_arguments['distance'] = distance
    if mask is not None:
        _dnakit_arguments['mask'] = mask
    return _call('splice_ai_predict_splice', _dnakit_arguments, verbose=verbose)

def structure_annotate_per_residue(*,
    operation: str | None = None,
    pdb_id: str | None | None = None,
    pdb_content: str | None | None = None,
    target_chain: str | None = None,
    partner_chains: list[str] | None = None,
    ligand_resnames: list[str] | None = None,
    distance_cutoff: float | None = None,
    core_rsa_cutoff: float | None = None,
    include_secondary_structure: bool | None = None,
    verbose: bool = True,
) -> _Object:
    'Per-residue structural annotation from a PDB structure: which residues sit at a binding interface (vs a partner chain), which line a ligand pocket, which are buried (core) vs solvent-exposed (surface). For each residue of the target chain returns: 1-letter amino acid, min side-chain-heavy-atom distance to the partner chain(s), min distance to ligand heavy atoms, relative SASA (computed with freesasa on the ISOLATED target chain — partner chains would bias the call by burying interface residues), a region label in {interface, ligand, both, other}, and an is_core boolean. Optionally fetches secondary structure (helix/strand/coil) from PDBe REST. This is the structural annotation track plotted under a DMS heatmap and the structural prior SAE feature drops are read against. Methodology adapted from an upstream research workflow.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if pdb_id is not None:
        _dnakit_arguments['pdb_id'] = pdb_id
    if pdb_content is not None:
        _dnakit_arguments['pdb_content'] = pdb_content
    if target_chain is not None:
        _dnakit_arguments['target_chain'] = target_chain
    if partner_chains is not None:
        _dnakit_arguments['partner_chains'] = partner_chains
    if ligand_resnames is not None:
        _dnakit_arguments['ligand_resnames'] = ligand_resnames
    if distance_cutoff is not None:
        _dnakit_arguments['distance_cutoff'] = distance_cutoff
    if core_rsa_cutoff is not None:
        _dnakit_arguments['core_rsa_cutoff'] = core_rsa_cutoff
    if include_secondary_structure is not None:
        _dnakit_arguments['include_secondary_structure'] = include_secondary_structure
    return _call('structure_annotate_per_residue', _dnakit_arguments, verbose=verbose)

def tpmi_phewas_by_variant(*,
    rsid: str | None = None,
    variant: str | None = None,
    limit: int | None = None,
    max_pval: float | None | None = None,
    verbose: bool = True,
) -> _Object:
    'Phenome-wide association study (PheWAS) lookup for a single variant in the Taiwan Precision Medicine Initiative (TPMI) PheWeb, a large Taiwanese (Han Chinese) biobank. Returns every phenotype associated with the variant, sorted by p-value. Accepts an rsID (auto-resolved to GRCh38 via Ensembl) or an explicit chr:pos:ref:alt variant in GRCh38. Use for cross-ancestry replication in an East Asian (Han Chinese) population.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if rsid is not None:
        _dnakit_arguments['rsid'] = rsid
    if variant is not None:
        _dnakit_arguments['variant'] = variant
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if max_pval is not None:
        _dnakit_arguments['max_pval'] = max_pval
    return _call('tpmi_phewas_by_variant', _dnakit_arguments, verbose=verbose)

def tark_get_mane_transcripts(*,
    gene: str | None | None = None,
    ensembl_id: str | None | None = None,
    refseq_id: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    'Look up the MANE (Matched Annotation from NCBI and EMBL-EBI) transcript(s) for a gene, mapping the Ensembl transcript (ENST) to its equivalent RefSeq transcript (NM_). Returns MANE Select and, where defined, MANE Plus Clinical pairings from the Ensembl Tark transcript archive. Query by HGNC gene symbol, an Ensembl transcript id, or a RefSeq transcript id (versions optional — they are matched on the base id). Use to find the clinically-recommended representative transcript for a gene and to convert between RefSeq and Ensembl transcript namespaces. No API key required.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    if ensembl_id is not None:
        _dnakit_arguments['ensembl_id'] = ensembl_id
    if refseq_id is not None:
        _dnakit_arguments['refseq_id'] = refseq_id
    return _call('tark_get_mane_transcripts', _dnakit_arguments, verbose=verbose)

def tark_get_transcript(*,
    stable_id: str,
    verbose: bool = True,
) -> _Object:
    'Get the archived transcript record for an Ensembl transcript id (ENST) from the Ensembl Tark transcript archive: assembly, biotype, genomic location (region/start/end/strand), the transcript sequence checksum, and the Ensembl release versions in which the transcript appears. Use to resolve transcript versions and verify transcript identity/equivalence across releases. No API key required.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['stable_id'] = stable_id
    return _call('tark_get_transcript', _dnakit_arguments, verbose=verbose)

def ucsc_get_encode_c_cres(*,
    genome: str | None = None,
    chrom: str,
    start: int,
    end: int,
    verbose: bool = True,
) -> _Object:
    'Get ENCODE4 candidate cis-Regulatory Elements (cCREs) for a genomic region from UCSC Genome Browser. cCREs are classified by epigenetic signatures into: Promoter-like (PLS, red), Proximal enhancer-like (pELS, orange), Distal enhancer-like (dELS, yellow), DNase-H3K4me3 (blue), CTCF-only/CTCF-bound (green). Each cCRE includes max Z-scores for DNase, H3K4me3, H3K27ac, and CTCF signals across biosamples. Use this to identify regulatory elements in a region of interest for functional genomics or GWAS variant interpretation.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if genome is not None:
        _dnakit_arguments['genome'] = genome
    _dnakit_arguments['chrom'] = chrom
    _dnakit_arguments['start'] = start
    _dnakit_arguments['end'] = end
    return _call('ucsc_get_encode_c_cres', _dnakit_arguments, verbose=verbose)

def ukbtopmed_phewas_by_variant(*,
    rsid: str | None = None,
    variant: str | None = None,
    limit: int | None = None,
    max_pval: float | None | None = None,
    verbose: bool = True,
) -> _Object:
    'Phenome-wide association study (PheWAS) lookup for a single variant in the UKB-TOPMed PheWeb (UK Biobank, ~400K European-ancestry samples, 1,400+ phecode phenotypes). Returns every phenotype associated with the variant, sorted by p-value. Accepts an rsID (auto-resolved to GRCh38 via Ensembl) or an explicit chr:pos:ref:alt variant in GRCh38. Use for variant pleiotropy assessment and replication in a European population.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if rsid is not None:
        _dnakit_arguments['rsid'] = rsid
    if variant is not None:
        _dnakit_arguments['variant'] = variant
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if max_pval is not None:
        _dnakit_arguments['max_pval'] = max_pval
    return _call('ukbtopmed_phewas_by_variant', _dnakit_arguments, verbose=verbose)

def usda_plants_get_characteristics(*,
    symbol: str,
    verbose: bool = True,
) -> _Object:
    "Get the USDA PLANTS morphology/physiology, growth-requirement, and reproduction characteristics for a plant by its PLANTS symbol (e.g. 'ABBA'). Returns a list of characteristic name/value/category entries (e.g. active growth period, shade tolerance, drought tolerance, mature height, bloom period). Resolves the symbol to the plant id automatically. No API key required.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['symbol'] = symbol
    return _call('usda_plants_get_characteristics', _dnakit_arguments, verbose=verbose)

def usda_plants_get_profile(*,
    symbol: str,
    verbose: bool = True,
) -> _Object:
    "Get a USDA PLANTS Database profile for a plant by its PLANTS symbol (e.g. 'ABBA' = balsam fir): scientific name, common name, taxonomic group/rank, duration (annual/perennial), growth habits (tree/shrub/forb), and synonyms. 'native_statuses' gives region-level (e.g. 'L48', 'CAN') native/introduced status; fine-grained per-state/county FIPS distribution codes are not returned by this endpoint (fips_distribution is included for forward-compatibility but is null in the current USDA API response). Authoritative US plant taxonomy. No API key required.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['symbol'] = symbol
    return _call('usda_plants_get_profile', _dnakit_arguments, verbose=verbose)

def uni_bind_get_dataset(*,
    operation: str | None = None,
    dataset_id: str,
    verbose: bool = True,
) -> _Object:
    "Retrieve full direct transcription factor (TF)-DNA binding-site detail for one UniBind dataset by its dataset_id (e.g. 'EXP030726.neural_stem_cells.SMAD3', obtained from UniBind_search_datasets). Returns the TF name, cell line/tissue, biological condition, source identifiers (GEO/ENCODE/PubMed), the JASPAR motif IDs used, and a list of TFBS model blocks. Each TFBS model gives the JASPAR motif (jaspar_id + version), total number of predicted binding sites (total_tfbs), the PWM score threshold, the centrality distance threshold, the CentriMo adjusted enrichment p-value, and direct download URLs for the binding sites in BED format, the underlying sequences in FASTA, and a summary plot. These are motif-anchored, base-pair-resolution, experimentally supported direct binding sites - distinct from raw ChIP-seq peaks. Public, no API key.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['dataset_id'] = dataset_id
    return _call('uni_bind_get_dataset', _dnakit_arguments, verbose=verbose)

def uni_bind_list_tfs(*,
    operation: str | None = None,
    search: str | None | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    "List the catalog of transcription factors (TFs) profiled in the UniBind database of direct TF-DNA binding sites. Optionally filter by a case-insensitive substring of the TF gene symbol via 'search' (e.g. 'SMAD' returns SMAD2, SMAD3, SMAD4...). Returns a sorted, de-duplicated list of TF gene symbols. Use this to discover valid tf_name values for UniBind_search_datasets. Public, no API key.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if search is not None:
        _dnakit_arguments['search'] = search
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('uni_bind_list_tfs', _dnakit_arguments, verbose=verbose)

def uni_bind_search_datasets(*,
    operation: str | None = None,
    tf_name: str | None | None = None,
    species: str | None | None = None,
    cell_line: str | None | None = None,
    collection: str | None | None = None,
    order: str | None | None = None,
    page: int | None = None,
    page_size: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Search the UniBind catalog of curated, experimentally derived direct transcription factor (TF)-DNA binding-site datasets (one dataset = one ChIP-seq experiment for one TF in one cell type/condition). Filter by TF gene symbol (tf_name), species scientific name (species, e.g. 'Homo sapiens'), cell type/tissue (cell_line), and collection ('Robust' or 'Permissive'). All filters are optional and compose. Returns a paginated list of matching datasets with their TF name, total peak count, and dataset_id (use UniBind_get_dataset for full binding-site detail). UniBind provides motif-anchored, base-pair-resolution direct binding sites derived from ChIP-seq via the DAMO/ChIP-eat pipeline with JASPAR motifs - distinct from JASPAR/HOCOMOCO motif matrices and from ReMap raw ChIP-seq peaks. Public, no API key.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if tf_name is not None:
        _dnakit_arguments['tf_name'] = tf_name
    if species is not None:
        _dnakit_arguments['species'] = species
    if cell_line is not None:
        _dnakit_arguments['cell_line'] = cell_line
    if collection is not None:
        _dnakit_arguments['collection'] = collection
    if order is not None:
        _dnakit_arguments['order'] = order
    if page is not None:
        _dnakit_arguments['page'] = page
    if page_size is not None:
        _dnakit_arguments['page_size'] = page_size
    return _call('uni_bind_search_datasets', _dnakit_arguments, verbose=verbose)

def uni_prot_get_disease_variants_by_accession(*,
    accession: str,
    verbose: bool = True,
) -> _Object:
    'Extract all variants (feature type = VARIANT) and their related annotations from UniProtKB entry.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['accession'] = accession
    return _call('uni_prot_get_disease_variants_by_accession', _dnakit_arguments, verbose=verbose)

def uni_prot_get_entry_by_accession(*,
    accession: str,
    compact: bool | None = None,
    verbose: bool = True,
) -> _Object:
    'Get the complete JSON entry for a specified UniProtKB accession. WARNING: This tool returns the complete UniProtKB entry, which can be extremely large (40,000+ lines for some entries) and may cause LLM context limits to be exceeded or crashes. For most use cases, consider using more specific extraction tools instead: UniProt_get_function_by_accession, UniProt_get_sequence_by_accession, UniProt_get_recommended_name_by_accession, UniProt_get_organism_by_accession, UniProt_get_subcellular_location_by_accession, UniProt_get_disease_variants_by_accession, or UniProt_get_ptm_processing_by_accession.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['accession'] = accession
    if compact is not None:
        _dnakit_arguments['compact'] = compact
    return _call('uni_prot_get_entry_by_accession', _dnakit_arguments, verbose=verbose)

def uni_prot_get_features_by_accession(*,
    accession: str,
    verbose: bool = True,
) -> _Object:
    "Extract all sequence features (domains, active/binding sites, motifs, modified residues, signal peptides, transmembrane regions, etc.) from a UniProtKB entry. Returns each feature's type, start/end positions, and description. Use to map functional/structural regions of a protein or to see what is annotated at a residue of interest. Example: P04637 (TP53) returns binding sites, the DNA-binding domain, modified residues, and more.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['accession'] = accession
    return _call('uni_prot_get_features_by_accession', _dnakit_arguments, verbose=verbose)

def uni_prot_get_function_by_accession(*,
    accession: str,
    verbose: bool = True,
) -> _Object:
    'Extract functional annotations from UniProtKB entry (Comment type = FUNCTION).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['accession'] = accession
    return _call('uni_prot_get_function_by_accession', _dnakit_arguments, verbose=verbose)

def uni_prot_get_sequence_by_accession(*,
    accession: str,
    verbose: bool = True,
) -> _Object:
    'Extract the canonical sequence from UniProtKB entry.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['accession'] = accession
    return _call('uni_prot_get_sequence_by_accession', _dnakit_arguments, verbose=verbose)

def uni_prot_search(*,
    query: str,
    organism: str | None = None,
    limit: int | None = None,
    min_length: int | None = None,
    max_length: int | None = None,
    fields: list[str] | None = None,
    verbose: bool = True,
) -> _Object:
    "Search UniProtKB database with flexible query syntax. Returns protein entries with accession numbers and metadata. Query syntax supports: field searches (gene:TP53, taxonomy_id:9606, reviewed:true), ranges (length:[100 TO 500], mass:[20000 TO 50000]), wildcards (gene:MEIOB*), boolean operators (AND/OR/NOT), and parentheses for grouping. Examples: 'gene:TP53 AND taxonomy_id:9606', 'length:[400 TO 500] AND reviewed:true', 'tissue:brain NOT taxonomy_id:10090'. Scope by taxon with taxonomy_id (lineage match, 'Taxonomy [OC]') rather than organism_id ('Organism [OS]', exact organism node only): most bacterial entries are annotated at strain level, so a species-level organism_id silently returns 0 or an undercount (gene:gyrA AND organism_id:1423 -> 0 hits; taxonomy_id:1423 -> 1). A caller-supplied organism_id:<taxid> is rewritten to taxonomy_id:<taxid> and the change is reported in data.normalization_note.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if min_length is not None:
        _dnakit_arguments['min_length'] = min_length
    if max_length is not None:
        _dnakit_arguments['max_length'] = max_length
    if fields is not None:
        _dnakit_arguments['fields'] = fields
    return _call('uni_prot_search', _dnakit_arguments, verbose=verbose)

def vcf_count_variants(*,
    operation: _Object | None = None,
    vcf_path: str,
    pass_only: bool | None | None = None,
    regions: str | None | None = None,
    min_qual: float | None | None = None,
    include_expr: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    'Count variants in a local VCF/BCF after applying filters, via bcftools. Returns total records plus SNP/indel/MNP breakdown. Supports PASS-only, region, min QUAL, and a custom bcftools include expression so SNP/indel counts under a filter are deterministic. Requires bcftools on PATH.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['vcf_path'] = vcf_path
    if pass_only is not None:
        _dnakit_arguments['pass_only'] = pass_only
    if regions is not None:
        _dnakit_arguments['regions'] = regions
    if min_qual is not None:
        _dnakit_arguments['min_qual'] = min_qual
    if include_expr is not None:
        _dnakit_arguments['include_expr'] = include_expr
    return _call('vcf_count_variants', _dnakit_arguments, verbose=verbose)

def vcf_normalize(*,
    operation: _Object | None = None,
    vcf_path: str,
    multiallelics: str | None | None = None,
    reference_fasta: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    'Normalize a local VCF/BCF with bcftools norm (split or join multiallelics; left-align indels against a reference FASTA when given) and report SNP/indel/record counts before vs after. Resolves the common cause of disagreeing indel counts: un-split multiallelics and non-left-aligned indels. Requires bcftools on PATH.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['vcf_path'] = vcf_path
    if multiallelics is not None:
        _dnakit_arguments['multiallelics'] = multiallelics
    if reference_fasta is not None:
        _dnakit_arguments['reference_fasta'] = reference_fasta
    return _call('vcf_normalize', _dnakit_arguments, verbose=verbose)

def vcf_summary_stats(*,
    operation: _Object | None = None,
    vcf_path: str,
    pass_only: bool | None | None = None,
    regions: str | None | None = None,
    min_qual: float | None | None = None,
    include_expr: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    'Canonical variant statistics for a local VCF/BCF file via bcftools stats. Returns records, SNPs, indels, MNPs, multiallelic sites, ts/tv ratio, and per-sample het/hom/missing counts. Optional PASS-only, region, QUAL, and bcftools-expression filters. Use this instead of writing an ad-hoc parser so counts are reproducible. Requires bcftools on PATH.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['vcf_path'] = vcf_path
    if pass_only is not None:
        _dnakit_arguments['pass_only'] = pass_only
    if regions is not None:
        _dnakit_arguments['regions'] = regions
    if min_qual is not None:
        _dnakit_arguments['min_qual'] = min_qual
    if include_expr is not None:
        _dnakit_arguments['include_expr'] = include_expr
    return _call('vcf_summary_stats', _dnakit_arguments, verbose=verbose)

def variant_validator_gene2transcripts(*,
    gene_symbol: str,
    gene: str | None = None,
    gene_name: str | None = None,
    transcript_set: str | None = None,
    genome_build: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get all transcripts for a gene using VariantValidator. Returns MANE Select and MANE Plus Clinical transcript annotations, Ensembl and RefSeq transcript identifiers, genomic coordinates, and cross-references. Useful for identifying the correct transcript for variant annotation. Supports filter for MANE transcripts only or all transcripts.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_symbol'] = gene_symbol
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    if gene_name is not None:
        _dnakit_arguments['gene_name'] = gene_name
    if transcript_set is not None:
        _dnakit_arguments['transcript_set'] = transcript_set
    if genome_build is not None:
        _dnakit_arguments['genome_build'] = genome_build
    return _call('variant_validator_gene2transcripts', _dnakit_arguments, verbose=verbose)

def variant_validator_validate_variant(*,
    genome_build: str,
    variant_description: str,
    select_transcripts: str,
    verbose: bool = True,
) -> _Object:
    'Validate and convert genetic variant descriptions using VariantValidator. Takes HGVS-formatted variant notation and returns validated variant descriptions with genomic context, transcript context, protein consequence, and cross-references. Validates variant nomenclature against current RefSeq sequences and LOVD syntax checker. Supports substitutions, deletions, insertions, duplications, inversions. Reference assemblies: GRCh37 (hg19) or GRCh38 (hg38).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['genome_build'] = genome_build
    _dnakit_arguments['variant_description'] = variant_description
    _dnakit_arguments['select_transcripts'] = select_transcripts
    return _call('variant_validator_validate_variant', _dnakit_arguments, verbose=verbose)

def wiki_pathways_search(*,
    query: str,
    organism: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search WikiPathways for pathways by free-text title match. Returns WPIDs, titles, and organisms. Backed by the sparql.wikipathways.org SPARQL endpoint.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('wiki_pathways_search', _dnakit_arguments, verbose=verbose)

def worm_base_get_expression(*,
    gene_id: str,
    verbose: bool = True,
) -> _Object:
    'Get expression data for a C. elegans gene from WormBase. Returns tissues where the gene is expressed, developmental stages, subcellular localization, and expression clusters. Expression data is curated from GFP reporters, antibody staining, RNA in situ, and microarray experiments. Example: WBGene00006763 (unc-26) is expressed in 38 tissues including neurons, pharynx, and body wall muscle.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_id'] = gene_id
    return _call('worm_base_get_expression', _dnakit_arguments, verbose=verbose)

def worm_base_get_gene(*,
    gene_id: str,
    verbose: bool = True,
) -> _Object:
    'Get detailed gene overview from WormBase for C. elegans (nematode worm) by WormBase gene ID. Returns gene name, sequence name, species, description, gene type, and status. C. elegans is one of the most important model organisms for genetics, neuroscience, and developmental biology. Example: WBGene00006763 returns unc-26 (synaptojanin, a polyphosphoinositide phosphatase required for synaptic vesicle recycling).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_id'] = gene_id
    return _call('worm_base_get_gene', _dnakit_arguments, verbose=verbose)

def worm_base_get_phenotypes(*,
    gene_id: str,
    verbose: bool = True,
) -> _Object:
    'Get phenotype annotations for a C. elegans gene from WormBase. Returns observed phenotypes (from RNAi, mutation, transgene experiments), phenotypes not observed, and total counts. Phenotype data uses the WormBase Phenotype Ontology (WPO). Example: WBGene00006763 (unc-26) has 125+ observed phenotypes including locomotion defects, synaptic vesicle recycling abnormalities, and pharyngeal pumping defects.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_id'] = gene_id
    return _call('worm_base_get_phenotypes', _dnakit_arguments, verbose=verbose)

def xenbase_get_gene(*,
    gene_id: str,
    verbose: bool = True,
) -> _Object:
    "Get detailed Xenopus gene information from Xenbase via the Alliance of Genome Resources. Returns gene symbol, name, synopsis, genomic location, synonyms, cross-references, and species details. Supports both Xenopus tropicalis and Xenopus laevis genes. Gene IDs use Xenbase format (e.g., 'XB-GENE-484286' or 'Xenbase:XB-GENE-484286').\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_id'] = gene_id
    return _call('xenbase_get_gene', _dnakit_arguments, verbose=verbose)

def xenbase_search_genes(*,
    query: str,
    species: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search for Xenopus (frog) genes in Xenbase via the Alliance of Genome Resources API. Returns gene IDs, symbols, names, and species for both Xenopus tropicalis and Xenopus laevis. Xenbase is the primary database for Xenopus genetics and developmental biology. Optionally filter by species (tropicalis or laevis). Useful for finding amphibian homologs of human disease genes.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if species is not None:
        _dnakit_arguments['species'] = species
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('xenbase_search_genes', _dnakit_arguments, verbose=verbose)

def zfin_get_gene(*,
    gene_id: str,
    verbose: bool = True,
) -> _Object:
    "Get detailed Danio rerio (zebrafish) gene information from ZFIN via the Alliance of Genome Resources. Returns gene symbol, name, synopsis, genomic location, synonyms, and cross-references. ZFIN is the primary database for zebrafish genetics. Example: 'ZFIN:ZDB-GENE-990415-8' returns pax2a (paired box 2a) with its role in nervous system development.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_id'] = gene_id
    return _call('zfin_get_gene', _dnakit_arguments, verbose=verbose)

def zfin_get_gene_expression(*,
    gene_id: str,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Get expression annotations for a zebrafish gene from the Alliance of Genome Resources. Returns per-annotation records of developmental stage (whenExpressed) and anatomical location (whereExpressed) for the gene, plus the total annotation count. Example: 'ZFIN:ZDB-GENE-980526-166' (shha) has hundreds of annotations across developmental stages. Zebrafish is a premier model for vertebrate development.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_id'] = gene_id
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('zfin_get_gene_expression', _dnakit_arguments, verbose=verbose)

def zfin_get_gene_phenotypes(*,
    gene_id: str,
    limit: int | None = None,
    page: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Get phenotype annotations for a zebrafish gene from ZFIN via the Alliance of Genome Resources. Returns phenotypes observed when the gene is mutated or knocked out in zebrafish, including affected anatomical structures and phenotype statements. Zebrafish is widely used for forward genetic screens and disease modeling. Example: 'ZFIN:ZDB-GENE-990415-8' (pax2a) returns phenotypes related to eye, kidney, and nervous system development.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_id'] = gene_id
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if page is not None:
        _dnakit_arguments['page'] = page
    return _call('zfin_get_gene_phenotypes', _dnakit_arguments, verbose=verbose)

def alphafold_get_prediction(*,
    qualifier: str | None = None,
    uniprot_id: str | None = None,
    uniprot_accession: str | None = None,
    accession: str | None = None,
    sequence_checksum: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Retrieve full AlphaFold 3D structure predictions for a given protein. IMPORTANT: The qualifier must be a UniProt ACCESSION (e.g., 'P69905' for HBA_HUMAN). Do NOT use UniProt entry names like 'HBA_HUMAN' - they will cause API errors. To find UniProt accession from a gene/protein name, use `UniProt_search` (e.g., query='gene:HBA' organism='human') or `UniProt_id_mapping` for ID conversion. Returns residue-level metadata including sequence, per-residue confidence scores (pLDDT), and structure download links (PDB, CIF, PAE). If you already have the accession and want UniProt details, call `UniProt_get_entry_by_accession`. For a quick overview, use `alphafold_get_summary`. For mutation/variant impact, see `alphafold_get_annotations`.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if qualifier is not None:
        _dnakit_arguments['qualifier'] = qualifier
    if uniprot_id is not None:
        _dnakit_arguments['uniprot_id'] = uniprot_id
    if uniprot_accession is not None:
        _dnakit_arguments['uniprot_accession'] = uniprot_accession
    if accession is not None:
        _dnakit_arguments['accession'] = accession
    if sequence_checksum is not None:
        _dnakit_arguments['sequence_checksum'] = sequence_checksum
    return _call('alphafold_get_prediction', _dnakit_arguments, verbose=verbose)

def c_bio_portal_get_cancer_studies(*,
    limit: int | None = None,
    offset: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get one page of the cBioPortal study catalogue. This returns a SLICE, not the whole catalogue: `limit` defaults to 20 while cBioPortal currently holds several hundred studies. The response always reports `total_available` (the full catalogue size) alongside `count` (studies actually returned), and sets `truncated: true` with a `truncation_note` whenever studies were left behind. Never conclude a study is missing from cBioPortal from a truncated page — raise `limit` to `total_available`, or page with `offset`, before drawing that conclusion.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    return _call('c_bio_portal_get_cancer_studies', _dnakit_arguments, verbose=verbose)

def c_bio_portal_get_gene_info(*,
    entrez_gene_id: int,
    verbose: bool = True,
) -> _Object:
    "Get basic information about a specific gene by Entrez Gene ID. Returns the HUGO gene symbol and gene type (e.g., 'protein-coding'); does not return aliases -- use HGNC_fetch_gene_by_symbol or NCBIGene_get_summary for gene aliases/chromosome location.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['entrez_gene_id'] = entrez_gene_id
    return _call('c_bio_portal_get_gene_info', _dnakit_arguments, verbose=verbose)

def c_bio_portal_get_molecular_profiles(*,
    study_id: str,
    verbose: bool = True,
) -> _Object:
    'Get molecular profiles for a cancer study. Molecular profiles include mutation data, copy number alterations, mRNA expression, etc.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['study_id'] = study_id
    return _call('c_bio_portal_get_molecular_profiles', _dnakit_arguments, verbose=verbose)

def c_bio_portal_get_mutations(*,
    study_id: str,
    gene_list: str,
    sample_list_id: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get mutation data for specific genes in a cancer study. This uses the updated cBioPortal API that requires molecular profiles.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['study_id'] = study_id
    _dnakit_arguments['gene_list'] = gene_list
    if sample_list_id is not None:
        _dnakit_arguments['sample_list_id'] = sample_list_id
    return _call('c_bio_portal_get_mutations', _dnakit_arguments, verbose=verbose)

def c_bio_portal_get_samples(*,
    study_id: str,
    page_size: int | None = None,
    page_number: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get all samples in a cancer study. Returns the full cohort by default; the response is a bare array with no truncation marker, so lower `page_size` only when you deliberately want a sample. Use this count as the denominator for mutation-frequency calculations.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['study_id'] = study_id
    if page_size is not None:
        _dnakit_arguments['page_size'] = page_size
    if page_number is not None:
        _dnakit_arguments['page_number'] = page_number
    return _call('c_bio_portal_get_samples', _dnakit_arguments, verbose=verbose)

def civic_get_molecular_profile(*,
    molecular_profile_id: int,
    evidence_limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Get detailed information about a specific molecular profile in CIViC database by molecular profile ID. Molecular profiles represent combinations of variants or features (e.g., 'BRAF V600E') that serve as biomarkers for clinical interpretation. Returns profile name, description, clinical significance score, and associated evidence items with diseases, therapies, and evidence levels.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['molecular_profile_id'] = molecular_profile_id
    if evidence_limit is not None:
        _dnakit_arguments['evidence_limit'] = evidence_limit
    return _call('civic_get_molecular_profile', _dnakit_arguments, verbose=verbose)

def civic_get_variant(*,
    variant_id: int,
    verbose: bool = True,
) -> _Object:
    'Get detailed information about a specific variant in CIViC database by variant ID. Variants represent specific genetic alterations (SNVs, indels, structural variants, etc.) with clinical significance in cancer. Returns variant name, coordinates, gene association, and related information.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['variant_id'] = variant_id
    return _call('civic_get_variant', _dnakit_arguments, verbose=verbose)

def civic_get_variants_by_gene(*,
    gene_id: int | None = None,
    gene_name: str | None = None,
    gene_symbol: str | None = None,
    gene: str | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Get all variants associated with a specific gene in CIViC database. Returns variant information including names, coordinates, and associated evidence. Use gene_name or gene_symbol (e.g., 'EGFR', 'BRAF') for convenience, or gene_id if already known. Uses cursor-based pagination to return up to 500 variants per call.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if gene_id is not None:
        _dnakit_arguments['gene_id'] = gene_id
    if gene_name is not None:
        _dnakit_arguments['gene_name'] = gene_name
    if gene_symbol is not None:
        _dnakit_arguments['gene_symbol'] = gene_symbol
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('civic_get_variants_by_gene', _dnakit_arguments, verbose=verbose)

def civic_search_assertions(*,
    limit: int | None = None,
    therapy: str | None = None,
    therapy_name: str | None = None,
    disease: str | None = None,
    disease_name: str | None = None,
    variant_name: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Search for assertions in CIViC database. Assertions are higher-level clinical interpretations that integrate multiple evidence items into formal statements about clinical actionability. Returns a list of assertions with descriptions and associated molecular profiles or variants. Filter by therapy (e.g., 'ponatinib'), disease (e.g., 'leukemia'), or variant name.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if therapy is not None:
        _dnakit_arguments['therapy'] = therapy
    if therapy_name is not None:
        _dnakit_arguments['therapy_name'] = therapy_name
    if disease is not None:
        _dnakit_arguments['disease'] = disease
    if disease_name is not None:
        _dnakit_arguments['disease_name'] = disease_name
    if variant_name is not None:
        _dnakit_arguments['variant_name'] = variant_name
    return _call('civic_search_assertions', _dnakit_arguments, verbose=verbose)

def civic_search_evidence_items(*,
    limit: int | None = None,
    status: str | None | None = None,
    therapy: str | None = None,
    therapy_name: str | None = None,
    disease: str | None = None,
    disease_name: str | None = None,
    molecular_profile: str | None = None,
    evidence_type: str | None | None = None,
    significance: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Search for evidence items in CIViC database. Evidence items are curated statements linking variants or molecular profiles to clinical outcomes, therapies, or disease states. Returns evidence with disease, therapies, significance, direction, and citation. Filter by therapy name (e.g., 'imatinib'), disease name (e.g., 'leukemia'), or molecular profile name (e.g., 'BRAF V600E').\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if status is not None:
        _dnakit_arguments['status'] = status
    if therapy is not None:
        _dnakit_arguments['therapy'] = therapy
    if therapy_name is not None:
        _dnakit_arguments['therapy_name'] = therapy_name
    if disease is not None:
        _dnakit_arguments['disease'] = disease
    if disease_name is not None:
        _dnakit_arguments['disease_name'] = disease_name
    if molecular_profile is not None:
        _dnakit_arguments['molecular_profile'] = molecular_profile
    if evidence_type is not None:
        _dnakit_arguments['evidence_type'] = evidence_type
    if significance is not None:
        _dnakit_arguments['significance'] = significance
    return _call('civic_search_evidence_items', _dnakit_arguments, verbose=verbose)

def civic_search_genes(*,
    limit: int | None = None,
    name: str | None = None,
    query: str | None = None,
    gene_name: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Search for genes in CIViC (Clinical Interpretation of Variants in Cancer) database. CIViC is a community knowledgebase for expert-curated interpretations of variants in cancer. Returns genes with their IDs, names, and descriptions. Use this to find genes of interest before querying variants or evidence.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if name is not None:
        _dnakit_arguments['name'] = name
    if query is not None:
        _dnakit_arguments['query'] = query
    if gene_name is not None:
        _dnakit_arguments['gene_name'] = gene_name
    return _call('civic_search_genes', _dnakit_arguments, verbose=verbose)

def civic_search_therapies(*,
    limit: int | None = None,
    name: str | None = None,
    query: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Search for therapies (drugs/treatments) in CIViC database. Returns a list of cancer therapies and drugs (with IDs and names) that are associated with variants and evidence in CIViC. Use 'query' or 'name' to filter by therapy name (e.g., 'imatinib', 'pembrolizumab').\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if name is not None:
        _dnakit_arguments['name'] = name
    if query is not None:
        _dnakit_arguments['query'] = query
    return _call('civic_search_therapies', _dnakit_arguments, verbose=verbose)

def coding_variant_fraction(*,
    file: str,
    vaf_threshold: float | None = None,
    annotation: str | None = None,
    header_rows: int | None = None,
    vaf_col: str | None | None = None,
    so_col: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Compute the fraction of CODING variants in an Excel/CSV variant table that match a given Sequence Ontology annotation and fall below a VAF threshold. Uses the CODING-only denominator convention (synonymous_variant, missense_variant, splice_region_variant, stop_gained, stop_lost, start_lost, frameshift_variant, inframe_insertion, inframe_deletion) - intronic, UTR, intergenic, and up/downstream records are excluded. Auto-detects 'variant allele frequency'/'VAF' and 'Sequence Ontology' columns if not specified. Supports flat headers (header_rows=1) and MultiIndex headers (header_rows=2) common in cancer variant reports. Returns the fraction with coding denominator plus the 'naive' fraction (over all VAF-filtered records) for comparison.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['file'] = file
    if vaf_threshold is not None:
        _dnakit_arguments['vaf_threshold'] = vaf_threshold
    if annotation is not None:
        _dnakit_arguments['annotation'] = annotation
    if header_rows is not None:
        _dnakit_arguments['header_rows'] = header_rows
    if vaf_col is not None:
        _dnakit_arguments['vaf_col'] = vaf_col
    if so_col is not None:
        _dnakit_arguments['so_col'] = so_col
    return _call('coding_variant_fraction', _dnakit_arguments, verbose=verbose)

def dbsnp_get_frequencies(*,
    rsid: str,
    verbose: bool = True,
) -> _Object:
    "Get allele frequencies for a variant from dbSNP. Returns population-specific allele frequency data. Each study reports an 'allele_count' -- the number of ALLELES (chromosomes) observed carrying the allele, NOT a number of individuals. Do not read it as study size: for a diploid cohort the number of people is at most half of it, and the denominator is the cohort's total allele count, which dbSNP does not report here.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['rsid'] = rsid
    return _call('dbsnp_get_frequencies', _dnakit_arguments, verbose=verbose)

def dbsnp_get_variant_by_rsid(*,
    rsid: str,
    verbose: bool = True,
) -> _Object:
    "Get variant information from dbSNP by rsID. Returns genomic coordinates, alleles, and basic variant information. Coordinates are reported on BOTH assemblies: 'position'/'position_grch38' are on GRCh38 (dbSNP's current assembly, named by the 'assembly' field) and 'position_grch37' is the same variant on GRCh37. Always pick the field matching your pipeline's build -- the two can differ by hundreds of kilobases.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['rsid'] = rsid
    return _call('dbsnp_get_variant_by_rsid', _dnakit_arguments, verbose=verbose)

def download_file(*,
    url: str,
    output_path: str | None = None,
    timeout: int | None = None,
    return_content: bool | None = None,
    chunk_size: int | None = None,
    follow_redirects: bool | None = None,
    verbose: bool = True,
) -> _Object:
    'Download files from HTTP/HTTPS URLs with cross-platform support (Windows, Mac, Linux). Similar to curl but platform-independent. Can save to specified path or temporary directory.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['url'] = url
    if output_path is not None:
        _dnakit_arguments['output_path'] = output_path
    if timeout is not None:
        _dnakit_arguments['timeout'] = timeout
    if return_content is not None:
        _dnakit_arguments['return_content'] = return_content
    if chunk_size is not None:
        _dnakit_arguments['chunk_size'] = chunk_size
    if follow_redirects is not None:
        _dnakit_arguments['follow_redirects'] = follow_redirects
    return _call('download_file', _dnakit_arguments, verbose=verbose)

def drugbank_get_drug_basic_info_by_drug_name_or_id(*,
    query: str,
    case_sensitive: bool | None = None,
    exact_match: bool | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get basic drug information including name, description, CAS number, and approval status by drug name or DrugBank ID.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if case_sensitive is not None:
        _dnakit_arguments['case_sensitive'] = case_sensitive
    if exact_match is not None:
        _dnakit_arguments['exact_match'] = exact_match
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('drugbank_get_drug_basic_info_by_drug_name_or_id', _dnakit_arguments, verbose=verbose)

def drugbank_get_pharmacology_by_drug_name_or_drugbank_id(*,
    query: str,
    case_sensitive: bool | None = None,
    exact_match: bool | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get drug pharmacodynamics, mechanism of action, and pharmacokinetics by drug name or Drugbank ID.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if case_sensitive is not None:
        _dnakit_arguments['case_sensitive'] = case_sensitive
    if exact_match is not None:
        _dnakit_arguments['exact_match'] = exact_match
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('drugbank_get_pharmacology_by_drug_name_or_drugbank_id', _dnakit_arguments, verbose=verbose)

def drugbank_get_targets_by_drug_name_or_drugbank_id(*,
    query: str,
    case_sensitive: bool | None = None,
    exact_match: bool | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get drug targets, enzymes, carriers, and transporters by drug name or DrugBank ID.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if case_sensitive is not None:
        _dnakit_arguments['case_sensitive'] = case_sensitive
    if exact_match is not None:
        _dnakit_arguments['exact_match'] = exact_match
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('drugbank_get_targets_by_drug_name_or_drugbank_id', _dnakit_arguments, verbose=verbose)

def drugbank_vocab_search(*,
    query: str,
    search_fields: list[str] | None = None,
    case_sensitive: bool | None = None,
    exact_match: bool | None = None,
    limit: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search the DrugBank vocabulary dataset for drugs by name, ID, synonyms, or other fields using text-based queries. Returns detailed drug information including DrugBank ID, common name, CAS number, UNII, and synonyms.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if search_fields is not None:
        _dnakit_arguments['search_fields'] = search_fields
    if case_sensitive is not None:
        _dnakit_arguments['case_sensitive'] = case_sensitive
    if exact_match is not None:
        _dnakit_arguments['exact_match'] = exact_match
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('drugbank_vocab_search', _dnakit_arguments, verbose=verbose)

def e_qtl_get_associations(*,
    dataset_id: str,
    size: int | None = None,
    gene_id: str | None | None = None,
    variant: str | None | None = None,
    molecular_trait_id: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    'Get eQTL associations (variant-gene expression correlations) from a specific dataset in the EBI eQTL Catalogue. Returns statistical associations between genetic variants and molecular traits (gene/exon/transcript expression or splicing). Each result includes beta effect size, p-value, allele frequency, and variant coordinates. Use eQTL_list_datasets first to find available dataset IDs. Can filter by gene_id (Ensembl gene ID) or variant (chr_pos_ref_alt format).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['dataset_id'] = dataset_id
    if size is not None:
        _dnakit_arguments['size'] = size
    if gene_id is not None:
        _dnakit_arguments['gene_id'] = gene_id
    if variant is not None:
        _dnakit_arguments['variant'] = variant
    if molecular_trait_id is not None:
        _dnakit_arguments['molecular_trait_id'] = molecular_trait_id
    return _call('e_qtl_get_associations', _dnakit_arguments, verbose=verbose)

def e_qtl_list_datasets(*,
    size: int | None = None,
    quant_method: str | None | None = None,
    tissue_id: str | None | None = None,
    verbose: bool = True,
) -> _Object:
    'List available eQTL (expression Quantitative Trait Loci) datasets from the EBI eQTL Catalogue. Each dataset represents a study-tissue-quantification combination with measured genetic associations to gene expression. Returns dataset IDs, study labels, tissue types, sample sizes, and quantification methods (ge=gene expression, exon, tx=transcript, txrev=transcript ratio, leafcutter=splicing). Use dataset IDs with eQTL_get_associations for specific variant-gene associations.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if size is not None:
        _dnakit_arguments['size'] = size
    if quant_method is not None:
        _dnakit_arguments['quant_method'] = quant_method
    if tissue_id is not None:
        _dnakit_arguments['tissue_id'] = tissue_id
    return _call('e_qtl_list_datasets', _dnakit_arguments, verbose=verbose)

def enrichr_gene_enrichment_analysis(*,
    gene_list: list[str],
    libs: list[str] | None = None,
    verbose: bool = True,
) -> _Object:
    'Perform gene enrichment analysis using Enrichr to find biological pathways, processes, and molecular functions associated with a list of genes. Returns connectivity paths between genes and enrichment terms.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_list'] = gene_list
    if libs is not None:
        _dnakit_arguments['libs'] = libs
    return _call('enrichr_gene_enrichment_analysis', _dnakit_arguments, verbose=verbose)

def ensembl_get_homology(*,
    species: str,
    symbol: str,
    target_species: str | None = None,
    target_taxon: str | None = None,
    type: str | None = None,
    sequence: str | None = None,
    aligned: bool | None = None,
    verbose: bool = True,
) -> _Object:
    'Get homology (orthologues and paralogues) for a gene by symbol across species. Returns evolutionary relationships with gene IDs, protein IDs, chromosomal locations, and homology types (1-to-1, 1-to-many, etc). Essential for comparative genomics studies.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['species'] = species
    _dnakit_arguments['symbol'] = symbol
    if target_species is not None:
        _dnakit_arguments['target_species'] = target_species
    if target_taxon is not None:
        _dnakit_arguments['target_taxon'] = target_taxon
    if type is not None:
        _dnakit_arguments['type'] = type
    if sequence is not None:
        _dnakit_arguments['sequence'] = sequence
    if aligned is not None:
        _dnakit_arguments['aligned'] = aligned
    return _call('ensembl_get_homology', _dnakit_arguments, verbose=verbose)

def ensembl_get_overlap_features(*,
    species: str,
    region: str,
    feature: str | None = None,
    biotype: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get genomic features (genes, transcripts, repeats) overlapping a region. Returns comprehensive annotations including coding genes, ncRNAs, pseudogenes, and structural features. Use this for genome annotation and feature discovery.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['species'] = species
    _dnakit_arguments['region'] = region
    if feature is not None:
        _dnakit_arguments['feature'] = feature
    if biotype is not None:
        _dnakit_arguments['biotype'] = biotype
    return _call('ensembl_get_overlap_features', _dnakit_arguments, verbose=verbose)

def ensembl_get_regulatory_features(*,
    region: str,
    species: str | None = None,
    feature: str,
    verbose: bool = True,
) -> _Object:
    'Get regulatory features (promoters, enhancers, transcription factor binding sites) in a genomic region. Returns regulatory annotations including feature types, positions, and evidence. Use this to identify regulatory elements that may control gene expression.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['region'] = region
    if species is not None:
        _dnakit_arguments['species'] = species
    _dnakit_arguments['feature'] = feature
    return _call('ensembl_get_regulatory_features', _dnakit_arguments, verbose=verbose)

def ensembl_get_sequence(*,
    id: str,
    type: str | None = None,
    species: str | None = None,
    multiple_sequences: bool | None = None,
    verbose: bool = True,
) -> _Object:
    'Get DNA or protein sequence for a gene, transcript, or genomic region. Returns sequence in JSON format with metadata. Supports multiple sequence types: genomic, CDS, cDNA, and peptide sequences.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['id'] = id
    if type is not None:
        _dnakit_arguments['type'] = type
    if species is not None:
        _dnakit_arguments['species'] = species
    if multiple_sequences is not None:
        _dnakit_arguments['multiple_sequences'] = multiple_sequences
    return _call('ensembl_get_sequence', _dnakit_arguments, verbose=verbose)

def ensembl_get_structural_variants(*,
    species: str,
    region: str | None = None,
    feature: str | None = None,
    chrom: str | None = None,
    start: int | None = None,
    end: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get structural variants overlapping a genomic region from Ensembl. Returns known SVs from DGVa (Database of Genomic Variants archive), ClinGen, and other sources. Each SV includes its ID (nsv/esv accession), coordinates, and feature type. Use ensembl_get_sv_detail to get clinical significance, variant class (CNV, deletion, duplication, inversion), and supporting evidence for specific variants. Region must be < 5Mb.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['species'] = species
    if region is not None:
        _dnakit_arguments['region'] = region
    if feature is not None:
        _dnakit_arguments['feature'] = feature
    if chrom is not None:
        _dnakit_arguments['chrom'] = chrom
    if start is not None:
        _dnakit_arguments['start'] = start
    if end is not None:
        _dnakit_arguments['end'] = end
    return _call('ensembl_get_structural_variants', _dnakit_arguments, verbose=verbose)

def ensembl_get_sv_detail(*,
    species: str,
    id: str | None = None,
    sv_id: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get detailed information for a specific structural variant from Ensembl by its accession (nsv/esv ID). Returns variant class (CNV, deletion, duplication, inversion, gain, loss), clinical significance (pathogenic, benign, etc.), source database, supporting evidence with individual sample observations, study information, and genomic mappings. Essential for CNV clinical interpretation.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['species'] = species
    if id is not None:
        _dnakit_arguments['id'] = id
    if sv_id is not None:
        _dnakit_arguments['sv_id'] = sv_id
    return _call('ensembl_get_sv_detail', _dnakit_arguments, verbose=verbose)

def ensembl_lookup_gene(*,
    gene_id: str | None = None,
    species: str | None = None,
    id: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Lookup gene information by Ensembl gene ID or symbol. Returns gene details including location, description, and external references. Supports both stable IDs (e.g., ENSG00000139618) and gene symbols (e.g., BRCA1). For gene symbols, species parameter is required.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if gene_id is not None:
        _dnakit_arguments['gene_id'] = gene_id
    if species is not None:
        _dnakit_arguments['species'] = species
    if id is not None:
        _dnakit_arguments['id'] = id
    return _call('ensembl_lookup_gene', _dnakit_arguments, verbose=verbose)

def ensembl_vep_region(*,
    species: str,
    region: str,
    allele: str | None = None,
    AncestralAllele: int | None = None,
    Blosum62: int | None = None,
    CADD: int | None = None,
    Conservation: int | None = None,
    DisGeNET: int | None = None,
    GeneSplicer: int | None = None,
    GO: int | None = None,
    LoF: int | None = None,
    MaxEntScan: int | None = None,
    Phenotypes: int | None = None,
    SIFT: str | None = None,
    PolyPhen: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Variant Effect Predictor (VEP) for genomic variants. Predicts functional consequences of variants including impact on genes, transcripts, and proteins. Returns SIFT/PolyPhen scores, consequence types, and affected transcripts. Essential for variant interpretation.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['species'] = species
    _dnakit_arguments['region'] = region
    if allele is not None:
        _dnakit_arguments['allele'] = allele
    if AncestralAllele is not None:
        _dnakit_arguments['AncestralAllele'] = AncestralAllele
    if Blosum62 is not None:
        _dnakit_arguments['Blosum62'] = Blosum62
    if CADD is not None:
        _dnakit_arguments['CADD'] = CADD
    if Conservation is not None:
        _dnakit_arguments['Conservation'] = Conservation
    if DisGeNET is not None:
        _dnakit_arguments['DisGeNET'] = DisGeNET
    if GeneSplicer is not None:
        _dnakit_arguments['GeneSplicer'] = GeneSplicer
    if GO is not None:
        _dnakit_arguments['GO'] = GO
    if LoF is not None:
        _dnakit_arguments['LoF'] = LoF
    if MaxEntScan is not None:
        _dnakit_arguments['MaxEntScan'] = MaxEntScan
    if Phenotypes is not None:
        _dnakit_arguments['Phenotypes'] = Phenotypes
    if SIFT is not None:
        _dnakit_arguments['SIFT'] = SIFT
    if PolyPhen is not None:
        _dnakit_arguments['PolyPhen'] = PolyPhen
    return _call('ensembl_vep_region', _dnakit_arguments, verbose=verbose)

def europepmc_disease_target_score(*,
    efoId: str,
    pageSize: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Extract disease-target association scores from Europe PMC literature. This includes literature-based evidence. Returns EVERY target with a Europe PMC score for the disease, sorted strongest-first -- it is not a paged or top-N view. Scanning is capped at 25 seconds; when the cap is hit the response sets 'truncated' and a 'note' saying how many of the disease's associated targets were scanned.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['efoId'] = efoId
    if pageSize is not None:
        _dnakit_arguments['pageSize'] = pageSize
    return _call('europepmc_disease_target_score', _dnakit_arguments, verbose=verbose)

def execute_tool(*,
    tool_name: str,
    arguments: _Object | None = None,
    verbose: bool = True,
) -> _Object:
    'Execute a ToolUniverse tool directly with custom arguments. This is the primary way to run any tool in the ToolUniverse system. You must get the tool definition first to figure out the arguments to pass.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['tool_name'] = tool_name
    if arguments is not None:
        _dnakit_arguments['arguments'] = arguments
    return _call('execute_tool', _dnakit_arguments, verbose=verbose)

def expression_atlas_disease_target_score(*,
    efoId: str,
    pageSize: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Extract disease-target association scores from Expression Atlas. This provides gene expression data.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['efoId'] = efoId
    if pageSize is not None:
        _dnakit_arguments['pageSize'] = pageSize
    return _call('expression_atlas_disease_target_score', _dnakit_arguments, verbose=verbose)

def get_joint_associated_diseases_by_hpo_id_list(*,
    HPO_ID_list: list[str],
    limit: int | None = None,
    offset: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Retrieve diseases associated with a list of phenotypes or symptoms by a list of HPO IDs.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['HPO_ID_list'] = HPO_ID_list
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    return _call('get_joint_associated_diseases_by_hpo_id_list', _dnakit_arguments, verbose=verbose)

def gnomad_get_gene_constraints(*,
    gene_symbol: str,
    reference_genome: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get gene constraint metrics from gnomAD for a gene symbol (GRCh38). Returns constraint fields under `exac_constraint` and `gnomad_constraint` (e.g., pLI, oe_lof/oe_mis/oe_syn, exp/obs counts). Mitochondrially-encoded genes (MT-*) are scored under `mitochondrial_constraint` instead -- their `gnomad_constraint` and `exac_constraint` are always null, which does not mean unconstrained. Note: does NOT include LOEUF, missense Z-score, or synonymous Z-score -- for the full constraint panel (pLI + LOEUF + mis_z + syn_z), use gnomad_get_constraint instead.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_symbol'] = gene_symbol
    if reference_genome is not None:
        _dnakit_arguments['reference_genome'] = reference_genome
    return _call('gnomad_get_gene_constraints', _dnakit_arguments, verbose=verbose)

def gnomad_get_sv_by_gene(*,
    gene_symbol: str,
    reference_genome: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get structural variants (SVs) from gnomAD for a gene. The SV callset follows the requested assembly: GRCh38 (default) uses the gnomAD v4 SV callset (gnomad_sv_r4), GRCh37 uses the gnomAD v2.1 SV callset (gnomad_sv_r2_1), so SV coordinates always match the gene coordinates returned. Returns deletions, duplications, inversions, translocations, and complex SVs with population allele frequencies, SV type, consequence, and homozygote counts (null for multi-allelic CNVs). Essential for interpreting CNVs in clinical genomics.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_symbol'] = gene_symbol
    if reference_genome is not None:
        _dnakit_arguments['reference_genome'] = reference_genome
    return _call('gnomad_get_sv_by_gene', _dnakit_arguments, verbose=verbose)

def gnomad_get_sv_by_region(*,
    chrom: str,
    start: int,
    stop: int | None = None,
    reference_genome: str | None = None,
    end: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get structural variants (SVs) from gnomAD overlapping a genomic region. The SV callset follows the requested assembly: GRCh38 (default) uses the gnomAD v4 SV callset (gnomad_sv_r4), GRCh37 uses the gnomAD v2.1 SV callset (gnomad_sv_r2_1), so SV coordinates are always in the same assembly as the requested region. Returns DEL/DUP/INV/BND/CPX variants with allele frequencies and consequences. Use for CNV interpretation in specific chromosomal regions. Region should be < 5Mb.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['chrom'] = chrom
    _dnakit_arguments['start'] = start
    if stop is not None:
        _dnakit_arguments['stop'] = stop
    if reference_genome is not None:
        _dnakit_arguments['reference_genome'] = reference_genome
    if end is not None:
        _dnakit_arguments['end'] = end
    return _call('gnomad_get_sv_by_region', _dnakit_arguments, verbose=verbose)

def gnomad_get_sv_detail(*,
    variant_id: str | None = None,
    dataset: str | None = None,
    sv_id: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Get detailed information for a specific gnomAD structural variant by its ID. Resolves IDs from both SV callsets, selected automatically from the ID format: GRCh38 IDs such as 'DEL_chr17_24e4872b' resolve against the gnomAD v4 SV callset (gnomad_sv_r4), while GRCh37 IDs such as 'DEL_1_1937' resolve against the gnomAD v2.1 SV callset (gnomad_sv_r2_1). Returns SV coordinates, type, length, population allele frequency, consequence, and filter status (consequence and homozygote/hemizygote counts are null for multi-allelic CNVs and intergenic SVs). Use gnomad_get_sv_by_gene or gnomad_get_sv_by_region to discover variant IDs first.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if variant_id is not None:
        _dnakit_arguments['variant_id'] = variant_id
    if dataset is not None:
        _dnakit_arguments['dataset'] = dataset
    if sv_id is not None:
        _dnakit_arguments['sv_id'] = sv_id
    return _call('gnomad_get_sv_detail', _dnakit_arguments, verbose=verbose)

def gnomad_get_variant(*,
    variant_id: str,
    dataset: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Get basic variant metadata from gnomAD by `variant_id` (format like '19-44908822-C-T'). Use `gnomad_search_variants` to discover valid `variant_id` values from an rsID query. `dataset` is required by gnomAD and defaults to `gnomad_r3`. `gnomad_r3` is a GENOMES-ONLY callset with no exome component, so under that default `exome` is null for every variant -- that means 'not present in this callset', NOT 'not observed in exomes'; responses carry a `dataset_note` saying so. Pass dataset='gnomad_r4' for exome allele counts (on GJB2 c.35delG that is a further 10,423 allele observations out of 1,461,600 alleles). Defaults deliberately differ across this gnomAD family: gnomad_get_variant, gnomad_search_variants and gnomad_get_region default to `gnomad_r3` while gnomad_get_variant_populations and gnomad_get_constraint default to `gnomad_r4`, so gnomad_get_variant and gnomad_get_variant_populations report different numbers for the same variant unless you pass the same `dataset` to both. Every response echoes `dataset` and `reference_genome` (gnomad_r2_1* and exac = GRCh37; gnomad_r3* and gnomad_r4* = GRCh38) so the callset that answered is never a guess.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['variant_id'] = variant_id
    if dataset is not None:
        _dnakit_arguments['dataset'] = dataset
    return _call('gnomad_get_variant', _dnakit_arguments, verbose=verbose)

def gnomad_get_variant_populations(*,
    variant_id: str,
    dataset: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Get per-ancestry (population-stratified) allele frequencies for a variant from gnomAD by `variant_id` (format like '1-55051215-G-GA'). Unlike `gnomad_get_variant` (aggregate AF only), this returns per-ancestry rows (`amr`, `nfe`, `fin`, `afr`, `eas`, `sas`, `asj`, `ami`, `mid`, `remaining`, plus `_XX`/`_XY` sex splits) with computed `af = ac / an` (af is null when an = 0), separated into `genome` and `exome` callsets. `dataset` defaults to `gnomad_r4`. Use `gnomad_search_variants` to find a `variant_id` from an rsID. `gnomad_r4` carries both genomes and exomes; `gnomad_r3` is genomes-only, so selecting it makes `exome` null for every variant (a property of the callset, not of the variant) and the response says so in `dataset_note`. Defaults deliberately differ across this gnomAD family: gnomad_get_variant, gnomad_search_variants and gnomad_get_region default to `gnomad_r3` while gnomad_get_variant_populations and gnomad_get_constraint default to `gnomad_r4`, so gnomad_get_variant and gnomad_get_variant_populations report different numbers for the same variant unless you pass the same `dataset` to both. Every response echoes `dataset` and `reference_genome` (gnomad_r2_1* and exac = GRCh37; gnomad_r3* and gnomad_r4* = GRCh38) so the callset that answered is never a guess.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['variant_id'] = variant_id
    if dataset is not None:
        _dnakit_arguments['dataset'] = dataset
    return _call('gnomad_get_variant_populations', _dnakit_arguments, verbose=verbose)

def gnomad_search_variants(*,
    query: str | None = None,
    dataset: str | None = None,
    variant: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Search for variants in gnomAD by free-text query (commonly an rsID like 'rs7412'). Returns matching `variant_id` strings. Note: gnomAD's `variant_search` response only exposes `variant_id`. Use the returned `variant_id` with `gnomad_get_variant`. `dataset` defaults to `gnomad_r3` (GRCh38, genomes only); which variants exist depends on the callset searched. Defaults deliberately differ across this gnomAD family: gnomad_get_variant, gnomad_search_variants and gnomad_get_region default to `gnomad_r3` while gnomad_get_variant_populations and gnomad_get_constraint default to `gnomad_r4`, so gnomad_get_variant and gnomad_get_variant_populations report different numbers for the same variant unless you pass the same `dataset` to both. Every response echoes `dataset` and `reference_genome` (gnomad_r2_1* and exac = GRCh37; gnomad_r3* and gnomad_r4* = GRCh38) so the callset that answered is never a guess.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if query is not None:
        _dnakit_arguments['query'] = query
    if dataset is not None:
        _dnakit_arguments['dataset'] = dataset
    if variant is not None:
        _dnakit_arguments['variant'] = variant
    return _call('gnomad_search_variants', _dnakit_arguments, verbose=verbose)

def gwas_get_associations_for_snp(*,
    rs_id: str,
    sort: str | None = None,
    direction: str | None = None,
    size: int | None = None,
    page: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get all associations for a specific SNP with optional sorting.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['rs_id'] = rs_id
    if sort is not None:
        _dnakit_arguments['sort'] = sort
    if direction is not None:
        _dnakit_arguments['direction'] = direction
    if size is not None:
        _dnakit_arguments['size'] = size
    if page is not None:
        _dnakit_arguments['page'] = page
    return _call('gwas_get_associations_for_snp', _dnakit_arguments, verbose=verbose)

def gwas_get_associations_for_study(*,
    accession_id: str,
    size: int | None = None,
    page: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get all associations for a specific study, sorted by p-value.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['accession_id'] = accession_id
    if size is not None:
        _dnakit_arguments['size'] = size
    if page is not None:
        _dnakit_arguments['page'] = page
    return _call('gwas_get_associations_for_study', _dnakit_arguments, verbose=verbose)

def gwas_get_associations_for_trait(*,
    disease_trait: str | None = None,
    efo_uri: str | None = None,
    efo_id: str | None = None,
    efo_trait: str | None = None,
    size: int | None = None,
    page: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Get all associations for a specific trait, sorted by p-value (most significant first). A free-text disease_trait is resolved to an EFO/MONDO term first, and that resolution is disclosed in the response: 'query_trait' (what you asked for), 'resolved_efo_id', 'resolved_efo_label' (the term's label in GWAS Catalog) and 'trait_resolution_source', plus a 'trait_resolution_note' whenever the resolved label is not an exact match for your query. Check those before trusting the results: e.g. 'cardiorespiratory fitness' resolves to EFO_0009184 'heart rate response to exercise', which is a different phenotype from VO2max (EFO_0004887 'maximal oxygen uptake measurement'). Pass efo_id explicitly to bypass the resolver.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if disease_trait is not None:
        _dnakit_arguments['disease_trait'] = disease_trait
    if efo_uri is not None:
        _dnakit_arguments['efo_uri'] = efo_uri
    if efo_id is not None:
        _dnakit_arguments['efo_id'] = efo_id
    if efo_trait is not None:
        _dnakit_arguments['efo_trait'] = efo_trait
    if size is not None:
        _dnakit_arguments['size'] = size
    if page is not None:
        _dnakit_arguments['page'] = page
    return _call('gwas_get_associations_for_trait', _dnakit_arguments, verbose=verbose)

def gwas_get_snp_by_id(*,
    rs_id: str,
    verbose: bool = True,
) -> _Object:
    'Get a specific GWAS SNP by its rs ID.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['rs_id'] = rs_id
    return _call('gwas_get_snp_by_id', _dnakit_arguments, verbose=verbose)

def gwas_get_snps_for_gene(*,
    size: int | None = None,
    page: int | None = None,
    gene_symbol: str | None = None,
    gene: str | None = None,
    mapped_gene: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Get all GWAS Catalog SNPs mapped to a specific gene by gene name (e.g., 'BRCA1', 'TP53'). Returns variant details including rs IDs, chromosomal positions, and functional annotations. Use this when you want to find all catalogued variants within or near a specific gene.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if size is not None:
        _dnakit_arguments['size'] = size
    if page is not None:
        _dnakit_arguments['page'] = page
    if gene_symbol is not None:
        _dnakit_arguments['gene_symbol'] = gene_symbol
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    if mapped_gene is not None:
        _dnakit_arguments['mapped_gene'] = mapped_gene
    return _call('gwas_get_snps_for_gene', _dnakit_arguments, verbose=verbose)

def gwas_get_studies_for_trait(*,
    disease_trait: str | None = None,
    efo_uri: str | None = None,
    efo_id: str | None = None,
    efo_trait: str | None = None,
    cohort: str | None = None,
    gxe: bool | None = None,
    full_pvalue_set: bool | None = None,
    size: int | None = None,
    page: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Get studies for a specific trait with optional filters for cohort, GxE interactions, and summary statistics.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if disease_trait is not None:
        _dnakit_arguments['disease_trait'] = disease_trait
    if efo_uri is not None:
        _dnakit_arguments['efo_uri'] = efo_uri
    if efo_id is not None:
        _dnakit_arguments['efo_id'] = efo_id
    if efo_trait is not None:
        _dnakit_arguments['efo_trait'] = efo_trait
    if cohort is not None:
        _dnakit_arguments['cohort'] = cohort
    if gxe is not None:
        _dnakit_arguments['gxe'] = gxe
    if full_pvalue_set is not None:
        _dnakit_arguments['full_pvalue_set'] = full_pvalue_set
    if size is not None:
        _dnakit_arguments['size'] = size
    if page is not None:
        _dnakit_arguments['page'] = page
    return _call('gwas_get_studies_for_trait', _dnakit_arguments, verbose=verbose)

def gwas_get_study_by_id(*,
    study_id: str | None = None,
    accession_id: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Get a specific GWAS study by its unique identifier. Accepts either 'study_id' or 'accession_id' (the same GCST identifier field name used by GWASAssociationsForStudy and returned on every association record) so IDs can be chained directly between gwas_* tools without renaming.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if study_id is not None:
        _dnakit_arguments['study_id'] = study_id
    if accession_id is not None:
        _dnakit_arguments['accession_id'] = accession_id
    return _call('gwas_get_study_by_id', _dnakit_arguments, verbose=verbose)

def gwas_get_variants_for_trait(*,
    disease_trait: str | None = None,
    trait: str | None = None,
    efo_uri: str | None = None,
    efo_id: str | None = None,
    efo_trait: str | None = None,
    size: int | None = None,
    limit: int | None = None,
    page: int | None = None,
    sort: str | None = None,
    direction: str | None = None,
    p_value_threshold: float | None = None,
    verbose: bool = True,
) -> _Object:
    "Search the GWAS Catalog for all genetic variants (SNPs) linked to a specific disease or trait. Accepts free-text disease/trait names (e.g., 'diabetes', 'breast cancer'), EFO ontology IDs, or EFO trait labels. Returns variant details including rsIDs, p-values, mapped genes, and genomic locations. Results are sorted by p-value ascending (most significant first) by default, so the first page holds the strongest loci; override with the 'sort' and 'direction' parameters. This is the best tool for finding all GWAS variants for a disease when you do NOT already have a specific variant ID. A free-text disease_trait is resolved to an EFO/MONDO term first, and that resolution is disclosed: 'query_trait', 'resolved_efo_id', 'resolved_efo_label' and 'trait_resolution_source' are returned, plus a 'trait_resolution_note' whenever the resolved label is not an exact match for your query - check them before trusting the results, and pass efo_id explicitly to bypass the resolver.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if disease_trait is not None:
        _dnakit_arguments['disease_trait'] = disease_trait
    if trait is not None:
        _dnakit_arguments['trait'] = trait
    if efo_uri is not None:
        _dnakit_arguments['efo_uri'] = efo_uri
    if efo_id is not None:
        _dnakit_arguments['efo_id'] = efo_id
    if efo_trait is not None:
        _dnakit_arguments['efo_trait'] = efo_trait
    if size is not None:
        _dnakit_arguments['size'] = size
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if page is not None:
        _dnakit_arguments['page'] = page
    if sort is not None:
        _dnakit_arguments['sort'] = sort
    if direction is not None:
        _dnakit_arguments['direction'] = direction
    if p_value_threshold is not None:
        _dnakit_arguments['p_value_threshold'] = p_value_threshold
    return _call('gwas_get_variants_for_trait', _dnakit_arguments, verbose=verbose)

def gwas_search_associations(*,
    disease_trait: str | None = None,
    query: str | None = None,
    efo_uri: str | None = None,
    efo_id: str | None = None,
    efo_trait: str | None = None,
    rs_id: str | None = None,
    accession_id: str | None = None,
    p_value: float | None = None,
    p_value_threshold: float | None = None,
    sort: str | None = None,
    direction: str | None = None,
    size: int | None = None,
    page: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Search the GWAS Catalog for genome-wide association study results by keyword. Accepts free-text disease/trait names (e.g., 'melanoma', 'diabetes', 'breast cancer'), SNP rs IDs, EFO ontology IDs, or study accession IDs. This is the primary tool for discovering GWAS associations when you do NOT already have a specific association ID. Returns SNPs, p-values, odds ratios, mapped genes, and study metadata with pagination support. Response envelope keys beyond 'data': 'metadata.pagination' (upstream page block - its totalElements/totalPages always describe the UNFILTERED result set); when a p_value threshold is given, 'metadata.p_value_filter', 'metadata.p_value_filter_scope' and 'metadata.p_value_filter_note' state that the threshold was applied CLIENT-SIDE (GWAS Catalog REST v2 has no server-side p-value filter), 'metadata.pagination_totals_are_prefilter' is true, 'metadata.filtered_total' / 'filtered_total_is_exact' / 'filtered_total_at_least' / 'filtered_total_note' give the post-filter denominator, and 'metadata.sort_override' records that the page was fetched p_value-ascending so the threshold could match anything. When a disease_trait is resolved to an ontology term, 'query_trait', 'resolved_efo_id', 'resolved_efo_label' and 'trait_resolution_source' are returned, plus a 'trait_resolution_note' whenever the resolved label is not an exact match for the trait you asked for. An empty 'data' comes with a 'note' that names the real cause: a client-side filter that removed every row is reported as such, and only a genuinely empty unfiltered result set suggests trying a trait synonym or a gene-based tool.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if disease_trait is not None:
        _dnakit_arguments['disease_trait'] = disease_trait
    if query is not None:
        _dnakit_arguments['query'] = query
    if efo_uri is not None:
        _dnakit_arguments['efo_uri'] = efo_uri
    if efo_id is not None:
        _dnakit_arguments['efo_id'] = efo_id
    if efo_trait is not None:
        _dnakit_arguments['efo_trait'] = efo_trait
    if rs_id is not None:
        _dnakit_arguments['rs_id'] = rs_id
    if accession_id is not None:
        _dnakit_arguments['accession_id'] = accession_id
    if p_value is not None:
        _dnakit_arguments['p_value'] = p_value
    if p_value_threshold is not None:
        _dnakit_arguments['p_value_threshold'] = p_value_threshold
    if sort is not None:
        _dnakit_arguments['sort'] = sort
    if direction is not None:
        _dnakit_arguments['direction'] = direction
    if size is not None:
        _dnakit_arguments['size'] = size
    if page is not None:
        _dnakit_arguments['page'] = page
    return _call('gwas_search_associations', _dnakit_arguments, verbose=verbose)

def gwas_search_snps(*,
    rs_id: str | None = None,
    mapped_gene: str | None = None,
    size: int | None = None,
    page: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Search the GWAS Catalog for single nucleotide polymorphisms (SNPs) by rs ID (e.g., 'rs7903146') or mapped gene name (e.g., 'TCF7L2'). Returns SNP details including chromosomal location, mapped genes, and functional class. Use this to look up specific variants or find all GWAS-catalogued SNPs for a gene.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if rs_id is not None:
        _dnakit_arguments['rs_id'] = rs_id
    if mapped_gene is not None:
        _dnakit_arguments['mapped_gene'] = mapped_gene
    if size is not None:
        _dnakit_arguments['size'] = size
    if page is not None:
        _dnakit_arguments['page'] = page
    return _call('gwas_search_snps', _dnakit_arguments, verbose=verbose)

def gwas_search_studies(*,
    disease_trait: str | None = None,
    efo_uri: str | None = None,
    efo_id: str | None = None,
    efo_trait: str | None = None,
    cohort: str | None = None,
    gxe: bool | None = None,
    full_pvalue_set: bool | None = None,
    size: int | None = None,
    page: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search for GWAS studies by various criteria including disease trait, cohort, GxE interactions, and summary statistics availability.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if disease_trait is not None:
        _dnakit_arguments['disease_trait'] = disease_trait
    if efo_uri is not None:
        _dnakit_arguments['efo_uri'] = efo_uri
    if efo_id is not None:
        _dnakit_arguments['efo_id'] = efo_id
    if efo_trait is not None:
        _dnakit_arguments['efo_trait'] = efo_trait
    if cohort is not None:
        _dnakit_arguments['cohort'] = cohort
    if gxe is not None:
        _dnakit_arguments['gxe'] = gxe
    if full_pvalue_set is not None:
        _dnakit_arguments['full_pvalue_set'] = full_pvalue_set
    if size is not None:
        _dnakit_arguments['size'] = size
    if page is not None:
        _dnakit_arguments['page'] = page
    return _call('gwas_search_studies', _dnakit_arguments, verbose=verbose)

def humanbase_ppi_analysis(*,
    gene_list: list[str] | None = None,
    tissue: str | None = None,
    max_node: int | None = None,
    interaction: str | None = None,
    string_mode: bool | None = None,
    genes: list[str] | None = None,
    top_n: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Retrieve tissue-specific protein-protein interactions and biological processes from HumanBase. Returns a NetworkX graph of tissue specific protein-protein interactions and a list of associated biological processes involeed by the given genes from Gene Ontology.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if gene_list is not None:
        _dnakit_arguments['gene_list'] = gene_list
    if tissue is not None:
        _dnakit_arguments['tissue'] = tissue
    if max_node is not None:
        _dnakit_arguments['max_node'] = max_node
    if interaction is not None:
        _dnakit_arguments['interaction'] = interaction
    if string_mode is not None:
        _dnakit_arguments['string_mode'] = string_mode
    if genes is not None:
        _dnakit_arguments['genes'] = genes
    if top_n is not None:
        _dnakit_arguments['top_n'] = top_n
    return _call('humanbase_ppi_analysis', _dnakit_arguments, verbose=verbose)

def iedb_get_epitope_antigens(*,
    structure_id: int,
    limit: int | None = None,
    offset: int | None = None,
    order: str | None = None,
    select: _Object | None = None,
    filters: dict[str, _Object] | None = None,
    verbose: bool = True,
) -> _Object:
    'Given an epitope `structure_id`, get linked antigens (IEDB Query API). Returned rows contain `parent_source_antigen_iri` (an antigen identifier) that you can use in other filters/searches.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['structure_id'] = structure_id
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    if order is not None:
        _dnakit_arguments['order'] = order
    if select is not None:
        _dnakit_arguments['select'] = select
    if filters is not None:
        _dnakit_arguments['filters'] = filters
    return _call('iedb_get_epitope_antigens', _dnakit_arguments, verbose=verbose)

def iedb_get_epitope_mhc(*,
    structure_id: int,
    limit: int | None = None,
    offset: int | None = None,
    order: str | None = None,
    select: _Object | None = None,
    filters: dict[str, _Object] | None = None,
    verbose: bool = True,
) -> _Object:
    'Given an epitope `structure_id`, get linked MHC restriction association rows (IEDB Query API). Returned rows include `elution_id` for the association.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['structure_id'] = structure_id
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    if order is not None:
        _dnakit_arguments['order'] = order
    if select is not None:
        _dnakit_arguments['select'] = select
    if filters is not None:
        _dnakit_arguments['filters'] = filters
    return _call('iedb_get_epitope_mhc', _dnakit_arguments, verbose=verbose)

def iedb_search_bcell(*,
    limit: int | None = None,
    offset: int | None = None,
    order: str | None = None,
    select: _Object | None = None,
    filters: dict[str, _Object] | None = None,
    verbose: bool = True,
) -> _Object:
    'Search B-cell assay rows (IEDB Query API). Use this to discover `bcell_id` values and related epitope `structure_id` values.\n\nTip: use `select` to request only a few fields (smaller + faster responses).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    if order is not None:
        _dnakit_arguments['order'] = order
    if select is not None:
        _dnakit_arguments['select'] = select
    if filters is not None:
        _dnakit_arguments['filters'] = filters
    return _call('iedb_search_bcell', _dnakit_arguments, verbose=verbose)

def iedb_search_epitopes(*,
    sequence_contains: str | None = None,
    structure_type: str | None = None,
    limit: int | None = None,
    offset: int | None = None,
    order: str | None = None,
    select: _Object | None = None,
    filters: dict[str, _Object] | None = None,
    verbose: bool = True,
) -> _Object:
    'Search immune epitopes (IEDB Query API). Use this to discover epitope `structure_id` values, then use the epitope linking tools (e.g., `iedb_get_epitope_references`) to retrieve related references, MHC, antigens, etc.\n\nTips:\n- Use `select` to request only a few fields (smaller + faster responses).\n- For advanced filtering, pass `filters` with values like `eq.123` or `ilike.*KVF*`.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if sequence_contains is not None:
        _dnakit_arguments['sequence_contains'] = sequence_contains
    if structure_type is not None:
        _dnakit_arguments['structure_type'] = structure_type
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    if order is not None:
        _dnakit_arguments['order'] = order
    if select is not None:
        _dnakit_arguments['select'] = select
    if filters is not None:
        _dnakit_arguments['filters'] = filters
    return _call('iedb_search_epitopes', _dnakit_arguments, verbose=verbose)

def iedb_search_mhc(*,
    limit: int | None = None,
    offset: int | None = None,
    order: str | None = None,
    select: _Object | None = None,
    filters: dict[str, _Object] | None = None,
    verbose: bool = True,
) -> _Object:
    'Search MHC-related assay rows (IEDB Query API). Use this to discover epitope `structure_id` values involved in MHC restriction.\n\nTip: use `select` to request only a few fields (smaller + faster responses).\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    if order is not None:
        _dnakit_arguments['order'] = order
    if select is not None:
        _dnakit_arguments['select'] = select
    if filters is not None:
        _dnakit_arguments['filters'] = filters
    return _call('iedb_search_mhc', _dnakit_arguments, verbose=verbose)

def iedb_search_tcell_assays(*,
    sequence: str | None = None,
    sequence_contains: str | None = None,
    mhc_class: str | None = None,
    qualitative_measure: str | None = None,
    limit: int | None = None,
    offset: int | None = None,
    select: _Object | None = None,
    filters: dict[str, _Object] | None = None,
    verbose: bool = True,
) -> _Object:
    'Search T-cell assay data from the IEDB. Returns T-cell response assays with epitope sequences, MHC restriction, source organism, and assay results.\n\nUseful for finding T-cell epitopes for vaccine design, immunotherapy, and understanding cellular immunity. Filter by peptide sequence, MHC class, organism, or qualitative result.\n\nTip: use `select` to request only key fields for faster responses.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if sequence is not None:
        _dnakit_arguments['sequence'] = sequence
    if sequence_contains is not None:
        _dnakit_arguments['sequence_contains'] = sequence_contains
    if mhc_class is not None:
        _dnakit_arguments['mhc_class'] = mhc_class
    if qualitative_measure is not None:
        _dnakit_arguments['qualitative_measure'] = qualitative_measure
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    if select is not None:
        _dnakit_arguments['select'] = select
    if filters is not None:
        _dnakit_arguments['filters'] = filters
    return _call('iedb_search_tcell_assays', _dnakit_arguments, verbose=verbose)

def intact_get_complex_details(*,
    complex_ac: str,
    format: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Get detailed information about a specific protein complex by its complex AC (Complex Accession, e.g., 'CPX-915'). Returns comprehensive complex information including: complex AC, name, synonyms, systematic name, functions, properties, species, ligands, complex assemblies, participants (with interactor details), and associated diseases. To find complex ACs, use intact_get_interactions_by_complex to search for complexes first.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['complex_ac'] = complex_ac
    if format is not None:
        _dnakit_arguments['format'] = format
    return _call('intact_get_complex_details', _dnakit_arguments, verbose=verbose)

def intact_get_interaction_network(*,
    identifier: str | None = None,
    uniprot_id: str | None = None,
    protein_id: str | None = None,
    gene_symbol: str | None = None,
    gene_name: str | None = None,
    depth: int | None = None,
    limit: int | None = None,
    size: int | None = None,
    format: str | None = None,
    protein_name: str | None = None,
    protein: str | None = None,
    gene: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Get interaction network centered on a specific interactor. Uses EBI Search API (IntAct domain) for reliable access. Accepts IntAct identifiers (e.g., 'EBI-1004115'), UniProt IDs (e.g., 'P04637'), or gene names (e.g., 'BRCA1'). Returns interactions with interaction IDs in the 'interaction_ids' field. These IDs can be used with intact_get_interaction_details for detailed information.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if identifier is not None:
        _dnakit_arguments['identifier'] = identifier
    if uniprot_id is not None:
        _dnakit_arguments['uniprot_id'] = uniprot_id
    if protein_id is not None:
        _dnakit_arguments['protein_id'] = protein_id
    if gene_symbol is not None:
        _dnakit_arguments['gene_symbol'] = gene_symbol
    if gene_name is not None:
        _dnakit_arguments['gene_name'] = gene_name
    if depth is not None:
        _dnakit_arguments['depth'] = depth
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if size is not None:
        _dnakit_arguments['size'] = size
    if format is not None:
        _dnakit_arguments['format'] = format
    if protein_name is not None:
        _dnakit_arguments['protein_name'] = protein_name
    if protein is not None:
        _dnakit_arguments['protein'] = protein
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    return _call('intact_get_interaction_network', _dnakit_arguments, verbose=verbose)

def intact_get_interactions(*,
    identifier: str,
    format: str | None = None,
    size: int | None | None = None,
    verbose: bool = True,
) -> _Object:
    "Get all interactions involving a specific interactor (protein or gene) by identifier. Uses EBI Search API (IntAct domain) for reliable access. Returns list of interactions with interaction IDs (e.g., 'EBI-366083-EBI-366083') in the 'interaction_ids' field. These IDs can be used with intact_get_interaction_details or intact_get_interaction_network to get detailed information. Accepts UniProt IDs, gene names, or IntAct identifiers. Returns 25 interactions per call by default; compare 'count' with 'hitCount' in the response metadata and raise 'size' to retrieve the rest.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['identifier'] = identifier
    if format is not None:
        _dnakit_arguments['format'] = format
    if size is not None:
        _dnakit_arguments['size'] = size
    return _call('intact_get_interactions', _dnakit_arguments, verbose=verbose)

def intact_search_interactions(*,
    query: str,
    first: int | None = None,
    max: int | None = None,
    format: str | None = None,
    species: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Search for interactions in IntAct database using EBI Search API (IntAct domain). Supports searching by protein names, gene names, interaction types, or other criteria. Use '*' as query to get all interactions (may be slow). Returns interaction entries with interaction IDs (e.g., 'EBI-366083-EBI-366083') in the 'interaction_ids' field. These IDs can be used with intact_get_interaction_details or intact_get_interaction_network to get detailed information.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if first is not None:
        _dnakit_arguments['first'] = first
    if max is not None:
        _dnakit_arguments['max'] = max
    if format is not None:
        _dnakit_arguments['format'] = format
    if species is not None:
        _dnakit_arguments['species'] = species
    return _call('intact_search_interactions', _dnakit_arguments, verbose=verbose)

def jaspar_get_matrix(*,
    matrix_id: str,
    verbose: bool = True,
) -> _Object:
    'Get full details for a single JASPAR matrix by `matrix_id` (e.g., MA0002.2). Returns metadata plus PFM/PWM. For reproducible analyses, use `jaspar_get_matrix_versions` with the corresponding `base_id` (e.g., MA0002) to see all versions and choose/cite the exact version you used.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['matrix_id'] = matrix_id
    return _call('jaspar_get_matrix', _dnakit_arguments, verbose=verbose)

def jaspar_get_matrix_versions(*,
    base_id: str,
    page: int | None = None,
    page_size: int | None = None,
    verbose: bool = True,
) -> _Object:
    'List all versions available for a JASPAR matrix base ID (e.g., base_id MA0002). This is useful for reproducibility (choose a specific matrix version like MA0002.2) and for tracking updates across releases.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['base_id'] = base_id
    if page is not None:
        _dnakit_arguments['page'] = page
    if page_size is not None:
        _dnakit_arguments['page_size'] = page_size
    return _call('jaspar_get_matrix_versions', _dnakit_arguments, verbose=verbose)

def jaspar_search_matrices(*,
    search: str | None = None,
    name: str | None = None,
    collection: str | None = None,
    tax_group: str | None = None,
    species: str | None = None,
    page: int | None = None,
    page_size: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Search/filter JASPAR matrices using common query parameters. Typical filters include `search` (free text), `name`, `collection`, `tax_group`, `species`, and paging via `page`/`page_size`.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if search is not None:
        _dnakit_arguments['search'] = search
    if name is not None:
        _dnakit_arguments['name'] = name
    if collection is not None:
        _dnakit_arguments['collection'] = collection
    if tax_group is not None:
        _dnakit_arguments['tax_group'] = tax_group
    if species is not None:
        _dnakit_arguments['species'] = species
    if page is not None:
        _dnakit_arguments['page'] = page
    if page_size is not None:
        _dnakit_arguments['page_size'] = page_size
    return _call('jaspar_search_matrices', _dnakit_arguments, verbose=verbose)

def kegg_find_genes(*,
    keyword: str | None = None,
    organism: str | None = None,
    query: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Find genes in KEGG database by keyword. Can search across all organisms or within a specific organism.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if keyword is not None:
        _dnakit_arguments['keyword'] = keyword
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if query is not None:
        _dnakit_arguments['query'] = query
    return _call('kegg_find_genes', _dnakit_arguments, verbose=verbose)

def kegg_get_gene_info(*,
    gene_id: str,
    verbose: bool = True,
) -> _Object:
    'Get detailed gene information from KEGG by gene ID. Returns gene data including sequence, function, and pathway associations.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['gene_id'] = gene_id
    return _call('kegg_get_gene_info', _dnakit_arguments, verbose=verbose)

def kegg_get_pathway_info(*,
    pathway_id: str,
    verbose: bool = True,
) -> _Object:
    'Get detailed pathway information from KEGG by pathway ID. Returns pathway data including genes, compounds, and reactions.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['pathway_id'] = pathway_id
    return _call('kegg_get_pathway_info', _dnakit_arguments, verbose=verbose)

def kegg_search_pathway(*,
    keyword: str,
    verbose: bool = True,
) -> _Object:
    'Search KEGG pathways by keyword. Returns pathway IDs and descriptions matching the search term.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['keyword'] = keyword
    return _call('kegg_search_pathway', _dnakit_arguments, verbose=verbose)

def ols_get_term_ancestors(*,
    operation: str | None = None,
    term_iri: str | None = None,
    term_id: str | None = None,
    ontology: str | None = None,
    include_obsolete: bool | None = None,
    size: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Get ancestor terms of a specific term in an ontology. Use `term_id` with CURIE format (e.g., 'HP:0001928') — the ontology is auto-inferred from the prefix. Or supply `term_iri` + `ontology` explicitly.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if term_iri is not None:
        _dnakit_arguments['term_iri'] = term_iri
    if term_id is not None:
        _dnakit_arguments['term_id'] = term_id
    if ontology is not None:
        _dnakit_arguments['ontology'] = ontology
    if include_obsolete is not None:
        _dnakit_arguments['include_obsolete'] = include_obsolete
    if size is not None:
        _dnakit_arguments['size'] = size
    return _call('ols_get_term_ancestors', _dnakit_arguments, verbose=verbose)

def ols_get_term_children(*,
    operation: str | None = None,
    term_iri: str | None = None,
    term_id: str | None = None,
    ontology: str | None = None,
    include_obsolete: bool | None = None,
    size: int | None = None,
    verbose: bool = True,
) -> _Object:
    "Get child terms of a specific term in an ontology. Use `term_id` with CURIE format (e.g., 'HP:0001928') — the ontology is auto-inferred from the prefix. Or supply `term_iri` + `ontology` explicitly.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if term_iri is not None:
        _dnakit_arguments['term_iri'] = term_iri
    if term_id is not None:
        _dnakit_arguments['term_id'] = term_id
    if ontology is not None:
        _dnakit_arguments['ontology'] = ontology
    if include_obsolete is not None:
        _dnakit_arguments['include_obsolete'] = include_obsolete
    if size is not None:
        _dnakit_arguments['size'] = size
    return _call('ols_get_term_children', _dnakit_arguments, verbose=verbose)

def ols_get_term_info(*,
    operation: str | None = None,
    id: str | None = None,
    term_id: str | None = None,
    term_iri: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get detailed information about a specific term in OLS\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    if id is not None:
        _dnakit_arguments['id'] = id
    if term_id is not None:
        _dnakit_arguments['term_id'] = term_id
    if term_iri is not None:
        _dnakit_arguments['term_iri'] = term_iri
    return _call('ols_get_term_info', _dnakit_arguments, verbose=verbose)

def ols_search_terms(*,
    operation: str | None = None,
    query: str,
    rows: int | None = None,
    size: int | None = None,
    ontology: str | None = None,
    exact_match: bool | None = None,
    include_obsolete: bool | None = None,
    limit: int | None | None = None,
    verbose: bool = True,
) -> _Object:
    'Search for terms in OLS (Ontology Lookup Service)\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if operation is not None:
        _dnakit_arguments['operation'] = operation
    _dnakit_arguments['query'] = query
    if rows is not None:
        _dnakit_arguments['rows'] = rows
    if size is not None:
        _dnakit_arguments['size'] = size
    if ontology is not None:
        _dnakit_arguments['ontology'] = ontology
    if exact_match is not None:
        _dnakit_arguments['exact_match'] = exact_match
    if include_obsolete is not None:
        _dnakit_arguments['include_obsolete'] = include_obsolete
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    return _call('ols_search_terms', _dnakit_arguments, verbose=verbose)

def openalex_get_author(*,
    author_id: str,
    mailto: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get a single OpenAlex author by Author ID (A...). You can pass either the short ID (e.g., "A5001226970") or the full URL (e.g., "https://openalex.org/A5001226970").\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['author_id'] = author_id
    if mailto is not None:
        _dnakit_arguments['mailto'] = mailto
    return _call('openalex_get_author', _dnakit_arguments, verbose=verbose)

def openalex_literature_search(*,
    search_keywords: str | None = None,
    query: str | None = None,
    max_results: int | None = None,
    limit: int | None = None,
    num_results: int | None = None,
    year_from: int | None = None,
    start_year: int | None = None,
    year_to: int | None = None,
    end_year: int | None = None,
    open_access: bool | None = None,
    require_has_fulltext: bool | None = None,
    fulltext_terms: list[str] | None = None,
    verbose: bool = True,
) -> _Object:
    'Search for academic literature using OpenAlex. Supports optional full-text-index filtering (has_fulltext + fulltext.search) when you need body-only term discovery, but note: OpenAlex full-text indexing is incomplete, and a work can be open-access elsewhere (e.g., PMC) while still having has_fulltext=false in OpenAlex. Returns paper metadata and (when available) OA links/content URLs.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if search_keywords is not None:
        _dnakit_arguments['search_keywords'] = search_keywords
    if query is not None:
        _dnakit_arguments['query'] = query
    if max_results is not None:
        _dnakit_arguments['max_results'] = max_results
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if num_results is not None:
        _dnakit_arguments['num_results'] = num_results
    if year_from is not None:
        _dnakit_arguments['year_from'] = year_from
    if start_year is not None:
        _dnakit_arguments['start_year'] = start_year
    if year_to is not None:
        _dnakit_arguments['year_to'] = year_to
    if end_year is not None:
        _dnakit_arguments['end_year'] = end_year
    if open_access is not None:
        _dnakit_arguments['open_access'] = open_access
    if require_has_fulltext is not None:
        _dnakit_arguments['require_has_fulltext'] = require_has_fulltext
    if fulltext_terms is not None:
        _dnakit_arguments['fulltext_terms'] = fulltext_terms
    return _call('openalex_literature_search', _dnakit_arguments, verbose=verbose)

def openalex_search_works(*,
    search: str | None = None,
    query: str | None = None,
    filter: str | None = None,
    require_has_fulltext: bool | None = None,
    fulltext_terms: list[str] | None = None,
    per_page: int | None = None,
    limit: int | None = None,
    page: int | None = None,
    sort: str | None = None,
    mailto: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Search OpenAlex works (papers) via the /works endpoint. Supports both general search and optional full-text-index filtering: set require_has_fulltext=true and/or pass fulltext_terms to add OpenAlex filters (has_fulltext:true and fulltext.search:<term>). Important: has_fulltext/fulltext.search only match works where OpenAlex has indexed full text; this can miss papers whose full text is available elsewhere (e.g., PMC) but not indexed by OpenAlex. Use this to discover Work IDs (W...) and DOIs, then follow up with openalex_get_work/openalex_get_work_by_doi to retrieve OA links and (sometimes) content_urls.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    if search is not None:
        _dnakit_arguments['search'] = search
    if query is not None:
        _dnakit_arguments['query'] = query
    if filter is not None:
        _dnakit_arguments['filter'] = filter
    if require_has_fulltext is not None:
        _dnakit_arguments['require_has_fulltext'] = require_has_fulltext
    if fulltext_terms is not None:
        _dnakit_arguments['fulltext_terms'] = fulltext_terms
    if per_page is not None:
        _dnakit_arguments['per_page'] = per_page
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if page is not None:
        _dnakit_arguments['page'] = page
    if sort is not None:
        _dnakit_arguments['sort'] = sort
    if mailto is not None:
        _dnakit_arguments['mailto'] = mailto
    return _call('openalex_search_works', _dnakit_arguments, verbose=verbose)

def proteins_api_get_features(*,
    accession: str | list[_Object],
    format: str | None = None,
    verbose: bool = True,
) -> _Object:
    'Get detailed sequence features for a protein, including domains, regions, motifs, modified residues, and other annotated features. Supports batch operations: pass a list of accessions or comma-separated string to retrieve multiple proteins in parallel. Note: This tool extracts features from the main protein endpoint (the separate features endpoint is not available). Returns all sequence feature annotations for the protein.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['accession'] = accession
    if format is not None:
        _dnakit_arguments['format'] = format
    return _call('proteins_api_get_features', _dnakit_arguments, verbose=verbose)

def proteins_api_search(*,
    query: str | None = None,
    size: int | None = None,
    offset: int | None = None,
    format: str | None = None,
    gene_symbol: str | None = None,
    gene: str | None = None,
    reviewed: bool | None = None,
    organism: str | None = None,
    taxid: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Search proteins in Proteins API by gene name, protein name, or accession. Defaults to human proteins (taxid=9606) for gene/name searches; use organism param to override. Gene names (e.g., 'BRCA1') search by gene parameter, accessions (e.g., 'P05067') search by accession parameter, and longer queries (e.g., 'insulin') search by protein parameter. Returns list of matching proteins with full annotations.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if query is not None:
        _dnakit_arguments['query'] = query
    if size is not None:
        _dnakit_arguments['size'] = size
    if offset is not None:
        _dnakit_arguments['offset'] = offset
    if format is not None:
        _dnakit_arguments['format'] = format
    if gene_symbol is not None:
        _dnakit_arguments['gene_symbol'] = gene_symbol
    if gene is not None:
        _dnakit_arguments['gene'] = gene
    if reviewed is not None:
        _dnakit_arguments['reviewed'] = reviewed
    if organism is not None:
        _dnakit_arguments['organism'] = organism
    if taxid is not None:
        _dnakit_arguments['taxid'] = taxid
    return _call('proteins_api_search', _dnakit_arguments, verbose=verbose)

def reactome_disease_target_score(*,
    efoId: str,
    pageSize: int | None = None,
    verbose: bool = True,
) -> _Object:
    'Extract disease-target association scores from Reactome pathway data. This includes pathway-based evidence.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0.'
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['efoId'] = efoId
    if pageSize is not None:
        _dnakit_arguments['pageSize'] = pageSize
    return _call('reactome_disease_target_score', _dnakit_arguments, verbose=verbose)

def read_executed_notebook(*,
    data_folder: str | None | None = None,
    notebook_path: str | None | None = None,
    search: str | None = None,
    regex: bool | None = None,
    max_output_chars: int | None = None,
    include_source: bool | None = None,
    verbose: bool = True,
) -> _Object:
    "Read a pre-executed Jupyter notebook (.ipynb, often named '*_executed.ipynb') and return its cells with printed outputs as structured JSON. Use this BEFORE reimplementing an analysis when the data folder ships with an executed notebook — the notebook already ran the authoritative pipeline with the exact package versions, filters, and thresholds that match the reference answers. Reimplementing with pydeseq2/scanpy/etc. produces different numbers. Accepts either a folder (auto-discovers the notebook) or a direct path. Optionally filters cells by a comma-separated search term list or a regex so the agent sees only the cells relevant to its question.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if data_folder is not None:
        _dnakit_arguments['data_folder'] = data_folder
    if notebook_path is not None:
        _dnakit_arguments['notebook_path'] = notebook_path
    if search is not None:
        _dnakit_arguments['search'] = search
    if regex is not None:
        _dnakit_arguments['regex'] = regex
    if max_output_chars is not None:
        _dnakit_arguments['max_output_chars'] = max_output_chars
    if include_source is not None:
        _dnakit_arguments['include_source'] = include_source
    return _call('read_executed_notebook', _dnakit_arguments, verbose=verbose)

def search_clinical_trials(*,
    condition: str | None = None,
    intervention: str | None = None,
    query_term: str | None = None,
    pageSize: int | None = None,
    pageToken: str | None = None,
    overall_status: list[str] | None = None,
    max_results: int | None = None,
    limit: int | None = None,
    status: _Object | None = None,
    keyword: str | None = None,
    verbose: bool = True,
) -> _Object:
    "Search ClinicalTrials.gov for clinical trials by disease/condition, drug/intervention, or keywords. This is the PRIMARY tool for finding clinical trials — use it when you need to discover trials for a drug (e.g., 'olaparib', 'pembrolizumab'), a disease (e.g., 'lung cancer', 'diabetes'), or any study criteria. Accepts condition names, intervention/treatment names, and free-text queries. Returns NCT IDs and brief descriptions with pagination. Use the returned NCT IDs with 'get_clinical_trial_*' tools for detailed protocol info, or 'extract_clinical_trial_*' tools for efficacy/safety results. Limited to trials beyond phase 1. Query rewriting: multi-word `intervention` values are matched as an EXACT PHRASE against each study's registered intervention names/synonyms, so a phrase no registrant recorded verbatim can match 0 studies even when many trials mention it. The query actually executed is returned in `executed_query` (detail in `query_rewrites`), and a zero-hit rewritten query reports in `relaxed_match_check` how many studies your original wording would have matched. Wrap a value in parentheses (e.g. '(cell therapy)') to send it verbatim.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    if condition is not None:
        _dnakit_arguments['condition'] = condition
    if intervention is not None:
        _dnakit_arguments['intervention'] = intervention
    if query_term is not None:
        _dnakit_arguments['query_term'] = query_term
    if pageSize is not None:
        _dnakit_arguments['pageSize'] = pageSize
    if pageToken is not None:
        _dnakit_arguments['pageToken'] = pageToken
    if overall_status is not None:
        _dnakit_arguments['overall_status'] = overall_status
    if max_results is not None:
        _dnakit_arguments['max_results'] = max_results
    if limit is not None:
        _dnakit_arguments['limit'] = limit
    if status is not None:
        _dnakit_arguments['status'] = status
    if keyword is not None:
        _dnakit_arguments['keyword'] = keyword
    return _call('search_clinical_trials', _dnakit_arguments, verbose=verbose)

def web_search(*,
    query: str,
    max_results: int | None = None,
    search_type: str | None = None,
    backend: str | None = None,
    region: str | None = None,
    safesearch: str | None = None,
    verbose: bool = True,
) -> _Object:
    "General web search using DDGS (Dux Distributed Global Search) supporting multiple search engines including Google, Bing, Brave, Yahoo, and DuckDuckGo. No API keys required. Parallel Search MCP and SerpBase are opt-in; selecting either sends the query to that hosted third-party service, so do not submit sensitive or patient-identifying information. SerpBase additionally requires SERPBASE_API_KEY and returns real Google SERP results (billed per request), useful when DDGS's Google backend is unavailable or rate-limited.\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
    _dnakit_arguments: dict[str, _Object] = {}
    _dnakit_arguments['query'] = query
    if max_results is not None:
        _dnakit_arguments['max_results'] = max_results
    if search_type is not None:
        _dnakit_arguments['search_type'] = search_type
    if backend is not None:
        _dnakit_arguments['backend'] = backend
    if region is not None:
        _dnakit_arguments['region'] = region
    if safesearch is not None:
        _dnakit_arguments['safesearch'] = safesearch
    return _call('web_search', _dnakit_arguments, verbose=verbose)

__all__ = ['alpha_genome_predict_interval', 'alpha_genome_score_variant', 'alpha_missense_get_protein_scores', 'alpha_missense_get_residue_scores', 'alpha_missense_get_variant_score', 'alphafold_get_prediction', 'ar_xiv_search_papers', 'bio_bank_japan_phewas_by_variant', 'bio_grid_get_interactions', 'bio_rxiv_get_preprint', 'bio_rxiv_list_recent_preprints', 'blast_protein_search', 'bvbrc_search_epitopes', 'c_bio_portal_get_cancer_studies', 'c_bio_portal_get_gene_info', 'c_bio_portal_get_molecular_profiles', 'c_bio_portal_get_mutations', 'c_bio_portal_get_samples', 'cadd_get_position_scores', 'cadd_get_range_scores', 'cadd_get_variant_score', 'cellx_gene_download_h5ad', 'cellx_gene_get_cell_metadata', 'cellx_gene_get_embeddings', 'cellx_gene_get_expression_data', 'ch_embl_get_drug_mechanisms', 'ch_embl_get_molecule', 'ch_embl_get_target_activities', 'ch_embl_search_drugs', 'ch_embl_search_targets', 'ch_ipatlas_enrichment_analysis', 'ch_ipatlas_get_experiments', 'ch_ipatlas_get_peak_data', 'ch_ipatlas_search_datasets', 'civic_get_molecular_profile', 'civic_get_variant', 'civic_get_variants_by_gene', 'civic_search_assertions', 'civic_search_evidence_items', 'civic_search_genes', 'civic_search_therapies', 'clin_gen_dosage_by_gene', 'clin_gen_dosage_region_search', 'clin_gen_get_actionability_adult', 'clin_gen_get_actionability_pediatric', 'clin_gen_get_dosage_sensitivity', 'clin_gen_get_gene_validity', 'clin_gen_get_variant_classifications', 'clin_gen_search_actionability', 'clin_gen_search_dosage_sensitivity', 'clin_gen_search_gene_validity', 'clin_var_get_clinical_significance', 'clin_var_get_variant_details', 'clin_var_search_variants', 'clinical_trials_search_studies', 'coding_variant_fraction', 'cosmic_get_mutations_by_gene', 'cosmic_search_mutations', 'crossref_search_works', 'ctd_get_chemical_gene_interactions', 'dbsnp_get_frequencies', 'dbsnp_get_variant_by_rsid', 'dep_map_get_gene_dependencies', 'dep_map_search_cell_lines', 'dgidb_get_drug_gene_interactions', 'dgidb_get_gene_druggability', 'dis_ge_net_get_disease_genes', 'dis_ge_net_get_gda', 'dis_ge_net_get_vda', 'dis_ge_net_search_disease', 'dis_ge_net_search_gene', 'dna_gibson_design', 'dna_golden_gate_assemble', 'dna_golden_gate_design', 'dna_translate_reading_frames', 'dna_virtual_digest', 'download_file', 'drugbank_get_drug_basic_info_by_drug_name_or_id', 'drugbank_get_pharmacology_by_drug_name_or_drugbank_id', 'drugbank_get_targets_by_drug_name_or_drugbank_id', 'drugbank_vocab_search', 'dyna_mut2_predict_stability', 'e_qtl_get_associations', 'e_qtl_list_datasets', 'enaportal_search_samples', 'enaportal_search_studies', 'encode_get_biosample', 'encode_get_chromatin_state', 'encode_get_experiment', 'encode_list_files', 'encode_search_annotations', 'encode_search_biosamples', 'encode_search_chromatin_accessibility', 'encode_search_experiments', 'encode_search_histone_experiments', 'encode_search_rnaseq_experiments', 'enrichr_enrich', 'enrichr_gene_enrichment_analysis', 'ensembl_compara_get_gene_tree', 'ensembl_compara_get_orthologues', 'ensembl_compara_get_paralogues', 'ensembl_get_homology', 'ensembl_get_overlap_features', 'ensembl_get_regulatory_features', 'ensembl_get_sequence', 'ensembl_get_structural_variants', 'ensembl_get_sv_detail', 'ensembl_lookup_gene', 'ensembl_seq_get_id_sequence', 'ensembl_seq_get_region_sequence', 'ensembl_var_get_variant_consequences', 'ensembl_vep_annotate_hgvs', 'ensembl_vep_annotate_rsid', 'ensembl_vep_region', 'ensembl_vep_variant_recoder', 'esm2_score_missense_variant', 'esm_describe_sae_feature', 'esm_explain_variant_mechanism', 'esm_get_region_sae_features', 'esm_get_sae_features', 'esm_score_sequence', 'esm_score_variant_sae_batch', 'esm_score_variant_sae_disruption', 'europe_pmc_search_articles', 'europepmc_disease_target_score', 'eve_get_gene_info', 'eve_get_variant_score', 'evo2_score_variant', 'execute_tool', 'expression_atlas_disease_target_score', 'expression_atlas_search_differential', 'expression_atlas_search_experiments', 'favor_annotate_variant', 'fda_get_active_ingredient_info_by_drug_name', 'fda_get_adverse_reactions_by_drug_name', 'fda_get_boxed_warning_info_by_drug_name', 'fda_get_clinical_studies_info_by_drug_name', 'fda_get_indications_by_drug_name', 'fda_get_mechanism_of_action_by_drug_name', 'finn_gen_get_variant_finemapping', 'fly_base_get_gene', 'fly_base_get_gene_alleles', 'fly_base_get_gene_disease_models', 'fly_base_get_gene_expression', 'fly_base_get_gene_interactions', 'fly_base_get_gene_orthologs', 'fly_mine_search', 'four_dn_search_data', 'gbif_get_species', 'gbif_search_species', 'gdc_get_clinical_data', 'gdc_get_mutation_frequency', 'gdc_get_ssm_by_gene', 'gdc_get_survival', 'gdc_list_files', 'gdc_list_projects', 'gdc_search_cases', 'gen_cc_get_classifications', 'gen_cc_search_disease', 'gen_cc_search_gene', 'genebass_gene_burden_phewas', 'geo_search_atacseq_datasets', 'geo_search_chipseq_datasets', 'geo_search_rnaseq_datasets', 'get_joint_associated_diseases_by_hpo_id_list', 'gnomad_get_gene_constraints', 'gnomad_get_sv_by_gene', 'gnomad_get_sv_by_region', 'gnomad_get_sv_detail', 'gnomad_get_variant', 'gnomad_get_variant_populations', 'gnomad_search_variants', 'go_get_annotations_for_gene', 'go_search_terms', 'gtdb_get_species', 'gtdb_get_taxon_info', 'gtdb_search_genomes', 'gtdb_search_taxon', 'gtex_calculate_eqtl', 'gtex_get_expression_summary', 'gtex_get_median_gene_expression', 'gtex_get_multi_tissue_eqtls', 'gtex_get_single_tissue_eqtls', 'gtex_get_tissue_sites', 'gtex_query_eqtl', 'gwas_get_associations_for_snp', 'gwas_get_associations_for_study', 'gwas_get_associations_for_trait', 'gwas_get_snp_by_id', 'gwas_get_snps_for_gene', 'gwas_get_studies_for_trait', 'gwas_get_study_by_id', 'gwas_get_variants_for_trait', 'gwas_search_associations', 'gwas_search_associations_by_gene', 'gwas_search_snps', 'gwas_search_studies', 'harmonizome_get_gene', 'hgnc_fetch_gene_by_symbol', 'hmdb_get_diseases', 'hmdb_get_metabolite', 'hmdb_search', 'hpa_get_cancer_prognostics_by_gene', 'hpa_get_protein_interactions_by_gene', 'hpa_get_rna_expression_by_source', 'hpa_get_rna_expression_in_specific_tissues', 'hpa_get_subcellular_location', 'hpa_search_genes_by_query', 'hpo_get_diseases_by_phenotype', 'hpo_get_genes_by_phenotype', 'hpo_get_term', 'hpo_search_terms', 'humanbase_ppi_analysis', 'iedb_get_epitope_antigens', 'iedb_get_epitope_mhc', 'iedb_search_bcell', 'iedb_search_epitopes', 'iedb_search_mhc', 'iedb_search_tcell_assays', 'igsr_list_data_collections', 'igsr_search_populations', 'igsr_search_samples', 'imgt_get_gene_info', 'imgt_search_genes', 'imm_port_search_studies', 'intact_get_complex_details', 'intact_get_interaction_network', 'intact_get_interactions', 'intact_search_interactions', 'inter_pro_get_entries_for_protein', 'inter_pro_get_protein_domains', 'jaspar_get_matrix', 'jaspar_get_matrix_versions', 'jaspar_get_transcription_factors', 'jaspar_search_matrices', 'kegg_find_genes', 'kegg_get_gene_info', 'kegg_get_gene_pathways', 'kegg_get_pathway_genes', 'kegg_get_pathway_info', 'kegg_search_pathway', 'ldlink_get_proxies', 'marrvel_get_gene', 'marrvel_get_omim_phenotypes', 'mave_db_get_effect_matrix', 'med_rxiv_get_preprint', 'mgi_get_gene', 'mgi_get_phenotypes', 'mgi_search_genes', 'mgnify_get_genome', 'mgnify_get_go_terms', 'mgnify_get_interpro', 'mgnify_get_study_detail', 'mgnify_get_taxonomy', 'mgnify_list_analyses', 'mgnify_list_biomes', 'mgnify_search_genomes', 'mgnify_search_studies', 'monarch_get_gene_diseases', 'monarch_get_gene_phenotypes', 'monarch_search', 'monarch_search_gene', 'monarch_v3_get_associations', 'monarch_v3_get_entity', 'monarch_v3_get_histopheno', 'monarch_v3_get_mappings', 'monarch_v3_phenotype_similarity_search', 'monarch_v3_search', 'my_gene_batch_query', 'my_gene_get_gene_annotation', 'my_gene_query_genes', 'my_variant_get_pathogenicity_scores', 'my_variant_get_variant_annotation', 'my_variant_query_variants', 'ncbi_fetch_accessions', 'ncbi_get_sequence', 'ncbi_search_nucleotide', 'ncbi_sra_search_runs', 'ncbidatasets_get_gene', 'ncbidatasets_get_gene_by_symbol', 'ncbidatasets_get_genome_assembly', 'ncbidatasets_get_orthologs', 'ncbidatasets_get_sequence_reports', 'ncbidatasets_get_taxonomy', 'ncbidatasets_list_genomes_by_taxon', 'ncbidatasets_suggest_taxonomy', 'ncbigene_get_summary', 'ncbigene_search', 'ncbiprotein_get_summary', 'nvidia_nim_alphafold2', 'nvidia_nim_esmfold', 'ols_get_term_ancestors', 'ols_get_term_children', 'ols_get_term_info', 'ols_search_terms', 'omim_get_clinical_synopsis', 'omim_get_entry', 'omim_get_gene_map', 'omim_search', 'onco_kb_annotate_variant', 'open_targets_get_associated_drugs_by_disease_efo_id', 'open_targets_get_associated_drugs_by_target_ensembl_id', 'open_targets_get_associated_phenotypes_by_disease_efo_id', 'open_targets_get_associated_targets_by_disease_efo_id', 'open_targets_get_credible_set_detail', 'open_targets_get_disease_ancestors_parents_by_efo_id', 'open_targets_get_disease_descendants_children_by_efo_id', 'open_targets_get_disease_description_by_efo_id', 'open_targets_get_disease_id_description_by_name', 'open_targets_get_disease_synonyms_by_efo_id', 'open_targets_get_disease_therapeutic_areas_by_efo_id', 'open_targets_get_diseases_phenotypes_by_target_ensembl', 'open_targets_get_drug_adverse_events_by_chembl_id', 'open_targets_get_drug_chemb_id_by_generic_name', 'open_targets_get_drug_indications_by_chembl_id', 'open_targets_get_drug_mechanisms_of_action_by_chembl_id', 'open_targets_get_drug_warnings_by_chembl_id', 'open_targets_get_evidence_by_datasource', 'open_targets_get_gwas_study', 'open_targets_get_publications_by_target_ensembl_id', 'open_targets_get_similar_entities_by_disease_efo_id', 'open_targets_get_study_credible_sets', 'open_targets_get_target_classes_by_ensembl_id', 'open_targets_get_target_gene_ontology_by_ensembl_id', 'open_targets_get_target_genomic_location_by_ensembl_id', 'open_targets_get_target_homologues_by_ensembl_id', 'open_targets_get_target_id_description_by_name', 'open_targets_get_target_interactions_by_ensembl_id', 'open_targets_get_target_safety_profile_by_ensembl_id', 'open_targets_get_target_tractability_by_ensembl_id', 'open_targets_get_variant_credible_sets', 'open_targets_get_variant_info', 'open_targets_map_any_disease_id_to_all_other_ids', 'open_targets_multi_entity_search_by_query_string', 'open_targets_search_gwas_studies_by_disease', 'open_targets_target_disease_evidence', 'openalex_get_author', 'openalex_literature_search', 'openalex_search_works', 'orphanet_get_classification', 'orphanet_get_disease', 'orphanet_get_epidemiology', 'orphanet_get_gene_diseases', 'orphanet_get_genes', 'orphanet_get_icd_mapping', 'orphanet_get_natural_history', 'orphanet_get_phenotypes', 'orphanet_search_by_name', 'orphanet_search_diseases', 'osl_get_efo_id_by_disease_name', 'panther_enrichment', 'panther_ortholog', 'pdbe_get_uniprot_mappings', 'pdbe_sifts_get_best_structures', 'pfam_get_protein_annotations', 'plant_reactome_get_pathway', 'plant_reactome_list_species', 'plant_reactome_search_pathways', 'powo_search_plants', 'progenetix_cnv_search', 'progenetix_list_cohorts', 'progenetix_list_filtering_terms', 'progenetix_search_biosamples', 'proteins_api_get_features', 'proteins_api_search', 'pub_chem_get_cid_by_compound_name', 'pub_chem_get_compound_properties_by_cid', 'pub_med_get_article', 'pub_med_get_cited_by', 'pub_med_get_related', 'pub_med_search_articles', 'quick_go_annotations_by_gene', 'rcsbdata_get_entry', 're_map_get_transcription_factor_binding', 'reactome_analysis_pathway_enrichment', 'reactome_content_search', 'reactome_disease_target_score', 'reactome_get_participants', 'reactome_get_pathway', 'reactome_get_pathway_reactions', 'reactome_map_uniprot_to_pathways', 'read_executed_notebook', 'regulome_db_query_variant', 'rgd_get_annotations', 'rgd_get_gene', 'rgd_get_orthologs', 'rgd_get_qtls_in_region', 'rgd_get_strain', 'rgd_get_strain_annotations', 'rgd_search_genes', 'rgd_search_strains', 'screen_get_regulatory_elements', 'search_clinical_trials', 'semantic_scholar_search_papers', 'sequence_count_residues', 'sequence_dn_ds', 'sequence_gc_content', 'sequence_reverse_complement', 'sequence_stats', 'sgd_get_gene', 'sgd_get_go_annotations', 'sgd_get_interactions', 'sgd_get_phenotypes', 'sgd_search', 'splice_ai_get_max_delta', 'splice_ai_predict_pangolin', 'splice_ai_predict_splice', 'sra_get_experiment', 'sra_search_experiments', 'string_functional_enrichment', 'string_get_functional_annotations', 'string_get_interaction_partners', 'string_get_network', 'string_ppi_enrichment', 'structure_annotate_per_residue', 'tark_get_mane_transcripts', 'tark_get_transcript', 'tpmi_phewas_by_variant', 'ucsc_get_encode_c_cres', 'ukbtopmed_phewas_by_variant', 'uni_bind_get_dataset', 'uni_bind_list_tfs', 'uni_bind_search_datasets', 'uni_prot_get_disease_variants_by_accession', 'uni_prot_get_entry_by_accession', 'uni_prot_get_features_by_accession', 'uni_prot_get_function_by_accession', 'uni_prot_get_sequence_by_accession', 'uni_prot_search', 'usda_plants_get_characteristics', 'usda_plants_get_profile', 'variant_validator_gene2transcripts', 'variant_validator_validate_variant', 'vcf_count_variants', 'vcf_normalize', 'vcf_summary_stats', 'web_search', 'wiki_pathways_search', 'worm_base_get_expression', 'worm_base_get_gene', 'worm_base_get_phenotypes', 'xenbase_get_gene', 'xenbase_search_genes', 'zfin_get_gene', 'zfin_get_gene_expression', 'zfin_get_gene_phenotypes']
