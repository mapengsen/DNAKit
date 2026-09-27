# Organism, taxonomy and microbiome queries: parameters and notes

See the [tutorial](../../api/features/22_database_query_organisms.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

## Retrieve Drosophila gene details {#dnakit.databases.fly_base_get_gene-parameters}

`dnakit.databases.fly_base_get_gene(*, gene_id, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_id` | string | Required | FlyBase gene ID with &#x27;FB:&#x27; prefix. Examples: &#x27;FB:FBgn0003996&#x27; (white/w), &#x27;FB:FBgn0000490&#x27; (decapentaplegic/dpp), &#x27;FB:FBgn0004647&#x27; (Notch/N), &#x27;FB:FBgn0000546&#x27; (eyeless/ey), &#x27;FB:FBgn0003731&#x27; (hedgehog/hh). |
| `verbose` | bool | True | Progress flag. |


## Retrieve Drosophila gene alleles and phenotypes {#dnakit.databases.fly_base_get_gene_alleles-parameters}

`dnakit.databases.fly_base_get_gene_alleles(*, gene_id, limit=None, page=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_id` | string | Required | FlyBase gene ID with &#x27;FB:&#x27; prefix. Examples: &#x27;FB:FBgn0003996&#x27; (white/w), &#x27;FB:FBgn0000490&#x27; (dpp), &#x27;FB:FBgn0004647&#x27; (Notch/N). |
| `limit` | integer | None | Maximum number of alleles to return (1-100). Default: 20. |
| `page` | integer | None | Page number for pagination. Default: 1. |
| `verbose` | bool | True | Progress flag. |


## Retrieve Drosophila gene disease models {#dnakit.databases.fly_base_get_gene_disease_models-parameters}

`dnakit.databases.fly_base_get_gene_disease_models(*, gene_id, limit=None, page=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_id` | string | Required | FlyBase gene ID with &#x27;FB:&#x27; prefix. Examples: &#x27;FB:FBgn0000490&#x27; (dpp), &#x27;FB:FBgn0004647&#x27; (Notch/N), &#x27;FB:FBgn0003996&#x27; (white/w). |
| `limit` | integer | None | Maximum number of disease models to return (1-100). Default: 20. |
| `page` | integer | None | Page number for pagination. Default: 1. |
| `verbose` | bool | True | Progress flag. |


## Retrieve Drosophila developmental and anatomical expression {#dnakit.databases.fly_base_get_gene_expression-parameters}

`dnakit.databases.fly_base_get_gene_expression(*, gene_id, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_id` | string | Required | FlyBase gene ID, with or without the &#x27;FB:&#x27; prefix. Examples: &#x27;FB:FBgn0000490&#x27; (dpp), &#x27;FB:FBgn0003996&#x27; (white/w), &#x27;FB:FBgn0004647&#x27; (Notch/N). |
| `limit` | integer | None | Max annotations to return (default 50, max 100). {&quot;default&quot;: 50} |
| `verbose` | bool | True | Progress flag. |


## Retrieve Drosophila molecular and genetic interactions {#dnakit.databases.fly_base_get_gene_interactions-parameters}

`dnakit.databases.fly_base_get_gene_interactions(*, gene_id, interaction_type=None, limit=None, page=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_id` | string | Required | FlyBase gene ID with &#x27;FB:&#x27; prefix. Examples: &#x27;FB:FBgn0003996&#x27; (white/w), &#x27;FB:FBgn0000490&#x27; (dpp), &#x27;FB:FBgn0004647&#x27; (Notch/N). |
| `interaction_type` | string | None | Type of interaction: &#x27;molecular&#x27; (physical PPI) or &#x27;genetic&#x27; (phenotype-based). Default: &#x27;molecular&#x27;. |
| `limit` | integer | None | Maximum number of interactions to return (1-100). Default: 20. |
| `page` | integer | None | Page number for pagination. Default: 1. |
| `verbose` | bool | True | Progress flag. |


## Retrieve cross-species orthologs of a Drosophila gene {#dnakit.databases.fly_base_get_gene_orthologs-parameters}

`dnakit.databases.fly_base_get_gene_orthologs(*, gene_id, limit=None, page=None, stringency=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_id` | string | Required | FlyBase gene ID with &#x27;FB:&#x27; prefix. Examples: &#x27;FB:FBgn0004647&#x27; (Notch/N), &#x27;FB:FBgn0003996&#x27; (white/w), &#x27;FB:FBgn0000546&#x27; (eyeless/ey). |
| `limit` | integer | None | Maximum number of orthologs to return (1-100). Default: 20. |
| `page` | integer | None | Page number for pagination. Default: 1. |
| `stringency` | string | None | Stringency filter for orthology calls: &#x27;stringent&#x27; (high confidence) or &#x27;moderate&#x27;. Default: &#x27;stringent&#x27;. |
| `verbose` | bool | True | Progress flag. |


## Search integrated Drosophila genomics records {#dnakit.databases.fly_mine_search-parameters}

`dnakit.databases.fly_mine_search(*, q, size=None, format=None, facet_Category=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `q` | string | Required | Search query: gene symbol (e.g., &#x27;dpp&#x27;, &#x27;Notch&#x27;, &#x27;hedgehog&#x27;), protein name, domain name, or free text. |
| `size` | integer | None | Number of results to return (default: 10, max: 100). {&quot;default&quot;: 10} |
| `format` | string | None | Response format. Must be &#x27;json&#x27;. {&quot;default&quot;: &quot;json&quot;} |
| `facet_Category` | string | None | Filter by entity category. Options: &#x27;Gene&#x27;, &#x27;Protein&#x27;, &#x27;Publication&#x27;, &#x27;ProteinDomain&#x27;, &#x27;CDNAClone&#x27;, &#x27;GOTerm&#x27;, &#x27;InteractionExperiment&#x27;, &#x27;RNAiScreen&#x27;, &#x27;UniProtFeature&#x27;. |
| `verbose` | bool | True | Progress flag. |


## Retrieve GBIF species details by taxon key {#dnakit.databases.gbif_get_species-parameters}

`dnakit.databases.gbif_get_species(*, speciesKey, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `speciesKey` | integer | Required | GBIF species/taxon key. Examples: 2435099 (Puma concolor / puma), 2878688 (Quercus robur / English oak), 5219404 (Canis lupus / wolf), 2480528 (Homo sapiens), 2495662 (Aquila chrysaetos / golden eagle). Use GBIF_search_species to find keys. |
| `verbose` | bool | True | Progress flag. |


## Search GBIF species names and taxonomy {#dnakit.databases.gbif_search_species-parameters}

`dnakit.databases.gbif_search_species(*, query, limit=None, offset=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required | Query keyword, e.g., Homo |
| `limit` | integer | None | Maximum results for the current request, not a database-wide total. {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 300} |
| `offset` | integer | None | Pagination offset within the current query. {&quot;default&quot;: 0, &quot;minimum&quot;: 0} |
| `verbose` | bool | True | Progress flag. |


## Retrieve GTDB species clusters and their genomes {#dnakit.databases.gtdb_get_species-parameters}

`dnakit.databases.gtdb_get_species(*, operation=None, species, max_genomes=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_species) {&quot;enum&quot;: &#91;&quot;get_species&quot;&#93;, &quot;default&quot;: &quot;get_species&quot;} |
| `species` | string | Required | Species name WITHOUT the s__ prefix (e.g., &#x27;Escherichia coli&#x27;, &#x27;Hydrogenothermus marinus&#x27;, &#x27;Bacillus subtilis&#x27;) |
| `max_genomes` | integer | None | Maximum number of genomes to return (default: 20). Species clusters can have thousands of genomes. |
| `verbose` | bool | True | Progress flag. |


## Retrieve GTDB taxonomy and release history {#dnakit.databases.gtdb_get_taxon_info-parameters}

`dnakit.databases.gtdb_get_taxon_info(*, operation=None, taxon, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_taxon_info) {&quot;enum&quot;: &#91;&quot;get_taxon_info&quot;&#93;, &quot;default&quot;: &quot;get_taxon_info&quot;} |
| `taxon` | string | Required | GTDB taxon name with rank prefix (e.g., &#x27;f__Lachnospiraceae&#x27;, &#x27;d__Bacteria&#x27;, &#x27;g__Escherichia&#x27;, &#x27;s__Escherichia coli&#x27;). Use GTDB_search_taxon to find valid names. |
| `verbose` | bool | True | Progress flag. |


## Search prokaryotic genomes by organism name {#dnakit.databases.gtdb_search_genomes-parameters}

`dnakit.databases.gtdb_search_genomes(*, operation=None, query, page=None, items_per_page=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search_genomes) {&quot;enum&quot;: &#91;&quot;search_genomes&quot;&#93;, &quot;default&quot;: &quot;search_genomes&quot;} |
| `query` | string | Required | Search query for organism name (e.g., &#x27;Escherichia&#x27;, &#x27;Staphylococcus aureus&#x27;, &#x27;Lactobacillus&#x27;) |
| `page` | integer | None | Page number for pagination (default: 1) |
| `items_per_page` | integer | None | Number of results per page (default: 10, max: 50) |
| `verbose` | bool | True | Progress flag. |


## Search prokaryotic taxa by partial name {#dnakit.databases.gtdb_search_taxon-parameters}

`dnakit.databases.gtdb_search_taxon(*, operation=None, query, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search_taxon) {&quot;enum&quot;: &#91;&quot;search_taxon&quot;&#93;, &quot;default&quot;: &quot;search_taxon&quot;} |
| `query` | string | Required | Partial taxon name to search for (e.g., &#x27;Lachnospiraceae&#x27;, &#x27;Escherichia&#x27;, &#x27;Bacillota&#x27;). Case-insensitive. |
| `limit` | integer | None | Maximum number of results per rank (default: 20, max: 100) |
| `verbose` | bool | True | Progress flag. |


## Retrieve mouse gene details by MGI identifier {#dnakit.databases.mgi_get_gene-parameters}

`dnakit.databases.mgi_get_gene(*, gene_id, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_id` | string | Required | MGI gene identifier with &#x27;MGI:&#x27; prefix. Examples: &#x27;MGI:98834&#x27; (Trp53/p53), &#x27;MGI:104537&#x27; (Brca1), &#x27;MGI:96677&#x27; (Kit), &#x27;MGI:88039&#x27; (Apoe), &#x27;MGI:88276&#x27; (Braf). |
| `verbose` | bool | True | Progress flag. |


## Retrieve mouse gene mutant phenotype annotations {#dnakit.databases.mgi_get_phenotypes-parameters}

`dnakit.databases.mgi_get_phenotypes(*, gene_id, limit=None, page=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_id` | string | Required | MGI gene identifier with &#x27;MGI:&#x27; prefix. Examples: &#x27;MGI:98834&#x27; (Trp53), &#x27;MGI:104537&#x27; (Brca1), &#x27;MGI:96677&#x27; (Kit). |
| `limit` | integer | None | Number of phenotype annotations to return per page (default: 20, max: 100). {&quot;default&quot;: 20} |
| `page` | integer | None | Page number for pagination (default: 1). {&quot;default&quot;: 1} |
| `verbose` | bool | True | Progress flag. |


## Search mouse genes by name or keyword {#dnakit.databases.mgi_search_genes-parameters}

`dnakit.databases.mgi_search_genes(*, query, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required | Gene search query: mouse gene symbol (e.g., &#x27;Trp53&#x27;, &#x27;Brca1&#x27;, &#x27;Sox9&#x27;), gene name keyword (e.g., &#x27;insulin&#x27;, &#x27;tumor&#x27;), or MGI ID (e.g., &#x27;MGI:98834&#x27;). |
| `limit` | integer | None | Maximum number of results to return (default: 10, max: 50). {&quot;default&quot;: 10} |
| `verbose` | bool | True | Progress flag. |


## Retrieve microbial genome quality and taxonomy {#dnakit.databases.mgnify_get_genome-parameters}

`dnakit.databases.mgnify_get_genome(*, genome_id, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `genome_id` | string | Required | MGnify genome accession. Examples: &#x27;MGYG000000001&#x27;, &#x27;MGYG000535630&#x27;. |
| `verbose` | bool | True | Progress flag. |


## Retrieve GO annotations from a metagenomics analysis {#dnakit.databases.mgnify_get_go_terms-parameters}

`dnakit.databases.mgnify_get_go_terms(*, analysis_id, page_size=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `analysis_id` | string | Required | MGnify analysis accession (e.g., &#x27;MGYA00585482&#x27;). Find IDs via MGnify_list_analyses. |
| `page_size` | integer | None | Number of GO terms to return per page (default 25, max 100). {&quot;default&quot;: 25} |
| `verbose` | bool | True | Progress flag. |


## Retrieve metagenomic protein domain annotations {#dnakit.databases.mgnify_get_interpro-parameters}

`dnakit.databases.mgnify_get_interpro(*, analysis_id, page_size=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `analysis_id` | string | Required | MGnify analysis accession (e.g., &#x27;MGYA00585482&#x27;). Find IDs via MGnify_list_analyses. |
| `page_size` | integer | None | Number of InterPro entries to return per page (default 25, max 100). {&quot;default&quot;: 25} |
| `verbose` | bool | True | Progress flag. |


## Retrieve MGnify metagenomics study details {#dnakit.databases.mgnify_get_study_detail-parameters}

`dnakit.databases.mgnify_get_study_detail(*, study_accession, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `study_accession` | string | Required | MGnify study accession. Examples: &#x27;MGYS00002008&#x27;, &#x27;MGYS00005292&#x27;. |
| `verbose` | bool | True | Progress flag. |


## Retrieve metagenomic taxonomic composition and abundance {#dnakit.databases.mgnify_get_taxonomy-parameters}

`dnakit.databases.mgnify_get_taxonomy(*, analysis_id, rna_type=None, page_size=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `analysis_id` | string | Required | MGnify analysis accession (e.g., &#x27;MGYA00585482&#x27;). Find IDs via MGnify_list_analyses. |
| `rna_type` | string | None | rRNA marker type: &#x27;ssu&#x27; for 16S/18S (default) or &#x27;lsu&#x27; for 23S/28S. {&quot;default&quot;: &quot;ssu&quot;} |
| `page_size` | integer | None | Number of taxa to return per page (default 25, max 100). {&quot;default&quot;: 25} |
| `verbose` | bool | True | Progress flag. |


## List analyses belonging to a metagenomics study {#dnakit.databases.mgnify_list_analyses-parameters}

`dnakit.databases.mgnify_list_analyses(*, study_accession, size=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `study_accession` | string | Required | MGnify study accession, e.g., &#x27;MGYS00000001&#x27; |
| `size` | integer | None | Limit for the returned page; does not automatically retrieve later pages. {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | Progress flag. |


## Browse the metagenomics biome hierarchy {#dnakit.databases.mgnify_list_biomes-parameters}

`dnakit.databases.mgnify_list_biomes(*, depth=None, page=None, page_size=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `depth` | integer | None | Hierarchy depth to filter by (1=root, 2=second level, etc.). |
| `page` | integer | None | Page number (default 1). |
| `page_size` | integer | None | Results per page (default 25, max 100). |
| `verbose` | bool | True | Progress flag. |


## Search metagenome-assembled and isolate genomes {#dnakit.databases.mgnify_search_genomes-parameters}

`dnakit.databases.mgnify_search_genomes(*, taxonomy=None, page=None, page_size=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `taxonomy` | string | None | GTDB taxonomy lineage filter, matched anywhere in the lineage. Examples: &#x27;Bacteroidota&#x27; (phylum), &#x27;Bacteroides&#x27; (genus), &#x27;Bacillota&#x27;, &#x27;Pseudomonadota&#x27;. NCBI-style phylum names such as &#x27;Bacteroidetes&#x27; or &#x27;Firmicutes&#x27; return few or no hits -- use the GTDB spelling. |
| `page` | integer | None | Page number (default 1). |
| `page_size` | integer | None | Results per page (default 25, max 100). |
| `verbose` | bool | True | Progress flag. |


## Search metagenomics studies by biome and keyword {#dnakit.databases.mgnify_search_studies-parameters}

`dnakit.databases.mgnify_search_studies(*, biome=None, search=None, size=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `biome` | string | None | Biome lineage identifier, e.g., &#x27;root:Host-associated&#x27; or &#x27;root:Host-associated:Human:Digestive system&#x27; |
| `search` | string | None | Keyword to search in study title/description |
| `size` | integer | None | Limit for the returned page; does not automatically retrieve later pages. {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | Progress flag. |


## Search global plant names and taxonomy {#dnakit.databases.powo_search_plants-parameters}

`dnakit.databases.powo_search_plants(*, q, f=None, perPage=None, cursor=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `q` | string | Required | Search query. Can be scientific name, common name, genus, or family. Examples: &#x27;rosa canina&#x27;, &#x27;orchid&#x27;, &#x27;Asteraceae&#x27;, &#x27;oak&#x27; |
| `f` | string &#124; null | None | Filter facets (comma-separated). Common filters: &#x27;species_f:true&#x27; (species only), &#x27;accepted_names:true&#x27; (only accepted names). Example: &#x27;species_f:true&#x27; |
| `perPage` | integer &#124; null | None | Results per page (default 24, max 500) |
| `cursor` | string &#124; null | None | Pagination cursor from previous response |
| `verbose` | bool | True | Progress flag. |


## Retrieve rat gene disease, phenotype and functional annotations {#dnakit.databases.rgd_get_annotations-parameters}

`dnakit.databases.rgd_get_annotations(*, rgd_id, aspect=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `rgd_id` | string | Required | RGD gene ID (e.g., &#x27;2218&#x27; for Brca1) |
| `aspect` | string | None | Filter by annotation type: &#x27;D&#x27; (Disease), &#x27;P&#x27; (Pathway), &#x27;W&#x27; (Phenotype), &#x27;F&#x27; (Function), &#x27;C&#x27; (Component), &#x27;E&#x27; (Expression). Leave empty for all. {&quot;default&quot;: &quot;&quot;} |
| `verbose` | bool | True | Progress flag. |


## Retrieve rat gene details by RGD identifier {#dnakit.databases.rgd_get_gene-parameters}

`dnakit.databases.rgd_get_gene(*, rgd_id, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `rgd_id` | string | Required | RGD gene ID (numeric, e.g., &#x27;2218&#x27; for Brca1). Can include &#x27;RGD:&#x27; prefix. |
| `verbose` | bool | True | Progress flag. |


## Retrieve cross-species orthologs of a rat gene {#dnakit.databases.rgd_get_orthologs-parameters}

`dnakit.databases.rgd_get_orthologs(*, rgd_id, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `rgd_id` | string | Required | RGD gene ID (e.g., &#x27;2218&#x27; for Brca1) |
| `verbose` | bool | True | Progress flag. |


## Retrieve rat QTLs overlapping a genomic interval {#dnakit.databases.rgd_get_qtls_in_region-parameters}

`dnakit.databases.rgd_get_qtls_in_region(*, chromosome, start, stop, map_key, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `chromosome` | string | Required | Chromosome name, e.g. &#x27;10&#x27; or &#x27;chr10&#x27;. |
| `start` | integer | Required | Region start coordinate (1-based). |
| `stop` | integer | Required | Region stop coordinate. |
| `map_key` | integer | Required | Assembly map key. 360=rat rn7/GRCr8, 372=rat mRatBN7.2, 38=human GRCh38, 17=human GRCh37, 35=mouse GRCm39, 18=mouse GRCm38. |
| `verbose` | bool | True | Progress flag. |


## Retrieve rat strain details by strain symbol {#dnakit.databases.rgd_get_strain-parameters}

`dnakit.databases.rgd_get_strain(*, symbol, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `symbol` | string | Required | Exact rat strain symbol (case-insensitive), e.g. &#x27;SHR&#x27; (spontaneously hypertensive rat), &#x27;BN&#x27; (Brown Norway), &#x27;SS&#x27; (Dahl salt-sensitive). |
| `verbose` | bool | True | Progress flag. |


## Retrieve rat strain disease and phenotype annotations {#dnakit.databases.rgd_get_strain_annotations-parameters}

`dnakit.databases.rgd_get_strain_annotations(*, symbol, category=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `symbol` | string | Required | Exact rat strain symbol (case-insensitive), e.g. &#x27;SHR&#x27;, &#x27;BN&#x27;, &#x27;SS&#x27;. |
| `category` | string &#124; null | None | Optional annotation category filter: &#x27;disease&#x27;, &#x27;phenotype&#x27;, &#x27;strain&#x27;, &#x27;pathway&#x27;, etc. Omit to return all categories. |
| `verbose` | bool | True | Progress flag. |


## Search RGD genes by name or keyword {#dnakit.databases.rgd_search_genes-parameters}

`dnakit.databases.rgd_search_genes(*, query, species=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required | Gene symbol or keyword (e.g., &#x27;Brca1&#x27;, &#x27;Tp53&#x27;, &#x27;insulin&#x27;) |
| `species` | string | None | Species to search: &#x27;rat&#x27; (default), &#x27;human&#x27;, &#x27;mouse&#x27; {&quot;default&quot;: &quot;rat&quot;} |
| `verbose` | bool | True | Progress flag. |


## Search rat strains by keyword and type {#dnakit.databases.rgd_search_strains-parameters}

`dnakit.databases.rgd_search_strains(*, query=None, strain_type=None, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string &#124; null | None | Keyword to match against strain symbol, name, genetics, and description (case-insensitive substring). Example: &#x27;hypertensive&#x27; or &#x27;Dahl&#x27;. |
| `strain_type` | string &#124; null | None | Filter by strain type, e.g. &#x27;inbred&#x27;, &#x27;congenic&#x27;, &#x27;transgenic&#x27;, &#x27;recombinant_inbred&#x27;, &#x27;consomic&#x27;. |
| `limit` | integer &#124; null | None | Maximum number of strains to return (default 25, max 200). |
| `verbose` | bool | True | Progress flag. |


## Retrieve yeast gene and locus overviews {#dnakit.databases.sgd_get_gene-parameters}

`dnakit.databases.sgd_get_gene(*, sgd_id, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `sgd_id` | string | Required | SGD locus identifier. Examples: &#x27;S000000259&#x27; (GAL4), &#x27;S000003219&#x27; (RMR1), &#x27;S000005564&#x27; (ACT1). |
| `verbose` | bool | True | Progress flag. |


## Retrieve yeast gene GO annotations and evidence {#dnakit.databases.sgd_get_go_annotations-parameters}

`dnakit.databases.sgd_get_go_annotations(*, sgd_id, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `sgd_id` | string | Required | SGD locus identifier. Examples: &#x27;S000003219&#x27; (RMR1), &#x27;S000005564&#x27; (ACT1), &#x27;S000000259&#x27; (GAL4). |
| `verbose` | bool | True | Progress flag. |


## Retrieve yeast genetic and physical interactions {#dnakit.databases.sgd_get_interactions-parameters}

`dnakit.databases.sgd_get_interactions(*, sgd_id, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `sgd_id` | string | Required | SGD locus identifier. Examples: &#x27;S000003219&#x27; (RMR1), &#x27;S000005564&#x27; (ACT1), &#x27;S000000259&#x27; (GAL4). |
| `verbose` | bool | True | Progress flag. |


## Retrieve yeast mutant phenotypes and experimental conditions {#dnakit.databases.sgd_get_phenotypes-parameters}

`dnakit.databases.sgd_get_phenotypes(*, sgd_id, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `sgd_id` | string | Required | SGD locus identifier. Examples: &#x27;S000003219&#x27; (RMR1), &#x27;S000000259&#x27; (GAL4), &#x27;S000005564&#x27; (ACT1). |
| `verbose` | bool | True | Progress flag. |


## Search the yeast database across data types {#dnakit.databases.sgd_search-parameters}

`dnakit.databases.sgd_search(*, query, category=None, limit=None, offset=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required | Search query string. Examples: &#x27;GAL4&#x27;, &#x27;cell cycle&#x27;, &#x27;sporulation&#x27;. |
| `category` | string | None | Filter by category. Options: &#x27;locus&#x27;, &#x27;reference&#x27;, &#x27;biological_process&#x27;, &#x27;molecular_function&#x27;, &#x27;cellular_component&#x27;, &#x27;phenotype&#x27;, &#x27;disease&#x27;, &#x27;allele&#x27;, &#x27;complex&#x27;, &#x27;pathway&#x27;, &#x27;chemical&#x27;, &#x27;strain&#x27;, &#x27;download&#x27;. |
| `limit` | integer | None | Maximum results to return (default 10, max 50). |
| `offset` | integer | None | Pagination offset (default 0). |
| `verbose` | bool | True | Progress flag. |


## Retrieve plant morphology, growth and reproduction characteristics {#dnakit.databases.usda_plants_get_characteristics-parameters}

`dnakit.databases.usda_plants_get_characteristics(*, symbol, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `symbol` | string | Required | USDA PLANTS symbol, e.g. &#x27;ABBA&#x27;, &#x27;ACRU&#x27;. |
| `verbose` | bool | True | Progress flag. |


## Retrieve plant names, taxonomy and growth habits {#dnakit.databases.usda_plants_get_profile-parameters}

`dnakit.databases.usda_plants_get_profile(*, symbol, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `symbol` | string | Required | USDA PLANTS symbol, e.g. &#x27;ABBA&#x27; (balsam fir), &#x27;ACRU&#x27; (red maple), &#x27;QUAL&#x27; (white oak). |
| `verbose` | bool | True | Progress flag. |


## Retrieve nematode tissue and developmental expression {#dnakit.databases.worm_base_get_expression-parameters}

`dnakit.databases.worm_base_get_expression(*, gene_id, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_id` | string | Required | WormBase gene identifier. Examples: &#x27;WBGene00006763&#x27; (unc-26), &#x27;WBGene00000898&#x27; (daf-2). |
| `verbose` | bool | True | Progress flag. |


## Retrieve C. elegans gene details {#dnakit.databases.worm_base_get_gene-parameters}

`dnakit.databases.worm_base_get_gene(*, gene_id, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_id` | string | Required | WormBase gene identifier (WBGene format). Examples: &#x27;WBGene00006763&#x27; (unc-26), &#x27;WBGene00006789&#x27; (unc-52), &#x27;WBGene00000898&#x27; (daf-2). |
| `verbose` | bool | True | Progress flag. |


## Retrieve experimentally observed nematode gene phenotypes {#dnakit.databases.worm_base_get_phenotypes-parameters}

`dnakit.databases.worm_base_get_phenotypes(*, gene_id, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_id` | string | Required | WormBase gene identifier. Examples: &#x27;WBGene00006763&#x27; (unc-26), &#x27;WBGene00000898&#x27; (daf-2). |
| `verbose` | bool | True | Progress flag. |


## Retrieve Xenopus gene details and cross-references {#dnakit.databases.xenbase_get_gene-parameters}

`dnakit.databases.xenbase_get_gene(*, gene_id, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_id` | string | Required | Xenbase gene ID. Can be with or without &#x27;Xenbase:&#x27; prefix. Examples: &#x27;XB-GENE-484286&#x27; (tp53, X. tropicalis), &#x27;XB-GENE-6252168&#x27; (tp53.L, X. laevis), &#x27;XB-GENE-487192&#x27; (pax6, X. tropicalis). |
| `verbose` | bool | True | Progress flag. |


## Search genes in Xenopus tropicalis and Xenopus laevis {#dnakit.databases.xenbase_search_genes-parameters}

`dnakit.databases.xenbase_search_genes(*, query, species=None, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required | Gene symbol or name to search. Examples: &#x27;tp53&#x27;, &#x27;pax6&#x27;, &#x27;sox2&#x27;, &#x27;bmp4&#x27;, &#x27;wnt3a&#x27;, &#x27;fgf8&#x27;. |
| `species` | string | None | Optional species filter: &#x27;tropicalis&#x27; for X. tropicalis, &#x27;laevis&#x27; for X. laevis. Leave empty for both species. {&quot;default&quot;: &quot;&quot;} |
| `limit` | integer | None | Maximum number of results to return. Default: 10. {&quot;default&quot;: 10} |
| `verbose` | bool | True | Progress flag. |


## Retrieve zebrafish gene details {#dnakit.databases.zfin_get_gene-parameters}

`dnakit.databases.zfin_get_gene(*, gene_id, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_id` | string | Required | ZFIN gene ID with &#x27;ZFIN:&#x27; prefix. Examples: &#x27;ZFIN:ZDB-GENE-990415-8&#x27; (pax2a), &#x27;ZFIN:ZDB-GENE-980526-166&#x27; (shha/sonic hedgehog a), &#x27;ZFIN:ZDB-GENE-990415-72&#x27; (tp53/tumor protein p53), &#x27;ZFIN:ZDB-GENE-030131-7740&#x27; (bmp4/bone morphogenetic protein 4). |
| `verbose` | bool | True | Progress flag. |


## Retrieve zebrafish developmental and anatomical expression {#dnakit.databases.zfin_get_gene_expression-parameters}

`dnakit.databases.zfin_get_gene_expression(*, gene_id, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_id` | string | Required | ZFIN gene ID, with or without the &#x27;ZFIN:&#x27; prefix. Examples: &#x27;ZFIN:ZDB-GENE-990415-8&#x27; (pax2a), &#x27;ZFIN:ZDB-GENE-980526-166&#x27; (shha), &#x27;ZFIN:ZDB-GENE-990415-72&#x27; (tp53). |
| `limit` | integer | None | Max annotations to return (default 50, max 100). {&quot;default&quot;: 50} |
| `verbose` | bool | True | Progress flag. |


## Retrieve zebrafish gene mutant phenotypes {#dnakit.databases.zfin_get_gene_phenotypes-parameters}

`dnakit.databases.zfin_get_gene_phenotypes(*, gene_id, limit=None, page=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_id` | string | Required | ZFIN gene ID with &#x27;ZFIN:&#x27; prefix. Examples: &#x27;ZFIN:ZDB-GENE-990415-8&#x27; (pax2a), &#x27;ZFIN:ZDB-GENE-990415-72&#x27; (tp53), &#x27;ZFIN:ZDB-GENE-980526-166&#x27; (shha). |
| `limit` | integer | None | Maximum number of phenotypes to return (1-100). Default: 20. |
| `page` | integer | None | Page number for pagination. Default: 1. |
| `verbose` | bool | True | Progress flag. |
