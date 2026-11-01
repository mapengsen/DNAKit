# P44 Public database operations

Execute real queries and downloads against fixed official endpoints. Ordinary tests and the replay example below use no network; they do not claim per-record live validation.

Downloads inherit byte limits, atomic writes, checksums and progress. Results are JSON serialisable and retain responses, retrieval dates and page information. Follow upstream page tokens explicitly; list release files before choosing an archived dataset.

## `datasets_gene_id_download_summary()`

- **Purpose:** Gene data package download summary (preview) by GeneID. NCBI Datasets v2; one page,
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_gene_id_download_summary(gene_ids: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_gene_id_expression_report()`

- **Purpose:** Gene expression report by GeneID. NCBI Datasets v2; one page, preserving next_page_token.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_gene_id_expression_report(gene_ids: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_gene_id_links()`

- **Purpose:** Gene external resource links by GeneID. NCBI Datasets v2; one page, preserving
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_gene_id_links(gene_ids: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_gene_product_report()`

- **Purpose:** Gene product report by GeneID, RefSeq accession, symbol+taxon or locus tag. NCBI
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_gene_product_report(identifier: 'str', *, by: 'str' = 'id', taxon: 'str | None' = None, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_gene_taxon_annotation_chromosome_summary()`

- **Purpose:** Gene counts per chromosome by taxon and annotation name. NCBI Datasets v2; one page,
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_gene_taxon_annotation_chromosome_summary(taxon: 'str', annotation_name: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_gene_taxon_counts()`

- **Purpose:** Gene counts by taxon. NCBI Datasets v2; one page, preserving next_page_token.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_gene_taxon_counts(taxon: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_gene_taxon_dataset_report()`

- **Purpose:** Gene data report by taxon (search). NCBI Datasets v2; one page, preserving
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_gene_taxon_dataset_report(taxon: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_gene_taxon_product_report()`

- **Purpose:** Gene product report by taxon (search). NCBI Datasets v2; one page, preserving
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_gene_taxon_product_report(taxon: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_genome_accession_annotation_report()`

- **Purpose:** Genome annotation report (gene records) by accession. NCBI Datasets v2; one page,
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_genome_accession_annotation_report(accession: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_genome_accession_annotation_report_download_summary()`

- **Purpose:** Genome annotation data package download summary (preview) by accession. NCBI Datasets
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_genome_accession_annotation_report_download_summary(accession: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_genome_accession_annotation_summary()`

- **Purpose:** Genome annotation report summary (facets) by accession. NCBI Datasets v2; one page,
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_genome_accession_annotation_summary(accession: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_genome_accession_check()`

- **Purpose:** Genome assembly accession validity check. NCBI Datasets v2; one page, preserving
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_genome_accession_check(accessions: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_genome_accession_download_summary()`

- **Purpose:** Genome data package download summary (preview) by accession. NCBI Datasets v2; one page,
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_genome_accession_download_summary(accessions: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_genome_accession_links()`

- **Purpose:** Genome assembly external resource links by accession. NCBI Datasets v2; one page,
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_genome_accession_links(accessions: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_genome_accession_revision_history()`

- **Purpose:** Genome assembly revision history by accession. NCBI Datasets v2; one page, preserving
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_genome_accession_revision_history(accession: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_genome_report()`

- **Purpose:** Genome assembly dataset report by identifier (accession, WGS, assembly name). NCBI
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_genome_report(identifier: 'str', *, by: 'str' = 'accession', parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_genome_sequence_accession_sequence_assemblies()`

- **Purpose:** Genome assembly accession for a nucleotide sequence accession. NCBI Datasets v2; one
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_genome_sequence_accession_sequence_assemblies(accession: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_genome_taxon_checkm_histogram()`

- **Purpose:** CheckM quality histogram of genome assemblies by species taxon. NCBI Datasets v2; one
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_genome_taxon_checkm_histogram(species_taxon: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_organelle_accession_download()`

- **Purpose:** Organelle genome data package download by nucleotide accession. NCBI Datasets v2; one
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_organelle_accession_download(accessions: 'str', *, output_path: 'str', overwrite: 'bool' = False, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_organelle_accessions_dataset_report()`

- **Purpose:** Organelle genome data report by nucleotide accession. NCBI Datasets v2; one page,
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_organelle_accessions_dataset_report(accessions: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_organelle_taxon_dataset_report()`

- **Purpose:** Organelle genome data report by taxon (search). NCBI Datasets v2; one page, preserving
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_organelle_taxon_dataset_report(taxons: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_protein_accession_download()`

- **Purpose:** Prokaryote gene data package download by RefSeq protein accession. NCBI Datasets v2; one
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_protein_accession_download(accessions: 'str', *, output_path: 'str', overwrite: 'bool' = False, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_virus_accession_annotation_report()`

- **Purpose:** Virus annotation report by nucleotide accession. NCBI Datasets v2; one page, preserving
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_virus_accession_annotation_report(accessions: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_virus_accession_check()`

- **Purpose:** Virus genome nucleotide accession validity check. NCBI Datasets v2; one page, preserving
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_virus_accession_check(accessions: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_virus_download_summary()`

- **Purpose:** Virus genome data package download summary (preview) by taxon or SARS-CoV-2 protein.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_virus_download_summary(identifier: 'str', *, by: 'str' = 'taxon', parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `datasets_virus_taxon_annotation_report()`

- **Purpose:** Virus annotation report by taxon (search). NCBI Datasets v2; one page, preserving
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `datasets_virus_taxon_annotation_report(taxon: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_ga4gh_beacon_query()`

- **Purpose:** GA4GH Beacon allele presence query. Official GET /ga4gh/beacon/query; upstream defaults
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_ga4gh_beacon_query(*, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_ga4gh_callsets_record()`

- **Purpose:** GA4GH callset by ID. Official GET /ga4gh/callsets/:id; upstream defaults retained.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_ga4gh_callsets_record(id: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_ga4gh_callsets_search()`

- **Purpose:** GA4GH genotype callsets search. Official POST /ga4gh/callsets/search; returns the
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_ga4gh_callsets_search(query: 'Mapping[str, Any]', *, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_ga4gh_datasets_record()`

- **Purpose:** GA4GH dataset by ID. Official GET /ga4gh/datasets/:id; upstream defaults retained.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_ga4gh_datasets_record(id: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_ga4gh_datasets_search()`

- **Purpose:** GA4GH datasets search. Official POST /ga4gh/datasets/search; returns the service page
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_ga4gh_datasets_search(query: 'Mapping[str, Any]', *, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_ga4gh_features_record()`

- **Purpose:** GA4GH feature record by ID. Official GET /ga4gh/features/:id; upstream defaults retained.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_ga4gh_features_record(id: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_ga4gh_features_search()`

- **Purpose:** GA4GH sequence features in a region. Official POST /ga4gh/features/search; returns the
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_ga4gh_features_search(query: 'Mapping[str, Any]', *, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_ga4gh_featuresets_record()`

- **Purpose:** GA4GH feature set by ID. Official GET /ga4gh/featuresets/:id; upstream defaults retained.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_ga4gh_featuresets_record(id: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_ga4gh_featuresets_search()`

- **Purpose:** GA4GH feature sets search. Official POST /ga4gh/featuresets/search; returns the service
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_ga4gh_featuresets_search(query: 'Mapping[str, Any]', *, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_ga4gh_references_record()`

- **Purpose:** GA4GH reference sequence by ID. Official GET /ga4gh/references/:id; upstream defaults
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_ga4gh_references_record(id: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_ga4gh_references_search()`

- **Purpose:** GA4GH reference sequences search. Official POST /ga4gh/references/search; returns the
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_ga4gh_references_search(query: 'Mapping[str, Any]', *, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_ga4gh_referencesets_record()`

- **Purpose:** GA4GH reference set by ID. Official GET /ga4gh/referencesets/:id; upstream defaults
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_ga4gh_referencesets_record(id: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_ga4gh_referencesets_search()`

- **Purpose:** GA4GH reference sets search. Official POST /ga4gh/referencesets/search; returns the
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_ga4gh_referencesets_search(query: 'Mapping[str, Any]', *, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_ga4gh_variantannotations_search()`

- **Purpose:** GA4GH variant annotations in a region. Official POST /ga4gh/variantannotations/search;
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_ga4gh_variantannotations_search(query: 'Mapping[str, Any]', *, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_ga4gh_variantannotationsets_record()`

- **Purpose:** GA4GH variant annotation set by ID. Official GET /ga4gh/variantannotationsets/:id;
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_ga4gh_variantannotationsets_record(id: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_ga4gh_variantannotationsets_search()`

- **Purpose:** GA4GH variant annotation sets search. Official POST /ga4gh/variantannotationsets/search;
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_ga4gh_variantannotationsets_search(query: 'Mapping[str, Any]', *, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_ga4gh_variants_record()`

- **Purpose:** GA4GH variant record by ID. Official GET /ga4gh/variants/:id; upstream defaults retained.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_ga4gh_variants_record(id: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_ga4gh_variants_search()`

- **Purpose:** GA4GH variant calls in a region. Official POST /ga4gh/variants/search; returns the
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_ga4gh_variants_search(query: 'Mapping[str, Any]', *, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_ga4gh_variantsets_record()`

- **Purpose:** GA4GH variant set by ID. Official GET /ga4gh/variantsets/:id; upstream defaults retained.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_ga4gh_variantsets_record(id: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_ga4gh_variantsets_search()`

- **Purpose:** GA4GH variant sets search. Official POST /ga4gh/variantsets/search; returns the service
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_ga4gh_variantsets_search(query: 'Mapping[str, Any]', *, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_info_analysis()`

- **Purpose:** List annotation analyses of a species. Official GET /info/analysis/:species; upstream
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_info_analysis(species: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_info_assembly()`

- **Purpose:** Return top-level assembly sequence metadata from GET info/assembly/:species.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_info_assembly(species: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_info_biotypes()`

- **Purpose:** List gene and transcript biotypes of a species. Official GET /info/biotypes/:species;
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_info_biotypes(species: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_info_biotypes_groups()`

- **Purpose:** List biotype groups. Official GET /info/biotypes/groups/:group/:object_type; upstream
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_info_biotypes_groups(group: 'str', object_type: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_info_biotypes_name()`

- **Purpose:** Biotype properties by name. Official GET /info/biotypes/name/:name/:object_type;
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_info_biotypes_name(name: 'str', object_type: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_info_compara_methods()`

- **Purpose:** List comparative genomics methods. Official GET /info/compara/methods; upstream defaults
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_info_compara_methods(*, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_info_compara_species_sets()`

- **Purpose:** List species sets analysed by a comparative method. Official GET
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_info_compara_species_sets(method: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_info_external_dbs()`

- **Purpose:** List external databases linked to a species. Official GET /info/external_dbs/:species;
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_info_external_dbs(species: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_info_genomes_accession()`

- **Purpose:** Genome record by INSDC accession. Official GET /info/genomes/accession/:accession;
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_info_genomes_accession(accession: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_info_genomes_assembly()`

- **Purpose:** Genome record by assembly ID. Official GET /info/genomes/assembly/:assembly_id; upstream
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_info_genomes_assembly(assembly_id: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_info_genomes_division()`

- **Purpose:** List genomes in a division. Official GET /info/genomes/division/:division_name; upstream
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_info_genomes_division(division_name: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_info_genomes_taxonomy()`

- **Purpose:** List genomes under a taxonomy node. Official GET /info/genomes/taxonomy/:taxon_name;
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_info_genomes_taxonomy(taxon_name: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_info_variation()`

- **Purpose:** List variation sources of a species. Official GET /info/variation/:species; upstream
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_info_variation(species: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_info_variation_consequence_types()`

- **Purpose:** List variant consequence types. Official GET /info/variation/consequence_types; upstream
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_info_variation_consequence_types(*, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_info_variation_populations()`

- **Purpose:** List individuals of a population. Official GET
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_info_variation_populations(species: 'str', population_name: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_info_variation_populations_list()`

- **Purpose:** List populations of a species. Official GET /info/variation/populations/:species;
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_info_variation_populations_list(species: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_map_cds()`

- **Purpose:** CDS to genomic coordinate conversion. Official GET /map/cds/:id/:region; upstream
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_map_cds(id: 'str', region: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_map_translation()`

- **Purpose:** Protein to genomic coordinate conversion. Official GET /map/translation/:id/:region;
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_map_translation(id: 'str', region: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_overlap_translation()`

- **Purpose:** Features of a translation (domains, variants). Official GET /overlap/translation/:id;
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_overlap_translation(id: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_phenotype_accession()`

- **Purpose:** Phenotype annotations by ontology accession. Official GET
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_phenotype_accession(species: 'str', accession: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_phenotype_gene()`

- **Purpose:** Phenotype annotations of a gene. Official GET /phenotype/gene/:species/:gene; upstream
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_phenotype_gene(species: 'str', gene: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_phenotype_term()`

- **Purpose:** Phenotype annotations by ontology term. Official GET /phenotype/term/:species/:term;
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_phenotype_term(species: 'str', term: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_variant_recoder()`

- **Purpose:** Variant ID, HGVS and SPDI recoding. Official GET /variant_recoder/:species/:id; upstream
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_variant_recoder(species: 'str', id: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_variation_publication()`

- **Purpose:** Variants in a publication by PMID or PMCID. Official GET
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_variation_publication(species: 'str', publication: 'str', *, publication_type: 'str' = 'pmcid', parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ensembl_vep_region()`

- **Purpose:** Variant consequences by genomic region and allele. Official GET
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ensembl_vep_region(species: 'str', region: 'str', allele: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## Offline input validation example

**Example code:**

```python
from dnakit.databases.public_apis_p44_functions import ensembl_info_analysis
from dnakit.exceptions import ConfigurationError
try:
    ensembl_info_analysis("", verbose=False)
except ConfigurationError:
    print("invalid input rejected before network access")
```

**Example output:**

```text
invalid input rejected before network access
```

## Upstream references

- [官方文档 1](https://rest.ensembl.org/)
- [官方文档 2](https://github.com/ncbi/datasets/blob/master/datasets.openapi.yaml)
- [官方文档 3](https://www.ebi.ac.uk/ena/browser/api/swagger-ui/index.html)
