# Gene names and intervals from genomic DNA: parameters and notes

See the [tutorial](../../api/features/26_gene_annotation.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

## Retrieve a cross-species gene family tree {#dnakit.databases.ensembl_compara_get_gene_tree-parameters}

`dnakit.databases.ensembl_compara_get_gene_tree(*, gene, species=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene` | string | Required | Gene symbol or Ensembl gene ID. Examples: &#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;ENSG00000141510&#x27;. Ensembl IDs are more reliable than symbols for gene tree retrieval. |
| `species` | string | None | Species (default: &#x27;human&#x27;). Used when providing a gene symbol. Examples: &#x27;human&#x27;, &#x27;mouse&#x27;. |
| `verbose` | bool | True | Progress flag. |


## Retrieve orthologous genes across species {#dnakit.databases.ensembl_compara_get_orthologues-parameters}

`dnakit.databases.ensembl_compara_get_orthologues(*, gene, species=None, target_species=None, target_taxon=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene` | string | Required | Gene symbol or Ensembl gene ID. Examples: &#x27;BRCA1&#x27;, &#x27;TP53&#x27;, &#x27;ENSG00000141510&#x27;, &#x27;ENSG00000012048&#x27;. |
| `species` | string | None | Source species (default: &#x27;human&#x27;). Examples: &#x27;human&#x27;, &#x27;mouse&#x27;, &#x27;zebrafish&#x27;, &#x27;drosophila_melanogaster&#x27;. |
| `target_species` | string | None | Target species to limit orthologue search. Examples: &#x27;mouse&#x27;, &#x27;zebrafish&#x27;, &#x27;rat&#x27;, &#x27;drosophila_melanogaster&#x27;. Default: all species. |
| `target_taxon` | integer | None | NCBI taxon ID to limit orthologue search. Examples: 10090 (mouse), 7955 (zebrafish), 7227 (fruit fly). Alternative to target_species. |
| `verbose` | bool | True | Progress flag. |


## Retrieve paralogous genes {#dnakit.databases.ensembl_compara_get_paralogues-parameters}

`dnakit.databases.ensembl_compara_get_paralogues(*, gene, species=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene` | string | Required | Gene symbol or Ensembl gene ID. Examples: &#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;ENSG00000141510&#x27;. |
| `species` | string | None | Species (default: &#x27;human&#x27;). Examples: &#x27;human&#x27;, &#x27;mouse&#x27;, &#x27;zebrafish&#x27;. |
| `verbose` | bool | True | Progress flag. |


## Retrieve gene homology by symbol {#dnakit.databases.ensembl_get_homology-parameters}

`dnakit.databases.ensembl_get_homology(*, species, symbol, target_species=None, target_taxon=None, type=None, sequence=None, aligned=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `species` | string | Required | Source species (e.g., &#x27;human&#x27;, &#x27;homo_sapiens&#x27;) |
| `symbol` | string | Required | Gene symbol (e.g., &#x27;BRCA2&#x27;, &#x27;TP53&#x27;) |
| `target_species` | string | None | Target species to find homologues in (optional, e.g., &#x27;mouse&#x27;, &#x27;mus_musculus&#x27;). Leave empty to search all species. |
| `target_taxon` | string | None | Target taxonomic group (optional, e.g., &#x27;Mammalia&#x27;, &#x27;Primates&#x27;) |
| `type` | string | None | Homology type: &#x27;orthologues&#x27; (different species), &#x27;paralogues&#x27; (same species), &#x27;all&#x27; (both) {&quot;enum&quot;: &#91;&quot;orthologues&quot;, &quot;paralogues&quot;, &quot;all&quot;&#93;, &quot;default&quot;: &quot;all&quot;} |
| `sequence` | string | None | Include sequences in response {&quot;enum&quot;: &#91;&quot;cdna&quot;, &quot;protein&quot;, &quot;none&quot;&#93;, &quot;default&quot;: &quot;none&quot;} |
| `aligned` | boolean | None | Return aligned sequences {&quot;default&quot;: false} |
| `verbose` | bool | True | Progress flag. |


## Retrieve annotations overlapping a genomic interval {#dnakit.databases.ensembl_get_overlap_features-parameters}

`dnakit.databases.ensembl_get_overlap_features(*, species, region, feature=None, biotype=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `species` | string | Required | Species (e.g., &#x27;human&#x27;, &#x27;homo_sapiens&#x27;) |
| `region` | string | Required | Genomic region (e.g., &#x27;7:140424943-140624564&#x27;) |
| `feature` | string | None | Feature type to retrieve {&quot;enum&quot;: &#91;&quot;gene&quot;, &quot;transcript&quot;, &quot;cds&quot;, &quot;exon&quot;, &quot;repeat&quot;, &quot;simple_feature&quot;, &quot;misc_feature&quot;, &quot;constrained_element&quot;&#93;, &quot;default&quot;: &quot;gene&quot;} |
| `biotype` | string | None | Filter by biotype (optional, e.g., &#x27;protein_coding&#x27;, &#x27;lincRNA&#x27;, &#x27;miRNA&#x27;) |
| `verbose` | bool | True | Progress flag. |


## Retrieve regulatory features in a genomic interval {#dnakit.databases.ensembl_get_regulatory_features-parameters}

`dnakit.databases.ensembl_get_regulatory_features(*, region, species=None, feature, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `region` | string | Required | Genomic region in format &#x27;chromosome:start-end&#x27; (e.g., &#x27;1:1000000-2000000&#x27;) |
| `species` | string | None | Species name (default &#x27;human&#x27;) {&quot;default&quot;: &quot;human&quot;} |
| `feature` | string | Required | Feature type (must be &#x27;regulatory&#x27; for this tool) |
| `verbose` | bool | True | Progress flag. |


## Look up Ensembl genes by identifier or symbol {#dnakit.databases.ensembl_lookup_gene-parameters}

`dnakit.databases.ensembl_lookup_gene(*, gene_id=None, species=None, id=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_id` | string | None | Ensembl gene ID or symbol (e.g., &#x27;ENSG00000139618&#x27; or &#x27;BRCA1&#x27;). If using a stable ID, the tool will automatically route to /lookup/id endpoint. |
| `species` | string | None | Species name required for gene symbols (default &#x27;homo_sapiens&#x27;). Examples: &#x27;homo_sapiens&#x27;, &#x27;mus_musculus&#x27;, &#x27;rattus_norvegicus&#x27; {&quot;default&quot;: &quot;homo_sapiens&quot;} |
| `id` | string | None | DNAKit optional alias for gene_id |
| `verbose` | bool | True | Progress flag. |


## Retrieve HGNC annotations by official gene symbol {#dnakit.databases.hgnc_fetch_gene_by_symbol-parameters}

`dnakit.databases.hgnc_fetch_gene_by_symbol(*, symbol, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `symbol` | string | Required | Official HGNC gene symbol. Examples: &#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;EGFR&#x27;, &#x27;INS&#x27;, &#x27;ACE2&#x27;. |
| `verbose` | bool | True | Progress flag. |


## Batch-query gene annotations and identifier mappings {#dnakit.databases.my_gene_batch_query-parameters}

`dnakit.databases.my_gene_batch_query(*, gene_ids, species=None, fields=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_ids` | array | Required | List of gene IDs to query. Can mix different ID types (Entrez, Ensembl, symbols). {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `species` | string | None | Species filter for the query. {&quot;default&quot;: &quot;human&quot;} |
| `fields` | string | None | Comma-separated list of fields to return for each gene. {&quot;default&quot;: &quot;symbol,name,entrezgene,ensembl.gene&quot;} |
| `verbose` | bool | True | Progress flag. |


## Retrieve comprehensive annotations for a gene identifier {#dnakit.databases.my_gene_get_gene_annotation-parameters}

`dnakit.databases.my_gene_get_gene_annotation(*, gene_id, fields=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_id` | string | Required | Gene ID to query. Can be Entrez Gene ID (e.g., &#x27;1017&#x27; for CDK2) or Ensembl ID (e.g., &#x27;ENSG00000123374&#x27;). |
| `fields` | string | None | Comma-separated list of fields to return. Available fields include: symbol, name, summary, go (gene ontology), pathway (KEGG/Reactome), interpro (protein domains), generif (gene references), pdb, uniprot. {&quot;default&quot;: &quot;symbol,name,entrezgene,ensembl,summary,go,pathway,interpro&quot;} |
| `verbose` | bool | True | Progress flag. |


## Search genes by keyword or identifier {#dnakit.databases.my_gene_query_genes-parameters}

`dnakit.databases.my_gene_query_genes(*, query=None, species=None, fields=None, size=None, q=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | None | Search query. Can be gene symbol (e.g., &#x27;CDK2&#x27;), name (&#x27;cyclin dependent kinase&#x27;), Entrez ID (&#x27;1017&#x27;), or Ensembl ID (&#x27;ENSG00000123374&#x27;). Supports wildcards (*) and boolean operators (AND, OR, NOT). |
| `species` | string | None | Species filter. Use common name or NCBI taxonomy ID. Examples: &#x27;human&#x27;, &#x27;mouse&#x27;, &#x27;9606&#x27; (human), &#x27;all&#x27;. {&quot;default&quot;: &quot;human&quot;} |
| `fields` | string | None | Comma-separated list of fields to return. Common fields: symbol, name, entrezgene, ensembl.gene, summary, go, pathway, interpro. {&quot;default&quot;: &quot;symbol,name,entrezgene,ensembl.gene,summary&quot;} |
| `size` | integer | None | Maximum number of results to return (1-100). {&quot;default&quot;: 10} |
| `q` | string | None | DNAKit optional alias for query |
| `verbose` | bool | True | Progress flag. |


## Retrieve gene details by NCBI Gene ID {#dnakit.databases.ncbidatasets_get_gene-parameters}

`dnakit.databases.ncbidatasets_get_gene(*, gene_id, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_id` | string | Required | NCBI Gene ID (numeric). Examples: &#x27;7157&#x27; (TP53), &#x27;672&#x27; (BRCA1), &#x27;3630&#x27; (INS), &#x27;22059&#x27; (mouse Trp53). {&quot;pattern&quot;: &quot;^&#91;1-9&#93;&#91;0-9&#93;*$&quot;} |
| `verbose` | bool | True | Progress flag. |


## Look up NCBI genes by symbol and organism {#dnakit.databases.ncbidatasets_get_gene_by_symbol-parameters}

`dnakit.databases.ncbidatasets_get_gene_by_symbol(*, symbol, taxon=None, organism=None, species=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `symbol` | string | Required | Gene symbol. Examples: &#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;INS&#x27;, &#x27;EGFR&#x27;. {&quot;minLength&quot;: 1, &quot;pattern&quot;: &quot;.*\\S.*&quot;} |
| `taxon` | string | None | Organism as common name or taxonomy ID. Examples: &#x27;human&#x27;, &#x27;mouse&#x27;, &#x27;rat&#x27;, &#x27;9606&#x27;, &#x27;10090&#x27;. Default: &#x27;human&#x27; (reported as metadata.query_taxon_defaulted). Synonyms `organism` and `species` are accepted. {&quot;minLength&quot;: 1, &quot;pattern&quot;: &quot;.*\\S.*&quot;} |
| `organism` | string | None | Synonym for `taxon`. Used when `taxon` is not supplied. {&quot;minLength&quot;: 1, &quot;pattern&quot;: &quot;.*\\S.*&quot;} |
| `species` | string | None | Synonym for `taxon`. Used when neither `taxon` nor `organism` is supplied. {&quot;minLength&quot;: 1, &quot;pattern&quot;: &quot;.*\\S.*&quot;} |
| `verbose` | bool | True | Progress flag. |


## Retrieve orthologs for an NCBI Gene ID {#dnakit.databases.ncbidatasets_get_orthologs-parameters}

`dnakit.databases.ncbidatasets_get_orthologs(*, gene_id, page_size=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_id` | string | Required | NCBI Gene ID (numeric). Examples: &#x27;7157&#x27; (TP53), &#x27;672&#x27; (BRCA1), &#x27;3630&#x27; (INS). {&quot;pattern&quot;: &quot;^&#91;1-9&#93;&#91;0-9&#93;*$&quot;} |
| `page_size` | integer | None | Maximum number of orthologs to return (1-100). Default: 20. {&quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | Progress flag. |


## Retrieve NCBI gene summaries and functional information {#dnakit.databases.ncbigene_get_summary-parameters}

`dnakit.databases.ncbigene_get_summary(*, id, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `id` | string | Required | NCBI Gene ID (Entrez Gene ID). Can be a single ID or comma-separated list. Examples: &#x27;7157&#x27; (TP53), &#x27;672&#x27; (BRCA1), &#x27;1956&#x27; (EGFR), &#x27;3845&#x27; (KRAS). Get IDs from NCBIGene_search. |
| `verbose` | bool | True | Progress flag. |


## Search NCBI Gene IDs by organism and keyword {#dnakit.databases.ncbigene_search-parameters}

`dnakit.databases.ncbigene_search(*, term, retmax=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `term` | string | Required | Search query. Use gene symbol with organism filter for best results. Format: &#x27;GENE&#91;Symbol&#93; AND Organism&#91;Organism&#93;&#x27;. Examples: &#x27;TP53&#91;Symbol&#93; AND Homo sapiens&#91;Organism&#93;&#x27;, &#x27;BRCA1&#91;Symbol&#93; AND Homo sapiens&#91;Organism&#93;&#x27;, &#x27;insulin AND Homo sapiens&#91;Organism&#93;&#x27;, &#x27;EGFR AND Mus musculus&#91;Organism&#93;&#x27;. |
| `retmax` | integer | None | Maximum number of results to return (default 10, max 100). |
| `verbose` | bool | True | Progress flag. |


## Retrieve MANE transcripts and RefSeq mappings {#dnakit.databases.tark_get_mane_transcripts-parameters}

`dnakit.databases.tark_get_mane_transcripts(*, gene=None, ensembl_id=None, refseq_id=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene` | string &#124; null | None | HGNC gene symbol, e.g. &#x27;BRCA2&#x27;, &#x27;CFTR&#x27;. |
| `ensembl_id` | string &#124; null | None | Ensembl transcript id, e.g. &#x27;ENST00000380152&#x27; (version optional). |
| `refseq_id` | string &#124; null | None | RefSeq transcript id, e.g. &#x27;NM_000059&#x27; (version optional). |
| `verbose` | bool | True | Progress flag. |


## Retrieve archived Ensembl transcript records {#dnakit.databases.tark_get_transcript-parameters}

`dnakit.databases.tark_get_transcript(*, stable_id, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `stable_id` | string | Required | Ensembl transcript stable id, e.g. &#x27;ENST00000380152&#x27; (version optional). |
| `verbose` | bool | True | Progress flag. |


## Retrieve gene transcripts and MANE annotations {#dnakit.databases.variant_validator_gene2transcripts-parameters}

`dnakit.databases.variant_validator_gene2transcripts(*, gene_symbol, gene=None, gene_name=None, transcript_set=None, genome_build=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_symbol` | string | Required | HGNC gene symbol (e.g., &#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;EGFR&#x27;). Aliases: gene, gene_name. |
| `gene` | string | None | Alias for gene_symbol. |
| `gene_name` | string | None | Alias for gene_symbol. |
| `transcript_set` | string | None | Transcript filter: &#x27;mane&#x27; for MANE Select/Plus Clinical only, &#x27;refseq&#x27; for RefSeq transcripts, &#x27;ensembl&#x27; for Ensembl, &#x27;all&#x27; for everything {&quot;default&quot;: &quot;mane&quot;} |
| `genome_build` | string | None | Reference genome assembly: &#x27;GRCh37&#x27; or &#x27;GRCh38&#x27; (default: GRCh38) {&quot;default&quot;: &quot;GRCh38&quot;} |
| `verbose` | bool | True | Progress flag. |
