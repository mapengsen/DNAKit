# Disease, phenotype and ontology queries: parameters and notes

See the [tutorial](../../api/features/22_database_query_phenotypes.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

## Retrieve disease-associated genes and association scores {#dnakit.databases.dis_ge_net_get_disease_genes-parameters}

`dnakit.databases.dis_ge_net_get_disease_genes(*, operation=None, disease, min_score=None, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_disease_genes) {&quot;const&quot;: &quot;get_disease_genes&quot;, &quot;default&quot;: &quot;get_disease_genes&quot;} |
| `disease` | string | Required | Disease ID (UMLS CUI) or disease name |
| `min_score` | number | None | Minimum association score (0-1). Default: no filter {&quot;minimum&quot;: 0, &quot;maximum&quot;: 1} |
| `limit` | integer | None | Maximum genes to return (default: 50) {&quot;default&quot;: 50} |
| `verbose` | bool | True | Progress flag. |

Configure before running: `DISGENET_API_KEY`.


## Filter gene–disease associations by source and score {#dnakit.databases.dis_ge_net_get_gda-parameters}

`dnakit.databases.dis_ge_net_get_gda(*, operation=None, gene=None, disease=None, source=None, min_score=None, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_gda) {&quot;const&quot;: &quot;get_gda&quot;, &quot;default&quot;: &quot;get_gda&quot;} |
| `gene` | string | None | Gene symbol (optional if disease provided) |
| `disease` | string | None | Disease ID/UMLS CUI (optional if gene provided) |
| `source` | string | None | Data source filter: CURATED, ANIMAL_MODELS, LITERATURE, etc. {&quot;enum&quot;: &#91;&quot;CURATED&quot;, &quot;ANIMAL_MODELS&quot;, &quot;LITERATURE&quot;, &quot;ALL&quot;&#93;} |
| `min_score` | number | None | Minimum GDA score (0-1). Higher = stronger evidence. {&quot;minimum&quot;: 0, &quot;maximum&quot;: 1} |
| `limit` | integer | None | Maximum results (default: 25) {&quot;default&quot;: 25} |
| `verbose` | bool | True | Progress flag. |

Configure before running: `DISGENET_API_KEY`.


## Retrieve variant–disease association evidence {#dnakit.databases.dis_ge_net_get_vda-parameters}

`dnakit.databases.dis_ge_net_get_vda(*, operation=None, variant=None, gene=None, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_vda) {&quot;const&quot;: &quot;get_vda&quot;, &quot;default&quot;: &quot;get_vda&quot;} |
| `variant` | string | None | Variant rsID (e.g., rs1234) - optional if gene provided |
| `gene` | string | None | Gene symbol to get all variants - optional if variant provided |
| `limit` | integer | None | Maximum results (default: 25) {&quot;default&quot;: 25} |
| `verbose` | bool | True | Progress flag. |

Configure before running: `DISGENET_API_KEY`.


## Find disease-associated genes by UMLS concept ID {#dnakit.databases.dis_ge_net_search_disease-parameters}

`dnakit.databases.dis_ge_net_search_disease(*, operation=None, disease, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search_disease) {&quot;const&quot;: &quot;search_disease&quot;, &quot;default&quot;: &quot;search_disease&quot;} |
| `disease` | string | Required | UMLS CUI for the disease, e.g. &#x27;C0006142&#x27; (breast cancer). Free-text disease names are NOT accepted -- resolve one via umls_search_concepts first. |
| `limit` | integer | None | Maximum associations to return (default: 10) {&quot;default&quot;: 10} |
| `verbose` | bool | True | Progress flag. |

Configure before running: `DISGENET_API_KEY`.


## Find diseases associated with a gene symbol {#dnakit.databases.dis_ge_net_search_gene-parameters}

`dnakit.databases.dis_ge_net_search_gene(*, operation=None, gene, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search_gene) {&quot;const&quot;: &quot;search_gene&quot;, &quot;default&quot;: &quot;search_gene&quot;} |
| `gene` | string | Required | Gene symbol (e.g., BRCA1, TP53, EGFR) |
| `limit` | integer | None | Maximum associations to return (default: 10) {&quot;default&quot;: 10} |
| `verbose` | bool | True | Progress flag. |

Configure before running: `DISGENET_API_KEY`.


## Summarize GenCC gene–disease validity classifications {#dnakit.databases.gen_cc_get_classifications-parameters}

`dnakit.databases.gen_cc_get_classifications(*, operation=None, gene_symbol=None, disease=None, submitter=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_classifications) {&quot;const&quot;: &quot;get_classifications&quot;, &quot;default&quot;: &quot;get_classifications&quot;} |
| `gene_symbol` | string | None | Optional HGNC gene symbol to filter classification stats (e.g., &#x27;GBA&#x27;, &#x27;BRCA2&#x27;) |
| `disease` | string | None | Optional disease name substring to filter classification stats (e.g., &#x27;Gaucher&#x27;, &#x27;breast cancer&#x27;) |
| `submitter` | string | None | Optional filter by submitting organization name (e.g., &#x27;ClinGen&#x27;, &#x27;Ambry&#x27;, &#x27;Genomics England&#x27;) {&quot;default&quot;: &quot;&quot;} |
| `verbose` | bool | True | Progress flag. |


## Search GenCC gene evidence for a disease {#dnakit.databases.gen_cc_search_disease-parameters}

`dnakit.databases.gen_cc_search_disease(*, operation=None, disease=None, disease_name=None, classification=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search_disease) {&quot;const&quot;: &quot;search_disease&quot;, &quot;default&quot;: &quot;search_disease&quot;} |
| `disease` | string | None | Disease name or ID to search (e.g., &#x27;Marfan syndrome&#x27;, &#x27;MONDO:0007947&#x27;, &#x27;breast cancer&#x27;) |
| `disease_name` | string | None | Alias for disease. Disease name or ID to search. |
| `classification` | string | None | Optional filter by classification level (e.g., Definitive, Strong) {&quot;default&quot;: &quot;&quot;} |
| `verbose` | bool | True | Progress flag. |


## Retrieve GenCC gene–disease classifications {#dnakit.databases.gen_cc_search_gene-parameters}

`dnakit.databases.gen_cc_search_gene(*, operation=None, gene_symbol, classification=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search_gene) {&quot;const&quot;: &quot;search_gene&quot;, &quot;default&quot;: &quot;search_gene&quot;} |
| `gene_symbol` | string | Required | HGNC gene symbol (e.g., BRCA2, TP53, FBN1, CFTR) |
| `classification` | string | None | Optional filter by classification level (e.g., Definitive, Strong, Moderate, Limited) {&quot;default&quot;: &quot;&quot;} |
| `verbose` | bool | True | Progress flag. |


## Retrieve diseases associated with multiple HPO phenotypes {#dnakit.databases.get_joint_associated_diseases_by_hpo_id_list-parameters}

`dnakit.databases.get_joint_associated_diseases_by_hpo_id_list(*, HPO_ID_list, limit=None, offset=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `HPO_ID_list` | array | Required | List of phenotypes or symptoms {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;, &quot;description&quot;: &quot;The HPO ID of the phenotype or symptom.&quot;}} |
| `limit` | integer | None | Number of entries to fetch. |
| `offset` | integer | None | Number of initial entries to skip. |
| `verbose` | bool | True | Progress flag. |


## Retrieve diseases associated with an HPO phenotype {#dnakit.databases.hpo_get_diseases_by_phenotype-parameters}

`dnakit.databases.hpo_get_diseases_by_phenotype(*, term_id=None, limit=None, offset=None, hp_id=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `term_id` | string | None | HPO term identifier (e.g., &#x27;HP:0001250&#x27;). Must start with &#x27;HP:&#x27; followed by 7 digits. |
| `limit` | integer &#124; null | None | Maximum number of results to return (1-500, default 50). Phenotype terms can be linked to thousands of diseases. |
| `offset` | integer &#124; null | None | Number of results to skip before returning &#x27;limit&#x27; results (default 0). Use with the response&#x27;s metadata.has_more/metadata.total to page through phenotypes with more diseases than fit in one call. |
| `hp_id` | string | None | DNAKit optional alias for term_id |
| `verbose` | bool | True | Progress flag. |


## Retrieve genes associated with an HPO phenotype {#dnakit.databases.hpo_get_genes_by_phenotype-parameters}

`dnakit.databases.hpo_get_genes_by_phenotype(*, term_id, limit=None, offset=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `term_id` | string | Required | HPO term identifier (e.g., &#x27;HP:0001250&#x27;). Must start with &#x27;HP:&#x27; followed by 7 digits. |
| `limit` | integer &#124; null | None | Maximum number of results to return (1-500, default 50). Phenotype terms can be linked to thousands of genes. |
| `offset` | integer &#124; null | None | Number of results to skip before returning &#x27;limit&#x27; results (default 0). Use with the response&#x27;s metadata.has_more/metadata.total to page through phenotypes with more genes than fit in one call. |
| `verbose` | bool | True | Progress flag. |


## Retrieve HPO term definitions and synonyms {#dnakit.databases.hpo_get_term-parameters}

`dnakit.databases.hpo_get_term(*, term_id, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `term_id` | string | Required | HPO term identifier. Must start with &#x27;HP:&#x27; followed by 7 digits. Examples: &#x27;HP:0001250&#x27; (Seizure), &#x27;HP:0001263&#x27; (Global developmental delay), &#x27;HP:0000252&#x27; (Microcephaly), &#x27;HP:0001249&#x27; (Intellectual disability). |
| `verbose` | bool | True | Progress flag. |


## Search HPO phenotype terms by keyword {#dnakit.databases.hpo_search_terms-parameters}

`dnakit.databases.hpo_search_terms(*, query, max_results=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required | Search keyword for phenotype terms. Examples: &#x27;seizure&#x27;, &#x27;intellectual disability&#x27;, &#x27;cardiomyopathy&#x27;, &#x27;scoliosis&#x27;, &#x27;microcephaly&#x27;. |
| `max_results` | integer &#124; null | None | Maximum number of results to return (default 10, max 50). |
| `verbose` | bool | True | Progress flag. |


## Aggregate human gene identifiers and annotations {#dnakit.databases.marrvel_get_gene-parameters}

`dnakit.databases.marrvel_get_gene(*, symbol=None, gene_symbol=None, gene=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `symbol` | string &#124; null | None | HGNC gene symbol, e.g. &#x27;CFTR&#x27;, &#x27;BRCA1&#x27;. |
| `gene_symbol` | string &#124; null | None | Alias for symbol. |
| `gene` | string &#124; null | None | Alias for symbol. |
| `verbose` | bool | True | Progress flag. |

Provide at least one of these parameter groups: `symbol`; `gene_symbol`; `gene`.


## Retrieve OMIM phenotypes and inheritance for a gene {#dnakit.databases.marrvel_get_omim_phenotypes-parameters}

`dnakit.databases.marrvel_get_omim_phenotypes(*, symbol=None, gene_symbol=None, gene=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `symbol` | string &#124; null | None | HGNC gene symbol, e.g. &#x27;CFTR&#x27;, &#x27;BRCA1&#x27;. |
| `gene_symbol` | string &#124; null | None | Alias for symbol. |
| `gene` | string &#124; null | None | Alias for symbol. |
| `verbose` | bool | True | Progress flag. |

Provide at least one of these parameter groups: `symbol`; `gene_symbol`; `gene`.


## Retrieve Monarch gene–disease associations {#dnakit.databases.monarch_get_gene_diseases-parameters}

`dnakit.databases.monarch_get_gene_diseases(*, subject, limit=None, offset=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `subject` | string | Required | Gene CURIE identifier (e.g., &#x27;HGNC:11998&#x27; for TP53). Use Monarch_search_gene to find this. |
| `limit` | integer | None | Maximum disease associations to return (default: 100) {&quot;default&quot;: 100} |
| `offset` | integer | None | Number of initial entries to skip. {&quot;default&quot;: 0} |
| `verbose` | bool | True | Progress flag. |


## Retrieve cross-species gene–phenotype associations {#dnakit.databases.monarch_get_gene_phenotypes-parameters}

`dnakit.databases.monarch_get_gene_phenotypes(*, subject, limit=None, offset=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `subject` | string | Required | Gene CURIE identifier (e.g., &#x27;HGNC:11998&#x27; for TP53, &#x27;HGNC:1100&#x27; for BRCA1). Use Monarch_search_gene to find this. |
| `limit` | integer | None | Maximum phenotype associations to return (default: 100) {&quot;default&quot;: 100} |
| `offset` | integer | None | Number of initial entries to skip. {&quot;default&quot;: 0} |
| `verbose` | bool | True | Progress flag. |


## Search Monarch biomedical entities by name {#dnakit.databases.monarch_search-parameters}

`dnakit.databases.monarch_search(*, query, category=None, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required | Search query: gene names, disease names, phenotype terms. Examples: &#x27;BRCA1&#x27;, &#x27;Alzheimer&#x27;, &#x27;seizure&#x27;, &#x27;insulin resistance&#x27;. |
| `category` | string | None | Filter by entity category. Options: &#x27;biolink:Gene&#x27;, &#x27;biolink:Disease&#x27;, &#x27;biolink:PhenotypicFeature&#x27;. Default: all categories. |
| `limit` | integer | None | Maximum results (default: 10, max: 50). |
| `verbose` | bool | True | Progress flag. |


## Search Monarch genes and cross-species identifiers {#dnakit.databases.monarch_search_gene-parameters}

`dnakit.databases.monarch_search_gene(*, query, limit=None, offset=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required | Gene symbol or name (e.g., &#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;WDR7&#x27;) |
| `limit` | integer | None | Number of results to return (default: 10) {&quot;default&quot;: 10} |
| `offset` | integer | None | Number of initial entries to skip. {&quot;default&quot;: 0} |
| `verbose` | bool | True | Progress flag. |


## Query cross-species biomedical associations {#dnakit.databases.monarch_v3_get_associations-parameters}

`dnakit.databases.monarch_v3_get_associations(*, subject=None, object=None, category, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `subject` | string | None | Subject-side entity CURIE, i.e. the entity on the LEFT of the association (the arrow&#x27;s source). For gene-&gt;disease categories (CausalGeneToDiseaseAssociation, CorrelatedGeneToDiseaseAssociation) the subject is the GENE; for gene-&gt;phenotype and gene-&gt;pathway the subject is the gene; for DiseaseToPhenotypicFeatureAssociation the subject is the disease. Examples: &#x27;HGNC:11998&#x27; (TP53), &#x27;MONDO:0005148&#x27; (type 2 diabetes), &#x27;HP:0001250&#x27; (seizure). Optional, but at least one of &#x27;subject&#x27; / &#x27;object&#x27; must be supplied. |
| `object` | string | None | Object-side entity CURIE, i.e. the entity on the RIGHT of the association (the arrow&#x27;s target). Use this to query the association backwards. For gene-&gt;disease categories the object is the DISEASE, so &#x27;which genes cause this disease?&#x27; is object=&#x27;MONDO:0006559&#x27; with category=&#x27;biolink:CausalGeneToDiseaseAssociation&#x27;; for gene/disease-&gt;phenotype categories the object is the HP phenotype term. Optional, but at least one of &#x27;subject&#x27; / &#x27;object&#x27; must be supplied. Supplying both narrows to associations linking exactly that pair. |
| `category` | string | Required | Biolink association category, written as subject-&gt;object. Options: &#x27;biolink:GeneToPhenotypicFeatureAssociation&#x27; (gene-&gt;phenotypes), &#x27;biolink:DiseaseToPhenotypicFeatureAssociation&#x27; (disease-&gt;phenotypes), &#x27;biolink:CorrelatedGeneToDiseaseAssociation&#x27; (gene-&gt;diseases, correlated), &#x27;biolink:CausalGeneToDiseaseAssociation&#x27; (gene-&gt;diseases, causal), &#x27;biolink:GeneToPathwayAssociation&#x27; (gene-&gt;pathways), &#x27;biolink:VariantToDiseaseAssociation&#x27; (variant-&gt;diseases), &#x27;biolink:GeneToExpressionSiteAssociation&#x27; (gene-&gt;expression sites/tissues). The left-hand side goes in &#x27;subject&#x27; and the right-hand side goes in &#x27;object&#x27;. |
| `limit` | integer | None | Maximum number of associations to return (default: 20, max: 200). |
| `verbose` | bool | True | Progress flag. |


## Retrieve biomedical entity details by CURIE {#dnakit.databases.monarch_v3_get_entity-parameters}

`dnakit.databases.monarch_v3_get_entity(*, entity_id, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `entity_id` | string | Required | Entity CURIE identifier. Examples: &#x27;HGNC:11998&#x27; (TP53), &#x27;HGNC:1100&#x27; (BRCA1), &#x27;MONDO:0005148&#x27; (type 2 diabetes), &#x27;HP:0001250&#x27; (seizure), &#x27;OMIM:191170&#x27; (Li-Fraumeni syndrome). |
| `verbose` | bool | True | Progress flag. |


## Summarize disease phenotypes by body system {#dnakit.databases.monarch_v3_get_histopheno-parameters}

`dnakit.databases.monarch_v3_get_histopheno(*, entity_id, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `entity_id` | string | Required | Disease CURIE identifier. Examples: &#x27;MONDO:0007947&#x27; (Marfan syndrome), &#x27;MONDO:0005148&#x27; (type 2 diabetes), &#x27;MONDO:0004975&#x27; (Alzheimer disease), &#x27;MONDO:0005015&#x27; (diabetes mellitus). |
| `verbose` | bool | True | Progress flag. |


## Retrieve cross-ontology identifier mappings {#dnakit.databases.monarch_v3_get_mappings-parameters}

`dnakit.databases.monarch_v3_get_mappings(*, entity_id, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `entity_id` | string | Required | Entity CURIE identifier to find mappings for. Examples: &#x27;MONDO:0005148&#x27; (type 2 diabetes), &#x27;MONDO:0007947&#x27; (Marfan syndrome), &#x27;HP:0001250&#x27; (seizure). |
| `limit` | integer | None | Maximum number of mappings to return (default: 20, max: 100). |
| `verbose` | bool | True | Progress flag. |


## Find diseases by phenotype semantic similarity {#dnakit.databases.monarch_v3_phenotype_similarity_search-parameters}

`dnakit.databases.monarch_v3_phenotype_similarity_search(*, phenotypes, group=None, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `phenotypes` | array | Required | List of HPO term CURIEs representing observed phenotypes. Examples: &#91;&#x27;HP:0001250&#x27;, &#x27;HP:0001249&#x27;, &#x27;HP:0002353&#x27;&#93; for seizure + intellectual disability + EEG abnormality. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `group` | string | None | Entity group to search against. Default: &#x27;Human Diseases&#x27;. Options: &#x27;Human Diseases&#x27;, &#x27;Mouse&#x27;, &#x27;Rat&#x27;, &#x27;Zebrafish&#x27;, &#x27;C. elegans&#x27;, &#x27;Drosophila&#x27;. |
| `limit` | integer | None | Maximum number of matching diseases to return (default: 10, max: 50). |
| `verbose` | bool | True | Progress flag. |


## Search Monarch knowledge graph entities and categories {#dnakit.databases.monarch_v3_search-parameters}

`dnakit.databases.monarch_v3_search(*, query=None, limit=None, category=None, q=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | None | Search query. Can be gene names, disease names, phenotype terms, etc. Examples: &#x27;BRCA1&#x27;, &#x27;breast cancer&#x27;, &#x27;seizure&#x27;, &#x27;insulin resistance&#x27;. |
| `limit` | integer | None | Maximum results to return (default: 10, max: 50). |
| `category` | string | None | Filter by entity category. Options: &#x27;biolink:Gene&#x27;, &#x27;biolink:Disease&#x27;, &#x27;biolink:PhenotypicFeature&#x27;, &#x27;biolink:NamedThing&#x27;. Default: all categories. |
| `q` | string | None | Alias for query parameter. |
| `verbose` | bool | True | Progress flag. |


## Retrieve ontology term ancestors {#dnakit.databases.ols_get_term_ancestors-parameters}

`dnakit.databases.ols_get_term_ancestors(*, operation=None, term_iri=None, term_id=None, ontology=None, include_obsolete=None, size=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | The operation to perform (get_term_ancestors) {&quot;const&quot;: &quot;get_term_ancestors&quot;, &quot;default&quot;: &quot;get_term_ancestors&quot;} |
| `term_iri` | string | None | The full IRI of the term (e.g., http://www.ebi.ac.uk/efo/EFO_0000408) |
| `term_id` | string | None | Short-form CURIE (e.g., &#x27;HP:0001928&#x27;, &#x27;GO:0008150&#x27;). Ontology is auto-inferred from prefix. |
| `ontology` | string | None | Ontology ID (e.g., &#x27;hp&#x27;, &#x27;go&#x27;, &#x27;efo&#x27;). Auto-inferred from term_id prefix when not supplied. |
| `include_obsolete` | boolean | None | Include obsolete terms (default: false) {&quot;default&quot;: false} |
| `size` | integer | None | Number of results to return (default: 20) {&quot;default&quot;: 20} |
| `verbose` | bool | True | Progress flag. |


## Retrieve ontology term children {#dnakit.databases.ols_get_term_children-parameters}

`dnakit.databases.ols_get_term_children(*, operation=None, term_iri=None, term_id=None, ontology=None, include_obsolete=None, size=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | The operation to perform (get_term_children) {&quot;const&quot;: &quot;get_term_children&quot;, &quot;default&quot;: &quot;get_term_children&quot;} |
| `term_iri` | string | None | The full IRI of the term (e.g., http://www.ebi.ac.uk/efo/EFO_0000408) |
| `term_id` | string | None | Short-form CURIE (e.g., &#x27;HP:0001928&#x27;, &#x27;GO:0008150&#x27;). Ontology is auto-inferred from prefix. |
| `ontology` | string | None | Ontology ID (e.g., &#x27;hp&#x27;, &#x27;go&#x27;, &#x27;efo&#x27;). Auto-inferred from term_id prefix when not supplied. |
| `include_obsolete` | boolean | None | Include obsolete terms (default: false) {&quot;default&quot;: false} |
| `size` | integer | None | Number of results to return (default: 20) {&quot;default&quot;: 20} |
| `verbose` | bool | True | Progress flag. |


## Retrieve ontology term details {#dnakit.databases.ols_get_term_info-parameters}

`dnakit.databases.ols_get_term_info(*, operation=None, id=None, term_id=None, term_iri=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | The operation to perform (get_term_info) {&quot;const&quot;: &quot;get_term_info&quot;, &quot;default&quot;: &quot;get_term_info&quot;} |
| `id` | string | None | Ontology term ID or IRI (e.g. HP:0001903, http://purl.obolibrary.org/obo/HP_0001903) |
| `term_id` | string | None | Alias for id — ontology term ID (e.g. HP:0001903) |
| `term_iri` | string | None | Alias for id. Full IRI or short-form ID (e.g., GO:0008150, EFO:0003971). Same as id. |
| `verbose` | bool | True | Progress flag. |


## Search ontology terms by keyword {#dnakit.databases.ols_search_terms-parameters}

`dnakit.databases.ols_search_terms(*, operation=None, query, rows=None, size=None, ontology=None, exact_match=None, include_obsolete=None, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | The operation to perform (search_terms) {&quot;const&quot;: &quot;search_terms&quot;, &quot;default&quot;: &quot;search_terms&quot;} |
| `query` | string | Required | The search query for terms |
| `rows` | integer | None | Number of results to return (default: 10). Alias: size. {&quot;default&quot;: 10} |
| `size` | integer | None | Alias for rows. Number of results to return (default: 10). {&quot;default&quot;: 10} |
| `ontology` | string | None | Filter by specific ontology (optional) |
| `exact_match` | boolean | None | Restrict results to exact matches (default: false). The query must equal a term&#x27;s label or one of its exact synonyms in full - substring and word-overlap hits are excluded, so &#x27;fibroblast&#x27; returns CL:0000057 but not &#x27;skin fibroblast&#x27;. If the query is instead an identifier (a CURIE such as &#x27;CL:0000084&#x27;, the OBO underscore form &#x27;CL_0000084&#x27;, or a full IRI), the exact match is applied to the obo_id/short_form/iri fields, which returns that identifier in every ontology that declares it - pass `ontology` to narrow it to one. {&quot;default&quot;: false} |
| `include_obsolete` | boolean | None | Include obsolete terms (default: false) {&quot;default&quot;: false} |
| `limit` | integer &#124; null | None | Alias for rows. Maximum number of results to return (default: 10, max: 100). |
| `verbose` | bool | True | Progress flag. |


## Retrieve an OMIM clinical phenotype synopsis {#dnakit.databases.omim_get_clinical_synopsis-parameters}

`dnakit.databases.omim_get_clinical_synopsis(*, operation=None, mim_number, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_clinical_synopsis) {&quot;const&quot;: &quot;get_clinical_synopsis&quot;, &quot;default&quot;: &quot;get_clinical_synopsis&quot;} |
| `mim_number` | string | Required | OMIM MIM number for phenotype (e.g., 219700 for cystic fibrosis) |
| `verbose` | bool | True | Progress flag. |

Configure before running: `OMIM_API_KEY`.


## Retrieve an OMIM entry by MIM number {#dnakit.databases.omim_get_entry-parameters}

`dnakit.databases.omim_get_entry(*, operation=None, mim_number, include=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_entry) {&quot;const&quot;: &quot;get_entry&quot;, &quot;default&quot;: &quot;get_entry&quot;} |
| `mim_number` | string | Required | OMIM MIM number (e.g., 164730, 219700). Can include OMIM: prefix. |
| `include` | string | None | Data sections to include (default: text,clinicalSynopsis,geneMap) {&quot;default&quot;: &quot;text,clinicalSynopsis,geneMap&quot;} |
| `verbose` | bool | True | Progress flag. |

Configure before running: `OMIM_API_KEY`.


## Retrieve OMIM gene–disease and chromosome mappings {#dnakit.databases.omim_get_gene_map-parameters}

`dnakit.databases.omim_get_gene_map(*, operation=None, mim_number=None, chromosome=None, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_gene_map) {&quot;const&quot;: &quot;get_gene_map&quot;, &quot;default&quot;: &quot;get_gene_map&quot;} |
| `mim_number` | string | None | OMIM MIM number for specific gene (optional if chromosome provided) |
| `chromosome` | string | None | Chromosome number (1-22, X, Y) to get all genes (optional if mim_number provided) |
| `limit` | integer | None | Maximum results for chromosome query (default: 50, max: 100) {&quot;default&quot;: 50} |
| `verbose` | bool | True | Progress flag. |

Configure before running: `OMIM_API_KEY`.


## Search Mendelian disorders and associated genes {#dnakit.databases.omim_search-parameters}

`dnakit.databases.omim_search(*, operation=None, query, limit=None, start=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search) {&quot;const&quot;: &quot;search&quot;, &quot;default&quot;: &quot;search&quot;} |
| `query` | string | Required | Search query - gene name (BRCA1), disease (cystic fibrosis), or phenotype (intellectual disability) |
| `limit` | integer | None | Maximum results (default: 10, max: 20) {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 20} |
| `start` | integer | None | Offset for pagination (default: 0) {&quot;default&quot;: 0} |
| `verbose` | bool | True | Progress flag. |

Configure before running: `OMIM_API_KEY`.


## Retrieve a rare disease classification hierarchy {#dnakit.databases.orphanet_get_classification-parameters}

`dnakit.databases.orphanet_get_classification(*, operation=None, orpha_code=None, lang=None, orphacode=None, orpha_id=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_classification) {&quot;const&quot;: &quot;get_classification&quot;, &quot;default&quot;: &quot;get_classification&quot;} |
| `orpha_code` | string &#124; integer | None | Orphanet ORPHA code |
| `lang` | string | None | Language code (default: en) {&quot;default&quot;: &quot;en&quot;} |
| `orphacode` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `orpha_id` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `verbose` | bool | True | Progress flag. |


## Retrieve rare disease details by ORPHA code {#dnakit.databases.orphanet_get_disease-parameters}

`dnakit.databases.orphanet_get_disease(*, operation=None, orpha_code=None, lang=None, orphacode=None, orpha_id=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_disease) {&quot;const&quot;: &quot;get_disease&quot;, &quot;default&quot;: &quot;get_disease&quot;} |
| `orpha_code` | string &#124; integer | None | Orphanet ORPHA code (e.g., 558, 166024). Can include ORPHA: prefix. |
| `lang` | string | None | Language code (default: en) {&quot;default&quot;: &quot;en&quot;} |
| `orphacode` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `orpha_id` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `verbose` | bool | True | Progress flag. |


## Retrieve reported gene associations for an ORPHA disorder {#dnakit.databases.orphanet_get_disease_genes-parameters}

`dnakit.databases.orphanet_get_disease_genes(*, orpha_code, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `orpha_code` | string | Required | Orphanet ORPHA disorder identifier. |
| `verbose` | bool | True | Progress flag. |


## Retrieve rare disease prevalence and incidence {#dnakit.databases.orphanet_get_epidemiology-parameters}

`dnakit.databases.orphanet_get_epidemiology(*, operation=None, orpha_code=None, orphacode=None, orpha_id=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_epidemiology) {&quot;const&quot;: &quot;get_epidemiology&quot;, &quot;default&quot;: &quot;get_epidemiology&quot;} |
| `orpha_code` | string &#124; integer | None | Orphanet ORPHA code (e.g., 558 for Marfan, 586 for Cystic Fibrosis) |
| `orphacode` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `orpha_id` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `verbose` | bool | True | Progress flag. |


## Find rare diseases associated with a gene {#dnakit.databases.orphanet_get_gene_diseases-parameters}

`dnakit.databases.orphanet_get_gene_diseases(*, operation=None, gene_name=None, gene_symbol=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_gene_diseases) {&quot;const&quot;: &quot;get_gene_diseases&quot;, &quot;default&quot;: &quot;get_gene_diseases&quot;} |
| `gene_name` | string | None | Gene symbol or name keyword (e.g., &#x27;FBN1&#x27;, &#x27;BRCA1&#x27;, &#x27;huntingtin&#x27;). Use gene_symbol for standard HGNC symbols. Note: Orphanet uses current HGNC symbols (e.g., &#x27;GBA1&#x27; not &#x27;GBA&#x27;). |
| `gene_symbol` | string &#124; null | None | Alternative to gene_name: standard HGNC gene symbol (e.g., &#x27;GBA1&#x27;, &#x27;FBN1&#x27;). Note: Orphanet uses current HGNC symbols (e.g., &#x27;GBA1&#x27; not &#x27;GBA&#x27;). |
| `verbose` | bool | True | Progress flag. |


## Retrieve rare disease genes and association types {#dnakit.databases.orphanet_get_genes-parameters}

`dnakit.databases.orphanet_get_genes(*, operation=None, orpha_code=None, orphacode=None, orpha_id=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_genes) {&quot;const&quot;: &quot;get_genes&quot;, &quot;default&quot;: &quot;get_genes&quot;} |
| `orpha_code` | string &#124; integer | None | Orphanet ORPHA code for the disease |
| `orphacode` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `orpha_id` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `verbose` | bool | True | Progress flag. |


## Map rare disease codes to other medical coding systems {#dnakit.databases.orphanet_get_icd_mapping-parameters}

`dnakit.databases.orphanet_get_icd_mapping(*, operation=None, orpha_code=None, coding_system=None, lang=None, orphacode=None, orpha_id=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_icd_mapping) {&quot;const&quot;: &quot;get_icd_mapping&quot;, &quot;default&quot;: &quot;get_icd_mapping&quot;} |
| `orpha_code` | string &#124; integer | None | Orphanet ORPHA code (e.g., 558 for Marfan, 586 for Cystic Fibrosis) |
| `coding_system` | string | None | Which coding system to retrieve: all, icd10, icd11, omim, snomed. Default: all {&quot;default&quot;: &quot;all&quot;} |
| `lang` | string | None | Language code (default: en) {&quot;default&quot;: &quot;en&quot;} |
| `orphacode` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `orpha_id` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `verbose` | bool | True | Progress flag. |


## Retrieve rare disease onset and inheritance information {#dnakit.databases.orphanet_get_natural_history-parameters}

`dnakit.databases.orphanet_get_natural_history(*, operation=None, orpha_code=None, orphacode=None, orpha_id=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_natural_history) {&quot;const&quot;: &quot;get_natural_history&quot;, &quot;default&quot;: &quot;get_natural_history&quot;} |
| `orpha_code` | string &#124; integer | None | Orphanet ORPHA code (e.g., 558 for Marfan, 399 for Huntington) |
| `orphacode` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `orpha_id` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `verbose` | bool | True | Progress flag. |


## Retrieve rare disease phenotypes and frequencies {#dnakit.databases.orphanet_get_phenotypes-parameters}

`dnakit.databases.orphanet_get_phenotypes(*, operation=None, orpha_code=None, orpha_id=None, orphacode=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_phenotypes) {&quot;const&quot;: &quot;get_phenotypes&quot;, &quot;default&quot;: &quot;get_phenotypes&quot;} |
| `orpha_code` | string &#124; integer | None | Orphanet ORPHA code (e.g., 558 for Marfan, 586 for Cystic Fibrosis) |
| `orpha_id` | string &#124; integer | None | Alias for orpha_code. Orphanet ORPHA code (e.g., 558 for Marfan) |
| `orphacode` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `verbose` | bool | True | Progress flag. |


## Search rare diseases by complete or partial name {#dnakit.databases.orphanet_search_by_name-parameters}

`dnakit.databases.orphanet_search_by_name(*, operation=None, name=None, exact=None, lang=None, limit=None, query=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search_by_name) {&quot;const&quot;: &quot;search_by_name&quot;, &quot;default&quot;: &quot;search_by_name&quot;} |
| `name` | string | None | Disease name to search (e.g., &#x27;Marfan syndrome&#x27;) |
| `exact` | boolean | None | If true, match exact name only. Default: false {&quot;default&quot;: false} |
| `lang` | string | None | Language code (default: en) {&quot;default&quot;: &quot;en&quot;} |
| `limit` | integer | None | Maximum number of results to return (default: 20, max: 200) |
| `query` | string | None | Alias for name. Disease name to search for in Orphanet. |
| `verbose` | bool | True | Progress flag. |


## Search rare diseases by keyword {#dnakit.databases.orphanet_search_diseases-parameters}

`dnakit.databases.orphanet_search_diseases(*, operation=None, query, limit=None, lang=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search_diseases) {&quot;const&quot;: &quot;search_diseases&quot;, &quot;default&quot;: &quot;search_diseases&quot;} |
| `query` | string | Required | Search query - disease name or keyword (e.g., &#x27;Marfan&#x27;, &#x27;muscular dystrophy&#x27;) |
| `limit` | integer | None | Maximum number of results to return (default: 20, max: 200). The API may return thousands of fuzzy matches; use this to control output size. {&quot;default&quot;: 20} |
| `lang` | string | None | Language code (en, fr, de, es, it, pt, pl, nl). Default: en {&quot;default&quot;: &quot;en&quot;} |
| `verbose` | bool | True | Progress flag. |


## Resolve a disease name to an EFO term {#dnakit.databases.osl_get_efo_id_by_disease_name-parameters}

`dnakit.databases.osl_get_efo_id_by_disease_name(*, disease, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `disease` | string | Required | Disease name or free-text query (e.g., &#x27;diabetes mellitus&#x27;). |
| `verbose` | bool | True | Progress flag. |
