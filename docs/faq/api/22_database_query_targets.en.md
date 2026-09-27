# Target, drug and association queries: parameters and notes

See the [tutorial](../../api/features/22_database_query_targets.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

## Retrieve drug mechanisms by name or identifier {#dnakit.databases.ch_embl_get_drug_mechanisms-parameters}

`dnakit.databases.ch_embl_get_drug_mechanisms(*, limit=None, offset=None, drug_chembl_id=None, drug_name=None, molecule_chembl_id=None, chembl_id=None, drug_chembl_id__exact=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `limit` | integer | None | Maximum results for the current request, not a database-wide total. {&quot;default&quot;: 20, &quot;maximum&quot;: 1000} |
| `offset` | integer | None | Pagination offset within the current query. {&quot;default&quot;: 0} |
| `drug_chembl_id` | string | None | ChEMBL drug/molecule ID (e.g., &quot;CHEMBL1201580&quot; for adalimumab, &quot;CHEMBL4535757&quot; for sotorasib). |
| `drug_name` | string | None | Drug name for automatic ChEMBL ID lookup (e.g., &quot;trastuzumab&quot;, &quot;lapatinib&quot;, &quot;aspirin&quot;). Case-insensitive. Triggers internal ChEMBL molecule search — use drug_chembl_id if you already know the ID. |
| `molecule_chembl_id` | string | None | Alias for drug_chembl_id. ChEMBL molecule ID (e.g., &quot;CHEMBL25&quot; for aspirin). |
| `chembl_id` | string | None | Alias for drug_chembl_id. ChEMBL ID (e.g., &quot;CHEMBL3301622&quot;). |
| `drug_chembl_id__exact` | string | None | DNAKit optional alias for drug_chembl_id |
| `verbose` | bool | True | Progress flag. |


## Retrieve ChEMBL molecular structures and properties {#dnakit.databases.ch_embl_get_molecule-parameters}

`dnakit.databases.ch_embl_get_molecule(*, chembl_id=None, format=None, molecule_chembl_id=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `chembl_id` | string | None | ChEMBL molecule ID, e.g., &#x27;CHEMBL25&#x27; |
| `format` | string | None | Response format {&quot;enum&quot;: &#91;&quot;json&quot;, &quot;xml&quot;, &quot;yaml&quot;&#93;, &quot;default&quot;: &quot;json&quot;} |
| `molecule_chembl_id` | string | None | Alias for chembl_id. ChEMBL molecule ID, e.g., CHEMBL1229517. |
| `verbose` | bool | True | Progress flag. |


## Retrieve bioactivity measurements for a target {#dnakit.databases.ch_embl_get_target_activities-parameters}

`dnakit.databases.ch_embl_get_target_activities(*, target_chembl_id__exact=None, limit=None, offset=None, target_chembl_id=None, chembl_id=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `target_chembl_id__exact` | string | None | ChEMBL target ID (e.g., &#x27;CHEMBL2074&#x27;). To find a target ID, use ChEMBL_search_targets with a target name or gene symbol. |
| `limit` | integer | None | Maximum results for the current request, not a database-wide total. {&quot;default&quot;: 20, &quot;maximum&quot;: 1000} |
| `offset` | integer | None | Pagination offset within the current query. {&quot;default&quot;: 0} |
| `target_chembl_id` | string | None | Alias for target_chembl_id__exact. ChEMBL target ID (e.g., CHEMBL213). |
| `chembl_id` | string | None | DNAKit optional alias for target_chembl_id |
| `verbose` | bool | True | Progress flag. |

Provide exactly one of these parameter groups: `target_chembl_id__exact`; `target_chembl_id`.


## Search drugs by name and development phase {#dnakit.databases.ch_embl_search_drugs-parameters}

`dnakit.databases.ch_embl_search_drugs(*, query=None, max_phase=None, limit=None, offset=None, molecule_chembl_id=None, q=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | None | Drug name to search for (partial match, case-insensitive). E.g., &quot;sotorasib&quot;, &quot;olaparib&quot;, &quot;imatinib&quot;. |
| `max_phase` | integer | None | Filter by maximum development phase (0-4) |
| `limit` | integer | None | Maximum results for the current request, not a database-wide total. {&quot;default&quot;: 20, &quot;maximum&quot;: 1000} |
| `offset` | integer | None | Pagination offset within the current query. {&quot;default&quot;: 0} |
| `molecule_chembl_id` | string | None | Filter by ChEMBL molecule ID (e.g., &quot;CHEMBL1201580&quot; for adalimumab). |
| `q` | string | None | DNAKit optional alias for query |
| `verbose` | bool | True | Progress flag. |


## Search drug targets by name, organism or type {#dnakit.databases.ch_embl_search_targets-parameters}

`dnakit.databases.ch_embl_search_targets(*, target_chembl_id=None, pref_name__contains=None, organism=None, target_type=None, fields=None, limit=None, offset=None, q=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `target_chembl_id` | string | None | Filter by target ChEMBL ID |
| `pref_name__contains` | string | None | Filter by target name (contains) |
| `organism` | string | None | Filter by organism (e.g., &#x27;Homo sapiens&#x27;) |
| `target_type` | string | None | Filter by target type (e.g., &#x27;SINGLE PROTEIN&#x27;, &#x27;PROTEIN COMPLEX&#x27;) |
| `fields` | array | None | Optional list of ChEMBL target fields to include in each returned target object (projection). ToolUniverse maps this to ChEMBL&#x27;s `only` query parameter (comma-separated). Supported fields: target_chembl_id, pref_name, organism, target_type, target_components. {&quot;uniqueItems&quot;: true, &quot;items&quot;: {&quot;type&quot;: &quot;string&quot;, &quot;enum&quot;: &#91;&quot;target_chembl_id&quot;, &quot;pref_name&quot;, &quot;organism&quot;, &quot;target_type&quot;, &quot;target_components&quot;&#93;}} |
| `limit` | integer | None | Maximum results for the current request, not a database-wide total. {&quot;default&quot;: 20, &quot;maximum&quot;: 1000} |
| `offset` | integer | None | Pagination offset within the current query. {&quot;default&quot;: 0} |
| `q` | string | None | DNAKit optional alias for pref_name__contains |
| `verbose` | bool | True | Progress flag. |


## Inspect the unavailable CTD chemical–gene interaction query {#dnakit.databases.ctd_get_chemical_gene_interactions-parameters}

`dnakit.databases.ctd_get_chemical_gene_interactions(*, input_terms, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `input_terms` | string | Required | Chemical name, CAS RN, or MeSH ID. Unused -- this tool always returns an error; see the tool description. |
| `verbose` | bool | True | Progress flag. |


## Retrieve drug–gene interaction evidence in batches {#dnakit.databases.dgidb_get_drug_gene_interactions-parameters}

`dnakit.databases.dgidb_get_drug_gene_interactions(*, genes=None, gene_name=None, gene=None, interaction_sources=None, interaction_types=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `genes` | array | None | List of gene symbols (e.g., &#91;&#x27;EGFR&#x27;, &#x27;BRAF&#x27;&#93;). Also accepts a single gene as string. Aliases: gene_name, gene. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `gene_name` | string | None | Alias for genes. Single gene symbol (e.g., &#x27;EGFR&#x27;). |
| `gene` | string | None | Alias for genes. Single gene symbol (e.g., &#x27;EGFR&#x27;). |
| `interaction_sources` | array | None | Optional filter by data sources (e.g., &#91;&#x27;DrugBank&#x27;, &#x27;ChEMBL&#x27;&#93;). {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `interaction_types` | array | None | Optional filter by interaction types (e.g., &#91;&#x27;inhibitor&#x27;, &#x27;antagonist&#x27;&#93;). {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `verbose` | bool | True | Progress flag. |


## Retrieve gene druggability categories {#dnakit.databases.dgidb_get_gene_druggability-parameters}

`dnakit.databases.dgidb_get_gene_druggability(*, genes=None, gene_name=None, gene=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `genes` | array | None | List of gene symbols to check druggability. Aliases: gene_name, gene. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `gene_name` | string | None | Alias for genes. Single gene symbol (e.g., &#x27;EGFR&#x27;). |
| `gene` | string | None | Alias for genes. Single gene symbol (e.g., &#x27;EGFR&#x27;). |
| `verbose` | bool | True | Progress flag. |


## Retrieve basic drug information and approval status {#dnakit.databases.drugbank_get_drug_basic_info_by_drug_name_or_id-parameters}

`dnakit.databases.drugbank_get_drug_basic_info_by_drug_name_or_id(*, query, case_sensitive=None, exact_match=None, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required | Drug name or DrugBank ID to search for {&quot;minLength&quot;: 1} |
| `case_sensitive` | boolean | None | Select True to perform a case-sensitive search {&quot;default&quot;: false} |
| `exact_match` | boolean | None | Select True to require an exact match with the queried name or ID {&quot;default&quot;: false} |
| `limit` | integer | None | Maximum number of results to return {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | Progress flag. |


## Retrieve pharmacodynamics, mechanisms and pharmacokinetics {#dnakit.databases.drugbank_get_pharmacology_by_drug_name_or_drugbank_id-parameters}

`dnakit.databases.drugbank_get_pharmacology_by_drug_name_or_drugbank_id(*, query, case_sensitive=None, exact_match=None, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required | Drug name or Drugbank ID to search for {&quot;minLength&quot;: 1} |
| `case_sensitive` | boolean | None | Select True to perform a case-sensitive search {&quot;default&quot;: false} |
| `exact_match` | boolean | None | Select True to require an exact match {&quot;default&quot;: false} |
| `limit` | integer | None | Maximum number of results to return {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | Progress flag. |


## Retrieve drug targets, enzymes, carriers and transporters {#dnakit.databases.drugbank_get_targets_by_drug_name_or_drugbank_id-parameters}

`dnakit.databases.drugbank_get_targets_by_drug_name_or_drugbank_id(*, query, case_sensitive=None, exact_match=None, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required | Drug name or Drugbank ID to search for {&quot;minLength&quot;: 1} |
| `case_sensitive` | boolean | None | Select True to perform a case-sensitive search {&quot;default&quot;: false} |
| `exact_match` | boolean | None | Select True to require an exact match {&quot;default&quot;: false} |
| `limit` | integer | None | Maximum number of results to return {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | Progress flag. |


## Search the DrugBank vocabulary by name or synonym {#dnakit.databases.drugbank_vocab_search-parameters}

`dnakit.databases.drugbank_vocab_search(*, query, search_fields=None, case_sensitive=None, exact_match=None, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required | Search query string. Can be drug name, synonym, DrugBank ID, or any text to search for. {&quot;minLength&quot;: 1} |
| `search_fields` | array | None | Fields to search in. Available fields: &#x27;DrugBank ID&#x27;, &#x27;Accession Numbers&#x27;, &#x27;Common name&#x27;, &#x27;CAS&#x27;, &#x27;UNII&#x27;, &#x27;Synonyms&#x27;, &#x27;Standard InChI Key&#x27;. Defaults to &#91;&#x27;Common name&#x27;, &#x27;Synonyms&#x27;, &#x27;DrugBank ID&#x27;&#93;. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `case_sensitive` | boolean | None | Whether the search should be case sensitive. {&quot;default&quot;: false} |
| `exact_match` | boolean | None | Whether to perform exact matching instead of substring matching. {&quot;default&quot;: false} |
| `limit` | integer | None | Maximum number of results to return. {&quot;minimum&quot;: 1, &quot;maximum&quot;: 1000, &quot;default&quot;: 10} |
| `verbose` | bool | True | Progress flag. |


## Retrieve literature-based disease–target association scores {#dnakit.databases.europepmc_disease_target_score-parameters}

`dnakit.databases.europepmc_disease_target_score(*, efoId, pageSize=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | Required | The EFO (Experimental Factor Ontology) ID of the disease, e.g., &#x27;MONDO_0011996&#x27; for chronic myelogenous leukemia |
| `pageSize` | integer | None | Upstream fetch batch size (max 100), NOT a limit on results returned. The tool pages through all of the disease&#x27;s associated targets regardless of this value -- &#x27;pageSize&#x27;: 10 on MONDO_0005011 (Crohn disease) still returns 390 scored targets. Lower it only to make each upstream request smaller. {&quot;default&quot;: 100} |
| `verbose` | bool | True | Progress flag. |


## Retrieve expression-based disease–target association scores {#dnakit.databases.expression_atlas_disease_target_score-parameters}

`dnakit.databases.expression_atlas_disease_target_score(*, efoId, pageSize=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | Required | The EFO (Experimental Factor Ontology) ID of the disease, e.g., &#x27;MONDO_0011996&#x27; for chronic myelogenous leukemia |
| `pageSize` | integer | None | Number of results per page (default: 100, max: 100) {&quot;default&quot;: 100} |
| `verbose` | bool | True | Progress flag. |


## Retrieve active ingredients in a drug product {#dnakit.databases.fda_get_active_ingredient_info_by_drug_name-parameters}

`dnakit.databases.fda_get_active_ingredient_info_by_drug_name(*, drug_name, limit=None, skip=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `drug_name` | string | Required | The name of the drug. |
| `limit` | integer | None | The number of records to return. |
| `skip` | integer | None | The number of records to skip. |
| `verbose` | bool | True | Progress flag. |


## Retrieve adverse reactions from drug labels {#dnakit.databases.fda_get_adverse_reactions_by_drug_name-parameters}

`dnakit.databases.fda_get_adverse_reactions_by_drug_name(*, drug_name, limit=None, skip=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `drug_name` | string | Required | The name of the drug. |
| `limit` | integer | None | The number of records to return. |
| `skip` | integer | None | The number of records to skip. |
| `verbose` | bool | True | Progress flag. |


## Retrieve drug boxed warnings and safety information {#dnakit.databases.fda_get_boxed_warning_info_by_drug_name-parameters}

`dnakit.databases.fda_get_boxed_warning_info_by_drug_name(*, drug_name, limit=None, skip=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `drug_name` | string | Required | The name of the drug. |
| `limit` | integer | None | The number of records to return. |
| `skip` | integer | None | The number of records to skip. |
| `verbose` | bool | True | Progress flag. |


## Retrieve clinical study information from drug labels {#dnakit.databases.fda_get_clinical_studies_info_by_drug_name-parameters}

`dnakit.databases.fda_get_clinical_studies_info_by_drug_name(*, drug_name, limit=None, skip=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `drug_name` | string | Required | The name of the drug. |
| `limit` | integer | None | The number of records to return. |
| `skip` | integer | None | The number of records to skip. |
| `verbose` | bool | True | Progress flag. |


## Retrieve drug indications and usage information {#dnakit.databases.fda_get_indications_by_drug_name-parameters}

`dnakit.databases.fda_get_indications_by_drug_name(*, drug_name, limit=None, skip=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `drug_name` | string | Required | The name of the drug. |
| `limit` | integer | None | The number of records to return. |
| `skip` | integer | None | The number of records to skip. |
| `verbose` | bool | True | Progress flag. |


## Retrieve mechanisms of action from drug labels {#dnakit.databases.fda_get_mechanism_of_action_by_drug_name-parameters}

`dnakit.databases.fda_get_mechanism_of_action_by_drug_name(*, drug_name, limit=None, skip=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `drug_name` | string | Required | The name of the drug. |
| `limit` | integer | None | The number of records to return. |
| `skip` | integer | None | The number of records to skip. |
| `verbose` | bool | True | Progress flag. |


## Retrieve metabolite–disease associations {#dnakit.databases.hmdb_get_diseases-parameters}

`dnakit.databases.hmdb_get_diseases(*, operation=None, hmdb_id=None, compound_name=None, pubchem_cid=None, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Provider operation mode; see this row&#x27;s default and allowed values. {&quot;const&quot;: &quot;get_diseases&quot;, &quot;default&quot;: &quot;get_diseases&quot;} |
| `hmdb_id` | string | None | HMDB ID (e.g., HMDB0000122) |
| `compound_name` | string | None | Compound name (e.g., glucose) |
| `pubchem_cid` | string | None | PubChem CID as an alternative identifier to hmdb_id or compound_name. |
| `limit` | integer | None | Maximum number of disease associations to return (default 50). {&quot;default&quot;: 50} |
| `verbose` | bool | True | Progress flag. |


## Retrieve metabolite information by identifier or name {#dnakit.databases.hmdb_get_metabolite-parameters}

`dnakit.databases.hmdb_get_metabolite(*, operation=None, hmdb_id=None, compound_name=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Provider operation mode; see this row&#x27;s default and allowed values. {&quot;const&quot;: &quot;get_info&quot;, &quot;default&quot;: &quot;get_info&quot;} |
| `hmdb_id` | string | None | HMDB ID (e.g., HMDB0000122) |
| `compound_name` | string | None | Compound name (e.g., glucose) |
| `verbose` | bool | True | Progress flag. |


## Search metabolites by name or molecular formula {#dnakit.databases.hmdb_search-parameters}

`dnakit.databases.hmdb_search(*, operation=None, query, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Provider operation mode; see this row&#x27;s default and allowed values. {&quot;const&quot;: &quot;search&quot;, &quot;default&quot;: &quot;search&quot;} |
| `query` | string | Required | Compound name or molecular formula to search |
| `verbose` | bool | True | Progress flag. |


## Retrieve known drugs associated with a disease {#dnakit.databases.open_targets_get_associated_drugs_by_disease_efo_id-parameters}

`dnakit.databases.open_targets_get_associated_drugs_by_disease_efo_id(*, efoId, size=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | Required | The EFO ID of the disease. |
| `size` | integer | None | Limit for the returned page; does not automatically retrieve later pages. {&quot;minimum&quot;: 1, &quot;maximum&quot;: 1000} |
| `verbose` | bool | True | Progress flag. |


## Retrieve target-associated drugs and development phases {#dnakit.databases.open_targets_get_associated_drugs_by_target_ensembl_id-parameters}

`dnakit.databases.open_targets_get_associated_drugs_by_target_ensembl_id(*, ensemblId, size=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `ensemblId` | string | Required | The Ensembl ID of the target. |
| `size` | integer | None | Limit for the returned page; does not automatically retrieve later pages. {&quot;minimum&quot;: 1, &quot;maximum&quot;: 1000} |
| `verbose` | bool | True | Progress flag. |


## Retrieve HPO phenotypes associated with a disease {#dnakit.databases.open_targets_get_associated_phenotypes_by_disease_efo_id-parameters}

`dnakit.databases.open_targets_get_associated_phenotypes_by_disease_efo_id(*, efoId, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | Required | The efoId of a disease or phenotype. |
| `verbose` | bool | True | Progress flag. |


## Retrieve disease-associated targets with pagination {#dnakit.databases.open_targets_get_associated_targets_by_disease_efo_id-parameters}

`dnakit.databases.open_targets_get_associated_targets_by_disease_efo_id(*, efoId, size=None, index=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | Required | The efoId of a disease or phenotype. |
| `size` | integer | None | Number of top-scored associated targets to return (default 50). A disease can have thousands of associated targets; raise this to retrieve more than the first page. {&quot;default&quot;: 50} |
| `index` | integer | None | Zero-based page index for paginating beyond the first `size` targets (default 0). |
| `verbose` | bool | True | Progress flag. |


## Retrieve GWAS credible set and locus-to-gene details {#dnakit.databases.open_targets_get_credible_set_detail-parameters}

`dnakit.databases.open_targets_get_credible_set_detail(*, studyLocusId, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `studyLocusId` | string | Required | Credible set identifier (studyLocusId), a 32-character hash (e.g., &#x27;b758d8fb10924f5338cbad8d27c7dee8&#x27;). Obtain from variant credible sets or study credible sets queries. |
| `verbose` | bool | True | Progress flag. |


## Retrieve disease ontology parents and ancestors {#dnakit.databases.open_targets_get_disease_ancestors_parents_by_efo_id-parameters}

`dnakit.databases.open_targets_get_disease_ancestors_parents_by_efo_id(*, efoId, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | Required | The EFO ID of the disease. |
| `verbose` | bool | True | Progress flag. |


## Retrieve disease ontology children and descendants {#dnakit.databases.open_targets_get_disease_descendants_children_by_efo_id-parameters}

`dnakit.databases.open_targets_get_disease_descendants_children_by_efo_id(*, efoId, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | Required | The EFO ID of the disease. |
| `verbose` | bool | True | Progress flag. |


## Retrieve disease descriptions and cross-references {#dnakit.databases.open_targets_get_disease_description_by_efo_id-parameters}

`dnakit.databases.open_targets_get_disease_description_by_efo_id(*, efoId, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | Required | The EFO ID of the disease. |
| `verbose` | bool | True | Progress flag. |


## Resolve disease names to ontology identifiers {#dnakit.databases.open_targets_get_disease_id_description_by_name-parameters}

`dnakit.databases.open_targets_get_disease_id_description_by_name(*, diseaseName=None, query=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `diseaseName` | string | None | The name of the disease to search for. |
| `query` | string | None | DNAKit optional alias for diseaseName |
| `verbose` | bool | True | Progress flag. |


## Retrieve disease synonyms {#dnakit.databases.open_targets_get_disease_synonyms_by_efo_id-parameters}

`dnakit.databases.open_targets_get_disease_synonyms_by_efo_id(*, efoId, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | Required | The EFO ID of the disease. |
| `verbose` | bool | True | Progress flag. |


## Retrieve therapeutic areas associated with a disease {#dnakit.databases.open_targets_get_disease_therapeutic_areas_by_efo_id-parameters}

`dnakit.databases.open_targets_get_disease_therapeutic_areas_by_efo_id(*, efoId, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | Required | The EFO ID of the disease. |
| `verbose` | bool | True | Progress flag. |


## Retrieve target-associated diseases and phenotypes with pagination {#dnakit.databases.open_targets_get_diseases_phenotypes_by_target_ensembl-parameters}

`dnakit.databases.open_targets_get_diseases_phenotypes_by_target_ensembl(*, ensemblId, page=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `ensemblId` | string | Required | The ensemblId of a target. |
| `page` | object | None | Pagination parameters. `count` in the response is the total across all pages; omit this to get the first page only. {&quot;properties&quot;: {&quot;index&quot;: {&quot;type&quot;: &quot;integer&quot;, &quot;description&quot;: &quot;The index of the page to retrieve (0-based).&quot;}, &quot;size&quot;: {&quot;type&quot;: &quot;integer&quot;, &quot;description&quot;: &quot;The number of associations per page.&quot;}}} |
| `verbose` | bool | True | Progress flag. |


## Retrieve significant drug adverse event records {#dnakit.databases.open_targets_get_drug_adverse_events_by_chembl_id-parameters}

`dnakit.databases.open_targets_get_drug_adverse_events_by_chembl_id(*, chemblId, page=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `chemblId` | string | Required | The ChEMBL ID of the drug. |
| `page` | object | None | Pagination settings. {&quot;properties&quot;: {&quot;index&quot;: {&quot;type&quot;: &quot;integer&quot;, &quot;description&quot;: &quot;Index of the page to fetch, starting from 0.&quot;}, &quot;size&quot;: {&quot;type&quot;: &quot;integer&quot;, &quot;description&quot;: &quot;Number of entries per page.&quot;}}} |
| `verbose` | bool | True | Progress flag. |


## Resolve a generic drug name to a ChEMBL identifier {#dnakit.databases.open_targets_get_drug_chemb_id_by_generic_name-parameters}

`dnakit.databases.open_targets_get_drug_chemb_id_by_generic_name(*, drugName=None, query=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `drugName` | string | None | The generic name of the drug for which the ID is required. |
| `query` | string | None | DNAKit optional alias for drugName |
| `verbose` | bool | True | Progress flag. |


## Retrieve drug indications by ChEMBL identifier {#dnakit.databases.open_targets_get_drug_indications_by_chembl_id-parameters}

`dnakit.databases.open_targets_get_drug_indications_by_chembl_id(*, chemblId, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `chemblId` | string | Required | The chemblId of the drug for which to retrieve treatable phenotypes information. |
| `verbose` | bool | True | Progress flag. |


## Retrieve drug mechanisms by ChEMBL identifier {#dnakit.databases.open_targets_get_drug_mechanisms_of_action_by_chembl_id-parameters}

`dnakit.databases.open_targets_get_drug_mechanisms_of_action_by_chembl_id(*, chemblId, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `chemblId` | string | Required | The ChEMBL ID of the drug. |
| `verbose` | bool | True | Progress flag. |


## Retrieve drug warnings by ChEMBL identifier {#dnakit.databases.open_targets_get_drug_warnings_by_chembl_id-parameters}

`dnakit.databases.open_targets_get_drug_warnings_by_chembl_id(*, chemblId, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `chemblId` | string | Required | The ChEMBL ID of the drug. |
| `verbose` | bool | True | Progress flag. |


## Filter target–disease evidence by data source {#dnakit.databases.open_targets_get_evidence_by_datasource-parameters}

`dnakit.databases.open_targets_get_evidence_by_datasource(*, efoId=None, ensemblId=None, gene_symbol=None, disease_name=None, datasourceIds=None, size=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | None | Disease EFO/MONDO ID (e.g., &#x27;MONDO_0005011&#x27; for Crohn disease). Alternative to disease_name. |
| `ensemblId` | string | None | Target Ensembl gene ID (e.g., &#x27;ENSG00000141510&#x27; for TP53). Alternative to gene_symbol. |
| `gene_symbol` | string | None | HGNC gene symbol (e.g., &#x27;TP53&#x27;, &#x27;BRCA1&#x27;). Auto-resolved to ensemblId. |
| `disease_name` | string | None | Disease or phenotype name (e.g., &#x27;Crohn disease&#x27;). Auto-resolved to efoId. |
| `datasourceIds` | array | None | List of datasource IDs to filter evidence. Examples: &#91;&#x27;gwas_credible_sets&#x27;, &#x27;europepmc&#x27;&#93;, &#91;&#x27;eva&#x27;, &#x27;clinical_precedence&#x27;&#93;, &#91;&#x27;intogen&#x27;, &#x27;cancer_gene_census&#x27;&#93;. For genetic (GWAS) association evidence use &#x27;gwas_credible_sets&#x27; — the retired &#x27;ot_genetics_portal&#x27; name is remapped automatically. Omit or pass empty array for all sources. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `size` | integer | None | Maximum evidence rows to return (default: 50) {&quot;default&quot;: 50} |
| `verbose` | bool | True | Progress flag. |


## Retrieve GWAS cohort and study metadata {#dnakit.databases.open_targets_get_gwas_study-parameters}

`dnakit.databases.open_targets_get_gwas_study(*, studyId, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `studyId` | string | Required | GWAS study accession ID (e.g., &#x27;GCST000392&#x27; for T1D Barrett et al., &#x27;GCST90029024&#x27; for T2D Loh et al.) |
| `verbose` | bool | True | Progress flag. |


## Retrieve publications associated with a target {#dnakit.databases.open_targets_get_publications_by_target_ensembl_id-parameters}

`dnakit.databases.open_targets_get_publications_by_target_ensembl_id(*, entityId, additionalIds=None, startYear=None, startMonth=None, endYear=None, endMonth=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `entityId` | string | Required | The ID of the entity (ensemblID). |
| `additionalIds` | array | None | List of additional IDs to include in the search. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `startYear` | integer | None | Year at the lower end of the filter. |
| `startMonth` | integer | None | Month at the lower end of the filter. |
| `endYear` | integer | None | Year at the higher end of the filter. |
| `endMonth` | integer | None | Month at the higher end of the filter. |
| `verbose` | bool | True | Progress flag. |


## Retrieve biomedical entities similar to a disease {#dnakit.databases.open_targets_get_similar_entities_by_disease_efo_id-parameters}

`dnakit.databases.open_targets_get_similar_entities_by_disease_efo_id(*, efoId, threshold, size, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | Required | The EFO ID of the disease. |
| `threshold` | number | Required | Threshold similarity between 0 and 1. Only results above threshold are returned. |
| `size` | integer | Required | Number of similar entities to fetch. |
| `verbose` | bool | True | Progress flag. |


## Retrieve credible sets from a GWAS study {#dnakit.databases.open_targets_get_study_credible_sets-parameters}

`dnakit.databases.open_targets_get_study_credible_sets(*, studyIds, size=None, index=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `studyIds` | array | Required | GWAS study accession IDs (e.g., &#91;&#x27;GCST000392&#x27;&#93;) {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `size` | integer | None | Number of credible sets to return (default 20) {&quot;default&quot;: 20} |
| `index` | integer | None | Page index for pagination (0-based, default 0) {&quot;default&quot;: 0} |
| `verbose` | bool | True | Progress flag. |


## Retrieve target functional classes {#dnakit.databases.open_targets_get_target_classes_by_ensembl_id-parameters}

`dnakit.databases.open_targets_get_target_classes_by_ensembl_id(*, ensemblId, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `ensemblId` | string | Required | The Ensembl ID of the target. |
| `verbose` | bool | True | Progress flag. |


## Retrieve GO annotations for a target {#dnakit.databases.open_targets_get_target_gene_ontology_by_ensembl_id-parameters}

`dnakit.databases.open_targets_get_target_gene_ontology_by_ensembl_id(*, ensemblId, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `ensemblId` | string | Required | The Ensembl ID of the target for which to retrieve Gene Ontology annotations. |
| `verbose` | bool | True | Progress flag. |


## Retrieve target genomic coordinates and strand {#dnakit.databases.open_targets_get_target_genomic_location_by_ensembl_id-parameters}

`dnakit.databases.open_targets_get_target_genomic_location_by_ensembl_id(*, ensemblId, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `ensemblId` | string | Required | The Ensembl ID of the target for which to retrieve genomic location information. |
| `verbose` | bool | True | Progress flag. |


## Retrieve homologues of a target {#dnakit.databases.open_targets_get_target_homologues_by_ensembl_id-parameters}

`dnakit.databases.open_targets_get_target_homologues_by_ensembl_id(*, ensemblId, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `ensemblId` | string | Required | The Ensembl ID of the target for which to retrieve homologues. |
| `verbose` | bool | True | Progress flag. |


## Resolve target names to identifiers and descriptions {#dnakit.databases.open_targets_get_target_id_description_by_name-parameters}

`dnakit.databases.open_targets_get_target_id_description_by_name(*, targetName, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `targetName` | string | Required | The name of the target for which the ID is required. |
| `verbose` | bool | True | Progress flag. |


## Retrieve target interaction partners and evidence with pagination {#dnakit.databases.open_targets_get_target_interactions_by_ensembl_id-parameters}

`dnakit.databases.open_targets_get_target_interactions_by_ensembl_id(*, ensemblId, page=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `ensemblId` | string | Required | The Ensembl ID of the target. |
| `page` | object | None | Pagination parameters. {&quot;properties&quot;: {&quot;index&quot;: {&quot;type&quot;: &quot;integer&quot;, &quot;description&quot;: &quot;The index of the page to retrieve.&quot;}, &quot;size&quot;: {&quot;type&quot;: &quot;integer&quot;, &quot;description&quot;: &quot;The number of items per page.&quot;}}} |
| `verbose` | bool | True | Progress flag. |


## Retrieve known target safety liabilities {#dnakit.databases.open_targets_get_target_safety_profile_by_ensembl_id-parameters}

`dnakit.databases.open_targets_get_target_safety_profile_by_ensembl_id(*, ensemblId, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `ensemblId` | string | Required | The Ensembl ID of the target for which to retrieve safety liabilities. |
| `verbose` | bool | True | Progress flag. |


## Retrieve target tractability across treatment modalities {#dnakit.databases.open_targets_get_target_tractability_by_ensembl_id-parameters}

`dnakit.databases.open_targets_get_target_tractability_by_ensembl_id(*, ensemblId, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `ensemblId` | string | Required | The Ensembl ID of the target. |
| `verbose` | bool | True | Progress flag. |


## Retrieve GWAS credible sets containing a variant {#dnakit.databases.open_targets_get_variant_credible_sets-parameters}

`dnakit.databases.open_targets_get_variant_credible_sets(*, variantId=None, size=None, variant_id=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `variantId` | string | None | Open Targets variant ID in format chr_position_ref_alt (e.g., &#x27;10_112998590_C_T&#x27; for rs7903146) |
| `size` | integer | None | Number of credible sets to return (default 10) {&quot;default&quot;: 10} |
| `variant_id` | string | None | DNAKit optional alias for variantId |
| `verbose` | bool | True | Progress flag. |


## Retrieve variant coordinates, frequencies and consequences {#dnakit.databases.open_targets_get_variant_info-parameters}

`dnakit.databases.open_targets_get_variant_info(*, variantId, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `variantId` | string | Required | Open Targets variant ID in format chr_position_ref_alt (e.g., &#x27;10_112998590_C_T&#x27; for rs7903146, &#x27;19_44908684_T_C&#x27; for rs429358). Use the search tool to find variant IDs from rsIDs. |
| `verbose` | bool | True | Progress flag. |


## Map disease and phenotype identifiers across databases {#dnakit.databases.open_targets_map_any_disease_id_to_all_other_ids-parameters}

`dnakit.databases.open_targets_map_any_disease_id_to_all_other_ids(*, inputId, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `inputId` | string | Required | Any known disease ID (e.g. OMIM:604302, UMLS:C0003873, ICD10:M05, MONDO_0008383, etc.) |
| `verbose` | bool | True | Progress flag. |


## Search disease, target and drug entities together {#dnakit.databases.open_targets_multi_entity_search_by_query_string-parameters}

`dnakit.databases.open_targets_multi_entity_search_by_query_string(*, queryString, entityNames=None, page=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `queryString` | string | Required | The search string for querying information. |
| `entityNames` | array | None | List of entity names to search for (e.g., target, disease, drug). {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `page` | object | None | Pagination settings with index and size. {&quot;properties&quot;: {&quot;index&quot;: {&quot;type&quot;: &quot;integer&quot;, &quot;description&quot;: &quot;Pagination index.&quot;}, &quot;size&quot;: {&quot;type&quot;: &quot;integer&quot;, &quot;description&quot;: &quot;Pagination size.&quot;}}} |
| `verbose` | bool | True | Progress flag. |


## Find GWAS studies by disease or trait {#dnakit.databases.open_targets_search_gwas_studies_by_disease-parameters}

`dnakit.databases.open_targets_search_gwas_studies_by_disease(*, diseaseIds=None, enableIndirect=None, disease_name=None, disease=None, trait=None, size=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `diseaseIds` | array | None | Disease ontology IDs (e.g., &#91;&#x27;MONDO_0005148&#x27;&#93; for type 2 diabetes, &#91;&#x27;MONDO_0004975&#x27;&#93; for Alzheimer disease) {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `enableIndirect` | boolean | None | Include studies for child disease terms (default true) {&quot;default&quot;: true} |
| `disease_name` | string | None | Disease or trait name for auto-resolution (e.g., &#x27;type 2 diabetes&#x27;). Resolved to diseaseIds automatically. |
| `disease` | string | None | Alias for disease_name. |
| `trait` | string | None | Alias for disease_name. |
| `size` | integer | None | Limit for the returned page; does not automatically retrieve later pages. {&quot;minimum&quot;: 1, &quot;maximum&quot;: 1000} |
| `verbose` | bool | True | Progress flag. |


## Retrieve IntOGen cancer driver association evidence {#dnakit.databases.open_targets_target_disease_evidence-parameters}

`dnakit.databases.open_targets_target_disease_evidence(*, efoId=None, ensemblId=None, gene_symbol=None, disease_name=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | None | EFO/MONDO disease ID (e.g., MONDO_0005011). Alternative to disease_name. |
| `ensemblId` | string | None | Ensembl gene ID (e.g., ENSG00000141510). Alternative to gene_symbol. |
| `gene_symbol` | string | None | HGNC gene symbol (e.g., &#x27;TP53&#x27;, &#x27;BRCA1&#x27;). Auto-resolved to ensemblId. |
| `disease_name` | string | None | Disease or phenotype name (e.g., &#x27;Crohn disease&#x27;, &#x27;breast carcinoma&#x27;). Auto-resolved to efoId. |
| `verbose` | bool | True | Progress flag. |


## Resolve compound names to PubChem CIDs {#dnakit.databases.pub_chem_get_cid_by_compound_name-parameters}

`dnakit.databases.pub_chem_get_cid_by_compound_name(*, name=None, compound_name=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `name` | string | None | Chemical compound name (e.g., &quot;Aspirin&quot;, &quot;Acetaminophen&quot;) or IUPAC name. Do not use disease names or medical conditions. |
| `compound_name` | string | None | Alias for name. The compound name to look up. |
| `verbose` | bool | True | Progress flag. |


## Retrieve molecular properties by PubChem CID {#dnakit.databases.pub_chem_get_compound_properties_by_cid-parameters}

`dnakit.databases.pub_chem_get_compound_properties_by_cid(*, cid, properties=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `cid` | integer | Required | PubChem compound ID to query, e.g., 2244 (Aspirin). |
| `properties` | array &#124; string &#124; null | None | Optional list of PubChem property names to retrieve (e.g. &#91;&#x27;MolecularFormula&#x27;,&#x27;InChIKey&#x27;,&#x27;XLogP&#x27;,&#x27;IUPACName&#x27;,&#x27;MolecularWeight&#x27;&#93;). If omitted, returns a default set (MolecularWeight, IUPACName, ConnectivitySMILES). For SMILES use the current names &#x27;ConnectivitySMILES&#x27; (was &#x27;CanonicalSMILES&#x27;) or &#x27;SMILES&#x27; (was &#x27;IsomericSMILES&#x27;) — a legacy name still returns data, but under the current key, not the one requested. See PubChem PUG-REST property tables for valid names. Accepts a list or a comma-separated string. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `verbose` | bool | True | Progress flag. |


## Retrieve pathway-based disease–target association scores {#dnakit.databases.reactome_disease_target_score-parameters}

`dnakit.databases.reactome_disease_target_score(*, efoId, pageSize=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `efoId` | string | Required | The EFO (Experimental Factor Ontology) ID of the disease, e.g., &#x27;MONDO_0011996&#x27; for chronic myelogenous leukemia |
| `pageSize` | integer | None | Number of results per page (default: 100, max: 100) {&quot;default&quot;: 100} |
| `verbose` | bool | True | Progress flag. |
