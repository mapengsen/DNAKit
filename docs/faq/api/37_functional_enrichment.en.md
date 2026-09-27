# Gene Set and Pathway Enrichment: parameters and notes

See the [tutorial](../../api/features/37_functional_enrichment.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

## Run gene set enrichment against an Enrichr library {#dnakit.databases.enrichr_enrich-parameters}

`dnakit.databases.enrichr_enrich(*, operation=None, gene_list, library=None, top_n=None, species=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: enrich) {&quot;const&quot;: &quot;enrich&quot;, &quot;default&quot;: &quot;enrich&quot;} |
| `gene_list` | array | Required | List of gene symbols (e.g., &#91;&#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;EGFR&#x27;&#93;). Minimum 1 gene. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `library` | string | None | Enrichr library name (default: GO_Biological_Process_2023). Common: GO_Biological_Process_2023, KEGG_2021_Human, Reactome_2022, WikiPathways_2024_Human, MSigDB_Hallmark_2020 {&quot;default&quot;: &quot;GO_Biological_Process_2023&quot;} |
| `top_n` | integer | None | Number of top enriched terms to return (default: 10) {&quot;default&quot;: 10} |
| `species` | string &#124; null | None | Organism instance to query: &#x27;human&#x27; (default), &#x27;mouse&#x27;, &#x27;fly&#x27; (FlyEnrichr), &#x27;worm&#x27; (WormEnrichr), &#x27;yeast&#x27; (YeastEnrichr), or &#x27;fish&#x27; (FishEnrichr). Gene set libraries differ per organism, and gene symbols must match that organism&#x27;s nomenclature. {&quot;enum&quot;: &#91;null, &quot;human&quot;, &quot;mouse&quot;, &quot;fly&quot;, &quot;worm&quot;, &quot;yeast&quot;, &quot;fish&quot;&#93;} |
| `verbose` | bool | True | Progress flag. |


## Retrieve connections between genes and enriched pathways {#dnakit.databases.enrichr_gene_enrichment_analysis-parameters}

`dnakit.databases.enrichr_gene_enrichment_analysis(*, gene_list, libs=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_list` | array | Required | List of gene names or symbols to analyze. At least 2 genes are required for path ranking analysis. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}, &quot;minItems&quot;: 2} |
| `libs` | array | None | List of Enrichr enrichment libraries to query (e.g. &#x27;KEGG_2021_Human&#x27;, &#x27;GO_Biological_Process_2023&#x27;). Omitting this or passing an empty list uses a broad default set of pathway/ontology libraries. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}, &quot;default&quot;: &#91;&quot;WikiPathways_2024_Human&quot;, &quot;Reactome_Pathways_2024&quot;, &quot;MSigDB_Hallmark_2020&quot;, &quot;GO_Molecular_Function_2023&quot;, &quot;GO_Biological_Process_2023&quot;&#93;} |
| `verbose` | bool | True | Progress flag. |


## Run PANTHER gene set overrepresentation analysis {#dnakit.databases.panther_enrichment-parameters}

`dnakit.databases.panther_enrichment(*, gene_list, organism=None, annotation_dataset=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_list` | string | Required | Comma-separated list of gene symbols, UniProt IDs, or Ensembl IDs. Example: &#x27;TP53,BRCA1,EGFR,KRAS,PTEN,APC,RB1&#x27;. |
| `organism` | integer &#124; null | None | NCBI taxonomy ID. Default: 9606 (human). Others: 10090 (mouse), 7227 (fly), 559292 (yeast). |
| `annotation_dataset` | string &#124; null | None | Annotation dataset for enrichment. Options: &#x27;GO:0008150&#x27; (biological_process, default), &#x27;GO:0003674&#x27; (molecular_function), &#x27;GO:0005575&#x27; (cellular_component), &#x27;ANNOT_TYPE_ID_PANTHER_PATHWAY&#x27; (PANTHER pathway), &#x27;ANNOT_TYPE_ID_REACTOME_PATHWAY&#x27; (Reactome pathway), &#x27;ANNOT_TYPE_ID_PANTHER_PC&#x27; (protein class). |
| `verbose` | bool | True | Progress flag. |


## Run Reactome pathway enrichment analysis {#dnakit.databases.reactome_analysis_pathway_enrichment-parameters}

`dnakit.databases.reactome_analysis_pathway_enrichment(*, identifiers, page_size=None, include_disease=None, projection=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `identifiers` | string | Required | Newline-separated list of gene/protein identifiers. Supports gene symbols (TP53, BRCA1), UniProt IDs (P04637, P38398), or Ensembl IDs. Example: &#x27;TP53\nBRCA1\nCDH1\nEGFR\nKRAS&#x27;. |
| `page_size` | integer &#124; null | None | Number of pathways to return (default 20, max 50). |
| `include_disease` | boolean &#124; null | None | Include disease pathways in results (default true). |
| `projection` | boolean &#124; null | None | Project identifiers to human Reactome pathways for cross-species analysis (default true). |
| `verbose` | bool | True | Progress flag. |


## Analyze functional and pathway enrichment for a protein set {#dnakit.databases.string_functional_enrichment-parameters}

`dnakit.databases.string_functional_enrichment(*, protein_ids, species=None, category=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `protein_ids` | array | Required | List of protein identifiers (UniProt IDs, gene names, Ensembl IDs). Minimum 3 proteins recommended for meaningful enrichment. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}, &quot;minItems&quot;: 2} |
| `species` | integer | None | NCBI taxonomy ID (default: 9606 for human) {&quot;default&quot;: 9606} |
| `category` | string | None | Enrichment category: &#x27;Process&#x27; (GO Biological Process), &#x27;Component&#x27; (GO Cellular Component), &#x27;Function&#x27; (GO Molecular Function), &#x27;KEGG&#x27;, &#x27;Reactome&#x27;, &#x27;WikiPathways&#x27;, &#x27;COMPARTMENTS&#x27;, &#x27;TISSUES&#x27;, &#x27;DISEASES&#x27; {&quot;enum&quot;: &#91;&quot;Process&quot;, &quot;Component&quot;, &quot;Function&quot;, &quot;KEGG&quot;, &quot;Reactome&quot;, &quot;WikiPathways&quot;, &quot;COMPARTMENTS&quot;, &quot;TISSUES&quot;, &quot;DISEASES&quot;&#93;, &quot;default&quot;: &quot;Process&quot;} |
| `verbose` | bool | True | Progress flag. |


## Run STRING enrichment from newline-separated protein identifiers {#dnakit.databases.string_get_functional_enrichment-parameters}

`dnakit.databases.string_get_functional_enrichment(*, identifiers, species=9606, category=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `identifiers` | string | Required | Protein identifiers separated by newline or carriage return; converted to a STRING identifier list. |
| `species` | integer | 9606 | Species name or taxonomy ID; type and supported values depend on the operation. {&quot;default&quot;: 9606} |
| `category` | string | None | Category or result class defined by this operation&#x27;s provider/schema. |
| `verbose` | bool | True | Progress flag. |
