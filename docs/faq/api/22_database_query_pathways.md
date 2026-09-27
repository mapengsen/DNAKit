# 通路、互作、蛋白与结构查询：参数与补充说明

用法和示例见[功能教程](../../api/features/22_database_query_pathways.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

## 获取已收录的 AlphaFold 蛋白质结构预测 {#dnakit.databases.alphafold_get_prediction-parameters}

`dnakit.databases.alphafold_get_prediction(*, qualifier=None, uniprot_id=None, uniprot_accession=None, accession=None, sequence_checksum=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `qualifier` | string | None | UniProt ACCESSION (e.g., &#x27;P69905&#x27;). Do NOT use entry names like &#x27;HBA_HUMAN&#x27;. Aliases: uniprot_id, uniprot_accession. |
| `uniprot_id` | string | None | Alias for qualifier: UniProt accession (e.g., &#x27;P69905&#x27;). |
| `uniprot_accession` | string | None | Alias for qualifier: UniProt accession (e.g., &#x27;P69905&#x27;). |
| `accession` | string | None | Alias for qualifier. UniProt accession (e.g., &#x27;P04637&#x27;) -- the same parameter name the UniProt_*_by_accession tools use. |
| `sequence_checksum` | string | None | Optional CRC64 checksum of the UniProt sequence. |
| `verbose` | bool | True | 运行提示开关。 |


## 查询实验验证的蛋白质与遗传互作 {#dnakit.databases.bio_grid_get_interactions-parameters}

`dnakit.databases.bio_grid_get_interactions(*, gene_names=None, gene_name=None, organism=None, interaction_type=None, evidence_types=None, limit=None, throughput=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_names` | array | None | List of gene names or protein identifiers (e.g., &#91;&#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;MYC&#x27;&#93;). Accepts official gene symbols. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}, &quot;minItems&quot;: 1} |
| `gene_name` | string | None | Alias for gene_names: single gene symbol (e.g., &#x27;TP53&#x27;). Converted to gene_names list automatically. |
| `organism` | string | None | Organism name (e.g., &#x27;Homo sapiens&#x27;, &#x27;Mus musculus&#x27;) or NCBI taxonomy ID (e.g., &#x27;9606&#x27; for human, &#x27;10090&#x27; for mouse). Default: 9606 (human) {&quot;default&quot;: &quot;9606&quot;} |
| `interaction_type` | string | None | Type of interaction: &#x27;physical&#x27; (protein-protein), &#x27;genetic&#x27; (epistasis, synthetic lethality), &#x27;both&#x27; (all interactions) {&quot;enum&quot;: &#91;&quot;physical&quot;, &quot;genetic&quot;, &quot;both&quot;&#93;, &quot;default&quot;: &quot;both&quot;} |
| `evidence_types` | array | None | Filter by evidence types (e.g., &#91;&#x27;Affinity Capture-MS&#x27;, &#x27;Two-hybrid&#x27;&#93; for physical, &#91;&#x27;Synthetic Lethality&#x27;, &#x27;Dosage Rescue&#x27;&#93; for genetic). Leave empty for all evidence types. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `limit` | integer | None | Maximum number of interactions to return (default: 100, max: 10000) {&quot;minimum&quot;: 1, &quot;maximum&quot;: 10000, &quot;default&quot;: 100} |
| `throughput` | string | None | Filter by throughput: &#x27;low&#x27; (low-throughput studies), &#x27;high&#x27; (high-throughput screens), or leave empty for all {&quot;enum&quot;: &#91;&quot;low&quot;, &quot;high&quot;, &quot;&quot;&#93;} |
| `verbose` | bool | True | 运行提示开关。 |

运行前配置：`BIOGRID_ACCESS_KEY`。


## 检索病原体 B 细胞与 T 细胞表位 {#dnakit.databases.bvbrc_search_epitopes-parameters}

`dnakit.databases.bvbrc_search_epitopes(*, taxon_id=None, protein_name=None, epitope_type=None, organism=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `taxon_id` | string &#124; null | None | NCBI Taxonomy ID for the pathogen organism. Examples: &#x27;2697049&#x27; (SARS-CoV-2), &#x27;11320&#x27; (Influenza A), &#x27;1773&#x27; (M. tuberculosis). |
| `protein_name` | string &#124; null | None | Target protein name to filter epitopes. Examples: &#x27;Spike glycoprotein&#x27;, &#x27;Nucleocapsid protein&#x27;, &#x27;Replicase polyprotein 1ab&#x27;. |
| `epitope_type` | string &#124; null | None | Type of epitope. Options: &#x27;Linear peptide&#x27;, &#x27;Discontinuous peptide&#x27;. Most epitopes are linear peptides. |
| `organism` | string &#124; null | None | Organism name keyword to search across epitope records. Example: &#x27;coronavirus&#x27;, &#x27;influenza&#x27;. |
| `limit` | integer &#124; null | None | Maximum number of results. Default: 25. Max: 100. |
| `verbose` | bool | True | 运行提示开关。 |


## 查询基因或蛋白质的 GO 注释 {#dnakit.databases.go_get_annotations_for_gene-parameters}

`dnakit.databases.go_get_annotations_for_gene(*, gene_id, rows=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填 | A gene identifier such as gene symbol (e.g., &#x27;TP53&#x27;) or database ID. |
| `rows` | integer | None | Maximum number of annotation rows to return. Default: 100. Raise it (e.g. 2000) for a heavily annotated gene such as TP53, which has 1812 annotation rows; lower it to keep responses small. {&quot;minimum&quot;: 1, &quot;default&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 按关键词检索 GO 术语 {#dnakit.databases.go_search_terms-parameters}

`dnakit.databases.go_search_terms(*, query, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填 | The keyword to search for, e.g., &#x27;apoptosis&#x27; or &#x27;kinase activity&#x27;. |
| `verbose` | bool | True | 运行提示开关。 |


## 获取组织特异性蛋白互作网络 {#dnakit.databases.humanbase_ppi_analysis-parameters}

`dnakit.databases.humanbase_ppi_analysis(*, gene_list=None, tissue=None, max_node=None, interaction=None, string_mode=None, genes=None, top_n=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_list` | array | None | List of gene names or symbols to analyze for protein-protein interactions. The gene name should be the official gene symbol, not the synonym. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}, &quot;minItems&quot;: 1} |
| `tissue` | string | None | Tissue type for tissue-specific interactions. Use simple names like &#x27;brain&#x27;, &#x27;heart&#x27;, &#x27;liver&#x27;, &#x27;kidney&#x27;, &#x27;lung&#x27;, &#x27;blood&#x27;. For v3 networks with more cell-type specificity, append &#x27;-v3&#x27; (e.g., &#x27;brain-v3&#x27;, &#x27;liver-v3&#x27;). {&quot;default&quot;: &quot;brain&quot;} |
| `max_node` | integer | None | Maximum number of nodes to retrieve in the interaction network. Warning: the more nodes, the more time it takes to retrieve the data. Default is 10 (~30 seconds). {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `interaction` | string | None | Specific interaction type to filter by. Available types: &#x27;co-expression&#x27;, &#x27;interaction&#x27;, &#x27;tf-binding&#x27;, &#x27;gsea-microrna-targets&#x27;, &#x27;gsea-perturbations&#x27;. If not specified, all types will be included. {&quot;enum&quot;: &#91;&quot;co-expression&quot;, &quot;interaction&quot;, &quot;tf-binding&quot;, &quot;gsea-microrna-targets&quot;, &quot;gsea-perturbations&quot;&#93;} |
| `string_mode` | boolean | None | Whether to return the result in string mode. If True, the result will be a string of the network graph and the biological processes. If False, the result will be a NetworkX graph and a list of biological processes. {&quot;default&quot;: true} |
| `genes` | array | None | Alias for gene_list. List of gene symbols to analyze. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `top_n` | integer | None | Alias for max_node. Maximum number of network nodes to retrieve. {&quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询免疫表位关联的抗原 {#dnakit.databases.iedb_get_epitope_antigens-parameters}

`dnakit.databases.iedb_get_epitope_antigens(*, structure_id, limit=None, offset=None, order=None, select=None, filters=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `structure_id` | integer | 必填 | IEDB epitope structure ID (from `iedb_search_epitopes`). |
| `limit` | integer | None | Maximum number of rows to return. {&quot;default&quot;: 10} |
| `offset` | integer | None | Pagination offset. {&quot;default&quot;: 0} |
| `order` | string | None | Optional PostgREST order clause. |
| `select` | object &#124; None | None | Optional projection list (comma-separated string or array of field names). {&quot;anyOf&quot;: &#91;{&quot;type&quot;: &quot;string&quot;}, {&quot;type&quot;: &quot;array&quot;, &quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}}&#93;} |
| `filters` | object | None | Advanced PostgREST filters mapping column -&gt; filter expression. {&quot;additionalProperties&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询免疫表位的 MHC 限制性关联 {#dnakit.databases.iedb_get_epitope_mhc-parameters}

`dnakit.databases.iedb_get_epitope_mhc(*, structure_id, limit=None, offset=None, order=None, select=None, filters=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `structure_id` | integer | 必填 | IEDB epitope structure ID (from `iedb_search_epitopes`). |
| `limit` | integer | None | Maximum number of rows to return. {&quot;default&quot;: 10} |
| `offset` | integer | None | Pagination offset. {&quot;default&quot;: 0} |
| `order` | string | None | Optional PostgREST order clause. |
| `select` | object &#124; None | None | Optional projection list (comma-separated string or array of field names). {&quot;anyOf&quot;: &#91;{&quot;type&quot;: &quot;string&quot;}, {&quot;type&quot;: &quot;array&quot;, &quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}}&#93;} |
| `filters` | object | None | Advanced PostgREST filters mapping column -&gt; filter expression. {&quot;additionalProperties&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `verbose` | bool | True | 运行提示开关。 |


## 检索 B 细胞表位实验记录 {#dnakit.databases.iedb_search_bcell-parameters}

`dnakit.databases.iedb_search_bcell(*, limit=None, offset=None, order=None, select=None, filters=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `limit` | integer | None | Maximum number of rows to return. {&quot;default&quot;: 10} |
| `offset` | integer | None | Pagination offset. {&quot;default&quot;: 0} |
| `order` | string | None | Optional PostgREST order clause. |
| `select` | object &#124; None | None | Optional projection list (comma-separated string or array of field names). {&quot;anyOf&quot;: &#91;{&quot;type&quot;: &quot;string&quot;}, {&quot;type&quot;: &quot;array&quot;, &quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}}&#93;} |
| `filters` | object | None | Advanced PostgREST filters mapping column -&gt; filter expression. {&quot;additionalProperties&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `verbose` | bool | True | 运行提示开关。 |


## 检索免疫表位及其结构编号 {#dnakit.databases.iedb_search_epitopes-parameters}

`dnakit.databases.iedb_search_epitopes(*, sequence_contains=None, structure_type=None, limit=None, offset=None, order=None, select=None, filters=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `sequence_contains` | string | None | Sequence fragment to search within `linear_sequence` using a contains match. Example: &#x27;KVF&#x27; finds epitopes containing that motif. |
| `structure_type` | string | None | Filter by epitope structure type. Example: &#x27;Linear peptide&#x27;. {&quot;enum&quot;: &#91;&quot;Linear peptide&quot;, &quot;Discontinuous peptide&quot;, &quot;Non-peptidic&quot;&#93;} |
| `limit` | integer | None | Maximum number of rows to return. {&quot;default&quot;: 10} |
| `offset` | integer | None | Pagination offset. {&quot;default&quot;: 0} |
| `order` | string | None | Optional PostgREST order clause (e.g., &#x27;structure_id.asc&#x27;). |
| `select` | object &#124; None | None | Optional projection list to reduce payload size. Provide a comma-separated string or an array of field names (e.g., &#91;&#x27;structure_id&#x27;,&#x27;linear_sequence&#x27;&#93;). {&quot;anyOf&quot;: &#91;{&quot;type&quot;: &quot;string&quot;}, {&quot;type&quot;: &quot;array&quot;, &quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}}&#93;} |
| `filters` | object | None | Advanced PostgREST filters mapping column -&gt; filter expression (e.g., {&quot;linear_sequence&quot;:&quot;ilike.*KVF*&quot;}). {&quot;additionalProperties&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `verbose` | bool | True | 运行提示开关。 |


## 检索 MHC 相关实验记录 {#dnakit.databases.iedb_search_mhc-parameters}

`dnakit.databases.iedb_search_mhc(*, limit=None, offset=None, order=None, select=None, filters=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `limit` | integer | None | Maximum number of rows to return. {&quot;default&quot;: 10} |
| `offset` | integer | None | Pagination offset. {&quot;default&quot;: 0} |
| `order` | string | None | Optional PostgREST order clause. |
| `select` | object &#124; None | None | Optional projection list (comma-separated string or array of field names). {&quot;anyOf&quot;: &#91;{&quot;type&quot;: &quot;string&quot;}, {&quot;type&quot;: &quot;array&quot;, &quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}}&#93;} |
| `filters` | object | None | Advanced PostgREST filters mapping column -&gt; filter expression. {&quot;additionalProperties&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `verbose` | bool | True | 运行提示开关。 |


## 检索 T 细胞免疫反应实验 {#dnakit.databases.iedb_search_tcell_assays-parameters}

`dnakit.databases.iedb_search_tcell_assays(*, sequence=None, sequence_contains=None, mhc_class=None, qualitative_measure=None, limit=None, offset=None, select=None, filters=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `sequence` | string | None | Exact peptide sequence to search (e.g., GILGFVFTL for influenza M1 epitope). |
| `sequence_contains` | string | None | Partial peptide sequence (substring match). Example: &#x27;SIINFEKL&#x27;. |
| `mhc_class` | string | None | MHC class restriction (I or II). {&quot;enum&quot;: &#91;&quot;I&quot;, &quot;II&quot;&#93;} |
| `qualitative_measure` | string | None | Filter by assay result. {&quot;enum&quot;: &#91;&quot;Positive&quot;, &quot;Negative&quot;, &quot;Positive-High&quot;, &quot;Positive-Intermediate&quot;, &quot;Positive-Low&quot;&#93;} |
| `limit` | integer | None | Maximum rows to return. {&quot;default&quot;: 10} |
| `offset` | integer | None | Pagination offset. {&quot;default&quot;: 0} |
| `select` | object &#124; None | None | Columns to return (array of strings). Example: &#91;&#x27;tcell_id&#x27;,&#x27;linear_sequence&#x27;,&#x27;mhc_restriction&#x27;,&#x27;qualitative_measure&#x27;&#93;. {&quot;oneOf&quot;: &#91;{&quot;type&quot;: &quot;string&quot;}, {&quot;type&quot;: &quot;array&quot;, &quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}}&#93;} |
| `filters` | object | None | Advanced PostgREST filters (e.g., {&quot;source_organism_iri&quot;:&quot;eq.NCBITaxon:11320&quot;}). {&quot;additionalProperties&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询免疫遗传学资源与基因命名规则 {#dnakit.databases.imgt_get_gene_info-parameters}

`dnakit.databases.imgt_get_gene_info(*, operation=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | 提供方操作模式；默认值与允许值见本行约束。 {&quot;const&quot;: &quot;get_gene_info&quot;, &quot;default&quot;: &quot;get_gene_info&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 定位免疫球蛋白与 TCR 胚系基因检索入口 {#dnakit.databases.imgt_search_genes-parameters}

`dnakit.databases.imgt_search_genes(*, operation=None, query=None, gene_type=None, species=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | 提供方操作模式；默认值与允许值见本行约束。 {&quot;const&quot;: &quot;search_genes&quot;, &quot;default&quot;: &quot;search_genes&quot;} |
| `query` | string | None | Search query (gene name) |
| `gene_type` | string | None | Gene type: IGHV, IGHD, IGHJ, IGKV, IGLV, TRAV, TRBV, etc. |
| `species` | string | None | Species name (default: Homo sapiens) {&quot;default&quot;: &quot;Homo sapiens&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 按疾病与实验方法检索免疫学研究 {#dnakit.databases.imm_port_search_studies-parameters}

`dnakit.databases.imm_port_search_studies(*, query, condition_or_disease=None, assay_method=None, research_focus=None, species=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填 | Search keywords (e.g., &#x27;influenza vaccine&#x27;, &#x27;COVID-19&#x27;, &#x27;T cell response&#x27;, &#x27;malaria&#x27;) |
| `condition_or_disease` | string | None | Filter by disease/condition (e.g., &#x27;influenza&#x27;, &#x27;malaria&#x27;, &#x27;HIV&#x27;, &#x27;dengue&#x27;, &#x27;tuberculosis&#x27;). Client-side substring match. |
| `assay_method` | string | None | Filter by assay method (e.g., &#x27;Flow Cytometry&#x27;, &#x27;ELISPOT&#x27;, &#x27;RNA sequencing&#x27;, &#x27;ELISA&#x27;, &#x27;CyTOF&#x27;, &#x27;Luminex xMAP&#x27;). Client-side substring match. |
| `research_focus` | string | None | Filter by research focus (e.g., &#x27;Vaccine Response&#x27;, &#x27;Infection Response&#x27;, &#x27;Immune Response&#x27;, &#x27;Transplantation&#x27;, &#x27;Autoimmune&#x27;). Client-side substring match. |
| `species` | string | None | Filter by species (e.g., &#x27;Homo sapiens&#x27;, &#x27;Mus musculus&#x27;). Applied server-side. |
| `limit` | integer | None | Maximum number of studies to return (1-100, default 10) {&quot;default&quot;: 10} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询蛋白质复合物组成与功能 {#dnakit.databases.intact_get_complex_details-parameters}

`dnakit.databases.intact_get_complex_details(*, complex_ac, format=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `complex_ac` | string | 必填 | Complex AC (Complex Accession) in format &#x27;CPX-XXXXX&#x27; (e.g., &#x27;CPX-915&#x27;, &#x27;CPX-10081&#x27;). Get complex ACs from intact_get_interactions_by_complex search results. |
| `format` | string | None | 返回格式；允许值与服务默认值见当前接口说明。 {&quot;enum&quot;: &#91;&quot;json&quot;, &quot;xml&quot;&#93;, &quot;default&quot;: &quot;json&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询以指定蛋白质为中心的互作网络 {#dnakit.databases.intact_get_interaction_network-parameters}

`dnakit.databases.intact_get_interaction_network(*, identifier=None, uniprot_id=None, protein_id=None, gene_symbol=None, gene_name=None, depth=None, limit=None, size=None, format=None, protein_name=None, protein=None, gene=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `identifier` | string | None | IntAct identifier, UniProt ID, or gene name. Aliases: uniprot_id, protein_id, gene_symbol, gene_name. |
| `uniprot_id` | string | None | Alias for identifier: UniProt accession (e.g., &#x27;P04637&#x27;). |
| `protein_id` | string | None | Alias for identifier: protein identifier. |
| `gene_symbol` | string | None | Alias for identifier: gene symbol (e.g., &#x27;BRCA1&#x27;). |
| `gene_name` | string | None | Alias for identifier: gene name. |
| `depth` | integer | None | Network depth: 1 for direct interactions only, 2 for 2-hop network, etc. (default: 1) {&quot;default&quot;: 1, &quot;minimum&quot;: 1, &quot;maximum&quot;: 3} |
| `limit` | integer | None | Maximum number of interactions to return (default: 50, max: 200). Alias: size. {&quot;default&quot;: 50} |
| `size` | integer | None | Alias for limit. Maximum number of interactions to return (default: 50). {&quot;default&quot;: 50} |
| `format` | string | None | 返回格式；允许值与服务默认值见当前接口说明。 {&quot;enum&quot;: &#91;&quot;json&quot;, &quot;xml&quot;&#93;, &quot;default&quot;: &quot;json&quot;} |
| `protein_name` | string | None | Alias for gene_symbol/identifier. Common protein name (e.g., MDM2, TP53). |
| `protein` | string | None | Alias for identifier. Gene symbol or protein name (e.g., &#x27;TP53&#x27;, &#x27;BRCA1&#x27;). |
| `gene` | string | None | DNAKit optional alias for gene_symbol |
| `verbose` | bool | True | 运行提示开关。 |


## 查询指定蛋白质或基因的互作记录 {#dnakit.databases.intact_get_interactions-parameters}

`dnakit.databases.intact_get_interactions(*, identifier, format=None, size=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `identifier` | string | 必填 | IntAct identifier, UniProt ID, or gene name |
| `format` | string | None | Response format. Only &#x27;json&#x27; is produced by this endpoint. {&quot;enum&quot;: &#91;&quot;json&quot;&#93;, &quot;default&quot;: &quot;json&quot;} |
| `size` | integer &#124; null | None | Maximum number of interactions to return (default 25). IntAct entries commonly have far more -- the response reports the true total as &#x27;hitCount&#x27;, so raise &#x27;size&#x27; to that value to retrieve them all. |
| `verbose` | bool | True | 运行提示开关。 |


## 按关键词检索 IntAct 分子互作 {#dnakit.databases.intact_search_interactions-parameters}

`dnakit.databases.intact_search_interactions(*, query, first=None, max=None, format=None, species=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填 | Search query string (e.g., protein name, gene name, interaction type) |
| `first` | integer | None | First result index for pagination (default: 0) {&quot;default&quot;: 0, &quot;minimum&quot;: 0} |
| `max` | integer | None | Maximum number of results to return (default: 25, max: 100) {&quot;default&quot;: 25, &quot;maximum&quot;: 100} |
| `format` | string | None | 返回格式；允许值与服务默认值见当前接口说明。 {&quot;enum&quot;: &#91;&quot;json&quot;, &quot;xml&quot;&#93;, &quot;default&quot;: &quot;json&quot;} |
| `species` | string | None | 物种名或 Taxonomy ID；类型与支持值取决于当前接口。 |
| `verbose` | bool | True | 运行提示开关。 |


## 查询蛋白质关联的结构域与家族条目 {#dnakit.databases.inter_pro_get_entries_for_protein-parameters}

`dnakit.databases.inter_pro_get_entries_for_protein(*, accession, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string | 必填 | UniProt protein accession. Examples: &#x27;P04637&#x27; (TP53, 9 entries), &#x27;P00533&#x27; (EGFR), &#x27;P01308&#x27; (Insulin), &#x27;P38398&#x27; (BRCA1). |
| `verbose` | bool | True | 运行提示开关。 |


## 查询蛋白质结构域与功能注释 {#dnakit.databases.inter_pro_get_protein_domains-parameters}

`dnakit.databases.inter_pro_get_protein_domains(*, protein_id=None, accession=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `protein_id` | string | None | UniProt protein ID (e.g., P05067, Q9Y6K9) |
| `accession` | string | None | DNAKit optional alias for protein_id |
| `verbose` | bool | True | 运行提示开关。 |


## 按关键词与物种检索 KEGG 基因 {#dnakit.databases.kegg_find_genes-parameters}

`dnakit.databases.kegg_find_genes(*, keyword=None, organism=None, query=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `keyword` | string | None | Search keyword for gene names or descriptions |
| `organism` | string | None | Organism code (e.g., &#x27;hsa&#x27; for human, &#x27;mmu&#x27; for mouse). Optional - searches all organisms if not specified {&quot;default&quot;: &quot;&quot;} |
| `query` | string | None | DNAKit optional alias for keyword |
| `verbose` | bool | True | 运行提示开关。 |


## 查询 KEGG 基因序列、功能与通路 {#dnakit.databases.kegg_get_gene_info-parameters}

`dnakit.databases.kegg_get_gene_info(*, gene_id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填 | KEGG gene identifier (e.g., &#x27;hsa:348&#x27;, &#x27;hsa:3480&#x27;) |
| `verbose` | bool | True | 运行提示开关。 |


## 查询基因参与的 KEGG 通路 {#dnakit.databases.kegg_get_gene_pathways-parameters}

`dnakit.databases.kegg_get_gene_pathways(*, gene_id=None, gene_symbol=None, gene=None, organism=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | None | KEGG gene identifier in organism:id format. Examples: &#x27;hsa:7157&#x27; (human TP53), &#x27;hsa:672&#x27; (human BRCA1), &#x27;mmu:22059&#x27; (mouse Trp53). Use organism prefix: hsa=human, mmu=mouse, rno=rat, dme=fly, sce=yeast. |
| `gene_symbol` | string | None | Alias for gene_id: human gene symbol (e.g. &#x27;TP53&#x27;, &#x27;BRCA1&#x27;). Auto-resolved to KEGG gene ID using organism=&#x27;hsa&#x27;. |
| `gene` | string | None | Alias for gene_id: human gene symbol (e.g. &#x27;TP53&#x27;, &#x27;EGFR&#x27;). |
| `organism` | string | None | KEGG organism code for gene symbol resolution (default: &#x27;hsa&#x27; for human). Used only when gene_symbol is provided instead of gene_id. {&quot;default&quot;: &quot;hsa&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 列出 KEGG 通路中的成员基因 {#dnakit.databases.kegg_get_pathway_genes-parameters}

`dnakit.databases.kegg_get_pathway_genes(*, pathway_id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `pathway_id` | string | 必填 | KEGG pathway identifier (e.g., &#x27;hsa04115&#x27; for p53 signaling pathway, &#x27;hsa04110&#x27; for cell cycle, &#x27;hsa04210&#x27; for apoptosis). Use organism prefix: hsa=human, mmu=mouse. |
| `verbose` | bool | True | 运行提示开关。 |


## 查询 KEGG 通路的基因、化合物与反应 {#dnakit.databases.kegg_get_pathway_info-parameters}

`dnakit.databases.kegg_get_pathway_info(*, pathway_id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `pathway_id` | string | 必填 | KEGG pathway identifier (e.g., &#x27;hsa00010&#x27;, &#x27;path:hsa00010&#x27;) |
| `verbose` | bool | True | 运行提示开关。 |


## 按关键词检索 KEGG 通路 {#dnakit.databases.kegg_search_pathway-parameters}

`dnakit.databases.kegg_search_pathway(*, keyword, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `keyword` | string | 必填 | Search keyword for pathway names or descriptions (e.g., &#x27;diabetes&#x27;, &#x27;metabolism&#x27;) |
| `verbose` | bool | True | 运行提示开关。 |


## 按编号查询 NCBI 蛋白质记录摘要 {#dnakit.databases.ncbiprotein_get_summary-parameters}

`dnakit.databases.ncbiprotein_get_summary(*, id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `id` | string | 必填 | NCBI Protein ID (GI number or accession). Can be a single ID or comma-separated list. Examples: &#x27;119395750&#x27; (keratin 1, NP_006112), &#x27;4507879&#x27; (TP53, NP_000537), &#x27;4503747&#x27; (EGFR, NP_005219). GI numbers can be obtained from NCBI Gene genomic info or from BLAST results. |
| `verbose` | bool | True | 运行提示开关。 |


## 查询 PANTHER 跨物种直系同源基因 {#dnakit.databases.panther_ortholog-parameters}

`dnakit.databases.panther_ortholog(*, gene_id, organism=None, target_organism=None, ortholog_type=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填 | Gene identifier - UniProt accession, Ensembl ID, or gene symbol. Examples: &#x27;P04637&#x27; (TP53), &#x27;P38398&#x27; (BRCA1). |
| `organism` | integer &#124; null | None | Source organism NCBI taxonomy ID. Default: 9606 (human). |
| `target_organism` | integer &#124; null | None | Target organism NCBI taxonomy ID. Default: 10090 (mouse). Others: 10116 (rat), 7227 (fly), 6239 (worm). |
| `ortholog_type` | string &#124; null | None | Ortholog type filter. Options: &#x27;LDO&#x27; (least diverged ortholog, default), &#x27;O&#x27; (all orthologs), &#x27;P&#x27; (paralogs), &#x27;X&#x27; (xenologs), &#x27;LDX&#x27; (least diverged xenolog), &#x27;all&#x27; (all types). |
| `verbose` | bool | True | 运行提示开关。 |


## 读取 PDB 结构坐标文件文本 {#dnakit.databases.pdb_get_structure-parameters}

`dnakit.databases.pdb_get_structure(*, pdb_id, format='cif', max_bytes=10000000, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `pdb_id` | string | 必填 | PDB 结构编号，支持 4 字符编号或 pdb_ 加 8 字符编号。 |
| `format` | string | &#x27;cif&#x27; | 返回格式；允许值与服务默认值见当前接口说明。 {&quot;default&quot;: &quot;cif&quot;} |
| `max_bytes` | integer | 10000000 | 坐标文件大小上限，1–100000000 字节。 {&quot;default&quot;: 10000000} |
| `verbose` | bool | True | 运行提示开关。 |


## 映射 PDB 链与 UniProt 残基区间 {#dnakit.databases.pdbe_get_uniprot_mappings-parameters}

`dnakit.databases.pdbe_get_uniprot_mappings(*, pdb_id=None, uniprot_id=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `pdb_id` | string | None | PDB identifier (4-character code). Examples: &#x27;4hhb&#x27; (hemoglobin), &#x27;1cbs&#x27; (CRABP2), &#x27;6lu7&#x27; (SARS-CoV-2 main protease). |
| `uniprot_id` | string | None | UniProt accession；反向结构映射时不要同时传 pdb_id。 |
| `limit` | integer | None | 当前请求的返回数量上限；不是全库命中总数。 {&quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 按覆盖度与分辨率筛选蛋白质结构 {#dnakit.databases.pdbe_sifts_get_best_structures-parameters}

`dnakit.databases.pdbe_sifts_get_best_structures(*, uniprot_accession, limit=None, offset=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `uniprot_accession` | string | 必填 | UniProt accession for the protein. Examples: &#x27;P04637&#x27; (TP53), &#x27;P00533&#x27; (EGFR), &#x27;P38398&#x27; (BRCA1), &#x27;P01308&#x27; (Insulin). |
| `limit` | integer &#124; null | None | Maximum number of entries to return (default 50, max 500). Proteins such as P04637 (TP53) have hundreds of structure-chain entries; raise this or use &#x27;offset&#x27; to reach past the first page. {&quot;default&quot;: 50} |
| `offset` | integer &#124; null | None | Zero-based index of the first entry to return, for paging through the ranked list (e.g. offset=50 returns entries 51-100). Default: 0. {&quot;default&quot;: 0} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询蛋白质 Pfam 结构域与位置 {#dnakit.databases.pfam_get_protein_annotations-parameters}

`dnakit.databases.pfam_get_protein_annotations(*, accession, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string | 必填 | UniProt protein accession. Examples: &#x27;P04637&#x27; (TP53, 4 domains), &#x27;P00533&#x27; (EGFR), &#x27;P68431&#x27; (Histone H3.1), &#x27;P38398&#x27; (BRCA1), &#x27;P01308&#x27; (Insulin). |
| `verbose` | bool | True | 运行提示开关。 |


## 查询植物通路的反应与注释 {#dnakit.databases.plant_reactome_get_pathway-parameters}

`dnakit.databases.plant_reactome_get_pathway(*, pathway_id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `pathway_id` | string | 必填 | Plant Reactome stable identifier. Format: R-{species_code}-{number}. Examples: &#x27;R-OSA-1119312&#x27; (Photorespiration in rice), &#x27;R-OSA-1119519&#x27; (Calvin cycle in rice). |
| `verbose` | bool | True | 运行提示开关。 |


## 列出植物通路数据库覆盖的物种 {#dnakit.databases.plant_reactome_list_species-parameters}

`dnakit.databases.plant_reactome_list_species(*, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `verbose` | bool | True | 运行提示开关。 |


## 检索植物代谢与调控通路 {#dnakit.databases.plant_reactome_search_pathways-parameters}

`dnakit.databases.plant_reactome_search_pathways(*, query, species=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填 | Search query for plant pathways. Examples: &#x27;photosynthesis&#x27;, &#x27;Calvin cycle&#x27;, &#x27;flavonoid biosynthesis&#x27;, &#x27;nitrogen fixation&#x27;, &#x27;auxin signaling&#x27;. |
| `species` | string &#124; null | None | Species name to filter results. Examples: &#x27;Oryza sativa&#x27;, &#x27;Arabidopsis thaliana&#x27;, &#x27;Zea mays&#x27;. Default: all species. |
| `verbose` | bool | True | 运行提示开关。 |


## 查询蛋白质序列特征与修饰位点 {#dnakit.databases.proteins_api_get_features-parameters}

`dnakit.databases.proteins_api_get_features(*, accession, format=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string &#124; array | 必填 | UniProt protein accession(s). Can be a single accession (e.g., &#x27;P05067&#x27;), comma-separated string (e.g., &#x27;P05067,P04637&#x27;), or list (e.g., &#91;&#x27;P05067&#x27;, &#x27;P04637&#x27;&#93;). Use UniProt search or ebi_search_domain with domain=&#x27;uniprot&#x27; to find valid accessions. |
| `format` | string | None | Response format. JSON is recommended for most use cases. {&quot;enum&quot;: &#91;&quot;json&quot;, &quot;xml&quot;&#93;, &quot;default&quot;: &quot;json&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 按基因、蛋白名称或编号检索蛋白质 {#dnakit.databases.proteins_api_search-parameters}

`dnakit.databases.proteins_api_search(*, query=None, size=None, offset=None, format=None, gene_symbol=None, gene=None, reviewed=None, organism=None, taxid=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | None | Search query. Can be: gene name (e.g., &#x27;BRCA1&#x27;, &#x27;TP53&#x27;), protein name (e.g., &#x27;insulin&#x27;, &#x27;hemoglobin&#x27;), or UniProt accession (e.g., &#x27;P05067&#x27;). The tool automatically determines the appropriate search parameter. |
| `size` | integer | None | Maximum number of results to return (default: 25, max: 100) {&quot;default&quot;: 25, &quot;maximum&quot;: 100} |
| `offset` | integer | None | Offset for pagination (default: 0) {&quot;default&quot;: 0} |
| `format` | string | None | Response format. JSON is recommended for most use cases. {&quot;enum&quot;: &#91;&quot;json&quot;, &quot;xml&quot;&#93;, &quot;default&quot;: &quot;json&quot;} |
| `gene_symbol` | string | None | Alias for query. Gene symbol to search (e.g., TP53, BRCA1). |
| `gene` | string | None | Alias for query. Gene name to search. |
| `reviewed` | boolean | None | If true, return only reviewed Swiss-Prot entries. If false or omitted, includes unreviewed TrEMBL entries. |
| `organism` | string | None | Restrict the search to an organism by NAME (e.g. &#x27;Klebsiella pneumoniae&#x27;, &#x27;Escherichia coli&#x27;); use &#x27;taxid&#x27; for numeric NCBI taxonomy ids. A purely numeric value here (e.g. &#x27;9606&#x27;) is auto-corrected to a taxid lookup and the response reports it in &#x27;normalization_note&#x27;. By default a gene/protein-name search is biased toward human (taxid 9606); set this for non-human proteins. If omitted and no human match is found, the search auto-widens to all organisms. |
| `taxid` | string | None | Restrict the search to an NCBI taxonomy id (e.g. &#x27;573&#x27; for Klebsiella pneumoniae, &#x27;562&#x27; for E. coli). Alternative to &#x27;organism&#x27; for precise scoping. |
| `verbose` | bool | True | 运行提示开关。 |


## 查询基因产物的 GO 注释与证据代码 {#dnakit.databases.quick_go_annotations_by_gene-parameters}

`dnakit.databases.quick_go_annotations_by_gene(*, gene_product_id, aspect=None, taxon_id=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_product_id` | string | 必填 | Gene product identifier, typically a UniProtKB accession. Format: &#x27;UniProtKB:P04637&#x27; (TP53), &#x27;UniProtKB:P38398&#x27; (BRCA1), &#x27;UniProtKB:P00533&#x27; (EGFR). |
| `aspect` | string &#124; null | None | Filter by GO aspect: &#x27;biological_process&#x27;, &#x27;molecular_function&#x27;, or &#x27;cellular_component&#x27;. Default: all aspects. |
| `taxon_id` | integer &#124; null | None | Filter by NCBI taxonomy ID. Example: 9606 (human), 10090 (mouse). |
| `limit` | integer | None | Maximum number of annotations to return. Default: 25, max: 100. {&quot;default&quot;: 25} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询 PDB 结构实验方法与元数据 {#dnakit.databases.rcsbdata_get_entry-parameters}

`dnakit.databases.rcsbdata_get_entry(*, pdb_id=None, entry_id=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `pdb_id` | string | None | PDB entry ID (4 characters). Examples: &#x27;4HHB&#x27; (hemoglobin), &#x27;1TUP&#x27; (p53-DNA complex), &#x27;1M17&#x27; (EGFR kinase), &#x27;6LU7&#x27; (SARS-CoV-2 main protease). |
| `entry_id` | string | None | Alias for pdb_id (the RCSB-native term; the REST endpoint is /core/entry/{id}). Provide either pdb_id or entry_id. |
| `verbose` | bool | True | 运行提示开关。 |

以下参数组至少提供一组：`pdb_id`; `entry_id`。


## 检索 Reactome 通路、反应与生物实体 {#dnakit.databases.reactome_content_search-parameters}

`dnakit.databases.reactome_content_search(*, query, species=None, types=None, rows=None, start=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填 | Search query. Examples: &#x27;apoptosis&#x27;, &#x27;TP53&#x27;, &#x27;cell cycle&#x27;, &#x27;DNA repair&#x27;, &#x27;insulin signaling&#x27;. |
| `species` | string | None | Species name for filtering results. Default: &#x27;Homo sapiens&#x27;. Other options: &#x27;Mus musculus&#x27;, &#x27;Rattus norvegicus&#x27;, &#x27;Danio rerio&#x27;, &#x27;Drosophila melanogaster&#x27;. {&quot;default&quot;: &quot;Homo sapiens&quot;} |
| `types` | string | None | Entity type to search. Options: &#x27;Pathway&#x27;, &#x27;Reaction&#x27;, &#x27;Complex&#x27;, &#x27;Protein&#x27;, &#x27;SmallMolecule&#x27;. Default: &#x27;Pathway&#x27;. {&quot;default&quot;: &quot;Pathway&quot;} |
| `rows` | integer &#124; null | None | Number of results to return per page (1-1000). Omit to use Reactome&#x27;s default page of 10. Queries often match far more than that (&#x27;DNA repair&#x27; in Homo sapiens matches 455 pathways), so raise this to see beyond the first page. Compare the returned &#x27;count&#x27; against &#x27;total_results&#x27; to know how much you are holding. |
| `start` | integer &#124; null | None | Zero-based offset of the first result to return, for paging through large result sets. Use together with &#x27;rows&#x27; (e.g. rows=50 with start=0, then start=50). Omit to start at the first result. |
| `verbose` | bool | True | 运行提示开关。 |


## 查询生物反应的参与分子与复合物 {#dnakit.databases.reactome_get_participants-parameters}

`dnakit.databases.reactome_get_participants(*, stId, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `stId` | string | 必填 | Event Stable ID (reaction or pathway, e.g., &#x27;R-HSA-111289&#x27;). To find reaction IDs, use Reactome_get_pathway_reactions with a pathway ID. To find pathway IDs, use Reactome_list_top_pathways. |
| `verbose` | bool | True | 运行提示开关。 |


## 查询 Reactome 通路完整信息 {#dnakit.databases.reactome_get_pathway-parameters}

`dnakit.databases.reactome_get_pathway(*, stId, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `stId` | string | 必填 | Pathway Stable ID (e.g., &#x27;R-HSA-73817&#x27;). To find pathway IDs, use Reactome_list_top_pathways to get top-level pathways, or Reactome_map_uniprot_to_pathways if you have a UniProt protein ID. |
| `verbose` | bool | True | 运行提示开关。 |


## 列出通路包含的反应与子通路 {#dnakit.databases.reactome_get_pathway_reactions-parameters}

`dnakit.databases.reactome_get_pathway_reactions(*, stId, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `stId` | string | 必填 | Pathway Stable ID (e.g., &#x27;R-HSA-73817&#x27;). To find pathway IDs, use Reactome_list_top_pathways to get top-level pathways, or Reactome_map_uniprot_to_pathways if you have a UniProt protein ID. |
| `verbose` | bool | True | 运行提示开关。 |


## 将 UniProt 蛋白映射到 Reactome 通路 {#dnakit.databases.reactome_map_uniprot_to_pathways-parameters}

`dnakit.databases.reactome_map_uniprot_to_pathways(*, uniprot_id=None, id=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `uniprot_id` | string | None | UniProt protein accession (e.g., &#x27;P04637&#x27; for TP53, &#x27;P00533&#x27; for EGFR) |
| `id` | string | None | DNAKit optional alias for uniprot_id |
| `verbose` | bool | True | 运行提示开关。 |


## 查询 STRING 蛋白质综合功能注释 {#dnakit.databases.string_get_functional_annotations-parameters}

`dnakit.databases.string_get_functional_annotations(*, identifiers, species=None, category=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `identifiers` | string | 必填 | Protein identifier (gene name or STRING ID). Examples: &#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;EGFR&#x27;, &#x27;9606.ENSP00000269305&#x27;. |
| `species` | integer | None | NCBI taxonomy ID (default 9606 for human). Examples: 9606 (human), 10090 (mouse), 7227 (fly), 6239 (worm). {&quot;default&quot;: 9606} |
| `category` | string | None | Filter by annotation category. Options: &#x27;Process&#x27; (GO BP), &#x27;Function&#x27; (GO MF), &#x27;Component&#x27; (GO CC), &#x27;KEGG&#x27;, &#x27;RCTM&#x27; (Reactome), &#x27;WikiPathways&#x27;, &#x27;InterPro&#x27;, &#x27;Pfam&#x27;, &#x27;HPO&#x27;, &#x27;DISEASES&#x27;, &#x27;TISSUES&#x27;, &#x27;COMPARTMENTS&#x27;, &#x27;Keyword&#x27;. Leave null for all categories. |
| `verbose` | bool | True | 运行提示开关。 |


## 查询蛋白质互作伙伴与置信度 {#dnakit.databases.string_get_interaction_partners-parameters}

`dnakit.databases.string_get_interaction_partners(*, identifiers=None, species=None, limit=None, required_score=None, protein_ids=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `identifiers` | string | None | Protein identifier (gene name or STRING ID). Examples: &#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;EGFR&#x27;, &#x27;9606.ENSP00000269305&#x27; |
| `species` | integer | None | NCBI taxonomy ID (9606=human, 10090=mouse, 7227=fly, 6239=worm) {&quot;default&quot;: 9606} |
| `limit` | integer | None | Maximum number of interaction partners to return {&quot;default&quot;: 10} |
| `required_score` | integer &#124; null | None | Minimum combined STRING score (0-1000). 400=medium, 700=high, 900=highest confidence. |
| `protein_ids` | array &#124; string | None | 蛋白标识；输入列表或字符串类型见本行类型说明。 {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `verbose` | bool | True | 运行提示开关。 |


## 构建查询蛋白及其邻居的互作网络 {#dnakit.databases.string_get_network-parameters}

`dnakit.databases.string_get_network(*, identifiers=None, species=None, limit=None, required_score=None, protein_ids=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `identifiers` | string | None | Protein identifier(s). For multiple proteins, separate with &#x27;\r&#x27; (carriage return character). Examples: &#x27;TP53&#x27;, &#x27;BRCA1\rBRCA2\rTP53&#x27;, &#x27;9606.ENSP00000269305&#x27; |
| `species` | integer | None | NCBI taxonomy ID. Examples: 9606 (human), 10090 (mouse), 7227 (Drosophila), 6239 (C. elegans) {&quot;default&quot;: 9606} |
| `limit` | integer | None | Maximum number of interaction partners to include (per protein). 0 = only interactions among input proteins. {&quot;default&quot;: 10} |
| `required_score` | integer &#124; null | None | Minimum combined STRING score (0-1000). 400=medium, 700=high, 900=highest confidence. Default: 400. |
| `protein_ids` | array &#124; string | None | 蛋白标识；输入列表或字符串类型见本行类型说明。 {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `verbose` | bool | True | 运行提示开关。 |


## 检验蛋白质集合的互作富集程度 {#dnakit.databases.string_ppi_enrichment-parameters}

`dnakit.databases.string_ppi_enrichment(*, protein_ids, species=None, confidence_score=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `protein_ids` | array | 必填 | List of protein identifiers (UniProt IDs, gene names, Ensembl IDs). Minimum 3 proteins required. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}, &quot;minItems&quot;: 3} |
| `species` | integer | None | NCBI taxonomy ID (default: 9606 for human) {&quot;default&quot;: 9606} |
| `confidence_score` | number | None | Minimum confidence score for counting interactions (0-1, default: 0.4) {&quot;minimum&quot;: 0, &quot;maximum&quot;: 1, &quot;default&quot;: 0.4} |
| `verbose` | bool | True | 运行提示开关。 |


## 逐残基标注结合界面、口袋与溶剂暴露 {#dnakit.databases.structure_annotate_per_residue-parameters}

`dnakit.databases.structure_annotate_per_residue(*, operation=None, pdb_id=None, pdb_content=None, target_chain=None, partner_chains=None, ligand_resnames=None, distance_cutoff=None, core_rsa_cutoff=None, include_secondary_structure=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;annotate_per_residue&quot;&#93;, &quot;default&quot;: &quot;annotate_per_residue&quot;} |
| `pdb_id` | string &#124; null | None | 4-character PDB ID (e.g. &#x27;6VJJ&#x27;). If provided, the structure is downloaded from RCSB. Exactly one of pdb_id or pdb_content must be supplied. |
| `pdb_content` | string &#124; null | None | Raw PDB file content (string). Use when the structure is local (e.g. AlphaFold-predicted, in-memory). |
| `target_chain` | string | None | Chain ID of the protein to annotate (e.g. &#x27;A&#x27;). {&quot;default&quot;: &quot;A&quot;} |
| `partner_chains` | array | None | Chain IDs that define the binding partner(s). Used to compute the interface. Empty list = no interface analysis (all dist_partner values null). {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}, &quot;default&quot;: &#91;&#93;} |
| `ligand_resnames` | array | None | Three-letter residue names of bound ligands to define the ligand pocket (e.g. &#91;&#x27;GNP&#x27;, &#x27;MG&#x27;&#93; for GTP-analogue + magnesium). Searched across all chains. Empty list = no ligand-pocket analysis. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}, &quot;default&quot;: &#91;&#93;} |
| `distance_cutoff` | number | None | Distance cutoff in angstroms for interface and ligand-pocket classification. {&quot;default&quot;: 5.0} |
| `core_rsa_cutoff` | number | None | Relative solvent accessibility threshold below which a residue is classified as &#x27;core&#x27;. {&quot;default&quot;: 0.25} |
| `include_secondary_structure` | boolean | None | If true and pdb_id is provided, fetch per-residue secondary structure (helix/strand/coil) from PDBe REST. Only meaningful for deposited PDB entries; ignored when pdb_content is supplied. {&quot;default&quot;: false} |
| `verbose` | bool | True | 运行提示开关。 |


## 提取 UniProt 蛋白质变异与相关注释 {#dnakit.databases.uni_prot_get_disease_variants_by_accession-parameters}

`dnakit.databases.uni_prot_get_disease_variants_by_accession(*, accession, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string | 必填 | UniProtKB accession, e.g., P05067. |
| `verbose` | bool | True | 运行提示开关。 |


## 按编号读取完整 UniProt 条目 {#dnakit.databases.uni_prot_get_entry_by_accession-parameters}

`dnakit.databases.uni_prot_get_entry_by_accession(*, accession, compact=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string | 必填 | UniProtKB entry accession, e.g., P05067. |
| `compact` | boolean | None | Return a bounded summary instead of the complete UniProtKB JSON entry. Defaults to true to avoid oversized LLM outputs. Set compact=false only when you explicitly need the raw UniProtKB JSON. {&quot;default&quot;: true} |
| `verbose` | bool | True | 运行提示开关。 |


## 提取蛋白质结构域、活性位点与序列特征 {#dnakit.databases.uni_prot_get_features_by_accession-parameters}

`dnakit.databases.uni_prot_get_features_by_accession(*, accession, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string | 必填 | UniProtKB accession, e.g., &#x27;P04637&#x27; (TP53) or &#x27;P05067&#x27; (APP). |
| `verbose` | bool | True | 运行提示开关。 |


## 提取 UniProt 蛋白质功能注释 {#dnakit.databases.uni_prot_get_function_by_accession-parameters}

`dnakit.databases.uni_prot_get_function_by_accession(*, accession, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string | 必填 | UniProtKB accession, e.g., P05067. |
| `verbose` | bool | True | 运行提示开关。 |


## 提取 UniProt 标准蛋白质序列 {#dnakit.databases.uni_prot_get_sequence_by_accession-parameters}

`dnakit.databases.uni_prot_get_sequence_by_accession(*, accession, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string | 必填 | UniProtKB accession, e.g., P05067. |
| `verbose` | bool | True | 运行提示开关。 |


## 按组合检索式查询 UniProt 蛋白质 {#dnakit.databases.uni_prot_search-parameters}

`dnakit.databases.uni_prot_search(*, query, organism=None, limit=None, min_length=None, max_length=None, fields=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填 | Search query using UniProt syntax. Simple: &#x27;MEIOB&#x27;, &#x27;insulin&#x27;. Field searches: &#x27;gene:TP53&#x27;, &#x27;protein_name:insulin&#x27;, &#x27;taxonomy_id:9606&#x27;, &#x27;reviewed:true&#x27;. Ranges: &#x27;length:&#91;100 TO 500&#93;&#x27;, &#x27;mass:&#91;20000 TO 50000&#93;&#x27;. Wildcards: &#x27;gene:MEIOB*&#x27;. Boolean: &#x27;gene:TP53 AND taxonomy_id:9606&#x27;, &#x27;tissue:brain OR tissue:liver&#x27;, &#x27;reviewed:true NOT fragment:true&#x27;. Use parentheses for grouping: &#x27;(taxonomy_id:9606 OR taxonomy_id:10090) AND gene:TP53&#x27;. Taxon scoping: use &#x27;taxonomy_id:&lt;taxid&gt;&#x27; (matches the whole lineage, so strain-level entries under a species taxid are included) or &#x27;organism_name:&quot;...&quot;&#x27;. A bare &#x27;organism:&#x27; auto-converts, and &#x27;organism_id:&lt;taxid&gt;&#x27; is rewritten to &#x27;taxonomy_id:&lt;taxid&gt;&#x27; with a note in data.normalization_note, because organism_id matches only the exact organism node and misses strain-level entries. |
| `organism` | string | None | Optional organism filter. Use common names (&#x27;human&#x27;, &#x27;mouse&#x27;, &#x27;rat&#x27;, &#x27;yeast&#x27;) or taxonomy ID (&#x27;9606&#x27;). Automatically combined with query using AND. Will not duplicate if organism is already in query. |
| `limit` | integer | None | Maximum number of results to return (default: 25, max: 500). Accepts string or integer. |
| `min_length` | integer | None | Minimum sequence length. Auto-converts to an open-ended length range query (min to unbounded). |
| `max_length` | integer | None | Maximum sequence length. Auto-converts to an open-ended length range query (unbounded to max). |
| `fields` | array | None | List of field names to return (e.g., &#91;&#x27;accession&#x27;,&#x27;gene_primary&#x27;,&#x27;length&#x27;,&#x27;organism_name&#x27;&#93;). When specified, returns raw API response with requested fields. Common fields: accession, id, gene_names, gene_primary, protein_name, organism_name, organism_id, length, mass, sequence, reviewed, cc_function. See UniProt API docs for full list. Default (no fields): returns formatted response with accession, id, protein_name, gene_names, organism, length. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `verbose` | bool | True | 运行提示开关。 |


## 按通路标题检索 WikiPathways {#dnakit.databases.wiki_pathways_search-parameters}

`dnakit.databases.wiki_pathways_search(*, query, organism=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填 | Text to search, e.g., p53 |
| `organism` | string | None | Optional organism filter, e.g., &#x27;Homo sapiens&#x27; |
| `limit` | integer | None | 当前请求的返回数量上限；不是全库命中总数。 {&quot;default&quot;: 20} |
| `verbose` | bool | True | 运行提示开关。 |
