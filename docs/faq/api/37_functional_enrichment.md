# 基因集与通路富集：参数与补充说明

用法和示例见[功能教程](../../api/features/37_functional_enrichment.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

## 对基因集合执行指定库的富集分析 {#dnakit.databases.enrichr_enrich-parameters}

`dnakit.databases.enrichr_enrich(*, operation=None, gene_list, library=None, top_n=None, species=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: enrich) {&quot;const&quot;: &quot;enrich&quot;, &quot;default&quot;: &quot;enrich&quot;} |
| `gene_list` | array | 必填 | List of gene symbols (e.g., &#91;&#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;EGFR&#x27;&#93;). Minimum 1 gene. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `library` | string | None | Enrichr library name (default: GO_Biological_Process_2023). Common: GO_Biological_Process_2023, KEGG_2021_Human, Reactome_2022, WikiPathways_2024_Human, MSigDB_Hallmark_2020 {&quot;default&quot;: &quot;GO_Biological_Process_2023&quot;} |
| `top_n` | integer | None | Number of top enriched terms to return (default: 10) {&quot;default&quot;: 10} |
| `species` | string &#124; null | None | Organism instance to query: &#x27;human&#x27; (default), &#x27;mouse&#x27;, &#x27;fly&#x27; (FlyEnrichr), &#x27;worm&#x27; (WormEnrichr), &#x27;yeast&#x27; (YeastEnrichr), or &#x27;fish&#x27; (FishEnrichr). Gene set libraries differ per organism, and gene symbols must match that organism&#x27;s nomenclature. {&quot;enum&quot;: &#91;null, &quot;human&quot;, &quot;mouse&quot;, &quot;fly&quot;, &quot;worm&quot;, &quot;yeast&quot;, &quot;fish&quot;&#93;} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询基因与富集通路的连接关系 {#dnakit.databases.enrichr_gene_enrichment_analysis-parameters}

`dnakit.databases.enrichr_gene_enrichment_analysis(*, gene_list, libs=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_list` | array | 必填 | List of gene names or symbols to analyze. At least 2 genes are required for path ranking analysis. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}, &quot;minItems&quot;: 2} |
| `libs` | array | None | List of Enrichr enrichment libraries to query (e.g. &#x27;KEGG_2021_Human&#x27;, &#x27;GO_Biological_Process_2023&#x27;). Omitting this or passing an empty list uses a broad default set of pathway/ontology libraries. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}, &quot;default&quot;: &#91;&quot;WikiPathways_2024_Human&quot;, &quot;Reactome_Pathways_2024&quot;, &quot;MSigDB_Hallmark_2020&quot;, &quot;GO_Molecular_Function_2023&quot;, &quot;GO_Biological_Process_2023&quot;&#93;} |
| `verbose` | bool | True | 运行提示开关。 |


## 执行 PANTHER 基因集合过度代表分析 {#dnakit.databases.panther_enrichment-parameters}

`dnakit.databases.panther_enrichment(*, gene_list, organism=None, annotation_dataset=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_list` | string | 必填 | Comma-separated list of gene symbols, UniProt IDs, or Ensembl IDs. Example: &#x27;TP53,BRCA1,EGFR,KRAS,PTEN,APC,RB1&#x27;. |
| `organism` | integer &#124; null | None | NCBI taxonomy ID. Default: 9606 (human). Others: 10090 (mouse), 7227 (fly), 559292 (yeast). |
| `annotation_dataset` | string &#124; null | None | Annotation dataset for enrichment. Options: &#x27;GO:0008150&#x27; (biological_process, default), &#x27;GO:0003674&#x27; (molecular_function), &#x27;GO:0005575&#x27; (cellular_component), &#x27;ANNOT_TYPE_ID_PANTHER_PATHWAY&#x27; (PANTHER pathway), &#x27;ANNOT_TYPE_ID_REACTOME_PATHWAY&#x27; (Reactome pathway), &#x27;ANNOT_TYPE_ID_PANTHER_PC&#x27; (protein class). |
| `verbose` | bool | True | 运行提示开关。 |


## 执行 Reactome 通路富集分析 {#dnakit.databases.reactome_analysis_pathway_enrichment-parameters}

`dnakit.databases.reactome_analysis_pathway_enrichment(*, identifiers, page_size=None, include_disease=None, projection=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `identifiers` | string | 必填 | Newline-separated list of gene/protein identifiers. Supports gene symbols (TP53, BRCA1), UniProt IDs (P04637, P38398), or Ensembl IDs. Example: &#x27;TP53\nBRCA1\nCDH1\nEGFR\nKRAS&#x27;. |
| `page_size` | integer &#124; null | None | Number of pathways to return (default 20, max 50). |
| `include_disease` | boolean &#124; null | None | Include disease pathways in results (default true). |
| `projection` | boolean &#124; null | None | Project identifiers to human Reactome pathways for cross-species analysis (default true). |
| `verbose` | bool | True | 运行提示开关。 |


## 分析蛋白质集合的功能与通路富集 {#dnakit.databases.string_functional_enrichment-parameters}

`dnakit.databases.string_functional_enrichment(*, protein_ids, species=None, category=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `protein_ids` | array | 必填 | List of protein identifiers (UniProt IDs, gene names, Ensembl IDs). Minimum 3 proteins recommended for meaningful enrichment. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}, &quot;minItems&quot;: 2} |
| `species` | integer | None | NCBI taxonomy ID (default: 9606 for human) {&quot;default&quot;: 9606} |
| `category` | string | None | Enrichment category: &#x27;Process&#x27; (GO Biological Process), &#x27;Component&#x27; (GO Cellular Component), &#x27;Function&#x27; (GO Molecular Function), &#x27;KEGG&#x27;, &#x27;Reactome&#x27;, &#x27;WikiPathways&#x27;, &#x27;COMPARTMENTS&#x27;, &#x27;TISSUES&#x27;, &#x27;DISEASES&#x27; {&quot;enum&quot;: &#91;&quot;Process&quot;, &quot;Component&quot;, &quot;Function&quot;, &quot;KEGG&quot;, &quot;Reactome&quot;, &quot;WikiPathways&quot;, &quot;COMPARTMENTS&quot;, &quot;TISSUES&quot;, &quot;DISEASES&quot;&#93;, &quot;default&quot;: &quot;Process&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 使用多行蛋白标识执行 STRING 富集分析 {#dnakit.databases.string_get_functional_enrichment-parameters}

`dnakit.databases.string_get_functional_enrichment(*, identifiers, species=9606, category=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `identifiers` | string | 必填 | 以换行或回车分隔的蛋白标识；将转换成 STRING 标识列表。 |
| `species` | integer | 9606 | 物种名或 Taxonomy ID；类型与支持值取决于当前接口。 {&quot;default&quot;: 9606} |
| `category` | string | None | 分类或结果类别，按当前接口的枚举/来源定义填写。 |
| `verbose` | bool | True | 运行提示开关。 |
