# 模式生物、物种与微生物组查询：参数与补充说明

用法和示例见[功能教程](../../api/features/22_database_query_organisms.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

## 查询果蝇基因详情 {#dnakit.databases.fly_base_get_gene-parameters}

`dnakit.databases.fly_base_get_gene(*, gene_id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填 | FlyBase gene ID with &#x27;FB:&#x27; prefix. Examples: &#x27;FB:FBgn0003996&#x27; (white/w), &#x27;FB:FBgn0000490&#x27; (decapentaplegic/dpp), &#x27;FB:FBgn0004647&#x27; (Notch/N), &#x27;FB:FBgn0000546&#x27; (eyeless/ey), &#x27;FB:FBgn0003731&#x27; (hedgehog/hh). |
| `verbose` | bool | True | 运行提示开关。 |


## 查询果蝇基因等位变异与表型 {#dnakit.databases.fly_base_get_gene_alleles-parameters}

`dnakit.databases.fly_base_get_gene_alleles(*, gene_id, limit=None, page=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填 | FlyBase gene ID with &#x27;FB:&#x27; prefix. Examples: &#x27;FB:FBgn0003996&#x27; (white/w), &#x27;FB:FBgn0000490&#x27; (dpp), &#x27;FB:FBgn0004647&#x27; (Notch/N). |
| `limit` | integer | None | Maximum number of alleles to return (1-100). Default: 20. |
| `page` | integer | None | Page number for pagination. Default: 1. |
| `verbose` | bool | True | 运行提示开关。 |


## 查询果蝇基因相关疾病模型 {#dnakit.databases.fly_base_get_gene_disease_models-parameters}

`dnakit.databases.fly_base_get_gene_disease_models(*, gene_id, limit=None, page=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填 | FlyBase gene ID with &#x27;FB:&#x27; prefix. Examples: &#x27;FB:FBgn0000490&#x27; (dpp), &#x27;FB:FBgn0004647&#x27; (Notch/N), &#x27;FB:FBgn0003996&#x27; (white/w). |
| `limit` | integer | None | Maximum number of disease models to return (1-100). Default: 20. |
| `page` | integer | None | Page number for pagination. Default: 1. |
| `verbose` | bool | True | 运行提示开关。 |


## 查询果蝇基因的发育阶段与组织表达 {#dnakit.databases.fly_base_get_gene_expression-parameters}

`dnakit.databases.fly_base_get_gene_expression(*, gene_id, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填 | FlyBase gene ID, with or without the &#x27;FB:&#x27; prefix. Examples: &#x27;FB:FBgn0000490&#x27; (dpp), &#x27;FB:FBgn0003996&#x27; (white/w), &#x27;FB:FBgn0004647&#x27; (Notch/N). |
| `limit` | integer | None | Max annotations to return (default 50, max 100). {&quot;default&quot;: 50} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询果蝇基因分子互作与遗传互作 {#dnakit.databases.fly_base_get_gene_interactions-parameters}

`dnakit.databases.fly_base_get_gene_interactions(*, gene_id, interaction_type=None, limit=None, page=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填 | FlyBase gene ID with &#x27;FB:&#x27; prefix. Examples: &#x27;FB:FBgn0003996&#x27; (white/w), &#x27;FB:FBgn0000490&#x27; (dpp), &#x27;FB:FBgn0004647&#x27; (Notch/N). |
| `interaction_type` | string | None | Type of interaction: &#x27;molecular&#x27; (physical PPI) or &#x27;genetic&#x27; (phenotype-based). Default: &#x27;molecular&#x27;. |
| `limit` | integer | None | Maximum number of interactions to return (1-100). Default: 20. |
| `page` | integer | None | Page number for pagination. Default: 1. |
| `verbose` | bool | True | 运行提示开关。 |


## 查询果蝇基因的跨物种直系同源关系 {#dnakit.databases.fly_base_get_gene_orthologs-parameters}

`dnakit.databases.fly_base_get_gene_orthologs(*, gene_id, limit=None, page=None, stringency=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填 | FlyBase gene ID with &#x27;FB:&#x27; prefix. Examples: &#x27;FB:FBgn0004647&#x27; (Notch/N), &#x27;FB:FBgn0003996&#x27; (white/w), &#x27;FB:FBgn0000546&#x27; (eyeless/ey). |
| `limit` | integer | None | Maximum number of orthologs to return (1-100). Default: 20. |
| `page` | integer | None | Page number for pagination. Default: 1. |
| `stringency` | string | None | Stringency filter for orthology calls: &#x27;stringent&#x27; (high confidence) or &#x27;moderate&#x27;. Default: &#x27;stringent&#x27;. |
| `verbose` | bool | True | 运行提示开关。 |


## 检索果蝇综合基因组学记录 {#dnakit.databases.fly_mine_search-parameters}

`dnakit.databases.fly_mine_search(*, q, size=None, format=None, facet_Category=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `q` | string | 必填 | Search query: gene symbol (e.g., &#x27;dpp&#x27;, &#x27;Notch&#x27;, &#x27;hedgehog&#x27;), protein name, domain name, or free text. |
| `size` | integer | None | Number of results to return (default: 10, max: 100). {&quot;default&quot;: 10} |
| `format` | string | None | Response format. Must be &#x27;json&#x27;. {&quot;default&quot;: &quot;json&quot;} |
| `facet_Category` | string | None | Filter by entity category. Options: &#x27;Gene&#x27;, &#x27;Protein&#x27;, &#x27;Publication&#x27;, &#x27;ProteinDomain&#x27;, &#x27;CDNAClone&#x27;, &#x27;GOTerm&#x27;, &#x27;InteractionExperiment&#x27;, &#x27;RNAiScreen&#x27;, &#x27;UniProtFeature&#x27;. |
| `verbose` | bool | True | 运行提示开关。 |


## 按分类编号查询 GBIF 物种详情 {#dnakit.databases.gbif_get_species-parameters}

`dnakit.databases.gbif_get_species(*, speciesKey, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `speciesKey` | integer | 必填 | GBIF species/taxon key. Examples: 2435099 (Puma concolor / puma), 2878688 (Quercus robur / English oak), 5219404 (Canis lupus / wolf), 2480528 (Homo sapiens), 2495662 (Aquila chrysaetos / golden eagle). Use GBIF_search_species to find keys. |
| `verbose` | bool | True | 运行提示开关。 |


## 检索 GBIF 物种名称与分类 {#dnakit.databases.gbif_search_species-parameters}

`dnakit.databases.gbif_search_species(*, query, limit=None, offset=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填 | Query keyword, e.g., Homo |
| `limit` | integer | None | 当前请求的返回数量上限；不是全库命中总数。 {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 300} |
| `offset` | integer | None | 分页起始偏移；对应当前查询结果。 {&quot;default&quot;: 0, &quot;minimum&quot;: 0} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询 GTDB 物种簇及其基因组 {#dnakit.databases.gtdb_get_species-parameters}

`dnakit.databases.gtdb_get_species(*, operation=None, species, max_genomes=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_species) {&quot;enum&quot;: &#91;&quot;get_species&quot;&#93;, &quot;default&quot;: &quot;get_species&quot;} |
| `species` | string | 必填 | Species name WITHOUT the s__ prefix (e.g., &#x27;Escherichia coli&#x27;, &#x27;Hydrogenothermus marinus&#x27;, &#x27;Bacillus subtilis&#x27;) |
| `max_genomes` | integer | None | Maximum number of genomes to return (default: 20). Species clusters can have thousands of genomes. |
| `verbose` | bool | True | 运行提示开关。 |


## 查询 GTDB 分类层级与版本历史 {#dnakit.databases.gtdb_get_taxon_info-parameters}

`dnakit.databases.gtdb_get_taxon_info(*, operation=None, taxon, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_taxon_info) {&quot;enum&quot;: &#91;&quot;get_taxon_info&quot;&#93;, &quot;default&quot;: &quot;get_taxon_info&quot;} |
| `taxon` | string | 必填 | GTDB taxon name with rank prefix (e.g., &#x27;f__Lachnospiraceae&#x27;, &#x27;d__Bacteria&#x27;, &#x27;g__Escherichia&#x27;, &#x27;s__Escherichia coli&#x27;). Use GTDB_search_taxon to find valid names. |
| `verbose` | bool | True | 运行提示开关。 |


## 按物种名称检索原核基因组 {#dnakit.databases.gtdb_search_genomes-parameters}

`dnakit.databases.gtdb_search_genomes(*, operation=None, query, page=None, items_per_page=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search_genomes) {&quot;enum&quot;: &#91;&quot;search_genomes&quot;&#93;, &quot;default&quot;: &quot;search_genomes&quot;} |
| `query` | string | 必填 | Search query for organism name (e.g., &#x27;Escherichia&#x27;, &#x27;Staphylococcus aureus&#x27;, &#x27;Lactobacillus&#x27;) |
| `page` | integer | None | Page number for pagination (default: 1) |
| `items_per_page` | integer | None | Number of results per page (default: 10, max: 50) |
| `verbose` | bool | True | 运行提示开关。 |


## 按部分名称检索原核分类单元 {#dnakit.databases.gtdb_search_taxon-parameters}

`dnakit.databases.gtdb_search_taxon(*, operation=None, query, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search_taxon) {&quot;enum&quot;: &#91;&quot;search_taxon&quot;&#93;, &quot;default&quot;: &quot;search_taxon&quot;} |
| `query` | string | 必填 | Partial taxon name to search for (e.g., &#x27;Lachnospiraceae&#x27;, &#x27;Escherichia&#x27;, &#x27;Bacillota&#x27;). Case-insensitive. |
| `limit` | integer | None | Maximum number of results per rank (default: 20, max: 100) |
| `verbose` | bool | True | 运行提示开关。 |


## 按 MGI 编号查询小鼠基因详情 {#dnakit.databases.mgi_get_gene-parameters}

`dnakit.databases.mgi_get_gene(*, gene_id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填 | MGI gene identifier with &#x27;MGI:&#x27; prefix. Examples: &#x27;MGI:98834&#x27; (Trp53/p53), &#x27;MGI:104537&#x27; (Brca1), &#x27;MGI:96677&#x27; (Kit), &#x27;MGI:88039&#x27; (Apoe), &#x27;MGI:88276&#x27; (Braf). |
| `verbose` | bool | True | 运行提示开关。 |


## 查询小鼠基因的突变表型注释 {#dnakit.databases.mgi_get_phenotypes-parameters}

`dnakit.databases.mgi_get_phenotypes(*, gene_id, limit=None, page=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填 | MGI gene identifier with &#x27;MGI:&#x27; prefix. Examples: &#x27;MGI:98834&#x27; (Trp53), &#x27;MGI:104537&#x27; (Brca1), &#x27;MGI:96677&#x27; (Kit). |
| `limit` | integer | None | Number of phenotype annotations to return per page (default: 20, max: 100). {&quot;default&quot;: 20} |
| `page` | integer | None | Page number for pagination (default: 1). {&quot;default&quot;: 1} |
| `verbose` | bool | True | 运行提示开关。 |


## 按名称或关键词检索小鼠基因 {#dnakit.databases.mgi_search_genes-parameters}

`dnakit.databases.mgi_search_genes(*, query, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填 | Gene search query: mouse gene symbol (e.g., &#x27;Trp53&#x27;, &#x27;Brca1&#x27;, &#x27;Sox9&#x27;), gene name keyword (e.g., &#x27;insulin&#x27;, &#x27;tumor&#x27;), or MGI ID (e.g., &#x27;MGI:98834&#x27;). |
| `limit` | integer | None | Maximum number of results to return (default: 10, max: 50). {&quot;default&quot;: 10} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询微生物基因组质量与分类信息 {#dnakit.databases.mgnify_get_genome-parameters}

`dnakit.databases.mgnify_get_genome(*, genome_id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `genome_id` | string | 必填 | MGnify genome accession. Examples: &#x27;MGYG000000001&#x27;, &#x27;MGYG000535630&#x27;. |
| `verbose` | bool | True | 运行提示开关。 |


## 查询宏基因组分析的 GO 功能注释 {#dnakit.databases.mgnify_get_go_terms-parameters}

`dnakit.databases.mgnify_get_go_terms(*, analysis_id, page_size=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `analysis_id` | string | 必填 | MGnify analysis accession (e.g., &#x27;MGYA00585482&#x27;). Find IDs via MGnify_list_analyses. |
| `page_size` | integer | None | Number of GO terms to return per page (default 25, max 100). {&quot;default&quot;: 25} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询宏基因组分析的蛋白结构域组成 {#dnakit.databases.mgnify_get_interpro-parameters}

`dnakit.databases.mgnify_get_interpro(*, analysis_id, page_size=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `analysis_id` | string | 必填 | MGnify analysis accession (e.g., &#x27;MGYA00585482&#x27;). Find IDs via MGnify_list_analyses. |
| `page_size` | integer | None | Number of InterPro entries to return per page (default 25, max 100). {&quot;default&quot;: 25} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询 MGnify 宏基因组研究详情 {#dnakit.databases.mgnify_get_study_detail-parameters}

`dnakit.databases.mgnify_get_study_detail(*, study_accession, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `study_accession` | string | 必填 | MGnify study accession. Examples: &#x27;MGYS00002008&#x27;, &#x27;MGYS00005292&#x27;. |
| `verbose` | bool | True | 运行提示开关。 |


## 查询宏基因组分类组成与丰度 {#dnakit.databases.mgnify_get_taxonomy-parameters}

`dnakit.databases.mgnify_get_taxonomy(*, analysis_id, rna_type=None, page_size=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `analysis_id` | string | 必填 | MGnify analysis accession (e.g., &#x27;MGYA00585482&#x27;). Find IDs via MGnify_list_analyses. |
| `rna_type` | string | None | rRNA marker type: &#x27;ssu&#x27; for 16S/18S (default) or &#x27;lsu&#x27; for 23S/28S. {&quot;default&quot;: &quot;ssu&quot;} |
| `page_size` | integer | None | Number of taxa to return per page (default 25, max 100). {&quot;default&quot;: 25} |
| `verbose` | bool | True | 运行提示开关。 |


## 列出指定宏基因组研究的分析记录 {#dnakit.databases.mgnify_list_analyses-parameters}

`dnakit.databases.mgnify_list_analyses(*, study_accession, size=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `study_accession` | string | 必填 | MGnify study accession, e.g., &#x27;MGYS00000001&#x27; |
| `size` | integer | None | 当前返回页的数量上限；不自动遍历后续页。 {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 浏览宏基因组生境分类 {#dnakit.databases.mgnify_list_biomes-parameters}

`dnakit.databases.mgnify_list_biomes(*, depth=None, page=None, page_size=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `depth` | integer | None | Hierarchy depth to filter by (1=root, 2=second level, etc.). |
| `page` | integer | None | Page number (default 1). |
| `page_size` | integer | None | Results per page (default 25, max 100). |
| `verbose` | bool | True | 运行提示开关。 |


## 检索宏基因组组装与分离株基因组 {#dnakit.databases.mgnify_search_genomes-parameters}

`dnakit.databases.mgnify_search_genomes(*, taxonomy=None, page=None, page_size=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `taxonomy` | string | None | GTDB taxonomy lineage filter, matched anywhere in the lineage. Examples: &#x27;Bacteroidota&#x27; (phylum), &#x27;Bacteroides&#x27; (genus), &#x27;Bacillota&#x27;, &#x27;Pseudomonadota&#x27;. NCBI-style phylum names such as &#x27;Bacteroidetes&#x27; or &#x27;Firmicutes&#x27; return few or no hits -- use the GTDB spelling. |
| `page` | integer | None | Page number (default 1). |
| `page_size` | integer | None | Results per page (default 25, max 100). |
| `verbose` | bool | True | 运行提示开关。 |


## 按生境与关键词检索宏基因组研究 {#dnakit.databases.mgnify_search_studies-parameters}

`dnakit.databases.mgnify_search_studies(*, biome=None, search=None, size=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `biome` | string | None | Biome lineage identifier, e.g., &#x27;root:Host-associated&#x27; or &#x27;root:Host-associated:Human:Digestive system&#x27; |
| `search` | string | None | Keyword to search in study title/description |
| `size` | integer | None | 当前返回页的数量上限；不自动遍历后续页。 {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 检索全球植物名称与分类记录 {#dnakit.databases.powo_search_plants-parameters}

`dnakit.databases.powo_search_plants(*, q, f=None, perPage=None, cursor=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `q` | string | 必填 | Search query. Can be scientific name, common name, genus, or family. Examples: &#x27;rosa canina&#x27;, &#x27;orchid&#x27;, &#x27;Asteraceae&#x27;, &#x27;oak&#x27; |
| `f` | string &#124; null | None | Filter facets (comma-separated). Common filters: &#x27;species_f:true&#x27; (species only), &#x27;accepted_names:true&#x27; (only accepted names). Example: &#x27;species_f:true&#x27; |
| `perPage` | integer &#124; null | None | Results per page (default 24, max 500) |
| `cursor` | string &#124; null | None | Pagination cursor from previous response |
| `verbose` | bool | True | 运行提示开关。 |


## 查询大鼠基因的疾病、表型与功能注释 {#dnakit.databases.rgd_get_annotations-parameters}

`dnakit.databases.rgd_get_annotations(*, rgd_id, aspect=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `rgd_id` | string | 必填 | RGD gene ID (e.g., &#x27;2218&#x27; for Brca1) |
| `aspect` | string | None | Filter by annotation type: &#x27;D&#x27; (Disease), &#x27;P&#x27; (Pathway), &#x27;W&#x27; (Phenotype), &#x27;F&#x27; (Function), &#x27;C&#x27; (Component), &#x27;E&#x27; (Expression). Leave empty for all. {&quot;default&quot;: &quot;&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 按 RGD 编号查询大鼠基因详情 {#dnakit.databases.rgd_get_gene-parameters}

`dnakit.databases.rgd_get_gene(*, rgd_id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `rgd_id` | string | 必填 | RGD gene ID (numeric, e.g., &#x27;2218&#x27; for Brca1). Can include &#x27;RGD:&#x27; prefix. |
| `verbose` | bool | True | 运行提示开关。 |


## 查询大鼠基因的跨物种直系同源关系 {#dnakit.databases.rgd_get_orthologs-parameters}

`dnakit.databases.rgd_get_orthologs(*, rgd_id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `rgd_id` | string | 必填 | RGD gene ID (e.g., &#x27;2218&#x27; for Brca1) |
| `verbose` | bool | True | 运行提示开关。 |


## 查询大鼠基因组区间的数量性状位点 {#dnakit.databases.rgd_get_qtls_in_region-parameters}

`dnakit.databases.rgd_get_qtls_in_region(*, chromosome, start, stop, map_key, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `chromosome` | string | 必填 | Chromosome name, e.g. &#x27;10&#x27; or &#x27;chr10&#x27;. |
| `start` | integer | 必填 | Region start coordinate (1-based). |
| `stop` | integer | 必填 | Region stop coordinate. |
| `map_key` | integer | 必填 | Assembly map key. 360=rat rn7/GRCr8, 372=rat mRatBN7.2, 38=human GRCh38, 17=human GRCh37, 35=mouse GRCm39, 18=mouse GRCm38. |
| `verbose` | bool | True | 运行提示开关。 |


## 按品系名称查询大鼠品系详情 {#dnakit.databases.rgd_get_strain-parameters}

`dnakit.databases.rgd_get_strain(*, symbol, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `symbol` | string | 必填 | Exact rat strain symbol (case-insensitive), e.g. &#x27;SHR&#x27; (spontaneously hypertensive rat), &#x27;BN&#x27; (Brown Norway), &#x27;SS&#x27; (Dahl salt-sensitive). |
| `verbose` | bool | True | 运行提示开关。 |


## 查询大鼠品系的疾病与表型注释 {#dnakit.databases.rgd_get_strain_annotations-parameters}

`dnakit.databases.rgd_get_strain_annotations(*, symbol, category=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `symbol` | string | 必填 | Exact rat strain symbol (case-insensitive), e.g. &#x27;SHR&#x27;, &#x27;BN&#x27;, &#x27;SS&#x27;. |
| `category` | string &#124; null | None | Optional annotation category filter: &#x27;disease&#x27;, &#x27;phenotype&#x27;, &#x27;strain&#x27;, &#x27;pathway&#x27;, etc. Omit to return all categories. |
| `verbose` | bool | True | 运行提示开关。 |


## 按名称或关键词检索 RGD 基因 {#dnakit.databases.rgd_search_genes-parameters}

`dnakit.databases.rgd_search_genes(*, query, species=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填 | Gene symbol or keyword (e.g., &#x27;Brca1&#x27;, &#x27;Tp53&#x27;, &#x27;insulin&#x27;) |
| `species` | string | None | Species to search: &#x27;rat&#x27; (default), &#x27;human&#x27;, &#x27;mouse&#x27; {&quot;default&quot;: &quot;rat&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 按关键词与类型检索大鼠品系 {#dnakit.databases.rgd_search_strains-parameters}

`dnakit.databases.rgd_search_strains(*, query=None, strain_type=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string &#124; null | None | Keyword to match against strain symbol, name, genetics, and description (case-insensitive substring). Example: &#x27;hypertensive&#x27; or &#x27;Dahl&#x27;. |
| `strain_type` | string &#124; null | None | Filter by strain type, e.g. &#x27;inbred&#x27;, &#x27;congenic&#x27;, &#x27;transgenic&#x27;, &#x27;recombinant_inbred&#x27;, &#x27;consomic&#x27;. |
| `limit` | integer &#124; null | None | Maximum number of strains to return (default 25, max 200). |
| `verbose` | bool | True | 运行提示开关。 |


## 查询酵母基因与基因座概况 {#dnakit.databases.sgd_get_gene-parameters}

`dnakit.databases.sgd_get_gene(*, sgd_id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `sgd_id` | string | 必填 | SGD locus identifier. Examples: &#x27;S000000259&#x27; (GAL4), &#x27;S000003219&#x27; (RMR1), &#x27;S000005564&#x27; (ACT1). |
| `verbose` | bool | True | 运行提示开关。 |


## 查询酵母基因的 GO 注释与证据 {#dnakit.databases.sgd_get_go_annotations-parameters}

`dnakit.databases.sgd_get_go_annotations(*, sgd_id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `sgd_id` | string | 必填 | SGD locus identifier. Examples: &#x27;S000003219&#x27; (RMR1), &#x27;S000005564&#x27; (ACT1), &#x27;S000000259&#x27; (GAL4). |
| `verbose` | bool | True | 运行提示开关。 |


## 查询酵母基因的遗传与物理互作 {#dnakit.databases.sgd_get_interactions-parameters}

`dnakit.databases.sgd_get_interactions(*, sgd_id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `sgd_id` | string | 必填 | SGD locus identifier. Examples: &#x27;S000003219&#x27; (RMR1), &#x27;S000005564&#x27; (ACT1), &#x27;S000000259&#x27; (GAL4). |
| `verbose` | bool | True | 运行提示开关。 |


## 查询酵母突变表型与实验条件 {#dnakit.databases.sgd_get_phenotypes-parameters}

`dnakit.databases.sgd_get_phenotypes(*, sgd_id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `sgd_id` | string | 必填 | SGD locus identifier. Examples: &#x27;S000003219&#x27; (RMR1), &#x27;S000000259&#x27; (GAL4), &#x27;S000005564&#x27; (ACT1). |
| `verbose` | bool | True | 运行提示开关。 |


## 跨数据类型检索酵母数据库 {#dnakit.databases.sgd_search-parameters}

`dnakit.databases.sgd_search(*, query, category=None, limit=None, offset=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填 | Search query string. Examples: &#x27;GAL4&#x27;, &#x27;cell cycle&#x27;, &#x27;sporulation&#x27;. |
| `category` | string | None | Filter by category. Options: &#x27;locus&#x27;, &#x27;reference&#x27;, &#x27;biological_process&#x27;, &#x27;molecular_function&#x27;, &#x27;cellular_component&#x27;, &#x27;phenotype&#x27;, &#x27;disease&#x27;, &#x27;allele&#x27;, &#x27;complex&#x27;, &#x27;pathway&#x27;, &#x27;chemical&#x27;, &#x27;strain&#x27;, &#x27;download&#x27;. |
| `limit` | integer | None | Maximum results to return (default 10, max 50). |
| `offset` | integer | None | Pagination offset (default 0). |
| `verbose` | bool | True | 运行提示开关。 |


## 查询植物形态、生长与繁殖特征 {#dnakit.databases.usda_plants_get_characteristics-parameters}

`dnakit.databases.usda_plants_get_characteristics(*, symbol, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `symbol` | string | 必填 | USDA PLANTS symbol, e.g. &#x27;ABBA&#x27;, &#x27;ACRU&#x27;. |
| `verbose` | bool | True | 运行提示开关。 |


## 查询植物名称、分类与生活型 {#dnakit.databases.usda_plants_get_profile-parameters}

`dnakit.databases.usda_plants_get_profile(*, symbol, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `symbol` | string | 必填 | USDA PLANTS symbol, e.g. &#x27;ABBA&#x27; (balsam fir), &#x27;ACRU&#x27; (red maple), &#x27;QUAL&#x27; (white oak). |
| `verbose` | bool | True | 运行提示开关。 |


## 查询线虫基因的组织与发育表达 {#dnakit.databases.worm_base_get_expression-parameters}

`dnakit.databases.worm_base_get_expression(*, gene_id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填 | WormBase gene identifier. Examples: &#x27;WBGene00006763&#x27; (unc-26), &#x27;WBGene00000898&#x27; (daf-2). |
| `verbose` | bool | True | 运行提示开关。 |


## 查询秀丽隐杆线虫基因详情 {#dnakit.databases.worm_base_get_gene-parameters}

`dnakit.databases.worm_base_get_gene(*, gene_id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填 | WormBase gene identifier (WBGene format). Examples: &#x27;WBGene00006763&#x27; (unc-26), &#x27;WBGene00006789&#x27; (unc-52), &#x27;WBGene00000898&#x27; (daf-2). |
| `verbose` | bool | True | 运行提示开关。 |


## 查询线虫基因的实验表型 {#dnakit.databases.worm_base_get_phenotypes-parameters}

`dnakit.databases.worm_base_get_phenotypes(*, gene_id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填 | WormBase gene identifier. Examples: &#x27;WBGene00006763&#x27; (unc-26), &#x27;WBGene00000898&#x27; (daf-2). |
| `verbose` | bool | True | 运行提示开关。 |


## 查询爪蟾基因详情与跨库标识 {#dnakit.databases.xenbase_get_gene-parameters}

`dnakit.databases.xenbase_get_gene(*, gene_id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填 | Xenbase gene ID. Can be with or without &#x27;Xenbase:&#x27; prefix. Examples: &#x27;XB-GENE-484286&#x27; (tp53, X. tropicalis), &#x27;XB-GENE-6252168&#x27; (tp53.L, X. laevis), &#x27;XB-GENE-487192&#x27; (pax6, X. tropicalis). |
| `verbose` | bool | True | 运行提示开关。 |


## 检索热带爪蟾与非洲爪蟾基因 {#dnakit.databases.xenbase_search_genes-parameters}

`dnakit.databases.xenbase_search_genes(*, query, species=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填 | Gene symbol or name to search. Examples: &#x27;tp53&#x27;, &#x27;pax6&#x27;, &#x27;sox2&#x27;, &#x27;bmp4&#x27;, &#x27;wnt3a&#x27;, &#x27;fgf8&#x27;. |
| `species` | string | None | Optional species filter: &#x27;tropicalis&#x27; for X. tropicalis, &#x27;laevis&#x27; for X. laevis. Leave empty for both species. {&quot;default&quot;: &quot;&quot;} |
| `limit` | integer | None | Maximum number of results to return. Default: 10. {&quot;default&quot;: 10} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询斑马鱼基因详情 {#dnakit.databases.zfin_get_gene-parameters}

`dnakit.databases.zfin_get_gene(*, gene_id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填 | ZFIN gene ID with &#x27;ZFIN:&#x27; prefix. Examples: &#x27;ZFIN:ZDB-GENE-990415-8&#x27; (pax2a), &#x27;ZFIN:ZDB-GENE-980526-166&#x27; (shha/sonic hedgehog a), &#x27;ZFIN:ZDB-GENE-990415-72&#x27; (tp53/tumor protein p53), &#x27;ZFIN:ZDB-GENE-030131-7740&#x27; (bmp4/bone morphogenetic protein 4). |
| `verbose` | bool | True | 运行提示开关。 |


## 查询斑马鱼基因的发育与组织表达 {#dnakit.databases.zfin_get_gene_expression-parameters}

`dnakit.databases.zfin_get_gene_expression(*, gene_id, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填 | ZFIN gene ID, with or without the &#x27;ZFIN:&#x27; prefix. Examples: &#x27;ZFIN:ZDB-GENE-990415-8&#x27; (pax2a), &#x27;ZFIN:ZDB-GENE-980526-166&#x27; (shha), &#x27;ZFIN:ZDB-GENE-990415-72&#x27; (tp53). |
| `limit` | integer | None | Max annotations to return (default 50, max 100). {&quot;default&quot;: 50} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询斑马鱼基因突变表型 {#dnakit.databases.zfin_get_gene_phenotypes-parameters}

`dnakit.databases.zfin_get_gene_phenotypes(*, gene_id, limit=None, page=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_id` | string | 必填 | ZFIN gene ID with &#x27;ZFIN:&#x27; prefix. Examples: &#x27;ZFIN:ZDB-GENE-990415-8&#x27; (pax2a), &#x27;ZFIN:ZDB-GENE-990415-72&#x27; (tp53), &#x27;ZFIN:ZDB-GENE-980526-166&#x27; (shha). |
| `limit` | integer | None | Maximum number of phenotypes to return (1-100). Default: 20. |
| `page` | integer | None | Page number for pagination. Default: 1. |
| `verbose` | bool | True | 运行提示开关。 |
