# 表达、调控与肿瘤数据查询：参数与补充说明

用法和示例见[功能教程](../../api/features/22_database_query_regulation.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

ChIP-Atlas 按 gene 查峰时，未提供 peak_file 会返回 requires_data 和检索计划；提供已核对组装的 BED 后才计算区间重叠。富集入口返回 requires_submission 提交说明，不能把它当成已计算的富集统计。

## 分页查询 cBioPortal 肿瘤研究目录 {#dnakit.databases.c_bio_portal_get_cancer_studies-parameters}

`dnakit.databases.c_bio_portal_get_cancer_studies(*, limit=None, offset=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `limit` | integer | None | Number of studies to return in this page. The default of 20 covers only a small fraction of the catalogue; set it to the `total_available` reported in the response to retrieve every study. {&quot;default&quot;: 20} |
| `offset` | integer | None | 0-based index of the first study to return, for paging through the catalogue. Applied client-side because the cBioPortal /studies endpoint ignores its own pageNumber parameter. {&quot;default&quot;: 0} |
| `verbose` | bool | True | 运行提示开关。 |


## 按 Entrez 编号查询肿瘤数据库基因信息 {#dnakit.databases.c_bio_portal_get_gene_info-parameters}

`dnakit.databases.c_bio_portal_get_gene_info(*, entrez_gene_id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `entrez_gene_id` | integer | 必填 | Entrez Gene ID (e.g., 672 for BRCA1) |
| `verbose` | bool | True | 运行提示开关。 |


## 查询肿瘤研究的分子数据类型 {#dnakit.databases.c_bio_portal_get_molecular_profiles-parameters}

`dnakit.databases.c_bio_portal_get_molecular_profiles(*, study_id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `study_id` | string | 必填 | Cancer study ID |
| `verbose` | bool | True | 运行提示开关。 |


## 查询肿瘤研究中指定基因的突变 {#dnakit.databases.c_bio_portal_get_mutations-parameters}

`dnakit.databases.c_bio_portal_get_mutations(*, study_id, gene_list, sample_list_id=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `study_id` | string | 必填 | Cancer study ID (e.g., &#x27;brca_tcga&#x27;) |
| `gene_list` | string | 必填 | Comma-separated gene symbols (e.g., &#x27;BRCA1,BRCA2&#x27;) |
| `sample_list_id` | string | None | Optional sample list ID. If not provided, uses all samples in the study. |
| `verbose` | bool | True | 运行提示开关。 |


## 查询肿瘤研究的样本集合 {#dnakit.databases.c_bio_portal_get_samples-parameters}

`dnakit.databases.c_bio_portal_get_samples(*, study_id, page_size=None, page_number=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `study_id` | string | 必填 | Cancer study ID |
| `page_size` | integer | None | Maximum number of samples to return. Defaults high enough to return every sample in a study; a smaller value silently truncates the cohort, which will corrupt any frequency denominator computed from it. {&quot;default&quot;: 100000} |
| `page_number` | integer | None | 0-based page index, used with `page_size` to page through very large studies. {&quot;default&quot;: 0} |
| `verbose` | bool | True | 运行提示开关。 |


## 获取 CELLxGENE 原始 H5AD 数据文件 {#dnakit.databases.cellx_gene_download_h5ad-parameters}

`dnakit.databases.cellx_gene_download_h5ad(*, operation=None, dataset_id, output_path=None, census_version=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;download_h5ad&quot;&#93;, &quot;default&quot;: &quot;download_h5ad&quot;} |
| `dataset_id` | string | 必填 | CELLxGENE dataset identifier (required) |
| `output_path` | string | None | Local path to save H5AD file (omit to get URI only) |
| `census_version` | string | None | Census version to query. &#x27;stable&#x27; (recommended, Long-Term Support release), &#x27;latest&#x27; (newest data, may change), or specific date &#x27;YYYY-MM-DD&#x27; (for reproducibility). Default &#x27;stable&#x27; is best for production analyses. {&quot;default&quot;: &quot;stable&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

运行前安装 `cellxgene-census`。


## 按筛选条件查询单细胞元数据 {#dnakit.databases.cellx_gene_get_cell_metadata-parameters}

`dnakit.databases.cellx_gene_get_cell_metadata(*, operation=None, organism=None, obs_value_filter=None, column_names=None, census_version=None, gene=None, tissue=None, max_cells=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;get_obs_metadata&quot;&#93;, &quot;default&quot;: &quot;get_obs_metadata&quot;} |
| `organism` | string | None | Organism name {&quot;enum&quot;: &#91;&quot;Homo sapiens&quot;, &quot;Mus musculus&quot;&#93;, &quot;default&quot;: &quot;Homo sapiens&quot;} |
| `obs_value_filter` | string | None | REQUIRED - filter cells using SQL-like syntax. Unfiltered queries timeout (50M+ cells). Format: &#x27;field == &quot;value&quot;&#x27;. Operators: ==, !=, in, &lt;, &gt;, &lt;=, &gt;=. Combine with &#x27;and&#x27;/&#x27;or&#x27;. Common fields and values: tissue_general (brain, heart, lung, liver, blood, kidney, skin, eye, intestine, pancreas, bone marrow), cell_type (T cell, B cell, neuron, macrophage, fibroblast, epithelial cell, endothelial cell), disease (normal, COVID-19, Alzheimer disease, lung adenocarcinoma), assay (10x 3&#x27; v3, Smart-seq2), sex (male, female). Examples: &#x27;tissue_general == &quot;lung&quot;&#x27;, &#x27;disease == &quot;COVID-19&quot; and tissue_general == &quot;blood&quot;&#x27;, &#x27;cell_type in &#91;&quot;T cell&quot;, &quot;B cell&quot;&#93; and disease == &quot;normal&quot;&#x27;. |
| `column_names` | array | None | Specific columns to return (default: all columns) {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `census_version` | string | None | Census version to query. &#x27;stable&#x27; (recommended, Long-Term Support release), &#x27;latest&#x27; (newest data, may change), or specific date &#x27;YYYY-MM-DD&#x27; (for reproducibility). Default &#x27;stable&#x27; is best for production analyses. {&quot;default&quot;: &quot;stable&quot;} |
| `gene` | string | None | 基因标识或基因列表；命名要求和类型见当前接口。 |
| `tissue` | string | None | 目标组织名称或标识，按当前数据源填写。 |
| `max_cells` | integer | None | 本次返回或处理的细胞数上限。 {&quot;minimum&quot;: 1, &quot;maximum&quot;: 10000} |
| `verbose` | bool | True | 运行提示开关。 |

运行前安装 `cellxgene-census`。


## 获取预计算的单细胞嵌入表示 {#dnakit.databases.cellx_gene_get_embeddings-parameters}

`dnakit.databases.cellx_gene_get_embeddings(*, operation=None, organism=None, embedding_name=None, census_version=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;get_embeddings&quot;&#93;, &quot;default&quot;: &quot;get_embeddings&quot;} |
| `organism` | string | None | Organism name {&quot;enum&quot;: &#91;&quot;Homo sapiens&quot;, &quot;Mus musculus&quot;&#93;, &quot;default&quot;: &quot;Homo sapiens&quot;} |
| `embedding_name` | string | None | Name of specific embedding to retrieve (omit to list available embeddings) |
| `census_version` | string | None | Census version to query. &#x27;stable&#x27; (recommended, Long-Term Support release), &#x27;latest&#x27; (newest data, may change), or specific date &#x27;YYYY-MM-DD&#x27; (for reproducibility). Default &#x27;stable&#x27; is best for production analyses. {&quot;default&quot;: &quot;stable&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

运行前安装 `cellxgene-census`。


## 按细胞与基因条件查询表达矩阵 {#dnakit.databases.cellx_gene_get_expression_data-parameters}

`dnakit.databases.cellx_gene_get_expression_data(*, operation=None, organism=None, obs_value_filter=None, var_value_filter=None, obs_column_names=None, var_column_names=None, census_version=None, gene=None, tissue=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;get_anndata&quot;&#93;, &quot;default&quot;: &quot;get_anndata&quot;} |
| `organism` | string | None | Organism name {&quot;enum&quot;: &#91;&quot;Homo sapiens&quot;, &quot;Mus musculus&quot;&#93;, &quot;default&quot;: &quot;Homo sapiens&quot;} |
| `obs_value_filter` | string | None | REQUIRED (or use var_value_filter) - Cell filter. Common values: tissue_general (brain, lung, heart, blood, liver), cell_type (T cell, B cell, neuron, macrophage), disease (normal, COVID-19, Alzheimer disease). Example: &#x27;tissue_general == &quot;lung&quot; and disease == &quot;normal&quot;&#x27; |
| `var_value_filter` | string | None | REQUIRED (or use obs_value_filter) - Gene filter by symbol or Ensembl ID. Examples: &#x27;feature_name in &#91;&quot;TP53&quot;, &quot;BRCA1&quot;, &quot;EGFR&quot;&#93;&#x27;, &#x27;feature_name == &quot;CD4&quot;&#x27;. Common genes: TP53, EGFR, BRCA1, CD4, CD8A, IL6, TNF, GAPDH. |
| `obs_column_names` | array | None | Cell metadata columns to include {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `var_column_names` | array | None | Gene metadata columns to include {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `census_version` | string | None | Census version to query. &#x27;stable&#x27; (recommended, Long-Term Support release), &#x27;latest&#x27; (newest data, may change), or specific date &#x27;YYYY-MM-DD&#x27; (for reproducibility). Default &#x27;stable&#x27; is best for production analyses. {&quot;default&quot;: &quot;stable&quot;} |
| `gene` | string | None | 基因标识或基因列表；命名要求和类型见当前接口。 |
| `tissue` | string | None | 目标组织名称或标识，按当前数据源填写。 |
| `verbose` | bool | True | 运行提示开关。 |

运行前安装 `cellxgene-census`。


## 生成 ChIP-Atlas 富集任务提交说明 {#dnakit.databases.ch_ipatlas_enrichment_analysis-parameters}

`dnakit.databases.ch_ipatlas_enrichment_analysis(*, operation=None, bed_data=None, motif=None, gene_list=None, genome=None, antigen_class=None, cell_type_class=None, threshold=None, distance=None, gene=None, cell_type=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | 提供方操作模式；默认值与允许值见本行约束。 {&quot;enum&quot;: &#91;&quot;enrichment_analysis&quot;&#93;, &quot;default&quot;: &quot;enrichment_analysis&quot;} |
| `bed_data` | string | None | **Option 1**: BED format genomic regions (tab-separated: chr, start, end). For finding proteins bound to specific genomic regions. Example: &#x27;chr1\t1000\t2000\nchr2\t3000\t4000&#x27;. |
| `motif` | string | None | **Option 2**: DNA sequence motif in IUPAC notation. Use: A/T/G/C (bases), W=A&#124;T, S=G&#124;C, M=A&#124;C, K=G&#124;T, R=A&#124;G, Y=C&#124;T, B=C&#124;G&#124;T, D=A&#124;G&#124;T, H=A&#124;C&#124;T, V=A&#124;C&#124;G, N=any. For finding proteins bound to specific DNA sequences. Example: &#x27;CANNTG&#x27; (E-box motif). |
| `gene_list` | array &#124; string | None | **Option 3**: Gene symbols (HGNC for human, MGI for mouse, RGD for rat, FlyBase, WormBase, SGD for yeast). Provide as array or single gene. For finding transcription factors regulating genes. Example: &#91;&#x27;TP53&#x27;, &#x27;MDM2&#x27;, &#x27;CDKN1A&#x27;&#93;. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `genome` | string | None | Genome assembly {&quot;enum&quot;: &#91;&quot;hg38&quot;, &quot;hg19&quot;, &quot;mm10&quot;, &quot;mm9&quot;, &quot;rn6&quot;, &quot;dm6&quot;, &quot;dm3&quot;, &quot;ce11&quot;, &quot;ce10&quot;, &quot;sacCer3&quot;&#93;, &quot;default&quot;: &quot;hg38&quot;} |
| `antigen_class` | string | None | Filter by antigen class (e.g., &#x27;TFs and others&#x27;, &#x27;Histone&#x27;, &#x27;RNA polymerase&#x27;) |
| `cell_type_class` | string | None | Filter by cell type class (e.g., &#x27;Blood&#x27;, &#x27;Liver&#x27;, &#x27;Brain&#x27;) |
| `threshold` | string | None | Peak calling stringency (MACS2 Q-value). Options: &#x27;05&#x27;=1e-5 (permissive, more peaks, broader features), &#x27;10&#x27;=1e-10 (moderate, balanced), &#x27;20&#x27;=1e-20 (strict, high confidence only, narrow peaks). Default &#x27;05&#x27; suitable for most analyses. Higher values = fewer but more confident peaks. {&quot;enum&quot;: &#91;&quot;05&quot;, &quot;10&quot;, &quot;20&quot;&#93;, &quot;default&quot;: &quot;05&quot;} |
| `distance` | string | None | Distance from Transcription Start Site (TSS) in base pairs for gene-TF association. Defines promoter region. Default 5000 (±5kb, captures typical promoters). Use 1000-2000 for narrow promoters, 10000+ for enhancer regions. {&quot;default&quot;: &quot;5000&quot;} |
| `gene` | array &#124; string | None | DNAKit optional alias for gene_list {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `cell_type` | string | None | DNAKit optional alias for cell_type_class |
| `verbose` | bool | True | 运行提示开关。 |


## 检索 ChIP-Atlas 实验元数据 {#dnakit.databases.ch_ipatlas_get_experiments-parameters}

`dnakit.databases.ch_ipatlas_get_experiments(*, operation=None, genome=None, antigen=None, cell_type=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | 提供方操作模式；默认值与允许值见本行约束。 {&quot;enum&quot;: &#91;&quot;get_experiment_list&quot;&#93;, &quot;default&quot;: &quot;get_experiment_list&quot;} |
| `genome` | string | None | Filter by genome assembly {&quot;enum&quot;: &#91;&quot;hg38&quot;, &quot;hg19&quot;, &quot;mm10&quot;, &quot;mm9&quot;, &quot;rn6&quot;, &quot;dm6&quot;, &quot;dm3&quot;, &quot;ce11&quot;, &quot;ce10&quot;, &quot;sacCer3&quot;&#93;} |
| `antigen` | string | None | Filter by antigen/protein name (e.g., &#x27;CTCF&#x27;, &#x27;H3K4me3&#x27;) |
| `cell_type` | string | None | Filter by cell type (e.g., &#x27;K-562&#x27;, &#x27;HepG2&#x27;, &#x27;CD4&#x27;) |
| `limit` | integer | None | Maximum results to return {&quot;default&quot;: 100, &quot;minimum&quot;: 1, &quot;maximum&quot;: 1000} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询 ChIP-Atlas 峰资源与本地区间重叠 {#dnakit.databases.ch_ipatlas_get_peak_data-parameters}

`dnakit.databases.ch_ipatlas_get_peak_data(*, operation=None, experiment_id=None, genome=None, format=None, threshold=None, gene=None, experiment_type=None, limit=None, peak_file=None, flank_bp=None, max_file_bytes=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | 提供方操作模式；默认值与允许值见本行约束。 {&quot;enum&quot;: &#91;&quot;get_peak_data&quot;&#93;, &quot;default&quot;: &quot;get_peak_data&quot;} |
| `experiment_id` | string | None | Experiment ID (SRX/ERX/DRX format, required) |
| `genome` | string | None | Genome assembly {&quot;enum&quot;: &#91;&quot;hg38&quot;, &quot;hg19&quot;, &quot;mm10&quot;, &quot;mm9&quot;, &quot;rn6&quot;, &quot;dm6&quot;, &quot;dm3&quot;, &quot;ce11&quot;, &quot;ce10&quot;, &quot;sacCer3&quot;&#93;, &quot;default&quot;: &quot;hg38&quot;} |
| `format` | string | None | Output format {&quot;enum&quot;: &#91;&quot;bigwig&quot;, &quot;bed&quot;, &quot;bigbed&quot;&#93;, &quot;default&quot;: &quot;bigwig&quot;} |
| `threshold` | string | None | Q-value threshold for BED/BigBed peak files. &#x27;05&#x27;=1e-5 (more peaks), &#x27;10&#x27;=1e-10 (moderate), &#x27;20&#x27;=1e-20 (high confidence only). Only applies to BED/BigBed formats. Default &#x27;05&#x27;. {&quot;enum&quot;: &#91;&quot;05&quot;, &quot;10&quot;, &quot;20&quot;&#93;, &quot;default&quot;: &quot;05&quot;} |
| `gene` | string | None | 基因标识或基因列表；命名要求和类型见当前接口。 |
| `experiment_type` | string | None | 实验类别，按数据源接受的类型筛选。 |
| `limit` | integer | None | 当前请求的返回数量上限；不是全库命中总数。 {&quot;minimum&quot;: 1, &quot;maximum&quot;: 10000} |
| `peak_file` | string | None | 本地 BED 峰文件；必须核对其基因组版本和轨道来源。 |
| `flank_bp` | integer | None | 基因区域两侧扩展的碱基数。 {&quot;minimum&quot;: 0, &quot;maximum&quot;: 1000000} |
| `max_file_bytes` | integer | None | 允许读取的峰文件字节上限。 {&quot;minimum&quot;: 1} |
| `verbose` | bool | True | 运行提示开关。 |


## 按抗原或细胞类型检索 ChIP-Atlas 数据集 {#dnakit.databases.ch_ipatlas_search_datasets-parameters}

`dnakit.databases.ch_ipatlas_search_datasets(*, operation=None, antigen=None, cell_type=None, genome=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | 提供方操作模式；默认值与允许值见本行约束。 {&quot;enum&quot;: &#91;&quot;search_datasets&quot;&#93;, &quot;default&quot;: &quot;search_datasets&quot;} |
| `antigen` | string | None | Search by antigen/protein name (e.g., &#x27;CTCF&#x27;, &#x27;H3K27ac&#x27;) |
| `cell_type` | string | None | Search by cell type (e.g., &#x27;K-562&#x27;, &#x27;HeLa&#x27;, &#x27;MCF-7&#x27;) |
| `genome` | string | None | Genome assembly {&quot;enum&quot;: &#91;&quot;hg38&quot;, &quot;hg19&quot;, &quot;mm10&quot;, &quot;mm9&quot;, &quot;rn6&quot;, &quot;dm6&quot;, &quot;dm3&quot;, &quot;ce11&quot;, &quot;ce10&quot;, &quot;sacCer3&quot;&#93;, &quot;default&quot;: &quot;hg38&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询 CRISPR 筛选的基因依赖信息 {#dnakit.databases.dep_map_get_gene_dependencies-parameters}

`dnakit.databases.dep_map_get_gene_dependencies(*, gene_symbol, model_id=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 必填 | Gene symbol (e.g., &#x27;EGFR&#x27;, &#x27;KRAS&#x27;, &#x27;TP53&#x27;) |
| `model_id` | string | None | Optional: Filter by specific cell line |
| `verbose` | bool | True | 运行提示开关。 |


## 按名称检索肿瘤细胞系 {#dnakit.databases.dep_map_search_cell_lines-parameters}

`dnakit.databases.dep_map_search_cell_lines(*, query, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填 | Search query (cell line name like &#x27;A549&#x27;, &#x27;MCF7&#x27;) |
| `verbose` | bool | True | 运行提示开关。 |


## 查询 eQTL 数据集中的变异表达关联 {#dnakit.databases.e_qtl_get_associations-parameters}

`dnakit.databases.e_qtl_get_associations(*, dataset_id, size=None, gene_id=None, variant=None, molecular_trait_id=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `dataset_id` | string | 必填 | Dataset identifier from eQTL_list_datasets (e.g., &#x27;QTD000001&#x27;, &#x27;QTD000584&#x27;) |
| `size` | integer | None | Number of associations to return (default 20) {&quot;default&quot;: 20} |
| `gene_id` | string &#124; null | None | Filter by Ensembl gene ID (e.g., &#x27;ENSG00000141510&#x27; for TP53). Returns all eQTLs affecting this gene. |
| `variant` | string &#124; null | None | Filter by variant in chr_pos_ref_alt format (e.g., &#x27;chr1_791100_G_GGGA&#x27;). Returns all genes affected by this variant. |
| `molecular_trait_id` | string &#124; null | None | Filter by molecular trait ID (usually same as gene_id for &#x27;ge&#x27; datasets) |
| `verbose` | bool | True | 运行提示开关。 |


## 列出可用的 eQTL 研究与组织数据集 {#dnakit.databases.e_qtl_list_datasets-parameters}

`dnakit.databases.e_qtl_list_datasets(*, size=None, quant_method=None, tissue_id=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `size` | integer | None | Number of datasets to return (default 20) {&quot;default&quot;: 20} |
| `quant_method` | string &#124; null | None | Filter by quantification method: &#x27;ge&#x27; (gene expression), &#x27;exon&#x27;, &#x27;tx&#x27; (transcript), &#x27;txrev&#x27; (transcript ratio), &#x27;leafcutter&#x27; (splicing) |
| `tissue_id` | string &#124; null | None | Filter by tissue ontology ID (e.g., &#x27;CL_0000235&#x27; for macrophage, &#x27;UBERON_0002107&#x27; for liver) |
| `verbose` | bool | True | 运行提示开关。 |


## 查询 ENCODE 生物样本完整元数据 {#dnakit.databases.encode_get_biosample-parameters}

`dnakit.databases.encode_get_biosample(*, accession, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string | 必填 | ENCODE biosample accession identifier (format: ENCBS######, e.g., &#x27;ENCBS000AAA&#x27;, &#x27;ENCBS130ENC&#x27;). Find accessions using ENCODE_search_biosamples or from experiment metadata. |
| `verbose` | bool | True | 运行提示开关。 |


## 检索染色质状态分割注释 {#dnakit.databases.encode_get_chromatin_state-parameters}

`dnakit.databases.encode_get_chromatin_state(*, biosample_term_name=None, organism=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `biosample_term_name` | string &#124; null | None | Biosample name (e.g., &#x27;K562&#x27;, &#x27;HepG2&#x27;, &#x27;GM12878&#x27;, &#x27;liver&#x27;). Leave empty for all. |
| `organism` | string | None | Organism scientific name. {&quot;default&quot;: &quot;Homo sapiens&quot;} |
| `limit` | integer | None | Maximum number of results. {&quot;default&quot;: 25, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询 ENCODE 实验元数据与分析信息 {#dnakit.databases.encode_get_experiment-parameters}

`dnakit.databases.encode_get_experiment(*, accession, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string | 必填 | ENCODE experiment accession identifier (format: ENCSR######, e.g., &#x27;ENCSR000AKS&#x27;, &#x27;ENCSR000CAG&#x27;). Find accessions using ENCODE_search_experiments. |
| `verbose` | bool | True | 运行提示开关。 |


## 按格式与实验类型列出 ENCODE 文件 {#dnakit.databases.encode_list_files-parameters}

`dnakit.databases.encode_list_files(*, file_type=None, assay_title=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `file_type` | string | None | File format filter (e.g., &#x27;fastq&#x27;, &#x27;bam&#x27;, &#x27;bigWig&#x27;, &#x27;bed&#x27;). Common formats: fastq (raw reads), bam (alignments), bigWig (coverage), bed (peaks/regions). |
| `assay_title` | string | None | Assay filter (e.g., &#x27;ChIP-seq&#x27;, &#x27;RNA-seq&#x27;). Filters files by the experimental assay that generated them. |
| `limit` | integer | None | Maximum number of results to return (1–100). {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 检索顺式调控元件与表观组注释 {#dnakit.databases.encode_search_annotations-parameters}

`dnakit.databases.encode_search_annotations(*, annotation_type=None, biosample_term_name=None, organism=None, assembly=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `annotation_type` | string | None | Annotation type filter. Options: &#x27;candidate Cis-Regulatory Elements&#x27; (cCREs), &#x27;chromatin state&#x27; (ChromHMM), &#x27;imputed signal&#x27; (predicted epigenetic marks). Default: cCREs. {&quot;default&quot;: &quot;candidate Cis-Regulatory Elements&quot;} |
| `biosample_term_name` | string &#124; null | None | Biosample filter (e.g., &#x27;K562&#x27;, &#x27;HepG2&#x27;). Leave empty for all. |
| `organism` | string | None | Organism scientific name. {&quot;default&quot;: &quot;Homo sapiens&quot;} |
| `assembly` | string | None | Genome assembly (e.g., &#x27;GRCh38&#x27;, &#x27;hg19&#x27;, &#x27;mm10&#x27;). {&quot;default&quot;: &quot;GRCh38&quot;} |
| `limit` | integer | None | Maximum number of results. {&quot;default&quot;: 25, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 按物种与处理条件检索 ENCODE 样本 {#dnakit.databases.encode_search_biosamples-parameters}

`dnakit.databases.encode_search_biosamples(*, organism=None, biosample_type=None, treatment=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `organism` | string | None | Organism filter (e.g., &#x27;Homo sapiens&#x27;, &#x27;Mus musculus&#x27;). |
| `biosample_type` | string | None | Biosample classification filter — maps to ENCODE&#x27;s biosample_ontology.classification field (e.g., &#x27;cell line&#x27;, &#x27;tissue&#x27;, &#x27;primary cell&#x27;, &#x27;in vitro differentiated cells&#x27;). |
| `treatment` | string | None | Treatment filter (e.g., &#x27;interferon gamma&#x27;, &#x27;ethanol&#x27;). Use &#x27;*&#x27; to find any treated samples. |
| `limit` | integer | None | Maximum number of results to return (1–100). {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 检索 ATAC-seq 与 DNase-seq 实验 {#dnakit.databases.encode_search_chromatin_accessibility-parameters}

`dnakit.databases.encode_search_chromatin_accessibility(*, assay_type=None, biosample_term_name=None, organism=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `assay_type` | string | None | Chromatin accessibility assay: &#x27;ATAC-seq&#x27; (Assay for Transposase-Accessible Chromatin, newer/preferred) or &#x27;DNase-seq&#x27; (DNase I hypersensitive sites, larger dataset). {&quot;default&quot;: &quot;ATAC-seq&quot;} |
| `biosample_term_name` | string &#124; null | None | Biosample name filter (e.g., &#x27;K562&#x27;, &#x27;GM12878&#x27;, &#x27;liver&#x27;). Leave empty for all. |
| `organism` | string | None | Organism scientific name. {&quot;default&quot;: &quot;Homo sapiens&quot;} |
| `limit` | integer | None | Maximum number of results to return. {&quot;default&quot;: 25, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 通过数据库模块检索 ENCODE 实验 {#dnakit.databases.encode_search_experiments-parameters}

`dnakit.databases.encode_search_experiments(*, assay_title=None, target=None, organism=None, status=None, limit=None, biosample=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `assay_title` | string | None | 实验方法名称，例如 ATAC-seq。 |
| `target` | string | None | 实验靶标或目标标识，按当前数据源填写。 |
| `organism` | string | None | 目标生物名称，按提供方命名规则填写。 |
| `status` | string | None | 记录的发布或处理状态筛选条件。 {&quot;default&quot;: &quot;released&quot;} |
| `limit` | integer | None | 当前请求的返回数量上限；不是全库命中总数。 {&quot;default&quot;: 10} |
| `biosample` | string | None | 生物样本名称或检索条件。 |
| `verbose` | bool | True | 运行提示开关。 |


## 按组蛋白修饰检索 ChIP-seq 实验 {#dnakit.databases.encode_search_histone_experiments-parameters}

`dnakit.databases.encode_search_histone_experiments(*, histone_mark=None, target=None, biosample_term_name=None, biosample=None, organism=None, limit=None, biosample_term=None, cell_type=None, tissue=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `histone_mark` | string &#124; null | None | Histone modification mark to filter by (e.g., &#x27;H3K4me3&#x27;, &#x27;H3K27ac&#x27;, &#x27;H3K27me3&#x27;, &#x27;H3K4me1&#x27;, &#x27;H3K36me3&#x27;, &#x27;H3K9me3&#x27;). Leave empty to search all histone marks. |
| `target` | string &#124; null | None | Alias for histone_mark. Histone modification mark (e.g., &#x27;H3K27ac&#x27;, &#x27;H3K4me3&#x27;). |
| `biosample_term_name` | string &#124; null | None | Biosample name from ENCODE ontology (cell lines or tissues, NOT disease names). Examples: &#x27;K562&#x27;, &#x27;HepG2&#x27;, &#x27;GM12878&#x27;, &#x27;liver&#x27;, &#x27;brain&#x27;. For disease-specific data, search by relevant cell lines (e.g., K562/MOLM-13 for AML). Leave empty to search all. |
| `biosample` | string &#124; null | None | Alias for biosample_term_name. Biosample (tissue or cell line, e.g., &#x27;liver&#x27;, &#x27;K562&#x27;). |
| `organism` | string | None | Organism scientific name (e.g., &#x27;Homo sapiens&#x27;, &#x27;Mus musculus&#x27;). {&quot;default&quot;: &quot;Homo sapiens&quot;} |
| `limit` | integer | None | Maximum number of results to return (1-100). {&quot;default&quot;: 25, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `biosample_term` | string &#124; null | None | Alias for biosample_term_name. Biosample tissue or cell line (e.g., &quot;breast epithelium&quot;, &quot;liver&quot;, &quot;K562&quot;). |
| `cell_type` | string &#124; null | None | Alias for biosample_term_name. Cell type or tissue (e.g., &#x27;GM12878&#x27;, &#x27;K562&#x27;, &#x27;liver&#x27;, &#x27;brain&#x27;). |
| `tissue` | string &#124; null | None | Alias for biosample_term_name. Tissue or cell line (e.g., &#x27;liver&#x27;, &#x27;brain&#x27;, &#x27;breast epithelium&#x27;, &#x27;MCF-7&#x27;). |
| `verbose` | bool | True | 运行提示开关。 |


## 按样本与实验类型检索 ENCODE RNA-seq {#dnakit.databases.encode_search_rnaseq_experiments-parameters}

`dnakit.databases.encode_search_rnaseq_experiments(*, biosample_term_name=None, biosample=None, cell_type=None, tissue=None, organism=None, assay_type=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `biosample_term_name` | string &#124; null | None | Biosample name from ENCODE ontology (e.g., &#x27;K562&#x27;, &#x27;HepG2&#x27;, &#x27;GM12878&#x27;, &#x27;liver&#x27;, &#x27;brain&#x27;). Leave empty to search all biosamples. |
| `biosample` | string &#124; null | None | Alias for biosample_term_name. Cell type or tissue (e.g., &#x27;K562&#x27;, &#x27;liver&#x27;). |
| `cell_type` | string &#124; null | None | Alias for biosample_term_name. Cell type (e.g., &#x27;K562&#x27;, &#x27;HepG2&#x27;, &#x27;GM12878&#x27;). |
| `tissue` | string &#124; null | None | Alias for biosample_term_name. Tissue type (e.g., &#x27;liver&#x27;, &#x27;brain&#x27;, &#x27;heart&#x27;). |
| `organism` | string | None | Organism scientific name (e.g., &#x27;Homo sapiens&#x27;, &#x27;Mus musculus&#x27;). {&quot;default&quot;: &quot;Homo sapiens&quot;} |
| `assay_type` | string | None | RNA-seq assay type: &#x27;total RNA-seq&#x27; (default, all RNA), &#x27;polyA plus RNA-seq&#x27; (mRNA), &#x27;small RNA-seq&#x27;, &#x27;microRNA-seq&#x27;. Aliases: &#x27;total&#x27;, &#x27;polya&#x27;, &#x27;mirna&#x27;, &#x27;small&#x27;. {&quot;default&quot;: &quot;total RNA-seq&quot;} |
| `limit` | integer | None | Maximum number of results to return (1-100). {&quot;default&quot;: 25, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 检索差异表达实验的项目记录 {#dnakit.databases.expression_atlas_search_differential-parameters}

`dnakit.databases.expression_atlas_search_differential(*, gene=None, condition=None, species=None, limit=None, offset=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | None | Gene symbol or Ensembl ID (e.g., &#x27;TP53&#x27;, &#x27;ENSG00000141510&#x27;) |
| `condition` | string | None | Condition/disease to filter by (e.g., &#x27;cancer&#x27;, &#x27;inflammation&#x27;, &#x27;breast&#x27;) |
| `species` | string | None | Species name (default: &#x27;homo sapiens&#x27;) {&quot;default&quot;: &quot;homo sapiens&quot;} |
| `limit` | integer &#124; null | None | Maximum number of experiments to return (default 50, max 500). Compare with the total reported in the response to see how many exist. {&quot;default&quot;: 50} |
| `offset` | integer &#124; null | None | Zero-based index of the first experiment to return, for paging through results beyond the first page (e.g. offset=50 returns 51-100). Default: 0. {&quot;default&quot;: 0} |
| `verbose` | bool | True | 运行提示开关。 |


## 检索基线与差异表达实验 {#dnakit.databases.expression_atlas_search_experiments-parameters}

`dnakit.databases.expression_atlas_search_experiments(*, gene=None, condition=None, species=None, limit=None, offset=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | None | Gene symbol or Ensembl ID |
| `condition` | string | None | Biological condition or tissue (e.g., &#x27;brain&#x27;, &#x27;cancer&#x27;, &#x27;kidney&#x27;) |
| `species` | string | None | Species name (optional, e.g., &#x27;homo sapiens&#x27;) |
| `limit` | integer &#124; null | None | Maximum number of experiments to return (default 50, max 500). Compare with the total reported in the response to see how many exist. {&quot;default&quot;: 50} |
| `offset` | integer &#124; null | None | Zero-based index of the first experiment to return, for paging through results beyond the first page (e.g. offset=50 returns 51-100). Default: 0. {&quot;default&quot;: 0} |
| `verbose` | bool | True | 运行提示开关。 |


## 检索 Hi-C 与染色质构象数据 {#dnakit.databases.four_dn_search_data-parameters}

`dnakit.databases.four_dn_search_data(*, operation=None, query=None, item_type=None, file_type=None, assay_title=None, biosource_name=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | 提供方操作模式；默认值与允许值见本行约束。 {&quot;enum&quot;: &#91;&quot;search&quot;&#93;, &quot;default&quot;: &quot;search&quot;} |
| `query` | string | None | Search query (default: &#x27;*&#x27; for all) {&quot;default&quot;: &quot;*&quot;} |
| `item_type` | string | None | Type of item to search {&quot;enum&quot;: &#91;&quot;File&quot;, &quot;Experiment&quot;, &quot;Biosource&quot;&#93;, &quot;default&quot;: &quot;File&quot;} |
| `file_type` | string | None | Filter by file type. Common types: &#x27;contact list&#x27; (processed contact matrices), &#x27;pairs&#x27; (aligned read pairs), &#x27;cooler&#x27; (multi-resolution contact matrices), &#x27;mcool&#x27; (multi-resolution cooler), &#x27;hic&#x27; (Juicer format). |
| `assay_title` | string | None | Filter by assay (e.g., &#x27;Hi-C&#x27;, &#x27;Micro-C&#x27;, &#x27;Capture Hi-C&#x27;) |
| `biosource_name` | string | None | Filter by cell type (e.g., &#x27;H1-hESC&#x27;, &#x27;HFFc6&#x27;, &#x27;GM12878&#x27;) |
| `limit` | integer | None | Maximum results to return {&quot;default&quot;: 25, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询 GDC 肿瘤病例的临床资料 {#dnakit.databases.gdc_get_clinical_data-parameters}

`dnakit.databases.gdc_get_clinical_data(*, project_id=None, primary_site=None, disease_type=None, vital_status=None, gender=None, size=None, offset=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `project_id` | string | None | GDC project identifier (e.g., &#x27;TCGA-BRCA&#x27;, &#x27;TCGA-LUAD&#x27;, &#x27;TARGET-AML&#x27;) |
| `primary_site` | string | None | Primary anatomical site (e.g., &#x27;Breast&#x27;, &#x27;Lung&#x27;, &#x27;Brain&#x27;) |
| `disease_type` | string | None | Disease type filter (e.g., &#x27;Ductal and Lobular Neoplasms&#x27;) |
| `vital_status` | string | None | Vital status filter: &#x27;Alive&#x27; or &#x27;Dead&#x27; {&quot;enum&quot;: &#91;&quot;Alive&quot;, &quot;Dead&quot;&#93;} |
| `gender` | string | None | Gender filter: &#x27;female&#x27; or &#x27;male&#x27; {&quot;enum&quot;: &#91;&quot;female&quot;, &quot;male&quot;&#93;} |
| `size` | integer | None | Number of cases to return (1-100) {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `offset` | integer | None | Pagination offset (0-based) {&quot;default&quot;: 0, &quot;minimum&quot;: 0} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询基因在泛癌队列中的突变频率 {#dnakit.databases.gdc_get_mutation_frequency-parameters}

`dnakit.databases.gdc_get_mutation_frequency(*, gene_symbol=None, gene=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_symbol` | string | None | Gene symbol (e.g., &#x27;TP53&#x27;, &#x27;KRAS&#x27;) |
| `gene` | string | None | Gene symbol alias — alternative to gene_symbol |
| `verbose` | bool | True | 运行提示开关。 |


## 按基因查询 GDC 体细胞突变 {#dnakit.databases.gdc_get_ssm_by_gene-parameters}

`dnakit.databases.gdc_get_ssm_by_gene(*, gene_symbol, project_id=None, size=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 必填 | Gene symbol (e.g., &#x27;TP53&#x27;, &#x27;EGFR&#x27;, &#x27;BRAF&#x27;) |
| `project_id` | string | None | Optional: Filter by project (e.g., &#x27;TCGA-BRCA&#x27;) |
| `size` | integer | None | Number of results (1–100) {&quot;default&quot;: 20, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询肿瘤队列的生存与删失数据 {#dnakit.databases.gdc_get_survival-parameters}

`dnakit.databases.gdc_get_survival(*, project_id, gene_symbol=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `project_id` | string | 必填 | GDC project identifier (e.g., &#x27;TCGA-BRCA&#x27;, &#x27;TCGA-LUAD&#x27;, &#x27;TCGA-GBM&#x27;) |
| `gene_symbol` | string | None | Optional: gene symbol to filter cases with mutations in this gene (e.g., &#x27;TP53&#x27;, &#x27;KRAS&#x27;) |
| `verbose` | bool | True | 运行提示开关。 |


## 按数据类型列出 GDC 文件 {#dnakit.databases.gdc_list_files-parameters}

`dnakit.databases.gdc_list_files(*, data_type=None, size=None, offset=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `data_type` | string | None | Data type filter (e.g., &#x27;Gene Expression Quantification&#x27;) |
| `size` | integer | None | Number of results (1–100) {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `offset` | integer | None | Offset for pagination (0-based) {&quot;default&quot;: 0, &quot;minimum&quot;: 0} |
| `verbose` | bool | True | 运行提示开关。 |


## 列出 GDC 肿瘤项目与汇总信息 {#dnakit.databases.gdc_list_projects-parameters}

`dnakit.databases.gdc_list_projects(*, program=None, size=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `program` | string | None | Filter by program (e.g., &#x27;TCGA&#x27;, &#x27;TARGET&#x27;) |
| `size` | integer | None | Number of results (1–100) {&quot;default&quot;: 20, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 按筛选条件检索 GDC 肿瘤病例 {#dnakit.databases.gdc_search_cases-parameters}

`dnakit.databases.gdc_search_cases(*, project_id=None, size=None, offset=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `project_id` | string | None | GDC project identifier (e.g., &#x27;TCGA-BRCA&#x27;) |
| `size` | integer | None | Number of results (1–100) {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `offset` | integer | None | Offset for pagination (0-based) {&quot;default&quot;: 0, &quot;minimum&quot;: 0} |
| `verbose` | bool | True | 运行提示开关。 |


## 检索 GEO 染色质可及性数据集 {#dnakit.databases.geo_search_atacseq_datasets-parameters}

`dnakit.databases.geo_search_atacseq_datasets(*, query=None, organism=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string &#124; null | None | Search query (e.g., &#x27;T cell differentiation&#x27;, &#x27;liver hepatocyte&#x27;, &#x27;hematopoiesis&#x27;). Leave empty to search broadly. |
| `organism` | string | None | Organism scientific name (e.g., &#x27;Homo sapiens&#x27;, &#x27;Mus musculus&#x27;). {&quot;default&quot;: &quot;Homo sapiens&quot;} |
| `limit` | integer | None | Maximum number of results to return (1-100). {&quot;default&quot;: 20, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 检索 GEO 蛋白质 DNA 结合数据集 {#dnakit.databases.geo_search_chipseq_datasets-parameters}

`dnakit.databases.geo_search_chipseq_datasets(*, query, organism=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填 | Search terms for ChIP-seq datasets (e.g., &#x27;H3K27ac liver&#x27;, &#x27;CTCF cancer&#x27;, &#x27;p53 ChIP-seq&#x27;). |
| `organism` | string | None | Organism filter. {&quot;default&quot;: &quot;Homo sapiens&quot;} |
| `limit` | integer | None | Maximum number of dataset IDs to return. {&quot;default&quot;: 20, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 检索 GEO 转录组表达数据集 {#dnakit.databases.geo_search_rnaseq_datasets-parameters}

`dnakit.databases.geo_search_rnaseq_datasets(*, query=None, organism=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string &#124; null | None | Search query (e.g., &#x27;BRCA1 breast cancer&#x27;, &#x27;T cell activation&#x27;, &#x27;Alzheimer disease&#x27;). Leave empty to search broadly. |
| `organism` | string | None | Organism scientific name (e.g., &#x27;Homo sapiens&#x27;, &#x27;Mus musculus&#x27;). {&quot;default&quot;: &quot;Homo sapiens&quot;} |
| `limit` | integer | None | Maximum number of results to return (1-100). {&quot;default&quot;: 20, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 计算指定组织中基因与变异的 eQTL 关联 {#dnakit.databases.gtex_calculate_eqtl-parameters}

`dnakit.databases.gtex_calculate_eqtl(*, operation=None, gencode_id, variant_id, tissue_site_detail_id, dataset_id=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;calculate_eqtl&quot;&#93;, &quot;default&quot;: &quot;calculate_eqtl&quot;} |
| `gencode_id` | string | 必填 | Required: Versioned GENCODE ID (e.g. ENSG00000141510.18) |
| `variant_id` | string | 必填 | Required: GTEx variant ID |
| `tissue_site_detail_id` | string | 必填 | Required: Tissue ID (e.g. &#x27;Whole_Blood&#x27;) |
| `dataset_id` | string | None | 提供方的数据集标识。 {&quot;enum&quot;: &#91;&quot;gtex_v8&quot;, &quot;gtex_v10&quot;, &quot;gtex_snrnaseq_pilot&quot;&#93;, &quot;default&quot;: &quot;gtex_v8&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询基因的 GTEx 表达摘要 {#dnakit.databases.gtex_get_expression_summary-parameters}

`dnakit.databases.gtex_get_expression_summary(*, gene_symbol=None, ensembl_gene_id=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_symbol` | string | None | Gene symbol (e.g., TP53, BRCA1). Auto-resolved to GENCODE ID. |
| `ensembl_gene_id` | string | None | Ensembl gene ID, e.g., ENSG00000141510 |
| `verbose` | bool | True | 运行提示开关。 |


## 查询基因在不同组织中的表达中位数 {#dnakit.databases.gtex_get_median_gene_expression-parameters}

`dnakit.databases.gtex_get_median_gene_expression(*, operation=None, gencode_id=None, tissue_site_detail_id=None, dataset_id=None, page=None, items_per_page=None, gene_symbol=None, gene=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;get_median_gene_expression&quot;&#93;, &quot;default&quot;: &quot;get_median_gene_expression&quot;} |
| `gencode_id` | string &#124; array | None | Gene identifier(s): gene symbol (e.g. &#x27;TP53&#x27;), unversioned Ensembl ID (e.g. &#x27;ENSG00000141510&#x27;), or versioned GENCODE ID (e.g. &#x27;ENSG00000141510.18&#x27;). Auto-resolved to versioned GENCODE ID. Can be single string or array. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `tissue_site_detail_id` | array | None | Optional: Tissue IDs to filter (e.g. &#91;&#x27;Liver&#x27;, &#x27;Brain_Cortex&#x27;&#93;). Omit for all tissues. See GTEx_get_tissue_sites for valid IDs {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `dataset_id` | string | None | GTEx dataset version (default: gtex_v8, the most widely cited release). gtex_v10 is the newer release (943 eQTL donors vs 838) and is fully queryable -- gene IDs are resolved to each dataset&#x27;s own GENCODE version automatically. {&quot;enum&quot;: &#91;&quot;gtex_v8&quot;, &quot;gtex_v10&quot;, &quot;gtex_snrnaseq_pilot&quot;&#93;, &quot;default&quot;: &quot;gtex_v8&quot;} |
| `page` | integer | None | Page number for pagination (0-based) {&quot;minimum&quot;: 0, &quot;default&quot;: 0} |
| `items_per_page` | integer | None | Results per page {&quot;minimum&quot;: 1, &quot;maximum&quot;: 100000, &quot;default&quot;: 250} |
| `gene_symbol` | string | None | Gene symbol alias for gencode_id (e.g., &quot;TP53&quot;, &quot;COL5A1&quot;) |
| `gene` | string | None | DNAKit optional alias for gene_symbol |
| `verbose` | bool | True | 运行提示开关。 |


## 查询跨组织 eQTL 荟萃分析结果 {#dnakit.databases.gtex_get_multi_tissue_eqtls-parameters}

`dnakit.databases.gtex_get_multi_tissue_eqtls(*, operation=None, gencode_id, variant_id=None, dataset_id=None, page=None, items_per_page=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;get_multi_tissue_eqtls&quot;&#93;, &quot;default&quot;: &quot;get_multi_tissue_eqtls&quot;} |
| `gencode_id` | string | 必填 | Required: Versioned GENCODE ID |
| `variant_id` | string | None | Optional: GTEx variant ID to filter specific variant |
| `dataset_id` | string | None | 提供方的数据集标识。 {&quot;enum&quot;: &#91;&quot;gtex_v8&quot;, &quot;gtex_v10&quot;, &quot;gtex_snrnaseq_pilot&quot;&#93;, &quot;default&quot;: &quot;gtex_v8&quot;} |
| `page` | integer | None | 结果页码；起始页约定按提供方参数说明。 {&quot;minimum&quot;: 0, &quot;default&quot;: 0} |
| `items_per_page` | integer | None | 每页返回的记录数。 {&quot;minimum&quot;: 1, &quot;maximum&quot;: 100000, &quot;default&quot;: 250} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询显著的单组织 eQTL 关联 {#dnakit.databases.gtex_get_single_tissue_eqtls-parameters}

`dnakit.databases.gtex_get_single_tissue_eqtls(*, operation=None, gencode_id=None, variant_id=None, tissue_site_detail_id=None, dataset_id=None, page=None, items_per_page=None, gene_symbol=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;get_single_tissue_eqtls&quot;&#93;, &quot;default&quot;: &quot;get_single_tissue_eqtls&quot;} |
| `gencode_id` | array | None | Optional: Versioned GENCODE ID(s) to query {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `variant_id` | array | None | Optional: GTEx variant ID(s) to query {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `tissue_site_detail_id` | array | None | Optional: Tissue ID(s) to filter. At least one of gencode_id, variant_id, or tissue_site_detail_id required {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `dataset_id` | string | None | 提供方的数据集标识。 {&quot;enum&quot;: &#91;&quot;gtex_v8&quot;, &quot;gtex_v10&quot;, &quot;gtex_snrnaseq_pilot&quot;&#93;, &quot;default&quot;: &quot;gtex_v8&quot;} |
| `page` | integer | None | 结果页码；起始页约定按提供方参数说明。 {&quot;minimum&quot;: 0, &quot;default&quot;: 0} |
| `items_per_page` | integer | None | 每页返回的记录数。 {&quot;minimum&quot;: 1, &quot;maximum&quot;: 100000, &quot;default&quot;: 250} |
| `gene_symbol` | string | None | 基因符号，例如 BRCA1。 |
| `verbose` | bool | True | 运行提示开关。 |


## 列出 GTEx 组织类型与样本统计 {#dnakit.databases.gtex_get_tissue_sites-parameters}

`dnakit.databases.gtex_get_tissue_sites(*, operation=None, dataset_id=None, page=None, items_per_page=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;get_tissue_sites&quot;&#93;, &quot;default&quot;: &quot;get_tissue_sites&quot;} |
| `dataset_id` | string | None | GTEx dataset version {&quot;enum&quot;: &#91;&quot;gtex_v8&quot;, &quot;gtex_v10&quot;, &quot;gtex_snrnaseq_pilot&quot;&#93;, &quot;default&quot;: &quot;gtex_v8&quot;} |
| `page` | integer | None | Page number {&quot;minimum&quot;: 0, &quot;default&quot;: 0} |
| `items_per_page` | integer | None | Results per page {&quot;minimum&quot;: 1, &quot;maximum&quot;: 100000, &quot;default&quot;: 250} |
| `verbose` | bool | True | 运行提示开关。 |


## 按条件查询 GTEx 单组织 eQTL {#dnakit.databases.gtex_query_eqtl-parameters}

`dnakit.databases.gtex_query_eqtl(*, gene_symbol=None, ensembl_gene_id=None, page=None, size=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_symbol` | string | None | Gene symbol (e.g., TP53, BRCA1). Auto-resolved to GENCODE ID. |
| `ensembl_gene_id` | string | None | Ensembl gene ID, e.g., ENSG00000141510 |
| `page` | integer | None | Page number (1-based) {&quot;default&quot;: 1, &quot;minimum&quot;: 1} |
| `size` | integer | None | Page size (1–100) {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 汇总多组学数据中的基因信息 {#dnakit.databases.harmonizome_get_gene-parameters}

`dnakit.databases.harmonizome_get_gene(*, gene_symbol, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 必填 | Official gene symbol (case-insensitive). Examples: &#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;EGFR&#x27;, &#x27;KRAS&#x27;, &#x27;CDK2&#x27;. |
| `verbose` | bool | True | 运行提示开关。 |


## 查询基因表达与肿瘤生存的关联 {#dnakit.databases.hpa_get_cancer_prognostics_by_gene-parameters}

`dnakit.databases.hpa_get_cancer_prognostics_by_gene(*, ensembl_id=None, gene_name=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `ensembl_id` | string | None | Ensembl Gene ID of the gene to check, e.g., &#x27;ENSG00000141510&#x27; for TP53, &#x27;ENSG00000012048&#x27; for BRCA1. |
| `gene_name` | string | None | DNAKit optional alias for ensembl_id |
| `verbose` | bool | True | 运行提示开关。 |


## 查看已停用的 HPA 蛋白互作查询 {#dnakit.databases.hpa_get_protein_interactions_by_gene-parameters}

`dnakit.databases.hpa_get_protein_interactions_by_gene(*, gene_name, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_name` | string | 必填 | Official gene symbol, e.g., &#x27;EGFR&#x27;, &#x27;TP53&#x27;, &#x27;BRCA1&#x27;, etc. |
| `verbose` | bool | True | 运行提示开关。 |


## 按组织或细胞来源查询 RNA 表达 {#dnakit.databases.hpa_get_rna_expression_by_source-parameters}

`dnakit.databases.hpa_get_rna_expression_by_source(*, gene_name, source_type, source_name, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_name` | string | 必填 | Gene name or gene symbol, e.g., &#x27;GFAP&#x27;, &#x27;TP53&#x27;, &#x27;BRCA1&#x27;, etc. |
| `source_type` | string | 必填 | The type of biological source. Choose from: &#x27;tissue&#x27; (51 consensus tissues, nTPM), &#x27;blood&#x27; (19 immune cell subsets, nTPM), &#x27;brain&#x27; (13 brain regions, nTPM), &#x27;single_cell&#x27; (~160 single cell types, nCPM). |
| `source_name` | string | 必填 | The specific source name, e.g. &#x27;liver&#x27;, &#x27;heart_muscle&#x27;, &#x27;skin&#x27;, &#x27;t_cell&#x27;, &#x27;hepatocyte&#x27;, &#x27;cerebellum&#x27;. Case-insensitive; spaces and hyphens are treated as underscores. Names HPA publishes only as subsets (blood &#x27;t_cell&#x27;, &#x27;b_cell&#x27;, &#x27;monocyte&#x27;, &#x27;dendritic_cell&#x27;) resolve to a representative subset column and the response reports which one in a &#x27;note&#x27;. An invalid name returns an error listing valid names -- it never returns a silent &#x27;N/A&#x27;. |
| `verbose` | bool | True | 运行提示开关。 |


## 查询基因在指定组织中的 RNA 表达 {#dnakit.databases.hpa_get_rna_expression_in_specific_tissues-parameters}

`dnakit.databases.hpa_get_rna_expression_in_specific_tissues(*, ensembl_id=None, tissue_names=None, tissues=None, gene_name=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `ensembl_id` | string | None | Ensembl Gene ID for the gene, e.g., &#x27;ENSG00000141510&#x27; for TP53. |
| `tissue_names` | array | None | List of tissue names to query, e.g., &#91;&#x27;skin&#x27;, &#x27;liver&#x27;, &#x27;heart muscle&#x27;, &#x27;kidney&#x27;&#93;. Case-insensitive; underscores, spaces and hyphens are interchangeable. HPA has no RNA column for &#x27;bronchus&#x27;, &#x27;nasopharynx&#x27;, &#x27;oral mucosa&#x27; or &#x27;soft tissue&#x27; (protein/IHC only). {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `tissues` | array | None | DNAKit optional alias for tissue_names {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `gene_name` | string | None | DNAKit optional alias for ensembl_id |
| `verbose` | bool | True | 运行提示开关。 |


## 查询蛋白质亚细胞定位注释 {#dnakit.databases.hpa_get_subcellular_location-parameters}

`dnakit.databases.hpa_get_subcellular_location(*, gene_name, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_name` | string | 必填 | Gene name or gene symbol, e.g., &#x27;CCNB1&#x27;, &#x27;TP53&#x27;, &#x27;EGFR&#x27;, etc. |
| `verbose` | bool | True | 运行提示开关。 |


## 检索 HPA 基因与 Ensembl 标识 {#dnakit.databases.hpa_search_genes_by_query-parameters}

`dnakit.databases.hpa_search_genes_by_query(*, search_query, max_results=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `search_query` | string | 必填 | Gene name, alias, keyword, or cell line name to search for, e.g., &#x27;EGFR&#x27;, &#x27;TP53&#x27;, or &#x27;MCF7&#x27;. |
| `max_results` | integer &#124; null | None | Maximum number of gene matches to return (default: 50). This search does a broad full-text match with no server-side limit, so a short/common query can otherwise return thousands of results; exact gene-symbol matches are ranked first. |
| `verbose` | bool | True | 运行提示开关。 |


## 获取转录因子结合矩阵与元数据 {#dnakit.databases.jaspar_get_matrix-parameters}

`dnakit.databases.jaspar_get_matrix(*, matrix_id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `matrix_id` | string | 必填 | JASPAR matrix ID (e.g., &#x27;MA0002.2&#x27;). You can discover IDs with `JASPAR_get_transcription_factors`. |
| `verbose` | bool | True | 运行提示开关。 |


## 列出 JASPAR 结合矩阵的历史版本 {#dnakit.databases.jaspar_get_matrix_versions-parameters}

`dnakit.databases.jaspar_get_matrix_versions(*, base_id, page=None, page_size=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `base_id` | string | 必填 | JASPAR base matrix ID without version (e.g., &#x27;MA0002&#x27;). You can get `base_id` from `jaspar_get_matrix`. |
| `page` | integer | None | Page number (1-based). {&quot;default&quot;: 1} |
| `page_size` | integer | None | Results per page. {&quot;default&quot;: 20} |
| `verbose` | bool | True | 运行提示开关。 |


## 分页列出 JASPAR 转录因子结合谱 {#dnakit.databases.jaspar_get_transcription_factors-parameters}

`dnakit.databases.jaspar_get_transcription_factors(*, collection=None, page=None, page_size=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `collection` | string | None | JASPAR collection (e.g., CORE). {&quot;default&quot;: &quot;CORE&quot;} |
| `page` | integer | None | Page number (1-based). {&quot;default&quot;: 1} |
| `page_size` | integer | None | Results per page (JASPAR `page_size`). {&quot;default&quot;: 20} |
| `verbose` | bool | True | 运行提示开关。 |


## 按名称、物种或分类检索结合矩阵 {#dnakit.databases.jaspar_search_matrices-parameters}

`dnakit.databases.jaspar_search_matrices(*, search=None, name=None, collection=None, tax_group=None, species=None, page=None, page_size=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `search` | string | None | Free-text search across matrices (e.g., &#x27;CTCF&#x27;). |
| `name` | string | None | Filter by TF name (e.g., &#x27;CTCF&#x27;). |
| `collection` | string | None | Filter by collection (e.g., &#x27;CORE&#x27;). |
| `tax_group` | string | None | Filter by taxonomic group (e.g., &#x27;vertebrates&#x27;). |
| `species` | string | None | Filter by NCBI taxonomy ID (e.g., &#x27;9606&#x27; for human). Sent to JASPAR as its `tax_id` query parameter. |
| `page` | integer | None | Page number (1-based). {&quot;default&quot;: 1} |
| `page_size` | integer | None | Results per page. {&quot;default&quot;: 20} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询基因与细胞类型对应的转录因子结合位点 {#dnakit.databases.re_map_get_transcription_factor_binding-parameters}

`dnakit.databases.re_map_get_transcription_factor_binding(*, gene_name, cell_type=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_name` | string | 必填 | Gene symbol (e.g., BRCA1, TP53, MYC) |
| `cell_type` | string | None | Cell type (e.g., HepG2, K562, MCF-7) {&quot;default&quot;: &quot;HepG2&quot;} |
| `limit` | integer | None | Number of results to return {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 50} |
| `verbose` | bool | True | 运行提示开关。 |


## 按 rsID 查询变异调控注释 {#dnakit.databases.regulome_db_query_variant-parameters}

`dnakit.databases.regulome_db_query_variant(*, rsid=None, variant=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `rsid` | string | None | dbSNP rsID (e.g., &#x27;rs4994&#x27;, &#x27;rs429358&#x27;) |
| `variant` | string | None | Alias for rsid. dbSNP rsID (e.g., &#x27;rs7903146&#x27;, &#x27;rs429358&#x27;). |
| `verbose` | bool | True | 运行提示开关。 |


## 查询基因组区间的候选顺式调控元件 {#dnakit.databases.screen_get_regulatory_elements-parameters}

`dnakit.databases.screen_get_regulatory_elements(*, gene_name=None, query=None, element_type=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_name` | string | None | Gene symbol to search for regulatory elements (e.g., BRCA1, TP53). Alias: &#x27;query&#x27;. |
| `query` | string | None | Gene symbol to search for regulatory elements (alias for gene_name). |
| `element_type` | string | None | Type of regulatory element (promoter, enhancer, insulator) {&quot;default&quot;: &quot;enhancer&quot;} |
| `limit` | integer | None | Number of results to return {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 50} |
| `verbose` | bool | True | 运行提示开关。 |


## 从 UCSC 获取 ENCODE 候选顺式调控元件 {#dnakit.databases.ucsc_get_encode_c_cres-parameters}

`dnakit.databases.ucsc_get_encode_c_cres(*, genome=None, chrom, start, end, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `genome` | string | None | Genome assembly (e.g., &#x27;hg38&#x27;, &#x27;mm10&#x27;). {&quot;default&quot;: &quot;hg38&quot;} |
| `chrom` | string | 必填 | Chromosome name (e.g., &#x27;chr17&#x27;). |
| `start` | integer | 必填 | Start position (0-based). |
| `end` | integer | 必填 | End position (0-based). |
| `verbose` | bool | True | 运行提示开关。 |


## 查询单个直接转录因子结合数据集 {#dnakit.databases.uni_bind_get_dataset-parameters}

`dnakit.databases.uni_bind_get_dataset(*, operation=None, dataset_id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation selector. Must be &#x27;get_dataset&#x27;. {&quot;enum&quot;: &#91;&quot;get_dataset&quot;&#93;, &quot;default&quot;: &quot;get_dataset&quot;} |
| `dataset_id` | string | 必填 | UniBind dataset identifier, e.g. &#x27;EXP030726.neural_stem_cells.SMAD3&#x27;. Obtain from UniBind_search_datasets (the &#x27;dataset_id&#x27; field). |
| `verbose` | bool | True | 运行提示开关。 |


## 列出 UniBind 收录的转录因子 {#dnakit.databases.uni_bind_list_tfs-parameters}

`dnakit.databases.uni_bind_list_tfs(*, operation=None, search=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation selector. Must be &#x27;list_tfs&#x27;. {&quot;enum&quot;: &#91;&quot;list_tfs&quot;&#93;, &quot;default&quot;: &quot;list_tfs&quot;} |
| `search` | string &#124; null | None | Case-insensitive substring filter applied client-side to TF gene symbols, e.g. &#x27;SMAD&#x27;, &#x27;FOX&#x27;, &#x27;GATA&#x27;. Optional; omit to list all TFs. |
| `limit` | integer | None | Maximum number of TF names to return after filtering. {&quot;default&quot;: 200, &quot;minimum&quot;: 1, &quot;maximum&quot;: 1000} |
| `verbose` | bool | True | 运行提示开关。 |


## 按转录因子与细胞条件检索结合数据集 {#dnakit.databases.uni_bind_search_datasets-parameters}

`dnakit.databases.uni_bind_search_datasets(*, operation=None, tf_name=None, species=None, cell_line=None, collection=None, order=None, page=None, page_size=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation selector. Must be &#x27;search_datasets&#x27;. {&quot;enum&quot;: &#91;&quot;search_datasets&quot;&#93;, &quot;default&quot;: &quot;search_datasets&quot;} |
| `tf_name` | string &#124; null | None | Transcription factor gene symbol to filter by, e.g. &#x27;CTCF&#x27;, &#x27;SMAD3&#x27;, &#x27;ESR1&#x27;. Exact match on UniBind&#x27;s TF name. Optional. |
| `species` | string &#124; null | None | Species scientific name to filter by, e.g. &#x27;Homo sapiens&#x27;, &#x27;Mus musculus&#x27;, &#x27;Drosophila melanogaster&#x27;. Optional. |
| `cell_line` | string &#124; null | None | Cell type / tissue / condition name to filter by, e.g. &#x27;neural stem cells&#x27;, &#x27;MCF7&#x27;. Optional. |
| `collection` | string &#124; null | None | UniBind collection: &#x27;Robust&#x27; (high-confidence, metadata-curated) or &#x27;Permissive&#x27; (broader). Optional. {&quot;enum&quot;: &#91;&quot;Robust&quot;, &quot;Permissive&quot;, null&#93;} |
| `order` | string &#124; null | None | Field to order results by, e.g. &#x27;tf_name&#x27; or &#x27;-tf_name&#x27; (prefix &#x27;-&#x27; for descending). Optional. |
| `page` | integer | None | 1-based page number for pagination. {&quot;default&quot;: 1, &quot;minimum&quot;: 1} |
| `page_size` | integer | None | Number of datasets per page (max 1000). {&quot;default&quot;: 25, &quot;minimum&quot;: 1, &quot;maximum&quot;: 1000} |
| `verbose` | bool | True | 运行提示开关。 |
