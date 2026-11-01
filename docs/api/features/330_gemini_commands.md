# GEMINI：查询、导出、统计、遗传模式与基因水平工具

`dnakit.variants` 中的这一组函数在 GEMINI 变异数据库（由 `dnakit.variants.load_gemini_database` 生成的 SQLite 文件）上运行 [GEMINI](https://github.com/arq5x/gemini) 的子命令：合并分块载入、修改样本信息、用 BED/VCF 注释、标记体细胞变异、按表型汇总携带者、表格导出、区域与基因查询、统计（SNP 替换类型计数、位点频率谱、MDS 距离、每样本变异数、汇总）、窗口化 Hardy-Weinberg 或核苷酸多样性、五种遗传模式（常染色体与 X 连锁）、KEGG 通路、LoF 优先级、基因 burden 与 C-alpha 检验、互作伙伴、基因融合、按基因聚合以及性别检查。每个函数用明确的参数列表（不经过 shell）运行一个固定子命令。输出表格的子命令把表格流式写入新文件；修改数据库的子命令默认在副本上工作，除非 `in_place=True`。结果含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`operation` 和 `provenance`；表格还有 `row_count` 和很短的 `rows` 预览。

上游接口的测试版本为 GEMINI 0.30.2（MIT 许可）；DNAKit 把它作为独立程序启动，不复制其代码。用 `conda install -c bioconda gemini` 安装并放到 PATH 上，或传入 `gemini_path`；载入数据库所需的注释文件（约 20 GB，`gemini update --dataonly`）由调用者自行负责，DNAKit 绝不下载。程序缺失时会在创建任何东西之前抛出 `BackendUnavailableError`。**验证状态：** 每个函数的参数映射依据 0.30.2 的命令行定义（`gemini_main.py`）编写，并用替身可执行文件测试；由于这里无法安装 GEMINI 的注释数据，尚未在真实的 GEMINI 数据库上运行，因此本页没有带真实输出的可运行示例。

## 示例数据

本页没有示例文件。

## 1) `dnakit.variants.gemini_amend()`

- **作用：** 用 PED 文件修改 GEMINI 数据库的样本或家系信息（`gemini amend`）。
- **处理方法：** `clear=True` 先把所提供各列的值置为 NULL；在 `output_database` 副本上修改，或 `in_place=True` 直接修改。
- **API：** `gemini_amend(database: 'str', sample_file: 'str | None' = None, *, output_database: 'str | None' = None, in_place: 'bool' = False, clear: 'bool' = False, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 GEMINI 程序（`conda install -c bioconda gemini`）和已有的 GEMINI 数据库；DNAKit 先看 `gemini_path`，再看 PATH。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
没有可运行示例：文档环境中没有安装 GEMINI 及其注释数据（见上面的验证状态）。

## 2) `dnakit.variants.gemini_annotate()`

- **作用：** 用带 tabix 索引的 BED 或 VCF 为变异添加注释列（`gemini annotate`）。
- **处理方法：** `annotation_type` 为 boolean、count 或 extract；extract 可用 `extract_columns` 与 `operations` 合并多个重叠注释。
- **API：** `gemini_annotate(database: 'str', annotation_file: 'str', column_names: 'str', *, output_database: 'str | None' = None, in_place: 'bool' = False, annotation_type: 'str' = 'extract', extract_columns: 'str | None' = None, column_types: 'str | None' = None, operations: 'str | None' = None, region_only: 'bool' = False, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 GEMINI 程序（`conda install -c bioconda gemini`）和已有的 GEMINI 数据库；DNAKit 先看 `gemini_path`，再看 PATH。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
没有可运行示例：文档环境中没有安装 GEMINI 及其注释数据（见上面的验证状态）。

## 3) `dnakit.variants.gemini_burden()`

- **作用：** 按样本（或病例/对照组）统计基因水平的变异 burden（`gemini burden`）。
- **处理方法：** 默认统计编码区的高影响变异；`nonsynonymous=True` 统计所有非同义变异；给出 `cases` 与 `controls` 时增加病例对照 burden 检验。
- **API：** `gemini_burden(database: 'str', output_path: 'str', *, nonsynonymous: 'bool' = False, cases: 'Sequence[str] | None' = None, controls: 'Sequence[str] | None' = None, min_aaf: 'float | None' = None, max_aaf: 'float | None' = None, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 GEMINI 程序（`conda install -c bioconda gemini`）和已有的 GEMINI 数据库；DNAKit 先看 `gemini_path`，再看 PATH。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
没有可运行示例：文档环境中没有安装 GEMINI 及其注释数据（见上面的验证状态）。

## 4) `dnakit.variants.gemini_calpha_test()`

- **作用：** 病例与对照之间的基因水平 C-alpha 关联检验（`gemini burden --calpha`）。
- **处理方法：** `cases` 与 `controls` 必填；`permutations` 为置换次数（0 只给渐近 p 值）。
- **API：** `gemini_calpha_test(database: 'str', output_path: 'str', cases: 'Sequence[str]', controls: 'Sequence[str]', *, permutations: 'int' = 0, nonsynonymous: 'bool' = False, min_aaf: 'float | None' = None, max_aaf: 'float | None' = None, save_tscores: 'bool' = False, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 GEMINI 程序（`conda install -c bioconda gemini`）和已有的 GEMINI 数据库；DNAKit 先看 `gemini_path`，再看 PATH。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
没有可运行示例：文档环境中没有安装 GEMINI 及其注释数据（见上面的验证状态）。

## 5) `dnakit.variants.gemini_carrier_summary_by_phenotype()`

- **作用：** 按样本表型列统计每个变异的携带者与非携带者数（`gemini query --carrier-summary-by-phenotype`）。
- **处理方法：** `query` 选择变异的 SQL；可选 `gt_filter` 与 `region`；结果写入 TSV。
- **API：** `gemini_carrier_summary_by_phenotype(database: 'str', phenotype_column: 'str', output_path: 'str', *, query: 'str' = 'SELECT chrom, start, end, ref, alt, gene FROM variants', gt_filter: 'str | None' = None, region: 'str | None' = None, header: 'bool' = True, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 GEMINI 程序（`conda install -c bioconda gemini`）和已有的 GEMINI 数据库；DNAKit 先看 `gemini_path`，再看 PATH。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
没有可运行示例：文档环境中没有安装 GEMINI 及其注释数据（见上面的验证状态）。

## 6) `dnakit.variants.gemini_dump()`

- **作用：** 把 GEMINI 数据库的表导出为带分隔符的文本（`gemini dump`）。
- **处理方法：** `table` 为 variants、genotypes、samples 或 tfam（TFAM 格式的样本信息）；`separator` 为列分隔符，默认制表符。
- **API：** `gemini_dump(database: 'str', output_path: 'str', *, table: 'str' = 'variants', header: 'bool' = False, separator: 'str' = '\t', annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 GEMINI 程序（`conda install -c bioconda gemini`）和已有的 GEMINI 数据库；DNAKit 先看 `gemini_path`，再看 PATH。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
没有可运行示例：文档环境中没有安装 GEMINI 及其注释数据（见上面的验证状态）。

## 7) `dnakit.variants.gemini_fusions()`

- **作用：** 从结构变异 breakend 检测体细胞基因融合（`gemini fusions`）。
- **处理方法：** 可选 `in_cosmic_census`、`min_quality` 与 `evidence_type`（PE、SR 或 PE,SR）。
- **API：** `gemini_fusions(database: 'str', output_path: 'str', *, in_cosmic_census: 'bool' = False, min_quality: 'float | None' = None, evidence_type: 'str | None' = None, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 GEMINI 程序（`conda install -c bioconda gemini`）和已有的 GEMINI 数据库；DNAKit 先看 `gemini_path`，再看 PATH。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
没有可运行示例：文档环境中没有安装 GEMINI 及其注释数据（见上面的验证状态）。

## 8) `dnakit.variants.gemini_gene_wise()`

- **作用：** 用基因型过滤器按基因聚合变异，报告通过足够多过滤器的基因（`gemini gene_wise`）。
- **处理方法：** `genotype_filters` 为逐变异评估的基因型表达式；`required_genotype_filters` 必须通过；`where` 为选择变异的 SQL。
- **API：** `gemini_gene_wise(database: 'str', output_path: 'str', *, genotype_filters: 'Sequence[str] | None' = None, required_genotype_filters: 'Sequence[str] | None' = None, min_filters: 'int | None' = None, where: 'str | None' = None, filter: 'str | None' = None, columns: 'str | None' = None, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 GEMINI 程序（`conda install -c bioconda gemini`）和已有的 GEMINI 数据库；DNAKit 先看 `gemini_path`，再看 PATH。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
没有可运行示例：文档环境中没有安装 GEMINI 及其注释数据（见上面的验证状态）。

## 9) `dnakit.variants.gemini_inheritance_model()`

- **作用：** 寻找符合孟德尔遗传模式的候选变异（`gemini autosomal_recessive`、`autosomal_dominant`、`x_linked_recessive`、`x_linked_dominant`、`x_linked_de_novo`）。
- **处理方法：** `model` 选择子命令；常染色体模式另有 `lenient` 与 `max_phred_likelihood`，X 连锁模式另有 `x_chromosomes`；其他选项（`min_kindreds`、`families`、`min_gq` 等）共用。
- **API：** `gemini_inheritance_model(database: 'str', output_path: 'str', *, model: 'str' = 'autosomal_recessive', columns: 'str | None' = None, filter: 'str | None' = None, min_kindreds: 'int | None' = None, families: 'str | None' = None, allow_unaffected: 'bool' = False, min_sample_depth: 'int | None' = None, min_gq: 'int | None' = None, lenient: 'bool' = False, max_phred_likelihood: 'int | None' = None, x_chromosomes: 'Sequence[str] | None' = None, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 GEMINI 程序（`conda install -c bioconda gemini`）和已有的 GEMINI 数据库；DNAKit 先看 `gemini_path`，再看 PATH。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
没有可运行示例：文档环境中没有安装 GEMINI 及其注释数据（见上面的验证状态）。

## 10) `dnakit.variants.gemini_interactions()`

- **作用：** 找出样本变异所在基因的互作伙伴（`gemini interactions` 与 `gemini lof_interactions`）。
- **处理方法：** `lof=True` 使用 `lof_interactions`（基因携带 LoF 变异）；`gene` 为 BFS 的根，`radius` 为搜索深度，`edges_file` 替换默认 HPRD 网络。
- **API：** `gemini_interactions(database: 'str', output_path: 'str', *, lof: 'bool' = False, gene: 'str | None' = None, radius: 'int | None' = None, edges_file: 'str | None' = None, variant_mode: 'bool' = False, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 GEMINI 程序（`conda install -c bioconda gemini`）和已有的 GEMINI 数据库；DNAKit 先看 `gemini_path`，再看 PATH。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
没有可运行示例：文档环境中没有安装 GEMINI 及其注释数据（见上面的验证状态）。

## 11) `dnakit.variants.gemini_lof_sieve()`

- **作用：** 按是真实 LoF 的可能性对功能丧失变异排序（`gemini lof_sieve`）。
- **处理方法：** 报告每个 LoF 变异在转录本中的位置、受影响比例以及是否为规范转录本。
- **API：** `gemini_lof_sieve(database: 'str', output_path: 'str', *, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 GEMINI 程序（`conda install -c bioconda gemini`）和已有的 GEMINI 数据库；DNAKit 先看 `gemini_path`，再看 PATH。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
没有可运行示例：文档环境中没有安装 GEMINI 及其注释数据（见上面的验证状态）。

## 12) `dnakit.variants.gemini_merge_chunks()`

- **作用：** 把分块并行载入产生的部分数据库合并成一个数据库（`gemini merge_chunks`）。
- **处理方法：** 输入分块数据库列表，输出新的合并数据库；`index=True` 在最后一次合并时建立索引。
- **API：** `gemini_merge_chunks(chunk_databases: 'Sequence[str]', output_database: 'str', *, vcf_path: 'str | None' = None, anno_type: 'str | None' = None, index: 'bool' = False, skip_pls: 'bool' = False, tempdir: 'str | None' = None, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 GEMINI 程序（`conda install -c bioconda gemini`）和已有的 GEMINI 数据库；DNAKit 先看 `gemini_path`，再看 PATH。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
没有可运行示例：文档环境中没有安装 GEMINI 及其注释数据（见上面的验证状态）。

## 13) `dnakit.variants.gemini_pathways()`

- **作用：** 把基因和变异映射到 KEGG 通路（`gemini pathways`）。
- **处理方法：** `ensembl_version` 为 Ensembl 基因版本（66–71，默认 68）；`lof=True` 报告携带 LoF 变异的个体/基因/位点的通路。
- **API：** `gemini_pathways(database: 'str', output_path: 'str', *, ensembl_version: 'str' = '68', lof: 'bool' = False, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 GEMINI 程序（`conda install -c bioconda gemini`）和已有的 GEMINI 数据库；DNAKit 先看 `gemini_path`，再看 PATH。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
没有可运行示例：文档环境中没有安装 GEMINI 及其注释数据（见上面的验证状态）。

## 14) `dnakit.variants.gemini_region()`

- **作用：** 列出某个基因（按基因名）或基因组区域内的变异（`gemini region --gene` / `--reg`）。
- **处理方法：** `gene` 与 `region` 给且只给一个；可选 `columns`、`filter`、`show_samples` 与 `output_format`。
- **API：** `gemini_region(database: 'str', output_path: 'str', *, gene: 'str | None' = None, region: 'str | None' = None, columns: 'str | None' = None, filter: 'str | None' = None, header: 'bool' = False, show_samples: 'bool' = False, output_format: 'str | None' = None, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 GEMINI 程序（`conda install -c bioconda gemini`）和已有的 GEMINI 数据库；DNAKit 先看 `gemini_path`，再看 PATH。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
没有可运行示例：文档环境中没有安装 GEMINI 及其注释数据（见上面的验证状态）。

## 15) `dnakit.variants.gemini_set_somatic()`

- **作用：** 通过比较肿瘤与正常样本标记体细胞变异 `is_somatic`（`gemini set_somatic`）。
- **处理方法：** 各阈值与 GEMINI 同名选项一致（默认均为 0）；`dry_run=True` 只报告将要标记的内容。
- **API：** `gemini_set_somatic(database: 'str', *, output_database: 'str | None' = None, in_place: 'bool' = False, min_depth: 'float | None' = None, min_quality: 'float | None' = None, min_somatic_score: 'float | None' = None, max_normal_alt_freq: 'float | None' = None, max_normal_alt_count: 'int | None' = None, min_normal_depth: 'int | None' = None, min_tumor_alt_freq: 'float | None' = None, min_tumor_alt_count: 'int | None' = None, min_tumor_depth: 'int | None' = None, chrom: 'str | None' = None, dry_run: 'bool' = False, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 GEMINI 程序（`conda install -c bioconda gemini`）和已有的 GEMINI 数据库；DNAKit 先看 `gemini_path`，再看 PATH。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
没有可运行示例：文档环境中没有安装 GEMINI 及其注释数据（见上面的验证状态）。

## 16) `dnakit.variants.gemini_sex_check()`

- **作用：** 把记录的样本性别与 X 染色体基因型推断的性别对照（`gemini qc --mode sex`）。
- **处理方法：** `chromosome` 为数据库中 X 染色体的名称（默认 `chrX`）。
- **API：** `gemini_sex_check(database: 'str', output_path: 'str', *, chromosome: 'str' = 'chrX', annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 GEMINI 程序（`conda install -c bioconda gemini`）和已有的 GEMINI 数据库；DNAKit 先看 `gemini_path`，再看 PATH。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
没有可运行示例：文档环境中没有安装 GEMINI 及其注释数据（见上面的验证状态）。

## 17) `dnakit.variants.gemini_stats()`

- **作用：** GEMINI 数据库的基本统计（`gemini stats`）。
- **处理方法：** `statistic` 为 snp_counts、sfs、mds、vars_by_sample、gts_by_sample、tstv、tstv_coding、tstv_noncoding 或 summarize（需要 SQL `query`）。
- **API：** `gemini_stats(database: 'str', output_path: 'str', *, statistic: 'str' = 'snp_counts', query: 'str | None' = None, gt_filter: 'str | None' = None, annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 GEMINI 程序（`conda install -c bioconda gemini`）和已有的 GEMINI 数据库；DNAKit 先看 `gemini_path`，再看 PATH。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
没有可运行示例：文档环境中没有安装 GEMINI 及其注释数据（见上面的验证状态）。

## 18) `dnakit.variants.gemini_windower()`

- **作用：** 沿基因组的窗口化群体遗传统计：Hardy-Weinberg 或核苷酸多样性（`gemini windower`）。
- **处理方法：** `analysis` 为 hwe 或 nucl_div；`window_size`、`step_size` 与 `operation`（mean、median、min、max、collapse）。
- **API：** `gemini_windower(database: 'str', output_path: 'str', *, analysis: 'str' = 'hwe', window_size: 'int' = 1000000, step_size: 'int' = 0, operation: 'str' = 'mean', annotation_dir: 'str | None' = None, gemini_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 GEMINI 程序（`conda install -c bioconda gemini`）和已有的 GEMINI 数据库；DNAKit 先看 `gemini_path`，再看 PATH。
- **输入：** 参数含义见上面的说明和 API 签名。
- **输出：** 通用结果字典，另有上面列出的键。
没有可运行示例：文档环境中没有安装 GEMINI 及其注释数据（见上面的验证状态）。
