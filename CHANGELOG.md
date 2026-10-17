# Changelog。

## [Unreleased]

- 将 615 个原先只返回输入摘要、状态却为 `ok` 的候选函数改为真实实现，涉及 `io_formats`、`variants`、`clinical`、`gwas`、`singlecell_spatial`、`transcriptomics`、`ml_models`、`knowledge_graph`、`sequencing_qc`、`assembly`、`microbiome` 等 29 个模块。约 518 个函数按各自文档定义计算结果或校验输入；97 个在不带选项的普通调用下返回 `status: "not_implemented"`，并附 `required_dependency` 和 `message` 说明缺少什么（外部数据库、训练好的模型、比对或变异检测软件、必要选项，或函数名没有确定含义），不再返回伪成功。启发式或近似算法（如基于规则的 CRISPR 打分、k-mer 近似 ANI、KING 亲缘系数）在函数说明里写明定义。
- 重命名三个名字与实际行为不符的函数：`variants.impute_missing_genotypes` 改为 `fill_missing_genotypes_with_mode`（按位点取最常见基因型），`variants.calculate_identity_by_descent` 改为 `calculate_kinship_king`（KING 亲缘系数，不是 IBD 片段），`ml_models.export_uncertainty_model` 改为 `export_conformal_uncertainty_model`（split-conformal 区间）。旧名在 0.2.0 中只是不能工作的占位函数，不保留别名。
- `design_crispr.calculate_primer_tm` 和 `evaluate_primer_pair` 改为真实现：前者长度不超过 13 nt 用 Wallace 规则（`ATGCATGCAT` 为 28 °C），更长的用 `dnakit.thermodynamics` 的最近邻法；后者给出正反向引物各自的 Tm、GC、GC 夹和两引物间的互补性，并按可调阈值给出标志。返回结构改为 `primers` / `forward` / `reverse`。`encode_data_into_dna` 对字节和文件按真实字节编码，不再编码成 JSON 文本。
- 脱敏函数改名并重做：`workflow_security.remove_phi_from_metadata` 改为 `pseudonymize_metadata_fields`，`knowledge_graph.deidentify_genomic_metadata` 改为 `pseudonymize_genomic_metadata_fields`。旧实现只处理最外层的几个键，并用不加盐的 SHA-256 前缀，姓名可被猜出。新实现递归处理嵌套字段，默认整体替换（`redact`），可选 `drop` 或带盐的 HMAC（`hash`，必须提供盐），并扫描文本中的邮箱、电话和 SSN；输入不是映射时报错；返回结果不含明文或明文的哈希，始终警告“只按键名和模式匹配，不保证真正脱敏”。
- `io_formats` 的 16 个 `compress_*` / `decompress_*` 改为对字节无损：文本、二进制、文件路径都能逐字节往返，zarr 目录先打成 tar 再 gzip，解压时拒绝绝对路径、`..` 和链接；gzip 结果固定（同样输入每次相同）；无效或损坏的载荷现在报错，不再返回 `ok`。`validate_map` / `validate_ped` / `validate_paf` 对每一行报告全部问题，不再泄露 Python 原始异常信息；`write_phylip` 按 `datatype` 校验字符；`index_map` 报告重复 ID；`convert_nexus` 转换矩阵并提示丢弃了树。
- 修复若干与说明不符的行为：`epigenome.calculate` 的区间重叠不再把一个区间重复计数；`genome3d.annotate` 默认报告 TAD 边缘域（`include_edge_domains`）；`detect_de_novo_variants` 的 `rejected_by` 取固定顺序的第一个失败检查；`infer_parentage` 新增 `not_assessed`；`reconstruct` 的 `log_likelihood` 取最终频率；`experimental_design.check` / `export` / `validate` 接受 `design()` 的完整返回值；`external_adapters.build_command` 新增 `subcommand`；`extract_barcodes` 的 `n_extracted` 只计提取到的读段并新增 `n_unmatched`；`build_cell_peak_matrix` 的字段 `n_genes` 改为 `n_peaks`；遗传相关系数超出 [-1, 1]、遗传力超出 [0, 1]、缺少人群频率、分不开的亲缘关系都会给出警告；一批函数说明改为与实际行为一致。
- 所有候选函数返回的 `input_summary` 按输入的真实格式统计：FASTQ 的质量字符、Newick 的标签和 `None` 不再被数成碱基。
- 工具目录登记新实现函数的真实副作用：会写文件的 142 个函数（仅在传入输出路径时写入，默认不覆盖）需要文件系统授权，调用外部程序的 9 个和加载或训练模型的 4 个需要相应授权，Agent 须设置 `allow_side_effects` 才能调用。
- 重写 `docs/api/features/candidate_catalog/` 中 29 个模块、约 640 个函数的中英文小节（含 615 个新实现的函数，以及压缩、解压、计数、引物 Tm、元数据脱敏等 24 个原本就可用但被写成“本地适配器”的函数）：给出真实的作用、方法、输入、输出，示例均实际运行并粘贴真实输出，并删去“只报告缺少依赖”的旧说明。其余约 200 个函数（如 `fasta_fastq`、`sequence`、`alignment`）原本就是真实现，小节保持原样。
- 新增图形界面：`pip install "dnakit[gui]"` 后运行 `dnakit gui` 打开桌面窗口（pywebview；无窗口组件时改用浏览器），也可用 `dnakit.gui.build_app(GuiSettings.public_demo())` 部署网页演示。两者共用同一套界面代码，包含单序列工具、数据集工具、全部 Agent 可调用函数的检索与运行，以及 Skills 检索；桌面版不限时，写文件等有副作用的函数经确认后执行，网页版限时且不执行这类函数。Linux 的 Qt 窗口后端为可选依赖组 `gui-qt`。
- `dnakit compare --help` 中的方法名由 `levenshtein` 更正为实际可用的 `edit`。
- MCP 目录新增 17 个第三方服务器（共 49 个），每个条目固定上游提交、记录许可证、环境变量和对外影响；12 个已通过 stdio 实际启动并核对工具名。未收录 10 个候选，原因见文档 110。

## [0.2.0] - 2026-10-13

- 新增约 1,200 个可选依赖封装：`dnakit.xxx()` 内部调用原包，未安装时给出 `dnakit[bio-...]` 安装提示；Python 公开函数由 1,826 个增至 3,042 个，覆盖变异、群体遗传、比对、注释、数据库、表观基因组、微生物组、分子生物学实验、序列、调控和预测等方向。
- 新增约 90 个 `bio-*` 可选依赖组（如 scikit-allel、parasail、pyhmmer、bioframe、hail、scanpy 等）。
- 新增 `parasail_dp_variant()` 与 `parasail_ssw_profile_align()`，按名称部件选择 parasail 的动态规划内核。
- 对已包装函数的别名、再导出和参数分支逐行复核，台账共 15,485 行，覆盖表见 `new_vision/2026-10-10/full_function_coverage.csv`（作者本地的审计文件，不在仓库中）。
- 补充第 48 到 109 页中英文 API 文档，每个函数给出作用、处理方法、API、输入、输出以及真实运行的示例。
- 导入可选依赖时不再改动调用方的日志设置；测试配置会在用例之间还原全局日志状态。

## [0.1.5] - 2026-09-27

- 修复已安装包中的工具质量检查，并确保 Git 检出后保留 Skills 资源的可执行权限。
- 新增 Agent 检查点、恢复、持久知识，以及工具生成、修复、优化和贡献流程。
- 扩展自定义训练后端与调度、配对序列分块、表格数据集与采样、评价绘图和参考一致性检查。
- 新增可通过 pip 安装的基因组学后端，更新中英文 API 教程。
- 新增统一 DNA／基因组学 Skills 目录、搜索与安装入口，并保留来源、许可证和资源校验信息。
- 新增数据库、ENCODE 和工作流接口，补充 484 个公开函数的中英文文档及源码驱动的覆盖校验。
- 扩展统一文件输入输出、序列处理并行执行、进度显示、模型下载与 Agent 工具调用。
- 更新中英文 README、文档导航、MCP 安装指南与可视化示例。

## [0.1.4] - 2026-09-09

- 新增基于参考基因组和 minimap2 的基因注释，支持参考资源准备及 Python/CLI 调用。
- 新增通用同步/异步 LLM 调用、提供方配置和有界文件输入。
- 新增传统机器学习性质建模与分类/回归标签分布可视化。
- 新增 `dnakit.bioinformatics` 生物信息分析 API、可选依赖及对应验证。
- 新增自定义 DNA 性质训练，支持基础模型线性探测、全量微调与文件结构分析。

- 将测试、示例、Notebook、MkDocs override 和本地验证资源统一迁移至 `others/`，并同步更新打包、CI、文档和工具配置。
- Diversity 和 Novelty 新增论文定义的原始 Levenshtein 距离计算方式，保留原有相似度方法，并提供资源上限、进度条、测试及中英文引用说明。
- 新增 NT Revised 18 项与 Genomic Benchmarks 9 项完整微调 checkpoint 分类预测，支持独立的统一下载说明页、本地可配置路径、安全加载、Agent 调用及逐功能文档。

## [0.1.3] - 2026-09-01

- 新增 AlphaGenome、Enformer、Evo 2、GENERator、LucaOneTasks 和 SegmentNT 的预训练性质预测接口，并对模型下载、远程代码和输出大小设置显式边界。
- 新增 Ensembl VEP、ClinVar、dbSNP、gnomAD、dN/dS 和 Golden Gate 等按功能命名的可选生物信息 API。
- 新增 `dnakit[agent]` 可选依赖、`dnakit-mcp` 服务和紧凑工具目录，使支持 MCP 的 Agent 可以搜索、查看并有界调用 DNAKit 公开功能。
- 新增中英文微信社区页面和文档导航入口。

## [0.1.2] - 2026-08-30

- 新增 Bioconda、GNU Guix 与 Galaxy Tool Shed 的发布适配文件。
- 可视化统一导出接口新增 JPG，支持通过 `image_type` 选择 PNG、SVG 或 JPG，无扩展名时默认导出 PNG；`SequencePlotConfig` 新增 `column_spacing`、`line_spacing` 和 `symbol_map`，可控制文字列间距、序列行间距及 DNA/IUPAC 显示字符；所有绘图改为无装饰性图内标题的正方形画布。
- GitHub 仓库首页和 PyPI 项目说明默认改为英文，并保留简体中文切换入口。
- 文档站提供完整中英文页面、导航和独立搜索结果，站内跳转保持当前语言。
- 精简不再对外展示的旧规划、演示和交付报告文档。

## [0.1.1] - 2026-08-29

- 发布首个不含开发版本后缀的 DNAKit 版本。
- 精简仓库根目录文件，将第三方声明合并到 `DISCLAIMER.md`，将引用信息移至 `README.md`。
- 更新 README、安装文档、打包配置和 GitHub 发布工作流。

## [0.1.0.dev0] - 2026-08-28

初始版本
