# Changelog。

## [Unreleased]

## [0.2.0] - 2026-10-13

- 新增约 1,200 个可选依赖封装：`dnakit.xxx()` 内部调用原包，未安装时给出 `dnakit[bio-...]` 安装提示；公开函数由 2,432 个增至 3,648 个，覆盖变异、群体遗传、比对、注释、数据库、表观基因组、微生物组、分子生物学实验、序列、调控和预测等方向。
- 新增约 90 个 `bio-*` 可选依赖组（如 scikit-allel、parasail、pyhmmer、bioframe、hail、scanpy 等）。
- 新增 `parasail_dp_variant()` 与 `parasail_ssw_profile_align()`，按名称部件选择 parasail 的动态规划内核。
- 对已包装函数的别名、再导出和参数分支逐行复核，台账共 15,485 行，覆盖表见 `new_vision/2026-10-10/full_function_coverage.csv`。
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
