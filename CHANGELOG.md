# Changelog。

## [Unreleased]

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
