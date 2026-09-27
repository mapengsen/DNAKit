# 生物信息学教程按任务归类

原综合页的 60 组用法迁入各功能教程；读写与绘图、双链切割与连接等组合条目按用途拆开，形成中英文各 64 个用法条目。原页列出的 94 个 API、类和方法名称全部保留。VCF 的 3 个生成接口与遗传模拟接口分别迁到变异处理和遗传模拟页面。

本次只调整文档组织及生成路由，不修改公开 Python API 的导入路径或行为。依赖、公式、参数约定继续放在 FAQ。

## 功能归位

| 目标页面（中英文同步） | 手写用法条目数 |
| --- | ---: |
| [文件读写](../../docs/api/features/02_io_data.md) | 6 |
| [OPS-001 序列补全](../../docs/api/features/04_ops_001_sequence_direction.md) | 1 |
| [OPS-003 子序列提取](../../docs/api/features/04_ops_003_subsequence_extraction.md) | 2 |
| [OPS-004 序列编辑](../../docs/api/features/04_ops_004_sequence_editing.md) | 1 |
| [OPS-006 序列拼接与组装](../../docs/api/features/04_ops_006_sequence_concatenation.md) | 4 |
| [序列功能搜索](../../docs/api/features/06_patterns.md) | 1 |
| [序列表征与特征处理](../../docs/api/features/08_fingerprints.md) | 2 |
| [序列距离与比对](../../docs/api/features/09_alignment.md) | 3 |
| [通用搜索](../../docs/api/features/09_search.md) | 2 |
| [聚类](../../docs/api/features/10_clustering.md) | 1 |
| [常用评价指标](../../docs/api/features/12_evaluation.md) | 10 |
| [可视化](../../docs/api/features/13_visualization.md) | 4 |
| [工程化与扩展](../../docs/api/features/14_engineering.md) | 1 |
| [测序数据整理](../../docs/api/features/28_read_processing.md) | 4 |
| [序列优化](../../docs/api/features/29_sequence_optimization.md) | 4 |
| [系统发育分析](../../docs/api/features/30_phylogeny.md) | 2 |
| [统计与差异分析](../../docs/api/features/31_statistical_analysis.md) | 10 |
| [测序数据分析](../../docs/api/features/32_sequencing_analysis.md) | 5 |
| [变异检测与 VCF 处理](../../docs/api/features/33_variant_calling.md) | 1 |

## 验证

- 文档生成一致性检查通过：每种语言均为 484 个唯一接口，迁移的 4 个生成接口具有正确教程与参数链接。
- 原有 57 组未拆分用法的代码与文本输出完整保留，拆分后的 7 组用法均有独立示例。
- 中英文各 64 个手写条目的六项结构、锚点、Python 语法及 94 个原 API 名称覆盖检查通过。
- 实际执行 21 个自包含的读写、序列操作、组装、优化、工作流和绘图示例，输出与文档一致。依赖真实文件或外部程序的示例未执行。
- 严格双语构建通过；检查 185 个 HTML 页面和 61731 个内部链接/图片路径，无失效引用。
- 生成路由代码经独立审查，Ruff 检查通过。

## 本次修改或新建的全部文件

- [README.md](../../README.md)
- [README_CN.md](../../README_CN.md)
- [docs/api/features/02_io_data.en.md](../../docs/api/features/02_io_data.en.md)
- [docs/api/features/02_io_data.md](../../docs/api/features/02_io_data.md)
- [docs/api/features/04_ops_001_sequence_direction.en.md](../../docs/api/features/04_ops_001_sequence_direction.en.md)
- [docs/api/features/04_ops_001_sequence_direction.md](../../docs/api/features/04_ops_001_sequence_direction.md)
- [docs/api/features/04_ops_003_subsequence_extraction.en.md](../../docs/api/features/04_ops_003_subsequence_extraction.en.md)
- [docs/api/features/04_ops_003_subsequence_extraction.md](../../docs/api/features/04_ops_003_subsequence_extraction.md)
- [docs/api/features/04_ops_004_sequence_editing.en.md](../../docs/api/features/04_ops_004_sequence_editing.en.md)
- [docs/api/features/04_ops_004_sequence_editing.md](../../docs/api/features/04_ops_004_sequence_editing.md)
- [docs/api/features/04_ops_006_sequence_concatenation.en.md](../../docs/api/features/04_ops_006_sequence_concatenation.en.md)
- [docs/api/features/04_ops_006_sequence_concatenation.md](../../docs/api/features/04_ops_006_sequence_concatenation.md)
- [docs/api/features/06_patterns.en.md](../../docs/api/features/06_patterns.en.md)
- [docs/api/features/06_patterns.md](../../docs/api/features/06_patterns.md)
- [docs/api/features/08_fingerprints.en.md](../../docs/api/features/08_fingerprints.en.md)
- [docs/api/features/08_fingerprints.md](../../docs/api/features/08_fingerprints.md)
- [docs/api/features/09_alignment.en.md](../../docs/api/features/09_alignment.en.md)
- [docs/api/features/09_alignment.md](../../docs/api/features/09_alignment.md)
- [docs/api/features/09_search.en.md](../../docs/api/features/09_search.en.md)
- [docs/api/features/09_search.md](../../docs/api/features/09_search.md)
- [docs/api/features/10_clustering.en.md](../../docs/api/features/10_clustering.en.md)
- [docs/api/features/10_clustering.md](../../docs/api/features/10_clustering.md)
- [docs/api/features/12_evaluation.en.md](../../docs/api/features/12_evaluation.en.md)
- [docs/api/features/12_evaluation.md](../../docs/api/features/12_evaluation.md)
- [docs/api/features/13_visualization.en.md](../../docs/api/features/13_visualization.en.md)
- [docs/api/features/13_visualization.md](../../docs/api/features/13_visualization.md)
- [docs/api/features/14_engineering.en.md](../../docs/api/features/14_engineering.en.md)
- [docs/api/features/14_engineering.md](../../docs/api/features/14_engineering.md)
- [docs/api/features/28_read_processing.en.md](../../docs/api/features/28_read_processing.en.md)
- [docs/api/features/28_read_processing.md](../../docs/api/features/28_read_processing.md)
- [docs/api/features/29_sequence_optimization.en.md](../../docs/api/features/29_sequence_optimization.en.md)
- [docs/api/features/29_sequence_optimization.md](../../docs/api/features/29_sequence_optimization.md)
- [docs/api/features/30_phylogeny.en.md](../../docs/api/features/30_phylogeny.en.md)
- [docs/api/features/30_phylogeny.md](../../docs/api/features/30_phylogeny.md)
- [docs/api/features/31_statistical_analysis.en.md](../../docs/api/features/31_statistical_analysis.en.md)
- [docs/api/features/31_statistical_analysis.md](../../docs/api/features/31_statistical_analysis.md)
- [docs/api/features/32_sequencing_analysis.en.md](../../docs/api/features/32_sequencing_analysis.en.md)
- [docs/api/features/32_sequencing_analysis.md](../../docs/api/features/32_sequencing_analysis.md)
- [docs/api/features/33_variant_calling.en.md](../../docs/api/features/33_variant_calling.en.md)
- [docs/api/features/33_variant_calling.md](../../docs/api/features/33_variant_calling.md)
- [docs/api/features/34_inheritance.en.md](../../docs/api/features/34_inheritance.en.md)
- [docs/api/features/34_inheritance.md](../../docs/api/features/34_inheritance.md)
- [docs/api/features/function_tree.en.md](../../docs/api/features/function_tree.en.md)
- [docs/api/features/function_tree.md](../../docs/api/features/function_tree.md)
- [docs/api/index.en.md](../../docs/api/index.en.md)
- [docs/api/index.md](../../docs/api/index.md)
- [docs/changelog.en.md](../../docs/changelog.en.md)
- [docs/changelog.md](../../docs/changelog.md)
- [docs/faq/api/04_ops_006_sequence_concatenation.en.md](../../docs/faq/api/04_ops_006_sequence_concatenation.en.md)
- [docs/faq/api/04_ops_006_sequence_concatenation.md](../../docs/faq/api/04_ops_006_sequence_concatenation.md)
- [docs/faq/api/33_variant_calling.en.md](../../docs/faq/api/33_variant_calling.en.md)
- [docs/faq/api/33_variant_calling.md](../../docs/faq/api/33_variant_calling.md)
- [docs/faq/api/34_inheritance.en.md](../../docs/faq/api/34_inheritance.en.md)
- [docs/faq/api/34_inheritance.md](../../docs/faq/api/34_inheritance.md)
- [docs/faq/api/index.en.md](../../docs/faq/api/index.en.md)
- [docs/faq/api/index.md](../../docs/faq/api/index.md)
- [docs/faq/bioinformatics.en.md](../../docs/faq/bioinformatics.en.md)
- [docs/faq/bioinformatics.md](../../docs/faq/bioinformatics.md)
- [docs/faq/tutorial_overview.en.md](../../docs/faq/tutorial_overview.en.md)
- [docs/faq/tutorial_overview.md](../../docs/faq/tutorial_overview.md)
- [docs/installation.en.md](../../docs/installation.en.md)
- [docs/installation.md](../../docs/installation.md)
- [mkdocs.yml](../../mkdocs.yml)
- [others/docs/bioinformatics_reorganization.md](bioinformatics_reorganization.md)
- [others/docs/tutorial_api_routes.py](tutorial_api_routes.py)
- [others/docs/tutorial_layout_changes.md](tutorial_layout_changes.md)
- [registry/native-api-documentation-validation.json](../../registry/native-api-documentation-validation.json)
- [registry/native-api-documentation.json](../../registry/native-api-documentation.json)

## 删除的页面

- `docs/api/features/27_bioinformatics.md`
- `docs/api/features/27_bioinformatics.en.md`
- `docs/faq/api/27_bioinformatics.md`
- `docs/faq/api/27_bioinformatics.en.md`

## 执行命令

按当前分类重新生成接口教程、FAQ 参数和索引：

```bash
/opt/conda/envs/admet_ai/bin/python others/docs/build_native_api_docs.py
```

检查生成文件是否与路由一致：

```bash
/opt/conda/envs/admet_ai/bin/python others/docs/build_native_api_docs.py --check
```

在本地预览中英文文档：

```bash
/opt/conda/envs/admet_ai/bin/python -m mkdocs serve --dev-addr 127.0.0.1:8080
```
