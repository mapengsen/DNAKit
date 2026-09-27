# 迁移功能标题编号修正

中英文各 64 个迁移功能标题补齐 `1)`、`2)` 等编号；每页按正文顺序续排。原有子层级编号（如 `1.1`、`2.1`）、功能锚点、示例代码和结果保留。子序列提取和序列编辑页原先没有标题的开头功能补为 `1)`。

生成脚本改为在生成块实际出现的位置续号，支持生成块前后混排手写功能；忽略代码围栏及小数形式子编号，并检查中英文生成编号一致。保留原有进度显示。

## 编号范围（中英文一致）

| 页面 | 本次迁移功能编号 |
| --- | --- |
| [文件读写](../../docs/api/features/02_io_data.md) | 5)–10) |
| [OPS-001 序列补全](../../docs/api/features/04_ops_001_sequence_direction.md) | 5) |
| [OPS-003 子序列提取](../../docs/api/features/04_ops_003_subsequence_extraction.md) | 2)–3) |
| [OPS-004 序列编辑](../../docs/api/features/04_ops_004_sequence_editing.md) | 3) |
| [OPS-006 序列拼接与组装](../../docs/api/features/04_ops_006_sequence_concatenation.md) | 6)–9) |
| [序列功能搜索](../../docs/api/features/06_patterns.md) | 11) |
| [序列表征与特征处理](../../docs/api/features/08_fingerprints.md) | 6)–7) |
| [序列距离与比对](../../docs/api/features/09_alignment.md) | 3)–5) |
| [通用搜索](../../docs/api/features/09_search.md) | 7)–8) |
| [聚类](../../docs/api/features/10_clustering.md) | 7) |
| [常用评价指标](../../docs/api/features/12_evaluation.md) | 10)–19) |
| [可视化](../../docs/api/features/13_visualization.md) | 10)–13) |
| [工程化与扩展](../../docs/api/features/14_engineering.md) | 16) |
| [测序数据整理](../../docs/api/features/28_read_processing.md) | 1)–4) |
| [序列优化](../../docs/api/features/29_sequence_optimization.md) | 1)–4) |
| [系统发育分析](../../docs/api/features/30_phylogeny.md) | 1)–2) |
| [统计与差异分析](../../docs/api/features/31_statistical_analysis.md) | 1)–10) |
| [测序数据分析](../../docs/api/features/32_sequencing_analysis.md) | 1)–5) |
| [变异检测与 VCF 处理](../../docs/api/features/33_variant_calling.md) | 1) |

## 验证

- 38 个双语页面的主编号均从 1 连续排列，无缺号、重号；128 个迁移标题及目录编号已在构建页面中核验。
- 484 个生成接口一致性检查、Ruff 和严格双语构建通过；功能锚点与示例未改变。
- 全站 185 个 HTML 页面、61351 个内部链接检查通过。
- 生成脚本经过独立代码审查；混排、新路由、双语差异、重复块与代码围栏等离线场景已验证。

## 本次修改或新建的全部文件

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
- [others/docs/build_native_api_docs.py](build_native_api_docs.py)
- [others/docs/tutorial_numbering_changes.md](tutorial_numbering_changes.md)

## 执行命令

检查生成文档与脚本的编号是否一致：

```bash
/opt/conda/envs/admet_ai/bin/python others/docs/build_native_api_docs.py --check
```

按当前规则重新生成文档：

```bash
/opt/conda/envs/admet_ai/bin/python others/docs/build_native_api_docs.py
```
