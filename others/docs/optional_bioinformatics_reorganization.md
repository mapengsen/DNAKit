# 可选生物信息功能按任务归位

删除原来的“可选生物信息功能”集中教程和参数页，将 104 个功能完整归入 15 个任务页面：新增 11 页，补充 4 个已有页。中英文同步，依赖说明留在 FAQ。

## 功能去向

| 目录位置 | 任务页 | 迁移条目数 | 处理 |
| --- | --- | ---: | --- |
| 数据准备 → 数据查询 → 遗传变异与人群关联查询 | [变异记录与群体频率查询](../../docs/api/features/22_database_query_variant_frequency.md) | 11 | 新建子页 |
| 数据准备 → 数据查询 → 遗传变异与人群关联查询 | [结构变异查询](../../docs/api/features/22_database_query_structural_variants.md) | 5 | 新建子页 |
| 数据准备 → 数据查询 → 遗传变异与人群关联查询 | [ClinVar 临床意义查询](../../docs/api/features/22_database_query_clinical_variants.md) | 6 | 新建子页 |
| 数据准备 → 数据查询 → 遗传变异与人群关联查询 | [遗传病关联、剂量与可干预性](../../docs/api/features/22_database_query_genetic_evidence.md) | 10 | 新建子页 |
| 数据准备 → 数据查询 → 遗传变异与人群关联查询 | [肿瘤变异与治疗证据查询](../../docs/api/features/22_database_query_cancer_variants.md) | 14 | 新建子页 |
| 数据准备 → 数据查询 → 遗传变异与人群关联查询 | [GWAS 研究与关联查询](../../docs/api/features/22_database_query_gwas.md) | 12 | 新建子页 |
| 数据准备 → 数据查询 → 遗传变异与人群关联查询 | [队列表型关联与精细定位查询](../../docs/api/features/22_database_query_phewas.md) | 5 | 新建子页 |
| 数据准备 → 数据查询 → 遗传变异与人群关联查询 | [人群、样本与连锁查询](../../docs/api/features/22_database_query_population_resources.md) | 4 | 新建子页 |
| 数据分析 → 变异检测、注释与评分 | [变异注释与标识转换](../../docs/api/features/35_variant_annotation.md) | 12 | 新建子页 |
| 数据分析 → 变异检测、注释与评分 | [变异效应评分与实验数据](../../docs/api/features/36_variant_effect_scores.md) | 11 | 新建子页 |
| 数据分析 → 性质预测 | [变异效应与蛋白质稳定性](../../docs/api/features/23_prediction_variant_effects.md) | 3 | 并入已有页 |
| 数据分析 → 系统发育与分子进化 | [系统发育与分子进化](../../docs/api/features/30_phylogeny.md) | 2 | 并入已有页 |
| 数据处理 → 序列内操作 → 序列拼接与组装 | [OPS-006 序列拼接与组装](../../docs/api/features/04_ops_006_sequence_concatenation.md) | 2 | 并入已有页 |
| 数据分析 → 变异检测、注释与评分 | [变异检测与统计](../../docs/api/features/33_variant_calling.md) | 1 | 并入已有页 |
| 数据分析 → 统计、差异与富集分析 | [基因集与通路富集](../../docs/api/features/37_functional_enrichment.md) | 6 | 新建子页 |

“数据查询”下的八个子页共接收 67 个查询条目。变异注释与标识转换 12 项、评分与实验数据 11 项独立成页；剪接影响预测 3 项并入变异预测；dN/dS 和 Golden Gate 各 2 项并入已有分析或组装页；编码变异比例 1 项并入变异统计；6 项富集分析进入新的富集子页。

同名服务的不同 Python 入口均予保留。例如 `dnakit.annotation` 与 `dnakit.databases` 的参数和返回结构不同，不能只按名称去重。`comparative`/`databases` 的 dN/dS，以及 `molbio`/`databases` 的 Golden Gate 入口也均保留。

## 脚本修改

- `tutorial_api_routes.py`：移除含糊的 `variants` 文档分组，将 90 个生成 API 显式分配到 14 个任务 section；其余 394 个 API 的路由不变。
- `build_native_api_docs.py`：更新 STRING 富集功能之间的双语交叉链接，重新生成后仍指向新位置。现有生成进度条保留。

## 验证

- 每种语言 104 个迁移条目的说明、代码、示例输出和函数锚点完整保留，90 个参数段与迁移前一致。
- 4 个已有目标页的原有功能完整保留；全部功能页编号连续。
- 全站中英文各 70 个功能页、825 个编号条目；484 个生成 API 及其 FAQ 参数锚点有效。
- 严格构建通过，245 个 HTML 页面中的 79,283 个站内链接与资源引用无失效。
- 中英文桌面与手机共 16 个场景通过：查询、变异分析、富集和 FAQ 目录顺序、缩进、编号与宽度正常。
- 两个修改过的 Python 模块已完成代码审查，Ruff、格式检查及生成一致性检查通过；`git diff --check` 通过。
- 未执行科学计算或远程数据库查询；本次检查针对文档迁移与呈现。

## 执行命令

重新生成各任务教程、FAQ 参数页与 API 路由清单：

```bash
/opt/conda/envs/admet_ai/bin/python others/docs/build_native_api_docs.py
```

检查生成文档是否与脚本和路由一致：

```bash
/opt/conda/envs/admet_ai/bin/python others/docs/build_native_api_docs.py --check
```

启动本地预览：

```bash
/opt/conda/envs/admet_ai/bin/python -m mkdocs serve --dev-addr 127.0.0.1:8080
```

## 本次修改文件

- [docs/api/features/04_ops_006_sequence_concatenation.en.md](../../docs/api/features/04_ops_006_sequence_concatenation.en.md)
- [docs/api/features/04_ops_006_sequence_concatenation.md](../../docs/api/features/04_ops_006_sequence_concatenation.md)
- [docs/api/features/23_prediction_variant_effects.en.md](../../docs/api/features/23_prediction_variant_effects.en.md)
- [docs/api/features/23_prediction_variant_effects.md](../../docs/api/features/23_prediction_variant_effects.md)
- [docs/api/features/30_phylogeny.en.md](../../docs/api/features/30_phylogeny.en.md)
- [docs/api/features/30_phylogeny.md](../../docs/api/features/30_phylogeny.md)
- [docs/api/features/33_variant_calling.en.md](../../docs/api/features/33_variant_calling.en.md)
- [docs/api/features/33_variant_calling.md](../../docs/api/features/33_variant_calling.md)
- [docs/api/features/function_tree.en.md](../../docs/api/features/function_tree.en.md)
- [docs/api/features/function_tree.md](../../docs/api/features/function_tree.md)
- [docs/api/index.en.md](../../docs/api/index.en.md)
- [docs/api/index.md](../../docs/api/index.md)
- [docs/faq/api/23_prediction_variant_effects.en.md](../../docs/faq/api/23_prediction_variant_effects.en.md)
- [docs/faq/api/23_prediction_variant_effects.md](../../docs/faq/api/23_prediction_variant_effects.md)
- [docs/faq/api/33_variant_calling.en.md](../../docs/faq/api/33_variant_calling.en.md)
- [docs/faq/api/33_variant_calling.md](../../docs/faq/api/33_variant_calling.md)
- [docs/faq/api/index.en.md](../../docs/faq/api/index.en.md)
- [docs/faq/api/index.md](../../docs/faq/api/index.md)
- [docs/faq/tutorial_basics.en.md](../../docs/faq/tutorial_basics.en.md)
- [docs/faq/tutorial_basics.md](../../docs/faq/tutorial_basics.md)
- [docs/faq/tutorial_integrations.en.md](../../docs/faq/tutorial_integrations.en.md)
- [docs/faq/tutorial_integrations.md](../../docs/faq/tutorial_integrations.md)
- [docs/llm_tools/skills.en.md](../../docs/llm_tools/skills.en.md)
- [docs/llm_tools/skills.md](../../docs/llm_tools/skills.md)
- [mkdocs.yml](../../mkdocs.yml)
- [others/docs/build_native_api_docs.py](build_native_api_docs.py)
- [others/docs/navigation_layout_changes.md](navigation_layout_changes.md)
- [others/docs/tutorial_api_routes.py](tutorial_api_routes.py)
- [others/docs/tutorial_layout_changes.md](tutorial_layout_changes.md)
- [registry/native-api-documentation-validation.json](../../registry/native-api-documentation-validation.json)
- [registry/native-api-documentation.json](../../registry/native-api-documentation.json)

## 本次新建文件

- [docs/api/features/22_database_query_cancer_variants.en.md](../../docs/api/features/22_database_query_cancer_variants.en.md)
- [docs/api/features/22_database_query_cancer_variants.md](../../docs/api/features/22_database_query_cancer_variants.md)
- [docs/api/features/22_database_query_clinical_variants.en.md](../../docs/api/features/22_database_query_clinical_variants.en.md)
- [docs/api/features/22_database_query_clinical_variants.md](../../docs/api/features/22_database_query_clinical_variants.md)
- [docs/api/features/22_database_query_genetic_evidence.en.md](../../docs/api/features/22_database_query_genetic_evidence.en.md)
- [docs/api/features/22_database_query_genetic_evidence.md](../../docs/api/features/22_database_query_genetic_evidence.md)
- [docs/api/features/22_database_query_gwas.en.md](../../docs/api/features/22_database_query_gwas.en.md)
- [docs/api/features/22_database_query_gwas.md](../../docs/api/features/22_database_query_gwas.md)
- [docs/api/features/22_database_query_phewas.en.md](../../docs/api/features/22_database_query_phewas.en.md)
- [docs/api/features/22_database_query_phewas.md](../../docs/api/features/22_database_query_phewas.md)
- [docs/api/features/22_database_query_population_resources.en.md](../../docs/api/features/22_database_query_population_resources.en.md)
- [docs/api/features/22_database_query_population_resources.md](../../docs/api/features/22_database_query_population_resources.md)
- [docs/api/features/22_database_query_structural_variants.en.md](../../docs/api/features/22_database_query_structural_variants.en.md)
- [docs/api/features/22_database_query_structural_variants.md](../../docs/api/features/22_database_query_structural_variants.md)
- [docs/api/features/22_database_query_variant_frequency.en.md](../../docs/api/features/22_database_query_variant_frequency.en.md)
- [docs/api/features/22_database_query_variant_frequency.md](../../docs/api/features/22_database_query_variant_frequency.md)
- [docs/api/features/35_variant_annotation.en.md](../../docs/api/features/35_variant_annotation.en.md)
- [docs/api/features/35_variant_annotation.md](../../docs/api/features/35_variant_annotation.md)
- [docs/api/features/36_variant_effect_scores.en.md](../../docs/api/features/36_variant_effect_scores.en.md)
- [docs/api/features/36_variant_effect_scores.md](../../docs/api/features/36_variant_effect_scores.md)
- [docs/api/features/37_functional_enrichment.en.md](../../docs/api/features/37_functional_enrichment.en.md)
- [docs/api/features/37_functional_enrichment.md](../../docs/api/features/37_functional_enrichment.md)
- [docs/faq/api/22_database_query_cancer_variants.en.md](../../docs/faq/api/22_database_query_cancer_variants.en.md)
- [docs/faq/api/22_database_query_cancer_variants.md](../../docs/faq/api/22_database_query_cancer_variants.md)
- [docs/faq/api/22_database_query_clinical_variants.en.md](../../docs/faq/api/22_database_query_clinical_variants.en.md)
- [docs/faq/api/22_database_query_clinical_variants.md](../../docs/faq/api/22_database_query_clinical_variants.md)
- [docs/faq/api/22_database_query_genetic_evidence.en.md](../../docs/faq/api/22_database_query_genetic_evidence.en.md)
- [docs/faq/api/22_database_query_genetic_evidence.md](../../docs/faq/api/22_database_query_genetic_evidence.md)
- [docs/faq/api/22_database_query_gwas.en.md](../../docs/faq/api/22_database_query_gwas.en.md)
- [docs/faq/api/22_database_query_gwas.md](../../docs/faq/api/22_database_query_gwas.md)
- [docs/faq/api/22_database_query_phewas.en.md](../../docs/faq/api/22_database_query_phewas.en.md)
- [docs/faq/api/22_database_query_phewas.md](../../docs/faq/api/22_database_query_phewas.md)
- [docs/faq/api/22_database_query_population_resources.en.md](../../docs/faq/api/22_database_query_population_resources.en.md)
- [docs/faq/api/22_database_query_population_resources.md](../../docs/faq/api/22_database_query_population_resources.md)
- [docs/faq/api/22_database_query_structural_variants.en.md](../../docs/faq/api/22_database_query_structural_variants.en.md)
- [docs/faq/api/22_database_query_structural_variants.md](../../docs/faq/api/22_database_query_structural_variants.md)
- [docs/faq/api/22_database_query_variant_frequency.en.md](../../docs/faq/api/22_database_query_variant_frequency.en.md)
- [docs/faq/api/22_database_query_variant_frequency.md](../../docs/faq/api/22_database_query_variant_frequency.md)
- [docs/faq/api/30_phylogeny.en.md](../../docs/faq/api/30_phylogeny.en.md)
- [docs/faq/api/30_phylogeny.md](../../docs/faq/api/30_phylogeny.md)
- [docs/faq/api/35_variant_annotation.en.md](../../docs/faq/api/35_variant_annotation.en.md)
- [docs/faq/api/35_variant_annotation.md](../../docs/faq/api/35_variant_annotation.md)
- [docs/faq/api/36_variant_effect_scores.en.md](../../docs/faq/api/36_variant_effect_scores.en.md)
- [docs/faq/api/36_variant_effect_scores.md](../../docs/faq/api/36_variant_effect_scores.md)
- [docs/faq/api/37_functional_enrichment.en.md](../../docs/faq/api/37_functional_enrichment.en.md)
- [docs/faq/api/37_functional_enrichment.md](../../docs/faq/api/37_functional_enrichment.md)
- [others/docs/optional_bioinformatics_reorganization.md](optional_bioinformatics_reorganization.md)

## 本次删除文件

- `docs/api/features/24_optional_bioinformatics.en.md`
- `docs/api/features/24_optional_bioinformatics.md`
- `docs/faq/api/24_optional_bioinformatics.en.md`
- `docs/faq/api/24_optional_bioinformatics.md`

