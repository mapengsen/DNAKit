# 性质预测文档拆分

将原有 70 个功能条目拆为 9 个任务子目录，删除原集中教程与参数页。中英文同步，每页独立从 `1)` 编号；导航按层级交替使用 `.` 与 `)`。没有新增中转总览页。

## 子目录

| 子目录 | 功能条目数 | 中文教程 | English |
| --- | ---: | --- | --- |
| 1. 自定义模型训练与预测 | 2 | [教程](../../docs/api/features/23_prediction_custom_training.md) | [Tutorial](../../docs/api/features/23_prediction_custom_training.en.md) |
| 2. 表达、调控与接触图预测 | 11 | [教程](../../docs/api/features/23_prediction_regulatory_tracks.md) | [Tutorial](../../docs/api/features/23_prediction_regulatory_tracks.en.md) |
| 3. 剪接与基因组区域预测 | 9 | [教程](../../docs/api/features/23_prediction_splicing_regions.md) | [Tutorial](../../docs/api/features/23_prediction_splicing_regions.en.md) |
| 4. 表观遗传与开放染色质判别 | 11 | [教程](../../docs/api/features/23_prediction_chromatin_classification.md) | [Tutorial](../../docs/api/features/23_prediction_chromatin_classification.en.md) |
| 5. 启动子与增强子判别 | 11 | [教程](../../docs/api/features/23_prediction_promoter_enhancer.md) | [Tutorial](../../docs/api/features/23_prediction_promoter_enhancer.en.md) |
| 6. 变异效应与蛋白质稳定性 | 10 | [教程](../../docs/api/features/23_prediction_variant_effects.md) | [Tutorial](../../docs/api/features/23_prediction_variant_effects.en.md) |
| 7. 物种与分子功能分类 | 6 | [教程](../../docs/api/features/23_prediction_sequence_classification.md) | [Tutorial](../../docs/api/features/23_prediction_sequence_classification.en.md) |
| 8. 分子关系与互作预测 | 4 | [教程](../../docs/api/features/23_prediction_molecular_interactions.md) | [Tutorial](../../docs/api/features/23_prediction_molecular_interactions.en.md) |
| 9. 蛋白质特征与结构预测 | 6 | [教程](../../docs/api/features/23_prediction_protein_features.md) | [Tutorial](../../docs/api/features/23_prediction_protein_features.en.md) |

详细依赖与训练配置仍在 FAQ。原集中参数页拆为调控轨道、变异效应、蛋白质特征三个参数页，保持教程与 FAQ 双向链接。

## 生成脚本变更

`tutorial_api_routes.py` 将原 `models` 文档路由拆为 `model-tracks`、`model-variants`、`model-proteins`，分别管理 1、7、6 个 API。原接口分类、函数实现与调用示例保持不变；重新生成时会写入相应子页面。

## 验证

- 中英文原有 70 个功能条目的说明、代码、输出及功能锚点完整保留，无重复或遗漏；14 个 API 参数段逐段一致。
- 全部 484 个自动生成 API 的路由和 FAQ 参数锚点有效，生成一致性检查通过。
- 全站中英文各 60 个功能页、825 个编号功能条目连续编号。
- 严格构建通过；203 个 HTML 页面、69,025 个站内链接与资源引用无失效。
- 中英文桌面和手机共 4 个浏览器场景通过，九个子目录顺序正确且无水平溢出。
- 路由模块代码审查、Ruff 检查、格式检查与 `git diff --check` 通过。
- 本次为文档拆分，未执行模型训练、预测或远程服务调用。

## 执行命令

重新生成文档及参数页：

```bash
/opt/conda/envs/admet_ai/bin/python others/docs/build_native_api_docs.py
```

检查生成结果是否与路由和模板一致：

```bash
/opt/conda/envs/admet_ai/bin/python others/docs/build_native_api_docs.py --check
```

启动本地预览：

```bash
/opt/conda/envs/admet_ai/bin/python -m mkdocs serve --dev-addr 127.0.0.1:8080
```

## 本次修改文件

- [docs/api/features/function_tree.en.md](../../docs/api/features/function_tree.en.md)
- [docs/api/features/function_tree.md](../../docs/api/features/function_tree.md)
- [docs/api/index.en.md](../../docs/api/index.en.md)
- [docs/api/index.md](../../docs/api/index.md)
- [docs/faq.en.md](../../docs/faq.en.md)
- [docs/faq.md](../../docs/faq.md)
- [docs/faq/api/index.en.md](../../docs/faq/api/index.en.md)
- [docs/faq/api/index.md](../../docs/faq/api/index.md)
- [docs/faq/tutorial_integrations.en.md](../../docs/faq/tutorial_integrations.en.md)
- [docs/faq/tutorial_integrations.md](../../docs/faq/tutorial_integrations.md)
- [docs/installation.en.md](../../docs/installation.en.md)
- [docs/installation.md](../../docs/installation.md)
- [docs/llm_tools/skills.en.md](../../docs/llm_tools/skills.en.md)
- [docs/llm_tools/skills.md](../../docs/llm_tools/skills.md)
- [mkdocs.yml](../../mkdocs.yml)
- [others/docs/navigation_layout_changes.md](navigation_layout_changes.md)
- [others/docs/tutorial_api_routes.py](tutorial_api_routes.py)
- [others/docs/tutorial_layout_changes.md](tutorial_layout_changes.md)
- [registry/native-api-documentation-validation.json](../../registry/native-api-documentation-validation.json)
- [registry/native-api-documentation.json](../../registry/native-api-documentation.json)

## 本次新建文件

- [docs/api/features/23_prediction_chromatin_classification.en.md](../../docs/api/features/23_prediction_chromatin_classification.en.md)
- [docs/api/features/23_prediction_chromatin_classification.md](../../docs/api/features/23_prediction_chromatin_classification.md)
- [docs/api/features/23_prediction_custom_training.en.md](../../docs/api/features/23_prediction_custom_training.en.md)
- [docs/api/features/23_prediction_custom_training.md](../../docs/api/features/23_prediction_custom_training.md)
- [docs/api/features/23_prediction_molecular_interactions.en.md](../../docs/api/features/23_prediction_molecular_interactions.en.md)
- [docs/api/features/23_prediction_molecular_interactions.md](../../docs/api/features/23_prediction_molecular_interactions.md)
- [docs/api/features/23_prediction_promoter_enhancer.en.md](../../docs/api/features/23_prediction_promoter_enhancer.en.md)
- [docs/api/features/23_prediction_promoter_enhancer.md](../../docs/api/features/23_prediction_promoter_enhancer.md)
- [docs/api/features/23_prediction_protein_features.en.md](../../docs/api/features/23_prediction_protein_features.en.md)
- [docs/api/features/23_prediction_protein_features.md](../../docs/api/features/23_prediction_protein_features.md)
- [docs/api/features/23_prediction_regulatory_tracks.en.md](../../docs/api/features/23_prediction_regulatory_tracks.en.md)
- [docs/api/features/23_prediction_regulatory_tracks.md](../../docs/api/features/23_prediction_regulatory_tracks.md)
- [docs/api/features/23_prediction_sequence_classification.en.md](../../docs/api/features/23_prediction_sequence_classification.en.md)
- [docs/api/features/23_prediction_sequence_classification.md](../../docs/api/features/23_prediction_sequence_classification.md)
- [docs/api/features/23_prediction_splicing_regions.en.md](../../docs/api/features/23_prediction_splicing_regions.en.md)
- [docs/api/features/23_prediction_splicing_regions.md](../../docs/api/features/23_prediction_splicing_regions.md)
- [docs/api/features/23_prediction_variant_effects.en.md](../../docs/api/features/23_prediction_variant_effects.en.md)
- [docs/api/features/23_prediction_variant_effects.md](../../docs/api/features/23_prediction_variant_effects.md)
- [docs/faq/api/23_prediction_protein_features.en.md](../../docs/faq/api/23_prediction_protein_features.en.md)
- [docs/faq/api/23_prediction_protein_features.md](../../docs/faq/api/23_prediction_protein_features.md)
- [docs/faq/api/23_prediction_regulatory_tracks.en.md](../../docs/faq/api/23_prediction_regulatory_tracks.en.md)
- [docs/faq/api/23_prediction_regulatory_tracks.md](../../docs/faq/api/23_prediction_regulatory_tracks.md)
- [docs/faq/api/23_prediction_variant_effects.en.md](../../docs/faq/api/23_prediction_variant_effects.en.md)
- [docs/faq/api/23_prediction_variant_effects.md](../../docs/faq/api/23_prediction_variant_effects.md)
- [others/docs/property_prediction_split.md](property_prediction_split.md)

## 本次删除文件

- `docs/api/features/23_deep_learning_property_prediction.en.md`
- `docs/api/features/23_deep_learning_property_prediction.md`
- `docs/faq/api/23_deep_learning_property_prediction.en.md`
- `docs/faq/api/23_deep_learning_property_prediction.md`

