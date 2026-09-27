# 全站导航与功能页整理

中英文导航的每级子目录逐层缩进，同级项目从 `1.` 连续编号；删除“查询总览与示例”页面，将其独有的 12 个查询函数归入相应功能页。功能页去除“接口列表”等分组标题，直接以 `1)`、`2)` 连续展示功能用法。原有函数与分组锚点保留。

生成脚本已同步采用二级功能标题和不可见的分组锚点，重新生成文档时会保留此结构。

## 验证

- 中英文各 90 个导航页面、90 个子项通过编号和路径一致性检查。
- 中英文各 52 个功能页、825 个编号功能条目通过标题层级与连续编号检查。
- 484 个自动生成 API 条目及对应 FAQ 参数锚点完整；54 个生成区块的示例代码和输出保持不变。
- 新整理的 24 个双语查询示例通过语法与函数签名检查；未执行联网查询，页面明确标注示例结果为返回结构。
- MkDocs 中英文严格构建通过；183 个 HTML 页面中的 62,824 个站内链接与资源引用无失效。
- 双语桌面与手机共 12 个浏览器场景通过缩进、编号、折叠和当前页状态检查。
- 生成脚本通过代码审查、Ruff 检查、格式检查及生成一致性检查；`git diff --check` 通过。

## 执行命令

重新生成自动管理的功能文档：

```bash
/opt/conda/envs/admet_ai/bin/python others/docs/build_native_api_docs.py
```

检查已生成文档与生成脚本是否一致：

```bash
/opt/conda/envs/admet_ai/bin/python others/docs/build_native_api_docs.py --check
```

启动本地预览：

```bash
/opt/conda/envs/admet_ai/bin/python -m mkdocs serve --dev-addr 127.0.0.1:8080
```

## 本次修改与新建的完整文件清单

- [docs/api/features/03_validation.en.md](../../docs/api/features/03_validation.en.md)
- [docs/api/features/03_validation.md](../../docs/api/features/03_validation.md)
- [docs/api/features/04_ops_001_sequence_direction.en.md](../../docs/api/features/04_ops_001_sequence_direction.en.md)
- [docs/api/features/04_ops_001_sequence_direction.md](../../docs/api/features/04_ops_001_sequence_direction.md)
- [docs/api/features/04_ops_002_transcription_translation.en.md](../../docs/api/features/04_ops_002_transcription_translation.en.md)
- [docs/api/features/04_ops_002_transcription_translation.md](../../docs/api/features/04_ops_002_transcription_translation.md)
- [docs/api/features/04_ops_004_sequence_editing.en.md](../../docs/api/features/04_ops_004_sequence_editing.en.md)
- [docs/api/features/04_ops_004_sequence_editing.md](../../docs/api/features/04_ops_004_sequence_editing.md)
- [docs/api/features/04_ops_005_mutation_generation.en.md](../../docs/api/features/04_ops_005_mutation_generation.en.md)
- [docs/api/features/04_ops_005_mutation_generation.md](../../docs/api/features/04_ops_005_mutation_generation.md)
- [docs/api/features/04_ops_006_sequence_concatenation.en.md](../../docs/api/features/04_ops_006_sequence_concatenation.en.md)
- [docs/api/features/04_ops_006_sequence_concatenation.md](../../docs/api/features/04_ops_006_sequence_concatenation.md)
- [docs/api/features/04_ops_007_masking.en.md](../../docs/api/features/04_ops_007_masking.en.md)
- [docs/api/features/04_ops_007_masking.md](../../docs/api/features/04_ops_007_masking.md)
- [docs/api/features/04_ops_007_trimming.en.md](../../docs/api/features/04_ops_007_trimming.en.md)
- [docs/api/features/04_ops_007_trimming.md](../../docs/api/features/04_ops_007_trimming.md)
- [docs/api/features/04_ops_008_circular_sequence_operations.en.md](../../docs/api/features/04_ops_008_circular_sequence_operations.en.md)
- [docs/api/features/04_ops_008_circular_sequence_operations.md](../../docs/api/features/04_ops_008_circular_sequence_operations.md)
- [docs/api/features/05_all_descriptors.en.md](../../docs/api/features/05_all_descriptors.en.md)
- [docs/api/features/05_all_descriptors.md](../../docs/api/features/05_all_descriptors.md)
- [docs/api/features/07_physicochemical.en.md](../../docs/api/features/07_physicochemical.en.md)
- [docs/api/features/07_physicochemical.md](../../docs/api/features/07_physicochemical.md)
- [docs/api/features/09_alignment.en.md](../../docs/api/features/09_alignment.en.md)
- [docs/api/features/09_alignment.md](../../docs/api/features/09_alignment.md)
- [docs/api/features/09_search.en.md](../../docs/api/features/09_search.en.md)
- [docs/api/features/09_search.md](../../docs/api/features/09_search.md)
- [docs/api/features/09_similarity_alignment.en.md](../../docs/api/features/09_similarity_alignment.en.md)
- [docs/api/features/09_similarity_alignment.md](../../docs/api/features/09_similarity_alignment.md)
- [docs/api/features/10_datasets.en.md](../../docs/api/features/10_datasets.en.md)
- [docs/api/features/10_datasets.md](../../docs/api/features/10_datasets.md)
- [docs/api/features/14_engineering.en.md](../../docs/api/features/14_engineering.en.md)
- [docs/api/features/14_engineering.md](../../docs/api/features/14_engineering.md)
- [docs/api/features/15_download.en.md](../../docs/api/features/15_download.en.md)
- [docs/api/features/15_download.md](../../docs/api/features/15_download.md)
- [docs/api/features/22_database_query_discovery.en.md](../../docs/api/features/22_database_query_discovery.en.md)
- [docs/api/features/22_database_query_discovery.md](../../docs/api/features/22_database_query_discovery.md)
- [docs/api/features/22_database_query_encode.en.md](../../docs/api/features/22_database_query_encode.en.md)
- [docs/api/features/22_database_query_encode.md](../../docs/api/features/22_database_query_encode.md)
- [docs/api/features/22_database_query_literature.en.md](../../docs/api/features/22_database_query_literature.en.md)
- [docs/api/features/22_database_query_literature.md](../../docs/api/features/22_database_query_literature.md)
- [docs/api/features/22_database_query_organisms.en.md](../../docs/api/features/22_database_query_organisms.en.md)
- [docs/api/features/22_database_query_organisms.md](../../docs/api/features/22_database_query_organisms.md)
- [docs/api/features/22_database_query_pathways.en.md](../../docs/api/features/22_database_query_pathways.en.md)
- [docs/api/features/22_database_query_pathways.md](../../docs/api/features/22_database_query_pathways.md)
- [docs/api/features/22_database_query_phenotypes.en.md](../../docs/api/features/22_database_query_phenotypes.en.md)
- [docs/api/features/22_database_query_phenotypes.md](../../docs/api/features/22_database_query_phenotypes.md)
- [docs/api/features/22_database_query_regulation.en.md](../../docs/api/features/22_database_query_regulation.en.md)
- [docs/api/features/22_database_query_regulation.md](../../docs/api/features/22_database_query_regulation.md)
- [docs/api/features/22_database_query_sequences.en.md](../../docs/api/features/22_database_query_sequences.en.md)
- [docs/api/features/22_database_query_sequences.md](../../docs/api/features/22_database_query_sequences.md)
- [docs/api/features/22_database_query_targets.en.md](../../docs/api/features/22_database_query_targets.en.md)
- [docs/api/features/22_database_query_targets.md](../../docs/api/features/22_database_query_targets.md)
- `docs/api/features/23_deep_learning_property_prediction.en.md`（已删除，按预测任务拆分）
- `docs/api/features/23_deep_learning_property_prediction.md`（已删除，按预测任务拆分）
- `docs/api/features/24_optional_bioinformatics.en.md`（已按任务拆分，原集中页面已删除）
- `docs/api/features/24_optional_bioinformatics.md`（已按任务拆分，原集中页面已删除）
- [docs/api/features/26_gene_annotation.en.md](../../docs/api/features/26_gene_annotation.en.md)
- [docs/api/features/26_gene_annotation.md](../../docs/api/features/26_gene_annotation.md)
- [docs/api/features/33_variant_calling.en.md](../../docs/api/features/33_variant_calling.en.md)
- [docs/api/features/33_variant_calling.md](../../docs/api/features/33_variant_calling.md)
- [docs/api/features/34_inheritance.en.md](../../docs/api/features/34_inheritance.en.md)
- [docs/api/features/34_inheritance.md](../../docs/api/features/34_inheritance.md)
- [docs/api/features/function_tree.en.md](../../docs/api/features/function_tree.en.md)
- [docs/api/features/function_tree.md](../../docs/api/features/function_tree.md)
- [docs/api/index.en.md](../../docs/api/index.en.md)
- [docs/api/index.md](../../docs/api/index.md)
- [docs/faq/tutorial_basics.en.md](../../docs/faq/tutorial_basics.en.md)
- [docs/faq/tutorial_basics.md](../../docs/faq/tutorial_basics.md)
- [docs/faq/tutorial_integrations.en.md](../../docs/faq/tutorial_integrations.en.md)
- [docs/faq/tutorial_integrations.md](../../docs/faq/tutorial_integrations.md)
- [docs/llm_tools/skills.en.md](../../docs/llm_tools/skills.en.md)
- [docs/llm_tools/skills.md](../../docs/llm_tools/skills.md)
- [docs/stylesheets/extra.css](../../docs/stylesheets/extra.css)
- [mkdocs.yml](../../mkdocs.yml)
- [others/docs/build_native_api_docs.py](build_native_api_docs.py)
- [others/docs/navigation_layout_changes.md](navigation_layout_changes.md)（新建）
- [others/docs/tutorial_layout_changes.md](tutorial_layout_changes.md)
- [registry/native-api-documentation-validation.json](../../registry/native-api-documentation-validation.json)

## 本次删除的文件

- `docs/api/features/22_database_query.en.md`
- `docs/api/features/22_database_query.md`
