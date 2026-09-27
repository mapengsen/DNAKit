# 文档功能编号清理

本次统一移除中英文页面标题、导航、目录、功能树、FAQ 表格及正文中的内部功能编号，保留原有序号和功能名称。二级及以下标题保留原锚点，以兼容现有链接。

Python 更新：`src/dnakit/evaluation/frechet.py` 的函数文档字符串改用功能名称，确保自动生成的 API 参考页同步清理。本次 Python 修改仅涉及说明文字。

验证范围：文档生成一致性、MkDocs 严格构建、生成页面的可见编号扫描、章节锚点与代码示例对比，以及独立代码审查。

## 修改文件

共 75 个文件，包含本清单。

- [docs/api/features/01_core_objects.en.md](../../docs/api/features/01_core_objects.en.md)
- [docs/api/features/01_core_objects.md](../../docs/api/features/01_core_objects.md)
- [docs/api/features/02_io_data.en.md](../../docs/api/features/02_io_data.en.md)
- [docs/api/features/02_io_data.md](../../docs/api/features/02_io_data.md)
- [docs/api/features/03_standardization.en.md](../../docs/api/features/03_standardization.en.md)
- [docs/api/features/03_standardization.md](../../docs/api/features/03_standardization.md)
- [docs/api/features/04_ops_001_sequence_direction.en.md](../../docs/api/features/04_ops_001_sequence_direction.en.md)
- [docs/api/features/04_ops_001_sequence_direction.md](../../docs/api/features/04_ops_001_sequence_direction.md)
- [docs/api/features/04_ops_002_transcription_translation.en.md](../../docs/api/features/04_ops_002_transcription_translation.en.md)
- [docs/api/features/04_ops_002_transcription_translation.md](../../docs/api/features/04_ops_002_transcription_translation.md)
- [docs/api/features/04_ops_003_subsequence_extraction.en.md](../../docs/api/features/04_ops_003_subsequence_extraction.en.md)
- [docs/api/features/04_ops_003_subsequence_extraction.md](../../docs/api/features/04_ops_003_subsequence_extraction.md)
- [docs/api/features/04_ops_004_sequence_editing.en.md](../../docs/api/features/04_ops_004_sequence_editing.en.md)
- [docs/api/features/04_ops_004_sequence_editing.md](../../docs/api/features/04_ops_004_sequence_editing.md)
- [docs/api/features/04_ops_005_mutation_generation.en.md](../../docs/api/features/04_ops_005_mutation_generation.en.md)
- [docs/api/features/04_ops_005_mutation_generation.md](../../docs/api/features/04_ops_005_mutation_generation.md)
- [docs/api/features/04_ops_006_sequence_concatenation.en.md](../../docs/api/features/04_ops_006_sequence_concatenation.en.md)
- [docs/api/features/04_ops_006_sequence_concatenation.md](../../docs/api/features/04_ops_006_sequence_concatenation.md)
- [docs/api/features/04_ops_008_circular_sequence_operations.en.md](../../docs/api/features/04_ops_008_circular_sequence_operations.en.md)
- [docs/api/features/04_ops_008_circular_sequence_operations.md](../../docs/api/features/04_ops_008_circular_sequence_operations.md)
- [docs/api/features/05_all_descriptors.en.md](../../docs/api/features/05_all_descriptors.en.md)
- [docs/api/features/05_all_descriptors.md](../../docs/api/features/05_all_descriptors.md)
- [docs/api/features/06_patterns.en.md](../../docs/api/features/06_patterns.en.md)
- [docs/api/features/06_patterns.md](../../docs/api/features/06_patterns.md)
- [docs/api/features/07_physicochemical.en.md](../../docs/api/features/07_physicochemical.en.md)
- [docs/api/features/07_physicochemical.md](../../docs/api/features/07_physicochemical.md)
- [docs/api/features/08_fingerprints.en.md](../../docs/api/features/08_fingerprints.en.md)
- [docs/api/features/08_fingerprints.md](../../docs/api/features/08_fingerprints.md)
- [docs/api/features/09_alignment.en.md](../../docs/api/features/09_alignment.en.md)
- [docs/api/features/09_alignment.md](../../docs/api/features/09_alignment.md)
- [docs/api/features/09_search.en.md](../../docs/api/features/09_search.en.md)
- [docs/api/features/09_search.md](../../docs/api/features/09_search.md)
- [docs/api/features/09_similarity_alignment.en.md](../../docs/api/features/09_similarity_alignment.en.md)
- [docs/api/features/09_similarity_alignment.md](../../docs/api/features/09_similarity_alignment.md)
- [docs/api/features/10_clustering.en.md](../../docs/api/features/10_clustering.en.md)
- [docs/api/features/10_clustering.md](../../docs/api/features/10_clustering.md)
- [docs/api/features/10_datasets.en.md](../../docs/api/features/10_datasets.en.md)
- [docs/api/features/10_datasets.md](../../docs/api/features/10_datasets.md)
- [docs/api/features/10_deduplication.en.md](../../docs/api/features/10_deduplication.en.md)
- [docs/api/features/10_deduplication.md](../../docs/api/features/10_deduplication.md)
- [docs/api/features/12_evaluation.en.md](../../docs/api/features/12_evaluation.en.md)
- [docs/api/features/12_evaluation.md](../../docs/api/features/12_evaluation.md)
- [docs/api/features/13_visualization.en.md](../../docs/api/features/13_visualization.en.md)
- [docs/api/features/13_visualization.md](../../docs/api/features/13_visualization.md)
- [docs/api/features/14_engineering.en.md](../../docs/api/features/14_engineering.en.md)
- [docs/api/features/14_engineering.md](../../docs/api/features/14_engineering.md)
- [docs/api/features/15_download.en.md](../../docs/api/features/15_download.en.md)
- [docs/api/features/15_download.md](../../docs/api/features/15_download.md)
- [docs/api/features/17_sequence_chunking.en.md](../../docs/api/features/17_sequence_chunking.en.md)
- [docs/api/features/17_sequence_chunking.md](../../docs/api/features/17_sequence_chunking.md)
- [docs/api/features/function_tree.en.md](../../docs/api/features/function_tree.en.md)
- [docs/api/features/function_tree.md](../../docs/api/features/function_tree.md)
- [docs/api/index.en.md](../../docs/api/index.en.md)
- [docs/api/index.md](../../docs/api/index.md)
- [docs/changelog.en.md](../../docs/changelog.en.md)
- [docs/changelog.md](../../docs/changelog.md)
- [docs/faq.en.md](../../docs/faq.en.md)
- [docs/faq.md](../../docs/faq.md)
- [docs/faq/api/04_ops_001_sequence_direction.en.md](../../docs/faq/api/04_ops_001_sequence_direction.en.md)
- [docs/faq/api/04_ops_001_sequence_direction.md](../../docs/faq/api/04_ops_001_sequence_direction.md)
- [docs/faq/api/04_ops_002_transcription_translation.en.md](../../docs/faq/api/04_ops_002_transcription_translation.en.md)
- [docs/faq/api/04_ops_002_transcription_translation.md](../../docs/faq/api/04_ops_002_transcription_translation.md)
- [docs/faq/api/04_ops_004_sequence_editing.en.md](../../docs/faq/api/04_ops_004_sequence_editing.en.md)
- [docs/faq/api/04_ops_004_sequence_editing.md](../../docs/faq/api/04_ops_004_sequence_editing.md)
- [docs/faq/api/04_ops_006_sequence_concatenation.en.md](../../docs/faq/api/04_ops_006_sequence_concatenation.en.md)
- [docs/faq/api/04_ops_006_sequence_concatenation.md](../../docs/faq/api/04_ops_006_sequence_concatenation.md)
- [docs/faq/api/index.en.md](../../docs/faq/api/index.en.md)
- [docs/faq/api/index.md](../../docs/faq/api/index.md)
- [docs/faq/tutorial_basics.en.md](../../docs/faq/tutorial_basics.en.md)
- [docs/faq/tutorial_basics.md](../../docs/faq/tutorial_basics.md)
- [docs/faq/tutorial_integrations.en.md](../../docs/faq/tutorial_integrations.en.md)
- [docs/faq/tutorial_integrations.md](../../docs/faq/tutorial_integrations.md)
- [mkdocs.yml](../../mkdocs.yml)
- [others/docs/documentation_label_changes.md](../../others/docs/documentation_label_changes.md)
- [src/dnakit/evaluation/frechet.py](../../src/dnakit/evaluation/frechet.py)

## 执行示例

以下命令在项目根目录执行，使用当前工作区可用的 Python 3.10 与文档依赖。

设置文档依赖与源码查找路径：

```bash
export PYTHONPATH="$PWD/.venv/lib/python3.10/site-packages:$PWD/src"
```

检查自动生成文档与当前源码一致：

```bash
/opt/conda/envs/qwen3-vl/bin/python others/docs/build_native_api_docs.py --check
```

严格构建中英文站点，输出到 `site/`：

```bash
/opt/conda/envs/qwen3-vl/bin/python -m mkdocs build --strict --config-file mkdocs.yml
```

