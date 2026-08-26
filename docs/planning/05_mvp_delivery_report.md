# 阶段 1～3 本地交付报告

!!! note "历史 MVP 快照"
    本页的 `planned` 计数和“后续工作”是阶段 3 封板时的历史记录，不代表当前状态。现行 184 项矩阵已无 `planned`；请以[阶段 4/5 最终报告](06_stage4_stage5_delivery_report.md)为准。

## 1. 交付结论

DNAKit 已完成需求审查、架构设计和阶段 3 MVP 基线的本地实现与验证。当前开发版本为 `0.1.0.dev0`，只在 Linux 和独立 Conda 环境中测试，没有推送 GitHub、部署网站或上传 TestPyPI/PyPI。

原始 CSV 共 184 项需求；当前追踪状态为：

| 状态 | 数量 |
|---|---:|
| `complete` | 72 |
| `partial` | 13 |
| `planned` | 91 |
| `conditional` | 8 |

需求 ID 全部唯一。原始需求文件未修改，SHA-256 为 `89453c851657134e04f80b7d040ba97a77f3e31ae34b67727351843bee7dd907`。

## 2. 已实现范围

- 统一不可变对象：`DNASequence`、`DNARecord`、`DNASet`、`DNAFeature`、显式 `Gap` 和坐标对象。
- FASTA、FASTQ、CSV、TSV、JSON、JSONL、gzip、流式读取、原子写入及有损写出保护。
- 标准化、IUPAC/质量/集合验证、结构化错误和参数校验。
- 反向、互补、反向互补、转录、翻译、切片、编辑、突变、拼接、trim 和 mask。
- 基础组成、GC/AT、skew、CpG、k-mer、熵、同聚物、滑窗和密码子描述符。
- 整数、one-hot、k-mer/canonical k-mer 指纹。
- exact/subsequence/反向互补搜索、Hamming、Levenshtein、k-mer/向量相似度和有界矩阵。
- exact/反向互补去重；随机、分层、group 和基础相似度划分。
- 有界串行批处理、稳定顺序、错误收集、进度回调。
- 纯 SVG 序列图、相似度矩阵和原子导出。
- Python API、CLI、本地 MkDocs Material 网站、固定输入只读演示和端到端 Notebook。

## 3. 本地环境与验证

- Conda 环境：`/home/mapengsen/anaconda3/envs/dnakit-dev`
- Python：`3.10.20`
- `PYTHONNOUSERSITE=1`，`site.ENABLE_USER_SITE=False`
- 全量测试：397 passed
- 分支覆盖率：86%
- Ruff lint/format：通过
- mypy strict：通过，104 个源码/测试文件无问题
- compileall、pip check：通过
- MkDocs strict、Notebook nbmake：通过
- wheel 和 sdist 构建、`twine check`、隔离安装及 CLI smoke：通过
- wheel：`dist/dnakit-0.1.0.dev0-py3-none-any.whl`
- sdist：`dist/dnakit-0.1.0.dev0.tar.gz`
- 两个归档均包含 `py.typed`，且不含原始 CSV、`dashing_similarity`、站点、缓存或字节码。

## 4. 代码审查结论

最终独立审查未发现阻断级问题。审查期间已修复路径穿越与多文件事务写出、结构化 I/O 的 feature/Gap 丢失、超长 CSV 边界、不可跨越 Gap、密码子 Gap 相位、生成器资源上限、随机状态与权重审计、metadata 类型碰撞、浮点累计不确定性和编辑距离回溯不一致等问题。

当前结果对象尚未全部统一到架构文档中的完整 `algorithm_version/provenance/parameters/issues` schema，因此 `CORE-008` 保持 `partial`，没有夸大为完成。

## 5. 未完成和条件功能

- 阶段 4：热力学后端、结构分析、近似去重/聚类、多约束划分、novelty/memorization、synthesis-risk、综合指纹和分子生物学模拟。
- 阶段 5：benchmark、论文复现实验、TestPyPI/PyPI 验证和正式发布。
- 批处理 resume/并行、分析类 CLI、统一 workflow manifest、PNG/TIFF/PDF/交互可视化仍为后续工作。
- 项目当前采用 MIT 许可证并提供 `LICENSE`；正式发布仍需完成后续版本和远程发布门禁。
- Primer3 只能作为用户自行安装的条件 adapter；尚未做 Tm/hairpin/dimer 差分验收。
- NUPACK 自动安装和真实差分验收受单独许可/下载条件限制；后续已新增被动探测与显式 adapter，但项目仍不自动安装、下载或打包，也不将其用作网站后端。真实验收需要用户在适用许可范围内独立安装。
- Dashing 本地 GPLv3 副本不进入 wheel/sdist；阶段 3 当时只计划外部 CLI adapter，现已提供要求显式可执行文件的有界 Jaccard/Top-k adapter，但尚无真实科学差分，最终边界以阶段 4/5 报告为准。
- 限制酶、成熟搜索/聚类工具对照和高级数值容差表尚待对应模块实现后完成。

完整状态见[需求追踪矩阵](requirements_traceability.csv)和[条件支持清单](04_conditional_and_unavailable_features.md)。

## 6. 所有新增或修改的持久化项目文件

以下清单不含原始需求 CSV、现有 `dashing_similarity/`、`DNA-Rdkit.md`、`goal.md`、IDE 设置，以及可重建的 `build/`、`site/`、缓存、字节码和 `src/dnakit.egg-info/`。

### 6.1 工程与发布配置

```text
.gitattributes
.gitignore
.github/workflows/ci.yml
CHANGELOG.md
CITATION.cff
CONTRIBUTING.md
MANIFEST.in
README.md
environment-dev.yml
mkdocs.yml
pyproject.toml
```

### 6.2 Python 包

```text
src/dnakit/__init__.py
src/dnakit/_version.py
src/dnakit/batch.py
src/dnakit/exceptions.py
src/dnakit/py.typed
src/dnakit/cli/__init__.py
src/dnakit/cli/__main__.py
src/dnakit/cli/app.py
src/dnakit/core/__init__.py
src/dnakit/core/_json.py
src/dnakit/core/backend_info.py
src/dnakit/core/collection.py
src/dnakit/core/coordinates.py
src/dnakit/core/enums.py
src/dnakit/core/feature.py
src/dnakit/core/gap.py
src/dnakit/core/issues.py
src/dnakit/core/provenance.py
src/dnakit/core/record.py
src/dnakit/core/results.py
src/dnakit/core/sequence.py
src/dnakit/datasets/__init__.py
src/dnakit/datasets/_metadata.py
src/dnakit/datasets/config.py
src/dnakit/datasets/deduplicate.py
src/dnakit/datasets/results.py
src/dnakit/datasets/split.py
src/dnakit/descriptors/__init__.py
src/dnakit/descriptors/_shared.py
src/dnakit/descriptors/basic.py
src/dnakit/descriptors/codon.py
src/dnakit/descriptors/entropy.py
src/dnakit/descriptors/homopolymer.py
src/dnakit/descriptors/kmer.py
src/dnakit/descriptors/results.py
src/dnakit/descriptors/window.py
src/dnakit/fingerprints/__init__.py
src/dnakit/fingerprints/_shared.py
src/dnakit/fingerprints/encoding.py
src/dnakit/fingerprints/kmer.py
src/dnakit/fingerprints/results.py
src/dnakit/io/__init__.py
src/dnakit/io/_formats.py
src/dnakit/io/api.py
src/dnakit/io/config.py
src/dnakit/io/results.py
src/dnakit/io/source.py
src/dnakit/ops/__init__.py
src/dnakit/ops/_common.py
src/dnakit/ops/concat.py
src/dnakit/ops/direction.py
src/dnakit/ops/edit.py
src/dnakit/ops/mutation.py
src/dnakit/ops/translation.py
src/dnakit/similarity/__init__.py
src/dnakit/similarity/_shared.py
src/dnakit/similarity/compare.py
src/dnakit/similarity/distance.py
src/dnakit/similarity/matrix.py
src/dnakit/similarity/results.py
src/dnakit/similarity/search.py
src/dnakit/similarity/vector.py
src/dnakit/standardize/__init__.py
src/dnakit/standardize/_shared.py
src/dnakit/standardize/config.py
src/dnakit/standardize/normalize.py
src/dnakit/standardize/results.py
src/dnakit/standardize/validate.py
src/dnakit/visualization/__init__.py
src/dnakit/visualization/_svg.py
src/dnakit/visualization/config.py
src/dnakit/visualization/export.py
src/dnakit/visualization/matrix.py
src/dnakit/visualization/results.py
src/dnakit/visualization/sequence.py
```

### 6.3 测试

```text
tests/test_cli_workflows.py
tests/test_import.py
tests/unit/test_batch.py
tests/unit/test_core_coordinates.py
tests/unit/test_core_records.py
tests/unit/test_core_results.py
tests/unit/test_core_sequence.py
tests/unit/test_datasets_deduplicate.py
tests/unit/test_datasets_split.py
tests/unit/test_descriptors_basic.py
tests/unit/test_descriptors_kmer_entropy.py
tests/unit/test_descriptors_runs_windows_codons.py
tests/unit/test_fingerprints_encoding.py
tests/unit/test_fingerprints_kmer.py
tests/unit/test_io_config_results.py
tests/unit/test_io_sequence_formats.py
tests/unit/test_io_source.py
tests/unit/test_io_structured.py
tests/unit/test_ops_direction_translation.py
tests/unit/test_ops_edit.py
tests/unit/test_ops_mutation_concat.py
tests/unit/test_similarity_distance.py
tests/unit/test_similarity_search.py
tests/unit/test_similarity_vector_matrix.py
tests/unit/test_standardize_normalize.py
tests/unit/test_standardize_validate.py
tests/unit/test_visualization_charts.py
tests/unit/test_visualization_export.py
tests/unit/test_visualization_sequence.py
```

### 6.4 文档、教程和固定示例

```text
docs/api/index.md
docs/acknowledgements.md
docs/demo/data/fixed_demo.json
docs/demo/index.md
docs/examples/index.md
docs/faq.md
docs/index.md
docs/installation.md
docs/planning/01_requirements_review.md
docs/planning/02_architecture_design.md
docs/planning/03_dependency_validation_development_plan.md
docs/planning/04_conditional_and_unavailable_features.md
docs/planning/05_mvp_delivery_report.md
docs/planning/README.md
docs/planning/requirements_traceability.csv
docs/quickstart.md
docs/stylesheets/extra.css
docs/tutorials/index.md
examples/README.md
examples/fixed_demo.fasta
examples/fixed_demo_expected.json
notebooks/00_skeleton_check.ipynb
```

### 6.5 本地生成工件

```text
dist/dnakit-0.1.0.dev0-py3-none-any.whl
dist/dnakit-0.1.0.dev0.tar.gz
build/
site/
src/dnakit.egg-info/
```

这些生成工件均可由源码重新构建，不应手工编辑。
