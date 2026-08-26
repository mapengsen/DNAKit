# DNAKit 功能树

本页说明功能分类的调整，并按用户实际使用顺序展示新的功能层级。功能编号和各功能页 URL 保持不变，仍可用于需求追踪和 API 查找。

## 分类调整说明

| 原来 | 现在 | 变化 |
| --- | --- | --- |
| 1. 准备数据 | 1. 数据准备 | 基本只是改名 |
| 2. 处理与变换序列 | 2. 数据处理 | 扩大范围，吸收原“整理 DNA 序列”中的去重、聚类和数据划分 |
| 3. 分析 DNA 序列 | 3. 数据分析 | 改名，同时吸收序列搜索 |
| 4. 整理 DNA 序列 | — | 取消一级分类；功能拆入“数据处理”和“数据分析” |
| 5. 评价 | 4. 数据评价 | 包含综合评价指标、相似度计算和序列比对 |
| 6. 可视化与报告 | 5. 可视化 | 取消报告子模块，可视化直接作为一级分类 |
| 7. 其他 | 6. 工程化与扩展 | 改名，明确后端、CLI、工作流、性能和可复现性 |

## 新的功能树

### 1. 数据准备

- [核心数据对象](01_core_objects.md)
  - 普通统一入口：`DNA(...)`（单条、多条、ID、拓扑、metadata、feature）
  - 选择后仍为 `DNA`：`data[index]`、`data[slice]`
  - 高级兼容模型：CORE-001～CORE-008
- [文件读写](02_io_data.md)
  - 统一入口：`read()`、`write()`
  - 普通读取：`read(..., mode="dna")`
  - 大文件读写：仍使用 `read(..., mode="stream")` 和 `write(...)`
  - 高级分块与索引：归入 IO-001；IO-005 仅保留为需求追踪编号
  - 兼容入口：`read_one()`、`read_set()`
- [下载（含常见物种名称）](15_download.md)
- [数据查询](22_database_query.md)
- [合法性检查](03_validation.md)：统一使用 `validate()`，`validate_set()` 仅兼容

### 2. 数据处理

- [DNA 序列标准化](03_standardization.md)
- **序列内操作**
  - 普通编辑、反向互补和环状操作统一使用无后缀名称；`*_record()` 仅用于详细审计
  - [OPS-001 序列补全](04_ops_001_sequence_direction.md)
  - [OPS-002 转录与翻译](04_ops_002_transcription_translation.md)
  - [OPS-003 子序列提取](04_ops_003_subsequence_extraction.md)
  - [OPS-004 序列编辑](04_ops_004_sequence_editing.md)
  - [OPS-005 序列生成](04_ops_005_mutation_generation.md)
  - [OPS-006 序列拼接](04_ops_006_sequence_concatenation.md)
  - **OPS-007**
    - [Trimming（修剪）](04_ops_007_trimming.md)
    - [Masking（掩蔽）](04_ops_007_masking.md)
  - [OPS-008 环状序列操作](04_ops_008_circular_sequence_operations.md)
- [序列切分](17_sequence_chunking.md)
- [序列去重](10_deduplication.md)
- [序列聚类](10_clustering.md)（含 DATA-027 神经网络聚类）
- [数据集整理与划分](10_datasets.md)

### 3. 数据分析

- **DNA 描述符、表征与指纹**
  - [DNA 描述符](05_all_descriptors.md)
  - [序列表征](08_fingerprints.md)
    - [神经网络表征](08_fingerprints.md#neural-representations)（DATA-027：11 种 DNA 基础模型 rep）
  - [DNA 指纹](08_feature_engineering.md)
- **序列搜索**
  - [通用搜索](09_search.md)
  - [序列模式搜索](06_patterns.md)
- [理化性质](07_physicochemical.md)
- [双链热力学扩展](19_duplex_thermodynamics.md)
- [二级结构性质](20_secondary_structure.md)
- [三维结构与力学性质](21_structure3d.md)
- [光学与浓度换算](18_optics_concentration.md)

### 4. 数据评价

- [DNA 综合评价体系](12_evaluation.md)
- [相似度计算](09_similarity_alignment.md)
- [序列距离与比对](09_alignment.md)

### 5. 可视化

- [可视化](13_visualization.md)

### 6. 工程化与扩展

- [后端、性能与可复现性](14_engineering.md)

## 导航原则

- **数据准备**解决“输入从哪里来、是否有效”。
- **数据处理**解决“如何整理和变换数据”。
- **数据分析**解决“如何搜索、描述和解释序列”。
- **数据评价**解决“如何比较序列，以及结果是否满足质量、唯一性和新颖性等指标”。
- **可视化**负责结果展示。
- **工程化与扩展**负责后端、CLI、工作流、缓存、性能和可复现性。
