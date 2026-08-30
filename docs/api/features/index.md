# DNAKit 功能模块索引

## 1) 数据准备

- [核心数据对象](01_core_objects.md)（8 项）
- [文件读写](02_io_data.md)（4 类；大文件能力作为 IO-001 的可选模式）
- [下载](15_download.md)（常见物种名称、参考基因组 FASTA 和公共数据库数据）
- [数据查询](22_database_query.md)（公共数据库查询能力）
- [合法性检查](03_validation.md)（统一检查入口）

## 2) 数据处理

- [DNA 序列标准化](03_standardization.md)（STD-001 字符标准化）
- **序列内操作**
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
  - [OPS-010 序列切分](17_sequence_chunking.md)
- [序列去重](10_deduplication.md)（DATA-001–006）
- [序列聚类](10_clustering.md)（DATA-007–011、DATA-027 神经网络聚类）
- [数据集整理与划分](10_datasets.md)（DATA-012–018、DATA-023，含自定义 label 划分）

## 3) 数据分析

- **DNA 描述符、表征与指纹**
  - [DNA 描述符](05_all_descriptors.md)
  - [序列表征](08_fingerprints.md)
    - [神经网络表征](08_fingerprints.md#neural-representations)（DATA-027：11 种 DNA 基础模型 rep）
  - [DNA 指纹](08_feature_engineering.md)
- **序列搜索**
  - [通用搜索](09_search.md)（SIM-001、SIM-002、SIM-004、SIM-005、SIM-014、SIM-015）
  - [序列功能搜索](06_patterns.md)（PAT-001、PAT-003–007、PAT-009–012）
- [理化性质](07_physicochemical.md)
- [双链热力学扩展](19_duplex_thermodynamics.md)
- [二级结构性质](20_secondary_structure.md)
- [三维结构与力学性质](21_structure3d.md)
- [换算](18_optics_concentration.md)

## 4) 数据评价

- [常用评价指标](12_evaluation.md)（EVAL-001、EVAL-002、EVAL-005–008、EVAL-016–018）
- [相似度计算](09_similarity_alignment.md)（SIM-010–013、SIM-016 和参考/分布相似度）
- [序列距离与比对](09_alignment.md)（SIM-006–008）

## 5) 可视化

- [可视化](13_visualization.md)（VIZ-001–009）

## 6) 工程化与扩展

- [后端、性能与可复现性](14_engineering.md)（ENG-001–017）
