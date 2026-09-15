# DNAKit 功能模块索引

> 更新日期：2026-09-13。按当前工作区公开 API 核对。完整层级见[功能树总览](function_tree.md)；本页同时列出功能页与公开 API 扩展的详情入口。条目按能力归类，跨分类链接不重复计数。

## 1) 数据准备

本轮新增：[基因、序列与下载](../native/sequences.md)（41 个函数）、[ENCODE 数据管理](../native/encode.md)（20）、[物种与微生物组](../native/organisms.md)（49）、[文献与研究检索](../native/literature.md)（17）。完整参数和返回值均已逐项列出。

- [核心数据对象](01_core_objects.md)（8 项）
- [文件读写](02_io_data.md)（4 类；大文件能力作为 IO-001 的可选模式）
- [索引与元数据 API](../index.md)（FASTA/FASTQ 随机提取、BED/GFF3/AGP、元数据与结果表）
- [下载](15_download.md)（常见物种名称、参考基因组 FASTA 和公共数据库数据）
- [数据查询](22_database_query.md)（公共数据库查询能力）
- [下载与外部索引 API](../index.md)（通用文件下载、BLAST/BWA/Bowtie2/minimap2/samtools 索引）
- [测序记录与区域](27_bioinformatics.md)（FASTX/gzip、配对、四行 FASTQ 恢复、排序、索引与 BAM 共识）
- [合法性检查](03_validation.md)（统一检查入口）

## 2) 数据处理

新增的序列统计、VCF 处理、酶切与组装接口见[序列计算、VCF、克隆与遗传模拟](../native/analysis.md)。该页也从分析栏目链接，函数不重复计数。

- [DNA 序列标准化](03_standardization.md)（STD-001 字符标准化）
- [Gap 与 AGP 标准化 API](../index.md#standardization)（连续大写 N 段转 Gap、AGP 组装与来源记录）
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
- [序列去重](10_deduplication.md)（DATA-001–005）
- [序列聚类](10_clustering.md)（DATA-007–011、DATA-027 神经网络聚类）
- [数据集整理与划分](10_datasets.md)（DATA-012–018、DATA-023，含自定义 label 划分）
- [数据集 API 扩展](../index.md#datasets)（时间/联合划分、元数据排除与划分质量评价）
- [分子生物学模拟与设计](../index.md#molecular-biology)（酶切、连接、PCR、引物设计、组装、CRISPR、密码子与规则优化、突变文库）
- [约束优化与双链克隆](27_bioinformatics.md)（原生规则/加权目标、圆形 DNA、克隆步骤、同源组装与 Sanger 共识）
- [NGS 外部程序流程](27_bioinformatics.md)（FastQC、比对/BAM、变异调用、SPAdes、DNA-seq/RNA-seq/ChIP-seq）

## 3) 数据分析

本轮新增：[变异与 GWAS](../native/variants.md)（80）、[疾病与表型](../native/phenotypes.md)（45）、[表达与调控](../native/regulation.md)（65）、[通路互作与结构](../native/pathways.md)（64）、[靶点与药物证据](../native/targets.md)（62）、[模型预测](../native/models.md)（14）、[序列计算与遗传模拟](../native/analysis.md)（15）。

- **DNA 描述符、表征与指纹**
  - [DNA 描述符](05_all_descriptors.md)
  - [序列表征](08_fingerprints.md)
    - [神经网络表征](08_fingerprints.md#neural-representations)（DATA-027：11 种 DNA 基础模型 rep）
  - [DNA 指纹](08_feature_engineering.md)
  - [高级指纹与预处理](../index.md#fingerprints)（motif、酶切、GC 空间分布、重复、编码区、热力学、多尺度、混合与训练集预处理）
  - [传统 DNA 编码与特征工程](27_bioinformatics.md)（44 种编码、选择、降维与聚类）
- **序列搜索**
  - [通用搜索](09_search.md)（SIM-001、SIM-002、SIM-004、SIM-005、SIM-014、SIM-015）
  - [序列功能搜索](06_patterns.md)（PAT-001、PAT-003–007、PAT-009–012）
- [理化性质](07_physicochemical.md)
- [双链热力学扩展](19_duplex_thermodynamics.md)
- [二级结构性质](20_secondary_structure.md)
- [三维结构与力学性质](21_structure3d.md)
- [换算](18_optics_concentration.md)
- [性质预测：深度学习与机器学习](23_deep_learning_property_prediction.md)（自定义分类/回归、37 类机器学习算法、模型重新加载，以及 54 个官方直接推理功能）
- [自定义训练 API](../index.md#deep-learning-property-prediction)（全量微调、线性探测、机器学习、模型目录与可选 LLM 字段分析）
- [可选生物信息功能](24_optional_bioinformatics.md)（变异注释、dN/dS 和 Golden Gate）
- [基因名称与区间注释](26_gene_annotation.md)（基因组 DNA 的参考比对、已有基因注释映射、重叠区间与缓存复用）
- [MSA、系统发育与概率 motif](27_bioinformatics.md)（多序列比对、BLAST/HMMER、NJ/UPGMA、树比较/共识、PWM/PSSM）
- [统计与计数推断](27_bioinformatics.md)（置换检验、排序、组成变换、ANCOM/ANCOM-BC、NB2 差异表达与 BH 校正）

## 4) 数据评价

- [常用评价指标](12_evaluation.md)（EVAL-001、EVAL-002、EVAL-005–008、EVAL-016–018）
- [相似度计算](09_similarity_alignment.md)（SIM-010–013、EVAL-011 参考相似度）
- [序列距离与比对](09_alignment.md)（SIM-006–008）
- [综合评价 API 扩展](../index.md#evaluation)（质量、复杂度、参考库、记忆/复制、分布相似度与评分卡）
- [生成集合与多目标评价](27_bioinformatics.md)（MMD/KID、PRDC/clipped、authenticity、性质 KL/命中率、conformity、Pareto 与 hypervolume）
- [群落多样性与距离](27_bioinformatics.md)（alpha/beta 多样性、weighted/unweighted UniFrac；区别于 DNA 生成集合多样性）

## 5) 可视化

- [可视化](13_visualization.md)（VIZ-001–008）
- [可视化 API 扩展](../index.md#visualization)（相似度矩阵、HTML 报告与图片导出）
- [生物信息学绘图](27_bioinformatics.md)（树图、dotplot、Sanger 色谱与定性凝胶图）

## 6) 工程化与扩展

新增[接口发现、工作流与报告](../native/execution.md)（12 个函数）。全部 **484 个新增函数**见[完整索引](../native/index.md)，使用方式见[调用示例](../native/usage.md)。

- [后端、性能与可复现性](14_engineering.md)（后端、批处理、有界线程、缓存、工作流与来源记录）
- [CLI 命令](../index.md#cli)（数据处理、分析、工作流、自定义性质训练/预测与基因注释）
- [大模型工具 → MCP](../../llm_tools/mcp.md)（ENG-013；服务配置、八个 MCP 入口、Schema、工具调用与执行授权）
- [大模型工具 → Python API](../../llm_tools/python_api.md)（无需安装 MCP；一句话驱动工具调用、通用大模型调用、同步/异步、自定义提示词、文本/JSON、本地文件、服务商预设与扩展参数）
- [Python DAG 与区间查询](27_bioinformatics.md)（依赖检查、失败传播与区间交/并/差）
- Galaxy 包装器（`galaxy/dnakit`）：convert、deduplicate、split、report，需 Galaxy 环境。

## 生物信息学 Python API

- [13 个子模块的完整树](function_tree.md#bioinformatics-extension)：alignment、phylogeny、reads、ngs、encoding、features、metrics、statistics、counts、optimization、motifs、cloning、workflows。
- [实现范围、依赖与示例](27_bioinformatics.md)：当前开发版仅提供 Python API，未加入默认 MCP 目录；原生算法、可选 Python 库与外部程序分别说明。
