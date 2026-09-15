# DNAKit 功能树

> 更新日期：2026-09-13。按当前工作区 `src/dnakit` 的公开导出、模型/编码目录、CLI、Agent 目录和功能页核对，包含尚未提交的实现。覆盖六大功能模块，并在下方展开生物信息学 Python API。编号沿用功能页或需求追踪编号，不表示连续的功能计数；这里按用户能力归类，不逐一列出配置类、结果类型、兼容别名或底层辅助函数。可选依赖、外部程序和服务要求见对应页面；列入树中表示已有接口与实现，不表示本机已安装或完成所有外部后端的实测。

详情入口：[功能模块索引](index.md) · 大模型工具：[MCP](../../llm_tools/mcp.md) / [Python API](../../llm_tools/python_api.md) · [自定义性质训练与预测](23_deep_learning_property_prediction.md) · [基因名称与区间注释](26_gene_annotation.md) · [生物信息学 API](27_bioinformatics.md) · [完整 API 参考](../index.md)。

## 新增公开函数分支 (2026-09-13) {#native-apis}

本轮新增 **484 个公开函数**，下列分类中均已逐项列出签名、参数和返回值。按模块统计为数据库 455、ENCODE 20、工作流 6、序列检索与下载 3。下表按用途分类，每个函数只计一次；函数数不等于独立科学算法数。

| 用途分支与详细文档 | 函数数 |
| --- | ---: |
| [基因、序列、测序数据与下载](../native/sequences.md) | 41 |
| [遗传变异、GWAS 与变异评分](../native/variants.md) | 80 |
| [疾病、表型与本体](../native/phenotypes.md) | 45 |
| [表达、表观调控与肿瘤数据](../native/regulation.md) | 65 |
| [模式生物、物种与微生物组](../native/organisms.md) | 49 |
| [通路、互作、蛋白与结构](../native/pathways.md) | 64 |
| [靶点、药物与关联证据](../native/targets.md) | 62 |
| [文献、预印本与临床研究检索](../native/literature.md) | 17 |
| [已训练模型的预测接口](../native/models.md) | 14 |
| [序列计算、VCF、克隆与遗传模拟](../native/analysis.md) | 15 |
| [ENCODE 查询、下载与数据管理](../native/encode.md) | 20 |
| [接口发现、工作流与报告](../native/execution.md) | 12 |
| **[484 / 484](../native/index.md)** | **484** |

[完整索引与调用示例](../native/usage.md)

## 核心功能树


```
DNAKit 功能总览（快速入门教程 + 独立工具与 API 扩展）
├── 功能树总览
├── 功能模块索引
├── 1. 数据准备
│   ├── （1）核心数据对象
│   │   ├── 普通用户统一入口：dnakit.DNA(...)
│   │   │   ├── 一条序列 + 可选 ID / topology / metadata / feature
│   │   │   ├── 多条序列或逐条 record mapping
│   │   │   └── 下标与切片仍返回 DNA
│   │   ├── 1. CORE-001 DNA序列对象
│   │   ├── 2. CORE-002 DNA记录对象
│   │   ├── 3. CORE-003 DNA数据集对象
│   │   ├── 4. CORE-004 特征对象
│   │   ├── 5. CORE-005 Gap对象
│   │   ├── 6. CORE-006 序列类型声明
│   │   └── 7. CORE-007 坐标系统
│   ├── （2）文件读写
│   │   ├── 1. IO-001 序列格式
│   │   │   ├── read(..., mode="dna")：普通单条/多条读取
│   │   │   ├── read(..., mode="stream")：大文件流式读取
│   │   │   ├── write(...)：普通数据或流式数据使用同一写出入口
│   │   │   ├── 批次迭代：iter_chunks（高级模式沿用 IO-005 追踪编号）
│   │   │   └── FASTA/FASTQ 索引构建、加载与随机提取：build_fasta_index / build_fastq_index / load_fasta_index / load_fastq_index / fetch
│   │   ├── 2. IO-002 注释格式：BED / GFF3 / AGP 读写与 GenBank feature
│   │   ├── 3. IO-003 表格格式
│   │   ├── 4. IO-004 压缩文件
│   │   └── 元数据与结果表：筛选、合并、字段选择、校验、CSV/TSV/JSON/Parquet 导入导出
│   ├── （3）下载
│   │   ├── 常见物种名称（下载输入辅助）
│   │   ├── 1. DBD-001 完整基因组
│   │   ├── 2. DBD-002 部分基因组数据包
│   │   ├── 3. DBD-003 区域序列
│   │   ├── 4. DBD-004 基因组注释
│   │   ├── 5. DBD-005 基因序列
│   │   ├── 6. DBD-006 蛋白序列
│   │   ├── 7. DBD-007 组装信息
│   │   ├── 8. DBD-008 分类信息
│   │   ├── 9. DBD-009 基因元数据
│   │   ├── 10. DBD-010 病毒数据包
│   │   ├── 11. DBD-011 原始测序数据
│   │   ├── 12. DBD-012 比对数据
│   │   ├── 13. DBD-013 测序元数据
│   │   ├── 14. DBD-014 变异数据
│   │   ├── 15. DBD-015 ClinVar 数据
│   │   ├── 16. DBD-016 表达数据
│   │   ├── 17. DBD-017 调控数据
│   │   ├── 18. DBD-018 重复序列数据
│   │   ├── 19. DBD-019 保守性数据
│   │   ├── 20. DBD-020 多物种比对
│   │   ├── 21. DBD-021 坐标转换文件
│   │   ├── 通用链接/文件下载、校验与来源清单：download_url / download_file / dataset
│   │   └── 外部索引构建：build_index（makeblastdb / BWA / Bowtie2 / minimap2 / samtools faidx、dict）
│   ├── （4）数据查询
│   │   ├── 1. 基础查询：基因、序列、accession、组装、分类、项目、样本、文献与数据库版本
│   │   ├── 2. 坐标、转录本与区域注释：附近基因、ID/坐标转换、Ensembl lookup
│   │   ├── 3. 测序、表达、调控、变异、病毒、同源/直系同源与比较基因组
│   │   ├── ENA / ENCODE / UCSC 专用查询与轨道、区域序列访问
│   │   └── 远程 BLAST：submit_blast / blast_status / wait_for_blast / blast_results；identify / novelty 高层查询
│   └── （5）合法性检查
│       └── 1. 统一序列、记录与数据集合法性检查
├── 2. 数据处理
│   ├── （1）DNA序列标准化
│   │   ├── 1. STD-001 字符标准化
│   │   ├── 连续大写 N 段转为已知长度 Gap：normalize_gaps
│   │   └── 按 AGP 组装序列及 Gap、来源记录：sequence_from_agp
│   ├── （2）序列内操作
│   │   ├── 统一无后缀入口：insert / delete / substitute / mask / trim / reverse_complement / rotate / canonical_origin
│   │   ├── OPS-001 序列补全
│   │   │   ├── 1. OPS-001.1 逆序
│   │   │   ├── 2. OPS-001.2 互补
│   │   │   └── 3. OPS-001.3 反向互补
│   │   ├── OPS-002 转录与翻译
│   │   │   ├── 1. OPS-002.1 转录
│   │   │   └── 2. OPS-002.2 翻译
│   │   ├── OPS-003 子序列提取
│   │   ├── OPS-004 序列编辑
│   │   ├── OPS-005 序列生成
│   │   │   ├── 1. OPS-005.1 突变生成
│   │   │   ├── 2. OPS-005.2 插入生成
│   │   │   ├── 3. OPS-005.3 删除生成
│   │   │   ├── 4. OPS-005.4 类进化算法操作
│   │   │   ├── 5. OPS-005.5 片段重排
│   │   │   ├── 6. OPS-005.6 k-mer 随机打乱
│   │   │   └── 7. OPS-005.7 多序列重组 / 交叉
│   │   ├── OPS-006 序列拼接
│   │   │   ├── 1. 普通序列拼接
│   │   │   └── 2. 去重叠拼接
│   │   ├── OPS-007
│   │   │   ├── 1. Trimming（修剪）
│   │   │   └── 2. Masking（掩蔽）
│   │   ├── OPS-008 环状序列操作
│   │   └── OPS-010 序列切分
│   │       ├── 1. 固定长度、不重叠
│   │       ├── 2. 固定长度、重叠滑窗
│   │       ├── 3. 随机区间
│   │       ├── 4. 多尺度
│   │       └── 5. 先短后长
│   ├── （3）序列去重
│   │   ├── 1. DATA-001 · 标准去重
│   │   ├── 2. DATA-002 · 反向互补去重
│   │   ├── 3. DATA-003 · 环状等价去重
│   │   ├── 4. DATA-004 · IUPAC-aware去重
│   │   └── 5. DATA-005 · 近似去重
│   ├── （4）序列聚类
│   │   ├── 1. DATA-007 · Identity聚类
│   │   ├── 2. DATA-008 · k-mer聚类
│   │   ├── 3. DATA-009 · 指纹聚类
│   │   ├── 4. DATA-010 · 层次聚类
│   │   ├── 5. DATA-011 · 代表序列选择
│   │   └── 6. DATA-027 · 神经网络聚类
│   ├── （5）数据集整理与划分
│   │   ├── 1. DATA-012 · 随机与稳定哈希划分
│   │   ├── 2. DATA-013 · 分层随机划分
│   │   ├── 3. DATA-014 · 相似度划分
│   │   ├── 4. DATA-015 · Cluster split
│   │   ├── 5. DATA-016 · 物种划分
│   │   ├── 6. DATA-017 · 染色体划分
│   │   ├── 7. DATA-018 · 个体划分
│   │   ├── 8. 按自定义 label 划分
│   │   ├── 9. DATA-023 · 泄漏检测
│   │   ├── 时间划分与多约束联合划分：temporal_split / joint_split
│   │   ├── 按物种、染色体或自定义元数据排除记录
│   │   └── 划分质量评价：evaluate_split_quality
│   └── （6）分子生物学模拟与设计（API 扩展）
│       ├── 限制酶消化、末端分类、末端兼容性检查与片段连接
│       ├── 引物性质、引物匹配与 PCR 模拟
│       ├── 引物设计请求与 Primer3 CLI 设计适配器
│       ├── 序列级组装模拟：simulate_assembly
│       ├── CRISPR 候选与脱靶扫描：scan_crispr_candidates
│       ├── 密码子优化、规则约束优化与突变文库生成
│       └── 更完整的约束优化、双链克隆与 Sanger 验证：见下方 bioinformatics 扩展
├── 3. 数据分析
│   ├── （1）DNA描述符+表征+指纹
│   │   ├── DNA描述符
│   │   │   ├── 1. STD-005 模糊碱基统计
│   │   │   ├── 2. DESC-001 长度特征
│   │   │   ├── 3. DESC-002 碱基组成
│   │   │   ├── 4. DESC-003 GC/AT特征
│   │   │   ├── 5. DESC-004 碱基偏斜
│   │   │   ├── 6. DESC-005 CpG特征
│   │   │   ├── 7. DESC-006 k-mer统计
│   │   │   ├── 8. DESC-007 序列熵
│   │   │   ├── 9. DESC-008 序列复杂度
│   │   │   ├── 10. DESC-009 Homopolymer
│   │   │   ├── 11. DESC-010 重复比例
│   │   │   ├── 12. DESC-011 窗口描述符
│   │   │   ├── 13. DESC-012 密码子统计
│   │   │   ├── 14. 全部描述符计算：all_descriptors（固定 schema 的 240 个输出字段）
│   │   │   └── 描述符 schema 查询与带来源的二核苷酸理化表加载
│   │   ├── DNA表征
│   │   │   ├── 1. FP-001 整数编码
│   │   │   ├── 2. FP-002 One-hot编码
│   │   │   ├── 3. FP-003 k-mer 特征
│   │   │   ├── 4. FP-005 k-mer Sketch（MinHash/FracMinHash）
│   │   │   └── 5. 神经网络表征
│   │   │       ├── 1. AlphaGenome
│   │   │       ├── 2. Caduceus
│   │   │       ├── 3. DNABERT-2
│   │   │       ├── 4. Enformer
│   │   │       ├── 5. Evo 2
│   │   │       ├── 6. GENERator
│   │   │       ├── 7. GROVER
│   │   │       ├── 8. HyenaDNA
│   │   │       ├── 9. JanusDNA
│   │   │       ├── 10. LucaOne（默认）
│   │   │       └── 11. Nucleotide Transformer v2
│   │   └── DNA指纹
│   │       ├── 1. Hashed k-mer 位指纹
│   │       ├── 2. Panel 存在性指纹
│   │       ├── motif / 限制酶位点：motif_fingerprint / restriction_fingerprint
│   │       ├── GC 空间分布 / 重复 / 编码区：gc_spatial_fingerprint / repeat_fingerprint / coding_fingerprint
│   │       ├── 热力学 / 多尺度 / 混合：thermodynamic_fingerprint / multiscale_fingerprint / hybrid_fingerprint
│   │       ├── 特征预处理：fit_preprocessor / transform（缺失值、标准化、归一化与低方差过滤）
│   │       └── 44 种传统 DNA 编码、特征选择、降维与聚类：见下方 encoding / features
│   ├── （2）序列搜索
│   │   ├── 通用搜索
│   │   │   ├── 1. SIM-001 · Exact search
│   │   │   ├── 2. SIM-002 · Subsequence search
│   │   │   ├── 3. SIM-004 · Approximate matching
│   │   │   ├── 4. SIM-005 · 反向互补搜索
│   │   │   ├── 5. SIM-014 · 最近邻搜索
│   │   │   ├── 6. SIM-015 · 数据库索引
│   │   │   └── Sketch 索引构建、保存与重载：build_sketch_index / save_sketch_index / load_sketch_index
│   │   └── 序列功能搜索
│   │       ├── 1. PAT-001 / SIM-003 motif搜索
│   │       ├── 2. PAT-003 起止密码子
│   │       ├── 3. PAT-004 启动子motif
│   │       ├── 4. PAT-005 TF motif
│   │       ├── 5. PAT-006 限制酶位点
│   │       ├── 6. PAT-007 CRISPR PAM
│   │       ├── 7. PAT-009 回文序列
│   │       ├── 8. PAT-010 倒置重复
│   │       ├── 9. PAT-011 串联重复
│   │       ├── 10. PAT-012 微卫星
│   │       ├── 开放阅读框：scan_orfs
│   │       └── CpG 岛与低复杂度区域：find_cpg_islands / find_low_complexity_regions
│   ├── （3）理化性质
│   │   ├── 1. THERMO-001 分子量
│   │   ├── 2. THERMO-014 260 nm 消光系数
│   │   ├── 3. THERMO-002 熔解温度Tm
│   │   ├── 4. THERMO-003 盐浓度修正
│   │   ├── 5. THERMO-012 局部熔解特征
│   │   ├── 6. EVAL-013 合成风险
│   │   └── 热力学性质
│   │       ├── 1. THERMO-004 热力学参数
│   │       ├── 2. THERMO-005 Nearest-neighbor
│   │       ├── 3. THERMO-006 Duplex stability
│   │       ├── 4. THERMO-007 碱基堆积
│   │       ├── 5. THERMO-008 Hairpin
│   │       ├── 6. THERMO-009 Self-dimer
│   │       └── 7. THERMO-010 Heterodimer
│   ├── （4）双链热力学扩展
│   │   ├── 1. 完全互补双链的 ΔH、ΔS、ΔG 和 Tm
│   │   ├── 2. 统一双链稳定性结果
│   │   ├── 3. 相邻碱基对步骤贡献
│   │   ├── 4. 条件与 Na⁺/K⁺ 单价盐
│   │   ├── 5. Ka、Kd 和双链比例
│   │   ├── 6. 理论熔解曲线
│   │   ├── 7. 5′/3′ 末端稳定性
│   │   ├── 8. DMSO 和甲酰胺经验修正
│   │   └── 9. Primer3 CLI 的 Mg²⁺、dNTP、mismatch 和 dangling end
│   ├── （5）二级结构性质
│   │   ├── 1. Dot-bracket 结构解析
│   │   ├── 2. 配对概率与窗口可接近性
│   │   ├── 3. 目标结构的集合缺陷
│   │   ├── 4. 目标结构热力学概率
│   │   ├── 5. NUPACK 被动可用性检查
│   │   ├── 6. NUPACK 单复合物集合分析
│   │   └── 7. NUPACK tube 多复合物平衡
│   ├── （6）三维结构与力学性质
│   │   ├── 1. 读取单个 PDB 模型
│   │   ├── 2. 读取 PDB 多模型集合
│   │   ├── 3. 显式坐标几何分析
│   │   ├── 4. NMR/多模型 RMSF 柔性
│   │   ├── 5. 3DNA bp_step.par 标准参数解析
│   │   └── 6. DSSR JSON 摘要解析
│   ├── （7）换算
│   │   ├── 1. 260 nm 单链摩尔消光系数
│   │   ├── 2. 单链/双链理论光学性质
│   │   ├── 3. 1 OD260 对应的 nmol 和质量
│   │   ├── 4. A260 到摩尔浓度、质量浓度和总量
│   │   ├── 5. 浓度、物质的量和质量互换
│   │   └── 6. 染料和修饰基团的显式校正
│   ├── （8）性质预测：深度学习与机器学习
│   │   ├── 用户自定义分类与回归
│   │   │   ├── 多格式文件检查、自动或显式字段映射与 DataReader 扩展
│   │   │   ├── 二分类、多分类、多标签、多输出分类与单/多目标回归
│   │   │   ├── 全量微调 full 与冻结表征线性探测 linear_probe
│   │   │   ├── 机器学习 machine_learning：37 类算法，序列/embedding/数值特征输入
│   │   │   │   ├── 分类与回归：ridge / sgd / svm / linear_svm / nu_svm / knn / decision_tree / extra_tree
│   │   │   │   ├── 分类与回归：random_forest / extra_trees / gradient_boosting / hist_gradient_boosting / adaboost / bagging
│   │   │   │   ├── 分类与回归：gaussian_process / mlp / dummy / xgboost / lightgbm / catboost
│   │   │   │   ├── 仅分类：logistic_regression / perceptron / gaussian_nb / multinomial_nb / bernoulli_nb / complement_nb / lda / qda
│   │   │   │   └── 仅回归：linear_regression / lasso / elastic_net / bayesian_ridge / ard / huber / ransac / theil_sen / kernel_ridge
│   │   │   ├── 数据划分、训练参数、分类/回归指标与训练产物导出
│   │   │   ├── 训练与重新加载推理：train_custom_property_model / predict_custom_properties
│   │   │   └── 可选 LLM 字段分析：服务商预设、结构/有限样本分享与声明式映射
│   │   ├── 官方任务 1–11. 基因组功能轨道及 SNV 轨道差值
│   │   │   ├── RNA-seq、CAGE、PRO-cap、ATAC-seq、DNase-seq
│   │   │   ├── 组蛋白修饰、转录因子结合
│   │   │   └── 剪接位点、剪接使用比例、剪接连接计数、三维 DNA 接触图
│   │   ├── 官方任务 12–13. 人/鼠调控轨道及 SNV 轨道差值
│   │   ├── 官方任务 14. 14 类单碱基基因组分割
│   │   │   ├── 14.1) 蛋白质编码基因（protein-coding gene）
│   │   │   ├── 14.2) 长链非编码 RNA（lncRNA）
│   │   │   ├── 14.3) 外显子（exon）
│   │   │   ├── 14.4) 内含子（intron）
│   │   │   ├── 14.5) 剪接供体位点（splice donor）
│   │   │   ├── 14.6) 剪接受体位点（splice acceptor）
│   │   │   ├── 14.7) 5′ 非翻译区（5′ UTR）
│   │   │   ├── 14.8) 3′ 非翻译区（3′ UTR）
│   │   │   ├── 14.9) CTCF 结合区域（CTCF-bound）
│   │   │   ├── 14.10) 多聚腺苷酸化信号（polyA signal）
│   │   │   ├── 14.11) 组织特异性增强子（tissue-specific enhancer）
│   │   │   ├── 14.12) 组织不变性增强子（tissue-invariant enhancer）
│   │   │   ├── 14.13) 组织特异性启动子（tissue-specific promoter）
│   │   │   └── 14.14) 组织不变性启动子（tissue-invariant promoter）
│   │   ├── 官方任务 15–16. 长上下文零样本变异评分与外显子概率
│   │   ├── 官方任务 17. 等位基因条件概率变异评分
│   │   ├── 官方任务 18–27. 10 个核酸/蛋白已训练下游任务
│   │   │   ├── 核酸–蛋白中心法则关系、SupKTax、GenusTax、物种分类
│   │   │   ├── 原核蛋白亚细胞定位、蛋白稳定性、ncRNA 家族
│   │   │   └── 甲型流感抗原关系、蛋白–蛋白互作、ncRNA–蛋白互作
│   │   ├── 官方任务 28–45. NT Revised 18 个完整微调分类任务
│   │   │   ├── H2A.Z、H3K27ac、H3K27me3、H3K36me3、H3K4me1/2/3、H3K9ac/me3、H4K20me1
│   │   │   ├── 增强子判别与类型分类；启动子、非 TATA / TATA 启动子
│   │   │   └── 剪接受体、剪接位点类型与剪接供体
│   │   └── 官方任务 46–54. Genomic Benchmarks 9 个完整微调分类任务
│   │       ├── 编码区/基因间区、人类/线虫、果蝇增强子、小鼠增强子
│   │       ├── 人类 Cohn / Ensembl 增强子、人类调控元件类型
│   │       └── 人类非 TATA 启动子、人类开放染色质区域
│   ├── （9）可选生物信息功能
│   │   ├── 1. Ensembl VEP 与变异标识转换
│   │   ├── 2. ClinVar、dbSNP 与 gnomAD 注释
│   │   ├── 3. Nei–Gojobori dN/dS
│   │   ├── 4. Golden Gate 设计与反应组装
│   │   └── MSA、系统发育、NGS、统计与分子设计：见下方 bioinformatics 扩展
│   └── （10）基因名称与区间注释
│       ├── 基因组 DNA + 物种 → annotate_genes / GeneAnnotationConfig
│       ├── 默认 RefSeq 组装解析、显式组装或本地 FASTA + GFF3
│       ├── minimap2 参考比对与已有基因注释映射（非从头基因预测）
│       ├── 基因名称、ID、方向、重叠基因与 gene/intergenic/unknown 分段
│       ├── Python 0-based 半开区间；JSON/CLI 1-based 闭区间
│       └── 组装与注释固定缓存、文件校验、离线复用及来源记录
├── 4. 数据评价
│   ├── （1）常用评价指标
│   │   ├── 1. EVAL-001 Validity
│   │   ├── 2. EVAL-005 Uniqueness
│   │   ├── 3. EVAL-006 Diversity（相似度 / 平均两两 Levenshtein）
│   │   ├── 4. EVAL-008 Novelty（相似度 / 平均最近参考 Levenshtein）
│   │   ├── 5. EVAL-016 Fréchet DNA distance
│   │   ├── 6. EVAL-017 Frag
│   │   ├── 7. EVAL-018 SNN
│   │   ├── 8. EVAL-002 Ambiguity
│   │   ├── 9. EVAL-007 Redundancy
│   │   ├── 复杂度与质量评价：evaluate_complexity / evaluate_quality
│   │   ├── 参考库构建、最近参考查询与记忆/复制评价
│   │   ├── 综合评分卡：evaluate_scorecard
│   │   └── MMD/KID、PRDC、authenticity、性质分布与多目标评价：见下方 metrics
│   ├── （2）相似度计算
│   │   ├── 1. SIM-010 · k-mer相似度
│   │   ├── 2. SIM-011 · 指纹相似度
│   │   ├── 3. SIM-012 · Sketch相似度
│   │   ├── 4. SIM-013 · Dashing sketch相似度
│   │   ├── 5. EVAL-011 Reference similarity
│   │   ├── 描述符分布相似度：evaluate_distribution_similarity
│   │   └── 向量与矩阵：kmer_vector_similarity / similarity_matrix / pairwise_matrix
│   └── （3）序列距离与比对
│       ├── 1. 序列距离
│       │   ├── 1.1 SIM-006 · Hamming distance
│       │   └── 1.2 SIM-007 · Edit distance
│       └── 2. 成对序列比对（SIM-008）
│           ├── 2.1 全局比对
│           ├── 2.2 局部比对
│           └── 2.3 半全局比对
├── 5. 可视化
│   ├── 1. VIZ-001 序列文字图
│   ├── 2. VIZ-002 位置高亮
│   ├── 3. VIZ-003 样式控制
│   ├── 4. VIZ-004 Gap 显示
│   ├── 5. VIZ-005 线性序列图
│   ├── 6. VIZ-006 环状 DNA 图
│   ├── 7. VIZ-007 自定义字符表示
│   ├── 8. VIZ-008 Alignment 图
│   ├── 相似度矩阵热图：plot_similarity_matrix
│   ├── 自包含、可筛选与展开的 HTML 报告：build_html_report / save_html_report
│   ├── 系统发育树、dotplot、Sanger 色谱与定性凝胶图：见下方 phylogeny / motifs / cloning
│   └── 图片导出：save_svg / save_image（PNG、SVG、JPG，兼容 TIFF/PDF）
└── 6. 工程化与扩展
    ├── 1. ENG-001 统一后端接口
    ├── 2. ENG-002 原生与外部标记
    ├── 3. ENG-003 Python API
    ├── 4. ENG-004 CLI
    │   ├── 环境与后端：info / backends
    │   ├── 数据准备：download-genome / convert / normalize / validate
    │   ├── 分析与比较：describe / fingerprint / search / orfs / compare
    │   ├── 数据集与报告：deduplicate / split / report / workflow
    │   ├── 自定义性质：property-models / inspect-property-data / train-property / predict-property
    │   └── 基因注释：annotate-genes
    ├── 5. ENG-005 配置工作流：load_workflow / run_workflow；产物 manifest 读写、校验与复用
    ├── 6. ENG-006 批量计算：iter_batch / run_batch；进度、取消与逐条错误
    ├── 7. ENG-007 并行计算：有界线程任务与稳定输出顺序
    ├── 8. ENG-008 分块与流式处理
    ├── 9. ENG-009 缓存
    ├── 10. ENG-010 随机种子
    ├── 11. ENG-011 版本追踪
    ├── 12. ENG-012 错误与警告
    ├── 13. ENG-013 Agent 与 MCP
    │   ├── 本地 stdio 服务：dnakit-mcp；可选依赖 dnakit[agent]
    │   ├── 支持 MCP 的客户端接入（Codex、Claude、Cursor、VS Code 等）
    │   ├── 工具统计：dnakit_catalog_stats
    │   ├── 功能分类：list_dnakit_categories
    │   ├── 分页目录：list_dnakit_tools
    │   ├── 按功能搜索：search_dnakit_tools
    │   ├── 参数 Schema 与执行属性：describe_dnakit_tool
    │   ├── JSON 参数调用：call_dnakit_tool
    │   ├── Python 工具注册表与服务构建：default_tool_registry / create_server
    │   └── 有界 JSON 返回、来源追踪、Python-only 标记与 allow_side_effects 授权
    ├── 14. ENG-014 一致性验证
    ├── 16. ENG-016 文档与教程
    ├── 17. ENG-017 Galaxy 包装器：convert / deduplicate / split / report（需 Galaxy 环境）
    ├── Python DAG 与区间集合工作流：见下方 workflows
    ├── 通用大模型调用（独立 Python API：dnakit.llm）
    │   ├── 配置与结果：LLMConfig / LLMResult
    │   ├── 同步/异步：call_llm / acall_llm；自定义提示词与当前步骤数据
    │   ├── 文本/JSON 输入，文本/JSON 对象输出与 token 用量
    │   ├── 本地文本文件：file_path；本地预览 read_llm_file
    │   ├── 编码、文件与请求大小限制；显式指定文件后发送完整文本
    │   ├── 服务商预设：DeepSeek / Gemini / Qwen / Kimi / OpenRouter / OpenAI / custom
    │   ├── 预设查询：available_llm_providers / get_llm_provider
    │   └── 自定义兼容地址、模型、密钥配置及 extra_body 扩展参数
    └── 开发验证：单元测试与一致性检查
```

## 生物信息学 Python API 扩展 {#bioinformatics-extension}

[实现方式、范围与示例](27_bioinformatics.md)。以下 13 个子模块均为当前开发版的 Python API，未加入默认 MCP 工具目录，也没有对应的专用 DNAKit CLI 命令。按子模块导入；原生领域算法可能需要 NumPy/SciPy，格式编解码、通用机器学习及专业组装使用下方注明的可选依赖，外部程序须单独安装。

```text
dnakit.bioinformatics
├── 1. alignment：多序列比对、BLAST 与 Profile HMM
│   ├── 原生 MSA 对象、等宽/ID 校验与 p-distance：MultipleAlignment.from_sequences / distances
│   ├── MSA 格式编解码：read_alignment / write_alignment（Biopython）
│   ├── 多序列比对：align_sequences（外部 MAFFT / MUSCLE 5 / Clustal Omega）
│   ├── 本地核酸 BLAST：build_blast_database / local_blast（外部 BLAST+）
│   ├── BLAST XML/tab 与蛋白 HMMER2/3 结果解析：parse_search_results（Biopython）
│   └── Profile HMM 构建/搜索与核酸表解析：build_profile / search_profile / parse_nhmmer_table（hmmbuild / nhmmer / hmmsearch）
├── 2. phylogeny：系统发育
│   ├── 原生树对象与建树：TreeNode / build_tree（NJ / 按簇大小加权 UPGMA）
│   ├── 树格式编解码：read_tree / write_tree（Newick / Nexus / phyloXML / NeXML；Biopython）
│   ├── 原生树比较、共识与剪枝：robinson_foulds / consensus_tree / prune_tree
│   ├── 原生群落系统发育距离：unifrac（weighted / unweighted，要求定根树）
│   └── 树图：plot_tree（Matplotlib）
├── 3. reads：测序记录、索引与区域
│   ├── FASTA / 多行 FASTQ / gzip 与 Phred 质量：Read / read_fastx / write_fastx
│   ├── 配对与损坏四行 FASTQ 记录恢复：pair_reads / repair_fastq（保留恢复/丢弃统计）
│   ├── 排序、打乱与公共记录：arrange_reads / common_reads
│   ├── 正则重命名与顺序敏感、ID 无关的序列摘要：rename_reads / sequence_digest
│   ├── 区间对象与 BED/GTF 读取：Region / read_regions
│   ├── FASTA 索引及区间提取：index_reference / IndexedReference.fetch（pysam）
│   └── BAM 汇总、覆盖度与参考坐标共识：alignment_summary / coverage_consensus（pysam；共识不写入 indel）
├── 4. ngs：外部程序适配与流程
│   ├── 原始 reads 质控：quality_control（FastQC）
│   ├── 单端/双端比对及私有索引：map_reads（BWA / Bowtie2）
│   ├── BAM 排序、mate 修复、可选重复标记、索引与统计：prepare_bam（samtools）
│   ├── 变异调用、参考规范化与 VCF 索引：call_variants（bcftools）
│   ├── 短读长/混合从头组装：assemble_reads（SPAdes，可选 Nanopore）
│   ├── DNA-seq：dna_seq_workflow（质控→比对→BAM→变异）
│   ├── 单样本 RNA-seq：rna_seq_workflow（HISAT2 / StringTie / featureCounts；差异分析另见 counts）
│   └── 已处理 BAM 的 ChIP-seq：chip_seq_workflow（MACS3 peaks + samtools coverage）
├── 5. encoding：44 种原生 DNA 编码（encode_dna / DNA_ENCODINGS）
│   ├── 词频与近似/非连续词：Kmer / RCKmer / Mismatch / Subsequence
│   ├── 碱基及位置：NAC / ANF / ENAC / binary / PS2 / PS3 / PS4
│   ├── 双核苷酸：CKSNAP / ASDC / DBE / LPDF
│   ├── 化学/EIIP/互信息：NCP / EIIP / PseEIIP / MMI
│   ├── Z-curve：Z_curve_9bit / Z_curve_12bit / Z_curve_36bit / Z_curve_48bit / Z_curve_144bit
│   ├── 理化性质位置编码：DPCP / TPCP
│   ├── 自/互协方差：DAC / DCC / DACC / TAC / TCC / TACC
│   ├── 自相关：NMBroto / Moran / Geary
│   ├── 伪核苷酸组成：PseDNC / PseKNC / PCPseDNC / PCPseTNC / SCPseDNC / SCPseTNC
│   └── 监督编码：PSTNPss / PSTNPds / KNN（显式 training_ids；理化编码需显式性质表）
├── 6. features：特征处理（scikit-learn；MCL 为原生算法）
│   ├── 拟合与持出数据复用：analyze_features / FeatureAnalysis.transform
│   ├── 缩放：standard / minmax / l1 / l2 / none
│   ├── 降维：pca / tsne / lda / latent_dirichlet（t-SNE 无持出 transform）
│   ├── 筛选：chi2 / f_classif / f_regression / mutual_info_classif / mutual_info_regression
│   ├── 筛选：information_gain / pearson / variance；none 仅预处理
│   └── 特征聚类：cluster_features
│       ├── kmeans / minibatch_kmeans / gmm / spectral / affinity_propagation
│       └── mean_shift / dbscan / agglomerative / birch / mcl
├── 7. metrics：生成集合、性质分布与多目标评价
│   ├── 自定义标量性质矩阵：sequence_properties
│   ├── 高斯 MMD² / 多项式 KID：kernel_distance
│   ├── Precision / Recall / Density / Coverage：manifold_metrics
│   ├── 校准后的 clipped density / coverage：clipped_manifold_metrics
│   ├── 相对生成器训练集的真实性：authenticity
│   ├── Pareto 层级与精确超体积：pareto_ranks / hypervolume
│   ├── 多性质区间命中率及直方图 KL：property_hit_rate / property_kl
│   └── KDE 分布一致性及标准误：conformity_score
├── 8. statistics：距离检验、排序与组成统计
│   ├── 距离矩阵校验与置换检验：distance_matrix / distance_test（PERMANOVA / ANOSIM / centroid PERMDISP）
│   ├── 矩阵相关与环境变量子集搜索：mantel / bioenv
│   ├── 主坐标与约束排序：pcoa / ordination（CA / CCA / RDA）
│   ├── 组成变换及逆变换：composition_transform（CLR / ILR / ALR）
│   ├── 原生组成差异分析：ancom
│   ├── alpha 多样性：alpha_diversity（observed / Shannon / Simpson / inverse Simpson）
│   └── beta 距离：beta_diversity（Bray-Curtis / Jaccard / Euclidean / cityblock / cosine）
├── 9. counts：计数数据推断
│   ├── 原始 ANCOM-BC：ancom_bc（OLS / sandwich 方差 / EM 采样偏差校正）
│   ├── 差异表达：differential_expression（未收缩 NB2 GLM / size factors / 显式 contrast / Wald）
│   └── 多重检验：adjust_pvalues（Benjamini–Hochberg）
├── 10. optimization：原生约束优化
│   ├── 规则构建与自定义评价：specification / Specification / SpecEvaluation
│   ├── GC、模式避免/计数、序列/翻译/长度、保留位置：gc / avoid_pattern / pattern_count / sequence / translation / length / avoid_changes
│   ├── 唯一 k-mer、终止/禁用/稀有密码子：unique_kmers / avoid_stop / avoid_codons / avoid_rare_codons
│   ├── 精确互补结构/匹配规则：hairpin / heterodimer / avoid_matches / regions_compatibility
│   ├── Tm 与密码子目标：tm / cai / codon_usage / harmonize
│   ├── 外部同源匹配约束：blast_match_specification（BLAST）
│   ├── 硬约束优先、加权目标、同义替换与圆形序列：optimize_sequence
│   ├── GenBank qualifier 规则导入：specifications_from_record
│   └── 修改、失败位置、搜索预算与 JSON 报告：OptimizationResult.write_report
├── 11. motifs：概率 motif 与 dotplot
│   ├── 原生 counts / PWM / PSSM、背景与简并共识：MotifModel / build_motif
│   ├── 保守 FPR 阈值与双链/圆形扫描：fpr_threshold / scan_motif_model
│   ├── JASPAR / MEME / TRANSFAC 读取：read_motifs（Biopython）
│   └── 正反链精确 k-mer dotplot：dotplot（Matplotlib）
├── 12. cloning：显式双链与克隆
│   ├── 双链、offset、黏/平端与拓扑：Duplex.from_sequence / ends / cut / ligate / reverse_complement
│   ├── 可执行克隆步骤与中间产物：CloningPlan.add / apply / steps
│   ├── Gibson / 同源重叠组装：assemble_duplexes（pydna）
│   ├── ABI 读取、Mott 修剪与色谱图：read_sanger / plot_chromatogram（Biopython / Matplotlib）
│   ├── 双向局部比对与参考坐标共识：sanger_consensus
│   ├── SnapGene / GenBank 构建文件：read_constructs（Biopython；宽容解析可选 pydna）
│   └── 定性凝胶模拟：simulate_gel（Matplotlib）
└── 13. workflows：原生 Python 工作流
    ├── DAG 定义、依赖/环检查与失败传播：Step / run_workflow
    └── 按染色体和链的区间交/并/差：query_regions（intersection / union / difference）
```

## 核对口径与数量

本次对照了公开模块导出、`bioinformatics` 各子模块的公开函数/类方法、`tools/catalog.py`、CLI 命令及 Galaxy 包装器。上方保留原来的六类组织方式；扩展树展开对应能力，不作为额外重复计数。

| 代码目录 | 本次核对结果 | 数量含义 |
| --- | --- | --- |
| `descriptors.DESCRIPTOR_NAMES_V1` | 240 | 固定 schema 的输出字段，包含不同统计量和理化表相关字段，不是 240 个独立 API |
| `bioinformatics.encoding.DNA_ENCODINGS` | 44 | 编码名称，各编码输出维数随参数/序列而异 |
| `representations.available_embedding_models()` | 11 | 表征模型，不代表所有模型都支持全量微调 |
| `predictions.available_machine_learning_models()` | 37 | 算法条目，分类/回归支持范围分别记录 |
| `predictions.available_prediction_tasks()` | 54 | 登记的直接预测任务，不含用户新训练的性质 |

新增扩展的公式、输入条件与限制以[生物信息学 API 说明](27_bioinformatics.md)和代码 docstring 为准。例如基因名称注释依赖参考比对；NB2 差异表达没有 DESeq2/edgeR 的完整收缩流程；约束优化采用有预算的搜索。这些能力均按实际实现范围列出。
