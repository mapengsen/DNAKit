# DNAKit 功能树

## 核心功能树

```
DNAKit 功能总览
├── 功能树总览
├── 1. 数据准备
│   ├── 1) 核心数据对象
│   │   ├── 普通用户统一入口：dnakit.DNA(...)
│   │   │   ├── 一条序列 + 可选 ID / topology / metadata / feature
│   │   │   ├── 多条序列或逐条 record mapping
│   │   │   └── 下标与切片仍返回 DNA
│   │   ├── 1. DNA序列对象
│   │   ├── 2. DNA记录对象
│   │   ├── 3. DNA数据集对象
│   │   ├── 4. 特征对象
│   │   ├── 5. Gap对象
│   │   ├── 6. 序列类型声明
│   │   └── 7. 坐标系统
│   ├── 2) 文件读写
│   │   ├── 1. 序列格式
│   │   │   ├── read(..., mode="dna")：普通单条/多条读取
│   │   │   ├── read(..., mode="stream")：大文件流式读取
│   │   │   ├── write(...)：普通数据或流式数据使用同一写出入口
│   │   │   ├── 批次迭代：iter_chunks
│   │   │   └── FASTA/FASTQ 索引构建、加载与随机提取：build_fasta_index / build_fastq_index / load_fasta_index / load_fastq_index / fetch
│   │   ├── 2. 注释格式：BED / GFF3 / AGP 读写与 GenBank feature
│   │   ├── 3. 表格格式
│   │   ├── 4. 压缩文件
│   │   ├── 元数据与结果表：筛选、合并、字段选择、校验、CSV/TSV/JSON/Parquet 导入导出
│   │   ├── MSA 格式编解码：read_alignment / write_alignment（Biopython）
│   │   ├── 树格式编解码：read_tree / write_tree（Newick / Nexus / phyloXML / NeXML；Biopython）
│   │   ├── FASTA / 多行 FASTQ / gzip 与 Phred 质量：Read / read_fastx / write_fastx
│   │   ├── JASPAR / MEME / TRANSFAC 读取：read_motifs（Biopython）
│   │   ├── SnapGene / GenBank 构建文件：read_constructs（Biopython；宽容解析可选 pydna）
│   │   └── ABI 文件读取与 Mott 修剪：read_sanger
│   ├── 3) 下载
│   │   ├── 常见物种名称（下载输入辅助）
│   │   ├── 1. 完整基因组
│   │   ├── 2. 部分基因组数据包
│   │   ├── 3. 区域序列
│   │   ├── 4. 基因组注释
│   │   ├── 5. 基因序列
│   │   ├── 6. 蛋白序列
│   │   ├── 7. 组装信息
│   │   ├── 8. 分类信息
│   │   ├── 9. 基因元数据
│   │   ├── 10. 病毒数据包
│   │   ├── 11. 原始测序数据
│   │   ├── 12. 比对数据
│   │   ├── 13. 测序元数据
│   │   ├── 14. 变异数据
│   │   ├── 15. ClinVar 数据
│   │   ├── 16. 表达数据
│   │   ├── 17. 调控数据
│   │   ├── 18. 重复序列数据
│   │   ├── 19. 保守性数据
│   │   ├── 20. 多物种比对
│   │   ├── 21. 坐标转换文件
│   │   ├── NCBI／ENA accession 序列下载与 ENCODE 批量下载
│   │   ├── 通用链接/文件下载、校验与来源清单：download_url / download_file / dataset
│   │   └── 外部索引构建：build_index（makeblastdb / BWA / Bowtie2 / minimap2 / samtools faidx、dict）
│   ├── 4) 数据查询
│   │   ├── 1. 查找与调用数据库接口
│   │   ├── 2. 序列、组装与测序数据查询
│   │   ├── 3. 疾病、表型与本体查询
│   │   ├── 4. 表达、调控与肿瘤数据查询
│   │   ├── 5. 模式生物、物种与微生物组查询
│   │   ├── 6. 通路、互作、蛋白与结构查询
│   │   ├── 7. 靶点、药物与关联证据查询
│   │   ├── 8. 文献、预印本与临床研究查询
│   │   ├── 9. ENCODE 实验与文件查询
│   │   └── 10. 遗传变异与人群关联查询
│   │       ├── 1) 变异记录与群体频率查询
│   │       ├── 2) 结构变异查询
│   │       ├── 3) ClinVar 临床意义查询
│   │       ├── 4) 遗传病关联、剂量与可干预性
│   │       ├── 5) 肿瘤变异与治疗证据查询
│   │       ├── 6) GWAS 研究与关联查询
│   │       ├── 7) 队列表型关联与精细定位查询
│   │       └── 8) 人群、样本与连锁查询
│   └── 5) DNA序列合法性检查
│       └── 1. 统一序列、记录与数据集合法性检查
├── 2. 数据处理
│   ├── 1) DNA序列标准化
│   │   ├── 1. 字符标准化
│   │   ├── 连续大写 N 段转为已知长度 Gap：normalize_gaps
│   │   └── 按 AGP 组装序列及 Gap、来源记录：sequence_from_agp
│   ├── 2) 序列内操作
│   │   ├── 统一无后缀入口：insert / delete / substitute / mask / trim / reverse_complement / rotate / canonical_origin
│   │   ├── 序列补全
│   │   │   ├── 1) 逆序
│   │   │   ├── 2) 互补
│   │   │   ├── 3) 反向互补
│   │   │   └── 显式双链反向互补：Duplex.reverse_complement
│   │   ├── 转录与翻译
│   │   │   ├── 1) 转录
│   │   │   └── 2) 翻译
│   │   ├── 子序列提取
│   │   │   ├── 区间对象与 BED/GTF 读取：Region / read_regions
│   │   │   ├── FASTA 索引及区间提取：index_reference / IndexedReference.fetch（pysam）
│   │   │   └── 按染色体和链的区间交/并/差：query_regions（intersection / union / difference）
│   │   ├── 序列编辑
│   │   │   ├── 显式双链、offset、拓扑、黏端/平端与酶切：Duplex.from_sequence / ends / cut
│   │   │   ├── Trimming（修剪）：trim
│   │   │   └── Masking（掩蔽）：mask
│   │   ├── 序列生成
│   │   │   ├── 1) 突变生成
│   │   │   ├── 2) 插入生成
│   │   │   ├── 3) 删除生成
│   │   │   ├── 4) 类进化算法操作
│   │   │   ├── 5) 片段重排
│   │   │   ├── 6) k-mer 随机打乱
│   │   │   ├── 7) 多序列重组 / 交叉
│   │   │   └── 突变文库生成
│   │   ├── 序列拼接与组装
│   │   │   ├── 1) 普通序列拼接
│   │   │   ├── 2) 去重叠拼接
│   │   │   ├── 短读长/混合从头组装：assemble_reads（SPAdes，可选 Nanopore）
│   │   │   ├── 可执行克隆步骤与中间产物：CloningPlan.add / apply / steps
│   │   │   ├── Gibson / 同源重叠组装：assemble_duplexes（pydna）
│   │   │   ├── Golden Gate 部件设计与反应组装：design_golden_gate / assemble_golden_gate
│   │   │   ├── 序列级组装模拟：simulate_assembly
│   │   │   └── 双链末端兼容性与连接：Duplex.ligate
│   │   ├── 环状序列操作
│   │   └── 序列切分
│   │       ├── 1) 固定长度、不重叠
│   │       ├── 2) 固定长度、重叠滑窗
│   │       ├── 3) 随机区间
│   │       ├── 4) 多尺度
│   │       └── 5) 先短后长
│   ├── 3) 序列去重
│   │   ├── 1. 标准去重
│   │   ├── 2. 反向互补去重
│   │   ├── 3. 环状等价去重
│   │   ├── 4. IUPAC-aware去重
│   │   └── 5. 近似去重
│   ├── 4) 序列聚类
│   │   ├── 1. Identity聚类
│   │   ├── 2. k-mer聚类
│   │   ├── 3. 指纹聚类
│   │   ├── 4. 层次聚类
│   │   ├── 5. 代表序列选择
│   │   ├── 6. 神经网络聚类
│   │   ├── 按标签进行多样性采样：diversity_sample
│   │   └── 特征聚类：cluster_features
│   │       ├── kmeans / minibatch_kmeans / gmm / spectral / affinity_propagation
│   │       └── mean_shift / dbscan / agglomerative / birch / mcl
│   ├── 5) 数据集整理与划分
│   │   ├── ENCODE 实验追踪、集合汇总、引用导出与文件来源记录
│   │   ├── 1. 随机与稳定哈希划分
│   │   ├── 2. 分层随机划分
│   │   ├── 3. 相似度划分
│   │   ├── 4. Cluster split
│   │   ├── 5. 物种划分
│   │   ├── 6. 染色体划分
│   │   ├── 7. 个体划分
│   │   ├── 8. 按自定义 label 划分
│   │   ├── 9. 泄漏检测
│   │   ├── 时间划分与多约束联合划分：temporal_split / joint_split
│   │   ├── 按物种、染色体或自定义元数据排除记录
│   │   └── 划分质量评价：evaluate_split_quality
│   ├── 6) 测序数据整理
│   │   ├── 配对与损坏四行 FASTQ 记录恢复：pair_reads / repair_fastq（保留恢复/丢弃统计）
│   │   ├── 排序、打乱与公共记录：arrange_reads / common_reads
│   │   ├── 正则重命名与顺序敏感、ID 无关的序列摘要：rename_reads / sequence_digest
│   │   └── BAM 排序、mate 修复、可选重复标记、索引与统计：prepare_bam（samtools）
│   └── 7) 序列优化
│       ├── 规则构建与自定义评价：specification / Specification / SpecEvaluation
│       ├── GC、模式避免/计数、序列/翻译/长度、保留位置：gc / avoid_pattern / pattern_count / sequence / translation / length / avoid_changes
│       ├── 唯一 k-mer、终止/禁用/稀有密码子：unique_kmers / avoid_stop / avoid_codons / avoid_rare_codons
│       ├── 精确互补结构/匹配规则：hairpin / heterodimer / avoid_matches / regions_compatibility
│       ├── Tm 与密码子目标：tm / cai / codon_usage / harmonize
│       ├── 外部同源匹配约束：blast_match_specification（BLAST）
│       ├── 硬约束优先、加权目标、同义替换与圆形序列：optimize_sequence
│       ├── GenBank qualifier 规则导入：specifications_from_record
│       └── 修改、失败位置、搜索预算与 JSON 报告：OptimizationResult.write_report
├── 3. 数据分析
│   ├── 1) DNA描述符+表征+指纹
│   │   ├── DNA描述符
│   │   │   ├── 1) 模糊碱基统计
│   │   │   ├── 2) 长度特征
│   │   │   ├── 3) 碱基组成
│   │   │   ├── 4) GC/AT特征
│   │   │   ├── 5) 碱基偏斜
│   │   │   ├── 6) CpG特征
│   │   │   ├── 7) k-mer统计
│   │   │   ├── 8) 序列熵
│   │   │   ├── 9) 序列复杂度
│   │   │   ├── 10) Homopolymer
│   │   │   ├── 11) 重复比例
│   │   │   ├── 12) 窗口描述符
│   │   │   ├── 13) 密码子统计
│   │   │   ├── 14) 全部描述符计算：all_descriptors（固定 schema 的 240 个输出字段）
│   │   │   └── 描述符 schema 查询与带来源的二核苷酸理化表加载
│   │   ├── 序列表征与特征处理
│   │   │   ├── 1) 整数编码
│   │   │   ├── 2) One-hot编码
│   │   │   ├── 3) k-mer 特征
│   │   │   ├── 4) k-mer Sketch（MinHash/FracMinHash）
│   │   │   ├── 5) 神经网络表征
│   │   │   │   ├── 1. AlphaGenome
│   │   │   │   ├── 2. Caduceus
│   │   │   │   ├── 3. DNABERT-2
│   │   │   │   ├── 4. Enformer
│   │   │   │   ├── 5. Evo 2
│   │   │   │   ├── 6. GENERator
│   │   │   │   ├── 7. GROVER
│   │   │   │   ├── 8. HyenaDNA
│   │   │   │   ├── 9. JanusDNA
│   │   │   │   ├── 10. LucaOne（默认）
│   │   │   │   └── 11. Nucleotide Transformer v2
│   │   │   ├── 44 种 DNA 编码：encode_dna / DNA_ENCODINGS
│   │   │   │   ├── 词频与近似/非连续词：Kmer / RCKmer / Mismatch / Subsequence
│   │   │   │   ├── 碱基及位置：NAC / ANF / ENAC / binary / PS2 / PS3 / PS4
│   │   │   │   ├── 双核苷酸：CKSNAP / ASDC / DBE / LPDF
│   │   │   │   ├── 化学/EIIP/互信息：NCP / EIIP / PseEIIP / MMI
│   │   │   │   ├── Z-curve：Z_curve_9bit / Z_curve_12bit / Z_curve_36bit / Z_curve_48bit / Z_curve_144bit
│   │   │   │   ├── 理化性质位置编码：DPCP / TPCP
│   │   │   │   ├── 自/互协方差：DAC / DCC / DACC / TAC / TCC / TACC
│   │   │   │   ├── 自相关：NMBroto / Moran / Geary
│   │   │   │   ├── 伪核苷酸组成：PseDNC / PseKNC / PCPseDNC / PCPseTNC / SCPseDNC / SCPseTNC
│   │   │   │   └── 监督编码：PSTNPss / PSTNPds / KNN（显式 training_ids；理化编码需显式性质表）
│   │   │   ├── 拟合与持出数据复用：analyze_features / FeatureAnalysis.transform
│   │   │   ├── 缩放：standard / minmax / l1 / l2 / none
│   │   │   ├── 降维：pca / tsne / lda / latent_dirichlet（t-SNE 无持出 transform）
│   │   │   ├── 筛选：chi2 / f_classif / f_regression / mutual_info_classif / mutual_info_regression
│   │   │   └── 筛选：information_gain / pearson / variance；none 仅预处理
│   │   └── DNA指纹
│   │       ├── 1) Hashed k-mer 位指纹
│   │       ├── 2) Panel 存在性指纹
│   │       ├── motif / 限制酶位点：motif_fingerprint / restriction_fingerprint
│   │       ├── GC 空间分布 / 重复 / 编码区：gc_spatial_fingerprint / repeat_fingerprint / coding_fingerprint
│   │       ├── 热力学 / 多尺度 / 混合：thermodynamic_fingerprint / multiscale_fingerprint / hybrid_fingerprint
│   │       └── 特征预处理：fit_preprocessor / transform（缺失值、标准化、归一化与低方差过滤）
│   ├── 2) 序列搜索
│   │   ├── 通用搜索
│   │   │   ├── 1) Exact search
│   │   │   ├── 2) Subsequence search
│   │   │   ├── 3) Approximate matching
│   │   │   ├── 4) 反向互补搜索
│   │   │   ├── 5) 最近邻搜索
│   │   │   ├── 6) 数据库索引
│   │   │   ├── Sketch 索引构建、保存与重载：build_sketch_index / save_sketch_index / load_sketch_index
│   │   │   ├── 本地核酸 BLAST：build_blast_database / local_blast（外部 BLAST+）
│   │   │   ├── BLAST XML/tab 与蛋白 HMMER2/3 结果解析：parse_search_results（Biopython）
│   │   │   └── Profile HMM 构建/搜索与核酸表解析：build_profile / search_profile / parse_nhmmer_table（hmmbuild / nhmmer / hmmsearch）
│   │   └── 序列功能搜索
│   │       ├── 1) motif搜索
│   │       ├── 2) 起止密码子
│   │       ├── 3) 启动子motif
│   │       ├── 4) TF motif
│   │       ├── 5) 限制酶位点
│   │       ├── 6) CRISPR PAM
│   │       ├── 7) 回文序列
│   │       ├── 8) 倒置重复
│   │       ├── 9) 串联重复
│   │       ├── 10) 微卫星
│   │       ├── 开放阅读框：scan_orfs
│   │       ├── CpG 岛与低复杂度区域：find_cpg_islands / find_low_complexity_regions
│   │       ├── 原生 counts / PWM / PSSM、背景与简并共识：MotifModel / build_motif
│   │       ├── 保守 FPR 阈值与双链/圆形扫描：fpr_threshold / scan_motif_model
│   │       ├── 引物性质、引物匹配与 PCR 模拟
│   │       ├── 引物设计请求与 Primer3 CLI 设计适配器
│   │       └── CRISPR 候选与脱靶扫描：scan_crispr_candidates
│   ├── 3) 理化性质
│   │   ├── 1. 分子量
│   │   ├── 2. 260 nm 消光系数
│   │   ├── 3. 熔解温度Tm
│   │   ├── 4. 盐浓度修正
│   │   ├── 5. 局部熔解特征
│   │   ├── 6. 合成风险
│   │   └── 热力学性质
│   │       ├── 1) 热力学参数
│   │       ├── 2) Nearest-neighbor
│   │       ├── 3) Duplex stability
│   │       ├── 4) 碱基堆积
│   │       ├── 5) Hairpin
│   │       ├── 6) Self-dimer
│   │       └── 7) Heterodimer
│   ├── 4) 双链热力学扩展
│   │   ├── 1. 完全互补双链的 ΔH、ΔS、ΔG 和 Tm
│   │   ├── 2. 统一双链稳定性结果
│   │   ├── 3. 相邻碱基对步骤贡献
│   │   ├── 4. 条件与 Na⁺/K⁺ 单价盐
│   │   ├── 5. Ka、Kd 和双链比例
│   │   ├── 6. 理论熔解曲线
│   │   ├── 7. 5′/3′ 末端稳定性
│   │   ├── 8. DMSO 和甲酰胺经验修正
│   │   └── 9. Primer3 CLI 的 Mg²⁺、dNTP、mismatch 和 dangling end
│   ├── 5) 二级结构性质
│   │   ├── 1. Dot-bracket 结构解析
│   │   ├── 2. 配对概率与窗口可接近性
│   │   ├── 3. 目标结构的集合缺陷
│   │   ├── 4. 目标结构热力学概率
│   │   ├── 5. NUPACK 被动可用性检查
│   │   ├── 6. NUPACK 单复合物集合分析
│   │   └── 7. NUPACK tube 多复合物平衡
│   ├── 6) 三维结构与力学性质
│   │   ├── 1. 读取单个 PDB 模型
│   │   ├── 2. 读取 PDB 多模型集合
│   │   ├── 3. 显式坐标几何分析
│   │   ├── 4. NMR/多模型 RMSF 柔性
│   │   ├── 5. 3DNA bp_step.par 标准参数解析
│   │   └── 6. DSSR JSON 摘要解析
│   ├── 7) 换算
│   │   ├── 1. 260 nm 单链摩尔消光系数
│   │   ├── 2. 单链/双链理论光学性质
│   │   ├── 3. 1 OD260 对应的 nmol 和质量
│   │   ├── 4. A260 到摩尔浓度、质量浓度和总量
│   │   ├── 5. 浓度、物质的量和质量互换
│   │   └── 6. 染料和修饰基团的显式校正
│   ├── 8) 性质预测：深度学习与机器学习
│   │   ├── 1. 自定义模型训练与预测
│   │   │   ├── 用户自定义分类与回归
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
│   │   ├── 2. 表达、调控与接触图预测
│   │   │   ├── AlphaGenome 基因组功能轨道：RNA-seq、CAGE、PRO-cap、ATAC-seq、DNase-seq
│   │   │   ├── 组蛋白修饰、转录因子结合与三维 DNA 接触图
│   │   │   └── Enformer 人类/小鼠调控轨道
│   │   ├── 3. 剪接与基因组区域预测
│   │   │   ├── 剪接位点、剪接使用比例与剪接连接计数
│   │   │   ├── 14 类单碱基基因组分割
│   │   │   │   ├── 1. 蛋白质编码基因（protein-coding gene）
│   │   │   │   ├── 2. 长链非编码 RNA（lncRNA）
│   │   │   │   ├── 3. 外显子（exon）
│   │   │   │   ├── 4. 内含子（intron）
│   │   │   │   ├── 5. 剪接供体位点（splice donor）
│   │   │   │   ├── 6. 剪接受体位点（splice acceptor）
│   │   │   │   ├── 7. 5′ 非翻译区（5′ UTR）
│   │   │   │   ├── 8. 3′ 非翻译区（3′ UTR）
│   │   │   │   ├── 9. CTCF 结合区域（CTCF-bound）
│   │   │   │   ├── 10. 多聚腺苷酸化信号（polyA signal）
│   │   │   │   ├── 11. 组织特异性增强子（tissue-specific enhancer）
│   │   │   │   ├── 12. 组织不变性增强子（tissue-invariant enhancer）
│   │   │   │   ├── 13. 组织特异性启动子（tissue-specific promoter）
│   │   │   │   └── 14. 组织不变性启动子（tissue-invariant promoter）
│   │   │   ├── 长上下文外显子概率
│   │   │   ├── NT Revised 完整微调分类任务：剪接受体、剪接位点类型与剪接供体
│   │   │   └── Genomic Benchmarks 完整微调分类任务：编码区/基因间区
│   │   ├── 4. 表观遗传与开放染色质判别
│   │   │   ├── NT Revised 完整微调分类任务：H2A.Z、H3K27ac、H3K27me3、H3K36me3、H3K4me1/2/3、H3K9ac/me3、H4K20me1
│   │   │   └── Genomic Benchmarks 完整微调分类任务：人类开放染色质区域
│   │   ├── 5. 启动子与增强子判别
│   │   │   ├── NT Revised 完整微调分类任务：增强子判别与类型分类；启动子、非 TATA / TATA 启动子
│   │   │   ├── Genomic Benchmarks 完整微调分类任务：果蝇增强子、小鼠增强子
│   │   │   └── 人类 Cohn / Ensembl 增强子、人类调控元件类型与人类非 TATA 启动子
│   │   ├── 6. 变异效应与蛋白质稳定性
│   │   │   ├── 基因组功能轨道及人/鼠调控轨道的 SNV 轨道差值
│   │   │   ├── Evo 2 长上下文零样本变异评分
│   │   │   ├── 等位基因条件概率变异评分
│   │   │   ├── 蛋白稳定性与 ESM 蛋白质变异效应
│   │   │   └── 剪接变异影响：SpliceAI / Pangolin / 最大 delta 评分
│   │   ├── 7. 物种与分子功能分类
│   │   │   ├── 已训练核酸/蛋白下游任务：SupKTax、GenusTax、物种分类
│   │   │   ├── 原核蛋白亚细胞定位、ncRNA 家族
│   │   │   └── Genomic Benchmarks 完整微调分类任务：人类/线虫
│   │   ├── 8. 分子关系与互作预测
│   │   │   ├── 已训练核酸/蛋白下游任务：核酸–蛋白中心法则关系
│   │   │   └── 甲型流感抗原关系、蛋白–蛋白互作、ncRNA–蛋白互作
│   │   └── 9. 蛋白质特征与结构预测
│   │       ├── ESM 蛋白质特征与序列评分
│   │       └── NVIDIA NIM：AlphaFold2 / ESMFold 蛋白质结构预测
│   ├── 9) 基因名称与区间注释
│   │   ├── 基因标识、转录本、同源关系与数据库注释查询
│   │   ├── 基因组 DNA + 物种 → annotate_genes / GeneAnnotationConfig
│   │   ├── 默认 RefSeq 组装解析、显式组装或本地 FASTA + GFF3
│   │   ├── minimap2 参考比对与已有基因注释映射（非从头基因预测）
│   │   ├── 基因名称、ID、方向、重叠基因与 gene/intergenic/unknown 分段
│   │   ├── Python 0-based 半开区间；JSON/CLI 1-based 闭区间
│   │   └── 组装与注释固定缓存、文件校验、离线复用及来源记录
│   ├── 10) 系统发育与分子进化
│   │   ├── 原生树对象与建树：TreeNode / build_tree（NJ / 按簇大小加权 UPGMA）
│   │   ├── 原生树比较、共识与剪枝：robinson_foulds / consensus_tree / prune_tree
│   │   ├── 原生群落系统发育距离：unifrac（weighted / unweighted，要求定根树）
│   │   └── Nei–Gojobori dN/dS：比较分析与数据库命名空间的选择压力接口
│   ├── 11) 统计、差异与富集分析
│   │   ├── 1. 统计与差异分析
│   │   │   ├── 距离矩阵校验与置换检验：distance_matrix / distance_test（PERMANOVA / ANOSIM / centroid PERMDISP）
│   │   │   ├── 矩阵相关与环境变量子集搜索：mantel / bioenv
│   │   │   ├── 主坐标与约束排序：pcoa / ordination（CA / CCA / RDA）
│   │   │   ├── 组成变换及逆变换：composition_transform（CLR / ILR / ALR）
│   │   │   ├── 原生组成差异分析：ancom
│   │   │   ├── alpha 多样性：alpha_diversity（observed / Shannon / Simpson / inverse Simpson）
│   │   │   ├── beta 距离：beta_diversity（Bray-Curtis / Jaccard / Euclidean / cityblock / cosine）
│   │   │   ├── 原始 ANCOM-BC：ancom_bc（OLS / sandwich 方差 / EM 采样偏差校正）
│   │   │   ├── 差异表达：differential_expression（未收缩 NB2 GLM / size factors / 显式 contrast / Wald）
│   │   │   └── 多重检验：adjust_pvalues（Benjamini–Hochberg）
│   │   └── 2. 基因集与通路富集
│   │       └── Enrichr / PANTHER / Reactome / STRING 富集分析
│   ├── 12) 测序数据分析
│   │   ├── BAM 汇总、覆盖度与参考坐标共识：alignment_summary / coverage_consensus（pysam；共识不写入 indel）
│   │   ├── DNA-seq：dna_seq_workflow（质控→比对→BAM→变异）
│   │   ├── 单样本 RNA-seq：rna_seq_workflow（HISAT2 / StringTie / featureCounts）
│   │   ├── 已处理 BAM 的 ChIP-seq：chip_seq_workflow（MACS3 peaks + samtools coverage）
│   │   └── 双向局部比对与参考坐标共识：sanger_consensus
│   ├── 13) 变异检测、注释与评分
│   │   ├── 1. 变异检测与统计
│   │   │   ├── 变异调用、参考规范化与 VCF 索引：call_variants（bcftools）
│   │   │   ├── 按过滤条件统计 VCF 变异数量：vcf_count_variants
│   │   │   ├── 规范化 VCF 表示与多等位位点：vcf_normalize
│   │   │   ├── 汇总变异类型与样本基因型：vcf_summary_stats
│   │   │   └── 编码区变异比例：coding_variant_fraction
│   │   ├── 2. 变异注释与标识转换
│   │   │   └── Ensembl VEP / MyVariant / FAVOR / VariantValidator 注释与标识转换
│   │   └── 3. 变异效应评分与实验数据
│   │       └── AlphaMissense / CADD / EVE / MaveDB 评分与实验效应数据
│   ├── 14) 遗传模拟
│   │   └── 模拟独立位点的等位基因遗传：simulate_inheritance
│   └── 15) 基因组与 DNA 工程
│       ├── 基因组关系与育种预测：genomic_relationship_matrix / genomic_prediction
│       ├── 原核 CDS 预测：predict_prokaryotic_genes（Pyrodigal）
│       ├── 甲基化汇总、DNA 存储与纳米结构设计
│       ├── 基因型 PCA：genotype_pca（NumPy；剂量矩阵输入）
│       └── 基因组 ANI：average_nucleotide_identity（PyFastANI；无需 FastANI 命令）
├── 4. 数据评价
│   ├── 1) 常用评价指标
│   │   ├── 1. Validity
│   │   │   ├── 1.1 语法合法性：evaluate_validity
│   │   │   └── 1.2 参考基因组一致性：evaluate_reference_consistency
│   │   ├── 2. Uniqueness
│   │   ├── 3. Diversity（平均两两原始 Levenshtein 距离）
│   │   ├── 4. Novelty（相似度 / 平均最近参考 Levenshtein）
│   │   ├── 5. Fréchet DNA distance
│   │   ├── 6. Frag
│   │   ├── 7. SNN
│   │   ├── 8. Ambiguity
│   │   ├── 9. Redundancy
│   │   ├── 复杂度与质量评价：evaluate_complexity / evaluate_quality
│   │   ├── 参考库构建、最近参考查询与记忆/复制评价
│   │   ├── 综合评分卡：evaluate_scorecard
│   │   ├── 原始 reads 质控：quality_control（FastQC）
│   │   ├── 自定义标量性质矩阵：sequence_properties
│   │   ├── 高斯 MMD² / 多项式 KID：kernel_distance
│   │   ├── Precision / Recall / Density / Coverage：manifold_metrics
│   │   ├── 校准后的 clipped density / coverage：clipped_manifold_metrics
│   │   ├── 相对生成器训练集的真实性：authenticity
│   │   ├── Pareto 层级与精确超体积：pareto_ranks / hypervolume
│   │   ├── 多性质区间命中率及直方图 KL：property_hit_rate / property_kl
│   │   └── KDE 分布一致性及标准误：conformity_score
│   ├── 2) 相似度计算
│   │   ├── 1. k-mer相似度
│   │   ├── 2. 指纹相似度
│   │   ├── 3. Sketch相似度
│   │   ├── 4. Dashing sketch相似度
│   │   ├── 5. Reference similarity
│   │   ├── 描述符分布相似度：evaluate_distribution_similarity
│   │   └── 向量与矩阵：kmer_vector_similarity / similarity_matrix / pairwise_matrix
│   └── 3) 序列距离与比对
│       ├── 1. 序列距离
│       │   ├── 1) Hamming distance
│       │   └── 2) Edit distance
│       ├── 2. 成对序列比对
│       │   ├── 2) 1 全局比对
│       │   ├── 2) 2 局部比对
│       │   └── 2) 3 半全局比对
│       ├── 原生 MSA 对象、等宽/ID 校验与 p-distance：MultipleAlignment.from_sequences / distances
│       ├── 多序列比对：align_sequences（外部 MAFFT / MUSCLE 5 / Clustal Omega）
│       └── 单端/双端比对及私有索引：map_reads（BWA / Bowtie2）
├── 5. 可视化
│   ├── 1) 序列文字图
│   ├── 2) 位置高亮
│   ├── 3) 样式控制
│   ├── 4) Gap 显示
│   ├── 5) 线性序列图
│   ├── 6) 环状 DNA 图
│   ├── 7) 自定义字符表示
│   ├── 8) Alignment 图
│   ├── 相似度矩阵热图：plot_similarity_matrix
│   ├── 自包含、可筛选与展开的 HTML 报告：build_html_report / save_html_report
│   ├── 图片导出：save_svg / save_image（PNG、SVG、JPG，兼容 TIFF/PDF）
│   ├── 树图：plot_tree（Matplotlib）
│   ├── 正反链精确 k-mer dotplot：dotplot（Matplotlib）
│   ├── Sanger 色谱图：plot_chromatogram
│   └── 定性凝胶模拟：simulate_gel（Matplotlib）
├── 6. 大语言模型工具
│   ├── Python API：通用大模型调用（dnakit.llm）
│   │   ├── 配置与结果：LLMConfig / LLMResult
│   │   ├── 同步/异步：call_llm / acall_llm；自定义提示词与当前步骤数据
│   │   ├── 文本/JSON 输入，文本/JSON 对象输出与 token 用量
│   │   ├── 本地文本文件：file_path；本地预览 read_llm_file
│   │   ├── 编码、文件与请求大小限制；显式指定文件后发送完整文本
│   │   ├── 服务商预设：DeepSeek / Gemini / Qwen / Kimi / OpenRouter / OpenAI / custom
│   │   ├── 预设查询：available_llm_providers / get_llm_provider
│   │   └── 自定义兼容地址、模型、密钥配置及 extra_body 扩展参数
│   ├── MCP：Agent 与 MCP
│   │   ├── 本地 stdio 服务：dnakit-mcp；可选依赖 dnakit[agent]
│   │   ├── 支持 MCP 的客户端接入（Codex、Claude、Cursor、VS Code 等）
│   │   ├── 工具统计：dnakit_catalog_stats
│   │   ├── 功能分类：list_dnakit_categories
│   │   ├── 分页目录：list_dnakit_tools
│   │   ├── 按功能搜索：search_dnakit_tools
│   │   ├── 参数 Schema 与执行属性：describe_dnakit_tool
│   │   ├── JSON 参数调用：call_dnakit_tool
│   │   ├── Python 工具注册表与服务构建：default_tool_registry / create_server
│   │   └── 有界 JSON 返回、来源追踪、Python-only 标记与 allow_side_effects 授权
│   ├── Skills：安装 DNAKit 技能与查找分析工作流
│   └── Agent：内置自然语言任务执行（run_agent）
│       ├── 自动检索工具、查看参数、调用与结果反馈
│       ├── 可选 Python 补充执行、失败恢复与运行日志
│       └── 关键词／语义／混合检索、执行权限与步数限制
└── FAQ
    └── 工程与扩展
        ├── 随包脚本执行与资源定位、ENCODE 凭据、Notebook 结果与 PDF 报告
        ├── 1. 统一后端接口
        ├── 2. 原生与外部标记
        ├── 3. Python API
        ├── 4. CLI
        │   ├── 环境与后端：info / backends
        │   ├── 数据准备：download-genome / convert / normalize / validate
        │   ├── 分析与比较：describe / fingerprint / search / orfs / compare
        │   ├── 数据集与报告：deduplicate / split / report / workflow
        │   ├── 自定义性质：property-models / inspect-property-data / train-property / predict-property
        │   └── 基因注释：annotate-genes
        ├── 5. 配置工作流：load_workflow / run_workflow；产物 manifest 读写、校验与复用
        ├── 6. 批量计算：iter_batch / run_batch；进度、取消与逐条错误
        ├── 7. 并行计算：有界线程任务与稳定输出顺序
        ├── 8. 分块与流式处理
        ├── 9. 缓存
        ├── 10. 随机种子
        ├── 11. 版本追踪
        ├── 12. 错误与警告
        ├── 14. 一致性验证
        ├── 16. 文档与教程
        ├── 17. Galaxy 包装器：convert / deduplicate / split / report（需 Galaxy 环境）
        ├── 开发验证：单元测试与一致性检查
        └── DAG 定义、依赖/环检查与失败传播：Step / run_workflow
```

目录核对范围与数量见 [FAQ](../../faq/tutorial_overview.md)。
