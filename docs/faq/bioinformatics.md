# 生物信息学补充说明

各函数的用法与示例已归入对应任务教程；本页补充依赖、参数清单、实现约定和验证方法。

`dnakit.bioinformatics` 补充传统生物信息学、DNA 特征工程、生成评价、统计分析和分子设计。以下接口属于开发版，**仅提供 Python API**。本页按真实实现范围列出能力，不表示与 Biopython、SeqKit、scikit-bio、UGENE、iLearnPlus、iFeatureOmega、seqme、pydna、DNA Chisel 的全部功能等价。

## 需要安装哪些依赖？ {#dependencies}

核心领域逻辑由 DNAKit 实现；NumPy/SciPy 提供数组、线性代数、距离和数值求解，scikit-learn 提供通用机器学习算法。没有调用 iFeatureOmega、iLearnPlus、seqme 或 DNA Chisel 的功能接口。

| 可选依赖 | 用途 | 安装包 |
| --- | --- | --- |
| `bio-analysis` | 原生编码、评价、统计；通用特征选择/降维/聚类 | NumPy、SciPy、scikit-learn |
| `bio-formats` | MSA、树、motif、BLAST/HMMER、ABI、SnapGene/GenBank 文件编解码 | Biopython |
| `bio-reads` | FASTA 索引与 BAM 随机访问 | pysam/htslib |
| `bio-plots` | 树、dotplot、色谱、凝胶绘图 | Matplotlib |
| `bio-assembly` | 复杂同源重叠组装、宽容构建文件解析 | pydna |
| `bio` | 上述依赖的合集 | 不安装外部可执行程序 |
| `bio-validation` | 测试中的独立数值对照；运行 API 不需要 | scikit-bio、statsmodels |

安装当前工作区的完整 Python 依赖：

```bash
python -m pip install -e ".[bio]"
```

只安装数组/统计/特征分析依赖：

```bash
python -m pip install -e ".[bio-analysis]"
```

按模块导入，缺少可选依赖时仅对应操作报错；不自动安装、下载数据或提交远程任务。格式读取返回的 Biopython 对象是明确的编解码边界；原生计算通常返回 dataclass、字典和 NumPy 数组。

## 比对、搜索与进化树有哪些实现约定？ {#alignment}

| 模块 | Python API | 行为 |
| --- | --- | --- |
| `alignment` | `MultipleAlignment.from_sequences`, `.distances` | 原生 MSA 对象；等宽、唯一 ID，p-distance 的 pairwise deletion / gap-as-state |
| `alignment` | `read_alignment`, `write_alignment` | 单个 FASTA/Clustal/Stockholm/PHYLIP 等 MSA 的 Biopython 编解码 |
| `alignment` | `align_sequences` | MAFFT `--auto`、MUSCLE **5**、Clustal Omega；还原原始 ID 并验证序列内容 |
| `alignment` | `build_blast_database`, `local_blast` | 本地核酸数据库与 BLAST XML 输出 |
| `alignment` | `parse_search_results` | 流式 Biopython SearchIO QueryResult；BLAST XML/tab、蛋白 HMMER2/3 格式 |
| `alignment` | `build_profile`, `search_profile`, `parse_nhmmer_table` | hmmbuild、核酸 nhmmer / 蛋白 hmmsearch；原生 nhmmer 表解析 |
| `phylogeny` | `TreeNode`, `build_tree` | 原生 NJ / 按簇大小加权的 UPGMA |
| `phylogeny` | `read_tree`, `write_tree`, `plot_tree` | Newick/Nexus/phyloXML/NeXML 编解码、树图 |
| `phylogeny` | `robinson_foulds`, `consensus_tree`, `prune_tree` | 原生 rooted/unrooted RF、严格多数/严格共识、剪枝 |

远程 BLAST 已有 `dnakit.search.blast` 的显式提交/轮询 API，可与本地结果解析配合使用，不重复实现网络服务。

NJ 保留负枝长估计，根位置任意；UPGMA 假定分子钟。RF 要求完全相同的叶节点集合，归一化分母为两树分裂总数。共识树要求相同的根解释，阈值为 `0.5` 时保留严格多数分支，`1` 时保留所有树共有分支。未提供 GME/BME、树重排搜索或最大似然建树。

## FASTX、坐标和外部测序流程有哪些约定？ {#reads}

| 模块 | Python API | 输入/输出与约定 |
| --- | --- | --- |
| `reads` | `Read`, `read_fastx`, `write_fastx` | FASTA、multiline FASTQ、gzip；质量为整数 Phred；写入返回记录数 |
| `reads` | `pair_reads` | 两个 read 迭代器，按去除 `/1`、`/2` 的 ID 配对；返回 ReadPair，孤儿一侧为 None |
| `reads` | `repair_fastq` | 只恢复合法四行记录，返回输出路径/恢复数/丢弃行数 |
| `reads` | `arrange_reads`, `common_reads`, `rename_reads`, `sequence_digest` | 排序/固定种子打乱、集合交集、正则重命名、有记录边界的 SHA256 |
| `reads` | `Region`, `read_regions` | BED 或指定 GTF feature；GTF 1-based 闭区间转换一次 |
| `reads` | `index_reference`, `IndexedReference.fetch` | 显式建立 `.fai`、上下文管理随机读取，支持反向互补 |
| `reads` | `alignment_summary`, `coverage_consensus` | pysam 解码、原生 flag 统计及区域 A/C/G/T 覆盖/多数共识 |
| `ngs` | `map_reads` | BWA MEM / Bowtie2 私有索引 → SAM |
| `ngs` | `prepare_bam` | samtools name-sort/fixmate/coordinate-sort/markdup/index/stats |
| `ngs` | `call_variants` | bcftools mpileup/call/norm → 索引 VCF；显式倍性 |
| `ngs` | `assemble_reads` | SPAdes 短读长/双端/含 Nanopore 的混合组装 |
| `ngs` | `quality_control`, `dna_seq_workflow` | FastQC；QC→BWA→BAM→变异调用，各阶段结果保留 |
| `ngs` | `rna_seq_workflow` | HISAT2→BAM→StringTie→featureCounts；链特异性 0/1/2 |
| `ngs` | `chip_seq_workflow` | 已准备的 treatment/control BAM → MACS3 peaks 与 coverage |

FASTA→FASTQ 必须显式传 `default_quality`，不推测质量。repair 不保证损坏文件的歧义重同步，multiline 文件应使用严格读取。配对的待匹配 payload 受 `max_pending` 限制，重复 ID 检查集合仍随记录总数增长。

Python 区间为 **0-based 半开区间**。BAM 共识保持参考正链坐标，`Region.strand` 不改变 BAM 共识方向；只计对齐 A/C/G/T，过滤 secondary/supplementary/duplicate/QC-failed 及低质量 read，平票或不足深度输出 N。插入/删除不写入共识，使用 VCF 获取这些事件。

外部程序须自行安装到 PATH。`ExecutionConfig(threads=4, timeout_seconds=3600)` 控制执行；不使用 shell。输出目录必须尚不存在，保留中间文件与诊断信息。`max_output_bytes` 限制进程日志，`max_artifact_bytes` 限制 stdout 或显式监控的输出文件，并非整个第三方工作目录的磁盘配额。原始输入不覆盖；索引 API 是明确例外，会在指定参考旁建立索引。

参考 FASTA 必须已有 `.fai`，变异调用需要已排序、索引 BAM。RNA 流程处理一个样本；跨样本生物学重复的差异分析见下文。ChIP 流程从已处理 BAM 开始，可先调用 mapping/BAM API。没有 assembly browser GUI。

## 支持哪些编码和特征分析方法？ {#encoding}

`encoding.encode_dna({id: sequence}, descriptor, parameters=..., property_tables=..., labels=..., training_ids=...)` 返回 `EncodingResult(ids, names, values, descriptor, parameters, schema_version)`。`values` 为只读 samples×features 数组；当前模式为 `dnakit.dna-encoding.v1`，接受 A/C/G/T，未知参数拒绝。

| 编码 | 参数/定义 |
| --- | --- |
| Kmer / RCKmer | `k=3, normalized=True`；重叠词频 / 按字典序合并反向互补词 |
| Mismatch / Subsequence | `k=3, normalized=True`；分别 `mismatches=1` 的 Hamming 邻域累计、`delta=0.5` 的非连续子序列跨度权重 |
| NAC / ANF / ENAC | 单核苷酸频率 / 当前碱基的前缀频率 / `window=5` 滑窗频率 |
| binary / PS2 / PS3 / PS4 | 按 A,C,G,T 顺序的逐位置 k-mer one-hot |
| CKSNAP / ASDC / LPDF / DBE | `gap=3, normalized=True` 间隔双核苷酸 / 所有有序位置对频率 / 局部双核苷酸前缀频率 / 双核苷酸 4-bit |
| NCP / EIIP / PseEIIP | 三组化学分类 / EIIP 位置值 / 三核苷酸频率×EIIP 和 |
| Z_curve_9bit/12bit/36bit/48bit/144bit | 相位与前缀条件下的 x/y/z 对比，分母为总重叠窗口数 |
| DPCP / TPCP | 显式理化性质表标准化后逐位置输出 |
| DAC / DCC / DACC、TAC / TCC / TACC | `lag=3`；居中 profile 的 auto/cross covariance；纯 cross 需要至少两张表 |
| NMBroto / Moran / Geary | `lag=3`；标准化性质乘积 / 方差归一自相关 / 平方差统计量 |
| PseDNC / PseKNC、PCPseDNC/TNC、SCPseDNC/TNC | `lag=3, weight=0.05`；PseKNC 另有 `k=3`；成分与滞后项共享分母 |
| MMI | 位置联合分布的互信息贡献按碱基多重集合聚合：10 项二元 MI + 20 项三元 co-information |
| PSTNPss / PSTNPds | 显式二分类训练集的逐位置三核苷酸频率差；ds 折叠 T→A、G→C |
| KNN | `neighbors=3`；按 Hamming 距离选训练邻居，输出正类比例，距离平票按 ID |

全部名称见 `encoding.DNA_ENCODINGS`。位置编码要求等长。理化性质表为 `{property_name: {word: finite_value}}`，必须覆盖全部 16 个二核苷酸或 64 个三核苷酸，并在完整词表上做总体标准差标准化；不捆绑未注明来源的理化表。`NMBroto`、Moran、Geary 使用二核苷酸表。序列过短、常数 profile、无定义分母报错。

Pseudo type-I 的滞后项为跨性质的平方差均值；SC 使用各性质的滞后乘积均值。Subsequence 的每项权重为 `delta ** span_length`，normalized 时除以总权重。MMI 使用一致的同长度窗口边缘分布，二元/三元特征分别求和得到总 MI/co-information；这是明确的 DNAKit 定义，**不承诺与 iFeatureOmega 的同名向量逐值相等**。PSTNP/KNN 只读取 `training_ids` 的标签，并逐样本排除自身，测试标签不参与计算。

`features.analyze_features` 支持 standard/minmax/L1/L2/none、PCA、t-SNE、监督 LDA、Latent Dirichlet Allocation、chi2、F-classification/regression、MI、Information Gain（训练集分位数分箱，`parameters={"bins": 10}`）、Pearson、variance、none。监督方法需要 `labels`；chi2/主题模型需要非负数据；t-SNE 无持出 transform。`cluster_features` 支持 kmeans、minibatch_kmeans、gmm、spectral、affinity_propagation、mean_shift、dbscan、agglomerative、birch，返回已拟合 sklearn 对象；`mcl` 是原生随机矩阵 expansion/inflation，输入对称非负 affinity，返回索引簇，未收敛报错。聚类不自动缩放。

## 生成评价指标如何解释？ {#metrics}

| `metrics` API | 定义与输入 |
| --- | --- |
| `kernel_distance` | embeddings；高斯 MMD² 或三次多项式 KID，默认无偏估计，不截断负估计值 |
| `manifold_metrics` | PRDC precision/recall/density/coverage，kNN 严格开球、自身排除 |
| `clipped_manifold_metrics` | Salvy 等的 median-clipped density、leave-one-out 密度校准及 beta-binomial coverage 校准；闭球、两组样本数相同 |
| `authenticity` | 生成点最近训练点距离 > 该训练点最近其他训练点距离的比例；返回 0–1 |
| `pareto_ranks` | 每列显式 `maximize`，非支配层从 0 开始，重复点同层 |
| `hypervolume` | 每列显式方向、较差 `reference_point`，原生精确递归；最多 5 维/1000 点并限制递归工作量 |
| `property_hit_rate` | 每列闭区间 `lower/upper`，返回联合命中率、逐列命中率和布尔 mask |
| `property_kl` | KL(generated∥reference)，共同直方图边界、显式 pseudocount、自然对数 |
| `conformity_score` | 固定带宽各向同性 Gaussian KDE，参考集 cross-fit，再按验证密度经验分布评分 |
| `sequence_properties` | 接受命名的 Python 标量 predictor，返回稳定 IDs/性质名/矩阵，供上述 API 使用 |

`authenticity` 必须使用**生成器训练集**；其他 reference 通常是独立参考/验证集。几何指标的 embedding 模型、缩放及 k 必须在比较前固定。PRDC 与 clipped 指标不能混用名称；KID 不是 FID。性质 KL 的直方图估计与 seqme KDE 估计不等价。Conformity 返回 `value`、`standard_error` 和 fold 分数；这是分布一致性指标，不是类别后验概率。

## 统计与计数推断有哪些适用条件？ {#statistics}

| 模块 | API | 范围 |
| --- | --- | --- |
| `statistics` | `distance_matrix`, `distance_test` | 距离矩阵校验；PERMANOVA、ANOSIM、centroid PERMDISP，固定 seed 置换 |
| `statistics` | `mantel`, `bioenv` | Pearson/Spearman Mantel；小规模变量子集 BIOENV |
| `statistics` | `pcoa`, `ordination` | PCoA、CA、CCA、RDA；PCoA 报告负特征值，CCA/RDA 仅约束轴 |
| `statistics` | `composition_transform` | CLR/ILR/ALR 及逆变换；ILR 使用正交 Helmert 基 |
| `statistics` | `ancom` | pairwise log-ratio、每 taxon Holm 校正、W；不擅自给 W 设显著性阈值 |
| `statistics` | `alpha_diversity`, `beta_diversity` | observed/Shannon/Simpson/inverse Simpson；Bray-Curtis/Jaccard/Euclidean/cityblock/cosine |
| `phylogeny` | `unifrac` | 原生 rooted weighted/unweighted UniFrac；taxa 列须与树叶完全一致 |
| `counts` | `ancom_bc` | 原始 ANCOM-BC：OLS、稳健 sandwich 方差、三高斯分量 EM 采样偏差校正 |
| `counts` | `differential_expression` | 原生未收缩 NB2 GLM、size factors、显式 contrast、Wald/BH |
| `counts` | `adjust_pvalues` | 一维 Benjamini–Hochberg 校正 |

置换 p 值采用 `(extreme+1)/(permutations+1)`，`permutations=0` 返回 `None`；不支持分层/重复测量置换。PERMDISP 当前要求欧氏距离。组成变换及 ANCOM-BC 的零值必须由调用者显式处理；ANCOM 可显式给 pseudocount。UniFrac 树必须按生物学意义定根、枝长非负。

ANCOM-BC 接受 samples×taxa、包含截距的满秩数值 design，至少 10 taxa，返回 feature×covariate 的**自然对数**效应、SE、Wald、p/q、bias、EM 迭代数；BH 按协变量分别校正。不是 ANCOM-BC2，不包括 structural zeros、global test 或 bias-variance adjustment。不收敛或方差退化直接报错。

NB2 接受 samples×genes 整数 counts、数值 design、显式 contrast。至少两个剩余自由度，需真实生物学重复。默认以全样本正计数基因计算 median-ratio size factors，默认离散度为设计内残差矩估计；也可显式提供 `size_factors` 和 `dispersions`。输出 log2 fold change、log2 SE、Wald、p/q、收敛标记和离散度；全零基因及拟合失败的推断为 NaN。**没有 DESeq2/edgeR 的经验贝叶斯离散度趋势、效应收缩和离群值策略**。

## 优化、motif 和克隆模型有哪些规则？ {#optimization}

`optimization.specification(kind, region=None, weight=1, **parameters)` 构造原生规则；`Specification(name, evaluator)` 允许自定义 `(sequence, circular) -> SpecEvaluation`。`SpecEvaluation` 分开表示非负 `violation`、待最大化 `utility`、失败位置与说明。

| kind | 参数 |
| --- | --- |
| `gc` | `minimum=0.4, maximum=0.6, window=整段长度` |
| `avoid_pattern`, `pattern_count` | `pattern`（IUPAC）、`both_strands=False`；后者 `count=1` |
| `sequence`, `translation`, `length` | `pattern` / `protein` / `minimum=1, maximum=∞` |
| `avoid_changes` | 原始完整 `sequence`、可选绝对 `positions`；不提供 positions 时锁定 region 或全部位置 |
| `unique_kmers`, `avoid_stop`, `avoid_codons` | `k=12, both_strands=True` / `allow_terminal=True` / `codons` |
| `avoid_rare_codons` | 全 64 项 `usage`、相对最高同义密码子频率 `threshold=0.1` |
| `hairpin` | `match_length=12, window=200, minimum_loop=3`，精确反向互补茎 |
| `heterodimer`, `avoid_matches` | `references, match_length=12`；后者 `both_strands=True` |
| `regions_compatibility` | 至少两段线性绝对 `regions`、`match_length=12`；避免区间间互补片段 |
| `tm` | `minimum, maximum`、可选原有热力学 `conditions` |
| `cai`, `codon_usage`, `harmonize` | 全 64 项 `usage`；harmonize 另需逐密码子 `relative_adaptiveness` |

结构规则是明确长度的精确互补模型，不冒充完整热力学能量模型。`blast_match_specification(database, minimum_identity=0.8, minimum_length=20)` 在专用 BLAST 后端检测同源 HSP，并按序列缓存；搜索为空不证明绝对不存在同源性。`avoid_matches` 是原生精确 k-mer，未提供 Bowtie 基因组约束包装。

硬约束先于加权 objectives；默认单碱基贪心坐标搜索，`preserve_translation=True` 使用完整同义密码子替换。固定长度；`length` 只检查，不插入/删除。圆形 DNA 的窗口/模式扫描覆盖原点，跨原点位置允许 `end > length`。`search_exhausted` 不是不可满足性证明；`constraints_satisfied` 也不是全局最优证明。

GenBank `dnakit_constraint` / `dnakit_objective` qualifier 可写 JSON，例如 `{"kind":"gc","minimum":0.4,"maximum":0.6}`，经 `specifications_from_record` 读取。只自动处理精确、正向简单 feature；反向、模糊、compound 位置要求显式 Python 声明。

| 模块 | API | 实现与范围 |
| --- | --- | --- |
| `motifs` | `MotifModel`, `build_motif`, `fpr_threshold`, `scan_motif_model` | 原生 counts/PWM/PSSM、背景、简并共识、量化保守 FPR、双链/圆形扫描 |
| `motifs` | `read_motifs`, `dotplot` | JASPAR/MEME/TRANSFAC 解码；原生正反链精确 k-mer dotplot |
| `cloning` | `Duplex.from_sequence`, `.cut`, `.ligate`, `.reverse_complement`, `.ends` | 原生显式两条 5′→3′ 链、offset、黏端/平端和拓扑 |
| `cloning` | `CloningPlan.add`, `.apply`, `.steps` | 可执行步骤和全部中间产物 |
| `cloning` | `assemble_duplexes` | pydna 复杂 Gibson/同源重叠组装 |
| `cloning` | `read_sanger`, `plot_chromatogram`, `sanger_consensus` | ABI 编解码、原生 Mott trim、色谱图、原生双向局部比对及参考坐标共识 |
| `cloning` | `read_constructs` | 严格 SnapGene/GenBank 编解码；显式 permissive 使用 pydna |
| `cloning` | `simulate_gel` | 返回 Matplotlib Axes；按 log 长度迁移的定性凝胶图 |

Duplex 的 `offset` 是 Crick 3′ 端相对 Watson 5′ 端的位置，配对区必须完全互补；不支持 nick、内部 gap 或错配。凝胶未模拟染料/缓冲液动力学，不应作为定量迁移预测。Sanger 共识只统计质量达标且 identity/覆盖度达标的比对，N 表示平票/不足覆盖，不插入 indel；ABI 信号采样坐标与碱基坐标分开。

## 如何运行工作流示例与验证？ {#workflows}

`workflows.Step(name, function, dependencies={parameter: step_name}, parameters=...)` 和 `run_workflow` 是原生 DAG，检查重复名称、依赖和环；失败时可停止，或收集失败并阻断下游。`query_regions(left, right, operation='intersection')` 支持 union/intersection/difference，按染色体和链组合。均不包含图形设计器。

运行离线 Python API 示例（只需 `bio-analysis`）：

```bash
python others/examples/bioinformatics_api.py
```

运行新增测试（真实外部工具测试在对应程序不在 PATH 时跳过）：

```bash
python -m pytest others/tests/unit/test_bioinformatics_native.py others/tests/unit/test_bioinformatics_io_workflows.py others/tests/unit/test_bioinformatics_external.py -q
```

实现参考：[Biopython 格式与 motif](https://biopython.org/docs/latest/Tutorial/chapter_motifs.html)、[PRDC](https://github.com/clovaai/generative-evaluation-prdc)、[clipped 指标论文](https://arxiv.org/abs/2507.01761)、[iLearnPlus 编码族](https://academic.oup.com/nar/article/49/10/e60/6154472)、[ANCOM-BC](https://doi.org/10.1038/s41467-020-17041-7)、[DNA Chisel specifications](https://edinburgh-genome-foundry.github.io/DnaChisel/ref/builtin_specifications.html)。公式约定以本页和 API docstring 为准。
