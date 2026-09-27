# 外部分析功能的安装简化与内置路线

审查日期：2026-09-21。范围是 `src/dnakit/integrations/analyses.json` 中的 **43 个分析 ID**；逐项检查了输入、实际命令和输出约定，并对照现有 Python 实现。本文是迁移评估；表中的“建议”不代表功能已经完成，也不代表对应 extra 已经发布到 PyPI。

目标是用户安装 DNAKit 或指定 extra 后即可计算。优先使用 Python 计算库和已有 API；复杂算法保留专业实现。由 pip 安装的程序仍可通过子进程运行，用户无需另行准备可执行文件，但这与进程内计算库是两种实现方式。

## 本轮已提供

当前源码已提供 `dnakit.genomics.genotype_pca()` 与 `average_nucleotide_identity()`，通过 `genomics` extra 安装 NumPy、Pyrodigal、pysam 与 PyFastANI。前者是边界明确的基因型矩阵 PCA，后者使用 PyFastANI 计算库；原有 `plink2-pca` 和 `fastani` 外部命令入口保留。安装与示例见[基因组教程](../api/features/38_genomics_extensions.md#pip-genomics)。这些新增能力尚未发布到 PyPI。

新接口的 80 项单元测试已通过，并在隔离 Python 3.11 环境使用真实 NumPy、PyFastANI、Pyrodigal 和 pysam 完成人工数据调用验证。仓库内运行 `python others/validation/pip_genomics.py` 可复核 PCA、ANI、CDS 预测与 BAM 统计；该脚本拦截 `subprocess.Popen`。验证范围是接口和人工数据，不代表真实样本的科学准确性或全部平台兼容性。

现有基因预测、VanRaden 关系矩阵与育种预测也可共用这些依赖。基因型数据、组装 FASTA 和参考资源仍由用户提供。

## 分类与优先级

每个分析 ID 只计入一个主类：**A：3 项，B：8 项，C：6 项，D：26 项，共 43 项**。

| 类别 | 含义 | 实施原则 |
| --- | --- | --- |
| A | 可评估通过 pip 计算库或程序包提供 | 优先验证 wheel、Python 版本、结果和安装后的程序定位 |
| B | 已有相关的 DNAKit API | 复用组件；明确模型、输入或输出的差异，不能直接宣称等价 |
| C | 适合新增边界明确的独立功能 | 使用清楚的新 API 名称，限定数据规模和算法语义 |
| D | 暂保留专业后端 | 先解决发行包与资源管理，不按相似名称重写算法 |

P0 表示优先实施；P1 表示后续限定范围实现或集成；P2 表示需专门的算法、打包或数据资源项目。优先级是工程建议，不是发布日期承诺。

## 43 项逐项评估

| 分析 ID | 类别／优先级 | 源码核对与建议边界 |
| --- | --- | --- |
| `plink2-pca` | B／P0 | 已提供独立 `genotype_pca()`，明确频率标准化、缺失值和不变位点规则；输入为剂量矩阵，协方差特征值按样本数缩放，保留原外部入口的 BED／eigenvec／eigenval 契约。通用 `analyze_features(method="pca")` 则继续使用通用特征缩放。 |
| `plink2-king` | C／P1 | 当前调用 `--make-king-table`。可独立实现限定于二倍体双等位常染色体硬调用的 KING-robust；需按样本对处理缺失、杂合与 IBS0，排除分母为零情形，并对照 PLINK。已有 VanRaden 矩阵不是 KING。 |
| `plink2-ld-prune` | C／P1 | 当前是按位点数指定窗口／步长的 `--indep-pairwise`。可提供独立的窗口内相关性剪枝；剪枝顺序、并列规则和缺失位点集合会改变保留集合，未经逐项对照不能复用 PLINK 的方法标识。 |
| `plink2-gwas` | C／P1 | 当前 `--glm` 根据表型编码选择回归。建议先提供有协变量的定量性状逐位点线性检验；二分类、完全分离／Firth、样本筛选和多表型规则另行实现验证，不能声称覆盖整个 `--glm`。 |
| `gcta-grm` | B／P1 | 已有 `genomic_relationship_matrix()`，源码是 VanRaden method 1，分母为所有位点 `2p(1-p)` 之和。GCTA 的标准化估计、样本对有效位点数与二进制输出需单独实现；现有矩阵不能直接冒充 GCTA GRM。 |
| `gcta-reml` | D／P2 | 当前读取 GCTA GRM、定量表型和协变量运行 GREML。已有固定惩罚的 `genomic_prediction()` 不估计 REML 或遗传力；保留后端，未来专门验证收敛、边界解和标准误。 |
| `minigraph-pangenome` | D／P2 | 当前是参考加一个组装的 `-cxggs` 建图。保留图构建算法；FASTA 合并或普通序列比对不能替代 GFA 图构建。 |
| `hifiasm-assembly` | D／P2 | 当前单样本 HiFi 默认组装、输出 GFA。保留专业组装器；后续可评估独立平台 wheel，不能用简化拼接替代。 |
| `flye-assembly` | D／P2 | 当前指定 `--nano-raw`。上游有 Python 安装入口，但源码安装仍涉及编译后端；先验证完整发行包，不能仅声明 Python 包就保证无系统构建依赖。[安装说明](https://github.com/mikolmogorov/Flye/blob/flye/docs/INSTALL.md) |
| `sniffles-sv` | A／P1 | 上游支持 pip 安装。可用 extra 自动安装 Sniffles 及其依赖，继续执行同一 `sniffles` 入口，保持长读长单样本 SV 语义；这是安装自动化，不是重新实现 SV 算法。需验证依赖版本和虚拟环境入口定位。[官方说明](https://github.com/fritzsedlazeck/Sniffles) |
| `modkit-pileup` | B／P2 | 已有 `read_bedmethyl()`／`summarize_methylation()`，只处理统计结果。当前任务从含 MM/ML 标签的 BAM 生成 pileup，涉及概率阈值、修饰代码和覆盖定义；继续保留 modkit，不把下游汇总当作替代。 |
| `mosdepth-coverage` | B／P1 | 已有 pysam 实现的 `coverage_consensus()`，限制区域长度，统计过滤后的 A/C/G/T 深度。可新增明确过滤规则的区间覆盖 API；mosdepth 的配对重叠、缺失／跳跃 CIGAR、全基因组分布和文件格式要分别验证。 |
| `fastani` | A／P0 | 已提供 `average_nucleotide_identity()`，通过 PyFastANI 直接使用 FastANI 内部算法，免 `fastANI` 可执行文件。记录 query/reference 方向、命中及总片段数；无命中返回 null ANI。[计算库](https://pyfastani.readthedocs.io/en/stable/)／[安装范围](https://pyfastani.readthedocs.io/en/stable/guide/install.html) |
| `kraken2-classify` | D／P2 | 当前对显式本地 Kraken2 库分类读段。保留分类器和数据库版本；通用 k-mer 相似度不能直接替代分类规则或物种丰度。 |
| `amrfinder-nucleotide` | D／P2 | 当前是组装核酸加显式数据库的 AMRFinder。保留参考数据及专业解释规则；一般相似序列搜索不能替代耐药注释流程。 |
| `prokka-annotation` | B／P2 | 已有 Pyrodigal CDS 预测，可供“仅基因预测”使用。Prokka 当前返回整套注释文件，包含超出 CDS 定位的步骤；完整功能注释仍保留原后端。 |
| `barrnap-rrna` | D／P2 | rRNA 识别需模型和命中过滤规则。先固定上游版本及模型；当前上游已包含 Infernal／Rfam 路线，不能仅用 PyHMMER 宣称复现全部版本。[上游说明](https://github.com/tseemann/barrnap) |
| `fimo-motifs` | A／P1 | 可评估 LightMotif 的 pip 计算库；已有 `scan_motif_model()` 可做基础 PSSM 扫描。现有 FIMO 契约按位点 p 值筛选，模型背景、伪计数、双链、p/q 值及输出格式均需对照。PyMEMEsuite 上游仍标为实验性并提示改用 LightMotif，其内含 MEME 组件有独立条款，不默认纳入。[LightMotif](https://lightmotif.readthedocs.io/en/stable/)／[PyMEMEsuite](https://github.com/althonos/pymemesuite) |
| `cas-offinder` | C／P1 | 可新增 CPU 上限定规模、仅错配且无 bulge 的 PAM 匹配／脱靶候选枚举。需明确双链与坐标，限制输入规模；不冒称完整 Cas-OFFinder 设备和规格文件兼容实现。 |
| `bakta-annotation` | B／P2 | Pyrodigal 只覆盖 CDS 预测组件。当前 Bakta 使用专用数据库，返回完整注释；功能转移、RNA 和数据库匹配流程仍需原后端。 |
| `busco-proteins` | D／P2 | 当前在显式 lineage 上执行离线蛋白评估。PyHMMER 可供 HMM 搜索组件使用，但 BUSCO 的分类、重复处理和谱系阈值仍需独立验证，不能把 HMM 命中比例叫 BUSCO 完整度。[PyHMMER](https://pyhmmer.readthedocs.io/en/stable/) |
| `busco-genome` | D／P2 | 除谱系数据还涉及基因预测路线。保留完整 BUSCO；原核 CDS 预测并不覆盖真核基因组模式。 |
| `augustus-genes` | D／P2 | 当前明确要求物种模型与配置目录。保留真核基因结构预测；不能由 ORF 检索或 Pyrodigal 替换含内含子的基因模型。 |
| `diamond-blastp` | D／P2 | 当前先从用户蛋白 FASTA 建 DIAMOND 库，再输出 12 列相似性命中。保留搜索引擎；全对全动态规划可以作为小规模新功能，不能宣称同等搜索性能和 E 值。 |
| `diamond-blastx` | D／P2 | 当前翻译核酸搜索蛋白数据库。还需保留阅读框、翻译和命中统计语义，不能用普通核酸比对代替。 |
| `eggnog-proteins` | D／P2 | 当前 emapper 蛋白模式依赖 DIAMOND 和版本匹配的 eggNOG 数据。保留注释流程；数据库／软件兼容版本必须一起记录。 |
| `miniprot-align` | D／P2 | 当前蛋白到基因组比对输出 GFF，支持剪接和移码。保留专业比对器；一般蛋白比对不具备这些基因组坐标语义。 |
| `liftoff-annotation` | D／P2 | 当前读取参考注释数据库，在目标组装上转移注释。Python 包本身不代表全部比对依赖已解决；保留迁移算法及 unmapped／低相似度标记。 |
| `checkm2-predict` | D／P2 | 当前输入 bins 和版本匹配的 DIAMOND 数据库，输出模型完整度／污染度估计。需管理搜索引擎、模型和参考数据；基础长度／GC 统计不能替代。 |
| `quast-assembly` | D／P2 | 当前是有参考组装评估。可另增纯 Python 的长度、N50/L50、GC 和 N 含量统计，但它们不覆盖参考比对、错组装和完整 QUAST 报告。 |
| `vsearch-dereplicate` | B／P0 | 已有 `deduplicate(equivalence="exact")`。可复用精确序列分组并补 FASTA 读写／计数；需明确大小写、输入 abundance、最短长度和输出排序，UC 文件和 VSEARCH 标签语义不能自动宣称兼容。 |
| `vsearch-cluster` | B／P1 | 已有 `cluster_sequences()`，源码是全对全阈值图的连通分量。当前 VSEARCH 是按长度排序的贪心聚类，二者会产生不同簇；保留不同的方法名，按需求新增独立的贪心模式。 |
| `mmseqs-protein-search` | D／P2 | 当前 `easy-search --search-type 1`。保留蛋白搜索、预过滤和显著性算法；DNA 相似度 API 不覆盖此任务。 |
| `mmseqs-protein-cluster` | D／P2 | 当前 `easy-cluster` 显式设置 identity、coverage 和 `--cov-mode 1`。保留蛋白覆盖定义与聚类策略；现有 DNA 连通分量聚类不等价。 |
| `abricate-screen` | D／P2 | 当前显式本地库的耐药／毒力序列匹配。保留数据库版本、覆盖和一致性语义；可复用结果解析，但结果解析不是筛查计算。 |
| `mlst-typing` | D／P2 | 当前匹配已安装方案及 BLAST 数据库并输出 ST。可以规划确定等位基因的精确匹配子集，但完整自动方案选择、部分命中和新等位基因处理仍保留专业后端。 |
| `genomad-mobile-elements` | D／P2 | 当前运行 end-to-end 且依赖显式数据库／模型。保留分类模型、标记证据和阈值；pip 安装 Python 代码不能取代模型资源配置。 |
| `repeatmasker-library` | D／P2 | 当前使用指定重复库和配置好的搜索引擎。保留重复注释流程；单纯低复杂度屏蔽只可作为独立子功能。 |
| `mageck-test` | D／P2 | 当前是两组非配对 RRA、median 归一化。保留完整 MAGeCK 检验流程；简单计数差异或均值比不是等价的基因显著性检验。 |
| `mageck-count` | C／P1 | 可提供单样本、固定长度 sgRNA、显式截取位置的精确计数 API。明确重复 guide、反向互补、低质量和未匹配读段；MAGeCK 的自动截取／多长度等行为需后续逐项验证。 |
| `crispresso2-amplicon` | D／P2 | 当前单扩增子、单 guide 的编辑定量。保留比对、切割窗口和 indel 归属规则；错配计数不能直接称作编辑效率。 |
| `drugz-unpaired` | C／P1 | 上游是可读的 Python 实现。可评估固定版本内置适配，或对照原实现独立实现非配对模式；需保持归一化、方差窗口、guide 汇总和 FDR。若复制源码，保留来源版本与原始通知。[上游源码](https://github.com/hart-lab/drugz/blob/master/drugz.py) |
| `gemma-lmm` | D／P2 | 当前 `-lmm 4` 输出 Wald／似然比／score 检验。固定 lambda 的育种预测或普通线性回归不覆盖亲缘校正的关联检验；需专门的混合模型验证。 |

## 已有可用组件与后续交付

审查开始时已存在的实现：

| 已有 Python API | 依赖入口 | 已核对的实现范围 |
| --- | --- | --- |
| `predict_prokaryotic_genes()` | `genomics-genes` | Pyrodigal 进程内 CDS 预测，不是完整功能注释 |
| `genomic_relationship_matrix()`／`genomic_prediction()` | `bio-analysis` | VanRaden method 1／显式固定 lambda 的 ridge 预测 |
| `analyze_features(method="pca")` | `bio-analysis` | 通用矩阵 PCA 和特征缩放 |
| `build_motif()`／`scan_motif_model()` | `bio-analysis` | PSSM 建模和按分数阈值扫描，FPR 阈值不是全基因组 q 值 |
| `read_bedmethyl()`／`summarize_methylation()` | 基础安装 | bedMethyl 读取、验证和计数汇总 |
| `coverage_consensus()` | `bio-reads`、`bio-analysis` | 指定区域的 A/C/G/T 覆盖和多数碱基共识 |
| `deduplicate()` | 基础安装 | 规范化 DNA 记录的精确／反向互补等去重与审计 |
| `cluster_sequences()` | 由所选相似度方法决定 | 有规模限制的全对全阈值图连通分量聚类 |

源码依据：`src/dnakit/genomics/{genes,breeding,methylation}.py`、`src/dnakit/bioinformatics/{features,motifs,reads}.py`、`src/dnakit/datasets/{deduplicate,clustering}.py`。使用示例见[基因组功能教程](../api/features/38_genomics_extensions.md)和[生物信息学接口说明](bioinformatics.md)。

新增实现应使用独立 API 或显式 backend 选择，返回计算方法、版本、输入处理规则和引用。保留现有外部 ID 的可复现命令，不根据机器上是否有工具而静默改变统计方法。只有安装、调用和验证完成后，才将对应项从“建议”更新为“已提供”。

## 判定“安装即可用”的条件

1. 在声明支持的 Python／操作系统／架构上做干净环境安装和真实调用；存在 sdist 不代表用户无需编译器。尤其 PyFastANI 上游明确说明部分平台需要从源码构建。
2. 对照上游结果，验证坐标、单位、缺失数据、边界、无命中和失败行为；PCA 比较子空间／符号不变性质，不能要求特征向量逐位相同。
3. 将人工小数据接口验证与科学准确性、真实数据性能分开报告；现有外部命令能执行不代表新实现已科学验证。
4. 数据库、模型和参考基因组单独记录版本、校验和、来源、大小与缓存位置。按需下载应由明确操作触发；不能把软件包安装成功写成数据库已就绪。
5. 纳入上游源码时记录固定 commit、改动、许可证与第三方组件来源，建立更新方式；成熟的可安装计算库优先于复制一份长期维护的源码。

本轮落地上述两个新 API 和组合安装选项；其余表项保留逐项建议和边界，不表示已经全部内置。
