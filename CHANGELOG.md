# Changelog

本项目遵循语义化版本。`0.1.0.dev0` 仍是本地开发快照，公开接口在首个正式版本前可能调整。

## [Unreleased]

### Added

- 面向普通用户的统一 `dnakit.DNA(...)` 门面：同一个入口接收单条/多条序列、逐条 record mapping、ID、拓扑、metadata、feature 和旧核心对象；下标与切片仍返回 `DNA`。
- `dnakit.representations` 可选择 11 种 DNA 基础模型提取固定长度 rep，checkpoint 默认下载到当前目录 `ckpt/` 并复用；`neural_cluster_sequences()` 完成 L2、可选 PCA、seed 固定的 k-means、中心最近代表序列和结果审计。
- `evaluate_frechet_distance()` 默认复用 DATA-027 的 LucaOne 表征和 L2 归一化，以样本空间等价算法计算两个 DNA 集合的 Fréchet 表征分布距离。
- `evaluate_fragment_similarity()` 和 `evaluate_snn()`：分别以 exact k-mer 出现次数余弦相似度和 hashed k-mer 指纹最近邻 Tanimoto 均值实现 MOSES Frag/SNN 的 DNA 适配。
- `dnakit.search` 公共数据库层：NCBI Datasets/Entrez/BLAST、Ensembl、ENA、ENCODE、UCSC 的有界查询、坐标标准化、凭据脱敏、分页和 provenance。
- `dnakit.download` 公共数据层：事务式多文件安装、checksum manifest、基因组/基因/病毒/ENA reads/ClinVar/dbSNP/GEO/ENCODE/UCSC 下载、查询元数据五格式导出及显式 BLAST/BWA/Bowtie2/minimap2/samtools 索引构建。
- 附件 29 项查询和 24 项下载的独立能力矩阵、API 教程与远程/受控数据边界。
- `dnakit.datasets.exclude_species()` 和 `exclude_chromosomes()`：按显式 metadata 精确排除一个或多个物种/染色体；缺少目标字段时结构化报错。
- MkDocs 文档站启用中文站内搜索框、搜索建议和结果关键词高亮。
- `ValidationConfig(sequence_length=...)`、`dnakit validate --sequence-length` 和工作流 `validate` 步骤现可按标准化后的精确 `symbol_length` 判定序列是否合法。
- 固定 `descriptor_schema_v1` 的 240 项 DNA 描述符、字段元数据、统一 `all_descriptors()` 结果、CLI 完整输出及 DiProDB 15 组公开参数汇总。
- 根目录新增 MIT 许可证文件，版权人标注为 Pengsen Ma，并将许可证纳入发行归档。
- 有界 GenBank 子集、GFF3、BED3–6、AGP 2.1、普通 FASTA/严格四行 FASTQ 持久索引、坐标提取、分块迭代和元数据管理。
- 以显式 `TableSchema` 有界读写 CSV/TSV/JSON/Parquet 的 `read_table()`、`export_table()` 和不可变审计结果。
- 复杂度/重复描述符，motif/PWM、ORF、限制酶、PAM、CpG island、回文、倒置重复、串联重复、STR 和低复杂度扫描。
- 分子量、基于公开 nearest-neighbor 吸光参数的未修饰 ssDNA 260 nm 理论消光系数，以及内部复现的 SantaLucia 1998 Tm、盐修正、nearest-neighbor、ΔG/ΔH/ΔS、堆积、完整互补 duplex 和局部 Tm；`duplex_stability()` 可显式使用 primer3-py adapter 评价 canonical mismatch/dangling heterodimer。
- 单/双链光学性质、OD260/A260/摩尔浓度/质量浓度/物质的量/质量换算、显式染料/修饰修正，以及 Ka/Kd、双链比例、理论熔解曲线、末端稳定性、Na⁺+K⁺ 和 DMSO/甲酰胺经验修正。
- dot-bracket 二级结构与配对概率派生指标、单独许可安装后才执行的显式 NUPACK 4 adapter，以及 PDB 坐标几何、SASA、体积、形状、骨架二面角、近似螺旋参数、NMR ensemble RMSF 和 3DNA/DSSR 解析。
- RCSB 1BNA、1AC7、139D 本地校验样本、SHA-256 清单、带进度条的结构分析脚本和机器可读结果。
- 显式可选的 `Primer3PyAdapter`，支持 Tm、hairpin、self-dimer 和 heterodimer，映射 Na/Mg/dNTP/strand concentration/temperature 条件，并记录后端版本、许可和 provenance；下游统一结果会校验序列、条件和结构选项绑定。
- MinHash/FracMinHash、motif/限制酶/GC/repeat/coding/固定 16 维热力学/混合/多尺度指纹及训练集拟合的预处理器；结构特征使用显式 adapter 或固定缺失策略。
- 有界近似匹配、global/local/semi-global 线性/仿射 gap 比对、sketch 相似度、持久索引和 Top-k 最近邻。
- 环状/IUPAC/近似去重、identity/edit/k-mer/fingerprint 聚类、层次聚类、代表序列、时间划分、多约束启发式划分、leakage 和划分质量评价。
- `SplitConfig(method="hash")` 顺序无关的稳定哈希划分：使用 SHA-256、seed 和唯一 `record.id`，并保持精确比例配额。
- validity、ambiguity、quality、complexity、uniqueness、diversity、redundancy、reference-scoped novelty/memorization、分布相似度、synthesis-risk 和规则型 scorecard。
- 酶切、末端/连接、PCR、引物匹配、可选结构引物属性、严格设计请求及显式 `Primer3PyDesignAdapter`、Gibson/LCR/Golden Gate/BioBrick 序列级组装、CRISPR 候选、规则优化、密码子优化和突变文库；Primer3 设计候选会反查模板坐标、左右引物序列和产物长度。
- 线性/环状 feature 图、alignment 图、相似度矩阵、自包含 HTML 报告，以及可选 PNG/TIFF/PDF 600 dpi 导出。
- 后端 registry、BLAST/MMseqs2/sourmash 被动 metadata/显式版本句柄、要求调用方提供可执行文件的严格 Dashing Jaccard/Top-k adapter、内容寻址 JSON 缓存、线程批处理、稳定 resume、JSON/YAML 多步骤工作流和可审计运行清单。
- `describe`、`fingerprint`、`search`、`orfs`、`compare`、`report` 和 `workflow` 等统一 CLI 能力，以及兼容开发入口 `python -m dnakit.cli.workflow run`。
- 本地正确性验证器、带 seed/环境/参数/逐样本记录的 microbenchmark、DNAKit/Biopython 对等任务比较、所选公开 callable 源码行口径、高级 Notebook、完整工作流示例和阶段 4/5 交付报告。
- TestPyPI 和正式 PyPI 的手动工作流配置；构建审查 job 无发布凭据，只有上传 job 具有 `id-token: write`，二者均带许可证、版本与显式确认门禁，本阶段未触发。

### Changed

- 文档将 FP-004 Canonical k-mer 并入 FP-003 k-mer 特征的 Canonical 模式，并将 FP-005 统一命名为 k-mer Sketch；公开 Python API 和稳定追踪编号保持兼容。
- 普通文件读取统一为 `read(..., mode="dna"|"stream")`，验证统一为 `validate()`；无后缀编辑、反向互补和环状操作对 `DNA` 输入同步注释并返回 `DNA`，旧 `read_one/read_set/validate_set/*_record` 入口继续兼容。
- DNA 基础模型表征和神经网络聚类的默认 checkpoint 从 GROVER 改为 `LucaGroup/LucaOne-gene-step36.8M`；默认缓存目录为 `ckpt/lucaone-gene-step36-8m/`，自定义 checkpoint 代码仍需显式授权。
- `normalize()` 新增 `keep_ambiguous`、`keep_u`和 `keep_other`；默认保留 IUPAC 模糊碱基，审计删除 `U` 和其他非 DNA 字符。
- 需求追踪矩阵改为逐项按当前源码、测试和本地后端状态核对，不再保留过期 `planned` 声明。
- 将独立的 identity、query coverage 和 target coverage 功能项合并到 `SIM-008` 序列比对结果，并删除重复编号。
- 批处理支持有界线程执行、稳定原始索引/派生 seed 和已完成 ID 跳过。
- 可视化范围收缩为序列、feature、alignment 和相似度矩阵，通用图片导出和 HTML 报告保持不变。
- 可视化导出扩展为原子多格式输出，并对像素、记录数和内嵌结果字节数设硬上限。
- README、安装、快速入门、API、FAQ、引用/许可证和文档站导航同步到阶段 4/5 真实状态。

### Validation

- 11 种模型 registry、checkpoint 复用/下载契约、IUPAC/Gap/矩阵边界和 rep → k-means/PCA 由单元测试覆盖；GROVER 使用真实本地 checkpoint 完成 4 条序列的 rep + k-means smoke，其余模型保持条件状态。
- 人工小序列、IUPAC、空序列、200,000 nt 输入和跨原点环状限制酶位点通过边界验证。
- 限制酶、分子量终端约定、全局比对、字面搜索和 single/complete/average linkage 与 Biopython 对照通过。
- DNAKit 内部复现的 SantaLucia Tm 与 primer3-py 在对齐条件下绝对差不超过 `0.02 °C`。
- `Primer3PyAdapter` 的 hairpin/self-dimer/heterodimer 与同一已安装 primer3-py 直接调用逐字段一致；这只验证 adapter 转换，不等价于独立算法验证。
- 独立集成审查另通过 Mg/dNTP 条件的 Primer3 Tm/primer properties、固定 16 维热力学指纹、mismatch duplex 和 5 个真实设计候选 smoke，并确认伪造条件、候选序列/坐标及错配 synthesis-risk 序列会被结构化拒绝；这些 smoke 不扩充 17 项正式差分报告。
- Dashing adapter 的固定命令、矩阵解析、失败/超时/输出/输入突变边界由受控替身覆盖，本地 `v1.0.2-4-g0635` 两序列 exact 文档示例 smoke 通过；尚未做科学差分。
- primer3-py 2.3.0 按当前发行 metadata 保守记录为 `GPL-2.0-only`；它是显式可选依赖，不随 DNAKit 打包。
- NUPACK adapter 的字段映射、边界和错误处理已用受控替身测试；当前环境没有 NUPACK，因此仍不存在真实 NUPACK 数值差分结论。

### Remaining boundaries

- DNAKit 自身许可证已确定为 MIT；正式发行仍需完成依赖兼容性复核、版本准备和远程发布授权。
- 未推送 GitHub、未部署 Pages、未上传 TestPyPI/PyPI，也未执行论文复现实验。
- NUPACK 不随项目安装/下载，当前环境不可用且真实差分仍受单独许可与安装条件约束；Dashing adapter 尚无真实科学差分；Dashing/BLAST/MMseqs2/sourmash/外部数据库没有随包分发。
- PyArrow 25.0.1 下的 DNAKit Parquet 写出/`read_table()` 往返已通过；其他引擎/版本交叉矩阵尚未执行。

## [0.1.0.dev0] - 2026-08-13

### Added

- `src` 布局、Python 3.10+、类型注解、`py.typed` 和语义化开发版本。
- 核心对象、结构化错误、标准化、基础 I/O/操作/描述符/指纹/相似度/去重/划分/可视化。
- pytest、mypy、Ruff、构建、Conda 环境、GitHub Actions 配置和本地 MkDocs Material 网站。
- 184 项需求审查、架构设计、依赖分析、阶段计划、固定输入只读演示和 MVP Notebook。

### Limitations

- 该版本号仍代表开发快照，不是已上传的软件发行版。
