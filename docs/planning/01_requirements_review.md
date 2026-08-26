# DNAKit 阶段 1：需求审查报告

!!! note "历史阶段快照"
    本页保留阶段 1 时的需求审查和当时状态语义，不是当前实现清单。当时的 183 项状态已扩展为当前 184 项，现行状态以[追踪矩阵](requirements_traceability.csv)为准；API/测试/文档映射见[共享证据索引](07_requirements_evidence_index.md)，运行结论见[阶段 4/5 最终报告](06_stage4_stage5_delivery_report.md)。

审查日期：2026-08-13  
需求基线：`DNAKit_完整功能与输入表.csv`  
基线 SHA-256：`89453c851657134e04f80b7d040ba97a77f3e31ae34b67727351843bee7dd907`

## 1. 结论

CSV 可以完整解析，共 **184 条逻辑需求、14 个一级模块、7 列**。所有字段均非空，没有畸形记录，也没有二级功能或“DNAKit具体功能”完全同名的重复项。文件有 378 个物理文本行，是因为“所需输入”单元格内包含换行，不能按物理行数统计需求。

当前 CSV 适合作为功能愿望清单，但还不能直接作为验收规范。主要问题不是缺功能，而是：

1. 多个跨模块功能共享同一底层算法，若逐行独立实现会重复和不一致。
2. identity、coverage、复杂度、Gap、IUPAC、热力学条件、指纹维度等关键口径尚未冻结。
3. “建议实现类型/论文定位”有 33 种自由文本值，尚未真正落实为 `native`、`adapter`、`reimplementation`、`novel`。
4. NUPACK、Dashing、Primer3、参考数据库等存在安装、平台或许可证约束，不能全部打进默认 wheel。
5. “重点创新/核心创新”目前只能视为候选研究方向，不能在缺少文献检索、数学定义、消融和 benchmark 时直接宣称创新。

本轮不实现功能，不修改原始 CSV，不发布 GitHub、TestPyPI 或正式 PyPI。

## 2. 已确认的产品边界

- DNAKit 只处理 DNA 序列及其确定性、统计性、热力学和参考库检索分析。
- 不集成启动子活性、表达量、TF 结合强度、CRISPR 编辑效率等任务型深度学习预测模型。
- 内部与外部算法都必须记录方法、参数、软件版本、数据库版本和警告。
- 当前只认证 Linux；不要求 Windows 和 macOS 测试。
- 网站演示使用固定输入和预生成结果，不提供用户输入、上传或在线计算后端。
- 当前只做本地开发、测试、文档构建和包构建；所有外部发布均设置为后续人工审批门。

## 3. 模块规模与稳定需求 ID

原 CSV 不修改。规划阶段按“模块缩写 + 模块内顺序”分配稳定 ID。184 项逐项阶段、实现类型、状态、共享内核和范围说明见 [`requirements_traceability.csv`](requirements_traceability.csv)；最终 API、测试和文档映射由[共享证据索引](07_requirements_evidence_index.md)提供。矩阵的 `csv_table_row` 是包含表头在内的 CSV 表格行号，因此取值为 2–185；它与 1–184 的需求序号不是同一字段。

矩阵的 `default_label` 表示默认方法对外公布的四类标签，`alternative_label` 表示同一功能可选择的另一实现方法。例如 Tm 的 Primer3 方法是 `adapter`，独立 nearest-neighbor 方法是 `reimplementation`；不能只在模块级笼统标一次。

| 模块 | 数量 | 计划 ID |
|---|---:|---|
| 1. 核心数据对象 | 8 | `CORE-001`–`CORE-008` |
| 2. 文件读写与数据管理 | 6 | `IO-001`–`IO-006` |
| 3. 标准化与合法性检查 | 9 | `STD-001`–`STD-009` |
| 4. 序列内操作 | 9 | `OPS-001`–`OPS-009` |
| 5. 基础序列描述符 | 12 | `DESC-001`–`DESC-012` |
| 6. 序列模式与结构注释 | 13 | `PAT-001`–`PAT-013` |
| 7. 理化与热力学性质 | 13 | `THERMO-001`–`THERMO-013` |
| 8. 序列表征与特征工程 | 14 | `FP-001`–`FP-014` |
| 9. 相似度、比对与搜索 | 16 | `SIM-001`–`SIM-016` |
| 10. 去重、聚类与数据划分 | 24 | `DATA-001`–`DATA-024` |
| 11. 分子生物学模拟与规则设计 | 12 | `MOLBIO-001`–`MOLBIO-012` |
| 12. DNA 综合评价体系 | 15 | `EVAL-001`–`EVAL-015` |
| 13. 可视化 | 11 | `VIZ-001`–`VIZ-006`、`VIZ-008`–`VIZ-009`、`VIZ-014`–`VIZ-016` |
| 14. 后端、性能与可复现性 | 17 | `ENG-001`–`ENG-017` |

## 4. 重复与共享实现审查

### 4.1 明确的语义重复或规格碰撞

| CSV 功能 | 问题 | 处理决定 |
|---|---|---|
| 通用 motif 搜索 / Motif search | IUPAC、正则、PWM、正负链和重叠匹配重复 | 一个扫描内核；模式注释和搜索 API 只做不同包装 |
| k-mer 统计 / k-mer 特征 | presence、count、frequency 重复 | 一个计数内核；数据集 vectorizer 增加词表与 TF-IDF |
| CRISPR PAM / CRISPR 候选扫描 | PAM/guide 扫描重复 | 拆为 `scan_guides()` 与可选 `search_off_targets()` |
| 表格格式写出 / 数据导出 | CSV、TSV、JSON、Parquet 重复 | 统一由 `io.tables` 实现，报告模块调用它 |
| Exact search / Subsequence search | Exact search 同时写“整条相等”和“返回位置”，语义冲突 | 前者只做整条记录相等；后者做精确子串定位 |

### 4.2 保留公共名称、共享底层结果

- 字母表检查、非法字符检测和 validity 共用 `ValidationReport`。
- 模糊碱基统计和 Ambiguity 评价共用 ambiguity kernel。
- 数据质量检查和 Sequence quality 共用基础 QC 指标。
- entropy、复杂度、低复杂度区域和 Complexity 评价共用版本化复杂度指标。
- 指定/随机突变和突变文库共用 mutation engine。
- 倒置重复只给字符串候选；hairpin 才给热力学结构，二者不能混写。
- Tm、盐修正、nearest-neighbor、ΔG/ΔH/ΔS 和 duplex stability 共用统一热力学参数与条件模型。
- 最近邻搜索只计算一次；Novelty、最近参考序列、Memorization、Reference similarity 消费同一检索结果。
- identity 近似去重和 identity 聚类共用分组内核，但代表选择和标签冲突策略不同。
- Uniqueness 与 Redundancy 共用同一份重复分组报告。
- 限制酶位点扫描和限制酶切共享同一酶数据库及切点定义。
- 标准化审计、热力学条件、参考库追踪和版本追踪统一进入 provenance 模型。

## 5. 实现前必须冻结的定义

### 5.1 核心对象、Gap 与坐标

- 内部坐标统一为 **0-based、半开区间 `[start, end)`**；负链仍保持 `start < end`。
- 环状跨原点 feature 使用复合区间，禁止用 `start > end` 暗示跨原点。
- `N` 默认是模糊碱基，不自动等同于 assembly gap；只有显式策略或 AGP 信息才能转换为 `Gap`。
- `[500 bp]` assembly gap、未知长度 gap 和 alignment 中的 `-` 是三种不同语义。
- 未知长度 gap 存在时，总长度和绝对坐标可能不可解析；相关操作必须返回未解析位置或抛出明确异常。
- `strict` 与 `iupac` 是字母表策略；`gapped` 是序列表示能力，不再作为第三种字母表。
- `strandedness=double` 表示保存一条参考链并隐式得到完全互补链；含 mismatch/overhang 的真实 duplex 使用高级结果对象。
- `DNASequence.__eq__` 只比较规范化 parts 和类型元数据；反向互补或环状旋转等价必须由显式去重策略处理。
- FASTQ 要求无 Gap，质量数组长度与序列符号数完全相同；N 和其他 IUPAC 字符也各自对应一个质量值。A/C/G/T 数、模糊字符数和含 Gap 的坐标跨度分别统计。

### 5.2 描述符和模式

- “有效碱基长度”、GC/AT、CpG O/E 必须固定分母以及 IUPAC 的 `error/ignore/mask/probabilistic` 策略。
- 滑动窗口必须固定尾部窗口的 `drop/partial/pad` 行为。
- complexity、重复比例、低复杂度必须指定算法、阈值、重叠合并和比例分母。
- ORF 必须固定完整/部分 ORF、嵌套 ORF、替代起始密码子、最短长度及环状跨原点规则。
- CpG island、倒置重复、串联重复和微卫星必须固定匹配、合并与排序规则。

### 5.3 热力学和结构

- 每次计算必须保存算法、参数集版本、单位、温度、链浓度、盐模型、标准态和 IUPAC 行为。
- hairpin/dimer 必须定义“最佳结构”、3′ 风险、并列结构和排序规则。
- 二级结构必须声明是否支持 pseudoknot、多链复合物以及 ensemble 条件。
- “局部熔解特征”当前实际是窗口 Tm 近似，不应命名为实验意义上的完整 melting profile。

### 5.4 指纹、相似度、聚类和划分

- 每种指纹必须有不可变的 schema/version：维度、词表或面板、空间分箱、哈希、聚合、缺失值和标准化方法。
- identity 与 query/target coverage 必须固定 terminal gap、alignment gap、IUPAC 和分母公式。
- Tanimoto 必须区分二值、非负计数和一般连续向量；带负值向量默认不允许 Tanimoto。
- sketch 必须记录 k、seed、scaled/num、canonical 和估计误差；最近邻结果必须标记 exact 或 approximate。
- IUPAC compatibility 不具传递性，不能直接用普通连通分量当作等价类。
- 近似去重和阈值聚类必须明确 single-linkage、complete-linkage 或 greedy representative。
- split 必须定义比例取整、小 strata、不可满足约束、超大连通分量和松弛策略。
- 使用近似索引做泄漏检测时必须明确可能漏报，不能声称严格无泄漏。

### 5.5 模拟、评价和输出

- ligation、PCR、assembly 必须定义方向、磷酸化、环化、多结合位点、多产物排序和最大枚举规模。
- 规则优化和 codon optimization 必须定义目标函数、约束优先级、不可行处理和确定性 tie-break。
- diversity、distribution similarity、synthesis risk 和综合评分都必须公开公式、范围、阈值和权重。
- synthesis risk 只能表示透明规则或统计风险，不能等同实验合成成功率。
- 600 dpi 只适用于 PNG/TIFF 等栅格输出；SVG/PDF 按矢量结构验收。

## 6. `native / adapter / reimplementation / novel` 判定

分类必须落实到具体方法，而不是只标整个模块。

| 类型 | 定义 | DNAKit 中的典型内容 |
|---|---|---|
| `native` | DNAKit 自己设计的框架、组合逻辑或简单确定性算法 | 核心对象、标准化流程、简单计数、exact/RC 去重、split orchestration、CLI、provenance |
| `adapter` | 调用外部 Python 库、CLI 或参考数据库，不复制其算法 | Biopython I/O、Primer3、BLAST+、MMseqs2、Dashing、sourmash、参考数据库 |
| `reimplementation` | 按公开标准、论文或已知算法独立实现 | IUPAC、坐标转换、Shannon entropy、ORF、DP alignment、nearest-neighbor 等 |
| `novel` | 数学定义和实现由 DNAKit 新提出，并有检索、消融和 benchmark 证据 | 当前没有已确认项 |

四类标签互斥，优先级固定：外部实现为 `adapter`；内部复现公开算法为 `reimplementation`；经证据确认的新方法为 `novel`；其余为 `native`。执行方式和研究来源另存为正交字段，`integration` 不是第五种公开标签。

混合指纹、多尺度指纹、多约束划分、泄漏评价和合成风险评分当前公开标签仍为 `native`，研究状态可记为 `candidate-novel`；论文实验完成前不得标为 `novel`。

## 7. 版本范围建议

### 7.1 MVP（计划版本 `0.1.x`）

- 核心：DNASequence、DNARecord、DNASet、DNAFeature、Gap、Location、MetricResult、Issue、Provenance。
- I/O：FASTA、FASTQ、CSV/TSV/JSON、gzip、流式读取和基础元数据；GenBank、GFF3/BED/AGP、Parquet 只做后续完整支持。
- 标准化：大小写/空白/U、strict/IUPAC、非法字符、ambiguity、基础 Gap、QC、审计。
- 操作：reverse/complement/RC、转录翻译、切片、编辑、突变、拼接、trim/mask；复杂环状与 feature 同步后移。
- 描述符：长度、组成、GC/AT、skew、CpG、k-mer、entropy、homopolymer、窗口和基础 codon 统计。
- 表征：整数、one-hot、k-mer、canonical k-mer 与固定 schema 的基础指纹。
- 相似度：exact/subsequence/RC、Hamming、edit distance、k-mer/基础指纹相似度、小规模矩阵。
- 数据集：exact/RC 去重、冲突报告、随机/分层/单字段 group/basic similarity split。
- 可视化：序列文字图、高亮、Gap 和相似度矩阵已有原生 SVG；PNG/TIFF/PDF、任意结果表和 Parquet 导出仍为部分/计划能力。
- 工程：Python API、CLI、批量、进度、流式、seed、错误、版本追踪、pytest、类型检查、README 和本地 MkDocs。

### 7.2 高级版本（计划版本 `0.2.x`–`0.4.x`）

- 完整注释格式、索引访问、环状与 feature 编辑同步。
- motif/PWM/ORF、限制酶、repeat、CRISPR 候选扫描。
- Primer3 热力学适配与独立 nearest-neighbor 实现。
- alignment、BLAST/MMseqs2/Dashing/sourmash 搜索、近似去重和聚类。
- 多约束 split、leakage、novelty、memorization、distribution similarity、synthesis risk。
- 综合/混合/多尺度指纹及其 schema、消融和 benchmark。
- 酶切、PCR、ligation、assembly、规则优化和突变文库。
- 交互 HTML 报告、静态只读网站演示、benchmark 和论文复现实验。

### 7.3 条件支持或无法按“单一纯 Python 包”承诺

- 在后端支持的长度、模型和资源范围内，可靠 hairpin、dimer、MFE 和配对概率需要外部后端；不承诺任意长度。
- 超大参考库搜索、identity 聚类、复杂 repeats 需要外部 CLI 和独立数据库。
- 任意规模的全相似度矩阵本身是 `O(n²)` 输出，必须限制规模或改用 Top-k。
- 多约束划分可能无可行解，不能保证全局最优，只能返回可行解、松弛结果或冲突报告。
- novelty 只能相对于明确版本和过滤条件的参考库定义，不能输出“绝对新颖”。
- PCR、assembly、codon optimization 和 synthesis risk 不能预测实验产率、表达效果或实际合成成功率。
- 不实现任务型深度学习预测，也不输出法律意义上的 `legal=True/False`。
- NUPACK 只提供被动探测与显式 adapter，不自动安装/下载；真实执行和对照受单独许可与独立安装条件限制，处理方式见依赖规划文档。

## 8. 验收规格需要补充的字段

正式编码前，每条需求至少补充以下字段：

- 稳定需求 ID、计划版本、状态：`planned/in_progress/complete/partial/conditional/unavailable/blocked`；每个单元格只允许一个状态，方法级差异写入范围说明。
- Python API、CLI 子命令、输入类型、输出 schema。
- 实现类型、默认算法、默认参数、单位和确定性规则。
- 必需/可选依赖、最低/最高兼容版本、许可证和引用。
- 后端缺失时的失败或降级行为。
- 边界、最大建议规模、复杂度和内存策略。
- 单元测试、差分测试、性能测试和明确容差。
- 已知限制、数据库版本和 provenance 要求。

184 条原需求能力全部进入追踪矩阵；语义重复项共享实现，不删除用户可见能力。三处输入契约做了显式细化：DNASequence 的名称/描述归 DNARecord，Gap 上下文由 DNASequence parts 表达，DNASet 的文件/迭代器输入使用显式工厂和 RecordSource。实现前需审查并冻结这些映射。
