# 条件支持、部分能力与阻断项

本文件记录不能由 DNAKit 默认安装包独立、无条件完成的功能。`conditional` 表示公共接口存在但依赖外部条件；`partial` 表示有可用子集且缺口明确；`blocked` 表示当前不能在授权范围内完成。

## 条件功能

| 功能 | 状态 | 当前可用部分 | 条件/缺口 |
| --- | --- | --- | --- |
| Primer3 Tm/hairpin/self-dimer/heterodimer/design | `conditional` | `Primer3CLIAdapter` 与 `Primer3CLIDesignAdapter` 已实现显式路径、CLI 协议、有界解析和错误处理 | 用户按实际许可单独安装 Primer3，并提供 `oligotm`、`ntthal` 或 `primer3_core` 路径；DNAKit 无 `primer3` extra、不自动发现/打包；当前只有受控替身契约测试，不是科学差分 |
| NUPACK 二级结构/多复合物平衡 | `conditional` | `NupackAdapter` 已实现 MFE、pfunc、pairs、subopt、sample、ensemble size、目标指标和 tube concentration 的有界映射；dot-bracket/概率派生指标可独立运行 | NUPACK 需用户按单独许可独立安装；当前环境不可用，只有受控替身契约测试，没有真实数值差分 |
| Dashing sketch Jaccard | `conditional` | `DashingAdapter` 可用调用方显式提供的可执行文件计算 exact/HLL Jaccard 矩阵和稳定 Top-k；固定命令协议、资源限制、SHA/provenance 及故障测试完成；本地 v1.0.2-4-g0635 exact 示例 smoke 通过 | 未做科学差分；GPLv3 第三方目录从 wheel/sdist 排除，DNAKit 不自动发现、安装或选择它 |
| 二核苷酸数值描述符 | `conditional` | 固定 240 字段 schema 保留；有界 JSON loader 验证 15×16 表并记录 SHA-256 | DNAKit 不内置 DiProDB 数值；默认后 60 项为 `None`，调用者须提供有权使用并可正确引用的表 |
| 外部数据库搜索/聚类 | `conditional` | DNAKit 有内部有界搜索、alignment、sketch 索引和参考库；BLAST/MMseqs2/sourmash 有被动 metadata/版本句柄 | 没有外部科学计算执行器；程序/数据库由用户合法提供 |
| DNA 基础模型 rep + k-means | `conditional` | 11 种模型 registry、checkpoint 下载/复用、专用 adapter、rep 提取、L2/PCA/k-means、结果审计和进度条已接入；GROVER 真实 checkpoint smoke 通过 | 其余模型尚未在同一环境逐项完成真实数值验证，并受 checkpoint 条款、远程代码、独立源码、CUDA/JAX、显存和存储条件限制 |

## 能力边界（含 partial 与 complete 的有界定义）

| 功能 | 已完成 | 未完成 |
| --- | --- | --- |
| 结构化结果表 I/O（已 complete，在此仅记边界） | `TableSchema`、`read_table`、`export_table` 有界读写 CSV/TSV/JSON/Parquet；PyArrow 25.0.1 DNAKit 往返复核 | 不推断/静默拓宽 schema；通用表读取会在资源上限内物化，不是流式 DataFrame 引擎 |
| 大文件随机访问（已 complete，在此仅记边界） | FASTA/FASTQ 流式；分块；本地未压缩普通 FASTA/严格四行 FASTQ 按 ID/坐标/strand 索引且质量同步 | gzip、remote/bgzip 明确不在当前索引定义域 |
| Gap 标准化（已 complete，在此仅记边界） | 核心 known/unknown Gap、结构化 I/O、AGP codec、显式 `normalize_gaps`、有界 `sequence_from_agp` | 不进行无配置 N-run 猜测或非 AGP 2.1/模糊组件组装 |
| Duplex stability（已 complete，在此仅记边界） | 完整 Watson–Crick 互补使用 `backend="native"` 的内部 NN 复现；显式 `backend="primer3-cli"` 和用户路径 adapter 覆盖 canonical mismatch/dangling heterodimer | 修饰、用户预设 alignment 和更宽复杂离子模型不在定义域 |
| 热力学指纹（已 complete，在此仅记边界） | 固定 16 维 v2 schema；内部 duplex 与显式 Primer3 hairpin/self/heterodimer；available 位及 zero/sentinel/error 缺失策略 | 缺少显式 adapter 时不会自动探测/执行，也不会把填充值声称为后端计算 |
| 引物性质/设计（已 complete，在此仅记边界） | native GC/Tm；显式 adapter 结构属性与严格请求绑定；后端中立请求及 `Primer3CLIDesignAdapter.design()` 有界执行，候选反查模板坐标/序列/产物长度 | 不自动安装、从 `PATH` 探测或执行 Primer3；只有调用方法才执行，实验成功率/特异性仍不由此保证 |
| Pairwise alignment（已 complete，在此仅记边界） | global/local/双端 free-end semi-global；linear/三状态 affine gap；identity/coverage | IUPAC 是 literal symbol；不包含多序列比对 |
| 统一外部后端 | registry、BackendInfo、Primer3/NUPACK 执行 adapter、严格显式 Dashing adapter，以及 BLAST/MMseqs2/sourmash 被动 metadata/版本句柄 | BLAST/MMseqs2/sourmash 科学计算 adapter 及统一领域协议仍未全部覆盖 |
| 正确性验证 | 人工边界、Biopython restriction/MW/alignment/search/linkage；Primer3/NUPACK/Dashing 协议与安全契约；GROVER rep + k-means 真实 smoke | 其余 10 种 DNA 模型、Primer3、NUPACK 和 CD-HIT/MMseqs2/Dashing/BLAST/sourmash 等真实数值或大库科学差分 |
| benchmark（已 complete，在此仅记边界） | 本机时间、Python allocator 峰值、规模扩展、seed/环境/参数；construct/GC/reverse-complement 与 Biopython 1.88 同输入比较；所选公开 callable 源码行口径 | 不声称跨机器、进程 RSS或未配置外部工具排名；源码行数不是质量指标；论文统计实验不属于该 microbenchmark |

## 阻断项 { #blocked-items }

| 功能 | 状态 | 原因 | 解除条件 |
| --- | --- | --- | --- |
| NUPACK 自动安装和真实数值对照 | `blocked` | NUPACK 具有单独许可/下载条件；项目只提供被动探测与显式 adapter，不自动安装/下载，当前环境未安装 | 用户在适用许可范围内独立安装并完成真实差分；此前不得声称 NUPACK 数值验证完成 |
| NUPACK 网站演示 | `blocked` | 当前许可边界与只读网站要求 | 网站继续只展示不依赖 NUPACK 的固定结果 |
| 用户输入型 widget/Web/Galaxy | `blocked-current-scope` | 用户明确要求网站演示不可输入 | 若未来产品要求改变，另行设计隔离、资源和许可边界 |
| TestPyPI/PyPI 正式发布 | `blocked` | 项目已采用 MIT，但开发版本、依赖复核、远程授权和正式上传仍未完成 | 完成依赖/发布审计并显式授权发布 |
| GitHub push/Pages/远程 CI 结果 | `blocked-current-scope` | 用户明确要求只做本地测试 | 用户显式授权并提供远程仓库/环境后执行 |
| 论文复现实验 | `blocked-by-specification` | 没有指定论文、数据集、金标准或实验方案 | 用户指定论文/数据/口径，并完成许可与资源审查 |

## 始终不在范围内

- 启动子活性、表达量、TF 结合强度、CRISPR 编辑效率等任务型深度学习预测；
- 实验合成成功率、PCR/assembly 产率和表达效果；
- 法律意义上的序列合法性；
- 不带参考库版本的“绝对 novelty”；
- 对不可满足的多约束 split 承诺全局最优；
- 任意规模 `O(n²)` 全相似度矩阵。

条件或部分状态不等于静默降级。所有公共函数必须返回真实能力、明确缺口、参数/版本/provenance 或结构化错误。
