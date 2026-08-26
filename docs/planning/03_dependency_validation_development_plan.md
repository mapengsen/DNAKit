# DNAKit 依赖、正确性验证与开发计划

!!! note "计划快照"
    本页保留开发前的依赖和分阶段计划，其中“当前状态”段落已由后续实现取代。当前依赖边界见[安装](../installation.md)，逐项状态见[追踪矩阵](requirements_traceability.csv)，API/测试/文档映射见[共享证据索引](07_requirements_evidence_index.md)，最终运行结论见[阶段 4/5 报告](06_stage4_stage5_delivery_report.md)。

日期：2026-08-13  
状态：开发计划快照；现行实现和最终门禁见阶段 4/5 报告，未发布任何内容。

## 1. 本地开发环境

已创建独立 conda 环境，未修改当前 `MDT2` 环境：

| 项目 | 当前状态 |
|---|---|
| 环境名 | `dnakit-dev` |
| 路径 | `/home/mapengsen/anaconda3/envs/dnakit-dev` |
| Python | `3.10.20` |
| pip | `26.1.2` |
| 已安装内容 | DNAKit 可编辑安装及 `dev`、`docs`、`io`、`validation`、`viz` extras；Primer3 不属于 Python extra |
| 用户 site 隔离 | 已固定 `PYTHONNOUSERSITE=1`，`site.ENABLE_USER_SITE=False` |

选择 Python 3.10 是为了从最低支持版本开始验证。当前环境已安装 DNAKit `0.1.0.dev0`，并用于 pytest、Ruff、mypy、MkDocs、构建和 Twine 的本地门禁；正式发布前仍需在全新环境中重新解析和验证依赖。

常用本地命令：

```bash
conda activate dnakit-dev
python --version
```

无需切换当前 shell 时：

```bash
conda run -n dnakit-dev python --version
```

仓库现有 `environment-dev.yml` 固定 Python 3.10、以可编辑方式安装 `.[dev]`，并设置 `PYTHONNOUSERSITE=1`。文档依赖可另外安装 `.[docs]`；Linux lock 文件仍是发布候选阶段的待办。所有测试、构建和文档命令必须在该环境或由该定义新建的干净环境中运行。

每次测试会话先执行隔离断言，防止 `~/.local` 用户包污染结果：

```bash
conda run -n dnakit-dev python -c "import site, sys; assert not site.ENABLE_USER_SITE; assert site.getusersitepackages() not in sys.path"
```

## 2. 依赖分组

当前直接依赖和 extras 以 `pyproject.toml` 为准；下表把已经配置的依赖与后续候选明确分开。

### 2.1 默认运行依赖

| 依赖 | 用途 | 原则 |
|---|---|---|
| PyYAML | 严格 YAML workflow 配置 | 只解析白名单 schema |
| Typer | CLI | CLI 层专用 |
| Rich | 错误、表格和进度条 | Python API 默认静默 |
| tomli | Python 3.10 TOML fallback | 仅在 Python <3.11 安装 |

当前核心对象、标准化、基础 I/O、序列操作、描述符、指纹、相似度和数据集算法主要使用 Python 标准库。NumPy、Biopython、Pydantic 和 platformdirs 未列入默认依赖；Biopython 仅用于 `validation` extra。

### 2.2 可选 extras

| Extra | 当前依赖 | 用途 |
|---|---|---|
| `docs` | MkDocs、MkDocs Material、mkdocstrings-python、nbmake | 本地文档网站和可执行 Notebook |
| `dev` | pytest、pytest-cov、Hypothesis、Ruff、mypy、build、twine、types-PyYAML | 测试、检查和本地构建 |
| `io` | PyArrow | Parquet 表读写 |
| `validation` | Biopython | 本地差分验证；不安装或执行 Primer3 |
| `viz` | CairoSVG、Pillow | PNG/TIFF/PDF 导出 |

当前未配置 mkdocs-jupyter，Notebook 通过 nbmake 独立执行而不嵌入文档站。alignment、cluster 和 molecular 当前使用内部实现，不需要额外运行依赖；基础 SVG 可视化仍只使用标准库。项目没有会自动安装受限后端的 `all` extra。

### 2.3 外部程序和数据库

以下均由用户自行安装并通过路径或环境发现，DNAKit 不自动下载：

- Primer3 的 `primer3_core`、`oligotm`、`ntthal` 和可选 thermodynamic parameter 目录；必须由用户单独安装并提供显式路径。
- BLAST+、MMseqs2、VSEARCH/CD-HIT。
- Dashing、sourmash CLI。
- NUPACK 或将来选定的合规二级结构后端；当前只提供被动探测和显式 adapter，不自动安装或下载。
- samtools/htslib、TRF、RepeatMasker、MEME/FIMO。
- REBASE、JASPAR、NCBI 或用户自定义参考库。

## 3. 许可证与发布约束

| 组件 | 已确认情况 | DNAKit 决定 |
|---|---|---|
| Biopython | 宽松 Biopython License | 当前只用于 `validation` extra；引用并记录版本 |
| BLAST+ | NCBI 说明为 public domain | 当前只有被动路径/显式有界版本句柄，无搜索执行器；数据库仍逐库记录许可和版本 |
| MMseqs2 | `GPL-3.0-or-later` | 当前只有被动路径/显式有界版本句柄，无搜索/聚类执行器；不捆绑二进制 |
| Dashing | 当前本地副本为 GPLv3；运行时记录 `GPL-3.0-only` | `dashing_similarity/` 全部排除 wheel/sdist；已有要求调用方显式路径的有界 Jaccard/Top-k adapter，不自动发现/安装；真实科学差分未完成 |
| Primer3 CLI | 官方源码头声明 GPL-2.0-or-later；仍应核对用户实际安装版本 | 不 vendor 源码/二进制、不设 Python extra、不自动按 `PATH` 发现；用户把显式路径交给有界 CLI adapter |
| sourmash | BSD-3-Clause；当前版本要求 Python 3.11+ | 当前只有被动路径/显式有界版本句柄，无 sketch/search 执行器；不作为默认依赖 |
| NUPACK 4 | 单独许可/订阅和下载条件；使用者须按当前条款复核 | adapter 标记 `conditional`；不自动安装/下载、不打包、不作为 Web 后端，当前环境无真实差分 |
| REBASE | 可访问不等于允许随包再分发 | 不捆绑数据库快照；使用 Biopython 随版本数据或用户导入数据 |
| JASPAR | CC BY 4.0 | 可选数据 adapter，必须保存版本和署名 |

关键官方依据：

- [NUPACK 下载与许可概览](https://www.nupack.org/download/overview)要求使用者按适用类别取得订阅/许可并单独下载。DNAKit 不再分发、不自动安装/下载，也不把它作为 Web 后端；adapter 仅映射公开 API 文档并要求用户显式调用。
- [NUPACK 4.1 安装说明](https://docs.nupack.org/4.1/start/)显示其由用户接受许可后单独下载，Windows 通过 WSL 使用。
- [Primer3 许可证](https://github.com/primer3-org/primer3/blob/main/LICENSE)与各 CLI 源文件头表明其 GPL 属性；DNAKit 只记录 `GPL-2.0-or-later` 提示，最终以安装版本为准。
- [BLAST 开发者说明](https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html)说明 BLAST 软件为 public domain；[本地 BLAST+ 文档](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)适用于批量任务。
- [MMseqs2 官方仓库](https://github.com/soedinglab/mmseqs2)的许可证文件标明 GPL-3.0；DNAKit metadata 使用 `GPL-3.0-or-later` 提示，发布前仍应按实际版本复核。
- [sourmash 当前包信息](https://pypi.org/pypi/sourmash/json)说明当前支持 Python 3.11 及以上，与 DNAKit 的 3.10 下限不一致。
- [Biopython 官方网站](https://biopython.org/)说明其宽松许可证。

NUPACK 自动差分验证只有两种可接受路径：取得允许本项目用法的书面许可，或由用户在其许可范围内独立人工验收后提供不含受限软件内容的结果摘要。否则文档明确列为未自动完成。

项目自身许可证已确定为 MIT。发布前仍必须结合可选 GPL adapter 的调用方式做一次法律/许可证复核。

## 4. 正确性验证策略

所有容差集中放在将来的版本化 `tests/differential/tolerances.yml`，不得散落在测试代码中。容差必须绑定算法、参数集和后端版本。以下是发布前必须达到的初始验收值；若基准集证明不合理，只能通过带依据的 ADR 修改，不能由单个测试放宽。

| 功能 | 对照 | 验收原则 |
|---|---|---|
| reverse/complement/RC、坐标、编辑 | 人工小序列 | 离散结果完全一致；验证代数不变量 |
| IUPAC | IUPAC 表和人工枚举 | 全字符 complement、兼容性和错误位置完全一致 |
| FASTA/FASTQ/GenBank | Biopython | 语义 round-trip 一致；不要求注释空白和字节完全一致 |
| 组成、k-mer、entropy、窗口 | 人工计数和小型参考实现 | 整数完全一致；无单位比例/entropy 绝对误差 `≤1e-12` |
| 分子量 | 固定 Biopython 版本和相同末端/链模型 | 绝对误差 `≤ max(1e-6 Da, 1e-10 × expected)` |
| Primer3 Tm/hairpin/dimer adapter | 用户许可环境中的固定 Primer3 版本、参数和原始输出 | 自动门禁只验证 CLI 协议/解析；真实差分应另存版本化摘要，CLI 数值误差不超过原始输出最末位的 0.5 单位 |
| nearest-neighbor reimplementation | 相同参数集的 Primer3/Biopython | Tm `≤0.1 °C`；ΔG `≤0.01 kcal/mol`；ΔH `≤0.1 kcal/mol`；ΔS `≤0.1 cal/(mol·K)` |
| 二级结构 | 同版本 NUPACK 或合规替代后端 | 许可允许时 MFE `≤1e-6 kcal/mol`、配对概率绝对误差 `≤1e-8`；MFE 并列时比较允许的最优配对集合，不强制唯一 dot-bracket |
| 限制酶 | 同版本 Biopython Restriction/REBASE 数据 | 酶名、链方向、切点和片段完全一致 |
| edit/alignment | edlib、parasail 或 EMBOSS | 相同打分下 score 一致；并列最优路径允许不同但成本相同 |
| sketch | 固定 k、seed、scaled、canonical 的 sourmash/Dashing | 同实现签名一致；跨工具比较误差和 recall@k，不要求 hash 表达一致 |
| 搜索 | BLAST/MMseqs2/人工小库 | exact 模式逐条一致；近似模式报告 recall、identity 和 coverage 口径 |
| 聚类 | CD-HIT/VSEARCH/MMseqs2 或人工图 | 不比较任意 cluster 编号，比较成对同簇关系和代表选择规则 |
| 划分 | 人工约束数据集 | 比例、组完整性、seed 重现、跨集合泄漏和不可行报告 |
| 图形 | SVG 结构检查和图片回归 | 当前 SVG 验 XML 结构、元素、顺序、转义和确定性 hash；未来 PDF 验结构，PNG/TIFF 验尺寸、dpi 和像素容差 |

必须覆盖的边界：

- 空序列、单碱基、全 N、所有 IUPAC 字符、非法 Unicode/不可见字符和 U。
- 已知 Gap、未知 Gap、相邻 Gap、N-run 不自动变 Gap。
- 环状旋转、跨原点切片和 CompoundLocation。
- 重复 ID、缺失 metadata、FASTQ quality 长度错误。
- k 大于序列长度、空 k-mer 集、空相似度集合、零向量。
- IUPAC compatibility 非传递的反例。
- 极小 strata、巨大相似度连通分量和不可满足 split。
- 超长序列使用流式/性能测试，不在普通单元测试中制造 OOM。
- 固定 seed 下，MVP 串行与 resume 结果完全一致；高级阶段启用并行后，确定性 native 离散算法的串行、并行和 resume 结果完全一致，浮点归约与外部后端按对应容差一致并记录线程数。

## 5. 软件质量门

每个功能合入前必须满足：

1. 类型注解完整，公共 API 有 docstring 和最小示例。
2. Ruff format/check、mypy 和 pytest 通过。
3. 至少一个正常测试、一个边界测试和一个错误提示测试。
4. 已知算法有引用；adapter 有后端探测和缺失测试。
5. 随机功能有固定 seed；启用并行的功能要求确定性 native 结果串并行完全一致，浮点/外部结果按冻结容差一致。
6. 结果包含 method、parameters、provenance 和 issues。
7. 大数据 API 有复杂度说明、流式路径和资源上限。
8. 功能矩阵状态、限制文档和 CHANGELOG 同步更新。
9. 所有教程 Notebook 使用 `python -m pytest --nbmake` 或等效方式从干净内核完整执行；任何 cell 失败都阻断本地发布候选。

## 6. 分阶段开发计划

### 阶段 A：需求冻结与骨架

当前状态：**在当前有界交付中完成**。需求矩阵、共享证据索引、核心对象、配置/结果基础类型、异常、provenance、`pyproject.toml`、src 布局、Linux CI、backend registry 和 workflow 均已落地；关键决定保存在架构与边界文档中，未另建独立 ADR 目录。

交付：

- 补全现有 184 项初始 traceability matrix 的 API、CLI、依赖、引用、测试和验收链接。
- 关键口径 ADR：Gap、坐标、identity/coverage、IUPAC、窗口、指纹 schema、seed。
- `pyproject.toml`、src 布局、异常、配置、provenance 和 backend protocol。
- Linux Python 3.10 的最小 CI 配置文件；当前只在本地检查 YAML，不推送 GitHub。

验收门：每项需求有唯一 ID、计划版本、实现类型、API/CLI、依赖、测试和限制。

### 阶段 B：MVP 核心与 I/O

当前状态：**在文档化定义域内完成**。核心对象、标准化/验证、记录与表 I/O、流式/chunk、FASTA/FASTQ 索引、metadata、GenBank 常用子集、GFF3/BED/AGP 以及基础/环状/record 同步操作均有公开 API 和自动测试；不支持的格式方言和索引来源会明确拒绝。

顺序：

1. 枚举、Gap、Location、DNASequence、DNAFeature、DNARecord、DNASet/RecordSource。
2. standardize、validate、audit 和结构化 issue。
3. FASTA/FASTQ、CSV/TSV/JSON、gzip 与流式读取。
4. 基础操作和坐标映射。

当前所选范围的验收门：人工小样、IUPAC、空序列、未知 Gap、FASTQ 和已支持格式的 round-trip 测试通过。环状操作属于后续高级需求，不能由当前阶段完成状态推断为可用。

### 阶段 C：MVP 分析、数据集与展示

当前状态：**在文档化定义域内完成**。描述符、指纹、相似度、去重/划分、基础 SVG/可选栅格导出、CLI、线程批处理/resume、workflow/manifest、README、Notebook 和本地 MkDocs 均已完成；每项边界以追踪矩阵为准。

顺序：

1. 基础描述符和窗口 kernel。
2. k-mer、canonical k-mer、整数/one-hot 和基础指纹。
3. exact/Hamming/edit/k-mer/fingerprint similarity。
4. exact/RC 去重、随机/分层/group/basic similarity split。
5. 基础图、CLI、批处理、进度、manifest、README、Notebook 和本地 MkDocs。

验收门：Python API 与 CLI 语义一致；MVP 串行批量保持行序；实现 resume 后，固定 seed 下串行与 resume 必须一致；其他结果按冻结容差一致；本地文档无断链。

### 阶段 D：高级模块

当前状态：**全部可行本地项已实现**；外部后端、开放式完整范围和许可缺口按追踪矩阵保持 `partial`、`conditional` 或 `blocked`。

- Primer3 adapter 与有边界的 nearest-neighbor reimplementation。
- 并行 BatchRunner、多进程/多线程执行和串并行一致性验证。
- motif/PWM/ORF、限制酶、repeat 和 CRISPR 候选。
- alignment、search、sketch、近似去重和聚类。
- motif/repeat/coding/thermo 指纹、混合指纹和多尺度指纹，以及 schema、消融和基线对比。
- 多约束 split、leakage、novelty、memorization 和 reference tracking。
- synthesis risk、分子生物学模拟和高级报告。
- NUPACK 仍受单独许可和安装条件控制；即使 adapter 已存在，也不会因其他高级模块完成而自动启用。

验收门：每个 adapter 有固定版本差分测试；每个综合指标公开公式；candidate novel 完成消融和成熟基线对比。

### 阶段 E：本地发布候选

当前状态：**本地候选门禁已完成，正式发布仍阻断**。wheel/sdist、Twine、归档审计、隔离安装、文档、Notebook、验证和 benchmark 已本地执行；项目已采用 MIT 许可证，但论文复现实验未指定，所有远程发布均未执行。

- Linux 测试 Python 3.10 到当时最新稳定小版本；Windows/macOS 不列入认证。
- 构建 wheel 和 sdist，执行 `twine check`。
- 建立两个全新 conda 环境：一个只从 wheel 安装，另一个只从 sdist 构建并安装；两者分别运行相同 smoke tests。
- 检查归档清单，确保没有 Dashing、NUPACK、数据库、测试大文件或本地缓存。
- 本地构建 MkDocs 和固定输入的只读演示。
- 从单一版本源生成包版本、`dnakit --version`、构建元数据和文档版本，并验证语义化版本一致。
- 完成根 README、安装、快速入门、每个主要模块 API 示例、可执行 Notebook 和完整工作流。
- 核对 LICENSE、CITATION.cff、CHANGELOG、CONTRIBUTING、引用/许可证清单和已完成/部分/不可用功能文档的一致性。
- 输出性能 benchmark 和论文复现实验的固定数据、配置、seed、命令与结果清单。
- 在本地审查 GitHub Actions 的 Linux 测试、构建、TestPyPI、正式 PyPI 和 Pages workflow；发布 job 保持手动审批/环境保护，当前不触发。

### 阶段 F：外部发布门（当前不执行）

1. 用户明确批准后才上传 TestPyPI。验证环境先从正式 PyPI/conda 安装锁定依赖，再以 `pip install --no-deps --index-url https://test.pypi.org/simple dnakit==VERSION` 仅安装 TestPyPI 上的 DNAKit，避免把 TestPyPI 当成不完整的依赖索引；随后运行 smoke tests。
2. TestPyPI 验证报告审查通过后，才可请求 GitHub 远程仓库、Actions/Pages 部署和正式 PyPI 的下一次独立批准。
3. 正式发布使用受保护环境和可信发布，版本 tag、wheel/sdist hash、CITATION、CHANGELOG 与文档版本必须一致。

当前阶段不会执行 `twine upload`、`git push`、创建远程仓库或部署网站。

## 7. 当前状态与下一开发入口

当前已完成全部可在本地、当前授权和许可范围内合理实现的阶段 A～E 工作；现行精确计数以阶段 4/5 报告和追踪矩阵为准。公开 API、测试和文档证据见共享证据索引，数值结果与容差见验证页，最终门禁与文件清单见阶段 4/5 报告。

后续入口不再是补齐旧 MVP 骨架，而是由所有者处理依赖法律复核与发布授权，并按矩阵缺口决定是否新增 BLAST/MMseqs2/sourmash 科学计算 adapter、Dashing/NUPACK 真实差分、更多 benchmark 对照或指定论文复现。任何真实 NUPACK 执行/差分仍须先满足许可和独立安装条件；项目只保留被动探测与显式 adapter。

补充发布风险：截至 2026-08-13，`https://pypi.org/project/dnakit/` 返回 404，名称看起来尚未被占用，但包名可用性会变化，真正 TestPyPI/PyPI 发布前必须再次实时核对。
