# 阶段 4/5 本地最终交付报告

!!! warning "2026-08-15 合规更新"
    本页保留阶段 4/5 的交付背景，但其中 2026-08-13 的测试数字属于历史快照。现行实现已移除内置 DiProDB 数值和 `primer3-py`，Primer3 改为用户显式路径的纯 CLI；当前门禁以[验证页](../validation.md)、[第三方声明](../acknowledgements.md#third-party-notices)和实际构建工件为准。

!!! note "2026-08-16 功能扩展"
    新增按显式 metadata 排除物种和染色体的 `DNASet` API；下方当前功能数量、哈希和状态统计已同步更新。

!!! note "2026-08-23 模型表征与聚类扩展"
    新增 11 种可选 DNA 基础模型 rep adapter 和 seed 固定的 k-means 聚类；当前 GROVER 真实 checkpoint smoke 通过，其余模型保持条件状态。下方数量、哈希和状态统计已同步更新。

## 1. 交付结论

DNAKit 已完成全部可在当前授权、许可和本地环境中合理实现的阶段 4/5 工作。当前版本为 `0.1.0.dev0`，只在独立 Conda Linux/WSL2 环境验证；没有推送 GitHub、部署 GitHub Pages、上传 TestPyPI/PyPI 或触发正式发布。

当前需求追踪矩阵为 184 项，SHA-256 为：

```text
e5f632295ccd1c17d576b7fd66e83a205a6a3174c4db70af5dd23559f5c740c0
```

最终逐项状态：

| 状态 | 数量 | 含义 |
| --- | ---: | --- |
| `complete` | 174 | 本地公开 API 可用且有自动测试/验证证据 |
| `partial` | 2 | 有可用且已测试子集，缺口逐项写明 |
| `conditional` | 7 | 依赖用户安装/许可/外部程序或模型运行环境；接口或当前边界明确 |
| `blocked` | 1 | 当前产品范围不允许完成 |

矩阵没有保留 `planned`，每一项 `scope_note` 均非空。`complete` 不表示无限输入、所有格式方言或实验有效性；它只表示原需求的有界、文档化实现可在本地使用和测试。

默认实现标签按当前公开方法与运行时 provenance（结果 schema 支持时）核对为：147 项 `native`、30 项 `reimplementation`、7 项 `adapter`、0 项 `novel`。具有可选实现的功能另在 `alternative_label` 中记录；没有把候选创新或被动外部句柄夸大为 `novel`/完整 adapter。

## 2. 阶段 4 已完成

### 2.1 高级 I/O 与核心操作

- 有界 GenBank 常用字段子集，严格拒绝模糊/远程/不支持 location；
- GFF3、BED3–6、AGP 2.1 codec；
- FASTA/FASTQ 流式处理、输入/序列/输出字节及 JSON 深度/节点硬上限、`iter_chunks()`、普通未压缩 FASTA/严格四行 FASTQ 的持久 ID/坐标索引；
- metadata 添加、合并、筛选、投影和 schema 校验；
- 显式 N-run 到已知长度 Gap 的有界转换，以及从 AGP 2.1 document 和组件表到已知/未知长度 Gap 序列的可审计组装；
- 任意结果表 CSV/TSV/JSON/Parquet 显式 schema 有界读写；PyArrow 25.0.1 DNAKit 往返验证；
- 环状 rotate/canonical origin/跨原点切片；record 编辑后的 feature/letter annotation 同步。

### 2.2 描述符、模式和热力学

- linguistic complexity、exact tandem-repeat union coverage；
- exact/IUPAC/regex/PWM、ORF/start/stop、启动子模式、TF PWM、限制酶、PAM/guide、CpG island、回文、倒置/串联重复、STR 和低复杂度；
- 单/双链分子量、Wallace/nearest-neighbor Tm、SantaLucia 单价盐修正、ΔH/ΔS/ΔG、stacking、局部 Tm，以及内部完整互补/显式 Primer3 CLI canonical mismatch/dangling duplex stability；
- `Primer3CLIAdapter` 的 Tm、hairpin、self-dimer、heterodimer，条件、单位、显式路径、许可和 provenance 完整保存。

内部热力学复现的范围固定为线性、无 Gap、标准 A/C/G/T、完整互补、2–60 nt 和单价 Na+；不会静默外推。canonical mismatch/dangling duplex 只有显式 `backend="primer3-cli"` 且调用方提供 `ntthal` 路径时交给外部 adapter；修饰、用户预设 alignment、Mg²+/dNTP 的内部模型不在定义域。

### 2.3 指纹、相似度、聚类和数据划分

- MinHash/FracMinHash、motif、限制酶、GC 空间、repeat、coding、固定 16 维热力学、混合、多尺度指纹；热力学结构项使用显式 adapter 或可审计缺失策略；
- 缺失值、standard/min-max/L1/L2 和低方差预处理；
- mismatch/indel 近似搜索、global/local/semi-global linear/affine pairwise alignment、identity/coverage；
- sketch similarity、持久 sketch index 和稳定 Top-k exact scan；
- circular/IUPAC/approximate deduplication；identity/edit/k-mer/fingerprint threshold graph、层次聚类和代表选择；
- 时间划分、联合 group/similarity/label 启发式划分、跨 split leakage 和质量评分。

### 2.4 综合评价与分子生物学

- validity、ambiguity、quality、complexity、uniqueness、diversity、redundancy；
- 带名称、版本、来源、日期、digest、filters 和 index parameters 的 `ReferenceLibrary`；
- reference-scoped nearest hit、novelty、memorization、reference/distribution similarity；
- 透明 synthesis-risk 和保留每项贡献的 scorecard；
- 限制酶切、末端分类、兼容性和连接；PCR/primer matching；
- 统一 primer properties、后端中立 primer design request 与显式 `Primer3CLIDesignAdapter` 有界执行；热力学结果严格绑定序列/条件/结构选项，设计候选反查模板坐标/序列/产物长度；Gibson/LCR/Golden Gate/BioBrick 序列级组装；
- CRISPR 候选/sequence-only off-target、规则型序列优化、codon optimization 和突变文库。

这些结果不预测表达、活性、结合强度、编辑效率、反应产率或实验成功率。

### 2.5 可视化

- 序列文字/高亮/Gap、线性/环状 feature map、alignment 和相似度矩阵；
- SVG 原子导出；CairoSVG/Pillow PNG/TIFF/PDF 和 600 dpi metadata；
- 完全自包含、HTML escaped、CSP 禁止网络、可搜索/展开的本地报告；
- CSV/TSV/JSON/Parquet 结果表有界读取和原子导出。

网站固定输入演示没有 DNA 序列输入框、上传、数据库查询或在线 backend。

## 3. 阶段 5 本地完成部分

- `src` 布局、Python 3.10+、类型注解、`py.typed`、pytest、Ruff、mypy 和 Linux CI 配置；
- 后端 registry、结构化后端错误、Primer3 执行 adapter、要求调用方提供可执行路径的严格 Dashing Jaccard/Top-k adapter，以及 BLAST/MMseqs2/sourmash 被动 metadata/版本句柄；
- 内容寻址 JSON 缓存、线程 batch、稳定 resume、显式 seed；
- 严格 YAML/JSON 多步骤 workflow、专用输出目录、artifact SHA-256 和 `RunManifest`；workflow 不使用通用 cache；
- 统一 CLI 覆盖标准化、验证、描述符、指纹、模式、ORF、比较、转换、去重、划分、HTML 报告和 `dnakit workflow`；
- README、安装、快速入门、主要模块 API 示例、两个真实 Notebook、固定 workflow、FAQ、引用/许可证矩阵和 MkDocs 网站；
- wheel/sdist、本地 `twine check` 与隔离安装门禁；
- TestPyPI/PyPI 手动 workflow 配置：构建审查 job 只读且无发布凭据，只有独立上传 job 具有 `id-token: write`；均含 `LICENSE`、版本和显式确认门禁，未触发；
- 本地正确性验证和 microbenchmark 机器可读 JSON。

## 4. 正确性证据

`validation/results/local_validation_report.json` 已于 2026-08-15 按合规更新后的源码重新生成，当前为 15/15 pass，0 fail、0 not-comparable、0 not-run：

- 人工 A/C/G/T、GC、k-mer、overlap search、reverse complement；
- 空序列、完整 IUPAC、200,000 nt、环状跨原点限制酶位点；
- 人工 threshold graph 聚类；
- Biopython 限制酶、分子量约定、global alignment、literal search；
- 同一 DNAKit identity distance 矩阵下，single/complete/average linkage 与 Biopython Bio.Cluster 合并距离。

主要容差：限制酶/search/人工边界为精确相等；alignment/linkage 为 `<=1e-12`；分子量按明确端基/质量表差异使用 `1.0 Da` 容差。

Primer3 不再进入正式验证器。CLI adapter 只用临时受控假可执行文件验证命令、解析、结果绑定、失败和资源边界；这不是实际 Primer3 数值对照。报告的禁止行为审计同时记录验证器没有自动发现、安装、导入或调用 Primer3，也没有安装、探测、导入或调用 NUPACK。

2026-08-13 报告中的 NUPACK 四项行为审计均为 `false`，只描述当时验证器行为。2026-08-14 新增被动探测与显式 adapter 后，受控替身契约测试已完成；当前仍没有真实 NUPACK MFE/base-pair probability 数值验证结论。

## 5. benchmark 证据

`benchmarks/results/local_benchmark_report.json` 于 `2026-08-13T11:36:49Z` 从最终稳定源码重新生成，记录 seed `20260813`、Python/DNAKit/Rich/Biopython 版本、平台、CPU 数、输入 SHA-256、任务参数、预热/重复次数、逐样本纳秒统计和 Python `tracemalloc` 峰值。

当前报告覆盖 100/1000 nt：DNAKit 运行 construct、normalize、GC、k-mer fingerprint、MinHash、subsequence search、reverse complement，Biopython 1.88 在同输入下运行三个具有直接公开对等入口的任务，共 20 个 case。报告另以统一 `inspect.getsourcelines` 口径记录所选公开 callable 的非空非注释行（DNAKit 637、Biopython 241）。这是一次本机 microbenchmark；`tracemalloc` 不是进程 RSS，源码行数不是质量指标，也不外推跨机器或未配置工具排名。

## 6. 条件、部分和阻断项

### 6.1 `partial`（2 项）

- `ENG-001`：registry、Primer3/NUPACK/Dashing 执行 adapter，以及 BLAST/MMseqs2/sourmash 被动 metadata/版本句柄完成；后三个候选外部工具无科学计算 adapter，统一执行层仍未覆盖全部候选后端；
- `ENG-014`：人工边界和 Biopython restriction/MW/alignment/search/linkage 差分通过；Primer3/NUPACK/Dashing 已通过受控协议、解析和安全测试，但没有当前真实科学差分；CD-HIT/MMseqs2/BLAST/sourmash 等大型外部搜索/聚类工具也未对照。

### 6.2 `conditional`（7 项）

- `THERMO-008..010`：依赖用户按实际许可单独安装 Primer3，并显式提供 `ntthal` 路径；DNAKit 记录 `GPL-2.0-or-later` 提示但不作最终许可结论；
- `THERMO-011`：dot-bracket/概率派生指标与显式 NUPACK adapter 已完成；真实 MFE/pairs/tube 计算仍需用户按单独许可独立安装 NUPACK，当前环境未验证；
- `SIM-013`：Dashing GPLv3 外部条件；严格显式 exact/HLL Jaccard 矩阵与稳定 Top-k adapter 已完成协议/解析/失败/资源契约测试，本地 Dashing `v1.0.2-4-g0635` 两序列 exact 文档示例 smoke 通过；仍没有科学差分，且不自动发现、安装、选择或打包第三方程序；
- `DATA-027`：11 种 DNA 基础模型 rep adapter 和 k-means 已接入；GROVER 真实 checkpoint smoke 通过，其余模型受独立依赖、源码、访问条款、远程代码和硬件条件限制。
- `EVAL-016`：Fréchet DNA distance 的统计量实现可用，真实 LucaOne 表征运行受 checkpoint、远程代码授权、依赖和硬件条件限制。

### 6.3 `blocked`（1 项）

- `ENG-017`：用户明确要求网站不可输入，因此 widget/输入型 Web/Galaxy 与当前产品范围冲突。

### 6.4 不属于当前 184 项算法实现、但阻止正式阶段 5 发布

- DNAKit 已采用 MIT `LICENSE`，但依赖兼容性复核和正式发布授权仍未完成；
- 没有远程 GitHub 仓库/发布授权；
- 用户明确要求本阶段不得推送 GitHub、部署 Pages 或上传 TestPyPI/PyPI；
- 没有指定要复现的论文、数据集、金标准、统计口径或许可，因此没有论文复现实验。

## 7. 本地环境与门禁

- Conda 环境：`/home/mapengsen/anaconda3/envs/dnakit-dev`；
- Python：3.10.20；
- `PYTHONNOUSERSITE=1`；
- 主要可选 Python 依赖：Biopython 1.88、CairoSVG 2.9.0、Pillow 12.3.0、PyArrow 25.0.1、NumPy 2.2.6、scikit-learn 1.7.2；完整 PyTorch/Transformers 模型栈在独立兼容环境按模型安装，Primer3 是用户单独安装的外部 CLI，不属于 extra；
- 外部示例 smoke：调用方显式提供的 Dashing `v1.0.2-4-g0635`；该程序不属于 Conda/Python 依赖，也不进入归档；
- 当前完整 pytest：`945 passed, 1 skipped`；Primer3 CLI 契约测试包含在完整套件中，未把替身测试表述为真实程序验证；
- 本地验证报告：`15 passed, 0 failed, 0 not run, 0 not comparable`；
- 全量 Ruff lint 通过，strict mypy 检查 170 个源码文件无问题，compileall 通过；
- MkDocs clean strict 构建通过；外部程序、数据库及真实科学差分只在用户显式提供并符合许可时运行；
- 最终 wheel/sdist 已从同一工作树重建；两包 `twine check`、归档内容审计、`src/dnakit` 逐字节一致性，以及分别在新的临时 venv 中强制安装归档后的 API/CLI help/info smoke 均通过。临时 venv 使用 `--system-site-packages` 离线复用当前已验证的 Conda 依赖，未联网重新解析依赖。

以上测试数和归档摘要均来自当前最终稳定工作树，不从旧的阶段 1～3 报告外推。

## 8. 本目标新增或修改的持久化文件

下列清单覆盖阶段 1～5 的 277 个本地持久化交付文件。花括号表示同一目录下多个真实文件的简写，不是字面文件名。清单不包含只读原始需求 `DNAKit_完整功能与输入表.csv`、用户已有的 `DNA-Rdkit.md`/`goal.md`、第三方 `dashing_similarity/`、IDE 设置，以及可重建的 `build/`、`site/`、缓存、coverage、字节码和 egg-info；两份 `dist/` 归档及其校验文件虽可重建，但作为本地发布门禁工件单独列出。

### 8.1 工程与发布配置

```text
.gitattributes
.gitignore
.github/workflows/{ci.yml,release.yml,testpypi.yml}
CHANGELOG.md
CITATION.cff
CONTRIBUTING.md
MANIFEST.in
README.md
environment-dev.yml
mkdocs.yml
pyproject.toml
```

### 8.2 Python 包

```text
src/dnakit/{__init__.py,_version.py,batch.py,exceptions.py,py.typed}
src/dnakit/alignment/{__init__.py,pairwise.py,results.py}
src/dnakit/backends/{__init__.py,external.py,registry.py}
src/dnakit/cache/{__init__.py,store.py}
src/dnakit/cli/{__init__.py,__main__.py,app.py,workflow.py}
src/dnakit/config/{__init__.py,loader.py}
src/dnakit/core/{__init__.py,_json.py,backend_info.py,collection.py,coordinates.py,enums.py,feature.py,gap.py,issues.py,provenance.py,record.py,results.py,sequence.py}
src/dnakit/datasets/{__init__.py,_advanced_shared.py,_metadata.py,advanced_deduplicate.py,advanced_split.py,clustering.py,config.py,deduplicate.py,evaluation.py,results.py,split.py}
src/dnakit/descriptors/{__init__.py,_shared.py,basic.py,codon.py,complexity.py,entropy.py,homopolymer.py,kmer.py,results.py,window.py}
src/dnakit/evaluation/{__init__.py,_shared.py,collection.py,config.py,distribution.py,reference.py,results.py,scorecard.py,sequence.py,synthesis.py}
src/dnakit/fingerprints/{__init__.py,_shared.py,advanced.py,encoding.py,kmer.py,results.py,sketch.py}
src/dnakit/io/{__init__.py,_advanced_common.py,_formats.py,annotations.py,api.py,config.py,genbank.py,indexed.py,metadata.py,parquet.py,results.py,source.py,tables.py}
src/dnakit/molbio/{__init__.py,_shared.py,assembly.py,crispr.py,optimization.py,primer3_design.py,primers.py,restriction.py,results.py}
src/dnakit/ops/{__init__.py,_common.py,circular.py,concat.py,direction.py,edit.py,mutation.py,records.py,translation.py}
src/dnakit/patterns/{__init__.py,_shared.py,coding.py,crispr.py,motif.py,regions.py,repeats.py,restriction.py,results.py}
src/dnakit/similarity/{__init__.py,_shared.py,approximate.py,compare.py,dashing.py,distance.py,index.py,matrix.py,results.py,search.py,sketch.py,vector.py}
src/dnakit/standardize/{__init__.py,_shared.py,config.py,gaps.py,normalize.py,results.py,validate.py}
src/dnakit/thermodynamics/{__init__.py,_shared.py,backends.py,calculations.py,config.py,parameters.py,results.py}
src/dnakit/visualization/{__init__.py,_svg.py,advanced.py,config.py,export.py,matrix.py,report.py,results.py,sequence.py}
src/dnakit/workflows/{__init__.py,manifest.py,runner.py,schema.py}
```

### 8.3 测试、验证与 benchmark

```text
tests/{test_benchmarks.py,test_cli_workflows.py,test_import.py,test_static_demo.py}
tests/unit/{test_alignment_pairwise.py,test_backends_registry.py,test_batch.py,test_cache_store.py,test_config_loader.py}
tests/unit/{test_core_coordinates.py,test_core_records.py,test_core_results.py,test_core_sequence.py}
tests/unit/{test_datasets_advanced_deduplicate.py,test_datasets_advanced_split_evaluation.py,test_datasets_clustering.py,test_datasets_deduplicate.py,test_datasets_split.py}
tests/unit/{test_descriptors_basic.py,test_descriptors_complexity.py,test_descriptors_kmer_entropy.py,test_descriptors_runs_windows_codons.py}
tests/unit/{test_evaluation_collection.py,test_evaluation_distribution_risk_scorecard.py,test_evaluation_reference.py,test_evaluation_sequence.py}
tests/unit/{test_external_backend_adapters.py,test_fingerprints_advanced.py,test_fingerprints_encoding.py,test_fingerprints_kmer.py,test_fingerprints_sketch.py}
tests/unit/{test_io_annotations_advanced.py,test_io_config_results.py,test_io_genbank_advanced.py,test_io_indexed_advanced.py,test_io_metadata_advanced.py,test_io_sequence_formats.py,test_io_source.py,test_io_structured.py,test_io_tables_advanced.py}
tests/unit/{test_molbio_optimization.py,test_molbio_primers_crispr.py,test_molbio_restriction_assembly.py}
tests/unit/{test_ops_circular.py,test_ops_direction_translation.py,test_ops_edit.py,test_ops_mutation_concat.py,test_ops_records.py}
tests/unit/{test_patterns_coding_restriction_crispr.py,test_patterns_motif.py,test_patterns_regions_repeats.py}
tests/unit/{test_similarity_approximate.py,test_similarity_dashing.py,test_similarity_distance.py,test_similarity_index.py,test_similarity_search.py,test_similarity_sketch.py,test_similarity_vector_matrix.py}
tests/unit/{test_standardize_gaps.py,test_standardize_normalize.py,test_standardize_validate.py}
tests/unit/{test_thermodynamics_backends.py,test_thermodynamics_calculations.py}
tests/unit/{test_visualization_advanced.py,test_visualization_charts.py,test_visualization_export.py,test_visualization_sequence.py}
tests/unit/{test_workflow_cli.py,test_workflow_manifest.py,test_workflow_runner.py}
tests/validation/test_validation_runner.py
validation/{__init__.py,CODE_REVIEW.md,README.md,run_validation.py}
validation/fixtures/README.md
validation/results/local_validation_report.json
benchmarks/{__init__.py,README.md,benchmark_core.py}
benchmarks/results/local_benchmark_report.json
```

### 8.4 文档、教程和固定示例

```text
docs/{index.md,installation.md,quickstart.md,faq.md,acknowledgements.md,validation.md}
docs/api/index.md
docs/demo/{index.md,data/fixed_demo.json}
docs/examples/index.md
docs/planning/{README.md,01_requirements_review.md,02_architecture_design.md,03_dependency_validation_development_plan.md,04_conditional_and_unavailable_features.md,05_mvp_delivery_report.md,06_stage4_stage5_delivery_report.md,07_requirements_evidence_index.md,requirements_traceability.csv}
docs/stylesheets/extra.css
docs/tutorials/index.md
examples/{README.md,advanced_workflow.yml,fixed_demo.fasta,fixed_demo_expected.json}
notebooks/{00_skeleton_check.ipynb,01_advanced_workflow.ipynb}
```

### 8.5 本地构建工件

```text
dist/dnakit-0.1.0.dev0-py3-none-any.whl
dist/dnakit-0.1.0.dev0.tar.gz
dist/SHA256SUMS
```

这两份归档已从最终稳定工作树重建：

| 工件 | 归档成员 | `src/dnakit` 文件 |
| --- | ---: | ---: |
| wheel | 156 | 151 |
| sdist | 316 | 151 |

两包的最终 SHA-256 保存在归档外的标准两列 `dist/SHA256SUMS`，可用 `sha256sum -c` 校验；最终字节数由 `stat` 或归档审计读取。这样避免了 sdist 内报告引用自身哈希造成不可收敛的自引用。两包均通过本地 `twine check`，且打包的 151 个 `src/dnakit` 文件与工作树逐字节一致。sdist 包含 `docs/planning/requirements_traceability.csv`；原始需求 CSV、`dashing_similarity/`、`site/`、缓存和 native binary 均被排除。wheel 与 sdist 分别在新的临时 venv 中强制安装归档，并通过 API、CLI help 和 `info` smoke；这些 venv 使用 `--system-site-packages` 离线复用当前已验证的 Conda 依赖，未联网重新解析依赖。它们仍只是本地工件，不代表已上传包索引。

## 9. 发布状态

| 动作 | 状态 |
| --- | --- |
| 本地 editable install | 已完成 |
| 本地测试/文档/验证/benchmark | 已完成 |
| 本地 wheel/sdist | 已完成并复核 |
| GitHub push/PR | 未执行 |
| GitHub Actions 远程运行 | 未执行 |
| GitHub Pages | 未部署 |
| TestPyPI upload | 未执行 |
| PyPI upload | 未执行 |
| 正式版本 tag/release | 未执行 |

## 10. 后续需要所有者或外部条件的动作

1. 复核 MIT `LICENSE`、`pyproject.toml`、CITATION 和引用页的一致性；
2. 对发布依赖和可选 GPL adapter 做最终法律/分发复核；
3. 若要真实 NUPACK 数值验证，由用户先在适用许可范围内独立安装并完成差分；项目不会自动安装或下载它；
4. 若要外部大库能力，在现有被动 metadata/版本句柄上逐个实现 BLAST/MMseqs2/sourmash 科学计算执行器和数据库 digest，并为 Dashing 及后续外部 adapter 增加真实科学差分；
5. 用户明确授权后，先上传 TestPyPI 并做 clean-environment install，再考虑 PyPI/GitHub/Pages；
6. 若要论文复现，先指定论文、数据集、版本、评价分母、统计方法和许可。

在这些条件满足前，当前最准确的描述是“功能完整度高、边界可审计的本地开发快照”，不是已正式发布的软件。
