# 184 项需求共享证据索引

本页把当前 CSV 的需求 ID 范围映射到当前公开 API、自动测试和文档边界。它与[184 项追踪矩阵](requirements_traceability.csv)共同使用：

- 矩阵逐项给出 `status`、实现类型和 `scope_note`；
- 本索引给出共享 kernel 的公开入口和测试路径；
- [验证与 benchmark](../validation.md)和[最终交付报告](06_stage4_stage5_delivery_report.md)给出运行结果、容差和尚未完成的外部对照。

路径均相对项目根目录。为避免表格过宽，`test_core_{sequence,records}.py` 这类写法使用 shell brace expansion 简写，表示两个真实文件 `test_core_sequence.py` 和 `test_core_records.py`，不是一个字面文件名；本页所有展开后路径均已用 `rg --files` 核对。同一测试文件可为多项共享 kernel 提供证据；不能由“文件存在”推断无限范围的 complete，仍必须阅读对应 `scope_note`。

## 核心、I/O、标准化和操作

| 需求 ID | 公开 API / 实现入口 | 自动测试 | 文档与边界 |
| --- | --- | --- | --- |
| `CORE-001..008` | 普通入口 `dnakit.DNA`；高级 `dnakit.core` 的 `DNASequence/DNARecord/DNASet/DNAFeature/Gap/Location/MetricResult` | `tests/unit/test_core_facade.py`、`test_core_{sequence,records,coordinates,results}.py` | [API：核心对象](../api/index.md#core-objects)、[快速入门](../quickstart.md#core-standardization)、[简化审计](08_api_simplification_audit.md) |
| `IO-001`, `IO-004` | `dnakit.io.read(mode="dna"|"stream")/write`；兼容 `read_one/read_set`；`ReadConfig/WriteConfig/RecordSource` | `tests/unit/test_core_facade.py`、`test_io_{sequence_formats,structured,source,config_results}.py` | [API：I/O](../api/index.md#io)、[安装](../installation.md) |
| `IO-002` | `read_gff3/write_gff3`、`read_bed/write_bed`、`read_agp/write_agp` | `tests/unit/test_io_annotations_advanced.py` | [API：I/O](../api/index.md#io)、[FAQ](../faq.md#annotation-format-scope) |
| `IO-003`, `VIZ-016` | `TableSchema`、`read_table`、`export_table`、`export_result` | `tests/unit/test_io_tables_advanced.py` | [API：I/O](../api/index.md#io)、[安装：Parquet](../installation.md#parquet) |
| `IO-005` | `iter_chunks`、FASTA/FASTQ 的 `build_*_index/load_*_index` 与 `FastaIndex.fetch/FastqIndex.fetch` | `tests/unit/test_io_indexed_advanced.py` | [FAQ：大型文件](../faq.md#large-files) |
| `IO-006` | `with_metadata/merge_metadata/filter_by_metadata/select_metadata/validate_metadata` | `tests/unit/test_io_metadata_advanced.py` | [API：I/O](../api/index.md#io) |
| `STD-001..006`, `STD-008..009` | `dnakit.standardize.normalize/validate`；`validate` 统一覆盖单条序列、记录和记录集合；相应 config/result 对象 | `tests/unit/test_standardize_{normalize,validate}.py` | [API：标准化](../api/index.md#standardization)、[快速入门](../quickstart.md#core-standardization) |
| `STD-007` | `normalize_gaps`、`sequence_from_agp`；`GapNormalizationConfig/AGPAssemblyConfig` | `tests/unit/test_standardize_gaps.py` | [API：标准化](../api/index.md#standardization)、[条件与范围](04_conditional_and_unavailable_features.md) |
| `OPS-001..009` | 普通用户使用无后缀编辑/反向互补/环状操作；`*_record` 保留详细审计；其余不同算法保持明确名称 | `tests/unit/test_core_facade.py`、`test_ops_{direction_translation,edit,mutation_concat,circular,records}.py`、`test_ops_evolution.py` | [API：序列操作](../api/index.md#sequence-operations)、[简化审计](08_api_simplification_audit.md) |

## 描述符、模式、热力学和指纹

| 需求 ID | 公开 API / 实现入口 | 自动测试 | 文档与边界 |
| --- | --- | --- | --- |
| `DESC-001..012` | `dnakit.descriptors` 的 basic/k-mer/entropy/complexity/homopolymer/window/codon API | `tests/unit/test_descriptors_{basic,kmer_entropy,complexity,runs_windows_codons}.py` | [API：描述符](../api/index.md#descriptors) |
| `PAT-001..013` | `dnakit.patterns` 的 motif/PWM、coding/ORF、restriction/PAM、region/repeat API | `tests/unit/test_patterns_{motif,coding_restriction_crispr,regions_repeats}.py` | [API：模式](../api/index.md#patterns) |
| `THERMO-001..007`, `THERMO-012..013` | `dnakit.thermodynamics` 的分子量、光学/浓度、NN、稳定性、平衡、熔解曲线和局部窗口 API | `tests/unit/test_thermodynamics_calculations.py`；`tests/unit/test_thermodynamics_optics_equilibrium.py`；`validation/run_validation.py` | [API：热力学](../api/index.md#thermodynamics)、[双链扩展](../api/features/19_duplex_thermodynamics.md) |
| `THERMO-008..010`, `MOLBIO-006` | `Primer3CLIAdapter`、`validate_primer3_result`、`probe_primer3`；`primer_properties(..., structure_adapter=...)`；序列/条件/选项绑定 | `tests/unit/test_thermodynamics_backends.py`；`tests/unit/test_molbio_primers_crispr.py`；`tests/validation/test_validation_runner.py` | [安装：Primer3](../installation.md#primer3)、[验证](../validation.md#primer3) |
| `THERMO-011` | `dnakit.secondary_structure` 的 dot-bracket/概率指标及显式 `NupackAdapter` | `tests/unit/test_secondary_structure.py`；`tests/validation/test_validation_runner.py` | [二级结构](../api/features/20_secondary_structure.md)、[NUPACK 审计](../validation.md#nupack) |
| `FP-001..014` | `dnakit.fingerprints` 的 encoding/k-mer/sketch/advanced/preprocessor API | `tests/unit/test_fingerprints_{encoding,kmer,sketch,advanced}.py` | [API：指纹](../api/index.md#fingerprints) |

当前 184 项功能清单之外还包括光学/浓度、双链平衡与三维结构扩展：测试分别见 `tests/unit/test_thermodynamics_optics_equilibrium.py`、`tests/unit/test_secondary_structure.py`、`tests/unit/test_structure3d.py` 和 `tests/integration/test_downloaded_dna_structures.py`，文档见本节新增功能页。

## 相似度、数据集和综合评价

| 需求 ID | 公开 API / 实现入口 | 自动测试 | 文档与边界 |
| --- | --- | --- | --- |
| `SIM-001..007`, `SIM-010..012`, `SIM-016` | `dnakit.similarity` 的 search/distance/compare/vector/sketch/matrix API | `tests/unit/test_similarity_{search,distance,approximate,sketch,vector_matrix}.py` | [API：相似度](../api/index.md#similarity-alignment) |
| `SIM-008` | `dnakit.alignment.align_pairwise`；`AlignmentConfig/AlignmentResult` | `tests/unit/test_alignment_pairwise.py`；`validation/run_validation.py` | [API：比对](../api/index.md#similarity-alignment)、[验证](../validation.md) |
| `SIM-013` | `DashingAdapter.matrix/top_k`；`DashingJaccardMatrixResult/DashingTopKResult`；registry 被动 metadata/版本句柄 | `tests/unit/test_similarity_dashing.py`；`tests/unit/test_external_backend_adapters.py` | [API：Dashing](../api/index.md#similarity-alignment)、[条件支持](04_conditional_and_unavailable_features.md) |
| `SIM-014..015` | `build_sketch_index/save_sketch_index/load_sketch_index/nearest_neighbors` | `tests/unit/test_similarity_index.py` | [API：相似度](../api/index.md#similarity-alignment) |
| `DATA-001..006` | `deduplicate/deduplicate_iupac/deduplicate_approximate` | `tests/unit/test_datasets_{deduplicate,advanced_deduplicate}.py` | [API：数据集](../api/index.md#datasets) |
| `DATA-007..011` | `cluster_sequences/hierarchical_cluster/select_representatives` | `tests/unit/test_datasets_clustering.py`；`validation/run_validation.py` | [API：数据集](../api/index.md#datasets)、[聚类对照](../validation.md) |
| `DATA-012..021` | `split(random/hash)/temporal_split`；`SplitConfig` | `tests/unit/test_datasets_split.py` | [快速入门：划分](../quickstart.md#split-reference-evaluation) |
| `DATA-022..024` | `joint_split/detect_leakage/evaluate_split_quality` | `tests/unit/test_datasets_advanced_split_evaluation.py` | [API：数据集](../api/index.md#datasets)、[FAQ：联合划分](../faq.md#joint-split) |
| `DATA-025..026` | `exclude_by_metadata/exclude_species/exclude_chromosomes` | `tests/unit/test_datasets_filter.py` | [API：数据集](../api/index.md#datasets)、[数据集整理与划分](../api/features/10_datasets.md) |
| `DATA-027` | `dnakit.representations.extract_representations`；`dnakit.datasets.neural_cluster_sequences` | `tests/unit/test_neural_clustering.py`；GROVER 本地 checkpoint smoke | [序列表征：神经网络表征](../api/features/08_fingerprints.md#neural-representations)、[神经网络聚类](../api/features/10_clustering.md#data-027-neural-clustering)、[安装](../installation.md#dna-rep-k-means) |
| `EVAL-001..007` | `evaluate_validity/ambiguity/quality/complexity/uniqueness/diversity/redundancy` | `tests/unit/test_evaluation_{sequence,collection}.py` | [API：综合评价](../api/index.md#evaluation) |
| `EVAL-008..012`, `EVAL-015` | `ReferenceLibrary`、`create_reference_library`、novelty/memorization/reference/distribution API | `tests/unit/test_evaluation_{reference,distribution_risk_scorecard}.py` | [快速入门：参考库](../quickstart.md#split-reference-evaluation)、[FAQ](../faq.md#novelty-memorization) |
| `EVAL-013..014` | `evaluate_synthesis_risk/evaluate_scorecard` | `tests/unit/test_evaluation_distribution_risk_scorecard.py` | [API：综合评价](../api/index.md#evaluation)、[只读演示](../demo/index.md) |
| `EVAL-016` | `evaluate_frechet_distance` | `tests/unit/test_evaluation_frechet.py` | [综合评价：Fréchet DNA distance](../api/features/12_evaluation.md#eval-016-frechet-dna-distance) |
| `EVAL-017..018` | `evaluate_fragment_similarity/evaluate_snn` | `tests/unit/test_evaluation_generative.py` | [综合评价：Frag](../api/features/12_evaluation.md#eval-017-frag)、[SNN](../api/features/12_evaluation.md#eval-018-snn) |

## 分子生物学、可视化和工程

| 需求 ID | 公开 API / 实现入口 | 自动测试 | 文档与边界 |
| --- | --- | --- | --- |
| `MOLBIO-001..003` | `digest_restriction/classify_restriction_end/check_end_compatibility/ligate_fragments` | `tests/unit/test_molbio_restriction_assembly.py` | [API：分子生物学](../api/index.md#molecular-biology) |
| `MOLBIO-004..007` | `simulate_pcr/match_primer/primer_properties/prepare_primer_design/Primer3CLIDesignAdapter.design`；候选模板坐标/序列/product-size 反查 | `tests/unit/test_molbio_primers_crispr.py` | [API：分子生物学](../api/index.md#molecular-biology)、[安装：Primer3](../installation.md#primer3) |
| `MOLBIO-008` | `simulate_assembly` | `tests/unit/test_molbio_restriction_assembly.py` | [API：分子生物学](../api/index.md#molecular-biology) |
| `MOLBIO-009..012` | `scan_crispr_candidates`、规则/密码子优化、`generate_mutation_library` | `tests/unit/test_molbio_{primers_crispr,optimization}.py` | [API：分子生物学](../api/index.md#molecular-biology) |
| `VIZ-001..006`、`VIZ-008..009` | `dnakit.visualization` 的 sequence/highlight/map/alignment/matrix API | `tests/unit/test_visualization_{sequence,charts,advanced}.py` | [API：可视化](../api/index.md#visualization) |
| `VIZ-014..015` | `save_svg/save_image/build_html_report/save_html_report` | `tests/unit/test_visualization_export.py`；`tests/test_cli_workflows.py` | [API：可视化](../api/index.md#visualization)、[快速入门：可视化](../quickstart.md#visualization) |
| `ENG-001..002` | `BackendRegistry/BackendInfo`；Primer3/Dashing 执行 adapter；BLAST/MMseqs2/sourmash 被动 metadata/版本句柄 | `tests/unit/test_{backends_registry,external_backend_adapters,similarity_dashing}.py` | [API：后端](../api/index.md#engineering)、[致谢与引用](../acknowledgements.md#citations-and-license) |
| `ENG-003..005` | 类型化 Python API；`dnakit` CLI；`dnakit.workflows`、`dnakit workflow` | `tests/test_{import,cli_workflows}.py`；`tests/unit/test_workflow_{runner,cli,manifest}.py` | [API](../api/index.md)、[工作流示例](../examples/index.md#yaml) |
| `ENG-006..010` | `run_batch/iter_batch`、`JSONCache/CacheKey`、config loader、固定 seed | `tests/unit/test_{batch,cache_store,config_loader}.py` 及 split/mutation 测试 | [API：批处理/缓存](../api/index.md#engineering) |
| `ENG-011..012` | `Provenance/BackendInfo/RunManifest`；结构化异常/Issue/CLI exit | `tests/unit/test_{core_results,workflow_manifest}.py`；`tests/test_cli_workflows.py` | [API](../api/index.md)、[FAQ](../faq.md) |
| `ENG-013` | pytest 测试集与 `pyproject.toml` 严格配置 | `tests/**` | [最终报告：门禁](06_stage4_stage5_delivery_report.md#7) |
| `ENG-014` | `validation.run_validation` 和机器可读报告 | `tests/validation/test_validation_runner.py` | [验证项与容差](../validation.md) |
| `ENG-015` | `benchmarks.benchmark_core`；DNAKit/Biopython 对等任务、规模/时间/`tracemalloc`/所选 callable 行数的机器可读报告 | `tests/test_benchmarks.py` | [benchmark 结果与限定](../validation.md#microbenchmark) |
| `ENG-016` | README、MkDocs、mkdocstrings、Notebook、workflow fixture | `python -m pytest --nbmake notebooks`；`mkdocs build --strict` | [教程](../tutorials/index.md)、[致谢与引用](../acknowledgements.md#citations-and-license) |
| `ENG-017` | 固定 JSON/FASTA 夹具和预生成静态页；无 DNA 输入 API | `tests/test_static_demo.py`；构建后 HTML 输入/外部字体审计 | [固定输入演示](../demo/index.md)、[阻断原因](04_conditional_and_unavailable_features.md#blocked-items) |

## 机器可读证据

- `validation/results/local_validation_report.json`：人工边界、Biopython/Bio.Cluster 和 Primer3 对照；
- `benchmarks/results/local_benchmark_report.json`：固定 seed/规模/预热/重复的本机 microbenchmark；
- `examples/fixed_demo_expected.json` 与 `docs/demo/data/fixed_demo.json`：固定演示预期值；
- workflow 每次运行的 `run-manifest.json`：resolved config、seed、版本、步骤状态和 artifact SHA-256。

NUPACK 当前没有真实成功结果工件；受控替身测试只证明 adapter 契约。`temp/dna_structures/analysis_results.json` 是 PDB/DBN/DSSR 本地分析，不是 NUPACK 输出。
