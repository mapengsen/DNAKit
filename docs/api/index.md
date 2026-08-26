# API 参考

本页列出当前源码中已实现、可导入并有测试覆盖的公开命名空间。所有坐标除格式转换处明确说明外，均为 0-based 半开区间。

如果希望按当前追踪矩阵的 186 项功能逐项查看作用、输入、示例、输出、状态和限制，请从[功能模块索引](features/index.md)进入。

## 顶层与核心对象 {#core-objects}

普通用户只需从 `dnakit` 顶层使用 `DNA`：一条或多条序列、ID、feature、metadata 和拓扑都由 `DNA(...)` 的输入及可选参数确定，返回类型不变。`DNASequence`、`DNARecord`、`DNASet` 保留为高级内部模型和旧代码兼容接口。

构造、读取、验证和序列/记录重复方法的逐项处理见[API 简化审计](../planning/08_api_simplification_audit.md)。

::: dnakit

核心结果 `MetricResult` 保存 value、unit、method、algorithm version、parameters、conditions、provenance 和 issues；领域结果可使用更专门的不可变 schema。

## I/O 与元数据

::: dnakit.io

主要边界：

- 普通读取统一使用 `read(..., mode="dna")` 并返回 `DNA`；大文件使用 `read(..., mode="stream")` 返回单次消费的 `RecordSource`；`read_one()`、`read_set()` 只保留兼容；
- `ReadConfig.max_sequence_symbols/max_input_bytes/max_json_depth/max_json_nodes` 和 `WriteConfig.max_output_bytes` 对记录 I/O 与 CSV/TSV/JSON/JSONL 内嵌 JSON 设硬上限；
- FASTA/FASTQ/CSV/TSV/JSON/JSONL 和 gzip 由统一 `read`/`write` 处理；
- GenBank 是明确的常用字段子集，不声称完整 INSDC 兼容；
- GFF3、BED3–6 和 AGP 2.1 使用独立严格 codec；
- `FastaIndex` 和 `FastqIndex` 只支持未压缩本地文件，并用源大小、mtime 和 SHA-256 检测陈旧索引；FASTQ 限定严格四行记录；
- `export_table()`/`read_table()` 以显式 `TableSchema` 有界读写 CSV/TSV/JSON/Parquet；读写分别限制 row/column/cell/input/output/decoded 字节，Parquet 依赖 `io` extra，并只接受 `ParquetCompression` 白名单中的压缩方式。

CSV/TSV 默认用字面量 `\N` 表示 null，空字符串仍是字符串；可通过成对设置 `null_value`/`missing_values` 改变该约定。

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit.io import TableSchema, export_table, iter_chunks, read_table

chunks = tuple(iter_chunks(range(5), chunk_size=2))
assert chunks == ((0, 1), (2, 3), (4,))

schema = TableSchema(
    ("id", "score"),
    column_types={"id": "string", "score": "number"},
    nullable=(),
)
with TemporaryDirectory() as directory:
    path = Path(directory) / "scores.json"
    export_table([{"id": "a", "score": 0.75}], path, format="json", schema=schema)
    loaded = read_table(path, format="json", schema=schema)
    assert loaded.rows[0]["score"] == 0.75
```

## 参考基因组下载

::: dnakit.references

该模块通过 NCBI Datasets v2 解析组装 accession 或物种名，流式下载 genomic FASTA，并保留 assembly report、MD5 校验和下载 provenance。完整下载流程见[下载](features/15_download.md)。

## 公共数据库查询

::: dnakit.search

`SearchConfig` 对超时、响应字节和记录数设硬上限。结果统一为不可变 `QueryResult`，包含脱敏 URL、原始提供方字段、分页信息和 provenance。Ensembl 与 UCSC 区域入口统一接收 0-based 半开坐标；NCBI BLAST 默认异步提交，不会隐式轮询。

完整查询示例和 29 项能力状态见[数据查询](features/22_database_query.md)。

## 公共数据下载

::: dnakit.download

新增下载先在目标文件系统内暂存，计算 MD5/SHA-256 并写 manifest，再事务安装；多文件失败会回滚。ENA/ENCODE 有提供方 MD5 时进行验证。`build_index()` 只执行显式可执行路径，`metadata()` 可把查询结果导出为 JSON/JSONL/CSV/TSV/XML。完整下载示例和 24 项能力状态见[下载](features/15_download.md)。

## 标准化与验证 {#standardization}

::: dnakit.standardize

`normalize()` 默认保留 IUPAC 模糊碱基，删除 `U` 和其他非 DNA 字符，并返回原始输入快照、各步骤和每个修改的原始位置。`normalize_gaps()` 只按 `GapNormalizationConfig` 将大写连续 N-run 转为已知长度 `Gap`；短 N、显式 Gap 和环状原点边界不会被猜测合并。`sequence_from_agp()` 从严格 AGP 2.1 document/entries 和调用方组件表组装序列；拒绝缺组件、不连续、越界、gapped/circular 组件和不支持方向。

`validate(value, config=...)` 是唯一推荐的合法性入口：单记录 `DNA` 返回单条报告，多记录 `DNA` 返回集合报告；集合只要有一条记录或一个集合规则失败，`is_valid` 就为 `False`。旧的 `validate_set()` 只保留兼容。

```python
from dnakit import DNA
from dnakit.standardize import GapNormalizationConfig, normalize_gaps

result = normalize_gaps(
    DNA("ACNNNNTG", alphabet="iupac"),
    config=GapNormalizationConfig(min_run_length=4),
)
assert result.sequence.is_gapped
assert result.changes[0].symbol_interval.start == 2
```

## 序列操作 {#sequence-operations}

::: dnakit.ops

普通用户始终调用无后缀名称：`insert()`、`delete()`、`substitute()`、`mask()`、`trim()`、`reverse_complement()`、`rotate()` 和 `canonical_origin()`。传入 `DNA` 时会同步 feature/逐碱基注释并返回新的 `DNA`。原有 `*_record()` 名称只作为需要详细变更审计的高级接口保留。

```python
from dnakit import DNA
from dnakit.ops import reverse_complement, rotate

linear = DNA("AACG")
circular = DNA("AACG", topology="circular")
assert reverse_complement(linear).symbols == "CGTT"
assert rotate(circular, 2).symbols == "CGAA"
```

## OPS-010 序列切分 {#sequence-chunking}

::: dnakit.chunking

`ChunkingConfig()` 默认使用固定长度 1024 bp、不重叠、`train` 标签。可通过
`strategy` 选择 `fixed`、`sliding`、`random`、`multiscale` 或 `curriculum`：

```python
from dnakit import ChunkingConfig, iter_fasta_chunks

config = ChunkingConfig(strategy="sliding", length=1024, step=512)
for chunk in iter_fasta_chunks("input.fa", config=config):
    record = chunk.to_record()
    print(record.id, chunk.source_start, chunk.source_end, chunk.split)
```

传入 `bed="regions.bed"` 后，BED 使用 0-based 半开区间；第 4 列作为
`train`、`valid`、`test` 等 split 标签。未传 BED 时，FASTA 中的每条记录全部按
`train` 处理。`LengthCurriculum((1024, 4096, 16384))` 只描述训练阶段，调用
`to_config()` 后即可生成对应的多阶段切分配置。

## 描述符 {#descriptors}

::: dnakit.descriptors

描述符覆盖长度、组成、GC/AT、skew、CpG、k-mer、Shannon entropy、homopolymer、窗口、密码子、linguistic complexity 和 exact tandem-repeat union coverage。`all_descriptors()` 按 `descriptor_schema_v1` 一次返回固定 240 字段，并为每个 `None` 值记录不可计算原因；DNAKit 不内置 DiProDB 数值，后 60 项只有在调用方显式加载有权使用的 15×16 表时才计算。[完整字段、公式、单位和来源见专页](features/05_all_descriptors.md)。

## 模式与注释 {#patterns}

::: dnakit.patterns

模式模块提供 exact/IUPAC/regex/PWM、六阅读框 ORF、start/stop、固定启动子模式、调用方提供的 TF PWM、限制酶、PAM/guide、CpG island、反向互补回文、倒置/串联重复、STR 和低复杂度区域。它只返回序列模式，不预测活性、结合强度或编辑效率。

## 热力学 {#thermodynamics}

::: dnakit.thermodynamics

内部模型复现版本化 SantaLucia 1998 参数，适用范围为线性、无 Gap、标准 A/C/G/T、完整互补、2–60 nt 和 Na⁺+K⁺ 总单价盐。`ThermodynamicConditions` 会随结果保存。`duplex_stability()` 默认的 `backend="native"` 使用该内部完整互补模型；用户提供含显式 `ntthal_path` 的 adapter 并选择 `backend="primer3-cli"` 时，可处理 canonical mismatch/dangling heterodimer。修饰和用户预设 alignment 不支持。

`optical_properties()` 在原有单链 ε260 上增加双链平均/显式 hypochromicity、每 OD260 的 nmol/质量及显式修饰修正；`concentration_from_a260()` 和 `convert_oligo_quantity()` 完成 Beer–Lambert 与单位安全的浓度/物质量/质量换算。`binding_equilibrium()`、`theoretical_melting_curve()`、`terminal_stability()` 和 `cosolvent_tm_correction()` 提供理想两态平衡、理论曲线、末端窗口及显式经验修正。

`Primer3CLIAdapter` 只执行用户明确给出的 `oligotm`/`ntthal` 路径。它不导入 Python binding、不从 `PATH` 搜索，也不自动安装或下载；统一 API 会核对返回结果的序列、完整条件、max-loop 和结构选项：

```python
from dnakit import DNASequence
from dnakit.thermodynamics import Primer3CLIAdapter, duplex_stability

adapter = Primer3CLIAdapter(
    ntthal_path="/opt/primer3/src/ntthal",
    thermodynamic_parameters_path="/opt/primer3/src/primer3_config",
)
result = adapter.self_dimer(DNASequence("GCGTTTTTCGC"))
print(result.delta_g_kcal_per_mol)
mismatch = duplex_stability(
    DNASequence("GTGCAT"),
    DNASequence("ATGCAA"),
    backend="primer3-cli",
    adapter=adapter,
)
assert not mismatch.fully_complementary
assert mismatch.provenance.implementation.label == "adapter"
```

## 二级结构

::: dnakit.secondary_structure

`analyze_dot_bracket()` 和概率派生指标不依赖外部后端。`probe_nupack()` 只被动定位；`NupackAdapter` 只在用户已获得适用许可并独立安装后由调用方显式执行，项目不自动安装、下载或静默调用。当前项目环境没有 NUPACK，因此没有真实 NUPACK 数值差分结论。

## 三维结构

::: dnakit.structure3d

PDB 原生坐标分析只报告从显式坐标可计算的几何；标准局部 12 参数通过 `read_3dna_bp_step()` 读取 3DNA 输出。沟槽、电荷、力学模量等缺少必要模型时保留为条件功能，不从普通序列猜值。

## 序列表征、指纹与预处理 {#fingerprints}

::: dnakit.fingerprints

提供整数/one-hot 编码、普通或 Canonical exact k-mer 特征、MinHash/FracMinHash k-mer Sketch，以及 motif、限制酶、GC 空间、repeat、coding、内部热力学、混合和多尺度指纹。`FeaturePreprocessor` 只用训练数据拟合，支持缺失值、standard/min-max/L1/L2 和低方差过滤。

热力学指纹使用固定 16 维 v2 schema：内部 Tm/ΔH/ΔS/ΔG 加 hairpin/self-dimer/heterodimer 的 available/found/Tm/ΔG。结构项只在调用方显式传入 adapter 时执行；否则按 `zero`、`sentinel` 或 `error` 缺失策略处理，available 位防止填充值被误读为真实计算。

## DNA 基础模型表征 {#model-representations}

::: dnakit.representations

`extract_representations()` 按所选模型为每条 `DNARecord` 返回一个只读 float32
rep。默认使用 LucaOne；缺少 checkpoint 时下载到运行目录的
`ckpt/lucaone-gene-step36-8m/`，已有完整 checkpoint 时复用。LucaOne 需要加载
checkpoint 自带代码，因此标准后端仍要求显式设置 `allow_remote_code=True`。
全部 11 种模型、checkpoint、依赖和远程代码边界见
[神经网络表征](features/08_fingerprints.md#neural-representations)；rep → k-means 见
[神经网络聚类](features/10_clustering.md#data-027-neural-clustering)。

## 相似度、搜索和比对 {#similarity-alignment}

::: dnakit.similarity

::: dnakit.alignment

`approximate_search()` 支持有界 mismatch/indel；`align_pairwise()` 支持 global、local 和双端 free-end `semi_global`，可选线性或三状态 affine gap，并返回 identity/coverage。IUPAC 当作 literal symbol 比较；DP 受 `max_cells` 约束。MinHash 索引是对内存 sketch 的确定性 exact scan，不等于 ANN 数据库。

`DashingAdapter` 是 opt-in 外部科学计算 adapter。它只接受调用方显式提供的可执行文件路径，不从 `PATH` 或项目第三方目录自动选择。`matrix()` 可执行 Dashing exact k-mer set 或 HLL-sketch Jaccard；`top_k()` 在同一次已验证矩阵上按分数和原索引稳定排序。内存序列必须线性、无 Gap 且长度不小于 k，也可传入显式 FASTA/FASTQ 路径。固定 `dist` 命令与 flag 白名单不经 shell，并限制项目数、输入/输出/捕获字节、sketch 内存、线程和超时；结果记录输入/原始输出 SHA-256、后端版本、许可和 provenance。当前没有真实 Dashing 科学差分，因此状态仍为 `conditional`。

```python
import os
from pathlib import Path

from dnakit import DNASequence
from dnakit.similarity import DashingAdapter

# 应用程序显式提供路径；DNAKit 不会自动寻找或安装 Dashing。
configured_path = os.environ.get("DNAKIT_DASHING_EXECUTABLE")
if configured_path is not None:
    dashing = DashingAdapter(Path(configured_path))
    matrix = dashing.matrix(
        (DNASequence("AACCGG"), DNASequence("AACCTT")),
        k=2,
        mode="exact",
    )
    assert matrix.provenance.backend is not None
    assert matrix.provenance.backend.name == "dashing"
```

```python
from dnakit import DNASequence
from dnakit.alignment import AlignmentConfig, align_pairwise

alignment = align_pairwise(
    DNASequence("ACGT"),
    DNASequence("ACCT"),
    config=AlignmentConfig(mode="global", gap_score=-1.0),
)
assert alignment.identity == 0.75
assert alignment.query_coverage == 1.0
```

## 数据集整理 {#datasets}

::: dnakit.datasets

功能包括 exact/reverse-complement/circular/IUPAC/近似去重、identity/edit/k-mer/fingerprint 阈值聚类、层次聚类、代表序列、随机/稳定哈希/分层/group/相似度/时间划分、联合约束启发式划分、跨集合 leakage 和划分质量。

`joint_split()` 不声称全局最优；不可满足时按 `infeasible_policy` 报错或返回带松弛标记的审计结果。

## 综合评价 {#evaluation}

::: dnakit.evaluation

reference-based 方法必须先调用 `create_reference_library()` 绑定名称、版本、来源、日期、filter、索引参数和内容 digest。novelty 定义为相对于该库的 `1 - nearest_similarity`；memorization 是 exact 或显式阈值近似复制。

`evaluate_synthesis_risk()` 的输出是透明规则和命中位置，不是供应商接单规则、结构预测或实验成功概率。`evaluate_scorecard()` 保留每一分项、归一化方向、权重、缺失策略和贡献。

## 分子生物学模拟 {#molecular-biology}

::: dnakit.molbio

酶切/末端/连接、PCR/引物匹配、组装、CRISPR 和序列优化均为确定性序列级模型。Golden Gate/BioBrick 要求预先带有已验证末端的片段；PCR 不模拟产率；CRISPR 不预测效率或生物学风险；codon optimization 需要调用方提供宿主表。

`prepare_primer_design()` 生成经过验证、后端中立的设计请求；只有调用方随后对含显式 `primer3_core_path` 的 `Primer3CLIDesignAdapter` 调用 `design()` 才执行外部程序。adapter 以 Boulder-IO 白名单映射参数，限制模板、候选数、结果键和文本长度，并返回候选坐标、Tm、GC、penalty、警告、许可和 provenance。每个候选的坐标、左右引物序列和 product size 会反查原模板，不一致结果按后端输出错误拒绝。

```python
from dnakit import DNASequence
from dnakit.molbio import (
    Primer3CLIDesignAdapter,
    digest_restriction,
    generate_mutation_library,
    prepare_primer_design,
)

digest = digest_restriction(DNASequence("TTTGAATTCAAA"), ("EcoRI",))
library = generate_mutation_library(
    DNASequence("ACGT"),
    {1: ("A", "G")},
    mode="single",
    seed=13,
)
assert len(digest.fragments) == 2
assert len(library.variants) == 2

request = prepare_primer_design(
    DNASequence("ACGT" * 100),
    target_start=150,
    target_end=200,
    candidate_count=2,
)
designer = Primer3CLIDesignAdapter(
    "/opt/primer3/src/primer3_core",
    thermodynamic_parameters_path="/opt/primer3/src/primer3_config",
)
if designer.info.available:
    designed = designer.design(request)
    assert designed.execution_performed
    assert designed.provenance.backend is not None
```

## 可视化 {#visualization}

::: dnakit.visualization

SVG 原生生成；`save_image()` 的 PNG/TIFF/PDF 需要 `viz` extra。`build_html_report()` 可将调用方显式提供的结果渲染为自包含、可筛选、可展开的只读 HTML。

## 批处理、缓存、后端和工作流 {#engineering}

::: dnakit.batch

::: dnakit.cache

::: dnakit.backends

::: dnakit.workflows

批处理支持串行或线程模式、稳定输入顺序、有界 in-flight、错误收集、每条记录派生 seed、进度回调和已完成 ID resume。`CacheKey` 自动纳入 key schema 和 DNAKit 版本；调用方必须把输入、参数及所用算法/后端版本显式放入 components。`JSONCache` 使用规范化内容 key、原子写入、payload 校验和 `max_entry_bytes` 读写上限；它没有 TTL 或自动淘汰。

内置 registry 有 6 个稳定 ID：`primer3-cli`、`nupack`、`blastn`、`mmseqs2`、`sourmash`、`dashing`。默认 `primer3-cli` 注册项不含路径，因此保持 unavailable；实际科学调用须直接构造显式路径 adapter。BLAST/MMseqs2/sourmash 的 `ExternalCLIAdapter` 只做被动路径定位，并在用户显式调用时执行有超时/输出上限的版本命令；它们不接收序列，也不是搜索、聚类或 sketch 执行器。registry 的 Dashing metadata/版本句柄同样保持被动；领域级 `DashingAdapter` 是另一条必须显式提供可执行路径、显式调用的有界科学计算入口。

YAML/JSON pipeline 使用 `dnakit-workflow-v1` 严格 schema，只允许 `normalize/validate/descriptors/fingerprint/deduplicate/split/write/report`。输出限制在配置文件下的专用目录；manifest 记录 resolved config、seed、版本、步骤状态和 artifact SHA-256。`--resume` 只跳过通过完整性校验的 `write/report`，不是通用 cache。工作流不加载任意 Python callable，不执行 shell/网络，也不自动下载后端/数据库。

```python
from pathlib import Path
from shutil import copy2
from tempfile import TemporaryDirectory

from dnakit import DNARecord, DNASequence
from dnakit.backends import backend_registry
from dnakit.batch import BatchConfig, run_batch
from dnakit.cache import CacheKey, JSONCache
from dnakit.workflows import load_workflow, run_workflow

tiny = [DNARecord(DNASequence("ACGT"), "a"), DNARecord(DNASequence("GGGG"), "b")]
batch = run_batch(
    tiny,
    lambda record, context: {"id": record.id, "seed": context.seed},
    name="api-example",
    config=BatchConfig(seed=7, jobs=2, execution_mode="thread", max_in_flight=2),
)
assert batch.success_count == 2
assert {
    "blastn",
    "dashing",
    "mmseqs2",
    "nupack",
    "primer3-cli",
    "sourmash",
}.issubset(set(backend_registry))

with TemporaryDirectory() as directory:
    root = Path(directory)
    cache = JSONCache(root / "cache")
    key = CacheKey.from_components(
        "example",
        {
            "sequence": "ACGT",
            "k": 2,
            "algorithm_version": "example-kmer-v1",
            "backend_version": None,
        },
    )
    cache.put(key, {"value": 1})
    assert cache.get(key) == {"value": 1}

    config_path = copy2("examples/advanced_workflow.yml", root)
    copy2("examples/fixed_demo.fasta", root)
    loaded = load_workflow(config_path)
    dry_run = run_workflow(config_path, dry_run=True)
    assert loaded.spec.schema_version == "dnakit-workflow-v1"
    assert dry_run.status == "dry-run"
```

## CLI

| 命令                                | 作用                           |
| ----------------------------------- | ------------------------------ |
| `info` / `backends`             | 报告运行环境与已注册后端       |
| `normalize` / `validate`        | 标准化，以及单序列或记录集合的统一验证 |
| `describe`                        | 基础组成、GC、复杂度和重复报告 |
| `fingerprint`                     | k-mer、MinHash 或 FracMinHash  |
| `search` / `orfs` / `compare` | 模式、ORF 和序列比较           |
| `convert`                         | 流式格式转换                   |
| `deduplicate` / `split`         | 数据集整理                     |
| `report`                          | 自包含只读 HTML 报告           |
| `workflow`                        | 严格 JSON/YAML 多步骤 workflow |

```bash
dnakit --help
dnakit workflow --help
```

`python -m dnakit.cli.workflow run CONFIG` 保留为兼容/开发入口。

逐项范围以[需求追踪矩阵](../planning/requirements_traceability.csv)为准。
