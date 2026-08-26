# DNAKit 阶段 2：项目架构设计

!!! note "设计基线"
    本页是阶段 2 设计基线，包含当时的拟议 API/CLI，不保证与最终公开签名逐字一致。当前可执行入口以 [API 参考](../api/index.md)、[示例](../examples/index.md)和[最终报告](06_stage4_stage5_delivery_report.md)为准。

本页的签名围栏均为不可执行的设计记法，因此使用 `text` 标记；它们不属于教程代码示例。

状态：设计基线，尚未创建包源码。  
目标平台：Linux，Python 3.10 及以上。  
设计原则：轻量顶层 API、不可变核心对象、显式算法口径、可选外部后端、流式批处理、完整 provenance。

## 1. 目录结构

```text
dnakit/
├── pyproject.toml
├── README.md
├── LICENSE
├── CITATION.cff
├── CHANGELOG.md
├── CONTRIBUTING.md
├── environment-dev.yml
├── mkdocs.yml
├── src/
│   └── dnakit/
│       ├── __init__.py
│       ├── _version.py
│       ├── py.typed
│       ├── exceptions.py
│       ├── core/
│       │   ├── enums.py
│       │   ├── gap.py
│       │   ├── coordinates.py
│       │   ├── sequence.py
│       │   ├── feature.py
│       │   ├── record.py
│       │   ├── collection.py
│       │   ├── issues.py
│       │   ├── backend_info.py
│       │   ├── results.py
│       │   └── provenance.py
│       ├── standardize/
│       ├── io/
│       ├── ops/
│       ├── descriptors/
│       ├── patterns/
│       ├── fingerprints/
│       ├── alignment/
│       ├── similarity/
│       ├── search/
│       ├── datasets/
│       │   ├── dedup.py
│       │   ├── cluster.py
│       │   ├── split.py
│       │   └── leakage.py
│       ├── thermo/
│       ├── structure/
│       ├── molecular_biology/
│       ├── evaluation/
│       ├── visualization/
│       ├── reporting/
│       ├── backends/
│       │   ├── registry.py
│       │   └── adapters/
│       ├── cache/
│       ├── batch/
│       ├── workflows/
│       ├── config/
│       └── cli/
├── tests/
│   ├── unit/
│   ├── integration/
│   ├── differential/
│   ├── property/
│   ├── performance/
│   └── data/
├── docs/
│   ├── index.md
│   ├── installation.md
│   ├── quickstart.md
│   ├── tutorials/
│   ├── api/
│   ├── examples/
│   ├── faq.md
│   ├── acknowledgements.md
│   └── demo/
├── notebooks/
├── examples/
├── benchmarks/
└── .github/workflows/
```

当前 `dashing_similarity/` 中的第三方源码、目标文件、静态库和二进制不属于 `src/dnakit`，未来也必须从 wheel/sdist 排除。

## 2. 依赖层次

只允许高层依赖低层：

```text
L0  exceptions、enums、JSON 基础类型
 ↓
L1  core 对象、Issue、通用 Result、BackendInfo、ImplementationInfo、Provenance
 ↓
L2  standardize、I/O codec、ops、config、cache 接口、backend registry 元数据
 ↓
L3  各领域 contracts：config、result、protocol；以及 descriptors、k-mer、patterns 等 kernel
 ↓
L4  基础 fingerprints、similarity、各领域外部 backend adapters
 ↓
L5  综合 fingerprints、search、dedup、cluster、split、simulation
 ↓
L6  evaluation、batch、workflows、reporting、visualization
 ↓
L7  CLI、本地文档站和静态演示
```

强制规则：

- `core` 不导入 I/O、算法、可视化、后端或 CLI。
- `DNASet` 不导入 dedup、split、cluster 或 evaluation；这些是独立服务函数。
- 可视化只消费对象和结果，不承担指标计算。
- CLI 只解析参数并调用 workflows，不包含算法。
- 后端和可选依赖延迟加载；`import dnakit` 不检查外部程序。
- 同级模块共享逻辑时下沉为更低层 kernel，禁止相互反向导入。
- 后端 registry 通过延迟 import 或 entry point 发现 adapter，避免注册时加载重依赖。
- `BackendInfo` 是 L1 的被动值对象，不包含领域方法；每个 backend protocol 由对应领域的 L3 `contracts` 模块拥有，可依赖 core 和本领域 config/result。L2 registry 只保存后端 ID、`BackendInfo`、entry point/延迟加载信息，不导入领域 protocol。
- 具体 adapter 位于 L4 并实现 L3 protocol；领域 kernel 和 contract 不反向导入 adapter。search/sketch 的 contract/result 也在 L3，L5 只放搜索编排、索引管理和综合工作流。
- 基础指纹只依赖 descriptors/k-mer；motif、repeat、thermo 和混合指纹属于 L5 的综合指纹。

## 3. 核心数据对象

核心对象采用不可变值语义；序列转换和编辑返回新对象。实现时优先使用冻结 dataclass，并复制传入的可变容器。

### 3.1 枚举

```text
DNAAlphabet = STRICT | IUPAC
Topology = LINEAR | CIRCULAR
Strandedness = SINGLE | DOUBLE
Strand = FORWARD | REVERSE | BOTH | UNKNOWN
```

`gapped` 不作为 alphabet。Gap 是 sequence parts 的一种显式元素。alignment 中的 `-` 只存在于 `AlignmentResult`，不写入普通 `DNASequence`。

### 3.2 `Gap`

概念接口：

```text
Gap(
    length: int | None,
    kind: GapKind,
    crossable: bool | None = None,
    evidence: tuple[str, ...] = (),
    metadata: Mapping[str, JSONValue] | None = None,
)
```

- `length=None` 表示未知长度。
- `crossable` 显式保存“算法是否允许跨越该 Gap”的策略；`None` 表示由调用配置决定。
- Gap 不复制保存上下游片段；上下游由其在 `DNASequence.parts` 中的位置确定。接受“Gap + 上下游片段”的便捷输入时，由 `DNASequence.from_fragments(...)` 校验并组合，而不是在 Gap 中重复存储序列。
- N-run 只有在显式策略或 AGP 依据下才转为 Gap，并写入标准化审计。
- 未知 Gap 不参与需要确定全局长度的静默估算。

### 3.3 用户门面 `DNA`

普通用户只有一个核心构造入口：

```text
DNA(
    data: str | bytes | DNASequence | DNARecord | DNASet | Iterable,
    *,
    id: str | None = None,
    description: str | None = None,
    features: Iterable[DNAFeature | Mapping] | None = None,
    metadata: Mapping | None = None,
    alphabet: DNAAlphabet | str | None = None,
    topology: Topology | str | None = None,
    strandedness: Strandedness | str | None = None,
    name/source/version/collection_metadata = None,
    max_records: int = 100_000,
) -> DNA
```

- 字符串表示一条记录，字符串列表表示多条记录，record mapping 列表可为每条记录设置不同可选信息。
- 包含显式 `Gap` 的字符串/Gap 列表表示一条 gapped 序列；无 Gap 的 multipart 序列使用 `{"parts": [...]}` 消除歧义。
- 原始输入通过现有 `normalize()` kernel 清理，审计保存在 `normalization/normalizations`。
- 单条和多条都返回 `DNA`；整数下标和切片也返回 `DNA`。
- 内部仍使用严格的 `DNASequence -> DNARecord -> DNASet` 分层，不把含糊类型判断扩散到科学算法。
- 构造迭代器时受 `max_records` 限制；大文件不经 `DNA(path)` 猜测，使用 `read(..., mode="stream")`。

### 3.4 `DNASequence`

概念接口：

```text
DNASequence(
    parts: str | bytes | Iterable[str | Gap],
    *,
    alphabet: DNAAlphabet = DNAAlphabet.STRICT,
    topology: Topology = Topology.LINEAR,
    strandedness: Strandedness = Strandedness.SINGLE,
)
```

内部规范为 `tuple[str | Gap, ...]`，合并相邻字符串。高级构造器只接受已经规范的序列；普通用户的原始文本由 `DNA(...)` 调用顶层 `normalize(...) -> NormalizationResult`，确保保留完整审计。

建议属性与规则：

- `parts`、`alphabet`、`topology`、`strandedness`。
- `symbol_length: int` 统计所有序列符号，包括 A/C/G/T 和 IUPAC 模糊字符，但不包含 Gap 长度。
- `canonical_base_count: int` 只统计 A/C/G/T；`ambiguity_count: int` 统计其他 IUPAC 符号。
- `coordinate_span: int | None` 等于 symbol 数加已知 Gap 长度；存在未知 Gap 时为 `None`。
- `length` 是 `coordinate_span` 的明确别名。
- `__len__()` 遇到未知长度抛 `UnknownLengthError`，不猜测长度。
- `is_gapped`、`has_unknown_length`。
- `reverse()`、`complement()`、`reverse_complement()` 返回新对象。
- `__eq__` 不隐式执行 U→T、RC 或环状旋转等价。
- 空线性序列可表示并由验证策略判定；空环状序列拒绝构造。
- FASTQ 不允许 Gap；quality 数组长度必须等于 `symbol_length`，因此 N 和其他 IUPAC 字符也各自对应一个质量值。

### 3.5 坐标与 `DNAFeature`

```text
Interval(start: int, end: int)  # 永远是内部 0-based、半开区间
CompoundLocation(parts: tuple[Interval, ...])
UnresolvedLocation(reason: str, anchors: ...)
ExternalInterval(start: int, end: int, system: CoordinateSystem, strand: Strand)

DNAFeature(
    type: str,
    location: Location,
    *,
    id: str | None = None,
    strand: Strand = Strand.UNKNOWN,
    label: str | None = None,
    score: float | None = None,
    phase: int | None = None,
    qualifiers: Mapping[str, JSONValue] | None = None,
    source: str | None = None,
)
```

- 内部一律 0-based、半开区间。
- `CompoundLocation` 表示跨环状原点或多段 feature。
- `UnresolvedLocation` 表示未知 Gap 造成的坐标不可解析。
- `DNAFeature` 不保存序列副本，也不自己执行 motif/ORF 扫描。

### 3.6 `DNARecord`

```text
DNARecord(
    sequence: DNASequence,
    id: str,
    *,
    description: str = "",
    features: tuple[DNAFeature, ...] = (),
    metadata: Mapping[str, JSONValue] | None = None,
    letter_annotations: Mapping[str, tuple[int | float, ...]] | None = None,
)
```

- 物种、染色体、个体、家系、批次等进入经过 schema 验证的 metadata。
- FASTQ quality 存于 `letter_annotations["phred_quality"]`。
- 复杂编辑使用 `edit_record()`，返回新 record、`CoordinateMap` 和 issues，不在对象上原地改坐标。

### 3.7 `DNASet` 与流式输入

```text
DNASet(
    records: Iterable[DNARecord],
    *,
    name: str | None = None,
    source: str | None = None,
    version: str | None = None,
    metadata: Mapping[str, JSONValue] | None = None,
)
```

MVP 中 `DNASet` 是可重复迭代的内存集合，只提供长度、迭代、索引、选择和集合元数据。大文件使用独立 `RecordSource`/`Iterable[DNARecord]`；禁止 `DNASet(path)` 隐式把大文件全部载入内存。

高级兼容入口继续保留：

- `DNASet.from_sequences(sequences, id_factory=...)` 将 DNASequence 物化为 DNARecord；默认 ID 按稳定输入序号生成。
- `DNASet.from_records(records)` 物化记录迭代器。
- 文件普通读取使用 `read(..., mode="dna") -> DNA`；大文件使用 `read(..., mode="stream") -> RecordSource`。旧的 `read_set(...)` 继续返回 DNASet。

`RecordSource` 的生命周期契约固定为：

- 它是单次消费的 `Iterator[DNARecord]`，`iter(source) is source`；耗尽后保持耗尽，不会自动重新打开文件。
- 支持 `with read(...) as source:`；`close()` 幂等，离开上下文时必须关闭底层文件、索引或外部流。
- 对手动关闭后的对象继续取值时抛 `RecordSourceClosedError`，而不是静默重新打开。
- `collect() -> DNASet` 只物化尚未消费的记录并在结束后关闭；普通用户需要完整集合时使用 `read(..., mode="dna")`。
- codec 可以声明 seek 能力，但 seek/restart 不属于通用 `RecordSource` API；可重复读取必须重新调用 `read(...)`。

### 3.8 不可变性、hash 与 CSV 契约映射

- 所有传入序列和容器都做防御性复制；metadata/qualifiers 递归限制为 JSONValue 并转换为只读结构。
- `DNASequence`、Gap、Interval 和 Location 可哈希；其 hash 覆盖所有参与结构相等的字段。
- DNARecord、DNAFeature 和 DNASet 采用不可变接口和结构相等，但因含丰富 metadata，MVP 不承诺可哈希。
- 序列名称和描述统一放入 DNARecord。CSV 对 DNASequence 的可选名称/描述属于字段位置调整，能力未删除；迁移映射记录在需求矩阵。
- CSV 对 Gap 的上下游片段通过 DNASequence parts 表达，`crossable` 保留为 Gap 字段；避免一份序列在 Gap 和容器中重复保存。
- CSV 对 DNASet 的多种输入通过有界、规则明确的 `DNA` 用户门面统一；内部模型仍保持显式类型，旧工厂继续兼容。

## 4. 公共 Python API

顶层只暴露稳定且轻量的对象和入口：

```text
from dnakit import DNA, Gap, read, validate, write

# 高级/兼容类型仍可导入：DNAFeature, DNARecord, DNASequence, DNASet
# 兼容函数仍可导入：normalize, read_one, read_set, validate_set
```

核心入口的返回语义固定为：

```text
normalize(
    raw: str | bytes | Iterable[str | Gap] | DNASequence,
    *,
    keep_ambiguous: bool = True,
    keep_u: bool = False,
    keep_other: bool = False,
    config: NormalizationConfig | None = None,
) -> NormalizationResult

validate(
    value: DNA | DNASequence | DNARecord | DNASet | Iterable[DNARecord],
    *,
    config: ValidationConfig | DatasetValidationConfig | None = None,
) -> ValidationReport | DatasetValidationReport

read(
    source: PathLike | TextIO | BinaryIO,
    *,
    format: str | None = None,
    config: ReadConfig | None = None,
    mode: Literal["dna", "stream"] = "stream",
) -> DNA | RecordSource

read_one(...) -> DNARecord       # 兼容入口；输入不是恰好一条时明确报错
read_set(...) -> DNASet          # 兼容入口；显式物化全部记录

write(
    records: DNA | DNASequence | DNARecord | Iterable[DNASequence | DNARecord],
    target: PathLike | TextIO | BinaryIO,
    *,
    format: str,
    config: WriteConfig | None = None,
) -> WriteResult
```

内存中的 `NormalizationResult` 保存完整原始输入快照、构造出的 `DNASequence`、每一步操作、issues 和修改坐标；`DNA(...)` 保存这份结果，但不增加返回类型含糊的 `DNASequence.from_raw()`。输入已经是 `DNASequence` 时，默认生成可审计的 no-op 结果并复用该不可变对象；只有配置明确请求转换时才返回新对象。为了避免敏感序列意外落盘，manifest/cache 默认只持久化原始输入的 hash 和摘要，只有显式配置才写出完整原文。

`write()` 可直接接收 `DNA`；遇到没有 record ID 的高级 `DNASequence` 时，按稳定输入顺序生成 `sequence_1`、`sequence_2` 等 ID，并把映射写入 `WriteResult.generated_ids`。统一的 `validate()` 根据 `DNA` 中的记录数完成单条或集合级检查，包括重复 ID、元数据和 letter annotation；旧的 `validate_set()` 仅作为兼容别名保留。

坐标 API：

```text
import_location(
    external: ExternalInterval | Sequence[ExternalInterval],
    *,
    sequence_length: int | None = None,
) -> Location

export_location(
    location: Location,
    *,
    target_system: CoordinateSystem,
    sequence_length: int | None = None,
) -> tuple[ExternalInterval, ...]

reverse_strand_location(
    location: Location,
    *,
    sequence_length: int,
) -> Location
```

`Location` 永远只表示内部 0-based 半开坐标，不能同时附加任意 `source_system`。`ExternalInterval` 才保存来源/目标的 origin、端点闭合方式和 strand reference；导入时先规范化为内部 `Location`，导出时每个复合片段各返回一个 DTO。缺少转换所需长度时抛 `CoordinateError`。

领域功能使用命名空间：

```text
dnakit.ops.reverse_complement(...)
dnakit.descriptors.compute(...)
dnakit.fingerprints.kmer(...)
dnakit.similarity.compare(...)
dnakit.datasets.deduplicate(...)
dnakit.datasets.split(...)
dnakit.thermo.tm(...)
dnakit.search.nearest(...)
dnakit.visualization.plot_sequence(...)
```

MVP 领域 API 的概念签名：

```text
descriptors.compute(
    sequence: DNASequence | DNARecord,
    names: Sequence[str],
    *,
    config: DescriptorConfig | None = None,
) -> DescriptorResult

fingerprints.kmer(
    sequence: DNASequence | DNARecord,
    *,
    k: int,
    canonical: bool = False,
    mode: Literal["presence", "count", "frequency"] = "count",
) -> FingerprintResult

similarity.compare(
    left: DNASequence | FingerprintResult,
    right: DNASequence | FingerprintResult,
    *,
    method: str,
    config: SimilarityConfig | None = None,
) -> SimilarityResult

datasets.deduplicate(
    records: Iterable[DNARecord],
    *,
    equivalence: Literal["exact", "reverse_complement"],
    config: DeduplicationConfig | None = None,
) -> DeduplicationResult

datasets.split(
    records: Iterable[DNARecord],
    *,
    config: SplitConfig,
) -> SplitResult
```

设计约束：

- 单序列函数统一接受 `DNASequence | DNARecord`，内部立即解析为明确对象。
- 批量函数接受 `Iterable[DNARecord]`，不强制创建 `DNASet`。
- 核心返回结构化结果；便捷数值属性可以提供，但不能丢失方法、条件和版本。
- 公共名称一旦进入稳定版本，按语义化版本管理弃用。

## 5. 结果对象、Issue 与异常

### 5.1 通用结果

```text
MetricResult[T](
    name: str,
    value: T,
    unit: str | None,
    method: str,
    algorithm_version: str,
    parameters: Mapping[str, JSONValue],
    conditions: Mapping[str, JSONValue],
    uncertainty: Uncertainty | None,
    provenance: Provenance,
    issues: tuple[Issue, ...],
)
```

实现来源、许可证、引用和置信信息使用固定 schema：

```text
ImplementationInfo(
    label: Literal["native", "adapter", "reimplementation", "novel"],
    execution_mode: Literal["internal", "external", "hybrid"],
    origin_class: Literal[
        "dnakit", "standard", "published_algorithm", "integration", "novel"
    ],
    license_expression: str | None,
    citations: tuple[Citation, ...],
)

Uncertainty(
    confidence_interval: tuple[float, float] | None,
    standard_error: float | None,
    method: str | None,
)
```

对外四类标签的判定优先级固定：调用外部实现为 `adapter`；内部复现公开算法为 `reimplementation`；经检索、定义和实验确认的新方法为 `novel`；其余 DNAKit 框架和简单逻辑为 `native`。`integration` 只属于 `origin_class`，不是第五个公开标签。结果对象具有 `provenance` 字段时，以运行时 `provenance.implementation` 为权威；部分早期/轻量结果 schema（例如部分 descriptor/sketch 结果）尚无该字段，其标签以追踪矩阵和模块文档为权威，统一 schema 属后续兼容性改进项。

对于具有 `provenance` 字段的结果，`provenance.implementation` 和 `provenance.backend` 分别是实现分类与后端信息的运行时权威位置；结果 schema 不应复制两份可能矛盾的真值。尚无该字段的早期/轻量结果继续以追踪矩阵和模块文档分类，后续统一时需保持序列化兼容。

领域结果至少包括：

- `NormalizationResult`、`ValidationReport`、`DatasetValidationReport`、`EditResult`、`WriteResult`。
- `DescriptorResult`、`FingerprintResult`。
- `AlignmentResult`、`SimilarityResult`、`SearchResult`、`SketchResult`。
- `DeduplicationResult`、`ClusterResult`、`SplitResult`、`LeakageReport`。
- `ThermoResult`、`StructureResult`、`RestrictionMapResult`、`SimulationResult`。
- `EvaluationReport`、`BatchResult`。

结果必须可序列化为 JSON 兼容结构；大型矩阵允许引用 NPZ/Parquet artifact，而不是嵌入 JSON。

`WriteResult` 至少保存 `format`、写出记录数、字节数（流不支持统计时为 `None`）、生成的匿名 ID 映射、目标 artifact、provenance 和 issues。`DatasetValidationReport` 至少保存记录数、ID 唯一性、集合级 metadata/quality 问题和逐记录报告引用，避免把所有明细强制常驻内存。

### 5.2 结构化问题

```text
Issue(
    code: str,
    severity: INFO | WARNING | ERROR,
    message: str,
    location: Location | None,
    details: Mapping[str, JSONValue],
)
```

普通序列质量问题进入 report；违反对象不变量、配置错误或后端执行失败才抛异常。

### 5.3 异常树

```text
DNAKitError
├── SequenceError
│   ├── InvalidAlphabetError
│   ├── UnknownLengthError
│   └── UnsupportedGapOperationError
├── CoordinateError
├── FeatureError
├── DuplicateIDError
├── InputFormatError
├── RecordSourceClosedError
├── ConfigurationError
├── CacheError
├── BackendError
│   ├── BackendUnavailableError
│   ├── BackendVersionError
│   ├── BackendTimeoutError
│   └── BackendExecutionError
└── BatchExecutionError
```

所有异常应包含稳定错误码、简短原因、相关参数和可执行的修复提示。

## 6. 配置系统

核心对象使用标准库类型；Pydantic 只用于配置和外部输入验证，避免把核心模型绑定到验证框架。

配置分组：

- `NormalizationConfig`
- `ExecutionConfig`
- `CacheConfig`
- `BackendConfig`
- `ProvenanceConfig`
- 各算法的版本化 config

优先级固定为：

```text
代码默认值 < YAML/JSON 配置文件 < DNAKIT_* 环境变量 < CLI 参数
```

每次运行保存最终解析后的配置；不能只保存用户原始配置。未知字段默认报错，防止拼写错误被忽略。

## 7. 缓存设计

- Python API 默认关闭磁盘缓存；CLI workflow 可显式开启。
- 内容寻址 key 至少包含：规范化序列/Gap 表示、算法和 schema 版本、完整参数、随机 seed/RNG 版本、后端及版本、数据库 checksum、DNAKit 版本。
- 不使用 pickle；元数据用 JSON/SQLite，数值用 NPZ/Parquet，图形和日志作为独立 artifact。
- 默认不缓存失败结果；超时和暂时性错误不得污染缓存。
- 输入或参数变化时 checkpoint/resume 必须拒绝错误恢复。
- manifest 默认记录输入 hash，而不是复制敏感序列；原始序列落盘必须显式选择。

缓存契约：

```text
Cache.get(key: CacheKey) -> ArtifactRef | None
Cache.put(key: CacheKey, value: SerializableResult) -> ArtifactRef
Cache.invalidate(key: CacheKey) -> bool
Cache.clear(namespace: str | None = None) -> CacheClearReport
```

`ArtifactRef` 至少含相对路径、media type、schema version、SHA-256、字节数和创建时间；读取时先验证 hash，再反序列化。

## 8. 后端接口

不建立一个万能后端，而是定义小型 capability protocol：

- `TmBackend`
- `DimerBackend`
- `HairpinBackend`
- `StructureBackend`
- `AlignmentBackend`
- `SearchBackend`
- `SketchBackend`
- `RestrictionDatabaseBackend`

共同能力：

```text
probe() -> BackendInfo
capabilities() -> frozenset[str]

TmBackend.calculate_tm(
    sequence: DNASequence,
    config: TmConfig,
) -> ThermoResult

DimerBackend.calculate_dimer(
    left: DNASequence,
    right: DNASequence,
    config: DimerConfig,
) -> ThermoResult

HairpinBackend.calculate_hairpin(
    sequence: DNASequence,
    config: HairpinConfig,
) -> ThermoResult

StructureBackend.fold(
    strands: Sequence[DNASequence],
    config: StructureConfig,
) -> StructureResult

AlignmentBackend.align(
    query: DNASequence,
    target: DNASequence,
    config: AlignmentConfig,
) -> AlignmentResult

SearchBackend.search(
    queries: Iterable[DNARecord],
    reference: ReferenceSpec,
    config: SearchConfig,
) -> SearchResult

SketchBackend.sketch(
    records: Iterable[DNARecord],
    config: SketchConfig,
) -> SketchResult

RestrictionDatabaseBackend.find_sites(
    sequence: DNASequence,
    config: RestrictionConfig,
) -> RestrictionMapResult
```

`BackendInfo` 保存名称、版本、可执行文件路径或包位置、许可证提示和已探测能力。adapter 要求：

- 延迟加载，缺少后端不影响导入 DNAKit。
- 外部进程使用参数列表而非 shell 字符串，并有超时、退出码、stdout/stderr 和临时目录管理。
- 缺失时报告探测命令和安装文档，不自动下载安装二进制或数据库。
- 输出统一转换为 DNAKit 结果对象，同时保留原始后端字段和单位转换记录。
- NUPACK 提供被动探测与显式 adapter，状态为 `conditional`；不自动安装/下载，也不能作为 Web 后端，真实执行仍需用户满足单独许可与安装条件。

## 9. 批处理、进度与可复现性

`BatchRunner` 设计要求：

- MVP 的流式批处理先固定 `jobs=1`；多进程/多线程和串并行一致性验收随 `ENG-007` 在高级阶段启用。
- 接受输入迭代器并按 chunk 处理，默认保持原输入顺序。
- `on_error=raise|skip|record`，默认批量任务使用 `record`。
- CPU 任务使用进程，I/O/外部进程协调可使用线程；确定性 native 整数/字符串算法要求串并行完全一致，浮点归约和外部后端按冻结容差一致并记录线程数。
- 提供 `jobs`、`chunk_size`、progress callback、result sink、checkpoint/resume。
- Python API 默认不打印；CLI 使用 Rich 进度条，数据写 stdout，进度和警告写 stderr。
- 随机行为只使用显式 RNG；子 seed 由主 seed 与稳定记录 ID/输入位置确定。

正式 provenance schema：

```text
Provenance(
    dnakit_version: str,
    python_version: str,
    platform: str,
    dependency_versions: Mapping[str, str],
    implementation: ImplementationInfo,
    backend: BackendInfo | None,
    reference: ReferenceInfo | None,
)

RunManifest(
    run_id: str,
    command: tuple[str, ...],
    resolved_config: Mapping[str, JSONValue],
    seed: int | None,
    seed_derivation: str | None,
    inputs: tuple[ArtifactRef, ...],
    outputs: tuple[ArtifactRef, ...],
    provenance: Provenance,
    started_at: str,
    finished_at: str | None,
    status: str,
    issues: tuple[Issue, ...],
)
```

`RunManifest` 至少记录：

- DNAKit、Python、平台、依赖版本。
- 完整解析配置、CLI 和运行 ID。
- 主 seed 与子 seed 派生规则。
- 后端名称、路径、版本、能力和参数。
- 参考库名称、版本、日期、checksum 和筛选条件。
- 输入文件 hash、大小、记录数；输出 artifact hash。
- 开始/结束时间、状态、失败项和 cache schema。

## 10. CLI 设计

MVP 子命令：

```text
dnakit info
dnakit backends
dnakit validate
dnakit normalize
dnakit convert
dnakit describe
dnakit fingerprint
dnakit similarity
dnakit dedup
dnakit split
dnakit plot
dnakit run
```

高级阶段增加 `thermo`、`structure`、`search`、`cluster`、`evaluate` 和 `simulate`。

统一参数语义：`--input`、`--output`、format、`--jobs`、`--chunk-size`、`--seed`、`--on-error`、`--progress`、`--cache-dir`、`--resume`、`--force`。相同 flag 在所有子命令中必须保持同一含义。

## 11. 文档站与只读演示

- 使用 MkDocs Material、mkdocstrings 和可执行 Notebook 文档。
- 页面包含首页、安装、快速入门、教程、API、完整工作流、FAQ、引用和许可证。
- `docs/demo/` 只展示固定 DNA 示例、预生成 JSON/CSV/SVG/PNG 和可展开说明。
- 不提供输入框、文件上传、数据库查询或 NUPACK/BLAST 在线运行。
- 本地使用 `mkdocs serve` 预览，当前不配置实际 GitHub Pages 部署动作。
