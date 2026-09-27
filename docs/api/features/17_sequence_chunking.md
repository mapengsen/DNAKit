# 序列切分

坐标统一为 0-based 半开区间；切分结果逐条惰性返回。

## 1) 固定长度、不重叠 <span id="2"></span> {#fixed-non-overlapping-chunks}

- **作用：** 从长序列生成互不重叠的固定长度片段。
- **计算方法：** 从起点每次前进 length；默认丢弃末尾不足长度的片段，可用 include_partial=True 保留。
- **API：** `dnakit.iter_sequence_chunks(value, config=ChunkingConfig(...))`
- **输入：** DNASequence 或 DNARecord，以及下例中的长度、步长或采样配置。
- **输出：** 返回 `SequenceChunk` 迭代器，每项包含片段序列和起止坐标。
- **示例代码：**

```python
from dnakit import DNASequence, ChunkingConfig, iter_sequence_chunks

sequence = DNASequence("AACCGGTT")
config = ChunkingConfig(length=4)
chunks = iter_sequence_chunks(sequence, config=config)
print([(c.source_start, c.source_end, c.sequence.symbols) for c in chunks])
```

- **示例结果：**

```text
[(0, 4, 'AACC'), (4, 8, 'GGTT')]
```

## 2) 固定长度、重叠滑窗 <span id="3"></span> {#overlapping-sliding-chunks}

- **作用：** 按固定步长生成可重叠窗口。
- **计算方法：** 每次前进 step；要求 1 ≤ step ≤ length。
- **API：** `dnakit.iter_sequence_chunks(value, config=ChunkingConfig(...))`
- **输入：** DNASequence 或 DNARecord，以及下例中的长度、步长或采样配置。
- **输出：** 返回 `SequenceChunk` 迭代器，每项包含重叠窗口的序列与坐标。
- **示例代码：**

```python
from dnakit import DNASequence, ChunkingConfig, iter_sequence_chunks

sequence = DNASequence("AACCGGTT")
config = ChunkingConfig(strategy="sliding", length=4, step=2)
chunks = iter_sequence_chunks(sequence, config=config)
print([(c.source_start, c.source_end, c.sequence.symbols) for c in chunks])
```

- **示例结果：**

```text
[(0, 4, 'AACC'), (2, 6, 'CCGG'), (4, 8, 'GGTT')]
```

## 3) 随机区间 <span id="4"></span> {#random-interval-chunks}

- **作用：** 抽取可复现的随机片段。
- **计算方法：** 用 seed 固定随机源，在指定长度范围内抽样 num_samples 次；结果可重复。
- **API：** `dnakit.iter_sequence_chunks(value, config=ChunkingConfig(...))`
- **输入：** DNASequence 或 DNARecord，以及下例中的长度、步长或采样配置。
- **输出：** 返回 `SequenceChunk` 迭代器，每项包含随机抽取片段的序列与坐标。
- **示例代码：**

```python
from dnakit import DNASequence, ChunkingConfig, iter_sequence_chunks

sequence = DNASequence("AACCGGTT")
config = ChunkingConfig(strategy="random", min_length=2, max_length=4, num_samples=3, seed=19)
chunks = iter_sequence_chunks(sequence, config=config)
print([(c.source_start, c.source_end, c.sequence.symbols) for c in chunks])
```

- **示例结果：**

```text
[(0, 4, 'AACC'), (0, 4, 'AACC'), (1, 5, 'ACCG')]
```

## 4) 多尺度 <span id="5"></span> {#multiscale-chunks}

- **作用：** 以多个长度分别切分同一序列。
- **计算方法：** 按 lengths 的顺序生成每个尺度的窗口；默认每个尺度不重叠。
- **API：** `dnakit.iter_sequence_chunks(value, config=ChunkingConfig(...))`
- **输入：** DNASequence 或 DNARecord，以及下例中的长度、步长或采样配置。
- **输出：** 返回 `SequenceChunk` 迭代器，每项包含窗口长度、序列和坐标。
- **示例代码：**

```python
from dnakit import DNASequence, ChunkingConfig, iter_sequence_chunks

sequence = DNASequence("AACCGGTT")
config = ChunkingConfig(strategy="multiscale", lengths=(2, 4))
chunks = iter_sequence_chunks(sequence, config=config)
print([(c.source_start, c.source_end, c.sequence.symbols) for c in chunks])
```

- **示例结果：**

```text
[(0, 2, 'AA'), (2, 4, 'CC'), (4, 6, 'GG'), (6, 8, 'TT'), (0, 4, 'AACC'), (4, 8, 'GGTT')]
```

## 5) 先短后长 <span id="6"></span> {#length-curriculum-chunks}

- **作用：** 创建先短后长的切分阶段配置。
- **计算方法：** 将长度和 stage_steps 组成计划，再用 to_config() 生成切分配置；不会自动启动训练。
- **API：** `dnakit.make_length_curriculum(lengths, stage_steps=...)` → `curriculum.to_config()` → `iter_sequence_chunks(...)`
- **输入：** DNASequence 或 DNARecord，以及下例中的长度、步长或采样配置。
- **输出：** 返回 `LengthCurriculum`，按阶段提供所选片段长度。
- **示例代码：**

```python
from dnakit import DNASequence, make_length_curriculum, iter_sequence_chunks

sequence = DNASequence("AACCGGTT")
curriculum = make_length_curriculum((2, 4), stage_steps=(10, 20))
config = curriculum.to_config()
chunks = iter_sequence_chunks(sequence, config=config)
print([(c.source_start, c.source_end, c.sequence.symbols) for c in chunks])
```

- **示例结果：**

```text
[(0, 2, 'AA'), (2, 4, 'CC'), (4, 6, 'GG'), (6, 8, 'TT'), (0, 4, 'AACC'), (4, 8, 'GGTT')]
```

<span id="1"></span><span id="7-fastabed"></span>[FAQ：补充说明](../../faq/tutorial_basics.md#chunking-details)

## 正反链配对分块与多格式导出 {#paired-chunk-export}

`dnakit.chunking.iter_paired_chunks(source, reverse_source=None, config=None, read_config=None, source_ids=None)` 读取正链 FASTA，逐对生成正链、反向互补链的行字典。每行包含 `id, sequence, pair_id, source_id, source_start, source_end, strand`，先正链、后反链。两条记录共用 `pair_id`，坐标均相对于正链，采用 0-based 半开区间。

提供 `reverse_source` 时，文件必须是正链 FASTA 各完整记录的反向互补，ID、顺序及长度一致；库检查这些结构条件，调用者负责提供正确的反向互补序列。正链 `[start:end)` 对应反链 `[L-end:L-start)`。省略该参数则直接对每个分块计算反向互补。`source_ids` 可筛选 FASTA ID。传入文件路径且未指定 `read_config` 时，读取上限会按较大的输入文件自动设置；传入流时可自行提供 `ReadConfig`。实现逐条处理染色体，内存需求与一对染色体大小有关，而不是与全部输出分块数有关。

`export_paired_chunks(source, target, ..., overwrite=False)` 将同一迭代流直接写到一个或多个目标文件，返回 `pair_count`、`record_count`、配置的 `length` 和 `coordinates`。后缀支持 CSV、TSV、JSON、JSONL、Parquet、FASTA；空结果也能写出。CSV 等表格保存全部字段，FASTA 保存 ID 和序列。

```python
from dnakit.chunking import ChunkingConfig, export_paired_chunks
from dnakit.io import save_json

result = export_paired_chunks(
    "normalized.fasta",
    ["output/chunks.csv", "output/chunks.fasta"],
    reverse_source="reverse_complement.fasta",
    config=ChunkingConfig(length=1000),
    overwrite=True,
)
save_json(result, "output/result.json", overwrite=True)
```

上例每对为两条 1000 bp 记录，共 2000 bp；不将两条链拼成一条序列。标准化、反向互补和分块可分别运行，分块脚本无需再导入 Biopython 或手工写 CSV／FASTA。文件路径的读取上限会自动适配 hg38；对流或特殊输入仍可显式传入 `read_config`。
