# OPS-010 序列切分

**作用：** 按固定窗口、重叠滑窗、随机区间或多尺度策略把长 DNA 切成带来源坐标的片段，并以迭代器逐条返回，用于模型输入、分区分析和大序列流式处理。

按照固定长度、滑动窗口、随机区间或多尺度策略切分 DNA 序列，从而生成带坐标的序列片段。

所有切分都使用 0-based、半开区间坐标，并惰性返回 `SequenceChunk`，不会一次性把全部结果读入内存。

<span id="1"></span>**序列切分总览**

API 入口如下：

- `dnakit.iter_sequence_chunks(value[必须], config[可选], start[可选], end[可选], source_id[可选], split[可选], region_index[可选], progress[可选])`：切分单条 `DNASequence` 或 `DNARecord`；
- `dnakit.iter_fasta_chunks(source[必须], config[可选], bed[可选], read_config[可选], progress[可选])`：逐条读取并切分 FASTA；
- `dnakit.ChunkingConfig(strategy[可选], length[可选], step[可选], min_length[可选], max_length[可选], num_samples[可选], lengths[可选], steps[可选], stage_steps[可选], include_partial[可选], seed[可选], split[可选], allow_gaps[可选])`：定义切分策略；
- `dnakit.make_length_curriculum(lengths[必须], stage_steps[可选])`：创建先短后长的长度训练计划；
- `dnakit.LengthCurriculum.to_config(window_step[可选], include_partial[可选], seed[可选], split[可选], allow_gaps[可选])`：把长度计划转换为切分配置；
- `dnakit.SequenceChunk.to_record()`：把切分结果转换为 `DNARecord`，无调用参数。

所有策略都返回惰性的 `SequenceChunk`，并保留来源 ID、起止坐标、split、策略和尺度编号等 provenance 信息。

## 1) 固定长度、不重叠 <span id="2"></span> {#fixed-non-overlapping-chunks}

**作用：** 从序列起点开始生成长度一致且互不重叠的片段，返回每段的来源坐标，适合建立无重复覆盖的训练样本或分析单元。

这是 `OPS-010 序列切分` 的默认策略。每个窗口长度由 `length` 指定，窗口之间不重叠，末尾不足一个完整窗口的部分默认丢弃。

```python
from dnakit import ChunkingConfig, iter_fasta_chunks

config = ChunkingConfig(length=1024)
for chunk in iter_fasta_chunks("genome.fa", config=config):
    print(chunk.id, chunk.source_start, chunk.source_end)
```

如果需要保留末尾不足长度的窗口，可设置 `include_partial=True`。

## 2) 固定长度、重叠滑窗 <span id="3"></span> {#overlapping-sliding-chunks}

**作用：** 按指定窗口长度和步长生成可重叠片段，使边界附近的碱基出现在多个窗口中，适合需要连续上下文覆盖的局部预测。

使用 `strategy="sliding"` 和 `step` 设置窗口起点之间的距离。`step` 小于 `length` 时窗口会重叠。

```python
from dnakit import ChunkingConfig, iter_sequence_chunks

config = ChunkingConfig(strategy="sliding", length=1024, step=512)
for chunk in iter_sequence_chunks(sequence, config=config):
    print(chunk.source_start, chunk.source_end)
```

`step` 不能大于窗口长度，否则会产生未覆盖的间隔。

## 3) 随机区间 <span id="4"></span> {#random-interval-chunks}

**作用：** 在允许的长度范围和坐标范围内抽取指定数量的随机片段，并通过 seed 保证复现，适合随机训练采样或抽样检查长序列。

随机策略在 `min_length` 和 `max_length` 范围内生成指定数量的随机区间。设置 `seed` 可以获得可复现结果。

```python
from dnakit import ChunkingConfig, iter_sequence_chunks

config = ChunkingConfig(
    strategy="random",
    min_length=512,
    max_length=2048,
    num_samples=100,
    seed=19,
)
chunks = iter_sequence_chunks(sequence, config=config)
```

## 4) 多尺度 <span id="5"></span> {#multiscale-chunks}

**作用：** 使用多个窗口长度分别生成短、中、长尺度片段并标记尺度编号，用于同时建模局部 motif 与较长范围上下文。

多尺度策略会按 `lengths` 中的每个长度生成窗口，适合同时提取局部和较长范围的上下文。

```python
from dnakit import ChunkingConfig, iter_fasta_chunks

config = ChunkingConfig(
    strategy="multiscale",
    lengths=(1024, 4096, 16384),
)
chunks = iter_fasta_chunks("genome.fa", config=config, bed="splits.bed")
for chunk in chunks:
    print(chunk.id, chunk.level_index, chunk.requested_length)
```

各尺度默认不重叠，也可以通过 `steps` 为每个尺度指定步长。

## 5) 先短后长 <span id="6"></span> {#length-curriculum-chunks}

**作用：** 定义训练阶段使用的序列长度及持续步数，再转换为对应切分配置，用于实现从短片段逐步增加到长片段的 curriculum 数据计划。

使用 `make_length_curriculum()` 创建逐步增加序列长度的训练计划，再通过 `to_config()` 转换为切分配置。

```python
from dnakit import iter_fasta_chunks, make_length_curriculum

curriculum = make_length_curriculum(
    (1024, 4096, 16384),
    stage_steps=(1000, 2000, 3000),
)
config = curriculum.to_config()
chunks = iter_fasta_chunks("genome.fa", config=config)
```

该计划只描述切分阶段，不会自动启动训练。

<span id="7-fastabed"></span>**FASTA、BED 与结果元数据**

`iter_sequence_chunks()` 用于切分单条 `DNASequence` 或 `DNARecord`；`iter_fasta_chunks()` 用于逐条读取并切分 FASTA。

不传 `bed` 时，所有 FASTA 记录默认使用 `train` 标签；传入 BED 时，要求 BED 第 1 列与 FASTA ID 完全相同，第 4 列可使用 `train`、`valid`、`test` 等标签。

每个 `SequenceChunk` 包含 `source_id`、`source_start`、`source_end`、`split`、`strategy` 和 `level_index`。调用 `to_record()` 后，来源 ID、起止坐标、split 和尺度编号会写入记录 metadata。

传入 `progress=callback` 可接收 `started`、`yielded`、`completed` 事件。当前 FASTA 读取器一次保留一条 FASTA 记录；处理染色体级大记录时，请根据实际内存调整 `ReadConfig` 的输入限制。
