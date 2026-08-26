# OPS-005 序列生成

通过突变、插入缺失、片段重排、随机打乱或序列重组生成新的 DNA 序列。

所有方法都返回新的 `DNASequence`，不修改输入；随机方法要求显式 `seed` 或 `random.Random`，以便复现。

## 1. OPS-005.1 突变生成

- **作用：** 根据指定位置生成单碱基替换，或按固定 seed 生成可复现的随机单点突变，同时记录原碱基和新碱基，用于突变模拟和数据增强。
- **API：** `dnakit.ops.mutate(sequence[必须], position[可选], replacement[可选], seed[可选], rng[可选], allowed_bases[可选])`；结果类型为 `MutationResult`。
- **输入：** 必填原 `DNASequence`；指定突变传 `position` 和 `replacement`，随机突变传 `seed` 或 `random.Random`，可选 `allowed_bases`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.ops import mutate

seq = DNASequence("AAAA")
specified = mutate(seq, position=1, replacement="T")
randomized = mutate(seq, seed=19)
print(specified.sequence.symbols)
print(randomized.sequence.symbols, randomized.seed)
```

- **示例结果：**

```text
ATAA
TAAA 19
```

- **限制：** 每次调用只生成一个 SNV，拒绝 no-op；含未知长度 Gap 的序列无法解析突变坐标。组合突变文库仍由分子生物学模块的 `generate_mutation_library()` 提供。

## 2. OPS-005.2 进化生成

`evolution_generate()` 将 EvoAug 的离散 DNA 增强操作适配为 `DNASequence` 级 API。候选操作先无放回抽样，再按论文规定的顺序执行：局部反向互补、缺失、易位、插入、整条反向互补、突变。

- **API：** `dnakit.ops.evolution_generate(sequence[必须], augmentations[可选], max_augmentations[可选], hard_aug[可选], seed[可选], rng[可选], mut_frac[可选], delete_min/delete_max[可选], insert_min/insert_max[可选], shift_min/shift_max[可选], invert_min/invert_max[可选], rc_prob[可选], pad_indels[可选])`。
- **结果：** `EvolutionGenerationResult.sequence` 是生成序列；`steps` 保存每个操作的坐标、长度、shift、是否实际生效和突变尝试数；结果同时保存随机源、随机状态和 `dnakit-evoaug-v1` 算法版本。
- **操作：**

| 操作名                 | 作用                                                                |
| ---------------------- | ------------------------------------------------------------------- |
| `mutation`           | 按`mut_frac` 对随机位置进行 A/C/G/T 替换                          |
| `deletion`           | 删除连续片段；`pad_indels=True` 时按 EvoAug 用随机 DNA 补齐       |
| `insertion`          | 在随机边界插入随机 DNA；`pad_indels=True` 时补齐到 `insert_max` |
| `translocation`      | 随机 roll 序列，交换断点两侧片段顺序                                |
| `inversion`          | 将连续片段替换为其反向互补序列                                      |
| `reverse_complement` | 按`rc_prob` 将整条序列替换为反向互补                              |

- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.ops import evolution_generate

seq = DNASequence("ACGT" * 10)
generated = evolution_generate(
    seq,
    augmentations=("insertion", "translocation", "mutation"),
    max_augmentations=2,
    seed=19,
    insert_min=1,
    insert_max=3,
    shift_min=1,
    shift_max=5,
    mut_frac=0.05,
    pad_indels=False,
)
print(generated.sequence.symbols)
print(generated.augmentations)
print(generated.steps)
```

- **抽样规则：** `hard_aug=True` 时每条序列恰好应用 `max_augmentations` 个候选操作；`False` 时从 1 到最大值随机选择。默认 `pad_indels=True` 与 EvoAug 的 one-hot 固定形状实现一致；设为 `False` 可得到真实变长的插入/缺失结果。

## 3. 插入 / 删除（Indel） <span id="4-indel"></span> {#indel}

`indel_generate()` 是对单个插入或删除的直接入口。默认返回自然变长结果；如果要模拟 EvoAug 固定长度 one-hot 的 padding 行为，可设置 `pad_indels=True`。

- **API：** `dnakit.ops.indel_generate(sequence[必须], operation[必须], min_length[可选], max_length[可选], seed[可选], rng[可选], pad_indels[可选])`。
- **操作：** `operation="insertion"` 在随机边界插入随机 A/C/G/T；`operation="deletion"` 删除随机连续片段。
- **结果：** 返回 `EvolutionGenerationResult`，实际区间和长度位于 `result.steps[0]`。

```python
from dnakit import DNASequence
from dnakit.ops import indel_generate

seq = DNASequence("ACGT" * 10)
inserted = indel_generate(
    seq,
    operation="insertion",
    min_length=2,
    max_length=6,
    seed=19,
)
deleted = indel_generate(
    seq,
    operation="deletion",
    min_length=2,
    max_length=6,
    seed=19,
)
print(inserted.sequence.symbols)
print(deleted.sequence.symbols)
```

## 4. 片段重排 <span id="5"></span> {#fragment-rearrangement}

`rearrange_generate()` 先用随机断点把序列切成 `segment_count` 个非空片段，再执行一种重排：

| `operation`   | 作用                                     |
| --------------- | ---------------------------------------- |
| `exchange`    | 随机交换片段顺序，长度不变               |
| `inversion`   | 选中一个片段并替换为其反向互补，长度不变 |
| `duplication` | 选中一个片段并在其后复制一份，长度增加   |

- **API：** `dnakit.ops.rearrange_generate(sequence[必须], operation[可选], segment_count[可选], seed[可选], rng[可选])`。
- **结果：** `RearrangementResult` 保存 `breakpoints`、交换 `permutation` 或 `selected_segment`，便于复核生成过程。

```python
from dnakit import DNASequence
from dnakit.ops import rearrange_generate

seq = DNASequence("AAGTCC")
result = rearrange_generate(
    seq,
    operation="exchange",
    segment_count=3,
    seed=19,
)
print(result.sequence.symbols)
print(result.breakpoints, result.permutation)
```

## 5. k-mer 保持的随机打乱 <span id="6-k-mer"></span> {#kmer-preserving-shuffle}

`kmer_shuffle()` 不直接独立打乱碱基，而是在原序列的 de Bruijn 多重图上随机采样 Euler 路径。因此会精确保持重叠 k-mer 的计数；例如 `k=2` 保持二联体计数，`k=3` 保持三联体计数。默认 `ensure_different=True`，会尝试生成与原序列不同的序列；若计数结构只有唯一重建结果，则抛出明确错误。

- **API：** `dnakit.ops.kmer_shuffle(sequence[必须], k[可选], seed[可选], rng[可选], ensure_different[可选], max_attempts[可选])`。
- **结果：** `KmerShuffleResult` 保存 `k`、原始 k-mer 计数、尝试次数和随机状态。

```python
from collections import Counter

from dnakit import DNASequence
from dnakit.ops import kmer_shuffle

seq = DNASequence("AAGATCGATCGGATC")
result = kmer_shuffle(seq, k=3, seed=19)
generated = result.sequence.symbols
assert Counter(
    generated[index : index + 3] for index in range(len(generated) - 3 + 1)
) == Counter(seq.symbols[index : index + 3] for index in range(len(seq) - 3 + 1))
print(generated)
```

该实现参考 k-mer 保持打乱的 Euler 路径方法（见 [uShuffle 论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC2375906/)）；它只保证指定 `k` 的重叠计数，不代表保持更高阶 k-mer、功能、表达量或实验性质。

## 6. 多序列重组 / 交叉（Crossover） <span id="7-crossover"></span> {#crossover}

`crossover()` 接收两条等长 DNA 序列和一个单点边界，生成“第一条序列前缀 + 第二条序列后缀”的子代。边界使用 0-based 半开切分位置；省略 `position` 时必须传 `seed` 或 `random.Random`，随机选择一个同时保留两侧碱基的边界。

- **API：** `dnakit.ops.crossover(first[必须], second[必须], position[可选], seed[可选], rng[可选])`。
- **结果：** `CrossoverResult.sequence` 是子代；`position`、两个亲本长度和随机状态保存在结果中。
- **示例：**

```python
from dnakit import DNASequence
from dnakit.ops import crossover

first = DNASequence("AAAACCCC")
second = DNASequence("GGGGTTTT")
result = crossover(first, second, position=4)
print(result.sequence.symbols)  # AAAATTTT
```

该 API 当前是两条等长序列的一点交叉；不同长度、环状、含 Gap、含 IUPAC 模糊碱基的输入会被拒绝。
