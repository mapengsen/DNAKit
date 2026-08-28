# OPS-005 序列生成

通过突变、插入、删除、片段重排、随机打乱或序列重组生成新的 DNA 序列。

## 1. OPS-005.1 突变生成

- **作用：** 对序列中的每个碱基独立判断是否发生替换突变。
- **API：** `dnakit.ops.evolution_generate(sequence[必须], augmentations=("mutation",), seed/rng[二选一], mut_frac[可选])`。
- **概率：** `mut_frac` 表示每个碱基发生突变的概率，范围为 `0.0～1.0`，默认值为 `0.05`。命中的碱基会等概率替换为另外三种 canonical 碱基之一，因此一定会改变。
- **补充：** `mutate()` 用于指定位置或随机位置的单次单碱基替换，不接收概率参数。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.ops import evolution_generate

seq = DNASequence("ACGTACGTACGT")
mutation_probability = 0.25
result = evolution_generate(
    seq,
    augmentations=("mutation",),
    seed=7,
    mut_frac=mutation_probability,  # 每个碱基有 25% 的独立突变概率
)
print(result.sequence.symbols)
print(result.steps[0].length)  # 实际发生突变的碱基数
```

- **示例结果：**

```text
ACATTCGACCTT
5
```

## 2. OPS-005.2 插入生成

- **作用：** 对序列中的每个碱基独立判断是否在其后插入随机 DNA。
- **API：** `dnakit.ops.evolution_generate(sequence[必须], augmentations=("insertion",), seed/rng[二选一], insert_frac[可选], insert_min[可选], insert_max[可选])`。
- **概率：** `insert_frac` 表示每个碱基位置发生插入的概率，范围为 `0.0～1.0`，默认值为 `0.05`。
- **插入长度：** `insert_min=1, insert_max=1` 时，每次命中只插入 1 个随机碱基；设置其他最小值和最大值时，每次命中会在闭区间 `[insert_min, insert_max]` 内随机选择片段长度。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.ops import evolution_generate

seq = DNASequence("ACGTACGTACGT")
insertion_probability = 0.25

single_base = evolution_generate(
    seq,
    augmentations=("insertion",),
    seed=7,
    insert_frac=insertion_probability,  # 每个位置有 25% 的独立插入概率
    insert_min=1,
    insert_max=1,
)
segment = evolution_generate(
    seq,
    augmentations=("insertion",),
    seed=7,
    insert_frac=insertion_probability,
    insert_min=2,
    insert_max=4,
)
print(single_base.sequence.symbols, single_base.steps[0].length)
print(segment.sequence.symbols, segment.steps[0].length)
```

- **示例结果：**

```text
ACGGTACTGTAACAGT 4
ACGAGACTTACAAACCGTTACAACCAGGT 17
```

如果只想在整条序列中随机选择一个位置，并插入一次连续片段，可使用 `indel_generate(operation="insertion", min_length=..., max_length=...)`。

## 3. OPS-005.3 删除生成

- **作用：** 对序列中的每个碱基独立判断是否删除该碱基。
- **API：** `dnakit.ops.evolution_generate(sequence[必须], augmentations=("deletion",), seed/rng[二选一], delete_frac[可选])`。
- **概率：** `delete_frac` 表示每个碱基位置发生删除的概率，范围为 `0.0～1.0`，默认值为 `0.05`。命中后只删除当前位置的一个碱基。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.ops import evolution_generate

seq = DNASequence("ACGTACGTACGT")
deletion_probability = 0.25
result = evolution_generate(
    seq,
    augmentations=("deletion",),
    seed=7,
    delete_frac=deletion_probability,  # 每个碱基有 25% 的独立删除概率
)
print(result.sequence.symbols)
print(result.steps[0].length)  # 实际删除的碱基数
```

- **示例结果：**

```text
ACTCGCT
5
```

如果只想在整条序列中随机选择一个连续片段并删除一次，可使用 `indel_generate(operation="deletion", min_length=..., max_length=...)`。

## 4. OPS-005.4 类进化算法操作

`evolution_generate()` 还支持局部反向互补、易位和整条反向互补。多个候选操作先无放回抽样，再按固定顺序执行：局部反向互补、删除、易位、插入、整条反向互补、突变。逐碱基概率针对各操作执行时的当前序列。

- **完整 API：** `dnakit.ops.evolution_generate(sequence[必须], augmentations[可选], max_augmentations[可选], hard_aug[可选], seed[可选], rng[可选], mut_frac[可选], insert_frac[可选], delete_frac[可选], insert_min/insert_max[可选], shift_min/shift_max[可选], invert_min/invert_max[可选], rc_prob[可选])`。
- **结果：** `EvolutionGenerationResult.sequence` 是生成序列；`steps` 保存每个操作是否生效及其实际影响长度；结果同时保存随机源、随机状态和 `dnakit-evoaug-v3` 算法版本。
- **抽样规则：** `hard_aug=True` 时每条序列恰好选择 `max_augmentations` 种操作；`False` 时从 1 到该最大值随机选择。选中某种操作不保证一定命中碱基，例如概率较低且序列较短时可能不产生变化。

| 操作名                 | 作用                                             |
| ---------------------- | ------------------------------------------------ |
| `mutation`           | 每个碱基独立以`mut_frac` 概率替换              |
| `insertion`          | 每个碱基后独立以`insert_frac` 概率插入随机 DNA |
| `deletion`           | 每个碱基独立以`delete_frac` 概率删除           |
| `translocation`      | 随机 roll 序列，交换断点两侧片段顺序             |
| `inversion`          | 将连续片段替换为其反向互补序列                   |
| `reverse_complement` | 按`rc_prob` 将整条序列替换为反向互补           |

## 5. OPS-005.5 片段重排 <span id="5"></span>

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

## 6. OPS-005.6 k-mer 随机打乱 <span id="6-k-mer"></span>

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

## 7. OPS-005.7 多序列重组 / 交叉<span id="7-crossover"></span>

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
