# 序列生成

通过突变、插入、删除、片段重排、随机打乱或序列重组生成新的 DNA 序列。


## 1) 突变生成 {#1-ops-0051}

- **作用：** 按碱基概率生成替换突变。
- **计算方法：** 逐位独立以 mut_frac 概率命中；命中后等概率选择其他三种碱基之一。
- **API：** `dnakit.ops.evolution_generate(sequence[必须], augmentations=("mutation",), seed/rng[二选一], mut_frac[可选])`。
- **输入：** canonical DNASequence；概率在 [0,1]，默认 0.05；seed/rng 二选一。
- **输出：** 返回 `EvolutionGenerationResult`，包含生成的 `sequence` 和突变步骤 `steps`。
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

## 2) 插入生成 {#2-ops-0052}

- **作用：** 按位置概率插入随机 DNA 片段。
- **计算方法：** 逐位独立以 insert_frac 概率命中，在该位置后插入长度为 [insert_min,insert_max] 的随机片段。
- **API：** `dnakit.ops.evolution_generate(sequence[必须], augmentations=("insertion",), seed/rng[二选一], insert_frac[可选], insert_min[可选], insert_max[可选])`。
- **输入：** canonical DNASequence；概率在 [0,1]，插入长度为正整数，seed/rng 二选一。
- **输出：** 返回 `EvolutionGenerationResult`，包含插入后的 `sequence` 和插入步骤 `steps`。
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

## 3) 删除生成 {#3-ops-0053}

- **作用：** 按碱基概率删除序列位置。
- **计算方法：** 逐位独立以 delete_frac 概率命中并删除该碱基。
- **API：** `dnakit.ops.evolution_generate(sequence[必须], augmentations=("deletion",), seed/rng[二选一], delete_frac[可选])`。
- **输入：** canonical DNASequence；概率在 [0,1]，默认 0.05；seed/rng 二选一。
- **输出：** 返回 `EvolutionGenerationResult`，包含删除后的 `sequence` 和删除步骤 `steps`。
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

## 4) 类进化算法操作 {#4-ops-0054}

- **作用：** 组合多个随机序列变换。
- **计算方法：** 无放回选择操作，依次执行局部反向互补、删除、易位、插入、整条反向互补、突变。
- **API：** `dnakit.ops.evolution_generate(sequence[必须], augmentations[可选], max_augmentations[可选], hard_aug[可选], seed[可选], rng[可选], mut_frac[可选], insert_frac[可选], delete_frac[可选], insert_min/insert_max[可选], shift_min/shift_max[可选], invert_min/invert_max[可选], rc_prob[可选])`。
- **输入：** canonical DNASequence、候选操作和 seed/rng；选中操作不保证一定改变序列。
- **输出：** 返回 `EvolutionGenerationResult`，包含生成的 `sequence` 和各步操作 `steps`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.ops import evolution_generate

result = evolution_generate(
    DNASequence("AAGTCC"),
    augmentations=("reverse_complement",),
    rc_prob=1.0,
    seed=19,
)
print(result.sequence.symbols)
```

- **示例结果：**

```text
GGACTT
```

## 5) 片段重排 <span id="5"></span> {#5-ops-0055}

- **作用：** 随机交换、反向互补或复制序列片段。
- **计算方法：** 随机切成非空片段，再按 operation 交换片段顺序、反向互补一个片段或复制一个片段。
- **API：** `dnakit.ops.rearrange_generate(sequence[必须], operation[可选], segment_count[可选], seed[可选], rng[可选])`。
- **输入：** canonical DNASequence、片段数、操作和 seed/rng。
- **输出：** 返回 `RearrangementResult`，示例展示 `sequence`、`breakpoints`、`permutation`。
- **示例代码：**

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

- **示例结果：**

```text
AGTCAC
(1, 5) (1, 0, 2)
```

## 6) k-mer 随机打乱 <span id="6-k-mer"></span> {#6-ops-0056-k-mer}

- **作用：** 随机打乱序列并保持指定 k 的重叠 k-mer 计数。
- **计算方法：** 在 de Bruijn 多重图上采样 Euler 路径，精确保留重叠 k-mer 计数。
- **API：** `dnakit.ops.kmer_shuffle(sequence[必须], k[可选], seed[可选], rng[可选], ensure_different[可选], max_attempts[可选])`。
- **输入：** canonical DNASequence、k 和 seed/rng；默认必须产生不同序列，无其他可行结果时报错。
- **输出：** 返回 `KmerShuffleResult`，示例展示 `sequence`。
- **示例代码：**

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

- **示例结果：**

```text
AAGATCGGATCGATC
```

## 7) 多序列重组 / 交叉<span id="7-crossover"></span> {#7-ops-0057}

- **作用：** 以单个边界拼接两个亲本，生成子代。
- **计算方法：** 在 0-based 边界 position 处取第一条的前缀和第二条的后缀。
- **API：** `dnakit.ops.crossover(first[必须], second[必须], position[可选], seed[可选], rng[可选])`。
- **输入：** 两条等长、线性、无 Gap 的 A/C/G/T 序列；省略 position 时必须提供 seed/rng。
- **输出：** 返回 `CrossoverResult`，示例展示 `sequence`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.ops import crossover

first = DNASequence("AAAACCCC")
second = DNASequence("GGGGTTTT")
result = crossover(first, second, position=4)
print(result.sequence.symbols)  # AAAATTTT
```

- **示例结果：**

```text
AAAATTTT
```

[FAQ：补充说明](../../faq/tutorial_basics.md#generation-details)
