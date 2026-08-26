# OPS-001 序列补全

对 DNA 序列进行逆序、互补或反向互补，从而得到方向转换后的新序列。

三种操作语义不同，所以保留三个名称。普通注释数据最常用的反向互补可直接接收并返回 `DNA`；纯 `reverse()`、`complement()` 仍是高级序列值操作。

## 1) OPS-001.1 逆序

- **作用：** 只把 DNA 序列的字符顺序倒过来而不替换碱基，用于检查或构造相反的位置顺序。
- **API：** `dnakit.ops.reverse(sequence[必须])`。
- **输入：** 必填一条 `DNASequence`；无额外必填参数。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.ops import reverse

seq = DNASequence("AATCGC")
print(reverse(seq).symbols)
```

- **示例结果：**

```text
CGCTAA
```

## 2) OPS-001.2 互补

- **作用：** 按碱基配对规则把每个 DNA 碱基替换为互补碱基，但保持原有位置顺序，用于构造对应互补链。
- **API：** `dnakit.ops.complement(sequence[必须])`。
- **输入：** 必填一条 `DNASequence`；无额外必填参数。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.ops import complement

seq = DNASequence("AATCGC")
print(complement(seq).symbols)
```

- **示例结果：**

```text
TTAGCG
```

## 3) OPS-001.3 反向互补

- **作用：** 同时反转序列顺序并替换为互补碱基，得到另一条链的 5′→3′ 表示，可用于反向链搜索、方向统一和引物分析。
- **API：** `dnakit.ops.reverse_complement(dna[必须], feature_policy[可选])`。
- **输入：** 普通用户传单记录 `DNA`；会同步 feature、strand 和逐碱基注释，并返回新的 `DNA`。
- **示例代码：**

```python
from dnakit import DNA
from dnakit.ops import reverse_complement

dna = DNA("AATCGC", id="seq-1")
print(reverse_complement(dna).symbols)
```

- **示例结果：**

```text
GCGATT
```

- **限制：** 三个函数都会保留 alphabet、topology 和 strandedness；`reverse()` 与 `reverse_complement()` 还会反转 multipart Gap 顺序。普通用户只需调用 `reverse_complement()`；`reverse_complement_record()` 仅在需要读取完整坐标变更审计时使用。
