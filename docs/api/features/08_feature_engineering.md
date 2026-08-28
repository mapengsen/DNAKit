# DNA指纹

DNAKit 目前提供两种 DNA 指纹计算方法。每一位只有 `0` 或 `1`，可直接用于 Tanimoto/Jaccard 相似度、检索、聚类和去重。

## 1) Hashed k-mer 位指纹

- **作用：** 提取序列中出现的 k-mer，并将其稳定哈希到固定长度 0/1 位向量，用于 Tanimoto/Jaccard 相似度、聚类和检索，同时避免显式 `4^k` 特征空间过大。
- **特点：** 功能上最接近分子中的 Morgan/ECFP 位指纹；无论 `k` 多大，输出长度始终为 `n_bits`。
- **API：** `dnakit.fingerprints.hashed_kmer_fingerprint(value[必须], k[必须], n_bits[可选], canonical[可选], seed[可选], representation[可选], ambiguity_policy[可选], overlapping[可选], cross_gaps[可选])`
- **常用参数：** 默认 `n_bits=2048`、`canonical=True`、`seed=0`；`representation="sparse"` 只保存值为 `1` 的位。

**示例：**

```python
from dnakit import DNASequence
from dnakit.fingerprints import hashed_kmer_fingerprint

result = hashed_kmer_fingerprint(
    DNASequence("ACGTAC"),
    k=3,
    n_bits=16,
    seed=7,
    representation="sparse",
)

print(result.dimension)
print(dict(result.values))
```

**输出：**

```text
16
{'bit:1': 1, 'bit:10': 1}
```

## 2) Panel 存在性指纹

- **作用：** 把用户定义的命名模式面板转换成可解释的 0/1 位向量，每一位表示对应 motif 或识别序列是否存在，便于按已知功能模式比较和筛选序列。
- **特点：** 类似分子中的 MACCS Keys，位含义明确，适合检测已知 motif、启动子模式或酶切识别序列。
- **API：** `dnakit.fingerprints.panel_fingerprint(value[必须], panel[必须], mode[可选], overlapping[可选], representation[可选], max_panel_size[可选], max_matches_per_pattern[可选])`
- **模式：** 默认 `mode="iupac"`，也可使用 `mode="exact"`；默认扫描正反两条链。

**示例：**

```python
from dnakit import DNASequence
from dnakit.fingerprints import panel_fingerprint

panel = {
    "start": "ATG",
    "EcoRI": "GAATTC",
    "TATA_box": "TATAWAWR",
}
result = panel_fingerprint(DNASequence("ATGGAATTC"), panel)

print(result.feature_names)
print(result.dense_values())
```

输出：

```text
('panel:EcoRI', 'panel:TATA_box', 'panel:start')
(1, 0, 1)
```

- **解释：** `EcoRI` 和 `start` 对应模式存在，因此相应位为 `1`；`TATA_box` 不存在，因此为 `0`。
