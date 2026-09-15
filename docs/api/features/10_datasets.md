# 数据划分

公共处理接口支持[统一的文件与文本输入](00_processing_inputs.md)，具体入口及输出规则见该说明。

按照随机、分层、相似度或元数据约束划分 DNA 数据集，从而得到训练、验证和测试等相互隔离的数据子集。

所有随机或启发式过程都记录 seed、配置、分组和资源上限。

## 1) DATA-012 · 随机与稳定哈希划分

**作用：** 按目标比例把记录分配到 `train / valid / test` 等命名子集，或用固定 seed 生成 K-fold 划分并选择一折作为 `valid`、其余折作为 `train`。返回各子集和分配清单。


**API：** `dnakit.datasets.split(records[必须], config[必须])`；`config` 使用 `dnakit.datasets.SplitConfig`，本项支持 `method="random"`、`method="hash"` 和 `method="kfold"`。

**输入：** 必填记录；比例划分通过 `ratios` 指定总和为 1 的比例，默认 `train=0.8, valid=0.1, test=0.1`。K-fold 通过 `n_splits` 指定折数、`fold_id` 指定验证折。默认 `seed=42`、`shuffle=True`。`hash` 模式要求每条记录具有稳定且唯一的 `record.id`。

`SplitConfig` 和 `dnakit split` 的默认 seed 已由 `0` 调整为 `42`；需要复现此前省略 seed 的比例划分时，请显式传入 `seed=0`。

**示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [DNARecord(DNASequence("AC"), f"r{i}") for i in range(10)]
result = split(
    records,
    config=SplitConfig(
        method="random", ratios={"train": 0.8, "test": 0.2}, seed=17
    ),
)
print(dict(result.counts))
```

**示例结果：**

```text
{'train': 8, 'test': 2}
```

<span id="1-hash"></span>**顺序无关的 `hash` 划分**

`hash` 模式在 `shuffle=True` 时使用版本化 SHA-256 计算 `seed + record.id` 的稳定排序键，在 `shuffle=False` 时按 `record.id` 的稳定字节序排列，再按比例分配记录。因此，同一批记录即使输入顺序变化，按 `record.id` 查看时分类仍完全一致；它不使用 Python 内置的进程随机化 `hash()`。

**示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [DNARecord(DNASequence("AC"), f"r{i}") for i in range(10)]
config = SplitConfig(
    method="hash",
    ratios={"train": 0.6, "valid": 0.2, "test": 0.2},
    seed=17,
    preserve_order=False,
)
first = split(records, config=config)
second = split(list(reversed(records)), config=config)
first_by_id = {item.record_id: item.split for item in first.assignments}
second_by_id = {item.record_id: item.split for item in second.assignments}
print(first_by_id == second_by_id)
print(first.get("train").ids == second.get("train").ids)
```

**示例结果：**

```text
True
True
```

### 固定 K-fold 划分：train / valid

设置 `method="kfold"`，行为对应 [scikit-learn KFold](https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.KFold.html)：

```python
KFold(n_splits=10, shuffle=True, random_state=42)
```

**同一随机种子 → 同一套固定的 K-fold 划分 → 依次选择 K 个验证折。** 保持输入记录及其顺序、`seed`、`shuffle`、`n_splits` 不变，只改变 `fold_id`。例如 10 折使用 `fold_id=0..9`，每个样本恰好作为验证样本一次；不要使用 `seed + fold_id`。每次调用会重建相同的划分并返回所选折。

| `SplitConfig` 参数 | 默认值 | 含义 |
| --- | --- | --- |
| `method` | `"random"` | 使用 `"kfold"` 启用 K-fold；`"random"` / `"hash"` 保留比例划分 |
| `n_splits` | K-fold 模式下为 `10` | 总折数，可改为 5、10 等；要求 `2 <= n_splits <= 样本数` |
| `fold_id` | K-fold 模式下为 `0` | 当前验证折，范围为 `0` 到 `n_splits - 1` |
| `seed` | `42` | 生成整套划分的固定随机种子，K-fold 要求 `0 <= seed < 2**32` |
| `shuffle` | `True` | 划分前打乱；`False` 时按输入顺序生成连续的验证折，seed 不影响分配 |
| `ratios` | 比例模式为 `0.8 / 0.1 / 0.1` | K-fold 模式请省略；`train / valid` 目标比例由折数推导 |

K-fold 返回的子集名称固定为 `train` 和 `valid`，每个子集内部保持原始输入顺序，与 scikit-learn 一致；`preserve_order=False` 也不会额外打乱子集。样本数不能整除 K 时，前 `样本数 % K` 个验证折多一个样本，各验证折大小相差至多 1。结果和导出的分配清单记录 `seed`、`shuffle`、`n_splits`、`fold_id` 及每条记录的归属。

此模式需要 scikit-learn；可安装 `python -m pip install -e ".[bio-analysis]"`。已安装 `.[ml]` 或 `.[training]` 的环境也包含所需依赖。普通比例划分不需要这项依赖。

**示例：先留出固定 test，再对 train 做 10-fold。**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [DNARecord(DNASequence("ACGT"), f"r{i}") for i in range(100)]
seed = 42
n_splits = 10  # 改为 5 即可运行 5 折

# 第一步只做一次，留出固定测试集。
outer = split(
    records,
    config=SplitConfig(method="random", ratios={"train": 0.8, "test": 0.2}, seed=seed),
)
training_pool = outer.get("train")
test = outer.get("test")

# 每轮只改变 fold_id；测试集不传入 K-fold。
for fold_id in range(n_splits):
    fold = split(
        training_pool,
        config=SplitConfig(
            method="kfold", n_splits=n_splits, fold_id=fold_id,
            shuffle=True, seed=seed,
        ),
    )
    train = fold.get("train")
    valid = fold.get("valid")
    print(f"fold={fold_id}: train={len(train)}, valid={len(valid)}, test={len(test)}")
    # 在这里创建并训练一个新的模型，再用 valid 评估。
```

上述示例的 `fold=0` 到 `fold=9` 均为 `train=72, valid=8, test=20`，验证样本逐折轮换，20 条测试样本保持固定。若已准备好训练数据，直接将其作为 `training_pool` 即可。`n_splits` 控制总折数，`fold_id` 选择某一折；遍历所有折才形成完整 K-fold 流程，`split()` 本身只划分数据，不会自动训练模型或汇总指标。

接入 `train_custom_property_model()` 时，需要在每一轮把当前 `train`、`valid` 和固定 `test` 写成训练接口支持的 split 列值，并用 `DataMapping.split_field` 映射；不能直接把 `fold_0` 等名称当作训练接口的 split 值。相关说明见[自定义模型训练](23_deep_learning_property_prediction.md)。

**CLI：** 输入已留出的训练池，输出当前验证折对应的两个子集和分配清单：

```bash
dnakit split train.fasta runs/fold_0 --method kfold --n-splits 10 --fold-id 0 --seed 42
```

运行其余折时保持输入和 seed 不变，修改 `--fold-id` 和输出目录；可用 `--no-shuffle` 关闭打乱。Agent/MCP 的 `dnakit.datasets.split` 使用相同的 `config` 参数；JSON/YAML workflow 的 `split.params` 也接受 `method`、`n_splits`、`fold_id`，随机种子使用 workflow 的全局 `seed`，请设为 `42`。

## 2) DATA-013 · 分层随机划分

**作用：** 在 train、valid、test 中尽量保持指定类别标签的原始比例，返回分层后的子集及标签统计，减少类别分布偏移。


**API：** `dnakit.datasets.split(records[必须], config[必须])`；`config` 使用 `dnakit.datasets.SplitConfig`，本项设置 `method="stratified"`。

**输入：** 必填记录、比例和 `metadata_key`；可选 seed、缺失字段策略和顺序策略。

**示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [
    DNARecord(DNASequence("AC"), f"p{i}", metadata={"label": "p"})
    for i in range(4)
] + [
    DNARecord(DNASequence("GT"), f"n{i}", metadata={"label": "n"})
    for i in range(4)
]
result = split(
    records,
    config=SplitConfig(
        method="stratified", ratios={"train": 0.5, "test": 0.5},
        metadata_key="label", seed=2,
    ),
)
print(dict(result.counts))
```

**示例结果：**

```text
{'train': 4, 'test': 4}
```

## 3) DATA-014 · 相似度划分

**作用：** 先按相似度阈值建立序列组，再以整组为单位分配 train、valid、test，防止近似序列跨集合形成数据泄漏。


**API：** `dnakit.datasets.split(records[必须], config[必须])`；`config` 使用 `dnakit.datasets.SplitConfig`，本项设置 `method="similarity"`。

**输入：** 必填记录和比例；可选 k、阈值、IUPAC/Gap 策略、seed 和规模上限。

**示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [
    DNARecord(DNASequence("AAAA"), "a"),
    DNARecord(DNASequence("AAAT"), "b"),
    DNARecord(DNASequence("CCCC"), "c"),
    DNARecord(DNASequence("CCCG"), "d"),
]
result = split(
    records,
    config=SplitConfig(
        method="similarity", ratios={"train": 0.5, "test": 0.5},
        similarity_k=2, similarity_threshold=0.5, seed=5,
    ),
)
print({item.record_id: item.split for item in result.assignments})
```

**示例结果：**

```text
{'a': 'train', 'b': 'train', 'c': 'test', 'd': 'test'}
```

## 4) DATA-015 · Cluster split

**作用：** 使用已有 cluster 标签作为不可拆分分组进行数据划分，保证同一 cluster 的全部记录只进入一个子集，并报告实际比例偏差。


**API：** `dnakit.datasets.split(records[必须], config[必须])`；`config` 使用 `dnakit.datasets.SplitConfig`，本项设置 `method="group", metadata_key="cluster"`。

**输入：** 必填记录、每条记录的 cluster metadata 和比例；可选 seed、缺失值策略。

**示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [
    DNARecord(DNASequence("AC"), "a1", metadata={"cluster": "A"}),
    DNARecord(DNASequence("AG"), "a2", metadata={"cluster": "A"}),
    DNARecord(DNASequence("GT"), "b1", metadata={"cluster": "B"}),
]
result = split(
    records,
    config=SplitConfig(
        method="group", ratios={"train": 0.5, "test": 0.5},
        metadata_key="cluster", seed=11,
    ),
)
print({item.record_id: item.split for item in result.assignments})
```

**示例结果：**

```text
{'a1': 'train', 'a2': 'train', 'b1': 'test'}
```

## 5) DATA-016 · 物种划分

**作用：** 按物种 metadata 将记录成组后整体划分，确保同一物种不同时出现在训练集和评估集，用于检验跨物种泛化。


**API：** `dnakit.datasets.split(records[必须], config[必须])`；使用 `dnakit.datasets.SplitConfig` 的 `group` 模式。

**输入：** 必填 species metadata key 和比例；可选 seed、缺失值策略。

**示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [
    DNARecord(DNASequence("AC"), "h1", metadata={"species": "human"}),
    DNARecord(DNASequence("AG"), "h2", metadata={"species": "human"}),
    DNARecord(DNASequence("GT"), "m1", metadata={"species": "mouse"}),
]
result = split(records, config=SplitConfig(
    method="group", ratios={"train": 0.5, "test": 0.5}, metadata_key="species"
))
print({item.record_id: item.split for item in result.assignments})
```

**示例结果：**

```text
{'h1': 'train', 'h2': 'train', 'm1': 'test'}
```

## 6) DATA-017 · 染色体划分

**作用：** 按染色体 metadata 将记录成组后整体划分，避免同一染色体的高度相关区域跨 train、valid、test 造成泄漏。


**API：** `dnakit.datasets.split(records[必须], config[必须])`；使用 `dnakit.datasets.SplitConfig` 的 `group` 模式。

**输入：** 必填 chromosome metadata key 和比例；可选 seed、缺失值策略。

**示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [
    DNARecord(DNASequence("AC"), "c1-a", metadata={"chromosome": "chr1"}),
    DNARecord(DNASequence("AG"), "c1-b", metadata={"chromosome": "chr1"}),
    DNARecord(DNASequence("GT"), "c2-a", metadata={"chromosome": "chr2"}),
]
result = split(records, config=SplitConfig(
    method="group", ratios={"train": 0.5, "test": 0.5}, metadata_key="chromosome"
))
print({item.record_id: item.split for item in result.assignments})
```

**示例结果：**

```text
{'c1-a': 'train', 'c1-b': 'train', 'c2-a': 'test'}
```

## 7) DATA-018 · 个体划分

**作用：** 按个体或 donor metadata 对样本分组并整体划分，确保同一个体不跨数据子集，用于避免个体特异信息泄漏。


**API：** `dnakit.datasets.split(records[必须], config[必须])`；使用 `dnakit.datasets.SplitConfig` 的 `group` 模式。

**输入：** 必填 individual metadata key 和比例；可选 seed、缺失值策略。

**示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [
    DNARecord(DNASequence("AC"), "p1-a", metadata={"individual": "p1"}),
    DNARecord(DNASequence("AG"), "p1-b", metadata={"individual": "p1"}),
    DNARecord(DNASequence("GT"), "p2-a", metadata={"individual": "p2"}),
]
result = split(records, config=SplitConfig(
    method="group", ratios={"train": 0.5, "test": 0.5}, metadata_key="individual"
))
print({item.record_id: item.split for item in result.assignments})
```

**示例结果：**

```text
{'p1-a': 'train', 'p1-b': 'train', 'p2-a': 'test'}
```

## 8) 按自定义 label 划分

**作用：** 按调用方指定的任意 metadata label（如 family、locus、batch）对记录分组，确保相同 label 值只进入一个子集，并返回分组分配结果。

**API：** `dnakit.datasets.split(records[必须], config[必须])`；使用 `dnakit.datasets.SplitConfig(method="group", metadata_key=label)`。

**输入：** 必填记录、划分比例和自定义 `label`。`label` 是每条记录中的 metadata 字段名，例如 `donor`、`family`、`locus`、`batch` 或其他调用方提供的字段。

**示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

label = "family"  # 也可以替换为 donor、locus 或其他 metadata 字段名
records = [
    DNARecord(DNASequence("AC"), "f1-a", metadata={label: "f1"}),
    DNARecord(DNASequence("AG"), "f1-b", metadata={label: "f1"}),
    DNARecord(DNASequence("GT"), "f2-a", metadata={label: "f2"}),
]
result = split(
    records,
    config=SplitConfig(
        method="group",
        ratios={"train": 0.5, "test": 0.5},
        metadata_key=label,
    ),
)
print({item.record_id: item.split for item in result.assignments})
```

**示例结果：**

```text
{'f1-a': 'train', 'f1-b': 'train', 'f2-a': 'test'}
```

## 9) DATA-023 · 泄漏检测

**作用：** 对 train、valid、test 等集合执行跨集合精确或近似比较，返回泄漏序列对、相似度、涉及子集和汇总比例，用于验证划分独立性。


**API：** `dnakit.datasets.detect_leakage(splits[必须], config[可选])`；`config` 使用 `dnakit.datasets.LeakageConfig`。

**输入：** 必填 split 名称到 `DNASet` 的映射；可选 identity/edit/k-mer/fingerprint、阈值和资源上限。

**示例代码：**

```python
from dnakit import DNARecord, DNASequence, DNASet
from dnakit.datasets import LeakageConfig, detect_leakage

splits = {
    "train": DNASet([DNARecord(DNASequence("AAAA"), "a")]),
    "test": DNASet([DNARecord(DNASequence("AAAA"), "same")]),
}
report = detect_leakage(
    splits, config=LeakageConfig(method="identity", threshold=0.9)
)
print(report.has_leakage, report.exact_event_count)
```

**示例结果：**

```text
True 1
```
