# 序列优化

## 1) 按约束优化 DNA 序列 {#bio-45}

- **作用：** 修改 DNA 以满足指定规则，必要时保留蛋白翻译。
- **计算方法：** 先满足硬约束，再优化加权目标；示例通过同义密码子替换去除 `GCT`。
- **API：** `dnakit.bioinformatics.optimization.specification(kind, region=None, weight=1, **parameters)`; `dnakit.bioinformatics.optimization.optimize_sequence(sequence, constraints=(), objectives=(), preserve_translation=False, seed=0)`
- **输入：** A/C/G/T 序列和规则；固定长度；[GC、模式、结构与密码子规则](../../faq/bioinformatics.md#optimization)。
- **输出：** 返回 `OptimizationResult`，包含优化序列、约束是否满足及编辑记录。
- **示例代码：**

```python
from dnakit.bioinformatics.optimization import optimize_sequence, specification

result = optimize_sequence(
    "GCTGCT", preserve_translation=True,
    constraints=[specification("avoid_pattern", pattern="GCT")], seed=0,
)
print(result.sequence, result.constraints_satisfied)
print(result.edits)
```

- **示例结果：**

```text
GCAGCA True
((2, 'T', 'A'), (5, 'T', 'A'))
```

[依赖与实现说明](../../faq/bioinformatics.md#optimization)。

## 2) 自定义序列约束与目标 {#bio-46}

- **作用：** 把自己的评价逻辑加入优化流程。
- **计算方法：** 函数返回 `SpecEvaluation`；约束最小化非负 `violation`，目标最大化 `utility`。
- **API：** `dnakit.bioinformatics.optimization.Specification(name, evaluator, weight=1)`; `dnakit.bioinformatics.optimization.SpecEvaluation(violation=0, utility=0, locations=(), message='')`
- **输入：** 接收 `(sequence, circular)` 的 Python 函数；返回有限数值。
- **输出：** 返回 `OptimizationResult`，包含优化序列和约束是否满足。
- **示例代码：**

```python
from dnakit.bioinformatics.optimization import Specification, SpecEvaluation, optimize_sequence

rule = Specification("starts_with_A",
                     lambda sequence, circular: SpecEvaluation(violation=float(sequence[0] != "A")))
result = optimize_sequence("CCGT", constraints=[rule])
print(result.sequence, result.constraints_satisfied)
```

- **示例结果：**

```text
ACGT True
```

[依赖与实现说明](../../faq/bioinformatics.md#optimization)。

## 3) 从 GenBank 注释读取优化规则 {#bio-47}

- **作用：** 把构建文件中的规则注释转为约束和目标。
- **计算方法：** 读取 feature 上 `dnakit_constraint` / `dnakit_objective` 的 JSON，并使用该 feature 区间。
- **API：** `dnakit.bioinformatics.optimization.specifications_from_record(record)`
- **输入：** Biopython SeqRecord，feature 必须为精确、正向简单区间；需要 `bio-formats`。
- **输出：** 返回约束和目标规则的元组，可直接用于序列优化。
- **示例代码：**

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation
from dnakit.bioinformatics.optimization import specifications_from_record

record = SeqRecord(Seq("ACGTACGT"), id="construct")
record.features = [SeqFeature(SimpleLocation(0, 8, strand=1),
    qualifiers={"dnakit_constraint": ['{"kind":"gc","minimum":0.4,"maximum":0.6}']})]
constraints, objectives = specifications_from_record(record)
print(len(constraints), len(objectives), constraints[0].evaluate(str(record.seq)).violation)
```

- **示例结果：**

```text
1 0 0.0
```

[依赖与实现说明](../../faq/bioinformatics.md#optimization)。

## 4) 用 BLAST 命中约束序列设计 {#bio-48}

- **作用：** 在优化中检测与指定数据库的同源 HSP。
- **计算方法：** 执行 BLAST 搜索，并按 identity 和命中长度阈值计算违反量；结果按序列缓存。
- **API：** `dnakit.bioinformatics.optimization.blast_match_specification(database, minimum_identity=0.8, minimum_length=20, evalue=1e-5)`
- **输入：** 已建立的本地 BLAST 核酸库前缀、blastn；搜索为空不证明绝对无同源性。
- **输出：** 返回 BLAST 匹配约束 `Specification`；评估时得到违规程度和命中区间。
- **示例代码：**

```python
from dnakit.bioinformatics.optimization import blast_match_specification

rule = blast_match_specification("local_db/database", minimum_identity=0.9, minimum_length=20)
evaluation = rule.evaluate("ACGTACGTACGTACGTACGTACGT")
print(evaluation.violation, evaluation.locations)
```

- **示例结果：**

返回 `SpecEvaluation`，包括违反量和命中位置；真实结果依赖数据库，此处未执行。

[依赖与实现说明](../../faq/bioinformatics.md#optimization)。
