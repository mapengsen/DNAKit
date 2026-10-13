# 无外部依赖的 VCF 记录操作

这一页记录 DNAKit 已实现的 8 个变异记录级函数。它们只处理显式提供的 VCF 风格文本或记录，不调用变异检测器、参考基因组、注释数据库，也不推断缺失基因型。每个函数都返回 `status`、`operation`、`result` 和 `provenance`，适合在论文方法中逐项记录输入格式和阈值。

## 1) 解析 VCF 记录

- **作用：** 将 VCF 文本解析为可序列化的表头、样本和记录映射。
- **处理方法：** 校验至少 8 个 VCF 列、1-based `POS`、`QUAL` 数值和 `ALT`；解析 INFO、FORMAT、样本字段及 GT，不进行参考序列校验。
- **API：** `dnakit.variants.parse_vcf_records(input[必须], options[可选])`。
- **输入：** VCF 字符串、UTF-8 bytes、VCF 行序列，或包含 `text`/`vcf`/`content` 的映射；`options.strict=False` 可跳过格式错误行并记录警告。
- **输出：** `result` 包含 `headers`、`samples`、`records`、`record_count` 和 `warnings`。
- **示例代码：**

```python
from dnakit.variants import parse_vcf_records

text = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\trs1\tA\tG\t50\tPASS\tDP=12"
result = parse_vcf_records(text, verbose=False)
print(result["result"]["records"][0]["alt"])
```

- **示例结果：**

```text
['G']
```

## 2) 拆分多等位记录

- **作用：** 将一个含多个 ALT 的记录拆为多个双等位记录。
- **处理方法：** 为每个 ALT 生成一条记录，写入 1-based `original_alt_index`；默认把 GT 中当前 ALT 映射为 `1`，其他 ALT 映射为 `0`，缺失仍为 `.`。
- **API：** `dnakit.variants.split_multiallelic(input[必须], options[可选])`。
- **输入：** 记录序列或包含 `variants`/`records` 的映射；每条记录需有 `chrom`、`pos`、`ref`、`alt`。设置 `remap_genotypes=False` 可保留原 GT 字符串。
- **输出：** `result.records` 为拆分后的记录，并报告输入/输出记录数。
- **示例代码：**

```python
from dnakit.variants import split_multiallelic

result = split_multiallelic({"variants": [{
    "chrom": "chr1", "pos": 10, "ref": "A", "alt": ["G", "T"],
    "genotypes": ["0/1", "1/2"],
}]}, verbose=False)
print([row["alt"] for row in result["result"]["records"]])
```

- **示例结果：**

```text
['G', 'T']
```

## 3) 统计记录和变异类型

- **作用：** 统计记录数、ALT 等位基因数、染色体分布及 SNP/MNP/插入/缺失类别。
- **处理方法：** 按显式 `CHROM/POS/REF/ALT` 计算；多 ALT 记录按 ALT 等位基因分别计数，不执行规范化。
- **API：** `dnakit.variants.count(input[必须], options[可选])`。
- **输入：** 记录序列或包含 `variants`/`records` 的映射。
- **输出：** `result` 包含 `record_count`、`alternate_allele_count`、`by_chromosome` 和 `allele_type_counts`。
- **示例代码：**

```python
from dnakit.variants import count

result = count({"variants": [
    {"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G"},
    {"chrom": "chr1", "pos": 20, "ref": "C", "alt": "CT"},
]}, verbose=False)
print(result["result"]["allele_type_counts"])
```

- **示例结果：**

```text
{'insertion': 1, 'snp': 1}
```

## 4) 精确比较两组变异

- **作用：** 按 `CHROM/POS/REF/ALT` 精确比较两个记录集合。
- **处理方法：** 将每个 ALT 展开为键，计算交集、左右差集和 Jaccard 指数；不执行左对齐、等价变异归一化或参考序列校正。
- **API：** `dnakit.variants.compare(input[必须], options[可选])`。
- **输入：** 包含 `left` 和 `right` 的映射；两侧分别是记录序列或包含 `variants`/`records` 的映射。
- **输出：** `result` 包含左右集合大小、`only_left`、`only_right`、`shared` 和 `jaccard`。
- **示例代码：**

```python
from dnakit.variants import compare

result = compare({
    "left": {"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G"}]},
    "right": {"variants": [
        {"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G"},
        {"chrom": "chr1", "pos": 20, "ref": "C", "alt": "T"},
    ]},
}, verbose=False)
print(result["result"]["jaccard"])
```

- **示例结果：**

```text
0.5
```

## 5) 计算次要等位基因频率

- **作用：** 从显式基因型计算每个位点的 MAF。
- **处理方法：** 统计已调用的参考和所有声明 ALT 等位基因；零计数 ALT 也纳入频率，因此全参考位点的 MAF 为 0。缺失 GT 不参与分母，不进行插补。频率并列时按等位基因索引最小者输出 `minor_allele`，不改变 MAF 数值；例如 `1/2` 且参考计数为 0 时，参考等位基因是确定性的最小频率等位基因。
- **API：** `dnakit.variants.calculate_minor_allele_frequency(input[必须], options[可选])`。
- **输入：** VCF 风格记录，必须提供 `genotypes`；不允许只凭 ALT 字符串推断频率。
- **输出：** 每个位点返回 `allele_frequencies`、`minor_allele`、`maf` 和 `called_alleles`，并返回 `mean_maf`。
- **示例代码：**

```python
from dnakit.variants import calculate_minor_allele_frequency

result = calculate_minor_allele_frequency({"variants": [{
    "chrom": "chr1", "pos": 10, "ref": "A", "alt": "G",
    "genotypes": ["0/1", "1/1", "0/0"],
}]}, verbose=False)
print(result["result"]["variants"][0]["maf"])
```

- **示例结果：**

```text
0.5
```

## 6) 按 QUAL 过滤变异

- **作用：** 保留 QUAL 不低于指定阈值的记录。
- **处理方法：** 只使用记录中的显式 `qual`/`QUAL`/`quality` 数值；缺失 QUAL 默认删除，设置 `keep_missing=True` 才保留。
- **API：** `dnakit.variants.filter_variants_by_quality(input[必须], options[必须])`。
- **输入：** 记录序列和 `options.min_quality`；可选 `keep_missing`。
- **输出：** `result.variants` 为保留记录，`metrics` 记录每条记录的质量值和保留决定。
- **示例代码：**

```python
from dnakit.variants import filter_variants_by_quality

result = filter_variants_by_quality({"variants": [
    {"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "qual": 50},
    {"chrom": "chr1", "pos": 20, "ref": "C", "alt": "T", "qual": 10},
]}, options={"min_quality": 20}, verbose=False)
print(result["result"]["output_record_count"])
```

- **示例结果：**

```text
1
```

## 7) 按基因型缺失率过滤

- **作用：** 删除缺失基因型比例高于阈值的记录。
- **处理方法：** 根据 GT 中的 `.`，或显式 `missing_count/sample_count` 计算位点缺失率；不进行缺失填补。
- **API：** `dnakit.variants.filter_variants_by_missingness(input[必须], options[必须])`。
- **输入：** 带基因型或缺失/样本计数字段的记录，以及 `options.max_missingness`（0 到 1）。
- **输出：** 返回保留记录、每个位点缺失率和过滤计数。
- **示例代码：**

```python
from dnakit.variants import filter_variants_by_missingness

result = filter_variants_by_missingness({"variants": [{
    "chrom": "chr1", "pos": 10, "ref": "A", "alt": "G",
    "genotypes": ["0/1", "./.", "0/0"],
}]}, options={"max_missingness": 0.5}, verbose=False)
print(result["result"]["metrics"][0]["missingness"])
```

- **示例结果：**

```text
0.3333333333333333
```

## 8) 按显式等位基因频率过滤

- **作用：** 按指定 ALT 等位基因的显式频率范围过滤记录。
- **处理方法：** 从 GT 或 `alt_count/allele_count` 计算频率；通过 `allele_index` 选择多 ALT 记录中的目标 ALT，不推断未调用样本。
- **API：** `dnakit.variants.filter_variants_by_allele_frequency(input[必须], options[可选])`。
- **输入：** 变异记录；可选 `min_frequency`、`max_frequency` 和正整数 `allele_index`，默认范围为 0 到 1。
- **输出：** 返回保留记录、每个位点频率、阈值和过滤计数。
- **示例代码：**

```python
from dnakit.variants import filter_variants_by_allele_frequency

result = filter_variants_by_allele_frequency({"variants": [
    {"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1"]},
    {"chrom": "chr1", "pos": 20, "ref": "C", "alt": "T", "genotypes": ["0/0", "0/0"]},
]}, options={"min_frequency": 0.25}, verbose=False)
print(result["result"]["output_record_count"])
```

- **示例结果：**

```text
1
```

这些函数的 `provenance.scientific_boundary` 明确表示它们是记录级确定性处理。论文仍需报告 DNAKit 版本、输入 VCF 版本、参考坐标约定、阈值和外部变异检测/注释软件版本。
