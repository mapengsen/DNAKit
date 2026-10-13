# G6b：PLINK 关系矩阵、bed-reader 基因型矩阵与 BGEN 期望和剂量

这一组函数是 `dnakit.variants` 中基于 [pandas-plink](https://github.com/limix/pandas-plink)、[bed-reader](https://github.com/fastlmm/bed-reader) 和 [bgen-reader](https://github.com/limix/bgen-reader) 的适配器。关系矩阵读取器读取 GCTA GRM 与 PLINK 2 REL 文件；bed-reader 函数标准化、取子集并从 `.bed` 文件计算标准化关系矩阵；bgen-reader 函数计算无相位 BGEN 位点的等位基因期望与剂量，并生成元数据文件。所有函数只处理显式给出的本地文件或列表，不访问网络，也不做关联检验或质控过滤。

## 示例数据

下面的关系矩阵示例使用同一个 3 个样本的小矩阵，只用于演示接口，数值不代表真实群体。BGEN 示例由 `bgen_fixture.py` 在运行时写出两个位点、三个样本的文件。

**示例文件：** `study.grm`
```text
1	1	10	1.0
2	1	10	0.5
2	2	10	2.0
3	1	10	0.25
3	2	10	-0.5
3	3	10	1.5
```

**示例文件：** `study.grm.id`
```text
fam1	ind1
fam1	ind2
fam2	ind3
```

**示例文件：** `study.rel`
```text
1.0
0.5	2.0
0.25	-0.5	1.5
```

**示例文件：** `study.rel.id`
```text
fam1	ind1
fam1	ind2
fam2	ind3
```

**示例文件：** `bgen_fixture.py`
```python
import struct
import zlib


def write_bgen(path, samples, variants):
    """Write a diallelic, unphased BGEN 1.2 file with zlib-compressed 8-bit probabilities."""
    block = struct.pack("<I", len(samples)) + b"".join(
        struct.pack("<H", len(s)) + s.encode() for s in samples)
    block = struct.pack("<I", len(block) + 4) + block
    head = struct.pack("<III", 20, len(variants), len(samples)) + b"bgen"
    head += struct.pack("<I", 0x80000009)
    out = bytearray(struct.pack("<I", len(head) + len(block)) + head + block)
    for vid, rsid, chrom, pos, (a, b), probs in variants:
        item = b"".join(struct.pack("<H", len(t)) + t.encode() for t in (vid, rsid, chrom))
        item += struct.pack("<IH", pos, 2) + b"".join(struct.pack("<I", len(x)) + x.encode() for x in (a, b))
        data = struct.pack("<IHBB", len(samples), 2, 2, 2) + bytes(0x82 if p is None else 2 for p in probs)
        data += bytes([0, 8]) + bytes(v for p in probs for v in ([0, 0] if p is None else [round(x * 255) for x in p[:2]]))
        packed = zlib.compress(data)
        out += item + struct.pack("<II", len(packed) + 4, len(data)) + packed
    open(path, "wb").write(bytes(out))
```

## 1) `dnakit.variants.read_gcta_relationship_matrix()`

- **作用：** 读取 GCTA 遗传关系矩阵（GRM），文本或二进制文件集均可。
- **处理方法：** 调用 pandas-plink 的 `read_grm`。样本从 `.id` 文件读取（制表符分隔的家系号与个体号），超过 `max_cells` 允许的样本数时在读取矩阵前报错。`id_path` 缺省时按 pandas-plink 的规则推断：二进制文件 `x.bin` 对应 `x.id`；其他文件去掉末尾的 `.gz`（若有）后加 `.id`，例如 `study.grm` 对应 `study.grm.id`。
- **API：** `dnakit.variants.read_gcta_relationship_matrix(path, *, id_path=None, n_snps_path=None, max_cells=5000000, verbose=True)`。
- **输入：** GRM 文件路径；`id_path` 与 `n_snps_path` 可选；`max_cells` 限制 n × n 矩阵的元素数。
- **输出：** 含 `samples`（`family_id`、`sample_id`）、`matrix`（n × n 列表）、`n_snps`（每个存储的样本对一个计数，按文件顺序）和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.variants import read_gcta_relationship_matrix

result = read_gcta_relationship_matrix("study.grm")
print(result["samples"][0])
print(result["matrix"][1])
print(result["n_snps"][:3])
```

- **示例输出：**

```text
{'family_id': 'fam1', 'sample_id': 'ind1'}
[0.5, 2.0, -0.5]
[10, 10, 10]
```

## 2) `dnakit.variants.read_plink_relationship_matrix()`

- **作用：** 读取 PLINK 2 的实现关系矩阵（`.rel` 或 `.rel.bin`）。
- **处理方法：** 调用 pandas-plink 的 `read_rel`。文本文件按下三角逐行读入并对称化，二进制文件按 n × n 的 float64 读入。`.id` 文件名为去掉 `.bin`、`.zst` 或 `.gz` 后加 `.id`；`.zst` 文件需要 `zstandard` 包。
- **API：** `dnakit.variants.read_plink_relationship_matrix(path, *, id_path=None, max_cells=5000000, verbose=True)`。
- **输入：** REL 文件路径；`id_path` 可选；`max_cells` 限制 n × n 矩阵的元素数。
- **输出：** 含 `samples`、`matrix` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.variants import read_plink_relationship_matrix

result = read_plink_relationship_matrix("study.rel")
print(result["matrix"][2])
print(len(result["samples"]))
```

- **示例输出：**

```text
[0.25, -0.5, 1.5]
3
```

## 3) `dnakit.variants.standardize_genotypes()`

- **作用：** 按列标准化基因型矩阵，使每个位点的均值为 0、总体标准差为 1。
- **处理方法：** 调用 bed-reader 的 `standardize_f64`；`dtype="float32"` 时调用 `standardize_f32`。标准差按总体公式（ddof=0）计算。缺失值与零方差位点输出 0；零方差位点的标准差报告为 `None`。全部缺失的位点会报错，因为 bed-reader 无法计算其均值。
- **API：** `dnakit.variants.standardize_genotypes(genotypes, *, dtype="float64", max_cells=5000000, verbose=True)`。
- **输入：** 样本为行、位点为列的列表；元素为 0、1、2 或 `None`（缺失）；`dtype` 为 `"float64"` 或 `"float32"`。
- **输出：** 含 `standardized`（标准化矩阵）、`variant_statistics`（每列的 `mean` 与 `std`）和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.variants import standardize_genotypes

result = standardize_genotypes([[0, 1], [1, 1], [2, 1], [None, 2]])
print([round(value, 4) for value in result["standardized"][0]])
print(result["variant_statistics"][0])
```

- **示例输出：**

```text
[-1.2247, -0.5774]
{'index': 0, 'mean': 1.0, 'std': 0.816496580927726}
```

## 4) `dnakit.variants.calculate_genotype_relationship_matrix()`

- **作用：** 从 PLINK 1 `.bed` 文件计算标准化基因型关系矩阵 `Z Z' / m`。
- **处理方法：** 用 bed-reader 的 `open_bed` 读入全部基因型，跳过全部缺失的位点，用 `standardize_f64` 标准化，再用 NumPy 计算 `Z Z'`。`m` 为标准差大于零的位点数，缺失基因型贡献 0。矩阵乘法由 NumPy 完成，不是 bed-reader 的函数。
- **API：** `dnakit.variants.calculate_genotype_relationship_matrix(path, *, max_cells=5000000, verbose=True)`。
- **输入：** `.bed` 文件路径；同名的 `.bim` 与 `.fam` 文件必须存在；`max_cells` 限制样本数乘以位点数。
- **输出：** 含 `samples`、`matrix`（n × n 列表）、`variant_count`、`observed_variants`、`polymorphic_variants` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.variants import calculate_genotype_relationship_matrix
from dnakit.variants.genotype_formats import write_plink_bed

write_plink_bed(
    "panel.bed",
    [[0, 1, 2], [1, 1, 0], [2, None, 1], [0, 0, 2]],
    sample_ids=["a", "b", "c", "d"],
    variant_ids=["s1", "s2", "s3"],
    chromosomes=["1", "1", "1"],
    positions=[1, 2, 3],
    allele_1=["A", "A", "A"],
    allele_2=["G", "G", "G"],
    overwrite=True,
)
result = calculate_genotype_relationship_matrix("panel.bed")
print(result["observed_variants"], result["polymorphic_variants"])
print([round(value, 4) for value in result["matrix"][0]])
```

- **示例输出：**

```text
3 3
[0.7121, -0.3788, -0.5455, 0.2121]
```

## 5) `dnakit.variants.subset_genotype_matrix()`

- **作用：** 按给定的样本与位点下标，从内存中的基因型矩阵取出子矩阵。
- **处理方法：** 调用 bed-reader 的 `subset_f64_f64`；`dtype="float32"` 时调用 `subset_f32_f32`。输出顺序与下标顺序一致，下标允许重复，缺失值保持为 `None`。
- **API：** `dnakit.variants.subset_genotype_matrix(genotypes, *, sample_indices, variant_indices, dtype="float64", max_cells=5000000, verbose=True)`。
- **输入：** 基因型矩阵；`sample_indices` 与 `variant_indices` 为非空的整数列表，超出范围时报错。
- **输出：** 含 `sample_indices`、`variant_indices`、`genotypes` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.variants import subset_genotype_matrix

result = subset_genotype_matrix(
    [[0, 1, 2], [None, 1, 0], [2, 2, 2]],
    sample_indices=[0, 2],
    variant_indices=[1, 2],
)
print(result["genotypes"])
```

- **示例输出：**

```text
[[1.0, 2.0], [2.0, 2.0]]
```

## 6) `dnakit.variants.read_bgen_allele_expectation()`

- **作用：** 计算无相位 BGEN 位点中每个样本的等位基因期望拷贝数。
- **处理方法：** 用 bgen-reader 打开文件，对选定位点调用 `allele_expectation`。样本期望为各基因型概率乘以该基因型中等位基因的拷贝数之和，例如 AG 对 A 和 G 各贡献 1。有相位的位点会被拒绝，缺失样本输出 `None`。
- **API：** `dnakit.variants.read_bgen_allele_expectation(path, *, variant_indices=None, max_cells=5000000, max_variants=100000, verbose=True)`。
- **输入：** BGEN 文件路径；`variant_indices` 缺省时读取全部位点；`max_variants` 与 `max_cells` 为显式上限。
- **输出：** 含 `samples`、`variants`、`expectations`（变异、样本、等位基因的顺序）和 `provenance` 的字典。
- **示例代码：**

```python
from bgen_fixture import write_bgen
from dnakit.variants import read_bgen_allele_expectation

write_bgen(
    "small.bgen",
    ["s1", "s2", "s3"],
    [
        ("v1", "rs1", "1", 100, ("A", "G"), [[51 / 255, 76 / 255, 128 / 255], [1, 0, 0], None]),
        ("v2", "rs2", "1", 200, ("C", "T"), [[0, 0, 1], [64 / 255, 128 / 255, 63 / 255], None]),
    ],
)
result = read_bgen_allele_expectation("small.bgen")
print(result["variants"][0]["rsid"], result["samples"])
print([round(value, 4) for value in result["expectations"][0][0]])
print(result["expectations"][0][2])
```

- **示例输出：**

```text
rs1 ['s1', 's2', 's3']
[0.698, 1.302]
None
```

## 7) `dnakit.variants.read_bgen_dosage()`

- **作用：** 计算每个样本中指定替代等位基因的剂量（期望拷贝数）。
- **处理方法：** 对 `read_bgen_allele_expectation` 的期望调用 `compute_dosage`。`alt_allele_index` 为 `None` 时使用 bgen-reader 的默认值，即最后一个等位基因列。注意：bgen-reader 的 `compute_dosage` 文档写的是次要等位基因（MAF），但代码实际取最后一列。
- **API：** `dnakit.variants.read_bgen_dosage(path, *, variant_indices=None, alt_allele_index=None, max_cells=5000000, max_variants=100000, verbose=True)`。
- **输入：** 与 `read_bgen_allele_expectation` 相同，另加 `alt_allele_index`（非负整数，超出等位基因数时报错）。
- **输出：** 含 `samples`、`variants`（每个位点附带 `alt_allele_index` 与 `alt_allele`）、`dosages` 和 `provenance` 的字典。
- **示例代码：**

```python
from bgen_fixture import write_bgen
from dnakit.variants import read_bgen_dosage

write_bgen(
    "small.bgen",
    ["s1", "s2", "s3"],
    [
        ("v1", "rs1", "1", 100, ("A", "G"), [[51 / 255, 76 / 255, 128 / 255], [1, 0, 0], None]),
        ("v2", "rs2", "1", 200, ("C", "T"), [[0, 0, 1], [64 / 255, 128 / 255, 63 / 255], None]),
    ],
)
result = read_bgen_dosage("small.bgen", alt_allele_index=0)
print([variant["alt_allele"] for variant in result["variants"]])
print([round(value, 4) for value in result["dosages"][0][:2]])
```

- **示例输出：**

```text
['A', 'C']
[0.698, 2.0]
```

## 8) `dnakit.variants.create_bgen_metafile()`

- **作用：** 为 BGEN 文件生成变异元数据文件，之后用 `read_bgen` 读取会更快。
- **处理方法：** 调用 bgen-reader 的 `create_metafile`。它不会覆盖已有文件，因此 `overwrite=True` 时先删除旧文件。生成后用 `read_bgen(metafile_filepath=...)` 读回变异表，确认文件可用并返回变异数。
- **API：** `dnakit.variants.create_bgen_metafile(path, metafile_path, *, overwrite=False, verbose=True)`。
- **输入：** BGEN 文件路径与元数据文件路径；元数据文件不能与 BGEN 文件相同，其目录必须存在。
- **输出：** 含 `path`、`metafile`、`variants`、`bytes` 和 `provenance` 的字典。
- **示例代码：**

```python
from bgen_fixture import write_bgen
from dnakit.variants import create_bgen_metafile

write_bgen(
    "small.bgen",
    ["s1", "s2", "s3"],
    [
        ("v1", "rs1", "1", 100, ("A", "G"), [[51 / 255, 76 / 255, 128 / 255], [1, 0, 0], None]),
        ("v2", "rs2", "1", 200, ("C", "T"), [[0, 0, 1], [64 / 255, 128 / 255, 63 / 255], None]),
    ],
)
result = create_bgen_metafile("small.bgen", "small.metafile", overwrite=True)
print(result["variants"], result["metafile"].endswith("small.metafile"))
```

- **示例输出：**

```text
2 True
```

## 安装可选依赖

```bash
python -m pip install "dnakit[bio-analysis]"      # numpy，本组所有函数都需要
python -m pip install "dnakit[bio-pandas-plink]"  # pandas-plink 与 xarray，关系矩阵读取
python -m pip install "dnakit[bio-bed-reader]"    # bed-reader，标准化、取子集与 GRM
python -m pip install "dnakit[bio-bgen-reader]"   # bgen-reader，BGEN 期望、剂量与元数据
```

各包的项目地址：[pandas-plink](https://github.com/limix/pandas-plink)、[bed-reader](https://github.com/fastlmm/bed-reader)、[bgen-reader](https://github.com/limix/bgen-reader)。

## 验证方式

测试文件是 `others/tests/unit/test_full_G6b.py`，共 20 项，主要检查：

- GRM 与 REL 的二进制及文本文件读回的矩阵与手写数值一致，并与直接调用 pandas-plink 的 `read_grm`、`read_rel` 结果相同。
- `standardize_genotypes` 与手算的均值、总体标准差一致，float32 结果与 float64 结果相差在 1e-6 以内；`subset_genotype_matrix` 与直接调用 bed-reader 的子集内核一致。
- `calculate_genotype_relationship_matrix` 与按公式逐项计算的结果一致，全部缺失和单态位点的处理符合文档。
- BGEN 期望与剂量与手算的 k/255 概率一致，并与 `bgen_reader.read_bgen` 加 `allele_expectation`、`compute_dosage` 的结果一致；有相位位点被拒绝；元数据文件可以用 `read_bgen` 读回。
- 缺少可选依赖时给出安装提示；超出 `max_cells` 或 `max_variants`、非法下标和越界的等位基因索引都会被拒绝。
- 所有 BGEN 与 PLINK 文件都在 `tmp_path` 中生成，测试不访问网络。

## 当前范围

这一批包装 pandas-plink 的 GRM 与 REL 读取，bed-reader 的标准化、取子集和标准化关系矩阵，以及 bgen-reader 的期望、剂量和元数据文件。未包装的上游函数如下：

- 底层 Rust 内核（`read_f32`、`read_f64`、`read_i8`、`write_*`、`encode1_*`、`file_*_piece` 等）：`read_*` 与 `write_*` 已由 `read_plink_bed` 与 `write_plink_bed` 覆盖，其余是内部线性代数或编码步骤。`subset_f32_f64` 为混合精度变体，未对外开放。
- 云端读取（`read_cloud_*`、`url_to_bytes`、`check_file_cloud`）与示例数据下载（`sample_file`、`sample_url`、`example_filepath`）：需要网络。
- 已由 `genotype_formats.py` 调用的上游函数（已包装）：`read_plink`、`read_plink1_bin`、`write_plink1_bin`、`to_bed`，以及类 `open_bed`、`open_bgen`。与现有函数功能重复的：`read_bgen`、`allele_frequency`，以及类 `create_bed`。
- 测试套件、conftest、`tmp_path`、`get_num_threads` 等辅助函数，以及 `Chunk` 参数类。

示例从 `dnakit.variants` 导入。
