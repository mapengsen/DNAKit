# 变异命名与基因型文件：varcode、hgvs、PLINK 与 BGEN

这一组函数是 `dnakit.variants` 中基于可选依赖的变异命名与基因型文件适配器。变异命名部分使用 [varcode](https://github.com/openvax/varcode) 给出变异的最小表示和类型，并使用 [hgvs](https://github.com/biocommons/hgvs) 解析、格式化和检查 HGVS 字符串。基因型文件部分使用 [bed-reader](https://github.com/fastlmm/bed-reader) 与 [pandas-plink](https://github.com/limix/pandas-plink) 读写 PLINK 1 的 `.bed`、`.bim` 和 `.fam`，使用 [bgen-reader](https://github.com/limix/bgen-reader) 读取 BGEN 1.2 文件的元数据、概率和等位基因频率。这些函数不计算基因功能效应，也不做需要参考序列的 HGVS 操作（序列获取、UTA 投影、校验和规范化）。依赖都是可选的，只有调用对应函数时才会导入。

## 示例数据

下面的 VCF 只有 5 条记录，用于演示 VCF 的读取。PLINK 与 BGEN 是二进制格式，因此示例代码自己生成文件：bgen-reader 没有写入器，它的示例文件需要从网络下载，所以 BGEN 示例文件按 BGEN 1.2 格式由代码写出，并由 bgen-reader 读回核对。

**示例文件：** `records.vcf`
```text
##fileformat=VCFv4.2
##contig=<ID=1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	101	rs1	A	G	50	PASS	DP=5
1	200	.	ACG	A	40	PASS	DP=6
1	300	.	A	AGT	.	q10	DP=3
1	400	.	AC	TG	.	.	DP=3
1	500	.	A	G	.	PASS	DP=3
```

## 1) `dnakit.variants.describe_variant()`

- **作用：** 根据染色体、1-based 起点、REF 与 ALT 给出变异的最小表示和类型。
- **处理方法：** 调用 `varcode.Variant`。它去掉 REF 与 ALT 的公共前缀和公共后缀，得到 `trimmed_*` 坐标；类型按 SNV、插入、缺失、等长替换（substitution）和其他复杂变化（complex）划分，只对 SNV 判断是否为转换。不读取参考序列。
- **API：** `dnakit.variants.describe_variant(contig, start, ref, alt)`。
- **输入：** 染色体名（不含空白）、1-based 起点、REF 与 ALT（只含 A、C、G、T、N，不区分大小写）。REF 与 ALT 不能同时为空，也不能相同。
- **输出：** 包含 `contig`、`start`、`ref`、`alt`（输入形式）、`trimmed_start`、`trimmed_end`、`trimmed_ref`、`trimmed_alt`、`kind`、`is_transition` 与 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.variants import describe_variant

result = describe_variant("1", 200, "ACG", "A")
print(result["kind"], result["trimmed_start"], result["trimmed_end"], repr(result["trimmed_ref"]), repr(result["trimmed_alt"]))
```

- **示例输出：**

```text
deletion 201 202 'CG' ''
```

## 2) `dnakit.variants.load_vcf_variants()`

- **作用：** 用 varcode 读取 VCF 记录，并给出每条记录的最小表示和类型。
- **处理方法：** 调用 `varcode.load_vcf`。`distinct=False` 保留所有记录，不合并相同位点。`only_passing=True` 时只保留 FILTER 为 PASS 或缺失的记录。`genome` 是 varcode 需要的参考名称，本函数不读取参考序列。varcode 在 `max_variants` 处会静默截断，因此本函数多请求一条记录，超过上限时报错。
- **API：** `dnakit.variants.load_vcf_variants(path, *, genome="GRCh38", only_passing=True, max_variants=100000)`。
- **输入：** VCF 文件路径；`max_variants` 为 1 到 1000000 之间的整数。
- **输出：** `count` 与 `variants` 列表。每项字段与 `describe_variant()` 相同，`start`、`ref`、`alt` 保留 VCF 中的输入形式。
- **示例代码：**

```python
from dnakit.variants import load_vcf_variants

result = load_vcf_variants("records.vcf")
print(result["count"])
print([(item["start"], item["kind"]) for item in result["variants"]])
```

- **示例输出：**

```text
4
[(101, 'snv'), (200, 'deletion'), (400, 'substitution'), (500, 'snv')]
```

## 3) `dnakit.variants.parse_hgvs_variant()`

- **作用：** 解析一条 HGVS 字符串（c.、g.、n.、m.、r. 或 p.），返回坐标、编辑类型和碱基。
- **处理方法：** 调用 `hgvs.parser.Parser().parse_hgvs_variant`，只做语法解析，不访问参考数据。语法错误时报错。区间起点晚于终点（按碱基和偏移比较）时也报错，因为 hgvs 本身不检查这一点。
- **API：** `dnakit.variants.parse_hgvs_variant(text)`。
- **输入：** 一条不含空白的 HGVS 字符串，长度不超过 1000 个字符。
- **输出：** `reference`（参考序列 ID）、`gene`、`coordinate_type`、`start` 与 `end`（含 `base` 和 `offset`，p. 变异的 `offset` 为 None）、`edit_type`、`reference_bases`、`alternate_bases`、`formatted` 与 `provenance`。
- **示例代码：**

```python
from dnakit.variants import parse_hgvs_variant

result = parse_hgvs_variant("NM_004006.2:c.4375C>T")
print(result["coordinate_type"], result["start"]["base"], result["edit_type"], result["reference_bases"], result["alternate_bases"])
print(parse_hgvs_variant("NP_003997.1:p.Val600Glu")["formatted"])
```

- **示例输出：**

```text
c 4375 sub C T
NP_003997.1:p.Val600Glu
```

## 4) `dnakit.variants.format_hgvs_variant()`

- **作用：** 把 HGVS 字符串格式化为 hgvs 的标准写法，并说明是否改变了原文。
- **处理方法：** 解析后调用 hgvs 的 `format()`，使用默认格式设置。格式化不需要参考数据。例如 `c.123delA` 会写成 `c.123del`，因为 hgvs 的默认输出不写出被删除的碱基。
- **API：** `dnakit.variants.format_hgvs_variant(text)`。
- **输入：** 同 `parse_hgvs_variant()`。
- **输出：** `input`、`formatted`、`changed` 与 `provenance`。
- **示例代码：**

```python
from dnakit.variants import format_hgvs_variant

result = format_hgvs_variant("NM_004006.2:c.123delA")
print(result["formatted"], result["changed"])
print(format_hgvs_variant("NM_004006.2:c.4375C>T")["changed"])
```

- **示例输出：**

```text
NM_004006.2:c.123del True
False
```

## 5) `dnakit.variants.check_hgvs_syntax()`

- **作用：** 批量检查 HGVS 字符串的语法，逐条报告是否有效、标准写法或错误信息。
- **处理方法：** 对每条字符串使用与 `parse_hgvs_variant()` 相同的解析流程。语法错误和反向区间都记为无效，不中断整批。非字符串条目会报错。
- **API：** `dnakit.variants.check_hgvs_syntax(texts, *, max_variants=10000)`。
- **输入：** 字符串列表；长度不能超过 `max_variants`（1 到 10000）。
- **输出：** `count`、`valid`、`invalid` 与 `items`。每项含 `input`、`valid`、`formatted` 和 `error`。
- **示例代码：**

```python
from dnakit.variants import check_hgvs_syntax

result = check_hgvs_syntax(["NM_004006.2:c.4375C>T", "NM_004006.2:c.123_122del", "not an hgvs string"])
print(result["count"], result["valid"], result["invalid"])
print([(item["input"], item["valid"]) for item in result["items"]])
```

- **示例输出：**

```text
3 1 2
[('NM_004006.2:c.4375C>T', True), ('NM_004006.2:c.123_122del', False), ('not an hgvs string', False)]
```

## 6) `dnakit.variants.read_plink_bed()`

- **作用：** 读取 PLINK 1 的 `.bed` 基因型矩阵，可按样本 ID 和变异 ID 选取子集。
- **处理方法：** `bed_reader.open_bed` 读取同名 `.fam` 和 `.bim`，再按索引取出基因型。`count_allele="allele_1"` 计数第一个等位基因（A1），与 bed-reader 的 `count_A1=True` 一致。缺失基因型返回 `None`。读取前检查 `.fam`、`.bim` 的记录数和矩阵单元数是否超过上限。
- **API：** `dnakit.variants.read_plink_bed(path, *, samples=None, variants=None, count_allele="allele_1", max_records=100000, max_cells=5000000)`。
- **输入：** `.bed` 文件路径；同目录需要同名的 `.bim` 与 `.fam`。`samples` 与 `variants` 为 ID 列表或 None。
- **输出：** `samples`、`variants`（每项含 `id`、`chromosome`、`position`、`allele_1`、`allele_2`）、`genotypes`（行为样本、列为变异，取值为 0、1、2 或 None）、`count_allele` 与 `provenance`。
- **示例代码：**

```python
from dnakit.variants import read_plink_bed, write_plink_bed

write_plink_bed(
    "panel.bed",
    [[0, 1, 2], [2, None, 1], [1, 1, 0]],
    sample_ids=["S1", "S2", "S3"],
    variant_ids=["rs1", "rs2", "rs3"],
    chromosomes=["1", "1", "2"],
    positions=[100, 200, 300],
    allele_1=["A", "C", "G"],
    allele_2=["G", "T", "A"],
    overwrite=True,
)
result = read_plink_bed("panel.bed")
print(result["samples"], [variant["id"] for variant in result["variants"]])
print(result["genotypes"])
```

- **示例输出：**

```text
['S1', 'S2', 'S3'] ['rs1', 'rs2', 'rs3']
[[0, 1, 2], [2, None, 1], [1, 1, 0]]
```

## 7) `dnakit.variants.write_plink_bed()`

- **作用：** 把基因型矩阵写成 PLINK 1 的 `.bed`、`.bim` 和 `.fam` 三个文件。
- **处理方法：** 调用 `bed_reader.to_bed` 写入。`count_allele` 说明矩阵中的值计数哪个等位基因。默认不覆盖已有文件，需要 `overwrite=True`。
- **API：** `dnakit.variants.write_plink_bed(path, genotypes, *, sample_ids, variant_ids, chromosomes, positions, allele_1, allele_2, count_allele="allele_1", overwrite=False, max_cells=5000000)`。
- **输入：** `genotypes` 为样本行、变异列的矩阵，取值为 0、1、2 或 None。其余列表的长度与矩阵维度一致；样本 ID 与变异 ID 必须唯一。
- **输出：** `path`、`files`、`samples`、`variants`、`count_allele` 与 `provenance`。
- **示例代码：**

```python
from pathlib import Path

from dnakit.variants import write_plink_bed

result = write_plink_bed(
    "written.bed",
    [[0, 1], [2, None]],
    sample_ids=["S1", "S2"],
    variant_ids=["rs1", "rs2"],
    chromosomes=["1", "1"],
    positions=[10, 20],
    allele_1=["A", "C"],
    allele_2=["G", "T"],
    overwrite=True,
)
print(result["operation"], result["samples"], result["variants"], [Path(name).name for name in result["files"]])
```

- **示例输出：**

```text
write_plink_bed 2 2 ['written.bed', 'written.bim', 'written.fam']
```

## 8) `dnakit.variants.summarize_plink_bed()`

- **作用：** 统计每个 PLINK 变异的有效基因型数、缺失数和被计数等位基因的频率。
- **处理方法：** 读取完整矩阵（受 `max_cells` 限制）。频率按 `被计数等位基因的拷贝数 / (2 × 有效基因型数)` 计算，因为 PLINK 基因型是二倍体。
- **API：** `dnakit.variants.summarize_plink_bed(path, *, count_allele="allele_1", max_records=100000, max_cells=5000000)`。
- **输入：** `.bed` 文件路径；同 `read_plink_bed()`。
- **输出：** `samples`、`variants`（每项含 `id`、`called`、`missing`、`allele_frequency`）、`missing_genotypes`、`count_allele` 与 `provenance`。
- **示例代码：**

```python
from dnakit.variants import summarize_plink_bed, write_plink_bed

write_plink_bed(
    "summary.bed",
    [[0, 2, 1], [0, 2, None], [0, None, 2]],
    sample_ids=["S1", "S2", "S3"],
    variant_ids=["rs1", "rs2", "rs3"],
    chromosomes=["1", "1", "2"],
    positions=[100, 200, 300],
    allele_1=["A", "C", "G"],
    allele_2=["G", "T", "A"],
    overwrite=True,
)
result = summarize_plink_bed("summary.bed")
print([(item["id"], item["called"], item["missing"], item["allele_frequency"]) for item in result["variants"]])
print(result["missing_genotypes"])
```

- **示例输出：**

```text
[('rs1', 3, 0, 0.0), ('rs2', 2, 1, 1.0), ('rs3', 2, 1, 0.75)]
2
```

## 9) `dnakit.variants.read_plink_metadata()`

- **作用：** 读取 `.fam` 样本表和 `.bim` 变异表，不读取基因型。
- **处理方法：** `pandas_plink.read_plink` 读取两张表，基因型部分不访问。读取前检查记录数上限。
- **API：** `dnakit.variants.read_plink_metadata(path, *, max_records=100000)`。
- **输入：** `.bed` 文件路径（同目录需要同名的 `.bim` 与 `.fam`）。
- **输出：** `samples`（含 `family_id`、`sample_id`、`father_id`、`mother_id`、`sex`、`phenotype`）、`variants`（含 `chromosome`、`id`、`genetic_position`、`position`、`allele_1`、`allele_2`）与 `provenance`。
- **示例代码：**

```python
from dnakit.variants import read_plink_metadata

result = read_plink_metadata("panel.bed")
print(result["samples"][0]["sample_id"], len(result["samples"]))
print(result["variants"][2])
```

- **示例输出：**

```text
S1 3
{'chromosome': '2', 'id': 'rs3', 'genetic_position': 0.0, 'position': 300, 'allele_1': 'G', 'allele_2': 'A'}
```

## 10) `dnakit.variants.read_plink_genotypes()`

- **作用：** 用 pandas-plink 读取 PLINK 1 基因型矩阵，返回与 `read_plink_bed()` 相同的结构。
- **处理方法：** `pandas_plink.read_plink1_bin` 按 `ref` 返回被计数的等位基因：`allele_1` 对应 `ref="a0"`，`allele_2` 对应 `ref="a1"`。读取前检查单元数上限。本函数不支持选取子集。
- **API：** `dnakit.variants.read_plink_genotypes(path, *, count_allele="allele_1", max_records=100000, max_cells=5000000)`。
- **输入：** `.bed` 文件路径；同目录需要同名的 `.bim` 与 `.fam`。
- **输出：** `samples`、`variants`、`genotypes`、`count_allele` 与 `provenance`。
- **示例代码：**

```python
from dnakit.variants import read_plink_genotypes

print(read_plink_genotypes("panel.bed")["genotypes"])
print(read_plink_genotypes("panel.bed", count_allele="allele_2")["genotypes"])
```

- **示例输出：**

```text
[[0, 1, 2], [2, None, 1], [1, 1, 0]]
[[2, 1, 0], [0, None, 1], [1, 1, 2]]
```

## 11) `dnakit.variants.write_plink_genotypes()`

- **作用：** 用 pandas-plink 的写入器写出 PLINK 1 的三个文件。
- **处理方法：** 把矩阵包装为 xarray `DataArray`，交给 `pandas_plink.write_plink1_bin`。该写入器保存第二个等位基因（a1）的计数，所以 `count_allele="allele_1"` 时先把每个值转换为 `2 - 值`。
- **API：** `dnakit.variants.write_plink_genotypes(path, genotypes, *, sample_ids, variant_ids, chromosomes, positions, allele_1, allele_2, count_allele="allele_1", overwrite=False, max_cells=5000000)`。
- **输入：** 同 `write_plink_bed()`。
- **输出：** 同 `write_plink_bed()`：`path`、`files`、`samples`、`variants`、`count_allele` 与 `provenance`。
- **示例代码：**

```python
from dnakit.variants import read_plink_bed, write_plink_genotypes

write_plink_genotypes(
    "pandas.bed",
    [[0, 1, 2], [2, None, 1], [1, 1, 0]],
    sample_ids=["S1", "S2", "S3"],
    variant_ids=["rs1", "rs2", "rs3"],
    chromosomes=["1", "1", "2"],
    positions=[100, 200, 300],
    allele_1=["A", "C", "G"],
    allele_2=["G", "T", "A"],
    overwrite=True,
)
print(read_plink_bed("pandas.bed")["genotypes"])
```

- **示例输出：**

```text
[[0, 1, 2], [2, None, 1], [1, 1, 0]]
```

## 12) `dnakit.variants.read_bgen_variants()`

- **作用：** 读取 BGEN 文件的样本数和变异元数据，不读取概率。
- **处理方法：** `bgen_reader.open_bgen` 打开文件，读取变异 ID、rsID、染色体、位置、等位基因、相位和基因型组合数。首次打开时，bgen-reader 会在 BGEN 文件旁生成 `.metadata2.mmm` 缓存文件。`max_variants` 限制变异数。本函数只在 zlib 压缩的 BGEN 1.2（布局 2）上测试过，未压缩文件未测试。
- **API：** `dnakit.variants.read_bgen_variants(path, *, max_variants=100000)`。
- **输入：** BGEN 1.2 文件路径。
- **输出：** `sample_count`、`variants`（每项含 `id`、`rsid`、`chromosome`、`position`、`alleles`、`phased`、`genotype_combinations`）与 `provenance`。
- **示例代码：**

```python
import struct
import zlib

from dnakit.variants import read_bgen_variants


def write_bgen(path, samples, variants):
    n = len(samples)
    flags = 0x80000000 | (2 << 2) | 1
    header = struct.pack("<III", 20, len(variants), n) + b"bgen" + struct.pack("<I", flags)
    names = b"".join(struct.pack("<H", len(s)) + s.encode() for s in samples)
    body = b""
    for ident, rsid, chrom, pos, alleles, probs in variants:
        for text in (ident, rsid, chrom):
            body += struct.pack("<H", len(text)) + text.encode()
        body += struct.pack("<IH", pos, len(alleles))
        body += b"".join(struct.pack("<I", len(a)) + a.encode() for a in alleles)
        raw = struct.pack("<IHBB", n, 2, 2, 2) + bytes([2] * n) + bytes([0, 8])
        raw += bytes(value for row in probs for value in row)
        packed = zlib.compress(raw)
        body += struct.pack("<II", 4 + len(packed), len(raw)) + packed
    sample_block = struct.pack("<II", 8 + len(names), n) + names
    with open(path, "wb") as handle:
        handle.write(struct.pack("<I", len(header) + len(sample_block)) + header + sample_block + body)


write_bgen(
    "two.bgen",
    ["S1", "S2", "S3"],
    [
        ("v1", "rs1", "1", 100, ["A", "G"], [[255, 0], [0, 255], [128, 127]]),
        ("v2", "rs2", "1", 200, ["C", "T"], [[0, 0], [255, 0], [0, 0]]),
    ],
)
result = read_bgen_variants("two.bgen")
print(result["sample_count"], [(v["rsid"], v["position"], v["alleles"]) for v in result["variants"]])
```

- **示例输出：**

```text
3 [('rs1', 100, ['A', 'G']), ('rs2', 200, ['C', 'T'])]
```

## 13) `dnakit.variants.read_bgen_probabilities()`

- **作用：** 读取指定变异的基因型概率。
- **处理方法：** `open_bgen(...).read(..., return_missings=True)` 返回概率与缺失标记，缺失样本返回 `None`。输出按变异、样本、基因型组合排列。对二倍体双等位位点，组合顺序为 AA、AB、BB。
- **API：** `dnakit.variants.read_bgen_probabilities(path, *, variant_indices=None, max_cells=5000000, max_variants=100000)`。
- **输入：** `variant_indices` 为变异的 0 基索引列表；为 None 时读取全部变异。
- **输出：** `samples`、`variants`、`probabilities`（变异 → 样本 → 概率列表，缺失为 None）与 `provenance`。
- **示例代码：**

```python
from dnakit.variants import read_bgen_probabilities

result = read_bgen_probabilities("two.bgen", variant_indices=[1])
print(result["variants"][0]["rsid"], result["samples"])
print(result["probabilities"][0])
```

- **示例输出：**

```text
rs2 ['S1', 'S2', 'S3']
[[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
```

## 14) `dnakit.variants.calculate_bgen_allele_frequency()`

- **作用：** 计算未分相 BGEN 变异的等位基因频率。
- **处理方法：** `allele_expectation` 给出每个样本的期望等位基因拷贝数（由概率得到）。剔除缺失样本后，每个等位基因的频率为其期望拷贝数除以总拷贝数。bgen-reader 自带的 `allele_frequency` 返回的是拷贝数之和除以等位基因数（结果之和等于样本数），因此本函数不直接使用它。分相变异会报错。
- **API：** `dnakit.variants.calculate_bgen_allele_frequency(path, *, variant_indices=None, max_cells=5000000, max_variants=100000)`。
- **输入：** 同 `read_bgen_probabilities()`。
- **输出：** `variants`（每项含 `id`、`rsid`、`chromosome`、`position`、`alleles`、`called_samples`、`allele_copies`、`allele_frequencies`）与 `provenance`。
- **示例代码：**

```python
from dnakit.variants import calculate_bgen_allele_frequency

result = calculate_bgen_allele_frequency("two.bgen")
print([(v["rsid"], [round(f, 4) for f in v["allele_frequencies"]]) for v in result["variants"]])
```

- **示例输出：**

```text
[('rs1', [0.7503, 0.2497]), ('rs2', [0.3333, 0.6667])]
```

## 验证方式

- 测试文件是 `others/tests/unit/test_variant_formats_backends.py`，共 32 项。大部分测试在本机真实安装的包上运行；缺包测试用 `sys.modules` 模拟包缺失。
- PLINK：bed-reader 与 pandas-plink 的读写结果相互核对，两种计数约定都检查；与直接调用 `bed_reader.open_bed` 和 `pandas_plink.read_plink` 的结果比较；汇总频率与手工计算比较。
- BGEN：测试文件按 bgen-reader 样例文件的字节结构，写出 zlib 压缩的 BGEN 1.2 布局 2 文件，并由 bgen-reader 读回。概率、缺失样本和等位基因频率与手工计算比较，频率还与 `allele_expectation` 的直接调用结果比较。
- varcode：`describe_variant()` 的修剪结果与 `varcode.Variant` 比较；`load_vcf_variants()` 与 `varcode.load_vcf` 比较。
- hgvs：解析和格式化结果与直接调用 `hgvs.parser.Parser()` 比较；语法错误、反向区间、非字符串条目和超出上限的批量都会被拒绝。
- 缺少可选依赖时，每个模块都会给出对应的安装提示（`bio-varcode`、`bio-hgvs`、`bio-bed-reader`、`bio-pandas-plink`、`bio-bgen-reader`）。

## 当前范围

- 已包装：varcode 的最小表示、变异类型和 VCF 读取；hgvs 的语法解析、格式化和批量语法检查；bed-reader 与 pandas-plink 的 PLINK 1 读写和汇总；bgen-reader 的元数据、概率和等位基因频率。
- 未包装 varcode 的基因功能效应（`Variant.effects()` 等）：需要 pyensembl 下载 Ensembl 参考数据，本环境未运行。
- 未包装 hgvs 的参考数据函数（序列获取、UTA 投影、validator、normalizer）：需要 seqrepo 和 psycopg2 等参考数据与数据库，本环境缺少这些依赖，因此这些函数没有运行，也没有测试。
- 未包装 pandas-plink 的 `read_grm` 和 `read_rel`：它们读取亲缘关系矩阵，不是 DNA 或基因层面的操作。
- 未测试 bgen-reader 的未压缩 BGEN 文件和 `allow_complex=True` 的文件。
- 未包装 sgkit：当前虚拟环境中无法导入，因此没有 sgkit 相关函数。
