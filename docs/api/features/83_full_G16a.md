# DnaChisel、DnaCauldron 与 pydna 序列函数

这一组函数是 `dnakit.molecular_lab` 中基于 DnaChisel、DnaCauldron 和 pydna 的模块级序列函数。每个函数只在调用时导入对应的包；缺少包时报 `BackendUnavailableError`，并给出 `dnakit[bio-dnachisel]`、`dnakit[bio-dnacauldron]` 或 `dnakit[bio-assembly]` 的安装提示。

函数只处理显式给出的序列或文件。`fetch_genbank_record()` 需要网络，示例用替身函数代替下载。下面的示例代码从模块文件导入；接入 `dnakit.molecular_lab` 导出后，也可以直接从包中导入。已有封装的函数（如 `reverse_translate_protein()`、`design_pcr_primers()`、`find_recombination_sites()` 和 `optimize_dnachisel_record()`）不在这里重复。

## 1) `dnakit.molecular_lab.generate_constrained_dna_sequence()`

- **作用：** 生成满足 DnaChisel 约束的随机 DNA 序列。
- **处理方法：** 调用 DnaChisel 的 `random_compatible_dna_sequence` 生成并反复尝试，直到约束满足；再用 `DnaOptimizationProblem` 逐项检查。约束名只能取 `AvoidPattern`、`AvoidStopCodons`、`AvoidHairpins`、`EnforceGCContent`。给定 `seed` 时结果可复现。
- **API：** `dnakit.molecular_lab.generate_constrained_dna_sequence(sequence_length, *, constraints=(), seed=None, max_random_iters=5000, verbose=True)`
- **输入：** 长度为 1 到 20000 的整数；`constraints` 是映射列表，每项含 `name` 和该约束的参数，例如 `{"name": "AvoidPattern", "pattern": "GGTCTC"}`。
- **输出：** `sequence`、`sequence_length`、`seed`、`constraint_names`、`all_constraints_pass` 和 `provenance`。
- **示例代码：**

```python
from dnakit.molecular_lab import generate_constrained_dna_sequence

result = generate_constrained_dna_sequence(
    30,
    constraints=[{"name": "AvoidPattern", "pattern": "GGTCTC"}],
    seed=1,
)
print(result["sequence"])
print(result["all_constraints_pass"])
```

- **示例输出：**

```text
TGAAGTGTGAATAGTACTCACTCAGACATC
True
```

## 2) `dnakit.molecular_lab.find_sequence_difference_segments()`

- **作用：** 找出两条等长 DNA 序列中不同的连续区段。
- **处理方法：** 调用 DnaChisel 的 `sequences_differences_segments` 逐位比较，相邻的不同位置合并为一段。函数不做比对，也不处理插入和缺失。
- **API：** `dnakit.molecular_lab.find_sequence_difference_segments(sequence_a, sequence_b, *, verbose=True)`
- **输入：** 两条只含 A、C、G、T 的等长序列；空白被去掉，大小写不敏感。
- **输出：** `segments`（`[start, end)` 列表）、`segment_count`、`differing_positions`、`sequence_length` 和 `provenance`。
- **示例代码：**

```python
from dnakit.molecular_lab import find_sequence_difference_segments

result = find_sequence_difference_segments("ACGTTACGTA", "ACGAAACGTC")
print(result["segments"])
print(result["differing_positions"])
```

- **示例输出：**

```text
[[3, 5], [9, 10]]
3
```

## 3) `dnakit.molecular_lab.list_common_restriction_enzymes()`

- **作用：** 列出 DnaChisel 内置酶表中识别位点长度相同、供应商数量足够的限制性内切酶。
- **处理方法：** 调用 DnaChisel 的 `list_common_enzymes`，按 `site_length` 和 `min_suppliers` 过滤，结果按名称排序。不访问网络。
- **API：** `dnakit.molecular_lab.list_common_restriction_enzymes(*, site_length=6, min_suppliers=1, verbose=True)`
- **输入：** `site_length` 为 4 到 12 的整数；`min_suppliers` 为 0 到 20 的整数。
- **输出：** `enzyme_count`、排序后的 `enzymes` 列表，以及 `site_length`、`min_suppliers` 和 `provenance`。
- **示例代码：**

```python
from dnakit.molecular_lab import list_common_restriction_enzymes

result = list_common_restriction_enzymes(site_length=6, min_suppliers=5)
print(result["enzyme_count"])
print(result["enzymes"][:5])
```

- **示例输出：**

```text
37
['AatII', 'AccI', 'ApaI', 'BalI', 'BamHI']
```

## 4) `dnakit.molecular_lab.group_nearby_positions()`

- **作用：** 把位置列表按距离分组，用于聚类相邻的变异或突变位点。
- **处理方法：** 调用 DnaChisel 的 `group_nearby_indices`。相邻位置的差值小于 `max_gap` 时并入同一组；`max_group_spread` 限制一组内首末位置之差。两个条件都可省略（为 `None` 时不限制），结果按位置排序。
- **API：** `dnakit.molecular_lab.group_nearby_positions(positions, *, max_gap=None, max_group_spread=None, verbose=True)`
- **输入：** 不重复的非负整数列表，最多 100000 个；`max_gap` 与 `max_group_spread` 为非负整数或 `None`。
- **输出：** `groups`（每组已排序）、`group_count`、`position_count`、`max_gap`、`max_group_spread` 和 `provenance`。
- **示例代码：**

```python
from dnakit.molecular_lab import group_nearby_positions

result = group_nearby_positions([30, 1, 11, 2, 10, 3], max_gap=2)
print(result["groups"])
```

- **示例输出：**

```text
[[1, 2, 3], [10, 11], [30]]
```

## 5) `dnakit.molecular_lab.autoselect_golden_gate_enzyme()`

- **作用：** 判断一组 Golden Gate 元件应当使用哪一种 Type IIS 酶。
- **处理方法：** 把每条序列做成线性记录，交给 DnaCauldron 的 `autoselect_enzyme`。它从候选酶中选出切割结构最接近两个位点的一种。若没有候选酶符合结构，DnaCauldron 仍返回一个最佳候选，调用者应自行检查切点。
- **API：** `dnakit.molecular_lab.autoselect_golden_gate_enzyme(parts, *, enzymes=("BsmBI", "BsaI", "BbsI", "AarI", "SapI"), verbose=True)`
- **输入：** 1 到 200 条 DNA 字符串，每条两端带酶切位点；`enzymes` 是候选酶名的子集，不能重复。
- **输出：** `enzyme`、`candidate_enzymes`、`part_count` 和 `provenance`。
- **示例代码：**

```python
parts = [
    "GGTCTCAATGCTTTAAACCCGGGTTTTGCTTGGAGACC",
    "GGTCTCAGCTTAAAACCCCGGGGAAAACGTAGGAGACC",
    "GGTCTCACGTACCCCAAAATTTTGGGGATGCGGAGACC",
]

from dnakit.molecular_lab import autoselect_golden_gate_enzyme

print(autoselect_golden_gate_enzyme(parts)["enzyme"])
```

- **示例输出：**

```text
BsaI
```

## 6) `dnakit.molecular_lab.list_type2s_digestion_overhangs()`

- **作用：** 列出 Type IIS 酶切后形成的粘性末端（overhang）。
- **处理方法：** 调用 DnaCauldron 的 `list_digestion_overhangs`。结果去重并排序；DnaCauldron 对缺失末端返回的占位字符串 `"None"` 会被去掉并单独计数。`parts_only` 为真时只保留元件片段上的 overhang。
- **API：** `dnakit.molecular_lab.list_type2s_digestion_overhangs(sequences, *, enzyme="auto", parts_only=True, verbose=True)`
- **输入：** 1 到 200 条 DNA 字符串；`enzyme` 为 `"auto"` 或 `BsmBI`、`BsaI`、`BbsI`、`AarI`、`SapI` 之一。
- **输出：** `overhangs`、`overhang_count`、`none_placeholder_count`、`enzyme`、`parts_only`、`sequence_count` 和 `provenance`。
- **示例代码：**

```python
parts = [
    "GGTCTCAATGCTTTAAACCCGGGTTTTGCTTGGAGACC",
    "GGTCTCAGCTTAAAACCCCGGGGAAAACGTAGGAGACC",
    "GGTCTCACGTACCCCAAAATTTTGGGGATGCGGAGACC",
]

from dnakit.molecular_lab import list_type2s_digestion_overhangs

result = list_type2s_digestion_overhangs(parts, enzyme="BsaI")
print(result["overhangs"])
print(result["overhang_count"])
```

- **示例输出：**

```text
['ATGC', 'CGTA', 'GCTT']
3
```

## 7) `dnakit.molecular_lab.crop_dna_record_with_features()`

- **作用：** 截取一段 DNA，并保留与截取窗口有重叠的注释。
- **处理方法：** DnaCauldron 的 `crop_record_with_saddling_features` 决定截取的序列和保留的注释。DnaCauldron 2.0.12 对部分落在窗口内的注释保留绝对坐标，因此本函数按窗口裁剪并平移这些注释的坐标；窗口外的注释被丢弃。
- **API：** `dnakit.molecular_lab.crop_dna_record_with_features(sequence, *, start, end, features=(), verbose=True)`
- **输入：** DNA 字符串；`start` 和 `end` 是从 0 开始的半开区间；`features` 是映射列表，含 `type`、`start`、`end`，可选 `label`。
- **输出：** `sequence`、`sequence_length`、`features`（坐标相对于新序列）、`feature_count`、`source_length` 和 `provenance`。
- **示例代码：**

```python
from dnakit.molecular_lab import crop_dna_record_with_features

result = crop_dna_record_with_features(
    "ACGT" * 50,
    start=10,
    end=120,
    features=[
        {"type": "CDS", "start": 0, "end": 40, "label": "a"},
        {"type": "misc_feature", "start": 60, "end": 100},
        {"type": "misc_feature", "start": 150, "end": 180},
    ],
)
print(result["sequence_length"])
print(result["features"])
```

- **示例输出：**

```text
110
[{'type': 'CDS', 'start': 0, 'end': 30, 'label': 'a'}, {'type': 'misc_feature', 'start': 50, 'end': 90, 'label': ''}]
```

## 8) `dnakit.molecular_lab.simulate_pcr_amplicon()`

- **作用：** 用一对引物模拟线性模板的 PCR 产物。
- **处理方法：** 调用 pydna 的 `amplify.pcr`。每条引物至少要有 `min_annealing` 个碱基与模板配对；产物从正向引物的 5' 端延伸到反向引物的 5' 端。找不到产物时报错。
- **API：** `dnakit.molecular_lab.simulate_pcr_amplicon(template, forward_primer, reverse_primer, *, min_annealing=14, verbose=True)`
- **输入：** 模板 DNA（最长 20000 bp）；两条引物序列（各最长 200 nt）；`min_annealing` 为 8 到 60 的整数。
- **输出：** `product_sequence`、`product_length`、`template_length`、`forward_primer`、`reverse_primer` 和 `provenance`。
- **示例代码：**

```python
from dnakit.molecular_lab import simulate_pcr_amplicon

template = "AAAAA" + "ACGTTGCAAGGCTTAGCAT" + "C" * 16 + "GCAAGGCTTAGCATTGGACC" + "TTTTT"
result = simulate_pcr_amplicon(template, "ACGTTGCAAGGCTTAGCAT", "GGTCCAATGCTAAGCCTTGC")
print(result["product_length"])
print(result["product_sequence"])
```

- **示例输出：**

```text
55
ACGTTGCAAGGCTTAGCATCCCCCCCCCCCCCCCCGCAAGGCTTAGCATTGGACC
```

## 9) `dnakit.molecular_lab.assemble_pydna_fragments()`

- **作用：** 枚举若干 DNA 片段可以组装成的产物。
- **处理方法：** 调用 pydna 的 `assembly2` 算法。重叠法（`gibson`、`in_fusion`、`fusion_pcr`、`in_vivo`）按 `min_overlap` 查找端部重叠；酶切法（`golden_gate`、`restriction_ligation`）用给定的 `Bio.Restriction` 酶切后连接。产物按长度和序列排序，最多列出 `max_products` 个。
- **API：** `dnakit.molecular_lab.assemble_pydna_fragments(fragments, *, method="gibson", min_overlap=25, enzymes=(), circular_only=False, max_products=50, verbose=True)`
- **输入：** 2 到 20 条 DNA 片段，总长不超过 200000 bp；酶切法需要 `enzymes`，例如 `["BsaI"]`；重叠法需要省略 `enzymes`。
- **输出：** `product_count`、`products`（每个产物含序列、长度和是否环状）、`method`、`min_overlap`、`enzymes` 和 `provenance`。
- **示例代码：**

```python
parts = [
    "GGTCTCAATGCTTTAAACCCGGGTTTTGCTTGGAGACC",
    "GGTCTCAGCTTAAAACCCCGGGGAAAACGTAGGAGACC",
    "GGTCTCACGTACCCCAAAATTTTGGGGATGCGGAGACC",
]

from dnakit.molecular_lab import assemble_pydna_fragments

result = assemble_pydna_fragments(
    parts, method="golden_gate", enzymes=["BsaI"], circular_only=True
)
print(result["product_count"], result["products"][0]["length"], result["products"][0]["circular"])
print(result["products"][0]["sequence"])
```

- **示例输出：**

```text
1 60 True
ATGCTTTAAACCCGGGTTTTGCTTAAAACCCCGGGGAAAACGTACCCCAAAATTTTGGGG
```

## 10) `dnakit.molecular_lab.find_common_substrings()`

- **作用：** 找出两条序列之间所有长度不小于 `min_length` 的公共子串。
- **处理方法：** 调用 pydna 的 `common_sub_strings`，每个匹配为 `(x_start, y_start, length)`，按位置排序，最多列出 1000 条。
- **API：** `dnakit.molecular_lab.find_common_substrings(sequence_x, sequence_y, *, min_length=25, verbose=True)`
- **输入：** 两条 DNA 序列（各最长 20000 bp）；`min_length` 为 1 到 20000 的整数。
- **输出：** `matches`、`match_count`、`min_length` 和 `provenance`。
- **示例代码：**

```python
from dnakit.molecular_lab import find_common_substrings

result = find_common_substrings("ACGTTGCAAGGCTTAGCAT", "TTGCAAGGCTTAGCATCC", min_length=5)
print(result["matches"])
```

- **示例输出：**

```text
[{'x_start': 3, 'y_start': 0, 'length': 16}]
```

## 11) `dnakit.molecular_lab.find_terminal_overlaps()`

- **作用：** 找出 `sequence_x` 末端与 `sequence_y` 起始端的重叠，用于 Gibson 类拼接。
- **处理方法：** 调用 pydna 的 `common_sub_strings.terminal_overlap`。返回的匹配应当贴合 `sequence_x` 的 3' 末端和 `sequence_y` 的 5' 起始端。
- **API：** `dnakit.molecular_lab.find_terminal_overlaps(sequence_x, sequence_y, *, min_length=15, verbose=True)`
- **输入：** 两条 DNA 序列（各最长 20000 bp）；`min_length` 为 1 到 20000 的整数。
- **输出：** `matches`、`match_count`、`min_length` 和 `provenance`。
- **示例代码：**

```python
from dnakit.molecular_lab import find_terminal_overlaps

result = find_terminal_overlaps("ACGTTGCAAGGCTTAGCAT", "TTGCAAGGCTTAGCATCC", min_length=5)
print(result["match_count"])
print(result["matches"])
```

- **示例输出：**

```text
1
[{'x_start': 3, 'y_start': 0, 'length': 16}]
```

## 12) `dnakit.molecular_lab.primer_pair_hybridization()`

- **作用：** 计算两条引物之间可能形成的双链产物，用于检查引物二聚体。
- **处理方法：** 调用 pydna 的 `oligonucleotide_hybridization` 和 `oligonucleotide_hybridization_overhangs`。两条引物在共享的互补区至少有 `minimal_annealing` 个碱基时配对；存在错配的配对会报错。
- **API：** `dnakit.molecular_lab.primer_pair_hybridization(forward_primer, reverse_primer, *, minimal_annealing=10, verbose=True)`
- **输入：** 两条 5'→3' 方向的引物序列（各最长 200 nt）；`minimal_annealing` 为 4 到 60 的整数。
- **输出：** `overhangs`（pydna 返回的可能重叠长度，符号沿用 pydna 的约定）、`product_count`、`products` 和 `provenance`。
- **示例代码：**

```python
from dnakit.molecular_lab import primer_pair_hybridization

result = primer_pair_hybridization("GGGGGACGTTGCAAGGCTTAGCATG", "CATGCTAAGCCTTGCAACGTCCCCC")
print(result["overhangs"])
print(result["products"])
```

- **示例输出：**

```text
[0]
[{'sequence': 'GGGGGACGTTGCAAGGCTTAGCATG', 'length': 25}]
```

## 13) `dnakit.molecular_lab.primer_melting_temperature()`

- **作用：** 用 pydna 默认的最近邻参数计算引物解链温度。
- **处理方法：** 调用 pydna 的 `tm.tm_default`，使用严格模式，拒绝非 ACGT 字符。结果保留三位小数。
- **API：** `dnakit.molecular_lab.primer_melting_temperature(sequence, *, verbose=True)`
- **输入：** 6 到 200 个碱基的 DNA 字符串。
- **输出：** `tm_celsius`、`length`、`sequence` 和 `provenance`。
- **示例代码：**

```python
from dnakit.molecular_lab import primer_melting_temperature

result = primer_melting_temperature("ACGTTGCAAGGCTTAGCATG")
print(result["tm_celsius"])
```

- **示例输出：**

```text
62.197
```

## 14) `dnakit.molecular_lab.parse_primer_sequences()`

- **作用：** 读取 FASTA 格式或纯文本中的引物序列。
- **处理方法：** 调用 pydna 的 `parsers.parse_primers`。`source` 可以是已存在的文件路径，也可以直接是文本；名称取自 FASTA 标题行，顺序与输入一致。
- **API：** `dnakit.molecular_lab.parse_primer_sequences(source, *, verbose=True)`
- **输入：** 文件路径或文本；最多 1000 条引物，文件或文本不超过 1 MB。
- **输出：** `primers`（每条含名称、序列和长度）、`primer_count`、`source` 和 `provenance`。
- **示例文件：** `primers.fa`

```text
>p1
ACGTACGTACGT
>p2
GGGGCCCCAAAA
```
- **示例代码：**

```python
from dnakit.molecular_lab import parse_primer_sequences

result = parse_primer_sequences("primers.fa")
print(result["primer_count"])
print(result["primers"])
```

- **示例输出：**

```text
2
[{'name': 'p1', 'sequence': 'ACGTACGTACGT', 'length': 12}, {'name': 'p2', 'sequence': 'GGGGCCCCAAAA', 'length': 12}]
```

## 15) `dnakit.molecular_lab.fetch_genbank_record()`

- **作用：** 按 accession 下载一条 GenBank 核酸记录，用于基因或序列级查询。
- **处理方法：** 调用 pydna 的 `genbank.genbank`，它通过 NCBI Entrez 下载。accession 必须形如 `CS570233.1` 或 `NC_000913.3`。`email` 是 NCBI 要求的联系地址，只传给客户端，不写入结果。本函数需要网络；示例用替身函数代替下载，默认测试只用 monkeypatch 模拟。
- **API：** `dnakit.molecular_lab.fetch_genbank_record(accession, *, email, verbose=True)`
- **输入：** accession 字符串；`email` 为含 `@` 的联系地址。
- **输出：** `accession`、`name`、`description`、`topology`、`length`、`sequence` 和 `provenance`。
- **示例代码：**

```python
from pydna.dseqrecord import Dseqrecord
import pydna.genbank
from dnakit.molecular_lab import fetch_genbank_record


def offline_stand_in(accession, *args, email=None, **kwargs):
    return Dseqrecord("ACGTACGTAA", name="DEMO")


# Replaces the NCBI download for this offline example only.
pydna.genbank.genbank = offline_stand_in
result = fetch_genbank_record("NC_000913.3", email="reader@example.org")
print(result["accession"], result["length"], result["topology"])
print(result["sequence"])
```

- **示例输出：**

```text
NC_000913.3 10 linear
ACGTACGTAA
```
