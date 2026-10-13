# H2 功能：pydna、DnaCauldron、DnaChisel 与 primer3 的 family-variant 核查

本页记录对 pydna、DnaCauldron、DnaChisel 和 primer3 共 101 个 `family-variant` 台账行的核查结果，以及由此新增的 10 个包装函数。它们位于 `dnakit.molecular_lab` 下的三个模块：`pydna_family_functions`、`dnacauldron_family_functions` 与 `dnachisel_family_functions`；primer3 没有新增函数，因此没有对应模块。

- 离线运行：示例只使用调用参数或示例文件，不访问网络。
- 核查方法：`duplicate-of` 用 `is` 比较两个上游对象，或比较同一段代码的结果；保留 `family-variant` 的行用已有包装函数的参数调用，并与直接调用上游函数的结果比较；其余行改为更准确的状态或新增包装。检查全部写在 `others/tests/unit/test_family_H2.py`。
- 示例从各自的模块导入，例如 `from dnakit.molecular_lab import simulate_pcr_assembly`。

## 台账状态（101 行）

| 包 | 原 family-variant | duplicate-of | family-variant（已核实） | wrapped | 其他状态 |
| --- | ---: | ---: | ---: | --- | --- |
| pydna | 59 | 34 | 13 | 7 行（6 个函数） | deprecated 3、broken-upstream 1、needs-network 1 |
| dnacauldron | 21 | 18 | 0 | 3 行（2 个函数） | 无 |
| dnachisel | 13 | 11 | 0 | 2 行（2 个函数） | 无 |
| primer3 | 8 | 8 | 0 | 0 | 无 |
| 合计 | 101 | 71 | 13 | 12 行（10 个函数） | 5 |

## 主要发现

- `pydna.readers.read` 等三行被标为 `parse_primer_sequences` 的变体，这是错的：`read` 返回 DNA 记录而不是引物，它就是 `parse` 加上“恰好一条记录”的检查，已由 `parse_sequence_text` 覆盖。
- `pydna.assembly2.common_sub_strings` 并不是字符串版本的别名：它接收 Dseqrecord 并能跨过环状分子的起点，例如示例中得到 `(4, 0, 11)` 而不是 `(4, 0, 7)`，所以新增了包装函数。
- `pydna.assembly2.pcr_assembly` 不同于 `amplify.pcr`：它允许引物错配、支持环状模板并可给出多个产物。
- `pydna.tm` 中的 `tm_dbd`、`tm_product`、`ta_default`、`ta_dbd`、`tmbresluc` 被挂在 `primer_melting_temperature` 名下，但该包装函数没有任何参数能到达它们，因此各自新增了包装。`tm_neb` 需要访问 NEB 网络服务（`needs-network`），`Q5` 的函数体只有 `raise NotImplementedError`（`broken-upstream`）。
- `homologous_recombination_excision`、`cre_lox_excision`、`recombinase_excision` 在 pydna 5.5.16 中已弃用（运行时给出 DeprecationWarning，建议改用 `*_excision_or_inversion`），改为 `deprecated`；`simulate_genome_edit(operation="excision")` 仍能到达它们。
- `pfu_sso7d_program` 只是 `dbd_program` 的别名，所以 `suggest_pcr_program` 的 `polymerase="dbd"` 与 `"pfu_sso7d"` 给出相同文本。
- 大多数其余行本来就是同一个对象的再导出（`snapgene_history_parser`、`all`、`design` 等模块），台账状态改为 `duplicate-of`。
- DnaCauldron 的 `load_records_from_files`、`load_records_from_file` 与 `string_to_records` 被说成已由 `read_dnacauldron_parts` 覆盖，其实后者按扩展名每个文件只读一条记录，不能读多记录文件、压缩包、目录或纯序列文本。
- `annotate_loxP_sites` 与 `annotate_gateway_sites` 只是把 `find_*_sites` 找到的位点写成特征，位置与 `find_recombination_sites` 的结果相同。
- primer3 的六个 `calc_*_tm` 就是对应 `calc_*(...).tm`；两个 `_b` 函数是 bytes 形式，`sanitize_sequence` 内部就调用 `sanitize_sequence_b`。

## 1) `find_common_substrings_with_topology`

- **作用：** 在两条 DNA 分子之间找共同子串，环状分子可跨过起点。
- **处理方法：** 调用 pydna.assembly2.common_sub_strings。环状分子被当作首尾相连处理，因此匹配可以从末端接到开头，长度不超过较短的分子；两条都是线状时结果与 `find_common_substrings` 相同。每个匹配为 `(x_start, y_start, length)`，起点从 0 开始，排序后最多列出 1000 行。
- **API：** `find_common_substrings_with_topology(sequence_x[必须], sequence_y[必须], min_length[可选], circular_x[可选], circular_y[可选], verbose[可选])`
- **输入：** 两条 A/C/G/T 序列、最小匹配长度（默认 25）以及两个布尔值 `circular_x`、`circular_y`。
- **输出：** 字典，含 `min_length`、`circular_x`、`circular_y`、`match_count`、`matches`（`x_start`、`y_start`、`length`）与 `provenance`。

**示例代码：**

```python
from dnakit.molecular_lab import find_common_substrings_with_topology

x = "ACGTTTGCAAC"
y = "TTGCAACACGTT"
circular = find_common_substrings_with_topology(x, y, min_length=6, circular_x=True, verbose=False)
linear = find_common_substrings_with_topology(x, y, min_length=6, verbose=False)
print(circular["matches"])
print(linear["matches"])
```

**示例输出：**

```text
[{'x_start': 4, 'y_start': 0, 'length': 11}]
[{'x_start': 4, 'y_start': 0, 'length': 7}]
```

## 2) `simulate_pcr_assembly`

- **作用：** 枚举一对引物的全部 PCR 产物，允许引物错配。
- **处理方法：** 调用 pydna.assembly2.pcr_assembly。每条引物至少要有 `min_annealing` 个碱基配对，`mismatches` 允许退火区内有若干错配，这是 `simulate_pcr_amplicon` 做不到的；环状模板可以得到跨过起点的产物，多个结合位点会给出多个产物。产物按长度与序列排序，没有产物时列表为空。
- **API：** `simulate_pcr_assembly(template[必须], forward_primer[必须], reverse_primer[必须], min_annealing[可选], mismatches[可选], template_circular[可选], max_products[可选], verbose[可选])`
- **输入：** 模板、正向引物、反向引物（A/C/G/T）、最小退火长度（默认 14）、允许错配数（0 到 5）、模板是否为环状与最多列出的产物数。
- **输出：** 字典，含 `template_length`、`template_circular`、`min_annealing`、`mismatches`、`product_count`、`products`（`sequence`、`length`、`circular`）与 `provenance`。

**示例代码：**

```python
from dnakit.molecular_lab import simulate_pcr_assembly

template = "CCGTAATGCCTTTCCCTAACAGAGTTTTTCGAACTCGTGTTGTCGAGCGACGGAATTAGATCAGTTAAATGGCAGAAAACTGGCAGGGCTTTTAGTCGTGGGATGATCAGTGGGTAAAGG"
forward = "AGAGTTTTTCGAACTCGTGTTG"
reverse = "CCCACGACTAAAAGCCCTGCCA"
result = simulate_pcr_assembly(template, forward, reverse, verbose=False)
print(result["product_count"], result["products"][0]["length"], result["products"][0]["sequence"] == template[20:102])
changed = forward[:10] + ("A" if forward[10] != "A" else "C") + forward[11:]
strict = simulate_pcr_assembly(template, changed, reverse, verbose=False)
tolerant = simulate_pcr_assembly(template, changed, reverse, mismatches=1, verbose=False)
print(strict["product_count"], tolerant["product_count"], [row["length"] for row in tolerant["products"]])
```

**示例输出：**

```text
1 82 True
0 1 [82]
```

## 3) `primer_melting_temperature_dbd`

- **作用：** 用适合 DNA 结合域聚合酶（如 Phusion、Pfu-Sso7d）的规则计算引物解链温度。
- **处理方法：** 调用 pydna.tm.tm_dbd：最近邻法，SantaLucia 1998 参数表，50 mM Na+、1.5 mM Mg2+，每条链 250 nM。它与使用 Taq 缓冲液的 `primer_melting_temperature` 结果不同。结果保留三位小数。
- **API：** `primer_melting_temperature_dbd(sequence[必须], verbose[可选])`
- **输入：** 长度 6 到 200 的 A/C/G/T 引物序列。
- **输出：** 字典，含 `sequence`、`length`、`tm_celsius` 与 `provenance`。

**示例代码：**

```python
from dnakit.molecular_lab import primer_melting_temperature
from dnakit.molecular_lab import primer_melting_temperature_dbd

primer = "AGAGTTTTTCGAACTCGTGTTG"
print(primer_melting_temperature(primer, verbose=False)["tm_celsius"])
print(primer_melting_temperature_dbd(primer, verbose=False)["tm_celsius"])
```

**示例输出：**

```text
59.445
55.982
```

## 4) `primer_melting_temperature_breslauer`

- **作用：** 用 Breslauer 参数计算引物解链温度。
- **处理方法：** 调用 pydna.tm.tmbresluc：用 Breslauer 最近邻焓和熵，并按引物浓度 `primer_nm`（nM，默认 500）与一价阳离子浓度 `salt_mm`（mM，默认 50）修正。盐浓度加倍使结果增加 16.6 乘 log10(2)，约 5.0 度。结果保留三位小数。
- **API：** `primer_melting_temperature_breslauer(sequence[必须], primer_nm[可选], salt_mm[可选], verbose[可选])`
- **输入：** 长度 6 到 200 的 A/C/G/T 引物序列、引物浓度与盐浓度（均为正数）。
- **输出：** 字典，含 `sequence`、`length`、`primer_nm`、`salt_mm`、`tm_celsius` 与 `provenance`。

**示例代码：**

```python
from dnakit.molecular_lab import primer_melting_temperature_breslauer

primer = "AGAGTTTTTCGAACTCGTGTTG"
base = primer_melting_temperature_breslauer(primer, verbose=False)
salty = primer_melting_temperature_breslauer(primer, salt_mm=100.0, verbose=False)
print(base["tm_celsius"], salty["tm_celsius"])
print(round(salty["tm_celsius"] - base["tm_celsius"], 2))
```

**示例输出：**

```text
59.082 64.079
5.0
```

## 5) `amplicon_melting_temperature`

- **作用：** 由长度与 GC 含量估算 PCR 产物的解链温度。
- **处理方法：** 调用 pydna.tm.tm_product，即 Rychlik 等人（1990）的公式：81.5 + 0.41 乘 GC 百分比 + 16.6 乘 log10(K) - 675 除以长度，K 为一价阳离子浓度 `potassium_m`（mol/L，默认 0.05）。结果保留三位小数。
- **API：** `amplicon_melting_temperature(sequence[必须], potassium_m[可选], verbose[可选])`
- **输入：** 至少 20 个碱基的 A/C/G/T 产物序列与阳离子浓度（正数）。
- **输出：** 字典，含 `length`、`gc_fraction`、`potassium_m`、`tm_celsius` 与 `provenance`。

**示例代码：**

```python
from dnakit.molecular_lab import amplicon_melting_temperature

product = "AGAGTTTTTCGAACTCGTGTTGTCGAGCGACGGAATTAGATCAGTTAAATGGCAGAAAACTGGCAGGGCTTTTAGTCGTGGG"
result = amplicon_melting_temperature(product, verbose=False)
print(result["length"], result["gc_fraction"], result["tm_celsius"])
print(amplicon_melting_temperature(product, potassium_m=0.5, verbose=False)["tm_celsius"])
```

**示例输出：**

```text
82 0.463415 70.671
87.271
```

## 6) `pcr_annealing_temperature`

- **作用：** 为一对引物及其产物建议 PCR 退火温度。
- **处理方法：** `polymerase="taq"`（默认）调用 pydna.tm.ta_default：0.3 乘较低的引物 Tm + 0.7 乘产物 Tm - 14.9；`polymerase="dbd"` 调用 pydna.tm.ta_dbd：较低的引物 Tm（DBD 规则）加 3 度。两条引物的 Tm 也一并给出，均保留三位小数。
- **API：** `pcr_annealing_temperature(forward_primer[必须], reverse_primer[必须], product[必须], polymerase[可选], verbose[可选])`
- **输入：** 两条引物的退火部分（6 到 200 个碱基）、至少 20 个碱基的产物序列与聚合酶类型。
- **输出：** 字典，含 `polymerase`、`forward_tm_celsius`、`reverse_tm_celsius`、`product_length`、`annealing_celsius` 与 `provenance`。

**示例代码：**

```python
from dnakit.molecular_lab import pcr_annealing_temperature

forward = "AGAGTTTTTCGAACTCGTGTTG"
reverse = "CCCACGACTAAAAGCCCTGCCA"
product = "AGAGTTTTTCGAACTCGTGTTGTCGAGCGACGGAATTAGATCAGTTAAATGGCAGAAAACTGGCAGGGCTTTTAGTCGTGGG"
taq = pcr_annealing_temperature(forward, reverse, product, verbose=False)
dbd = pcr_annealing_temperature(forward, reverse, product, polymerase="dbd", verbose=False)
print(taq["forward_tm_celsius"], taq["reverse_tm_celsius"], taq["annealing_celsius"])
print(dbd["forward_tm_celsius"], dbd["reverse_tm_celsius"], dbd["annealing_celsius"])
```

**示例输出：**

```text
59.445 68.485 52.403
55.982 64.16 58.982
```

## 7) `read_dnacauldron_records`

- **作用：** 用 DnaCauldron 读取本地序列文件中的全部记录。
- **处理方法：** 调用 dnacauldron.biotools.sequence_io 的 `load_records_from_files` 与 `load_records_from_file`。与 `read_dnacauldron_parts` 不同，一个文件可含多条记录，格式按内容识别，可读 `.zip` 压缩包，也可用 `folder` 读取整个目录（按路径排序）。没有 ID 的记录用文件名命名，一个文件含多条记录时加编号；`use_file_names_as_ids=True` 时，只含一条记录的文件用文件名作 ID。
- **API：** `read_dnacauldron_records(sources[可选], folder[可选], use_file_names_as_ids[可选], max_records[可选], verbose[可选])`
- **输入：** 文件路径列表 `sources` 或目录 `folder`（至少给出一个）、是否用文件名作 ID 与最多记录数（默认 500）。
- **输出：** 字典，含 `file_count`、`record_count`、`files`（`file`、`format`、`record_count`，压缩包的格式为 `zip`）、`records`（`id`、`file`、`length`、`topology`、`sequence`）与 `provenance`。

**示例文件：** `parts.fa`

```text
>partA
ATGCGTACGTTAGC
>partB
GGCATCGATCGGA
```

**示例文件：** `backbone.fasta`

```text
>
ACGTACGTAGCTAGCTAGC
```

**示例代码：**

```python
from dnakit.molecular_lab import read_dnacauldron_records

result = read_dnacauldron_records(["parts.fa", "backbone.fasta"], verbose=False)
print(result["file_count"], result["record_count"])
print([(item["file"], item["format"], item["record_count"]) for item in result["files"]])
print([(row["id"], row["length"]) for row in result["records"]])
```

**示例输出：**

```text
2 3
[('parts.fa', 'fasta', 2), ('backbone.fasta', 'fasta', 1)]
[('partA', 14), ('partB', 13), ('backbone', 19)]
```

## 8) `parse_dnacauldron_sequence_text`

- **作用：** 识别序列文本的格式并读取其中的记录。
- **处理方法：** 调用 dnacauldron.biotools.sequence_io.string_to_records。文本可以是纯 ACGT 字符串（格式 `ATGC`，一条无名记录）、FASTA 或 GenBank。与 `parse_sequence_text` 不同，纯序列也能读取，并且会报告识别到的格式；其他格式抛出 `ConfigurationError`。
- **API：** `parse_dnacauldron_sequence_text(text[必须], verbose[可选])`
- **输入：** 非空文本，不超过 1 MB。
- **输出：** 字典，含 `format`、`record_count`、`records`（`id`、`length`、`topology`、`sequence`）与 `provenance`。

**示例代码：**

```python
from dnakit.molecular_lab import parse_dnacauldron_sequence_text

bare = parse_dnacauldron_sequence_text("ACGTACGT", verbose=False)
print(bare["format"], bare["records"][0]["sequence"])
fasta = parse_dnacauldron_sequence_text(">x\nACGT\n>y\nGGCC\n", verbose=False)
print(fasta["format"], [row["id"] for row in fasta["records"]])
```

**示例输出：**

```text
ATGC ACGTACGT
fasta ['x', 'y']
```

## 9) `find_sequence_difference_positions`

- **作用：** 返回两条等长 DNA 序列逐位置的差异向量。
- **处理方法：** 调用 dnachisel.biotools.sequences_differences.sequences_differences_array。`difference_array` 每个位置一个值：与参考序列不同为 1，相同为 0；`positions` 列出不同位置（从 0 开始）。与 `find_sequence_difference_segments` 不同，它不合并成片段。不做比对。
- **API：** `find_sequence_difference_positions(sequence[必须], reference[必须], verbose[可选])`
- **输入：** 两条长度相同的 A/C/G/T 序列。
- **输出：** 字典，含 `sequence_length`、`difference_count`、`positions`、`difference_array` 与 `provenance`。

**示例代码：**

```python
from dnakit.molecular_lab import find_sequence_difference_positions

result = find_sequence_difference_positions("ACGTACGTAC", "ACCTACGAAC", verbose=False)
print(result["difference_count"], result["positions"])
print(result["difference_array"])
```

**示例输出：**

```text
2 [2, 7]
[0, 0, 1, 0, 0, 0, 0, 1, 0, 0]
```

## 10) `annotate_sequence_differences`

- **作用：** 把序列与参考序列的差异标注为带标签的特征。
- **处理方法：** 调用 dnachisel.biotools.annotate_differences。每一段相邻的差异位置成为一个特征（`start` 含、`end` 不含），标签为 `prefix` 加上参考序列在该处的碱基，例如 `#A` 或 `#TT`，这是 `find_sequence_difference_segments` 不提供的。两条序列必须等长；完全相同时得到 0 个特征（上游在此情况下会抛出 IndexError）。
- **API：** `annotate_sequence_differences(sequence[必须], reference[必须], feature_type[可选], prefix[可选], verbose[可选])`
- **输入：** 两条长度相同的 A/C/G/T 序列、特征类型（默认 `misc_feature`）与标签前缀（默认 `#`，最长 20 个字符）。
- **输出：** 字典，含 `sequence_length`、`feature_count`、`features`（`start`、`end`、`type`、`label`）与 `provenance`。

**示例代码：**

```python
from dnakit.molecular_lab import annotate_sequence_differences

result = annotate_sequence_differences("AATTGGCCAC", "AACCGGCCGC", verbose=False)
for row in result["features"]:
    print(row["start"], row["end"], row["type"], row["label"])
same = annotate_sequence_differences("ACGT", "ACGT", verbose=False)
print(same["feature_count"])
```

**示例输出：**

```text
2 4 misc_feature #CC
8 9 misc_feature #G
0
```
