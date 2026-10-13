# 85 G17a：alphagenome、grelu 与 enformer-pytorch 的离线 DNA 工具

本文档介绍 G17a 组封装的 32 个函数。它们来自 alphagenome、grelu 与 enformer-pytorch 的 DNA 级公开函数，全部使用本地输入离线运行；模型推理、权重下载、基因组下载和绘图没有封装。示例中的输出均由实际运行生成。封装函数位于 `dnakit.predictions` 的 `alphagenome_functions`、`grelu_functions` 与 `enformer_pytorch_functions` 模块中，示例直接从模块导入；API 行给出的是包级导出名，需在 `dnakit.predictions` 的 `__init__` 中接入后才能从包顶层导入。

## 覆盖统计

| 包 | 状态 | 数量 |
| --- | --- | ---: |
| alphagenome | internal | 18 |
| alphagenome | plot | 12 |
| alphagenome | wrapped | 9 |
| alphagenome | other | 8 |
| alphagenome | family-variant | 7 |
| alphagenome | not-dna | 6 |
| alphagenome | already-wrapped | 2 |
| alphagenome | needs-network | 2 |
| alphagenome | duplicate-of | 1 |
| grelu | internal | 38 |
| grelu | needs-network | 33 |
| grelu | other | 31 |
| grelu | wrapped | 20 |
| grelu | duplicate-of | 17 |
| grelu | plot | 14 |
| grelu | family-variant | 10 |
| grelu | already-wrapped | 8 |
| grelu | deprecated | 6 |
| grelu | cap | 1 |
| grelu | needs-binary | 1 |
| enformer_pytorch | internal | 23 |
| enformer_pytorch | other | 10 |
| enformer_pytorch | already-wrapped | 3 |
| enformer_pytorch | wrapped | 3 |
| enformer_pytorch | duplicate-of | 2 |
| enformer_pytorch | needs-network | 2 |
| enformer_pytorch | family-variant | 1 |

状态说明：`wrapped` 为本文档中的函数；`already-wrapped`、`duplicate-of` 指向 dnakit 已有函数；`family-variant` 由封装函数的参数选择；`needs-network` 需要下载；`other` 中以 needs-data 开头的原因需要训练模型或外部数据。

## 1) 左对齐 indel（alphagenome normalize_variant）

- **作用：** 把插入或缺失左对齐到最左的等价表示，便于比较不同来源的 indel。
- **处理方法：** 用本地参考窗口提供参考碱基，调用 alphagenome 的 normalize_variant；SNV 原样返回。
- **API：** `dnakit.predictions.normalize_alphagenome_variant(variant, *, reference_sequence, reference_start)`
- **输入：** variant 为 `chr:pos:ref>alt`（1 基）；reference_sequence 为从 reference_start（0 基）开始的参考序列，需覆盖变异前一个碱基及参考跨度。
- **输出：** 返回 input、normalized、changed 三个字段，位置为 1 基。

- **示例代码：**

```python
from dnakit.predictions import normalize_alphagenome_variant

result = normalize_alphagenome_variant(
    "chr1:4:AA>A", reference_sequence="ACAAAG", reference_start=0
)
print(result["normalized"])
print(result["changed"])
```

- **示例输出：**

```text
{'chromosome': 'chr1', 'position': 2, 'reference': 'CA', 'alternate': 'C'}
True
```

## 2) 区间交集与并集（alphagenome interval algebra）

- **作用：** 对两组不重叠的基因组区间求交集或并集。
- **处理方法：** 调用 intersect_intervals 或 union_intervals；输入每侧先检查同一染色体内是否重叠。
- **API：** `dnakit.predictions.combine_alphagenome_intervals(left, right, *, operation='intersect')`
- **输入：** left、right 为 `chr:start-end` 字符串列表（0 基半开）；operation 为 intersect 或 union。
- **输出：** 返回按染色体和坐标排序的 intervals 列表与计数。

- **示例代码：**

```python
from dnakit.predictions import combine_alphagenome_intervals

left = ["chr1:0-10", "chr1:20-30"]
right = ["chr1:5-25"]
for mode in ("intersect", "union"):
    result = combine_alphagenome_intervals(left, right, operation=mode)
    print(mode, [(item["start"], item["end"]) for item in result["intervals"]])
```

- **示例输出：**

```text
intersect [(5, 10), (20, 25)]
union [(0, 30)]
```

## 3) 从 GTF 提取转录起始位点（TSS）

- **作用：** 从 GENCODE GTF 得到每个转录本或基因的 5' 端位点。
- **处理方法：** 读取 transcript 或 gene 行，调用 alphagenome 的 extract_tss；+ 链取起点，- 链取终点（开区间）。
- **API：** `dnakit.predictions.extract_alphagenome_tss(gtf_path, *, feature='transcript', max_records=1000000)`
- **输入：** gtf_path 为本地 GTF 文件（可 gzip）；feature 为 transcript 或 gene；max_records 限制读取行数。
- **输出：** 返回 tss 列表，每项含染色体、位置（0 基点）、链、基因和转录本 ID。
- **示例文件：** `small.gtf`

```text
##description: small annotation for G17a tests
chr1	HAVANA	gene	1001	2000	.	+	.	gene_id "ENSG1.1"; gene_type "protein_coding"; gene_name "GENEA";
chr1	HAVANA	transcript	1001	1800	.	+	.	gene_id "ENSG1.1"; transcript_id "ENST1.1"; gene_type "protein_coding"; gene_name "GENEA"; transcript_type "protein_coding"; transcript_support_level "1"; tag "basic MANE_Select";
chr1	HAVANA	transcript	1101	2000	.	+	.	gene_id "ENSG1.1"; transcript_id "ENST2.1"; gene_type "protein_coding"; gene_name "GENEA"; transcript_type "protein_coding"; transcript_support_level "3"; tag "basic";
chr1	HAVANA	gene	5001	5500	.	-	.	gene_id "ENSG2.1"; gene_type "lncRNA"; gene_name "GENEB";
chr1	HAVANA	transcript	5001	5500	.	-	.	gene_id "ENSG2.1"; transcript_id "ENST3.1"; gene_type "lncRNA"; gene_name "GENEB"; transcript_type "lncRNA"; transcript_support_level "NA"; tag "basic";
```

- **示例代码：**

```python
from dnakit.predictions import extract_alphagenome_tss

result = extract_alphagenome_tss("small.gtf")
print([item["position"] for item in result["tss"]])
print([item["strand"] for item in result["tss"]])
```

- **示例输出：**

```text
[1000, 1100, 5500]
['+', '+', '-']
```

## 4) 按转录本规则筛选 GTF（protein_coding、最长、MANE、TSL、类型）

- **作用：** 按编码、最长转录本、MANE Select、转录本支持等级或转录本类型筛选 GENCODE 注释。
- **处理方法：** 调用 alphagenome 的对应筛选函数；method 选择规则，protein_coding 为默认。
- **API：** `dnakit.predictions.filter_alphagenome_transcripts(gtf_path, *, method='protein_coding', transcript_types=None, support_levels=None, include_gene_entries=False, max_records=1000000)`
- **输入：** gtf_path 为本地 GTF；method 取 protein_coding、transcript_type、longest、mane_select、support_level 之一；transcript_types 与 support_levels 按所选规则给出。
- **输出：** 返回保留的行数、各 Feature 计数、基因 ID 与转录本 ID（排序）。
- **示例文件：** `small.gtf`

```text
##description: small annotation for G17a tests
chr1	HAVANA	gene	1001	2000	.	+	.	gene_id "ENSG1.1"; gene_type "protein_coding"; gene_name "GENEA";
chr1	HAVANA	transcript	1001	1800	.	+	.	gene_id "ENSG1.1"; transcript_id "ENST1.1"; gene_type "protein_coding"; gene_name "GENEA"; transcript_type "protein_coding"; transcript_support_level "1"; tag "basic MANE_Select";
chr1	HAVANA	transcript	1101	2000	.	+	.	gene_id "ENSG1.1"; transcript_id "ENST2.1"; gene_type "protein_coding"; gene_name "GENEA"; transcript_type "protein_coding"; transcript_support_level "3"; tag "basic";
chr1	HAVANA	gene	5001	5500	.	-	.	gene_id "ENSG2.1"; gene_type "lncRNA"; gene_name "GENEB";
chr1	HAVANA	transcript	5001	5500	.	-	.	gene_id "ENSG2.1"; transcript_id "ENST3.1"; gene_type "lncRNA"; gene_name "GENEB"; transcript_type "lncRNA"; transcript_support_level "NA"; tag "basic";
```

- **示例代码：**

```python
from dnakit.predictions import filter_alphagenome_transcripts

for method in ("protein_coding", "longest", "mane_select"):
    result = filter_alphagenome_transcripts("small.gtf", method=method)
    print(method, result["transcript_ids"])
```

- **示例输出：**

```text
protein_coding ['ENST1.1', 'ENST2.1']
longest ['ENST2.1', 'ENST3.1']
mane_select ['ENST1.1']
```

## 5) Ensembl 转录本/基因 ID 版本升级

- **作用：** 把旧版本的 Ensembl ID 映射到新版本，例如 ENST00010.1 到 ENST00010.3。
- **处理方法：** 去掉版本号后按无版本 ID 匹配，调用 alphagenome 的 upgrade_annotation_ids；无匹配的映射为 None。
- **API：** `dnakit.predictions.upgrade_alphagenome_annotation_ids(old_ids, new_ids, *, patchless=False)`
- **输入：** old_ids 与 new_ids 为非空字符串列表；新旧 ID 需带版本号，旧 ID 无版本时设 patchless=True。
- **输出：** 返回 mappings（保持输入顺序）、匹配数与总数。

- **示例代码：**

```python
from dnakit.predictions import upgrade_alphagenome_annotation_ids

result = upgrade_alphagenome_annotation_ids(
    ["ENST00010.1", "ENST00011.1"], ["ENST00010.3", "ENST00012.1"]
)
print(result["mappings"])
print(result["matched_count"])
```

- **示例输出：**

```text
[{'old_id': 'ENST00010.1', 'new_id': 'ENST00010.3'}, {'old_id': 'ENST00011.1', 'new_id': None}]
1
```

## 6) 列出 ISM 单碱基替换变异

- **作用：** 枚举一个区间内所有单碱基替换，用于体内突变扫描（ISM）。
- **处理方法：** 调用 alphagenome 的 ism_variants；跳过参考碱基，按位置再按 vocabulary 顺序输出。
- **API：** `dnakit.predictions.enumerate_alphagenome_ism_variants(interval, sequence, *, vocabulary='ACGT', skip_n=False)`
- **输入：** interval 为 `chr:start-end`；sequence 为与区间等长的 DNA；vocabulary 默认 ACGT；skip_n 控制是否跳过 N。
- **输出：** 返回 variants 列表（1 基位置）、数量与区间信息。

- **示例代码：**

```python
from dnakit.predictions import enumerate_alphagenome_ism_variants

result = enumerate_alphagenome_ism_variants("chr1:100-103", "ACG")
print(result["count"])
print(result["variants"][0])
```

- **示例输出：**

```text
9
{'chromosome': 'chr1', 'position': 101, 'reference': 'A', 'alternate': 'C'}
```

## 7) 由 ISM 分数构建位置-碱基矩阵

- **作用：** 把每个 ISM 变异的分数整理成 (位置, 碱基) 矩阵，供可视化或汇总。
- **处理方法：** 按 enumerate 的顺序生成变异，调用 alphagenome 的 ism_matrix；每个位置减去该位置的均值。
- **API：** `dnakit.predictions.build_alphagenome_ism_matrix(interval, sequence, scores, *, vocabulary='ACGT', multiply_by_sequence=True)`
- **输入：** interval、sequence（仅 ACGT）、scores（与 ISM 变异一一对应，长度为 3 × 区间宽度）；multiply_by_sequence 为 True 时只保留参考碱基位置。
- **输出：** 返回矩阵（行为位置，列为 A、C、G、T）。

- **示例代码：**

```python
from dnakit.predictions import build_alphagenome_ism_matrix

result = build_alphagenome_ism_matrix("chr1:100-101", "A", [0.5, 0.25, 0.75])
print(result["matrix"])
full = build_alphagenome_ism_matrix(
    "chr1:100-101", "A", [0.5, 0.25, 0.75], multiply_by_sequence=False
)
print(full["matrix"])
```

- **示例输出：**

```text
[[-0.5, 0.0, 0.0, 0.0]]
[[-0.5, 0.0, -0.25, 0.25]]
```

## 8) 校验 AlphaGenome 支持的序列长度

- **作用：** 检查一个输入长度是否属于模型支持的长度（16 KB、100 KB、500 KB、1 MB）。
- **处理方法：** 调用 alphagenome 的 validate_sequence_length；不加载模型，也不发送请求。
- **API：** `dnakit.predictions.validate_alphagenome_sequence_length(sequence_length)`
- **输入：** sequence_length 为正整数。
- **输出：** 通过时返回 supported_lengths；不支持的长度抛出 ConfigurationError。

- **示例代码：**

```python
from dnakit.predictions import validate_alphagenome_sequence_length

result = validate_alphagenome_sequence_length(16384)
print(result["supported_lengths"])
```

- **示例输出：**

```text
[16384, 131072, 524288, 1048576]
```

## 9) 读取 AlphaGenome 数据折（fold）区间

- **作用：** 从本地区域文件取出某个模型版本、某个子集（训练、验证、测试）的区间。
- **处理方法：** 调用 alphagenome 的 get_fold_intervals 并传入本地文件，不下载默认区域文件。
- **API：** `dnakit.predictions.get_alphagenome_fold_intervals(regions_path, *, model_version='fold_0', subset='test', organism='homo_sapiens')`
- **输入：** regions_path 为无表头、制表符分隔的 BED 文件（染色体、起点、终点、fold）；model_version 为 all_folds 或 fold_0 至 fold_3；subset 为 train、valid、test；organism 默认 homo_sapiens。
- **输出：** 返回 intervals、fold_names 与数量。
- **示例文件：** `regions.bed`

```text
chr1	0	100	fold0
chr2	200	300	fold1
chr3	400	500	fold2
chr1	600	700	fold0
```

- **示例代码：**

```python
from dnakit.predictions import get_alphagenome_fold_intervals

result = get_alphagenome_fold_intervals("regions.bed", model_version="fold_0", subset="valid")
print([(item["chromosome"], item["start"]) for item in result["intervals"]])
print(result["fold_names"])
```

- **示例输出：**

```text
[('chr1', 0), ('chr1', 600)]
['fold0']
```

## 10) DNA 编码互转（indices、one-hot、字符串）

- **作用：** 在 gReLU 的整数索引、one-hot 矩阵与 DNA 字符串之间转换一批序列。
- **处理方法：** 按 source 和 target 调用 gReLU 的 format 函数；索引 A=0、C=1、G=2、T=3、N=4，N 的 one-hot 列全为 0。
- **API：** `dnakit.predictions.convert_dna_encoding(values, *, source, target)`
- **输入：** values 为批量列表；one-hot 为 4 × L 矩阵（行为 A、C、G、T）；同一批序列长度需相同。字符串输入请用 dnakit 的 one_hot_encode_dna。
- **输出：** 返回 values（字符串、索引或 one-hot 列表）与 count。

- **示例代码：**

```python
from dnakit.predictions import convert_dna_encoding

one_hot = [[[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0]]]
print(convert_dna_encoding(one_hot, source="one_hot", target="strings")["values"])
print(convert_dna_encoding([[0, 3]], source="indices", target="one_hot")["values"])
```

- **示例输出：**

```text
['ACGTN']
[[[1.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 1.0]]]
```

## 11) 检查 DNA 编码是否有效

- **作用：** 判断一批字符串、索引或 one-hot 矩阵是否是合法的 DNA 编码。
- **处理方法：** 调用 gReLU 的检查函数，并补充 one-hot 每列至多一个 1、索引范围 0 至 4 的规则；字符串先转大写。
- **API：** `dnakit.predictions.check_dna_encoding(values, *, encoding='strings')`
- **输入：** values 为批量输入；encoding 取 strings、indices、one_hot。结构错误（不等长、非整数）会抛出 ConfigurationError。
- **输出：** 返回 valid（布尔）、encoding 与 count。

- **示例代码：**

```python
from dnakit.predictions import check_dna_encoding

print(check_dna_encoding(["ACGTN", "acgt"], encoding="strings")["valid"])
print(check_dna_encoding(["ACGX"], encoding="strings")["valid"])
print(check_dna_encoding([[[1, 1], [1, 0], [0, 0], [0, 0]]], encoding="one_hot")["valid"])
```

- **示例输出：**

```text
True
False
False
```

## 12) GC 含量与 GC 直方图

- **作用：** 计算每条序列的 GC 比例，或统计一组序列的 GC 分布。
- **处理方法：** mode=fraction 调用 gReLU 的 gc；mode=distribution 调用 gc_distribution，按 binwidth 分箱（从 0 开始）。
- **API：** `dnakit.predictions.gc_content_distribution(sequences, *, mode='distribution', binwidth=0.1, normalize=False)`
- **输入：** sequences 为 DNA 字符串或列表；binwidth 在 (0, 1] 内；normalize 为 True 时计数和为 1。
- **输出：** fraction 模式返回 gc 列表；distribution 模式返回 bin_edges 与 counts。

- **示例代码：**

```python
from dnakit.predictions import gc_content_distribution

seqs = ["ACGT", "GGGG", "AAAA"]
print(gc_content_distribution(seqs, mode="fraction")["gc"])
hist = gc_content_distribution(seqs, mode="distribution", binwidth=0.25)
print(hist["bin_edges"], hist["counts"])
```

- **示例输出：**

```text
[0.5, 1.0, 0.0]
[0.0, 0.25, 0.5, 0.75] [1.0, 0.0, 1.0, 1.0]
```

## 13) 序列长度汇总

- **作用：** 给出每条序列的长度，并判断是否等长。
- **处理方法：** 调用 gReLU 的 get_lengths、check_equal_lengths 与 get_unique_length（仅等长时）。
- **API：** `dnakit.predictions.summarize_sequence_lengths(sequences)`
- **输入：** sequences 为 DNA 字符串或列表。
- **输出：** 返回 lengths、all_equal、unique_length（不等长时为 None）。

- **示例代码：**

```python
from dnakit.predictions import summarize_sequence_lengths

print(summarize_sequence_lengths(["ACG", "TTT"]))
print(summarize_sequence_lengths(["ACG", "TT"])["unique_length"])
```

- **示例输出：**

```text
{'status': 'completed', 'operation': 'summarize_sequence_lengths', 'count': 2, 'lengths': [3, 3], 'all_equal': True, 'unique_length': 3, 'provenance': {'package': 'grelu', 'version': '1.1.0', 'length_unit': 'bases'}}
None
```

## 14) 两条等长序列的差异位点

- **作用：** 列出两条等长 DNA 序列不同的位置及两侧碱基。
- **处理方法：** 调用 gReLU 的 seq_differences（不打印），再补充参考碱基与替代碱基。
- **API：** `dnakit.predictions.dna_sequence_differences(sequence_a, sequence_b)`
- **输入：** sequence_a 与 sequence_b 为等长 DNA 字符串。
- **输出：** 返回 differences 列表（0 基位置、参考碱基、替代碱基）与 count。

- **示例代码：**

```python
from dnakit.predictions import dna_sequence_differences

print(dna_sequence_differences("ACGT", "ACCT")["differences"])
```

- **示例输出：**

```text
[{'position': 2, 'reference': 'G', 'alternate': 'C'}]
```

## 15) 按summit扩展为固定长度区间

- **作用：** 以峰的 summit 为中心，生成固定长度的区间（如 narrowPeak）。
- **处理方法：** 调用 gReLU 的 extend_from_coord，起点为 start + summit - seq_len // 2；summit 是相对 start 的偏移。
- **API：** `dnakit.predictions.extend_intervals_from_centers(intervals, *, seq_len, center_col='summit')`
- **输入：** intervals 为含 chrom、start、end 与 summit 列的映射列表；seq_len 为正整数；扩展后起点小于 0 时报错。
- **输出：** 返回 intervals（chrom、start、end）与 count。

- **示例代码：**

```python
from dnakit.predictions import extend_intervals_from_centers

rows = [{"chrom": "chr1", "start": 100, "end": 200, "summit": 50}]
print(extend_intervals_from_centers(rows, seq_len=10)["intervals"])
```

- **示例输出：**

```text
[{'chrom': 'chr1', 'start': 145, 'end': 155}]
```

## 16) 区间裁剪到边界

- **作用：** 把区间起点提高到 start、终点降低到 end。
- **处理方法：** 调用 gReLU 的 clip_intervals；裁剪后起点不小于终点的行保留并计入 invalid_count。
- **API：** `dnakit.predictions.clip_intervals_to_bounds(intervals, *, start=None, end=None)`
- **输入：** intervals 为映射列表；start 与 end 可选。
- **输出：** 返回 intervals、count 与 invalid_count。

- **示例代码：**

```python
from dnakit.predictions import clip_intervals_to_bounds

rows = [{"chrom": "chr1", "start": 10, "end": 50}, {"chrom": "chr2", "start": 200, "end": 300}]
result = clip_intervals_to_bounds(rows, start=0, end=150)
print(result["intervals"])
print(result["invalid_count"])
```

- **示例输出：**

```text
[{'chrom': 'chr1', 'start': 10, 'end': 50}, {'chrom': 'chr2', 'start': 200, 'end': 150}]
1
```

## 17) 按分组列合并区间

- **作用：** 把同一分组值的区间合并为一个跨度（最小起点到最大终点）。
- **处理方法：** 调用 gReLU 的 merge_intervals_by_column；同一分组跨染色体时报错。
- **API：** `dnakit.predictions.merge_intervals_by_group(intervals, *, group_col)`
- **输入：** intervals 为映射列表，需包含 group_col 列。
- **输出：** 返回每个分组一条区间及 count。

- **示例代码：**

```python
from dnakit.predictions import merge_intervals_by_group

rows = [
    {"chrom": "chr1", "start": 0, "end": 10, "grp": "a"},
    {"chrom": "chr1", "start": 5, "end": 20, "grp": "a"},
    {"chrom": "chr2", "start": 0, "end": 5, "grp": "b"},
]
print(merge_intervals_by_group(rows, group_col="grp")["intervals"])
```

- **示例输出：**

```text
[{'chrom': 'chr1', 'start': 0, 'end': 20, 'grp': 'a'}, {'chrom': 'chr2', 'start': 0, 'end': 5, 'grp': 'b'}]
```

## 18) 保留或剔除与参考区间重叠的区间

- **作用：** 按与参考区间集合的重叠关系筛选区间。
- **处理方法：** 调用 gReLU 的 filter_overlapping；window 会使参考区间两侧各扩展相应碱基数。
- **API：** `dnakit.predictions.filter_intervals_overlapping_reference(intervals, reference, *, window=0, invert=False, method='any')`
- **输入：** intervals 与 reference 为映射列表；method 为 any 或 all；invert 为 True 时反向保留。
- **输出：** 返回保留的 intervals、计数与 messages（上游进度信息）。

- **示例代码：**

```python
from dnakit.predictions import filter_intervals_overlapping_reference

intervals = [{"chrom": "chr1", "start": 0, "end": 50}, {"chrom": "chr1", "start": 100, "end": 150}]
reference = [{"chrom": "chr1", "start": 40, "end": 60}]
print(filter_intervals_overlapping_reference(intervals, reference)["intervals"])
```

- **示例输出：**

```text
[{'chrom': 'chr1', 'start': 0, 'end': 50}]
```

## 19) 按染色体划分训练、验证、测试集

- **作用：** 按染色体名单划分区间，默认验证集为 chr10、测试集为 chr11。
- **处理方法：** 调用 gReLU 的 split；train_chroms 为 None 时训练集为其余全部染色体；三组名单不得重叠。
- **API：** `dnakit.predictions.split_intervals_by_chromosome(intervals, *, train_chroms=None, val_chroms=('chr10',), test_chroms=('chr11',))`
- **输入：** intervals 为映射列表；train_chroms、val_chroms、test_chroms 为染色体名列表。
- **输出：** 返回 train、validation、test 三组区间与计数。

- **示例代码：**

```python
from dnakit.predictions import split_intervals_by_chromosome

rows = [
    {"chrom": "chr1", "start": 0, "end": 50},
    {"chrom": "chr10", "start": 0, "end": 50},
    {"chrom": "chr11", "start": 0, "end": 50},
]
result = split_intervals_by_chromosome(rows)
print(result["counts"])
print([item["chrom"] for item in result["test"]])
```

- **示例输出：**

```text
{'train': 1, 'validation': 1, 'test': 1}
['chr11']
```

## 20) 剔除 ENCODE 黑名单区域

- **作用：** 去掉与 ENCODE 黑名单重叠的区间（hg19、hg38、mm10 随 gReLU 一起提供，无需下载）。
- **处理方法：** 调用 gReLU 的 filter_blacklist，使用内置黑名单文件；window 为两侧扩展的碱基数。
- **API：** `dnakit.predictions.filter_intervals_blacklist(intervals, *, genome='hg38', window=0)`
- **输入：** intervals 为映射列表；genome 取 hg19、hg38 或 mm10。
- **输出：** 返回保留的 intervals、input_count、count 与 removed_count。

- **示例代码：**

```python
from dnakit.predictions import filter_intervals_blacklist

rows = [
    {"chrom": "chr1", "start": 0, "end": 500},
    {"chrom": "chr1", "start": 1000, "end": 1500},
    {"chrom": "chr2", "start": 5000, "end": 5500},
]
result = filter_intervals_blacklist(rows, genome="hg38")
print(result["intervals"])
print(result["removed_count"])
```

- **示例输出：**

```text
[{'chrom': 'chr2', 'start': 5000, 'end': 5500}]
2
```

## 21) 按染色体包含或排除区间

- **作用：** 只保留指定染色体的区间，或排除指定染色体。
- **处理方法：** 调用 gReLU 的 filter_chromosomes；include 与 exclude 至少给一个。
- **API：** `dnakit.predictions.filter_intervals_by_chromosome(intervals, *, include=None, exclude=None)`
- **输入：** intervals 为映射列表；include、exclude 为染色体名列表。
- **输出：** 返回保留的 intervals 与计数。

- **示例代码：**

```python
from dnakit.predictions import filter_intervals_by_chromosome

rows = [{"chrom": "chr1", "start": 0, "end": 10}, {"chrom": "chr2", "start": 0, "end": 10}]
print(filter_intervals_by_chromosome(rows, include=["chr2"])["intervals"])
```

- **示例输出：**

```text
[{'chrom': 'chr2', 'start': 0, 'end': 10}]
```

## 22) PWM 的共识序列

- **作用：** 把 4 × L 的位置权重矩阵转为共识 DNA 字符串。
- **处理方法：** 调用 gReLU 的 motifs_to_strings，每个位置取权重最高的碱基；并列时取 A、C、G、T 中靠前的碱基。
- **API：** `dnakit.predictions.pwm_to_consensus(pwm)`
- **输入：** pwm 为 4 行（A、C、G、T）的数值列表。
- **输出：** 返回 consensus 与 length。

- **示例代码：**

```python
from dnakit.predictions import pwm_to_consensus

pwm = [
    [0.9, 0.05, 0.02],
    [0.05, 0.9, 0.03],
    [0.03, 0.02, 0.9],
    [0.02, 0.03, 0.05],
]
print(pwm_to_consensus(pwm)["consensus"])
```

- **示例输出：**

```text
ACG
```

## 23) 修剪 PWM 两端低信息位置

- **作用：** 去掉 PWM 两端权重较小的列。
- **处理方法：** 调用 gReLU 的 trim_pwm：每列得分为绝对值之和，保留得分不低于最大得分乘以 trim_threshold 的列。
- **API：** `dnakit.predictions.trim_motif_pwm(pwm, *, trim_threshold=0.3)`
- **输入：** pwm 为 4 × L 数值矩阵，应使用对数似然比或中心化后的矩阵；概率矩阵每列得分均为 1，不会被修剪。
- **输出：** 返回 start、end（0 基半开）、修剪后的 pwm 与 input_length。

- **示例代码：**

```python
from dnakit.predictions import trim_motif_pwm

pwm = [
    [0.0, 0.65, -0.2, -0.2, 0.0],
    [0.0, -0.2, 0.65, -0.2, 0.0],
    [0.0, -0.2, -0.2, 0.65, 0.0],
    [0.0, -0.25, -0.25, -0.25, 0.0],
]
result = trim_motif_pwm(pwm)
print(result["start"], result["end"])
```

- **示例输出：**

```text
1 4
```

## 24) 多 PWM 扫描 DNA 序列（FIMO 风格）

- **作用：** 在两条链上扫描 DNA，找出满足 p 值阈值的 motif 位点。
- **处理方法：** 调用 gReLU 的 scan_sequences，按 pthresh 过滤；结果按序列、起点、链与 motif 排序。
- **API：** `dnakit.predictions.scan_motifs_with_pwms(sequences, motifs, *, pthresh=0.001, rc=True)`
- **输入：** sequences 为 ID 到序列的映射或序列列表（ID 为 seq0、seq1 ...）；motifs 为名称到 4 × L 矩阵的映射；rc 控制是否扫描反向互补链。
- **输出：** 返回 hits（motif、sequence、start、end、strand、score、p_value、matched_sequence）与计数。

- **示例代码：**

```python
from dnakit.predictions import scan_motifs_with_pwms

motif = [[0.9, 0.05, 0.03], [0.05, 0.9, 0.03], [0.03, 0.03, 0.9], [0.02, 0.02, 0.04]]
result = scan_motifs_with_pwms({"s1": "TTTTACGTTT"}, {"M1": motif}, pthresh=0.01)
for hit in result["hits"]:
    print(hit["start"], hit["strand"], hit["matched_sequence"])
```

- **示例输出：**

```text
4 + ACG
5 - CGT
```

## 25) 比较参考与替代等位基因上的 PWM 位点

- **作用：** 找出 motif 在参考序列与替代序列之间的得分变化。
- **处理方法：** 调用 gReLU 的 compare_motifs，每个位点给出参考与替代的 p 值和得分，并计算替代减参考。
- **API：** `dnakit.predictions.compare_motif_pwm_alleles(reference, alternate, motifs, *, pthresh=0.001, rc=True)`
- **输入：** reference 与 alternate 为 DNA 字符串；motifs 为名称到 4 × L 矩阵的映射。
- **输出：** 返回 sites 列表与计数；缺失值为 None。

- **示例代码：**

```python
from dnakit.predictions import compare_motif_pwm_alleles

motif = [[0.9, 0.05, 0.03], [0.05, 0.9, 0.03], [0.03, 0.03, 0.9], [0.02, 0.02, 0.04]]
result = compare_motif_pwm_alleles("TTTTACGTTT", "TTTTATGTTT", {"M1": motif}, pthresh=0.01)
for site in result["sites"]:
    print(site["start"], site["strand"], round(site["score_difference"], 4))
```

- **示例输出：**

```text
4 + -5.4848
5 - -14.8755
```

## 26) 按等位基因类型与 indel 长度过滤变异

- **作用：** 保留标准碱基的变异，并限制插入与缺失的长度。
- **处理方法：** 调用 gReLU 的 filter_variants；每行需含 chrom、pos（1 基）、ref、alt。
- **API：** `dnakit.predictions.filter_variants_by_length(variants, *, standard_bases=True, max_insert_len=0, max_del_len=0)`
- **输入：** variants 为映射列表；standard_bases 为 True 时只保留 A、C、G、T、N 等位基因；max_insert_len、max_del_len 为非负整数。
- **输出：** 返回保留的 variants、input_count、count 与 removed_count。

- **示例代码：**

```python
from dnakit.predictions import filter_variants_by_length

variants = [
    {"chrom": "chr1", "pos": 10, "ref": "A", "alt": "T"},
    {"chrom": "chr1", "pos": 20, "ref": "AC", "alt": "A"},
    {"chrom": "chr2", "pos": 30, "ref": "G", "alt": "GTT"},
]
result = filter_variants_by_length(variants)
print(result["variants"])
print(result["removed_count"])
```

- **示例输出：**

```text
[{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'T'}]
2
```

## 27) 变异周围的固定长度区间

- **作用：** 为每个变异生成以其位置为中心的固定长度区间，用于取序列窗口。
- **处理方法：** 调用 gReLU 的 variants_to_intervals：start = pos - ceil(seq_len / 2)，end = start + seq_len。
- **API：** `dnakit.predictions.variants_to_sequence_intervals(variants, *, seq_len=1)`
- **输入：** variants 为映射列表（chrom、pos、ref、alt）；seq_len 为正整数，默认 1。
- **输出：** 返回 intervals（chrom、start、end）与 count。

- **示例代码：**

```python
from dnakit.predictions import variants_to_sequence_intervals

variants = [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "T"}]
print(variants_to_sequence_intervals(variants, seq_len=10)["intervals"])
```

- **示例输出：**

```text
[{'chrom': 'chr1', 'start': 5, 'end': 15}]
```

## 28) 读取 FASTA 序列（gReLU 读取器）

- **作用：** 读取本地 FASTA 或 gzip FASTA 中的序列，按文件顺序返回。
- **处理方法：** 调用 gReLU 的 read_fasta；gReLU 不返回表头，因此结果只含序列与长度。
- **API：** `dnakit.predictions.read_fasta_sequences(fasta_path, *, max_sequences=1000)`
- **输入：** fasta_path 为本地文件；max_sequences 限制序列条数，超出时报错。
- **输出：** 返回 sequences、lengths、count 与 source。
- **示例文件：** `toy.fa`

```text
>a
ACGTTA
>b
GGCC
```

- **示例代码：**

```python
from dnakit.predictions import read_fasta_sequences

result = read_fasta_sequences("toy.fa")
print(result["sequences"])
print(result["lengths"])
```

- **示例输出：**

```text
['ACGTTA', 'GGCC']
[6, 4]
```

## 29) 按种子随机抽取区间

- **作用：** 不放回地随机抽取 n 个区间，结果可复现。
- **处理方法：** 调用 gReLU 的 filter_random 并在副本上原地执行（上游无 inplace 时返回 None）；输出按坐标排序。
- **API：** `dnakit.predictions.sample_intervals(intervals, *, n, seed=None)`
- **输入：** intervals 为映射列表；n 为不大于区间数的正整数；seed 为可选非负整数。
- **输出：** 返回 intervals（已排序）、input_count 与 seed。

- **示例代码：**

```python
from dnakit.predictions import sample_intervals

rows = [{"chrom": "chr1", "start": i * 10, "end": i * 10 + 5} for i in range(6)]
result = sample_intervals(rows, n=3, seed=1)
print([item["start"] for item in result["intervals"]])
```

- **示例输出：**

```text
[10, 20, 40]
```

## 30) DNA 转整数索引（enformer-pytorch）

- **作用：** 把 DNA 字符串转为 enformer-pytorch 使用的整数索引。
- **处理方法：** 调用 enformer_pytorch.data.str_to_seq_indices；A=0、C=1、G=2、T=3、N=4。
- **API：** `dnakit.predictions.encode_enformer_indices(sequence, *, max_length=10000)`
- **输入：** sequence 为 ACGTN 字符串；max_length 限制长度。
- **输出：** 返回 indices 列表、length 与索引约定。

- **示例代码：**

```python
from dnakit.predictions import encode_enformer_indices

print(encode_enformer_indices("ACGTN")["indices"])
```

- **示例输出：**

```text
[0, 1, 2, 3, 4]
```

## 31) 整数索引转 one-hot（支持填充）

- **作用：** 把整数索引转为 (L, 4) 的 one-hot 矩阵，并处理填充位置。
- **处理方法：** 调用 enformer_pytorch.data.seq_indices_to_one_hot；N（4）全为 0，填充位置（默认 -1）各通道为 0.25。
- **API：** `dnakit.predictions.indices_to_one_hot_enformer(indices, *, padding=-1)`
- **输入：** indices 为 0 至 4 或填充值的整数列表；padding 为负整数，默认 -1。
- **输出：** 返回 one_hot（每行一个位置）、length 与填充约定。

- **示例代码：**

```python
from dnakit.predictions import indices_to_one_hot_enformer

result = indices_to_one_hot_enformer([0, 1, 4, -1])
for row in result["one_hot"]:
    print(row)
```

- **示例输出：**

```text
[1.0, 0.0, 0.0, 0.0]
[0.0, 1.0, 0.0, 0.0]
[0.0, 0.0, 0.0, 0.0]
[0.25, 0.25, 0.25, 0.25]
```

## 32) enformer 编码的反向互补

- **作用：** 对整数索引或 one-hot 形式的 DNA 做反向互补。
- **处理方法：** 调用 seq_indices_reverse_complement（indices）或 one_hot_reverse_complement（one_hot）。
- **API：** `dnakit.predictions.reverse_complement_enformer(sequence, *, representation='indices')`
- **输入：** sequence 为 ACGTN 字符串；representation 取 indices 或 one_hot。
- **输出：** 返回 reverse_complement（索引列表或 one-hot 行）与 length。

- **示例代码：**

```python
from dnakit.predictions import reverse_complement_enformer

print(reverse_complement_enformer("ACGTN")["reverse_complement"])
```

- **示例输出：**

```text
[4, 0, 1, 2, 3]
```
