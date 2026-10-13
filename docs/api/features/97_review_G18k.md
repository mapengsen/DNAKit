# 审查 G18k：varcode 与 scikit-allel 变异工具

本页记录 varcode 与 scikit-allel 中状态为 `other` 或 `cap` 的函数的处理结果。新增 10 个 varcode 包装函数与 31 个 scikit-allel 包装函数，全部离线运行并只接受本地输入。模块为 `dnakit.variants` 与 `dnakit.variants`，尚未接入 `dnakit.variants` 的导出列表，由主线统一接入。依赖：`bio-varcode`（varcode）与 `bio-scikit-allel`（scikit-allel，HDF5 归档需要 h5py）。

## 1) 判断单个核苷酸是否为嘌呤

- **作用：** 判断一个单字母核苷酸码是否为嘌呤（A 或 G）。
- **计算方法：** 调用 varcode.nucleotides.is_purine；非标准码在未设置 allow_extended 时报错。
- **API：** `is_purine_base(nucleotide, *, allow_extended=False, verbose=True)`
- **输入：** 单个字母字符串 nucleotide；可选 allow_extended。
- **输出：** 返回字典，含大写的 nucleotide、is_purine 与 provenance。
- **示例代码：**

```python
from dnakit.variants import is_purine_base

result = is_purine_base("A", verbose=False)
print(result["nucleotide"], result["is_purine"])
```

- **示例输出：**

```text
A True
```

## 2) 检查序列是否只含 A、C、G、T

- **作用：** 检查碱基串是否只含标准碱基，并给出非标准碱基的 1 起始位置。
- **计算方法：** 调用 varcode.nucleotides.all_standard_nucleotides；逐位列出不在 A、C、G、T 中的位置。
- **API：** `has_only_standard_nucleotides(nucleotides, *, verbose=True)`
- **输入：** 非空字母字符串 nucleotides，不含空白，最多 1,000,000 个字符。
- **输出：** 返回 length、standard（布尔）和 non_standard_positions（1 起始位置列表）。
- **示例代码：**

```python
from dnakit.variants import has_only_standard_nucleotides

result = has_only_standard_nucleotides("ACGTNACG", verbose=False)
print(result["standard"], result["non_standard_positions"])
```

- **示例输出：**

```text
False [5]
```

## 3) 判断染色体名是否指线粒体

- **作用：** 判断 chrM、MT、mtDNA 等染色体名是否指线粒体基因组，用于选择遗传密码表。
- **计算方法：** 调用 varcode.effects.codon_tables.is_mitochondrial_contig；大小写不敏感，None 返回 False。
- **API：** `is_mitochondrial_contig(contig, *, verbose=True)`
- **输入：** contig 字符串或 None。
- **输出：** 返回 contig 与 mitochondrial 布尔值。
- **示例代码：**

```python
from dnakit.variants import is_mitochondrial_contig

print(is_mitochondrial_contig("chrM", verbose=False)["mitochondrial"], is_mitochondrial_contig("chr1", verbose=False)["mitochondrial"])
```

- **示例输出：**

```text
True False
```

## 4) 判断基因组名是否为 UCSC 组装名

- **作用：** 判断基因组名是否为 varcode 收录的 UCSC 组装名，例如 hg19、mm10。
- **计算方法：** 调用 varcode.reference.is_ucsc_reference_name；名称先做大小写和分隔符规范化。
- **API：** `is_ucsc_reference_name(name, *, verbose=True)`
- **输入：** 非空字符串 name。
- **输出：** 返回 name 与 ucsc 布尔值。
- **示例代码：**

```python
from dnakit.variants import is_ucsc_reference_name

print(is_ucsc_reference_name("hg38", verbose=False)["ucsc"], is_ucsc_reference_name("GRCh38", verbose=False)["ucsc"])
```

- **示例输出：**

```text
True False
```

## 5) 从参考名或路径推断组装名

- **作用：** 从参考名或 FASTA 路径推断规范组装名，不读取文件内容。
- **计算方法：** 调用 varcode.reference.infer_reference_name，在组装名表中匹配文件名；无匹配时报 ConfigurationError。
- **API：** `infer_reference_assembly(reference, *, verbose=True)`
- **输入：** 参考名或路径字符串 reference。
- **输出：** 返回 reference 与 assembly。
- **示例代码：**

```python
from dnakit.variants import infer_reference_assembly

print(infer_reference_assembly("/data/refs/GRCh37.fa", verbose=False)["assembly"])
```

- **示例输出：**

```text
GRCh37
```

## 6) 计算断点连接对应的 DEL、DUP 或 INV 区间

- **作用：** 给出同一染色体上一条新连接所描述的 DEL、DUP 或 INV 区间。
- **计算方法：** 调用 varcode.structural_variant.junction_span_details；两端需在同一染色体且保留侧已知，否则 span 为 None。
- **API：** `breakend_junction_span(near, far, *, verbose=True)`
- **输入：** near、far 两个映射，含 contig、position（1 起始）和 keeps（"left"、"right" 或 None）。
- **输出：** span 为 None 或含 sv_type、start、end、affected_start、affected_end 的字典。
- **示例代码：**

```python
from dnakit.variants import breakend_junction_span

span = breakend_junction_span(
    {"contig": "chr1", "position": 100, "keeps": "left"},
    {"contig": "chr1", "position": 200, "keeps": "right"},
    verbose=False,
)["span"]
print(span["sv_type"], span["start"], span["end"], span["affected_start"], span["affected_end"])
```

- **示例输出：**

```text
DEL 100 199 101 199
```

## 7) 按断点拼接两段局部布局

- **作用：** 把两段参考区间各保留断点一侧，拼成带有新连接的有向布局。
- **计算方法：** 调用 varcode.genomic_layout.join_breakends；必要时反向拼接段的顺序，不读取序列。
- **API：** `join_breakend_layouts(first_segments, first_breakend, second_segments, second_breakend, *, verbose=True)`
- **输入：** 两组区段（contig、start、end，可选 strand），以及两个断点（contig、position、keeps 为 "left" 或 "right"）。
- **输出：** 返回 segments（contig、start、end、strand、length）与总长度 length。
- **示例代码：**

```python
from dnakit.variants import join_breakend_layouts

result = join_breakend_layouts(
    [{"contig": "chr1", "start": 100, "end": 110}],
    {"contig": "chr1", "position": 105, "keeps": "left"},
    [{"contig": "chr1", "start": 300, "end": 310}],
    {"contig": "chr1", "position": 302, "keeps": "right"},
    verbose=False,
)
print([(item["start"], item["end"]) for item in result["segments"]], result["length"])
```

- **示例输出：**

```text
[(100, 105), (302, 310)] 15
```

## 8) 把参考区间投影到局部布局上

- **作用：** 把参考区间投影到局部布局上，得到每个连续匹配片段。
- **计算方法：** 调用 varcode.genomic_layout.project_intervals；只按坐标求交，不展开碱基；倒位段不参与匹配。
- **API：** `project_intervals_onto_layout(segments, intervals, *, verbose=True)`
- **输入：** 布局区段 segments 与参考区间 intervals（contig、start、end，1 起始闭区间）。
- **输出：** 返回 count 与 pieces；每个片段含 segments 与 length。
- **示例代码：**

```python
from dnakit.variants import project_intervals_onto_layout

result = project_intervals_onto_layout(
    [{"contig": "chr1", "start": 100, "end": 110}],
    [{"contig": "chr1", "start": 103, "end": 120}],
    verbose=False,
)
piece = result["pieces"][0]["segments"][0]
print(result["count"], piece["start"], piece["end"], result["pieces"][0]["length"])
```

- **示例输出：**

```text
1 103 110 8
```

## 9) 合并本地 VCF 中 MATEID 配对的 BND 记录

- **作用：** 读取本地 VCF，把 MATEID 互指的 BND 记录合并成每个事件一条组合记录。
- **计算方法：** 用 varcode.load_vcf（解析结构变异）读取文件，再调用 varcode.transforms.pair_breakends；genome 只用于推断染色体命名。
- **API：** `pair_breakend_mates(path, *, genome="GRCh38", verbose=True)`
- **输入：** 本地 VCF 路径 path，可选 genome（默认 GRCh38）。
- **输出：** 返回 input_count、output_count、paired_count 与按位置排序的 records（contig、start、ref、alt、source_starts）。
- **示例代码：**

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit.variants import pair_breakend_mates

ROWS = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=1000>",
    '##ALT=<ID=DEL,Description="Deletion">',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
    '##INFO=<ID=MATEID,Number=1,Type=String,Description="Mate">',
    '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Length">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "chr1\t100\tbnd_a\tA\tA]chr1:300]\t50\tPASS\tSVTYPE=BND;MATEID=bnd_b\tGT\t0/1",
    "chr1\t300\tbnd_b\tT\t[chr1:100[T\t50\tPASS\tSVTYPE=BND;MATEID=bnd_a\tGT\t0/1",
    "chr1\t500\tdel1\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=600;SVLEN=-100\tGT\t0/1",
]
with TemporaryDirectory() as folder:
    path = Path(folder) / "bnd.vcf"
    path.write_text("\n".join(ROWS) + "\n", encoding="utf-8")
    result = pair_breakend_mates(str(path), verbose=False)
print(result["input_count"], result["output_count"], result["paired_count"])
for row in result["records"]:
    print(row["contig"], row["start"], row["alt"], row["source_starts"])
```

- **示例输出：**

```text
3 2 1
chr1 100 A]chr1:300] [100, 300]
chr1 500 <DEL> []
```

## 10) 用 varcode 规则重写本地 VCF

- **作用：** 读取本地 VCF，按 varcode 的 INFO 与 FORMAT 规则写出一个新 VCF。
- **计算方法：** 调用 varcode.vcf_output.variants_to_vcf；输出必须是新文件，写入失败时删除未完成的文件。
- **API：** `rewrite_vcf_records(path, output, *, genome="GRCh38", verbose=True)`
- **输入：** 源 VCF 路径 path、新输出路径 output（父目录需存在）、可选 genome。
- **输出：** 返回 output 路径与 variants 数量。
- **示例代码：**

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit.variants import rewrite_vcf_records

ROWS = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=1000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "chr1\t10\trs1\tA\tG\t50\tPASS\tDP=10\tGT\t0|1\t0/0\t1|1",
    "chr1\t20\trs2\tC\tT\t50\tPASS\tDP=12\tGT\t0/1\t1/1\t0/0",
    "chr1\t35\trs3\tG\tA\t50\tPASS\tDP=8\tGT\t1/1\t0/1\t0|0",
]
with TemporaryDirectory() as folder:
    source = Path(folder) / "small.vcf"
    source.write_text("\n".join(ROWS) + "\n", encoding="utf-8")
    output = Path(folder) / "rewritten.vcf"
    result = rewrite_vcf_records(str(source), str(output), verbose=False)
    last = output.read_text(encoding="utf-8").splitlines()[-1]
print(result["variants"], " ".join(last.split("\t")))
```

- **示例输出：**

```text
3 chr1 35 rs3 G A 50.0 PASS DP=8 GT 1/1 0/1 0|0
```

## 11) 找出只在一个群体中出现的等位基因

- **作用：** 对两到五十个群体的等位基因计数，找出只存在于一个群体中的等位基因。
- **计算方法：** 调用 allel.locate_private_alleles；某等位基因在恰好一个群体中计数大于 0 时记为私有。
- **API：** `locate_population_private_alleles(allele_counts, *, verbose=True)`
- **输入：** 按群体排列的计数列表；每个群体为 (变异数, 等位基因数) 矩阵，变异数须一致。
- **输出：** 返回 private_alleles（variant_index、allele_index）与 private_variant_indices。
- **示例代码：**

```python
from dnakit.variants import locate_population_private_alleles

counts = [
    [[4, 0], [2, 2], [0, 4]],
    [[0, 4], [2, 2], [4, 0]],
    [[4, 0], [2, 2], [1, 3]],
]
result = locate_population_private_alleles(counts, verbose=False)
print(result["private_variant_indices"], result["private_alleles"])
```

- **示例输出：**

```text
[0] [{'variant_index': 0, 'allele_index': 1}]
```

## 12) 找出两群体间的固定差异位点

- **作用：** 找出两个群体之间没有共享等位基因的变异（固定差异）。
- **计算方法：** 调用 allel.locate_fixed_differences。
- **API：** `locate_population_fixed_differences(first_counts, second_counts, *, verbose=True)`
- **输入：** 两个群体的等位基因计数矩阵，行数相同。
- **输出：** 返回 fixed_indices 与 fixed_count。
- **示例代码：**

```python
from dnakit.variants import locate_population_fixed_differences

result = locate_population_fixed_differences([[4, 0], [2, 2], [0, 4]], [[0, 4], [2, 2], [3, 1]], verbose=False)
print(result["fixed_indices"], result["fixed_count"])
```

- **示例输出：**

```text
[0] 1
```

## 13) 统计滑动窗口内的变异数

- **作用：** 统计一条染色体上每个固定长度窗口内的变异数。
- **计算方法：** 调用 allel.windowed_count；窗口为 1 起始闭区间，默认步长等于窗口长度。
- **API：** `count_variants_in_windows(positions, *, window_size, start=None, stop=None, step=None, verbose=True)`
- **输入：** 升序的 1 起始位置列表 positions、窗口长度 window_size，以及可选 start、stop、step。
- **输出：** 返回 windows，每项含 start、stop 与 variants。
- **示例代码：**

```python
from dnakit.variants import count_variants_in_windows

result = count_variants_in_windows([10, 20, 30, 40, 50, 60], window_size=20, verbose=False)
print([(w["start"], w["stop"], w["variants"]) for w in result["windows"]])
```

- **示例输出：**

```text
[(10, 29, 2), (30, 49, 2), (50, 60, 2)]
```

## 14) 将窗口统计量除以可调用碱基数

- **作用：** 把每个窗口的统计量除以窗口内的碱基数（或可访问碱基数），得到每碱基数值。
- **计算方法：** 调用 allel.per_base；提供 is_accessible 时只计数掩码为 True 的碱基。
- **API：** `windowed_values_per_base(values, windows, *, is_accessible=None, verbose=True)`
- **输入：** 每个窗口一个数值 values；窗口 (start, stop) 为 1 起始闭区间；可选从第 1 位起的布尔掩码。
- **输出：** 返回 windows，每项含 start、stop、callable_bases 与 value_per_base。
- **示例代码：**

```python
from dnakit.variants import windowed_values_per_base

result = windowed_values_per_base([1.0, 2.0], [[1, 3], [4, 6]], verbose=False)
print([round(w["value_per_base"], 4) for w in result["windows"]])
```

- **示例输出：**

```text
[0.3333, 0.6667]
```

## 15) 构造可访问碱基数相同的窗口

- **作用：** 从可访问性掩码构造窗口，每个窗口包含相同数量的可访问碱基。
- **计算方法：** 调用 allel.equally_accessible_windows；window_size 以可访问碱基计。
- **API：** `build_equally_accessible_windows(is_accessible, *, window_size, step=None, verbose=True)`
- **输入：** 布尔掩码 is_accessible（从第 1 位起）、window_size 与可选 step。
- **输出：** 返回 windows，每项含 start、stop 与 accessible_bases。
- **示例代码：**

```python
from dnakit.variants import build_equally_accessible_windows

mask = [True, True, False, True, True, True]
result = build_equally_accessible_windows(mask, window_size=2, verbose=False)
print([(w["start"], w["stop"]) for w in result["windows"]])
```

- **示例输出：**

```text
[(1, 2), (4, 5)]
```

## 16) 计算两群体间的平均成对差异

- **作用：** 计算每个变异在两群体之间的平均成对差异（pi_between）。
- **计算方法：** 调用 allel.mean_pairwise_difference_between，按等位基因计数计算。
- **API：** `calculate_between_population_pi(first_counts, second_counts, *, verbose=True)`
- **输入：** 两个群体的等位基因计数矩阵。
- **输出：** 返回 values，每个变异一个值；未定义时为 None。
- **示例代码：**

```python
from dnakit.variants import calculate_between_population_pi

first = [[4, 0], [2, 2], [3, 1]]
second = [[0, 4], [2, 2], [1, 3]]
result = calculate_between_population_pi(first, second, verbose=False)
print(result["values"])
```

- **示例输出：**

```text
[1.0, 0.5, 0.625]
```

## 17) 计算窗口内固定差异密度

- **作用：** 计算两群体之间固定差异在每个窗口内的密度。
- **计算方法：** 调用 allel.windowed_df；密度为固定差异数除以窗口内可调用碱基数。
- **API：** `calculate_windowed_fixed_difference_density(positions, first_counts, second_counts, *, window_size, start=None, stop=None, step=None, is_accessible=None, verbose=True)`
- **输入：** 位置、两群体计数矩阵、窗口参数，可选可访问掩码。
- **输出：** 返回 windows，每项含 start、stop、sites、callable_bases 与 density。
- **示例代码：**

```python
from dnakit.variants import calculate_windowed_fixed_difference_density

positions = [10, 20, 30, 40]
first = [[4, 0], [2, 2], [0, 4], [3, 1]]
second = [[0, 4], [2, 2], [4, 0], [1, 3]]
result = calculate_windowed_fixed_difference_density(positions, first, second, window_size=20, verbose=False)
print([(w["start"], w["stop"], w["sites"], w["callable_bases"], w["density"]) for w in result["windows"]])
```

- **示例输出：**

```text
[(10, 29, 2, 20, 0.05), (30, 40, 2, 11, 0.0909090909)]
```

## 18) 计算 Patterson Fst

- **作用：** 计算两群体的 Patterson Fst，给出每个变异的分子、分母以及总比值。
- **计算方法：** 调用 allel.patterson_fst；总体 Fst 为分子之和除以分母之和。
- **API：** `calculate_patterson_fst(first_counts, second_counts, *, verbose=True)`
- **输入：** 两个群体的等位基因计数矩阵。
- **输出：** 返回 numerator、denominator 与 fst。
- **示例代码：**

```python
from dnakit.variants import calculate_patterson_fst

first = [[4, 0], [2, 2], [0, 4], [3, 1]]
second = [[0, 4], [2, 2], [4, 0], [1, 3]]
result = calculate_patterson_fst(first, second, verbose=False)
print(result["numerator"], result["denominator"], result["fst"])
```

- **示例输出：**

```text
[1.0, -0.1666666667, 1.0, 0.125] [1.0, 0.5, 1.0, 0.625] 0.6266666667
```

## 19) 计算窗口 Patterson Fst

- **作用：** 计算每个碱基窗口内的 Patterson Fst。
- **计算方法：** 调用 allel.windowed_patterson_fst。
- **API：** `calculate_windowed_patterson_fst(positions, first_counts, second_counts, *, window_size, start=None, stop=None, step=None, verbose=True)`
- **输入：** 位置、两群体计数矩阵与窗口参数。
- **输出：** 返回 windows，每项含 start、stop、sites 与 fst。
- **示例代码：**

```python
from dnakit.variants import calculate_windowed_patterson_fst

positions = [10, 20, 30, 40]
first = [[4, 0], [2, 2], [0, 4], [3, 1]]
second = [[0, 4], [2, 2], [4, 0], [1, 3]]
result = calculate_windowed_patterson_fst(positions, first, second, window_size=20, verbose=False)
print([(w["start"], w["stop"], w["sites"], w["fst"]) for w in result["windows"]])
```

- **示例输出：**

```text
[(10, 29, 2, 0.5555555556), (30, 40, 2, 0.6923076923)]
```

## 20) 计算分块 jackknife 的 Patterson Fst

- **作用：** 计算 Patterson Fst 及其块 jackknife 标准误。
- **计算方法：** 调用 allel.blockwise_patterson_fst；block_length 按变异数计。
- **API：** `calculate_blockwise_patterson_fst(first_counts, second_counts, *, block_length, verbose=True)`
- **输入：** 两群体计数矩阵与 block_length。
- **输出：** 返回 fst、standard_error、block_values 与 jackknife_values。
- **示例代码：**

```python
from dnakit.variants import calculate_blockwise_patterson_fst

first = [[4, 0], [2, 2], [0, 4], [3, 1]]
second = [[0, 4], [2, 2], [4, 0], [1, 3]]
result = calculate_blockwise_patterson_fst(first, second, block_length=2, verbose=False)
print(result["fst"], result["standard_error"], result["block_values"])
```

- **示例输出：**

```text
0.6266666667 0.0683760684 [0.5555555556, 0.6923076923]
```

## 21) 计算滑动变异窗口的 Patterson Fst

- **作用：** 在连续变异组成的滑动窗口中计算 Patterson Fst。
- **计算方法：** 调用 allel.moving_patterson_fst；窗口按变异序号划分，而非碱基位置。
- **API：** `calculate_moving_patterson_fst(first_counts, second_counts, *, window_size, start=0, stop=None, step=None, verbose=True)`
- **输入：** 两群体计数矩阵与窗口参数。
- **输出：** 返回 window_size、step 与 values。
- **示例代码：**

```python
from dnakit.variants import calculate_moving_patterson_fst

first = [[4, 0], [2, 2], [0, 4], [3, 1]]
second = [[0, 4], [2, 2], [4, 0], [1, 3]]
result = calculate_moving_patterson_fst(first, second, window_size=2, verbose=False)
print(result["values"])
```

- **示例输出：**

```text
[0.5555555556, 0.6923076923]
```

## 22) 计算样本两两 dxy 矩阵

- **作用：** 计算样本两两之间的核苷酸分歧（dxy），输出对称方阵。
- **计算方法：** 调用 allel.pairwise_dxy 得到压缩矩阵，再展开为方阵。
- **API：** `calculate_pairwise_dxy(positions, genotype_allele_counts, *, start=None, stop=None, is_accessible=None, verbose=True)`
- **输入：** 位置、形如 (变异, 样本, 等位基因) 的计数，可选区域与可访问掩码。
- **输出：** 返回 samples 与 matrix。
- **示例代码：**

```python
from dnakit.variants import calculate_pairwise_dxy

positions = [10, 20, 30]
genotypes = [
    [[4, 0], [2, 2], [0, 4]],
    [[3, 1], [1, 3], [2, 2]],
    [[4, 0], [4, 0], [0, 4]],
]
result = calculate_pairwise_dxy(positions, genotypes, verbose=False)
print(result["samples"], result["matrix"][0])
```

- **示例输出：**

```text
3 [0.0, 0.0535714286, 0.119047619]
```

## 23) 计算 Patterson f2

- **作用：** 计算两群体的 Patterson f2（分支长度估计），给出每个变异的值与均值。
- **计算方法：** 调用 allel.patterson_f2；均值只计算有定义的变异。
- **API：** `calculate_patterson_f2(first_counts, second_counts, *, verbose=True)`
- **输入：** 两群体计数矩阵。
- **输出：** 返回 values 与 mean。
- **示例代码：**

```python
from dnakit.variants import calculate_patterson_f2

first = [[4, 0], [2, 2], [0, 4], [3, 1]]
second = [[0, 4], [2, 2], [4, 0], [1, 3]]
result = calculate_patterson_f2(first, second, verbose=False)
print(result["values"], result["mean"])
```

- **示例输出：**

```text
[1.0, -0.1666666667, 1.0, 0.125] 0.4895833333
```

## 24) 计算 Patterson f3 与 f3*

- **作用：** 计算三群体检验 f3(C; A, B) 的分子、杂合度之和，以及 f3* 比值。
- **计算方法：** 调用 allel.patterson_f3；f3* 为分子之和除以杂合度之和。
- **API：** `calculate_patterson_f3(target_counts, first_counts, second_counts, *, verbose=True)`
- **输入：** 目标群体 C 与群体 A、B 的计数矩阵。
- **输出：** 返回 numerator、heterozygosity、numerator_sum、heterozygosity_sum 与 f3_star。
- **示例代码：**

```python
from dnakit.variants import calculate_patterson_f3

target = [[4, 0], [2, 2], [1, 3], [3, 1]]
first = [[4, 0], [2, 2], [0, 4], [3, 1]]
second = [[0, 4], [2, 2], [4, 0], [1, 3]]
result = calculate_patterson_f3(target, first, second, verbose=False)
print(result["numerator_sum"], result["heterozygosity_sum"], result["f3_star"])
```

- **示例输出：**

```text
-0.3958333333 1.6666666667 -0.2375
```

## 25) 计算分块 jackknife 的 Patterson f3

- **作用：** 计算 f3 的块 jackknife 标准误与 z 值。
- **计算方法：** 调用 allel.blockwise_patterson_f3；无定义时为 None。
- **API：** `calculate_blockwise_patterson_f3(target_counts, first_counts, second_counts, *, block_length, normed=True, verbose=True)`
- **输入：** 三组计数矩阵、block_length，可选 normed。
- **输出：** 返回 f3、standard_error、z_score、block_values 与 jackknife_values。
- **示例代码：**

```python
from dnakit.variants import calculate_blockwise_patterson_f3

target = [[4, 0], [2, 2], [1, 3], [3, 1]]
first = [[4, 0], [2, 2], [0, 4], [3, 1]]
second = [[0, 4], [2, 2], [4, 0], [1, 3]]
result = calculate_blockwise_patterson_f3(target, first, second, block_length=2, verbose=False)
print(result["f3"], result["standard_error"], result["z_score"])
```

- **示例输出：**

```text
-0.2375 0.09375 -2.5333333333
```

## 26) 计算滑动变异窗口的 Patterson f3

- **作用：** 在连续变异窗口中计算 f3。
- **计算方法：** 调用 allel.moving_patterson_f3。
- **API：** `calculate_moving_patterson_f3(target_counts, first_counts, second_counts, *, window_size, start=0, stop=None, step=None, normed=True, verbose=True)`
- **输入：** 三组计数矩阵与窗口参数。
- **输出：** 返回 window_size、step、normed 与 values。
- **示例代码：**

```python
from dnakit.variants import calculate_moving_patterson_f3

target = [[4, 0], [2, 2], [1, 3], [3, 1]]
first = [[4, 0], [2, 2], [0, 4], [3, 1]]
second = [[0, 4], [2, 2], [4, 0], [1, 3]]
result = calculate_moving_patterson_f3(target, first, second, window_size=2, verbose=False)
print(result["values"])
```

- **示例输出：**

```text
[-0.125, -0.3125]
```

## 27) 给共享前缀的单倍型染色

- **作用：** 按前缀排序单倍型，并为每段共享前缀赋予同一整数。
- **计算方法：** 调用 allel.voight_painting；只接受无缺失的双等位（0/1）矩阵。
- **API：** `paint_shared_haplotype_prefixes(haplotypes, *, verbose=True)`
- **输入：** 形如 (变异, 单倍型) 的 0/1 矩阵。
- **输出：** 返回 painting 矩阵与 haplotype_order。
- **示例代码：**

```python
from dnakit.variants import paint_shared_haplotype_prefixes

haplotypes = [[0, 0, 1, 1], [0, 0, 1, 0], [1, 1, 0, 0], [0, 1, 0, 1]]
result = paint_shared_haplotype_prefixes(haplotypes, verbose=False)
print(result["haplotype_order"], result["painting"][0])
```

- **示例输出：**

```text
[0, 1, 3, 2] [1, 1, 2, 2]
```

## 28) 计算单倍型多样性

- **作用：** 计算单倍型多样性：1 减去单倍型频率平方和，并做 n/(n-1) 校正。
- **计算方法：** 调用 allel.haplotype_diversity。
- **API：** `calculate_haplotype_diversity(haplotypes, *, verbose=True)`
- **输入：** 形如 (变异, 单倍型) 的 0/1 矩阵，至少两条单倍型。
- **输出：** 返回 diversity、variants 与 haplotypes。
- **示例代码：**

```python
from dnakit.variants import calculate_haplotype_diversity

haplotypes = [[0, 0, 1, 1], [0, 0, 1, 0], [1, 1, 0, 0], [0, 1, 0, 1]]
result = calculate_haplotype_diversity(haplotypes, verbose=False)
print(result["haplotypes"], result["diversity"])
```

- **示例输出：**

```text
4 1.0
```

## 29) 计算滑动窗口单倍型多样性

- **作用：** 在连续变异窗口中计算单倍型多样性。
- **计算方法：** 调用 allel.moving_haplotype_diversity。
- **API：** `calculate_moving_haplotype_diversity(haplotypes, *, window_size, start=0, stop=None, step=None, verbose=True)`
- **输入：** 0/1 单倍型矩阵与窗口参数。
- **输出：** 返回 window_size、step 与 values。
- **示例代码：**

```python
from dnakit.variants import calculate_moving_haplotype_diversity

haplotypes = [[0, 0, 1, 1], [0, 0, 1, 0], [1, 1, 0, 0], [0, 1, 0, 1], [1, 0, 1, 0], [1, 1, 1, 0]]
result = calculate_moving_haplotype_diversity(haplotypes, window_size=3, verbose=False)
print(result["values"])
```

- **示例输出：**

```text
[0.8333333333, 0.8333333333]
```

## 30) 计算跨群体 nSL

- **作用：** 计算两组单倍型之间的跨群体 nSL 分数（未标准化）。
- **计算方法：** 调用 allel.xpnsl，并关闭线程以保证结果确定。
- **API：** `calculate_xpnsl(haplotypes_a, haplotypes_b, *, verbose=True)`
- **输入：** 两组形如 (变异, 单倍型) 的 0/1 矩阵，变异数一致。
- **输出：** 返回 variants 与 xpnsl 列表。
- **示例代码：**

```python
from dnakit.variants import calculate_xpnsl

haplotypes_a = [[0, 0, 1, 1], [0, 0, 1, 0], [1, 1, 0, 0], [0, 1, 0, 1], [1, 0, 1, 0], [1, 1, 1, 0]]
haplotypes_b = [[0, 1, 1, 0], [0, 0, 1, 1], [1, 1, 0, 0], [1, 0, 0, 1], [0, 1, 0, 1], [1, 1, 1, 0]]
result = calculate_xpnsl(haplotypes_a, haplotypes_b, verbose=False)
print(result["xpnsl"])
```

- **示例输出：**

```text
[0.4054651081, 0.2876820725, 0.0, 0.5596157879, 0.3364722366, 0.1335313926]
```

## 31) 按等位基因计数分箱标准化选择评分

- **作用：** 把 iHS、nSL 等评分在等位基因计数分箱内标准化。
- **计算方法：** 调用 allel.standardize_by_allele_count，关闭绘图诊断；bins 必须严格递增。
- **API：** `standardize_score_by_allele_count(score, allele_counts, *, n_bins=None, bins=None, verbose=True)`
- **输入：** 评分列表 score、每个位点的等位基因计数 allele_counts，可选 n_bins 或 bins。
- **输出：** 返回 values 与实际使用的 bins。
- **示例代码：**

```python
from dnakit.variants import standardize_score_by_allele_count

score = [0.1, -0.2, 0.3, 0.5, -0.4, 0.2, 0.0, 1.0]
counts = [1, 2, 3, 4, 5, 6, 7, 8]
result = standardize_score_by_allele_count(score, counts, n_bins=2, verbose=False)
print(result["bins"], [round(v, 4) for v in result["values"][:3]])
```

- **示例输出：**

```text
[1, 8] [-0.2163, -0.9581, 0.2781]
```

## 32) 计算两群体 Tajima D 差值的滑动窗口

- **作用：** 在连续变异窗口中计算两群体 Tajima D 之差。
- **计算方法：** 调用 allel.moving_delta_tajima_d；位点不足的窗口返回 None。
- **API：** `calculate_moving_delta_tajimas_d(first_counts, second_counts, *, window_size, start=0, stop=None, step=None, verbose=True)`
- **输入：** 两群体计数矩阵与窗口参数。
- **输出：** 返回 window_size、step 与 values。
- **示例代码：**

```python
from dnakit.variants import calculate_moving_delta_tajimas_d

first = [[4, 0], [3, 1], [2, 2], [1, 3], [4, 0], [0, 4], [3, 1], [2, 2], [1, 3], [4, 0], [2, 2], [3, 1]]
second = [[0, 4], [2, 2], [1, 3], [3, 1], [2, 2], [4, 0], [1, 3], [3, 1], [2, 2], [0, 4], [4, 0], [2, 2]]
result = calculate_moving_delta_tajimas_d(first, second, window_size=6, verbose=False)
print(result["values"])
```

- **示例输出：**

```text
[-1.0, 1.0]
```

## 33) 折叠双群体联合位点频谱

- **作用：** 把两群体的联合位点频谱折叠到次要等位基因一侧。
- **计算方法：** 调用 allel.fold_joint_sfs；每一维的类别数不得超过样本量加一。
- **API：** `fold_joint_site_frequency_spectrum(spectrum, *, first_sample_size, second_sample_size, verbose=True)`
- **输入：** 二维整数频谱与两组样本量。
- **输出：** 返回 folded 矩阵及两组样本量。
- **示例代码：**

```python
from dnakit.variants import fold_joint_site_frequency_spectrum

spectrum = [[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11]]
result = fold_joint_site_frequency_spectrum(spectrum, first_sample_size=2, second_sample_size=3, verbose=False)
print(result["folded"])
```

- **示例输出：**

```text
[[14, 14]]
```

## 34) 缩放位点频谱

- **作用：** 把频谱第 i 类乘以 i，用于 i·ξi 形式的统计量。
- **计算方法：** 调用 allel.scale_sfs。
- **API：** `scale_site_frequency_spectrum(spectrum, *, verbose=True)`
- **输入：** 一维整数频谱。
- **输出：** 返回 classes 与 scaled。
- **示例代码：**

```python
from dnakit.variants import scale_site_frequency_spectrum

result = scale_site_frequency_spectrum([5, 3, 2, 1], verbose=False)
print(result["scaled"])
```

- **示例输出：**

```text
[0, 3, 4, 3]
```

## 35) 缩放折叠位点频谱

- **作用：** 按 k(n−k)/n 缩放折叠频谱第 k 类。
- **计算方法：** 调用 allel.scale_sfs_folded；类别数不得超过样本量加一。
- **API：** `scale_folded_site_frequency_spectrum(spectrum, *, sample_size, verbose=True)`
- **输入：** 一维折叠频谱与染色体数 sample_size。
- **输出：** 返回 sample_size 与 scaled。
- **示例代码：**

```python
from dnakit.variants import scale_folded_site_frequency_spectrum

result = scale_folded_site_frequency_spectrum([5, 3, 2], sample_size=6, verbose=False)
print(result["scaled"])
```

- **示例输出：**

```text
[0.0, 2.5, 2.6666666667]
```

## 36) 缩放联合位点频谱

- **作用：** 把联合频谱每个格子 (i, j) 乘以 i·j。
- **计算方法：** 调用 allel.scale_joint_sfs。
- **API：** `scale_joint_site_frequency_spectrum(spectrum, *, verbose=True)`
- **输入：** 二维整数频谱。
- **输出：** 返回 shape 与 scaled。
- **示例代码：**

```python
from dnakit.variants import scale_joint_site_frequency_spectrum

result = scale_joint_site_frequency_spectrum([[0, 1, 2], [3, 4, 5]], verbose=False)
print(result["scaled"])
```

- **示例输出：**

```text
[[0, 0, 0], [0, 4, 10]]
```

## 37) 缩放折叠联合位点频谱

- **作用：** 按 i·j·(n1−i)·(n2−j) 缩放折叠联合频谱的每个格子。
- **计算方法：** 调用 allel.scale_joint_sfs_folded。
- **API：** `scale_folded_joint_site_frequency_spectrum(spectrum, *, first_sample_size, second_sample_size, verbose=True)`
- **输入：** 二维频谱与两组样本量。
- **输出：** 返回样本量与 scaled。
- **示例代码：**

```python
from dnakit.variants import scale_folded_joint_site_frequency_spectrum

result = scale_folded_joint_site_frequency_spectrum(
    [[0, 1, 2], [3, 4, 5]], first_sample_size=2, second_sample_size=3, verbose=False
)
print(result["scaled"])
```

- **示例输出：**

```text
[[0, 0, 0], [0, 8, 10]]
```

## 38) 按亲本遗传给后代单倍型染色

- **作用：** 对后代每个等位基因标注它来自亲本哪条单倍型，或属于非亲本等位基因。
- **计算方法：** 调用 allel.paint_transmission；codes 字段给出编码含义。
- **API：** `paint_parental_transmission(parent_haplotypes, progeny_haplotypes, *, verbose=True)`
- **输入：** 亲本两条单倍型 (变异, 2) 与后代单倍型 (变异, 后代数)；-1 表示缺失。
- **输出：** 返回 painting 矩阵与 codes 编码表。
- **示例代码：**

```python
from dnakit.variants import paint_parental_transmission

parent = [[0, 1], [1, 1]]
progeny = [[0, 1, 0], [1, 0, 1]]
result = paint_parental_transmission(parent, progeny, verbose=False)
print(result["painting"])
```

- **示例输出：**

```text
[[1, 2, 1], [4, 5, 4]]
```

## 39) 读取本地 VCF 头部摘要

- **作用：** 读取本地 VCF 头部中的样本名、INFO 与 FORMAT 字段、过滤器和 contig 长度。
- **计算方法：** 调用 allel.read_vcf_headers，并解析 ##contig 行。
- **API：** `read_vcf_header_summary(path, *, verbose=True)`
- **输入：** 本地 VCF 路径 path。
- **输出：** 返回 samples、info_fields、format_fields、filters 与 contigs。
- **示例代码：**

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit.variants import read_vcf_header_summary

ROWS = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=1000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "chr1\t10\trs1\tA\tG\t50\tPASS\tDP=10\tGT\t0|1\t0/0\t1|1",
    "chr1\t20\trs2\tC\tT\t50\tPASS\tDP=12\tGT\t0/1\t1/1\t0/0",
    "chr1\t35\trs3\tG\tA\t50\tPASS\tDP=8\tGT\t1/1\t0/1\t0|0",
]
with TemporaryDirectory() as folder:
    source = Path(folder) / "small.vcf"
    source.write_text("\n".join(ROWS) + "\n", encoding="utf-8")
    result = read_vcf_header_summary(str(source), verbose=False)
print(result["samples"], result["info_fields"], result["contigs"])
```

- **示例输出：**

```text
['S1', 'S2', 'S3'] ['DP'] [{'id': 'chr1', 'length': 1000}]
```

## 40) 将本地 VCF 转为 npz 数组归档

- **作用：** 把本地 VCF 读成 NumPy 数组并保存为新的压缩 .npz 文件。
- **计算方法：** 调用 allel.vcf_to_npz（overwrite=False）；输出路径必须是新的 .npz 文件。
- **API：** `vcf_to_npz_archive(source, output, *, fields=None, verbose=True)`
- **输入：** 源 VCF、新的 .npz 路径、可选字段列表 fields。
- **输出：** 返回 output 路径与每个数组的形状 arrays。
- **示例代码：**

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit.variants import vcf_to_npz_archive

ROWS = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=1000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "chr1\t10\trs1\tA\tG\t50\tPASS\tDP=10\tGT\t0|1\t0/0\t1|1",
    "chr1\t20\trs2\tC\tT\t50\tPASS\tDP=12\tGT\t0/1\t1/1\t0/0",
    "chr1\t35\trs3\tG\tA\t50\tPASS\tDP=8\tGT\t1/1\t0/1\t0|0",
]
with TemporaryDirectory() as folder:
    source = Path(folder) / "small.vcf"
    source.write_text("\n".join(ROWS) + "\n", encoding="utf-8")
    result = vcf_to_npz_archive(
        str(source), str(Path(folder) / "small.npz"),
        fields=["variants/POS", "calldata/GT"], verbose=False,
    )
print(result["arrays"])
```

- **示例输出：**

```text
{'calldata/GT': [3, 3, 2], 'variants/POS': [3]}
```

## 41) 将本地 VCF 转为 HDF5 归档

- **作用：** 把本地 VCF 读成数组并写入新的 HDF5 文件。
- **计算方法：** 调用 allel.vcf_to_hdf5（需要 h5py）；完成后读取各数据集的形状。
- **API：** `vcf_to_hdf5_archive(source, output, *, fields=None, verbose=True)`
- **输入：** 源 VCF、新输出路径与可选字段列表 fields。
- **输出：** 返回 output 路径与 datasets（数据集路径到形状）。
- **示例代码：**

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit.variants import vcf_to_hdf5_archive

ROWS = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=1000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "chr1\t10\trs1\tA\tG\t50\tPASS\tDP=10\tGT\t0|1\t0/0\t1|1",
    "chr1\t20\trs2\tC\tT\t50\tPASS\tDP=12\tGT\t0/1\t1/1\t0/0",
    "chr1\t35\trs3\tG\tA\t50\tPASS\tDP=8\tGT\t1/1\t0/1\t0|0",
]
with TemporaryDirectory() as folder:
    source = Path(folder) / "small.vcf"
    source.write_text("\n".join(ROWS) + "\n", encoding="utf-8")
    result = vcf_to_hdf5_archive(
        str(source), str(Path(folder) / "small.h5"),
        fields=["variants/POS", "calldata/GT"], verbose=False,
    )
print(result["datasets"])
```

- **示例输出：**

```text
{'calldata/GT': [3, 3, 2], 'variants/POS': [3]}
```

## 未包装的函数

以下函数不在本页的包装中，原因写入账本的 `reason` 字段。

- **重复或已覆盖：** 共 18 个。varcode 的 compare_samples 与 query_in_cis，scikit-allel 的 average_patterson_fst、average_patterson_f3、fold_sfs（population 的 dadi 版本覆盖）与 write_fasta（dnakit.io.write 覆盖）。
- **同族变体：** 共 32 个。DataFrame、CSV 或 recarray 形式的读取函数，以及 allel.write_vcf；对应的记录级函数已由 dnakit.variants.read_vcf 或 rewrite_vcf_records 覆盖。
- **需要 pyensembl 转录本、效应或 ExACTO 文件：** 共 7 个 varcode 函数；allel.vcf_to_zarr 需要 zarr 2（当前环境为 zarr 3.4.1），共 1 个。
- **内部辅助：** 共 712 个，包括注册表与全局配置、分组与排序键、CSV 元数据、窗口与索引构造、通用滑动统计、HDF5 存储工具，以及只产生中间对象的函数。
