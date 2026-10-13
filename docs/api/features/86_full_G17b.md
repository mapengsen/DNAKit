# G17b 功能：selene-sdk、genomic-benchmarks 与 kipoi

本页记录 `predictions` 区域第 G17b 组的包装函数。每个函数位于对应模块中，可直接导入，例如 `from dnakit.predictions import decode_one_hot_sequence`。

- 本地函数：全部 13 个包装函数只读取本地输入，例如 one-hot 列表、DNA 字符串、本地 FASTA 与 VCF、本地分数表，不需要网络。
- kipoi 在本组没有可离线运行的 DNA 级函数，因此 `kipoi_functions.py` 不导出任何函数。kipoi 的列表与下载函数只记录在台账中；现有的 `list_kipoi_models` 与 `describe_kipoi_model` 在单元测试中使用替身模块，不会克隆模型库。
- 示例中的输入文件以 `示例文件` 标签给出；其余示例直接运行。

## 已有包装与未包装的函数

- 已由现有模块包装：selene-sdk 的 `sequence_to_encoding`（由 `Genome.sequence_to_encoding` 调用，即 `encode_one_hot_selene`）；genomic-benchmarks 的 `info`、`is_downloaded` 与 `list_datasets`（即 `describe_genomic_benchmark`、`check_genomic_benchmark_cache` 与 `list_genomic_benchmark_datasets`）；kipoi 的 `list_models` 与 `get_model_descr`（即 `list_kipoi_models` 与 `describe_kipoi_model`）。
- 重复：同一函数的模块别名共 21 个，其中 selene 15 个、genomic-benchmarks 3 个、kipoi 3 个。selene 的 `reverse_complement_sequence` 与 `reverse_complement_with_grelu` 做同样的 DNA 字符串反向互补，因此标记为重复。
- selene-sdk 未包装的 75 个键中：绘图函数 10 个；命令行与训练入口 7 个；日志、配置、运行目录、模型推理与输出写入等 27 个（其中 `get_reverse_encoding` 的 2 个键因上游索引错误而不包装，对 one-hot 输入抛出 IndexError）；模型权重、评估指标与标签等与 DNA 无关的函数 13 个；已弃用类 2 个；内部辅助 3 个。
- genomic-benchmarks：数据集下载与云端缓存 14 个（需要网络）；缓存删除与 `labels_in_order` 3 个（有副作用或依赖文件系统顺序）。
- kipoi：命令行 88 个；需要 singularity 可执行文件 3 个；需要网络或克隆模型库 37 个；通用数组、字典与文件工具 52 个；内部类与字段 23 个；插件、环境、配置、文件与路径等辅助函数 86 个。完整状态见 G17b 台账。

## 1) `decode_one_hot_sequence`

- **作用：** 把 (L, 4) 的 one-hot 矩阵（列为 A、C、G、T）还原为 DNA 字符串。
- **处理方法：** 每一行必须是 one-hot（一个 1、其余为 0）或均匀的 0.25。合法行交给 selene-sdk 的 `sequence.encoding_to_sequence` 解码，均匀 0.25 行输出 `unk_base`；其他行直接报错，避免静默解码。
- **API：** `decode_one_hot_sequence(encoding[必须], unk_base[可选], verbose[可选])`
- **输入：** 行列表，每行 4 个有限数值；`unk_base` 取 A、C、G、T 或 N，默认 N。
- **输出：** 字典，含 `length`、`sequence`、`unk_base` 与 `provenance`。

**示例代码：**

```python
from dnakit.predictions import decode_one_hot_sequence

rows = [[1, 0, 0, 0], [0, 0, 0, 1], [0.25, 0.25, 0.25, 0.25], [0, 1, 0, 0]]
result = decode_one_hot_sequence(rows, verbose=False)
print(result["sequence"], result["length"])
```

**示例输出：**

```text
ATNC 4
```

## 2) `reverse_complement_one_hot`

- **作用：** 对 (L, 4) 的 one-hot 矩阵做反向互补。
- **处理方法：** 调用 selene-sdk 的 `model_predict.get_reverse_complement_encoding`：先按 A-T、C-G 交换列，再把行顺序反转。N 行（均匀 0.25）保持不变。
- **API：** `reverse_complement_one_hot(encoding[必须], verbose[可选])`
- **输入：** 同 `decode_one_hot_sequence`，行必须是 one-hot 或均匀 0.25。
- **输出：** 字典，含 `length`、`encoding`（反向互补后的矩阵）与 `provenance`。

**示例代码：**

```python
from dnakit.predictions import (
    decode_one_hot_sequence,
    reverse_complement_one_hot,
)

rows = [[1, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]]
rc = reverse_complement_one_hot(rows, verbose=False)["encoding"]
print(decode_one_hot_sequence(rc, verbose=False)["sequence"])
```

**示例输出：**

```text
CGTT
```

## 3) `reverse_complement_allele`

- **作用：** 对一个变异等位基因字符串做反向互补。
- **处理方法：** 调用 selene-sdk 的 `model_predict.get_reverse_complement`。`-` 与 `*`（缺失或缺失型等位基因）单独出现时返回 `*`；其他情况只允许 A、C、G、T、N，并在输入中混入 `-` 或 `*` 时报错。
- **API：** `reverse_complement_allele(allele[必须], verbose[可选])`
- **输入：** 等位基因字符串，长度 1 至 256；大小写不敏感。
- **输出：** 字典，含 `allele`（大写）、`reverse_complement` 与 `provenance`。

**示例代码：**

```python
from dnakit.predictions import reverse_complement_allele

for allele in ["AC", "GATTACA", "-"]:
    print(allele, reverse_complement_allele(allele, verbose=False)["reverse_complement"])
```

**示例输出：**

```text
AC GT
GATTACA TGTAATC
- *
```

## 4) `apply_one_hot_mutations`

- **作用：** 把 one-hot 矩阵中指定的 0 起始位置改为新的碱基。
- **处理方法：** 调用 selene-sdk 的 `model_predict.mutate_sequence`，把该位置整行清零后置为新碱基的 one-hot。同一位置不能重复出现；位置超出范围或碱基不是 A、C、G、T 时报错。
- **API：** `apply_one_hot_mutations(encoding[必须], mutations[必须], verbose[可选])`
- **输入：** one-hot 矩阵；`mutations` 为 `[位置, 碱基]` 的列表，最多 1000 项。
- **输出：** 字典，含 `length`、`mutation_count`、按位置排序的 `mutations`、突变后的 `encoding` 与 `provenance`。

**示例代码：**

```python
from dnakit.predictions import (
    apply_one_hot_mutations,
    decode_one_hot_sequence,
)

rows = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
result = apply_one_hot_mutations(rows, [[0, "T"], [2, "A"]], verbose=False)
print(result["mutation_count"], decode_one_hot_sequence(result["encoding"], verbose=False)["sequence"])
```

**示例输出：**

```text
2 TCAT
```

## 5) `list_in_silico_mutants`

- **作用：** 列出 DNA 字符串指定窗口内，恰好改变 n 个碱基的全部突变体（in silico mutagenesis）。
- **处理方法：** 调用 selene-sdk 的 `model_predict.in_silico_mutagenesis_sequences`。ACGT 位置有 3 种替换，N 位置有 4 种替换。生成前先按组合公式计算突变数，超过 `max_mutants` 时报错，不会先生成再截断。位置是整条序列上的 0 起始坐标。
- **API：** `list_in_silico_mutants(sequence[必须], mutate_n_bases[可选], start_position[可选], end_position[可选], max_mutants[可选], verbose[可选])`
- **输入：** 大写 DNA 字符串，长度至多 10000；`mutate_n_bases` 为 1 至 4；`max_mutants` 默认 100000。
- **输出：** 字典，含 `count`、`mutants`（每个突变体为 `[位置, 碱基]` 列表）与 `provenance`。顺序为先按位置组合、再按碱基顺序。

**示例代码：**

```python
from dnakit.predictions import list_in_silico_mutants

result = list_in_silico_mutants("ACG", verbose=False)
print(result["count"], result["mutants"][0], result["mutants"][-1])
```

**示例输出：**

```text
9 [[0, 'C']] [[2, 'T']]
```

## 6) `read_selene_vcf`

- **作用：** 读取 VCF 记录，并对照本地 FASTA 检查染色体、坐标范围与 REF 碱基。
- **处理方法：** 用 selene-sdk 的 `Genome` 打开 FASTA（pyfaidx 会在 FASTA 旁写入 `.fai` 索引），再调用 `model_predict.read_vcf_file`。上游不检查 REF 和坐标，因此包装器逐条补充 `in_bounds` 与 `ref_matches_fasta`；VCF 中出现 FASTA 没有的染色体时报错。
- **API：** `read_selene_vcf(vcf_path[必须], fasta_path[必须], max_variants[可选], verbose[可选])`
- **输入：** VCF 文件（前五列须为 `#CHROM POS ID REF ALT`，文件不超过 50 MB）；FASTA 文件；`max_variants` 默认 100000。
- **输出：** 字典，含 `variant_count` 与 `variants`。每条记录包含 `chromosome`、`position`（VCF 的 1 起始 POS）、`id`、`ref`、`alt`、`strand`、`in_bounds` 与 `ref_matches_fasta`（无法比较时为 null）。

**示例文件：** `ref.fa`

```text
>chr1
ACGTACGTAC
>chr2
GGGGTTTTAA
```

**示例文件：** `variants.vcf`

```text
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	3	rs1	G	T	.	.	.
chr1	4	.	A	G	.	.	.
chr2	5	.	T	C	.	.	.
```

**示例代码：**

```python
from dnakit.predictions import read_selene_vcf

result = read_selene_vcf("variants.vcf", "ref.fa", verbose=False)
for v in result["variants"]:
    print(v["chromosome"], v["position"], v["ref"], v["alt"], v["ref_matches_fasta"])
```

**示例输出：**

```text
chr1 3 G T True
chr1 4 A G False
chr2 5 T C True
```

## 7) `unpack_bit_packed_sequence`

- **作用：** 还原 selene-sdk 的位压缩 one-hot 矩阵（每个字节存 8 个位置，每个通道一列）。
- **处理方法：** 调用 selene-sdk 的 `samplers.dataloader.unpackbits_sequence`。四个通道全为 1 的行是 N，还原为均匀 0.25；其他行保持 0 或 1。输入行数必须至少为 `ceil(sequence_length / 8)`。
- **API：** `unpack_bit_packed_sequence(packed_rows[必须], sequence_length[必须], verbose[可选])`
- **输入：** 行列表，每行 4 个 0 至 255 的整数；`sequence_length` 为正整数。
- **输出：** 字典，含 `length`、`encoding`（前 L 个位置）与 `provenance`。

**示例代码：**

```python
from dnakit.predictions import (
    decode_one_hot_sequence,
    unpack_bit_packed_sequence,
)

result = unpack_bit_packed_sequence([[136, 72, 40, 24]], 5, verbose=False)
print(decode_one_hot_sequence(result["encoding"], verbose=False)["sequence"])
```

**示例输出：**

```text
ACGTN
```

## 8) `interval_length_probabilities`

- **作用：** 计算基因组区间按长度加权的抽样概率。
- **处理方法：** 调用 selene-sdk 的 `utils.get_indices_and_probabilities`，权重为区间长度除以所选区间长度之和。长度为 0 的区间权重为 0，会被去掉并重新计算；上游输出被去掉的区间由 `dropped_indices` 列出。
- **API：** `interval_length_probabilities(interval_lengths[必须], indices[可选], verbose[可选])`
- **输入：** 非负整数的区间长度列表（至多 1000000 项）；`indices` 为不重复的区间下标，缺省为全部；所选区间中至少一个长度为正。
- **输出：** 字典，含 `interval_count`、`indices`、`probabilities`、`dropped_indices` 与 `provenance`。

**示例代码：**

```python
from dnakit.predictions import interval_length_probabilities

result = interval_length_probabilities([100, 300, 0, 600], verbose=False)
print(result["indices"], [round(p, 3) for p in result["probabilities"]], result["dropped_indices"])
```

**示例输出：**

```text
[0, 1, 3] [0.1, 0.3, 0.6] [2]
```

## 9) `rescale_attribution_scores`

- **作用：** 对 (L, 4) 的归因分数矩阵做逐碱基与逐位置的缩放。
- **处理方法：** 调用 selene-sdk 的 `interpret.rescale_score_matrix`。`max_effect` 在碱基方向减去每个位置的最小值，位置方向按总和归一化后乘以该位置的极差；`probability` 在碱基方向是恒等变换。若某个位置的除数为 0，结果不是有限值，包装器报错而不输出 null。
- **API：** `rescale_attribution_scores(score_matrix[必须], base_scaling[可选], position_scaling[可选], verbose[可选])`
- **输入：** 行列表，每行 4 个有限数值（行为位置，列为 A、C、G、T）；两个缩放参数取 identity、probability 或 max_effect，默认 identity。
- **输出：** 字典，含 `length`、两个缩放参数、`scores` 与 `provenance`。

**示例代码：**

```python
from dnakit.predictions import rescale_attribution_scores

scores = [[1, 2, 3, 4], [0, 1, 0, 3]]
result = rescale_attribution_scores(
    scores, base_scaling="max_effect", position_scaling="max_effect", verbose=False
)
print([[round(v, 2) for v in row] for row in result["scores"]])
```

**示例输出：**

```text
[[0.0, 0.5, 1.0, 1.5], [0.0, 0.75, 0.0, 2.25]]
```

## 10) `load_variant_diff_scores`

- **作用：** 读取变异与特征的分数表（selene-sdk 变异效应分析的输出格式）。
- **处理方法：** 调用 selene-sdk 的 `interpret.load_variant_abs_diff_scores`。前五列为变异（染色体、位置、名称、REF、ALT），其余列为特征。位置必须是整数，否则报错；没有特征列时报错。
- **API：** `load_variant_diff_scores(scores_path[必须], max_rows[可选], verbose[可选])`
- **输入：** 制表符分隔的文本文件，不超过 100 MB；`max_rows` 默认 100000。
- **输出：** 字典，含 `variant_count`、`feature_count`、`features`、`variants`（每项为 `chromosome`、`position`、`name`、`ref`、`alt`）、`scores` 与 `provenance`。

**示例文件：** `scores.tsv`

```text
chrom	pos	name	ref	alt	f1	f2
chr1	10	rs1	A	G	0.5	0.25
chr2	20	rs2	C	T	0.1	0.0
```

**示例代码：**

```python
from dnakit.predictions import load_variant_diff_scores

result = load_variant_diff_scores("scores.tsv", verbose=False)
print(result["variant_count"], result["features"], result["scores"][0])
```

**示例输出：**

```text
2 ['f1', 'f2'] [0.5, 0.25]
```

## 11) `order_variant_labels`

- **作用：** 按染色体、位置、REF 与 ALT 给变异排序，并返回每个变异在输入中的原始下标。
- **处理方法：** 调用 selene-sdk 的 `interpret.ordered_variants_and_indices`，染色体顺序由 `sort_standard_chrs` 决定（chr1 至 chr22、chrX、chrY、chrM，其余在最后）。包装器先按同一规则与名称、输入下标对输入排序，因此并列项的顺序与输入顺序无关。上游会丢弃变异名称，但 `index` 保留了原始位置。
- **API：** `order_variant_labels(labels[必须], verbose[可选])`
- **输入：** 列表，每项为 `[chrom, pos, name, ref, alt]`；`chrom` 必须以 `chr` 开头，`pos` 为非负整数，最多 200000 项。
- **输出：** 字典，含 `count`、`order`（排序后各项的原始下标）、`variants`（每项为 `index`、`chromosome`、`position`、`ref`、`alt`）与 `provenance`。

**示例代码：**

```python
from dnakit.predictions import order_variant_labels

labels = [["chr2", 30, "r1", "A", "C"], ["chr1", 20, "r2", "G", "T"], ["chr1", 10, "r3", "C", "A"]]
result = order_variant_labels(labels, verbose=False)
print(result["order"], [v["chromosome"] + ":" + str(v["position"]) for v in result["variants"]])
```

**示例输出：**

```text
[2, 1, 0] ['chr1:10', 'chr1:20', 'chr2:30']
```

## 12) `sort_standard_chromosomes`

- **作用：** 按 selene-sdk 的染色体顺序排列一组 `chr` 名称。
- **处理方法：** 用 selene-sdk 的 `interpret.sort_standard_chrs` 作为排序键：chr1 至 chr22 按数字，chrX 为 23，chrY 为 24，chrM 为 25，其余为 26；同一等级按名称排序。名称必须以 `chr` 开头，且不能重复。
- **API：** `sort_standard_chromosomes(chromosomes[必须], verbose[可选])`
- **输入：** 不重复的 `chr` 名称列表，最多 10000 项。
- **输出：** 字典，含 `count` 与 `chromosomes`（每项为 `chromosome` 与 `rank`）。

**示例代码：**

```python
from dnakit.predictions import sort_standard_chromosomes

result = sort_standard_chromosomes(["chr10", "chr2", "chrX", "chr1", "chrM", "chrUn_a"], verbose=False)
print([item["chromosome"] for item in result["chromosomes"]])
```

**示例输出：**

```text
['chr1', 'chr2', 'chr10', 'chrX', 'chrM', 'chrUn_a']
```

## 13) `locate_sequences_in_reference`

- **作用：** 在参考染色体中查找 FASTA 中每条序列的精确匹配位置（正链或反向互补链）。
- **处理方法：** 调用 genomic-benchmarks 的 `seq2loc.fasta2loc`，只做离线比较，不需要网络。上游每个序列名只记录一个位置，即从左到右扫描时最后找到的匹配。包装器另外统计每条序列的全部精确匹配数（`occurrences`），大于 1 时标记 `ambiguous`。包装器还拒绝两种会互相覆盖的输入：两条序列相同，或互为反向互补。回文序列在上游中报告为 `-` 链。
- **API：** `locate_sequences_in_reference(fasta_path[必须], reference[必须], max_reference_length[可选], max_sequences[可选], max_query_length[可选], verbose[可选])`
- **输入：** FASTA 文件，序列名不重复、只含大写 A、C、G、T、N；`reference` 为染色体名到 DNA 字符串的映射，总长默认上限 2000000 bp；`max_sequences` 默认 100，`max_query_length` 默认 100000。
- **输出：** 字典，含 `query_count`、`contig_count`、`matched_count`、`matches`（每项含 `sequence`、`chromosome`、`start`、`end`、`strand`、`occurrences`、`ambiguous`）、`unmatched` 与 `messages`。坐标为 0 起始、左闭右开；`-` 表示查询序列的反向互补出现在该区间。

**示例文件：** `queries.fa`

```text
>unique
GGGG
>repeat
GTAC
>absent
AAAAA
```

**示例代码：**

```python
from dnakit.predictions import locate_sequences_in_reference

reference = {"chr1": "ACGTACGTAC", "chr2": "GGGGTTTTAA"}
result = locate_sequences_in_reference("queries.fa", reference, verbose=False)
for match in result["matches"]:
    print(match["sequence"], match["chromosome"], match["start"], match["end"],
          match["strand"], match["occurrences"])
print(result["unmatched"])
```

**示例输出：**

```text
repeat chr1 6 10 - 2
unique chr2 0 4 + 1
['absent']
```
